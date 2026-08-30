import os
import sys
import csv
import json
import glob
import math
import argparse
from datetime import datetime

# =============================================================================
# MINIMAL BOOTSTRAP -- locates the shared phase1_common module.
# See Phase 1/_common/phase1_common.py for why this logic is centralized.
# =============================================================================
def _bootstrap_find_research_root(script_file):
    current = os.path.dirname(os.path.abspath(script_file))
    while os.path.basename(current) != "Research":
        parent = os.path.dirname(current)
        if parent == current:
            print(f"\n[FATAL ERROR] Could not locate a 'Research' anchor folder above: {script_file}")
            sys.exit(1)
        current = parent
    return current

_PROJECT_ROOT = _bootstrap_find_research_root(__file__)
_COMMON_DIR = os.path.join(_PROJECT_ROOT, "Phase 1", "_common")
if _COMMON_DIR not in sys.path:
    sys.path.insert(0, _COMMON_DIR)
_STEPD_DIR = os.path.dirname(os.path.abspath(__file__))
if _STEPD_DIR not in sys.path:
    sys.path.insert(0, _STEPD_DIR)

import phase1_common as common

# =============================================================================
# PHASE 1De: CONFORMATIONAL B-CELL CORROBORATION WITH SEMA-3D (LOCAL)
#
# Sec. I.D: B-cell candidates from BepiPred-2.0 are cross-referenced against
# SEMA 2.0 conformational-patch predictions on the native (unassembled) antigen
# fold. >=50% residue overlap with a patch -> "structurally corroborated",
# PRIORITIZED. No overlap -> retained as lower-priority linear-only, FLAGGED.
# This step never excludes a candidate.
#
# ---------------------------------------------------------------------------
# WHAT CHANGED, AND WHY
#
# 1. THE HOSTED API IS GONE FROM THIS STEP. sema.airi.net's prediction
#    endpoints sit behind a ServicePipe anti-bot WAF that returns the SPA shell
#    to every non-browser client (deviation #20), so this step produced nothing
#    and all seven targets stayed UNSCREENED. SEMA is open source, so the model
#    now runs LOCALLY (see sema3d_local.py). The blocker disappears rather than
#    being worked around.
#
# 2. SEMA-3D, NOT SEMA-1D. The previous implementation posted esm_switch="1d"
#    -- the sequence-only ESM-2 model -- against the raw FASTA. Sec. I.D
#    specifies the structure-based method: conformational epitope propensity
#    "directly from AlphaFold-predicted or experimentally solved structures".
#    Predicting *conformational* epitopes from sequence alone was the wrong
#    tool for the claim being made. This step now consumes a native fold per
#    antigen (see Phase1De_foldPrep.py).
#
# 3. THE THRESHOLD WAS MEANINGLESS. The old code treated ">= 0.5" as "in an
#    epitope patch", commenting that "SEMA scores are typically probability-
#    like in [0,1]". They are not. SEMA outputs the LOG-SCALED EXPECTED NUMBER
#    OF ANTIBODY CONTACTS -- an unbounded regression score. In SEMA 2.0's own
#    training data the target is ln(1 + raw_contacts): ln(2)=0.6931 is one
#    contact, ln(4)=1.3863 is three. A 0.5 cut has no interpretation at all.
#    EPITOPE_SCORE_CUT below is derived from SEMA's own labelled test set --
#    see its comment.
#
# 4. IT NO LONGER WRITES BACK INTO THE PRODUCTION POOL BY DEFAULT. See the
#    guardrail note on run_step1de_conformational().
# ---------------------------------------------------------------------------

REQUIRED_TARGETS = ["Mpox_L1R", "Mpox_B5R", "Mpox_A35R",
                    "HIV_gp120", "HIV_gp41", "HIV_p24", "HIV_p17"]

# Sec. I.D's own rule: a linear candidate counts as structurally corroborated
# when at least half of its residues fall inside a conformational patch.
OVERLAP_THRESHOLD = 0.50

# Per-residue score above which a residue is "in a conformational patch".
#
# DERIVED, NOT ASSUMED. Calibrated on SEMA 2.0's own bundled test set
# (external_tools/SEMAi/epitopes_prediction/data/sema_2.0/test_set.csv, 11,302
# labelled residues) using the same Youden criterion the upstream
# SEMA-3D_inference notebook uses to pick its operating point
# (np.argmax(tpr - fpr)). That optimum lands on 1.386294 = ln(4), i.e. three or
# more expected antibody contacts, giving TPR 0.662 at FPR 0.072 against the
# curated contact_number_binary epitope label (AUC 0.779).
#
# CAVEAT TO STATE IN THE WRITE-UP: this calibrates the cut in the model's
# TARGET space (ground-truth contact number vs. the binary epitope label), not
# on this ensemble's predictions over that test set, which would additionally
# require SEMA's processed PDB archive. The model regresses directly onto that
# target, so the cut transfers, but it is a calibration on the label definition
# rather than on measured model error. Raw per-residue scores are written
# alongside every call so the cut can be moved without re-running anything.
EPITOPE_SCORE_CUT = math.log(4.0)   # 1.386294

MODELS_DIR = os.path.join(_PROJECT_ROOT, "external_tools", "SEMAi", "models")
HF_CACHE = os.path.join(_PROJECT_ROOT, "external_tools", "SEMAi", "hf_cache")


def _reference_sequence(target):
    """
    The MATURE SUBUNIT sequence for this target -- what the fold represents and
    what epitope offsets are measured against.

    Phase 1A stores POLYPROTEINS for the four HIV targets (gp120/gp41 share one
    856-aa Env record; p17/p24 share one 1437-aa Gag-Pol record), so returning
    the raw Var_01 record here would give gp120 and gp41 the same reference,
    the same fold and the same SEMA scores. HIV_SUBUNIT_RANGES slices each to
    its own mature chain; Mpox targets are already mature proteins.
    """
    folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1A")
    if not os.path.isdir(folder):
        return None
    matches = sorted(f for f in os.listdir(folder) if f.startswith(f"{target}_Var_01_"))
    if not matches:
        return None
    with open(os.path.join(folder, matches[0])) as f:
        seq = "".join(l.strip() for l in f if not l.startswith(">")).upper()
    info = common.HIV_SUBUNIT_RANGES.get(target)
    if info is None:
        return seq
    return seq[info["start"] - 1:info["end"]]


def _find_fold(target):
    """
    Returns the path to this target's native structure, or None.
    Prefers .pdb (HADDOCK/foldseek-friendly, and what Phase 3A also emits);
    AlphaFold Server hands out .cif, so both are accepted and .cif is converted
    on the fly by _ensure_pdb().
    """
    folds = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D",
                         "Phase1De", "Folds", target)
    if not os.path.isdir(folds):
        return None
    pdbs = sorted(glob.glob(os.path.join(folds, "*.pdb")))
    if pdbs:
        return pdbs[0]
    cifs = sorted(glob.glob(os.path.join(folds, "*.cif")))
    return cifs[0] if cifs else None


def _ensure_pdb(path):
    """
    foldseek reads mmCIF, but the rest of this project standardised on legacy
    PDB (Phase 3A converts for exactly this reason), and chain-naming differs
    between the two in ways that would silently change which chain is scored.
    Converting up front keeps one code path.
    """
    if path.lower().endswith(".pdb"):
        return path
    out = os.path.splitext(path)[0] + "_converted.pdb"
    if os.path.isfile(out):
        return out
    from Bio.PDB import MMCIFParser, PDBIO
    structure = MMCIFParser(QUIET=True).get_structure("model", path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(out)
    return out



def _model_plddt(model_path):
    """
    Per-residue pLDDT for the model backing this target, or None if it cannot
    be read. Reuses Phase1De_foldIngest's parser so there is one definition of
    how pLDDT comes out of an AlphaFold file.
    """
    try:
        from Phase1De_foldIngest import _per_residue_plddt
    except Exception:
        return None
    cif = model_path
    if cif.lower().endswith(".pdb"):
        stem = cif[:-len("_converted.pdb")] if cif.endswith("_converted.pdb") else cif[:-4]
        for cand in (stem + ".cif", cif):
            if os.path.isfile(cand) and cand.lower().endswith(".cif"):
                cif = cand
                break
    if not cif.lower().endswith(".cif") or not os.path.isfile(cif):
        return None
    try:
        vals = _per_residue_plddt(cif)
        return vals or None
    except Exception:
        return None


def _map_scores_to_reference(model_seq, model_scores, ref_seq):
    """
    Returns a list of per-residue scores indexed to ref_seq, using None where
    the structure does not cover a reference residue.

    WHY AN ALIGNMENT AND NOT A DIRECT INDEX. The scores come back indexed to
    the residues foldseek actually saw, i.e. the MODELLED residues. Any
    unmodelled loop, cleaved signal peptide or terminal disorder shifts every
    subsequent offset. Indexing epitopes straight into that array would score
    the wrong window and never announce it.
    """
    from Bio import Align
    aligner = Align.PairwiseAligner(
        mode="global", open_gap_score=-11, extend_gap_score=-1,
        substitution_matrix=Align.substitution_matrices.load("BLOSUM62"))
    aln = aligner.align(ref_seq, model_seq)[0]
    mapped = [None] * len(ref_seq)
    for (rs, re_), (ms, me) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(re_ - rs):
            mapped[rs + k] = float(model_scores[ms + k])
    covered = sum(1 for s in mapped if s is not None)
    return mapped, covered


def _score_epitope(mapped_scores, start, length):
    """
    Sec. I.D's rule: percentage of the epitope's residues sitting inside a
    conformational patch. Residues not covered by the structure are excluded
    from the denominator rather than counted as misses -- an unmodelled residue
    is unknown, not negative.
    """
    window = mapped_scores[start:start + length]
    scored = [s for s in window if s is not None]
    if not scored:
        return None, 0
    hits = sum(1 for s in scored if s >= EPITOPE_SCORE_CUT)
    return 100.0 * hits / len(scored), len(scored)


def run_step1de_conformational(write_back=False, targets=None, device="auto"):
    """
    GUARDRAIL -- READ THIS BEFORE SETTING write_back=True.

    The previous version of this step wrote its enriched CSV straight back into
    Step_Outputs/Phase1/Phase1D/Phase1Db/ with a fresh timestamp. Phase 1G
    selects its input via common.latest_file(), and SEMA_Corroborated feeds
    n_prioritized in 1G's objective -- so merely RUNNING this step with real
    data silently armed a different construct on the next 1G run.

    Vax_Final_6f34b53e is frozen: it has been through all of Phase 2 and has
    been docked. So the default is now a sensitivity output that nothing reads
    automatically. Phase1G_sensitivity.py picks it up deliberately and reports
    what WOULD change, without changing it.

    write_back=True restores the old behaviour and should only be used if a new
    construct is actually intended.
    """
    common.print_banner("PHASE 1De: SEMA-3D CONFORMATIONAL B-CELL CORROBORATION (LOCAL)")

    input_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1Db")
    output_dir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1De")
    sens_dir = os.path.join(output_dir, "Sensitivity")
    tool_runs_dir = os.path.join(output_dir, "_tool_runs")
    os.makedirs(sens_dir, exist_ok=True)
    os.makedirs(tool_runs_dir, exist_ok=True)

    latest_db = common.latest_file(input_folder, suffix=".csv")
    if latest_db is None:
        print(f"[ERROR] No Phase 1Db output found at: {input_folder}")
        return

    with open(latest_db) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        original_fields = reader.fieldnames

    bcell_by_target = {}
    for r in rows:
        if r.get("Type") == "B-cell":
            bcell_by_target.setdefault(r["Target"], set()).add(r["Peptide"])

    wanted = targets or REQUIRED_TARGETS
    folds = {t: _find_fold(t) for t in wanted}
    missing = [t for t in wanted if folds[t] is None]
    if missing:
        print(f"[ACTION AVAILABLE, NOT BLOCKING] No native structure for: {missing}")
        print( "[INFO] Run Phase1De_foldPrep.py and follow its HOW_TO_SUBMIT.txt. Those "
               "targets stay SEMA_Corroborated=UNSCREENED.")
    runnable = [t for t in wanted if folds[t] is not None]
    if not runnable:
        print("[INFO] Nothing to score yet -- no structures present. Stopping without writing.")
        return

    os.environ.setdefault("HF_HOME", HF_CACHE)
    import sema3d_local as sema

    dev = sema.pick_device(device)
    print(f"[INFO] SEMA-3D device: {dev} | epitope score cut: {EPITOPE_SCORE_CUT:.6f} "
          f"(= ln 4, >=3 expected antibody contacts)")

    # ---- resolve inputs, reuse cache, and decide what still needs scoring ----
    # The ensemble is run via predict_many(), which streams checkpoints one at a
    # time. Holding all five resident costs 2.43 GB each = 12.2 GB, measured, and
    # this is a 16 GB machine -- load_ensemble(n=5) would thrash or die mid-run.
    # Streaming keeps the peak at ONE model while still reading each checkpoint
    # from disk exactly once, and produces bit-identical scores.
    resolved = {}
    cached = {}
    for target in runnable:
        ref_seq = _reference_sequence(target)
        if ref_seq is None:
            print(f"[WARNING] No Phase 1A reference for {target} -- skipping.")
            continue
        pdb_path = _ensure_pdb(folds[target])
        chains = sema.get_struc_seq(pdb_path)
        chain = sorted(chains)[0]
        if len(chains) > 1:
            print(f"[WARNING] {target}: {len(chains)} chains in {os.path.basename(pdb_path)} "
                  f"({sorted(chains)}); scoring chain {chain!r}. A native antigen model should "
                  f"be a monomer -- check this is the intended chain.")
        cache = os.path.join(tool_runs_dir, f"sema3d_{target}.json")
        if os.path.isfile(cache):
            blob = json.load(open(cache))
            cached[target] = (blob["sequence"], blob["scores"])
            print(f"[INFO] {target}: cached SEMA-3D result ({len(blob['scores'])} residues).")
        else:
            # Per-residue pLDDT from the AlphaFold model, used to mask 3Di
            # tokens in regions the fold does not determine. Without it, a
            # disordered tail contributes invented structure to the epitope
            # prediction -- see get_struc_seq's docstring and deviation #21.
            plddt = _model_plddt(folds[target])
            resolved[target] = (pdb_path, chain, plddt)

    if resolved:
        print(f"[INFO] Scoring {len(resolved)} antigen(s) with the {sema.N_ENSEMBLE}-model "
              f"ensemble (streamed, one checkpoint resident at a time)...")
        predicted = sema.predict_many(MODELS_DIR, resolved, dev, n=sema.N_ENSEMBLE)
        for target, (model_seq, scores) in predicted.items():
            model_scores = [float(x) for x in scores]
            cached[target] = (model_seq, model_scores)
            # resolved[target] is (pdb_path, chain, plddt) -- index rather than
            # unpack, so adding a fourth element later cannot break this again.
            spec = resolved[target]
            pdb_path, chain = spec[0], spec[1]
            n_masked = sum(1 for v in (spec[2] or []) if v < sema.PLDDT_MASK_THRESHOLD)
            with open(os.path.join(tool_runs_dir, f"sema3d_{target}.json"), "w") as fh:
                json.dump({"pdb": os.path.basename(pdb_path), "chain": chain,
                           "plddt_masked_residues": n_masked,
                           "plddt_mask_threshold": sema.PLDDT_MASK_THRESHOLD,
                           "sequence": model_seq, "scores": model_scores}, fh)

    corroboration = {}
    per_residue_out = {}
    for target in runnable:
        if target not in cached:
            continue
        peptides = bcell_by_target.get(target, set())
        ref_seq = _reference_sequence(target)
        model_seq, model_scores = cached[target]

        mapped, covered = _map_scores_to_reference(model_seq, model_scores, ref_seq)
        pct_cov = 100.0 * covered / len(ref_seq)
        print(f"       reference {len(ref_seq)} aa | structure {len(model_seq)} aa | "
              f"coverage {pct_cov:.1f}%")
        if pct_cov < 50.0:
            print(f"[WARNING] {target}: the structure covers under half the reference sequence. "
                  f"Epitopes outside the modelled region cannot be corroborated.")
        per_residue_out[target] = {"reference_length": len(ref_seq),
                                   "structure_length": len(model_seq),
                                   "coverage_pct": round(pct_cov, 1),
                                   "scores": mapped}

        for pep in peptides:
            start = ref_seq.find(pep)
            if start == -1:
                # Not in this target's own mature subunit -- expected, because
                # the Target label records which Entrez query fetched the
                # peptide's parent record, not which protein it lies in
                # (deviation #26). Scoring it against the labelled antigen's
                # fold would be scoring the wrong protein, so it is deferred to
                # _score_against_true_owner() below, which looks for the fold of
                # the protein the peptide is ACTUALLY in. Only if no scored
                # antigen contains it does it stay UNSCREENED.
                corroboration.setdefault(pep, (None, None, "NOT_IN_SUBUNIT"))
                continue
            pct, n_scored = _score_epitope(mapped, start, len(pep))
            if pct is None:
                corroboration[pep] = (None, None, "NOT_MODELLED")
                continue
            corroboration[pep] = (round(pct, 1), n_scored,
                                  "YES" if pct >= OVERLAP_THRESHOLD * 100 else "NO")

    # ---- rescue mislabelled epitopes against the fold they really belong to --
    # Without this, deviation #26's mislabelling silently costs SEMA coverage:
    # an epitope labelled HIV_p17 that actually lies in p24 is searched only in
    # the p17 subunit, is not found, and is reported UNSCREENED even though the
    # p24 fold was scored in this very run. The peptide is located by sequence
    # across every antigen that WAS scored, and attributed to that one.
    rescued = 0
    for pep, result in list(corroboration.items()):
        if result[2] != "NOT_IN_SUBUNIT":
            continue
        for other in runnable:
            if other not in per_residue_out or other not in cached:
                continue
            other_ref = _reference_sequence(other)
            if other_ref is None:
                continue
            start = other_ref.find(pep)
            if start == -1:
                continue
            pct, n_scored = _score_epitope(per_residue_out[other]["scores"], start, len(pep))
            if pct is None:
                continue
            corroboration[pep] = (round(pct, 1), n_scored,
                                  "YES" if pct >= OVERLAP_THRESHOLD * 100 else "NO")
            rescued += 1
            print(f"[INFO] {pep}: not in its labelled target's subunit; scored against "
                  f"{other}, which actually contains it ({pct:.1f}% overlap).")
            break
    if rescued:
        print(f"[INFO] {rescued} mislabelled epitope(s) recovered by true-provenance lookup.")

    # ---- emit the enriched table --------------------------------------------
    fieldnames = original_fields + ["SEMA_Overlap_Pct", "SEMA_Corroborated", "SEMA_Residues_Scored"]
    out_rows = []
    counts = {"YES": 0, "NO": 0, "UNSCREENED": 0}
    for r in rows:
        clean = {k: r[k] for k in original_fields}
        if r.get("Type") == "B-cell":
            result = corroboration.get(r["Peptide"])
            if result is None or result[2] in ("NOT_IN_SUBUNIT", "NOT_MODELLED"):
                clean["SEMA_Overlap_Pct"] = ""
                clean["SEMA_Corroborated"] = "UNSCREENED"
                clean["SEMA_Residues_Scored"] = ""
                counts["UNSCREENED"] += 1
            else:
                pct, n_scored, call = result
                clean["SEMA_Overlap_Pct"] = pct
                clean["SEMA_Corroborated"] = call
                clean["SEMA_Residues_Scored"] = n_scored
                counts[call] += 1
        else:
            clean["SEMA_Overlap_Pct"] = ""
            clean["SEMA_Corroborated"] = ""
            clean["SEMA_Residues_Scored"] = ""
        out_rows.append(clean)

    ts = datetime.now().strftime("%Y%m%d_%H%M")
    if write_back:
        # NOTE: Phase1Db/1Dc/1Da all use "%Y%m%d_%H%M" (no dashes), and
        # Phase1Dc_benchmark.py picks its input by a plain string sort, not by
        # ctime. A dash-formatted timestamp would sort BEFORE same-day no-dash
        # names ('-' < '0' in ASCII) and be silently skipped as "not the
        # latest" -- so this must match the sibling files' format exactly.
        out_path = os.path.join(input_folder, f"Phase1Db_Elite_{ts}.csv")
        print("\n[WARNING] write_back=True: this REPLACES the pool Phase 1G reads and will "
              "change the construct on the next 1G run.")
    else:
        out_path = os.path.join(sens_dir, f"Phase1De_SEMA3D_{ts}.csv")

    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out_rows)

    residues_path = os.path.join(sens_dir, f"Phase1De_SEMA3D_perResidue_{ts}.json")
    with open(residues_path, "w") as f:
        json.dump({"epitope_score_cut": EPITOPE_SCORE_CUT,
                   "score_units": "log-scaled expected number of antibody contacts, ln(1+n)",
                   "targets": per_residue_out}, f)

    print(f"\n[INFO] B-cell SEMA-3D corroboration -- corroborated: {counts['YES']} | "
          f"not corroborated: {counts['NO']} | unscreened: {counts['UNSCREENED']}")
    print(f"[INFO] Enriched table   : {out_path}")
    print(f"[INFO] Per-residue maps : {residues_path}")
    if not write_back:
        print("[INFO] SENSITIVITY OUTPUT -- the Phase 1G production pool was NOT modified. "
              "Run Phase1G_sensitivity.py to see whether this would change the construct.")
    print("=" * 90 + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="SEMA-3D conformational B-cell corroboration.")
    ap.add_argument("--write-back", action="store_true",
                    help="Write the enriched pool back into Phase1Db/ (CHANGES THE CONSTRUCT on "
                         "the next Phase 1G run). Off by default.")
    ap.add_argument("--targets", nargs="*", default=None, help="Subset of targets to score.")
    ap.add_argument("--device", default="auto", help="auto | mps | cpu")
    args = ap.parse_args()
    run_step1de_conformational(write_back=args.write_back, targets=args.targets, device=args.device)
