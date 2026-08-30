import os
import sys
import csv
import glob
from datetime import datetime
from collections import defaultdict, Counter

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

import phase1_common as common

# =============================================================================
# PHASE 1A PROVENANCE CORRECTION -- REPORT ONLY, CHANGES NOTHING
#
# THE DEFECT. Phase1A_sequenceRetrival.py defines its four HIV targets as four
# Entrez QUERY STRINGS:
#       "HIV_gp120": 'gp120 AND "HIV-1"[Organism] AND CRF01_AE'
#       "HIV_gp41" : 'gp41  AND "HIV-1"[Organism] AND CRF01_AE'
#       "HIV_p24"  : 'p24   AND "HIV-1"[Organism] AND CRF01_AE'
#       "HIV_p17"  : 'p17   AND "HIV-1"[Organism] AND CRF01_AE'
# and stores whatever record comes back WHOLE, labelled with the query name.
# For CRF01_AE, NCBI returns POLYPROTEINS, not mature subunits:
#   HIV_gp120_Var_01 and HIV_gp41_Var_01 are the SAME 856-aa Env  (YES72107.1)
#   HIV_p17_Var_01   and HIV_p24_Var_01   are the SAME 1437-aa Gag-Pol (YES72110.1)
# Phase 1Da then slid a k-mer window over each stored record and tagged every
# peptide with that record's target label. So an epitope's Target says which
# QUERY retrieved its parent record -- NOT which protein the peptide is in.
#
# WHY IT MATTERS. Three separate claims depend on the label being real:
#   (1) the per-antigen representation table in the manuscript,
#   (2) MAX_EPITOPES_PER_TARGET_PER_CLASS = 2 in Phase 1G, which is enforced
#       on the label and is therefore enforceable only if the label is true,
#   (3) the "which antigen contributes no B-cell epitopes" finding.
#
# WHAT THIS STEP DOES *NOT* DO. It does not rewrite the Phase 1F pool, the
# Phase 1G construct, or any Phase 1A FASTA. Phase 1G selects its input via
# common.latest_file() over Step_Outputs/.../Phase1F/Filtered -- writing a
# corrected pool there would silently change the construct on the next 1G run.
# The construct is frozen (it has been through all of Phase 2 and been docked),
# so this emits a standalone report and nothing else. Whether to re-run 1F/1G
# under corrected labels is a separate decision that WOULD change the design.
# =============================================================================

# Where the subunit boundaries come from, and how they were derived, is
# documented on HIV_SUBUNIT_RANGES in phase1_common.py. Summary: the GenBank
# records carry no mat_peptide features, so the ranges come from a global
# BLOSUM62 alignment of each Var_01 record against its annotated HXB2
# reference (P04578 Env, P04585 Gag-Pol) with the reference's UniProt CHAIN
# features mapped through the alignment -- corroborated independently by
# sequence landmarks (Env furin site + fusion peptide; Gag MA/CA junction).
HIV_TARGETS = tuple(common.HIV_SUBUNIT_RANGES.keys())

# Var_01 is the canonical isolate per target -- the same reference Phase 1De
# uses. An epitope selected from a NON-Var_01 variant will not be found here;
# that is reported as UNRESOLVED rather than guessed at.
def _var01_record(target):
    folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1A")
    if not os.path.isdir(folder):
        return None, None
    matches = sorted(f for f in os.listdir(folder) if f.startswith(f"{target}_Var_01_"))
    if not matches:
        return None, None
    path = os.path.join(folder, matches[0])
    with open(path) as f:
        seq = "".join(l.strip() for l in f if not l.startswith(">")).upper()
    return seq, matches[0]


def _variant_source(peptide):
    """
    Phase 1F records which VARIANT record each surviving peptide was drawn
    from. Epitopes selected from a non-Var_01 isolate are absent from the
    Var_01 reference, so they must be resolved against their own source record
    instead of being reported as unknown.
    """
    folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1F", "Filtered")
    latest = common.latest_file(folder, suffix=".csv")
    if latest is None:
        return None
    with open(latest) as f:
        for row in csv.DictReader(f):
            if row.get("Peptide") == peptide and row.get("Variant"):
                return row["Variant"]
    return None


def _load_variant(filename):
    path = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1A", filename)
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        return "".join(l.strip() for l in f if not l.startswith(">")).upper()


def _true_target(target, peptide, parent_seq):
    """
    Returns (true_target, position_1based, status, source).

    status is one of:
      OK           -- located in the Var_01 parent, subunit from the exact,
                      alignment-derived HIV_SUBUNIT_RANGES
      OK_VARIANT   -- absent from Var_01, resolved instead against the variant
                      record Phase 1F says it came from, using the cleavage
                      landmarks (see subunit_boundaries_by_landmark). Variant
                      records differ in length (Env here spans 854-868 aa), so
                      Var_01's fixed offsets do not transfer to them.
      UNRESOLVED   -- not locatable in either; reported, never guessed.
    """
    if parent_seq is not None:
        idx = parent_seq.find(peptide)
        if idx != -1:
            pos = idx + 1
            return common.subunit_of(target, pos), pos, "OK", "Var_01"

    variant = _variant_source(peptide)
    if variant:
        vseq = _load_variant(variant)
        if vseq:
            idx = vseq.find(peptide)
            if idx != -1:
                pos = idx + 1
                parent = common.HIV_SUBUNIT_RANGES[target]["parent"]
                bounds = common.subunit_boundaries_by_landmark(vseq, parent)
                if bounds:
                    for name, (s, e) in bounds.items():
                        if s <= pos <= e:
                            return name, pos, "OK_VARIANT", variant
                    if parent == "Gag-Pol":
                        return "Gag_downstream", pos, "OK_VARIANT", variant
    return None, None, "UNRESOLVED", ""


def _load_construct_epitopes():
    """
    Reads Epitope_Provenance ("PEPTIDE:Target;PEPTIDE:Target;...") and
    Boundary_Map ("...;MHC-I:PEP[start-end];...") from the latest Phase 1G
    construct CSV, so each epitope carries both its labelled target and its
    epitope CLASS (needed for the per-target-per-class cap check).
    """
    folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1G")
    latest = common.latest_file(folder, suffix=".csv")
    if latest is None:
        return None, None, None
    with open(latest) as f:
        row = next(csv.DictReader(f))

    labelled = {}
    for item in row["Epitope_Provenance"].split(";"):
        if ":" in item:
            pep, target = item.rsplit(":", 1)
            labelled[pep.strip()] = target.strip()

    # Boundary_Map segments look like "MHC-I:SMVISLLSM[50-59]" for epitopes and
    # "AAY[59-62]" / "adjuvant[0-45]" for everything else -- only the former
    # carry a class prefix, which is exactly the filter we want.
    klass = {}
    for seg in row["Boundary_Map"].split(";"):
        seg = seg.strip()
        for cls in ("MHC-I", "MHC-II", "B-cell"):
            prefix = cls + ":"
            if seg.startswith(prefix):
                pep = seg[len(prefix):].split("[")[0]
                klass[pep] = cls
                break
    return labelled, klass, os.path.basename(latest)



# =============================================================================
# ANTIGEN IDENTITY CHECK -- is each target the protein the study MEANT to use?
#
# The subunit check above assumes the retrieved RECORD is the right protein and
# only asks which part of it an epitope sits in. That assumption needs testing
# on its own, because Phase 1A selects antigens by literal gene-name string
# ('B5R AND "Monkeypox virus"[Organism]'), and orthologous poxvirus genes do
# NOT share names across species. Each Var_01 record is therefore aligned
# against a panel of candidate UniProt references and the BEST match reported,
# so a wrong antigen surfaces as a RESULT rather than as an assumption.
#
# The intended Mpox antigens are the classic neutralising trio, named by their
# vaccinia genes: A33R (EEV glycoprotein), L1R (IMV membrane protein) and B5R
# (EEV glycoprotein). The MPXV orthologs are OPG161 (=A35R), OPG095 (=M1R) and
# OPG190 (=B6R). Note the last: MPXV has its OWN gene named B5R and it is NOT
# the vaccinia-B5R ortholog -- it is OPG189, an ankyrin repeat protein.
# =============================================================================
ANTIGEN_REFERENCES = {
    "Mpox_A35R": {
        "intended": "OPG161 -- VACV A33R ortholog, EEV envelope glycoprotein",
        "candidates": {"A0A7H0DND2": "OPG161 (EEV glycoprotein, VACV A33R ortholog)"},
    },
    "Mpox_L1R": {
        "intended": "OPG095 / M1R -- VACV L1R ortholog, IMV membrane protein",
        "candidates": {"M1LBP0": "OPG095 / M1R (IMV membrane protein, VACV L1R ortholog)"},
    },
    "Mpox_B5R": {
        "intended": "OPG190 / B6R -- VACV B5R ortholog, EEV envelope glycoprotein",
        "candidates": {
            "P0DTN2": "OPG190 / B6R (EEV glycoprotein, VACV B5R ortholog)",
            "A0A7H0DNF1": "OPG189 (ankyrin repeat protein -- MPXV's own gene named B5R)",
        },
    },
    "HIV_gp120": {"intended": "Env gp160 polyprotein (HXB2 reference)",
                  "candidates": {"P04578": "Env gp160 (HXB2)"}},
    "HIV_gp41":  {"intended": "Env gp160 polyprotein (HXB2 reference)",
                  "candidates": {"P04578": "Env gp160 (HXB2)"}},
    "HIV_p17":   {"intended": "Gag-Pol polyprotein (HXB2 reference)",
                  "candidates": {"P04585": "Gag-Pol (HXB2)"}},
    "HIV_p24":   {"intended": "Gag-Pol polyprotein (HXB2 reference)",
                  "candidates": {"P04585": "Gag-Pol (HXB2)"}},
}

# Below this a target is not the protein it claims to be. Orthologous poxvirus
# surface antigens across Orthopoxvirus species sit well above 90%; anything in
# the 20-30% band is shared-fold background, not orthology.
IDENTITY_FLOOR_PCT = 50.0


def _global_identity(query, reference):
    """Percent identity over the shorter sequence, BLOSUM62 global alignment."""
    from Bio import Align
    aligner = Align.PairwiseAligner(
        mode="global", open_gap_score=-11, extend_gap_score=-1,
        substitution_matrix=Align.substitution_matrices.load("BLOSUM62"))
    aln = aligner.align(query, reference)[0]
    a, b = str(aln[0]), str(aln[1])
    matches = sum(1 for x, y in zip(a, b) if x == y and x != "-")
    return 100.0 * matches / min(len(query), len(reference))


def check_antigen_identity():
    """Aligns each Var_01 record against its candidate references; best match wins."""
    ref_dir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1A", "Reference_Antigens")
    rows = []
    print("\n" + "-" * 100)
    print("ANTIGEN IDENTITY -- is each retrieved record the protein the study intended?")
    for target, spec in ANTIGEN_REFERENCES.items():
        seq, fname = _var01_record(target)
        if seq is None:
            continue
        scored = []
        for acc, desc in spec["candidates"].items():
            path = os.path.join(ref_dir, f"{acc}.fasta")
            if not os.path.isfile(path):
                continue
            ref = "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()
            scored.append((_global_identity(seq, ref), acc, desc, len(ref)))
        if not scored:
            continue
        scored.sort(reverse=True)
        best_pct, best_acc, best_desc, _bl = scored[0]
        intended_acc = next(iter(spec["candidates"]))
        intended_pct = next((p for p, a, _d, _l in scored if a == intended_acc), 0.0)
        ok = (best_acc == intended_acc and best_pct >= IDENTITY_FLOOR_PCT)
        print(f"\n  {target}  ({len(seq)} aa, {fname})")
        print(f"    intended: {spec['intended']}")
        for pct, acc, desc, rlen in scored:
            mark = "  <-- BEST MATCH" if acc == best_acc else ""
            print(f"      {pct:5.1f}% vs {acc} ({rlen} aa)  {desc}{mark}")
        print(f"    VERDICT: {'CORRECT' if ok else 'WRONG PROTEIN'}")
        rows.append({"Target": target, "Var01_Length": len(seq),
                     "Intended": spec["intended"], "Best_Match": best_desc,
                     "Best_Match_Accession": best_acc, "Best_Identity_Pct": round(best_pct, 1),
                     "Intended_Identity_Pct": round(intended_pct, 1),
                     "Verdict": "CORRECT" if ok else "WRONG_PROTEIN"})
    return rows


# =============================================================================
# VARIANT-POOL HOMOGENEITY -- is every "variant" of a target the same protein?
#
# The identity check above tests only Var_01, the canonical isolate. But every
# target carries up to 30 variant records, and CONSERVANCY is computed across
# that whole pool. If a pool mixes two different genes, its conservancy figures
# are comparing unrelated proteins and are meaningless -- and epitopes can be
# selected from a record that is not the antigen at all.
#
# That is not hypothetical. Phase 1A's query 'L1R AND "Monkeypox virus"' also
# matches MPXV OPG053, a DIFFERENT entry-fusion-complex protein that submitters
# annotate with similar wording. Half the Mpox_L1R pool is OPG053, and both
# L1R B-cell epitopes in the construct come from it rather than from L1R.
#
# Each variant is aligned to its own Var_01. Below IDENTITY_FLOOR_PCT the
# record is a different protein, not a strain variant.
# =============================================================================
_STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def _sanitise(seq):
    """
    Drops non-standard residue codes (X, B, Z, U ...). Some GenBank records
    carry them and BLOSUM62 has no entry for them, so the aligner raises.
    """
    return "".join(c for c in seq if c in _STANDARD_AA)


def check_pool_homogeneity():
    import glob
    from collections import Counter
    folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1A")
    rows = []
    print("\n" + "-" * 100)
    print("VARIANT-POOL HOMOGENEITY -- are all variants of a target the same protein?")
    for target in list(ANTIGEN_REFERENCES.keys()):
        files = sorted(glob.glob(os.path.join(folder, f"{target}_Var_*.fasta")))
        if not files:
            continue
        seqs = []
        for f in files:
            with open(f) as fh:
                seqs.append((os.path.basename(f),
                             "".join(l.strip() for l in fh if not l.startswith(">")).upper()))
        ref = next((s for n, s in seqs if "_Var_01_" in n), None)
        if ref is None:
            continue
        ids = [(n, _global_identity(_sanitise(s), _sanitise(ref))) for n, s in seqs]
        outliers = [(n, i) for n, i in ids if i < IDENTITY_FLOOR_PCT]
        lens = Counter(len(s) for _n, s in seqs)
        verdict = "MIXED_PROTEINS" if outliers else "HOMOGENEOUS"
        print(f"\n  {target}: {len(seqs)} records | lengths "
              f"{', '.join(f'{L}({c})' for L, c in sorted(lens.items()))[:60]}")
        print(f"    min identity to Var_01: {min(i for _n, i in ids):.1f}%   VERDICT: {verdict}")
        for n, i in outliers[:3]:
            print(f"      DIFFERENT PROTEIN: {n} ({i:.1f}%)")
        if len(outliers) > 3:
            print(f"      ... and {len(outliers) - 3} more")
        rows.append({"Target": target, "N_Records": len(seqs),
                     "N_Different_Protein": len(outliers),
                     "Min_Identity_To_Var01_Pct": round(min(i for _n, i in ids), 1),
                     "Distinct_Lengths": len(lens), "Verdict": verdict,
                     "Example_Outlier": outliers[0][0] if outliers else ""})
    return rows

def run_provenance_correction():
    common.print_banner("PHASE 1A PROVENANCE CORRECTION: TRUE HIV SUBUNIT PER EPITOPE (REPORT ONLY)")

    labelled, klass, src = _load_construct_epitopes()
    if labelled is None:
        print("[ERROR] No Phase 1G construct CSV found -- nothing to check.")
        return

    parents = {t: _var01_record(t) for t in HIV_TARGETS}
    for t, (seq, fname) in parents.items():
        if seq is None:
            print(f"[WARNING] No Var_01 record found for {t}.")
        else:
            info = common.HIV_SUBUNIT_RANGES[t]
            print(f"[INFO] {t:<10} parent={info['parent']:<8} {fname}  ({len(seq)} aa) "
                  f"-> subunit {info['start']}-{info['end']}")

    print(f"\n[INFO] Construct source: {src}")
    print("-" * 100)

    rows = []
    mismatches = 0
    unresolved = 0
    for pep, lab in sorted(labelled.items(), key=lambda kv: (kv[1], kv[0])):
        if not lab.startswith("HIV"):
            # Mpox targets are genuine single proteins -- their Var_01 file IS
            # the mature protein, so the label cannot be wrong this way.
            rows.append({"Peptide": pep, "Class": klass.get(pep, ""), "Labelled_Target": lab,
                         "True_Target": lab, "Parent_Record": "", "Position_In_Parent": "",
                         "Status": "NOT_APPLICABLE", "Mislabelled": "NO",
                         "Resolved_Against": ""})
            continue

        parent_seq, _fname = parents.get(lab, (None, None))
        true_t, pos, status, source = _true_target(lab, pep, parent_seq)
        if status == "UNRESOLVED":
            unresolved += 1
            rows.append({"Peptide": pep, "Class": klass.get(pep, ""), "Labelled_Target": lab,
                         "True_Target": "", "Parent_Record": common.HIV_SUBUNIT_RANGES[lab]["accession"],
                         "Position_In_Parent": "", "Status": "UNRESOLVED", "Mislabelled": "UNKNOWN",
                         "Resolved_Against": ""})
            continue

        bad = (true_t != lab)
        mismatches += bad
        rows.append({"Peptide": pep, "Class": klass.get(pep, ""), "Labelled_Target": lab,
                     "True_Target": true_t or "", "Parent_Record": common.HIV_SUBUNIT_RANGES[lab]["accession"],
                     "Position_In_Parent": pos, "Status": status,
                     "Mislabelled": "YES" if bad else "NO",
                     "Resolved_Against": source})

    print(f"{'labelled':<11} {'class':<8} {'epitope':<18} {'pos':<7} {'true target':<16} flag")
    for r in rows:
        if r["Status"] == "NOT_APPLICABLE":
            continue
        flag = {"YES": "MISLABELLED", "UNKNOWN": "UNRESOLVED"}.get(r["Mislabelled"], "")
        if r["Status"] == "OK_VARIANT":
            flag = (flag + "  [via " + r["Resolved_Against"] + "]").strip()
        print(f"{r['Labelled_Target']:<11} {r['Class']:<8} {r['Peptide']:<18} "
              f"{str(r['Position_In_Parent']):<7} {r['True_Target']:<16} {flag}")

    # ---- Corrected representation -------------------------------------------
    lab_counts, true_counts = Counter(), Counter()
    for r in rows:
        if not r["Labelled_Target"].startswith("HIV"):
            continue
        lab_counts[r["Labelled_Target"]] += 1
        true_counts[r["True_Target"] or "UNRESOLVED"] += 1

    print("\n" + "-" * 100)
    print("PER-ANTIGEN REPRESENTATION (HIV only)")
    all_t = sorted(set(lab_counts) | set(true_counts))
    print(f"  {'target':<18} {'as labelled':>12} {'actual':>8}")
    for t in all_t:
        print(f"  {t:<18} {lab_counts.get(t,0):>12} {true_counts.get(t,0):>8}")

    # ---- Cap violations ------------------------------------------------------
    # Phase 1G enforces MAX_EPITOPES_PER_TARGET_PER_CLASS = 2 on the LABEL.
    # Recount on the true target to see whether the real constraint held.
    by_true = defaultdict(list)
    for r in rows:
        if r["True_Target"] and r["Class"]:
            by_true[(r["True_Target"], r["Class"])].append(r["Peptide"])
    violations = {k: v for k, v in by_true.items() if len(v) > 2}
    print("\n" + "-" * 100)
    if violations:
        print("MAX_EPITOPES_PER_TARGET_PER_CLASS = 2 VIOLATED ON TRUE TARGETS:")
        for (t, c), peps in sorted(violations.items()):
            print(f"  {t} / {c}: {len(peps)} epitopes -> {peps}")
    else:
        print("MAX_EPITOPES_PER_TARGET_PER_CLASS = 2 holds on true targets.")

    # ---- B-cell representation gap ------------------------------------------
    bcell_true = {r["True_Target"] for r in rows if r["Class"] == "B-cell" and r["True_Target"]}
    bcell_lab = {r["Labelled_Target"] for r in rows if r["Class"] == "B-cell"}
    missing_true = [t for t in HIV_TARGETS if t not in bcell_true]
    missing_lab = [t for t in HIV_TARGETS if t not in bcell_lab]
    print("\n" + "-" * 100)
    print(f"HIV antigens with NO B-cell epitope -- as labelled: {missing_lab or 'none'}")
    print(f"HIV antigens with NO B-cell epitope -- actual     : {missing_true or 'none'}")

    # ---- Write the report ----------------------------------------------------
    out_dir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1A", "Provenance_Correction")
    os.makedirs(out_dir, exist_ok=True)
    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    out_path = os.path.join(out_dir, f"Phase1A_ProvenanceCorrection_{ts}.csv")
    fields = ["Peptide", "Class", "Labelled_Target", "True_Target", "Parent_Record",
              "Position_In_Parent", "Status", "Mislabelled", "Resolved_Against"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    identity_rows = check_antigen_identity()
    pool_rows = check_pool_homogeneity()
    if pool_rows:
        pool_path = os.path.join(out_dir, f"Phase1A_PoolHomogeneity_{ts}.csv")
        with open(pool_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(pool_rows[0].keys()))
            w.writeheader()
            w.writerows(pool_rows)
        mixed = [r["Target"] for r in pool_rows if r["Verdict"] == "MIXED_PROTEINS"]
        print(f"\n[INFO] Pool homogeneity report: {pool_path}")
        if mixed:
            print(f"[WARNING] {len(mixed)} target(s) have variant pools mixing DIFFERENT "
                  f"proteins: {mixed}")
            print( "[WARNING] Conservancy for those targets is computed across unrelated "
                   "sequences and cannot be interpreted.")
    if identity_rows:
        id_path = os.path.join(out_dir, f"Phase1A_AntigenIdentity_{ts}.csv")
        with open(id_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(identity_rows[0].keys()))
            w.writeheader()
            w.writerows(identity_rows)
        wrong = [r["Target"] for r in identity_rows if r["Verdict"] == "WRONG_PROTEIN"]
        print(f"\n[INFO] Antigen identity report: {id_path}")
        if wrong:
            print(f"[WARNING] {len(wrong)} target(s) are NOT the intended protein: {wrong}")

    print("\n" + "-" * 100)
    n_variant = sum(1 for r in rows if r["Status"] == "OK_VARIANT")
    print(f"[INFO] {n_variant} epitope(s) resolved against their own variant record "
          f"(absent from Var_01).")
    print(f"[INFO] {mismatches} mislabelled, {unresolved} unresolved, "
          f"{sum(1 for r in rows if r['Mislabelled']=='NO' and r['Status']=='OK')} correct "
          f"(of {sum(1 for r in rows if r['Labelled_Target'].startswith('HIV'))} HIV epitopes).")
    print(f"[INFO] Report written to: {out_path}")
    print("[INFO] REPORT ONLY -- no pool, construct or FASTA was modified.")
    print("=" * 100 + "\n")


if __name__ == "__main__":
    run_provenance_correction()
