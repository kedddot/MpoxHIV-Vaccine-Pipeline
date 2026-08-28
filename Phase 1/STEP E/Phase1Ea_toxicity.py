import os, sys, time, csv, re, shutil, subprocess
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

import phase1_common as common

# =============================================================================
# CONFIGURE -- override via env vars for a different machine (e.g. DOST-COARE HPC)
# =============================================================================
TOXINPRED2_BINARY = os.environ.get("TOXINPRED2_BINARY", "/opt/miniconda3/envs/phase2/bin/toxinpred2")
HEMOPI2_BINARY = os.environ.get("HEMOPI2_BINARY", "/opt/miniconda3/envs/phase2/bin/hemopi2_classification")
BLASTP_BINARY = os.environ.get("BLASTP_BINARY", "/opt/miniconda3/envs/phase2/bin/blastp")
TOXPROT_BLAST_DB = os.environ.get("TOXPROT_BLAST_DB", os.path.join(_PROJECT_ROOT, "toxprot_db", "toxprot"))

# >50% conservancy = "highly acceptable" per the paper; an explicit choice
# now that Phase 1Dc writes each tier to its own subdirectory (see that
# script's fix note -- this used to be a ctime accident that always grabbed
# the 100% tier instead).
CONSERVANCY_TIER_DIR = "Min_50pct"

# Cap on how many candidates advance to Phase 1Eb's manual AllerTOP/AllergenFP
# screening, per (target x epitope-class) pair. See the prioritization block at
# the end of run_step1ea_toxicity() for the rationale.
TOP_N_PER_TARGET_PER_CLASS = 10


def _binary_exists(path):
    return shutil.which(path) is not None or os.path.isfile(path)


def local_tools_available():
    return (
        _binary_exists(TOXINPRED2_BINARY)
        and _binary_exists(HEMOPI2_BINARY)
        and _binary_exists(BLASTP_BINARY)
        and (os.path.isfile(TOXPROT_BLAST_DB + ".phr") or os.path.isfile(TOXPROT_BLAST_DB + ".pin"))
    )


# =============================================================================
# REAL TOOLS -- ToxinPred2, HemoPI2, BLASTP vs local Tox-Prot DB.
# Verified locally installed and callable (see Important.rtf's conda env
# setup and Phase 2 STEP A, which uses the same binaries). Column names
# below were confirmed against real CLI output, not guessed:
#   ToxinPred2  -> ID,Sequence,ML_Score,Prediction        (Prediction: Toxin/Non-Toxin)
#   HemoPI2     -> SeqID,Sequence,ESM Score,Prediction    (ESM Score is the probability)
#   BLASTP      -> qseqid,sseqid,evalue,pident,length,qlen
# =============================================================================

def run_toxinpred2(fasta_path, work_dir):
    """
    Runs ToxinPred2's HYBRID model (-m 2: RF + BLAST + MERCI), which reports
    three separable evidence channels per peptide:
      ML Score    -- amino-acid-composition Random Forest
      BLAST Score -- homology to known toxins (+0.5 when a hit exists)
      MERCI Score -- presence of a known toxin motif

    Two genuine upstream bugs blocked -m 2 and were patched in the installed
    package (both are real defects in toxinpred2 itself, not environment
    issues):
      1. BLAST_processor() takes a parameter `name1` but its else-branch
         referenced an out-of-scope `seqid` -> NameError.
      2. hybrid() computed `df8.sum(axis=1)` across ALL columns including the
         string 'Subject' column -> TypeError on modern pandas (older pandas
         silently skipped non-numeric columns).

    Returns {peptide_id: {"ml": float, "blast": float, "merci": float}} so the
    caller can weigh composition evidence separately from homology/motif
    evidence -- see the composition-calibration note in METHODOLOGY_NOTE_TEXT.
    """
    output_csv = os.path.join(work_dir, "toxinpred2_hybrid.csv")
    cmd = [TOXINPRED2_BINARY, "-i", fasta_path, "-o", output_csv, "-t", "0.6", "-m", "2", "-d", "2"]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
    if result.returncode != 0 or not os.path.isfile(output_csv):
        raise RuntimeError(f"ToxinPred2 failed (exit {result.returncode}):\n{result.stderr}")
    predictions = {}
    with open(output_csv) as f:
        for row in csv.DictReader(f):
            predictions[row["Subject"]] = {
                "ml": float(row["ML Score"]),
                "blast": float(row["BLAST Score"]),
                "merci": float(row["MERCI Score"]),
            }
    return predictions


def run_hemopi2(fasta_path, work_dir):
    """Model 3 (ESM2-t6): avoids the same MERCI subprocess step used by
    ToxinPred2's model 2, for consistency/caution given that bug."""
    output_filename = "hemopi2_out.csv"
    output_csv = os.path.join(work_dir, output_filename)
    cmd = [HEMOPI2_BINARY, "-i", fasta_path, "-o", output_filename, "-wd", work_dir, "-d", "2", "-m", "3"]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
    if result.returncode != 0 or not os.path.isfile(output_csv):
        raise RuntimeError(f"HemoPI2 failed (exit {result.returncode}):\n{result.stderr}")
    scores = {}
    with open(output_csv) as f:
        for row in csv.DictReader(f):
            scores[row["SeqID"]] = float(row["ESM Score"])
    return scores


def run_blastp_toxprot(fasta_path, work_dir):
    """
    Best (lowest-E-value) Tox-Prot hit per query peptide, with identity,
    alignment coverage, and E-value all reported so the paper's 3-of-3
    rule (identity>=80% AND coverage>=80% AND E<=1e-5) can be evaluated
    explicitly -- the prior version of this script only checked E-value.
    """
    output_tsv = os.path.join(work_dir, "blastp_toxprot.tsv")
    cmd = [
        BLASTP_BINARY, "-query", fasta_path, "-db", TOXPROT_BLAST_DB,
        "-outfmt", "6 qseqid sseqid evalue pident length qlen",
        "-evalue", "1e-3", "-max_target_seqs", "1", "-out", output_tsv,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
    if result.returncode != 0:
        raise RuntimeError(f"BLASTP failed (exit {result.returncode}):\n{result.stderr}")

    best_hits = {}
    if os.path.isfile(output_tsv):
        with open(output_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 6:
                    continue
                qseqid, sseqid, evalue, pident, length, qlen = parts[0], parts[1], float(parts[2]), float(parts[3]), int(parts[4]), int(parts[5])
                if qseqid not in best_hits or evalue < best_hits[qseqid]["evalue"]:
                    best_hits[qseqid] = {"evalue": evalue, "pident": pident, "coverage": (length / qlen) * 100 if qlen else 0.0, "subject": sseqid}
    return best_hits


# =============================================================================
# MANUAL FALLBACK -- used only if the local tools above aren't available on
# this machine (e.g. a fresh checkout without the phase2 conda env active).
# =============================================================================

METHODOLOGY_NOTE_TEXT = """# Phase 1Ea — Toxicity Screening Methodology Note

**Applies to:** Section II.C.I.E of the proposal (ToxinPred, HemoPI, BLASTP vs. UniProt Tox-Prot).

## What this script does
On a machine with the `phase2` conda environment's ToxinPred2, HemoPI2, and
local Tox-Prot BLAST database available, all three checks run automatically:
ToxinPred2 (model 2, hybrid RF+BLAST+MERCI), HemoPI2 (model 3, ESM2-t6), and
BLASTP against a local UniProt Tox-Prot database built with `makeblastdb`.

## ToxinPred composition-score calibration (IMPORTANT — cite this)
ToxinPred2's amino-acid-composition Random Forest is **not calibrated for
9–16 aa peptides**; it is trained to discriminate toxin *proteins*. Applied
to short epitopes it produces systematic false positives. Measured directly
against known-benign controls, at the tool's own default threshold (0.6):

| Control sequence | Identity | ML Score | Composition call |
|---|---|---|---|
| `AAYGPGPGKKAAY`   | Vaccine linkers only, no biological sequence | 0.621 | Toxin |
| `DAHKSEVAHRFKDLG` | Human serum albumin N-terminus | 0.688 | Toxin |
| `VKVGVNGFGRIGRLV` | Human GAPDH (housekeeping protein) | 0.670–0.721 | Toxin |
| `GIINTLQKYYCRVRG` | β-defensin-3 — **this study's own adjuvant** | 0.781 | Toxin |

That is a 100% false-positive rate on benign controls, including the
construct's own adjuvant. On the real candidate set the composition score
flagged 92.5% of peptides as toxins at threshold 0.5 (65.1% at 0.6).

By contrast, across all 1,025 real candidate peptides the two
**evidence-based** channels returned:
- MERCI toxin-motif hits: **0**
- ToxinPred internal BLAST hits to known toxins: **0**
- Independent Tox-Prot BLASTP 3-of-3 hits: **0**

Three independent homology/motif tests agree these viral epitopes bear no
resemblance to known toxins; only the miscalibrated composition heuristic
disagrees.

**Decision rule adopted:** a peptide is flagged TOXIC by ToxinPred only when
there is real toxin evidence (a BLAST homology hit or a MERCI motif hit).
The composition-only ML score is still computed and reported in every output
row (`ToxinPred_ML_Score`, with `COMPOSITION-FLAG` in the `ToxinPred` column
when ML >= 0.6), but is advisory and does not exclude a candidate. All other
toxicity criteria in Section II.C.I.E are applied unchanged as hard filters.

This is a documented, evidence-based deviation from a literal reading of the
methods section, in the same spirit as the Phase 1Dc conservancy note, and
should be disclosed in the manuscript rather than left implicit.

If those binaries/database are not found, this script falls back to a manual
workflow: it emits a query FASTA and a `Manual_Toxicity_Results.csv`
template, and a human runs the real web tools
(https://webs.iiitd.edu.in/raghava/toxinpred2/, HemoPI2, and a BLASTP search
against Tox-Prot) and pastes the results back in before re-running.

## Handling of incomplete data
A peptide is only classified NON-TOXIC if every required field is resolved.
Any peptide with a blank/unfilled manual field (fallback mode only) is
classified UNRESOLVED and routed to `Needs_Review`, never silently counted
as passing.

## BLASTP 3-of-3 rule
A peptide is TOXIC via Tox-Prot homology only if its best hit satisfies ALL
THREE: percent identity >= 80%, alignment coverage (alignment_length /
peptide_length) >= 80%, AND E-value <= 1e-5 -- not E-value alone.
"""


def run_step1ea_toxicity():
    start_time = time.time()

    input_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1Dc", "Filtered_Benchmarks", CONSERVANCY_TIER_DIR)
    output_base = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1E", "Phase1Ea")
    raw_dir = os.path.join(output_base, "Raw")
    filt_dir = os.path.join(output_base, "Filtered")
    review_dir = os.path.join(output_base, "Needs_Review")
    prio_dir = os.path.join(output_base, "Prioritized")
    tool_runs_dir = os.path.join(output_base, "_tool_runs")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(filt_dir, exist_ok=True)
    os.makedirs(review_dir, exist_ok=True)
    os.makedirs(tool_runs_dir, exist_ok=True)

    print("\n" + "="*80 + "\nPHASE 1Ea: TOXICITY SCREENING\n" + "="*80)
    print(f"[INFO] Conservancy tier    : {CONSERVANCY_TIER_DIR} (>50%, \"highly acceptable\" per methodology)")

    if not os.path.isdir(input_folder):
        print(f"[ERROR] Conservancy tier directory not found at: {input_folder}")
        print("[ERROR] Run Phase 1Dc first.")
        return
    csv_files = [f for f in os.listdir(input_folder) if f.endswith(".csv")]
    if not csv_files:
        print("[ERROR] No input CSV found.")
        return
    latest_csv = os.path.join(input_folder, sorted(csv_files)[-1])

    with open(latest_csv, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        original_fields = reader.fieldnames

    unique_peptides = sorted(set(row['Peptide'] for row in rows))
    print(f"[INFO] {len(rows)} candidate rows | {len(unique_peptides)} unique peptides")

    tox_predictions, hemo_scores, blast_hits = {}, {}, {}
    use_local_tools = local_tools_available()

    if use_local_tools:
        print("[INFO] Local ToxinPred2/HemoPI2/BLASTP detected -- running automated screen.")
        query_fasta = os.path.join(tool_runs_dir, "toxicity_query.fasta")
        pep_id_map = {}
        with open(query_fasta, "w") as f:
            for i, pep in enumerate(unique_peptides):
                pep_id = f"Pep_{i}"
                pep_id_map[pep_id] = pep
                f.write(f">{pep_id}\n{pep}\n")

        tp_raw = run_toxinpred2(query_fasta, tool_runs_dir)
        hemo_raw = run_hemopi2(query_fasta, tool_runs_dir)
        blast_raw = run_blastp_toxprot(query_fasta, tool_runs_dir)

        tox_predictions = {pep_id_map[k]: v for k, v in tp_raw.items() if k in pep_id_map}
        hemo_scores = {pep_id_map[k]: v for k, v in hemo_raw.items() if k in pep_id_map}
        blast_hits = {pep_id_map[k]: v for k, v in blast_raw.items() if k in pep_id_map}
    else:
        manual_results_file = os.path.join(output_base, "Manual_Toxicity_Results.csv")
        fasta_export = os.path.join(output_base, "Toxicity_Query.fasta")
        if not os.path.exists(manual_results_file):
            with open(fasta_export, "w") as f_fasta, open(manual_results_file, "w", newline='') as f_csv:
                writer = csv.writer(f_csv)
                writer.writerow(["Peptide", "ToxinPred_Result", "Hemolysis_Prob", "Pident", "Coverage_Pct", "Evalue"])
                for i, pep in enumerate(unique_peptides):
                    f_fasta.write(f">Pep_{i}\n{pep}\n")
                    writer.writerow([pep, "", "", "", "", ""])
            print(f"\n[ACTION REQUIRED] Local ToxinPred2/HemoPI2/BLASTP binaries or Tox-Prot DB not found.")
            print(f"1. FASTA export ready at: {os.path.basename(fasta_export)}")
            print(f"2. Query ToxinPred2, HemoPI2 (https://webs.iiitd.edu.in/raghava/) and BLASTP vs UniProt Tox-Prot manually.")
            print(f"3. Fill in every column of: {os.path.basename(manual_results_file)}")
            print(f"4. Re-run this script.\n")
            sys.exit(0)

        with open(manual_results_file, 'r') as f:
            for row in csv.DictReader(f):
                pep = row['Peptide']
                tp_raw = row.get('ToxinPred_Result', '').strip().upper()
                if tp_raw:
                    # Manual fallback has no separable evidence channels; a
                    # human-entered "TOXIN" is taken as an evidence-backed call.
                    tox_predictions[pep] = {
                        "ml": 1.0 if tp_raw == "TOXIN" else 0.0,
                        "blast": 0.5 if tp_raw == "TOXIN" else 0.0,
                        "merci": 0.0,
                    }
                hemo = row.get('Hemolysis_Prob', '').strip()
                if hemo:
                    try: hemo_scores[pep] = float(hemo)
                    except ValueError: pass
                pident, cov, ev = row.get('Pident', '').strip(), row.get('Coverage_Pct', '').strip(), row.get('Evalue', '').strip()
                if pident and cov and ev:
                    try:
                        blast_hits[pep] = {"pident": float(pident), "coverage": float(cov), "evalue": float(ev)}
                    except ValueError:
                        pass

    # ---------------------------------------------------------
    # Apply the toxicity matrix
    # ---------------------------------------------------------
    fieldnames = original_fields + ["Hydro_Fraction", "Cys_Count", "ToxinPred",
                                     "ToxinPred_ML_Score", "ToxinPred_BLAST", "ToxinPred_MERCI", "Hemolysis_Prob",
                                     "ToxProt_Pident", "ToxProt_Coverage", "ToxProt_Evalue",
                                     "Toxicity_Status", "Exclusion_Reason"]
    raw_data, filtered_data, review_data = [], [], []

    for i, row in enumerate(rows):
        pep = row['Peptide']
        hydro_ratio = sum(pep.count(aa) for aa in "AVILMFWY") / len(pep)
        c_count = pep.count('C')

        unresolved_reasons = []

        # ToxinPred: a TOXIC call requires actual toxin EVIDENCE (BLAST
        # homology to a known toxin, or a MERCI toxin motif). The
        # composition-only ML score is recorded and reported but is NOT a
        # hard exclusion -- it is demonstrably miscalibrated at this peptide
        # length (see METHODOLOGY_NOTE_TEXT for the control-experiment data).
        tp = tox_predictions.get(pep)
        if tp is None:
            unresolved_reasons.append("ToxinPred result missing")
            tp_ml, tp_blast, tp_merci = "", "", ""
            fail_toxinpred = False
            toxinpred_result = "UNKNOWN"
        else:
            tp_ml, tp_blast, tp_merci = tp["ml"], tp["blast"], tp["merci"]
            has_evidence = (tp_blast > 0) or (tp_merci > 0)
            fail_toxinpred = has_evidence
            if has_evidence:
                toxinpred_result = "TOXIN (homology/motif evidence)"
            elif tp_ml >= 0.6:
                toxinpred_result = "COMPOSITION-FLAG (advisory, no toxin evidence)"
            else:
                toxinpred_result = "NON-TOXIN"

        hemo_prob = hemo_scores.get(pep)
        if hemo_prob is None:
            unresolved_reasons.append("Hemolysis_Prob not resolved")
            fail_hemo = False
        else:
            fail_hemo = (hemo_prob >= 0.50) and (hydro_ratio > 0.60)

        hit = blast_hits.get(pep)
        if hit is None:
            # No significant Tox-Prot hit at all is a genuine, resolved "not toxic
            # by homology" result -- not the same as "unresolved" (unlike a blank
            # manual field). Only the manual-fallback path can produce a true
            # unresolved BLASTP state, handled by blast_hits simply omitting pep.
            fail_blast = False
            pident_val, cov_val, evalue_val = "", "", ""
        else:
            fail_blast = (hit["pident"] >= 80.0) and (hit["coverage"] >= 80.0) and (hit["evalue"] <= 1e-5)
            pident_val, cov_val, evalue_val = round(hit["pident"], 2), round(hit["coverage"], 2), hit["evalue"]

        fail_hydro = hydro_ratio > 0.80
        fail_cys = c_count > 2
        is_toxic = fail_hydro or fail_cys or fail_toxinpred or fail_hemo or fail_blast
        is_unresolved = len(unresolved_reasons) > 0

        clean_row = {k: row[k] for k in original_fields}
        clean_row.update({
            "Hydro_Fraction": round(hydro_ratio, 2), "Cys_Count": c_count,
            "ToxinPred": toxinpred_result,
            "ToxinPred_ML_Score": tp_ml, "ToxinPred_BLAST": tp_blast, "ToxinPred_MERCI": tp_merci,
            "Hemolysis_Prob": hemo_prob if hemo_prob is not None else "",
            "ToxProt_Pident": pident_val, "ToxProt_Coverage": cov_val, "ToxProt_Evalue": evalue_val,
        })

        if is_unresolved:
            clean_row["Toxicity_Status"] = "UNRESOLVED"
            clean_row["Exclusion_Reason"] = "; ".join(unresolved_reasons)
            review_data.append(clean_row)
        elif is_toxic:
            clean_row["Toxicity_Status"] = "TOXIC"
            reasons = []
            if fail_hydro: reasons.append("Hydrophobic>80%")
            if fail_cys: reasons.append("Cys>2")
            if fail_toxinpred: reasons.append("ToxinPred toxin homology/motif")
            if fail_hemo: reasons.append("Hemolytic")
            if fail_blast: reasons.append("Tox-Prot homology (3-of-3)")
            clean_row["Exclusion_Reason"] = "; ".join(reasons)
        else:
            clean_row["Toxicity_Status"] = "NON-TOXIC"
            clean_row["Exclusion_Reason"] = "N/A"
            filtered_data.append(clean_row)

        raw_data.append(clean_row)
        sys.stdout.write(f"\r[ PROCESS ] {i+1:03d}/{len(rows):03d} | Filtering Toxicity Matrix")
        sys.stdout.flush()

    print(f"\n[INFO] Screening complete. NON-TOXIC survivors: {len(filtered_data)}/{len(rows)} | Needs review: {len(review_data)}")
    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    if raw_data:
        with open(os.path.join(raw_dir, f"Phase1Ea_Raw_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(raw_data)
    if filtered_data:
        with open(os.path.join(filt_dir, f"Phase1Ea_Filtered_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(filtered_data)

        # ------------------------------------------------------------------
        # PRIORITIZED SUBSET for the manual allergenicity stage (Phase 1Eb).
        #
        # AllerTOP and AllergenFP have no API and must be queried by hand, one
        # peptide at a time, so the full NON-TOXIC pool is not manually
        # screenable. Rather than let that practical limit silently truncate
        # the candidate set, select the top TOP_N_PER_TARGET_PER_CLASS
        # candidates for EACH (target x epitope-class) pair. This caps manual
        # workload while structurally guaranteeing that every antigen and all
        # three epitope classes (MHC-I / MHC-II / B-cell) reach the construct
        # -- the previous run collapsed to B-cell-only from 3 of 7 antigens.
        #
        # Ranking: conservancy desc (primary; the paper's own prioritization
        # metric), then percentile rank asc (stronger binder first), then
        # peptide asc for deterministic tie-breaking across reruns.
        # ------------------------------------------------------------------
        def _rank_key(r):
            try:
                cons = -float(r.get("Conservancy") or 0.0)
            except ValueError:
                cons = 0.0
            try:
                rank = float(r.get("Percentile_Rank") or 999.0)
            except ValueError:
                rank = 999.0
            return (cons, rank, r.get("Peptide", ""))

        groups = {}
        for r in filtered_data:
            groups.setdefault((r.get("Target", ""), r.get("Type", "")), []).append(r)

        prioritized = []
        for key in sorted(groups):
            picked = sorted(groups[key], key=_rank_key)[:TOP_N_PER_TARGET_PER_CLASS]
            prioritized.extend(picked)

        os.makedirs(prio_dir, exist_ok=True)
        with open(os.path.join(prio_dir, f"Phase1Ea_Prioritized_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(prioritized)

        n_unique_prio = len({r["Peptide"] for r in prioritized})
        print(f"[INFO] Prioritized for manual screening: {len(prioritized)} rows "
              f"({n_unique_prio} unique peptides) "
              f"= top {TOP_N_PER_TARGET_PER_CLASS} per target x class across {len(groups)} groups")
        for key in sorted(groups):
            n_picked = min(len(groups[key]), TOP_N_PER_TARGET_PER_CLASS)
            print(f"          {key[0]:12s} {key[1]:8s} : {n_picked} of {len(groups[key])} available")

    if review_data:
        with open(os.path.join(review_dir, f"Phase1Ea_NeedsReview_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(review_data)
        print(f"[ACTION REQUIRED] {len(review_data)} peptides need manual completion -- see {review_dir}")

    note_path = os.path.join(output_base, "METHODOLOGY_NOTE.md")
    with open(note_path, "w") as f:
        f.write(METHODOLOGY_NOTE_TEXT)
    print(f"[INFO] Methodology note written to {note_path}")

if __name__ == "__main__":
    run_step1ea_toxicity()
