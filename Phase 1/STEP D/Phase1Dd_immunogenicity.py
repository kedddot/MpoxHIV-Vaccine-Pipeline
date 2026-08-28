import os, sys, csv, subprocess, tempfile
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
# Sec. I.D (current manuscript): "MHC Class I candidates surviving the
# top-1% binding-percentile threshold were additionally screened using
# BigMHC...to estimate T-cell recognition likelihood independent of binding
# affinity alone. Peptides receiving a high BigMHC immunogenicity score were
# retained as prioritized candidates; peptides with a negative score were
# deprioritized but not automatically excluded." Applied ONLY to MHC-I
# (the IEDB Immunogenicity tool's MHC-II counterpart does not exist, and
# the manuscript restricts this step to MHC-I on that basis too).
#
# PRIORITIZATION ONLY -- this script must never remove a peptide.
# =============================================================================
BIGMHC_ROOT = os.environ.get("BIGMHC_ROOT", os.path.join(_PROJECT_ROOT, "external_tools", "bigmhc"))
BIGMHC_PYTHON = os.environ.get("BIGMHC_PYTHON", "/opt/miniconda3/envs/phase2/bin/python")

# BigMHC IM's own output is a bounded [0,1] probability (a sigmoid output),
# not a raw pre-activation score -- so "negative score" from the manuscript
# is read as "below the tool's own neutral midpoint", the natural cutoff for
# a bounded probability. Documented explicitly rather than left implicit.
PRIORITY_CUTOFF = 0.5


def bigmhc_available():
    return (
        os.path.isfile(os.path.join(BIGMHC_ROOT, "src", "predict.py"))
        and os.path.isdir(os.path.join(BIGMHC_ROOT, "models", "bat512", "im"))
    )


def run_bigmhc_im(pairs, work_dir):
    """
    pairs: list of (allele, peptide) tuples.
    Returns {(allele, peptide): score}. Runs BigMHC's 7-model IM ensemble
    on CPU (no GPU required per the tool's own docs for inference/transfer
    learning; only large-batch training needs one).
    """
    input_csv = os.path.join(work_dir, "bigmhc_im_input.csv")
    output_csv = os.path.join(work_dir, "bigmhc_im_output.csv")
    with open(input_csv, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["mhc", "pep", "tgt"])
        for allele, pep in pairs:
            writer.writerow([allele, pep, 0])

    cmd = [
        BIGMHC_PYTHON, os.path.join(BIGMHC_ROOT, "src", "predict.py"),
        f"-i={input_csv}", "-m=im", "-t=2", "-a=0", "-p=1", "-c=1",
        "-d=cpu", f"-o={output_csv}",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.join(BIGMHC_ROOT, "src"))
    if result.returncode != 0 or not os.path.isfile(output_csv):
        raise RuntimeError(f"BigMHC IM failed (exit {result.returncode}):\n{result.stderr[-2000:]}")

    scores = {}
    with open(output_csv) as f:
        for row in csv.DictReader(f):
            scores[(row["mhc"], row["pep"])] = float(row["BigMHC_IM"])
    return scores


def run_step1dd_immunogenicity():
    common.print_banner("PHASE 1Dd: BigMHC IMMUNOGENICITY PRIORITIZATION (MHC-I ONLY)")

    input_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1Db")
    tool_runs_dir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1Dd", "_tool_runs")
    os.makedirs(tool_runs_dir, exist_ok=True)

    latest_db = common.latest_file(input_folder, suffix=".csv")
    if latest_db is None:
        print(f"[ERROR] No Phase 1Db output found at: {input_folder}")
        return

    with open(latest_db, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        original_fields = reader.fieldnames

    if not bigmhc_available():
        print(f"[ERROR] BigMHC not found at {BIGMHC_ROOT} (expected src/predict.py + models/bat512/im/).")
        print("[ERROR] Install: git clone https://github.com/KarchinLab/bigmhc.git into external_tools/")
        return

    # Build the (allele, peptide) query set from every MHC-I row's
    # Binding_Alleles column (already the top-1%-rank survivors from 1Db).
    pairs = set()
    for r in rows:
        if r.get("Type") != "MHC-I":
            continue
        for allele in (r.get("Binding_Alleles") or "").split(";"):
            allele = allele.strip()
            if allele:
                pairs.add((allele, r["Peptide"]))
    pairs = sorted(pairs)
    print(f"[INFO] {len(pairs)} unique (allele, peptide) MHC-I pairs to score.")

    if not pairs:
        print("[WARNING] No MHC-I (allele, peptide) pairs found -- nothing to score.")
        scores = {}
    else:
        print("[INFO] Running BigMHC IM (7-model ensemble, CPU)...")
        scores = run_bigmhc_im(pairs, tool_runs_dir)
        print(f"[SUCCESS] BigMHC IM scored {len(scores)} pairs.")

    fieldnames = original_fields + ["BigMHC_IM_Score", "Immunogenicity_Priority"]
    out_rows = []
    n_prio = n_deprio = n_na = 0
    for r in rows:
        clean = {k: r[k] for k in original_fields}
        if r.get("Type") == "MHC-I":
            alleles = [a.strip() for a in (r.get("Binding_Alleles") or "").split(";") if a.strip()]
            per_allele_scores = [scores[(a, r["Peptide"])] for a in alleles if (a, r["Peptide"]) in scores]
            if per_allele_scores:
                best = max(per_allele_scores)
                clean["BigMHC_IM_Score"] = round(best, 6)
                clean["Immunogenicity_Priority"] = "PRIORITIZED" if best >= PRIORITY_CUTOFF else "DEPRIORITIZED"
                n_prio += best >= PRIORITY_CUTOFF
                n_deprio += best < PRIORITY_CUTOFF
            else:
                clean["BigMHC_IM_Score"] = ""
                clean["Immunogenicity_Priority"] = "NOT_SCORED"
                n_na += 1
        else:
            clean["BigMHC_IM_Score"] = ""
            clean["Immunogenicity_Priority"] = ""
        out_rows.append(clean)

    ts = datetime.now().strftime("%Y%m%d_%H%M")  # match Phase1Db/1Dc/1De's format -- see 1De's note
    out_path = os.path.join(input_folder, f"Phase1Db_Elite_{ts}.csv")
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader(); writer.writerows(out_rows)

    print(f"\n[INFO] MHC-I immunogenicity -- PRIORITIZED: {n_prio} | DEPRIORITIZED: {n_deprio} | NOT_SCORED: {n_na}")
    print("[INFO] Deprioritized peptides are NOT excluded -- retained per the manuscript's own wording.")
    print(f"[INFO] Enriched Phase1Db_Elite CSV written to: {out_path}")
    print("=" * 90 + "\n")


if __name__ == "__main__":
    run_step1dd_immunogenicity()
