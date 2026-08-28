import os, sys, csv, shutil, subprocess
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
# Sec. I.E (current manuscript): "all candidate T-cell and B-cell epitopes
# surviving toxicity and allergenicity screening were additionally screened
# via BLASTP against the reviewed human reference proteome (UniProt Homo
# sapiens, Swiss-Prot subset)." This runs AFTER 1Eb (allergenicity), which is
# why it is Phase1Ec rather than folded into 1Ea.
#
# Exclude ("self-like") only when ALL THREE hold: pident >= 70%,
# coverage (alignment_length/peptide_length) >= 70%, E-value <= 1e-5.
# Anything with a hit below that bar is RETAINED but flagged for manual
# review -- the paper is explicit that partial homology should not be a
# hard exclusion ("even sub-threshold similarity to self-proteins can, in
# some cases, still influence immunogenicity through partial tolerance
# effects").
# =============================================================================
BLASTP_BINARY = os.environ.get("BLASTP_BINARY", "/opt/miniconda3/envs/phase2/bin/blastp")
HUMAN_SWISSPROT_DB = os.environ.get(
    "HUMAN_SWISSPROT_DB",
    os.path.join(_PROJECT_ROOT, "human_swissprot_db", "human_swissprot"),
)

PIDENT_CUTOFF = 70.0
COVERAGE_CUTOFF = 70.0
EVALUE_CUTOFF = 1e-5


def _binary_exists(path):
    return shutil.which(path) is not None or os.path.isfile(path)


def local_tools_available():
    return _binary_exists(BLASTP_BINARY) and (
        os.path.isfile(HUMAN_SWISSPROT_DB + ".phr") or os.path.isfile(HUMAN_SWISSPROT_DB + ".pin")
    )


def run_blastp_human_swissprot(fasta_path, work_dir):
    """
    Best (lowest-E-value) reviewed-human hit per query peptide. Same
    outfmt/columns as Phase1Ea's run_blastp_toxprot -- only the target DB
    and the acceptance thresholds (70/70/1e-5 vs Tox-Prot's 80/80/1e-5)
    differ, per the paper's own distinct wording for each screen.
    """
    output_tsv = os.path.join(work_dir, "blastp_human_swissprot.tsv")
    cmd = [
        BLASTP_BINARY, "-query", fasta_path, "-db", HUMAN_SWISSPROT_DB,
        "-outfmt", "6 qseqid sseqid evalue pident length qlen",
        "-evalue", "1e-3", "-max_target_seqs", "1", "-out", output_tsv,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
    if result.returncode != 0:
        raise RuntimeError(f"BLASTP (human Swiss-Prot) failed (exit {result.returncode}):\n{result.stderr}")

    best_hits = {}
    if os.path.isfile(output_tsv):
        with open(output_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 6:
                    continue
                qseqid, sseqid, evalue, pident, length, qlen = (
                    parts[0], parts[1], float(parts[2]), float(parts[3]), int(parts[4]), int(parts[5]),
                )
                if qseqid not in best_hits or evalue < best_hits[qseqid]["evalue"]:
                    best_hits[qseqid] = {
                        "evalue": evalue, "pident": pident,
                        "coverage": (length / qlen) * 100 if qlen else 0.0,
                        "subject": sseqid,
                    }
    return best_hits


METHODOLOGY_NOTE_TEXT = """# Phase 1Ec — Human Self-Homology Screening Methodology Note

**Applies to:** Sec. I.E of the manuscript (BLASTP against the reviewed human
reference proteome, UniProt Homo sapiens Swiss-Prot subset).

## What this script does
Runs AFTER Phase 1Eb (allergenicity), matching the paper's stated order
("all candidate T-cell and B-cell epitopes surviving toxicity and
allergenicity screening were additionally screened..."). Every peptide in
1Eb's Filtered output is BLASTP'd against a local database built from
`reviewed:true AND organism_id:9606` (UniProt REST, Swiss-Prot subset only
-- not the much larger unreviewed TrEMBL set).

## Decision rule
A peptide is flagged "self-like" and EXCLUDED only if its best hit meets
ALL THREE: percent identity >= 70%, alignment coverage
(alignment_length / peptide_length) >= 70%, E-value <= 1e-5.

Any hit below that bar is RETAINED but flagged `Self_Homology_Status =
PARTIAL` for manual review, per the paper's own text: "Epitopes with
partial homology below this threshold were retained but flagged for
manual review, given that even sub-threshold similarity to self-proteins
can, in some cases, still influence immunogenicity through partial
tolerance effects."

## Database provenance
`human_swissprot_db/human_swissprot.fasta`, fetched from
`rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&format=fasta`
-- 20,431 reviewed human sequences, indexed with `makeblastdb -dbtype prot`.
"""


def run_step1ec_self_homology():
    input_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1E", "Phase1Eb", "Filtered")
    output_base = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1E", "Phase1Ec")
    raw_dir = os.path.join(output_base, "Raw")
    filt_dir = os.path.join(output_base, "Filtered")
    tool_runs_dir = os.path.join(output_base, "_tool_runs")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(filt_dir, exist_ok=True)
    os.makedirs(tool_runs_dir, exist_ok=True)

    common.print_banner("PHASE 1Ec: HUMAN SELF-HOMOLOGY SCREENING (BLASTP vs Swiss-Prot Human)")

    latest_csv = common.latest_file(input_folder, suffix=".csv")
    if latest_csv is None:
        print(f"[ERROR] No Phase 1Eb Filtered output found at: {input_folder}")
        print("[ERROR] Run Phase 1Eb first.")
        return

    with open(latest_csv, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        original_fields = reader.fieldnames

    unique_peptides = sorted(set(row['Peptide'] for row in rows))
    print(f"[INFO] {len(rows)} candidate rows | {len(unique_peptides)} unique peptides")

    if not local_tools_available():
        print(f"[ERROR] blastp or the human Swiss-Prot DB not found.")
        print(f"[ERROR] BLASTP_BINARY = {BLASTP_BINARY}")
        print(f"[ERROR] HUMAN_SWISSPROT_DB = {HUMAN_SWISSPROT_DB}")
        return

    query_fasta = os.path.join(tool_runs_dir, "self_homology_query.fasta")
    pep_id_map = {}
    with open(query_fasta, "w") as f:
        for i, pep in enumerate(unique_peptides):
            pep_id = f"Pep_{i}"
            pep_id_map[pep_id] = pep
            f.write(f">{pep_id}\n{pep}\n")

    print("[INFO] Running BLASTP against local human Swiss-Prot database...")
    blast_raw = run_blastp_human_swissprot(query_fasta, tool_runs_dir)
    blast_hits = {pep_id_map[k]: v for k, v in blast_raw.items() if k in pep_id_map}

    fieldnames = original_fields + [
        "SelfHomology_Pident", "SelfHomology_Coverage", "SelfHomology_Evalue",
        "SelfHomology_Subject", "Self_Homology_Status",
    ]
    raw_data, filtered_data = [], []
    n_excluded, n_partial, n_clean = 0, 0, 0

    for row in rows:
        pep = row['Peptide']
        hit = blast_hits.get(pep)
        clean_row = {k: row[k] for k in original_fields}

        if hit is None:
            clean_row.update({
                "SelfHomology_Pident": "", "SelfHomology_Coverage": "",
                "SelfHomology_Evalue": "", "SelfHomology_Subject": "",
                "Self_Homology_Status": "NONE",
            })
            n_clean += 1
            filtered_data.append(clean_row)
        else:
            is_self_like = (
                hit["pident"] >= PIDENT_CUTOFF
                and hit["coverage"] >= COVERAGE_CUTOFF
                and hit["evalue"] <= EVALUE_CUTOFF
            )
            clean_row.update({
                "SelfHomology_Pident": round(hit["pident"], 2),
                "SelfHomology_Coverage": round(hit["coverage"], 2),
                "SelfHomology_Evalue": hit["evalue"],
                "SelfHomology_Subject": hit["subject"],
                "Self_Homology_Status": "SELF-LIKE" if is_self_like else "PARTIAL",
            })
            if is_self_like:
                n_excluded += 1
            else:
                n_partial += 1
                filtered_data.append(clean_row)

        raw_data.append(clean_row)

    print(f"[INFO] Self-like (excluded): {n_excluded} | Partial (retained+flagged): {n_partial} | "
          f"No hit (clean): {n_clean}")
    print(f"[INFO] Retained for Phase 1F: {len(filtered_data)}/{len(rows)} rows")

    # Per (target x class) survivor counts -- surfaces any group emptied by
    # this screen, per the plan's I-1 contingency (default: accept and
    # document; only re-screen deeper if a whole antigen loses all 3 classes).
    groups = {}
    for r in filtered_data:
        groups.setdefault((r.get("Target", ""), r.get("Type", "")), 0)
        groups[(r.get("Target", ""), r.get("Type", ""))] += 1
    for key in sorted(groups):
        print(f"          {key[0]:12s} {key[1]:8s} : {groups[key]} survivors")
    empty_targets = set()
    for row in rows:
        key_target = row.get("Target", "")
        empty_targets.add(key_target)
    for t in sorted(empty_targets):
        classes_present = {k[1] for k in groups if k[0] == t}
        if len(classes_present) < 3:
            missing = {"MHC-I", "MHC-II", "B-cell"} - classes_present
            if missing:
                print(f"[WARNING] {t}: no survivors for class(es) {sorted(missing)} after self-homology screen")

    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    if raw_data:
        with open(os.path.join(raw_dir, f"Phase1Ec_Raw_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(raw_data)
    if filtered_data:
        with open(os.path.join(filt_dir, f"Phase1Ec_Filtered_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(filtered_data)

    note_path = os.path.join(output_base, "METHODOLOGY_NOTE.md")
    with open(note_path, "w") as f:
        f.write(METHODOLOGY_NOTE_TEXT)
    print(f"[INFO] Methodology note written to {note_path}")


if __name__ == "__main__":
    run_step1ec_self_homology()
