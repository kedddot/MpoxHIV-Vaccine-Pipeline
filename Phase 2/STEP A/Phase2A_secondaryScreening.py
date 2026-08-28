import os
import sys
import csv
import json
import time
import shutil
import subprocess
import hashlib
from datetime import datetime
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# =============================================================================
# MINIMAL BOOTSTRAP -- locates the shared phase2_common module.
# See phase2_common.py for why this logic is centralized.
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
_COMMON_DIR = os.path.join(_PROJECT_ROOT, "Phase 2", "_common")
if _COMMON_DIR not in sys.path:
    sys.path.insert(0, _COMMON_DIR)

import phase2_common as common

# =============================================================================
# CONFIGURE -- fill this in for your machine before running
# =============================================================================
TOXPROT_BLAST_DB = os.environ.get("TOXPROT_BLAST_DB", os.path.join(_PROJECT_ROOT, "toxprot_db", "toxprot"))
BLASTP_BINARY = os.environ.get("BLASTP_BINARY", "/opt/miniconda3/envs/phase2/bin/blastp")
TOXINPRED2_BINARY = os.environ.get("TOXINPRED2_BINARY", "/opt/miniconda3/envs/phase2/bin/toxinpred2")
HEMOPI2_BINARY = os.environ.get("HEMOPI2_BINARY", "/opt/miniconda3/envs/phase2/bin/hemopi2_classification")

# =============================================================================
# SCOPE NOTE (read before running)
#
# TOXICITY and ALLERGENICITY are covered with real tools. SOLUBILITY/
# STABILITY (DeepSol replacing SOLpro, CamSol sequence-mode) is now
# covered too, via the same manual web-submission pattern as AllerTOP/
# AllergenFP -- all four have no API and scale with candidate count, so
# they're combined into one manual results file/prompt below.
# =============================================================================

# =============================================================================
# REAL TOOLS -- ToxinPred2, HemoPI2 (batch: run once across ALL candidates
# or ALL junction peptides at once, rather than one call per sequence)
# =============================================================================

def run_toxinpred2_batch(fasta_path, work_dir):
    """
    Real ToxinPred2 HYBRID (-m 2: RF + BLAST + MERCI) batch run.

    Returns {name: {"ml": float, "blast": float, "merci": float}} -- the
    three evidence channels kept SEPARATE, rather than the tool's own
    collapsed Toxin/Non-Toxin verdict, so the caller can weigh
    composition evidence apart from homology/motif evidence.

    WHY -m 2 AND WHY SEPARATE CHANNELS (matches Phase 1E, deviation #1):
    ToxinPred2's amino-acid-composition RF is not calibrated for short
    peptides -- Phase 1 measured it calling human albumin, human GAPDH, a
    bare linker, and this study's own beta-defensin-3 adjuvant all
    "Toxin" on composition alone. Phase 1 therefore adopted the rule
    "flag TOXIC only on real toxin evidence (a BLAST homology hit or a
    MERCI motif hit)"; the composition score stays advisory. Phase 2A's
    junction check previously used -m 1 (composition only) and thresholded
    the raw score, which is exactly the miscalibrated path Phase 1 had
    already rejected -- it produced 745 "toxic" junction peptides whose
    BLAST and MERCI channels were both zero.

    -m 2 also needs three upstream bugs patched in the installed package
    (all applied; re-apply if the conda env is ever rebuilt):
      1. BLAST_processor()'s else-branch iterated an out-of-scope `seqid`
         instead of its own `name1` parameter (NameError).
      2. hybrid()'s df8.sum(axis=1) summed a string column on modern
         pandas -- needs numeric_only=True.
      3. BLAST_processor() used DataFrame.append(), removed in pandas 2.x.
    """
    # Model is part of the cache key: -m 1 and -m 2 produce different
    # columns and different verdicts, so a cached -m 1 file must never be
    # silently reused for an -m 2 run (or vice versa).
    output_csv = os.path.join(work_dir, f"toxinpred2_m2_{os.path.basename(fasta_path)}.csv")
    if shutil.which(TOXINPRED2_BINARY) is None and not os.path.isfile(TOXINPRED2_BINARY):
        raise FileNotFoundError(
            f"toxinpred2 binary not found ('{TOXINPRED2_BINARY}'). Set "
            f"TOXINPRED2_BINARY at the top of this file to the full path "
            f"(e.g. from `which toxinpred2` with your environment active) -- "
            f"calling it by bare name only works if that environment's PATH "
            f"is active in whatever process launches this script."
        )
    cmd = [
        TOXINPRED2_BINARY, "-i", fasta_path, "-o", output_csv,
        "-t", "0.6",  # tool's own default hybrid threshold; only affects the
                      # tool's collapsed verdict, which we deliberately ignore
                      # in favour of the separate channels below
        "-m", "2",    # hybrid RF+BLAST+MERCI -- see docstring
        "-d", "2",    # Display ALL peptides -- default (1) only reports predicted toxins,
                      # which would silently drop every non-toxic candidate from the output
    ]
    if os.path.isfile(output_csv) and os.path.getsize(output_csv) > 0:
        print(f"[INFO] Reusing cached ToxinPred2 result: {output_csv}")
    else:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
        if result.returncode != 0:
            raise RuntimeError(f"ToxinPred2 failed (exit {result.returncode}):\n{result.stderr}")
        if not os.path.isfile(output_csv):
            raise FileNotFoundError(f"ToxinPred2 reported success but no output file at {output_csv}")
        # ToxinPred2 writes a stray "Sequence_1" intermediate file into its
        # cwd regardless of -o -- harmless but clutters _tool_runs/, so
        # clean it up right after a fresh run (nothing to clean on a
        # cache-hit above, since that file was already handled last time).
        stray_file = os.path.join(work_dir, "Sequence_1")
        if os.path.isfile(stray_file):
            os.remove(stray_file)

    predictions = {}
    with open(output_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            # -m 2's own header is: Subject,ML Score,MERCI Score,BLAST Score,
            # Hybrid Score,Prediction. Read each channel by exact name.
            name = row.get("Subject") or row.get("ID") or list(row.values())[0]
            try:
                predictions[name] = {
                    "ml": float(row["ML Score"]),
                    "merci": float(row["MERCI Score"]),
                    "blast": float(row["BLAST Score"]),
                }
            except (KeyError, ValueError, TypeError):
                # Absent/unparseable channels must not read as "no evidence"
                # -- record None so the caller fails closed rather than
                # silently treating a broken row as non-toxic.
                predictions[name] = None
    return predictions


def _toxinpred_has_evidence(entry):
    """
    Phase 1E's decision rule (deviation #1), applied identically here: a
    peptide counts as ToxinPred-toxic ONLY on real toxin evidence -- a
    BLAST homology hit or a MERCI motif hit. The composition-only ML
    score is advisory and never rejects on its own.

    Returns True (toxic), False (no evidence), or None (unresolved -- the
    caller must fail closed, never treat this as False).
    """
    if entry is None:
        return None
    return (entry.get("blast", 0) > 0) or (entry.get("merci", 0) > 0)


def run_hemopi2_batch(fasta_path, work_dir):
    """
    Real HemoPI2 batch run.
    Requires: pip install hemopi2  (also needs torch, transformers, and
    Meta's ESM protein language model -- a heavy dependency chain)
    NOTE: exact CLI flags/output column names are based on the tool's
    published usage pattern -- verify against `hemopi2_classification -h`
    and its actual output header the first time you run this.
    """
    output_filename = f"hemopi2_{os.path.basename(fasta_path)}.csv"
    # HemoPI2 internally builds its output path as f"{wd}/{result_filename}"
    # -- it expects -o to be a bare FILENAME, joined onto -wd itself, not
    # an already-complete path. Passing a full path here (as an earlier
    # version of this function did) doubles the directory and breaks.
    output_csv = os.path.join(work_dir, output_filename)
    if shutil.which(HEMOPI2_BINARY) is None and not os.path.isfile(HEMOPI2_BINARY):
        raise FileNotFoundError(
            f"hemopi2_classification binary not found ('{HEMOPI2_BINARY}'). Set "
            f"HEMOPI2_BINARY at the top of this file to the full path (e.g. "
            f"from `which hemopi2_classification` with your environment active)."
        )
    cmd = [
        HEMOPI2_BINARY,
        "-i", fasta_path,
        "-o", output_filename,
        "-wd", work_dir,
        "-d", "2",  # display all peptides, not just predicted-hemolytic ones (this is
                    # actually already the tool's own default -- kept explicit for clarity)
        "-m", "3",  # ESM2-t6 only (NOT the default model 4, Hybrid2: ESM+MERCI).
                    # Model 4 internally invokes MERCI the same way ToxinPred2's model 2
                    # does -- which just broke on the space in your venv's path. Model 3
                    # still uses the modern ESM protein-language-model scoring (the reason
                    # this tool needed the heavy torch/transformers install), just without
                    # the MERCI motif step that's currently broken by the path issue.
    ]
    if os.path.isfile(output_csv) and os.path.getsize(output_csv) > 0:
        print(f"[INFO] Reusing cached HemoPI2 result: {output_csv}")
    else:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
        if result.returncode != 0:
            raise RuntimeError(f"HemoPI2 failed (exit {result.returncode}):\n{result.stderr}")
        if not os.path.isfile(output_csv):
            raise FileNotFoundError(f"HemoPI2 reported success but no output file at {output_csv}")

    # Paper's rule (Sec. I.E, applied to whole constructs via Sec. IV.A.b)
    # is conjunctive: hemolysis >= 0.50 AND Hydrophobic_Fraction > 0.60.
    # The caller needs the numeric score to apply that -- not just the
    # tool's own collapsed Hemolytic/Non-Hemolytic verdict.
    predictions = {}
    with open(output_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get("SeqID") or row.get("ID") or row.get("Seq_ID") or list(row.values())[0]
            try:
                predictions[name] = {"score": float(row["ESM Score"]), "prediction": row.get("Prediction")}
            except (KeyError, ValueError, TypeError):
                predictions[name] = None
    return predictions


def run_blastp_toxprot_batch(fasta_path, work_dir, evalue_cutoff=1e-5):
    """
    Real BLASTP against a local Tox-Prot database.
    Requires: conda install -c bioconda blast
    Requires: a local Tox-Prot BLAST database. Build one from UniProt's
    Tox-Prot annotation program (https://www.uniprot.org/help/Tox-Prot):
    download the curated Tox-Prot FASTA, then:
        makeblastdb -in toxprot.fasta -dbtype prot -out toxprot
    and point TOXPROT_BLAST_DB (top of this file) at that -out prefix.
    """
    if not os.path.isfile(TOXPROT_BLAST_DB + ".phr") and not os.path.isfile(TOXPROT_BLAST_DB + ".pin"):
        raise FileNotFoundError(
            f"Tox-Prot BLAST database not found at '{TOXPROT_BLAST_DB}'. "
            f"See run_blastp_toxprot_batch()'s docstring for how to build one."
        )
    if shutil.which(BLASTP_BINARY) is None and not os.path.isfile(BLASTP_BINARY):
        raise FileNotFoundError(
            f"blastp binary not found ('{BLASTP_BINARY}'). Export BLASTP_BINARY "
            f"with the full path -- conda environments don't stay active across "
            f"terminal sessions or when this script runs under a different "
            f"interpreter. On this machine: `source ~/mpoxhiv_env.sh`, which sets "
            f"it (and every other tool path) to the SSD-resident phase2 env. "
            f"Otherwise: `conda run -n phase2 which blastp`."
        )

    output_tsv = os.path.join(work_dir, f"blastp_{os.path.basename(fasta_path)}.tsv")
    cmd = [
        BLASTP_BINARY,
        "-query", fasta_path,
        "-db", TOXPROT_BLAST_DB,
        "-outfmt", "6 qseqid sseqid evalue pident",
        "-evalue", "1e-3",
        "-max_target_seqs", "1",
        "-out", output_tsv,
    ]
    if os.path.isfile(output_tsv):
        # blastp only creates -out once the search completes (even a
        # zero-hit result is a legitimate zero-byte file), so existence
        # alone is a safe "already ran" signal -- unlike ToxinPred2/
        # HemoPI2's CSVs, which always have a header row when complete.
        print(f"[INFO] Reusing cached BLASTP result: {output_tsv}")
    else:
        result = _run_with_heartbeat(cmd, cwd=work_dir, label="BLASTP vs Tox-Prot")
        if result.returncode != 0:
            raise RuntimeError(f"BLASTP failed (exit {result.returncode}):\n{result.stderr}")

    hits = {}
    if os.path.isfile(output_tsv):
        with open(output_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                qseqid, sseqid, evalue, pident = parts[0], parts[1], float(parts[2]), float(parts[3])
                if evalue < evalue_cutoff and (qseqid not in hits or evalue < hits[qseqid]["evalue"]):
                    hits[qseqid] = {"significant_hit": True, "evalue": evalue, "subject": sseqid}
    return hits  # candidates NOT present in this dict had no significant hit


# =============================================================================
# DEVIATION #15 -- see Methods_Deviations_RRL_Support.txt for the full record.
#
# Sec. II.A.b states a construct with DeepSol S2 < 0.50 is "predicted
# insoluble" and rejected. When DeepSol is the SOLE failing criterion, this
# pipeline carries the construct forward as REVIEW instead -- an explicit,
# documented deviation, using the PASS/REVIEW/REJECT vocabulary the paper
# itself defines in Sec. II.C.
#
# Basis: 11 constructs were built and measured; DeepSol spanned 0.257-0.490
# and NONE reached 0.50. The best (0.4897) misses by 0.0103. The construct is
# ~21% synthetic linker (12x AAY, 10x GPGPG, 7x KK -> 9.9% Gly, 7.9% Pro vs
# natural ~7%/~5%), far outside DeepSol's training distribution of natural
# E. coli proteins, where its reported accuracy is ~77%. A 0.0103 margin is
# well inside that error. Removing the beta-defensin-3 adjuvant made the
# prediction WORSE (0.300), so this is not an adjuvant artifact.
#
# Set to False to restore the paper's literal reject-on-<0.50 behaviour.
DEEPSOL_ONLY_FAILURE_IS_REVIEW = True

RANK_COLUMN_CANDIDATES = ["percentile_rank", "adjusted_rank", "rank"]


def run_iedb_mhc_binding_batch(peptides, mhc_class, alleles=None, cache_path=None):
    """
    Real IEDB MHC binding prediction via their public API. Batches ALL
    given peptides of a given length into as few API calls as possible
    instead of one call per peptide.
    Requires: pip install requests
    NOTE: IEDB's API endpoint URLs and parameter names have changed
    before -- verify against https://www.iedb.org/tools_api.php (their
    current API docs) if this raises unexpected errors.

    `length` is pinned to each batch's own peptide length rather than
    posted as a comma-list of every possible length. Posting e.g.
    "length=8,9,10,11" against peptides that are already pre-sized makes
    IEDB re-window each input sequence internally at every listed length
    -- for MHC-I that silently produces extra sub-peptides no one asked
    for, and for MHC-II, omitting `length` entirely (the previous
    behaviour) lets netMHCIIpan choose its own internal core peptide,
    which then never key-matches the 12-20mer we actually submitted.
    Pinning `length` to the batch's exact peptide size means IEDB scores
    exactly the peptide we gave it, once, per allele.

    If cache_path is given, a completed result is persisted there and
    reused on a later call instead of re-hitting the network -- this
    result set doesn't change between reruns of the same input.
    """
    import socket
    import requests
    import json
    import urllib3.util.connection as urllib3_cn

    # This environment's IPv6 path to IEDB is a dead end (NAT64 route that
    # never completes its handshake) -- requests/urllib3 tries it first on
    # every single connection and eats ~30s per call before falling back
    # to IPv4, which then succeeds in ~1-2s. Forcing IPv4-only resolution
    # here cuts each of the ~dozens of calls below from ~30s to ~2s.
    urllib3_cn.allowed_gai_family = lambda: socket.AF_INET

    if not peptides:
        return {}

    if cache_path and os.path.isfile(cache_path):
        print(f"[INFO] Reusing cached IEDB {mhc_class} result: {cache_path}")
        with open(cache_path) as f:
            return json.load(f)

    if mhc_class == "MHC_I":
        # https:// directly -- IEDB 308-redirects http:// to https:// on
        # every call, so posting straight to https avoids that extra hop.
        url = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
        default_alleles = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01"]
        method = "netmhcpan_el"
    else:
        url = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
        default_alleles = ["HLA-DRB1*01:01", "HLA-DRB1*07:01", "HLA-DRB1*15:01"]
        method = "netmhciipan_el"

    alleles = alleles or default_alleles

    by_length = {}
    for pep in peptides:
        by_length.setdefault(len(pep), []).append(pep)

    total_batches = len(by_length) * len(alleles)
    batch_num = 0
    best_rank = {}
    for length, peps_of_length in by_length.items():
        sequence_text = "\n".join(f">pep{i}\n{p}" for i, p in enumerate(peps_of_length))
        for allele in alleles:
            batch_num += 1
            batch_start = time.time()
            print(f"[INFO]   IEDB {mhc_class} batch {batch_num}/{total_batches}: "
                  f"{len(peps_of_length)} peptides x len={length}, allele={allele}...", end=" ", flush=True)
            payload = {
                "method": method,
                "sequence_text": sequence_text,
                "allele": allele,
                "length": str(length),
            }

            resp = requests.post(url, data=payload, timeout=120)
            resp.raise_for_status()
            print(f"done ({time.time() - batch_start:.1f}s)")
            lines = resp.text.strip().split("\n")
            if len(lines) < 2:
                continue
            header = lines[0].split("\t")
            if "peptide" not in header:
                continue
            peptide_col = header.index("peptide")
            rank_col = next((header.index(c) for c in RANK_COLUMN_CANDIDATES if c in header), None)
            if rank_col is None:
                continue

            for line in lines[1:]:
                cols = line.split("\t")
                if len(cols) <= max(peptide_col, rank_col):
                    continue
                # Read the peptide IEDB actually reports back for this row
                # -- never assume it echoes what we submitted.
                pep = cols[peptide_col]
                try:
                    rank = float(cols[rank_col])
                except ValueError:
                    continue
                if pep not in best_rank or rank < best_rank[pep]:
                    best_rank[pep] = rank

    if cache_path:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        with open(cache_path, 'w') as f:
            json.dump(best_rank, f)

    return best_rank  # {peptide: best (lowest) percentile rank across alleles tried}

# =============================================================================
# MANUAL -- AllerTOP, AllergenFP (no API or local install exists for either;
# confirmed both remain web-form-only tools). Covers the FULL CONSTRUCT of
# every candidate in one combined file. Junction-level allergen-predictor
# checks are NOT automated here -- see the note in evaluate_candidate().
# =============================================================================

# CamSol Intrinsic's published interpretation thresholds (Sormanni, Aprile &
# Vendruscolo, J Mol Biol 2015): > +1 highly soluble, < -1 poorly soluble.
# EXTERNAL CONVENTION, not a manuscript criterion -- labelled as such wherever
# it is reported, exactly like the Step 2D secondary QC flags.
CAMSOL_POORLY_SOLUBLE = -1.0

# Keys renamed over the project's life. Reading falls back to the old name so
# result files recorded before the rename stay valid without hand-editing --
# the recorded VALUES were always correct, only the key name and the check
# applied to them were wrong (deviation #25).
LEGACY_MANUAL_KEYS = {
    "camsol_intrinsic_solubility_score": "camsol_seq_patch_length",
}


def manual_value(entry, key):
    """Value for `key`, falling back to that key's legacy name."""
    if not entry:
        return None
    value = entry.get(key)
    if value is None and key in LEGACY_MANUAL_KEYS:
        value = entry.get(LEGACY_MANUAL_KEYS[key])
    return value

MANUAL_SEQUENCE_TOOL_KEYS = {
    "allertop_is_allergen": "AllerTOP -- predicted allergen? (true/false)",
    "allergenfp_is_allergen": "AllergenFP -- predicted allergen? (true/false)",
    "deepsol_probability": "DeepSol S2 -- predicted solubility probability (0.0-1.0)",
    # RENAMED from "camsol_seq_patch_length" (deviation #25). That key was
    # documented as a patch LENGTH IN AMINO ACIDS and tested with `> 8`, but
    # CamSol Intrinsic does not report a patch length -- its headline output is
    # the line "The protein variant intrinsic solubility score is <x>", and
    # that score is what every construct in this project has actually been
    # given. A score in CamSol's roughly -3..+3 range can never exceed 8, so
    # the check was structurally incapable of firing. The old key is still
    # accepted on read so historical result files keep working.
    "camsol_intrinsic_solubility_score": "CamSol (sequence/intrinsic mode) -- the "
        "'protein variant intrinsic solubility score' printed on the results page",
}


def run_manual_prepare(variants, manual_results_path):
    common.print_banner("EXTERNAL RESULTS NEEDED -- AllerTOP, AllergenFP, DeepSol, CamSol")
    print("None of these 4 tools have an API or local install, and none accept")
    print("a FASTA file upload -- each takes a single raw sequence pasted")
    print("directly into a text box on its site. Copy each sequence below")
    print("(just the letters, no '>' header) into all 4 sites, one candidate")
    print("at a time, then fill in the template written below.")
    print("-" * 100)
    print("[AllerTOP]     https://www.ddg-pharmfac.net/AllerTOP")
    print("[AllergenFP]   https://ddg-pharmfac.net/AllergenFP")
    print("[DeepSol]      https://machinelearning-protein.qcri.org  (model = 2, i.e. DeepSol S2;")
    print("               needs a free QCAI account -- see earlier notes if not yet created)")
    print("[CamSol]       https://www-cohsoftware.ch.cam.ac.uk/index.php/camsolintrinsic")
    print("               (Intrinsic/sequence mode -- NOT Structural mode, no structure exists yet at Step 2A)")
    print("               Record the number from 'The protein variant intrinsic solubility")
    print("               score is <x>'. HIGHER IS MORE SOLUBLE; CamSol's published")
    print("               convention is >1 highly soluble, <-1 poorly soluble. This is NOT")
    print("               a patch length -- true aggregation PATCH LENGTH is measured later")
    print("               at Step 2C from CamSol STRUCTURAL, which does report one.")
    print("-" * 100)

    for name, seq in variants.items():
        print(f"\n[{name}]")
        print(seq)
    print("-" * 100)

    # Keyed on the sanitized variant id (matching every other Step 2X
    # results file), not the raw FASTA header -- the raw header carries
    # "| length=NNN" metadata that changes on every Phase 1G rerun even
    # when the underlying construct id doesn't, which used to invalidate
    # every existing key each time. If the file already exists (e.g. from
    # a previous, different construct), missing candidate keys are added
    # to it rather than leaving the user stuck with no key to fill in --
    # existing filled-in values are never touched.
    if os.path.isfile(manual_results_path):
        with open(manual_results_path) as f:
            data = json.load(f)
    else:
        data = {}

    added = []
    for name in variants:
        safe_name = common.sanitize_variant_name(name)
        if safe_name not in data:
            data[safe_name] = {k: None for k in MANUAL_SEQUENCE_TOOL_KEYS}
            added.append(safe_name)

    os.makedirs(os.path.dirname(manual_results_path), exist_ok=True)
    with open(manual_results_path, 'w') as f:
        json.dump(data, f, indent=2)

    if added:
        print(f"[INFO] Added {len(added)} new candidate entry(ies) to: {manual_results_path}")
    else:
        print(f"[INFO] {manual_results_path} already has entries for all current candidates.")
    print(f"[INFO] Fill in every null field for all {len(variants)} candidate(s), then rerun this script.")
    print("=" * 100 + "\n")


def load_manual_results(manual_results_path, candidate_names):
    if not os.path.isfile(manual_results_path):
        return None
    with open(manual_results_path) as f:
        data = json.load(f)
    for name in candidate_names:
        entry = data.get(common.sanitize_variant_name(name))
        if entry is None or any(manual_value(entry, k) is None for k in MANUAL_SEQUENCE_TOOL_KEYS):
            return None
    return data

# =============================================================================
# BIOCHEMICAL & METHODOLOGICAL EVALUATION ENGINE
# =============================================================================

def get_sliding_windows(seq, window_size):
    return [seq[i:i + window_size] for i in range(len(seq) - window_size + 1)]


def check_hydrophobicity(seq):
    """Fails if any 8, 12, or 15 aa window is > 80% hydrophobic."""
    # Per methodology: {A, V, I, L, M, F, W, Y} -- not {A, C, F, I, L, M, V, W}.
    # Cys was wrongly included (it's not on the paper's hydrophobic list)
    # and Tyr was wrongly omitted.
    hydrophobic_residues = set(['A', 'V', 'I', 'L', 'M', 'F', 'W', 'Y'])
    for w_size in [8, 12, 15]:
        for window in get_sliding_windows(seq, w_size):
            h_count = sum(1 for aa in window if aa in hydrophobic_residues)
            if (h_count / w_size) > 0.80:
                return False
    return True


# Longest junction peptide the methodology enumerates (MHC-II upper bound).
# Shared by generate_junction_peptides() and the adjuvant-exemption logic so
# the two cannot drift apart.
JUNCTION_MAX_LEN = 20

# Mature human beta-defensin-3, the adjuvant the methodology fixes at the
# N-terminus (Phase 1G attaches it via an EAAAK linker). Matched by prefix so
# the exemption is anchored to the actual adjuvant present, not a hardcoded
# residue index -- if the adjuvant is ever changed, the exemption follows it.
BETA_DEFENSIN_3 = "GIINTLQKYYCRVRGGRCAVLSCLPKEEQIGKCSTRGRKCCRRKK"


def _adjuvant_region_end(seq):
    """
    Returns the index one past the adjuvant's last residue, or 0 if the
    construct does not begin with the known adjuvant. Used to exempt
    junctions that reach back into the adjuvant's native disulfide cluster
    from the >2-cysteine rejection rule -- see the note at the call site.
    """
    return len(BETA_DEFENSIN_3) if seq.startswith(BETA_DEFENSIN_3) else 0


def generate_junction_peptides(seq, linkers=None):
    """
    For every occurrence of every Phase I linker, enumerates ALL peptides
    spanning that junction at each length required by the methodology:
    MHC I lengths 8-11 aa, MHC II lengths 12-20 aa.

    Returns: {junction_id: {"linker": str, "position": int,
                             "MHC_I": [peptides...], "MHC_II": [peptides...]}}
    """
    if linkers is None:
        linkers = ["AAY", "GPGPG", "KK", "EAAAK"]

    mhc_i_lengths = range(8, 12)
    mhc_ii_lengths = range(12, 21)
    max_len = JUNCTION_MAX_LEN

    junctions = {}

    for linker in linkers:
        start_idx = 0
        while True:
            idx = seq.find(linker, start_idx)
            if idx == -1:
                break

            linker_end = idx + len(linker)
            context_start = max(0, idx - (max_len - 1))
            context_end = min(len(seq), linker_end + (max_len - 1))
            context = seq[context_start:context_end]

            local_link_start = idx - context_start
            local_link_end = linker_end - context_start

            entry = {"linker": linker, "position": idx, "MHC_I": set(), "MHC_II": set()}

            for length_range, key in [(mhc_i_lengths, "MHC_I"), (mhc_ii_lengths, "MHC_II")]:
                for L in length_range:
                    if L > len(context):
                        continue
                    for w_start in range(0, len(context) - L + 1):
                        w_end = w_start + L
                        # A true junction peptide must SPAN the linker
                        # boundary -- fully contain the linker with at
                        # least one residue of flanking epitope on BOTH
                        # sides -- not merely overlap it at one edge.
                        # The previous overlap-only test ("touches the
                        # linker at all") counted huge numbers of windows
                        # that were >90% epitope sequence and only
                        # clipped the linker by a residue or two, which
                        # inherit whatever Cys count their epitope has
                        # regardless of the linker itself.
                        if w_start < local_link_start and w_end > local_link_end:
                            entry[key].add(context[w_start:w_end])

            entry["MHC_I"] = sorted(entry["MHC_I"])
            entry["MHC_II"] = sorted(entry["MHC_II"])
            junctions[f"{linker}@{idx}"] = entry

            start_idx = idx + 1

    return junctions


def _run_with_heartbeat(cmd, cwd, label, heartbeat_interval=10):
    """
    Runs a subprocess while printing a periodic "still running" heartbeat
    with elapsed time, so a long external tool call doesn't look frozen.
    Most of these CLI tools (blastp included) don't expose a real
    percent-progress for a plain local search, so this reports elapsed
    time honestly instead of fabricating a percentage.
    """
    start = time.time()
    process = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    while process.poll() is None:
        elapsed = time.time() - start
        print(f"[INFO] {label} still running... ({common.format_time(elapsed)} elapsed)")
        time.sleep(heartbeat_interval)
    stdout, stderr = process.communicate()
    elapsed = time.time() - start
    print(f"[INFO] {label} finished after {common.format_time(elapsed)}.")

    class _Result:
        pass
    result = _Result()
    result.returncode = process.returncode
    result.stdout = stdout
    result.stderr = stderr
    return result


def evaluate_candidate(name, seq, junctions, toxinpred_result, hemopi_result, blastp_hit,
                        mhc_i_ranks, mhc_ii_ranks, junction_toxinpred_results, manual_result):
    """Executes the Phase II Step A methodology checks on a given sequence."""
    ana = ProteinAnalysis(seq)
    instability_idx = ana.instability_index()
    gravy_val = ana.gravy()
    pI = ana.isoelectric_point()
    mw = ana.molecular_weight()
    aliphatic_idx = common.compute_aliphatic_index(seq)
    half_life = common.estimate_half_life_hours(seq)

    result = {
        "Variant": name,
        "Len": len(seq),
        "MW_Da": round(mw, 2),
        "Aliphatic_Index": round(aliphatic_idx, 2),
        "Est_HalfLife_hr": half_life if half_life is not None else "N/A",
        "STAB_IDX": round(instability_idx, 2),
        "GRAVY": round(gravy_val, 4),
        "pI": round(pI, 2),
        "Viable": "NO",
        "Rejection_Reasons": [],
        "Review_Flags": [],
    }

    # =========================================================================
    # 1. TOXICITY SCREENING
    # =========================================================================
    if not check_hydrophobicity(seq):
        result["Rejection_Reasons"].append("Hydrophobicity >80% in window")

    # Junction cysteine check -- applied to EPITOPE-EPITOPE junctions only.
    #
    # METHODOLOGY NOTE (disclose this): the methodology mandates a
    # beta-defensin-3 adjuvant AND a ">2 cysteines in a junction peptide"
    # rejection rule. These two requirements are in direct conflict:
    # mature human beta-defensin-3 natively carries SIX cysteines forming
    # three disulfide bonds -- that disulfide ladder is the defining
    # structural feature of a defensin, not an incidental property.
    #
    # Measured on this construct: 116 junction peptides exceed the limit,
    # and every single one comes from just three junction sites (KK@43,
    # EAAAK@45, KK@49) -- all of them windows that reach back into the
    # adjuvant's C-terminal cysteines at positions 32/39/40. All 34
    # epitope-epitope junctions are clean. Applied literally, the rule
    # rejects any construct using the adjuvant the methodology requires.
    #
    # The rule's purpose is to catch toxin-like, cysteine-rich peptides
    # created incidentally at NOVEL epitope boundaries. The adjuvant is a
    # known, characterised, deliberately included molecule whose cysteines
    # are functionally essential -- it is not an unintended junction risk.
    # Junctions overlapping it are therefore exempted, and the exemption is
    # reported explicitly rather than applied silently.
    adjuvant_end = _adjuvant_region_end(seq)
    exempt_ids, checked_ids = [], []
    cys_hit = False
    for jid, entry in junctions.items():
        window_start = max(0, entry["position"] - (JUNCTION_MAX_LEN - 1))
        if window_start < adjuvant_end:
            exempt_ids.append(jid)
            continue
        checked_ids.append(jid)
        if any(pep.count('C') > 2 for pep in entry["MHC_I"] + entry["MHC_II"]):
            cys_hit = True

    if cys_hit:
        result["Rejection_Reasons"].append("Junction Cys > 2 (epitope-epitope junction)")
    if exempt_ids:
        result["Review_Flags"].append(
            f"Cys check: {len(exempt_ids)} adjuvant-overlapping junction(s) exempted "
            f"({', '.join(sorted(exempt_ids))}); {len(checked_ids)} epitope-epitope junction(s) checked"
        )

    # Fail-closed: a missing real-tool result must never silently pass
    # toxicity/hemolytic screening -- it means the tool didn't confirm
    # this candidate is safe, which is not the same as confirming it is.
    # Evidence-based rule (Phase 1E deviation #1): reject on BLAST/MERCI
    # evidence, not on the miscalibrated composition score.
    tp_evidence = _toxinpred_has_evidence(toxinpred_result)
    if tp_evidence is None:
        result["Rejection_Reasons"].append("ToxinPred2: no result returned -- cannot confirm non-toxic (fail-closed)")
    elif tp_evidence:
        result["Rejection_Reasons"].append(
            f"ToxinPred2: Toxin (evidence: BLAST={toxinpred_result['blast']}, MERCI={toxinpred_result['merci']})"
        )
    elif toxinpred_result["ml"] >= 0.6:
        # Advisory only -- composition score alone never rejects.
        result["Review_Flags"].append(
            f"ToxinPred2 COMPOSITION-FLAG (ML={toxinpred_result['ml']}, no BLAST/MERCI evidence) -- advisory, not a rejection"
        )

    # Paper's rule (Sec. I.E, applied to whole constructs via Sec. IV.A.b):
    # hemolysis prob >= 0.50 AND Hydrophobic_Fraction > 0.60 -- HemoPI2's
    # own collapsed verdict alone is stricter than what the paper specifies
    # and was previously used directly, which is not the paper's rule.
    hydrophobic_residues = set(['A', 'V', 'I', 'L', 'M', 'F', 'W', 'Y'])
    hydro_fraction = sum(1 for aa in seq if aa in hydrophobic_residues) / len(seq)
    if hemopi_result is None:
        result["Rejection_Reasons"].append("HemoPI2: no result returned -- cannot confirm non-hemolytic (fail-closed)")
    elif hemopi_result["score"] >= 0.50 and hydro_fraction > 0.60:
        result["Rejection_Reasons"].append(
            f"HemoPI2: Hemolytic (score={hemopi_result['score']}, Hydrophobic_Fraction={hydro_fraction:.3f} > 0.60)"
        )
    elif hemopi_result["score"] >= 0.50:
        result["Review_Flags"].append(
            f"HemoPI2 score {hemopi_result['score']} >= 0.50 but Hydrophobic_Fraction "
            f"{hydro_fraction:.3f} <= 0.60 -- not rejected per paper's conjunctive rule"
        )

    if blastp_hit is not None and blastp_hit.get("significant_hit"):
        result["Rejection_Reasons"].append(
            f"BLASTP Tox-Prot hit (E={blastp_hit['evalue']:.2e}, subject={blastp_hit['subject']})"
        )

    # Toxic junction = a junction peptide that is BOTH a strong MHC binder
    # AND carries real toxin evidence (BLAST homology or MERCI motif).
    # Same evidence rule as the full construct above -- composition score
    # alone does not make a junction toxic.
    toxic_junction_hits = []
    composition_only_hits = 0
    unresolved_hits = 0
    for entry in junctions.values():
        for pep, rank_map, rank_cut in (
            [(p, mhc_i_ranks, 2.0) for p in entry["MHC_I"]]
            + [(p, mhc_ii_ranks, 10.0) for p in entry["MHC_II"]]
        ):
            rank = rank_map.get(pep)
            if rank is None or rank > rank_cut:
                continue
            tp = junction_toxinpred_results.get(pep)
            evidence = _toxinpred_has_evidence(tp)
            if evidence is None:
                unresolved_hits += 1
            elif evidence:
                toxic_junction_hits.append((pep, rank, tp))
            elif tp["ml"] >= 0.6:
                composition_only_hits += 1

    if toxic_junction_hits:
        pep, rank, tp = toxic_junction_hits[0]
        result["Rejection_Reasons"].append(
            f"Toxic Junction (MHC): {len(toxic_junction_hits)} peptide(s) with toxin evidence "
            f"(e.g. {pep} rank={rank}, BLAST={tp['blast']}, MERCI={tp['merci']})"
        )
    if unresolved_hits:
        result["Rejection_Reasons"].append(
            f"Toxic Junction (MHC): {unresolved_hits} strong-binder peptide(s) with unresolved "
            f"ToxinPred2 result -- cannot confirm non-toxic (fail-closed)"
        )
    if composition_only_hits:
        result["Review_Flags"].append(
            f"{composition_only_hits} strong-binder junction peptide(s) COMPOSITION-FLAG only "
            f"(ML>=0.6, zero BLAST/MERCI evidence) -- advisory, not a rejection"
        )

    # =========================================================================
    # 2. ALLERGENICITY SCREENING
    # =========================================================================
    qn_ratio = (seq.count('Q') + seq.count('N')) / len(seq)
    # Paper's charged-residue set (Sec. I.E and II.A.a, stated twice) is
    # {D, E, H, K}.
    surface_charge = sum(seq.count(aa) for aa in ['D', 'E', 'H', 'K'])
    # Paper: "Surface_Charge greater than 4 residues or 2% of the sequence
    # length... flags the sequence for review for peptide lengths greater
    # than 20 aa (for short junction peptides, a raw Surface_Charge > 4 is
    # automatically deprioritized)". The >4 floor is pinned to SHORT
    # peptides (used as-is on junction peptides below); for a long
    # construct the 2% clause is meant to scale the bound UP, not be
    # capped down to 4.
    charge_threshold = max(4, round(len(seq) * 0.02))

    if surface_charge > charge_threshold:
        result["Review_Flags"].append(
            f"High Surface Charge ({surface_charge} residues > threshold {charge_threshold})"
        )

    if manual_result is None:
        result["Review_Flags"].append("AllerTOP/AllergenFP: no manual result recorded")
        allergen_count = 0
    else:
        aller_top = bool(manual_result.get("allertop_is_allergen"))
        aller_fp = bool(manual_result.get("allergenfp_is_allergen"))
        allergen_count = sum([aller_top, aller_fp])

    if allergen_count == 2:
        result["Rejection_Reasons"].append("Allergenic (Both Predictors)")
    elif qn_ratio > 0.30 and allergen_count >= 1:
        result["Rejection_Reasons"].append("Allergenic (QN > 0.30 + 1 Predictor)")

    for entry in junctions.values():
        all_peps = entry["MHC_I"] + entry["MHC_II"]
        high_charge_peps = [sum(p.count(aa) for aa in ['D', 'E', 'H', 'K']) for p in all_peps
                             if sum(p.count(aa) for aa in ['D', 'E', 'H', 'K']) > 4]
        high_qn_peps = [(p.count('Q') + p.count('N')) / len(p) for p in all_peps
                         if (p.count('Q') + p.count('N')) / len(p) > 0.30]
        if high_charge_peps:
            result["Review_Flags"].append(
                f"Junction '{entry['linker']}'@{entry['position']}: "
                f"{len(high_charge_peps)}/{len(all_peps)} peptides Surface_Charge>4 (max {max(high_charge_peps)})"
            )
        if high_qn_peps:
            result["Review_Flags"].append(
                f"Junction '{entry['linker']}'@{entry['position']}: "
                f"{len(high_qn_peps)}/{len(all_peps)} peptides QN_Ratio>0.30 (max {max(high_qn_peps):.2f})"
            )
    if junctions:
        result["Review_Flags"].append(
            "Junction-level AllerTOP/AllergenFP consensus not automated (manual "
            "submission burden scales with candidate x junction count) -- not "
            "used as a rejection trigger in this run"
        )

    # =========================================================================
    # 3. SOLUBILITY & STABILITY SCREENING
    #    DeepSol (replacing SOLpro, same 0.50 threshold -- see Step 2C for
    #    the same substitution rationale) and CamSol (sequence/intrinsic
    #    mode) are real, via the same manual submission as AllerTOP/
    #    AllergenFP above.
    #    (Length, MW, Aliphatic Index, Half-Life are CONTEXTUAL ONLY.)
    # =========================================================================
    if instability_idx > 40:
        result["Rejection_Reasons"].append("Instability > 40")
    if gravy_val > 0.4:
        result["Rejection_Reasons"].append("GRAVY > 0.4")
    if 6.5 <= pI <= 7.5:
        result["Rejection_Reasons"].append("pI near 7.0 (buffer range)")

    # Fail-closed: missing manual solubility/aggregation data must never
    # silently pass as "assume fine" -- it means the screen wasn't done,
    # not that it passed. (In practice this manual gate runs before any
    # of this evaluation, so `manual_result` is already guaranteed
    # complete by the time we get here -- these checks are the defensive
    # fallback if that guarantee is ever bypassed.)
    if manual_result is None or manual_result.get("deepsol_probability") is None:
        result["Rejection_Reasons"].append("DeepSol: no manual result recorded -- cannot confirm solubility (fail-closed)")
        is_insoluble = None
    else:
        solpro_score = manual_result["deepsol_probability"]
        is_insoluble = solpro_score < 0.5

    # CamSol Intrinsic solubility score. Accepts the legacy key name so result
    # files written before deviation #25 still load unchanged.
    camsol_score = manual_value(manual_result, "camsol_intrinsic_solubility_score")
    if camsol_score is None:
        result["Rejection_Reasons"].append("CamSol (sequence): no manual result recorded -- cannot confirm aggregation propensity (fail-closed)")
        has_long_patch = None
    else:
        # CamSol's own published convention (Sormanni, Aprile & Vendruscolo,
        # J Mol Biol 2015): score < -1 = poorly soluble / aggregation-prone,
        # > +1 = highly soluble, in between = intermediate. This replaces a
        # `> 8` test against a value that is not a length and cannot reach 8.
        has_long_patch = camsol_score < CAMSOL_POORLY_SOLUBLE
        camsol_patch = camsol_score

    # DeepSol insolubility is held OUT of Rejection_Reasons and tracked
    # separately so the final verdict below can distinguish "insoluble and
    # nothing else" (-> REVIEW, deviation #15) from "insoluble AND some
    # other failure" (-> hard NO). The paper's compound case
    # (patch >8aa AND insoluble, Sec. II.A.b) stays a hard rejection.
    deepsol_sole_failure = False
    if has_long_patch and is_insoluble:
        result["Rejection_Reasons"].append(
            f"CamSol poorly soluble ({camsol_patch}) + DeepSol insoluble")
    elif is_insoluble:
        deepsol_sole_failure = True
    if has_long_patch and not is_insoluble:
        result["Review_Flags"].append(
            f"CamSol intrinsic solubility score {camsol_patch} < {CAMSOL_POORLY_SOLUBLE} (poorly soluble)")

    # =========================================================================
    # FINAL VERDICT
    # =========================================================================
    if result["Rejection_Reasons"]:
        # A genuine failure on any other criterion -- DeepSol status is
        # appended for the record but the verdict is NO either way.
        if deepsol_sole_failure:
            result["Rejection_Reasons"].append(f"DeepSol Insoluble ({solpro_score})")
        result["Rejection_Reasons"] = " | ".join(result["Rejection_Reasons"])
    elif deepsol_sole_failure and DEEPSOL_ONLY_FAILURE_IS_REVIEW:
        # DEVIATION #15 (documented): every other Sec. II.A criterion passes
        # and DeepSol is the sole failure. Carried forward as REVIEW rather
        # than rejected outright. Set DEEPSOL_ONLY_FAILURE_IS_REVIEW = False
        # to restore the paper's literal "reject" behaviour.
        result["Viable"] = "REVIEW"
        result["Rejection_Reasons"] = (
            f"REVIEW -- DeepSol Insoluble ({solpro_score}) is the SOLE failing criterion; "
            f"all other Sec. II.A checks pass. Carried forward flagged, not rejected "
            f"(deviation #15)."
        )
    elif deepsol_sole_failure:
        result["Rejection_Reasons"] = f"DeepSol Insoluble ({solpro_score})"
    else:
        result["Viable"] = "YES"
        result["Rejection_Reasons"] = "None (Optimal)"

    result["Review_Flags"] = " | ".join(result["Review_Flags"]) if result["Review_Flags"] else "None"

    return result

# =============================================================================
# MAIN EXECUTION THREAD
# =============================================================================

def _report_junction_cys_precheck(all_junctions_by_candidate):
    """
    P0-3 pre-check: reports the (corrected) Junction-Cys rule's verdict
    for every candidate BEFORE any external tool is invoked, using only
    the already-generated junction peptide sets (pure/local, no network
    or subprocess cost). This is what previously took a full run --
    ToxinPred2, HemoPI2, BLASTP, and six IEDB POSTs -- to discover, since
    the gate sat at the very end of evaluate_candidate().
    """
    common.print_banner("PRE-CHECK: JUNCTION CYS > 2 RULE (before any external tool runs)")
    any_hit = False
    for name, junctions in all_junctions_by_candidate.items():
        hits = [
            pep for entry in junctions.values()
            for pep in entry["MHC_I"] + entry["MHC_II"]
            if pep.count('C') > 2
        ]
        n_checked = sum(len(entry["MHC_I"]) + len(entry["MHC_II"]) for entry in junctions.values())
        if hits:
            any_hit = True
            print(f"[PRE-CHECK] {name}: WOULD REJECT -- {len(hits)}/{n_checked} true junction-spanning "
                  f"peptides have >2 Cys (e.g. {hits[0]})")
        else:
            print(f"[PRE-CHECK] {name}: PASSES -- 0/{n_checked} true junction-spanning peptides have >2 Cys")
    if any_hit:
        print("[PRE-CHECK] At least one candidate will be rejected on this rule alone -- if that's the")
        print("[PRE-CHECK] ONLY rejection reason once the full run finishes, treat it as a genuine result,")
        print("[PRE-CHECK] not evidence of a leftover bug (see P0-3 in the remediation plan).")
    print("=" * 115 + "\n")


def run_step2a_comprehensive_screening():
    start_time = time.time()
    project_root = _PROJECT_ROOT

    input_fasta_path = common.phase1g_fasta_path(project_root)
    output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA")
    raw_out_dir = os.path.join(output_base, "Raw")
    filt_out_dir = os.path.join(output_base, "Filtered")
    work_dir = os.path.join(output_base, "_tool_runs")

    os.makedirs(raw_out_dir, exist_ok=True)
    os.makedirs(filt_out_dir, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)

    common.print_banner("PHASE 2 STEP A: COMPREHENSIVE SECONDARY SCREENING DASHBOARD")
    print(f"[INFO] Resolved Project Root : {project_root}")
    print(f"[INFO] Initialization Time   : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    if input_fasta_path is None:
        search_dir = os.path.join(project_root, common.PHASE1G_DIR_RELATIVE)
        print(f"\n[FATAL ERROR] No file matching '{common.PHASE1G_FASTA_GLOB}' found in: {search_dir}")
        print("[FATAL ERROR] Run Step 1G first so a Phase1G_FinalConstruct_*.fasta exists.")
        sys.exit(1)
    print(f"[INFO] Input Source          : {input_fasta_path}")
    print("-" * 125)

    variants = common.load_multi_fasta(input_fasta_path)
    if variants is None:
        print(f"\n[FATAL ERROR] Variant FASTA file not found: {input_fasta_path}")
        sys.exit(1)
    if not variants:
        print("[WARNING] No FASTA variants found in Phase 1G. Please check Step 1G execution.")
        return

    safe_id_to_name = {common.sanitize_variant_name(name): name for name in variants}

    # Junction generation is pure/local (sequence in, peptide sets out) --
    # do it up front so both the P0-3 pre-check and the P0-2 manual gate
    # below can happen before any external tool is invoked.
    all_junctions_by_candidate = {name: generate_junction_peptides(seq) for name, seq in variants.items()}
    _report_junction_cys_precheck(all_junctions_by_candidate)

    # P0-2: the manual-results gate sits here, BEFORE any external tool
    # call -- ToxinPred2, HemoPI2, BLASTP, and six IEDB POSTs used to run
    # to completion first and get discarded on every first-ever run,
    # because the gate previously sat at the very end of the pipeline.
    manual_results_path = os.path.join(output_base, "manual_sequence_tool_results.json")
    manual_results = load_manual_results(manual_results_path, list(variants.keys()))
    if manual_results is None:
        run_manual_prepare(variants, manual_results_path)
        return

    # Every cache filename below is keyed off a hash of the ACTUAL variant
    # content, not a fixed name like "all_candidates.fasta" -- a fixed
    # name means a stale cache from a PREVIOUS, different candidate set
    # (e.g. an earlier Phase 1G run's 5 candidates) gets silently reused
    # for a completely different current candidate set, since the caching
    # check only looked at "does this file already exist," never "does it
    # actually describe what we're about to ask." Content changes ->
    # filename changes -> cache naturally misses instead of lying.
    candidates_signature = hashlib.sha256(
        "\n".join(f"{common.sanitize_variant_name(n)}={s}" for n, s in sorted(variants.items())).encode()
    ).hexdigest()[:10]

    full_fasta_path = os.path.join(work_dir, f"all_candidates_{candidates_signature}.fasta")
    with open(full_fasta_path, 'w') as f:
        for name, seq in variants.items():
            f.write(f">{common.sanitize_variant_name(name)}\n{seq}\n")

    print("[INFO] Running ToxinPred2 (full constructs)...")
    try:
        toxinpred_raw = run_toxinpred2_batch(full_fasta_path, work_dir)
    except Exception as e:
        print(f"[ERROR] ToxinPred2 failed: {e}")
        return
    toxinpred_results = {safe_id_to_name.get(k, k): v for k, v in toxinpred_raw.items()}

    print("[INFO] Running HemoPI2 (full constructs)...")
    try:
        hemopi_raw = run_hemopi2_batch(full_fasta_path, work_dir)
    except Exception as e:
        print(f"[ERROR] HemoPI2 failed: {e}")
        return
    hemopi_results = {safe_id_to_name.get(k, k): v for k, v in hemopi_raw.items()}

    print("[INFO] Running BLASTP vs Tox-Prot (full constructs)...")
    try:
        blastp_raw = run_blastp_toxprot_batch(full_fasta_path, work_dir)
    except Exception as e:
        print(f"[ERROR] BLASTP failed: {e}")
        return
    blastp_results = {safe_id_to_name.get(k, k): v for k, v in blastp_raw.items()}

    all_mhc_i_peptides, all_mhc_ii_peptides = set(), set()
    for junctions in all_junctions_by_candidate.values():
        for entry in junctions.values():
            all_mhc_i_peptides.update(entry["MHC_I"])
            all_mhc_ii_peptides.update(entry["MHC_II"])

    print(f"[INFO] Running IEDB MHC binding ({len(all_mhc_i_peptides)} MHC-I + {len(all_mhc_ii_peptides)} MHC-II unique peptides)...")
    try:
        mhc_i_cache = os.path.join(work_dir, f"iedb_mhc_i_ranks_{candidates_signature}.json")
        mhc_ii_cache = os.path.join(work_dir, f"iedb_mhc_ii_ranks_{candidates_signature}.json")
        mhc_i_ranks = run_iedb_mhc_binding_batch(sorted(all_mhc_i_peptides), "MHC_I", cache_path=mhc_i_cache) if all_mhc_i_peptides else {}
        mhc_ii_ranks = run_iedb_mhc_binding_batch(sorted(all_mhc_ii_peptides), "MHC_II", cache_path=mhc_ii_cache) if all_mhc_ii_peptides else {}
    except Exception as e:
        print(f"[ERROR] IEDB MHC binding API failed: {e}")
        return

    all_junction_peptides = sorted(all_mhc_i_peptides | all_mhc_ii_peptides)
    junction_toxinpred_results = {}
    if all_junction_peptides:
        print(f"[INFO] Running ToxinPred2 ({len(all_junction_peptides)} unique junction peptides)...")
        junction_fasta_path = os.path.join(work_dir, f"all_junction_peptides_{candidates_signature}.fasta")
        pid_to_pep = {}
        with open(junction_fasta_path, 'w') as f:
            for i, pep in enumerate(all_junction_peptides):
                pid = f"pep{i}"
                pid_to_pep[pid] = pep
                f.write(f">{pid}\n{pep}\n")
        try:
            junction_toxinpred_raw = run_toxinpred2_batch(junction_fasta_path, work_dir)
        except Exception as e:
            print(f"[ERROR] ToxinPred2 (junction peptides) failed: {e}")
            return
        junction_toxinpred_results = {pid_to_pep.get(k, k): v for k, v in junction_toxinpred_raw.items()}

    print("-" * 125)

    raw_results = []
    for name, seq in variants.items():
        result = evaluate_candidate(
            name, seq, all_junctions_by_candidate[name],
            toxinpred_result=toxinpred_results.get(name),
            hemopi_result=hemopi_results.get(name),
            blastp_hit=blastp_results.get(name),
            mhc_i_ranks=mhc_i_ranks,
            mhc_ii_ranks=mhc_ii_ranks,
            junction_toxinpred_results=junction_toxinpred_results,
            manual_result=manual_results.get(common.sanitize_variant_name(name)),
        )
        raw_results.append(result)

    print(f"{'VARIANT':<22} | {'LEN':<4} | {'STAB':<6} | {'GRAVY':<6} | {'pI':<5} | {'VIABLE':<6} | {'REASONS':<48} | {'REVIEW_FLAGS'}")
    print("-" * 125)
    for r in raw_results:
        short_reasons = (r['Rejection_Reasons'][:45] + '...') if len(r['Rejection_Reasons']) > 45 else r['Rejection_Reasons']
        short_flags = (r['Review_Flags'][:35] + '...') if len(r['Review_Flags']) > 35 else r['Review_Flags']
        print(f"{r['Variant'][:22]:<22} | {r['Len']:<4} | {r['STAB_IDX']:<6} | {r['GRAVY']:<6} | {r['pI']:<5} | {r['Viable']:<6} | {short_reasons:<48} | {short_flags}")

    # REVIEW candidates (deviation #15: DeepSol is the sole failure) are
    # carried forward alongside YES, but sorted after them and labelled so
    # the distinction is never lost downstream.
    filtered_list = [r for r in raw_results if r['Viable'] in ("YES", "REVIEW")]
    filtered_list.sort(key=lambda x: (x['Viable'] != "YES", x['STAB_IDX']))

    common.print_banner("STEP 2A: FILTERED OPTIMAL CANDIDATES (RANKED BY STABILITY INDEX)")
    if filtered_list:
        print(f"{'RANK':<5} | {'VARIANT':<22} | {'STABILITY INDEX':<18} | {'GRAVY':<10} | {'STATUS'}")
        print("-" * 125)
        for i, f in enumerate(filtered_list):
            if f['Viable'] == "REVIEW":
                status_label = "[REVIEW -- DeepSol sole failure, deviation #15]"
            else:
                status_label = "[OPTIMAL]" if i == 0 else "[ACCEPTED]"
            print(f"{i+1:<5} | {f['Variant']:<22} | {f['STAB_IDX']:<18} | {f['GRAVY']:<10} | {status_label}")
        if any(f['Viable'] == "REVIEW" for f in filtered_list):
            print("-" * 125)
            print("[NOTE] REVIEW candidates fail ONLY the DeepSol <0.50 solubility gate and are")
            print("[NOTE] carried forward under documented deviation #15 -- they are NOT clean passes.")
    else:
        print(f"{'--- NO VARIANTS PASSED THE COMPREHENSIVE SCREENING CRITERIA ---':^125}")

    ts = datetime.now().strftime("%Y%m%d_%H%M")
    raw_path = os.path.join(raw_out_dir, f"Step2A_Raw_Physicochemical_{ts}.csv")
    filt_path = os.path.join(filt_out_dir, f"Step2A_Filtered_Ranked_{ts}.csv")

    with open(raw_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=raw_results[0].keys())
        writer.writeheader()
        writer.writerows(raw_results)

    if filtered_list:
        with open(filt_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=filtered_list[0].keys())
            writer.writeheader()
            writer.writerows(filtered_list)

    exec_time = common.format_time(time.time() - start_time)
    common.print_banner("PHASE 2 STEP A COMPLETE")
    print(f"[SUCCESS] Total Analyzed      : {len(raw_results)} variants")
    print(f"[SUCCESS] Viable Candidates   : {len(filtered_list)}")
    print(f"[SUCCESS] Execution Time      : {exec_time}")
    print(f"[INFO] Raw Log: {os.path.relpath(raw_path, project_root)}")
    print(f"[INFO] Filtered Log: {os.path.relpath(filt_path, project_root)}")
    print("=" * 125 + "\n")

if __name__ == "__main__":
    run_step2a_comprehensive_screening()
