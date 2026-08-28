import os, sys, time, re, csv, requests
from datetime import datetime
from collections import defaultdict

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
# ALLELE PANEL -- MHC-I/II alleles queried against the IEDB recommended-method
# tools for percentile-rank binding prediction.
# =============================================================================
MHCI_ALLELES = ["HLA-A*24:02", "HLA-B*15:02", "HLA-B*40:01", "HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02", "HLA-B*08:01", "HLA-C*07:01", "HLA-C*04:01"]

# MHC-II panel: DRB1 + DQ + DP, matching the 21-allele panel Phase1F scores
# coverage against. Previously this was DRB1-only (5 alleles), which meant
# 11 of Phase1F's 21 panel alleles -- including HLA-DPB1*05:01, the single
# highest-frequency allele in the whole MHC-II panel -- could never appear
# in a candidate's Binding_Alleles no matter how well it actually bound,
# capping cumulative MHC-II coverage at ~84% regardless of epitope choice.
#
# IEDB requires paired alpha/beta notation for DQ/DP (verified live --
# "HLA-DPB1*05:01" alone returns "Invalid allele name"). MHCnuggets (used
# in Phase1F) instead expects single-chain, no-asterisk notation for the
# SAME allele ("HLA-DPB105:01"). MHCII_ALLELE_TO_SINGLE_CHAIN maps this
# script's IEDB-format keys to Phase1F's ALLELE_FREQ keys so both sides
# agree on which real-world allele each entry represents.
MHCII_ALLELES = [
    "HLA-DRB1*15:01", "HLA-DRB1*12:02", "HLA-DRB1*09:01", "HLA-DRB1*04:05",
    "HLA-DRB1*08:03", "HLA-DRB1*03:01", "HLA-DRB1*04:03", "HLA-DRB1*07:01", "HLA-DRB1*14:01",
    "HLA-DQA1*01:02/DQB1*06:02", "HLA-DQA1*05:01/DQB1*03:01", "HLA-DQA1*01:01/DQB1*05:01",
    "HLA-DQA1*05:01/DQB1*02:01", "HLA-DQA1*03:01/DQB1*04:02", "HLA-DQA1*04:01/DQB1*03:02",
    "HLA-DPA1*01:03/DPB1*05:01", "HLA-DPA1*01:03/DPB1*02:01", "HLA-DPA1*01:03/DPB1*04:01",
    "HLA-DPA1*02:01/DPB1*09:01", "HLA-DPA1*01:03/DPB1*14:01", "HLA-DPA1*02:01/DPB1*11:01",
]
# Maps this script's IEDB paired-notation allele names to the beta-chain key
# used by Phase1F's ALLELE_FREQ table (which uses standard asterisk notation).
# Phase1F recomputes its own binding set via MHCflurry/MHCnuggets rather than
# consuming this script's Binding_Alleles column, so nothing depends on this
# map today -- it exists so the Binding_Alleles values written to this step's
# CSV remain traceable to a frequency-table entry for reporting and audit.
# Target format MUST stay asterisk-notation to match ALLELE_FREQ; a no-asterisk
# form would silently resolve to frequency 0 via ALLELE_FREQ.get's default.
MHCII_ALLELE_TO_SINGLE_CHAIN = {
    "HLA-DQA1*01:02/DQB1*06:02": "HLA-DQB1*06:02", "HLA-DQA1*05:01/DQB1*03:01": "HLA-DQB1*03:01",
    "HLA-DQA1*01:01/DQB1*05:01": "HLA-DQB1*05:01", "HLA-DQA1*05:01/DQB1*02:01": "HLA-DQB1*02:01",
    "HLA-DQA1*03:01/DQB1*04:02": "HLA-DQB1*04:02", "HLA-DQA1*04:01/DQB1*03:02": "HLA-DQB1*03:02",
    "HLA-DPA1*01:03/DPB1*05:01": "HLA-DPB1*05:01", "HLA-DPA1*01:03/DPB1*02:01": "HLA-DPB1*02:01",
    "HLA-DPA1*01:03/DPB1*04:01": "HLA-DPB1*04:01", "HLA-DPA1*02:01/DPB1*09:01": "HLA-DPB1*09:01",
    "HLA-DPA1*01:03/DPB1*14:01": "HLA-DPB1*14:01", "HLA-DPA1*02:01/DPB1*11:01": "HLA-DPB1*11:01",
}

MHCI_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
MHCII_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
BCELL_URL = "https://tools-cluster-interface.iedb.org/tools_api/bcell/"


def print_banner(text): print(f"\n{'='*80}\n{text:^80}\n{'='*80}")

def get_gravy(pep):
    hydro = {'A': 1.8, 'L': 3.8, 'I': 4.5, 'V': 4.2, 'F': 2.8, 'M': 1.9, 'C': 2.5, 'G': -0.4,
             'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'N': -3.5, 'Q': -3.5,
             'D': -3.5, 'E': -3.5, 'K': -3.9, 'R': -4.5, 'H': -3.2}
    return sum(hydro.get(aa, 0) for aa in pep) / len(pep)

def classify_bcell_tier(mean_score, pct_above):
    if mean_score >= 0.60 and pct_above >= 75.0: return "High"
    elif mean_score >= 0.50 and pct_above >= 50.0: return "Medium"
    elif mean_score >= 0.45 and pct_above >= 37.5: return "Deprioritized"
    return "Excluded"

def find_column(header, keywords):
    for kw in keywords:
        for i, col in enumerate(header):
            if kw.lower() in col.lower(): return i
    return None

def find_column_exact(header, names):
    """
    Exact (case-insensitive) column match, tried in order. Needed for the
    peptide column specifically: IEDB's MHC-II response header is
    "... core_peptide  peptide  score  rank", and find_column's substring
    match returns core_peptide (the 9-mer binding core) since it's the
    first column containing "peptide" -- silently substituting a 9-mer for
    the intended 15-mer epitope on every MHC-II row. Falls back to
    find_column's substring behaviour only if no exact match exists, so
    this is safe to use even if a future IEDB response format differs.
    """
    lower_header = [c.lower() for c in header]
    for name in names:
        if name.lower() in lower_header:
            return lower_header.index(name.lower())
    return find_column(header, names)

def print_status(current, total, target, start_time, kept, action):
    elapsed = time.time() - start_time
    msg = f"[PROCESS] {current:02d}/{total:02d} | Target: {target:<12} | Action: {action:<20} | Elapsed: {elapsed:5.1f}s | Kept: {kept}"
    sys.stdout.write("\r" + " " * 100)
    sys.stdout.write(f"\r{msg}")
    sys.stdout.flush()


def run_step1db_optimized(mhcii_only=False):
    """
    mhcii_only=True re-runs ONLY the MHC-II predictions and merges them with
    the MHC-I and B-cell rows from the most recent prior Phase1Db_Elite CSV.

    This exists because the peptide-column bug fixed by find_column_exact()
    affected MHC-II exclusively -- MHC-I and B-cell rows from the prior run
    are provably correct (verified lengths 9/10 and 16 respectively) and are
    not touched by either that bug or the MHC-II allele-panel expansion.
    Re-running everything costs ~2h; MHC-II alone is a fraction of that, with
    identical output for the reused classes.
    """
    start_time = time.time()

    fasta_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1C", "Filtered_Antigenicity")
    identification_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1Da")
    output_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1Db")
    os.makedirs(output_folder, exist_ok=True)

    print_banner("PHASE 1Db: IEDB AFFINITY & SOLUBILITY FILTER")

    # In mhcii_only mode, load the prior run's MHC-I / B-cell rows up front so
    # we fail fast if they're missing rather than after a long IEDB run.
    carried_rows = []
    if mhcii_only:
        prior_files = sorted([f for f in os.listdir(output_folder) if f.endswith(".csv")])
        if not prior_files:
            print(f"\n[ERROR] --mhcii-only needs a prior Phase1Db_Elite CSV to merge with; none found in {output_folder}")
            return
        prior_path = os.path.join(output_folder, prior_files[-1])
        with open(prior_path, 'r', newline='') as f:
            for r in csv.DictReader(f):
                if r.get("Type") != "MHC-II":
                    carried_rows.append(r)
        n_i = sum(1 for r in carried_rows if r.get("Type") == "MHC-I")
        n_b = sum(1 for r in carried_rows if r.get("Type") == "B-cell")
        print(f"[MODE] MHC-II ONLY -- reusing {n_i} MHC-I and {n_b} B-cell rows from {os.path.basename(prior_path)}")
        print(f"[INFO] Querying {len(MHCII_ALLELES)} MHC-II alleles per peptide (MHC-I and B-cell skipped).")
    else:
        print(f"[INFO] Querying {len(MHCI_ALLELES)} MHC-I and {len(MHCII_ALLELES)} MHC-II alleles per peptide — runtime will scale accordingly.")

    if not os.path.exists(identification_folder):
        print(f"\n[ERROR] Phase 1Da directory not found at: {identification_folder}")
        return
    da_files = [f for f in os.listdir(identification_folder) if f.endswith(".csv")]
    if not da_files: return
    latest_da = os.path.join(identification_folder, sorted(da_files)[-1])

    candidates_by_variant = defaultdict(list)
    with open(latest_da, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader: candidates_by_variant[row["Variant"]].append(row)

    fasta_files = sorted([f for f in os.listdir(fasta_folder) if f.endswith(".fasta")])

    results_map = {} # Maps peptide string to data dict for deduplication
    skipped = {"MHC-I": 0, "MHC-II": 0, "B-cell": 0}
    gravy_deprioritized = {"MHC-I": 0, "MHC-II": 0, "B-cell": 0}
    api_errors = 0
    session = requests.Session()

    # Response cache keyed by (endpoint, allele/method, length, sequence).
    #
    # PERFORMANCE (this dominates total runtime -- read before changing):
    # IEDB accepts a COMMA-SEPARATED allele list in a single request and
    # returns one row per (allele, peptide), with the allele named in the
    # response's own first column. Verified live: all 21 MHC-II alleles in
    # one call returns 1513 rows in ~3.6s, versus ~12s for each of 21
    # separate calls under sustained load (IEDB throttles per request).
    #
    # Combined with iterating over UNIQUE SEQUENCES rather than files (167
    # input files collapse to 57 distinct sequences -- NCBI surveillance
    # dumps carry byte-identical isolates across many accessions), this cuts
    # the MHC-II pass from 3507 requests (~11.7 h) to 57 (~4 min). Same API,
    # same method, same thresholds -- purely how requests are packed. Fewer
    # requests also means far less exposure to throttle-induced errors.
    response_cache = {}
    cache_hits = [0]  # mutable counter closed over by cached_post

    def cached_post(url, payload, cache_key):
        if cache_key in response_cache:
            cache_hits[0] += 1
            return response_cache[cache_key]
        try:
            response = session.post(url, data=payload, timeout=180)
            text = response.text if response.status_code == 200 else None
        except Exception:
            text = None
        response_cache[cache_key] = text
        return text

    # Group input files by their cleaned sequence so each distinct sequence is
    # submitted to IEDB exactly once; results are then written back to every
    # file sharing it. Preserves per-file Target/Variant attribution.
    seq_groups = defaultdict(list)   # clean_seq -> [(f_name, target, candidates)]
    for f_name in fasta_files:
        my_candidates = candidates_by_variant.get(f_name, [])
        if not my_candidates:
            continue
        with open(os.path.join(fasta_folder, f_name), "r") as f:
            raw = "".join([l.strip() for l in f if not l.startswith(">")])
        clean_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', raw.upper())
        if len(clean_seq) < 9:
            continue
        seq_groups[clean_seq].append((f_name, f_name.split('_Var')[0], my_candidates))
    print(f"[INFO] {len(fasta_files)} input files collapse to {len(seq_groups)} distinct sequences "
          f"-- IEDB is queried once per distinct sequence, with all alleles batched per request.")

    for i, (clean_seq, members) in enumerate(seq_groups.items()):
        # Every file sharing this sequence gets the same IEDB result; the
        # first member supplies the display target for progress output.
        _, disp_target, _ = members[0]
        bcell_by_member = [(fn, tgt, [c for c in cands if c["Type"] == "B-cell"])
                           for fn, tgt, cands in members]

        # ---------------- MHC-I ----------------
        # All alleles batched into ONE request per length. The response names
        # its own allele per row, so Binding_Alleles is read from the row --
        # NOT from a loop variable (which is what the per-allele version did).
        if not mhcii_only:
            for length in (9, 10):
                print_status(i+1, len(seq_groups), disp_target, start_time, len(results_map), f"MHC-I all-alleles ({length}m)")
                allele_str = ",".join(MHCI_ALLELES)
                cache_key = ("mhci_batch", allele_str, length, clean_seq)
                was_cached = cache_key in response_cache
                payload = {'method': 'recommended', 'sequence_text': clean_seq,
                           'allele': allele_str, 'length': ",".join([str(length)] * len(MHCI_ALLELES))}
                text = cached_post(MHCI_URL, payload, cache_key)
                if text is not None:
                    lines = text.strip().split('\n')
                    if len(lines) > 1:
                        header = lines[0].split('\t')
                        r_idx = find_column(header, ['percentile_rank', 'rank'])
                        p_idx = find_column_exact(header, ['peptide'])
                        a_idx = find_column_exact(header, ['allele'])
                        if r_idx is not None and p_idx is not None and a_idx is not None:
                            for line in lines[1:]:
                                cols = line.split('\t')
                                if len(cols) <= max(r_idx, p_idx, a_idx): continue
                                pep = cols[p_idx]; row_allele = cols[a_idx]
                                try: rank = float(cols[r_idx])
                                except: continue

                                if rank <= 1.0:
                                    gravy = get_gravy(pep)
                                    is_deprioritized = gravy >= 0.2
                                    for f_name, target, _ in members:
                                        pep_key = f"{target}_{pep}_MHCI"
                                        if pep_key not in results_map:
                                            if is_deprioritized:
                                                gravy_deprioritized["MHC-I"] += 1
                                            results_map[pep_key] = {"Target": target, "Variant": f_name, "Type": "MHC-I", "Length": length, "Peptide": pep, "GRAVY": round(gravy, 3), "GRAVY_Deprioritized": is_deprioritized, "Percentile_Rank": rank, "mean_BepiPred": "", "pct_above": "", "Bcell_Tier": "", "Binding_Alleles": set()}
                                        results_map[pep_key]["Binding_Alleles"].add(row_allele)
                                        results_map[pep_key]["Percentile_Rank"] = min(results_map[pep_key]["Percentile_Rank"], rank)
                                else: skipped["MHC-I"] += 1
                else:
                    api_errors += 1
                if not was_cached:
                    time.sleep(1.0)

        # ---------------- MHC-II ----------------
        # All 21 alleles in ONE request (see the performance note above).
        print_status(i+1, len(seq_groups), disp_target, start_time, len(results_map), f"MHC-II all {len(MHCII_ALLELES)} alleles")
        allele_str = ",".join(MHCII_ALLELES)
        cache_key = ("mhcii_batch", allele_str, 15, clean_seq)
        was_cached = cache_key in response_cache
        payload = {'method': 'recommended', 'sequence_text': clean_seq,
                   'allele': allele_str, 'length': ",".join(['15'] * len(MHCII_ALLELES))}
        text = cached_post(MHCII_URL, payload, cache_key)
        if text is not None:
            lines = text.strip().split('\n')
            if len(lines) > 1:
                header = lines[0].split('\t')
                r_idx = find_column(header, ['percentile_rank', 'adjusted_rank', 'rank'])
                p_idx = find_column_exact(header, ['peptide'])
                a_idx = find_column_exact(header, ['allele'])
                if r_idx is not None and p_idx is not None and a_idx is not None:
                    for line in lines[1:]:
                        cols = line.split('\t')
                        if len(cols) <= max(r_idx, p_idx, a_idx): continue
                        pep = cols[p_idx]; row_allele = cols[a_idx]
                        try: rank = float(cols[r_idx])
                        except: continue

                        if rank <= 10.0:
                            gravy = get_gravy(pep)
                            is_deprioritized = gravy >= 0.2
                            for f_name, target, _ in members:
                                pep_key = f"{target}_{pep}_MHCII"
                                if pep_key not in results_map:
                                    if is_deprioritized:
                                        gravy_deprioritized["MHC-II"] += 1
                                    results_map[pep_key] = {"Target": target, "Variant": f_name, "Type": "MHC-II", "Length": 15, "Peptide": pep, "GRAVY": round(gravy, 3), "GRAVY_Deprioritized": is_deprioritized, "Percentile_Rank": rank, "mean_BepiPred": "", "pct_above": "", "Bcell_Tier": "", "Binding_Alleles": set()}
                                results_map[pep_key]["Binding_Alleles"].add(row_allele)
                                results_map[pep_key]["Percentile_Rank"] = min(results_map[pep_key]["Percentile_Rank"], rank)
                        else: skipped["MHC-II"] += 1
        else:
            api_errors += 1
        if not was_cached:
            time.sleep(1.0)

        # ---------------- B-cell ----------------
        # BepiPred scores the whole sequence once; per-residue scores are then
        # sliced per candidate. Since every member of this group shares the
        # sequence, one call serves them all -- but each member keeps its own
        # candidate list and Target attribution.
        any_bcell = any(bc for _, _, bc in bcell_by_member)
        if any_bcell and not mhcii_only:
            print_status(i+1, len(seq_groups), disp_target, start_time, len(results_map), "B-cell")
            cache_key = ("bcell", "Bepipred-2.0", None, clean_seq)
            was_cached = cache_key in response_cache
            payload = {'method': 'Bepipred-2.0', 'sequence_text': clean_seq}
            text = cached_post(BCELL_URL, payload, cache_key)
            if text is not None:
                lines = text.strip().split('\n')
                if len(lines) > 1:
                    header = lines[0].split('\t')
                    score_idx = find_column(header, ['score'])
                    if score_idx is not None:
                        residue_scores = [float(line.split('\t')[score_idx]) for line in lines[1:] if line.strip()]
                        for f_name, target, bcell_candidates in bcell_by_member:
                            for cand in bcell_candidates:
                                start = int(cand["Start_Position"])
                                window = residue_scores[start:start + 16]
                                if len(window) < 16: continue
                                mean_score = sum(window) / len(window)
                                pct_above = 100.0 * sum(1 for s in window if s >= 0.50) / len(window)
                                tier = classify_bcell_tier(mean_score, pct_above)

                                if tier in ("High", "Medium", "Deprioritized"):
                                    gravy = get_gravy(cand["Peptide"])
                                    is_deprioritized = gravy >= 0.2
                                    pep_key = f"{target}_{cand['Peptide']}_BCell"
                                    if pep_key not in results_map:
                                        if is_deprioritized:
                                            gravy_deprioritized["B-cell"] += 1
                                        results_map[pep_key] = {"Target": target, "Variant": f_name, "Type": "B-cell", "Length": 16, "Peptide": cand["Peptide"], "GRAVY": round(gravy, 3), "GRAVY_Deprioritized": is_deprioritized, "Percentile_Rank": "", "mean_BepiPred": round(mean_score, 3), "pct_above": round(pct_above, 2), "Bcell_Tier": tier, "Binding_Alleles": set()}
                                else: skipped["B-cell"] += 1
            else:
                api_errors += 1
            if not was_cached:
                time.sleep(1.0)

    print()

    final_results = []
    for val in results_map.values():
        val["Binding_Alleles"] = ";".join(sorted(list(val["Binding_Alleles"])))
        final_results.append(val)

    # Guard against a silent regression of the core_peptide bug: IEDB returns
    # both a 9-mer core_peptide and the real 15-mer peptide for MHC-II, and
    # picking the wrong column is invisible in the output unless length is
    # checked explicitly. Abort rather than write a corrupt candidate set.
    bad_len = [r for r in final_results if r["Type"] == "MHC-II" and len(r["Peptide"]) != 15]
    if bad_len:
        print(f"\n[FATAL] {len(bad_len)} MHC-II rows are not 15 aa (found lengths: "
              f"{sorted({len(r['Peptide']) for r in bad_len})}).")
        print("[FATAL] This is the core_peptide column bug -- check find_column_exact(). Aborting without writing.")
        sys.exit(1)

    n_new_mhcii = len(final_results)
    if mhcii_only:
        final_results = carried_rows + final_results

    out_file = os.path.join(output_folder, f"Phase1Db_Elite_{datetime.now().strftime('%Y%m%d_%H%M')}.csv")
    fieldnames = ["Target", "Variant", "Type", "Length", "Peptide", "GRAVY", "GRAVY_Deprioritized", "Percentile_Rank", "mean_BepiPred", "pct_above", "Bcell_Tier", "Binding_Alleles"]

    with open(out_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(final_results)

    print_banner("FILTRATION COMPLETE")
    print("[OK] All MHC-II peptides verified at 15 aa.")
    if mhcii_only:
        print(f"[MERGE] {n_new_mhcii} new MHC-II rows + {len(carried_rows)} carried MHC-I/B-cell rows")
    print(f"[SUCCESS] Candidates retained : {len(final_results)}")
    print(f"[INFO] Below-threshold (rank/tier) skipped : MHC-I={skipped['MHC-I']} MHC-II={skipped['MHC-II']} B-cell={skipped['B-cell']}")
    print(f"[INFO] Retained but GRAVY-deprioritized (>=0.2) : MHC-I={gravy_deprioritized['MHC-I']} MHC-II={gravy_deprioritized['MHC-II']} B-cell={gravy_deprioritized['B-cell']}")
    print(f"[INFO] API errors / non-200 responses : {api_errors}")
    print(f"[INFO] API calls avoided via duplicate-sequence cache : {cache_hits[0]}")

if __name__ == "__main__":
    run_step1db_optimized(mhcii_only="--mhcii-only" in sys.argv)
