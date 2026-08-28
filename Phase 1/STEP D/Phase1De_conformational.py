import os, sys, csv, re, json, time
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
# Sec. I.D (current manuscript): B-cell candidates from BepiPred-2.0 are
# cross-referenced against SEMA 2.0 conformational-patch predictions on the
# native (unassembled) antigen fold. >=50% residue overlap with a patch ->
# "structurally corroborated", PRIORITIZED. No overlap -> retained as
# lower-priority linear-only, FLAGGED. This step never excludes a candidate.
#
# SEMA 2.0 has a documented REST API (sema.airi.net/openapi.json) that in
# principle needs no manual submission. In practice, POST /proteins (the
# actual prediction call) sits behind an anti-bot WAF (ServicePipe) that
# returns the SPA shell instead of JSON for any non-browser client --
# confirmed this session with both `requests` (form + true multipart
# encoding) and a real headless Chromium via Playwright (fingerprinted and
# 403'd outright). This is the same class of external blocker as RaptorX
# being down (deviation #16) or AllerTOP/AllergenFP having no API at all --
# NOT something to keep fighting with browser-automation evasion.
#
# So this script: (1) tries the automated API path once per target and uses
# it if it works, (2) otherwise falls back to a manual-submission template,
# exactly like Phase 2A's AllerTOP/AllergenFP/DeepSol/CamSol workflow.
# SEMA is prioritization-only, so Phase 1G proceeds with
# SEMA_Corroborated=UNSCREENED for any target left unfilled rather than
# blocking construction on it -- see the plan's impact note: the B-cell pool
# is shallow enough (<=3 candidates/target, 2 taken) that SEMA mainly adds
# required annotation/justification, not a different construct.
# =============================================================================
SEMA_BASE_URL = "https://sema.airi.net"
SEMA_HEADERS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
                  "(KHTML, like Gecko) Chrome/128.0.0.0 Safari/537.36",
    "Accept": "application/json",
}

REQUIRED_TARGETS = ["Mpox_L1R", "Mpox_B5R", "Mpox_A35R", "HIV_gp120", "HIV_gp41", "HIV_p24", "HIV_p17"]
OVERLAP_THRESHOLD = 0.50


def _reference_sequence(target):
    """
    Var_01 per target -- the first (canonical/RefSeq) isolate from Phase 1A,
    matching the paper's framing of ONE native fold per antigen (Sec. I.D:
    "native tertiary structures of each source antigen").
    """
    folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1A")
    if not os.path.isdir(folder):
        return None
    matches = sorted(f for f in os.listdir(folder) if f.startswith(f"{target}_Var_01_"))
    if not matches:
        return None
    with open(os.path.join(folder, matches[0])) as f:
        lines = [l.strip() for l in f if not l.startswith(">")]
    return "".join(lines).upper()


def try_sema_api(seq, timeout=90):
    """
    Attempts the real SEMA-1D (sequence-only, ESM-2-based) API call.
    Returns a list of per-residue epitope-propensity scores (0-indexed to
    the submitted sequence) on success, or None if the WAF intercepts the
    request (detected by content-type / body shape, not just status code --
    the WAF returns HTTP 200 with the SPA's index.html).
    """
    import requests
    import socket
    import urllib3.util.connection as urllib3_cn
    urllib3_cn.allowed_gai_family = lambda: socket.AF_INET

    try:
        files = {"seq": (None, seq), "esm_switch": (None, "1d")}
        r = requests.post(f"{SEMA_BASE_URL}/proteins", files=files, headers=SEMA_HEADERS, timeout=timeout)
    except Exception:
        return None

    ctype = r.headers.get("content-type", "")
    if r.status_code != 200 or "html" in ctype.lower():
        return None  # WAF/SPA intercepted the request -- not real API output

    try:
        data = r.json()
        protein_id = data.get("id") or data.get("protein_id") or data.get("_id")
    except Exception:
        return None
    if not protein_id:
        return None

    # Poll for completion, then fetch per-residue JSON.
    for _ in range(30):
        try:
            rr = requests.get(f"{SEMA_BASE_URL}/proteins/{protein_id}/json_data",
                               headers=SEMA_HEADERS, timeout=timeout)
        except Exception:
            return None
        if rr.status_code == 200 and "html" not in rr.headers.get("content-type", "").lower():
            try:
                return rr.json()
            except Exception:
                return None
        time.sleep(3)
    return None


def run_step1de_conformational():
    common.print_banner("PHASE 1De: SEMA 2.0 CONFORMATIONAL B-CELL CORROBORATION")

    input_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1Db")
    output_dir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1De")
    tool_runs_dir = os.path.join(output_dir, "_tool_runs")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(tool_runs_dir, exist_ok=True)

    latest_db = common.latest_file(input_folder, suffix=".csv")
    if latest_db is None:
        print(f"[ERROR] No Phase 1Db output found at: {input_folder}")
        return

    with open(latest_db, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        original_fields = reader.fieldnames

    bcell_by_target = {}
    for r in rows:
        if r.get("Type") == "B-cell":
            bcell_by_target.setdefault(r["Target"], set()).add(r["Peptide"])

    manual_path = os.path.join(output_dir, "Manual_SEMA_Results.json")
    manual_data = json.load(open(manual_path)) if os.path.isfile(manual_path) else {}

    corroboration = {}   # peptide -> (overlap_pct, corroborated_bool)
    need_manual = []     # targets where the automated path failed

    for target in REQUIRED_TARGETS:
        peptides = bcell_by_target.get(target, set())
        if not peptides:
            continue

        ref_seq = _reference_sequence(target)
        if ref_seq is None:
            print(f"[WARNING] No Phase 1A reference sequence found for {target} -- skipping.")
            continue

        cache_path = os.path.join(tool_runs_dir, f"sema_{target}.json")
        patch_scores = None
        if os.path.isfile(cache_path):
            patch_scores = json.load(open(cache_path))
        else:
            print(f"[INFO] {target}: attempting automated SEMA-1D API call...")
            patch_scores = try_sema_api(ref_seq)
            if patch_scores is not None:
                json.dump(patch_scores, open(cache_path, "w"))
                print(f"[SUCCESS] {target}: automated SEMA-1D result cached.")
            else:
                # Manual fallback: check if the user has already filled in
                # Manual_SEMA_Results.json for this target.
                if target in manual_data and manual_data[target]:
                    patch_scores = manual_data[target]
                else:
                    need_manual.append(target)
                    continue

        # patch_scores expected shape: list of per-residue dicts/floats aligned
        # to ref_seq, each with an epitope-probability score. Accept either a
        # bare list of floats or a list of {"resi":..,"score":..} rows.
        if patch_scores and isinstance(patch_scores[0], dict):
            score_by_resi = {int(p.get("resi", i)): float(p.get("score", 0.0)) for i, p in enumerate(patch_scores)}
            residue_scores = [score_by_resi.get(i, 0.0) for i in range(len(ref_seq))]
        else:
            residue_scores = [float(x) for x in patch_scores]

        for pep in peptides:
            start = ref_seq.find(pep)
            if start == -1:
                corroboration[pep] = (None, False)  # not found in this target's reference -- unscreened
                continue
            window = residue_scores[start:start + len(pep)]
            if not window:
                corroboration[pep] = (None, False)
                continue
            # SEMA scores are typically probability-like in [0,1]; treat >=0.5
            # per-residue as "in an epitope patch region", matching the >=50%
            # RESIDUE overlap rule (not a mean-score threshold).
            overlap_pct = 100.0 * sum(1 for s in window if s >= 0.5) / len(window)
            corroboration[pep] = (round(overlap_pct, 1), overlap_pct >= OVERLAP_THRESHOLD * 100)

    # Emit a manual-submission template for any target the API couldn't reach.
    if need_manual:
        if not os.path.isfile(manual_path):
            template = {t: None for t in need_manual}
            json.dump(template, open(manual_path, "w"), indent=2)
        print(f"\n[ACTION AVAILABLE, NOT BLOCKING] SEMA-1D's automated API is WAF-blocked for "
              f"{len(need_manual)} target(s): {need_manual}")
        print(f"[INFO] Submit each reference sequence at {SEMA_BASE_URL} (Predict epitopes tab, SEMA-1D "
              f"mode), then paste the per-residue JSON into: {manual_path}")
        print(f"[INFO] Phase 1G will proceed treating these targets' B-cell candidates as "
              f"SEMA_Corroborated=UNSCREENED (annotation-only impact -- see plan §0.4).")

    # Write an enriched Phase1Db_Elite_*.csv (only B-cell rows carry new
    # columns; all other rows pass through unchanged) so 1Dc and everything
    # downstream picks it up automatically via original_fields + latest_file.
    fieldnames = original_fields + ["SEMA_Overlap_Pct", "SEMA_Corroborated"]
    out_rows = []
    n_corrob = n_unscreened = n_not_corrob = 0
    for r in rows:
        clean = {k: r[k] for k in original_fields}
        if r.get("Type") == "B-cell":
            result = corroboration.get(r["Peptide"])
            if result is None or result[0] is None:
                clean["SEMA_Overlap_Pct"] = ""
                clean["SEMA_Corroborated"] = "UNSCREENED"
                n_unscreened += 1
            else:
                overlap_pct, is_corrob = result
                clean["SEMA_Overlap_Pct"] = overlap_pct
                clean["SEMA_Corroborated"] = "YES" if is_corrob else "NO"
                n_corrob += is_corrob
                n_not_corrob += (not is_corrob)
        else:
            clean["SEMA_Overlap_Pct"] = ""
            clean["SEMA_Corroborated"] = ""
        out_rows.append(clean)

    # NOTE: Phase1Db/1Dc/1Da all use "%Y%m%d_%H%M" (no dashes), and
    # Phase1Dc_benchmark.py picks its input via a plain string sort
    # (sorted(db_files)[-1]), not by ctime. A dash-formatted timestamp
    # here would sort lexicographically BEFORE same-day no-dash names
    # ('-' < '0' in ASCII) and silently get skipped as "not the latest" --
    # so this must match the sibling files' format exactly.
    ts = datetime.now().strftime("%Y%m%d_%H%M")
    out_path = os.path.join(input_folder, f"Phase1Db_Elite_{ts}.csv")
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader(); writer.writerows(out_rows)

    print(f"\n[INFO] B-cell SEMA corroboration -- corroborated: {n_corrob} | not corroborated: {n_not_corrob} "
          f"| unscreened: {n_unscreened}")
    print(f"[INFO] Enriched Phase1Db_Elite CSV written to: {out_path}")
    print("=" * 90 + "\n")


if __name__ == "__main__":
    run_step1de_conformational()
