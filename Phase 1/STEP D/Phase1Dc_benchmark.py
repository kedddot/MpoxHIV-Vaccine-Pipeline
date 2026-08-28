import os
import csv
import re
import time
import sys
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
# METHODOLOGY NOTE: IEDB Conservancy Tool vs. Exact Substring Match
# =============================================================================
# The paper cites the "IEDB Epitope Conservancy Analysis Tool".
# However, IEDB does not provide a documented public REST API for this specific
# tool, and local standalone binaries have highly variable/unstandardized CLI
# interfaces depending on the installation environment.
#
# Because the IEDB tool set to 100% sequence identity (the standard for exact
# epitope matching) operates mathematically as a strict substring search, this
# script utilizes exact substring matching as a same-quality, deterministic
# substitute. This guarantees pipeline stability without hallucinating unverified
# CLI parameters or relying on fragile web-scraping endpoints.
#
# This substitution is written out as METHODOLOGY_NOTE.md alongside every run's
# output so the manuscript's methods section can be updated to describe it
# explicitly, rather than the code silently diverging from what the paper claims.
# =============================================================================

METHODOLOGY_NOTE_TEXT = """# Phase 1Dc — Conservancy Methodology Note

**Applies to:** Section II.C.I.D of the proposal ("IEDB Conservancy Tool").

## What the manuscript currently says
The methods section states that epitope conservancy was assessed using the
IEDB Epitope Conservancy Analysis Tool.

## What this script actually does
IEDB does not expose a documented public REST API for the Conservancy
Analysis tool (unlike its MHC-I/MHC-II binding APIs), and standalone local
distributions of the tool have inconsistent, environment-dependent CLI
interfaces. Rather than depend on an unverified binary or silently fall back
to something undocumented, this script computes conservancy directly:

- For each target (e.g. `Mpox_L1R`, `HIV_gp120`), every surviving variant
  FASTA file from Phase 1C is loaded into a per-target pool of full-length
  sequences.
- For each candidate peptide, conservancy is calculated as the percentage of
  that target's variant pool containing an **exact, 100%-identity substring
  match** of the peptide: `Conservancy = (hits / total_variants) * 100`.

## Why this is equivalent, not just convenient
The IEDB Conservancy tool's identity threshold controls how much mismatch is
tolerated when searching for the epitope within each sequence. At a 100%
identity threshold specifically, no mismatches or gaps are permitted — which
is mathematically the same operation as an exact substring search. This
script always evaluates at 100% identity, so its output is equivalent to
running the IEDB tool at that same threshold on the same sequence set.

## What this does NOT reproduce
This substitution only holds at the 100% identity threshold. It does not
reproduce IEDB's behavior at lower identity thresholds (e.g. 80% or 90%
conservancy allowing partial mismatches), which the current script does not
attempt to compute. If partial-identity conservancy is ever required, this
approach would need to be replaced with an actual alignment-based method.

## Recommended manuscript update
Section II.C.I.D should be revised to state that conservancy was computed via
exact substring matching against each target's full variant pool at 100%
sequence identity, functionally equivalent to the IEDB Conservancy Tool at
that threshold, rather than citing the tool as if its CLI/API were invoked
directly.
"""

def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

def run_step1dc_conservancy_benchmark():
    start_time = time.time()

    input_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1Db")
    fasta_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1C", "Filtered_Antigenicity")
    output_base = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1Dc")

    raw_dir = os.path.join(output_base, "Raw_Conservancy")
    filt_dir = os.path.join(output_base, "Filtered_Benchmarks")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(filt_dir, exist_ok=True)

    print("\n" + "="*80)
    print(f"{'PHASE 1Dc: CONSERVANCY & VARIANT BENCHMARKING':^80}")
    print("="*80)

    if not os.path.isdir(input_folder):
        print(f"[ERROR] Phase 1Db output directory not found at: {input_folder}")
        return
    db_files = [f for f in os.listdir(input_folder) if f.endswith(".csv")]
    if not db_files:
        print("[ERROR] No input CSV found in Phase 1Db.")
        return
    latest_db = os.path.join(input_folder, sorted(db_files)[-1])

    # Pre-load ALL variants per target into a comprehensive library
    print(f"[PROCESS] Loading full viral variant library from {os.path.basename(fasta_folder)}...")
    library = {}
    for f in os.listdir(fasta_folder):
        if f.endswith(".fasta"):
            target = f.split('_Var')[0]
            with open(os.path.join(fasta_folder, f), "r") as file:
                seq = "".join([l.strip() for l in file if not l.startswith(">")])
                clean_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq.upper())

                if target not in library:
                    library[target] = []
                library[target].append(clean_seq)

    total_variants = sum(len(v) for v in library.values())
    print(f"[INFO] Loaded {total_variants} total variant sequences across {len(library)} targets.")

    print(f"[PROCESS] Analyzing conservancy using exact-match substitute (100% identity threshold)...")
    final_results = []

    with open(latest_db, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        total_rows = len(rows)

        for i, row in enumerate(rows):
            pep = row['Peptide']
            target = row['Target']

            # Retrieve the full pool of variant sequences for this specific target
            relevant_vars = library.get(target, [])
            total_relevant = len(relevant_vars)

            if total_relevant == 0:
                row['Conservancy'] = 0.0
                row['Hit_Ratio'] = "0/0"
                final_results.append(row)
                continue

            # Exact substring match evaluation against the entire variant pool
            hits = sum(1 for seq in relevant_vars if pep in seq)

            row['Conservancy'] = round((hits / total_relevant) * 100, 2)
            row['Hit_Ratio'] = f"{hits}/{total_relevant}"
            final_results.append(row)

            # Progress tracker
            if i % 10 == 0 or i == total_rows - 1:
                elapsed = format_time(time.time() - start_time)
                sys.stdout.write(f"\r[ PROCESS ] {i+1:03d}/{total_rows:03d} | Calculating Conservancy | Elapsed: {elapsed}")
                sys.stdout.flush()

    print("\n[INFO] Conservancy analysis complete.")

    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    raw_path = os.path.join(raw_dir, f"Phase1Dc_Raw_Full_{ts}.csv")

    if final_results:
        with open(raw_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=final_results[0].keys())
            writer.writeheader()
            writer.writerows(final_results)

    print("-" * 80)
    print(f"{'BENCHMARK THRESHOLD RESULTS':^80}")

    # Strict thresholds matching the methodology wording (> 20%, > 50%, >= 100% ceiling)
    thresholds = [
        (20, "Minimal (> 20%)", lambda c: c > 20),
        (50, "Highly Acceptable (> 50%)", lambda c: c > 50),
        (100, "Ideal (100%)", lambda c: c >= 100),
    ]

    for t, label, test in thresholds:
        filtered = [r for r in final_results if test(float(r['Conservancy']))]
        if filtered:
            # Each tier gets its own subdirectory (Min_20pct/, Min_50pct/,
            # Min_100pct/) rather than a shared flat folder. Previously all
            # three landed in Filtered_Benchmarks/ together, and Phase 1Ea's
            # ctime-based "pick the latest file" always grabbed the 100%
            # tier (written last) regardless of which tier was intended --
            # silently discarding the >20% and >50% pools with no warning.
            tier_dir = os.path.join(filt_dir, f"Min_{t}pct")
            os.makedirs(tier_dir, exist_ok=True)
            f_path = os.path.join(tier_dir, f"Phase1Dc_Min_{t}pct_{ts}.csv")
            with open(f_path, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=final_results[0].keys())
                writer.writeheader()
                writer.writerows(filtered)
            print(f"[SUCCESS] {label}: {len(filtered)} epitopes found.")

            if "Type" in filtered[0]:
                type_counts = {}
                for r in filtered:
                    type_counts[r["Type"]] = type_counts.get(r["Type"], 0) + 1
                breakdown = " | ".join(f"{k}: {v}" for k, v in sorted(type_counts.items()))
                print(f"           -> {breakdown}")

    # Write the methodology note alongside this run's outputs so the
    # paper/code gap is closed at the source rather than left implicit.
    note_path = os.path.join(output_base, "METHODOLOGY_NOTE.md")
    with open(note_path, "w") as f:
        f.write(METHODOLOGY_NOTE_TEXT)
    print(f"[INFO] Methodology note written to {note_path}")

    print("="*80 + "\n")

if __name__ == "__main__":
    run_step1dc_conservancy_benchmark()
