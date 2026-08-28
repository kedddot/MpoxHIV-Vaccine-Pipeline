import os
import sys
import time
import re
import csv
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
# UTILITY FUNCTIONS
# =============================================================================

def format_time(seconds):
    """Formats time in seconds to a readable MM:SS format."""
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

def extract_sliding_window(sequence, window_size):
    """Generates overlapping k-mers from a sequence."""
    return [sequence[i:i + window_size] for i in range(len(sequence) - window_size + 1)]

# =============================================================================
# CORE PROCESSING FUNCTION
# =============================================================================

def run_phase1da_identification():
    start_time = time.time()

    input_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1C", "Filtered_Antigenicity")
    output_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1Da")
    os.makedirs(output_folder, exist_ok=True)

    if not os.path.exists(input_folder):
        print(f"\n[ERROR] Phase 1A directory not found at: {input_folder}")
        return

    # 2. EXPERIMENTAL PREVIEW & HEADER
    print("\n" + "="*80)
    print(f"{'PHASE 1Da: MULTI-TYPE EPITOPE LIBRARY GENERATION':^80}")
    print("="*80)
    print(f"[INFO] Initialization Time : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[INFO] Input Directory     : {input_folder}")
    print(f"[INFO] Target Protocol     : MHC-I (9/10-mer), MHC-II (15-mer), B-cell (16-mer)")
    print("-" * 80)

    fasta_files = sorted([f for f in os.listdir(input_folder) if f.endswith(".fasta")])
    total_files = len(fasta_files)
    
    all_peptides = []
    
    # Benchmarking counters
    stats = {"MHC-I": 0, "MHC-II": 0, "B-cell": 0}

    # 3. SLIDING WINDOW EXECUTION
    for i, filename in enumerate(fasta_files):
        target_group = filename.split('_Var')[0]
        
        # Dynamic terminal status
        elapsed = format_time(time.time() - start_time)
        status_str = f"[ PROCESS ] Record {i+1:03d}/{total_files:03d} | Target: {target_group:<12} | Elapsed: {elapsed}"
        sys.stdout.write(f"\r{status_str:<80}")
        sys.stdout.flush()

        with open(os.path.join(input_folder, filename), "r") as f:
            seq = "".join([line.strip() for line in f if not line.startswith(">")])
            clean_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq.upper())

        # Generate Peptides
        # NOTE: each loop below now uses enumerate() over the sliding-window output
        # to record the 0-based Start_Position of every peptide within clean_seq.
        # This is required by Phase1Db, which submits the *whole* sequence to the
        # BepiPred-2.0 B-cell endpoint and needs each 16-mer's exact position to
        # slice out the matching per-residue scores.

        # MHC-I
        for length in [9, 10]:
            for start, pep in enumerate(extract_sliding_window(clean_seq, length)):
                all_peptides.append({"Target": target_group, "Variant": filename, "Type": "MHC-I", "Length": length, "Start_Position": start, "Peptide": pep})
                stats["MHC-I"] += 1
        
        # MHC-II
        for start, pep in enumerate(extract_sliding_window(clean_seq, 15)):
            all_peptides.append({"Target": target_group, "Variant": filename, "Type": "MHC-II", "Length": 15, "Start_Position": start, "Peptide": pep})
            stats["MHC-II"] += 1
            
        # B-cell
        for start, pep in enumerate(extract_sliding_window(clean_seq, 16)):
            all_peptides.append({"Target": target_group, "Variant": filename, "Type": "B-cell", "Length": 16, "Start_Position": start, "Peptide": pep})
            stats["B-cell"] += 1

    print() # New line after status bar
    print("-" * 80)

    # 4. DATA EXPORT & BENCHMARK REPORTING
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M")
    csv_path = os.path.join(output_folder, f"Phase1Da_Master_Library_{timestamp}.csv")
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Target", "Variant", "Type", "Length", "Start_Position", "Peptide"])
        writer.writeheader()
        writer.writerows(all_peptides)

    # Final Academic Summary
    total_time = format_time(time.time() - start_time)
    print(f"{'EPITOPE GENERATION COMPLETE':^80}")
    print("="*80)
    print(f"[SUCCESS] Total Peptides Generated : {len(all_peptides)}")
    print(f"[SUCCESS] Breakdown:")
    print(f"          - MHC-I   : {stats['MHC-I']}")
    print(f"          - MHC-II  : {stats['MHC-II']}")
    print(f"          - B-cell  : {stats['B-cell']}")
    print(f"[SUCCESS] Total Execution Time     : {total_time}")
    print(f"[INFO] Master Library Path : {os.path.basename(csv_path)}")
    print("="*80 + "\n")

if __name__ == "__main__":
    run_phase1da_identification()
