"""
phase2_common.py

Shared utilities for the Phase II vaccine-candidate screening pipeline
(Steps 2A-2D).

WHY THIS FILE EXISTS:
Every bug found while debugging Steps 2A-2D this session came from the
same root cause: path-resolution and FASTA-parsing logic was copy-pasted
into each script, and the copies drifted out of sync with each other
(one script got fixed, the others didn't). Centralizing that logic here
means it only needs to be correct once. If you need to change how the
project root is found, or how FASTA headers are parsed, change it here
and all four steps pick it up automatically.

Each Step_2X script keeps a small (~10 line) local bootstrap that finds
this file and adds it to sys.path -- that bootstrap is intentionally
kept in each script rather than moved here, since a module can't be
used to locate itself.
"""

import os
import re
import sys
import csv
import glob

RESEARCH_ANCHOR = "Research"
PHASE1G_DIR_RELATIVE = os.path.join("Step_Outputs", "Phase1", "Phase1G")
PHASE1G_FASTA_GLOB = "Phase1G_FinalConstruct_*.fasta"
PHASE1G_CSV_GLOB = "Phase1G_FinalConstruct_*.csv"

# NOTE: the per-script `_bootstrap_find_research_root()` copies in each
# Step_2X file are intentionally NOT replaced by a call to this function --
# a script needs to know where this module lives before it can import it,
# so that bootstrap has to stay local. This is the version everything
# AFTER import uses.


def latest_file(folder, pattern="*"):
    """
    Returns the path to the most recently created file matching `pattern`
    (glob syntax) in `folder`, or None if the folder is missing or no file
    matches. Callers must check for None -- this never raises.
    """
    if not os.path.isdir(folder):
        return None
    candidates = glob.glob(os.path.join(folder, pattern))
    if not candidates:
        return None
    candidates.sort(key=os.path.getctime)
    return candidates[-1]


def phase1g_fasta_path(project_root):
    """
    Absolute path to the newest Phase 1G final-construct FASTA, or None if
    none exists yet. Every Step 1G rerun writes a fresh timestamped file
    (Phase1G_FinalConstruct_<timestamp>.fasta) directly under
    Step_Outputs/Phase1G/ -- there is no "Phase 1/" prefix and no fixed
    filename, so this always picks up the latest run instead of a
    hardcoded name that goes stale the moment Step 1G reruns.
    """
    folder = os.path.join(project_root, PHASE1G_DIR_RELATIVE)
    return latest_file(folder, PHASE1G_FASTA_GLOB)


def phase1g_boundary_map(project_root, construct_name=None):
    """
    Returns the Boundary_Map string for a Phase 1G construct, or None.

    Boundary_Map records where every architectural element sits in the final
    sequence -- "adjuvant[0-45];EAAAK[45-50];MHC-I:MIVGGLIGL[50-59];..." in
    0-based half-open coordinates. Reading it lets downstream steps derive
    the construct's regions instead of hardcoding residue ranges that go
    stale the moment Step 1G reruns with a different epitope set.

    construct_name, if given, is matched against the Construct_ID column
    (substring, either direction, so it tolerates the "| length=484" suffix
    Step 2A carries in its Variant column). With no match, or no name given,
    the first row is returned -- 1G writes one construct per file.
    """
    path = latest_file(os.path.join(project_root, PHASE1G_DIR_RELATIVE), PHASE1G_CSV_GLOB)
    if path is None:
        return None
    try:
        with open(path, newline="") as fh:
            rows = [r for r in csv.DictReader(fh) if r.get("Boundary_Map")]
    except OSError:
        return None
    if not rows:
        return None
    if construct_name:
        key = construct_name.split("|")[0].strip()
        for r in rows:
            cid = (r.get("Construct_ID") or "").strip()
            if cid and (cid in key or key in cid):
                return r["Boundary_Map"]
    return rows[0]["Boundary_Map"]


def load_multi_fasta(fasta_path):
    """
    Parses a multi-sequence FASTA file into an ordered dict of
    {header_text: sequence}. Header text is everything after '>' on the
    header line, stripped of whitespace and otherwise left untouched
    (metadata like "| length=138" is preserved as-is so it matches
    exactly whatever Step 2A wrote into its CSV's Variant column).

    Returns None if fasta_path is None or the file doesn't exist (caller
    decides how to report that), or an (possibly empty) dict otherwise.
    """
    if not fasta_path or not os.path.isfile(fasta_path):
        return None

    records = {}
    current_name = None
    current_seq_lines = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records[current_name] = "".join(current_seq_lines).upper()
                current_name = line[1:].strip()
                current_seq_lines = []
            else:
                current_seq_lines.append(line)
        if current_name is not None:
            records[current_name] = "".join(current_seq_lines).upper()

    return records


def sanitize_variant_name(name):
    """
    Variant names carry length metadata and spaces, which is fragile to
    put directly into filenames. This strips the metadata and swaps
    spaces for underscores, producing one stable filename-safe token
    that Steps 2B, 2C, and 2D all build filenames from identically.

    Two header shapes occur in practice and BOTH must reduce to the same
    token:
        "Vax_Var1_fb3faeb6 | length=138"   (pipe-separated, older runs)
        "Vax_Final_a6aaa7d2 length=496"    (no pipe -- Phase 1G writes
                                            SeqRecord.description like this)
    The pipe-only split silently left "length=496" in the name for the
    second shape, so filenames carried an "=". That breaks external tools
    that parse `key=value` command lines: phenix.geometry_minimization
    died with 'improper definition name', and APBS with 'Ignoring
    undefined keyword 496_pH7.0.pqr'. Any "=" is therefore removed here
    rather than patched per-tool, since every new tool would hit it.
    """
    token = name.split("|")[0].strip()
    # Drop a trailing "length=NNN" (with or without a preceding pipe).
    token = re.sub(r"\s*length\s*=\s*\d+\s*$", "", token, flags=re.IGNORECASE).strip()
    # Safety net: no "=" may survive into a filename, whatever its source.
    return token.replace(" ", "_").replace("=", "_")


def lookup_sequence(variants, full_name):
    """
    Looks up a sequence by trying the exact recorded name first, then
    falling back to the name with any '| metadata' stripped, in case
    the FASTA header and a CSV's Variant column ever drift slightly.
    """
    if variants is None:
        return None
    clean_name = full_name.split("|")[0].strip()
    return variants.get(full_name) or variants.get(clean_name)


def get_winner_from_filtered_csv(filtered_dir):
    """
    Finds the most recent Filtered CSV written by Step 2A and returns
    (winner_row_dict, csv_path) for the rank-1 (first) row -- Step 2A
    always writes this file sorted best-first by Stability Index.

    Returns (None, None) after printing a diagnostic message if
    anything is missing, so every downstream step reports failures
    the same way instead of each reinventing this check slightly
    differently (which is how Steps 2B/2C/2D ended up with three
    subtly different error-handling styles before this rewrite).
    """
    if not os.path.isdir(filtered_dir):
        print(f"[ERROR] Filtered output folder does not exist: {filtered_dir}")
        print("[ERROR] Run Step 2A first so this folder gets created.")
        return None, None

    csv_files = sorted([f for f in os.listdir(filtered_dir) if f.endswith(".csv")])
    if not csv_files:
        print(f"[ERROR] No filtered candidate CSVs found in: {filtered_dir}")
        print("[ERROR] Step 2A ran but produced zero viable candidates -- check Rejection_Reasons in its Raw log.")
        return None, None

    csv_path = os.path.join(filtered_dir, csv_files[-1])
    with open(csv_path, 'r') as f:
        reader = list(csv.DictReader(f))

    if not reader:
        print(f"[ERROR] Filtered CSV is empty: {csv_path}")
        return None, None

    return reader[0], csv_path


def print_banner(text, width=115):
    print("\n" + "=" * width)
    print(f"{text:^{width}}")
    print("=" * width)


def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"


# ---------------------------------------------------------------------------
# Contextual-only physicochemical helpers not provided directly by Biopython
# ---------------------------------------------------------------------------

def compute_aliphatic_index(seq):
    """
    Ikai (1980) aliphatic index -- relative volume occupied by aliphatic
    side chains (Ala, Val, Ile, Leu). Biopython's ProteinAnalysis does
    NOT implement this, so it's computed by hand here.

    Per methodology this is a CONTEXTUAL metric only -- it must never be
    used as a rejection criterion.
    """
    length = len(seq)
    if length == 0:
        return 0.0
    ala = seq.count('A') / length * 100
    val = seq.count('V') / length * 100
    ile = seq.count('I') / length * 100
    leu = seq.count('L') / length * 100
    return ala + 2.9 * val + 3.9 * (ile + leu)


# ExPASy ProtParam's mammalian (in vitro) N-end rule half-life table,
# in hours, keyed by N-terminal residue (Bachmair et al. 1986;
# Rogers et al. 1986).
_MAMMALIAN_HALF_LIFE_HOURS = {
    'A': 4.4, 'R': 1.0, 'N': 1.4, 'D': 1.1, 'C': 1.2,
    'Q': 0.8, 'E': 1.0, 'G': 30.0, 'H': 3.5, 'I': 20.0,
    'L': 5.5, 'K': 1.3, 'M': 30.0, 'F': 1.1, 'P': 20.0,
    'S': 1.9, 'T': 7.2, 'W': 2.8, 'Y': 2.8, 'V': 100.0,
}


def estimate_half_life_hours(seq):
    """
    Rough N-end-rule half-life estimate (mammalian, in vitro) based on
    the sequence's N-terminal residue, matching the convention used by
    ExPASy ProtParam. This is a approximation, not a real assay result.

    Per methodology this is a CONTEXTUAL metric only -- it must never be
    used as a rejection criterion.
    """
    if not seq:
        return None
    return _MAMMALIAN_HALF_LIFE_HOURS.get(seq[0])
