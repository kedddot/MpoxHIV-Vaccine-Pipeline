"""
phase1_common.py

Shared utilities for the Phase I epitope-selection pipeline (Steps 1A-1G).

WHY THIS FILE EXISTS:
Every Phase 1 script computed its project root as
os.path.join(script_dir, "..", "..", "..") -- three fixed hops. From
Research/Phase 1/STEP X/ that resolves to the folder ABOVE Research
itself, so every step's output was silently written outside the repo.
Phase 2 solved the same problem with an anchor-folder walk instead of a
fixed hop count (see Phase 2/_common/phase2_common.py); this file ports
that fix to Phase 1 so both phases resolve the same way.
"""

import os
import sys

RESEARCH_ANCHOR = "Research"


def resolve_project_root(script_file):
    """
    Walk upward from a script's own location until a folder literally
    named "Research" is found. Anchoring on a named folder (rather than
    a fixed hop count) keeps working regardless of how deep a given
    STEP script sits -- which is exactly what broke the fixed "../../.."
    version of this logic.
    """
    script_dir = os.path.dirname(os.path.abspath(script_file))
    current = script_dir
    while os.path.basename(current) != RESEARCH_ANCHOR:
        parent = os.path.dirname(current)
        if parent == current:
            print(f"\n[FATAL ERROR] Could not locate a '{RESEARCH_ANCHOR}' anchor folder above: {script_dir}")
            sys.exit(1)
        current = parent
    return current


def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"


def print_banner(text, width=90):
    print("\n" + "=" * width)
    print(f"{text:^{width}}")
    print("=" * width)


def latest_file(folder, suffix=".csv"):
    """
    Returns the path to the most recently created file with the given
    suffix in `folder`, or None if the folder is missing or empty.
    Callers must check for None -- this never raises on a missing dir.
    """
    if not os.path.isdir(folder):
        return None
    candidates = [f for f in os.listdir(folder) if f.endswith(suffix)]
    if not candidates:
        return None
    candidates.sort(key=lambda f: os.path.getctime(os.path.join(folder, f)))
    return os.path.join(folder, candidates[-1])
