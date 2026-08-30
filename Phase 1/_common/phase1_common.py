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


# =============================================================================
# HIV SUBUNIT RANGES -- one definition, read by 1A's provenance check and 1De.
#
# WHY THIS EXISTS: Phase 1A defines its four HIV targets as four Entrez QUERY
# STRINGS ("gp120 AND ...", "gp41 AND ...", "p17 AND ...", "p24 AND ..."), and
# stores whatever record NCBI returns whole, labelled with the query name. For
# CRF01_AE, NCBI returns POLYPROTEINS: HIV_gp120_Var_01 and HIV_gp41_Var_01 are
# the same 856-aa Env (YES72107.1), and HIV_p17_Var_01 and HIV_p24_Var_01 are
# the same 1437-aa Gag-Pol (YES72110.1). Every epitope sliced from them
# therefore inherited its query's label regardless of where in the polyprotein
# it actually lies -- see the Phase 1A provenance-correction report.
#
# PROVENANCE OF THE NUMBERS: the GenBank records carry NO mat_peptide features,
# so boundaries come from a global BLOSUM62 alignment of each Var_01 sequence
# against the annotated HXB2 reference, and the reference's own UniProt CHAIN
# features mapped through that alignment:
#   Env      YES72107.1 (856 aa)  vs  P04578 ENV_HV1H2 (856 aa)
#            CHAIN 33..511  "Surface protein gp120"      -> Var_01  32..504
#            CHAIN 512..856 "Transmembrane protein gp41" -> Var_01 505..856
#   Gag-Pol  YES72110.1 (1437 aa) vs  P04585 POL_HV1H2 (1435 aa)
#            CHAIN 2..132   "Matrix protein p17"         -> Var_01   2..135
#            CHAIN 133..363 "Capsid protein p24"         -> Var_01 136..366
# Independently corroborated by sequence landmarks in the stored Var_01 files:
# the Env furin site + fusion peptide "...RVVERPKR | AVGIGAMIFGF" puts gp41 at
# 505, and the Gag MA/CA junction "...SQNY | PIVQ" puts p24 at 136, ending at
# "...KARVL" (366). Both methods agree exactly.
#
# 1-BASED AND INCLUSIVE, matching UniProt/GenBank convention. Slice with
# seq[start-1:end].
# =============================================================================
HIV_SUBUNIT_RANGES = {
    "HIV_gp120": {"parent": "Env",     "accession": "YES72107.1", "start": 32,  "end": 504},
    "HIV_gp41":  {"parent": "Env",     "accession": "YES72107.1", "start": 505, "end": 856},
    "HIV_p17":   {"parent": "Gag-Pol", "accession": "YES72110.1", "start": 2,   "end": 135},
    "HIV_p24":   {"parent": "Gag-Pol", "accession": "YES72110.1", "start": 136, "end": 366},
}

# Targets whose Var_01 file IS the mature protein -- no slicing needed.
SINGLE_PROTEIN_TARGETS = ("Mpox_A35R", "Mpox_B5R", "Mpox_L1R")


def subunit_of(target, position_1based):
    """
    Given a 1-based offset into a target's Var_01 parent record, returns the
    name of the HIV subunit that offset actually falls in, or None if the
    target is not one of the polyprotein-derived HIV four (Mpox targets are
    already mature proteins, so their label is always correct).

    Returns "Gag_downstream" for Gag-Pol offsets past p24 (p2/p7/p6/pol) --
    those are real regions, just not one of this study's four HIV antigens.
    """
    if target in SINGLE_PROTEIN_TARGETS:
        return target
    info = HIV_SUBUNIT_RANGES.get(target)
    if info is None:
        return None
    siblings = [(n, d) for n, d in HIV_SUBUNIT_RANGES.items() if d["parent"] == info["parent"]]
    for name, d in siblings:
        if d["start"] <= position_1based <= d["end"]:
            return name
    if info["parent"] == "Gag-Pol" and position_1based > 366:
        return "Gag_downstream"
    return None


# =============================================================================
# LANDMARK-BASED SUBUNIT BOUNDARIES FOR NON-Var_01 RECORDS.
#
# HIV_SUBUNIT_RANGES above is exact but applies ONLY to the Var_01 records it
# was aligned against. Epitopes selected from other variants live in records of
# different length (Env variants here run 854-868 aa), so those fixed offsets
# do not transfer. Rather than align every variant, locate the SAME cleavage
# landmarks the alignment confirmed on Var_01:
#
#   Env      gp120 | gp41 at the host-furin site immediately followed by the
#            gp41 fusion peptide: "...RVVERPKR | AVGIGAMIFGF". On Var_01 this
#            puts gp41 at 505 -- identical to the alignment-derived boundary.
#   Gag-Pol  p17 | p24 at the MA/CA junction "...SQNY | PIVQ", and p24's C
#            terminus at "...KARVL". On Var_01 these give 136 and 366 --
#            again identical to the alignment-derived boundaries.
#
# Both landmarks are the actual protease/furin recognition sites, so they track
# indels correctly in a way fixed offsets cannot. Returns None when a landmark
# is absent, and callers must treat that as UNRESOLVED rather than guessing.
# =============================================================================
def subunit_boundaries_by_landmark(seq, parent):
    """
    Returns a dict of {subunit_name: (start, end)} 1-based inclusive for the
    given record, or None if the defining landmark cannot be found.
    """
    if parent == "Env":
        fp = seq.find("AVGIG")          # gp41 fusion peptide, first residue of gp41
        if fp == -1:
            return None
        return {"HIV_gp120": (1, fp), "HIV_gp41": (fp + 1, len(seq))}

    if parent == "Gag-Pol":
        j = seq.find("SQNYPIVQ")        # MA/CA junction: p24 starts at the P
        if j == -1:
            return None
        p24_start = j + 5               # 1-based position of 'P' in PIVQ
        ca_end = seq.find("KARVL", p24_start)
        p24_end = (ca_end + 5) if ca_end != -1 else None   # p24 ends at the L
        if p24_end is None:
            return None
        return {"HIV_p17": (1, p24_start - 1), "HIV_p24": (p24_start, p24_end)}

    return None
