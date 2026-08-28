import os
import re
import sys
import csv
import json
import time
import shutil
import zipfile
import tempfile
import glob
import argparse
import itertools
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

# Step 2C already has a tested AlphaFold mmCIF pLDDT parser (it reads
# `_ma_qa_metric_local` first and falls back to the CA B-factor column).
# Import it rather than writing a second one that can drift out of sync.
_STEPC_DIR = os.path.join(_PROJECT_ROOT, "Phase 2", "STEP C")
if _STEPC_DIR not in sys.path:
    sys.path.insert(0, _STEPC_DIR)
from Phase2C_solubilityAnalysis import parse_plddt_per_residue, region_mean_plddt

# =============================================================================
# ALPHAFOLD SERVER INTEGRATION (real, but web-submission-based -- see notes)
#
# AlphaFold3 requires an NVIDIA/CUDA GPU and does not run on a Mac. AlphaFold
# Server (alphafoldserver.com) is the official free DeepMind web version and
# has NO public API -- submission is browser-only. Per Sec. II.B only the
# VACCINE MONOMER is predicted here; the TLR-2/TLR-4 structures are taken
# from RCSB and the complexes are built by HADDOCK in Phase III (Sec. III.A).
# So this step is split into two phases you run separately:
#
#   1) python Phase2B_secondaryAndTertiary.py --prepare
#      Downloads the TLR-2 (6NIG) and TLR-4 (8WTA) EXPERIMENTAL structures
#      from RCSB automatically (Sec. II.B: they are "taken from the RCSB
#      Protein Data Bank", not predicted), then writes out the sequences to
#      paste into AlphaFold Server (vaccine monomer only), PsiPred, and
#      NetSurfP-3.0 -- 3 manual jobs.
#
#   2) [ you manually submit the jobs on alphafoldserver.com /
#        bioinf.cs.ucl.ac.uk/psipred / services.healthtech.dtu.dk NetSurfP-3.0,
#        wait for them to finish, and download each result ]
#
#   3) python Phase2B_secondaryAndTertiary.py --import monomer <path_to_download>
#      python Phase2B_secondaryAndTertiary.py --import psipred <path_to_.ss2>
#      python Phase2B_secondaryAndTertiary.py --import netsurfp <path_to_.csv>
#      Unzips each AlphaFold download, locates the top-ranked model (.cif),
#      and places it where Step 2C/2D expect it. PsiPred/NetSurfP results
#      are copied in as-is and parsed when Step 2B runs normally.
#
#   4) python Phase2B_secondaryAndTertiary.py
#      Runs normally once all 3 imports are done: computes the Biopython
#      Chou-Fasman baseline AND the real PsiPred/NetSurfP-3.0 secondary
#      structure calls, and finishes the step.
#
# NOTE: AlphaFold Server's real output format is mmCIF (.cif), not legacy
# PDB. Step 2C and 2D need one matching edit each to expect .cif -- see the
# comments at the top of those files.
# =============================================================================


def fetch_rcsb_entity_sequence(pdb_id, entity_id="1", residue_range=None):
    """
    Fetches a polymer entity's sequence directly from RCSB PDB's REST API,
    instead of a hand-typed literal -- avoids both transcription errors
    and, for chimeric crystallization constructs, accidentally including a
    fusion/chaperone partner's sequence in what's meant to be a pure
    receptor sequence.

    residue_range, if given, is an inclusive (start, end) 1-based slice
    into the ENTITY's own sequence (not UniProt numbering), used to trim
    off a fusion partner. Example: 6NIG's resolved TLR2 entity is a
    TLR2-ectodomain/VLR-B chimera used for crystallization -- only entity
    residues 1-507 map to TLR2 (UniProt O60603); 509-576 are the VLR-B
    chaperone (UniProt Q2YE02), confirmed via RCSB's own SIFTS alignment.
    """
    import requests

    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    seq = data["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", "").strip()
    if not seq:
        raise ValueError(f"RCSB returned an empty sequence for {pdb_id} entity {entity_id}")
    if residue_range:
        start, end = residue_range
        seq = seq[start - 1:end]
    return seq


def fetch_tlr2_ectodomain_sequence():
    """
    RETAINED FOR PHASE III DOCKING PREP -- not used by --prepare anymore
    (Sec. II.B takes the TLR structures from RCSB rather than predicting
    them). Kept because the 6NIG residue range below is the record of
    which part of that chimeric entity is actually TLR2, which is needed
    when stripping the VLR-B fusion partner before docking.
    
    Human TLR2 ectodomain, from RCSB 6NIG (Crystal structure of the human
    TLR2-Diprovocim complex), entity 1, residues 1-507 -- the TLR2 portion
    of that entity's TLR2/VLR-B crystallization chimera (see
    fetch_rcsb_entity_sequence's docstring). This is the resolved
    ectodomain only, not the full-length UniProt O60603 sequence (which
    also includes the transmembrane and TIR domains, irrelevant to a
    surface-docking model).
    """
    return fetch_rcsb_entity_sequence("6NIG", "1", residue_range=(1, 507))


def fetch_tlr4_ectodomain_sequence():
    """
    Human TLR4 ectodomain, from RCSB 8WTA (Cryo-EM Structure of Human
    TLR4/MD-2/DLAM3 Complex), entity 1 -- confirmed via RCSB's own SIFTS
    alignment to be a clean (non-chimeric) TLR4 sequence (UniProt O00206,
    residues 27-631), so no trimming is needed here the way 6NIG's TLR2
    entity requires.
    """
    return fetch_rcsb_entity_sequence("8WTA", "1")


def _parse_psipred_ss2(path):
    """
    Parses PsiPred's `.ss2` output: one header comment line, then one row
    per residue as whitespace-separated `pos aa ss coil helix sheet`.
    Tallies the single-letter SS call (C/H/E) per residue into fractions.
    Verify against your actual downloaded .ss2 the first time you use
    this -- the 3-column layout has been PsiPred's stable format for
    years, but hasn't been tested here against a live file.
    """
    counts = {'C': 0, 'H': 0, 'E': 0}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 3:
                continue
            ss = cols[2].upper()
            if ss in counts:
                counts[ss] += 1
    total = sum(counts.values())
    if total == 0:
        raise ValueError(f"No parseable SS calls found in {path}")
    return {
        "helix_fraction": counts['H'] / total,
        "sheet_fraction": counts['E'] / total,
        "coil_fraction": counts['C'] / total,
        "n_residues": total,
    }


def _parse_raptorx_horiz(path):
    """
    Parses RaptorX's PsiPred-style `.horiz` secondary-structure output:
    blocks of aligned lines, one of which per block starts with "Pred:"
    holding the 3-state (H/E/C) call for that block. Concatenates all
    "Pred:" lines into one full-length SS string and tallies fractions.
    Verify against your actual downloaded .horiz the first time you use
    this -- RaptorX has more than one output flavor; this targets the
    PsiPred-compatible one.
    """
    ss_chunks = []
    with open(path) as f:
        for line in f:
            if line.startswith("Pred:"):
                ss_chunks.append(line[len("Pred:"):].strip())
    full_ss = "".join(ss_chunks).upper()
    if not full_ss:
        raise ValueError(f"No 'Pred:' lines found in {path}")
    total = len(full_ss)
    return {
        "helix_fraction": full_ss.count('H') / total,
        "sheet_fraction": full_ss.count('E') / total,
        "coil_fraction": full_ss.count('C') / total,
        "n_residues": total,
    }


def _parse_netsurfp_csv(path):
    """
    Parses NetSurfP-3.0's per-residue CSV export.

    SUBSTITUTED FOR RAPTORX (deviation #16): raptorx.uchicago.edu returned
    HTTP 500 on every endpoint including its root -- its whole Django app
    was down, leaking a debug traceback. NetSurfP-3.0 (DTU, 2022) was
    chosen over the other live options because it is built on ESM-1b
    embeddings and is genuinely independent of the PSIPRED lineage, which
    is the entire point of running a second predictor (S4Pred, the
    convenient alternative, comes from the same group as PSIPRED).

    NetSurfP-3.0's CSV carries one row per residue with a 3-state call in
    a `q3` column. Column naming has varied between releases, so this
    matches the q3 column case-insensitively and falls back to any column
    whose values are drawn only from {H,E,C}. Verify against your actual
    download the first time -- same caution as the PsiPred/RaptorX parsers.
    """
    import csv as _csv
    with open(path, newline="") as f:
        rows = list(_csv.DictReader(f))
    if not rows:
        raise ValueError(f"No rows parsed from {path}")

    fieldnames = rows[0].keys()
    q3_col = next((c for c in fieldnames if c and c.strip().lower() == "q3"), None)
    if q3_col is None:
        # Fallback: first column whose values are purely 3-state SS letters.
        for c in fieldnames:
            vals = {(r.get(c) or "").strip().upper() for r in rows}
            vals.discard("")
            if vals and vals <= {"H", "E", "C"}:
                q3_col = c
                break
    if q3_col is None:
        raise ValueError(
            f"Could not find a 3-state (H/E/C) column in {path}. "
            f"Columns present: {list(fieldnames)}"
        )

    counts = {"H": 0, "E": 0, "C": 0}
    for r in rows:
        ss = (r.get(q3_col) or "").strip().upper()
        if ss in counts:
            counts[ss] += 1
    total = sum(counts.values())
    if total == 0:
        raise ValueError(f"No parseable H/E/C calls in column '{q3_col}' of {path}")
    return {
        "helix_fraction": counts["H"] / total,
        "sheet_fraction": counts["E"] / total,
        "coil_fraction": counts["C"] / total,
        "n_residues": total,
    }


def _dssp_from_structure(cif_path):
    """
    Assigns secondary structure directly from the AlphaFold model's 3D
    coordinates using DSSP (pydssp -- pure Python, no mkdssp binary
    needed).

    This is an ASSIGNMENT from predicted coordinates, NOT an independent
    sequence-based prediction: it inherits whatever AlphaFold got wrong,
    so it does not substitute for the second predictor Sec. II.B asks
    for. It is reported as a third, structure-derived column because it
    is self-consistent with the very structure Steps 2C/2D go on to use,
    and it costs nothing extra.

    Returns None (with a printed reason) rather than raising, so a DSSP
    problem can never block Step 2B -- it is supplementary, not required.
    """
    try:
        import numpy as np
        import pydssp
        from Bio.PDB.MMCIFParser import MMCIFParser

        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", cif_path)

        # pydssp needs backbone N, CA, C, O per residue, in that order.
        coords = []
        for residue in structure[0].get_residues():
            try:
                coords.append([residue["N"].coord, residue["CA"].coord,
                               residue["C"].coord, residue["O"].coord])
            except KeyError:
                continue  # skip residues missing backbone atoms
        if not coords:
            print("[WARN] DSSP: no complete backbone residues found -- skipping.")
            return None

        ss = pydssp.assign(np.array(coords, dtype=float), out_type="c3")
        ss_str = "".join(ss).upper().replace("-", "C").replace(" ", "C")
        total = len(ss_str)
        if total == 0:
            return None
        return {
            "helix_fraction": ss_str.count("H") / total,
            "sheet_fraction": ss_str.count("E") / total,
            "coil_fraction": ss_str.count("C") / total,
            "n_residues": total,
        }
    except Exception as e:
        print(f"[WARN] DSSP assignment unavailable ({e}) -- reported as n/a, not a failure.")
        return None


# =============================================================================
# ALPHAFOLD3 FOLD-CONFIDENCE ANALYSIS
#
# Sec. II.B only asks for the 3D model, so an earlier version of this step
# kept model_0 and deleted the rest of the AlphaFold Server download. That
# threw away every confidence signal AlphaFold produces: the other four
# ranked models, all five summary_confidences_*.json (pTM, iptm,
# fraction_disordered, ranking_score) and all five 484x484 PAE matrices.
#
# That matters because Step 2D's MolProbity verdict scores LOCAL
# STEREOCHEMISTRY -- bond geometry, clashes, rotamers, Ramachandran. A
# clean-geometry random coil passes it. Nothing downstream was looking at
# whether the fold itself was determined, which is the actual precondition
# for docking in Phase III (Sec. III.A): docking an undetermined region
# means docking one arbitrary conformer, and the result does not reproduce
# when a different ranked model is used.
#
# These functions read those files (they are archived under
# AlphaFold_Raw/<job>/ now) and report four independent confidence axes:
#   pLDDT      -- per-residue local confidence, from each model's mmCIF
#   pTM        -- whole-fold reliability, from the summary JSONs
#   PAE        -- confidence in the RELATIVE placement of two regions
#   RMSD       -- how far a region actually moves between ranked models
# Agreement between four independent axes is what makes the conclusion
# defensible; any one of them alone is arguable.
# =============================================================================

# Which cassette each Boundary_Map segment belongs to. Linkers are grouped
# with the cassette they introduce (AAY joins MHC-I epitopes, GPGPG joins
# MHC-II, KK joins B-cell), which is what makes the leading GPGPG/KK fall on
# the correct side of a cassette boundary.
_SEGMENT_TO_REGION = {
    "adjuvant": "Adjuvant_bDefensin3",
    "EAAAK": "Linker_EAAAK",
    "AAY": "Cassette_MHC_I",
    "MHC-I": "Cassette_MHC_I",
    "GPGPG": "Cassette_MHC_II",
    "MHC-II": "Cassette_MHC_II",
    "KK": "Cassette_Bcell",
    "B-cell": "Cassette_Bcell",
}
_REGION_ORDER = ["Adjuvant_bDefensin3", "Linker_EAAAK", "Cassette_MHC_I",
                 "Cassette_MHC_II", "Cassette_Bcell"]


def parse_boundary_map(boundary_map):
    """
    Turns Phase 1G's Boundary_Map column into contiguous functional regions.

    Input is 0-based half-open, semicolon-separated, e.g.
        "adjuvant[0-45];EAAAK[45-50];MHC-I:MIVGGLIGL[50-59];AAY[59-62];..."

    Returns [(region_name, start, end), ...] as 1-based INCLUSIVE residue
    ranges, matching the numbering used by parse_plddt_per_residue().

    Deriving the regions from 1G instead of hardcoding them means this
    follows the construct automatically if the architecture ever changes.
    """
    spans = {}
    for seg in boundary_map.split(";"):
        seg = seg.strip()
        if not seg:
            continue
        m = re.match(r"^([^\[]+)\[(\d+)-(\d+)\]$", seg)
        if not m:
            raise ValueError(f"Unparseable Boundary_Map segment: {seg!r}")
        label, start0, end0 = m.group(1), int(m.group(2)), int(m.group(3))
        key = label.split(":", 1)[0]
        region = _SEGMENT_TO_REGION.get(key)
        if region is None:
            raise ValueError(
                f"Boundary_Map segment {label!r} has no region mapping. Add it to "
                f"_SEGMENT_TO_REGION -- do not let it fall through silently.")
        lo, hi = start0 + 1, end0          # 0-based half-open -> 1-based inclusive
        if region in spans:
            spans[region] = (min(spans[region][0], lo), max(spans[region][1], hi))
        else:
            spans[region] = (lo, hi)

    regions = [(name, spans[name][0], spans[name][1])
               for name in _REGION_ORDER if name in spans]
    for name in spans:
        if name not in _REGION_ORDER:
            regions.append((name, spans[name][0], spans[name][1]))
    return regions


def longest_confident_run(plddt_table, threshold=70.0):
    """
    Longest contiguous stretch of residues at pLDDT >= threshold.

    Same shape as Step 2C's Aggrescan3D patch logic, and the same 70 cut --
    AlphaFold's own "confident" band, already the gate 2C uses before it
    trusts an aggregation patch.

    Returns (start, end, length) 1-based inclusive, or (None, None, 0).
    """
    best = (None, None, 0)
    run_start = None
    prev = None
    for resi in sorted(plddt_table):
        ok = plddt_table[resi] >= threshold
        if ok and (run_start is None or prev is None or resi != prev + 1):
            run_start = resi
        if ok:
            length = resi - run_start + 1
            if length > best[2]:
                best = (run_start, resi, length)
        else:
            run_start = None
        prev = resi
    return best


def _af3_model_files(raw_dir):
    """Groups an AlphaFold Server download into per-model (cif, summary, full_data)."""
    models = {}
    for path in sorted(glob.glob(os.path.join(raw_dir, "*"))):
        base = os.path.basename(path)
        m = re.search(r"_model_(\d+)\.cif$", base)
        if m:
            models.setdefault(int(m.group(1)), {})["cif"] = path
            continue
        m = re.search(r"_summary_confidences_(\d+)\.json$", base)
        if m:
            models.setdefault(int(m.group(1)), {})["summary"] = path
            continue
        m = re.search(r"_full_data_(\d+)\.json$", base)
        if m:
            models.setdefault(int(m.group(1)), {})["full"] = path
    return dict(sorted(models.items()))


def _pae_block_mean(pae, res_index, rows, cols):
    """
    Mean PAE over the block linking two residue ranges.

    PAE is asymmetric (pae[i][j] is the error in residue j's position when
    the frame is aligned on residue i), so for a cross-region figure both
    off-diagonal blocks are averaged together -- reporting only one direction
    would understate or overstate the coupling depending on which was picked.
    """
    ri = [res_index[r] for r in rows if r in res_index]
    ci = [res_index[c] for c in cols if c in res_index]
    if not ri or not ci:
        return None
    total = 0.0
    for i in ri:
        row = pae[i]
        for j in ci:
            total += row[j]
    n = len(ri) * len(ci)
    if rows is not cols and ri != ci:
        for i in ci:
            row = pae[i]
            for j in ri:
                total += row[j]
        n *= 2
    return total / n


def _region_ca_coords(cif_path, start, end):
    """CA coordinates for a 1-based inclusive residue range, ordered by residue id."""
    from Bio.PDB import MMCIFParser
    import warnings
    from Bio import BiopythonWarning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        structure = MMCIFParser(QUIET=True).get_structure("af3", cif_path)
    atoms = {}
    for chain in structure[0]:
        for residue in chain:
            resi = residue.get_id()[1]
            if start <= resi <= end and "CA" in residue:
                atoms[resi] = residue["CA"]
    return [atoms[k] for k in sorted(atoms)]


def _cross_model_rmsd(model_paths, start, end):
    """
    Pairwise CA RMSD over a residue range, after optimal superposition on
    THAT RANGE ONLY, across every pair of ranked models.

    Uses Bio.PDB.Superimposer deliberately. A hand-rolled Kabsch that applies
    the rotation on the wrong side (A @ R.T instead of A @ R) does not fail
    loudly -- it returns plausible-looking wrong RMSDs, which is exactly the
    kind of error this analysis exists to rule out.

    Returns (mean, min, max, n_pairs) in Angstroms, or (None, None, None, 0).
    """
    from Bio.PDB import Superimposer
    coords = {}
    for idx, path in model_paths.items():
        ca = _region_ca_coords(path, start, end)
        if ca:
            coords[idx] = ca

    values = []
    for a, b in itertools.combinations(sorted(coords), 2):
        fixed, moving = coords[a], coords[b]
        if len(fixed) != len(moving) or not fixed:
            continue
        sup = Superimposer()
        sup.set_atoms(fixed, list(moving))
        values.append(sup.rms)
    if not values:
        return None, None, None, 0
    return sum(values) / len(values), min(values), max(values), len(values)


def analyze_af3_confidence(raw_dir, regions, plddt_threshold=70.0):
    """
    Reads a whole AlphaFold Server download and returns the four confidence
    axes, per model and per region. Never raises on a partial download --
    missing files just leave their fields None, because this is diagnostic
    reporting and must not be able to block the step.
    """
    models = _af3_model_files(raw_dir)
    if not models:
        return None

    result = {"raw_dir": raw_dir, "n_models": len(models), "models": [],
              "regions": [], "region_pairs": [], "whole": {}}

    # ---- per-model global confidence -------------------------------------
    for idx, files in models.items():
        entry = {"model": idx, "ptm": None, "iptm": None,
                 "fraction_disordered": None, "ranking_score": None, "has_clash": None}
        if "summary" in files:
            try:
                s = json.load(open(files["summary"]))
                entry.update({k: s.get(k) for k in
                              ("ptm", "iptm", "fraction_disordered", "ranking_score", "has_clash")})
            except Exception as e:
                print(f"[WARN] Could not read {os.path.basename(files['summary'])}: {e}")
        result["models"].append(entry)

    # ---- pLDDT on the top-ranked model (model_0) -------------------------
    top = models.get(0) or models[min(models)]
    plddt = None
    if "cif" in top:
        try:
            plddt = parse_plddt_per_residue(top["cif"])
        except Exception as e:
            print(f"[WARN] Could not parse pLDDT from the top-ranked model: {e}")

    if plddt:
        values = list(plddt.values())
        n_conf = sum(1 for v in values if v >= plddt_threshold)
        run_start, run_end, run_len = longest_confident_run(plddt, plddt_threshold)
        result["whole"] = {
            "n_residues": len(values),
            "mean_plddt": sum(values) / len(values),
            "n_plddt_ge_threshold": n_conf,
            "pct_plddt_ge_threshold": 100.0 * n_conf / len(values),
            "pct_plddt_lt_50": 100.0 * sum(1 for v in values if v < 50.0) / len(values),
            "max_plddt": max(values),
            "longest_confident_run": (run_start, run_end),
            "longest_confident_run_len": run_len,
        }

    # ---- PAE blocks from the top-ranked model ----------------------------
    pae = res_index = None
    if "full" in top:
        try:
            full = json.load(open(top["full"]))
            pae = full.get("pae")
            token_res = full.get("token_res_ids")
            if pae and token_res:
                # First token per residue id -- for a monomer this is 1:1, but
                # do not assume it.
                res_index = {}
                for i, r in enumerate(token_res):
                    res_index.setdefault(int(r), i)
        except Exception as e:
            print(f"[WARN] Could not read PAE from {os.path.basename(top['full'])}: {e}")

    if pae:
        flat = [v for row in pae for v in row]
        result["whole"]["max_pae"] = max(flat)
        result["whole"]["mean_pae"] = sum(flat) / len(flat)

    model_cifs = {i: f["cif"] for i, f in models.items() if "cif" in f}

    # ---- per-region -------------------------------------------------------
    for name, start, end in regions:
        row = {"region": name, "start": start, "end": end, "length": end - start + 1,
               "mean_plddt": None, "pct_plddt_ge_threshold": None,
               "intra_pae": None, "rmsd_mean": None, "rmsd_min": None,
               "rmsd_max": None, "rmsd_pairs": 0}
        if plddt:
            try:
                row["mean_plddt"] = region_mean_plddt(plddt, start, end)
            except ValueError:
                pass
            in_region = [v for r, v in plddt.items() if start <= r <= end]
            if in_region:
                row["pct_plddt_ge_threshold"] = (
                    100.0 * sum(1 for v in in_region if v >= plddt_threshold) / len(in_region))
        if pae and res_index:
            rng = list(range(start, end + 1))
            row["intra_pae"] = _pae_block_mean(pae, res_index, rng, rng)
        mean, lo, hi, n = _cross_model_rmsd(model_cifs, start, end)
        row.update({"rmsd_mean": mean, "rmsd_min": lo, "rmsd_max": hi, "rmsd_pairs": n})
        result["regions"].append(row)

    # ---- cross-region PAE -------------------------------------------------
    if pae and res_index:
        for (n1, s1, e1), (n2, s2, e2) in itertools.combinations(regions, 2):
            result["region_pairs"].append({
                "region_a": n1, "region_b": n2,
                "inter_pae": _pae_block_mean(pae, res_index,
                                             list(range(s1, e1 + 1)), list(range(s2, e2 + 1))),
            })
    return result


def classify_docking_readiness(conf, plddt_threshold=70.0, ptm_threshold=0.5,
                               rmsd_threshold=2.0, min_domain_len=30,
                               whole_construct_pct=70.0):
    """
    Turns the four axes into the one statement Phase III needs: what, if
    anything, in this model can legitimately be docked.

    THRESHOLDS ARE EXTERNAL CONVENTION, NOT FROM THE MANUSCRIPT -- the paper
    specifies no fold-confidence criterion at all. Recorded here exactly as
    Step 2D's secondary QC flags are, so a reader can see the provenance:
      pLDDT >= 70   AlphaFold's own "confident" band; already Step 2C's gate
      pTM   >= 0.5  the usual "the global fold is probably right" line
      RMSD  <= 2.0 A  conventional "same structure" agreement between models

    Returns (readiness, regions_string, notes).
      WHOLE_CONSTRUCT -- global fold determined; dock the construct as-is
      DOMAIN_ONLY     -- only listed sub-domains are determined; scope to them
      NOT_READY       -- nothing in the model supports a docking claim
    """
    if not conf:
        return "NOT_ASSESSED", "", "no AlphaFold confidence data available"

    whole = conf.get("whole") or {}
    ptms = [m["ptm"] for m in conf.get("models", []) if m.get("ptm") is not None]
    max_ptm = max(ptms) if ptms else None
    pct_conf = whole.get("pct_plddt_ge_threshold")

    if (max_ptm is not None and max_ptm >= ptm_threshold
            and pct_conf is not None and pct_conf >= whole_construct_pct):
        return ("WHOLE_CONSTRUCT", f"1-{whole.get('n_residues')}",
                f"pTM {max_ptm:.2f} >= {ptm_threshold} and {pct_conf:.1f}% of residues "
                f"at pLDDT >= {plddt_threshold:.0f}")

    dockable, why = [], []
    for row in conf.get("regions", []):
        if row["length"] < min_domain_len:
            continue
        if row["mean_plddt"] is None or row["mean_plddt"] < plddt_threshold:
            continue
        if row["rmsd_mean"] is not None and row["rmsd_mean"] > rmsd_threshold:
            why.append(f"{row['region']} rejected: cross-model RMSD "
                       f"{row['rmsd_mean']:.2f} A > {rmsd_threshold} A")
            continue
        dockable.append(row)

    if dockable:
        spans = ";".join(f"{r['start']}-{r['end']}" for r in dockable)
        note = "; ".join(
            f"{r['region']} {r['start']}-{r['end']}: mean pLDDT {r['mean_plddt']:.1f}, "
            f"cross-model RMSD {r['rmsd_mean']:.2f} A" for r in dockable if r["rmsd_mean"] is not None)
        if max_ptm is not None:
            note += f" (whole-model pTM {max_ptm:.2f} -- global fold NOT determined)"
        return "DOMAIN_ONLY", spans, note

    reason = f"no region >= {min_domain_len} aa reaches mean pLDDT {plddt_threshold:.0f}"
    if why:
        reason += " | " + "; ".join(why)
    return "NOT_READY", "", reason


def write_confidence_report(conf, readiness, regions_str, notes, output_base,
                            project_root, variant_name):
    """Writes Model_Confidence/Step2B_AF3_Confidence_<ts>.csv (long format)."""
    out_dir = os.path.join(output_base, "Model_Confidence")
    os.makedirs(out_dir, exist_ok=True)
    # Sibling steps select their input by lexicographic sort, and '-' sorts
    # before '0' in ASCII -- a dash-formatted same-day name would be silently
    # skipped as "not the latest". Match the no-dash sibling format exactly.
    ts = datetime.now().strftime("%Y%m%d_%H%M")
    path = os.path.join(out_dir, f"Step2B_AF3_Confidence_{ts}.csv")

    def f(v, nd=2):
        return "" if v is None else (round(v, nd) if isinstance(v, float) else v)

    rows = []
    whole = conf.get("whole") or {}
    run = whole.get("longest_confident_run") or (None, None)
    rows.append({"Scope": "WHOLE_MODEL", "Name": variant_name, "Start": 1,
                 "End": whole.get("n_residues", ""), "Length": whole.get("n_residues", ""),
                 "Mean_pLDDT": f(whole.get("mean_plddt")),
                 "Pct_pLDDT_ge70": f(whole.get("pct_plddt_ge_threshold")),
                 "Pct_pLDDT_lt50": f(whole.get("pct_plddt_lt_50")),
                 "Max_pLDDT": f(whole.get("max_plddt")),
                 "Intra_PAE": f(whole.get("mean_pae")),
                 "Longest_Confident_Run": (f"{run[0]}-{run[1]}" if run[0] else ""),
                 "Longest_Confident_Run_Len": whole.get("longest_confident_run_len", ""),
                 "Notes": f"max PAE {f(whole.get('max_pae'))} A"})

    for m in conf.get("models", []):
        rows.append({"Scope": "MODEL", "Name": f"model_{m['model']}",
                     "pTM": m.get("ptm"), "iPTM": m.get("iptm"),
                     "Fraction_Disordered": m.get("fraction_disordered"),
                     "Ranking_Score": m.get("ranking_score"),
                     "Has_Clash": m.get("has_clash")})

    for r in conf.get("regions", []):
        rows.append({"Scope": "REGION", "Name": r["region"], "Start": r["start"],
                     "End": r["end"], "Length": r["length"],
                     "Mean_pLDDT": f(r["mean_plddt"]),
                     "Pct_pLDDT_ge70": f(r["pct_plddt_ge_threshold"]),
                     "Intra_PAE": f(r["intra_pae"]),
                     "CrossModel_RMSD_Mean": f(r["rmsd_mean"]),
                     "CrossModel_RMSD_Min": f(r["rmsd_min"]),
                     "CrossModel_RMSD_Max": f(r["rmsd_max"]),
                     "CrossModel_RMSD_Pairs": r["rmsd_pairs"]})

    for p in conf.get("region_pairs", []):
        rows.append({"Scope": "REGION_PAIR", "Name": f"{p['region_a']} <-> {p['region_b']}",
                     "Inter_PAE": f(p["inter_pae"])})

    rows.append({"Scope": "VERDICT", "Name": "Docking_Readiness",
                 "Docking_Readiness": readiness, "Docking_Ready_Regions": regions_str,
                 "Notes": notes})

    fieldnames = ["Scope", "Name", "Start", "End", "Length", "Mean_pLDDT",
                  "Pct_pLDDT_ge70", "Pct_pLDDT_lt50", "Max_pLDDT", "pTM", "iPTM",
                  "Fraction_Disordered", "Ranking_Score", "Has_Clash", "Intra_PAE",
                  "Inter_PAE", "CrossModel_RMSD_Mean", "CrossModel_RMSD_Min",
                  "CrossModel_RMSD_Max", "CrossModel_RMSD_Pairs",
                  "Longest_Confident_Run", "Longest_Confident_Run_Len",
                  "Docking_Readiness", "Docking_Ready_Regions", "Notes"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})
    return path


def _pick_top_model(cif_files, source):
    """model_0 is AlphaFold Server's top-ranked model; fall back to sort order."""
    if not cif_files:
        raise FileNotFoundError(f"No .cif files found in {source}")
    model_0 = [f for f in cif_files if "_model_0" in os.path.basename(f)]
    return model_0[0] if model_0 else sorted(cif_files)[0]


def import_alphafold_output(downloaded_path, target_cif_path, archive_dir=None):
    """
    Takes a .zip downloaded from AlphaFold Server, the UNZIPPED FOLDER, or a
    bare .cif, and copies the top-ranked model into target_cif_path.

    If archive_dir is given, the ENTIRE download is preserved there first --
    all five ranked models, all five summary_confidences_*.json and all five
    full_data_*.json (which carry the PAE matrices). An earlier version kept
    only model_0 and discarded the rest, so every confidence signal AlphaFold
    produced was destroyed at import and the archive step had to be redone by
    hand later. Keep the whole download; it is small next to what it proves.

    Returns (target_cif_path, archived_dir_or_None).

    AlphaFold Server's output naming has changed before and may change again,
    so verify against your actual download the first time you use this.
    """
    if not os.path.exists(downloaded_path):
        raise FileNotFoundError(f"Downloaded file not found: {downloaded_path}")

    os.makedirs(os.path.dirname(target_cif_path), exist_ok=True)

    def _archive(src_dir, job_name):
        if archive_dir is None:
            return None
        dest = os.path.join(archive_dir, job_name)
        if os.path.abspath(src_dir) == os.path.abspath(dest):
            return dest            # already sitting in the archive; nothing to do
        os.makedirs(archive_dir, exist_ok=True)
        if os.path.isdir(dest):
            shutil.rmtree(dest)
        shutil.copytree(src_dir, dest)
        return dest

    # --- unzipped folder (what AlphaFold Server actually hands you on a Mac)
    if os.path.isdir(downloaded_path):
        cif_files = glob.glob(os.path.join(downloaded_path, "**", "*.cif"), recursive=True)
        chosen = _pick_top_model(cif_files, downloaded_path)
        archived = _archive(downloaded_path, os.path.basename(downloaded_path.rstrip(os.sep)))
        shutil.copy(chosen, target_cif_path)
        return target_cif_path, archived

    if downloaded_path.endswith(".cif"):
        # A bare .cif carries no confidence files -- nothing to archive, and
        # analyze_af3_confidence() will correctly report NOT_ASSESSED.
        shutil.copy(downloaded_path, target_cif_path)
        return target_cif_path, None

    if downloaded_path.endswith(".zip"):
        with tempfile.TemporaryDirectory() as tmp:
            with zipfile.ZipFile(downloaded_path, 'r') as z:
                z.extractall(tmp)
            cif_files = glob.glob(os.path.join(tmp, "**", "*.cif"), recursive=True)
            chosen = _pick_top_model(cif_files, downloaded_path)
            # A server zip may wrap everything in one folder, or not.
            entries = [os.path.join(tmp, e) for e in os.listdir(tmp)]
            src_dir = entries[0] if len(entries) == 1 and os.path.isdir(entries[0]) else tmp
            job_name = os.path.splitext(os.path.basename(downloaded_path))[0]
            archived = _archive(src_dir, job_name)
            shutil.copy(chosen, target_cif_path)
            return target_cif_path, archived

    raise ValueError(f"Unrecognized file type for AlphaFold Server output: {downloaded_path}")


def _resolve_paths():
    project_root = _PROJECT_ROOT
    input_csv = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA", "Filtered")
    variant_fasta_path = common.phase1g_fasta_path(project_root)
    output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepB")
    out_sec_dir = os.path.join(output_base, "Secondary_Structure")
    out_ter_dir = os.path.join(output_base, "Tertiary_Structure")
    return project_root, input_csv, variant_fasta_path, output_base, out_sec_dir, out_ter_dir


def _get_target(input_csv, variant_fasta_path):
    winner_row, _ = common.get_winner_from_filtered_csv(input_csv)
    if winner_row is None:
        return None, None, None
    winner_name = winner_row["Variant"]

    variants = common.load_multi_fasta(variant_fasta_path)
    if variants is None:
        print(f"[ERROR] Variant FASTA file not found: {variant_fasta_path}")
        return None, None, None

    vax_sequence = common.lookup_sequence(variants, winner_name)
    if vax_sequence is None:
        print(f"[ERROR] Could not find a matching header for '{winner_name}' in: {variant_fasta_path}")
        available = list(variants.keys())
        print(f"[ERROR] Headers actually present ({len(available)}): {available[:10]}{' ...' if len(available) > 10 else ''}")
        return None, None, None

    safe_name = common.sanitize_variant_name(winner_name)
    return winner_name, vax_sequence, safe_name


def run_prepare():
    project_root, input_csv, variant_fasta_path, output_base, out_sec_dir, out_ter_dir = _resolve_paths()
    os.makedirs(out_sec_dir, exist_ok=True)
    os.makedirs(out_ter_dir, exist_ok=True)

    common.print_banner("PHASE 2 STEP B -- PREPARE STRUCTURE PREDICTION SUBMISSIONS")
    winner_name, vax_sequence, safe_name = _get_target(input_csv, variant_fasta_path)
    if winner_name is None:
        return

    print(f"[INFO] Target Variant : {winner_name}")
    print("-" * 100)

    # Sec. II.B: "AlphaFold3 would be incorporated to construct a Tertiary
    # (3D) Structure Prediction ... of how the amino acids of THE VACCINE
    # fold. Human Toll-Like Receptor (TLR) models, specifically TLR-2
    # (6NIG) and TLR-4 (8WTA), were TAKEN FROM the RCSB Protein Data Bank."
    #
    # So AlphaFold predicts the vaccine MONOMER only. The TLR structures are
    # downloaded experimental coordinates, not predictions -- and the
    # vaccine-TLR COMPLEXES are built later by HADDOCK docking (Sec. III.A:
    # "utilized Haddock v2.4 by taking the AlphaFold3 model of the vaccine
    # and TLR-2 and TLR-4 receptors and docking it"). An earlier version of
    # this script asked AlphaFold to co-fold vaccine+TLR complexes, which is
    # both off-spec and wasteful: it burned two extra AlphaFold Server jobs
    # and replaced experimental structures with predicted ones.
    print("[INFO] Downloading TLR2/TLR4 experimental structures from RCSB (6NIG, 8WTA)...")
    tlr_targets = {
        "TLR2": ("6NIG", os.path.join(out_ter_dir, "RCSB_TLR2_6NIG.cif")),
        "TLR4": ("8WTA", os.path.join(out_ter_dir, "RCSB_TLR4_8WTA.cif")),
    }
    try:
        # requests + forced IPv4, NOT urllib.request. This environment's
        # IPv6 (NAT64) route to several hosts never completes its handshake:
        # urllib.request.urlretrieve hung indefinitely on files.rcsb.org
        # (>4 min, no bytes) where curl fetched the same file in <2s. The
        # same stall was diagnosed on IEDB in Step 2A and fixed the same
        # way. urllib.request does not go through urllib3, so the Step 2A
        # patch does not cover it -- this needs its own.
        import socket
        import requests
        import urllib3.util.connection as urllib3_cn
        urllib3_cn.allowed_gai_family = lambda: socket.AF_INET
        for label, (pdb_id, dest) in tlr_targets.items():
            if os.path.isfile(dest) and os.path.getsize(dest) > 0:
                print(f"[INFO] Reusing cached {label} structure: {os.path.basename(dest)}")
                continue
            resp = requests.get(f"https://files.rcsb.org/download/{pdb_id}.cif", timeout=120)
            resp.raise_for_status()
            with open(dest, "wb") as fh:
                fh.write(resp.content)
            print(f"[SUCCESS] {label} ({pdb_id}) -> {os.path.basename(dest)} "
                  f"({os.path.getsize(dest):,} bytes)")
    except Exception as e:
        print(f"[ERROR] Failed to download TLR structures from RCSB: {e}")
        return
    print("[NOTE] 6NIG is a TLR2/VLR-B crystallization chimera -- the VLR-B fusion")
    print("[NOTE] partner must be stripped before docking in Phase III. Flagged here")
    print("[NOTE] so it is not carried into HADDOCK unnoticed.")
    print("-" * 100)

    print("Submit this 1 job at https://alphafoldserver.com. After it finishes,")
    print("download the result and run the --import command shown for it.")
    print("-" * 100)

    jobs = [
        ("monomer", f"Monomer_{safe_name}", [("Vaccine", vax_sequence)]),
    ]

    for job_key, job_name, entities in jobs:
        print(f"\n[JOB: {job_key}]  AlphaFold Server job name suggestion: {job_name}")
        for entity_label, seq in entities:
            print(f"   Entity ({entity_label}):")
            print(f"   {seq}")
        print(f"   --> After downloading the result, run:")
        print(f"       python {os.path.basename(__file__)} --import {job_key} /path/to/downloaded_file.zip")

    print("\n" + "-" * 100)
    print("Submit the monomer sequence below to the two secondary-structure")
    print("predictors (separate from AlphaFold's tertiary model).")
    print("-" * 100)
    psipred_fasta = os.path.join(out_sec_dir, f"PsiPred_{safe_name}.fasta")
    netsurfp_fasta = os.path.join(out_sec_dir, f"NetSurfP_{safe_name}.fasta")
    with open(psipred_fasta, 'w') as f:
        f.write(f">PsiPred_{winner_name}\n{vax_sequence}\n")
    with open(netsurfp_fasta, 'w') as f:
        f.write(f">NetSurfP_{winner_name}\n{vax_sequence}\n")
    print(f"[PsiPred]   http://bioinf.cs.ucl.ac.uk/psipred  -- paste sequence; tick PSIPRED 4.0 only")
    print(f"            FASTA also written to: {psipred_fasta}")
    print(f"            --> python {os.path.basename(__file__)} --import psipred /path/to/downloaded_result.ss2")
    print(f"[NetSurfP]  https://services.healthtech.dtu.dk/services/NetSurfP-3.0/")
    print(f"            PASTE the raw sequence (the form takes a sequence, not a file);")
    print(f"            FASTA also written to: {netsurfp_fasta}")
    print(f"            Download results as CSV.")
    print(f"            SUBSTITUTE FOR RAPTORX (deviation #16): raptorx.uchicago.edu returns")
    print(f"            HTTP 500 on all endpoints. NetSurfP-3.0 is ESM-1b based and, unlike")
    print(f"            S4Pred, independent of the PSIPRED lineage.")
    print(f"            --> python {os.path.basename(__file__)} --import netsurfp /path/to/downloaded_result.csv")

    print("\n" + "-" * 100)
    print("Once all 3 are imported, run this script with no arguments to finish Step 2B.")
    print("=" * 100 + "\n")


def run_import(job_key, downloaded_path):
    project_root, input_csv, variant_fasta_path, output_base, out_sec_dir, out_ter_dir = _resolve_paths()
    winner_name, vax_sequence, safe_name = _get_target(input_csv, variant_fasta_path)
    if winner_name is None:
        return

    target_map = {
        "monomer": os.path.join(out_ter_dir, f"AF3_Target_{safe_name}.cif"),
        "psipred": os.path.join(out_sec_dir, f"PsiPred_{safe_name}.ss2"),
        "netsurfp": os.path.join(out_sec_dir, f"NetSurfP_{safe_name}.csv"),
    }

    if job_key not in target_map:
        print(f"[ERROR] Unknown job key '{job_key}'. Expected one of: {list(target_map.keys())}")
        return

    target_path = target_map[job_key]
    if job_key == "monomer":
        archive_dir = os.path.join(output_base, "AlphaFold_Raw")
        try:
            result, archived = import_alphafold_output(downloaded_path, target_path, archive_dir)
        except (FileNotFoundError, ValueError) as e:
            print(f"[ERROR] {e}")
            return
        if archived:
            print(f"[SUCCESS] Full AlphaFold download archived -> "
                  f"{os.path.relpath(archived, project_root)}")
            run_confidence_analysis(archived)
        else:
            print("[WARN] Only a bare .cif was imported -- no ranked models, summary "
                  "confidences or PAE data, so fold confidence cannot be assessed. "
                  "Re-import the full AlphaFold Server download to enable it.")
    else:
        # PsiPred/NetSurfP results are plain text files downloaded as-is
        # (not zips of ranked models like AlphaFold Server) -- just copy
        # them into place; parsing happens when Step 2B runs normally.
        if not os.path.isfile(downloaded_path):
            print(f"[ERROR] Downloaded file not found: {downloaded_path}")
            return
        os.makedirs(os.path.dirname(target_path), exist_ok=True)
        shutil.copy(downloaded_path, target_path)
        result = target_path

    print(f"[SUCCESS] Imported '{job_key}' model -> {os.path.relpath(result, project_root)}")


def _latest_af3_raw_dir(output_base):
    """Newest AlphaFold_Raw/<job>/ that actually contains ranked models."""
    root = os.path.join(output_base, "AlphaFold_Raw")
    if not os.path.isdir(root):
        return None
    candidates = [os.path.join(root, d) for d in sorted(os.listdir(root))
                  if os.path.isdir(os.path.join(root, d)) and not d.startswith("_")]
    candidates = [c for c in candidates if _af3_model_files(c)]
    if not candidates:
        return None
    return max(candidates, key=os.path.getmtime)


def run_confidence_analysis(raw_dir=None):
    """
    Reports AlphaFold fold confidence for the current construct and states
    what Phase III may dock.

    This is DIAGNOSTIC AND NON-BLOCKING by design: it never changes the model
    files and never fails the step. If the confidence files are missing it
    says so and returns -- the same tolerant pattern Phase 1De uses for
    SEMA's UNSCREENED.
    """
    project_root, input_csv, variant_fasta_path, output_base, _, _ = _resolve_paths()
    winner_name, vax_sequence, safe_name = _get_target(input_csv, variant_fasta_path)
    if winner_name is None:
        return None

    if raw_dir is None:
        raw_dir = _latest_af3_raw_dir(output_base)
    if raw_dir is None:
        print("[WARN] No archived AlphaFold download found under StepB/AlphaFold_Raw/ -- "
              "fold confidence NOT_ASSESSED. Re-run: --import monomer <download folder>")
        return None

    boundary_map = common.phase1g_boundary_map(project_root, winner_name)
    if boundary_map:
        regions = parse_boundary_map(boundary_map)
    else:
        print("[WARN] Phase 1G Boundary_Map unavailable -- reporting whole-model "
              "confidence only, without per-region breakdown.")
        regions = []

    conf = analyze_af3_confidence(raw_dir, regions)
    if conf is None:
        print(f"[WARN] {os.path.relpath(raw_dir, project_root)} holds no AlphaFold ranked "
              f"models -- fold confidence NOT_ASSESSED.")
        return None

    readiness, regions_str, notes = classify_docking_readiness(conf)

    common.print_banner("PHASE 2 STEP B: ALPHAFOLD3 FOLD CONFIDENCE")
    print(f"[INFO] Source download : {os.path.relpath(raw_dir, project_root)}")
    print(f"[INFO] Ranked models   : {conf['n_models']}")
    print("-" * 110)

    whole = conf.get("whole") or {}
    if whole:
        run = whole.get("longest_confident_run") or (None, None)
        print(f"Mean pLDDT (model_0)           | {whole['mean_plddt']:.1f}")
        print(f"Residues at pLDDT >= 70        | {whole['n_plddt_ge_threshold']}/"
              f"{whole['n_residues']} ({whole['pct_plddt_ge_threshold']:.1f}%)")
        print(f"Residues below pLDDT 50        | {whole['pct_plddt_lt_50']:.1f}%")
        print(f"Longest confident run          | "
              f"{f'{run[0]}-{run[1]}' if run[0] else 'none'} "
              f"({whole['longest_confident_run_len']} aa)")
    ptms = [m["ptm"] for m in conf["models"] if m.get("ptm") is not None]
    if ptms:
        print(f"pTM across ranked models       | "
              f"{', '.join(f'{p:.2f}' for p in ptms)}")
    fds = [m["fraction_disordered"] for m in conf["models"] if m.get("fraction_disordered") is not None]
    if fds:
        print(f"fraction_disordered            | {', '.join(f'{v:.2f}' for v in fds)}")

    if conf["regions"]:
        print("-" * 110)
        print(f"{'REGION':<24}{'RANGE':<12}{'LEN':>5}{'meanpLDDT':>11}"
              f"{'%>=70':>8}{'intraPAE':>10}{'RMSD mean':>11}{'RMSD range':>16}")
        for r in conf["regions"]:
            def _n(value, spec):
                return format(value, spec) if value is not None else "n/a"
            rng = ("n/a" if r["rmsd_min"] is None
                   else f"{r['rmsd_min']:.2f}-{r['rmsd_max']:.2f}")
            print(f"{r['region']:<24}"
                  f"{str(r['start']) + '-' + str(r['end']):<12}"
                  f"{r['length']:>5}"
                  f"{_n(r['mean_plddt'], '.1f'):>11}"
                  f"{_n(r['pct_plddt_ge_threshold'], '.0f'):>8}"
                  f"{_n(r['intra_pae'], '.2f'):>10}"
                  f"{_n(r['rmsd_mean'], '.2f'):>11}"
                  f"{rng:>16}")

    if conf["region_pairs"]:
        print("-" * 110)
        print("Inter-region PAE (confidence in RELATIVE placement of two regions):")
        for p in sorted(conf["region_pairs"], key=lambda x: -(x["inter_pae"] or 0)):
            if p["inter_pae"] is None:
                continue
            print(f"  {p['region_a']:<24} <-> {p['region_b']:<24} {p['inter_pae']:>8.2f} A")

    print("-" * 110)
    print(f"DOCKING READINESS : [{readiness}]")
    if regions_str:
        print(f"  Docking_Ready_Regions : {regions_str}")
    print(f"  Basis                 : {notes}")
    print("  NOTE: pLDDT>=70 / pTM>=0.5 / RMSD<=2.0 A are external convention, not")
    print("        manuscript criteria -- see deviation #21.")
    print("-" * 110)

    path = write_confidence_report(conf, readiness, regions_str, notes,
                                   output_base, project_root, winner_name)
    print(f"[INFO] Confidence report saved : {os.path.relpath(path, project_root)}")
    print("=" * 110 + "\n")
    return conf


def run_step2b_structure_prediction():
    start_time = time.time()
    project_root, input_csv, variant_fasta_path, output_base, out_sec_dir, out_ter_dir = _resolve_paths()
    os.makedirs(out_sec_dir, exist_ok=True)
    os.makedirs(out_ter_dir, exist_ok=True)

    common.print_banner("PHASE 2 STEP B: SECONDARY AND TERTIARY STRUCTURE PREDICTION")
    print(f"[INFO] Resolved Project Root : {project_root}")
    print(f"[INFO] Looking for Filtered CSVs in : {input_csv}")

    winner_name, vax_sequence, safe_name = _get_target(input_csv, variant_fasta_path)
    if winner_name is None:
        return

    print(f"[INFO] Target Variant    : {winner_name} (Rank 1 from Step 2A)")
    print("[INFO] Methodology       : PsiPred + NetSurfP-3.0 (secondary; NetSurfP substitutes RaptorX, dev #16); AlphaFold Server (vaccine monomer, Sec. II.B);")
    print("[INFO]                     TLR-2/TLR-4 experimental structures downloaded from RCSB (6NIG, 8WTA)")
    print("-" * 110)

    # Sec. II.B requires the AlphaFold model of the VACCINE only; the TLR
    # structures are downloaded experimental coordinates, and the complexes
    # are built by HADDOCK in Phase III (Sec. III.A) -- not predicted here.
    # So 3 manual imports are required, not 5.
    required_files = {
        "monomer": os.path.join(out_ter_dir, f"AF3_Target_{safe_name}.cif"),
        "psipred": os.path.join(out_sec_dir, f"PsiPred_{safe_name}.ss2"),
        "netsurfp": os.path.join(out_sec_dir, f"NetSurfP_{safe_name}.csv"),
    }
    # Downloaded automatically by --prepare; verified here so a missing or
    # truncated download is caught now rather than at docking time.
    rcsb_files = {
        "TLR2 (6NIG)": os.path.join(out_ter_dir, "RCSB_TLR2_6NIG.cif"),
        "TLR4 (8WTA)": os.path.join(out_ter_dir, "RCSB_TLR4_8WTA.cif"),
    }
    missing = [k for k, p in required_files.items() if not os.path.isfile(p)]
    if missing:
        print(f"[ERROR] Missing results for: {missing}")
        print(f"[ERROR] Run: python {os.path.basename(__file__)} --prepare")
        print("[ERROR] then submit the jobs and import each result.")
        return
    missing_rcsb = [k for k, p in rcsb_files.items() if not os.path.isfile(p)]
    if missing_rcsb:
        print(f"[ERROR] Missing RCSB TLR structures: {missing_rcsb}")
        print(f"[ERROR] Re-run: python {os.path.basename(__file__)} --prepare")
        return
    for label, p in rcsb_files.items():
        print(f"[INFO] {label} structure: {os.path.relpath(p, project_root)} ({os.path.getsize(p):,} bytes)")

    print("[INFO] All 3 required manual results found:")
    for k, p in required_files.items():
        print(f"       {k}: {os.path.relpath(p, project_root)}")
    print("-" * 110)

    try:
        psipred_result = _parse_psipred_ss2(required_files["psipred"])
        netsurfp_result = _parse_netsurfp_csv(required_files["netsurfp"])
    except (ValueError, OSError) as e:
        print(f"[ERROR] Failed to parse PsiPred/NetSurfP output: {e}")
        return

    # BASELINE (Biopython Chou-Fasman) -- a quick sequence-only estimate,
    # kept and clearly labeled as a baseline ONLY. It must never be
    # confused with, or substituted for, PsiPred/NetSurfP's real output.
    analysis = ProteinAnalysis(vax_sequence)
    b_helix, b_turn, b_sheet = analysis.secondary_structure_fraction()
    b_coil = 1.0 - (b_helix + b_sheet + b_turn)

    # Third column: DSSP assigned from the AlphaFold model's coordinates.
    # Supplementary only -- an assignment, not an independent prediction
    # (see _dssp_from_structure). Returns None if unavailable; never fatal.
    dssp_result = _dssp_from_structure(required_files["monomer"])

    def _fmt(res, key):
        return f"{res[key]*100:>10.2f}%" if res else f"{'n/a':>11}"

    print(f"{'STRUCTURAL ELEMENT':<26} | {'BASELINE (Chou-Fasman)':<22} | {'PsiPred':<11} | {'NetSurfP-3.0':<12} | {'DSSP(AF3)'}")
    print("-" * 110)
    print(f"{'Alpha-Helix (H)':<26} | {b_helix*100:>20.2f}% | {psipred_result['helix_fraction']*100:>9.2f}% | {netsurfp_result['helix_fraction']*100:>10.2f}% | {_fmt(dssp_result,'helix_fraction')}")
    print(f"{'Beta-Sheet (E)':<26} | {b_sheet*100:>20.2f}% | {psipred_result['sheet_fraction']*100:>9.2f}% | {netsurfp_result['sheet_fraction']*100:>10.2f}% | {_fmt(dssp_result,'sheet_fraction')}")
    print(f"{'Random Coils (C)':<26} | {b_coil*100:>20.2f}% | {psipred_result['coil_fraction']*100:>9.2f}% | {netsurfp_result['coil_fraction']*100:>10.2f}% | {_fmt(dssp_result,'coil_fraction')}")
    print(f"{'Turns (T, baseline only)':<26} | {b_turn*100:>20.2f}% | {'n/a':>10} | {'n/a':>11} | {'n/a':>11}")
    print("-" * 110)
    print("[NOTE] NetSurfP-3.0 substitutes RaptorX (deviation #16: RaptorX server down, HTTP 500).")
    print("[NOTE] DSSP(AF3) is assigned from the predicted 3D model, NOT an independent")
    print("[NOTE] sequence predictor -- it inherits AlphaFold's errors. Supplementary only.")
    print("[NOTE] Weight that caveat heavily here: --confidence measures this model's fold")
    print("[NOTE] as undetermined outside residues 1-45 (pTM 0.17), so the DSSP column is")
    print("[NOTE] reading secondary structure off coordinates AlphaFold is not confident in.")
    print("[NOTE] PsiPred and NetSurfP-3.0 are the load-bearing calls. See deviation #21.")

    report_path = os.path.join(out_sec_dir, "Step2B_Secondary_Structure_Comparison.csv")
    with open(report_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Element", "Baseline_ChouFasman", "PsiPred", "NetSurfP_3.0", "DSSP_from_AF3_model"])
        def _c(res, key):
            return f"{res[key]*100:.2f}%" if res else "n/a"
        writer.writerow(["Alpha-Helix", f"{b_helix*100:.2f}%", f"{psipred_result['helix_fraction']*100:.2f}%", f"{netsurfp_result['helix_fraction']*100:.2f}%", _c(dssp_result,'helix_fraction')])
        writer.writerow(["Beta-Sheet", f"{b_sheet*100:.2f}%", f"{psipred_result['sheet_fraction']*100:.2f}%", f"{netsurfp_result['sheet_fraction']*100:.2f}%", _c(dssp_result,'sheet_fraction')])
        writer.writerow(["Random Coils", f"{b_coil*100:.2f}%", f"{psipred_result['coil_fraction']*100:.2f}%", f"{netsurfp_result['coil_fraction']*100:.2f}%", _c(dssp_result,'coil_fraction')])
        writer.writerow(["Turns", f"{b_turn*100:.2f}%", "n/a", "n/a", "n/a"])

    total_time = common.format_time(time.time() - start_time)
    common.print_banner("STEP 2B COMPLETE")
    print("[SUCCESS] Real PsiPred/NetSurfP-3.0 secondary structure parsed; AlphaFold vaccine model in place.")
    for label, p in rcsb_files.items():
        print(f"[INFO] {label} (experimental, for Phase III docking): {os.path.relpath(p, project_root)}")
    print("[INFO] Vaccine-TLR complexes are NOT built here -- Sec. III.A builds them by HADDOCK docking.")
    print(f"[SUCCESS] Total Execution Time : {total_time}")
    print(f"[INFO] Outputs routed to       : {os.path.relpath(output_base, project_root)}")
    print("=" * 110 + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 2B: Secondary + Tertiary Structure Prediction")
    parser.add_argument("--prepare", action="store_true", help="Print sequences to submit to AlphaFold Server / PsiPred / RaptorX")
    parser.add_argument("--import", dest="import_args", nargs=2, metavar=("JOB_KEY", "DOWNLOADED_PATH"),
                         help="Import a downloaded result. JOB_KEY is monomer, psipred, or netsurfp.")
    parser.add_argument("--confidence", nargs="?", const=True, metavar="RAW_DIR",
                         help="Report AlphaFold fold confidence and Phase III docking readiness. "
                              "Defaults to the newest archived download under StepB/AlphaFold_Raw/.")
    args = parser.parse_args()

    if args.prepare:
        run_prepare()
    elif args.import_args:
        run_import(args.import_args[0], args.import_args[1])
    elif args.confidence:
        run_confidence_analysis(None if args.confidence is True else args.confidence)
    else:
        run_step2b_structure_prediction()
