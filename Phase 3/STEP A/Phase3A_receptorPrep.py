import os
import re
import sys
import csv
import glob
import time
import argparse
from datetime import datetime

# =============================================================================
# MINIMAL BOOTSTRAP -- locates the shared phase2_common module.
# Phase III reuses Phase II's common helpers rather than forking a third copy;
# every bug in Steps 2A-2D came from copy-pasted path logic drifting out of
# sync, and a "phase3_common" that duplicates 90% of phase2_common would
# recreate exactly that. See phase2_common.py's own header.
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
# PHASE 3 STEP A: RECEPTOR AND LIGAND PREPARATION FOR HADDOCK
#
# Sec. III.A: "utilized Haddock v2.4 by taking the AlphaFold3 model of the
# vaccine and TLR-2 and TLR-4 receptors and docking it." Sec. II.B already
# established that the TLR structures are TAKEN FROM RCSB (6NIG, 8WTA), not
# predicted -- Step 2B downloads them. This step turns those raw depositions
# and our own model into docking-ready inputs, and does NOT dock.
#
# Three things in the raw files would quietly corrupt a docking run:
#
# 1. 6NIG IS A CRYSTALLIZATION CHIMERA. Its "TLR2" entity is a TLR2
#    ectodomain fused to variable lymphocyte receptor B (VLR-B) from the
#    hagfish Eptatretus stoutii, used as a crystallization chaperone. From
#    the deposition's own _struct_ref_seq records: TLR2 (UniProt O60603) =
#    author residues 1-507; VLR-B (UniProt Q2YE02) = author residues 509-576.
#    Residue 508 is an engineered linker and IS modelled. Docking against the
#    unstripped chain would let the vaccine dock onto hagfish protein and
#    score it as a human TLR2 interaction.
#
# 2. 6NIG HAS AN AGONIST SITTING IN THE TARGET SITE. Ligand KQD is
#    Diprovocim, a synthetic small-molecule TLR2 agonist -- it occupies the
#    exact pocket we want the vaccine to dock into. Left in place it blocks
#    the site; HADDOCK would find a different, meaningless interface.
#
# 3. HADDOCK v2.4 DOES NOT ACCEPT mmCIF. RCSB and AlphaFold Server both
#    hand out .cif now. Everything this step emits is legacy .pdb.
#
# 8WTA (TLR4/MD-2) needs none of that surgery -- it is a clean two-protein
# complex -- but MD-2's bound lipid-A acyl chains are a real decision rather
# than obvious debris, so they are flagged explicitly (see LIGAND_POLICY).
# =============================================================================

# ---- Receptor definitions ---------------------------------------------------
# Author (auth_seq_id) numbering throughout -- that is what the coordinate
# records and any HADDOCK restraint file will use. Ranges verified against
# each entry's own _struct_ref_seq, not assumed from the UniProt entry.
RECEPTORS = {
    "TLR2": {
        "pdb_id": "6NIG",
        "source_cif": "RCSB_TLR2_6NIG.cif",
        "keep_chains": ["A"],
        "residue_range": (27, 507),
        "range_reason": (
            "6NIG entity 1 is a TLR2/VLR-B crystallization chimera: TLR2 (O60603) "
            "is auth 1-507, VLR-B (Q2YE02, Eptatretus stoutii) is auth 509-576, "
            "and 508 is engineered linker. Residues below 27 are the signal "
            "peptide and are not modelled anyway."),
        "drop_hetero": ["KQD", "HOH"],
        "drop_reason": {
            "KQD": "Diprovocim -- a synthetic TLR2 agonist occupying the docking site",
            "HOH": "crystallographic waters (HADDOCK adds its own solvation)",
        },
        "flag_hetero": ["NAG"],
        "flag_reason": {
            "NAG": ("N-linked glycans. KEPT by default: they are covalently attached "
                    "to the real receptor surface and occlude parts of it. Drop them "
                    "with --strip-glycans if the docking protocol assumes an "
                    "unglycosylated receptor -- but say which you did."),
        },
    },
    "TLR4": {
        "pdb_id": "8WTA",
        "source_cif": "RCSB_TLR4_8WTA.cif",
        # TLR4 and MD-2 must BOTH be present: TLR4 does not engage a ligand
        # without MD-2, so docking bare TLR4 models a receptor that does not
        # exist. 8WTA holds two copies of the complex (A/B = TLR4, C/D =
        # MD-2); one copy is enough and two would give HADDOCK a spurious
        # second interface.
        "keep_chains": ["A", "C"],
        "residue_range": None,
        "range_reason": (
            "8WTA entity 1 is clean, non-chimeric TLR4 (O00206) -- RCSB's own SIFTS "
            "alignment confirms no fusion partner, so no trimming is needed. Chain A "
            "= TLR4 ectodomain, chain C = MD-2 (Q9Y6Y9)."),
        "drop_hetero": ["HOH"],
        "drop_reason": {"HOH": "crystallographic waters"},
        "flag_hetero": ["0IL", "2IL", "GP4", "XIQ", "BMA", "NAG", "MAN"],
        "flag_reason": {
            "_default": ("MD-2 hydrophobic-pocket lipid and glycans. KEPT by default: "
                         "that pocket is FUNCTIONALLY OCCUPIED in the activated "
                         "receptor, and emptying it would expose a large artificial "
                         "cavity for the vaccine to dock into. This is a deliberate "
                         "modelling decision, not leftover debris -- state it."),
        },
    },
}

# Cross-check: if a deposition ever changes and one of these disappears,
# that is a fact worth reporting rather than silently succeeding.
EXPECTED = {
    "TLR2": {"n_residues": 478, "first": 27, "last": 507},
}


# =============================================================================
# mmCIF -> PDB, without a full structure library
#
# Bio.PDB can do this, but its PDBIO writer renumbers/renames atoms when the
# mmCIF uses names longer than the legacy format allows, and it silently drops
# chains whose IDs are multi-character. Reading _atom_site directly and
# writing the fixed-width records ourselves keeps auth numbering EXACTLY as
# deposited -- which matters because HADDOCK restraint files reference these
# residue numbers, and an off-by-one there is invisible until the results are
# nonsense.
# =============================================================================

def parse_atom_site(cif_path):
    """Yields _atom_site rows from an mmCIF as dicts keyed by column name."""
    with open(cif_path) as fh:
        lines = fh.readlines()

    i = 0
    n = len(lines)
    while i < n:
        if lines[i].strip() == "loop_":
            j = i + 1
            columns = []
            while j < n and lines[j].lstrip().startswith("_"):
                columns.append(lines[j].strip())
                j += 1
            if columns and columns[0].startswith("_atom_site."):
                names = [c.split(".", 1)[1] for c in columns]
                while j < n:
                    line = lines[j]
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#") or stripped.startswith("loop_") \
                            or stripped.startswith("_") or stripped.startswith("data_"):
                        break
                    values = _split_cif_row(stripped)
                    if len(values) == len(names):
                        yield dict(zip(names, values))
                    j += 1
                return
            i = j
        else:
            i += 1


def _split_cif_row(row):
    """Whitespace split that respects mmCIF's single/double-quoted values."""
    return [m.group(1) or m.group(2) or m.group(3)
            for m in re.finditer(r"'([^']*)'|\"([^\"]*)\"|(\S+)", row)]


def _pdb_atom_line(serial, row):
    """
    One fixed-width ATOM/HETATM record.

    Columns are positional in the legacy format, so this is built by slice
    position rather than by f-string concatenation -- a single character of
    drift shifts every coordinate field and produces a file that parses but
    is wrong.
    """
    name = row.get("auth_atom_id") or row["label_atom_id"]
    name = name.strip('"')
    # Atom names are left-justified from column 14 unless the element symbol
    # is one character, in which case they start at 15. Getting this wrong is
    # the classic cause of "HADDOCK thinks every CA is a calcium ion".
    element = (row.get("type_symbol") or "").strip()
    if len(name) < 4 and len(element) == 1:
        name = f" {name}"
    altloc = row.get("label_alt_id", ".")
    altloc = " " if altloc in (".", "?") else altloc[:1]
    resname = (row.get("auth_comp_id") or row["label_comp_id"])[:3]
    chain = (row.get("auth_asym_id") or row["label_asym_id"])[:1]
    resseq = int(row.get("auth_seq_id") or row["label_seq_id"])
    icode = row.get("pdbx_PDB_ins_code", "?")
    icode = " " if icode in (".", "?") else icode[:1]
    x, y, z = float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"])
    occ = float(row.get("occupancy", 1.0) or 1.0)
    bfac = float(row.get("B_iso_or_equiv", 0.0) or 0.0)
    record = "HETATM" if row.get("group_PDB") == "HETATM" else "ATOM  "
    return (f"{record}{serial:>5d} {name:<4s}{altloc}{resname:>3s} {chain}"
            f"{resseq:>4d}{icode}   {x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"{occ:>6.2f}{bfac:>6.2f}          {element:>2s}\n")


def write_pdb(rows, out_path, header_lines=()):
    """Writes ATOM/HETATM rows as a legacy PDB, TERing between chains."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    serial = 0
    written = 0
    with open(out_path, "w") as fh:
        for line in header_lines:
            fh.write(f"REMARK 220 {line}\n"[:81].rstrip() + "\n")
        prev_chain = None
        for row in rows:
            chain = (row.get("auth_asym_id") or row["label_asym_id"])[:1]
            if prev_chain is not None and chain != prev_chain:
                serial += 1
                fh.write(f"TER   {serial:>5d}\n")
            serial += 1
            fh.write(_pdb_atom_line(serial, row))
            written += 1
            prev_chain = chain
        if written:
            serial += 1
            fh.write(f"TER   {serial:>5d}\n")
        fh.write("END\n")
    return written


# =============================================================================
# Receptor preparation
# =============================================================================

def prepare_receptor(name, spec, source_dir, out_dir, strip_glycans=False):
    """
    Applies one receptor's chain / residue-range / heteroatom policy and
    writes the docking-ready .pdb. Returns a report dict; never raises on a
    policy miss, because a surprising deposition is something to REPORT, not
    something to crash on.
    """
    cif_path = os.path.join(source_dir, spec["source_cif"])
    report = {
        "Receptor": name, "PDB_ID": spec["pdb_id"],
        "Source_File": spec["source_cif"], "Status": "OK",
        "Chains_Kept": ",".join(spec["keep_chains"]),
        "Residue_Range": ("full" if not spec["residue_range"]
                          else f"{spec['residue_range'][0]}-{spec['residue_range'][1]}"),
        "Range_Rationale": spec["range_reason"],
        "N_Protein_Residues": 0, "N_Atoms_Written": 0,
        "First_Residue": "", "Last_Residue": "",
        "Hetero_Removed": "", "Hetero_Kept": "", "Hetero_Decisions": "",
        "Output_File": "", "Warnings": "",
    }
    if not os.path.isfile(cif_path):
        report["Status"] = "MISSING SOURCE"
        report["Warnings"] = (f"{spec['source_cif']} not found in "
                              f"{os.path.basename(source_dir)} -- run Step 2B --prepare "
                              f"to download it from RCSB.")
        return report, None

    keep_chains = set(spec["keep_chains"])
    drop_het = set(spec["drop_hetero"])
    flag_het = set(spec["flag_hetero"])
    if strip_glycans:
        drop_het |= {"NAG", "MAN", "BMA", "FUC", "GAL"}

    rows = []
    residues = {}
    removed, kept_het = {}, {}
    warnings = []

    for row in parse_atom_site(cif_path):
        chain = row.get("auth_asym_id") or row["label_asym_id"]
        if chain not in keep_chains:
            continue

        comp = (row.get("auth_comp_id") or row["label_comp_id"]).strip()
        is_het = row.get("group_PDB") == "HETATM"

        # Alternate conformations: keep only the first. HADDOCK expects one
        # coordinate per atom; leaving both gives it overlapping copies.
        alt = row.get("label_alt_id", ".")
        if alt not in (".", "?", "", "A"):
            continue

        if is_het:
            if comp in drop_het:
                removed[comp] = removed.get(comp, 0) + 1
                continue
            kept_het[comp] = kept_het.get(comp, 0) + 1
            if comp not in flag_het:
                warnings.append(f"unlisted heteroatom {comp} KEPT -- decide explicitly")
            rows.append(row)
            continue

        try:
            resseq = int(row.get("auth_seq_id") or row["label_seq_id"])
        except (TypeError, ValueError):
            continue
        if spec["residue_range"]:
            lo, hi = spec["residue_range"]
            if not (lo <= resseq <= hi):
                continue
        residues.setdefault(chain, set()).add(resseq)
        rows.append(row)

    all_res = sorted(r for s in residues.values() for r in s)
    report["N_Protein_Residues"] = len(all_res)
    report["First_Residue"] = all_res[0] if all_res else ""
    report["Last_Residue"] = all_res[-1] if all_res else ""
    report["Hetero_Removed"] = "; ".join(
        f"{k} x{v} atoms ({spec['drop_reason'].get(k, 'policy')})"
        for k, v in sorted(removed.items())) or "none"
    report["Hetero_Kept"] = "; ".join(f"{k} x{v} atoms" for k, v in sorted(kept_het.items())) or "none"
    report["Hetero_Decisions"] = " | ".join(
        f"{k}: {spec['flag_reason'].get(k) or spec['flag_reason'].get('_default', 'kept')}"
        for k in sorted(kept_het)) or "n/a"

    expected = EXPECTED.get(name)
    if expected:
        if len(all_res) != expected["n_residues"]:
            warnings.append(f"expected {expected['n_residues']} residues, got {len(all_res)} "
                            f"-- the deposition may have changed; re-verify _struct_ref_seq")
        if all_res and (all_res[0] != expected["first"] or all_res[-1] != expected["last"]):
            warnings.append(f"expected residues {expected['first']}-{expected['last']}, "
                            f"got {all_res[0]}-{all_res[-1]}")

    out_path = os.path.join(out_dir, f"Phase3A_{name}_{spec['pdb_id']}_receptor.pdb")
    header = [
        f"{name} receptor prepared for HADDOCK v2.4 from RCSB {spec['pdb_id']}.",
        f"Chains kept: {report['Chains_Kept']}. Residue range: {report['Residue_Range']}.",
        spec["range_reason"],
        f"Heteroatoms removed: {report['Hetero_Removed']}",
        f"Heteroatoms kept: {report['Hetero_Kept']}",
    ]
    report["N_Atoms_Written"] = write_pdb(rows, out_path, header)
    report["Output_File"] = os.path.relpath(out_path, _PROJECT_ROOT)
    report["Warnings"] = " | ".join(warnings) if warnings else "none"
    if warnings:
        report["Status"] = "OK (with warnings)"
    return report, out_path


# =============================================================================
# Ligand (vaccine) preparation
# =============================================================================

def _parse_docking_regions(regions_str):
    """"1-45" or "1-45;120-160" -> [(1,45), (120,160)]."""
    out = []
    for chunk in (regions_str or "").split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        m = re.match(r"^(\d+)\s*-\s*(\d+)$", chunk)
        if m:
            out.append((int(m.group(1)), int(m.group(2))))
    return out


def read_step2d_docking_scope(project_root):
    """
    Reads Docking_Readiness / Docking_Ready_Regions from the newest Step 2D
    report.

    NEVER HARDCODE THE REGION HERE. It is measured per model in Step 2B (see
    deviation #21) and would go stale the moment the construct or the fold
    changes -- which is exactly the failure this whole check exists to
    prevent.
    """
    folder = os.path.join(project_root, "Step_Outputs", "Phase2", "StepD")
    reports = sorted(glob.glob(os.path.join(folder, "Step2D_Validation_Report_*.csv")))
    if not reports:
        return None
    try:
        with open(reports[-1], newline="") as fh:
            rows = list(csv.DictReader(fh))
    except OSError:
        return None
    if not rows:
        return None
    r = rows[0]
    return {
        "report": os.path.basename(reports[-1]),
        "readiness": r.get("Docking_Readiness", "UNKNOWN"),
        "regions": r.get("Docking_Ready_Regions", ""),
        "basis": r.get("Docking_Readiness_Basis", ""),
        "overall_status": r.get("Overall_Status", ""),
        "variant": r.get("Variant", ""),
    }


def prepare_ligand(project_root, scope, out_dir):
    """
    Slices the Phase III docking scope out of the Phenix-refined model.

    The refined pass1 PDB is the same file Step 2D validated -- using the raw
    AlphaFold .cif instead would dock coordinates that no MolProbity number
    in this project describes.
    """
    src_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepC", "Supplementary_Archive")
    candidates = sorted(glob.glob(os.path.join(src_dir, "*_phenix_pass1.pdb")))
    report = {"Receptor": "VACCINE (ligand)", "PDB_ID": "-", "Status": "OK",
              "Source_File": "", "Chains_Kept": "A", "Residue_Range": "",
              "Range_Rationale": "", "N_Protein_Residues": 0, "N_Atoms_Written": 0,
              "First_Residue": "", "Last_Residue": "", "Hetero_Removed": "n/a",
              "Hetero_Kept": "n/a", "Hetero_Decisions": "n/a", "Output_File": "",
              "Warnings": ""}
    if not candidates:
        report["Status"] = "MISSING SOURCE"
        report["Warnings"] = f"No *_phenix_pass1.pdb in {os.path.relpath(src_dir, project_root)} -- run Step 2C."
        return report, None

    src = candidates[-1]
    report["Source_File"] = os.path.basename(src)

    ranges = _parse_docking_regions(scope["regions"]) if scope else []
    if scope and scope["readiness"] == "WHOLE_CONSTRUCT":
        ranges = []          # whole model is dockable; no slicing
        report["Range_Rationale"] = "Step 2D: WHOLE_CONSTRUCT -- global fold determined."
    elif ranges:
        report["Range_Rationale"] = (
            f"Step 2D: {scope['readiness']} -- {scope['basis']}")
    else:
        report["Status"] = "BLOCKED"
        report["Warnings"] = (
            f"Step 2D reports Docking_Readiness="
            f"{scope['readiness'] if scope else 'NOT_ASSESSED'} with no dockable region. "
            f"Nothing in this model supports a docking claim -- do NOT proceed to "
            f"HADDOCK. Run Step 2B --confidence first if this is unexpected.")
        return report, None

    kept, residues = [], set()
    with open(src) as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            if ranges and not any(lo <= resseq <= hi for lo, hi in ranges):
                continue
            kept.append(line)
            residues.add(resseq)

    report["Residue_Range"] = scope["regions"] if ranges else "full"
    report["N_Protein_Residues"] = len(residues)
    report["N_Atoms_Written"] = len(kept)
    if residues:
        report["First_Residue"], report["Last_Residue"] = min(residues), max(residues)

    out_path = os.path.join(out_dir, "Phase3A_Vaccine_ligand.pdb")
    os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w") as fh:
        for line in [
            f"{scope['variant'] if scope else 'vaccine'} -- Phase III docking scope.",
            f"Sliced from {os.path.basename(src)} (the Phenix-refined model Step 2D validated).",
            f"Step 2D Docking_Readiness: {scope['readiness'] if scope else 'NOT_ASSESSED'}; "
            f"regions {report['Residue_Range']}.",
            "See deviation #21 for why the scope is not the whole construct.",
        ]:
            fh.write(f"REMARK 220 {line}\n"[:81].rstrip() + "\n")
        fh.writelines(kept)
        fh.write("TER\nEND\n")
    report["Output_File"] = os.path.relpath(out_path, project_root)
    report["Warnings"] = "none"
    return report, out_path


def verify_pdb(path):
    """Re-parses an emitted file with Bio.PDB -- catches column-drift bugs."""
    try:
        import warnings as _w
        from Bio.PDB import PDBParser
        from Bio import BiopythonWarning
        with _w.catch_warnings():
            _w.simplefilter("ignore", BiopythonWarning)
            s = PDBParser(QUIET=True).get_structure("x", path)
        chains = [c.id for c in s[0]]
        n = sum(1 for _ in s.get_residues())
        return True, f"parses OK -- chains {','.join(chains)}, {n} residues"
    except Exception as e:
        return False, f"FAILED to parse: {e}"


def run_step3a_receptor_prep(strip_glycans=False):
    start = time.time()
    project_root = _PROJECT_ROOT
    source_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepB", "Tertiary_Structure")
    out_dir = os.path.join(project_root, "Step_Outputs", "Phase3", "StepA")
    os.makedirs(out_dir, exist_ok=True)

    common.print_banner("PHASE 3 STEP A: RECEPTOR & LIGAND PREPARATION FOR HADDOCK")
    print(f"[INFO] Resolved Project Root : {project_root}")
    print(f"[INFO] Source structures     : {os.path.relpath(source_dir, project_root)}")
    print("[INFO] Output format         : legacy .pdb (HADDOCK v2.4 does not read mmCIF)")
    print("-" * 110)

    scope = read_step2d_docking_scope(project_root)
    if scope:
        print(f"[INFO] Step 2D report        : {scope['report']}")
        print(f"[INFO] Stereochemistry       : {scope['overall_status']}")
        print(f"[INFO] Docking readiness     : {scope['readiness']}  regions={scope['regions'] or 'none'}")
    else:
        print("[WARN] No Step 2D validation report found -- docking scope unknown. "
              "Run Step 2D first.")
    print("-" * 110)

    reports = []
    for name, spec in RECEPTORS.items():
        rep, path = prepare_receptor(name, spec, source_dir, out_dir, strip_glycans)
        if path:
            ok, msg = verify_pdb(path)
            rep["Warnings"] = rep["Warnings"] if ok else f"{rep['Warnings']} | {msg}"
            if not ok:
                rep["Status"] = "FAILED"
            print(f"[{'SUCCESS' if ok else 'ERROR'}] {name:<6} {spec['pdb_id']} -> "
                  f"{os.path.basename(path)}")
            print(f"          {rep['N_Protein_Residues']} residues "
                  f"({rep['First_Residue']}-{rep['Last_Residue']}), "
                  f"{rep['N_Atoms_Written']} atoms; {msg}")
            print(f"          removed: {rep['Hetero_Removed']}")
            print(f"          kept   : {rep['Hetero_Kept']}")
        else:
            print(f"[ERROR] {name}: {rep['Warnings']}")
        if rep["Warnings"] not in ("none", ""):
            print(f"          WARN   : {rep['Warnings']}")
        reports.append(rep)

    lig_rep, lig_path = prepare_ligand(project_root, scope, out_dir)
    if lig_path:
        ok, msg = verify_pdb(lig_path)
        if not ok:
            lig_rep["Status"] = "FAILED"
            lig_rep["Warnings"] = msg
        print(f"[{'SUCCESS' if ok else 'ERROR'}] VACCINE -> {os.path.basename(lig_path)}")
        print(f"          {lig_rep['N_Protein_Residues']} residues "
              f"({lig_rep['First_Residue']}-{lig_rep['Last_Residue']}), "
              f"{lig_rep['N_Atoms_Written']} atoms; {msg}")
        print(f"          scope  : {lig_rep['Range_Rationale']}")
    else:
        print(f"[BLOCKED] VACCINE: {lig_rep['Warnings']}")
    reports.append(lig_rep)

    print("-" * 110)
    blocked = [r for r in reports if r["Status"] in ("BLOCKED", "FAILED", "MISSING SOURCE")]
    if blocked:
        print(f"FINAL DECISION : [NOT READY] -- {len(blocked)} input(s) unresolved")
    else:
        print("FINAL DECISION : [INPUTS PREPARED]")
        print("  NOTE: this step PREPARES inputs only. It does not dock -- HADDOCK v2.4")
        print("        is run separately (Sec. III.A).")
    print("-" * 110)

    ts = datetime.now().strftime("%Y%m%d_%H%M")
    report_path = os.path.join(out_dir, f"Step3A_ReceptorPrep_Report_{ts}.csv")
    fields = ["Receptor", "PDB_ID", "Status", "Source_File", "Chains_Kept", "Residue_Range",
              "Range_Rationale", "N_Protein_Residues", "First_Residue", "Last_Residue",
              "N_Atoms_Written", "Hetero_Removed", "Hetero_Kept", "Hetero_Decisions",
              "Output_File", "Warnings"]
    with open(report_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in reports:
            w.writerow({k: r.get(k, "") for k in fields})

    common.print_banner("PHASE 3 STEP A COMPLETE")
    print(f"[SUCCESS] Execution Time : {common.format_time(time.time() - start)}")
    print(f"[INFO] Report Saved      : {os.path.relpath(report_path, project_root)}")
    print("=" * 110 + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase 3 Step A: Receptor & ligand prep for HADDOCK")
    parser.add_argument("--strip-glycans", action="store_true",
                        help="Also remove N-linked glycans (NAG/MAN/BMA/FUC/GAL). "
                             "Off by default -- glycans occlude real receptor surface.")
    args = parser.parse_args()
    run_step3a_receptor_prep(strip_glycans=args.strip_glycans)
