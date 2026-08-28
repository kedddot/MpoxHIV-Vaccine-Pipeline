import os
import sys
import csv
import json
import time
import glob
import shutil
import tempfile
import subprocess
import argparse
from datetime import datetime

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
# CONFIGURE -- fill these in for your machine before running
# =============================================================================
FOLDX_BINARY = os.environ.get("FOLDX_BINARY", os.path.join(_PROJECT_ROOT, "FoldX5", "foldx5.1_20261231"))
APBS_BINARY = os.environ.get("APBS_BINARY", "/opt/miniconda3/envs/apbs_env/bin/apbs")

# =============================================================================
# PER-RESIDUE pLDDT -- must be parsed from the ORIGINAL AlphaFold Server
# mmCIF, before any repair step. PDBFixer -> OpenMM -> Phenix each rewrite
# the B-factor column downstream, so pLDDT is unrecoverable once any of
# them has touched the file.
# =============================================================================

def _parse_mmcif_loop(lines, category_prefix):
    """
    Generic minimal mmCIF loop-block parser: given a category prefix like
    '_atom_site.' or '_ma_qa_metric_local.', returns a list of
    {column_name: value_str} dict rows for the first matching loop_ block
    -- reading column NAMES from the loop's own header instead of
    assuming a fixed column order/count, since mmCIF writers don't
    guarantee either (same "read the identifier from the response"
    principle applied elsewhere in this pipeline to IEDB output).
    """
    i, n = 0, len(lines)
    while i < n:
        if lines[i].strip() == "loop_":
            j = i + 1
            columns = []
            while j < n and lines[j].strip().startswith(category_prefix):
                columns.append(lines[j].strip()[len(category_prefix):])
                j += 1
            if columns:
                rows = []
                while j < n:
                    row_line = lines[j].strip()
                    if not row_line or row_line == "loop_" or row_line.startswith("_") or row_line.startswith("#"):
                        break
                    values = row_line.split()
                    if len(values) >= len(columns):
                        rows.append(dict(zip(columns, values)))
                    j += 1
                return rows
        i += 1
    return []


def parse_plddt_per_residue(cif_path):
    """
    Extracts real per-residue pLDDT confidence from an AlphaFold Server
    mmCIF. Tries the `_ma_qa_metric_local` category first (AlphaFold
    Server's dedicated per-residue confidence table); falls back to the
    `_atom_site` CA-atom B-factor column (AlphaFold's convention: pLDDT
    written into B-factor) if that category isn't present.

    Returns: {residue_seq_id (int): plddt (float)}.
    """
    with open(cif_path) as f:
        lines = f.readlines()

    plddt_by_residue = {}
    for row in _parse_mmcif_loop(lines, "_ma_qa_metric_local."):
        try:
            plddt_by_residue[int(row["label_seq_id"])] = float(row["metric_value"])
        except (KeyError, ValueError):
            continue
    if plddt_by_residue:
        return plddt_by_residue

    for row in _parse_mmcif_loop(lines, "_atom_site."):
        try:
            if row.get("label_atom_id") != "CA":
                continue
            plddt_by_residue[int(row["label_seq_id"])] = float(row["B_iso_or_equiv"])
        except (KeyError, ValueError):
            continue

    if not plddt_by_residue:
        raise ValueError(
            f"Could not parse per-residue pLDDT from {cif_path} via either "
            f"'_ma_qa_metric_local' or '_atom_site.B_iso_or_equiv' -- inspect "
            f"the file's actual mmCIF categories and adjust this parser."
        )
    return plddt_by_residue


def region_mean_plddt(plddt_table, start, end):
    """Mean pLDDT over an inclusive 1-based [start, end] residue range."""
    values = [v for seq_id, v in plddt_table.items() if start <= seq_id <= end]
    if not values:
        raise ValueError(f"No residues found in the pLDDT table for range {start}-{end}")
    return sum(values) / len(values)


# =============================================================================
# REAL TOOLS -- repair, protonation, SASA, FoldX, APBS
# =============================================================================

def repair_structure_pdbfixer(input_path, output_path):
    """
    Structure repair (methodology prep step 1, part A) via OpenMM PDBFixer.
    Adds missing atoms/residues, removes heterogens.

    NOTE: this step alone does NOT relax/minimize the structure -- it only
    patches in missing atoms. See minimize_structure_openmm() below for the
    actual "minimal relaxation" half of prep step 1; both must run in
    sequence, or FoldX (and to a lesser extent SASA) can report wildly
    unrealistic values off of un-relaxed steric clashes.

    Requires: pip install pdbfixer openmm   (or via conda-forge if pip
    gives you trouble with the OpenMM C-extension dependency)
    """
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    fixer = PDBFixer(filename=input_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.removeHeterogens(keepWater=False)

    with open(output_path, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    return output_path


def minimize_structure_openmm(input_pdb_path, output_pdb_path, max_iterations=500):
    """
    Real energy minimization (methodology prep step 1, part B) via OpenMM,
    using the AMBER14 forcefield with GBn2 implicit solvent. This resolves
    local steric clashes left over from PDBFixer's atom-repair alone --
    AlphaFold-predicted multi-domain constructs joined by short linkers are
    especially prone to minor clashes right at the junctions, which
    downstream tools (particularly FoldX) are very sensitive to: an
    un-minimized structure can report wildly unrealistic large-positive
    ddG values that reflect clash artifacts, not genuine instability.

    Requires: pip install openmm (already installed alongside pdbfixer)
    """
    from openmm.app import PDBFile, ForceField, Modeller, Simulation, HBonds, NoCutoff
    from openmm import LangevinMiddleIntegrator
    from openmm.unit import kelvin, picosecond, picoseconds

    pdb = PDBFile(input_pdb_path)
    forcefield = ForceField('amber14-all.xml', 'implicit/gbn2.xml')

    modeller = Modeller(pdb.topology, pdb.positions)
    # PDBFixer's repair step only adds missing HEAVY atoms -- it never adds
    # hydrogens, so the forcefield below (which needs a fully complete
    # molecule to match its templates) can't build a System without them
    # yet. This is a rough, generic hydrogen placement purely so
    # minimization has something physically complete to work with -- it is
    # NOT the rigorous, pKa-aware pH 7.0 protonation your methodology
    # calls for; that happens later, more accurately, via PDB2PQR/PROPKA.
    modeller.addHydrogens(forcefield, pH=7.0)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
    )
    integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    simulation.minimizeEnergy(maxIterations=max_iterations)
    positions = simulation.context.getState(getPositions=True).getPositions()

    # Strip hydrogens back out before writing the final file. They were
    # necessary for the forcefield's physics during minimization, but
    # every downstream consumer of this file -- our own FreeSASA call,
    # CamSol, Aggrescan3D, and PDB2PQR's own re-protonation -- expects a
    # heavy-atom-only structure, matching virtually any real PDB/
    # AlphaFold structure. Leaving hydrogens in caused external tools
    # (e.g. CamSol's own FreeSASA-based backend) to crash outright, and
    # likely silently skewed our own local SASA numbers too.
    minimized_modeller = Modeller(modeller.topology, positions)
    hydrogens = [a for a in minimized_modeller.topology.atoms() if a.element is not None and a.element.symbol == 'H']
    minimized_modeller.delete(hydrogens)

    with open(output_pdb_path, 'w') as f:
        PDBFile.writeFile(minimized_modeller.topology, minimized_modeller.positions, f)
    return output_pdb_path


def refine_structure_phenix(input_pdb_path, output_dir, safe_name, second_pass=False):
    """
    Real Phenix geometry refinement -- resolves bad angles, CaBLAM outliers,
    and improves the Ramachandran distribution beyond what OpenMM's generic
    physics-based minimization alone reliably achieves. OpenMM minimizes
    total energy; it has no concept of "Ramachandran favored region" or
    "CaBLAM outlier" (statistical categories from real solved structures),
    so it can land in a real physical energy minimum that's still
    statistically unusual by those standards. Phenix's tools restrain
    directly toward known-good conformational statistics.

    Requires: Phenix installed with its command-line tools on PATH
    (`phenix.geometry_minimization` runnable directly in your terminal).

    Pass 1 (always run): phenix.geometry_minimization -- fixes bond/angle
    strain and backbone geometry, no density map needed (real-space mode).
    This alone took a real run from 94.85% -> 99.26% Ramachandran favored,
    0.27% -> 0% bad angles, and 4.5% -> 1.5% CaBLAM outliers.

    Pass 2 (optional, second_pass=True): phenix.model_idealization -- a
    more targeted second pass that specifically re-fixes remaining rotamer
    outliers (additionally_fix_rotamer_outliers=True) alongside further
    Ramachandran restraints, for squeezing out residual outliers pass 1
    doesn't fully resolve (e.g. CaBLAM still above its <1.0% goal).
    """
    # Phenix parses EVERY command-line argument as a `name=value` PHIL
    # parameter definition, so any path containing "=" is misread as a
    # parameter assignment and the run dies with
    # 'improper definition name ...'. Our variant names carry the FASTA
    # header's "length=496", so the sanitized filename contains "=" and
    # trips this. Rather than change sanitize_variant_name() globally
    # (which would orphan the .cif/.ss2/.csv/manual-results files already
    # written under the "=" form), Phenix alone is run inside a scratch
    # directory on "="-free copies, and the outputs are copied back under
    # the canonical names.
    with tempfile.TemporaryDirectory() as scratch:
        scratch_input = os.path.join(scratch, "model_in.pdb")
        shutil.copy(input_pdb_path, scratch_input)

        # Parameter is `output_file_name_prefix`, NOT `output.file_name`
        # (confirmed against `phenix.geometry_minimization --show-defaults`;
        # `file_name` in that listing is the INPUT model, not the output).
        #
        # Ramachandran restraints are OFF by default in Phenix but are
        # enabled here. They guide backbone dihedrals toward the regions
        # actually observed in experimental structures -- standard practice
        # in refinement, and squarely within Sec. II.D's "structure
        # minimization and repair ... to adjust the atomic coordinates."
        # Measured on this construct (AlphaFold model, mean pLDDT 37):
        #     restraints OFF -> MolProbity 1.26, Rama favored 91.09%, FAILED
        #     restraints ON  -> MolProbity 0.71, Rama favored 98.99%, VALIDATED
        # RMS(bonds) 0.0032 and RMS(angles) 0.86 are IDENTICAL either way,
        # so this is not trading geometric strain for a better Rama score --
        # the backbone is guided, not distorted.
        scratch_p1 = os.path.join(scratch, "model_p1.pdb")
        cmd1 = [
            "phenix.geometry_minimization",
            scratch_input,
            "output_file_name_prefix=model_p1",
            "pdb_interpretation.ramachandran_plot_restraints.enabled=True",
        ]
        result = subprocess.run(cmd1, capture_output=True, text=True, cwd=scratch)
        if result.returncode != 0:
            raise RuntimeError(f"phenix.geometry_minimization failed (exit {result.returncode}):\n{result.stderr}")
        if not os.path.isfile(scratch_p1):
            # Phenix sometimes ignores output.file_name and derives its own
            # name from the input; fall back to whatever .pdb it produced.
            produced = [f for f in glob.glob(os.path.join(scratch, "*.pdb"))
                        if os.path.abspath(f) != os.path.abspath(scratch_input)]
            if not produced:
                raise FileNotFoundError(
                    f"phenix.geometry_minimization reported success but wrote no .pdb in {scratch}"
                )
            scratch_p1 = max(produced, key=os.path.getmtime)

        pass1_output = os.path.join(output_dir, f"{safe_name}_phenix_pass1.pdb")
        shutil.copy(scratch_p1, pass1_output)

        if not second_pass:
            return pass1_output

        cmd2 = [
            "phenix.model_idealization",
            scratch_p1,
            "output_prefix=model_p2",
        ]
        result2 = subprocess.run(cmd2, capture_output=True, text=True, cwd=scratch)
        if result2.returncode != 0:
            raise RuntimeError(f"phenix.model_idealization failed (exit {result2.returncode}):\n{result2.stderr}")
        produced2 = [f for f in glob.glob(os.path.join(scratch, "model_p2*.pdb"))]
        if not produced2:
            raise FileNotFoundError(
                f"phenix.model_idealization reported success but wrote no model_p2*.pdb in {scratch}"
            )
        pass2_output = os.path.join(output_dir, f"{safe_name}_phenix_pass2.pdb")
        shutil.copy(max(produced2, key=os.path.getmtime), pass2_output)
        return pass2_output


def protonate_pdb2pqr(input_pdb_path, output_pqr_path, ph=7.0):
    """
    Real protonation at target pH (methodology prep step 2) via PDB2PQR.
    Requires: pip install pdb2pqr
    """
    cmd = [
        "pdb2pqr30",
        "--ff=AMBER",
        f"--with-ph={ph}",
        "--titration-state-method=propka",
        input_pdb_path,
        output_pqr_path,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"PDB2PQR failed (exit {result.returncode}):\n{result.stderr}")
    return output_pqr_path


# Aligned to the paper's own hydrophobic-residue set (Sec. I.E / II.A.a:
# A, V, I, L, M, F, W, Y) -- this previously used a different, undocumented
# set (…PRO, CYS instead of TYR), which is internally inconsistent with the
# sequence-level hydrophobicity gate applied everywhere else in the
# pipeline. Measured effect of the correction on the current construct:
# hydrophobic surface fraction 0.4163 -> 0.4371 (see plan §0.3/II-1) --
# reported explicitly rather than left as a silent discrepancy.
HYDROPHOBIC_RESIDUES = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR'}
EXPOSURE_THRESHOLD_A2 = 15.0  # per-residue SASA above this counts as "exposed"


def analyze_sasa_freesasa(pdb_path):
    """
    Real SASA + hydrophobic surface fraction via FreeSASA.
    Requires: pip install freesasa

    NOTE on "exposed_patch_area": true spatial patch detection (like
    Aggrescan3D does) clusters residues by 3D proximity. This computes a
    simpler, sequence-adjacency-based approximation -- the longest run of
    consecutive, exposed, hydrophobic residues along the chain. That is a
    real, defensible metric, but it is NOT equivalent to a true 3D
    spatial-clustering patch. Flagged here so you don't mistake this for
    Aggrescan3D-grade analysis; Aggrescan3D itself is still handled
    manually (see run_step2c_manual_prepare / manual_results.json).
    """
    import freesasa

    structure = freesasa.Structure(pdb_path)
    result = freesasa.calc(structure)
    residue_areas = result.residueAreas()

    total_sasa = result.totalArea()
    hydrophobic_sasa = 0.0
    max_patch_area = 0.0
    current_patch = 0.0

    for chain_id in residue_areas:
        for res_num in sorted(residue_areas[chain_id], key=lambda x: int(x)):
            area = residue_areas[chain_id][res_num]
            is_hydrophobic = area.residueType in HYDROPHOBIC_RESIDUES
            is_exposed = area.total > EXPOSURE_THRESHOLD_A2

            if is_hydrophobic:
                hydrophobic_sasa += area.total

            if is_hydrophobic and is_exposed:
                current_patch += area.total
                max_patch_area = max(max_patch_area, current_patch)
            else:
                current_patch = 0.0

    hydrophobic_fraction = hydrophobic_sasa / total_sasa if total_sasa > 0 else 0.0

    return {
        "hydrophobic_fraction": hydrophobic_fraction,
        "exposed_patch_area": max_patch_area,
        "total_sasa": total_sasa,
    }


def run_foldx(pdb_path, output_dir):
    """
    Real FoldX Stability command -> Total Folding Energy.

    NOTE ON NAMING: this is FoldX's raw `--command=Stability` total folding
    energy, NOT a true ddG -- ddG requires a mutant vs. wild-type pair, and
    this pipeline has no designed mutants (a non-mutational vaccine
    construct). The paper's ddG > +1.5 kcal/mol gate is therefore N/A here
    and is not applied; see Methods_Deviations_RRL_Support.txt item #10.
    Kept as a CONTEXTUAL total-energy value only.

    Requires: FoldX binary, free academic registration at foldxsuite.crg.eu
    Set FOLDX_BINARY at the top of this file to your install path.

    NOTE: FoldX's exact CLI flags and output filename have changed across
    versions -- verify the "_0_ST.fxout" output name matches your version
    the first time you run this.
    """
    if not os.path.isfile(FOLDX_BINARY):
        raise FileNotFoundError(
            f"FoldX binary not found at '{FOLDX_BINARY}'. Set FOLDX_BINARY "
            f"at the top of this file to your actual FoldX install path."
        )

    pdb_dir = os.path.dirname(pdb_path)
    pdb_name = os.path.basename(pdb_path)

    # FoldX reads rotabase.txt from its OWN working directory at runtime,
    # but we run it with cwd set to the PDB's folder (below) -- so copy
    # rotabase.txt alongside the PDB each time rather than relying on you
    # to remember to place it there manually. It must live next to the
    # FoldX binary itself (that's the standard FoldX distribution layout).
    rotabase_src = os.path.join(os.path.dirname(FOLDX_BINARY), "rotabase.txt")
    rotabase_dst = os.path.join(pdb_dir, "rotabase.txt")
    if os.path.isfile(rotabase_src) and not os.path.isfile(rotabase_dst):
        shutil.copy(rotabase_src, rotabase_dst)
    elif not os.path.isfile(rotabase_src):
        raise FileNotFoundError(
            f"rotabase.txt not found next to FOLDX_BINARY at '{rotabase_src}'. "
            f"FoldX needs this file to run -- it should have come with your "
            f"FoldX download; place it in the same folder as the binary."
        )

    cmd = [
        FOLDX_BINARY,
        "--command=Stability",
        f"--pdb={pdb_name}",
        f"--output-dir={output_dir}",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=pdb_dir)
    if result.returncode != 0:
        raise RuntimeError(f"FoldX failed (exit {result.returncode}):\n{result.stderr}")

    stability_file = os.path.join(output_dir, f"{os.path.splitext(pdb_name)[0]}_0_ST.fxout")
    if not os.path.isfile(stability_file):
        raise FileNotFoundError(
            f"Expected FoldX output not found: {stability_file}\n"
            f"Your FoldX version may name its output differently -- check {output_dir} manually."
        )

    # Some FoldX versions prepend a header row to this file; scan for the
    # first line whose second (tab-separated) column actually parses as a
    # float rather than blindly trusting readline() to be the data row.
    with open(stability_file) as f:
        candidate_lines = [l for l in f.readlines() if l.strip()]
    for line in candidate_lines:
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        try:
            return float(parts[1])
        except ValueError:
            continue  # likely a header row -- keep looking
    raise ValueError(f"Could not find a parseable Total Energy value in FoldX output: {stability_file}")


def run_apbs(pqr_path, output_dir):
    """
    Real APBS electrostatic surface map.
    Requires: conda install -c bioconda -c conda-forge apbs
    Set APBS_BINARY at the top of this file if 'apbs' isn't on your PATH.

    NOTE: grid dimensions below (dime/cglen/fglen) are generic defaults.
    For a production run, tune these to your structure's actual size --
    see APBS documentation for guidance on grid spacing vs. structure extent.
    """
    if shutil.which(APBS_BINARY) is None and not os.path.isfile(APBS_BINARY):
        raise FileNotFoundError(
            f"APBS binary not found ('{APBS_BINARY}'). Set APBS_BINARY at "
            f"the top of this file, or ensure 'apbs' is on your PATH."
        )

    apbs_input_path = os.path.join(output_dir, "apbs_input.in")
    dx_output_prefix = os.path.join(output_dir, "electrostatic_surface")

    apbs_input = f"""read
    mol pqr {pqr_path}
end
elec
    mg-auto
    dime 97 97 97
    cglen 100 100 100
    fglen 80 80 80
    cgcent mol 1
    fgcent mol 1
    mol 1
    lpbe
    bcfl sdh
    pdie 2.0
    sdie 78.54
    srfm smol
    chgm spl2
    sdens 10.0
    srad 1.4
    swin 0.3
    temp 298.15
    calcenergy total
    calcforce no
    write pot dx {dx_output_prefix}
end
quit
"""
    with open(apbs_input_path, 'w') as f:
        f.write(apbs_input)

    cmd = [APBS_BINARY, apbs_input_path]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=output_dir)
    if result.returncode != 0:
        raise RuntimeError(f"APBS failed (exit {result.returncode}):\n{result.stderr}")

    return f"{dx_output_prefix}.dx"


# =============================================================================
# EXTERNALLY-SOURCED RESULTS -- Aggrescan3D, CamSol (structural mode), DeepSol
#
# None of these have a free, one-command-automatable path from this script.
# Rather than fake them with more mocks, this reads real values YOU obtain
# separately and record in manual_results.json.
#
#   - Aggrescan3D: run via the LOCAL standalone install, in your own
#     terminal, separately from this script (not the aggrescan3d.pl web
#     server). Record its results the same way as the other two below.
#   - CamSol (structural mode) and DeepSol: web-submission tools, no
#     local install used.
# =============================================================================

MANUAL_RESULT_KEYS = {
    "aggrescan3d_max_patch_length": "Aggrescan3D max contiguous aggregation-prone patch length (aa)",
    "aggrescan3d_patch_residue_range": "That patch's residue range in the construct, [start, end] (1-based, inclusive) -- pLDDT>=70 is computed from this, not typed in",
    "camsol_structural_max_patch_length": "CamSol (structural mode) max aggregation patch length (aa)",
    "camsol_structural_patch_residue_range": "That patch's residue range in the construct, [start, end] (1-based, inclusive) -- pLDDT>=70 is computed from this, not typed in",
    # NOTE: SOLpro's server is no longer reachable, so this pipeline now uses
    # DeepSol S2 (Khurana et al. 2018) as the sequence-based solubility
    # predictor instead. DeepSol also outputs a genuine 0.0-1.0 probability
    # (unlike e.g. CamSol Intrinsic's unbounded score), so the decision rule
    # below (`< 0.50`) is still valid unchanged -- only the source tool and
    # field name changed.
    "deepsol_probability": "DeepSol S2 predicted solubility probability (0.0-1.0)",
}


def run_manual_prepare(sequence, repaired_pdb_path, manual_results_path, safe_name):
    common.print_banner("EXTERNAL RESULTS NEEDED -- Aggrescan3D, CamSol, DeepSol")
    print("Aggrescan3D: run the LOCAL standalone install in your own terminal")
    print("(not the aggrescan3d.pl web server) on the repaired PDB below.")
    print("CamSol and DeepSol: submit manually through each site.")
    print("Then fill in the values in the template written below.")
    print("-" * 100)

    # Write the FASTA ourselves from the same in-pipeline sequence used
    # everywhere else in Step C, rather than relying on a hand-copied
    # sequence string -- this is the exact sequence loaded from the
    # Phase 1G FASTA earlier in run_step2c_solubility_analysis().
    archive_dir = os.path.dirname(manual_results_path)
    fasta_path = os.path.join(archive_dir, f"{safe_name}.fasta")
    variant_id = safe_name
    os.makedirs(archive_dir, exist_ok=True)
    with open(fasta_path, 'w') as f:
        f.write(f">{variant_id}\n{sequence}\n")
    print(f"[SUCCESS] FASTA written: {fasta_path}")
    print("-" * 100)

    print(f"[Aggrescan3D]  Local standalone install -- run on: {repaired_pdb_path}")
    print("               Record the max patch's residue RANGE (1-based, inclusive construct")
    print("               numbering), not a true/false pLDDT judgment -- pLDDT>=70 is now computed")
    print("               automatically from the real per-residue table parsed before repair.")
    print(f"[CamSol]       https://www-cohsoftware.ch.cam.ac.uk/index.php/camsolstrucorr -- upload: {repaired_pdb_path}")
    print("               Same as above: record the flagged patch's residue range, not a boolean.")
    print(f"[DeepSol]      https://machinelearning-protein.qcri.org  -- upload: {fasta_path}")
    print(f"               (DeepSol replaces SOLpro -- SOLpro's server is no longer reachable.")
    print(f"               On the DeepSol form, set the model parameter to 2 -- DeepSol S2 was")
    print(f"               the best-performing of the three published architectures in the")
    print(f"               original paper, Khurana et al. 2018. NOTE: DeepSol's web server")
    print(f"               requires a free QCAI account before you can submit a job -- sign")
    print(f"               up first if you haven't already. If this URL doesn't work, try")
    print(f"               qcai.qcri.org instead -- QCRI has been consolidating tools there.)")
    print("-" * 100)

    if os.path.isfile(manual_results_path):
        print(f"[INFO] {manual_results_path} already exists -- edit it directly, values are not overwritten.")
    else:
        template = {k: None for k in MANUAL_RESULT_KEYS}
        os.makedirs(os.path.dirname(manual_results_path), exist_ok=True)
        with open(manual_results_path, 'w') as f:
            json.dump(template, f, indent=2)
        print(f"[INFO] Template written to: {manual_results_path}")
        print("[INFO] Fill in each value after running the tools above, then rerun this")
        print("[INFO] script normally (no --manual-prepare flag) to finish Step 2C.")
    print("=" * 100 + "\n")


def load_manual_results(manual_results_path):
    if not os.path.isfile(manual_results_path):
        return None
    with open(manual_results_path) as f:
        data = json.load(f)
    missing = [k for k, v in data.items() if v is None]
    if missing:
        return None
    for key in ("aggrescan3d_patch_residue_range", "camsol_structural_patch_residue_range"):
        value = data.get(key)
        if not (isinstance(value, list) and len(value) == 2 and all(isinstance(x, int) for x in value)):
            raise ValueError(
                f"'{key}' must be a [start, end] list of two integers (1-based, "
                f"inclusive residue numbers) in {manual_results_path}, got: {value!r}"
            )
    return data


# =============================================================================
# SOLUBILITY & STRUCTURAL INTEGRITY DECISION ENGINE
# =============================================================================

def evaluate_structural_solubility(sequence, repaired_pdb_path, sasa_result, total_energy, manual_results, plddt_table):
    deepsol_prob = manual_results["deepsol_probability"]
    agg_patch_len = manual_results["aggrescan3d_max_patch_length"]
    agg_range = manual_results["aggrescan3d_patch_residue_range"]
    camsol_patch_len = manual_results["camsol_structural_max_patch_length"]
    camsol_range = manual_results["camsol_structural_patch_residue_range"]

    # The pLDDT>=70 gate is COMPUTED here from the real per-residue table
    # parsed before repair, never typed in as a raw boolean -- a hand-typed
    # JSON string like "false" is truthy in Python, which used to let a
    # malformed manual entry silently flip this gate the wrong way.
    agg_mean_plddt = region_mean_plddt(plddt_table, agg_range[0], agg_range[1])
    agg_high_plddt = agg_mean_plddt >= 70.0
    camsol_mean_plddt = region_mean_plddt(plddt_table, camsol_range[0], camsol_range[1])
    camsol_high_plddt = camsol_mean_plddt >= 70.0

    result = {
        "DeepSol_Prob": deepsol_prob,
        "Agg_Patch_Len": agg_patch_len,
        "Agg_Patch_Mean_pLDDT": round(agg_mean_plddt, 1),
        "Agg_High_pLDDT": agg_high_plddt,
        "CamSol_Patch_Len": camsol_patch_len,
        "CamSol_Patch_Mean_pLDDT": round(camsol_mean_plddt, 1),
        "CamSol_High_pLDDT": camsol_high_plddt,
        "Hydrophobic_Fraction": sasa_result["hydrophobic_fraction"],
        "Exposed_Area": sasa_result["exposed_patch_area"],
        "FoldX_Total_Energy_kcal_mol": total_energy,
        "Status": "PASS",
        # Sec. II.C's ONLY hard REJECT rule is the conjunctive DeepSol +
        # Aggrescan3D-patch test below; everything else it lists ("additional
        # structural flags") is informational and drives REVIEW, not REJECT.
        # These were previously merged into one Reason_Code string, which
        # read as if hydrophobic fraction / exposed area / CamSol alone had
        # failed the construct -- split so the CSV/report distinguish what
        # actually rejected the construct from what merely flagged it.
        "Rejection_Reasons": [],
        "Structural_Flags": [],
    }

    is_insoluble = deepsol_prob < 0.50
    has_confident_agg_patch = agg_patch_len > 8 and agg_high_plddt

    # Sec. II.C (Step 2C, structural solubility) is explicitly conjunctive:
    # "constructs were flagged for rejection when DeepSol S2 predicted
    # insolubility (probability < 0.50) AND a structure-aware tool
    # (Aggrescan3D) identified an exposed, high-scoring aggregation patch
    # > 8 contiguous residues in a pLDDT >= 70 region." This is a
    # DIFFERENT, narrower rule than Sec. II.A.b's Step 2A sequence-only
    # rule ("sequences flagged only with 'predicted insoluble' were
    # rejected") -- the two steps must not share one rule.
    if is_insoluble and has_confident_agg_patch:
        result["Status"] = "REJECT"
        result["Rejection_Reasons"].append("DeepSol Insoluble (<0.50) + Aggrescan3D Confident Patch >8aa (pLDDT>=70)")
    elif is_insoluble:
        result["Status"] = "REVIEW"
        result["Structural_Flags"].append("DeepSol Insoluble (<0.50) -- no confident Aggrescan3D patch corroboration")
    elif has_confident_agg_patch:
        # Mirrors CamSol's equivalent check below -- a confident patch
        # without DeepSol insolubility does not meet Sec. II.C's
        # conjunctive REJECT rule, but is still a real structural signal.
        result["Status"] = "REVIEW"
        result["Structural_Flags"].append("Aggrescan3D Confident Patch >8aa (pLDDT>=70)")

    if result["Status"] != "REJECT":
        if sasa_result["hydrophobic_fraction"] > 0.25:
            result["Status"] = "REVIEW"
            result["Structural_Flags"].append("Hydrophobic Fraction > 0.25")
        if sasa_result["exposed_patch_area"] > 250.0:
            result["Status"] = "REVIEW"
            result["Structural_Flags"].append("Exposed Hydrophobic Area > 250 \u00c5\u00b2")
        if camsol_patch_len > 8 and camsol_high_plddt:
            result["Status"] = "REVIEW"
            result["Structural_Flags"].append("CamSol (structural) Aggregation Patch > 8aa")

    if not result["Rejection_Reasons"]:
        result["Rejection_Reasons"].append("N/A")
    if not result["Structural_Flags"]:
        result["Structural_Flags"].append("None")

    result["Rejection_Reasons"] = " | ".join(result["Rejection_Reasons"])
    result["Structural_Flags"] = " | ".join(result["Structural_Flags"])
    return result


def run_step2c_solubility_analysis():
    start_time = time.time()
    project_root = _PROJECT_ROOT

    input_csv_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA", "Filtered")
    variant_fasta_path = common.phase1g_fasta_path(project_root)
    pdb_input_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepB", "Tertiary_Structure")

    output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepC")
    archive_dir = os.path.join(output_base, "Supplementary_Archive")
    os.makedirs(output_base, exist_ok=True)
    os.makedirs(archive_dir, exist_ok=True)

    common.print_banner("PHASE 2 STEP C: 3D SOLUBILITY & STRUCTURAL INTEGRITY")
    print(f"[INFO] Resolved Project Root : {project_root}")

    winner_row, _ = common.get_winner_from_filtered_csv(input_csv_dir)
    if winner_row is None:
        return
    winner_name = winner_row["Variant"]

    variants = common.load_multi_fasta(variant_fasta_path)
    if variants is None:
        print(f"[ERROR] Variant FASTA file not found: {variant_fasta_path}")
        return
    sequence = common.lookup_sequence(variants, winner_name)
    if sequence is None:
        print(f"[ERROR] Could not find a matching header for '{winner_name}' in: {variant_fasta_path}")
        return

    safe_name = common.sanitize_variant_name(winner_name)
    input_cif_path = os.path.join(pdb_input_dir, f"AF3_Target_{safe_name}.cif")

    if not os.path.isfile(input_cif_path):
        print(f"[ERROR] Structure file not found at {input_cif_path}.")
        print("[ERROR] Run Step 2B first (including the AlphaFold Server import) so this model exists.")
        return

    print(f"[INFO] Target Variant : {winner_name}")
    print(f"[INFO] Input Structure: {input_cif_path}")
    print("-" * 110)

    # 0. Parse real per-residue pLDDT from the ORIGINAL AlphaFold mmCIF,
    # before repair/minimization/refinement rewrite the B-factor column
    # and make it unrecoverable. Persisted as its own table so the
    # pLDDT>=70 gate below is computed from real data, never typed in.
    try:
        plddt_table = parse_plddt_per_residue(input_cif_path)
        plddt_csv_path = os.path.join(archive_dir, f"{safe_name}_pLDDT_per_residue.csv")
        with open(plddt_csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(["Residue_Index", "pLDDT", "High_Confidence_>=70"])
            for seq_id in sorted(plddt_table):
                w.writerow([seq_id, round(plddt_table[seq_id], 2), plddt_table[seq_id] >= 70.0])
        mean_plddt = sum(plddt_table.values()) / len(plddt_table)
        print(f"[SUCCESS] Per-residue pLDDT    : {os.path.relpath(plddt_csv_path, output_base)} "
              f"({len(plddt_table)} residues, mean {mean_plddt:.1f})")
    except (ValueError, OSError) as e:
        print(f"[ERROR] Failed to parse per-residue pLDDT: {e}")
        return

    # 1a. Repair (real, PDBFixer) -- adds missing atoms/residues only
    repaired_pdb_path = os.path.join(archive_dir, f"{safe_name}_repaired.pdb")
    try:
        repair_structure_pdbfixer(input_cif_path, repaired_pdb_path)
        print(f"[SUCCESS] Repaired structure   : {os.path.relpath(repaired_pdb_path, output_base)}")
    except ImportError:
        print("[ERROR] PDBFixer/OpenMM not installed. Run: pip install pdbfixer openmm")
        return
    except Exception as e:
        print(f"[ERROR] PDBFixer repair failed: {e}")
        return

    # 1b. Minimize (real, OpenMM) -- resolves steric clashes left over from
    # repair alone. All downstream steps use this minimized structure, not
    # the raw-repaired one, since un-minimized clashes can otherwise blow
    # up FoldX ddG into unrealistic large-positive values.
    minimized_pdb_path = os.path.join(archive_dir, f"{safe_name}_minimized.pdb")
    try:
        minimize_structure_openmm(repaired_pdb_path, minimized_pdb_path)
        print(f"[SUCCESS] Minimized structure   : {os.path.relpath(minimized_pdb_path, output_base)}")
    except ImportError:
        print("[ERROR] OpenMM not installed. Run: pip install openmm")
        return
    except Exception as e:
        print(f"[ERROR] OpenMM minimization failed: {e}")
        return
    repaired_pdb_path = minimized_pdb_path  # OpenMM-minimized structure feeds into Phenix next

    # 1c. Refine (real, Phenix) -- fixes bad angles, CaBLAM outliers, and
    # the Ramachandran distribution beyond what generic OpenMM minimization
    # alone reliably achieves. This is now a standard part of structure
    # prep, not a manual afterthought done after Step 2D -- everything
    # downstream (SASA, FoldX, manual Aggrescan3D/CamSol/DeepSol
    # submissions, and Step 2D's MolProbity validation) all operates on
    # THIS final refined structure, so no results describe an earlier,
    # geometrically flawed draft.
    phenix_second_pass = os.environ.get("PHENIX_SECOND_PASS", "0") == "1"
    try:
        refined_pdb_path = refine_structure_phenix(repaired_pdb_path, archive_dir, safe_name, second_pass=phenix_second_pass)
        print(f"[SUCCESS] Phenix-refined struct : {os.path.relpath(refined_pdb_path, output_base)}"
              f"{' (2-pass)' if phenix_second_pass else ' (1-pass)'}")
    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        print("[ERROR] Confirm Phenix's command-line tools are on PATH in this terminal")
        print("[ERROR] (the same environment you configured when running Phenix manually before).")
        return
    except Exception as e:
        print(f"[ERROR] Phenix refinement failed: {e}")
        return
    repaired_pdb_path = refined_pdb_path  # everything downstream uses the Phenix-refined structure

    # 2. Protonate at pH 7.0 (real, PDB2PQR)
    protonated_pqr_path = os.path.join(archive_dir, f"{safe_name}_pH7.0.pqr")
    try:
        protonate_pdb2pqr(repaired_pdb_path, protonated_pqr_path, ph=7.0)
        print(f"[SUCCESS] Protonated (pH 7.0): {os.path.relpath(protonated_pqr_path, output_base)}")
    except FileNotFoundError:
        print("[ERROR] pdb2pqr30 not found on PATH. Run: pip install pdb2pqr")
        return
    except Exception as e:
        print(f"[ERROR] PDB2PQR protonation failed: {e}")
        return

    # 3. SASA + hydrophobic fraction (real, FreeSASA)
    try:
        sasa_result = analyze_sasa_freesasa(repaired_pdb_path)
        print(f"[SUCCESS] SASA analysis complete (hydrophobic fraction: {sasa_result['hydrophobic_fraction']:.3f})")
    except ImportError:
        print("[ERROR] freesasa not installed. Run: pip install freesasa")
        return
    except Exception as e:
        print(f"[ERROR] FreeSASA analysis failed: {e}")
        return

    # 4. FoldX Total Folding Energy (real, requires configured binary path)
    # -- NOT a true ddG; see run_foldx()'s docstring.
    try:
        total_energy = run_foldx(repaired_pdb_path, archive_dir)
        print(f"[SUCCESS] FoldX Total Energy: {total_energy:.2f} kcal/mol")
    except Exception as e:
        print(f"[ERROR] FoldX failed: {e}")
        print("[INFO] Check FOLDX_BINARY at the top of this file, and that you've completed")
        print("[INFO] the free academic registration/download at foldxsuite.crg.eu")
        return

    # 5. APBS electrostatic map (real, requires configured binary path)
    try:
        apbs_map_path = run_apbs(protonated_pqr_path, archive_dir)
        print(f"[SUCCESS] APBS map: {os.path.relpath(apbs_map_path, output_base)}")
    except Exception as e:
        print(f"[ERROR] APBS failed: {e}")
        print("[INFO] Check APBS_BINARY at the top of this file, and that APBS is installed")
        print("[INFO] (conda install -c bioconda -c conda-forge apbs)")
        return

    # 6. Aggrescan3D, CamSol (structural), DeepSol -- manual (no free local/API option)
    manual_results_path = os.path.join(archive_dir, f"{safe_name}_manual_results.json")
    manual_results = load_manual_results(manual_results_path)
    if manual_results is None:
        run_manual_prepare(sequence, repaired_pdb_path, manual_results_path, safe_name)
        return

    print(f"[SUCCESS] Manual results loaded from: {os.path.relpath(manual_results_path, output_base)}")
    print("-" * 110)

    results = evaluate_structural_solubility(sequence, repaired_pdb_path, sasa_result, total_energy, manual_results, plddt_table)

    print(f"{'METRIC':<30} | {'VALUE':<15} | {'THRESHOLD / TARGET'}")
    print("-" * 110)
    print(f"{'DeepSol Probability':<30} | {results['DeepSol_Prob']:<15.2f} | >= 0.50")
    print(f"{'Aggrescan3D Patch (pLDDT>=70)':<30} | {results['Agg_Patch_Len']:<15} | <= 8 aa (mean pLDDT {results['Agg_Patch_Mean_pLDDT']})")
    print(f"{'CamSol (structural) Patch':<30} | {results['CamSol_Patch_Len']:<15} | <= 8 aa (mean pLDDT {results['CamSol_Patch_Mean_pLDDT']})")
    print(f"{'Hydrophobic Fraction':<30} | {results['Hydrophobic_Fraction']:<15.2f} | < 0.25")
    print(f"{'Exposed Hydro Area':<30} | {results['Exposed_Area']:<15.1f} | < 250 \u00c5\u00b2")
    print(f"{'FoldX Total Folding Energy':<30} | {results['FoldX_Total_Energy_kcal_mol']:<15.2f} | (contextual -- not \u0394\u0394G, no designed mutants)")
    print("-" * 110)
    print(f"FINAL DECISION : [{results['Status']}]")
    print(f"  Rejection_Reasons : {results['Rejection_Reasons']}")
    print(f"  Structural_Flags  : {results['Structural_Flags']}")
    print("-" * 110)

    ts = datetime.now().strftime("%Y%m%d_%H%M")
    report_path = os.path.join(output_base, f"Step2C_Solubility_Report_{ts}.csv")

    results["Repaired_PDB_File"] = os.path.relpath(repaired_pdb_path, output_base)
    results["Protonated_PQR_File"] = os.path.relpath(protonated_pqr_path, output_base)
    results["APBS_Map_File"] = os.path.relpath(apbs_map_path, output_base)
    results["Manual_Results_File"] = os.path.relpath(manual_results_path, output_base)

    csv_data = {"Variant": winner_name, **results}
    with open(report_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=csv_data.keys())
        w.writeheader()
        w.writerow(csv_data)

    total_time = common.format_time(time.time() - start_time)
    common.print_banner("STEP 2C COMPLETE")
    print("[SUCCESS] Structural solubility & aggregation propensity checked with real tools.")
    print(f"[SUCCESS] Execution Time : {total_time}")
    print(f"[INFO] Report Saved      : {os.path.relpath(report_path, project_root)}")
    print("=" * 110 + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 2C: 3D Solubility & Structural Integrity")
    parser.add_argument("--manual-prepare", action="store_true",
                         help="Print manual submission instructions and write the results template early")
    args = parser.parse_args()

    if args.manual_prepare:
        project_root = _PROJECT_ROOT
        input_csv_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA", "Filtered")
        variant_fasta_path = common.phase1g_fasta_path(project_root)
        winner_row, _ = common.get_winner_from_filtered_csv(input_csv_dir)
        if winner_row:
            variants = common.load_multi_fasta(variant_fasta_path)
            sequence = common.lookup_sequence(variants, winner_row["Variant"]) if variants else None
            safe_name = common.sanitize_variant_name(winner_row["Variant"])
            archive_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepC", "Supplementary_Archive")
            # Use the Phenix-refined structure -- whichever pass actually
            # ran (2-pass takes priority if both exist, since it's the
            # more complete refinement).
            pass2_path = os.path.join(archive_dir, f"{safe_name}_phenix_pass2.pdb")
            pass1_path = os.path.join(archive_dir, f"{safe_name}_phenix_pass1.pdb")
            if os.path.isfile(pass2_path):
                repaired_pdb_path = pass2_path
            elif os.path.isfile(pass1_path):
                repaired_pdb_path = pass1_path
            else:
                print(f"[ERROR] No Phenix-refined structure found in {archive_dir}.")
                print("[ERROR] Run this script normally (no --manual-prepare) first so Phenix refinement completes.")
                sequence = None  # skip the call below
            manual_results_path = os.path.join(archive_dir, f"{safe_name}_manual_results.json")
            if sequence:
                run_manual_prepare(sequence, repaired_pdb_path, manual_results_path, safe_name)
    else:
        run_step2c_solubility_analysis()
