import os
import sys
from datetime import datetime

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
# PHASE 1De FOLD PREPARATION -- native antigen structures for SEMA-3D.
#
# Sec. I.D: "native tertiary structures of each source antigen were obtained or
# predicted prior to chimeric assembly ... Discontinuous B-cell epitopes were
# predicted on each native fold using SEMA 2.0 ... directly from AlphaFold-
# predicted or experimentally solved structures." SEMA-3D therefore needs ONE
# STRUCTURE PER ANTIGEN, and this step produces the inputs for them.
#
# WHY THE SEQUENCES ARE SLICED. Phase 1A's four HIV "targets" are four Entrez
# query strings, and NCBI returned POLYPROTEINS for all of them: gp120 and gp41
# share one 856-aa Env record, p17 and p24 share one 1437-aa Gag-Pol record.
# Folding those whole would be wrong twice over -- gp120 and gp41 would get
# identical structures (and therefore identical SEMA scores), and a 1437-aa
# Gag-Pol never folds as one unit in vivo, so its model would be the same kind
# of undetermined coil that deviation #21 documents for the construct itself.
# Each HIV target is therefore cut to its MATURE SUBUNIT using the ranges in
# phase1_common.HIV_SUBUNIT_RANGES. The three Mpox targets are already single
# mature proteins and are emitted whole.
#
# WHY ALPHAFOLD AND NOT A SOLVED STRUCTURE. The paper prefers RCSB "where
# available". For Env the obvious candidate, 8TTW, is a BG505 SOSIP.664
# stabilised trimer -- clade A, not the CRF01_AE target -- carrying temsavir
# and two Fabs on exactly the surfaces SEMA scores, the same contamination
# Phase 3A had to strip from 6NIG. Reading "available" as available FOR THIS
# ANTIGEN rather than for any homolog, an AlphaFold model of the actual
# CRF01_AE subunit is the more defensible input. Recorded as a deviation.
#
# THIS STEP DOES NOT FOLD ANYTHING. AlphaFold Server is a manual submission,
# exactly as in Phase 2A/2B. It writes the FASTAs and the instructions, then
# ingests whatever comes back.
# =============================================================================

TARGETS = ("Mpox_L1R", "Mpox_B5R", "Mpox_A35R",
           "HIV_gp120", "HIV_gp41", "HIV_p24", "HIV_p17")


def _var01(target):
    folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1A")
    matches = sorted(f for f in os.listdir(folder) if f.startswith(f"{target}_Var_01_"))
    if not matches:
        return None, None
    with open(os.path.join(folder, matches[0])) as f:
        seq = "".join(l.strip() for l in f if not l.startswith(">")).upper()
    return seq, matches[0]


def subunit_sequence(target):
    """
    Returns (sequence_to_fold, source_description). For HIV targets this is the
    mature subunit slice; for Mpox it is the whole record.
    """
    parent_seq, fname = _var01(target)
    if parent_seq is None:
        return None, None
    info = common.HIV_SUBUNIT_RANGES.get(target)
    if info is None:
        return parent_seq, f"{fname} (whole record -- already a mature protein)"
    s, e = info["start"], info["end"]
    return parent_seq[s - 1:e], f"{fname} residues {s}-{e} ({info['parent']} subunit)"


def run_fold_prep():
    common.print_banner("PHASE 1De FOLD PREP: NATIVE ANTIGEN STRUCTURES FOR SEMA-3D")

    folds_dir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1De", "Folds")
    sub_dir = os.path.join(folds_dir, "_submissions")
    os.makedirs(sub_dir, exist_ok=True)

    prepared, missing_structs = [], []
    for target in TARGETS:
        seq, source = subunit_sequence(target)
        if seq is None:
            print(f"[WARNING] No Var_01 record for {target} -- skipped.")
            continue
        fasta_path = os.path.join(sub_dir, f"{target}_subunit.fasta")
        with open(fasta_path, "w") as f:
            f.write(f">{target} | {source} | {len(seq)} aa\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")

        # Has a structure already been ingested for this target?
        tgt_dir = os.path.join(folds_dir, target)
        have = []
        if os.path.isdir(tgt_dir):
            have = [f for f in os.listdir(tgt_dir) if f.lower().endswith((".pdb", ".cif"))]
        status = "READY" if have else "AWAITING STRUCTURE"
        if not have:
            missing_structs.append(target)
        prepared.append((target, len(seq), source, status))
        print(f"  {target:<11} {len(seq):>5} aa  {status:<19} {source}")

    # ---- submission instructions, in the house style of Phase 2A ------------
    howto = os.path.join(folds_dir, "HOW_TO_SUBMIT.txt")
    with open(howto, "w") as f:
        f.write("HOW TO PRODUCE THE NATIVE ANTIGEN STRUCTURES FOR SEMA-3D\n")
        f.write("=" * 60 + "\n\n")
        f.write("WHY: Sec. I.D requires conformational B-cell prediction on the native fold of\n")
        f.write("each source antigen. SEMA-3D takes a structure, not a sequence.\n\n")
        f.write("WHAT TO SUBMIT: one AlphaFold Server monomer job per FASTA in _submissions/.\n")
        f.write("These are MATURE SUBUNITS, not the raw Phase 1A records -- gp120/gp41 share one\n")
        f.write("Env polyprotein and p17/p24 share one Gag-Pol polyprotein, so folding the raw\n")
        f.write("records would give two antigens the same structure and the same SEMA scores.\n\n")
        for target, n, source, _status in prepared:
            f.write(f"  {target:<11} {n:>5} aa   {source}\n")
        f.write("\nWHERE: https://alphafoldserver.com/  (same tool and workflow as Phase 2B)\n\n")
        f.write("WHAT TO DO WITH EACH RESULT\n")
        f.write("---------------------------\n")
        f.write("Download the job's ZIP, and put the model file for each antigen at:\n\n")
        f.write("    Step_Outputs/Phase1/Phase1D/Phase1De/Folds/<TARGET>/\n\n")
        f.write("Keep the whole AlphaFold download (all five ranked models plus the confidence\n")
        f.write("JSONs), not just the top model -- deviation #21 exists because Step 2B kept only\n")
        f.write("model_0 and the fold-confidence problem stayed invisible for a whole phase.\n\n")
        f.write("Phase1De_conformational.py picks these up automatically. Any antigen without a\n")
        f.write("structure stays SEMA_Corroborated=UNSCREENED rather than blocking the step.\n\n")
        f.write("NOTE ON CONFIDENCE: record each model's mean pLDDT. A low-confidence fold\n")
        f.write("produces confident-looking SEMA patches that mean nothing, which is the exact\n")
        f.write("trap deviation #21 documents for the construct's own model.\n")

    print(f"\n[INFO] {len(prepared)} subunit FASTA(s) written to: {sub_dir}")
    print(f"[INFO] Submission instructions: {howto}")
    if missing_structs:
        print(f"[ACTION REQUIRED, NOT BLOCKING] {len(missing_structs)} antigen(s) have no structure yet: "
              f"{missing_structs}")
    else:
        print("[INFO] All antigens have a structure -- Phase 1De can run SEMA-3D.")
    print("=" * 90 + "\n")


if __name__ == "__main__":
    run_fold_prep()
