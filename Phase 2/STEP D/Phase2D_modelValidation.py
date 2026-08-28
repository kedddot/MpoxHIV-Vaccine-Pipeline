import os
import sys
import csv
import json
import time
import glob
import shutil
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
# MOLPROBITY -- run LOCALLY via phenix.molprobity, web server as fallback
#
# Phenix 2.2 is installed on this machine and bundles the Richardson Lab
# MolProbity code, so `phenix.molprobity` produces the same six validation
# numbers the web server reports -- same tool, run locally. That removes
# the last manual web submission in Phase 2 AND the transcription step
# where a value could be mistyped.
#
# The manual molprobity.biochem.duke.edu workflow is RETAINED as an
# automatic fallback: if Phenix is missing or its run fails, the script
# writes the template and asks for the values as before. Sec. II.D names
# "MolProbity" without mandating the web interface, so neither path is a
# deviation.
#
# NOTE: Rampage is intentionally NOT used here. MolProbity's own summary
# report already includes Ramachandran favored/allowed/outlier percentages
# directly, so running a second, separate tool for the same numbers would
# be redundant.
# =============================================================================

MANUAL_RESULT_KEYS = {
    "molprobity_score": "MolProbity score (from the summary table)",
    "clashscore": "All-atom clashscore",
    "poor_rotamers_percent": "Poor rotamers (%)",
    "rama_favored_percent": "Ramachandran favored (%)",
    "rama_allowed_percent": "Ramachandran allowed (%)",
    "rama_outlier_percent": "Ramachandran outliers (%)",
}


MOLPROBITY_BINARY = os.environ.get("MOLPROBITY_BINARY", "phenix.molprobity")


def run_phenix_molprobity(structure_path):
    """
    Runs MolProbity LOCALLY via phenix.molprobity, which ships with a
    Phenix install (confirmed present: Phenix 2.2-6143).

    This replaces the manual molprobity.biochem.duke.edu web submission
    with the same underlying tool -- MolProbity IS the Richardson Lab
    code that Phenix bundles, so this is not a substitution of one
    program for another, it is the same validation run locally. Sec. II.D
    names "MolProbity" without mandating the web interface.

    Returns a dict with the six values Sec. II.D needs, or None if Phenix
    isn't installed / the run fails -- in which case the caller falls
    back to the manual web workflow rather than blocking.
    """
    import tempfile
    import shutil as _shutil
    import re as _re

    if shutil.which(MOLPROBITY_BINARY) is None:
        print(f"[INFO] {MOLPROBITY_BINARY} not on PATH -- falling back to manual MolProbity entry.")
        return None

    with tempfile.TemporaryDirectory() as scratch:
        # Copy in under a plain name: Phenix parses every argument as a
        # `key=value` PHIL parameter, so any "=" in the path is misread
        # (this is the same class of bug that broke Step 2C's Phenix and
        # APBS calls -- see sanitize_variant_name's docstring).
        local_model = os.path.join(scratch, "model.pdb")
        _shutil.copy(structure_path, local_model)
        try:
            proc = subprocess.run(
                [MOLPROBITY_BINARY, "model.pdb"],
                capture_output=True, text=True, cwd=scratch, timeout=3600,
            )
        except (OSError, subprocess.TimeoutExpired) as e:
            print(f"[WARN] phenix.molprobity could not be run ({e}) -- falling back to manual entry.")
            return None
        if proc.returncode != 0:
            print(f"[WARN] phenix.molprobity failed (exit {proc.returncode}):\n{proc.stderr[:400]}")
            return None

        out_path = os.path.join(scratch, "molprobity.out")
        text = open(out_path).read() if os.path.isfile(out_path) else proc.stdout
        if not text.strip():
            print("[WARN] phenix.molprobity produced no output -- falling back to manual entry.")
            return None

        def _num(pattern):
            m = _re.search(pattern, text)
            return float(m.group(1)) if m else None

        # Rotamer/Ramachandran both have "Outliers :" lines, so the
        # Ramachandran block is isolated first rather than matched loosely.
        rama_block = ""
        m = _re.search(r"Ramachandran Plot:(.*?)Rotamer:", text, _re.S)
        if m:
            rama_block = m.group(1)

        def _in_block(block, label):
            mm = _re.search(rf"{label}\s*:\s*([\d.]+)\s*%", block)
            return float(mm.group(1)) if mm else None

        results = {
            "molprobity_score": _num(r"MolProbity score\s*=\s*([\d.]+)"),
            "clashscore": _num(r"Clashscore\s*=\s*([\d.]+)"),
            "poor_rotamers_percent": _num(r"Rotamer outliers\s*=\s*([\d.]+)\s*%"),
            "rama_favored_percent": _in_block(rama_block, "Favored"),
            "rama_allowed_percent": _in_block(rama_block, "Allowed"),
            "rama_outlier_percent": _in_block(rama_block, "Outliers"),
        }

        missing = [k for k, v in results.items() if v is None]
        if missing:
            print(f"[WARN] phenix.molprobity output missing {missing} -- falling back to manual entry.")
            return None

        # Archive the full report next to the structure for the record.
        archive = os.path.join(os.path.dirname(structure_path),
                               os.path.basename(structure_path).replace(".pdb", "") + "_molprobity.out")
        try:
            _shutil.copy(out_path, archive)
            print(f"[INFO] Full MolProbity report archived: {os.path.basename(archive)}")
        except OSError:
            pass
        return results


def run_manual_prepare(structure_path, manual_results_path):
    common.print_banner("EXTERNAL RESULTS NEEDED -- MolProbity")
    print("Submit the structure below to the MolProbity web server, then fill")
    print("in the values in the template written below.")
    print("-" * 100)
    print("[MolProbity]  https://molprobity.biochem.duke.edu")
    print(f"              Upload: {structure_path}")
    print("              Run the full validation (default options are fine).")
    print("              Read the 6 values below off its summary table.")
    print("-" * 100)

    if os.path.isfile(manual_results_path):
        print(f"[INFO] {manual_results_path} already exists -- edit it directly, values are not overwritten.")
    else:
        template = {k: None for k in MANUAL_RESULT_KEYS}
        os.makedirs(os.path.dirname(manual_results_path), exist_ok=True)
        with open(manual_results_path, 'w') as f:
            json.dump(template, f, indent=2)
        print(f"[INFO] Template written to: {manual_results_path}")
        print("[INFO] Fill in each value after running MolProbity, then rerun this")
        print("[INFO] script normally (no --manual-prepare flag) to finish Step 2D.")
    print("=" * 100 + "\n")


def load_manual_results(manual_results_path):
    if not os.path.isfile(manual_results_path):
        return None
    with open(manual_results_path) as f:
        data = json.load(f)
    missing = [k for k, v in data.items() if v is None]
    if missing:
        return None
    return data


def _resolve_validation_target(project_root):
    """
    Shared by both the normal run and --manual-prepare (previously each
    re-derived this inline, byte-for-byte duplicated). Returns
    (winner_name, safe_name, structure_path, pass_label), or None if no
    filtered candidate exists yet.

    pass_label ("pass1"/"pass2") is folded into the manual results
    filename downstream so that switching which Phenix pass produced the
    structure being validated naturally invalidates any manual MolProbity
    entry recorded against the OTHER pass, instead of silently reporting
    pass1 numbers against a pass2 structure path.
    """
    input_csv_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA", "Filtered")
    stepc_archive_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepC", "Supplementary_Archive")

    winner_row, _ = common.get_winner_from_filtered_csv(input_csv_dir)
    if winner_row is None:
        return None
    winner_name = winner_row["Variant"]
    safe_name = common.sanitize_variant_name(winner_name)

    pass2_path = os.path.join(stepc_archive_dir, f"{safe_name}_phenix_pass2.pdb")
    pass1_path = os.path.join(stepc_archive_dir, f"{safe_name}_phenix_pass1.pdb")
    if os.path.isfile(pass2_path):
        structure_path, pass_label = pass2_path, "pass2"
    else:
        structure_path, pass_label = pass1_path, "pass1"

    return winner_name, safe_name, structure_path, pass_label

# =============================================================================
# FOLD CONFIDENCE (from Step 2B) -- a SECOND, INDEPENDENT AXIS
#
# MolProbity answers "is the local geometry right?" -- bond lengths, angles,
# clashes, rotamers, Ramachandran. It does NOT answer "is the fold right?".
# A random coil with clean geometry scores well on every metric below, which
# is exactly what happened here: this construct is VALIDATED at MolProbity
# 0.64 while 88.6% of its residues sit below pLDDT 50 and its pTM is 0.17.
#
# Phase III (Sec. III.A) docks this model into TLR-2/TLR-4. That needs the
# FOLD to be determined, not just the geometry -- docking an undetermined
# region means docking one arbitrary conformer, and the answer changes when
# a different ranked AlphaFold model is used. So the two claims are reported
# as two separate verdicts that cannot be read as one:
#
#   Overall_Status    -- stereochemistry (methodology-defined, UNCHANGED)
#   Docking_Readiness -- fold confidence (external convention, see #21)
#
# NON-BLOCKING: if Step 2B has not produced a confidence report, every field
# reads NOT_ASSESSED/UNKNOWN and the step completes normally -- the same
# tolerant pattern Phase 1De uses for SEMA's UNSCREENED. Absence of the
# analysis must never look like a failed analysis.
# =============================================================================

def _fmt_conf(value, spec, unit=""):
    """
    Confidence values are legitimately absent when Step 2B has not run.
    The unit is dropped in that case -- "NOT_ASSESSED A" reads like a number.
    """
    if not isinstance(value, (int, float)):
        return "NOT_ASSESSED"
    return format(value, spec) + unit


def load_fold_confidence(project_root):
    """
    Reads the newest Step2B_AF3_Confidence_*.csv, or returns None.

    Returns a flat dict of the fields Step 2D reports. Never raises: a
    missing, partial or unparseable report degrades to None so the
    stereochemical validation still runs.
    """
    folder = os.path.join(project_root, "Step_Outputs", "Phase2", "StepB", "Model_Confidence")
    if not os.path.isdir(folder):
        return None
    reports = sorted(glob.glob(os.path.join(folder, "Step2B_AF3_Confidence_*.csv")))
    if not reports:
        return None
    path = reports[-1]

    def _f(value):
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    try:
        with open(path, newline="") as fh:
            rows = list(csv.DictReader(fh))
    except OSError:
        return None
    if not rows:
        return None

    out = {"Fold_Confidence_Report": os.path.basename(path)}
    ptms, rmsds = [], []
    for r in rows:
        scope = r.get("Scope")
        if scope == "WHOLE_MODEL":
            out["Mean_pLDDT"] = _f(r.get("Mean_pLDDT"))
            out["Pct_Residues_pLDDT_ge70"] = _f(r.get("Pct_pLDDT_ge70"))
            out["Longest_Confident_Run"] = r.get("Longest_Confident_Run") or ""
        elif scope == "MODEL":
            ptm = _f(r.get("pTM"))
            if ptm is not None:
                ptms.append(ptm)
        elif scope == "REGION":
            # Only regions that are actually confident tell us anything about
            # docking; a coil's cross-model RMSD is large by definition and
            # would swamp the maximum if it were included here.
            mean_plddt = _f(r.get("Mean_pLDDT"))
            rmsd = _f(r.get("CrossModel_RMSD_Mean"))
            if mean_plddt is not None and mean_plddt >= 70.0 and rmsd is not None:
                rmsds.append(rmsd)
        elif scope == "VERDICT":
            out["Docking_Readiness"] = r.get("Docking_Readiness") or "UNKNOWN"
            out["Docking_Ready_Regions"] = r.get("Docking_Ready_Regions") or ""
            out["Docking_Readiness_Basis"] = r.get("Notes") or ""

    out["pTM"] = max(ptms) if ptms else None
    out["Max_CrossModel_RMSD_Confident_Region"] = max(rmsds) if rmsds else None
    return out


# =============================================================================
# CORE STEREOCHEMICAL VALIDATION DECISION ENGINE
# =============================================================================

def evaluate_model_quality(manual_results, fold_confidence=None):
    """
    Evaluates stereochemical quality based on Phase II Step D rules.

    fold_confidence, if supplied by load_fold_confidence(), is reported
    alongside as a SEPARATE axis. It deliberately does NOT feed into
    Overall_Status: the methodology defines that verdict as stereochemistry,
    and quietly redefining it here would misrepresent the paper's own
    criterion. See the Docking_Readiness fields instead.
    """

    result = {
        "MolProbity_Score": manual_results["molprobity_score"],
        "MolProbity_Clashscore": manual_results["clashscore"],
        "MolProbity_Poor_Rotamers_Pct": manual_results["poor_rotamers_percent"],
        "Rama_Favored": manual_results["rama_favored_percent"],
        "Rama_Allowed": manual_results["rama_allowed_percent"],
        "Rama_Outlier": manual_results["rama_outlier_percent"],
        "MP_Quality": "UNKNOWN",
        "Rama_Label": "UNKNOWN",
        "Rama_Meets_Target": False,
        "Overall_Status": "REVIEW",
    }

    # 1. MolProbity score bands, per methodology: 0.5-1.5 excellent,
    #    <2 acceptable, >3 poor.
    #    NOTE: the methodology text does not define the 2.0-3.0 range.
    #    Rather than silently folding that gap into either ACCEPTABLE or
    #    POOR, it is surfaced explicitly as MARGINAL so a real candidate
    #    landing there gets flagged for manual review instead of an
    #    unstated assumption deciding its fate.
    score = result["MolProbity_Score"]
    if score <= 1.5:
        # A score below the stated 0.5 floor is not worse than "excellent"
        # -- it's a BETTER (lower) score than the excellent band's own
        # upper reference point (0 = a perfect structure). The previous
        # `0.5 <= score` lower bound meant a genuinely great score like
        # 0.3 fell through to the next branch and was mislabeled ACCEPTABLE.
        result["MP_Quality"] = "EXCELLENT"
    elif score < 2.0:
        result["MP_Quality"] = "ACCEPTABLE"
    elif score <= 3.0:
        result["MP_Quality"] = "MARGINAL (2.0-3.0 undefined by methodology)"
    else:
        result["MP_Quality"] = "POOR"

    # 2. Ramachandran: methodology's stated preferred target is 95-98%
    #    favored residues.
    favored = result["Rama_Favored"]
    result["Rama_Meets_Target"] = favored >= 95.0
    if 95.0 <= favored <= 98.0:
        result["Rama_Label"] = "MEETS TARGET (95-98%)"
    elif favored > 98.0:
        result["Rama_Label"] = "EXCEEDS TARGET (>98%)"
    else:
        result["Rama_Label"] = "BELOW TARGET (<95%)"

    # 2b. Clashscore, poor rotamers %, and Rama outliers % were previously
    #     computed and printed as "(contextual)" but never affected
    #     Overall_Status at all -- a structure with e.g. a clashscore of 80
    #     could still be marked VALIDATED. The methodology gives no
    #     explicit thresholds for these (unlike MolProbity score and Rama
    #     favored), so widely-used MolProbity community conventions (NOT
    #     from the paper -- labeled as such) are applied as defensive
    #     REVIEW triggers: a clearly bad value on any of these is now at
    #     least surfaced, never silently ignored. They can only downgrade
    #     VALIDATED to REVIEW, never force FAILED on their own.
    secondary_flags = []
    clashscore = result["MolProbity_Clashscore"]
    poor_rotamers = result["MolProbity_Poor_Rotamers_Pct"]
    rama_outlier = result["Rama_Outlier"]
    if clashscore > 20:
        secondary_flags.append(f"Clashscore {clashscore:.1f} > 20 (community default, not methodology-specified)")
    if poor_rotamers > 3.0:
        secondary_flags.append(f"Poor Rotamers {poor_rotamers:.1f}% > 3% (community default, not methodology-specified)")
    if rama_outlier > 0.5:
        secondary_flags.append(f"Rama Outliers {rama_outlier:.1f}% > 0.5% (community default, not methodology-specified)")
    result["Secondary_QC_Flags"] = " | ".join(secondary_flags) if secondary_flags else "None"

    # 3. Final verdict.
    #    NOTE: methodology grades MolProbity and Ramachandran independently
    #    but does not state how to combine them into one verdict. This
    #    requires BOTH to be satisfactory to call VALIDATED, treats a
    #    MARGINAL MolProbity score as REVIEW (not an automatic fail)
    #    rather than silently passing or failing it, and otherwise FAILS.
    if result["MP_Quality"] in ("EXCELLENT", "ACCEPTABLE") and result["Rama_Meets_Target"]:
        result["Overall_Status"] = "REVIEW" if secondary_flags else "VALIDATED"
    elif "MARGINAL" in result["MP_Quality"]:
        result["Overall_Status"] = "REVIEW"
    else:
        result["Overall_Status"] = "FAILED - REQUIRES REFINEMENT"

    # 4. Fold confidence -- reported, never folded into Overall_Status.
    result.update({
        "Mean_pLDDT": None,
        "Pct_Residues_pLDDT_ge70": None,
        "pTM": None,
        "Max_CrossModel_RMSD_Confident_Region": None,
        "Longest_Confident_Run": "",
        "Fold_Confidence": "NOT_ASSESSED",
        "Docking_Readiness": "NOT_ASSESSED",
        "Docking_Ready_Regions": "",
        "Docking_Readiness_Basis": "Step 2B fold-confidence report not found",
        "Fold_Confidence_Report": "",
    })
    if fold_confidence:
        for key in ("Mean_pLDDT", "Pct_Residues_pLDDT_ge70", "pTM",
                    "Max_CrossModel_RMSD_Confident_Region", "Longest_Confident_Run",
                    "Docking_Readiness", "Docking_Ready_Regions",
                    "Docking_Readiness_Basis", "Fold_Confidence_Report"):
            if fold_confidence.get(key) is not None:
                result[key] = fold_confidence[key]
        # Thresholds below are EXTERNAL CONVENTION, not methodology-specified
        # -- exactly like the secondary QC flags above, and labeled the same
        # way so the provenance is visible in the report itself.
        ptm = result["pTM"]
        pct = result["Pct_Residues_pLDDT_ge70"]
        if ptm is None or pct is None:
            result["Fold_Confidence"] = "UNKNOWN"
        elif ptm >= 0.5 and pct >= 70.0:
            result["Fold_Confidence"] = "GLOBAL FOLD DETERMINED (pTM >= 0.5, >=70% at pLDDT >= 70)"
        elif result["Docking_Readiness"] == "DOMAIN_ONLY":
            result["Fold_Confidence"] = "DOMAIN-LEVEL ONLY (global fold undetermined; pTM < 0.5)"
        else:
            result["Fold_Confidence"] = "UNDETERMINED (pTM < 0.5, no confident domain)"

    return result


def run_step2d_model_validation():
    start_time = time.time()
    project_root = _PROJECT_ROOT

    output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepD")
    os.makedirs(output_base, exist_ok=True)

    common.print_banner("PHASE 2 STEP D: STEREOCHEMICAL MODEL VALIDATION")
    print(f"[INFO] Resolved Project Root : {project_root}")
    print("[INFO] Methodology : MolProbity (phenix.molprobity local; web server as fallback)")
    print("[INFO] Targets     : MolProbity 0.5-1.5 excellent / <2 acceptable / >3 poor | Favored 95-98%")
    print("-" * 110)

    target = _resolve_validation_target(project_root)
    if target is None:
        return
    winner_name, safe_name, structure_path, pass_label = target

    if not os.path.isfile(structure_path):
        print(f"[ERROR] Structure file not found at {structure_path}.")
        print("[ERROR] Run Step 2C first so the Phenix-refined structure exists.")
        return

    print(f"[INFO] Validating Model : {os.path.relpath(structure_path, project_root)} ({pass_label})")
    print("-" * 110)

    # Manual results filename encodes which Phenix pass produced the
    # structure being validated -- if Step 2C later produces a pass2
    # structure where only pass1 MolProbity data was ever recorded, this
    # naturally requires a fresh MolProbity submission instead of silently
    # reporting pass1 numbers against the pass2 structure.
    # Prefer the LOCAL phenix.molprobity run -- same tool as the web
    # server, no manual submission, and no transcription step where a
    # value can be mistyped. Falls back to the manual web workflow only
    # if Phenix is unavailable or its run fails.
    manual_results_path = os.path.join(output_base, f"{safe_name}_{pass_label}_manual_results.json")
    manual_results = run_phenix_molprobity(structure_path)
    if manual_results is not None:
        print("[INFO] MolProbity source : phenix.molprobity (local run)")
        with open(manual_results_path, 'w') as f:
            json.dump(manual_results, f, indent=2)
        print(f"[INFO] Values recorded   : {os.path.relpath(manual_results_path, project_root)}")
    else:
        manual_results = load_manual_results(manual_results_path)
        if manual_results is None:
            run_manual_prepare(structure_path, manual_results_path)
            return
        print("[INFO] MolProbity source : manual web-server entry")

    fold_confidence = load_fold_confidence(project_root)
    results = evaluate_model_quality(manual_results, fold_confidence)

    print(f"{'VALIDATION METRIC':<30} | {'VALUE':<20} | {'ASSESSMENT'}")
    print("-" * 110)
    print(f"{'MolProbity Score':<30} | {results['MolProbity_Score']:<20.2f} | {results['MP_Quality']}")
    print(f"{'MolProbity Clashscore':<30} | {results['MolProbity_Clashscore']:<20.2f} | (contextual; see Secondary_QC_Flags)")
    print(f"{'MolProbity Poor Rotamers %':<30} | {results['MolProbity_Poor_Rotamers_Pct']:<19.2f}% | (contextual; see Secondary_QC_Flags)")
    print(f"{'Rama: Favored Region':<30} | {results['Rama_Favored']:<19.2f}% | {results['Rama_Label']}")
    print(f"{'Rama: Allowed Region':<30} | {results['Rama_Allowed']:<19.2f}% | (contextual)")
    print(f"{'Rama: Outlier Region':<30} | {results['Rama_Outlier']:<19.2f}% | (contextual; see Secondary_QC_Flags)")
    print(f"{'Secondary QC Flags':<30} | {results['Secondary_QC_Flags']}")
    print("-" * 110)
    # Labeled explicitly: this verdict is about STEREOCHEMISTRY only. It was
    # previously printed as a bare "FINAL DECISION", which reads as a verdict
    # on the model as a whole -- and a VALIDATED here was taken as clearance
    # to dock, which it never was.
    print(f"FINAL DECISION (stereochemistry) : [{results['Overall_Status']}]")
    print("-" * 110)
    print(f"{'Mean pLDDT (AlphaFold model)':<34} | {_fmt_conf(results['Mean_pLDDT'], '.1f')}")
    print(f"{'Residues at pLDDT >= 70':<34} | {_fmt_conf(results['Pct_Residues_pLDDT_ge70'], '.1f', '%')}")
    print(f"{'pTM (best ranked model)':<34} | {_fmt_conf(results['pTM'], '.2f')}")
    print(f"{'Longest confident run':<34} | {results['Longest_Confident_Run'] or 'n/a'}")
    print(f"{'Max cross-model RMSD (confident)':<34} | "
          f"{_fmt_conf(results['Max_CrossModel_RMSD_Confident_Region'], '.2f', ' A')}")
    print(f"{'Fold confidence':<34} | {results['Fold_Confidence']}")
    print("-" * 110)
    print(f"DOCKING READINESS (Phase III)    : [{results['Docking_Readiness']}]")
    if results["Docking_Ready_Regions"]:
        print(f"  Docking_Ready_Regions : {results['Docking_Ready_Regions']}")
    print(f"  Basis                 : {results['Docking_Readiness_Basis']}")
    if results["Docking_Readiness"] not in ("WHOLE_CONSTRUCT", "NOT_ASSESSED"):
        print("  NOTE: a VALIDATED stereochemistry verdict does NOT clear the whole")
        print("        construct for docking -- these are independent claims. pLDDT>=70 /")
        print("        pTM>=0.5 / RMSD<=2.0 A are external convention, see deviation #21.")
    print("-" * 110)

    ts = datetime.now().strftime("%Y%m%d_%H%M")
    report_path = os.path.join(output_base, f"Step2D_Validation_Report_{ts}.csv")

    csv_data = {
        "Variant": winner_name,
        "Structure_File": os.path.relpath(structure_path, output_base),
        **results,
        "Manual_Results_File": os.path.relpath(manual_results_path, output_base),
    }

    with open(report_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=csv_data.keys())
        w.writeheader()
        w.writerow(csv_data)

    total_time = common.format_time(time.time() - start_time)
    common.print_banner("PHASE 2 STEP D COMPLETE")
    print("[SUCCESS] Stereochemical quality verified using real MolProbity results.")
    print(f"[SUCCESS] Execution Time : {total_time}")
    print(f"[INFO] Report Saved      : {os.path.relpath(report_path, project_root)}")
    print("=" * 110 + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 2D: Stereochemical Model Validation")
    parser.add_argument("--manual-prepare", action="store_true",
                         help="Print MolProbity submission instructions and write the results template early")
    args = parser.parse_args()

    if args.manual_prepare:
        project_root = _PROJECT_ROOT
        target = _resolve_validation_target(project_root)
        if target is None:
            print("[ERROR] No filtered candidate found -- run Step 2A first.")
        else:
            winner_name, safe_name, structure_path, pass_label = target
            if not os.path.isfile(structure_path):
                print(f"[ERROR] Structure file not found at {structure_path}.")
                print("[ERROR] Run Step 2C first so the Phenix-refined structure exists.")
            else:
                output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepD")
                manual_results_path = os.path.join(output_base, f"{safe_name}_{pass_label}_manual_results.json")
                run_manual_prepare(structure_path, manual_results_path)
    else:
        run_step2d_model_validation()
