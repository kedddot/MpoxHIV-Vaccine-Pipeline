import os, sys, time, csv
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
# METHODOLOGY NOTE
# =============================================================================
# AllerTOP and AllergenFP have no public API and no local package -- verified
# (no CLI, no PyPI/conda package). This script uses the same manual
# web-submission + result-import pattern as Phase 1Ea's HemoPI/Tox-Prot BLAST
# checks: emit a query FASTA + CSV template, exit, and pick the filled-in
# results back up on the next run.
#
# DECISION RULE (a gap in the proposal, resolved here): the proposal defines
# an explicit "2 of {QN,AllerTop,AllergenFP} positive" consensus rule only
# for POST-CONSTRUCT allergenicity re-screening of the whole assembled
# vaccine (Section II.A.a) -- not for screening individual candidate
# peptides here in Section I.E. For individual peptides, the wording is
# "QN_Ratio > 0.30 flags the sequence as allergenic" (unconditional) and
# "Surface_Charge > 4 ... deprioritizes" (soft flag, not an exclusion --
# same semantics as the GRAVY >= 0.2 deprioritization used elsewhere in this
# pipeline).
#
# AllerTOP/AllergenFP now require CONSENSUS to exclude (deviation #24). An
# earlier version treated each as an independent exclusion signal (ANY
# positive -> ALLERGEN). Measured on the real pool the two predictors
# disagree on 41 of the 70 peptides where both can run (59%) and agree on
# ALLERGEN for only 14, so a single positive was excluding peptides on close
# to a coin flip -- and it cost HIV_gp41 every one of its B-cell candidates.
# Where only one predictor can run (AllergenFP cannot go below 16 aa, which
# in this pool is every MHC-I and MHC-II peptide) that predictor still
# decides alone, since consensus is undefined from a single opinion.
#
# QN_Ratio is unchanged and still excludes unconditionally: the proposal
# states it that way, and it is a biophysical criterion rather than one half
# of a disagreeing predictor pair. A still earlier version required 2-of-3
# positives INCLUDING QN_Ratio, which is the post-construct rule applied one
# stage too early and weakened a criterion the proposal states outright.
METHODOLOGY_NOTE_TEXT = """# Phase 1Eb — Allergenicity Screening Methodology Note

**Applies to:** Section II.C.I.E of the proposal (Q/N Fraction, Charged Residue Count, AllerTop, AllergenFP).

## Decision rule
- QN_Ratio = (count(Q)+count(N))/length > 0.30 -> ALLERGEN (excluded). Per the
  proposal this flags the sequence unconditionally, not as one vote among three.
- Surface_Charge = count(D+E+H+K) > 4 -> DEPRIORITIZED (flagged, retained),
  matching the proposal's "deprioritizes" wording and the same soft-flag
  semantics already used for GRAVY elsewhere in this pipeline. Not a hard
  exclusion.
- AllerTOP AND AllergenFP both positive -> ALLERGEN (excluded). Consensus is
  required because the two disagree on 59% of the peptides where both can run
  (41/70 measured), so either one alone is close to a coin flip. Where only
  one can run -- AllergenFP's ACC fingerprint fails below 16 aa, which covers
  every MHC-I and MHC-II peptide in this pool -- that one decides alone.
  Peptides retained on a split verdict are labelled as such in
  Exclusion_Reason and must not be reported as two clean negatives.
- A peptide is excluded if QN_Ratio>0.30 OR the applicable predictors agree on
  ALLERGEN. The proposal's explicit "2 of 3 predictors"
  consensus rule applies to POST-CONSTRUCT re-screening of the whole assembled
  vaccine (Section II.A.a), not to individual-peptide screening here; applying
  it at this stage would have silently weakened the unconditional QN_Ratio rule.

## Tooling
AllerTOP and AllergenFP have no public API and no local package. Both are
queried via a manual web-submission + result-import workflow, the same
pattern used for HemoPI/Tox-Prot BLAST in Phase 1Ea.

## AllergenFP length constraint (confirmed empirically)
AllergenFP's ACC (auto cross-covariance) fingerprint cannot be computed for
peptides shorter than 16 residues -- verified directly against the live
tool: AllerTOP accepted a 9-mer MHC-I candidate that AllergenFP rejected.
This affects every MHC-I 9/10-mer candidate (AllergenFP was only ever
usable for the 15/16-mer MHC-II and B-cell candidates). For peptides below
16 aa, `AllergenFP_Result` is pre-filled `N/A` and allergenicity is judged
on QN_Ratio and AllerTOP alone -- this is a resolved state (the tool
genuinely cannot score these), not a missing-data gap awaiting manual entry.
"""


def run_step1eb_allergenicity():
    start_time = time.time()

    # Reads Phase 1Ea's PRIORITIZED subset (top N per target x epitope class),
    # not its full Filtered output: AllerTOP/AllergenFP must be queried by hand
    # one peptide at a time, so the full pool is not manually screenable. The
    # prioritized subset caps that workload while guaranteeing every antigen and
    # all three epitope classes are represented. See Phase1Ea's prioritization
    # block for the selection rule.
    input_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1E", "Phase1Ea", "Prioritized")
    output_base = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1E", "Phase1Eb")
    raw_dir = os.path.join(output_base, "Raw")
    filt_dir = os.path.join(output_base, "Filtered")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(filt_dir, exist_ok=True)

    print("\n" + "="*80 + "\nPHASE 1Eb: ALLERGENICITY SCREENING\n" + "="*80)

    if not os.path.isdir(input_folder):
        print(f"[ERROR] Phase 1Ea Prioritized directory not found at: {input_folder}")
        print("[ERROR] Run Phase 1Ea first.")
        return
    csv_files = [f for f in os.listdir(input_folder) if f.endswith(".csv")]
    if not csv_files:
        print("[ERROR] No NON-TOXIC candidates found from Phase 1Ea.")
        return
    latest_csv = os.path.join(input_folder, sorted(csv_files)[-1])

    with open(latest_csv, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        original_fields = reader.fieldnames

    unique_peptides = sorted(set(row['Peptide'] for row in rows))
    print(f"[INFO] {len(rows)} candidate rows | {len(unique_peptides)} unique peptides")

    manual_results_file = os.path.join(output_base, "Manual_Allergenicity_Results.csv")
    fasta_export = os.path.join(output_base, "Allergenicity_Query.fasta")

    # Carry forward any peptide already answered in a PRIOR version of this
    # file (e.g. from an earlier, narrower TOP_N_PER_TARGET_PER_CLASS run).
    # Only genuinely new peptides -- ones with no prior valid answer -- are
    # written to the query FASTA / checklist, so widening the candidate pool
    # never asks for peptides that were already manually screened.
    manual_data = {}
    if os.path.exists(manual_results_file):
        with open(manual_results_file, 'r') as f:
            for row in csv.DictReader(f):
                manual_data[row['Peptide']] = row

    def _is_answered(pep):
        m = manual_data.get(pep)
        if not m:
            return False
        allertop_ok = m.get('AllerTOP_Result', '').strip().upper() in ("ALLERGEN", "NON-ALLERGEN")
        fp_val = m.get('AllergenFP_Result', '').strip().upper()
        allergenfp_ok = fp_val in ("ALLERGEN", "NON-ALLERGEN", "N/A", "NA", "TOO_SHORT")
        return allertop_ok and allergenfp_ok

    new_peptides = [p for p in unique_peptides if not _is_answered(p)]

    if new_peptides:
        # Rewrite the manual results file: keep every already-answered
        # peptide as-is, add blank rows only for the new ones.
        with open(fasta_export, "w") as f_fasta:
            for i, pep in enumerate(new_peptides):
                f_fasta.write(f">Pep_{i}\n{pep}\n")

        n_too_short = 0
        with open(manual_results_file, "w", newline='') as f_csv:
            writer = csv.writer(f_csv)
            writer.writerow(["Peptide", "AllerTOP_Result", "AllergenFP_Result"])
            for pep in unique_peptides:
                if pep in manual_data and _is_answered(pep):
                    m = manual_data[pep]
                    writer.writerow([pep, m.get('AllerTOP_Result', ''), m.get('AllergenFP_Result', '')])
                elif len(pep) < 16:
                    # AllergenFP's ACC fingerprint cannot be computed below 16
                    # residues (confirmed empirically -- AllerTOP accepts
                    # short peptides, AllergenFP rejects them).
                    writer.writerow([pep, "", "N/A"])
                    n_too_short += 1
                else:
                    writer.writerow([pep, "", ""])

        n_carried = len(unique_peptides) - len(new_peptides)
        print(f"\n[ACTION REQUIRED] {len(new_peptides)} NEW peptide(s) need screening "
              f"({n_carried} already-answered peptide(s) carried forward, no resubmission needed).")
        print(f"1. FASTA export ready at: {os.path.basename(fasta_export)} -- contains ONLY the {len(new_peptides)} new peptides")
        print(f"2. Query AllerTOP (https://ddg-pharmfac.net/AllerTOP/) for every new peptide, and AllergenFP (https://ddg-pharmfac.net/AllergenFP/) for new peptides >= 16 aa only")
        print(f"   ({n_too_short} of the new peptides are below 16 aa and pre-filled N/A for AllergenFP)")
        print(f"3. Fill in the blank rows (the new peptides) in: {os.path.basename(manual_results_file)} -- already-answered rows are untouched")
        print(f"4. Re-run this script to execute the filtering matrix.\n")
        sys.exit(0)

    fieldnames = original_fields + ["QN_Ratio", "Surface_Charge", "Surface_Charge_Deprioritized",
                                     "AllerTOP", "AllergenFP", "Allergenicity_Status", "Exclusion_Reason"]
    raw_data, filtered_data, review_data = [], [], []

    for i, row in enumerate(rows):
        pep = row['Peptide']
        qn_ratio = (pep.count('Q') + pep.count('N')) / len(pep)
        surface_charge = sum(pep.count(aa) for aa in "DEHK")

        m_res = manual_data.get(pep, {})
        unresolved_reasons = []

        allertop = m_res.get('AllerTOP_Result', '').strip().upper()
        if allertop not in ("ALLERGEN", "NON-ALLERGEN"):
            unresolved_reasons.append("AllerTOP_Result not filled in")

        allergenfp = m_res.get('AllergenFP_Result', '').strip().upper()
        # AllergenFP's ACC fingerprint genuinely cannot be computed below 16
        # residues (confirmed empirically: AllerTOP accepts short peptides,
        # AllergenFP rejects them). That's a permanent tool constraint for
        # every MHC-I 9/10-mer, not a temporary data-entry gap -- treating it
        # as "unresolved" would permanently strand these peptides in
        # Needs_Review since a valid AllergenFP result can never exist for
        # them. N/A is a resolved state: allergenicity is judged on QN_Ratio
        # and AllerTOP alone for these peptides.
        allergenfp_not_applicable = allergenfp in ("N/A", "NA", "TOO_SHORT", "N/A_TOO_SHORT")
        if not allergenfp_not_applicable and allergenfp not in ("ALLERGEN", "NON-ALLERGEN"):
            unresolved_reasons.append("AllergenFP_Result not filled in")

        fail_qn = qn_ratio > 0.30
        fail_allertop = (allertop == "ALLERGEN")
        fail_allergenfp = (allergenfp == "ALLERGEN")  # False when N/A, correctly contributing nothing
        is_deprioritized = surface_charge > 4

        # -------------------------------------------------------------------
        # ALLERGEN COMBINATION RULE: CONSENSUS, NOT "EITHER TOOL"  (dev. #24)
        #
        # This was OR -- either predictor alone excluded a peptide. Measured on
        # the real 210-peptide pool, AllerTOP and AllergenFP DISAGREE on 41 of
        # the 70 peptides where both can run (59%), and agree on ALLERGEN for
        # only 14. Using a single call from a pair that agrees less than half
        # the time as a HARD EXCLUSION deleted 41 peptides on what is close to
        # a coin flip.
        #
        # The consequence was concrete and costly: all 10 HIV_gp41 B-cell
        # candidates were excluded, 6 of them on split verdicts, leaving gp41 --
        # the HIV antigen whose MPER carries the best-characterised broadly
        # neutralising antibody epitopes -- with ZERO B-cell representation in
        # the construct.
        #
        # So exclusion now requires the two predictors to AGREE. This also
        # brings the step into line with the rest of the pipeline, where
        # stability (#2), population coverage (#9), GRAVY and Surface_Charge
        # were all already converted from hard filters to prioritise-not-exclude
        # signals; the allergenicity OR rule was the last over-strict filter.
        #
        # WHERE ONLY ONE TOOL CAN RUN, THAT TOOL STILL DECIDES ALONE. Consensus
        # is undefined with a single opinion, and silently downgrading it to
        # "cannot exclude" would discard AllerTOP's judgement on every short
        # peptide for no reason. AllergenFP's ACC fingerprint cannot be computed
        # below 16 aa, and in this pool EVERY MHC-I and MHC-II peptide is below
        # 16 aa -- so this rule change affects the B-cell pool ONLY (15 -> 56
        # survivors) and leaves both MHC pools untouched.
        #
        # QN_Ratio is unchanged and still excludes on its own: it is a
        # methodology-specified biophysical criterion (Sec. I.E), not one half
        # of a disagreeing predictor pair.
        # -------------------------------------------------------------------
        both_tools_ran = not allergenfp_not_applicable
        if both_tools_ran:
            fail_tools = fail_allertop and fail_allergenfp
            tools_disagree = fail_allertop != fail_allergenfp
        else:
            fail_tools = fail_allertop
            tools_disagree = False
        is_allergen = fail_qn or fail_tools
        is_unresolved = len(unresolved_reasons) > 0

        clean_row = {k: row[k] for k in original_fields}
        clean_row.update({
            "QN_Ratio": round(qn_ratio, 3), "Surface_Charge": surface_charge,
            "Surface_Charge_Deprioritized": is_deprioritized,
            "AllerTOP": allertop or "UNKNOWN", "AllergenFP": allergenfp or "UNKNOWN",
        })

        if is_unresolved:
            clean_row["Allergenicity_Status"] = "UNRESOLVED"
            clean_row["Exclusion_Reason"] = "; ".join(unresolved_reasons)
            review_data.append(clean_row)
        elif is_allergen:
            clean_row["Allergenicity_Status"] = "ALLERGEN"
            reasons = []
            if fail_qn: reasons.append("QN_Ratio>0.30")
            if fail_tools and both_tools_ran:
                reasons.append("AllerTOP=ALLERGEN AND AllergenFP=ALLERGEN (consensus)")
            elif fail_tools:
                reasons.append("AllerTOP=ALLERGEN (sole applicable tool; peptide <16 aa "
                               "so AllergenFP cannot run)")
            clean_row["Exclusion_Reason"] = "; ".join(reasons)
        else:
            clean_row["Allergenicity_Status"] = "NON-ALLERGEN"
            # A split verdict is retained, but it is NOT the same as two clean
            # negatives and must not be reported as one.
            clean_row["Exclusion_Reason"] = (
                "N/A -- RETAINED ON SPLIT VERDICT (AllerTOP and AllergenFP disagree; "
                "consensus required to exclude)" if tools_disagree else "N/A")
            filtered_data.append(clean_row)

        raw_data.append(clean_row)
        sys.stdout.write(f"\r[ PROCESS ] {i+1:03d}/{len(rows):03d} | Filtering Allergenicity Matrix")
        sys.stdout.flush()

    print(f"\n[INFO] Screening complete. NON-ALLERGEN survivors: {len(filtered_data)}/{len(rows)} | Needs review: {len(review_data)}")
    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    if raw_data:
        with open(os.path.join(raw_dir, f"Phase1Eb_Raw_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(raw_data)
    if filtered_data:
        with open(os.path.join(filt_dir, f"Phase1Eb_Filtered_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(filtered_data)
    if review_data:
        review_dir = os.path.join(output_base, "Needs_Review")
        os.makedirs(review_dir, exist_ok=True)
        with open(os.path.join(review_dir, f"Phase1Eb_NeedsReview_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(review_data)
        print(f"[ACTION REQUIRED] {len(review_data)} peptides need manual completion in {os.path.basename(manual_results_file)}")

    note_path = os.path.join(output_base, "METHODOLOGY_NOTE.md")
    with open(note_path, "w") as f:
        f.write(METHODOLOGY_NOTE_TEXT)
    print(f"[INFO] Methodology note written to {note_path}")

if __name__ == "__main__":
    run_step1eb_allergenicity()
