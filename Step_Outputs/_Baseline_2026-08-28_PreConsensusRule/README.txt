BASELINE SNAPSHOT -- STATE BEFORE THE ALLERGENICITY CONSENSUS-RULE CHANGE
=========================================================================
Taken 2026-08-28, immediately before Phase 1Eb's allergenicity combination
rule was changed from OR to CONSENSUS.

THIS IS NOT "SUPERSEDED" MATERIAL. Step_Outputs/_Superseded/ holds outputs
that were replaced by a newer run of the SAME method. This folder holds a
complete, internally consistent result set produced by a DIFFERENT method,
kept so the two rules can be compared directly and so every number already
written into the manuscript remains reproducible.

WHAT CHANGED AND WHY
--------------------
Phase 1Eb screens each peptide with AllerTOP v2.0 and AllergenFP v1.0. The
rule here was OR: a peptide was rejected if EITHER tool called it an
allergen. Measured on the real 210-peptide pool, the two tools DISAGREE on
96 of 210 (46%) and agree on ALLERGEN for only 14. So 96 peptides were being
deleted on a single-tool call from a pair that agrees less than half the
time.

The consequence was concrete: all 10 HIV_gp41 B-cell candidates were
rejected, 6 of them on split verdicts, leaving gp41 with ZERO B-cell
representation in the construct -- in an HIV vaccine, the antigen whose MPER
carries the best-characterised broadly neutralising antibody epitopes.

The new rule requires CONSENSUS (both tools) to reject. Survivors 98 -> 196.

CONSTRUCT IN THIS SNAPSHOT
--------------------------
  Vax_Final_0c83e66a  (484 aa)
  12 MHC-I, 11 MHC-II, 4 B-cell epitopes; B-cell coverage 3/7 antigens
  Phase 2A REVIEW (DeepSol 0.4131329 sole failing criterion)
  Phase 2C REVIEW (A3D patch 16 aa, CamSol 6 aa, hydrophobic 0.4036)
  Phase 2D VALIDATED stereochemistry (MolProbity 0.64, Rama 98.13% favored)
           DOMAIN_ONLY docking readiness, regions 1-45 (pTM 0.17)
  Phase 1H 9/27 experimentally confirmed, 12/27 overlapping, 6/27 novel
  Phase 1F cumulative coverage: MHC-I 99.9949%, MHC-II >99.99%

CONTENTS
--------
  Phase1/Phase1E/ toxicity, allergenicity (OR rule), self-homology
  Phase1/Phase1F/ population coverage, incl. the cumulative addendum
  Phase1/Phase1G/ the 484 aa construct + Boundary_Map
  Phase1/Phase1H/ IEDB corroboration (the _tool_runs/ IEDB cache is NOT
                  copied -- it is source data, not a result, and is reused
                  in place so a re-run does not re-download 19k records)
  Phase2/         2A-2D, including AlphaFold_Raw and all Supplementary_Archive
  Phase3/         prepared HADDOCK receptor + ligand inputs
  _project_docs/  Phase_I_II_Methods_Status.csv and
                  Methods_Deviations_RRL_Support.txt as they stood

RESTORING
---------
Paths mirror their original location under Step_Outputs/, so any entry can be
restored by copying it back along the same relative path. NOTE: Phase I outputs
were regrouped under Step_Outputs/Phase1/ on 2026-08-28 and this archive was
restructured to match -- see ARCHIVED_ITEMS.txt in the project root.

NOT READ BY THE PIPELINE
------------------------
Every input selector in Phases 1-3 uses os.listdir()/glob.glob() WITHOUT
recursion and filters on ".csv"/".fasta". This directory has no such
extension, so it is invisible to all of them -- the same property that makes
_Superseded/ safe.
