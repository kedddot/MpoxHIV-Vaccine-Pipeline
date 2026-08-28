# Phase 1Eb — Allergenicity Screening Methodology Note

**Applies to:** Section II.C.I.E of the proposal (Q/N Fraction, Charged Residue Count, AllerTop, AllergenFP).

## Decision rule
- QN_Ratio = (count(Q)+count(N))/length > 0.30 -> ALLERGEN (excluded). Per the
  proposal this flags the sequence unconditionally, not as one vote among three.
- Surface_Charge = count(D+E+H+K) > 4 -> DEPRIORITIZED (flagged, retained),
  matching the proposal's "deprioritizes" wording and the same soft-flag
  semantics already used for GRAVY elsewhere in this pipeline. Not a hard
  exclusion.
- AllerTOP positive OR AllergenFP positive -> ALLERGEN (excluded).
- A peptide is excluded if QN_Ratio>0.30 OR AllerTOP positive OR AllergenFP
  positive (independent signals). The proposal's explicit "2 of 3 predictors"
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
