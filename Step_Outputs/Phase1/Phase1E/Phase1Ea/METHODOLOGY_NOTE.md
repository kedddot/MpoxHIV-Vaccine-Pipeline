# Phase 1Ea — Toxicity Screening Methodology Note

**Applies to:** Section II.C.I.E of the proposal (ToxinPred, HemoPI, BLASTP vs. UniProt Tox-Prot).

## What this script does
On a machine with the `phase2` conda environment's ToxinPred2, HemoPI2, and
local Tox-Prot BLAST database available, all three checks run automatically:
ToxinPred2 (model 2, hybrid RF+BLAST+MERCI), HemoPI2 (model 3, ESM2-t6), and
BLASTP against a local UniProt Tox-Prot database built with `makeblastdb`.

## ToxinPred composition-score calibration (IMPORTANT — cite this)
ToxinPred2's amino-acid-composition Random Forest is **not calibrated for
9–16 aa peptides**; it is trained to discriminate toxin *proteins*. Applied
to short epitopes it produces systematic false positives. Measured directly
against known-benign controls, at the tool's own default threshold (0.6):

| Control sequence | Identity | ML Score | Composition call |
|---|---|---|---|
| `AAYGPGPGKKAAY`   | Vaccine linkers only, no biological sequence | 0.621 | Toxin |
| `DAHKSEVAHRFKDLG` | Human serum albumin N-terminus | 0.688 | Toxin |
| `VKVGVNGFGRIGRLV` | Human GAPDH (housekeeping protein) | 0.670–0.721 | Toxin |
| `GIINTLQKYYCRVRG` | β-defensin-3 — **this study's own adjuvant** | 0.781 | Toxin |

That is a 100% false-positive rate on benign controls, including the
construct's own adjuvant. On the real candidate set the composition score
flagged 92.5% of peptides as toxins at threshold 0.5 (65.1% at 0.6).

By contrast, across all 1,025 real candidate peptides the two
**evidence-based** channels returned:
- MERCI toxin-motif hits: **0**
- ToxinPred internal BLAST hits to known toxins: **0**
- Independent Tox-Prot BLASTP 3-of-3 hits: **0**

Three independent homology/motif tests agree these viral epitopes bear no
resemblance to known toxins; only the miscalibrated composition heuristic
disagrees.

**Decision rule adopted:** a peptide is flagged TOXIC by ToxinPred only when
there is real toxin evidence (a BLAST homology hit or a MERCI motif hit).
The composition-only ML score is still computed and reported in every output
row (`ToxinPred_ML_Score`, with `COMPOSITION-FLAG` in the `ToxinPred` column
when ML >= 0.6), but is advisory and does not exclude a candidate. All other
toxicity criteria in Section II.C.I.E are applied unchanged as hard filters.

This is a documented, evidence-based deviation from a literal reading of the
methods section, in the same spirit as the Phase 1Dc conservancy note, and
should be disclosed in the manuscript rather than left implicit.

If those binaries/database are not found, this script falls back to a manual
workflow: it emits a query FASTA and a `Manual_Toxicity_Results.csv`
template, and a human runs the real web tools
(https://webs.iiitd.edu.in/raghava/toxinpred2/, HemoPI2, and a BLASTP search
against Tox-Prot) and pastes the results back in before re-running.

## Handling of incomplete data
A peptide is only classified NON-TOXIC if every required field is resolved.
Any peptide with a blank/unfilled manual field (fallback mode only) is
classified UNRESOLVED and routed to `Needs_Review`, never silently counted
as passing.

## BLASTP 3-of-3 rule
A peptide is TOXIC via Tox-Prot homology only if its best hit satisfies ALL
THREE: percent identity >= 80%, alignment coverage (alignment_length /
peptide_length) >= 80%, AND E-value <= 1e-5 -- not E-value alone.
