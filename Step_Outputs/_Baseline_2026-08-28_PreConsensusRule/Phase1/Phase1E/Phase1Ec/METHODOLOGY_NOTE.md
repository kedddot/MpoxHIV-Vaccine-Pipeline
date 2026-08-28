# Phase 1Ec — Human Self-Homology Screening Methodology Note

**Applies to:** Sec. I.E of the manuscript (BLASTP against the reviewed human
reference proteome, UniProt Homo sapiens Swiss-Prot subset).

## What this script does
Runs AFTER Phase 1Eb (allergenicity), matching the paper's stated order
("all candidate T-cell and B-cell epitopes surviving toxicity and
allergenicity screening were additionally screened..."). Every peptide in
1Eb's Filtered output is BLASTP'd against a local database built from
`reviewed:true AND organism_id:9606` (UniProt REST, Swiss-Prot subset only
-- not the much larger unreviewed TrEMBL set).

## Decision rule
A peptide is flagged "self-like" and EXCLUDED only if its best hit meets
ALL THREE: percent identity >= 70%, alignment coverage
(alignment_length / peptide_length) >= 70%, E-value <= 1e-5.

Any hit below that bar is RETAINED but flagged `Self_Homology_Status =
PARTIAL` for manual review, per the paper's own text: "Epitopes with
partial homology below this threshold were retained but flagged for
manual review, given that even sub-threshold similarity to self-proteins
can, in some cases, still influence immunogenicity through partial
tolerance effects."

## Database provenance
`human_swissprot_db/human_swissprot.fasta`, fetched from
`rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&format=fasta`
-- 20,431 reviewed human sequences, indexed with `makeblastdb -dbtype prot`.
