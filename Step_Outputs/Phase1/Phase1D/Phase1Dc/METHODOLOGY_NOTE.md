# Phase 1Dc — Conservancy Methodology Note

**Applies to:** Section II.C.I.D of the proposal ("IEDB Conservancy Tool").

## What the manuscript currently says
The methods section states that epitope conservancy was assessed using the
IEDB Epitope Conservancy Analysis Tool.

## What this script actually does
IEDB does not expose a documented public REST API for the Conservancy
Analysis tool (unlike its MHC-I/MHC-II binding APIs), and standalone local
distributions of the tool have inconsistent, environment-dependent CLI
interfaces. Rather than depend on an unverified binary or silently fall back
to something undocumented, this script computes conservancy directly:

- For each target (e.g. `Mpox_L1R`, `HIV_gp120`), every surviving variant
  FASTA file from Phase 1C is loaded into a per-target pool of full-length
  sequences.
- For each candidate peptide, conservancy is calculated as the percentage of
  that target's variant pool containing an **exact, 100%-identity substring
  match** of the peptide: `Conservancy = (hits / total_variants) * 100`.

## Why this is equivalent, not just convenient
The IEDB Conservancy tool's identity threshold controls how much mismatch is
tolerated when searching for the epitope within each sequence. At a 100%
identity threshold specifically, no mismatches or gaps are permitted — which
is mathematically the same operation as an exact substring search. This
script always evaluates at 100% identity, so its output is equivalent to
running the IEDB tool at that same threshold on the same sequence set.

## What this does NOT reproduce
This substitution only holds at the 100% identity threshold. It does not
reproduce IEDB's behavior at lower identity thresholds (e.g. 80% or 90%
conservancy allowing partial mismatches), which the current script does not
attempt to compute. If partial-identity conservancy is ever required, this
approach would need to be replaced with an actual alignment-based method.

## Recommended manuscript update
Section II.C.I.D should be revised to state that conservancy was computed via
exact substring matching against each target's full variant pool at 100%
sequence identity, functionally equivalent to the IEDB Conservancy Tool at
that threshold, rather than citing the tool as if its CLI/API were invoked
directly.
