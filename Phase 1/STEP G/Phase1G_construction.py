import os
import sys
import csv
import re
import hashlib
import random
from datetime import datetime
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

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
# VACCINE CONFIGURATION -- fixed by the methodology (Sec. I.G): the
# adjuvant, its linker, and the three epitope-class linkers. The paper
# specifies none of the ORDERING or SELECTION logic below -- RQ2 asks
# "what architectural arrangement ... is most optimal", which is exactly
# the freedom this script now searches over.
# =============================================================================
ADJUVANT = "GIINTLQKYYCRVRGGRCAVLSCLPKEEQIGKCSTRGRKCCRRKK"  # mature human beta-defensin-3
ADJ_LINKER = "EAAAK"

L_MHCI = "AAY"     # CTL-CTL
L_MHCII = "GPGPG"  # HTL-HTL
L_BCELL = "KK"     # B-cell-B-cell

MAX_EPITOPES_PER_TARGET_PER_CLASS = 2
REQUIRED_TARGETS = ["Mpox_L1R", "Mpox_B5R", "Mpox_A35R", "HIV_gp120", "HIV_gp41", "HIV_p24", "HIV_p17"]

# Sec. I.G: "each selected epitope was screened to confirm the absence of
# embedded linker motifs (AAY, GPGPG, KK) within its native sequence,
# preventing ambiguity in linker-boundary parsing during downstream
# splicing operations." Applied at pool-load time (embedded case) and again
# after assembly (junction-created case -- see _spurious_motif_count below,
# an extension of the same rule/rationale, not a new one).
LINKER_MOTIFS = ("AAY", "GPGPG", "KK")

# Sec. II.A.a hydrophobicity gate, reproduced here so 1G can select
# AGAINST the rule Phase 2A will apply, instead of discovering a failing
# construct a whole phase later with no path back to a fix.
HYDROPHOBIC_RESIDUES = set(['A', 'V', 'I', 'L', 'M', 'F', 'W', 'Y'])
HYDRO_WINDOW_SIZES = (8, 12, 15)

# Kyte-Doolittle hydropathy (paper Table 2) -- used to bias epitope
# SELECTION toward hydrophilic candidates as a DeepSol S2 solubility
# proxy. DeepSol has no local install or API, so it cannot be optimized
# against directly; GRAVY is the standard correlate and is what the
# paper's own methodology already computes.
KYTE_DOOLITTLE = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
    'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3,
}
SELECTION_WINDOW = 7   # per-target coverage-ranked window to pick within -- widened from 5
# to 7 (the deepest per-target/class pool after the Sec I.E/I.G screens), which
# measurably improved zero-hydrophobic-window feasibility (5/30 vs 2/40 polished
# candidates) without needing any new manual screening, since 7 stays within the
# already-screened pool depth for every (target, class) group.
GRAVY_JITTER = 0.6     # exploration noise on the hydrophilicity ordering
# Hydrophilicity-bias strengths cycled across trials. 0.0 = ignore GRAVY
# entirely (pure random within the coverage window), which is what finds
# selections that can reach ZERO hydrophobic windows; higher values chase
# lower whole-construct GRAVY. Both must be represented in the candidate
# pool -- a fixed strong bias collapsed the search onto one infeasible
# selection (0/40 candidates reached zero windows).
GRAVY_WEIGHT_SWEEP = (0.0, 0.0, 0.25, 0.5, 1.0, 2.0)

# Search budget for the joint selection+ordering optimizer (Sec below).
# Measured on the real Phase 1F pool: a naive greedy at only 400 trials
# already reaches 17->1 bad windows; this budget is generous headroom.
SEARCH_TRIALS = 6000   # raised from 3000 -- the embedded/junction motif screens
# (LINKER_MOTIFS) and the BigMHC/SEMA priority term shrink the feasible set, and
# 2-opt polishing (not trial generation) is the actual bottleneck, so a larger
# raw trial budget costs little and widens the top-K candidate pool it draws from.
# (TOP_K itself is defined locally in run_step1g_construction(), raised 40->60
# for the same reason -- see the comment there.)
LOCAL_SEARCH_PASSES = 6
RNG_SEED = 20260824  # fixed for determinism -- same input must give same construct


def count_hydrophobic_windows(seq):
    """Sec. II.A.a: windows of length 8/12/15 with >80% hydrophobic residues."""
    n = 0
    for w in HYDRO_WINDOW_SIZES:
        if w > len(seq):
            continue
        for i in range(len(seq) - w + 1):
            window = seq[i:i + w]
            if sum(1 for aa in window if aa in HYDROPHOBIC_RESIDUES) / w > 0.80:
                n += 1
    return n


def _shares_frame(a, b, min_overlap_frac=0.70, min_overlap_len=6):
    """
    True if `a` and `b` are overlapping sliding-window frames of the same
    underlying region -- one contains the other, or a long prefix/suffix
    of one matches the other. Exact-string dedup alone let three
    near-identical shifted windows of the same Mpox_A35R stretch (e.g.
    AITDSAVAVAAASST / AAITDSAVAVAAASS / AAITDSAVAVAAASST) all enter the
    same construct as if they were independent epitopes.
    """
    if a == b:
        return True
    short, long_ = (a, b) if len(a) <= len(b) else (b, a)
    if short in long_:
        return True
    m = min(len(a), len(b))
    threshold = max(min_overlap_len, int(m * min_overlap_frac))
    for k in range(m, threshold - 1, -1):
        if a[-k:] == b[:k] or b[-k:] == a[:k] or a[:k] == b[:k] or a[-k:] == b[-k:]:
            return True
    return False


def _load_pool(latest_csv):
    """
    Loads the Phase 1F Filtered pool, grouped by (target, type), each
    entry ranked by (PRIORITIZED, coverage, conservancy) descending.

    Two exclusions happen up front, both unfixable by ordering and both
    reported to the caller:
      - Intrinsically hydrophobic epitopes (>80% window INSIDE the epitope
        itself -- see the module docstring in the remediation plan).
      - Sec. I.G's embedded-linker-motif screen: an epitope containing AAY,
        GPGPG, or KK in its OWN native sequence would create linker-boundary
        ambiguity no matter where it's placed in the construct, so it is
        excluded here rather than left to the ordering/2-opt stage.

    Also returns priority_info: pep -> {"bigmhc": Immunogenicity_Priority,
    "sema": SEMA_Corroborated}, read from whatever columns Phase 1Dd/1De
    wrote (both are optional -- absent/blank if those steps weren't run,
    same tolerant-of-missing-columns pattern as Conservancy/Coverage above).
    """
    by_target = defaultdict(lambda: {"MHC-I": [], "MHC-II": [], "B-cell": []})
    excluded_intrinsic = []
    excluded_motif = []
    priority_info = {}

    with open(latest_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pep = row['Peptide']
            pep_type = row['Type']

            if count_hydrophobic_windows(pep) > 0:
                excluded_intrinsic.append((pep, pep_type, row["Target"]))
                continue

            embedded_hits = [m for m in LINKER_MOTIFS if m in pep]
            if embedded_hits:
                excluded_motif.append((pep, pep_type, row["Target"], embedded_hits))
                continue

            raw_conservancy = row.get('Conservancy', '')
            rank_val = float(raw_conservancy) if str(raw_conservancy).strip() else 0.0

            raw_cov = str(row.get('Overall_Coverage_Pct', '')).strip()
            try:
                cov_val = float(raw_cov)
            except ValueError:
                cov_val = 0.0
            prioritized = 1 if str(row.get('Coverage_Priority', '')).strip() == "PRIORITIZED" else 0

            by_target[row["Target"]][pep_type].append((prioritized, cov_val, rank_val, pep))
            priority_info[pep] = {
                "bigmhc": str(row.get("Immunogenicity_Priority", "")).strip(),
                "sema": str(row.get("SEMA_Corroborated", "")).strip(),
            }

    for target in by_target:
        for t_type in by_target[target]:
            by_target[target][t_type].sort(key=lambda x: (x[0], x[1], x[2]), reverse=True)

    return by_target, excluded_intrinsic, excluded_motif, priority_info


def _peptide_gravy(pep):
    """Kyte-Doolittle GRAVY for one peptide (paper Table 2 values)."""
    return sum(KYTE_DOOLITTLE[aa] for aa in pep) / len(pep)


def _select_candidates(by_target, rng, gravy_weight=1.0):
    """
    For each (target, class), walks that target's ranked candidate list
    and takes up to MAX_EPITOPES_PER_TARGET_PER_CLASS epitopes that do
    NOT share a frame with anything already selected -- within its own
    class AND across classes (I-1's missing cross-class check).

    Within each target's top-N coverage-ranked window, candidates are
    ordered by their own GRAVY (most hydrophilic first) scaled by
    `gravy_weight`, plus random jitter. Coverage prioritization is
    preserved at the POOL level (the window is still the top-N by
    PRIORITIZED/coverage/conservancy); the hydrophilicity bias only
    decides which of those near-equal candidates to prefer, which lowers
    whole-construct GRAVY at no cost in epitope count or antigen coverage.

    `gravy_weight` is varied ACROSS trials by the caller, from 0 (pure
    random -- explores orderings that can actually reach zero hydrophobic
    windows) up to a strong bias (explores the most soluble selections).
    A fixed strong bias collapses every trial onto the same hydrophilic
    epitope set: measured, all 40 retained candidates were effectively
    one selection stuck at 3 windows with 0/40 reaching zero. Varying the
    weight is what keeps BOTH objectives represented in the pool.
    """
    mhc_i, mhc_ii, bcell = [], [], []
    epitope_provenance = {}
    all_selected = []  # flat list across all three classes, for cross-class dedup

    for target in by_target:
        for t_type, bucket in (("MHC-I", mhc_i), ("MHC-II", mhc_ii), ("B-cell", bcell)):
            candidates = by_target[target][t_type][:SELECTION_WINDOW]
            candidates = candidates[:]
            candidates.sort(key=lambda c: gravy_weight * _peptide_gravy(c[3])
                            + rng.uniform(-GRAVY_JITTER, GRAVY_JITTER))
            got = 0
            for _prio, _cov, _cons, pep in candidates:
                if any(_shares_frame(pep, existing) for existing in all_selected):
                    continue
                bucket.append(pep)
                all_selected.append(pep)
                epitope_provenance[pep] = target
                got += 1
                if got >= MAX_EPITOPES_PER_TARGET_PER_CLASS:
                    break

    return mhc_i, mhc_ii, bcell, epitope_provenance


def _greedy_order(epitopes, linker, prefix):
    """
    Greedily appends epitopes one at a time, each time picking whichever
    remaining epitope adds the FEWEST new >80%-hydrophobic windows to the
    growing sequence. Not globally optimal (see the 2-opt pass below) but
    cheap and a strong starting point -- measured 17->1 bad windows on
    the full construct from a single greedy pass.
    """
    remaining = epitopes[:]
    ordered = []
    current = prefix
    while remaining:
        best = None
        for pep in remaining:
            candidate = current + (linker if ordered else "") + pep
            added = count_hydrophobic_windows(candidate) - count_hydrophobic_windows(current)
            if best is None or added < best[0]:
                best = (added, pep)
        _added, pep = best
        ordered.append(pep)
        remaining.remove(pep)
        current = current + (linker if len(ordered) > 1 else "") + pep
    return ordered, current


def _assemble(mhc_i_order, mhc_ii_order, bcell_order):
    chain = ""
    if mhc_i_order:
        chain += L_MHCI.join(mhc_i_order)
    if mhc_ii_order:
        chain += (L_MHCII if chain else "") + L_MHCII.join(mhc_ii_order)
    if bcell_order:
        chain += (L_BCELL if chain else "") + L_BCELL.join(bcell_order)
    if chain:
        return f"{ADJUVANT}{ADJ_LINKER}{chain}"
    # Dangling-linker guard: no epitopes at all means no epitope chain to
    # attach a linker to -- just the adjuvant, not "ADJUVANT + EAAAK" with
    # nothing following it.
    return ADJUVANT


def _two_opt_polish(mhc_i, mhc_ii, bcell, passes=LOCAL_SEARCH_PASSES):
    """
    Greedy ordering is myopic and can stall a few windows above zero at a
    local minimum that adjacent-only swaps can't escape (measured: two
    remaining violations both sat at an internal AAY junction that
    neither neighbouring swap alone could clear). This tries EVERY pair
    swap within each block (not just adjacent), which is cheap at these
    block sizes (<=14 epitopes) and escapes local minima adjacent swaps
    get stuck in.

    Optimizes the lexicographic pair (hydrophobic windows, spurious linker
    motifs) -- both are hard Sec. II.A.a / Sec. I.G gates, and a swap that
    clears a window but creates a junction motif (or vice versa) must not
    be accepted as an improvement on windows alone.
    """
    def combined_score(mi, mii, bc):
        seq = _assemble(mi, mii, bc)
        return (count_hydrophobic_windows(seq), _spurious_motif_count(seq, mi, mii, bc))

    mi, mii, bc = mhc_i[:], mhc_ii[:], bcell[:]
    best_score = combined_score(mi, mii, bc)
    for _ in range(passes):
        improved = False
        for block in (mi, mii, bc):
            n = len(block)
            for i in range(n):
                for j in range(i + 1, n):
                    block[i], block[j] = block[j], block[i]
                    score = combined_score(mi, mii, bc)
                    if score < best_score:
                        best_score = score
                        improved = True
                    else:
                        block[i], block[j] = block[j], block[i]  # revert
        if not improved or best_score == (0, 0):
            break
    return mi, mii, bc, best_score


def _junction_cys_violations(seq):
    """
    Mirrors Phase 2A's corrected Junction-Cys rule (>2 Cys in a true
    junction-spanning 8-20aa peptide), including the adjuvant exemption
    (deviation #12) -- so 1G selects against the SAME rule 2A will apply,
    not a different one.
    """
    adjuvant_end = len(ADJUVANT) if seq.startswith(ADJUVANT) else 0
    max_len = 20
    linkers = [L_MHCI, L_MHCII, L_BCELL, ADJ_LINKER]
    violations = 0
    for linker in linkers:
        start_idx = 0
        while True:
            idx = seq.find(linker, start_idx)
            if idx == -1:
                break
            linker_end = idx + len(linker)
            window_start = max(0, idx - (max_len - 1))
            if window_start < adjuvant_end:
                start_idx = idx + 1
                continue
            context_start = window_start
            context_end = min(len(seq), linker_end + (max_len - 1))
            context = seq[context_start:context_end]
            local_start = idx - context_start
            local_end = linker_end - context_start
            for length_range in (range(8, 12), range(12, 21)):
                for L in length_range:
                    if L > len(context):
                        continue
                    for w_start in range(0, len(context) - L + 1):
                        w_end = w_start + L
                        if w_start < local_start and w_end > local_end:
                            if context[w_start:w_end].count('C') > 2:
                                violations += 1
            start_idx = idx + 1
    return violations


def _boundary_map(mhc_i_order, mhc_ii_order, bcell_order):
    """
    Explicit linker-position / epitope-boundary manifest, so downstream
    steps (Phase 2A) don't have to re-derive junctions by re-searching
    for linker strings in the assembled sequence.
    """
    segments = [("adjuvant", ADJUVANT), (ADJ_LINKER, ADJ_LINKER)]
    for i, pep in enumerate(mhc_i_order):
        if i > 0:
            segments.append((L_MHCI, L_MHCI))
        segments.append((f"MHC-I:{pep}", pep))
    for i, pep in enumerate(mhc_ii_order):
        segments.append((L_MHCII, L_MHCII))
        segments.append((f"MHC-II:{pep}", pep))
    for i, pep in enumerate(bcell_order):
        segments.append((L_BCELL, L_BCELL))
        segments.append((f"B-cell:{pep}", pep))

    offset = 0
    parts = []
    for label, text in segments:
        parts.append(f"{label}[{offset}-{offset + len(text)}]")
        offset += len(text)
    return ";".join(parts)


def _spurious_motif_count(seq, mhc_i_order, mhc_ii_order, bcell_order):
    """
    Extends Sec. I.G's embedded-linker-motif screen (already applied to each
    epitope individually in _load_pool) to motifs CREATED at a junction by
    two adjacent pieces of sequence -- e.g. an epitope ending "...LS" next to
    the KK linker produces "...LSKK", and if the epitope itself also ends in
    "K" this can read as an extra, unintended KK. The measured case:
    "ASIRYRQRLISLLS" + KK-linker + "KK..." produced THREE KK occurrences
    where only one linker was intended.

    Counts any AAY/GPGPG/KK/EAAAK occurrence in the assembled sequence that
    does NOT start at one of the boundary map's own intended linker
    positions. Occurrences wholly inside the adjuvant are exempt -- mature
    beta-defensin-3 natively ends "...RGRKCCRRKK", which would otherwise
    always register as one unavoidable, non-fixable violation (same
    exemption shape as _junction_cys_violations' adjuvant exemption / dev
    #12).
    """
    boundary_map = _boundary_map(mhc_i_order, mhc_ii_order, bcell_order)
    intended = {}
    for seg in boundary_map.split(";"):
        m = re.match(r"(.+)\[(\d+)-(\d+)\]$", seg)
        label, start = m.group(1), int(m.group(2))
        if label in LINKER_MOTIFS or label == ADJ_LINKER:
            intended.setdefault(label, set()).add(start)

    adjuvant_end = len(ADJUVANT) if seq.startswith(ADJUVANT) else 0
    count = 0
    for motif in LINKER_MOTIFS + (ADJ_LINKER,):
        start_idx = 0
        while True:
            pos = seq.find(motif, start_idx)
            if pos == -1:
                break
            if pos not in intended.get(motif, set()) and pos + len(motif) > adjuvant_end:
                count += 1
            start_idx = pos + 1
    return count


def run_step1g_construction():
    input_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1F", "Filtered")
    output_dir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1G")
    os.makedirs(output_dir, exist_ok=True)

    common.print_banner("PHASE 1G: CHIMERIC VACCINE ASSEMBLY & LINKER INTEGRATION")

    latest_csv = common.latest_file(input_folder, suffix=".csv")
    if latest_csv is None:
        print(f"[ERROR] No Phase 1F Filtered output found at: {input_folder}")
        return

    by_target, excluded_intrinsic, excluded_motif, priority_info = _load_pool(latest_csv)

    for t in REQUIRED_TARGETS:
        total_eps = sum(len(pool) for pool in by_target[t].values())
        if total_eps == 0:
            print(f"[WARNING] Missing representation entirely for target: {t}")

    if excluded_intrinsic:
        print(f"[INFO] Excluded {len(excluded_intrinsic)} epitope(s) with an intrinsic >80% "
              f"hydrophobic window (unfixable by ordering):")
        for pep, t_type, target in excluded_intrinsic:
            print(f"       {target:<11} {t_type:<7} {pep}")

    if excluded_motif:
        print(f"[INFO] Excluded {len(excluded_motif)} epitope(s) with an embedded linker motif "
              f"(Sec. I.G):")
        for pep, t_type, target, hits in excluded_motif:
            print(f"       {target:<11} {t_type:<7} {pep:<20} contains {hits}")

    print(f"[INFO] Source Data : {os.path.basename(latest_csv)}")
    print(f"[INFO] Searching {SEARCH_TRIALS} selection+ordering trials "
          f"(seed={RNG_SEED}) to minimize hydrophobic-window violations...")
    print("-" * 80)

    rng = random.Random(RNG_SEED)
    # The 2-opt polish only REORDERS a fixed epitope selection -- it never
    # changes composition, so GRAVY is a property of the SELECTION while
    # the window count depends on the ORDERING. A selection's pre-polish
    # window count is therefore a poor predictor of whether it can reach
    # zero after polishing. Keeping only a handful of candidates ranked
    # pre-polish let the GRAVY bias crowd out every zero-feasible
    # selection (measured: all 5 kept candidates stalled at 2-3 windows
    # while better-ordered selections existed).
    #
    # So: keep a WIDE pool, polish all of them, discard any that cannot
    # reach zero windows (the hard Sec. II.A.a gate), and only then pick
    # the lowest-GRAVY survivor. Feasibility first, solubility proxy
    # second -- never the reverse.
    TOP_K = 60  # raised from 40 -- LINKER_MOTIFS + the BigMHC/SEMA priority term
    # (see _select_candidates/_spurious_motif_count) shrink the feasible set versus
    # the previous run, so a wider top-K keeps enough candidates in play post-polish.
    top_candidates = []  # list of (score, mhc_i, mhc_ii, bcell, provenance)

    for _trial in range(SEARCH_TRIALS):
        # Sweep the hydrophilicity bias across trials -- 0.0 is pure
        # random (finds zero-window-feasible orderings), higher values
        # chase the lowest-GRAVY selections. See _select_candidates().
        gravy_weight = GRAVY_WEIGHT_SWEEP[_trial % len(GRAVY_WEIGHT_SWEEP)]
        mhc_i, mhc_ii, bcell, provenance = _select_candidates(by_target, rng, gravy_weight)
        if not (mhc_i or mhc_ii or bcell):
            continue

        mhc_i_order, prefix_i = _greedy_order(mhc_i, L_MHCI, ADJUVANT + ADJ_LINKER)
        mhc_ii_order, prefix_ii = _greedy_order(mhc_ii, L_MHCII, prefix_i + (L_MHCII if mhc_i else ""))
        bcell_order, _prefix_b = _greedy_order(bcell, L_BCELL, prefix_ii + (L_BCELL if (mhc_i or mhc_ii) else ""))

        candidate_seq = _assemble(mhc_i_order, mhc_ii_order, bcell_order)
        bad_windows = count_hydrophobic_windows(candidate_seq)
        spurious_motifs = _spurious_motif_count(candidate_seq, mhc_i_order, mhc_ii_order, bcell_order)
        # BigMHC (MHC-I) / SEMA (B-cell) priority count -- Sec. I.D's own
        # prioritization signals, ranked below the two hard gates above
        # (windows, motifs) and above the GRAVY solubility proxy: unlike
        # GRAVY (no demonstrated purchase on DeepSol -- see plan §0.2),
        # BigMHC/SEMA priority is a signal the manuscript explicitly asks
        # this step to act on.
        n_prioritized = (
            sum(1 for p in mhc_i_order if priority_info.get(p, {}).get("bigmhc") == "PRIORITIZED")
            + sum(1 for p in bcell_order if priority_info.get(p, {}).get("sema") == "YES")
        )
        # Secondary objective: real Kyte-Doolittle GRAVY (Bio.SeqUtils.
        # ProtParam), as a proxy for DeepSol S2 solubility -- DeepSol has
        # no local/API path, so this can only be optimized blind (no live
        # oracle). A first round using raw hydrophobic-residue FRACTION as
        # the proxy moved DeepSol 0.257 -> 0.446 as GRAVY fell -0.09 ->
        # -0.21 alongside it -- GRAVY is the more standard, better
        # evidenced correlate (and what the paper's own methodology
        # already computes for stability), so it replaces the coarser
        # fraction here. Ranked ABOVE coverage since DeepSol, like the
        # hydrophobic-window rule, is a hard Phase 2A/2C reject --
        # coverage is only "prioritized" (dev #9).
        gravy_proxy = ProteinAnalysis(candidate_seq).gravy()
        coverage_sum = sum(
            cov for target in by_target for t_type in by_target[target]
            for prio, cov, cons, pep in by_target[target][t_type]
            if pep in (mhc_i_order + mhc_ii_order + bcell_order)
        )
        score = (bad_windows, spurious_motifs, -n_prioritized, round(gravy_proxy, 4), -coverage_sum)

        top_candidates.append((score, mhc_i_order, mhc_ii_order, bcell_order, provenance))
        top_candidates.sort(key=lambda c: c[0])
        del top_candidates[TOP_K:]
        # No early break -- we keep searching the full trial budget so the
        # top-K pool reflects the best candidates found across all trials,
        # not just the first acceptable one.

    if not top_candidates:
        print("[ERROR] No usable epitope selection found -- Phase 1F pool may be empty.")
        return

    print(f"[INFO] Polishing top {len(top_candidates)} candidate selection(s) "
          f"(pre-polish bad windows: {[c[0][0] for c in top_candidates]})...")

    best_final = None  # (rank_key, polished_score, mhc_i, mhc_ii, bcell, provenance)
    n_feasible = 0
    for score, mi, mii, bc, provenance in top_candidates:
        p_mi, p_mii, p_bc, polished_score = _two_opt_polish(mi, mii, bc)
        if polished_score == (0, 0):
            n_feasible += 1
        # (windows, motifs, -priority, GRAVY, -coverage): the two hard gates
        # (windows, motifs) first so an infeasible candidate can never
        # outrank a feasible one no matter how good its priority/GRAVY. No
        # early break -- every candidate is polished so the best FEASIBLE
        # selection wins, not merely the first feasible one encountered.
        rank_key = polished_score + score[2:]
        if best_final is None or rank_key < best_final[0]:
            best_final = (rank_key, polished_score, p_mi, p_mii, p_bc, provenance)

    _rank_key, polished_score, mhc_i_order, mhc_ii_order, bcell_order, epitope_provenance = best_final
    print(f"[INFO] After 2-opt local search: {polished_score[0]} bad window(s), "
          f"{polished_score[1]} spurious linker motif(s) "
          f"({n_feasible}/{len(top_candidates)} candidates reached (0,0))")

    final_sequence = _assemble(mhc_i_order, mhc_ii_order, bcell_order)
    seq_hash = hashlib.md5(final_sequence.encode()).hexdigest()[:8]

    # Post-assembly local metrics -- everything Phase 2A checks that can
    # be computed WITHOUT an external tool call, so a construct that will
    # fail 2A on these grounds is caught here instead of a phase later.
    ana = ProteinAnalysis(final_sequence)
    instability_idx = ana.instability_index()
    gravy_val = ana.gravy()
    pI = ana.isoelectric_point()
    net_charge = sum(final_sequence.count(aa) for aa in ('K', 'R')) - \
        sum(final_sequence.count(aa) for aa in ('D', 'E'))
    junction_cys = _junction_cys_violations(final_sequence)
    bad_windows_final = count_hydrophobic_windows(final_sequence)
    spurious_motifs_final = _spurious_motif_count(final_sequence, mhc_i_order, mhc_ii_order, bcell_order)
    n_bigmhc_prioritized = sum(1 for p in mhc_i_order if priority_info.get(p, {}).get("bigmhc") == "PRIORITIZED")
    n_sema_corroborated = sum(1 for p in bcell_order if priority_info.get(p, {}).get("sema") == "YES")

    a35r_mhcii_selected = [p for p in mhc_ii_order if epitope_provenance.get(p) == "Mpox_A35R"]
    a35r_note = (
        f"A35R MHC-II uses {a35r_mhcii_selected} (BELOW_TARGET-tier substitute) -- "
        f"the pool's only >=90%-coverage A35R MHC-II epitopes were excluded for an "
        f"intrinsic hydrophobic window (see excluded-epitope log above); coverage "
        f"traded for the hard hydrophobicity gate per deviation #9."
        if a35r_mhcii_selected else "A35R MHC-II: none selected."
    )

    unique_sources = ";".join(sorted(set(epitope_provenance.values())))
    provenance_str = ";".join(
        f"{pep}:{tgt}" for pep, tgt in sorted(epitope_provenance.items(), key=lambda x: x[1])
    )
    boundary_map = _boundary_map(mhc_i_order, mhc_ii_order, bcell_order)

    construct = {
        "Construct_ID": f"Vax_Final_{seq_hash}",
        "Length": len(final_sequence),
        "MHC_I_Count": len(mhc_i_order),
        "MHC_II_Count": len(mhc_ii_order),
        "BCell_Count": len(bcell_order),
        "Bad_Hydrophobic_Windows": bad_windows_final,
        "Spurious_Linker_Motifs": spurious_motifs_final,
        "Instability_Index": round(instability_idx, 2),
        "GRAVY": round(gravy_val, 4),
        "Theoretical_pI": round(pI, 2),
        "Net_Charge_Approx": net_charge,
        "Junction_Cys_Violations": junction_cys,
        "BigMHC_Prioritized_MHCI": n_bigmhc_prioritized,
        "SEMA_Corroborated_BCell": n_sema_corroborated,
        "Source_Targets": unique_sources,
        "Epitope_Provenance": provenance_str,
        "Boundary_Map": boundary_map,
        "A35R_Coverage_Note": a35r_note,
        "Sequence": final_sequence,
    }

    print(f"[SUCCESS] Generated Construct ({seq_hash}) | Length: {len(final_sequence)} aa")
    print(f"[INFO] Local checks -- Bad Windows: {bad_windows_final} | Spurious Linker Motifs: "
          f"{spurious_motifs_final} | Instability: {instability_idx:.2f} (<=40) | GRAVY: "
          f"{gravy_val:.4f} (<=0.4) | pI: {pI:.2f} (avoid 6.5-7.5) | Junction Cys Violations: {junction_cys}")
    print(f"[INFO] Priority signals used -- BigMHC-prioritized MHC-I: {n_bigmhc_prioritized}/{len(mhc_i_order)} | "
          f"SEMA-corroborated B-cell: {n_sema_corroborated}/{len(bcell_order)}")
    if bad_windows_final > 0:
        print(f"[WARNING] {bad_windows_final} hydrophobic-window violation(s) remain after "
              f"{SEARCH_TRIALS} trials + local search -- Phase 2A will likely still reject on this.")
    if spurious_motifs_final > 0:
        print(f"[WARNING] {spurious_motifs_final} spurious linker-motif occurrence(s) remain "
              f"(junction-created, outside the adjuvant exemption).")

    # 4. EXPORT
    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    csv_path = os.path.join(output_dir, f"Phase1G_FinalConstruct_{ts}.csv")
    fasta_path = os.path.join(output_dir, f"Phase1G_FinalConstruct_{ts}.fasta")

    with open(csv_path, 'w', newline='') as f:
        fieldnames = list(construct.keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(construct)

    # Vaccine construction via Biopython SeqRecord, per methodology, rather
    # than hand-wrapping the FASTA text.
    record = SeqRecord(
        Seq(construct['Sequence']),
        id=construct['Construct_ID'],
        description=f"length={construct['Length']}",
    )
    with open(fasta_path, 'w') as f:
        SeqIO.write(record, f, "fasta")

    print("-" * 80)
    print(f"[INFO] Exported Matrix : {os.path.basename(csv_path)}")
    print(f"[INFO] Exported FASTA  : {os.path.basename(fasta_path)}")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    run_step1g_construction()
