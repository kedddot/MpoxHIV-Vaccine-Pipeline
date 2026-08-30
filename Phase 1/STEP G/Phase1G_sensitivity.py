import os
import sys
import csv
import random
import hashlib
import argparse
import importlib.util
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
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# =============================================================================
# PHASE 1G SENSITIVITY AUDIT -- would SEMA have changed the construct?
#
# THE QUESTION. Deviation #20 justifies never running SEMA with: "every B-cell
# slot in the final construct is FORCED, so SEMA could not have changed the
# selection even with complete data." That claim is false. Re-derived from the
# shipped Phase 1F pool under Phase 1G's own screens, Mpox_B5R has SEVEN
# genuinely distinct B-cell candidates competing for two slots and Mpox_L1R has
# three for two; only gp120/gp41/p17/p24/A35R are actually forced. Since
# SEMA_Corroborated=="YES" feeds n_prioritized -- the third key of 1G's
# lexicographic objective -- SEMA can decide between selections tied at the two
# hard gates. So the honest answer has to be measured, not argued.
#
# HOW IT STAYS HONEST. This module IMPORTS Phase1G_construction rather than
# copying its objective, so the two cannot drift apart. It replays the same
# 6000-trial search with the same RNG_SEED, changing ONLY the SEMA labels, and
# it verifies up front that the untouched pool reproduces Vax_Final_6f34b53e
# byte-for-byte -- an audit whose baseline does not reproduce the shipped
# construct is measuring its own bugs.
#
# WHAT IT NEVER DOES. It never writes to Step_Outputs/Phase1/Phase1G/. Phase 1G
# picks its input and its own history via common.latest_file(), so dropping a
# construct CSV there would silently redefine "the construct". Everything lands
# under Phase1G/Sensitivity/.
#
# MEASURED ENVELOPE (real pool, RNG_SEED=20260824):
#   baseline  all UNSCREENED            -> 6f34b53e, 566 aa  (== shipped)
#   best      only selected corroborated-> 6f34b53e, 566 aa, 100.0% identity
#   worst     only unselected corroborated-> 24d8424e, 587 aa, 27.2% identity
# The worst case moves BOTH T-cell blocks, not just the B-cell cassette,
# because n_prioritized sums BigMHC MHC-I priority and SEMA B-cell
# corroboration into ONE scalar -- so two corroborated B-cell epitopes can
# outrank two BigMHC-prioritised CTL epitopes. Those are not commensurable
# quantities and should be separate keys in the sort tuple; recorded as a
# finding rather than silently fixed, since changing it would alter the
# construct.
# =============================================================================

_G_PATH = os.path.join(_PROJECT_ROOT, "Phase 1", "STEP G", "Phase1G_construction.py")
_spec = importlib.util.spec_from_file_location("phase1g_construction", _G_PATH)
G = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(G)

# Phase 1G defines TOP_K inside run_step1g_construction(), so it is not
# importable. Mirrored here with the same value and the same rationale; if it
# changes there it must change here, which is why it is called out loudly.
TOP_K = 60


def _shipped_sequence():
    folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1G")
    fastas = [f for f in os.listdir(folder) if f.endswith(".fasta")] if os.path.isdir(folder) else []
    if not fastas:
        return None, None
    fastas.sort(key=lambda f: os.path.getctime(os.path.join(folder, f)))
    path = os.path.join(folder, fastas[-1])
    with open(path) as f:
        seq = "".join(l.strip() for l in f if not l.startswith(">"))
    return seq, os.path.basename(path)


def replay(by_target, priority_info, label_fn):
    """
    Re-runs Phase 1G's selection+ordering search with SEMA labels replaced by
    label_fn(peptide). Every other term is Phase 1G's own code.
    """
    for pep, info in priority_info.items():
        info["sema"] = label_fn(pep)

    rng = random.Random(G.RNG_SEED)
    top = []
    for trial in range(G.SEARCH_TRIALS):
        gravy_weight = G.GRAVY_WEIGHT_SWEEP[trial % len(G.GRAVY_WEIGHT_SWEEP)]
        mhc_i, mhc_ii, bcell, provenance = G._select_candidates(by_target, rng, gravy_weight)
        if not (mhc_i or mhc_ii or bcell):
            continue
        mi, p1 = G._greedy_order(mhc_i, G.L_MHCI, G.ADJUVANT + G.ADJ_LINKER)
        mii, p2 = G._greedy_order(mhc_ii, G.L_MHCII, p1 + (G.L_MHCII if mhc_i else ""))
        bc, _p3 = G._greedy_order(bcell, G.L_BCELL, p2 + (G.L_BCELL if (mhc_i or mhc_ii) else ""))
        seq = G._assemble(mi, mii, bc)
        n_prioritized = (
            sum(1 for p in mi if priority_info.get(p, {}).get("bigmhc") == "PRIORITIZED")
            + sum(1 for p in bc if priority_info.get(p, {}).get("sema") == "YES"))
        coverage = sum(cov for t in by_target for tt in by_target[t]
                       for _pr, cov, _cn, pep in by_target[t][tt]
                       if pep in (mi + mii + bc))
        score = (G.count_hydrophobic_windows(seq),
                 G._spurious_motif_count(seq, mi, mii, bc),
                 -n_prioritized,
                 round(ProteinAnalysis(seq).gravy(), 4),
                 -coverage)
        top.append((score, mi, mii, bc, provenance))
        top.sort(key=lambda c: c[0])
        del top[TOP_K:]

    best = None
    for score, mi, mii, bc, provenance in top:
        p_mi, p_mii, p_bc, polished = G._two_opt_polish(mi, mii, bc)
        rank_key = polished + score[2:]
        if best is None or rank_key < best[0]:
            best = (rank_key, polished, p_mi, p_mii, p_bc, provenance)

    _rk, polished, mi, mii, bc, provenance = best
    seq = G._assemble(mi, mii, bc)
    return {"sequence": seq, "hash": hashlib.md5(seq.encode()).hexdigest()[:8],
            "length": len(seq), "windows": polished[0], "motifs": polished[1],
            "gravy": round(ProteinAnalysis(seq).gravy(), 4),
            "mhc_i": mi, "mhc_ii": mii, "bcell": bc}


def _identity(a, b):
    import difflib
    return 100.0 * difflib.SequenceMatcher(None, a, b).ratio()


def _load_real_labels():
    """
    Reads SEMA_Corroborated from the Phase 1De sensitivity output, if it has
    been produced. Returns None when SEMA has not been run -- in which case
    only the bounding cases are meaningful.
    """
    folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D",
                          "Phase1De", "Sensitivity")
    latest = common.latest_file(folder, suffix=".csv")
    if latest is None:
        return None, None
    labels = {}
    with open(latest) as f:
        for row in csv.DictReader(f):
            call = (row.get("SEMA_Corroborated") or "").strip()
            if row.get("Peptide") and call:
                labels[row["Peptide"]] = call
    return (labels or None), os.path.basename(latest)


def run_sensitivity():
    common.print_banner("PHASE 1G SENSITIVITY: WOULD SEMA HAVE CHANGED THE CONSTRUCT?")

    pool_dir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1F", "Filtered")
    latest_csv = common.latest_file(pool_dir, suffix=".csv")
    if latest_csv is None:
        print(f"[ERROR] No Phase 1F pool at: {pool_dir}")
        return
    by_target, _exc_i, _exc_m, priority_info = G._load_pool(latest_csv)
    shipped, shipped_name = _shipped_sequence()
    if shipped is None:
        print("[ERROR] No shipped construct FASTA found.")
        return

    print(f"[INFO] Pool     : {os.path.basename(latest_csv)}")
    print(f"[INFO] Shipped  : {shipped_name} ({len(shipped)} aa)")
    print(f"[INFO] Replaying {G.SEARCH_TRIALS} trials per case (seed={G.RNG_SEED}, TOP_K={TOP_K})...")

    # ---- baseline: the audit must reproduce the shipped construct exactly ----
    base = replay(by_target, priority_info, lambda p: "UNSCREENED")
    reproduces = (base["sequence"] == shipped)
    print(f"\n[BASELINE] all UNSCREENED -> {base['hash']} ({base['length']} aa) | "
          f"reproduces shipped construct: {reproduces}")
    if not reproduces:
        print("[FATAL] Baseline does not reproduce the shipped construct. The audit would be "
              "measuring its own drift, not SEMA's effect. Refusing to report diffs.")
        return

    cases = [("BASELINE_all_unscreened", base)]

    # ---- bounding cases ------------------------------------------------------
    selected_bcell = set(base["bcell"])
    contestable = set()
    for target in by_target:
        window = [c[3] for c in by_target[target]["B-cell"][:G.SELECTION_WINDOW]]
        for pep in window:
            if pep not in selected_bcell:
                contestable.add(pep)

    cases.append(("WORST_only_unselected_corroborated",
                  replay(by_target, priority_info,
                         lambda p: "YES" if p in contestable else "NO")))
    cases.append(("BEST_only_selected_corroborated",
                  replay(by_target, priority_info,
                         lambda p: "YES" if p in selected_bcell else "NO")))

    # ---- the real thing, if SEMA has actually run ---------------------------
    real_labels, real_src = _load_real_labels()
    if real_labels:
        print(f"[INFO] Real SEMA labels loaded from: {real_src} ({len(real_labels)} peptides)")
        cases.append(("REAL_sema3d",
                      replay(by_target, priority_info,
                             lambda p: real_labels.get(p, "UNSCREENED"))))
    else:
        print("[INFO] No Phase 1De SEMA output yet -- reporting the bounding envelope only. "
              "Re-run this after SEMA-3D to get the REAL verdict.")

    # ---- report --------------------------------------------------------------
    out_dir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1G", "Sensitivity")
    os.makedirs(out_dir, exist_ok=True)
    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    out_path = os.path.join(out_dir, f"Phase1G_Sensitivity_{ts}.csv")

    rows = []
    print("\n" + "-" * 100)
    for name, res in cases:
        changed = (res["sequence"] != shipped)
        ident = _identity(res["sequence"], shipped)
        added = sorted(set(res["bcell"]) - set(base["bcell"]))
        removed = sorted(set(base["bcell"]) - set(res["bcell"]))
        tcell_changed = (set(res["mhc_i"]) != set(base["mhc_i"])
                         or set(res["mhc_ii"]) != set(base["mhc_ii"]))
        verdict = "CHANGED" if changed else "UNCHANGED"
        print(f"{name}")
        print(f"   {verdict:<9} hash={res['hash']} len={res['length']} identity={ident:.1f}% "
              f"windows={res['windows']} motifs={res['motifs']} GRAVY={res['gravy']}")
        if changed:
            print(f"   T-cell blocks changed: {tcell_changed}")
            print(f"   B-cell added   : {added}")
            print(f"   B-cell removed : {removed}")
        rows.append({"Case": name, "Verdict": verdict, "Construct_Hash": res["hash"],
                     "Length": res["length"], "Identity_To_Shipped_Pct": round(ident, 1),
                     "Bad_Windows": res["windows"], "Spurious_Motifs": res["motifs"],
                     "GRAVY": res["gravy"], "TCell_Blocks_Changed": tcell_changed,
                     "BCell_Added": ";".join(added), "BCell_Removed": ";".join(removed),
                     "Sequence": res["sequence"]})

    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print("\n" + "-" * 100)
    if real_labels:
        real = next(r for r in rows if r["Case"] == "REAL_sema3d")
        print(f"[VERDICT] With real SEMA-3D data the construct would be {real['Verdict']}.")
        if real["Verdict"] == "CHANGED":
            print(f"          Substituted B-cell epitopes -- added: {real['BCell_Added'] or 'none'} | "
                  f"removed: {real['BCell_Removed'] or 'none'}")
            print("          Deviation #20 must be restated: SEMA became available after the "
                  "construct was fixed, and would have altered it.")
        else:
            print("          Deviation #20 can be retired: SEMA was run and the construct is "
                  "unchanged under its real labels.")
    else:
        print("[VERDICT] PENDING -- bounding envelope only; SEMA-3D has not been run yet.")
    print(f"[INFO] Report written to: {out_path}")
    print("[INFO] Nothing under Step_Outputs/Phase1/Phase1G/ itself was modified.")
    print("=" * 100 + "\n")


if __name__ == "__main__":
    argparse.ArgumentParser(description="Phase 1G SEMA sensitivity audit (read-only).").parse_args()
    run_sensitivity()
