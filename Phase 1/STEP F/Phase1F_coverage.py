import os, sys, time, csv, re, glob, subprocess, tempfile
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
# CONFIGURE -- override via env vars for a different machine
# =============================================================================
# MHCflurry and MHCnuggets are installed in a DEDICATED conda env
# ("mhcpredict"), not the shared "phase2" env -- both packages pull in a
# modern numpy/scikit-learn stack that breaks HemoPI2's pinned
# scikit-learn==1.3.1 when installed alongside it (confirmed empirically:
# installing them into phase2 upgraded scikit-learn to 1.9.0 and broke
# hemopi2_classification). Keeping them isolated avoids that conflict.
MHCFLURRY_PREDICT_BINARY = os.environ.get("MHCFLURRY_PREDICT_BINARY", "/opt/miniconda3/envs/mhcpredict/bin/mhcflurry-predict")
MHCPREDICT_PYTHON = os.environ.get("MHCPREDICT_PYTHON", "/opt/miniconda3/envs/mhcpredict/bin/python3")

# =============================================================================
# ALLELE PANEL & FREQUENCY TABLE
# =============================================================================
# Source: Allele Frequency Net Database (AFND), allelefrequencies.net
# Population: "Singapore Riau Malay" (AFND pop_id=2042, n=132), retrieved 2026-08-23.
# Loci A/B/C/DRB1/DQB1/DPB1, 2-field resolution, filtered to alleles supported
# by MHCflurry (class I) / MHCnuggets (class II).
#
# WHY A PROXY POPULATION (disclose this in the manuscript):
# AFND holds no usable Philippine table for this purpose. Verified directly:
#   - Philippines National Capital Region (n=51) -- Metro Manila, the paper's
#     target population -- is 1-field only (A*02, not A*02:01) and cannot be used.
#   - The only 2-field Philippine data is "Philippines Ivatan" (n=50), an
#     indigenous Batanes population not representative of lowland Filipinos.
#   - HLA-DQB1 and HLA-DPB1 have NO Philippine data at any resolution, and those
#     loci carry most of the MHC-II coverage.
# Singapore Riau Malay was selected as the closest usable proxy: Austronesian,
# the same lineage as lowland Filipino populations, with complete 2-field data
# across all six loci. This is a DOCUMENTED SUBSTITUTION -- results must be
# described as Southeast Asian (Austronesian) proxy coverage, NOT as Philippine
# population coverage, and the paper's "downloaded from IEDB" wording needs
# amending to cite AFND.
ALLELE_FREQ = {
    # HLA-A (22 alleles, freq sum = 0.816)
    "HLA-A*11:01": 0.177, "HLA-A*24:07": 0.157, "HLA-A*33:03": 0.109,
    "HLA-A*02:01": 0.065, "HLA-A*34:01": 0.06, "HLA-A*02:03": 0.052,
    "HLA-A*24:10": 0.04, "HLA-A*01:01": 0.032, "HLA-A*02:06": 0.02,
    "HLA-A*02:07": 0.02, "HLA-A*26:01": 0.02, "HLA-A*02:11": 0.012,
    "HLA-A*11:04": 0.012, "HLA-A*24:02": 0.008, "HLA-A*02:24": 0.004,
    "HLA-A*02:36": 0.004, "HLA-A*24:06": 0.004, "HLA-A*24:08": 0.004,
    "HLA-A*30:01": 0.004, "HLA-A*32:01": 0.004, "HLA-A*68:02": 0.004,
    "HLA-A*74:01": 0.004,
    # HLA-B (33 alleles, freq sum = 0.898)
    "HLA-B*18:01": 0.099, "HLA-B*15:02": 0.084, "HLA-B*15:13": 0.069,
    "HLA-B*35:05": 0.069, "HLA-B*44:03": 0.064, "HLA-B*40:01": 0.059,
    "HLA-B*13:01": 0.054, "HLA-B*58:01": 0.05, "HLA-B*38:02": 0.045,
    "HLA-B*15:21": 0.04, "HLA-B*40:06": 0.04, "HLA-B*35:01": 0.035,
    "HLA-B*07:05": 0.02, "HLA-B*57:01": 0.02, "HLA-B*18:03": 0.015,
    "HLA-B*46:01": 0.015, "HLA-B*07:02": 0.01, "HLA-B*15:25": 0.01,
    "HLA-B*27:04": 0.01, "HLA-B*27:06": 0.01, "HLA-B*35:02": 0.01,
    "HLA-B*40:10": 0.01, "HLA-B*55:01": 0.01, "HLA-B*13:02": 0.005,
    "HLA-B*35:03": 0.005, "HLA-B*35:17": 0.005, "HLA-B*39:15": 0.005,
    "HLA-B*40:02": 0.005, "HLA-B*42:02": 0.005, "HLA-B*51:06": 0.005,
    "HLA-B*54:01": 0.005, "HLA-B*56:02": 0.005, "HLA-B*56:04": 0.005,
    # HLA-C (19 alleles, freq sum = 0.945)
    "HLA-C*08:01": 0.187, "HLA-C*04:01": 0.122, "HLA-C*07:03": 0.112,
    "HLA-C*04:03": 0.098, "HLA-C*07:04": 0.084, "HLA-C*07:01": 0.075,
    "HLA-C*03:02": 0.051, "HLA-C*14:02": 0.047, "HLA-C*01:02": 0.042,
    "HLA-C*12:02": 0.037, "HLA-C*15:02": 0.023, "HLA-C*06:02": 0.019,
    "HLA-C*03:03": 0.014, "HLA-C*04:06": 0.009, "HLA-C*05:01": 0.005,
    "HLA-C*05:04": 0.005, "HLA-C*08:05": 0.005, "HLA-C*12:03": 0.005,
    "HLA-C*15:05": 0.005,
    # HLA-DRB1 (18 alleles, freq sum = 0.932)
    "HLA-DRB1*12:02": 0.324, "HLA-DRB1*15:01": 0.12, "HLA-DRB1*15:02": 0.12,
    "HLA-DRB1*07:01": 0.083, "HLA-DRB1*03:01": 0.046, "HLA-DRB1*09:01": 0.046,
    "HLA-DRB1*11:01": 0.046, "HLA-DRB1*12:01": 0.028, "HLA-DRB1*16:02": 0.028,
    "HLA-DRB1*13:02": 0.019, "HLA-DRB1*01:01": 0.009, "HLA-DRB1*04:03": 0.009,
    "HLA-DRB1*04:05": 0.009, "HLA-DRB1*08:02": 0.009, "HLA-DRB1*10:01": 0.009,
    "HLA-DRB1*11:04": 0.009, "HLA-DRB1*13:01": 0.009, "HLA-DRB1*14:01": 0.009,
    # HLA-DQB1 (7 alleles, freq sum = 0.800)
    "HLA-DQB1*03:01": 0.355, "HLA-DQB1*05:02": 0.118, "HLA-DQB1*02:01": 0.109,
    "HLA-DQB1*05:01": 0.1, "HLA-DQB1*06:02": 0.055, "HLA-DQB1*05:03": 0.045,
    "HLA-DQB1*03:02": 0.018,
    # HLA-DPB1 (8 alleles, freq sum = 0.639)
    "HLA-DPB1*05:01": 0.248, "HLA-DPB1*04:01": 0.171, "HLA-DPB1*03:01": 0.076,
    "HLA-DPB1*14:01": 0.038, "HLA-DPB1*01:01": 0.029, "HLA-DPB1*02:01": 0.029,
    "HLA-DPB1*04:02": 0.029, "HLA-DPB1*09:01": 0.019,
}
MHCI_ALLELES = [a for a in ALLELE_FREQ if a.startswith("HLA-A") or a.startswith("HLA-B") or a.startswith("HLA-C")]
MHCII_ALLELES = [a for a in ALLELE_FREQ if a.startswith("HLA-DRB1") or a.startswith("HLA-DQB1") or a.startswith("HLA-DPB1")]


def locus_of(allele):
    """HLA-A*24:02 -> 'A', HLA-DRB1*15:01 -> 'DRB1'"""
    m = re.match(r"HLA-([A-Z0-9]+)\*", allele)
    return m.group(1) if m else allele


def to_two_field(allele_raw):
    """
    Normalizes an allele string to 2-field notation (e.g. "HLA-A*02:01").
    Handles 4-field input (HLA-A*02:01:01:01 -> HLA-A*02:01) and bare
    "A*02:01" (adds the HLA- prefix). Previously this conversion step was
    entirely missing -- it only worked by accident because every allele
    already in the panel was hardcoded in 2-field form.
    """
    a = allele_raw.strip()
    if not a.startswith("HLA-"):
        a = "HLA-" + a
    parts = a.split(":")
    if len(parts) > 2:
        a = ":".join(parts[:2])
    return a


def compute_overall_coverage(binding_alleles, ndigits=2):
    """
    s = sum(freq[a]) per locus, capped at 1.0
    coverage_locus = 1 - (1-s)^2
    overall = 1 - prod_loci(1 - coverage_locus)
    Missing allele frequencies are treated as 0 (ALLELE_FREQ.get default).
    """
    by_locus = {}
    for allele in binding_alleles:
        locus = locus_of(allele)
        by_locus.setdefault(locus, 0.0)
        by_locus[locus] += ALLELE_FREQ.get(allele, 0.0)

    complement_product = 1.0
    for locus, s in by_locus.items():
        s = min(s, 1.0)
        coverage_locus = 1 - (1 - s) ** 2
        complement_product *= (1 - coverage_locus)
    # ndigits defaults to 2 so every existing per-epitope caller is unchanged.
    # The cumulative caller asks for more precision, because a set of epitopes
    # lands at 99.997% and rounding that to "100.00%" overstates the claim.
    return round((1 - complement_product) * 100, ndigits)


# =============================================================================
# MHCflurry (MHC-I) -- real percentile-rank prediction, batched across all
# peptide x allele combinations in one process call.
# =============================================================================

def predict_mhci_ranks(peptides, work_dir):
    """Returns {peptide: {allele: percentile_rank}}."""
    if not peptides:
        return {}
    batch_csv = os.path.join(work_dir, "mhci_batch_input.csv")
    out_csv = os.path.join(work_dir, "mhci_batch_output.csv")
    with open(batch_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["allele", "peptide"])
        for pep in peptides:
            for allele in MHCI_ALLELES:
                writer.writerow([allele, pep])

    cmd = [MHCFLURRY_PREDICT_BINARY, batch_csv, "--out", out_csv, "--no-flanking"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not os.path.isfile(out_csv):
        raise RuntimeError(f"mhcflurry-predict failed (exit {result.returncode}):\n{result.stderr}")

    ranks = {}
    skipped_alleles = set()
    with open(out_csv) as f:
        for row in csv.DictReader(f):
            pep, allele = row["peptide"], row["allele"]
            raw = (row.get("mhcflurry_affinity_percentile") or "").strip()
            if not raw:
                # MHCflurry lists some alleles as "supported" but returns an
                # empty percentile for them (no percentile-calibration data
                # behind that allele) -- observed for HLA-C*15:05. Treat as
                # "no prediction available" rather than crashing: the allele
                # simply contributes nothing to that peptide's binding set,
                # which is the same handling as an allele that fails the rank
                # cutoff. Reported below so it is never silently invisible.
                skipped_alleles.add(allele)
                continue
            try:
                ranks.setdefault(pep, {})[allele] = float(raw)
            except ValueError:
                skipped_alleles.add(allele)
    if skipped_alleles:
        print(f"\n[WARNING] MHCflurry returned no percentile for {len(skipped_alleles)} allele(s): "
              f"{', '.join(sorted(skipped_alleles))}")
        print("[WARNING] These contribute 0 to MHC-I coverage. Their ALLELE_FREQ mass is")
        print("          therefore unreachable -- note this if MHC-I coverage looks capped.")
    return ranks


# =============================================================================
# MHCnuggets (MHC-II) -- real percentile-rank prediction via its built-in
# rank_output mode (ranks IC50 against a human-proteome reference peptide
# set). Called once per allele (mhcnuggets' predict() API takes a single
# allele per call, unlike MHCflurry's batch-CSV mode).
# =============================================================================

_MHCNUGGETS_DRIVER = """
import sys, csv
from mhcnuggets.src.predict import predict
peptides_path, allele, output_path = sys.argv[1], sys.argv[2], sys.argv[3]
predict(class_='II', peptides_path=peptides_path, mhc=allele, output=output_path, rank_output=True)
"""


def predict_mhcii_ranks(peptides, work_dir):
    """Returns {peptide: {allele: percentile_rank}}. Percentile is
    human_proteome_rank * 100 (mhcnuggets reports it as a 0-1 fraction)."""
    if not peptides:
        return {}
    driver_path = os.path.join(work_dir, "_mhcnuggets_driver.py")
    with open(driver_path, "w") as f:
        f.write(_MHCNUGGETS_DRIVER)

    peptides_path = os.path.join(work_dir, "mhcii_peptides.txt")
    with open(peptides_path, "w") as f:
        f.write("\n".join(peptides))

    ranks = {}
    for allele in MHCII_ALLELES:
        safe_allele = re.sub(r"[^A-Za-z0-9]", "_", allele)
        output_path = os.path.join(work_dir, f"mhcii_out_{safe_allele}.csv")
        # MHCnuggets' own trained-model keys have no "*" (e.g. "HLA-DRB112:02"),
        # while our ALLELE_FREQ table uses standard notation with one
        # ("HLA-DRB1*12:02"). Passing the asterisk made every allele fail
        # closest_mhcII()'s exact match and silently fall back to the same
        # default model (HLA-DRB1*01:01) for all 9 alleles -- confirmed by
        # all 9 output files being byte-identical. Verified all 9 of our
        # alleles have real distinct trained models once the "*" is stripped.
        mhcnuggets_allele = allele.replace("*", "")
        cmd = [MHCPREDICT_PYTHON, driver_path, peptides_path, mhcnuggets_allele, output_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"\n[WARNING] MHCnuggets failed for {allele} (exit {result.returncode}): {result.stderr[-500:]}")
            continue

        # mhcnuggets has a real naming bug: it builds the rank filename as
        # f"{name}_ranks.{ext}" from os.path.splitext, whose ext already
        # includes the leading dot -- producing "..csv" (double dot), not
        # ".csv". Glob for either form so a future upstream fix doesn't
        # silently break this.
        base, ext = os.path.splitext(output_path)
        candidates = glob.glob(f"{base}_ranks.{ext.lstrip('.')}") + glob.glob(f"{base}_ranks*{ext}")
        if not candidates:
            print(f"\n[WARNING] MHCnuggets ran for {allele} but no rank output file was found.")
            continue

        with open(candidates[0]) as f:
            for row in csv.DictReader(f):
                pep = row["peptide"]
                pct = float(row["human_proteome_rank"]) * 100
                ranks.setdefault(pep, {})[allele] = pct
    return ranks


# =============================================================================
# CUMULATIVE (WHOLE-VACCINE) POPULATION COVERAGE  -- ADDITION, see below
#
# compute_overall_coverage() above answers "what fraction of the population
# has an allele that binds THIS ONE peptide?". That is a per-epitope
# prioritization signal, and it is what the >=90% gate uses.
#
# It is NOT the number that belongs in the abstract. A vaccine does not
# deliver one epitope -- it delivers all of them at once, and a person is
# covered if ANY epitope in the construct binds one of their alleles. The
# right figure is therefore computed over the UNION of binding alleles across
# every epitope in the construct, which is how IEDB's own Population Coverage
# tool reports a multi-epitope set.
#
# The difference is large and it runs in the vaccine's favour: individually
# only 4 of 23 HLA-restricted epitopes clear 90%, which reads badly and is
# the wrong comparison to make.
#
# Reported per class (MHC-I alone, MHC-II alone) and combined, because CD8
# and CD4 coverage are separate immunological claims and the manuscript
# discusses them separately.
# =============================================================================

def _fmt_cov(pct):
    """
    Never print a bare "100.0000%". Across six loci the complement product
    underflows to ~1e-9, so the formula returns a value indistinguishable
    from 100 -- which is a property of the arithmetic, not evidence that
    every single person is covered. Report it as a bound instead.
    """
    return ">99.99%" if pct >= 99.995 else f"{pct:.4f}%"


def compute_cumulative_coverage(allele_sets):
    """
    Population coverage of a SET of epitopes.

    allele_sets: iterable of per-epitope binding-allele lists. Their UNION is
    what the construct actually presents, so overlapping epitopes correctly
    stop double-counting the same allele.

    Uses the identical per-locus formula as compute_overall_coverage() --
    same Hardy-Weinberg two-chromosome assumption, same AFND frequencies --
    so the two numbers are directly comparable and no second methodology is
    introduced.
    """
    union = set()
    for alleles in allele_sets:
        union.update(alleles)
    union = sorted(union)

    # Per-locus breakdown. This is the informative layer: the across-loci
    # formula multiplies complements, so three loci at >85% already force the
    # overall figure above 99.9% almost regardless of the fourth. Reporting
    # only the combined number hides which locus is actually weak.
    per_locus = {}
    for allele in union:
        L = locus_of(allele)
        per_locus.setdefault(L, {"n": 0, "freq": 0.0})
        per_locus[L]["n"] += 1
        per_locus[L]["freq"] += ALLELE_FREQ.get(allele, 0.0)
    for L, d in per_locus.items():
        s = min(d["freq"], 1.0)
        d["coverage_pct"] = round((1 - (1 - s) ** 2) * 100, 2)
        d["freq_sum"] = s
        d["panel_n"] = sum(1 for a in ALLELE_FREQ if locus_of(a) == L)
        d["saturated"] = s >= 1.0

    return compute_overall_coverage(union, ndigits=4), union, per_locus


def _parse_alleles(cell):
    return [a for a in (cell or "").replace(",", ";").split(";") if a.strip()]


def run_construct_coverage():
    """
    Cumulative coverage of the epitopes actually in the Phase 1G construct.
    Read-only: prints and writes its own report, never touches 1F's outputs.
    """
    import glob, re, json
    from datetime import datetime as _dt

    g = sorted(glob.glob(os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1G",
                                      "Phase1G_FinalConstruct_*.csv")))
    f = sorted(glob.glob(os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1F",
                                      "Filtered", "*.csv")))
    if not g or not f:
        print("[ERROR] Need both a Phase 1G construct and a Phase 1F pool.")
        return
    with open(g[-1], newline="") as fh:
        crow = list(csv.DictReader(fh))[0]
    with open(f[-1], newline="") as fh:
        pool = {r["Peptide"]: r for r in csv.DictReader(fh)}

    shipped = []
    for seg in crow["Boundary_Map"].split(";"):
        m = re.match(r"^(MHC-I|MHC-II|B-cell):([A-Z]+)\[", seg.strip())
        if m:
            shipped.append((m.group(1), m.group(2)))

    print("\n" + "=" * 80)
    print("PHASE 1F (ADDENDUM): CUMULATIVE POPULATION COVERAGE OF THE CONSTRUCT")
    print("=" * 80)
    print(f"[INFO] Construct : {crow['Construct_ID']}")
    print(f"[INFO] Panel     : {len(ALLELE_FREQ)} alleles, AFND Singapore Riau Malay (proxy -- see deviation #7)")
    print("-" * 80)

    out_rows, summary = [], {}
    for cls in ("MHC-I", "MHC-II"):
        peps = [p for c, p in shipped if c == cls]
        sets = []
        print(f"\n{cls} epitopes in the construct: {len(peps)}")
        print(f"  {'PEPTIDE':<18}{'individual cov %':>18}{'binding alleles':>18}")
        for pep in peps:
            row = pool.get(pep)
            alleles = _parse_alleles(row.get("Binding_Alleles_Recomputed") if row else "")
            sets.append(alleles)
            ind = row.get("Overall_Coverage_Pct", "n/a") if row else "n/a"
            print(f"  {pep:<18}{str(ind):>18}{len(alleles):>18}")
        cum, union, per_locus = compute_cumulative_coverage(sets)
        summary[cls] = (cum, union, per_locus)
        print(f"  {'-'*70}")
        print(f"  CUMULATIVE {cls} COVERAGE : {_fmt_cov(cum)}   "
              f"(union of {len(union)} alleles across {len(peps)} epitopes)")
        print(f"  {'locus':<8}{'bound/panel':>14}{'freq sum':>11}{'locus coverage':>17}")
        for L, d in sorted(per_locus.items()):
            bound_of_panel = f"{d['n']}/{d['panel_n']}"
            flag = "  <-- SATURATED (freq sum >= 1.0)" if d["saturated"] else ""
            print(f"  {L:<8}{bound_of_panel:>14}{d['freq']:>11.3f}"
                  f"{d['coverage_pct']:>16.2f}%{flag}")

        out_rows.append({"Scope": cls, "N_Epitopes": len(peps),
                         "N_Union_Alleles": len(union), "Cumulative_Coverage_Pct": cum,
                         "Loci": ";".join(f"{k}:{d['n']}/{d['panel_n']}"
                                          f"@{d['coverage_pct']:.2f}%"
                                          for k, d in sorted(per_locus.items())),
                         "Union_Alleles": ";".join(union)})

    combined_sets = [summary[c][1] for c in summary]
    comb, comb_union, comb_loci = compute_cumulative_coverage(combined_sets)
    print("\n" + "-" * 80)
    print(f"COMBINED (MHC-I + MHC-II) CUMULATIVE COVERAGE : {_fmt_cov(comb)}")
    print(f"  union of {len(comb_union)} alleles across all HLA-restricted epitopes")
    print("-" * 80)
    print("NOTE: B-cell epitopes are excluded -- they are not HLA-restricted, so")
    print("      population coverage is not defined for them.")
    print("NOTE: report the PER-LOCUS figures alongside the combined one. The across-loci")
    print("      formula multiplies complements, so three loci above ~85% force the")
    print("      combined value over 99.9% almost regardless of the rest -- the combined")
    print("      number is real but nearly uninformative on its own. HLA-A at 57.10% is")
    print("      the genuine weak spot here and is what a reviewer should be shown.")
    print("NOTE: this is SOUTHEAST ASIAN (AUSTRONESIAN) PROXY coverage, not Philippine")
    print("      -- see deviation #7. Report it with that wording.")
    out_rows.append({"Scope": "COMBINED", "N_Epitopes": sum(r["N_Epitopes"] for r in out_rows),
                     "N_Union_Alleles": len(comb_union), "Cumulative_Coverage_Pct": comb,
                     "Loci": "", "Union_Alleles": ";".join(comb_union)})

    out_dir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1F", "Cumulative")
    os.makedirs(out_dir, exist_ok=True)
    ts = _dt.now().strftime("%Y%m%d_%H%M")
    path = os.path.join(out_dir, f"Phase1F_CumulativeCoverage_{ts}.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["Scope", "N_Epitopes", "N_Union_Alleles",
                                           "Cumulative_Coverage_Pct", "Loci", "Union_Alleles"])
        w.writeheader(); w.writerows(out_rows)
    print(f"\n[INFO] Report Saved : {os.path.relpath(path, _PROJECT_ROOT)}")
    print("=" * 80 + "\n")
    return out_rows


def run_step1f_population_coverage():
    start_time = time.time()

    # Phase 1Ec (human self-homology screen) now runs after 1Eb, per the
    # manuscript's stated order ("...surviving toxicity AND allergenicity
    # screening were additionally screened via BLASTP against...Swiss-Prot").
    input_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1E", "Phase1Ec", "Filtered")
    raw_dir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1F", "Raw")
    filt_dir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1F", "Filtered")
    tool_runs_dir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1F", "_tool_runs")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(filt_dir, exist_ok=True)
    os.makedirs(tool_runs_dir, exist_ok=True)

    print("\n" + "="*80 + "\nPHASE 1F: POPULATION COVERAGE ANALYSIS\n" + "="*80)
    print(f"[INFO] MHC-I alleles ({len(MHCI_ALLELES)}) via MHCflurry, MHC-II alleles ({len(MHCII_ALLELES)}) via MHCnuggets")

    if not os.path.isdir(input_folder):
        print(f"[ERROR] Phase 1Ec Filtered directory not found at: {input_folder}")
        return
    csv_files = [f for f in os.listdir(input_folder) if f.endswith(".csv")]
    if not csv_files:
        print("[ERROR] No candidates found from Phase 1Ec.")
        return
    latest_csv = os.path.join(input_folder, sorted(csv_files)[-1])

    with open(latest_csv, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        original_fields = reader.fieldnames

    mhci_peptides = sorted(set(r["Peptide"] for r in rows if r.get("Type") == "MHC-I"))
    mhcii_peptides = sorted(set(r["Peptide"] for r in rows if r.get("Type") == "MHC-II"))
    print(f"[INFO] {len(mhci_peptides)} unique MHC-I peptides, {len(mhcii_peptides)} unique MHC-II peptides, "
          f"{sum(1 for r in rows if r.get('Type') == 'B-cell')} B-cell rows (not HLA-restricted -- excluded from this gate)")

    print("[PROCESS] Running MHCflurry for MHC-I percentile ranks...")
    mhci_ranks = predict_mhci_ranks(mhci_peptides, tool_runs_dir)
    print("[PROCESS] Running MHCnuggets for MHC-II percentile ranks (one pass per allele)...")
    mhcii_ranks = predict_mhcii_ranks(mhcii_peptides, tool_runs_dir)

    raw_results, filtered_results = [], []
    fieldnames = original_fields + ["Binding_Alleles_Recomputed", "Overall_Coverage_Pct",
                                     "Coverage_Status", "Coverage_Priority"]

    # METHODOLOGY NOTE -- the 90% criterion is a PRIORITIZATION, not an exclusion.
    #
    # Section II.C.I.F states: "Peptides with overall_coverage >= 90% were
    # PRIORITIZED as population-covering candidates" -- prioritized, not
    # required. This is the same wording pattern already corrected elsewhere
    # in this pipeline for GRAVY ("deprioritized") and Surface_Charge
    # ("deprioritizes"), both of which had likewise been implemented as hard
    # drops. A prior version excluded every MHC-I/MHC-II row below 90%.
    #
    # That exclusion was not survivable for MHC-II: measured cumulative
    # coverage across ALL MHC-II candidates saturates well below 90%, so a
    # hard gate produced a construct with zero HTL epitopes -- no GPGPG
    # linkers and no CD4+ help, contradicting the construct design in
    # Section I.G and making the Phase III helper-T-cell research questions
    # unanswerable. Rows are now retained and ranked; Phase 1G consumes
    # Coverage_Priority to give >=90% peptides genuine precedence.
    n_priority = {"MHC-I": 0, "MHC-II": 0}
    n_carried = {"MHC-I": 0, "MHC-II": 0}

    for i, row in enumerate(rows):
        pep, ptype = row["Peptide"], row.get("Type", "")
        clean_row = {k: row[k] for k in original_fields}

        if ptype in ("MHC-I", "MHC-II"):
            if ptype == "MHC-I":
                per_allele = mhci_ranks.get(pep, {})
                rank_cutoff = 2.0
            else:
                per_allele = mhcii_ranks.get(pep, {})
                rank_cutoff = 10.0
            binding = [to_two_field(a) for a, rank in per_allele.items() if rank <= rank_cutoff]
            coverage = compute_overall_coverage(binding)
            meets = coverage >= 90.0
            clean_row["Binding_Alleles_Recomputed"] = ";".join(sorted(binding))
            clean_row["Overall_Coverage_Pct"] = coverage
            clean_row["Coverage_Status"] = "PASS" if meets else "FAIL"
            clean_row["Coverage_Priority"] = "PRIORITIZED" if meets else "BELOW_TARGET"
            if meets:
                n_priority[ptype] += 1
            else:
                n_carried[ptype] += 1
            filtered_results.append(clean_row)
        else:
            # B-cell: population coverage is an HLA-restriction concept and
            # does not apply to B-cell (antibody) epitopes. A prior version
            # of this script hardcoded these to a fabricated 100.0 "pass" --
            # they are instead marked not-applicable and passed through
            # unfiltered by this specific gate, not silently scored.
            clean_row["Binding_Alleles_Recomputed"] = "N/A"
            clean_row["Overall_Coverage_Pct"] = "N/A"
            clean_row["Coverage_Status"] = "NOT_APPLICABLE (not HLA-restricted)"
            clean_row["Coverage_Priority"] = "NOT_APPLICABLE"
            filtered_results.append(clean_row)

        raw_results.append(clean_row)
        if i % 5 == 0 or i == len(rows) - 1:
            elapsed = common.format_time(time.time() - start_time)
            sys.stdout.write(f"\r[ PROCESS ] {i+1:03d}/{len(rows):03d} | Elapsed: {elapsed}")
            sys.stdout.flush()

    print()
    ts = datetime.now().strftime("%Y%m%d_%H%M")
    if raw_results:
        with open(os.path.join(raw_dir, f"Phase1F_Raw_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(raw_results)
    if filtered_results:
        with open(os.path.join(filt_dir, f"Phase1F_Elite_Vaccine_Candidates_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(filtered_results)

    n_mhci = sum(1 for r in filtered_results if r.get("Type") == "MHC-I")
    n_mhcii = sum(1 for r in filtered_results if r.get("Type") == "MHC-II")
    n_bcell = sum(1 for r in filtered_results if r.get("Type") == "B-cell")
    print(f"[SUCCESS] Candidates carried forward: {len(filtered_results)}/{len(rows)} "
          f"(MHC-I: {n_mhci}, MHC-II: {n_mhcii}, B-cell: {n_bcell})")
    print(f"[INFO] Coverage >=90% (PRIORITIZED) : MHC-I={n_priority['MHC-I']} MHC-II={n_priority['MHC-II']}")
    print(f"[INFO] Coverage <90% (BELOW_TARGET, retained & ranked lower) : "
          f"MHC-I={n_carried['MHC-I']} MHC-II={n_carried['MHC-II']}")
    print("[INFO] Per the methodology, >=90% is a prioritization criterion, not an")
    print("       exclusion -- Phase 1G ranks by coverage so prioritized peptides win.")

if __name__ == "__main__":
    import argparse as _ap
    _p = _ap.ArgumentParser(description="Step 1F: Population Coverage Analysis")
    _p.add_argument("--construct-coverage", action="store_true",
                    help="Report CUMULATIVE coverage of the Phase 1G construct's epitope set. "
                         "Read-only -- does not re-run or overwrite the 1F pool.")
    _a = _p.parse_args()
    if _a.construct_coverage:
        run_construct_coverage()
    else:
        run_step1f_population_coverage()
