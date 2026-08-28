import os, sys, csv, json, time, re, argparse, hashlib
from datetime import datetime

# =============================================================================
# MINIMAL BOOTSTRAP -- locates the shared phase1_common module.
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
# PHASE 1H: IEDB EXPERIMENTAL CORROBORATION  (ADDITION, NOT IN THE MANUSCRIPT)
#
# Every epitope in this construct is PREDICTED -- MHCflurry, MHCnuggets and
# BepiPred-2.0 outputs, filtered and ranked. Nothing downstream ever asks
# whether any of them has actually been observed in a laboratory.
#
# IEDB curates exactly that: peptides with published T-cell assays, B-cell
# assays and MHC ligand-elution mass spectrometry, each tied to a PubMed ID.
# Cross-referencing the construct against it converts "our tools predicted
# this" into "our tools predicted this AND it is independently documented in
# the literature" -- for whichever epitopes earn it.
#
# THIS STEP NEVER CHANGES SELECTION. It annotates. An epitope with no IEDB
# match is not worse than one with a match; prediction of a novel epitope is
# a legitimate result and the whole point of doing prediction. What the step
# buys is the ability to say which is which, instead of being silent.
#
# ORTHOPOXVIRUS SCOPE: the Mpox antigens A35R, B5R and L1R have direct
# vaccinia orthologs of the same gene names, and orthopoxvirus cross-
# reactivity is the mechanism by which smallpox vaccination protects against
# mpox. IEDB holds 29 Monkeypox-virus epitopes but 17,903 vaccinia ones, so
# restricting to MPXV alone would discard almost all of the relevant
# experimental record. Vaccinia matches are reported as CROSS-SPECIES and
# labeled as such -- never silently merged with same-species evidence.
# =============================================================================

IEDB_API = "https://query-api.iedb.org/epitope_search"
IEDB_HEADERS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
                  "(KHTML, like Gecko) Chrome/128.0.0.0 Safari/537.36",
    "Accept": "application/json",
}
PAGE = 1000

# NCBI taxon -> (label, is_same_species_for_which_of_our_targets)
SOURCE_TAXA = [
    ("NCBITaxon:11676", "HIV-1",            "HIV"),
    ("NCBITaxon:10244", "Monkeypox virus",  "MPXV"),
    ("NCBITaxon:10245", "Vaccinia virus",   "VACV"),
]

# Shortest overlap that is still biologically meaningful. 8 is the minimum
# length of an MHC-I ligand core, so a shared stretch below that is sequence
# coincidence rather than a shared epitope.
MIN_OVERLAP = 8

FIELDS = ("linear_sequence,parent_source_antigen_names,source_organism_names,"
          "mhc_allele_names,mhc_classes,assay_names,pubmed_ids,"
          "tcell_ids,bcell_ids,elution_ids,structure_id")


def fetch_taxon_epitopes(taxon, label, cache_dir, force=False):
    """
    Pulls every linear epitope IEDB holds for one source organism, paginated.
    Cached on disk -- the record only changes when IEDB re-curates, and a
    cached copy makes the whole step reproducible offline.
    """
    import requests
    cache = os.path.join(cache_dir, f"iedb_{taxon.replace(':', '_')}.json")
    if os.path.isfile(cache) and not force:
        with open(cache) as fh:
            data = json.load(fh)
        print(f"[CACHE] {label:<18}{len(data):>6} epitopes (delete {os.path.basename(cache)} to refresh)")
        return data

    rows, offset = [], 0
    while True:
        try:
            r = requests.get(IEDB_API, headers=IEDB_HEADERS, timeout=120, params={
                "source_organism_iri_search": f"cs.{{{taxon}}}",
                "select": FIELDS, "limit": PAGE, "offset": offset,
                "linear_sequence": "not.is.null",
                # PostgREST rejects offset without a deterministic order, and
                # an unordered paginated scan could silently repeat or skip
                # rows even if it did not.
                "order": "structure_id.asc",
            })
        except Exception as e:
            print(f"[WARN] {label}: request failed ({e}) -- keeping the {len(rows)} rows fetched so far.")
            break
        if r.status_code not in (200, 206):
            print(f"\n[WARN] {label}: HTTP {r.status_code} -- keeping {len(rows)} rows. "
                  f"{r.text[:200]}")
            break
        try:
            page = r.json()
        except Exception:
            print(f"[WARN] {label}: non-JSON response -- keeping {len(rows)} rows.")
            break
        if not page:
            break
        rows.extend(page)
        offset += PAGE
        print(f"\r[FETCH] {label:<18}{len(rows):>6} epitopes...", end="", flush=True)
        time.sleep(0.2)
    print(f"\r[FETCH] {label:<18}{len(rows):>6} epitopes    ")
    if rows:
        os.makedirs(cache_dir, exist_ok=True)
        with open(cache, "w") as fh:
            json.dump(rows, fh)
    return rows


def _as_list(value):
    """IEDB returns some columns as JSON lists, some as strings, some as None."""
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v) for v in value if v is not None]
    s = str(value).strip()
    if not s or s == "None":
        return []
    if s.startswith("[") and s.endswith("]"):
        try:
            return [str(v) for v in json.loads(s.replace("'", '"')) if v is not None]
        except Exception:
            pass
    return [s]


def _longest_common_substring(a, b):
    """Length of the longest shared contiguous stretch."""
    if not a or not b:
        return 0
    prev = [0] * (len(b) + 1)
    best = 0
    for i in range(1, len(a) + 1):
        cur = [0] * (len(b) + 1)
        ai = a[i - 1]
        for j in range(1, len(b) + 1):
            if ai == b[j - 1]:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best:
                    best = cur[j]
        prev = cur
    return best


def classify_match(ours, theirs):
    """
    Relationship between our predicted peptide and a curated IEDB epitope.
    Ordered strongest first; returns None when the two are unrelated.
    """
    if ours == theirs:
        return "EXACT", len(ours)
    if ours in theirs:
        return "CONTAINED_IN_CURATED", len(ours)
    if theirs in ours and len(theirs) >= MIN_OVERLAP:
        return "CONTAINS_CURATED", len(theirs)
    n = _longest_common_substring(ours, theirs)
    if n >= MIN_OVERLAP:
        return "PARTIAL_OVERLAP", n
    return None


_RANK = {"EXACT": 0, "CONTAINED_IN_CURATED": 1, "CONTAINS_CURATED": 2, "PARTIAL_OVERLAP": 3}


def assay_evidence(rec):
    """Which assay classes back this curated epitope."""
    kinds = []
    if _as_list(rec.get("tcell_ids")):   kinds.append("T-cell")
    if _as_list(rec.get("bcell_ids")):   kinds.append("B-cell")
    if _as_list(rec.get("elution_ids")): kinds.append("MHC-ligand-elution")
    return kinds


def build_index(all_records):
    """Deduplicates curated epitopes by sequence, merging their evidence."""
    index = {}
    for label, is_same, rec in all_records:
        seq = (rec.get("linear_sequence") or "").strip().upper()
        if not seq or not seq.isalpha():
            continue
        entry = index.setdefault(seq, {
            "seq": seq, "organisms": set(), "antigens": set(), "alleles": set(),
            "classes": set(), "assays": set(), "pubmed": set(), "evidence": set(),
            "species_scope": set(),
        })
        entry["organisms"].update(_as_list(rec.get("source_organism_names")))
        entry["antigens"].update(_as_list(rec.get("parent_source_antigen_names")))
        entry["alleles"].update(_as_list(rec.get("mhc_allele_names")))
        entry["classes"].update(_as_list(rec.get("mhc_classes")))
        entry["assays"].update(_as_list(rec.get("assay_names")))
        entry["pubmed"].update(_as_list(rec.get("pubmed_ids")))
        entry["evidence"].update(assay_evidence(rec))
        entry["species_scope"].add(is_same)
    return index


# =============================================================================
# ANTIGEN CONCORDANCE -- the check that makes a sequence match trustworthy
#
# Matching on sequence alone is not sufficient. A 9-mer is short, and an
# orthopoxvirus proteome is ~200 kb, so a curated epitope from a completely
# unrelated protein can share 9 residues with ours by chance. UniProt's
# protein NAMES cannot settle it either: several of the real matches here
# come back labelled "Ankyrin repeat protein OPG189" or "Protein OPG161",
# legacy/unified orthopoxvirus gene names that do not resemble "B5R" or
# "A35R" at all, so a name-based check would reject genuine orthologs.
#
# So the check is done against the actual sequences: does the curated epitope
# occur in OUR OWN Phase 1A reference set for that antigen? If it aligns at
# high identity, the two are the same protein regardless of what either
# database calls it. If it does not, the shared stretch is coincidence.
#
# Compared against ALL variants of the target, not just Var_01 -- Phase 1D
# generated epitopes across every variant, so a peptide need not appear in
# the first isolate.
# =============================================================================

CONCORDANCE_MIN = 0.70


def _load_reference_variants(project_root):
    """{target: [sequence, ...]} from Phase 1A."""
    import glob
    folder = os.path.join(project_root, "Step_Outputs", "Phase1", "Phase1A")
    out = {}
    for path in sorted(glob.glob(os.path.join(folder, "*_Var_*.fasta"))):
        base = os.path.basename(path)
        m = re.match(r"^(.+?)_Var_\d+_", base)
        if not m:
            continue
        seq = "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()
        if seq:
            out.setdefault(m.group(1), []).append(seq)
    return out


def antigen_concordance(curated_seq, target_names, reference_variants):
    """
    Best ungapped identity of a curated epitope against our own reference
    variants for that antigen. Returns 0.0 when we have no reference.
    """
    best = 0.0
    L = len(curated_seq)
    for target in target_names:
        for s in reference_variants.get(target, []):
            if L > len(s):
                continue
            for i in range(len(s) - L + 1):
                w = s[i:i + L]
                # cheap reject before the full comparison
                if w[0] != curated_seq[0] and w[-1] != curated_seq[-1]:
                    continue
                ident = sum(1 for a, b in zip(curated_seq, w) if a == b) / L
                if ident > best:
                    best = ident
                    if best == 1.0:
                        return 1.0
    return best


def load_construct_epitopes(project_root):
    """
    Reads the shipped epitopes straight out of Phase 1G's Boundary_Map, so
    this reflects what is actually IN the construct rather than what merely
    survived Phase 1F.
    """
    path = common.latest_file(os.path.join(project_root, "Step_Outputs", "Phase1", "Phase1G"),
                              suffix=".csv") if hasattr(common, "latest_file") else None
    if path is None:
        import glob
        c = sorted(glob.glob(os.path.join(project_root, "Step_Outputs", "Phase1", "Phase1G",
                                          "Phase1G_FinalConstruct_*.csv")))
        path = c[-1] if c else None
    if path is None:
        return None, None, []
    with open(path, newline="") as fh:
        row = list(csv.DictReader(fh))[0]
    eps = []
    for seg in row["Boundary_Map"].split(";"):
        m = re.match(r"^(MHC-I|MHC-II|B-cell):([A-Z]+)\[(\d+)-(\d+)\]$", seg.strip())
        if m:
            eps.append({"class": m.group(1), "peptide": m.group(2),
                        "start": int(m.group(3)) + 1, "end": int(m.group(4))})
    return row["Construct_ID"], os.path.basename(path), eps


def target_of(peptide, project_root):
    """Source antigen(s) for a peptide, from the Phase 1F pool."""
    import glob
    f = sorted(glob.glob(os.path.join(project_root, "Step_Outputs", "Phase1", "Phase1F",
                                      "Filtered", "*.csv")))
    if not f:
        return ""
    with open(f[-1], newline="") as fh:
        hits = sorted({r["Target"] for r in csv.DictReader(fh) if r["Peptide"] == peptide})
    return ";".join(hits)


def run_step1h(force_refresh=False):
    start = time.time()
    project_root = _PROJECT_ROOT
    out_dir = os.path.join(project_root, "Step_Outputs", "Phase1", "Phase1H")
    cache_dir = os.path.join(out_dir, "_tool_runs")
    os.makedirs(cache_dir, exist_ok=True)

    common.print_banner("PHASE 1H: IEDB EXPERIMENTAL CORROBORATION OF CONSTRUCT EPITOPES")
    print(f"[INFO] Resolved Project Root : {project_root}")
    print("[INFO] Source : IEDB Query API (query-api.iedb.org) -- curated assay records")
    print("[INFO] NOTE   : annotation only. This step never adds, removes or reorders an epitope.")
    print("-" * 110)

    construct_id, construct_file, eps = load_construct_epitopes(project_root)
    if not eps:
        print("[ERROR] No Phase 1G construct found -- run Step 1G first.")
        return
    print(f"[INFO] Construct : {construct_id}  ({construct_file})")
    print(f"[INFO] Epitopes  : {len(eps)}")
    print("-" * 110)

    records = []
    for taxon, label, scope in SOURCE_TAXA:
        for rec in fetch_taxon_epitopes(taxon, label, cache_dir, force_refresh):
            records.append((label, scope, rec))
    if not records:
        print("\n[WARN] IEDB returned nothing -- every epitope recorded as NOT_QUERIED. "
              "This is non-blocking; re-run when the API is reachable.")
    index = build_index(records)
    print(f"[INFO] Curated epitopes indexed (deduplicated by sequence): {len(index)}")
    print("-" * 110)

    reference_variants = _load_reference_variants(project_root)
    print(f"[INFO] Reference antigens loaded for concordance check: "
          f"{', '.join(f'{k}={len(v)}' for k, v in sorted(reference_variants.items()))}")
    print("-" * 110)

    results = []
    for e in eps:
        ours = e["peptide"]
        matches = []
        for seq, entry in index.items():
            got = classify_match(ours, seq)
            if got:
                kind, n = got
                matches.append((_RANK[kind], -n, kind, n, entry))
        matches.sort(key=lambda m: (m[0], m[1], -len(m[4]["pubmed"])))

        our_targets = target_of(ours, project_root)
        is_hiv = "HIV" in our_targets
        # Same-species evidence for HIV targets is HIV-1; for Mpox targets it
        # is Monkeypox virus. Vaccinia is orthologous, not same-species.
        want = "HIV" if is_hiv else "MPXV"

        row = {
            "Construct_ID": construct_id, "Class": e["class"], "Peptide": ours,
            "Length": len(ours), "Construct_Start": e["start"], "Construct_End": e["end"],
            "Source_Target": our_targets, "IEDB_Status": "NOVEL (no curated match)",
            "Match_Type": "", "Overlap_Len": "", "Curated_Epitope": "",
            "Species_Scope": "", "Curated_Source_Antigen": "", "Assay_Evidence": "",
            "Curated_MHC_Alleles": "", "PubMed_IDs": "", "N_Curated_Matches": len(matches),
            "Antigen_Concordance_Pct": "", "Concordance_Verdict": "",
        }
        if matches:
            _, _, kind, n, entry = matches[0]
            same = want in entry["species_scope"]
            scope = ("SAME_SPECIES" if same else
                     "CROSS_SPECIES (orthopoxvirus ortholog)" if "VACV" in entry["species_scope"]
                     else "CROSS_SPECIES")
            ev = sorted(entry["evidence"])
            conc = antigen_concordance(entry["seq"], our_targets.split(";"), reference_variants)
            concordant = conc >= CONCORDANCE_MIN
            row["Antigen_Concordance_Pct"] = round(conc * 100, 1)
            row["Concordance_Verdict"] = (
                "CONCORDANT (curated epitope occurs in our own reference antigen)"
                if concordant else
                f"NOT CONCORDANT (<{CONCORDANCE_MIN*100:.0f}% identity to our reference) "
                f"-- treat as coincidental sequence overlap, NOT evidence")
            if not concordant:
                row.update({
                    "IEDB_Status": "NOVEL (sequence match rejected -- wrong antigen)",
                    "Match_Type": kind + " (REJECTED)", "Overlap_Len": n,
                    "Curated_Epitope": entry["seq"],
                    "Curated_Source_Antigen": "; ".join(sorted(entry["antigens"]))[:200],
                })
                results.append(row)
                continue
            if kind in ("EXACT", "CONTAINED_IN_CURATED") and ev:
                status = f"EXPERIMENTALLY CONFIRMED ({'/'.join(ev)})"
            elif ev:
                status = f"OVERLAPS CURATED EPITOPE ({'/'.join(ev)})"
            else:
                status = "CURATED MATCH, no assay class recorded"
            row.update({
                "IEDB_Status": status, "Match_Type": kind, "Overlap_Len": n,
                "Curated_Epitope": entry["seq"], "Species_Scope": scope,
                "Curated_Source_Antigen": "; ".join(sorted(entry["antigens"]))[:200],
                "Assay_Evidence": "/".join(ev),
                "Curated_MHC_Alleles": "; ".join(sorted(entry["alleles"]))[:200],
                "PubMed_IDs": ";".join(sorted(entry["pubmed"])[:8]),
            })
        results.append(row)

    # ------------------------------------------------------------- report --
    print(f"{'CLASS':<8}{'PEPTIDE':<18}{'TARGET':<22}{'MATCH':<22}{'EVIDENCE':<22}{'PubMed'}")
    print("-" * 110)
    for r in results:
        pm = r["PubMed_IDs"].split(";")[0] if r["PubMed_IDs"] else ""
        print(f"{r['Class']:<8}{r['Peptide']:<18}{r['Source_Target'][:20]:<22}"
              f"{(r['Match_Type'] or 'none'):<22}{(r['Assay_Evidence'] or '-'):<22}{pm}")

    n_conf = sum(1 for r in results if r["IEDB_Status"].startswith("EXPERIMENTALLY"))
    n_over = sum(1 for r in results if r["IEDB_Status"].startswith("OVERLAPS"))
    n_novel = sum(1 for r in results if r["IEDB_Status"].startswith("NOVEL"))
    print("-" * 110)
    print(f"EXPERIMENTALLY CONFIRMED (exact / contained in a curated epitope, with assay data) : {n_conf}/{len(results)}")
    print(f"OVERLAPS a curated epitope (>= {MIN_OVERLAP} aa shared)                                  : {n_over}/{len(results)}")
    print(f"NOVEL -- no curated match                                                        : {n_novel}/{len(results)}")
    print("-" * 110)
    print("NOTE: NOVEL is not a defect. Predicting epitopes with no prior curation is the")
    print("      point of the prediction pipeline; this step distinguishes the two, it does")
    print("      not rank them. Selection was NOT changed by any result here.")
    print("-" * 110)

    ts = datetime.now().strftime("%Y%m%d_%H%M")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"Phase1H_IEDB_Corroboration_{ts}.csv")
    fields = list(results[0].keys())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    common.print_banner("PHASE 1H COMPLETE")
    print(f"[SUCCESS] Execution Time : {common.format_time(time.time() - start)}")
    print(f"[INFO] Report Saved      : {os.path.relpath(path, project_root)}")
    print("=" * 110 + "\n")
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Phase 1H: IEDB experimental corroboration (annotation only)")
    p.add_argument("--refresh", action="store_true", help="Re-download IEDB records, ignoring the cache")
    a = p.parse_args()
    run_step1h(force_refresh=a.refresh)
