import os
import sys
import csv
import time
from datetime import datetime
from Bio import Entrez
from Bio import SeqIO

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
# EXPERIMENTAL CONFIGURATION
# =============================================================================

# NCBI Identification (Required for API compliance)
# Uses an environment variable for safety/portability, falling back to your email
Entrez.email = os.getenv("NCBI_EMAIL", "enzoleonor.3309@gmail.com")

# Target viral proteomes for Chimeric Vaccine construct.
# HIV-1 CRF01_AE recombinant-form tagging is reliable in NCBI annotation,
# so those four queries are unchanged. The three Mpox queries are two-tier:
# CONFIRMED_QUERY asks for records whose own annotation mentions "clade IIa"
# or "clade IIb" text (verified against real GenBank /note fields -- see
# Clade tagging note below); FALLBACK_QUERY drops that requirement so a gene
# with sparse clade annotation (B5R, A35R) still yields sequences instead of
# zero results, at the cost of an honestly-reported "Unconfirmed" clade tag.
MPOX_TARGETS = ["L1R", "B5R", "A35R"]
HIV_TARGETS = {
    "HIV_gp120": 'gp120 AND "HIV-1"[Organism] AND CRF01_AE',
    "HIV_gp41": 'gp41 AND "HIV-1"[Organism] AND CRF01_AE',
    "HIV_p24": 'p24 AND "HIV-1"[Organism] AND CRF01_AE',
    "HIV_p17": 'p17 AND "HIV-1"[Organism] AND CRF01_AE',
}

# Sampling parameters
VARIANTS_PER_TARGET = 30
TOTAL_TARGETS = len(MPOX_TARGETS) + len(HIV_TARGETS)
TOTAL_PROJECTED = TOTAL_TARGETS * VARIANTS_PER_TARGET

# =============================================================================
# CLADE TAGGING (Mpox only)
# =============================================================================
# NCBI has no separate taxonomy node for Clade IIa/IIb (Monkeypox virus is a
# single taxid, 10244) -- clade membership is only recoverable from strain-level
# annotation text. Verified empirically: B5R and A35R records in NCBI's protein
# database (25 and 27 total, respectively) carry NO clade annotation at all --
# this is a genuine data-availability gap, not a retrieval bug. Requiring clade
# text in the query (the old behavior) silently zeroed both targets. Instead we
# retrieve confirmed-clade records first, then top up from the unconfirmed pool,
# and record which is which so the paper's clade-purity claim stays auditable.
CLADE_PATTERN_HINTS = ("clade iia", "clade iib", "clade ii", "west african", "genotype: wa", "genotype:wa")


def mpox_queries(gene):
    field_query = f'("{gene}"[Protein Name] OR "{gene}"[Title]) AND "Monkeypox virus"[Organism]'
    confirmed_query = f'{field_query} AND ("clade IIa" OR "clade IIb")'
    return confirmed_query, field_query


def classify_clade(gb_record):
    """Inspects a parsed GenBank record's source feature + description for
    clade signal text. Returns (tag, confirmed: bool)."""
    haystacks = [gb_record.description.lower()]
    for feature in gb_record.features:
        if feature.type == "source":
            for key in ("note", "isolate", "strain"):
                for val in feature.qualifiers.get(key, []):
                    haystacks.append(val.lower())
    blob = " | ".join(haystacks)
    for hint in CLADE_PATTERN_HINTS[:3]:
        if hint in blob:
            return "Clade " + hint.replace("clade ", "").upper(), True
    for hint in CLADE_PATTERN_HINTS[3:]:
        if hint in blob:
            return "Clade II (legacy West African designation)", True
    return "Unconfirmed", False

# =============================================================================
# CORE PROCESSING FUNCTION
# =============================================================================

def fetch_id_pool(query, retmax):
    handle = Entrez.esearch(db="protein", term=query, retmax=retmax)
    result = Entrez.read(handle)
    handle.close()
    return result["IdList"]


def fetch_records(id_list):
    """Fetches full GenBank records (not just FASTA) so source-feature
    qualifiers are available for clade classification."""
    if not id_list:
        return []
    handle = Entrez.efetch(db="protein", id=",".join(id_list), rettype="gb", retmode="text")
    records = list(SeqIO.parse(handle, "genbank"))
    handle.close()
    return records


def run_high_density_retrieval():
    start_time = time.time()

    phase1a_path = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1A")
    os.makedirs(phase1a_path, exist_ok=True)

    print("\n" + "=" * 80)
    print(f"{'PHASE 1A: VIRAL PROTEOME SEQUENCE RETRIEVAL (NCBI ENTREZ)':^80}")
    print("=" * 80)
    print(f"[INFO] Initialization Time : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[INFO] Target Antigens     : {TOTAL_TARGETS} defined")
    print(f"[INFO] Variants per Target : {VARIANTS_PER_TARGET} (ceiling -- see per-target availability below)")
    print(f"[INFO] Projected Yield     : {TOTAL_PROJECTED} FASTA sequences")
    print(f"[INFO] Output Directory    : {phase1a_path}")
    print("-" * 80)

    sys.stdout.write("[PROCESS] Purging prior sequence data to ensure experimental integrity...")
    sys.stdout.flush()
    purged_count = 0
    for f in os.listdir(phase1a_path):
        if f.endswith(".fasta"):
            os.remove(os.path.join(phase1a_path, f))
            purged_count += 1
    print(f" Done. ({purged_count} files removed)")
    print("-" * 80)

    successful_downloads = 0
    clade_audit_rows = []

    # --- Mpox targets: two-tier confirmed/fallback retrieval + clade tagging ---
    for gene in MPOX_TARGETS:
        label = f"Mpox_{gene}"
        print(f"\n[INFO] Establishing NCBI connection for target: {label}")
        try:
            confirmed_query, field_query = mpox_queries(gene)

            confirmed_ids = fetch_id_pool(confirmed_query, VARIANTS_PER_TARGET)
            remaining = VARIANTS_PER_TARGET - len(confirmed_ids)

            fallback_ids = []
            if remaining > 0:
                # Overfetch so we can drop the already-confirmed IDs and still
                # have enough left to fill the remaining slots.
                pool_ids = fetch_id_pool(field_query, VARIANTS_PER_TARGET + remaining + len(confirmed_ids))
                fallback_ids = [i for i in pool_ids if i not in confirmed_ids][:remaining]

            all_ids = confirmed_ids + fallback_ids
            if not all_ids:
                print(f"[WARNING] No records found for target: {label}")
                continue

            elapsed_search = common.format_time(time.time() - start_time)
            sys.stdout.write(f"\r[ FETCH ] {label} | Found {len(all_ids):02d} IDs "
                              f"({len(confirmed_ids)} confirmed-clade, {len(fallback_ids)} unconfirmed) "
                              f"| Downloading... | Elapsed: {elapsed_search}")
            sys.stdout.flush()

            records = fetch_records(all_ids)
            confirmed_id_set = set(confirmed_ids)

            for i, record in enumerate(records):
                clade_tag, is_confirmed = classify_clade(record)
                status = "confirmed" if record.id.split(".")[0] in confirmed_id_set or is_confirmed else "unconfirmed"
                record.description = f"{record.description} | Clade_Tag={clade_tag} ({status})"
                file_name = f"{label}_Var_{i+1:02d}_{record.id}.fasta"
                with open(os.path.join(phase1a_path, file_name), "w") as f:
                    SeqIO.write(record, f, "fasta")
                successful_downloads += 1
                clade_audit_rows.append({
                    "Target": label, "Accession": record.id,
                    "Clade_Tag": clade_tag, "Confirmed": status == "confirmed",
                })

            elapsed_done = common.format_time(time.time() - start_time)
            n_confirmed_saved = sum(1 for r in clade_audit_rows if r["Target"] == label and r["Confirmed"])
            sys.stdout.write(f"\r[ FETCH ] {label} | Saved {len(records):02d} records "
                              f"({n_confirmed_saved} clade-confirmed) | Elapsed: {elapsed_done:<15}\n")
            sys.stdout.flush()
            if len(records) < VARIANTS_PER_TARGET:
                print(f"[NOTICE] {label}: only {len(records)}/{VARIANTS_PER_TARGET} sequences exist in NCBI's "
                      f"protein database for this query -- this is the full available pool, not a retrieval failure.")

            time.sleep(1.0)

        except Exception as e:
            print(f"\n[ERROR] Protocol failure during {label} acquisition. Reason: {e}")

    # --- HIV targets: unchanged single-tier retrieval (CRF01_AE tags reliably) ---
    for label, query in HIV_TARGETS.items():
        print(f"\n[INFO] Establishing NCBI connection for target: {label}")
        try:
            id_list = fetch_id_pool(query, VARIANTS_PER_TARGET)
            if not id_list:
                print(f"[WARNING] No records found for query: '{query}'")
                continue

            elapsed_search = common.format_time(time.time() - start_time)
            sys.stdout.write(f"\r[ FETCH ] {label} | Found {len(id_list):02d} IDs | Downloading batch... "
                              f"| Elapsed: {elapsed_search}")
            sys.stdout.flush()

            fetch_handle = Entrez.efetch(db="protein", id=",".join(id_list), rettype="fasta", retmode="text")
            records = list(SeqIO.parse(fetch_handle, "fasta"))
            fetch_handle.close()

            for i, record in enumerate(records):
                file_name = f"{label}_Var_{i+1:02d}_{record.id}.fasta"
                with open(os.path.join(phase1a_path, file_name), "w") as f:
                    SeqIO.write(record, f, "fasta")
                successful_downloads += 1

            elapsed_done = common.format_time(time.time() - start_time)
            sys.stdout.write(f"\r[ FETCH ] {label} | Successfully saved {len(records):02d} records "
                              f"| Elapsed: {elapsed_done:<15}\n")
            sys.stdout.flush()

            time.sleep(1.0)

        except Exception as e:
            print(f"\n[ERROR] Protocol failure during {label} acquisition. Reason: {e}")

    # --- Clade audit CSV (Mpox only) ---
    if clade_audit_rows:
        audit_path = os.path.join(phase1a_path, "Phase1A_Clade_Audit.csv")
        with open(audit_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["Target", "Accession", "Clade_Tag", "Confirmed"])
            writer.writeheader()
            writer.writerows(clade_audit_rows)
        print(f"\n[INFO] Clade audit written : {audit_path}")

    total_time = common.format_time(time.time() - start_time)
    print("\n" + "=" * 80)
    print(f"{'ACQUISITION PROTOCOL COMPLETE':^80}")
    print("=" * 80)
    print(f"[SUCCESS] Total Sequences Retrieved : {successful_downloads}/{TOTAL_PROJECTED}")
    print(f"[SUCCESS] Total Execution Time      : {total_time}")
    print("[INFO] Data formatting complete. Proceed to Phase 1B for stability thresholding.")
    print("=" * 80 + "\n")

if __name__ == "__main__":
    run_high_density_retrieval()
