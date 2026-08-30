"""
Local SEMA-3D inference -- replaces the WAF-blocked sema.airi.net API.

WHY THIS EXISTS
    Sec. I.D specifies SEMA 2.0 for conformational B-cell corroboration. The
    hosted service's prediction endpoints sit behind a ServicePipe anti-bot WAF
    (deviation #20), so no data was ever produced. SEMA is open source
    (github.com/AIRI-Institute/SEMAi), so the model runs locally instead and the
    blocker disappears. This module is a faithful, device-agnostic port of the
    project's own SEMA-3D_inference.ipynb -- the upstream notebook hardcodes
    .cuda(), which no Apple-silicon machine can run.

WHAT SEMA-3D IS
    An ensemble of five SaProt-650M models fine-tuned to regress, per residue,
    the LOG-SCALED EXPECTED NUMBER OF ANTIBODY CONTACTS. SaProt consumes a
    "combined sequence": each residue's amino acid interleaved with its Foldseek
    3Di structural token, which is what makes it structure-aware.

    THE OUTPUT IS NOT A PROBABILITY. It is an unbounded regression score on a
    log-contact scale, and it is routinely negative. Any threshold of the form
    "score >= 0.5 means epitope" is meaningless against it -- see
    calibrate_threshold() and the note in Phase1De_conformational.py.

REQUIREMENTS
    - foldseek on PATH (macOS/arm64: `brew install brewsci/bio/foldseek`). The
      binary vendored in SEMAi/saprot_utils/bin is a Linux x86-64 ELF and will
      not execute here.
    - The five fine-tuned checkpoints (sema_3d_0..4.pth) and the SaProt-650M
      base model. See TOOLS_INSTALL.txt.
"""

import os
import shutil
import subprocess
import tempfile

import numpy as np
import torch
from torch import nn
from transformers import EsmTokenizer, EsmForMaskedLM

SAPROT_MODEL = "westlake-repl/SaProt_650M_PDB"
N_ENSEMBLE = 5

# The fine-tuned head is Linear(446, 2) over the MLM logit space -- 446 is
# SaProt's combined (amino-acid x 3Di) vocabulary size. Taken from the
# upstream notebook; changing it silently breaks checkpoint loading.
_CLASSIFIER_IN = 446
_CLASSIFIER_OUT = 2


def _foldseek_binary():
    """Prefers a real host binary; the vendored one is Linux-only."""
    found = shutil.which("foldseek")
    if found:
        return found
    raise RuntimeError(
        "foldseek not found on PATH. Install it (macOS: `brew install "
        "brewsci/bio/foldseek`). The copy in SEMAi/saprot_utils/bin is a Linux "
        "x86-64 binary and cannot run on this machine.")


PLDDT_MASK_THRESHOLD = 70.0


def get_struc_seq(pdb_path, chains=None, plddt=None, plddt_threshold=PLDDT_MASK_THRESHOLD):
    """
    Runs Foldseek's structureto3didescriptor and returns
    {chain: (aa_seq, tridi_seq, combined_seq)}.

    Ported from SEMAi/saprot_utils/foldseek_util.py, with three changes: the
    temporary files are written to a private temp directory rather than the
    current working directory (the original leaks them into wherever it was
    invoked), failures raise instead of being swallowed by os.system, and the
    pLDDT mask takes per-residue values directly instead of a JSON path.

    pLDDT MASKING -- why it is on by default for predicted structures.
    SaProt reads a 3Di structural token per residue. In a low-confidence region
    that token is not "unknown", it is a CONFIDENTLY WRONG description of a
    fold AlphaFold could not determine, and SaProt has no way to tell the
    difference. Upstream SEMA masks those positions with "#" for exactly this
    reason. It matters here: HIV_gp41's model has pTM 0.40 with only 38% of
    residues at pLDDT >= 70, because of its transmembrane and cytoplasmic
    tail -- while the B-cell epitope region itself sits at pLDDT 71.3. Masking
    keeps the trustworthy part and stops the disordered tail from contributing
    invented structure. Deviation #21 is the record of what ignoring fold
    confidence costs.
    """
    foldseek = _foldseek_binary()
    pdb_path = os.path.abspath(pdb_path)
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(pdb_path)

    with tempfile.TemporaryDirectory() as tmpdir:
        out_tsv = os.path.join(tmpdir, "struc_seq.tsv")
        proc = subprocess.run(
            [foldseek, "structureto3didescriptor", "-v", "0", "--threads", "1",
             "--chain-name-mode", "1", pdb_path, out_tsv],
            capture_output=True, text=True)
        if proc.returncode != 0 or not os.path.isfile(out_tsv):
            raise RuntimeError(
                f"foldseek failed on {pdb_path}:\n{proc.stdout}\n{proc.stderr}")

        name = os.path.basename(pdb_path)
        seq_dict = {}
        with open(out_tsv) as fh:
            for line in fh:
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                desc, seq, tridi = parts[0], parts[1], parts[2]
                name_chain = desc.split(" ")[0]
                chain = name_chain.replace(name, "").split("_")[-1]
                if chains is not None and chain not in chains:
                    continue
                if chain in seq_dict:
                    continue
                if plddt is not None:
                    if len(plddt) != len(tridi):
                        raise ValueError(
                            f"pLDDT has {len(plddt)} values for {len(tridi)} modelled "
                            f"residues in {pdb_path}: refusing to mask on a mismatched index.")
                    tridi = "".join("#" if p < plddt_threshold else c
                                    for c, p in zip(tridi, plddt))
                combined = "".join(a + b.lower() for a, b in zip(seq, tridi))
                seq_dict[chain] = (seq, tridi, combined)
    if not seq_dict:
        raise RuntimeError(f"foldseek returned no chains for {pdb_path}")
    return seq_dict


class SaProtForTokenClassification(nn.Module):
    """Exactly the upstream architecture -- the checkpoints depend on it."""

    def __init__(self, num_labels=_CLASSIFIER_OUT):
        super().__init__()
        self.num_labels = num_labels
        self.encoder = EsmForMaskedLM.from_pretrained(SAPROT_MODEL)
        self.classifier = nn.Linear(_CLASSIFIER_IN, self.num_labels)

    def forward(self, token_ids):
        logits = self.encoder(input_ids=token_ids)["logits"]
        logits = logits[:, 1:-1, :]      # drop <cls> and <eos>
        return self.classifier(logits)


def pick_device(prefer="auto"):
    """
    MPS is used when available. It is not merely a speed choice: these are
    650M-parameter models and CPU inference on a 500-residue chain is minutes
    per checkpoint, times five checkpoints, times seven antigens.
    """
    if prefer != "auto":
        return torch.device(prefer)
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def load_ensemble(models_dir, device, n=N_ENSEMBLE):
    models = []
    for seed in range(n):
        path = os.path.join(models_dir, f"sema_3d_{seed}.pth")
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"Missing SEMA-3D checkpoint: {path}. See TOOLS_INSTALL.txt.")
        model = SaProtForTokenClassification()
        state = torch.load(path, map_location="cpu", weights_only=False)
        model.load_state_dict(state)
        model.eval().to(device)
        models.append(model)
    return models


@torch.no_grad()
def predict_combined_seq(models, combined_seq, device, tokenizer=None):
    """
    Returns the ensemble-mean per-residue score for one combined sequence.
    Score = logit channel 1, matching the upstream notebook's
    `model.forward(...)[0][0][:, 1]`.
    """
    tokenizer = tokenizer or EsmTokenizer.from_pretrained(SAPROT_MODEL)
    token_ids = tokenizer(combined_seq, padding=False,
                          return_tensors="pt")["input_ids"].to(device)
    preds = []
    for model in models:
        out = model(token_ids)[0][:, 1]
        preds.append(out.float().cpu().numpy())
    return np.mean(np.stack(preds), axis=0)


def predict_pdb(models, pdb_path, chain, device, tokenizer=None):
    """
    Returns (aa_seq, per_residue_scores) for one chain of one PDB.

    The returned sequence is Foldseek's view of the MODELLED residues, which is
    what the scores are indexed to. Callers must align it to their reference
    sequence rather than assuming the two share a numbering -- a structure with
    unmodelled loops silently shifts every downstream offset otherwise.
    """
    chains = get_struc_seq(pdb_path)
    key = None
    for candidate in (chain, chain.upper(), chain.lower()):
        if candidate in chains:
            key = candidate
            break
    if key is None:
        raise KeyError(
            f"Chain {chain!r} not in {pdb_path} (found: {sorted(chains)})")
    aa_seq, _tridi, combined = chains[key]
    scores = predict_combined_seq(models, combined, device, tokenizer)
    if len(scores) != len(aa_seq):
        raise RuntimeError(
            f"SEMA-3D returned {len(scores)} scores for {len(aa_seq)} residues "
            f"in {pdb_path}:{key}")
    return aa_seq, scores

# =============================================================================
# STREAMING ENSEMBLE -- the memory-safe way to run all five checkpoints.
#
# MEASURED ON THIS MACHINE: each fine-tuned SaProt-650M occupies 2.43 GB of the
# Apple unified memory pool, exactly linearly (1 model 2.43 GB, 2 models 4.86
# GB), on top of ~5 GB of CPU-side peak RSS for the load itself. Five resident
# at once is 12.2 GB of a 16 GB pool -- it does not fit, and load_ensemble(n=5)
# would thrash or die partway through a long run.
#
# The fix inverts the loops. Foldseek preprocessing needs no model at all, so
# every antigen's combined sequence is computed FIRST. Then the checkpoint loop
# is the OUTER one: load seed k, score every antigen with it, accumulate into a
# running sum, free it, move on. Peak memory is one model (2.43 GB) regardless
# of how many antigens or checkpoints there are, and each checkpoint is still
# read from disk exactly once -- the same 5 loads load_ensemble would do, not
# 5 x n_antigens.
#
# The arithmetic is identical to the upstream notebook's
# np.mean(np.stack(ensembl_res), axis=0): a mean of five per-residue score
# vectors. Only the order of accumulation differs.
# =============================================================================

def predict_many(models_dir, jobs, device, n=N_ENSEMBLE, tokenizer=None, verbose=True):
    """
    Scores several structures with the full ensemble, holding ONE model in
    memory at a time.

        jobs: {name: (pdb_path, chain)} or {name: (pdb_path, chain, plddt_list)}
        -> {name: (aa_seq, ensemble_mean_scores)}

    Raises before loading anything if a checkpoint is missing, so a long run
    cannot die half way through having already spent the foldseek work.
    """
    for seed in range(n):
        path = os.path.join(models_dir, f"sema_3d_{seed}.pth")
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"Missing SEMA-3D checkpoint: {path}. See TOOLS_INSTALL.txt (ENV 5).")

    tokenizer = tokenizer or EsmTokenizer.from_pretrained(SAPROT_MODEL)

    # ---- 1. foldseek pass: no model needed, so do it once up front ----------
    prepared = {}
    for name, spec in jobs.items():
        pdb_path, chain = spec[0], spec[1]
        job_plddt = spec[2] if len(spec) > 2 else None
        chains = get_struc_seq(pdb_path, plddt=job_plddt)
        key = next((c for c in (chain, chain.upper(), chain.lower()) if c in chains), None)
        if key is None:
            raise KeyError(f"Chain {chain!r} not in {pdb_path} (found: {sorted(chains)})")
        aa_seq, tridi, combined = chains[key]
        prepared[name] = (aa_seq, combined)
        if verbose:
            masked = tridi.count("#")
            note = (f", {masked} masked (pLDDT < {PLDDT_MASK_THRESHOLD:.0f})"
                    if job_plddt is not None else "")
            print(f"       [foldseek] {name}: {len(aa_seq)} residues{note}")

    totals = {name: None for name in jobs}

    # ---- 2. checkpoint loop on the OUTSIDE ----------------------------------
    for seed in range(n):
        if verbose:
            print(f"       [ensemble] checkpoint {seed + 1}/{n} ...")
        model = SaProtForTokenClassification()
        state = torch.load(os.path.join(models_dir, f"sema_3d_{seed}.pth"),
                           map_location="cpu", weights_only=False)
        model.load_state_dict(state)
        del state
        model.eval().to(device)

        with torch.no_grad():
            for name, (_aa, combined) in prepared.items():
                token_ids = tokenizer(combined, padding=False,
                                      return_tensors="pt")["input_ids"].to(device)
                out = model(token_ids)[0][:, 1].float().cpu().numpy()
                totals[name] = out if totals[name] is None else totals[name] + out

        # Free before the next checkpoint -- without the cache empty, MPS holds
        # the previous model's blocks and the peak is 2x, defeating the point.
        del model
        if device.type == "mps":
            torch.mps.empty_cache()
        elif device.type == "cuda":
            torch.cuda.empty_cache()

    results = {}
    for name, (aa_seq, _combined) in prepared.items():
        scores = totals[name] / float(n)
        if len(scores) != len(aa_seq):
            raise RuntimeError(
                f"SEMA-3D returned {len(scores)} scores for {len(aa_seq)} residues ({name})")
        results[name] = (aa_seq, scores)
    return results
