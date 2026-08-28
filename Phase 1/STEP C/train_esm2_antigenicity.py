import os
import inspect
import torch
import numpy as np
import pandas as pd
from datetime import datetime
from sklearn.metrics import roc_auc_score, confusion_matrix
from sklearn.model_selection import train_test_split
from transformers import AutoTokenizer, EsmForSequenceClassification, Trainer, TrainingArguments, EvalPrediction
from datasets import Dataset, DatasetDict

# =============================================================================
# OFFLINE TRAINING JOB - run this once, manually, before Phase1C_antigenicity.py
# is ever executed as part of the pipeline. This script is intentionally NOT
# named Phase1C_antigenicity.py: that name is reserved for the per-run
# inference/screening script that reads Phase1B output and writes
# Filtered_Antigenicity. Overwriting that file with this one would silently
# remove antigenicity screening from the pipeline.
# =============================================================================

MODEL_NAME = "facebook/esm2_t12_35M_UR50D"
DATASET_PATH = "data/benchmark_dataset.csv"  # Format: sequence, label (1/0), source_id
OUTPUT_DIR = "models/esm2_antigenicity_finetuned"
MAX_LENGTH = 1024


def print_banner(text):
    print(f"\n{'='*80}\n{text:^80}\n{'='*80}")


def make_training_args(**kwargs):
    """
    transformers renamed TrainingArguments' `evaluation_strategy` to
    `eval_strategy` in v4.46+, and later versions (5.x, confirmed against
    the installed 5.15.1) dropped `logging_dir` entirely. Detect what the
    installed version actually accepts so this script doesn't hard-fail on
    a version mismatch.
    """
    sig = inspect.signature(TrainingArguments.__init__)
    strategy_key = "eval_strategy" if "eval_strategy" in sig.parameters else "evaluation_strategy"
    kwargs[strategy_key] = kwargs.pop("_strategy_value")
    if "logging_dir" in kwargs and "logging_dir" not in sig.parameters:
        kwargs.pop("logging_dir")
    return TrainingArguments(**kwargs)


# =============================================================================
# METRIC COMPUTATION & VALIDATION GATES
# =============================================================================
def compute_metrics(p: EvalPrediction):
    """Calculates all 5 metrics required by the Phase 1C Remediation Plan."""
    preds = torch.softmax(torch.tensor(p.predictions), dim=1)[:, 1].numpy()
    labels = p.label_ids

    # Paper threshold: >= 0.50 -> antigenic
    preds_binary = (preds >= 0.50).astype(int)

    auc = roc_auc_score(labels, preds)
    # labels=[0, 1] forces a full 2x2 matrix even if a batch happens to be
    # single-class, avoiding an unpack error on .ravel()
    tn, fp, fn, tp = confusion_matrix(labels, preds_binary, labels=[0, 1]).ravel()

    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0.0

    return {
        "auc": auc,
        "sensitivity": sens,
        "specificity": spec,
        "ppv": ppv,
        "npv": npv,
    }


def run_finetuning():
    print_banner("ESM-2 ANTIGENICITY FINE-TUNING & VALIDATION")

    # 1. LOAD & PREPARE DATA (Step 2)
    if not os.path.exists(DATASET_PATH):
        raise FileNotFoundError(
            f"Benchmark dataset missing at {DATASET_PATH}. "
            f"Please compile 500 IEDB antigenic / 500 non-antigenic (e.g. "
            f"UniProt housekeeping) sequences first."
        )

    df = pd.read_csv(DATASET_PATH)
    print(f"[INFO] Loaded {len(df)} sequences. Class balance:\n{df['label'].value_counts()}")

    # Stratified 70/15/15 split
    train_df, temp_df = train_test_split(df, test_size=0.30, stratify=df["label"], random_state=42)
    val_df, test_df = train_test_split(temp_df, test_size=0.50, stratify=temp_df["label"], random_state=42)

    print(f"[INFO] Split sizes - Train: {len(train_df)} | Val: {len(val_df)} | Test: {len(test_df)}")

    # 2. TOKENIZATION (Step 3)
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)

    def tokenize_data(data):
        dataset = Dataset.from_pandas(data)
        return dataset.map(
            lambda x: tokenizer(x["sequence"], padding="max_length", truncation=True, max_length=MAX_LENGTH),
            batched=True,
        )

    ds = DatasetDict({
        "train": tokenize_data(train_df),
        "val": tokenize_data(val_df),
        "test": tokenize_data(test_df),
    })

    # 3. MODEL SETUP & TRAINING (Step 3)
    model = EsmForSequenceClassification.from_pretrained(MODEL_NAME, num_labels=2)

    training_args = make_training_args(
        output_dir="./training_checkpoints",
        learning_rate=2e-5,
        per_device_train_batch_size=8,
        per_device_eval_batch_size=8,
        num_train_epochs=5,
        _strategy_value="epoch",
        save_strategy="epoch",
        load_best_model_at_end=True,
        metric_for_best_model="auc",
        logging_dir="./logs",
    )

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=ds["train"],
        eval_dataset=ds["val"],
        compute_metrics=compute_metrics,
    )

    print("\n[PROCESS] Initiating training protocol...")
    trainer.train()

    # 4. STRICT TEST SET VALIDATION (Step 4)
    print_banner("TEST SET VALIDATION & PASS/FAIL GATE")
    test_results = trainer.evaluate(ds["test"])

    m_auc = test_results["eval_auc"]
    m_sens = test_results["eval_sensitivity"]
    m_spec = test_results["eval_specificity"]
    m_ppv = test_results["eval_ppv"]
    m_npv = test_results["eval_npv"]

    print(f"[TEST METRICS] AUC: {m_auc:.3f} | Sens: {m_sens:.3f} | Spec: {m_spec:.3f} | PPV: {m_ppv:.3f} | NPV: {m_npv:.3f}")

    # The explicit pass/fail gate — do not relax these to force a pass
    passed = (m_auc > 0.85) and (m_sens >= 0.85) and (m_spec >= 0.80) and (m_ppv >= 0.85) and (m_npv >= 0.80)

    if not passed:
        raise RuntimeError(
            "\n[BLOCKED] Model failed to meet the paper's claimed thresholds on the test set.\n"
            "Do not ship this model. Adjust hyperparameters, unfreeze layers, or improve dataset quality."
        )

    # 5. ARTIFACT GENERATION (Step 5)
    print("\n[SUCCESS] Model passed all validation gates. Saving artifacts...")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    model.save_pretrained(OUTPUT_DIR)
    tokenizer.save_pretrained(OUTPUT_DIR)

    report_path = os.path.join(OUTPUT_DIR, "VALIDATION_REPORT.md")
    with open(report_path, "w") as f:
        f.write("# ESM-2 Antigenicity Finetuning Validation Report\n")
        f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("## Dataset Composition\n")
        f.write(f"- Total Sequences: {len(df)}\n")
        f.write(f"- Train Split: {len(train_df)}\n")
        f.write(f"- Validation Split: {len(val_df)}\n")
        f.write(f"- Held-Out Test Split: {len(test_df)}\n\n")
        f.write("## Achieved Test Metrics (Threshold = 0.50)\n")
        f.write(f"- **AUC:** {m_auc:.4f} (Required: > 0.85)\n")
        f.write(f"- **Sensitivity:** {m_sens:.4f} (Required: >= 0.85)\n")
        f.write(f"- **Specificity:** {m_spec:.4f} (Required: >= 0.80)\n")
        f.write(f"- **PPV:** {m_ppv:.4f} (Required: >= 0.85)\n")
        f.write(f"- **NPV:** {m_npv:.4f} (Required: >= 0.80)\n\n")
        f.write("> **Status:** CLEAR. This model meets all parameters claimed in Section II.C.I.C.\n")

    print(f"[INFO] Validation report written to {report_path}")
    print("[INFO] Phase 1C inference script is now unblocked.")


if __name__ == "__main__":
    run_finetuning()
