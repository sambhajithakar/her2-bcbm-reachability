#!/usr/bin/env python3
"""Tier 3: physicochemical prediction of barrier permeability, and an honest
measurement of how much it can be trusted.

The protocol allows prediction only where Tiers 1 and 2 are empty, so the
question that matters is not whether a predictor can be fitted but how well it
does on compounds it has not seen. This script fits two predictors on B3DB and
evaluates both on a held-out split:

  1. A transparent CNS-MPO-style desirability rule, which uses no training data.
  2. Logistic regression on the same descriptors.

Both are evaluated on the same held-out compounds, so the comparison says
whether fitting buys anything over the published rule.

Usage: tier3_predictor.py B3DB_CLASSIFICATION_TSV OUT_MD
"""
import csv
import sys

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from sklearn.calibration import calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")
SEED = 20261013
DESCRIPTORS = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "ArRings"]


def featurize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return [
        Descriptors.MolWt(mol),
        Crippen.MolLogP(mol),
        rdMolDescriptors.CalcTPSA(mol),
        rdMolDescriptors.CalcNumHBD(mol),
        rdMolDescriptors.CalcNumHBA(mol),
        rdMolDescriptors.CalcNumRotatableBonds(mol),
        rdMolDescriptors.CalcNumAromaticRings(mol),
    ]


def _ramp(value, good, bad):
    """1.0 at or beyond `good`, 0.0 at or beyond `bad`, linear between."""
    if good < bad:
        if value <= good:
            return 1.0
        if value >= bad:
            return 0.0
        return (bad - value) / (bad - good)
    if value >= good:
        return 1.0
    if value <= bad:
        return 0.0
    return (value - bad) / (good - bad)


def cns_mpo_like(features):
    """CNS-MPO desirability over the six properties available without pKa.

    The published score uses six properties including most basic pKa, which
    needs a predictor of its own; this variant drops it and renormalises, so it
    is a reimplementation in the spirit of the rule and not the rule itself.
    Scored 0-1 rather than 0-6.
    """
    mw, clogp, tpsa, hbd, hba, rotb, ar = features
    desirabilities = [
        _ramp(mw, 360.0, 500.0),
        _ramp(clogp, 3.0, 5.0),
        _ramp(tpsa, 40.0, 120.0) * _ramp(tpsa, 90.0, 20.0),  # favours 40-90
        _ramp(hbd, 0.5, 3.5),
        _ramp(hba, 4.0, 10.0),
        _ramp(rotb, 4.0, 10.0),
    ]
    return float(np.mean(desirabilities))


def load(path):
    smiles, labels, names = [], [], []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            smi = (row.get("SMILES") or "").strip()
            label = (row.get("BBB+/BBB-") or "").strip()
            if smi and label in {"BBB+", "BBB-"}:
                smiles.append(smi)
                labels.append(1 if label == "BBB+" else 0)
                names.append((row.get("compound_name") or "").strip())
    return smiles, labels, names


def query_skeletons(path):
    """InChIKey skeletons of the agents the predictor will be applied to.

    Passed so the held-out evaluation can exclude them, which makes the
    reported discrimination belong to the same model that produces the
    applied predictions rather than to a model fitted on a superset.
    """
    from rdkit.Chem.MolStandardize import rdMolStandardize
    chooser = rdMolStandardize.LargestFragmentChooser()
    uncharger = rdMolStandardize.Uncharger()
    skeletons = set()
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            mol = Chem.MolFromSmiles((row.get("smiles") or "").strip())
            if mol is None:
                continue
            key = Chem.MolToInchiKey(uncharger.uncharge(chooser.choose(mol)))
            if key:
                skeletons.add(key.split("-")[0])
    return skeletons


def skeleton_of(smiles):
    from rdkit.Chem.MolStandardize import rdMolStandardize
    chooser = rdMolStandardize.LargestFragmentChooser()
    uncharger = rdMolStandardize.Uncharger()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    key = Chem.MolToInchiKey(uncharger.uncharge(chooser.choose(mol)))
    return key.split("-")[0] if key else None


def main():
    data_path, out_path = sys.argv[1:3]
    structures_path = sys.argv[3] if len(sys.argv) > 3 else None
    smiles, labels, names = load(data_path)

    removed = 0
    if structures_path:
        queries = query_skeletons(structures_path)
        kept = [i for i, smi in enumerate(smiles)
                if skeleton_of(smi) not in queries]
        removed = len(smiles) - len(kept)
        smiles = [smiles[i] for i in kept]
        labels = [labels[i] for i in kept]
        names = [names[i] for i in kept]

    features, keep = [], []
    for i, smi in enumerate(smiles):
        vector = featurize(smi)
        if vector is not None:
            features.append(vector)
            keep.append(i)
    X = np.asarray(features)
    y = np.asarray([labels[i] for i in keep])
    parsed, total = len(X), len(smiles)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=SEED, stratify=y)

    rule_test = np.array([cns_mpo_like(row) for row in X_test])
    rule_auc = roc_auc_score(y_test, rule_test)

    model = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=5000, random_state=SEED))
    model.fit(X_train, y_train)
    proba = model.predict_proba(X_test)[:, 1]
    model_auc = roc_auc_score(y_test, proba)
    brier = brier_score_loss(y_test, proba)
    accuracy = ((proba >= 0.5).astype(int) == y_test).mean()
    base_rate = y_test.mean()

    # Bootstrap the held-out AUC difference.
    rng = np.random.default_rng(SEED)
    diffs = []
    for _ in range(2000):
        idx = rng.integers(0, len(y_test), len(y_test))
        if len(np.unique(y_test[idx])) < 2:
            continue
        diffs.append(roc_auc_score(y_test[idx], proba[idx])
                     - roc_auc_score(y_test[idx], rule_test[idx]))
    lo, hi = np.quantile(diffs, [0.025, 0.975])

    frac_pos, mean_pred = calibration_curve(y_test, proba, n_bins=8, strategy="quantile")
    coefficients = dict(zip(
        DESCRIPTORS, model.named_steps["logisticregression"].coef_[0].round(3)))

    lines = [
        "# Tier 3 predictor: fitted, and measured on held-out compounds\n",
        f"Source: B3DB classification set. {total} labelled records, "
        f"{parsed} parsed by RDKit ({parsed / total:.1%}). "
        f"Stratified 70/30 split, seed {SEED}. Held-out n = {len(y_test)}, "
        f"BBB+ base rate {base_rate:.1%}."
        + (f" {removed} record(s) whose InChIKey skeleton matches an agent the "
           f"predictor is applied to were removed before the split, so the "
           f"discrimination below belongs to the same model that produces the "
           f"applied predictions."
           if structures_path else
           " No query compounds were withheld; see the v2 run for the applied "
           "model's own discrimination.") + "\n",
        "| predictor | held-out ROC AUC |",
        "|---|---|",
        f"| CNS-MPO-style rule (no training data) | {rule_auc:.3f} |",
        f"| Logistic regression on the same 7 descriptors | {model_auc:.3f} |",
        "",
        f"Difference in AUC (regression minus rule), bootstrap 95% CI: "
        f"{model_auc - rule_auc:+.3f} [{lo:+.3f}, {hi:+.3f}].",
        "",
        f"Logistic regression also: accuracy {accuracy:.3f} at a 0.5 threshold, "
        f"Brier score {brier:.3f}.",
        "",
        "Standardised coefficients (positive favours BBB+): "
        + ", ".join(f"{k} {v:+.3f}" for k, v in coefficients.items()) + ".",
        "",
        "Calibration on the held-out split, by quantile bin:",
        "",
        "| mean predicted | observed BBB+ fraction |",
        "|---|---|",
    ]
    for predicted, observed in zip(mean_pred, frac_pos):
        lines.append(f"| {predicted:.3f} | {observed:.3f} |")

    lines += [
        "",
        "## What this predictor is, and what it is not\n",
        "Tier 3 is the last resort in the registered order: physicochemical "
        "prediction, used only for small molecules and only where clinical "
        "pharmacokinetics and B3DB are both empty. The agents it is applied to, "
        "and their predicted probabilities, are in `tier3_applied.md`.",
        "",
        "Two limits survive having the structures. The descriptors are "
        "tumour-agnostic and say nothing about efflux liability, which is the "
        "mechanism that demoted the CDK4/6 inhibitors on clinical evidence, so "
        "a compound can score well here and still be pumped out at the barrier. "
        "And they are undefined for antibodies and antibody-drug conjugates, so "
        "the four unresolved biologics are outside the predictor's domain "
        "regardless.",
    ]

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
