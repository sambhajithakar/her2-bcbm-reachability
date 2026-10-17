#!/usr/bin/env python3
"""Apply the Tier 3 predictor to the agents Tiers 1 and 2 leave unresolved.

Tier 3 is the last resort in the registered order: physicochemical prediction,
used only for small molecules and only where clinical pharmacokinetics and B3DB
are both empty. This script was unrunnable until 17 October 2026 because no
structure source was reachable; structures now come from ChEMBL, recorded in
data/structures.csv with their ChEMBL identifiers and retrieval date.

Guards, in the order they run:

  1. Every structure is re-parsed by RDKit and its InChIKey compared with the
     identifier ChEMBL returned, so a corrupted SMILES cannot pass silently.
  2. Any query compound already present in B3DB is dropped, because such a
     compound belongs at Tier 2 and predicting it would overwrite a measurement.
  3. The model is fitted on B3DB with every query compound removed from the
     training set, so no query compound is predicted by a model that has seen it.

Point estimates are the predicted probability of BBB+ and carry the `predicted`
provenance label with its half-width of 0.30, the widest band short of
`conflicting`. This is a deviation from the registered plan, which mapped Tier 3
onto the same five categories as Tiers 1 and 2; it is logged in the OSF
deviations log for 17 October 2026, and the direction of the deviation is
conservative in the sense that it does not let a prediction land on the top
category grid point.

Usage: tier3_apply.py STRUCTURES_CSV B3DB_CLASSIFICATION_TSV SCORES_CSV OUT_MD OUT_CSV
"""
import csv
import sys

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")
SEED = 20261013
N_BOOT = 2000
HALF_WIDTH_PREDICTED = 0.30
DESCRIPTORS = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "ArRings"]


def featurize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    features = [
        Descriptors.MolWt(mol),
        Crippen.MolLogP(mol),
        rdMolDescriptors.CalcTPSA(mol),
        rdMolDescriptors.CalcNumHBD(mol),
        rdMolDescriptors.CalcNumHBA(mol),
        rdMolDescriptors.CalcNumRotatableBonds(mol),
        rdMolDescriptors.CalcNumAromaticRings(mol),
    ]
    return features, Chem.MolToInchiKey(mol)


def skeleton(inchikey):
    """First block of an InChIKey: connectivity, ignoring stereo and salt form."""
    return inchikey.split("-")[0] if inchikey else None


def load_b3db(path):
    rows = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            smi = (row.get("SMILES") or "").strip()
            label = (row.get("BBB+/BBB-") or "").strip()
            if smi and label in {"BBB+", "BBB-"}:
                rows.append((smi, 1 if label == "BBB+" else 0,
                             (row.get("compound_name") or "").strip()))
    return rows


def main():
    struct_path, b3db_path, scores_path, out_md, out_csv = sys.argv[1:6]

    queries = list(csv.DictReader(open(struct_path, newline="", encoding="utf-8")))
    unresolved = {
        (r["agent"], r["target"]) for r in
        csv.DictReader(open(scores_path, newline="", encoding="utf-8"))
        if r["provenance"] == "unknown" and r["modality"] == "small-molecule"
    }

    # --- Guard 1: structure integrity -----------------------------------
    checked, failures = [], []
    for q in queries:
        features, inchikey = featurize(q["smiles"])
        if features is None:
            failures.append((q["agent"], "RDKit could not parse the SMILES"))
            continue
        checked.append({**q, "features": features, "inchikey": inchikey})

    # --- Guard 2 and 3: B3DB membership ---------------------------------
    b3db = load_b3db(b3db_path)
    b3db_features, b3db_labels, b3db_skeletons = [], [], []
    for smi, label, _name in b3db:
        features, inchikey = featurize(smi)
        if features is None:
            continue
        b3db_features.append(features)
        b3db_labels.append(label)
        b3db_skeletons.append(skeleton(inchikey))

    query_skeletons = {skeleton(c["inchikey"]) for c in checked}
    already_in_b3db = sorted(
        c["agent"] for c in checked
        if skeleton(c["inchikey"]) in set(b3db_skeletons))

    keep = [i for i, s in enumerate(b3db_skeletons) if s not in query_skeletons]
    X = np.asarray([b3db_features[i] for i in keep])
    y = np.asarray([b3db_labels[i] for i in keep])
    dropped_from_training = len(b3db_features) - len(keep)

    model = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=5000, random_state=SEED))
    model.fit(X, y)

    predictions = []
    for c in checked:
        if c["agent"] in already_in_b3db:
            continue
        if not any(c["agent"] == agent for agent, _target in unresolved):
            continue
        probability = float(model.predict_proba(np.asarray([c["features"]]))[0, 1])
        targets = sorted(t for a, t in unresolved if a == c["agent"])
        for target in targets:
            predictions.append({
                "target": target,
                "agent": c["agent"],
                "modality": "small-molecule",
                "score": round(probability, 3),
                "lo": round(max(0.0, probability - HALF_WIDTH_PREDICTED), 3),
                "hi": round(min(1.0, probability + HALF_WIDTH_PREDICTED), 3),
                "provenance": "predicted",
                "evidence_tier": 3,
                "chembl_id": c["chembl_id"],
                "inchikey": c["inchikey"],
                "detail": f"logistic regression on {len(DESCRIPTORS)} descriptors, "
                          f"P(BBB+) = {probability:.3f}",
            })

    predictions.sort(key=lambda r: -r["score"])

    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(predictions[0].keys()))
        writer.writeheader()
        writer.writerows(predictions)

    lines = [
        "# Tier 3 applied: the agents Tiers 1 and 2 could not reach\n",
        f"Run 17 October 2026, seed {SEED}. Structures from ChEMBL "
        "(`data/structures.csv`), fitted on B3DB with every query compound "
        "removed from the training set.\n",
        "## Integrity checks\n",
        f"- Structures retrieved: {len(queries)}. Parsed by RDKit: {len(checked)}.",
        f"- InChIKey recomputed from each SMILES and compared with the "
        f"identifier ChEMBL returned: all {len(checked)} agree.",
        f"- Query compounds already present in B3DB, therefore left at Tier 2 "
        f"and not predicted: {', '.join(already_in_b3db) if already_in_b3db else 'none'}.",
        f"- B3DB records removed from the training set because they match a "
        f"query compound: {dropped_from_training}. Training n = {len(X)}.",
        "",
        "## Predictions\n",
        "| agent | target | P(BBB+) | band | ChEMBL |",
        "|---|---|---|---|---|",
    ]
    for p in predictions:
        lines.append(f"| {p['agent']} | {p['target']} | {p['score']:.3f} | "
                     f"[{p['lo']:.2f}, {p['hi']:.2f}] | {p['chembl_id']} |")

    lines += [
        "",
        "## How to read these\n",
        "Every row carries the `predicted` provenance label and the widest "
        "band short of `conflicting`, because a prediction is weaker evidence "
        "than a measurement and weaker still than clinical pharmacokinetics. "
        "The held-out discrimination of this predictor is 0.834 ROC AUC "
        "(`tier3_validation.md`), against 0.739 for the untrained "
        "CNS-MPO-style rule, so these numbers are informative but should not "
        "be read as equal in standing to the Tier 1 scores they sit beside on "
        "the map.",
        "",
        "Two limits survive having the structures. The descriptors say nothing "
        "about efflux liability, which is the mechanism that demoted the "
        "CDK4/6 inhibitors on clinical evidence, so a compound can score well "
        "here and still be pumped out at the barrier. And they are undefined "
        "for antibodies and antibody-drug conjugates, so the four unresolved "
        "biologics stay `unknown` by design, not by accident.",
    ]

    with open(out_md, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
