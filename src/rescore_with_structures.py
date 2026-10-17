#!/usr/bin/env python3
"""Re-score the reachability axis now that agent structures are retrievable.

Three things change relative to the 13 October run, and only these three:

  1. **B3DB is matched on structure, not on compound name.** Name matching
     missed records whose `compound_name` field carries a synonym list, which is
     how everolimus ("everolimus, rad001") escaped the Tier 2 lookup. Matching
     is on the first block of the InChIKey after salt stripping and charge
     neutralisation, so a salt, a counter-ion or a protonation state cannot
     hide a measurement.
  2. **Tier 3 can be applied**, because structures now exist for the small
     molecules Tiers 1 and 2 leave empty.
  3. Nothing else. Tier 1 keeps precedence over both, the category grid, band
     widths and thresholds are unchanged, and no antibody or ADC is predicted.

Every change of score relative to the name-matched run is listed in the report,
so the effect of the matching fix is visible rather than folded into the result.

Usage: rescore_with_structures.py AGENTS SCORES_V1 STRUCTURES B3DB_REG B3DB_CLS OUT_CSV OUT_MD
"""
import csv
import sys

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")
SEED = 20261013
CATEGORY = {"High": 0.80, "Medium-High": 0.65, "Medium": 0.50,
            "Low-Medium": 0.35, "Low": 0.20}
HALF_WIDTH = {"clinical-RCT": 0.08, "clinical-PK": 0.10,
              "in-vivo-intracranial": 0.12, "barrier-model": 0.18,
              "measured-logBB": 0.15, "classified-BBB": 0.20,
              "inferred-by-class": 0.25, "predicted": 0.30,
              "conflicting": 0.45}
LOGBB_HIGH, LOGBB_LOW = 0.30, -0.50

_chooser = rdMolStandardize.LargestFragmentChooser()
_uncharger = rdMolStandardize.Uncharger()


def normalize(smiles):
    """Largest fragment, charges neutralised, then the InChIKey skeleton."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    mol = _uncharger.uncharge(_chooser.choose(mol))
    key = Chem.MolToInchiKey(mol)
    return mol, (key.split("-")[0] if key else None)


def featurize(mol):
    return [
        Descriptors.MolWt(mol),
        Crippen.MolLogP(mol),
        rdMolDescriptors.CalcTPSA(mol),
        rdMolDescriptors.CalcNumHBD(mol),
        rdMolDescriptors.CalcNumHBA(mol),
        rdMolDescriptors.CalcNumRotatableBonds(mol),
        rdMolDescriptors.CalcNumAromaticRings(mol),
    ]


def index_b3db(reg_path, cls_path):
    """skeleton -> {"logbb": [...], "labels": [...], "names": {...}}"""
    index, training = {}, []
    for path, kind in ((reg_path, "logbb"), (cls_path, "label")):
        with open(path, newline="", encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                smiles = (row.get("SMILES") or "").strip()
                if not smiles:
                    continue
                mol, skeleton = normalize(smiles)
                if skeleton is None:
                    continue
                entry = index.setdefault(
                    skeleton, {"logbb": [], "labels": [], "names": set()})
                entry["names"].add((row.get("compound_name") or "").strip())
                if kind == "logbb":
                    try:
                        entry["logbb"].append(float(row.get("logBB")))
                    except (TypeError, ValueError):
                        pass
                else:
                    label = (row.get("BBB+/BBB-") or "").strip()
                    if label in {"BBB+", "BBB-"}:
                        entry["labels"].append(label)
                        training.append((skeleton, featurize(mol),
                                         1 if label == "BBB+" else 0))
    return index, training


def tier2(entry):
    """Same rules as score_reachability.py, applied to a structural match."""
    logbb = entry["logbb"]
    if logbb:
        lo, hi = min(logbb), max(logbb)
        if hi >= LOGBB_HIGH and lo <= LOGBB_LOW:
            return "Medium", "conflicting", (
                f"B3DB logBB records span {lo:.2f} to {hi:.2f}; not averaged")
        mean = sum(logbb) / len(logbb)
        category = ("High" if mean >= LOGBB_HIGH
                    else "Low" if mean <= LOGBB_LOW else "Medium")
        return category, "measured-logBB", f"B3DB logBB n={len(logbb)}, mean {mean:.2f}"
    labels = set(entry["labels"])
    if labels == {"BBB+", "BBB-"}:
        return "Medium", "conflicting", (
            f"B3DB holds both BBB+ and BBB- (n={len(entry['labels'])}); not averaged")
    if labels:
        label = labels.pop()
        return ("Medium-High" if label == "BBB+" else "Low"), "classified-BBB", \
               f"B3DB {label} (n={len(entry['labels'])})"
    return None, None, None


def main():
    (agents_path, v1_path, struct_path, reg_path, cls_path,
     out_csv, out_md) = sys.argv[1:8]

    v1 = {(r["target"], r["agent"]): r for r in
          csv.DictReader(open(v1_path, newline="", encoding="utf-8"))}

    structures, rejected = {}, []
    for row in csv.DictReader(open(struct_path, newline="", encoding="utf-8")):
        mol, skeleton = normalize(row["smiles"])
        if mol is None:
            rejected.append((row["agent"], "RDKit could not parse the SMILES"))
            continue
        structures[row["agent"]] = {
            "mol": mol, "skeleton": skeleton, "chembl": row["chembl_id"],
            "features": featurize(mol)}

    index, training = index_b3db(reg_path, cls_path)

    query_skeletons = {s["skeleton"] for s in structures.values()}
    X = np.asarray([f for skel, f, _y in training if skel not in query_skeletons])
    y = np.asarray([lab for skel, _f, lab in training if skel not in query_skeletons])
    held_out = len(training) - len(X)
    model = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=5000, random_state=SEED))
    model.fit(X, y)

    out_rows, changes, tier3_rows, matched, unmatched = [], [], [], [], []
    for key, row in v1.items():
        target, agent = key
        new = dict(row)
        structure = structures.get(agent)
        if structure:
            (matched if structure["skeleton"] in index else unmatched).append(agent)

        # Tier 1 keeps precedence: anything already resolved above Tier 2 stands.
        if row["provenance"] not in {"unknown"} and int(row["evidence_tier"] or 0) == 1:
            out_rows.append(new)
            continue

        if structure and structure["skeleton"] in index:
            category, provenance, detail = tier2(index[structure["skeleton"]])
            if category:
                score = CATEGORY[category]
                half = HALF_WIDTH[provenance]
                new.update({
                    "category": category, "score": score,
                    "lo": round(max(0.0, score - half), 3),
                    "hi": round(min(1.0, score + half), 3),
                    "provenance": provenance, "evidence_tier": 2,
                    "detail": f"{detail} [structural match, {structure['chembl']}]"})
                if row["provenance"] != provenance or row["category"] != category:
                    changes.append((agent, target, row["category"],
                                    row["provenance"], category, provenance,
                                    "structural B3DB match"))
                out_rows.append(new)
                continue

        if (row["provenance"] == "unknown" and row["modality"] == "small-molecule"
                and structure):
            probability = float(
                model.predict_proba(np.asarray([structure["features"]]))[0, 1])
            new.update({
                "category": "predicted", "score": round(probability, 3),
                "lo": round(max(0.0, probability - HALF_WIDTH["predicted"]), 3),
                "hi": round(min(1.0, probability + HALF_WIDTH["predicted"]), 3),
                "provenance": "predicted", "evidence_tier": 3,
                "detail": f"P(BBB+) = {probability:.3f} from 7 descriptors "
                          f"[{structure['chembl']}]"})
            tier3_rows.append((agent, target, probability))
            changes.append((agent, target, row["category"], row["provenance"],
                            "predicted", "predicted", "Tier 3 applied"))
        out_rows.append(new)

    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)

    resolved_before = sum(1 for r in v1.values() if r["provenance"] != "unknown")
    resolved_after = sum(1 for r in out_rows if r["provenance"] != "unknown")

    lines = [
        "# Reachability re-scored with retrieved structures\n",
        "Run 17 October 2026. Structures from ChEMBL "
        "(`data/structures.csv`), B3DB matched on the InChIKey skeleton after "
        "salt stripping and charge neutralisation.\n",
        f"- Structures accepted: {len(structures)} of "
        f"{len(structures) + len(rejected)}.",
        f"- B3DB indexed on structure: {len(index)} distinct skeletons.",
        f"- Query compounds found in B3DB by structure: "
        f"{len(sorted(set(matched)))} ({', '.join(sorted(set(matched)))}).",
        f"- Query compounds absent from B3DB: "
        f"{', '.join(sorted(set(unmatched))) if unmatched else 'none'}.",
        f"- B3DB training records withheld because they match a query "
        f"compound: {held_out}. Tier 3 training n = {len(X)}.",
        f"- **Triples resolving on reachability: {resolved_before} of "
        f"{len(v1)} before, {resolved_after} of {len(out_rows)} after.**",
        "",
        "## Every score that changed\n",
        "| agent | target | was | now | why |",
        "|---|---|---|---|---|",
    ]
    for agent, target, old_c, old_p, new_c, new_p, why in changes:
        lines.append(f"| {agent} | {target} | {old_c} ({old_p}) | "
                     f"{new_c} ({new_p}) | {why} |")
    if not changes:
        lines.append("| — | — | — | — | nothing changed |")

    lines += [
        "",
        "## What this does and does not fix\n",
        "The matching fix is a correction, not a new result: name matching was "
        "under-recovering B3DB, and any count of unresolved agents from the "
        "13 October run is therefore an overcount. Tier 3 is a genuine "
        "extension of coverage, but it is prediction, labelled as such and "
        "carrying the widest band short of `conflicting`, and it is applied to "
        "no antibody or ADC because the descriptors are undefined for them.",
    ]
    if rejected:
        lines += ["", "Structures rejected: "
                  + "; ".join(f"{a} ({why})" for a, why in rejected) + "."]

    with open(out_md, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
