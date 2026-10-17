#!/usr/bin/env python3
"""Axis 2: tiered, modality- and compartment-aware CNS reachability scoring.

Implements the pre-registered rules (OSF decisions 9-13, 31-32) without
deviation. Nothing is imputed: an agent that no evidence tier resolves is
labelled `unknown` and carries no score.

Usage:
    score_reachability.py AGENTS_CSV B3DB_REGRESSION_TSV B3DB_CLASSIFICATION_TSV OUT_CSV
"""
import csv
import sys

import pandas as pd

# Category midpoints, pre-specified. Reachability is reported on 0-1.
CATEGORY = {
    "High": 0.80,
    "Medium-High": 0.65,
    "Medium": 0.50,
    "Low-Medium": 0.35,
    "Low": 0.20,
}

# Half-width of the uncertainty band, set by evidence provenance, not by
# category. Weaker provenance means a wider band; `conflicting` is widest.
HALF_WIDTH = {
    "clinical-RCT": 0.08,
    "clinical-PK": 0.10,
    "in-vivo-intracranial": 0.12,
    "barrier-model": 0.18,
    "measured-logBB": 0.15,
    "classified-BBB": 0.20,
    "inferred-by-class": 0.25,
    "predicted": 0.30,
    "conflicting": 0.45,
}

# Tier 1: direct clinical or in vivo intracranial evidence.
TIER1 = {
    "intracranial-rct": ("High", "clinical-RCT"),
    "intracranial-clinical": ("High", "clinical-PK"),
    "in-vivo-intracranial": ("High", "in-vivo-intracranial"),
    "in-vivo-no-intracranial-response": ("Low", "in-vivo-intracranial"),
    "clinical-pk-not-efflux-limited": ("Medium-High", "clinical-PK"),
    "clinical-pk-efflux-limited": ("Low-Medium", "clinical-PK"),
    "barrier-model-assessed": ("Medium", "barrier-model"),
}

# logBB thresholds for Tier 2 quantitative records.
LOGBB_HIGH, LOGBB_LOW = 0.30, -0.50


def load_b3db(reg_path, cls_path):
    """Index B3DB by lowercased compound name. Both files are untrusted input."""
    reg, cls = {}, {}
    with open(reg_path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = (row.get("compound_name") or "").strip().lower()
            try:
                value = float(row.get("logBB"))
            except (TypeError, ValueError):
                continue
            if name:
                reg.setdefault(name, []).append(value)
    with open(cls_path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = (row.get("compound_name") or "").strip().lower()
            label = (row.get("BBB+/BBB-") or "").strip()
            if name and label:
                cls.setdefault(name, []).append(label)
    return reg, cls


def lookup(index, agent):
    """Exact name match first, then a contained-name match (B3DB carries salts
    and stereo-qualified synonyms). Returns every record found, so that
    contradictions stay visible instead of being silently resolved."""
    key = agent.strip().lower()
    if key in index:
        return index[key]
    hits = []
    for name, values in index.items():
        if key in name.split() or key == name.replace(" ", ""):
            hits.extend(values)
    return hits


def score_agent(row, reg, cls):
    """Apply the tiers in the pre-registered order and stop at the first that
    resolves. Returns (category, provenance, tier, detail)."""
    agent = row["agent"]
    evidence = (row.get("t1_evidence") or "none").strip()
    modality = (row.get("modality") or "").strip()
    compartment = (row.get("compartment") or "").strip()

    if evidence == "no-clinical-agent":
        return None, "unknown", 0, "no clinical-grade agent exists for this target"

    # --- Tier 1 -----------------------------------------------------------
    if evidence in TIER1:
        category, provenance = TIER1[evidence]
        return category, provenance, 1, row.get("t1_note", "")

    # Modality prior, conditioned on target compartment (OSF decision 32).
    # The low prior for naked antibodies applies only where the target cell
    # sits behind an intact barrier.
    if evidence == "class-prior-low":
        if modality == "naked-mAb" and compartment != "tumour-intrinsic":
            return (None, "unknown", 0,
                    "naked-mAb prior withheld: target compartment is "
                    f"{compartment}, so the prior does not transfer")
        return "Low", "inferred-by-class", 1, "modality prior, tumour-intrinsic compartment"

    # --- Tier 2: B3DB -----------------------------------------------------
    logbb = lookup(reg, agent)
    labels = lookup(cls, agent)

    if logbb:
        lo, hi = min(logbb), max(logbb)
        if hi >= LOGBB_HIGH and lo <= LOGBB_LOW:
            return "Medium", "conflicting", 2, (
                f"B3DB logBB records span {lo:.2f} to {hi:.2f}; not averaged")
        mean = sum(logbb) / len(logbb)
        if mean >= LOGBB_HIGH:
            category = "High"
        elif mean <= LOGBB_LOW:
            category = "Low"
        else:
            category = "Medium"
        return category, "measured-logBB", 2, (
            f"B3DB logBB n={len(logbb)}, mean {mean:.2f}")

    if labels:
        unique = set(labels)
        if unique == {"BBB+", "BBB-"}:
            return "Medium", "conflicting", 2, (
                f"B3DB holds both BBB+ and BBB- for this compound (n={len(labels)}); "
                "not averaged")
        label = unique.pop()
        category = "Medium-High" if label == "BBB+" else "Low"
        return category, "classified-BBB", 2, f"B3DB {label} (n={len(labels)})"

    # --- Tier 3 -----------------------------------------------------------
    # Physicochemical prediction requires descriptors this environment cannot
    # compute (no rdkit). Per OSF decision 12 nothing is imputed.
    return None, "unknown", 0, "not resolved at Tier 1 or Tier 2; Tier 3 not run"


def main():
    agents_path, reg_path, cls_path, out_path = sys.argv[1:5]
    agents = pd.read_csv(agents_path).fillna("")
    reg, cls = load_b3db(reg_path, cls_path)
    print(f"B3DB indexed: {len(reg)} compounds with logBB, "
          f"{len(cls)} with a classification", file=sys.stderr)

    rows = []
    for _, row in agents.iterrows():
        category, provenance, tier, detail = score_agent(row, reg, cls)
        score = CATEGORY.get(category) if category else None
        half = HALF_WIDTH.get(provenance)
        rows.append({
            "target": row["target"],
            "agent": row["agent"],
            "modality": row["modality"],
            "compartment": row["compartment"],
            "category": category or "unknown",
            "score": score,
            "lo": max(0.0, score - half) if score is not None else None,
            "hi": min(1.0, score + half) if score is not None else None,
            "provenance": provenance,
            "evidence_tier": tier,
            "detail": detail,
            "source": row.get("t1_source", ""),
        })

    out = pd.DataFrame(rows)
    out.to_csv(out_path, index=False)

    scorable = out[out.category != "unknown"]
    print(f"\nscored {len(scorable)} of {len(out)} target-modality-agent triples", file=sys.stderr)
    print("\nresolvability by evidence tier:", file=sys.stderr)
    print(out.evidence_tier.value_counts().sort_index().to_string(), file=sys.stderr)
    print("\nprovenance:", file=sys.stderr)
    print(out.provenance.value_counts().to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
