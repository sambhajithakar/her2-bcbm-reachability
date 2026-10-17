#!/usr/bin/env python3
"""Axis 1 (replication tier) x Axis 2 (CNS reachability), Pareto front, and
bootstrap uncertainty over both axes.

The axes are never collapsed into a single primary score (OSF decision 14).
Tiers are recomputed here from the census by the pre-registered rule rather
than taken from the provisional column.

Usage:
    two_axis_map.py CENSUS_CSV TARGET_MAP_CSV SCORES_CSV AGENT_MAP_CSV OUT_PREFIX
"""
import sys

import numpy as np
import pandas as pd

RNG = np.random.default_rng(20261013)
N_BOOT = 2000
TIER_VALUE = {"A": 3, "B": 2, "C": 1}
# Nodes that are not single targets; carried in the census but kept off the
# prioritization map.
NOT_A_TARGET = {"PATHWAY-LEVEL", "SPATIAL-SIGNATURE", "THY1-PANEL"}


def assign_tier(group):
    """Pre-registered tiering (OSF decision 6).

    Tier A: >=2 non-overlapping cohorts AND >=2 distinct modalities AND >=1
            functional validation.
    Tier B: >=2 independent cohorts with a single modality, OR one cohort with
            in vivo functional validation.
    Tier C: single cohort, association only.

    Cohort independence is approximated by distinct source identifiers. Where
    two rows share a source the cohort is counted once, which scores the node
    at the lower tier - the pre-registered response to overlap ambiguity
    (decision 7).
    """
    cohorts = group.source_id.nunique()
    modalities = set()
    for entry in group.modality:
        modalities.update(part.strip() for part in str(entry).split("+"))
    functional = group.evidence_level.isin(["iv", "vv"]).any()
    in_vivo = (group.evidence_level == "vv").any()

    if cohorts >= 2 and len(modalities) >= 2 and functional:
        return "A", cohorts, len(modalities)
    if (cohorts >= 2) or in_vivo:
        return "B", cohorts, len(modalities)
    return "C", cohorts, len(modalities)


def build_nodes(census, target_map, scores, agent_map):
    census = census.merge(target_map, left_on="target", right_on="census_target")
    scores = scores.merge(agent_map, left_on="target", right_on="agents_target")

    rows = []
    for node, group in census.groupby("canonical_target"):
        if node in NOT_A_TARGET:
            continue
        tier, cohorts, modalities = assign_tier(group)
        agents = scores[scores.canonical_target == node]
        scored = agents[agents.score.notna()]
        if len(scored):
            best = scored.loc[scored.score.idxmax()]
            reach, lo, hi = best.score, best.lo, best.hi
            best_agent, provenance = best.agent, best.provenance
        else:
            reach = lo = hi = np.nan
            best_agent = agents.agent.iloc[0] if len(agents) else "none-identified"
            provenance = "unknown"
        rows.append({
            "node": node,
            "tier": tier,
            "tier_value": TIER_VALUE[tier],
            "n_cohorts": cohorts,
            "n_modalities": modalities,
            "n_nominations": len(group),
            "her2_stratified": (group.her2_stratified == "yes").any(),
            "reach": reach, "reach_lo": lo, "reach_hi": hi,
            "best_agent": best_agent, "provenance": provenance,
            "n_agents_scored": len(scored), "n_agents": len(agents),
        })
    return pd.DataFrame(rows).sort_values(["tier_value", "reach"],
                                          ascending=[False, False])


def pareto(points):
    """Indices of points not dominated on both axes (higher is better)."""
    front = []
    for i, (xi, yi) in enumerate(points):
        if np.isnan(yi):
            continue
        dominated = any(
            (xj >= xi and yj >= yi) and (xj > xi or yj > yi)
            for j, (xj, yj) in enumerate(points)
            if j != i and not np.isnan(yj))
        if not dominated:
            front.append(i)
    return front


def bootstrap_front(census, target_map, scores, agent_map, nodes):
    """Resample the census (tier uncertainty) and sample within each
    reachability band (score uncertainty); record Pareto membership."""
    counts = {node: 0 for node in nodes.node}
    merged = census.merge(target_map, left_on="target", right_on="census_target")
    merged = merged[~merged.canonical_target.isin(NOT_A_TARGET)]
    index = merged.index.to_numpy()

    for _ in range(N_BOOT):
        draw = merged.loc[RNG.choice(index, size=len(index), replace=True)]
        tiers = {}
        for node, group in draw.groupby("canonical_target"):
            tiers[node] = TIER_VALUE[assign_tier(group)[0]]
        points, names = [], []
        for _, row in nodes.iterrows():
            if np.isnan(row.reach):
                continue
            tier_value = tiers.get(row.node)
            if tier_value is None:        # node absent from this resample
                continue
            points.append((tier_value, RNG.uniform(row.reach_lo, row.reach_hi)))
            names.append(row.node)
        for i in pareto(points):
            counts[names[i]] += 1
    return {node: count / N_BOOT for node, count in counts.items()}


def main():
    census_p, map_p, scores_p, agent_map_p, prefix = sys.argv[1:6]
    census = pd.read_csv(census_p)
    target_map = pd.read_csv(map_p)
    scores = pd.read_csv(scores_p)
    agent_map = pd.read_csv(agent_map_p)

    nodes = build_nodes(census, target_map, scores, agent_map)
    points = list(zip(nodes.tier_value, nodes.reach))
    front = set(nodes.index[i] for i in pareto(points))
    nodes["on_pareto_front"] = nodes.index.isin(front)

    freq = bootstrap_front(census, target_map, scores, agent_map, nodes)
    nodes["pareto_frequency"] = nodes.node.map(freq)
    nodes.to_csv(f"{prefix}_nodes.csv", index=False)

    scorable = nodes[nodes.reach.notna()]
    print(f"nodes on the map: {len(nodes)}  "
          f"(with a scored agent: {len(scorable)}, "
          f"unresolved: {len(nodes) - len(scorable)})")
    print(f"tier distribution: {nodes.tier.value_counts().sort_index().to_dict()}")
    print(f"HER2-stratified somewhere: {int(nodes.her2_stratified.sum())} of {len(nodes)}")
    print("\nPareto front (point estimate), with bootstrap membership frequency:")
    cols = ["node", "tier", "n_cohorts", "n_modalities", "reach",
            "best_agent", "provenance", "pareto_frequency"]
    print(nodes[nodes.on_pareto_front][cols].to_string(index=False))
    print("\nAll nodes with a scored agent:")
    print(scorable[cols + ["on_pareto_front"]].to_string(index=False))
    print("\nNodes with no agent resolvable at any evidence tier "
          "(reported, never imputed):")
    print(nodes[nodes.reach.isna()][["node", "tier", "n_agents", "best_agent"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
