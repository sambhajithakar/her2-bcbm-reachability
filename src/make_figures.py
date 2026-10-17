#!/usr/bin/env python3
"""Figures 2 and 3. Usage: make_figures.py SCORES_CSV NODES_CSV OUT_DIR"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROV_COLOR = {
    "clinical-RCT": "#1b5e20", "clinical-PK": "#2e7d32",
    "in-vivo-intracranial": "#558b2f", "barrier-model": "#f9a825",
    "measured-logBB": "#1565c0", "classified-BBB": "#42a5f5",
    "inferred-by-class": "#8e24aa", "conflicting": "#c62828",
    "predicted": "#9e9e9e",
}
MODALITY_ORDER = ["small-molecule", "ADC", "peptide", "naked-mAb"]


def figure2(scores, out_dir):
    data = scores[scores.score.notna()].copy()
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    yticks, ylabels, y = [], [], 0
    for modality in MODALITY_ORDER:
        block = data[data.modality == modality].sort_values("score")
        if block.empty:
            continue
        for _, row in block.iterrows():
            ax.errorbar(row.score, y,
                        xerr=[[row.score - row.lo], [row.hi - row.score]],
                        fmt="o", ms=7, capsize=3, lw=1.4,
                        color=PROV_COLOR.get(row.provenance, "#555555"))
            label = f"{row.agent} — {row.target}  ({row.compartment})"
            yticks.append(y)
            ylabels.append(label)
            y += 1
        y += 0.8
        ax.text(1.02, y - 1.4, modality, transform=ax.get_yaxis_transform(),
                fontsize=9, style="italic", va="center")
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlim(-0.03, 1.03)
    ax.set_xlabel("CNS reachability score (0–1), with provenance-set uncertainty band")
    ax.axvspan(0, 0.40, color="#f5f5f5", zorder=0)
    ax.axvspan(0.60, 1.0, color="#f0f7f0", zorder=0)
    ax.set_title("Figure 2. CNS reachability of target–modality–agent triples",
                 fontsize=11, loc="left")
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=p)
               for p, c in PROV_COLOR.items() if p in set(data.provenance)]
    ax.legend(handles=handles, fontsize=7.5, loc="lower right",
              title="evidence provenance", title_fontsize=8, framealpha=0.95)
    fig.tight_layout()
    fig.savefig(f"{out_dir}/figure2_reachability.png", dpi=300)
    plt.close(fig)


def figure3(nodes, out_dir):
    data = nodes[nodes.reach.notna()].copy()
    fig, ax = plt.subplots(figsize=(9.0, 6.0))

    # Spread points vertically *within* their tier so labels do not collide.
    positions = {}
    for tier_value, block in data.groupby("tier_value"):
        block = block.sort_values("reach")
        n = len(block)
        offsets = np.linspace(-0.30, 0.30, n) if n > 1 else np.array([0.0])
        for offset, idx in zip(offsets, block.index):
            positions[idx] = tier_value + offset

    for idx, row in data.iterrows():
        yy = positions[idx]
        on_front = bool(row.on_pareto_front)
        ax.errorbar(row.reach, yy,
                    xerr=[[row.reach - row.reach_lo], [row.reach_hi - row.reach]],
                    fmt="D" if on_front else "o",
                    ms=9 if on_front else 6,
                    color="#b71c1c" if on_front else "#5c6bc0",
                    capsize=2.5, lw=1.3, zorder=3 if on_front else 2)
        ax.annotate(f"{row.node} · {row.best_agent} ({row.pareto_frequency:.0%})",
                    (row.reach_hi, yy), textcoords="offset points",
                    xytext=(7, 0), fontsize=7.6, va="center",
                    fontweight="bold" if on_front else "normal",
                    color="#b71c1c" if on_front else "#333333")

    ax.set_yticks([1, 2, 3])
    ax.set_yticklabels(["Tier C\nsingle cohort,\nassociation only",
                        "Tier B\n2+ cohorts or\nin vivo validation",
                        "Tier A\n2+ cohorts, 2+ modalities,\nfunctional validation"],
                       fontsize=8)
    ax.set_xlim(0, 1.85)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_ylim(0.4, 3.6)
    ax.axvline(0.60, ls=":", lw=1, color="#777777")
    ax.text(0.605, 0.46, "required region floor", fontsize=7,
            color="#777777", rotation=90, va="bottom")
    ax.set_xlabel("CNS reachability of the best-scoring agent for the node "
                  "(bars are provenance-set uncertainty)")
    ax.set_title("Figure 3. Two-axis prioritization map with Pareto front",
                 fontsize=11, loc="left", pad=16)
    ax.text(0, 1.02,
            "Red diamonds are non-dominated. The percentage is bootstrap Pareto "
            "membership over 2000 resamples of the census and the bands.",
            transform=ax.transAxes, fontsize=8, color="#555555")
    ax.grid(axis="x", ls=":", alpha=0.35)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{out_dir}/figure3_two_axis.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    scores_p, nodes_p, out_dir = sys.argv[1:4]
    figure2(pd.read_csv(scores_p), out_dir)
    figure3(pd.read_csv(nodes_p), out_dir)
    print(f"wrote figure2_reachability.png and figure3_two_axis.png to {out_dir}")
