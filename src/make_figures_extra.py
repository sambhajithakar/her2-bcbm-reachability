#!/usr/bin/env python3
"""Figures 1, 4 and 5: census flow, collapsing sensitivity, Tier 3 validation.

Usage: make_figures_extra.py CENSUS COLLAPSE_MD B3DB_CLS OUT_DIR
"""
import csv
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from sklearn.calibration import calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")
SEED = 20261013
plt.rcParams.update({"font.size": 11, "figure.dpi": 300,
                     "axes.spines.top": False, "axes.spines.right": False})


def figure1(census_path, out):
    rows = list(csv.DictReader(open(census_path, newline="", encoding="utf-8")))
    nodes = 23
    stages = [
        ("Records identified\nPubMed + bioRxiv/medRxiv/arXiv", "search complete;\nScopus and Web of Science outstanding"),
        ("Records screened on title and abstract", "excluded: no human BCBM;\nno molecular nomination"),
        ("Full texts assessed for eligibility", "excluded: review or case report;\nimaging- or radiotherapy-only;\nprimary CNS tumour"),
        (f"Nominations included\nn = {len(rows)} from {len({r['study'] for r in rows})} studies", "mapped onto canonical targets"),
        (f"Canonical target nodes\nn = {nodes}", ""),
    ]
    fig, ax = plt.subplots(figsize=(8.5, 7.5))
    ax.set_xlim(0, 10); ax.set_ylim(0, len(stages) * 2 + 0.5); ax.axis("off")
    for i, (main, side) in enumerate(stages):
        y = (len(stages) - i) * 2 - 1
        ax.add_patch(plt.Rectangle((0.3, y - 0.55), 4.6, 1.1, fill=True,
                                   facecolor="#eef2f8", edgecolor="#39506e", lw=1.2))
        ax.text(2.6, y, main, ha="center", va="center", fontsize=10.5)
        if side:
            ax.add_patch(plt.Rectangle((5.6, y - 0.5), 4.1, 1.0, fill=True,
                                       facecolor="#f7f7f7", edgecolor="#b0b0b0",
                                       lw=0.9, linestyle="--"))
            ax.text(7.65, y, side, ha="center", va="center", fontsize=8.5, color="#444")
            ax.annotate("", xy=(5.55, y), xytext=(4.95, y),
                        arrowprops=dict(arrowstyle="->", color="#888", lw=1))
        if i < len(stages) - 1:
            ax.annotate("", xy=(2.6, y - 1.45), xytext=(2.6, y - 0.6),
                        arrowprops=dict(arrowstyle="->", color="#39506e", lw=1.4))
    ax.set_title("Figure 1. Flow of records through the nomination census",
                 fontsize=12.5, loc="left", pad=14)
    fig.tight_layout(); fig.savefig(f"{out}/figure1_flow.png"); plt.close(fig)


def figure4(collapse_md, out):
    text = open(collapse_md, encoding="utf-8").read()
    rows = []
    for line in text.splitlines():
        m = re.match(r"\|\s*`?([a-zA-Z_0-9]+)`?\s*\|.*\|\s*([A-Za-z0-9/,\-\s]+)\s*\|$", line)
        if m and m.group(1) not in {"scenario"}:
            members = [x.strip().replace("/", "-") for x in m.group(2).split(",")]
            rows.append((m.group(1), members))
    nodes = sorted({n for _s, members in rows for n in members})
    grid = np.array([[1 if n in members else 0 for n in nodes] for _s, members in rows])
    fig, ax = plt.subplots(figsize=(1.6 + 1.5 * len(nodes), 0.7 + 0.55 * len(rows)))
    ax.imshow(grid, cmap="Blues", vmin=0, vmax=1.6, aspect="auto")
    ax.set_xticks(range(len(nodes))); ax.set_xticklabels(nodes, fontsize=10)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([s.replace("_", " ") for s, _m in rows], fontsize=10)
    for i in range(len(rows)):
        for j in range(len(nodes)):
            if grid[i, j]:
                ax.text(j, i, "on front", ha="center", va="center",
                        fontsize=8.5, color="#14365c")
    ax.set_xticks(np.arange(-.5, len(nodes), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", length=0)
    ax.set_title("Figure 4. Pareto-front membership under each target-collapsing scenario",
                 fontsize=12, loc="left", pad=12)
    fig.tight_layout(); fig.savefig(f"{out}/figure4_collapse.png"); plt.close(fig)


def _query_skeletons(path):
    """InChIKey skeletons of the agents Tier 3 is applied to.

    Withheld before the split so the curves describe the same model that
    produces the applied predictions (see src/tier3_predictor.py).
    """
    from rdkit.Chem.MolStandardize import rdMolStandardize
    chooser = rdMolStandardize.LargestFragmentChooser()
    uncharger = rdMolStandardize.Uncharger()
    out = set()
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            mol = Chem.MolFromSmiles((row.get("smiles") or "").strip())
            if mol is None:
                continue
            key = Chem.MolToInchiKey(uncharger.uncharge(chooser.choose(mol)))
            if key:
                out.add(key.split("-")[0])
    return out


def _skeleton(smiles):
    from rdkit.Chem.MolStandardize import rdMolStandardize
    chooser = rdMolStandardize.LargestFragmentChooser()
    uncharger = rdMolStandardize.Uncharger()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    key = Chem.MolToInchiKey(uncharger.uncharge(chooser.choose(mol)))
    return key.split("-")[0] if key else None


def figure5(b3db_path, out, structures_path=None):
    smiles, labels = [], []
    with open(b3db_path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            smi = (row.get("SMILES") or "").strip()
            lab = (row.get("BBB+/BBB-") or "").strip()
            if smi and lab in {"BBB+", "BBB-"}:
                smiles.append(smi); labels.append(1 if lab == "BBB+" else 0)
    if structures_path:
        queries = _query_skeletons(structures_path)
        kept = [i for i, smi in enumerate(smiles)
                if _skeleton(smi) not in queries]
        smiles = [smiles[i] for i in kept]
        labels = [labels[i] for i in kept]
    X, y = [], []
    for smi, lab in zip(smiles, labels):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        X.append([Descriptors.MolWt(mol), Crippen.MolLogP(mol),
                  rdMolDescriptors.CalcTPSA(mol), rdMolDescriptors.CalcNumHBD(mol),
                  rdMolDescriptors.CalcNumHBA(mol),
                  rdMolDescriptors.CalcNumRotatableBonds(mol),
                  rdMolDescriptors.CalcNumAromaticRings(mol)])
        y.append(lab)
    X, y = np.asarray(X), np.asarray(y)
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.3,
                                          random_state=SEED, stratify=y)

    def ramp(v, good, bad):
        if good < bad:
            return 1.0 if v <= good else 0.0 if v >= bad else (bad - v) / (bad - good)
        return 1.0 if v >= good else 0.0 if v <= bad else (v - bad) / (good - bad)

    def rule(f):
        mw, clogp, tpsa, hbd, hba, rotb, _ar = f
        return float(np.mean([ramp(mw, 360., 500.), ramp(clogp, 3., 5.),
                              ramp(tpsa, 40., 120.) * ramp(tpsa, 90., 20.),
                              ramp(hbd, .5, 3.5), ramp(hba, 4., 10.),
                              ramp(rotb, 4., 10.)]))

    rule_scores = np.array([rule(r) for r in Xte])
    model = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=5000, random_state=SEED))
    model.fit(Xtr, ytr)
    proba = model.predict_proba(Xte)[:, 1]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
    for scores, label, colour in ((proba, "logistic regression", "#1f4e9c"),
                                  (rule_scores, "CNS-MPO-style rule", "#b03a2e")):
        fpr, tpr, _ = roc_curve(yte, scores)
        ax1.plot(fpr, tpr, color=colour, lw=2,
                 label=f"{label} (AUC {roc_auc_score(yte, scores):.3f})")
    ax1.plot([0, 1], [0, 1], "--", color="#999", lw=1, label="chance")
    ax1.set_xlabel("false positive rate"); ax1.set_ylabel("true positive rate")
    ax1.set_title("Held-out discrimination", fontsize=11); ax1.legend(frameon=False, fontsize=9.5)

    frac, mean_pred = calibration_curve(yte, proba, n_bins=8, strategy="quantile")
    ax2.plot([0, 1], [0, 1], "--", color="#999", lw=1)
    ax2.plot(mean_pred, frac, "o-", color="#1f4e9c", lw=1.8)
    ax2.set_xlabel("mean predicted P(BBB+)"); ax2.set_ylabel("observed BBB+ fraction")
    ax2.set_title("Calibration, eight quantile bins", fontsize=11)
    fig.suptitle("Figure 5. Tier 3 predictor on compounds held out of B3DB",
                 fontsize=12.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(f"{out}/figure5_tier3.png"); plt.close(fig)


if __name__ == "__main__":
    census, collapse, b3db, outdir = sys.argv[1:5]
    structures = sys.argv[5] if len(sys.argv) > 5 else None
    figure1(census, outdir)
    figure4(collapse, outdir)
    figure5(b3db, outdir, structures)
    print("wrote figures 1, 4 and 5 to", outdir)
