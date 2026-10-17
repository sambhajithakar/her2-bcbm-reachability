# Pilot run of the prioritization pipeline

**Run date:** 13 October 2026 · **Seed:** 20261013 · **Bootstrap resamples:** 2000
**Reproduce:** `snakemake --cores 1 --use-conda` from `/mnt/project-files/analysis`

> ⛔ **These are pilot numbers, not the manuscript's Results.** Three conditions in the
> pre-registered protocol are not met: the census is not closed (Scopus and Web of
> Science are unsearched), extraction and tiering were done by one operator rather than
> two independent reviewers (OSF decisions 5 and 8), and the reachability curation was
> not blinded to the tiering. They are reported here to show the pipeline runs end to
> end and to expose what it will and will not be able to say. Nothing below belongs in
> the manuscript until the protocol conditions are satisfied.

## What ran

| Step | Script | Output |
|---|---|---|
| Reachability scoring, three evidence tiers | `src/score_reachability.py` | `reachability_scores.csv` |
| Pre-registered control test | `src/validate_controls.py` | `control_test.md` |
| Tiering, two-axis map, Pareto front, bootstrap | `src/two_axis_map.py` | `two_axis_nodes.csv` |
| Figures 2 and 3 | `src/make_figures.py` | two PNGs at 300 dpi |
| Collapse sensitivity | `src/collapse_sensitivity.py` | `collapse_sensitivity.md` |
| Tier 3 predictor, fitted and held-out tested | `src/tier3_predictor.py` | `tier3_validation.md` |

Inputs: a 31-row nomination census (`data/census.csv`), a 42-row target–modality–agent
table (`data/agents.csv`), two declared target-collapsing maps, and the two B3DB files
downloaded on 8 October 2026 (960 compounds with usable logBB, 5093 with a
classification, after parsing).

## 1. The control test passes, and the strict pair separates

All five pre-registered controls land where the protocol requires. The strict test —
lapatinib against tucatinib and against poziotinib, same target family and same
modality — separates with non-overlapping uncertainty bands: lapatinib [0.08, 0.32]
against tucatinib [0.72, 0.88] and poziotinib [0.68, 0.92].

The honest qualification is in `control_test.md` and is repeated here. Four of the five
controls were curated by hand at Tier 1 before the test ran, so for those the test
checks propagation and band arithmetic rather than independent discovery. Only
lapatinib is placed independently, by B3DB's barrier-negative classification at Tier 2.

## 2. The reachability axis can only be populated for a minority of targets

Of 42 target–modality–agent triples, **21 resolved** and 21 did not. Of the 21 that did
not, 8 are targets with no clinical-grade agent at all; the remaining 13 are real agents
that neither clinical pharmacokinetics nor B3DB covers — neratinib, paxalisib,
everolimus, capivasertib, olaparib, selumetinib, bemcentinib, selpercatinib,
galunisertib, two ADCs, and pembrolizumab twice. Tier 3 prediction was not run, so these
are labelled `unknown` and carry no score, per OSF decision 12.

At the level of target nodes the bound is sharper: **9 of 23** nodes have any agent whose
reachability resolves. This is the number that limits what the axis can presently claim,
and it is a result in its own right rather than a caveat.

## 3. The Pareto front is small, and CDK4/6 falls off it

Three nodes are non-dominated: **PI3K–MTOR** (Tier A, reachability 0.65 via temsirolimus,
bootstrap Pareto membership 74%), **ERBB3–NRG1** (Tier B, 0.80 via poziotinib, 40%), and
**IL34–CSF1R** (Tier B, 0.80 via anti-CSF1R antibody, 39%).

The informative exclusion is **CDK4/6**. It reaches Tier A on the census — nominated by
Yang 2026 on organoid evidence and again by Rastad 2025 from public-data reanalysis —
and is then dominated, because all three approved inhibitors score 0.35 on
pharmacokinetic evidence of efflux limitation. Bootstrap Pareto membership 13%. This is
the pre-specified Discussion point arriving as a computed result rather than an
assertion: a target can be among the best replicated and still not be the best bet, and
only the second axis shows it.

The same mechanism separates two PARP routes. Niraparib scores 0.65 on clinical
pharmacokinetics; talazoparib scores 0.20 on B3DB's barrier-negative classification.
"Use a PARP inhibitor" and "use niraparib" are not the same recommendation.

## 4. The compartment rule has a visible consequence

The naked-antibody prior was withheld twice, for pembrolizumab against both immune
nodes, because those targets are not tumour-intrinsic. Those nodes are therefore
`unknown` rather than scored Low. Meanwhile the anti-CSF1R antibody — also a naked
antibody, but against a myeloid target at the invading edge — scores 0.80 on in vivo
intracranial evidence. Had the modality prior been applied unconditionally, as the
earlier draft of the scheme did, IL34–CSF1R would have been scored Low and would not
have appeared on the Pareto front at all.

## 5. One input error the run exposed

The first execution gave ERBB3–NRG1 no agents, because poziotinib and lapatinib were
listed only under ERBB2. Ippolitov's nomination is the resistance route, not the
receptor, so both agents were added under ERBB3–NRG1 and the pipeline re-run. The node
then entered the Pareto front. Recorded because it changed a headline output.

## 6. The Pareto front is only partly robust to the collapsing decisions

`results/collapse_sensitivity.md` re-runs tiering and the front with each
contested collapse undone. **ERBB3–NRG1 and IL34–CSF1R stay on the front in every
scenario, including the one with no collapsing at all. PI3K–MTOR does not.**

Splitting mTOR back out from the PI3K node drops PI3K–MTOR off the front, because
the merged node's Tier A rests on counting Ni 2016 and Yang 2026 as cohorts for one
node. CDK4/6 then enters the front in its place — not because it became more
reachable, but because the higher-tier competitor left. Separating PIK3CA, the two
immune nodes, or HRR from HRD changes nothing.

So the two in vivo-anchored nodes are findings and PI3K–MTOR is, for now, an
artefact of a judgement call. That has to be reported that way.

## 7. Tier 3 is implemented and characterised, but cannot be applied here

`results/tier3_validation.md`. Fitting on B3DB with a stratified 70/30 split
(7805 of 7807 records parsed by RDKit, held-out n = 2342, BBB+ base rate 63.5%):

| predictor | held-out ROC AUC |
|---|---|
| CNS-MPO-style rule, no training data | 0.739 |
| Logistic regression on 7 descriptors | 0.834 |

Difference +0.095, bootstrap 95% CI [+0.076, +0.114]; accuracy 0.795, Brier 0.150,
and calibration close to the diagonal across all eight quantile bins. Polar surface
area dominates (standardised coefficient −1.87).

This is useful mainly as a measurement of what Tier 3 is worth: clearly informative,
clearly worse than clinical pharmacokinetics, which is what the wide Tier 3
uncertainty band already encodes.

**It still cannot be applied to the 13 unresolved agents**, because a predictor needs
a structure and no structure source is reachable: PubChem, ChEMBL and the NCBI
E-utilities are all denied by this environment's network policy, and only
raw.githubusercontent.com resolves. Writing SMILES from memory would risk a confident
wrong score, so those agents stay `unknown`.

## Honest limits of this pilot

1. The census is not closed. Scopus and Web of Science are unsearched.
2. One operator did extraction, tiering and reachability curation. The protocol requires
   two independent reviewers with the tiering blind to reachability.
3. Target collapsing moves the result, and now measurably so: the sensitivity analysis
   in section 6 shows one of the three front members is an artefact of it.
4. The 0–1 category midpoints and provenance-set band widths are a defensible encoding,
   not a measurement. Results should be read as ordinal.
5. Tier 3 is implemented and validated but not applicable here for want of a reachable
   structure source, so 13 agents stay unresolved that a connected environment could score.
6. The benchmark has not been run: the metrics module is complete and self-tested, but
   the agentic workflow has not been executed and the literature-volume counts cannot be
   collected while the E-utilities are unreachable.
7. The transcriptomic component has not been started. GEO is denied by the same network
   policy, so no cohort has been downloaded.
