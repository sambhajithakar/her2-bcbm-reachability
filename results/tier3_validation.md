# Tier 3 predictor: fitted, and measured on held-out compounds

Source: B3DB classification set. 7807 labelled records, 7805 parsed by RDKit (100.0%). Stratified 70/30 split, seed 20261013. Held-out n = 2342, BBB+ base rate 63.5%.

| predictor | held-out ROC AUC |
|---|---|
| CNS-MPO-style rule (no training data) | 0.739 |
| Logistic regression on the same 7 descriptors | 0.834 |

Difference in AUC (regression minus rule), bootstrap 95% CI: +0.095 [+0.076, +0.114].

Logistic regression also: accuracy 0.795 at a 0.5 threshold, Brier score 0.150.

Standardised coefficients (positive favours BBB+): MW +0.167, cLogP +0.133, TPSA -1.868, HBD -0.442, HBA +0.441, RotB -0.391, ArRings -0.279.

Calibration on the held-out split, by quantile bin:

| mean predicted | observed BBB+ fraction |
|---|---|
| 0.075 | 0.088 |
| 0.337 | 0.284 |
| 0.572 | 0.497 |
| 0.694 | 0.775 |
| 0.774 | 0.809 |
| 0.833 | 0.832 |
| 0.887 | 0.884 |
| 0.935 | 0.911 |

## Why this is implemented but not applied

Tier 3 exists in the protocol for agents that Tiers 1 and 2 leave unresolved. In this run that is 13 real agents: neratinib, paxalisib, everolimus, capivasertib, olaparib, selumetinib, bemcentinib, selpercatinib, galunisertib, two ADCs and pembrolizumab twice.

**None of them can be scored here, because a predictor needs a structure and no structure source is reachable from this environment.** PubChem, ChEMBL and the NCBI E-utilities are all denied by the network policy; only raw.githubusercontent.com resolves. Writing SMILES strings from memory is not an option: a silently wrong structure would produce a confident and wrong reachability score, which is precisely the failure the provenance labels exist to prevent. Those agents therefore stay `unknown`.

Two further limits apply even once structures are available. The descriptors are tumour-agnostic and say nothing about efflux liability, which is the mechanism that demoted the CDK4/6 inhibitors on clinical evidence. And they are undefined for antibodies and ADCs, so the four unresolved biologics are outside the predictor's domain regardless.
