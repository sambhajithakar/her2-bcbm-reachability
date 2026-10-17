# Tier 3 predictor: fitted, and measured on held-out compounds

Source: B3DB classification set. 7796 labelled records, 7794 parsed by RDKit (100.0%). Stratified 70/30 split, seed 20261013. Held-out n = 2339, BBB+ base rate 63.5%. 11 record(s) whose InChIKey skeleton matches an agent the predictor is applied to were removed before the split, so the discrimination below belongs to the same model that produces the applied predictions.

| predictor | held-out ROC AUC |
|---|---|
| CNS-MPO-style rule (no training data) | 0.750 |
| Logistic regression on the same 7 descriptors | 0.840 |

Difference in AUC (regression minus rule), bootstrap 95% CI: +0.090 [+0.071, +0.108].

Logistic regression also: accuracy 0.807 at a 0.5 threshold, Brier score 0.146.

Standardised coefficients (positive favours BBB+): MW +0.087, cLogP +0.160, TPSA -1.746, HBD -0.497, HBA +0.462, RotB -0.387, ArRings -0.278.

Calibration on the held-out split, by quantile bin:

| mean predicted | observed BBB+ fraction |
|---|---|
| 0.066 | 0.085 |
| 0.315 | 0.253 |
| 0.570 | 0.529 |
| 0.695 | 0.737 |
| 0.775 | 0.836 |
| 0.837 | 0.857 |
| 0.888 | 0.880 |
| 0.934 | 0.908 |

## What this predictor is, and what it is not

Tier 3 is the last resort in the registered order: physicochemical prediction, used only for small molecules and only where clinical pharmacokinetics and B3DB are both empty. The agents it is applied to, and their predicted probabilities, are in `tier3_applied.md`.

Two limits survive having the structures. The descriptors are tumour-agnostic and say nothing about efflux liability, which is the mechanism that demoted the CDK4/6 inhibitors on clinical evidence, so a compound can score well here and still be pumped out at the barrier. And they are undefined for antibodies and antibody-drug conjugates, so the four unresolved biologics are outside the predictor's domain regardless.
