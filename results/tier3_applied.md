# Tier 3 applied: the agents Tiers 1 and 2 could not reach

Run 17 October 2026, seed 20261013. Structures from ChEMBL (`data/structures.csv`), fitted on B3DB with every query compound removed from the training set.

## Integrity checks

- Structures retrieved: 9. Parsed by RDKit: 9.
- InChIKey recomputed from each SMILES and compared with the identifier ChEMBL returned: all 9 agree.
- Query compounds already present in B3DB, therefore left at Tier 2 and not predicted: everolimus.
- B3DB records removed from the training set because they match a query compound: 1. Training n = 7804.

## Predictions

| agent | target | P(BBB+) | band | ChEMBL |
|---|---|---|---|---|
| olaparib | HRR/PARP | 0.546 | [0.25, 0.85] | CHEMBL521686 |
| bemcentinib | AXL | 0.523 | [0.22, 0.82] | CHEMBL3809489 |
| paxalisib | PIK3CA | 0.516 | [0.22, 0.82] | CHEMBL3813842 |
| galunisertib | TGFB1 macrophage axis | 0.516 | [0.22, 0.82] | CHEMBL2364611 |
| selumetinib | MEK/ERK | 0.455 | [0.15, 0.76] | CHEMBL1614701 |
| selpercatinib | RET | 0.403 | [0.10, 0.70] | CHEMBL4559134 |
| capivasertib | PTEN | 0.341 | [0.04, 0.64] | CHEMBL2325741 |
| neratinib | ERBB2 | 0.322 | [0.02, 0.62] | CHEMBL180022 |

## How to read these

Every row carries the `predicted` provenance label and the widest band short of `conflicting`, because a prediction is weaker evidence than a measurement and weaker still than clinical pharmacokinetics. The held-out discrimination of this predictor is 0.834 ROC AUC (`tier3_validation.md`), against 0.739 for the untrained CNS-MPO-style rule, so these numbers are informative but should not be read as equal in standing to the Tier 1 scores they sit beside on the map.

Two limits survive having the structures. The descriptors say nothing about efflux liability, which is the mechanism that demoted the CDK4/6 inhibitors on clinical evidence, so a compound can score well here and still be pumped out at the barrier. And they are undefined for antibodies and antibody-drug conjugates, so the four unresolved biologics stay `unknown` by design, not by accident.
