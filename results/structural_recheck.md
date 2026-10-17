# Reachability re-scored with retrieved structures

Run 17 October 2026. Structures from ChEMBL (`data/structures.csv`), B3DB matched on the InChIKey skeleton after salt stripping and charge neutralisation.

- Structures accepted: 22 of 22.
- B3DB indexed on structure: 4019 distinct skeletons.
- Query compounds found in B3DB by structure: 8 (alpelisib, atorvastatin, everolimus, lapatinib, rosuvastatin, talazoparib, temsirolimus, trametinib).
- Query compounds absent from B3DB: abemaciclib, bemcentinib, capivasertib, galunisertib, neratinib, niraparib, olaparib, palbociclib, paxalisib, poziotinib, ribociclib, selpercatinib, selumetinib, tucatinib.
- B3DB training records withheld because they match a query compound: 11. Tier 3 training n = 7794.
- **Triples resolving on reachability: 21 of 42 before, 30 of 42 after.**

## Every score that changed

| agent | target | was | now | why |
|---|---|---|---|---|
| neratinib | ERBB2 | unknown (unknown) | predicted (predicted) | Tier 3 applied |
| paxalisib | PIK3CA | unknown (unknown) | predicted (predicted) | Tier 3 applied |
| everolimus | MTOR | unknown (unknown) | Medium-High (classified-BBB) | structural B3DB match |
| capivasertib | PTEN | unknown (unknown) | predicted (predicted) | Tier 3 applied |
| olaparib | HRR/PARP | unknown (unknown) | predicted (predicted) | Tier 3 applied |
| selumetinib | MEK/ERK | unknown (unknown) | predicted (predicted) | Tier 3 applied |
| bemcentinib | AXL | unknown (unknown) | predicted (predicted) | Tier 3 applied |
| selpercatinib | RET | unknown (unknown) | predicted (predicted) | Tier 3 applied |
| galunisertib | TGFB1 macrophage axis | unknown (unknown) | predicted (predicted) | Tier 3 applied |

## What this does and does not fix

The matching fix is a correction, not a new result: name matching was under-recovering B3DB, and any count of unresolved agents from the 13 October run is therefore an overcount. Tier 3 is a genuine extension of coverage, but it is prediction, labelled as such and carrying the widest band short of `conflicting`, and it is applied to no antibody or ADC because the descriptors are undefined for them.
