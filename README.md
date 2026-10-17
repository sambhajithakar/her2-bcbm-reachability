# CNS-reachability-aware target prioritization in HER2-positive breast cancer brain metastasis

Analysis archive for the manuscript **"Machine-Learning-Assisted, CNS-Reachability-Aware
Prioritization of Therapeutic Targets in HER2-Positive Breast Cancer Brain Metastasis:
A Replication Census and Two-Axis Evidence Map"**, submitted to the *Journal of the
American Medical Informatics Association*, Research and Applications.

Every number reported in the manuscript comes from this pipeline. Nothing in the paper
was computed outside it.

## Reproduce everything

```
conda env create -f environment.yml      # or: mamba env create -f environment.yml
snakemake --cores 1 --use-conda
```

The archive is self-contained: the B3DB barrier dataset is bundled (it is CC0), chemical
structures are bundled as retrieved, and no network access is needed. Outputs are written
to `results/` and `results/v2/`; `results/v2/` holds the five figures at 300 dpi as they
appear in the manuscript.

## What is here

| Path | What it is |
|---|---|
| `Snakefile` | The pipeline. Every reported result is a rule in it. |
| `environment.yml` | Pinned conda environment. |
| `src/` | The analysis code: scoring, two-axis mapping, controls, collapsing sensitivity, the Tier 3 predictor and the figure scripts. |
| `data/census.csv` | The nomination census: 31 therapeutic-target nominations for breast cancer brain metastasis, 2015–2026, one row each, with cohort size, modality, HER2 stratification, evidence level, provisional tier and target compartment. |
| `data/agents.csv` | The agent table: 42 target–modality–agent triples with their Tier 1 evidence and source. Scoring happens at the triple, not the target. |
| `data/agent_target_map.csv` | Maps the 25 agent-table target labels onto canonical targets. |
| `data/target_map.csv` | The target-collapsing map, 30 rows with a note per decision. The sensitivity analysis contests each one across six scenarios. |
| `data/structures.csv` | ChEMBL identifiers and SMILES for 22 agents, as retrieved on 17 October 2026. |
| `data/b3db-raw/` | B3DB, bundled under CC0. See `LICENSE`. |
| `results/` | Every scoring output and log. Files ending `_v2` are the 17 October 2026 run reported in the manuscript; the unsuffixed files are the 13 October pilot, kept so the difference is inspectable. |
| `results/v2/` | The five manuscript figures, 300 dpi. |
| `PROTOCOL.md` | The study protocol, with its dated deviations log. **It was never deposited in a public registry** — see the note at the top of that file and the manuscript's Limitations. |
| `.zenodo.json` | Metadata for the Zenodo record. |

## Two things to know before you reuse this

**Tier 2 matches on chemical structure, not on compound name.** B3DB stores synonym lists
in its name field, so an exact-name lookup silently misses records — `"everolimus, rad001"`
was invisible to a search for `everolimus`. Matching on InChIKey skeleton after salt
stripping and charge neutralisation raised resolution from 21 to 30 of 42 triples. Any
pipeline joining curated chemistry to a target list on free-text names has this bug.

**The Tier 3 predictor's reported discrimination belongs to the model that made the
predictions.** `src/tier3_predictor.py` takes `data/structures.csv` and withholds every
B3DB record matching a query compound before the train/test split. An earlier version
reported discrimination from a model fitted on all of B3DB while the predictions came from
a model with those records withheld; that was corrected on 17 October 2026 and the
superseded figures are in the deviations log in `PROTOCOL.md` rather than deleted.

## Citing

See `CITATION.cff`. Please cite both this archive and the paper it supports.

## Licence

MIT for the code and for the data files we authored; bundled third-party data keep their
own terms. See `LICENSE`.
