# Pre-registration — Open Science Framework

> ⚠️ **This document was never deposited on the Open Science Framework.** It is the study protocol, written before the analysis and kept here unchanged apart from the deviations log in section 9. No registration, DOI or time stamp exists for it. The manuscript used to claim an OSF registration with a placeholder DOI; that claim was removed on 17 October 2026, and this document now ships with the submission as **Supplementary Methods S6** instead. If it is ever deposited, it must be described as a protocol archived after analysis, never as a pre-registration.


**Title**
Benchmarking Agentic Artificial Intelligence for Drug-Target Nomination in HER2-Positive Breast Cancer Brain Metastasis: A Leakage-Controlled, CNS-Reachability-Aware Evaluation

**Registered by** Sambhaji Balaso Thakar and Pradnya Namdev Ghorpade (Advanced Data Science and Gen AI Unit, Aditya BioNova Analytics, Pune, India) · **Date** `[[registration date]]` · **Registration type** Preregistration (OSF standard)
**Status at registration** No analysis has been run. No result exists.

> 🇮🇳 **हे का महत्त्वाचे:** हा दस्तऐवज OSF वर नोंदवल्यानंतर timestamp मिळतो. त्यानंतर कोणीही म्हणू शकत नाही की "निकाल पाहून गृहीतक बदलले." Benchmark-प्रकारच्या पेपरसाठी reviewer हेच पहिल्यांदा तपासतो. **Analysis सुरू करण्यापूर्वी हे नोंदवा.**

---

## 1. Study type

Retrospective computational study with two components: (a) a systematic census and two-axis prioritization of published therapeutic-target nominations, and (b) a benchmark of agentic LLM nomination against a temporal hold-out. No new human-subjects data. No wet-lab experiment.

## 2. Hypotheses, stated before any analysis

**H1 (replication).** A minority of targets nominated for breast cancer brain metastasis (BCBM) between 2015 and 2026 are supported by two or more independent cohorts *and* two or more measurement modalities.

**H2 (HER2 stratification).** A minority of BCBM target nominations report a HER2-positive subset size, so most nominations currently applied to HER2-positive disease rest on evidence that was not HER2-stratified.

**H3 (reachability re-ranking).** Scoring CNS reachability as an explicit axis materially changes the ranking of nominated targets relative to replication alone. Pre-specified concrete instance: targets nominated on organoid or in vitro evidence without pharmacokinetic consideration will include at least one whose agents have documented efflux-limited brain penetration.

**H4 (benchmark).** An agentic LLM workflow restricted to pre-2024 evidence will recover 2024–2026 nominations at a rate that does **not** exceed a literature-volume baseline by a margin that is both statistically and practically meaningful.
*Direction is pre-specified as a null-leaning hypothesis deliberately.* We expect such systems to recover well-studied targets. If they beat the baseline, that is a positive finding; if they do not, that is the finding. Either outcome is reported.

**H5 (memorization).** Performance under gene-name masking will be lower than unmasked performance, and the gap is our estimate of the contribution of memorization rather than reasoning.

## 3. Locked analytic decisions

These cannot be changed after registration without a documented, dated deviation note.

| # | Decision | Locked value |
|---|---|---|
| 1 | Census window | 1 Jan 2015 to freeze date |
| 2 | Databases | PubMed, bioRxiv, medRxiv, arXiv, Scopus, Web of Science |
| 3 | Eligibility | Human BCBM, or brain metastasis including breast primaries with breast-specific results; nominates ≥1 molecular entity as therapeutic target or vulnerability |
| 4 | Exclusions | Reviews without primary nomination; case reports without molecular nomination; imaging/radiotherapy-only; primary CNS tumours |
| 5 | Extraction | Two independent reviewers, fixed schema, Krippendorff α reported. **Named at registration: S.B.T. and P.N.G.** The pilot run of 13 October 2026 was single-operator and is therefore explicitly not a protocol-conformant extraction |
| 6 | Replication tiers | Tier A: ≥2 non-overlapping cohorts AND ≥2 modalities AND ≥1 functional validation. Tier B: ≥2 cohorts single modality, OR 1 cohort with in vivo validation. Tier C: single cohort, association only |
| 7 | Overlap ambiguity | Scored at the **lower** tier, ambiguity recorded |
| 8 | Blinding | Tiering performed blind to reachability, by a different reviewer |
| 9 | Reachability unit | **target–modality–agent triple**, not target and not target–modality |
| 10 | Reachability evidence tiers | 1 = clinical/CSF PK; 2 = B3DB (measured logBB, else categorical); 3 = prediction (physicochemical, CNS MPO) |
| 11 | Provenance labels | `clinical-PK`, `measured-logBB`, `classified-BBB`, `predicted`, `inferred-by-class`, `conflicting`, `unknown` |
| 12 | Unknown handling | Reported as `unknown`; **never imputed** |
| 13 | Contradictory source records | Labelled `conflicting`, widest uncertainty band; **never averaged** |
| 14 | Axis combination | Two-dimensional map with **Pareto front**; no single collapsed score as primary. Weighted score is secondary with the weight swept across its full range |
| 15 | Positive controls | Tucatinib (small molecule), trastuzumab deruxtecan (ADC) and **poziotinib** (small molecule, intracranial remission on MRI in a patient-derived HER2+ BCBM model, Ippolitov 2024) must score high |
| 16 | Negative controls | Trastuzumab alone (naked mAb) and **lapatinib** (small molecule, B3DB barrier-negative, *and* no intracranial remission in the same model in which poziotinib succeeded) must score low. Lapatinib is the strict test: same target family and modality as tucatinib and poziotinib, opposite required region |
| 17 | Control failure | Reported as scheme failure. **Not** corrected by re-weighting |
| 18 | Expression cohorts | GSE43837 (discovery, HER2-positive stratified), GSE184869, GSE125989, GSE248830, GSE110590, GSE193103, GSE209998 (replication); GSE186344, GSE179572 (compartment attribution) |
| 19 | Excluded cohorts | TCGA-BRCA, METABRIC, CPTAC-BRCA — no brain-metastasis samples |
| 20 | Leakage control | Normalization, batch correction and feature selection fitted on discovery only; leave-one-cohort-out; replication cohorts not inspected before lock |
| 21 | Benchmark hold-out | Ground truth = nominations first published ≥1 Jan 2024. System evidence restricted to ≤31 Dec 2023 |
| 22 | Primary endpoint | **precision@20 of the agentic workflow versus the literature-volume baseline** |
| 23 | Baselines | Literature-volume (primary comparator), knowledge-base (Open Targets, pre-2024, version-pinned), random |
| 24 | Secondary metrics | precision@k, recall@k, nDCG@k, MRR; Brier score; expected calibration error |
| 25 | Inference | Bootstrap 95% CIs; paired bootstrap and permutation for pipeline-vs-baseline |
| 26 | Multiplicity | Benjamini–Hochberg q<0.05 over the enumerated secondary family |
| 27 | Ablations | Gene-name masking; reachability removed; retrieval removed; self-consistency across n runs |
| 28 | Citation audit | Every emitted citation checked against PubMed for existence and claim support; verifiable-claim and hallucinated-citation rates reported |
| 29 | Model disclosure | Identifier, version, temperature, call date and full prompt logged for every call and released |
| 30 | Subgroups | Only those listed here are confirmatory. HER2-positive pooled with HER2-low is a **secondary** analysis |
| 31 | Target compartment recorded | Every scored entry records the compartment of the target cell: tumour-intrinsic, myeloid or stromal, vascular, or extracellular matrix |
| 32 | Modality prior conditioned on compartment | The low-CNS-exposure prior for naked antibodies applies to tumour-intrinsic targets behind an intact barrier. It is **not** applied unchanged to myeloid or stromal targets at the invading edge (Cheng 2025: anti-CSF1R antibody shrank intracranial tumours in an immunocompetent HER2+ cerebellar metastasis model). This refinement is pre-specified, not validated, and is declared as a limitation |
| 34 | Target collapsing | Every collapse of a census target onto a canonical node is recorded with its reason in the released map. The analysis is re-run with each contested collapse undone and with no collapsing at all; Pareto-front membership is reported per scenario, and a node non-dominated only under the primary map is reported as an artefact |
| 35 | Tier 3 structures | Structures come from a named public resource and are never written from memory. An agent whose structure cannot be retrieved is reported `unknown`, not predicted |
| 36 | Tier 3 must be measured before it is used | The predictor is evaluated on compounds held out of B3DB; held-out discrimination and the gain over the untrained rule are reported alongside any score it produces |
| 33 | Benchmark differentiation declared | Relation to BioStudyBench (arXiv 2610.07614), DBench-Bio (arXiv 2603.03322), HindSight (arXiv 2603.15164) and LongevityBench is tabulated in the manuscript. Our distinguishing features are the ranked **target-nomination** task, the **literature-volume primary baseline**, the **gene-name masking** ablation, and the coupling to a deliverability axis |

## 4. Sample size and power

The census size is not under our control; it is whatever the search returns, and we report it rather than target it. For the benchmark, precision@20 comparisons are powered by the size of the 2024–2026 ground-truth set, which is likewise fixed by the literature. We therefore pre-commit to reporting confidence intervals and effect sizes rather than relying on significance thresholds, and to stating explicitly where the ground-truth set is too small to separate the pipeline from a baseline.

## 5. Known threats, and what we will do about them

| Threat | Pre-committed response |
|---|---|
| **Model knowledge cutoff after 1 Jan 2024** contaminates the temporal split | Report every model's cutoff; prefer open-weight models with documented cutoffs for the primary analysis; treat the masking ablation as the demonstrable bound; state the limitation **in the Abstract** |
| Ranking merely tracks publication volume | The literature-volume baseline is the primary comparator, not an afterthought |
| Reachability evidence is heterogeneous | Provenance labels on every score; report the proportion resolvable at each tier as a result |
| B3DB does not cover modern agents | Already established by pilot lookup; Tier 1 hand curation is planned, not improvised |
| Published signatures unavailable in full | Sub-analysis is **dropped and reported as dropped**, not approximated |
| Outcome-dependent reframing | This document. Hypothesis directions are fixed above, including H4's null-leaning direction |

## 6. What would make us say the method failed

Pre-specified, so that failure is recognizable rather than reinterpreted:

1. Either positive control scores low, or either negative control scores high — in particular if **lapatinib and tucatinib cannot be separated**.
2. Inter-rater agreement on tiering is poor (Krippendorff α with a confidence interval overlapping values conventionally regarded as unacceptable), indicating the tier definitions are not operable.
3. Masked performance equals unmasked performance at chance level for both, indicating the benchmark cannot discriminate at all.
4. Fewer than a reportable minimum of candidate agents resolve at Tiers 1 or 2, indicating the reachability axis cannot currently be populated for this disease.

In each case we report the failure and its diagnosis rather than adjusting the scheme until a control passes.

## 7. Outstanding dependency, declared at registration

The preprint and arXiv pass has now been executed and closed (8-12 October 2026). It produced four material corrections, all folded into the protocol and the manuscript before registration:

1. **Prior art for the PARP-inhibitor idea.** Diossy et al., *Ann Oncol* 2018, PMID 29917049. H3's PARP-related instance is re-scoped to the HER2+/PTEN route and to agent selection within the class.
2. **The three post-cutoff LLM benchmarks have been read**, and none of them performs our task. BioStudyBench (arXiv 2610.07614) evaluates agents on 25 long-horizon data-analysis tasks drawn from post-cutoff studies, with a no-data/no-tools prior-knowledge baseline, no literature-volume comparator and no masking ablation; target nomination is not one of its tasks. DBench-Bio (arXiv 2603.03322) scores hypothesis question answering with no named baselines. HindSight (arXiv 2603.15164) matches generated AI/ML research ideas to later papers; its finding that LLM-judged novelty is *negatively* correlated with later impact is independent support for decision 23's choice of a literature-volume comparator over system self-assessment. Decision 33 records the differentiation.
3. **A fifth control agent.** Ippolitov et al., bioRxiv 2024.02.19.581073, supplies poziotinib as a positive control and raises lapatinib's negative-control status from a database classification to in vivo HER2+ BCBM evidence. Decisions 15 and 16 are updated.
4. **A modality-rule refinement.** Cheng et al., bioRxiv 2025.06.17.660224, shows an anti-CSF1R antibody shrinking intracranial tumours, which forces the naked-antibody prior to be conditioned on target compartment. Decisions 31 and 32 are new.

⚠️ **Still outstanding: Scopus and Web of Science (decision 2) have not been searched**, and the preprint pass was keyword-based rather than PRISMA-conformant. Both are part of the locked protocol and must be completed before the census is closed. Until then census completeness is not established. This is declared here rather than discovered later.

## 8. Deviations log

| Date | Decision changed | Reason |
|---|---|---|
| 2026-10-08 | None (clarification before registration) | First preprint pass identified Diossy et al. 2018 as prior art for the PARP-inhibitor idea; the PTEN-specific and agent-selection components of that hypothesis were re-scoped in the manuscript *before* registration, so this is a pre-registration clarification rather than a deviation. |
| 2026-10-12 | Decisions 15, 16, 18 amended; decisions 31, 32, 33 added (all before registration) | Completion of the preprint and arXiv pass. Poziotinib added as a positive control and lapatinib's negative-control basis strengthened (Ippolitov 2024); three replication cohorts added; target compartment recorded and the naked-antibody modality prior conditioned on it (Cheng 2025); benchmark differentiation tabulated after reading BioStudyBench, DBench-Bio and HindSight. No hypothesis direction was changed and no analysis had been run, so these are pre-registration amendments rather than deviations. |
| 2026-10-13 | Decisions 34, 35, 36 added (before registration) | A pilot run of the pipeline showed that Pareto-front membership is not invariant to the target-collapsing decisions, so the sensitivity analysis is now part of the locked protocol rather than an optional check. The Tier 3 rules were tightened at the same time, after the pilot could not retrieve structures for 13 agents and the honest response was to leave them unresolved. No hypothesis direction changed and no confirmatory analysis had been run. |
| 2026-10-17 | **Deviations, now that the analysis has run** | Four, each reported in the manuscript where it bites. **(a) Single-operator extraction, tiering and reachability curation.** Decisions 5 and 8 require two independent reviewers with tiering blind to reachability; one operator did all three. No Krippendorff α is reported, the tiering is provisional, and failure criterion 2 is therefore untestable in this analysis rather than passed. **(b) Tier 2 is matched on chemical structure, not on compound name.** The registered text did not specify a matching rule; name matching proved to under-recover B3DB (it missed "everolimus, rad001"), so matching is by InChIKey skeleton after salt stripping and charge neutralisation. The change is a correction and is one-directional: it can only resolve agents that name matching left unresolved. **(c) Tier 3 point estimates are the predicted probability, not a category midpoint.** Decision 12 mapped Tier 3 onto the same five categories as Tiers 1 and 2; probabilities are used instead, with the registered `predicted` band of ±0.30, so no prediction can land on a top grid point. **(d) Two registered components were not run and are reported as not run:** the HER2-stratified transcriptomic replication (no cohort was downloaded) and the leakage-controlled benchmark (no available model has a documented training cutoff preceding 1 January 2024, so running it would contaminate the split by construction). Failure criterion 3 is likewise untestable rather than passed. Failure criteria 1 and 4 were testable and both were passed: all five controls landed as required with the strict pair separating on non-overlapping bands, and 30 of 42 triples resolve at Tiers 1–3. |

## 9. What the first full run actually found (17 October 2026)

Recorded here so that the registered protocol and its outcome sit in one place.

- **Controls: pass.** lapatinib [0.08, 0.32] against tucatinib [0.72, 0.88] and poziotinib [0.68, 0.92].
- **Coverage:** 30 of 42 target–modality–agent triples and 12 of 23 nodes resolve. The 12 unresolved triples are 8 targets with no clinical-grade agent and 4 biologics outside the Tier 3 domain.
- **Pareto front:** PI3K–MTOR (Tier A, 0.65, bootstrap membership 73%), IL34–CSF1R (Tier B, 0.80, 39%), ERBB3–NRG1 (Tier B, 0.80, 36%). CDK4/6 reaches Tier A and is dominated at 0.35 (13%).
- **Collapsing sensitivity:** ERBB3–NRG1 and IL34–CSF1R hold the front in every scenario; PI3K–MTOR holds under each individual split but not when no collapsing is applied, where CDK4/6 replaces it.
- **Tier 3 predictor:** held-out ROC AUC 0.840 against 0.750 for the untrained rule (difference +0.090, 95% CI +0.071 to +0.108), accuracy 0.807, Brier 0.146. These figures supersede the ones first recorded here (0.834 against 0.739, difference +0.095 [+0.076, +0.114], accuracy 0.795, Brier 0.150), which came from a model fitted on the whole of B3DB including the eleven records that match agents Tier 3 is applied to. The pre-submission audit of 17 October 2026 caught the mismatch: the discrimination reported must belong to the model that produced the predictions, so `src/tier3_predictor.py` now takes the structures file and withholds those records before the split (`rule tier3_validation_v2`, `results/tier3_validation_v2.md`), and `src/make_figures_extra.py` withholds the same records before drawing Figure 5. The correction moves the point estimates slightly up and the difference slightly down; it changes no conclusion.
