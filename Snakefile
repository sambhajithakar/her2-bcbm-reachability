# Two-axis prioritization in HER2+ breast cancer brain metastasis.
# Reproduces every computed result in the manuscript from the committed inputs.
# Run:  snakemake --cores 1 --use-conda

B3DB = "data/b3db-raw"   # self-contained in this archive; B3DB is CC0

rule all:
    input:
        "results/reachability_scores.csv",
        "results/control_test.md",
        "results/two_axis_nodes.csv",
        "results/figure2_reachability.png",
        "results/figure3_two_axis.png",
        "results/collapse_sensitivity.md",
        "results/tier3_validation.md",
        "results/reachability_scores_v2.csv",
        "results/structural_recheck.md",
        "results/two_axis_v2_nodes.csv",
        "results/control_test_v2.md",
        "results/collapse_sensitivity_v2.md",
        "results/v2/figure1_flow.png",
        "results/v2/figure2_reachability.png",
        "results/v2/figure3_two_axis.png",
        "results/v2/figure4_collapse.png",
        "results/v2/figure5_tier3.png",
        "results/tier3_validation_v2.md",

rule score_reachability:
    input:
        agents="data/agents.csv",
        reg=f"{B3DB}/B3DB_regression.tsv",
        cls=f"{B3DB}/B3DB_classification.tsv",
    output:
        "results/reachability_scores.csv"
    conda:
        "environment.yml"
    shell:
        "python3 -I src/score_reachability.py {input.agents} {input.reg} {input.cls} {output}"

rule validate_controls:
    input:
        "results/reachability_scores.csv"
    output:
        "results/control_test.md"
    conda:
        "environment.yml"
    shell:
        "python3 -I src/validate_controls.py {input} {output}"

rule two_axis_map:
    input:
        census="data/census.csv",
        target_map="data/target_map.csv",
        scores="results/reachability_scores.csv",
        agent_map="data/agent_target_map.csv",
    output:
        "results/two_axis_nodes.csv"
    conda:
        "environment.yml"
    shell:
        "python3 -I src/two_axis_map.py {input.census} {input.target_map} "
        "{input.scores} {input.agent_map} results/two_axis > results/two_axis_log.txt"

rule figures:
    input:
        scores="results/reachability_scores.csv",
        nodes="results/two_axis_nodes.csv",
    output:
        "results/figure2_reachability.png",
        "results/figure3_two_axis.png",
    conda:
        "environment.yml"
    shell:
        "python3 -I src/make_figures.py {input.scores} {input.nodes} results"

rule collapse_sensitivity:
    input:
        census="data/census.csv",
        target_map="data/target_map.csv",
        scores="results/reachability_scores.csv",
        agent_map="data/agent_target_map.csv",
    output:
        "results/collapse_sensitivity.md"
    conda:
        "environment.yml"
    shell:
        "python3 src/collapse_sensitivity.py {input.census} {input.target_map} "
        "{input.scores} {input.agent_map} {output}"

rule tier3_validation:
    input:
        f"{B3DB}/B3DB_classification.tsv"
    output:
        "results/tier3_validation.md"
    conda:
        "environment.yml"
    shell:
        "python3 -I src/tier3_predictor.py {input} {output}"

# --- 17 October 2026: structures retrievable, so Tier 2 is matched on
# --- chemical structure and Tier 3 can be applied. ---------------------

rule rescore_with_structures:
    input:
        agents="data/agents.csv",
        v1="results/reachability_scores.csv",
        structures="data/structures.csv",
        reg=f"{B3DB}/B3DB_regression.tsv",
        cls=f"{B3DB}/B3DB_classification.tsv",
    output:
        csv="results/reachability_scores_v2.csv",
        md="results/structural_recheck.md",
    conda:
        "environment.yml"
    shell:
        "python3 -I src/rescore_with_structures.py {input.agents} {input.v1} "
        "{input.structures} {input.reg} {input.cls} {output.csv} {output.md}"

rule two_axis_map_v2:
    input:
        census="data/census.csv",
        target_map="data/target_map.csv",
        scores="results/reachability_scores_v2.csv",
        agent_map="data/agent_target_map.csv",
    output:
        "results/two_axis_v2_nodes.csv"
    conda:
        "environment.yml"
    shell:
        "PYTHONPATH=src python3 src/two_axis_map.py {input.census} {input.target_map} "
        "{input.scores} {input.agent_map} results/two_axis_v2 > results/two_axis_v2_log.txt"

rule controls_v2:
    input: "results/reachability_scores_v2.csv"
    output: "results/control_test_v2.md"
    conda: "environment.yml"
    shell: "python3 -I src/validate_controls.py {input} {output}"

rule collapse_sensitivity_v2:
    input:
        census="data/census.csv",
        target_map="data/target_map.csv",
        scores="results/reachability_scores_v2.csv",
        agent_map="data/agent_target_map.csv",
    output: "results/collapse_sensitivity_v2.md"
    conda: "environment.yml"
    shell:
        "PYTHONPATH=src python3 src/collapse_sensitivity.py {input.census} "
        "{input.target_map} {input.scores} {input.agent_map} {output}"

rule figures_v2:
    input:
        scores="results/reachability_scores_v2.csv",
        nodes="results/two_axis_v2_nodes.csv",
        census="data/census.csv",
        collapse="results/collapse_sensitivity_v2.md",
        cls=f"{B3DB}/B3DB_classification.tsv",
        structures="data/structures.csv",
    output:
        "results/v2/figure1_flow.png",
        "results/v2/figure2_reachability.png",
        "results/v2/figure3_two_axis.png",
        "results/v2/figure4_collapse.png",
        "results/v2/figure5_tier3.png",
    conda:
        "environment.yml"
    shell:
        "mkdir -p results/v2 && "
        "PYTHONPATH=src python3 src/make_figures.py {input.scores} {input.nodes} results/v2 && "
        "PYTHONPATH=src python3 src/make_figures_extra.py {input.census} {input.collapse} "
        "{input.cls} results/v2 {input.structures}"

# The 13 October validation fitted on the whole of B3DB, including the records
# that match agents Tier 3 is applied to. This rule withholds those records, so
# the reported discrimination belongs to the model that made the predictions.
rule tier3_validation_v2:
    input:
        cls=f"{B3DB}/B3DB_classification.tsv",
        structures="data/structures.csv",
    output:
        "results/tier3_validation_v2.md"
    conda:
        "environment.yml"
    shell:
        "python3 -I src/tier3_predictor.py {input.cls} {output} {input.structures}"
