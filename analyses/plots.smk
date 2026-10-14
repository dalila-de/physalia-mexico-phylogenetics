from pathlib import Path
HERE = Path(workflow.basedir).resolve()
GENES = ["16S", "18S", "CO1", "ITS", "mitogenome_submitted"]
rule plot:
    input:
        script=str(HERE / "scripts/plot_genetic_groups.py"),
        tree=str(HERE / "results/figures/{gene}.rooted.nwk"),
        tips=str(HERE / "inputs/{gene}.tips.tsv"),
        annotations=str(HERE / "inputs/genetic_annotations.tsv")
    output:
        expand(str(HERE / "results/genetic_groups/{{gene}}{view}.{ext}"), view=["", ".focal"], ext=["pdf", "png", "svg"])
    shell:
        "python {input.script:q} {wildcards.gene:q}"
rule all:
    input:
        expand(str(HERE / "results/genetic_groups/{gene}.pdf"), gene=GENES)
    default_target: True
