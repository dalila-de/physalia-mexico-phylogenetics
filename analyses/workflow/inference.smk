"""Align all five panels de novo; infer, validate, root and plot trees."""
rule align_sequences_mafft:
    input:
        fasta=lambda w: str(HERE / "inputs" / ("combined" if w.gene in {"16S", "CO1"} else "") / f"{w.gene}.fasta"),
        validator=str(HERE / "scripts/validate_alignment.py")
    output:
        OUT + "/alignments/{gene}.fasta"
    params:
        raw=WORK + "/alignments/{gene}.raw.fasta"
    threads: 8
    resources: mem_mb=16000
    log: WORK + "/logs/mafft/{gene}.log"
    shell:
        "mkdir -p $(dirname {params.raw:q}) $(dirname {log:q}); "
        "mafft --auto --adjustdirectionaccurately --thread {threads} {input.fasta:q} > {params.raw:q} 2> {log:q}; "
        "{PYTHON:q} {input.validator:q} {params.raw:q} {input.fasta:q} {output:q}"

rule infer_phylogeny_iqtree:
    input:
        alignment=OUT + "/alignments/{gene}.fasta"
    output:
        tree=OUT + "/trees/{gene}.treefile",
        report=OUT + "/trees/{gene}.iqtree",
        consensus=OUT + "/trees/{gene}.contree"
    params:
        prefix=OUT + "/trees/{gene}",
        seed=config["seed"],
        runs=lambda w: config.get("search_runs_by_marker", {}).get(w.gene, config["search_runs"])
    threads: 8
    resources: mem_mb=16000
    log: WORK + "/logs/iqtree/{gene}.log"
    shell:
        "mkdir -p $(dirname {params.prefix:q}) $(dirname {log:q}); "
        "{IQTREE:q} -s {input.alignment:q} -m MFP -B 1000 --runs {params.runs} "
        "-nt {threads} -seed {params.seed} -pre {params.prefix:q} > {log:q} 2>&1"

rule report_phylogeny:
    input:
        tree=OUT + "/trees/{gene}.treefile",
        inference=OUT + "/trees/{gene}.iqtree",
        consensus=OUT + "/trees/{gene}.contree",
        alignment=OUT + "/alignments/{gene}.fasta",
        tips=lambda w: str(HERE / "inputs" / ("combined" if w.gene in {"16S", "CO1"} else "") / f"{w.gene}.tips.tsv"),
        script=str(HERE / "scripts/report_tree.py")
    output:
        report=OUT + "/tree_reports/{gene}.json",
        rooted=OUT + "/figures/{gene}.rooted.nwk"
    shell:
        "{PYTHON:q} {input.script:q} {wildcards.gene:q}"

rule plot_genetic_groups:
    input:
        report=OUT + "/tree_reports/{gene}.json",
        tree=OUT + "/figures/{gene}.rooted.nwk",
        tips=lambda w: str(HERE / "inputs" / ("combined" if w.gene in {"16S", "CO1"} else "") / f"{w.gene}.tips.tsv"),
        annotation=INPUTS + "/genetic_annotations.tsv",
        script=str(HERE / "scripts/plot_genetic_groups.py")
    output:
        pdf=OUT + "/genetic_groups/{gene}.pdf",
        png=OUT + "/genetic_groups/{gene}.png",
        table=OUT + "/genetic_groups/{gene}.annotations.tsv",
        focal_pdf=OUT + "/genetic_groups/{gene}.focal.pdf",
        focal_png=OUT + "/genetic_groups/{gene}.focal.png"
    shell:
        "{PYTHON:q} {input.script:q} {wildcards.gene:q}"
