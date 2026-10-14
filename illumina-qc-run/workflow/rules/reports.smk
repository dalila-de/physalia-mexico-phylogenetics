"""Validate completed tools, retain compact products, and summarize the cohort.

report_sample preserves FASTAs, reports, logs and checksums in project storage.
report_cohort requires a valid summary for every selected specimen.
"""

rule report_sample:
    input:
        fastqc=WORK + "/fastqc/{sample}",
        kraken=WORK + "/kraken/{sample}.report.txt",
        human=WORK + "/human/{sample}.mapped_paired",
        mito=WORK + "/getorganelle/{sample}_mtgenome",
        shark=WORK + "/shark_results/{sample}",
        rego=WORK + "/rego_results/{sample}",
        trim=WORK + "/logs/trimmomatic/{sample}.log",
        counts=WORK + "/reads/{sample}.json",
        script=str(REPO / "scripts/retain_results.py"),
        reporter=str(REPO / "scripts/qc_report.py")
    output:
        RESULTS + "/{sample}/summary.json"
    params:
        work=WORK,
        destination=RESULTS,
        maximum=config["max_pairs"]
    resources:
        mem_mb=2000
    shell:
        "python {input.script:q} sample --work {params.work:q} --destination {params.destination:q} "
        "--sample {wildcards.sample:q} --maximum {params.maximum}"

rule report_cohort:
    input:
        reports=expand(RESULTS + "/{sample}/summary.json", sample=SAMPLES),
        script=str(REPO / "scripts/retain_results.py")
    output:
        RESULTS + "/sample_summary.tsv"
    params: samples=str(REPO / config["samples"])
    resources:
        mem_mb=2000
    shell:
        "python {input.script:q} cohort --manifest {params.samples:q} --output {output:q} "
        "--reports {input.reports:q}"
