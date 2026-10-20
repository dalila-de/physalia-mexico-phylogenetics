"""Download NCBI records and prepare reference plus local sequence panels."""
rule download_ncbi_references:
    input:
        manifest=REFERENCE + "/accessions.tsv",
        script=str(HERE / "scripts/download_references.py")
    output:
        REFERENCE + "/records.gb.gz"
    shell:
        "{PYTHON:q} {input.script:q}"

rule prepare_tree_sequences:
    input:
        records=REFERENCE + "/records.gb.gz",
        manifest=REFERENCE + "/accessions.tsv",
        script=str(HERE / "scripts/prepare_tree_inputs.py"),
        extractor=str(HERE / "scripts/ncbi_references.py"),
        common=str(HERE / "scripts/common.py"),
        cohort=str(HERE / "inputs/cohort.tsv"),
        local=sorted(str(p) for p in (ROOT / "illumina-preparation-qc-mitogenome-assembly/results" / config["qc_run"]).glob("*/*/*.fasta")),
        annotations=sorted(str(p) for p in (ROOT / "mitogenome-annotation/results").glob("*/result.fas"))
    output:
        fasta=INPUTS + "/{gene}.fasta",
        tips=INPUTS + "/{gene}.tips.tsv"
    shell:
        "{PYTHON:q} {input.script:q} {wildcards.gene:q}"

rule combine_sanger:
    input:
        fasta=expand(INPUTS + "/{gene}.fasta", gene=["CO1", "16S"]),
        tips=expand(INPUTS + "/{gene}.tips.tsv", gene=["CO1", "16S"]),
        sanger=expand(str(ROOT / "sanger-2026-10-14/inputs/{gene}.sanger.fasta"), gene=["CO1", "16S"]),
        script=str(HERE / "scripts/combine_sanger.py"),
        common=str(HERE / "scripts/common.py")
    output:
        fasta=expand(INPUTS + "/combined/{gene}.fasta", gene=["CO1", "16S"]),
        tips=expand(INPUTS + "/combined/{gene}.tips.tsv", gene=["CO1", "16S"])
    shell:
        "{PYTHON:q} {input.script:q}"
