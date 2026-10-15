"""Validate/subset paired inputs, then trim adapters and low-quality bases.

subset_paired_reads -> trim_reads_trimmomatic -> quality and marker rules.
Two-run specimens contribute an equal maximum from each run.
"""

rule subset_paired_reads:
    input:
        reads=lambda w: [r[k] for r in READS[w.sample] for k in ("R1", "R2")],
        script=str(REPO / "scripts/subset_reads.py")
    output:
        r1=temp(WORK + "/reads/{sample}_sub_R1.fastq"),
        r2=temp(WORK + "/reads/{sample}_sub_R2.fastq"),
        counts=WORK + "/reads/{sample}.json"
    params:
        manifest=str(REPO / config["reads"]),
        maximum=config["max_pairs"]
    resources:
        mem_mb=2000
    shell:
        "python {input.script:q} --manifest {params.manifest:q} --sample {wildcards.sample:q} "
        "--maximum {params.maximum} --r1 {output.r1:q} --r2 {output.r2:q} --counts {output.counts:q}"

rule trim_reads_trimmomatic:
    input:
        r1=WORK + "/reads/{sample}_sub_R1.fastq",
        r2=WORK + "/reads/{sample}_sub_R2.fastq",
        adapters=str(REPO / "config/adapters.fa")
    output:
        r1=temp(WORK + "/trimmed/{sample}.R1.trimmed.fastq"),
        r2=temp(WORK + "/trimmed/{sample}.R2.trimmed.fastq"),
        u1=temp(WORK + "/trimmed/{sample}.R1.unpaired.fastq"),
        u2=temp(WORK + "/trimmed/{sample}.R2.unpaired.fastq"),
        report=WORK + "/logs/trimmomatic/{sample}.log"
    threads: 8
    resources:
        mem_mb=16000
    shell:
        "java -Xmx12g -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads {threads} "
        "{input.r1:q} {input.r2:q} {output.r1:q} {output.u1:q} {output.r2:q} {output.u2:q} "
        "ILLUMINACLIP:{input.adapters}:2:30:10:2:True LEADING:15 TRAILING:15 MINLEN:50 "
        "2> {output.report:q}"
