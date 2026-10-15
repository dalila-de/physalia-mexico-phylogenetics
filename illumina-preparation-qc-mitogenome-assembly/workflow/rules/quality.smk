"""Assess trimmed reads; screening does not remove reads from marker assembly.

assess_reads_fastqc, classify_reads_kraken2, screen_human_bwa -> reports.
"""

rule assess_reads_fastqc:
    input:
        r1=WORK + "/trimmed/{sample}.R1.trimmed.fastq",
        r2=WORK + "/trimmed/{sample}.R2.trimmed.fastq"
    output:
        directory(WORK + "/fastqc/{sample}")
    threads: 2
    resources:
        mem_mb=2000
    log: WORK + "/logs/fastqc/{sample}.log"
    shell:
        "mkdir -p {output:q}; fastqc -t {threads} {input.r1:q} {input.r2:q} "
        "--outdir {output:q} > {log:q} 2>&1"

rule classify_reads_kraken2:
    input:
        r1=WORK + "/trimmed/{sample}.R1.trimmed.fastq",
        r2=WORK + "/trimmed/{sample}.R2.trimmed.fastq",
        db=config["kraken_db"] + "/hash.k2d"
    output:
        WORK + "/kraken/{sample}.report.txt"
    params: db=config["kraken_db"]
    threads: 8
    resources:
        mem_mb=100000
    log: WORK + "/logs/kraken/{sample}.log"
    shell:
        "kraken2 --threads {threads} --db {params.db:q} --output /dev/null "
        "--report {output:q} --use-names --report-minimizer-data "
        "--paired {input.r1:q} {input.r2:q} 2> {log:q}"

rule screen_human_bwa:
    input:
        r1=WORK + "/trimmed/{sample}.R1.trimmed.fastq",
        r2=WORK + "/trimmed/{sample}.R2.trimmed.fastq",
        reference=config["human_bwa"]
    output:
        WORK + "/human/{sample}.mapped_paired"
    threads: 8
    resources:
        mem_mb=16000
    log: WORK + "/logs/bwa_human/{sample}.log"
    shell:
        r"""
        bwa mem -t {threads} {input.reference:q} {input.r1:q} {input.r2:q} 2> {log:q} |
          samtools view -q 30 -F 0x900 -f 0x0002 -h - |
          awk 'substr($0,1,1)=="@" || ($9>=200 && $9<=1000) || ($9<=-200 && $9>=-1000)' |
          samtools view -c - > {output:q}
        """
