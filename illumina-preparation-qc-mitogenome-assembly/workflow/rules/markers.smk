"""Recover mitochondrial assemblies and primer-defined marker sequences.

assemble_mitogenome_getorganelle and recover_markers_sharkmer feed report_sample.
Both consume trimmed reads.
"""

rule assemble_mitogenome_getorganelle:
    input:
        r1=WORK + "/trimmed/{sample}.R1.trimmed.fastq",
        r2=WORK + "/trimmed/{sample}.R2.trimmed.fastq",
        seed=config["mito_seed"]
    output:
        directory(WORK + "/getorganelle/{sample}_mtgenome")
    params: db=config["getorganelle_db"]
    threads: 4
    resources:
        mem_mb=32000
    log: WORK + "/logs/getorganelle/{sample}.log"
    shell:
        "GETORG_PATH={params.db:q} get_organelle_from_reads.py "
        "-1 {input.r1:q} -2 {input.r2:q} -t {threads} -R 10 -k 21,45,65,85,105 "
        "-F animal_mt -s {input.seed:q} --overwrite --reduce-reads-for-coverage inf "
        "--max-reads inf -o {output:q} > {log:q} 2>&1"

rule recover_markers_sharkmer:
    input:
        r1=WORK + "/trimmed/{sample}.R1.trimmed.fastq",
        tool=config["sharkmer"]
    output:
        directory(WORK + "/shark_results/{sample}")
    threads: 8
    resources:
        mem_mb=160000
    log:
        out=WORK + "/logs/sharkmer/{sample}.out",
        err=WORK + "/logs/sharkmer/{sample}.err"
    shell:
        "mkdir -p {output:q}; {input.tool:q} {input.r1:q} -k 31 -n 100 -t {threads} "
        "-s {wildcards.sample:q} -o {output:q} --pcr cnidaria --pcr metazoa --pcr bacteria "
        "> {log.out:q} 2> {log.err:q}"
