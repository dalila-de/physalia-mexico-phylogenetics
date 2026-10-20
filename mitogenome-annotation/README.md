# Mexican mitogenome annotations

Five specimens: YPM-IZ-116693 through YPM-IZ-116697. `inputs/` holds the
GetOrganelle assemblies used by MITOS2; `results/<sample>/` contains the original
MITOS nucleotide feature FASTA, GFF3 coordinates, and gene-order output.
`annotation_summary.tsv` summarizes the retained records.

MITOS2 2.1.9 used refseq89m and genetic code 4, with `--orih 0 --oril 0 --trna 0
--intron 0 --linear`. tRNAs and replication origins were not annotated in this
run. Each assembly has 15,038 bp. These are automated annotations, including split
gene fragments, not manually curated submission-ready annotations. A recovered
assembly is not by itself evidence of a complete/circular mitochondrial genome.

The annotation rule is `annotate_mitogenome_mitos` in `Snakefile`.
Install the separate MITOS environment and reference database, then configure
`mitos_env`, `mitos_refdir`, and `mitos_refname` in `../analyses/config/config.yaml`.
From the repository root, preview with
`snakemake all --snakefile mitogenome-annotation/Snakefile --cores 5 --dry-run`.
The NCBI tree rerun consumes the retained annotations without rerunning MITOS2.
The analysis accesses these outputs through the relative
`analyses/results/annotations` symlink. Historical environment/tool hashes are
in `../provenance/mitos-annotation/`.
