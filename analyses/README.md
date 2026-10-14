# Mexican and published specimen phylogenetics

The four marker analyses use 35 published Atlantic specimens and five Mexican
specimens. Six published P. utriculus specimens provide comparison; Rhizophysa
is the outgroup wherever the source panel provides it. Each specimen appears
once. The mitochondrial tree excludes YPM-IZ-110972 because its repository
mitogenome has unresolved alternative paths and was not deposited.

| Analysis | Atlantic published | Mexico | utriculus | Rhizophysa | Total |
|---|---:|---:|---:|---:|---:|
|16S|35|5|6|13|59|
|18S|35|5|6|2|48|
|CO1|35|5|6|10|56|
|ITS|35|5|6|0|46|
|Mitogenome|34|5|6|3|48|

MAFFT 7.526 aligns the marker sequences de novo with direction adjustment;
mitogenomes are added to the paper's aligned backbone. IQ-TREE 3.0.1 uses
ModelFinder/BIC, 1000 ultrafast bootstraps, five searches, and seed 20261012.
For Mexican 16S/CO1, select the longest individual MITOS rrnL/cox1 fragment or
canonical sharkmer product, with MITOS winning ties. 18S/ITS select the longest
canonical sharkmer product. Fragmented genes are never concatenated for length.
Candidate selection is recorded in `inputs/*.candidates.tsv`.

## Figures

`results/genetic_groups/` contains PDF, SVG and PNG figures. Full plots retain
all taxa; `.focal` plots prune utriculus/Rhizophysa only for display, preserving
focal branch lengths and topology. The original inference is unchanged.
Filled internal-node dots indicate bootstrap ≥95; hollow dots indicate 70–94;
values below 70 are unmarked. Exact values remain in Newick trees.
Rhizophysa roots the reported non-ITS trees; ITS is midpoint-rooted for display.

Blue/orange E/W groups come from the **prior genome-wide NGSadmix K2 results**,
not these gene trees or collection geography. Mexico is pink and unassigned here.
The three vagrants have a dagger:110844 W-like,111012 E-like,111018 predominantly
E-like (94.73% eastern component in the whole-assembly analysis). Frozen evidence
and annotations are in `provenance/genetic_groups/` and
`inputs/genetic_annotations.tsv`. Coefficients are not calibrated historical
ancestry proportions. Colored specimen labels and collection locations sit directly beside branch tips.

## Reproduce

From the repository root, use the top-level environment.yml:

```bash
snakemake --snakefile analyses/plots.smk --cores 2 --forceall
```

This needs only the files shipped here and recreates the figures. To rerun the
full annotation/alignment/inference chain, configure the MITOS environment and
refseq89m database in `config/config.yaml`, set `iqtree` to an executable path,
and run `snakemake all --snakefile analyses/Snakefile --cores 24`. Full runs
consume bundled reference panels and QC assemblies; no raw reads are needed.
Plotting does not depend on access to raw FASTQs or MITOS databases.

## Public-sequence audit

See [availability findings](results/README.md),
[specimen-by-marker accession table](results/specimen_marker_availability.tsv),
and [all sharkmer matches](results/all_sharkmer_product_matches.tsv).
The dated NCBI record snapshot is retained as `provenance/ncbi_records.gb.gz`;
`python analyses/scripts/audit_availability.py` rereads that frozen snapshot.
An absent voucher match does not prove a sequence is novel. Five Mexican
mitogenomes lack exact matches in the queried datasets; confirmation requires
assembly/read review. Sequence-identity relationships and specimen matches are
kept separate. The availability audit covers the 35 published specimens and five Mexican
specimens; local sequence comparisons use only the five Mexican assemblies.
`inputs/cohort.tsv` selects the phylogenetic cohort independently of the
five-specimen Illumina manifest.
