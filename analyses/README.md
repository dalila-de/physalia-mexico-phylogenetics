# Mexican and published specimen phylogenetics

The four marker analyses use 35 published Atlantic specimens and five Mexican
specimens. Six published P. utriculus specimens provide comparison; Rhizophysa
is the outgroup wherever the source panel provides it. Each specimen appears
once. COI and 16S also include one Caribbean Sanger specimen (CH887).
These are the only retained COI and 16S tree analyses. The mitochondrial tree excludes YPM-IZ-110972 because its mitogenome was not
included in the [reference repository](https://github.com/shchurch/Physalia_spatiotemporal_patterns)
following concerns about its assembly.

| Analysis | Atlantic published | Mexico | utriculus | Rhizophysa | Sanger | Total |
|---|---:|---:|---:|---:|---:|---:|
|16S|35|5|6|13|1|60|
|18S|35|5|6|2|0|48|
|COI|35|5|6|10|1|57|
|ITS|35|5|6|0|0|46|
|Mitogenome|34|5|6|3|0|48|

MAFFT 7.526 aligns the marker sequences de novo with direction adjustment;
mitogenomes are added to the reference repository's aligned backbone. IQ-TREE 3.0.1 uses
ModelFinder/BIC, 1000 ultrafast bootstraps, five searches, and base seed 20261012 (incremented by 1000 for each subsequent search).
For Mexican 16S/CO1, select the longest individual MITOS rrnL/cox1 fragment or
canonical sharkmer product, with MITOS winning ties. 18S/ITS select the longest
canonical sharkmer product. Fragmented genes are never concatenated for length.
Candidate selection is recorded in `inputs/*.candidates.tsv`.
In the retained selection, MITOS2 supplied COI and 16S for all five Progreso
specimens; canonical `cnidaria_18S` and `cnidaria_ITS` sharkmer products supplied
18S and ITS. Other recovered primer products were not added as duplicate tips.
The MITOS2 annotations used version 2.1.9 with `refseq89m`, code 4, and
`--linear --orih 0 --oril 0 --trna 0 --intron 0`.
The WGS/reference FASTAs in `inputs/` are intermediate panels. Final COI/16S
inputs in `inputs/combined/` add the quality-controlled Sanger sequences; all
alignments, inferred trees and figures are together under `results/`.
See [Sanger processing](../sanger-2026-10-14/README.md) for trace evidence.

Marker alignment uses `--auto --adjustdirectionaccurately`; mitogenome addition
uses `--add --adjustdirectionaccurately`. Validation normalizes direction prefixes,
checks unique IDs, exact input membership and the nucleotide/IUPAC alphabet,
and removes columns containing only `-` or `?`. None were removed from the
five retained alignments. ModelFinder selected models under the Bayesian
information criterion (BIC). Outgroups are *Rhizophysa filiformis* and
*R. eysenhardtii*; the table counts accession records, not verified individuals.

## Figures

`results/genetic_groups/` contains PDF, SVG and PNG figures of rooted bootstrap
consensus trees (`.contree`), rather than the separate maximum-likelihood trees. Full plots retain
all taxa; `.focal` plots prune utriculus/Rhizophysa only for display, preserving
focal branch lengths and topology. The existing combined COI/16S inference and the other three inferences are preserved.
Full and focal figures are two views of each single analysis, not independent trees.
Filled internal-node dots indicate bootstrap ≥95; hollow dots indicate 70–94;
values below 70 are unmarked. Exact values remain in Newick trees.
Rhizophysa roots the reported non-ITS trees; ITS is midpoint-rooted for display.

Blue/orange E/W groups come from the **prior genome-wide NGSadmix K2 results**,
not these gene trees or collection geography. The five Gulf of Mexico (Progreso) specimens are pink and unassigned here;
CH887, Caribbean Sea (Puerto Morelos), is purple and also unassigned.
The three vagrants have a dagger:110844 W-like,111012 E-like,111018 predominantly
E-like (94.73% eastern component in the whole-assembly analysis). Frozen evidence
and annotations are in `provenance/genetic_groups/` and
`inputs/genetic_annotations.tsv`. Coefficients are not calibrated historical
ancestry proportions. Colored specimen labels and collection locations sit directly beside branch tips.

## Interpretation of the retained placements

The new specimens group with published *P. physalis* references, consistent with
that identification. This is not definitive species confirmation or a genome-wide
E/W assignment. In the rooted bootstrap consensus trees, the clade containing all
published *P. physalis* references and the new specimens has support of 36% (COI),
58% (16S), 75% (ITS), and 77% (mitogenomes). The corresponding 18S node has no
recorded support value; a missing value is not a zero-support estimate.
CH887's immediate groupings have 14% COI and 30% 16S support. These values describe
particular branches, not probabilities of species identity. Its COI sequence is
partial and forward-only. Short or near-zero branches and weak groupings limit
within-species resolution. Some mitochondrial subclades have stronger support,
but mitochondrial placement alone does not establish genome-wide ancestry.

Support values are retained in `results/trees/*.contree` and their rooted display
copies in `results/figures/*.rooted.nwk`. Nearest-reference reports use the separate
maximum-likelihood trees and should not be substituted for consensus branch support.
The five mitochondrial assemblies have automated, fragmented MITOS2 annotations;
manual annotation validation and genome completeness have not been established.

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
