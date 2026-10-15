# Mexican Physalia in a published phylogenetic context

Illumina read QC, mitochondrial annotation, and marker/mitogenome phylogenies
for five Mexican specimens alongside published Physalia sequences.

The five new specimens are **YPM-IZ-116693–116697**. The Illumina workflow and bundled QC outputs cover only these five specimens.
The trees use published sequences for those 35 Atlantic specimens, the five Mexican
assemblies, six published P. utriculus comparators, and available Rhizophysa
outgroups. COI and 16S also include the Caribbean Sanger specimen CH887.
There is one sequence per specimen per analysis and one retained tree analysis
per marker; the full and focal figures are views of that same tree.

![Mitogenome tree, focal view](analyses/results/genetic_groups/mitogenome_submitted.focal.png)

## Explore the results

| Tree | Full tree with references/outgroups | Focal view |
|---|---|---|
|16S|[PDF](analyses/results/genetic_groups/16S.pdf)|[PDF](analyses/results/genetic_groups/16S.focal.pdf)|
|18S|[PDF](analyses/results/genetic_groups/18S.pdf)|[PDF](analyses/results/genetic_groups/18S.focal.pdf)|
|COI|[PDF](analyses/results/genetic_groups/CO1.pdf)|[PDF](analyses/results/genetic_groups/CO1.focal.pdf)|
|ITS|[PDF](analyses/results/genetic_groups/ITS.pdf)|[PDF](analyses/results/genetic_groups/ITS.focal.pdf)|
|Mitogenome|[PDF](analyses/results/genetic_groups/mitogenome_submitted.pdf)|[PDF](analyses/results/genetic_groups/mitogenome_submitted.focal.pdf)|

PNG and editable SVG versions are alongside each PDF. E/W colors represent
**prior genome-wide genetic assignments**, not classifications inferred from
these trees. A dagger marks the three known vagrants; Mexican specimens are
pink and labeled **Gulf of Mexico (Progreso)**. CH887 is labeled
**Caribbean Sea (Puerto Morelos)** and shown in purple, without an E/W assignment. Filled/hollow internal-node dots show bootstrap
support ≥95/70–94. All original support values remain in the Newick files.

ITS lacks a Rhizophysa sequence in the reference panel and is midpoint-rooted
for display. The mitochondrial tree excludes YPM-IZ-110972 because its published
repository assembly has unresolved alternative paths and no deposited mitogenome.
Focal views only prune comparison/outgroup taxa for display; they do not rerun
inference or change the retained subtree's branch lengths.

## Repository layout

| Directory | Contents |
|---|---|
|[illumina-qc-run](illumina-qc-run/)|Five-specimen Mexican workflow, manifests, environment records, compact QC reports, sharkmer products, and mitochondrial assemblies|
|[mitogenome-annotation](mitogenome-annotation/)|Five Mexican assemblies, original MITOS2 FASTA/GFF3/gene-order annotations, and a summary|
|[analyses](analyses/)|Selected sequences, alignments, IQ-TREE outputs, genetic-group annotations, figures, and public-sequence audit|
|[references/church](references/church/)|Pinned published reference panels, specimen metadata, and upstream license/attribution|
|[sanger-2026-10-14](sanger-2026-10-14/)|Original chromatograms, trace QC and processed CH887 sequences; combined trees are in analyses|
|[provenance](provenance/)|Historical execution settings, tool/environment records, input identities, and artifact checksums|

The Illumina manifests and outputs contain only the five Mexican specimens.
Their results were selected from the completed 40-specimen run without rerunning
assembly. Historical execution records retain their original scope; the complete
pre-cleanup snapshot is preserved outside this repository, as recorded in
[selection provenance](provenance/mexico-only-qc.json).

This is a curated snapshot of the completed analyses. Raw reads, large screening
reports, reference databases, software environments, job working directories,
and superseded exploratory tree runs are not bundled. The separate
Super-waster project is outside this repository.

## Recreate the figures

```bash
conda env create -f environment.yml
conda activate physalia-trees
snakemake --snakefile analyses/plots.smk --cores 2 --forceall
python scripts/validate_release.py
```

All files needed for plotting are included; raw-read storage and MITOS references
are not needed. The environment recipe is supplied; a clean environment rebuild
has not been tested. The figures were checked with the recorded existing Python,
Biopython and matplotlib environment.

For full inference, see [analysis methods](analyses/README.md). For QC, see
[the QC workflow](illumina-qc-run/README.md); full raw-read reruns require access
to the FASTQs and installation/configuration of the site-specific tools and
reference databases. Historical absolute paths in manifests/provenance identify
what was used; update execution configs when running elsewhere.

## Data availability and interpretation

The [public-sequence audit](analyses/results/README.md) and
[accession table](analyses/results/specimen_marker_availability.tsv) distinguish
specimen identity from sequence matches. Existing marker/mitogenome records were
checked against the pinned repository and a dated NCBI snapshot. Missing voucher
matches or non-identical assemblies alone do not establish biological novelty.
MITOS outputs are automated annotations, including fragmented features, and
have not been manually curated for GenBank submission.

Published reference data are from
[Church's Physalia repository](https://github.com/shchurch/Physalia_spatiotemporal_patterns),
commit `d08f0edf4d819e21c85705e4965de71673808dec`. See
[third-party attribution](references/church/README.md) and its preserved LICENSE.


## Caribbean Sanger specimen

The [CH887 Sanger analysis](sanger-2026-10-14/README.md) adds COI and 16S from one
Caribbean Sea specimen, with raw chromatograms, quality reports and reproducible
sequence processing. The two amplicons initially supplied as COI were identified
as different mitochondrial markers. CH887 is shown separately from genetic E/W.
The single COI and 16S trees linked above include CH887, the published reference
panels and all five Mexican specimens. Superseded WGS-only COI/16S results were
archived outside this repository; [consolidation provenance](provenance/tree-consolidation.json)
records their replacement by the previously completed combined analyses.

## Validation

With the analysis environment activated, run from the repository root:

```bash
python -m pytest -q
ruff check analyses/scripts illumina-qc-run/scripts sanger-2026-10-14/scripts
python scripts/validate_release.py
python sanger-2026-10-14/scripts/validate_sanger.py
```

The repository tracks README.md documentation, workflow code, selected scientific
inputs and results, environment recipes, and provenance. Other Markdown files,
local configuration, caches, bulk read files, temporary files, scheduler output,
and dry-run transcripts are excluded. Full execution requires the site resources
described in each workflow README.
