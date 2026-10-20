# Mitochondrial seed used for the Mexican Illumina assemblies

`siph_mitochondrion_ncbi.fasta` is an unchanged copy of the nine-sequence
GetOrganelle seed used for the five Progreso specimens. It is a small seed
reference set supplied through `-s`, separate from the GetOrganelle database
configured through `GETORG_PATH`.

The source recorded in `../../provenance/illumina-qc/references.tsv` is:

```text
/gpfs/gibbs/pi/dunn/sc2962/20241018_IlluminaQC/siph_mitochondrion_ncbi.fasta
```

The copied file is 134,470 bytes, with the same SHA256 as the historical
execution record:

```text
4d6bae1baae57862995925ba4515b9cccba467f7f1f0be11c322bee618436b78
```

## Contents

| Original record ID | Taxon or sample label in the supplied file | Length (bp) |
|---|---|---:|
| NC_080952.1 | *Craseoa lathetica* | 15,081 |
| Nanomia_bijuga_NC | *Nanomia bijuga* | 14,625 |
| Nanomia_sp_California | *Nanomia* sp., California label | 15,414 |
| Agalma_clausi | *Agalma clausi* | 14,582 |
| Agalma_elegans_Atlantic | *Agalma elegans*, Atlantic label | 15,105 |
| Agalma_okeni | *Agalma okenii* | 14,560 |
| Agalma_elegans_Pacific | *Agalma elegans*, Pacific label | 15,054 |
| KT809335.1 | *Rhizophysa eysenhardtii* | 13,915 |
| KT809328.1 | *Physalia physalis* | 14,825 |

These are nine sequence records, not nine distinct species. Geographic labels
are reproduced from record IDs; they have not been independently verified.
Only three headers provide explicit accession identifiers. The other six
records retain their original descriptive IDs; accession mappings and the
original compilation method have not been established. The filename does not
establish accession provenance for every sequence.

The supplied headers describe NC_080952.1 as a complete mitochondrial genome
and KT809335.1 and KT809328.1 as partial mitochondrial genomes. Completeness of
the remaining records has not been established here. All sequences, headers,
and line wrapping are preserved. `sequences.tsv` records their original headers,
lengths, and sequence hashes.

## Reuse and verification

For a rerun, set `mito_seed` in the execution configuration to the absolute path
of this bundled FASTA. Historical configuration and checksum manifests retain
the original cluster path. The configured GetOrganelle database is still a
separate dependency.

From this directory, verify the exact file:

```bash
sha256sum --check sha256sum.txt
```

The spatiotemporal repository at commit
`d08f0edf4d819e21c85705e4965de71673808dec` instead configures
`siph_NCBI_RefSeq_mtgenome.fasta`. The file currently at that upstream cluster
path contains 22 records and differs from this seed. Do not substitute it when
reproducing the Mexican assemblies.
