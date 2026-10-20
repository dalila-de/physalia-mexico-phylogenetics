# Reference identities for the Mexican Illumina analysis

This audit supplements the immutable execution manifest
[`../illumina-qc/references.tsv`](../illumina-qc/references.tsv).
It identifies the references behind the configuration without changing
sequences, historical paths, or recorded hashes.

## GetOrganelle

The animal mitochondrial **LabelDatabase and SeedDatabase both report version
0.0.1**. Copies of their `VERSION` records are included here. Their version-file
and FASTA SHA256 values match the historical execution manifest, and each
`VERSION` record also gives the matching FASTA checksum.

The custom nine-sequence siphonophore seed supplied with `-s` is separate from
these installed databases. Its exact original bytes are now bundled in
[`../../references/mitochondrial-seed/`](../../references/mitochondrial-seed/).
Original accession and compilation-provenance limits are documented there.
No newer upstream seed was substituted.

The executed command used:

```text
-F animal_mt -R 10 -k 21,45,65,85,105
--reduce-reads-for-coverage inf --max-reads inf
```

Both paired trimmed-read files were supplied, with the custom FASTA through
`-s` and the installation's separate database directory through `GETORG_PATH`.
The existing configuration retains its historical cluster seed path; to run
elsewhere, point `mito_seed` at the absolute path of the bundled copy.

## Kraken2

The configured site path `/gpfs/gibbs/data/kraken2/standard` currently resolves
to `/gpfs/gibbs/data/db/kraken2/k2_standard_20240904`. Thus the identifiable
release label is **20240904 (4 September 2024)**, not the file modification date
or a newly inferred exact build timestamp.

The `opts.k2d` and `taxo.k2d` checksums were recomputed and match the execution
manifest. The 85.9-GB `hash.k2d` table was not reread during this audit: its
recorded execution SHA256 remains the identity evidence, and its current size
and modification time match the historical record. The release-label
association is based on the current alias target plus those checks; it is not
a separately archived historical symlink observation.

The complete component identities and verification scope are in
`reference-identities.json`. These large external database files are not
distributed with this repository. An installation carrying the same release
label should still be checked against the recorded hashes before reuse.

## Human screening reference

The executed `genome.fa` and its BWA `.ann` file both match their historical
SHA256 values. The accession remains **GCA_000001405.15 (GRCh38)**, with the
following important specification: the local FASTA represents the primary
assembly and mitochondrial sequence, plus an **EBV decoy**; it is not the
entire assembly package including alternate loci.

All primary-assembly and non-nuclear accession/length pairs in the official
NCBI report are represented in the local BWA index. The extra record is
`chrEBV`, accession `AJ507799.2`, labeled `Human_herpesvirus_4` and `decoy`
in the original header. The supplementary reference inventory records each
sequence and its assembly unit.

The NCBI assembly report is preserved as `GRCh38-assembly-report.txt`, retrieved
from the [official assembly directory](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_report.txt).
`human-genome.fa.ann` is the matching local index metadata. This audit verifies
historical file identity plus accession/length membership; it does not claim
a fresh base-by-base comparison against a downloaded NCBI genome.

The screening command counted qualifying alignments against this entire
reference, including its decoy. The recorded count is zero in all five
specimens, so the presence of the decoy does not change that observed result.

## Verification

`reference-identities.json` records source paths, hashes, and the scope of each
check. `checksums.tsv` covers the new audit files and bundled seed without
rewriting the original execution or release manifests. Paths in this companion
manifest are relative to the repository root.

The database release labels do not replace the immutable checksums. No
database, human genome, index, or raw-read bulk files are added to git.
