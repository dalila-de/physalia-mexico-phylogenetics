# NCBI reference sequences

`accessions.tsv` explicitly selects **140 distinct accession.version records**
for five analyses. `records.gb.gz` contains the direct NCBI GenBank responses;
`download.json` records download time, request URLs and response/file hashes.
[Specimen accession key](specimen_accessions.tsv) and
[Rhizophysa accessions by marker](rhizophysa_accessions.tsv) provide compact
views of the input manifest. Empty cells indicate that a specimen is not used
for that mitochondrial marker analysis.
No reference alignment is supplied or consumed. All five panels are aligned
de novo with MAFFT after adding the locally generated Mexican sequences.

| Analysis | P. physalis | P. utriculus | Rhizophysa accession records |
|---|---:|---:|---:|
| COI | 34 | 6 | 8 |
| 16S | 34 | 6 | 11 |
| 18S | 35 | 6 | 2 |
| ITS | 35 | 6 | 0 |
| Mitogenome | 34 | 6 | 3 |

The selected *Physalia* records are subsets of GenBank accessions
PZ224317–PZ224484 (mitochondria, 168 deposited records), PZ802573–PZ802753
(18S, 181), and PZ804143–PZ804328 (ITS, 186). The full range totals are not
the numbers of records used in these trees. Every selected voucher/isolate
and species is checked against the source feature of its GenBank record.

For whole mitogenomes and stand-alone nuclear/outgroup marker records, use the
entire deposited sequence. For COI/16S contained in mitochondrial records,
extract the longest non-pseudogene annotated cox1 CDS or 16S rRNA feature,
respecting strand and the deposited feature coordinates. Separate fragments
are not concatenated. Fuzzy/partial boundaries are retained as deposited.
Feature coordinates and sequence hashes are in each `analyses/inputs/*.tips.tsv`.
This uses the deposited 1,552-bp cox1 feature of YPM-IZ-104460 rather than the
1,549-bp sequence in the superseded analysis.

`excluded.tsv` records deliberate exclusions. YPM-IZ-110972 lacks a deposited
mitogenome in the specified data; it remains in the nuclear analyses only.
NC_080941.1 and NC_080942.1 are RefSeq copies of OQ957199.1 and OQ957206.1,
respectively, according to their NCBI comments. The original accessions are
used once per marker, avoiding repeated representations of the same record.
Different vouchers are retained even when their marker sequences are identical.
Outgroup accession records are not assumed to represent distinct individuals.

## Retrieve and verify

From the repository root:

```bash
python analyses/scripts/download_references.py
python analyses/scripts/audit_availability.py
```

The downloader reuses the bundled snapshot only after checking its hash and
accession membership. If the snapshot is missing it requests the versioned
accessions; an unexpected change from the recorded response hash requires review.
The audit verifies every extracted sequence and voucher against the records.
It does not infer specimen identity from sequence similarity alone.
