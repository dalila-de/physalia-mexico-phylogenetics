# Public sequence availability for the 40 October specimens

Checked 2026-10-12 against upstream commit
`d08f0edf4d819e21c85705e4965de71673808dec` and a live NCBI nucleotide search.

| Cohort | Specimens | 16S and CO1 in repository | 18S and ITS in repository / NCBI | Mitogenomes in repository / NCBI |
|---|---:|---:|---:|---:|
| Previously published | 35 | 35 | 35 / 35 | 35 / 34 |
| Mexico, 116693–116697 | 5 | 0 | 0 / 0 | 0 / 0 |

For 34 published specimens, 16S and CO1 are available **inside deposited
mitogenomes**, rather than necessarily as separate marker records. The complete
accession mapping is in [specimen_marker_availability.tsv](specimen_marker_availability.tsv).
The 35 nuclear 18S records and 35 ITS records were verified from downloaded
GenBank specimen/isolate metadata, not inferred from accession ranges.

**YPM-IZ-110972 is the exception.** Its 16S, CO1 and mitochondrial assembly are
available in the repository, but no corresponding mitochondrial GenBank record
was found. Its 18S and ITS are public in GenBank. The authors explicitly withheld
its mitochondrial submission because the long assembly/annotation differed from
the other samples. The repository retains two mitochondrial path variants for
this specimen. This specimen is excluded from the mitochondrial tree.

## Mexico specimens

No records matched these five specimen IDs in either the pinned reference panel
or NCBI. A second NCBI search across **all taxa**, using hyphenated, colon-separated
and space-separated voucher forms, returned zero records. Thus these appear to
be newly represented specimens relative to the resources searched. Records under
unknown aliases cannot be excluded, and absence of an assembled marker record
does not establish absence of raw reads in SRA.

None of the five whole mitochondrial sequences has an exact whole-sequence match
in the repository/NCBI mitochondrial collection, including reverse complements
and equivalent rotations. Their highest-scoring single BLAST alignments to
either resource are:

| Mexico specimen | Identity | Mismatches | Gap bases | Aligned columns |
|---|---:|---:|---:|---:|
| YPM-IZ-116693 | 99.874% | 19 | 0 | 15,038 |
| YPM-IZ-116694 | 99.661% | 47 | 4 | 15,033 |
| YPM-IZ-116695 | 99.761% | 36 | 0 | 15,035 |
| YPM-IZ-116696 | 99.721% | 42 | 0 | 15,038 |
| YPM-IZ-116697 | 99.880% | 18 | 0 | 15,038 |

All five query assemblies are 15,038 bp. These are local alignment comparisons,
not global edit distances, and do not establish sequencing accuracy or novel
species. They support treating the assemblies as new candidate sequences pending
phylogenetic and assembly review. Individual short marker sequences may already
occur in other specimens; that does not make the Mexico specimen records duplicates.
See [nearest mitochondrial matches](mitogenome_nearest_matches.tsv) and
[whole-sequence matches](mitogenome_haplotype_matches.tsv).

## Scope and reproducibility

The genus-wide query was `Physalia[Organism] AND 100:100000[SLEN]`, yielding 729
records. Four WGS/TSA master records lacked sequence content; their metadata is
inventoried, but they were excluded from sequence identity checks. The timestamped
query responses, download URLs and hashes are preserved in `../provenance/`.
The original GenBank responses are bundled as `../provenance/ncbi_records.gb.gz`.
Availability covers all 40 focal specimens; local sequence comparisons use only
the five Mexican assemblies. The original published-specimen reassembly
comparisons are preserved in the external historical snapshot.

[all_sharkmer_product_matches.tsv](all_sharkmer_product_matches.tsv) covers the five Mexican specimens
only (131 products), including alternative primer
sets, 28S, EF1A and non-host primer products. A missing exact match in that table
is **not** a comprehensive novelty determination for those additional markers.
Bacterial-primer products are not automatically Physalia genes. Sequence identity
and specimen identity are kept separate throughout.

Sources: [upstream phylogenetics](https://github.com/shchurch/Physalia_spatiotemporal_patterns/tree/d08f0edf4d819e21c85705e4965de71673808dec/phylogenetics),
[mitochondrial submission scope](https://github.com/shchurch/Physalia_spatiotemporal_patterns/blob/d08f0edf4d819e21c85705e4965de71673808dec/phylogenetics/submission/build_gene_submission.py),
[example mitochondrial record PZ224317](https://www.ncbi.nlm.nih.gov/nuccore/PZ224317.1).
