# Current NCBI-based analysis results

All five alignments and phylogenies are regenerated from the versioned NCBI
reference records in [references/ncbi](../../references/ncbi/README.md), retained
Mexican marker/mitogenome sequences, and CH887 COI/16S Sanger sequences.
No pre-existing alignment backbone is used.

[ncbi_reference_audit.tsv](ncbi_reference_audit.tsv) records the exact accession,
voucher, feature coordinates, length and hash of every reference input.
[ncbi_reference_summary.json](ncbi_reference_summary.json) summarizes the audit.
Use `python analyses/scripts/audit_availability.py` to repeat these checks.
The audit covers the explicitly selected records, not all NCBI sequences or aliases.

The earlier availability/novelty tables and previous tree results are historical,
not inputs to this rerun. Their identities and archive locations are recorded in
[rerun provenance](../../provenance/ncbi-rerun.json) and
[retired reference material](../../provenance/retired-reference-material.json).
The new specimens' mitochondrial assemblies retain automated MITOS2 annotations;
this analysis does not establish manually validated genome completeness.
