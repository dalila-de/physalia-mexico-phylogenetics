#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
python analyses/scripts/compare_mitogenomes.py prepare
blastn \
    -query analyses/inputs/mitogenome_audit_queries.fasta \
    -subject analyses/inputs/mitogenome_audit_references.fasta \
    -task megablast -dust no -max_target_seqs 1000 -max_hsps 1 \
    -outfmt '6 qseqid sseqid qlen slen length pident mismatch gaps qstart qend sstart send qcovhsp evalue bitscore' \
    -out analyses/results/mitogenome_blast.tsv
python analyses/scripts/compare_mitogenomes.py summarize
