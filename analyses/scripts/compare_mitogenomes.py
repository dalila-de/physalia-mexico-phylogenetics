"""Prepare full-genome comparisons and report the best single BLAST HSP per source."""

import argparse
import csv
import gzip
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from common import ANALYSIS, QC, QC_IDS, RESULTS, UPSTREAM, records, ungap, write_table


def prepare() -> None:
    queries, references = [], []
    for sid in QC_IDS:
        recs = [r for p in (QC / sid / "mitogenomes").glob("*.fasta") for r in records(p)]
        if len(recs) != 1:
            raise ValueError(f"Expected one October path for {sid}, found {len(recs)}")
        queries.append(SeqRecord(recs[0].seq, id=sid, description=""))
    with gzip.open(ANALYSIS / "provenance/ncbi_records.gb.gz", "rt") as stream:
        for rec in SeqIO.parse(stream, "genbank"):
            if rec.seq.defined and len(rec) > 10000 and "mitochondrion" in rec.description.lower():
                references.append(SeqRecord(rec.seq, id="NCBI__" + rec.id, description=""))
    for rec in records(UPSTREAM / "mitogenome_tree/identification.aln.fa"):
        references.append(SeqRecord(Seq(ungap(rec.seq)), id="repo__" + rec.id, description=""))
    SeqIO.write(queries, ANALYSIS / "inputs/mitogenome_audit_queries.fasta", "fasta")
    SeqIO.write(references, ANALYSIS / "inputs/mitogenome_audit_references.fasta", "fasta")


def summarize() -> None:
    groups = defaultdict(list)
    with (RESULTS / "mitogenome_blast.tsv").open() as stream:
        for row in csv.reader(stream, delimiter="\t"):
            groups[(row[0], row[1].split("__")[0])].append(row)
    output = []
    for (sid, source), rows in sorted(groups.items()):
        score = max(float(r[14]) for r in rows)
        tied = [r for r in rows if float(r[14]) == score]
        row = sorted(tied, key=lambda r: r[1])[0]
        output.append(
            {
                "sample": sid,
                "reference_source": source,
                "best_reference": row[1],
                "equally_scoring_references": ";".join(r[1] for r in tied),
                "query_length": row[2],
                "reference_length": row[3],
                "aligned_columns": row[4],
                "percent_identity": row[5],
                "mismatches": row[6],
                "gap_bases": row[7],
                "query_start": row[8],
                "query_end": row[9],
                "reference_start": row[10],
                "reference_end": row[11],
                "query_coverage_percent_rounded": row[12],
                "bitscore": row[14],
            }
        )
    write_table(RESULTS / "mitogenome_nearest_matches.tsv", output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["prepare", "summarize"])
    if parser.parse_args().command == "prepare":
        prepare()
    else:
        summarize()
