"""Verify every retained reference input against its versioned NCBI record."""

import hashlib
import json
from collections import Counter

from Bio import SeqIO
from common import ANALYSIS, ROOT, table, write_table
from ncbi_references import extract_sequence, load_records, verify_identity


def main():
    reference = ROOT / "references/ncbi"
    records = load_records(reference / "records.gb.gz")
    manifest = table(reference / "accessions.tsv")
    expected = {(r["marker"], r["sample"]): r for r in manifest}
    specimen_rows = {}
    for row in manifest:
        if row["species"].startswith("Physalia"):
            sid = row["sample"]
            specimen_rows.setdefault(
                sid,
                dict(
                    specimen_id=sid,
                    species=row["species"],
                    COI="",
                    **{"16S": "", "18S": "", "ITS": "", "mitogenome": ""},
                ),
            )
            marker = {"CO1": "COI", "mitogenome_submitted": "mitogenome"}.get(
                row["marker"], row["marker"]
            )
            specimen_rows[sid][marker] = row["accession"]
    write_table(
        reference / "specimen_accessions.tsv", [specimen_rows[s] for s in sorted(specimen_rows)]
    )
    write_table(
        reference / "rhizophysa_accessions.tsv",
        [r for r in manifest if r["species"].startswith("Rhizophysa")],
    )
    results = []
    for gene in ("CO1", "16S", "18S", "ITS", "mitogenome_submitted"):
        panel = ANALYSIS / "inputs" / ("combined" if gene in {"CO1", "16S"} else "")
        sequences = {r.id: str(r.seq) for r in SeqIO.parse(panel / f"{gene}.fasta", "fasta")}
        for tip in table(panel / f"{gene}.tips.tsv"):
            q = sequences[tip["tip"]]
            if hashlib.sha256(q.encode()).hexdigest() != tip["sequence_sha256"]:
                raise ValueError(f"Input hash differs: {tip['tip']}")
            if tip["source"] != "NCBI_GenBank":
                continue
            entry = expected[(gene, tip["sample"])]
            record = records[entry["accession"]]
            verify_identity(record, entry)
            sequence, feature = extract_sequence(record, entry["extraction"])
            if str(sequence).upper() != q:
                raise ValueError(f"Input does not equal NCBI extraction: {tip['tip']}")
            results.append(
                dict(
                    marker=gene,
                    sample=tip["sample"],
                    species=entry["species"],
                    accession=record.id,
                    tip=tip["tip"],
                    extraction=entry["extraction"],
                    feature=feature,
                    length=len(q),
                    sequence_sha256=tip["sequence_sha256"],
                    verification="exact_feature_or_record_match",
                )
            )
    if len(results) != len(manifest):
        raise ValueError("Missing reference manifest rows")
    write_table(ANALYSIS / "results/ncbi_reference_audit.tsv", results)
    summary = dict(
        distinct_accessions=len(records),
        reference_input_occurrences=len(results),
        by_marker=dict(Counter(r["marker"] for r in results)),
        identity="versioned accession and voucher/isolate checked; exact feature/record sequence match",
    )
    (ANALYSIS / "results/ncbi_reference_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
