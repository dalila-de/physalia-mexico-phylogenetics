"""Audit public records by specimen identity separately from sequence identity."""

import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict

from Bio import SeqIO
from common import (
    ANALYSIS,
    GENES,
    IDS,
    QC,
    RESULTS,
    SAMPLES,
    UPSTREAM,
    records,
    sequence_match,
    specimen,
    table,
    ungap,
    write_table,
)


def normalize(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", value.upper())


def classify(record) -> str:
    title = record.description.lower()
    features = " ".join(
        " ".join(f.qualifiers.get("product", []) + f.qualifiers.get("gene", []))
        for f in record.features
    ).lower()
    if "mitochondrion" in title and len(record) > 10000:
        return "mitogenome"
    if "internal transcribed spacer" in title or "internal transcribed spacer" in features:
        return "ITS"
    if "small subunit ribosomal RNA".lower() in features and not any(
        f.type == "source" and "organelle" in f.qualifiers for f in record.features
    ):
        return "18S"
    for gene, pattern in [
        ("18S", "18s"),
        ("28S", "28s"),
        ("16S", "16s"),
        ("CO1", "cytochrome c oxidase subunit i"),
        ("CO1", "cytochrome oxidase subunit i"),
        ("EF1A", "elongation factor"),
    ]:
        if pattern in title or pattern in features:
            return gene
    return "other"


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    with gzip.open(ANALYSIS / "provenance/ncbi_records.gb.gz", "rt") as stream:
        downloaded = list(SeqIO.parse(stream, "genbank"))
    download = json.loads((ANALYSIS / "provenance/ncbi_download.json").read_text())
    if len(downloaded) != download["count"]:
        raise ValueError("NCBI retrieval incomplete")
    ncbi = []
    mapping = defaultdict(list)
    for record in downloaded:
        sources = [f for f in record.features if f.type == "source"]
        qualifiers = {k: ";".join(v) for f in sources for k, v in f.qualifiers.items()}
        identity_fields = [
            qualifiers.get(k, "") for k in ("isolate", "specimen_voucher", "strain", "clone")
        ]
        # Exact canonicalized voucher/isolate IDs, allowing submission suffixes.
        matches = [
            s
            for s in IDS
            if any(normalize(v).removesuffix("ANNOTATION") == normalize(s) for v in identity_fields)
        ]
        if len(matches) > 1:
            raise ValueError(f"Ambiguous specimen mapping: {record.id}")
        row = {
            "accession": record.id,
            "sample": matches[0] if matches else "",
            "marker": classify(record),
            "length": len(record),
            "organism": record.annotations.get("organism", ""),
            "isolate": qualifiers.get("isolate", ""),
            "specimen_voucher": qualifiers.get("specimen_voucher", ""),
            "description": record.description,
            "sequence_sha256": hashlib.sha256(ungap(record.seq).encode()).hexdigest()
            if record.seq.defined
            else "sequence_not_resolved",
        }
        ncbi.append(row)
        if matches and record.seq.defined:
            mapping[matches[0]].append((row, ungap(record.seq)))
    write_table(RESULTS / "ncbi_record_inventory.tsv", ncbi)
    repo = {
        g: {r.id: ungap(r.seq) for r in records(UPSTREAM / "gene_trees" / f"{g}.all.fasta")}
        for g in GENES
    }
    mitos = defaultdict(list)
    for r in records(UPSTREAM / "mitogenome_tree/identification.aln.fa"):
        mitos[specimen(r.id)].append((r.id, ungap(r.seq)))
    accession_table = {
        r["Isolate"]: r["GenBank_accession"]
        for r in table(
            UPSTREAM / "mitogenome_submission/Physalia_mitogenomes_sourcemodifiers_accessioned.tsv"
        )
    }
    specimen_rows, comparisons, all_markers = [], [], []
    for metadata in SAMPLES:
        sample = metadata["sample"]
        for gene in [*GENES, "mitogenome"]:
            references = (
                mitos[sample]
                if gene == "mitogenome"
                else [(sample, repo[gene][sample])]
                if sample in repo[gene]
                else []
            )
            public = [
                (row, seq)
                for row, seq in mapping[sample]
                if row["marker"] == gene
                or (gene in {"16S", "CO1"} and row["marker"] == "mitogenome")
            ]
            paths = (
                sorted((QC / sample / "mitogenomes").glob("*.fasta"))
                if gene == "mitogenome"
                else [QC / sample / "markers" / f"{sample}_{GENES[gene]}.fasta"]
            )
            local = [(p, r) for p in paths if p.exists() for r in records(p)]
            if gene != "mitogenome" and local:
                local = [max(local, key=lambda pair: len(pair[1]))]
            specimen_rows.append(
                {
                    "sample": sample,
                    "cohort": metadata["cohort"],
                    "marker": gene,
                    "local_records": len(local),
                    "repo_records": len(references),
                    "ncbi_accessions": ";".join(row["accession"] for row, _ in public),
                    "ncbi_record_type": ";".join(sorted({row["marker"] for row, _ in public})),
                    "availability": "NCBI_and_repo"
                    if public and references
                    else "repo_only"
                    if references
                    else "NCBI_only"
                    if public
                    else "not_found_in_search",
                    "repository_listed_mito_accession": accession_table.get(sample, "")
                    if gene == "mitogenome"
                    else "",
                }
            )
            for path, query in local:
                for origin, refs in [
                    ("repository", references),
                    ("NCBI", [(row["accession"], seq) for row, seq in public]),
                ]:
                    for ref_id, sequence in refs:
                        comparisons.append(
                            {
                                "sample": sample,
                                "marker": gene,
                                "query_file": str(path.relative_to(ANALYSIS.parent)),
                                "query_header": query.description,
                                "query_length": len(query),
                                "reference_source": origin,
                                "reference_id": ref_id,
                                "reference_length": len(sequence),
                                "relationship": sequence_match(
                                    ungap(query.seq), sequence, gene == "mitogenome"
                                ),
                            }
                        )
        # Inventory every sharkmer product, including alternative and non-host primers.
        for path in sorted((QC / sample / "markers").glob("*.fasta")):
            primer = path.stem.removeprefix(sample + "_")
            for index, query in enumerate(records(path)):
                exact = []
                for row, seq in mapping[sample]:
                    relation = sequence_match(ungap(query.seq), seq)
                    if relation != "not_exact_or_contained":
                        exact.append(f"{row['accession']}:{relation}")
                other_matches = []
                for record, row in zip(downloaded, ncbi, strict=True):
                    if row["sample"] != sample and record.seq.defined:
                        relation = sequence_match(ungap(query.seq), ungap(record.seq))
                        if relation in {
                            "identical",
                            "identical_reverse_complement",
                            "query_contained_in_reference",
                        }:
                            other_matches.append(record.id)
                all_markers.append(
                    {
                        "sample": sample,
                        "primer": primer,
                        "product_index": index,
                        "length": len(query),
                        "same_specimen_exact_or_contained": ";".join(exact),
                        "other_specimen_or_unassigned_exact_matches": ";".join(other_matches),
                        "interpretation": "non_host_primer_do_not_infer_specimen_novelty"
                        if primer.startswith("bacteria_")
                        else "sequence_match_is_not_specimen_identity",
                    }
                )
    write_table(RESULTS / "specimen_marker_availability.tsv", specimen_rows)
    # Newly sampled vouchers can have no same-specimen public comparisons.
    with (RESULTS / "same_specimen_sequence_comparisons.tsv").open("w") as stream:
        writer = csv.DictWriter(
            stream,
            delimiter="\t",
            fieldnames=[
                "sample",
                "marker",
                "query_file",
                "query_header",
                "query_length",
                "reference_source",
                "reference_id",
                "reference_length",
                "relationship",
            ],
        )
        writer.writeheader()
        writer.writerows(comparisons)
    write_table(RESULTS / "all_sharkmer_product_matches.tsv", all_markers)
    haplotypes = []
    for sid in IDS:
        for path in sorted((QC / sid / "mitogenomes").glob("*.fasta")):
            for query in records(path):
                ncbi_matches, repo_matches = [], []
                for record, row in zip(downloaded, ncbi, strict=True):
                    if row["marker"] == "mitogenome" and record.seq.defined:
                        relation = sequence_match(ungap(query.seq), ungap(record.seq), True)
                        if relation in {
                            "identical",
                            "identical_reverse_complement",
                            "rotation_equivalent",
                        }:
                            ncbi_matches.append(record.id)
                for ref_sample, refs in mitos.items():
                    for ref_id, seq in refs:
                        if sequence_match(ungap(query.seq), seq, True) in {
                            "identical",
                            "identical_reverse_complement",
                            "rotation_equivalent",
                        }:
                            repo_matches.append(ref_id)
                haplotypes.append(
                    {
                        "sample": sid,
                        "length": len(query),
                        "ncbi_whole_sequence_matches": ";".join(ncbi_matches),
                        "repository_whole_sequence_matches": ";".join(repo_matches),
                        "interpretation": "Sequence identity does not imply the specimen was previously deposited.",
                    }
                )
    write_table(RESULTS / "mitogenome_haplotype_matches.tsv", haplotypes)
    summary = {
        "ncbi_records_retrieved": len(ncbi),
        "specimen_marker_cells": len(specimen_rows),
        "availability": dict(Counter(r["availability"] for r in specimen_rows)),
        "same_specimen_relationships": dict(Counter(r["relationship"] for r in comparisons)),
        "local_sharkmer_records": len(all_markers),
        "query": download["query"],
        "query_utc": download["utc"],
    }
    (RESULTS / "availability_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
