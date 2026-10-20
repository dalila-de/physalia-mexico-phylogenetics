"""Extract reference sequences from GenBank features with explicit identity checks."""

import gzip
import re

from Bio import SeqIO


def normalize_voucher(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", value.upper()).removesuffix("ANNOTATION")


def extract_sequence(record, mode: str):
    """Select one annotated feature; never concatenate separate gene fragments."""
    if mode == "whole_record":
        return record.seq, "whole_record"
    options = []
    for feature in record.features:
        genes = {s.lower() for s in feature.qualifiers.get("gene", [])}
        products = " ".join(feature.qualifiers.get("product", [])).lower()
        if "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers:
            continue
        if mode == "cox1_CDS":
            selected = feature.type == "CDS" and bool(genes & {"cox1", "coi", "co1"})
        elif mode == "16S_rRNA":
            selected = feature.type == "rRNA" and (
                bool(genes & {"rrnl", "rnl", "16s"})
                or "16s" in products
                or "large subunit ribosomal" in products
            )
        else:
            raise ValueError(f"Unknown extraction mode: {mode}")
        if selected:
            options.append(feature)
    if not options:
        raise ValueError(f"Missing annotated {mode} in {record.id}")
    feature = max(options, key=len)
    return feature.extract(record.seq), f"{feature.type}:{feature.location}"


def load_records(path):
    with gzip.open(path, "rt") as stream:
        records = list(SeqIO.parse(stream, "genbank"))
    lookup = {r.id: r for r in records}
    if len(lookup) != len(records):
        raise ValueError("Duplicate accession in GenBank snapshot")
    return lookup


def verify_identity(record, row):
    if record.id != row["accession"] or record.annotations["organism"] != row["species"]:
        raise ValueError(f"Accession/species mismatch: {row}")
    if row["species"].startswith("Physalia"):
        identifiers = {
            normalize_voucher(value)
            for feature in record.features
            if feature.type == "source"
            for key in ("isolate", "specimen_voucher")
            for value in feature.qualifiers.get(key, [])
        }
        if normalize_voucher(row["sample"]) not in identifiers:
            raise ValueError(f"Voucher/isolate mismatch: {row['sample']} / {record.id}")
