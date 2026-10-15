"""Validate retained hashes, cohort identity and agreement among tree artifacts."""

import csv
import hashlib
import json
from pathlib import Path

from Bio import Phylo, SeqIO

ROOT = Path(__file__).resolve().parents[1]


def main():
    with (ROOT / "provenance/artifact_checksums.tsv").open() as stream:
        hashes = list(csv.DictReader(stream, delimiter="\t"))
    for row in hashes:
        path = ROOT / row["file"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
            raise ValueError(f"Changed artifact: {path}")
    consolidation = json.loads((ROOT / "provenance/tree-consolidation.json").read_text())
    for row in consolidation["relocated"]:
        if "/genetic_groups/" not in row["to"]:
            if hashlib.sha256((ROOT / row["to"]).read_bytes()).hexdigest() != row["sha256"]:
                raise ValueError(f"Consolidation changed a scientific artifact: {row['to']}")
    if list((ROOT / "sanger-2026-10-14/results").rglob("*.treefile")):
        raise ValueError("Duplicate Sanger tree results remain outside analyses")
    mexico = {f"YPM-IZ-{i}" for i in range(116693, 116698)}
    for name in ("samples", "reads"):
        with (ROOT / f"illumina-qc-run/metadata/{name}.tsv").open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        if len(rows) != 5 or {r["sample"] for r in rows} != mexico:
            raise ValueError(f"Unexpected QC scope: {name}")
    qc = ROOT / "illumina-qc-run/results/cohort-2026-10-09"
    if {p.name for p in qc.iterdir() if p.is_dir()} != mexico:
        raise ValueError("Unexpected bundled QC specimens")
    expected = {"16S": 60, "18S": 48, "CO1": 57, "ITS": 46, "mitogenome_submitted": 48}
    for gene, count in expected.items():
        panel = (
            ROOT / "analyses/inputs/combined"
            if gene in {"16S", "CO1"}
            else ROOT / "analyses/inputs"
        )
        with (panel / f"{gene}.tips.tsv").open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        tips = {r["tip"] for r in rows}
        if gene in {"16S", "CO1"} and "sanger__CH887" not in tips:
            raise ValueError(f"Missing CH887 from combined panel: {gene}")
        if len(rows) != count or len(tips) != count:
            raise ValueError(f"Unexpected number of tips: {gene}")
        focal = [
            r["sample"]
            for r in rows
            if r["source"].startswith("October") or r["species"] == "Physalia physalis"
        ]
        if len(focal) != len(set(focal)):
            raise ValueError(f"Duplicate focal specimen: {gene}")
        for ext in ("treefile", "contree"):
            tree = Phylo.read(ROOT / f"analyses/results/trees/{gene}.{ext}", "newick")
            if {t.name for t in tree.get_terminals()} != tips:
                raise ValueError(f"Tree/input mismatch: {gene}")
        alignment = list(SeqIO.parse(ROOT / f"analyses/results/alignments/{gene}.fasta", "fasta"))
        if {r.id for r in alignment} != tips or len({len(r) for r in alignment}) != 1:
            raise ValueError(f"Alignment/input mismatch: {gene}")
    with (ROOT / "analyses/inputs/genetic_annotations.tsv").open() as stream:
        annotations = {r["sample"]: r for r in csv.DictReader(stream, delimiter="\t")}
    if any(annotations[s]["location"] != "Gulf of Mexico (Progreso)" for s in mexico):
        raise ValueError("Mexican locality labels differ")
    if annotations["CH887"]["location"] != "Caribbean Sea (Puerto Morelos)":
        raise ValueError("CH887 locality differs")
    for i in range(116693, 116698):
        sid = f"YPM-IZ-{i}"
        records = list(SeqIO.parse(ROOT / f"mitogenome-annotation/inputs/{sid}.fasta", "fasta"))
        if len(records) != 1 or len(records[0]) != 15038:
            raise ValueError(f"Unexpected Mexican assembly: {sid}")
        for ext in ("fas", "gff", "geneorder"):
            if not (ROOT / f"mitogenome-annotation/results/{sid}/result.{ext}").is_file():
                raise ValueError(f"Missing annotation: {sid}")
    print(
        f"Validated {len(hashes)} retained file hashes, five tree/alignment tip sets, and five Mexican mitogenome annotations."
    )


if __name__ == "__main__":
    main()
