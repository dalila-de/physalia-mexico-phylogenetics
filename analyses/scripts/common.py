"""Paths, sequence handling and tabular output for the pinned analysis."""

import csv
import re
from pathlib import Path

import yaml
from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "analyses"
CONFIG = yaml.safe_load((ANALYSIS / "config/config.yaml").read_text())
SCRATCH = ROOT / CONFIG["scratch"]
RESULTS = ANALYSIS / "results"
QC = ROOT / "illumina-preparation-qc-mitogenome-assembly/results" / CONFIG["qc_run"]
GENES = {"16S": "cnidaria_16S", "18S": "cnidaria_18S", "CO1": "cnidaria_CO1", "ITS": "cnidaria_ITS"}


def table(path: Path) -> list[dict]:
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


SAMPLES = table(ANALYSIS / "inputs/cohort.tsv")
QC_IDS = [row["sample"] for row in table(ROOT / "illumina-preparation-qc-mitogenome-assembly/metadata/samples.tsv")]
IDS = [row["sample"] for row in SAMPLES]


def write_table(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def ungap(sequence) -> str:
    return str(sequence).upper().replace("-", "").replace("?", "N")


def specimen(identifier: str) -> str:
    identifier = identifier.removeprefix("_R_")
    return re.split(r"_(?:mtgenome|mito|annotation)(?:_|$)", identifier)[0]


def records(path: Path) -> list:
    seqs = list(SeqIO.parse(path, "fasta"))
    if not seqs or any(
        not x.seq or set(str(x.seq).upper()) - set("ACGTRYSWKMBDHVN-?") for x in seqs
    ):
        raise ValueError(f"Empty or malformed nucleotide FASTA: {path}")
    return seqs


def sequence_match(query: str, reference: str, genome: bool = False) -> str:
    """Exact string evidence only; absence of a match does not imply novelty."""
    if query == reference:
        return "identical"
    reverse = str(Seq(query).reverse_complement())
    if reverse == reference:
        return "identical_reverse_complement"
    if query in reference or reverse in reference:
        return "query_contained_in_reference"
    if reference in query or reference in reverse:
        return "reference_contained_in_query"
    if (
        genome
        and len(query) == len(reference)
        and (query in reference + reference or reverse in reference + reference)
    ):
        return "rotation_equivalent"
    return "not_exact_or_contained"
