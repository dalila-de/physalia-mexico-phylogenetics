"""Retain validated compact outputs and explicit recovery status on PI storage."""

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
from pathlib import Path

from Bio import SeqIO
from qc_report import write_report


def fasta_records(path: Path) -> list:
    records = list(SeqIO.parse(path, "fasta"))
    if not records or any(
        not r.seq or set(str(r.seq).upper()) - set("ACGTRYSWKMBDHVN-") for r in records
    ):
        raise ValueError(f"Empty or invalid nucleotide FASTA: {path}")
    return records


def retain_sample(work: Path, destination: Path, sample: str, maximum: int) -> None:
    work, destination = work.resolve(), destination.resolve()
    final = destination / sample
    staging = destination / f".{sample}.tmp"
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True, exist_ok=True)
    # Reports use the established output layout, now rooted entirely in scratch.
    os.chdir(work)
    write_report(sample, staging / "qc_report.txt", maximum)
    trim = (work / "logs/trimmomatic" / f"{sample}.log").read_text()
    counts = re.findall(
        r"Input Read Pairs: (\d+) Both Surviving: (\d+).*?Forward Only Surviving: (\d+)"
        r".*?Reverse Only Surviving: (\d+).*?Dropped: (\d+)",
        trim,
    )
    if len(counts) != 1:
        raise ValueError(f"Missing or ambiguous trimming counts: {sample}")
    input_pairs, surviving, forward, reverse, dropped = map(int, counts[0])
    selected = json.loads((work / "reads" / f"{sample}.json").read_text())
    if input_pairs != sum(r["selected_pairs"] for r in selected):
        raise ValueError(f"Subsampling/trimming count mismatch: {sample}")
    if surviving + forward + reverse + dropped != input_pairs:
        raise ValueError(f"Trimming categories do not sum to input: {sample}")
    mapped = int((work / "human" / f"{sample}.mapped_paired").read_text())
    markers = sorted((work / "shark_results" / sample).glob("*.fasta"))
    assemblies = sorted(
        (work / "getorganelle" / f"{sample}_mtgenome").glob("*.path_sequence.fasta")
    )
    lengths = {p.name: [len(r.seq) for r in fasta_records(p)] for p in markers + assemblies}
    getorg_log = work / "getorganelle" / f"{sample}_mtgenome" / "get_org.log.txt"
    errors = [line for line in getorg_log.read_text().splitlines() if "ERROR" in line]
    summary = {
        "sample": sample,
        "workflow_status": "complete",
        "input_pairs": input_pairs,
        "retained_pairs": surviving,
        "forward_only_reads": forward,
        "reverse_only_reads": reverse,
        "dropped_pairs": dropped,
        "human_mapped_records": mapped,
        "human_percent_reads": 100 * mapped / (2 * surviving),
        "marker_files": len(markers),
        "mitogenome_status": "paths_recovered" if assemblies else "no_path_recovered",
        "mitogenome_files": len(assemblies),
        "getorganelle_error_lines": len(errors),
        "source_work_directory": str(work),
    }
    if errors:
        summary["mitogenome_status"] += "_with_logged_errors"
    (staging / "sequence_lengths.json").write_text(json.dumps(lengths, indent=2) + "\n")
    for label, paths in [("markers", markers), ("mitogenomes", assemblies)]:
        folder = staging / label
        folder.mkdir(exist_ok=True)
        for path in paths:
            shutil.copy2(path, folder / path.name)
    # Keep reports and small evidence outside scratch's purge policy.
    evidence = staging / "evidence"
    evidence.mkdir(exist_ok=True)
    artifacts = [
        work / "reads" / f"{sample}.json",
        work / "logs/trimmomatic" / f"{sample}.log",
        work / "kraken" / f"{sample}.report.txt",
        work / "human" / f"{sample}.mapped_paired",
        work / "shark_results" / sample / f"{sample}.stats",
        getorg_log,
    ]
    artifacts.extend((work / "fastqc" / sample).glob("*"))
    artifacts.extend((work / "rego_results" / sample).glob("*.txt"))
    for path in artifacts:
        if path.is_file():
            shutil.copy2(path, evidence / path.name)
    inventory = []
    for path in sorted(staging.rglob("*")):
        if path.is_file() and path.name not in {"checksums.tsv", "summary.json"}:
            inventory.append(
                f"{hashlib.sha256(path.read_bytes()).hexdigest()}\t{path.relative_to(staging)}"
            )
    (staging / "checksums.tsv").write_text("sha256\tfile\n" + "\n".join(inventory) + "\n")
    (staging / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    final.mkdir(parents=True, exist_ok=True)
    # A rerun may recover fewer markers; remove obsolete files from these
    # workflow-owned result folders rather than mixing separate attempts.
    for label in ("markers", "mitogenomes", "evidence"):
        if (final / label).exists():
            shutil.rmtree(final / label)
    for path in staging.iterdir():
        if path.is_dir():
            shutil.copytree(path, final / path.name, dirs_exist_ok=True)
        elif path.name != "summary.json":
            shutil.copy2(path, final / path.name)
    # Completion is written last, after all retained files have been copied.
    (staging / "summary.json").replace(final / "summary.json")
    shutil.rmtree(staging)


def cohort(manifest: Path, reports: list[Path], output: Path) -> None:
    with manifest.open() as stream:
        samples = list(csv.DictReader(stream, delimiter="\t"))
    rows = {json.loads(p.read_text())["sample"]: json.loads(p.read_text()) for p in reports}
    if len(rows) != len(reports) or set(rows) != {r["sample"] for r in samples}:
        raise ValueError("Summary must contain exactly the approved specimens")
    combined = [{**sample, **rows[sample["sample"]]} for sample in samples]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(combined[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(combined)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    sample = subparsers.add_parser("sample")
    sample.add_argument("--work", type=Path, required=True)
    sample.add_argument("--destination", type=Path, required=True)
    sample.add_argument("--sample", required=True)
    sample.add_argument("--maximum", type=int, required=True)
    aggregate = subparsers.add_parser("cohort")
    aggregate.add_argument("--manifest", type=Path, required=True)
    aggregate.add_argument("--reports", type=Path, nargs="+", required=True)
    aggregate.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "sample":
        retain_sample(args.work, args.destination, args.sample, args.maximum)
    else:
        cohort(args.manifest, args.reports, args.output)


if __name__ == "__main__":
    main()
