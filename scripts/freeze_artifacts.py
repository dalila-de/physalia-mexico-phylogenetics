"""Record checksums of scientific artifacts; do not run to conceal unexpected changes."""

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
paths = [
    "illumina-qc-run/results",
    "mitogenome-annotation/inputs",
    "mitogenome-annotation/results",
    "analyses/inputs",
    "analyses/results/alignments",
    "analyses/results/trees",
    "references/church",
]
files = sorted({p for name in paths for p in (ROOT / name).rglob("*") if p.is_file()})
with (ROOT / "provenance/artifact_checksums.tsv").open("w") as stream:
    writer = csv.DictWriter(stream, fieldnames=["file", "sha256", "bytes"], delimiter="\t")
    writer.writeheader()
    for path in files:
        writer.writerow(
            {
                "file": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "bytes": path.stat().st_size,
            }
        )
