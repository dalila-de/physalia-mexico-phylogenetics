"""Select paired FASTQ prefixes, validating every consumed record and mate ID."""

import argparse
import csv
import gzip
import hashlib
import json
from contextlib import ExitStack
from pathlib import Path


def record_id(header: bytes) -> bytes:
    """Normalize only the conventional /1 or /2 mate suffix."""
    name = header.split()[0]
    return name[:-2] if name.endswith((b"/1", b"/2")) else name


def read_record(stream, path: str) -> list[bytes] | None:
    """Read one four-line Illumina record; reject malformed consumed input."""
    header = stream.readline()
    if not header:
        return None
    lines = [header, stream.readline(), stream.readline(), stream.readline()]
    seq, qual = lines[1].rstrip(b"\r\n"), lines[3].rstrip(b"\r\n")
    if (
        not header.startswith(b"@")
        or not lines[2].startswith(b"+")
        or not seq
        or len(seq) != len(qual)
        or set(seq.upper()) - set(b"ACGTN")
        or any(x < 33 or x > 126 for x in qual)
    ):
        raise ValueError(f"Malformed FASTQ record in {path}: {header[:120]!r}")
    return lines


def select_pairs(rows: list[dict], maximum: int, r1: Path, r2: Path) -> list[dict]:
    """Take an equal capped prefix per sequencing unit without separating mates."""
    if not rows or maximum < len(rows):
        raise ValueError("At least one pair per sequencing unit is required")
    records = []
    r1.parent.mkdir(parents=True, exist_ok=True)
    r2.parent.mkdir(parents=True, exist_ok=True)
    with r1.open("wb") as out1, r2.open("wb") as out2:
        for index, row in enumerate(rows):
            limit = maximum // len(rows) + (index < maximum % len(rows))
            hashes = [hashlib.sha256(), hashlib.sha256()]
            count = 0
            with ExitStack() as stack:
                streams = [
                    stack.enter_context(
                        gzip.open(row[mate], "rb")
                        if row[mate].endswith(".gz")
                        else open(row[mate], "rb")
                    )
                    for mate in ("R1", "R2")
                ]
                for _ in range(limit):
                    left = read_record(streams[0], row["R1"])
                    right = read_record(streams[1], row["R2"])
                    if left is None and right is None:
                        break
                    if left is None or right is None:
                        raise ValueError(f"Unequal mate counts in {row['sample']} {row['unit']}")
                    if record_id(left[0]) != record_id(right[0]):
                        raise ValueError(f"Mismatched mates: {left[0]!r}, {right[0]!r}")
                    for record, out, digest in zip(
                        (left, right), (out1, out2), hashes, strict=True
                    ):
                        chunk = b"".join(record)
                        out.write(chunk)
                        digest.update(chunk)
                    count += 1
            if count == 0:
                raise ValueError(f"No read pairs in {row['sample']} {row['unit']}")
            records.append(
                {
                    **row,
                    "requested_pairs": limit,
                    "selected_pairs": count,
                    "R1_selected_sha256": hashes[0].hexdigest(),
                    "R2_selected_sha256": hashes[1].hexdigest(),
                    "input_file_stats": {
                        mate: {
                            "bytes": Path(row[mate]).stat().st_size,
                            "mtime_ns": Path(row[mate]).stat().st_mtime_ns,
                        }
                        for mate in ("R1", "R2")
                    },
                }
            )
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--maximum", type=int, required=True)
    parser.add_argument("--r1", type=Path, required=True)
    parser.add_argument("--r2", type=Path, required=True)
    parser.add_argument("--counts", type=Path, required=True)
    args = parser.parse_args()
    with args.manifest.open() as stream:
        rows = [r for r in csv.DictReader(stream, delimiter="\t") if r["sample"] == args.sample]
    records = select_pairs(rows, args.maximum, args.r1, args.r2)
    args.counts.write_text(json.dumps(records, indent=2) + "\n")


if __name__ == "__main__":
    main()
