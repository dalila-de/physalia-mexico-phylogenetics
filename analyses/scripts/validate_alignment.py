"""Normalize MAFFT direction labels, check scope, and remove all-gap columns."""

import argparse
import json
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq


def validate(source: Path, expected: Path, destination: Path) -> None:
    alignment = AlignIO.read(source, "fasta")
    expected_ids = {r.id for r in SeqIO.parse(expected, "fasta")}
    for r in alignment:
        r.id = r.id.removeprefix("_R_")
        r.description = ""
        if set(str(r.seq).upper()) - set("ACGTRYSWKMBDHVN-?"):
            raise ValueError(f"Invalid alphabet for {r.id}")
    names = [r.id for r in alignment]
    if len(names) != len(set(names)) or set(names) != expected_ids:
        raise ValueError("Alignment and input tip sets differ")
    length = alignment.get_alignment_length()
    keep = [i for i in range(length) if any(r.seq[i] not in "-?" for r in alignment)]
    if not keep:
        raise ValueError("Empty alignment")
    for r in alignment:
        r.seq = Seq("".join(r.seq[i] for i in keep))
    destination.parent.mkdir(parents=True, exist_ok=True)
    AlignIO.write(alignment, destination, "fasta")
    destination.with_suffix(".json").write_text(
        json.dumps(
            {"tips": len(names), "sites": len(keep), "all_gap_columns_removed": length - len(keep)},
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("expected", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    validate(args.source, args.expected, args.destination)
