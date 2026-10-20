"""Record current accession-defined references and actual tree input identities."""

import argparse
import hashlib
from pathlib import Path

from common import ANALYSIS, ROOT, write_table


def main(output: Path) -> None:
    paths = list((ROOT / "references/ncbi").glob("*"))
    paths += list((ANALYSIS / "inputs").glob("*.fasta"))
    paths += list((ANALYSIS / "inputs").glob("*.tips.tsv"))
    paths += list((ANALYSIS / "inputs/combined").glob("*"))
    paths += [ANALYSIS / "config/config.yaml", ANALYSIS / "inputs/genetic_annotations.tsv"]
    paths += list((ROOT / "sanger-2026-10-14/inputs").glob("*.sanger.fasta"))
    rows = [
        dict(
            path=str(p.relative_to(ROOT)),
            bytes=p.stat().st_size,
            sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
        )
        for p in sorted(set(paths))
        if p.is_file()
    ]
    write_table(output, rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args().output)
