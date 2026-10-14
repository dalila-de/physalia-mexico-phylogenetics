"""Record identities of consumed QC, upstream, NCBI, and annotation resources."""

import argparse
import hashlib
from pathlib import Path

from common import CONFIG, IDS, QC, SCRATCH, UPSTREAM, write_table


def main(output: Path) -> None:
    paths = [p for sid in IDS for p in (QC / sid).rglob("*.fasta")]
    paths += list((UPSTREAM / "gene_trees").glob("*.all.fasta"))
    paths += [UPSTREAM / "tip_metadata.tsv"]
    paths += [
        UPSTREAM / "mitogenome_tree" / name
        for name in [
            "identification.aln.fa",
            "submitted_rooted.aln.fa",
            "rhizophysa_outgroups.fasta",
        ]
    ]
    paths += list((SCRATCH / "ncbi").glob("*.gb"))
    paths += [
        p
        for p in (Path(CONFIG["mitos_refdir"]) / CONFIG["mitos_refname"]).rglob("*")
        if p.is_file()
    ]
    paths += [Path(CONFIG["iqtree"]), Path(CONFIG["mitos_env"]) / "bin/runmitos.py"]
    rows = []
    for p in paths:
        digest = hashlib.sha256()
        with p.open("rb") as stream:
            for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                digest.update(block)
        rows.append({"path": str(p), "bytes": p.stat().st_size, "sha256": digest.hexdigest()})
    write_table(output, rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args().output)
