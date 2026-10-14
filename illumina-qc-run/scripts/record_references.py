"""Hash the actual references and auxiliary tools once inside the allocation."""

import argparse
import csv
import hashlib
from pathlib import Path

import yaml


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = yaml.safe_load(Path("config/config.yaml").read_text())
    paths = [Path(config[k]) for k in ("human_bwa", "mito_seed", "sharkmer", "rego")]
    paths += [Path(config["human_bwa"] + ext) for ext in (".amb", ".ann", ".bwt", ".pac", ".sa")]
    paths += [Path(config["kraken_db"]) / name for name in ("hash.k2d", "opts.k2d", "taxo.k2d")]
    paths += [p for p in Path(config["getorganelle_db"]).rglob("*") if p.is_file()]
    with args.output.open("w") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["path", "bytes", "mtime_ns", "sha256"])
        for path in paths:
            before = path.stat()
            digest = hashlib.sha256()
            with path.open("rb") as stream:
                for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                    digest.update(block)
            after = path.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise ValueError(f"Reference changed while hashing: {path}")
            writer.writerow([str(path), before.st_size, before.st_mtime_ns, digest.hexdigest()])


if __name__ == "__main__":
    main()
