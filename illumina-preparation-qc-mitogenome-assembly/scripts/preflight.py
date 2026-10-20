"""Check manifest joins, accessible inputs, paired prefixes and references."""

import argparse
import csv
import gzip
import json
from pathlib import Path

import yaml
from subset_reads import read_record, record_id


def check(root: Path) -> dict:
    config = yaml.safe_load((root / "config/config.yaml").read_text())
    with (root / config["samples"]).open() as stream:
        samples = list(csv.DictReader(stream, delimiter="\t"))
    with (root / config["reads"]).open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    ids = [r["sample"] for r in samples]
    if (
        len(ids) != 5
        or set(ids) != {f"YPM-IZ-{i}" for i in range(116693, 116698)}
        or len(rows) != 5
    ):
        raise ValueError("Approved scope is five Mexican specimens and five sequencing inputs")
    if set(ids) != {r["sample"] for r in rows}:
        raise ValueError("Sample/read manifest mismatch")
    for sample in samples:
        if int(sample["input_pairs"]) != sum(r["sample"] == sample["sample"] for r in rows):
            raise ValueError(f"Unexpected number of inputs: {sample['sample']}")
    paths = [r[k] for r in rows for k in ("R1", "R2")]
    if len(set(paths)) != len(paths):
        raise ValueError("A FASTQ is assigned more than once")
    for row in rows:
        for mate in ("R1", "R2"):
            if not Path(row[mate]).is_file() or not Path(row[mate]).stat().st_size:
                raise ValueError(f"Missing or empty input: {row[mate]}")
        with gzip.open(row["R1"], "rb") as left, gzip.open(row["R2"], "rb") as right:
            for _ in range(100):
                a, b = read_record(left, row["R1"]), read_record(right, row["R2"])
                if a is None or b is None or record_id(a[0]) != record_id(b[0]):
                    raise ValueError(f"Paired prefix validation failed: {row['sample']}")
    refs = [config[k] for k in ("human_bwa", "mito_seed", "sharkmer")]
    refs += [config["human_bwa"] + x for x in (".amb", ".ann", ".bwt", ".pac", ".sa")]
    refs += [str(Path(config["kraken_db"]) / x) for x in ("hash.k2d", "opts.k2d", "taxo.k2d")]
    for path in refs:
        if not Path(path).is_file() or not Path(path).stat().st_size:
            raise ValueError(f"Missing or empty reference/tool: {path}")
    return {
        "specimens": len(ids),
        "sequencing_inputs": len(rows),
        "fastq_files": len(paths),
        "validated_prefix_pairs_per_input": 100,
        "full_fastq_integrity_checked": False,
        "compressed_input_bytes": sum(Path(p).stat().st_size for p in paths),
        "reference_paths": refs,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = check(args.root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
