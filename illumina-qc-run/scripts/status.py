"""List every selected specimen, including samples without completed reports."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import yaml


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-pilot", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    config = yaml.safe_load((root / "config/config.yaml").read_text())
    with (root / config["samples"]).open() as stream:
        samples = list(csv.DictReader(stream, delimiter="\t"))
    print("sample\tstatus\tmitogenome_status\tmarker_files")
    incomplete_pilot = []
    for row in samples:
        sample = row["sample"]
        folder = root / "results" / config["run_id"] / sample
        path = folder / "summary.json"
        if path.exists():
            summary = json.loads(path.read_text())
            if (
                args.require_pilot
                and sample in config["pilot"]
                and summary["getorganelle_error_lines"]
            ):
                raise ValueError(f"Pilot has GetOrganelle errors requiring inspection: {sample}")
            with (folder / "checksums.tsv").open() as stream:
                for artifact in csv.DictReader(stream, delimiter="\t"):
                    digest = hashlib.sha256((folder / artifact["file"]).read_bytes()).hexdigest()
                    if digest != artifact["sha256"]:
                        raise ValueError(
                            f"Retained result checksum mismatch: {sample}/{artifact['file']}"
                        )
            print(
                f"{sample}\t{summary['workflow_status']}\t{summary['mitogenome_status']}\t{summary['marker_files']}"
            )
        else:
            print(f"{sample}\tnot_complete\tunassessed\tNA")
            if sample in config["pilot"]:
                incomplete_pilot.append(sample)
    if args.require_pilot and incomplete_pilot:
        raise ValueError(f"Pilot incomplete: {incomplete_pilot}")


if __name__ == "__main__":
    main()
