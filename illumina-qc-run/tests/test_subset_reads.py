"""Checks for consequential mate matching and multi-run allocation behavior."""

import gzip
from pathlib import Path

import pytest
from subset_reads import select_pairs


def make_unit(tmp_path: Path, unit: str, n: int, bad_mate: bool = False) -> dict:
    row = {"sample": "test", "unit": unit}
    for mate in (1, 2):
        path = tmp_path / f"{unit}_R{mate}.fq.gz"
        with gzip.open(path, "wt") as stream:
            for i in range(n):
                name = f"{unit}-{i}" if not (bad_mate and mate == 2) else "wrong"
                stream.write(f"@{name}/{mate}\nACGT\n+\nIIII\n")
        row[f"R{mate}"] = str(path)
    return row


def test_two_runs_share_the_per_specimen_cap(tmp_path):
    units = [make_unit(tmp_path, x, 10) for x in ("a", "b")]
    counts = select_pairs(units, 6, tmp_path / "out1", tmp_path / "out2")
    assert [r["selected_pairs"] for r in counts] == [3, 3]
    assert (tmp_path / "out1").read_text().count("\n+") == 6
    assert "@b-0/1" in (tmp_path / "out1").read_text()


def test_rejects_mismatched_mates(tmp_path):
    unit = make_unit(tmp_path, "a", 2, bad_mate=True)
    with pytest.raises(ValueError, match="Mismatched mates"):
        select_pairs([unit], 2, tmp_path / "out1", tmp_path / "out2")


def test_short_input_records_actual_count(tmp_path):
    unit = make_unit(tmp_path, "a", 2)
    counts = select_pairs([unit], 10, tmp_path / "out1", tmp_path / "out2")
    assert counts[0]["selected_pairs"] == 2


def test_rejects_truncated_record(tmp_path):
    unit = make_unit(tmp_path, "a", 2)
    with gzip.open(unit["R2"], "wt") as stream:
        stream.write("@a-0/2\nACGT\n+\nI\n")
    with pytest.raises(ValueError, match="Malformed FASTQ"):
        select_pairs([unit], 2, tmp_path / "out1", tmp_path / "out2")
