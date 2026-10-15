"""Check denominators, explicit missing assemblies, and complete cohort joins."""

import json
import zipfile

import pytest
from retain_results import cohort, retain_sample


def make_outputs(tmp_path):
    sample = "test"
    work = tmp_path / "work"

    def write(name, contents):
        path = work / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(contents)

    write("reads/test.json", json.dumps([{"selected_pairs": 10}]))
    write(
        "logs/trimmomatic/test.log",
        "Input Read Pairs: 10 Both Surviving: 8 (80%) Forward Only Surviving: 1 (10%) Reverse Only Surviving: 0 (0%) Dropped: 1 (10%)\n",
    )
    write("human/test.mapped_paired", "2\n")
    write("kraken/test.report.txt", "80.0\t8\t8\tU\t0\tunclassified\n")
    write("getorganelle/test_mtgenome/get_org.log.txt", "No assembly paths\n")
    write("shark_results/test/test.stats", "test statistics\n")
    write("shark_results/test/test_cnidaria_16S.fasta", ">test\nACGT\nACGT\n")
    for primer in ["16sF", "16sR", "co1F", "co1R", "18sF", "18sR", "ITSfullF"]:
        write(f"rego_results/test/test_{primer}.txt", "Hits found: 0\n")
    for mate in ["R1", "R2"]:
        path = work / "fastqc/test" / f"test.{mate}.trimmed_fastqc.zip"
        path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("test/summary.txt", "PASS\tBasic Statistics\ttest.fastq\n")
    return work, sample


def test_retains_multiline_markers_and_no_assembly_status(tmp_path, monkeypatch):
    work, sample = make_outputs(tmp_path)
    monkeypatch.chdir(tmp_path)
    out = tmp_path / "results"
    retain_sample(work, out, sample, 10)
    result = json.loads((out / sample / "summary.json").read_text())
    assert result["human_percent_reads"] == 12.5  # 2 records / (2 * 8 pairs)
    assert result["mitogenome_status"] == "no_path_recovered"
    assert result["marker_files"] == 1
    assert json.loads((out / sample / "sequence_lengths.json").read_text())[
        "test_cnidaria_16S.fasta"
    ] == [8]


def test_rejects_inconsistent_selected_count(tmp_path, monkeypatch):
    work, sample = make_outputs(tmp_path)
    (work / "reads/test.json").write_text(json.dumps([{"selected_pairs": 9}]))
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="Subsampling/trimming count mismatch"):
        retain_sample(work, tmp_path / "results", sample, 10)
    assert not (tmp_path / "results/test/summary.json").exists()


def test_cohort_does_not_silently_drop_missing_sample(tmp_path):
    manifest = tmp_path / "samples.tsv"
    manifest.write_text("sample\ntest\nmissing\n")
    report = tmp_path / "report.json"
    report.write_text(json.dumps({"sample": "test"}))
    with pytest.raises(ValueError, match="exactly the approved specimens"):
        cohort(manifest, [report], tmp_path / "summary.tsv")
