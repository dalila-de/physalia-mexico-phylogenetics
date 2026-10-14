"""Selection must choose one record, not concatenate split mitochondrial genes."""

import prepare_tree_inputs as prepare
from Bio import SeqIO


def test_longest_mitos_fragment_competes_with_sharkmer(tmp_path, monkeypatch):
    upstream = tmp_path / "upstream"
    qc = tmp_path / "qc"
    result = tmp_path / "results"
    for p in [upstream / "gene_trees", qc / "test/markers", result / "annotations/test"]:
        p.mkdir(parents=True)
    (upstream / "tip_metadata.tsv").write_text("seq\tspecies\nreference\tPhysalia physalis\n")
    (upstream / "gene_trees/CO1.all.fasta").write_text(">reference\nACGTACGT\n")
    (result / "annotations/test/result.fas").write_text(
        ">scaffold; 1-6; +; cox1_0\nACGTAC\n"
        ">scaffold; 7-11; +; cox1_1\nACGTA\n"
        ">scaffold; 12-19; +; cox2\nACGTACGT\n"
    )
    (qc / "test/markers/test_cnidaria_CO1.fasta").write_text(">test product 0\nACGT\n")
    for name, value in [
        ("CONFIG", {}),
        ("UPSTREAM", upstream),
        ("QC", qc),
        ("RESULTS", result),
        ("ANALYSIS", tmp_path),
        ("IDS", ["test"]),
    ]:
        monkeypatch.setattr(prepare, name, value)
    prepare.prepare("CO1")
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(tmp_path / "inputs/CO1.fasta", "fasta")}
    assert seqs["local__test"] == "ACGTAC"
