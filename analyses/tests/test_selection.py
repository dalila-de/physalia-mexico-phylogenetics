"""Selection must choose one record, not concatenate split mitochondrial genes."""

import gzip

import prepare_tree_inputs as prepare
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def test_longest_mitos_fragment_competes_with_sharkmer(tmp_path, monkeypatch):
    upstream = tmp_path / "upstream"
    qc = tmp_path / "qc"
    result = tmp_path / "results"
    for p in [upstream, qc / "test/markers", result / "annotations/test"]:
        p.mkdir(parents=True)
    (upstream / "accessions.tsv").write_text(
        "marker\tsample\tspecies\taccession\textraction\nCO1\treference\tPhysalia physalis\tTEST.1\twhole_record\n"
    )
    record = SeqRecord(Seq("ACGTACGT"), id="TEST.1", name="TEST", description="test")
    record.annotations.update(molecule_type="DNA", organism="Physalia physalis")
    record.features = [
        SeqFeature(FeatureLocation(0, 8), type="source", qualifiers={"isolate": ["reference"]})
    ]
    with gzip.open(upstream / "records.gb.gz", "wt") as f:
        SeqIO.write(record, f, "genbank")
    (result / "annotations/test/result.fas").write_text(
        ">scaffold; 1-6; +; cox1_0\nACGTAC\n"
        ">scaffold; 7-11; +; cox1_1\nACGTA\n"
        ">scaffold; 12-19; +; cox2\nACGTACGT\n"
    )
    (qc / "test/markers/test_cnidaria_CO1.fasta").write_text(">test product 0\nACGT\n")
    for name, value in [
        ("REFERENCE", upstream),
        ("QC", qc),
        ("RESULTS", result),
        ("ANALYSIS", tmp_path),
        ("QC_IDS", ["test"]),
    ]:
        monkeypatch.setattr(prepare, name, value)
    prepare.prepare("CO1")
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(tmp_path / "inputs/CO1.fasta", "fasta")}
    assert seqs["local__test"] == "ACGTAC"
