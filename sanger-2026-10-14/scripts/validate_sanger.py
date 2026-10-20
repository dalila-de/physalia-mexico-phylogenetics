"""Verify retained trace identities, expanded input panels, and final tree scope."""

import csv
import hashlib
import json
from pathlib import Path

from Bio import Phylo, SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent


def main():
    with (ROOT / "provenance/raw_checksums.tsv").open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            assert hashlib.sha256((ROOT / row["file"]).read_bytes()).hexdigest() == row["sha256"], (
                row["file"]
            )
    for gene in ["CO1", "16S"]:
        original = {
            r.id: str(r.seq)
            for r in SeqIO.parse(REPO / "analyses/inputs" / f"{gene}.fasta", "fasta")
        }
        records = list(SeqIO.parse(REPO / "analyses/inputs/combined" / f"{gene}.fasta", "fasta"))
        new = {r.id: str(r.seq) for r in records}
        assert len(new) == len(records) == len(original) + 1
        assert {k: v for k, v in new.items() if k != "sanger__CH887"} == original
        seq = SeqIO.read(ROOT / "inputs" / f"{gene}.sanger.fasta", "fasta")
        assert new[seq.id] == str(seq.seq)
        qc = json.loads((ROOT / "results/qc" / f"{gene}.consensus.json").read_text())
        assert (
            len(seq) == qc["length"]
            and sum(b not in "ACGT" for b in str(seq.seq)) == qc["ambiguous_bases"]
        )
        aln = list(SeqIO.parse(REPO / "analyses/results/alignments" / f"{gene}.fasta", "fasta"))
        assert len({len(r) for r in aln}) == 1 and {r.id for r in aln} == set(new)
        for r in aln:
            assert str(r.seq).replace("-", "").upper() in (
                new[r.id].upper(),
                str(Seq(new[r.id]).reverse_complement()).upper(),
            )
        for suffix in ["treefile", "contree"]:
            tree = Phylo.read(REPO / "analyses/results/trees" / f"{gene}.{suffix}", "newick")
            names = [t.name for t in tree.get_terminals()]
            assert len(names) == len(set(names)) and set(names) == set(new)
    print(
        "Verified 8 original files, consistent expanded reference panels, consensus sequences, alignments and both tree tip sets."
    )


if __name__ == "__main__":
    main()
