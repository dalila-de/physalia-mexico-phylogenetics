"""Validate expanded marker trees and root consensus trees for shared plotting."""

import csv
import json
from pathlib import Path

from Bio import Phylo, SeqIO

ROOT = Path(__file__).resolve().parents[2] / "analyses"


def main():
    rows = []
    for gene in ["CO1", "16S"]:
        with (ROOT / "inputs/combined" / f"{gene}.tips.tsv").open() as stream:
            tips = {r["tip"]: r for r in csv.DictReader(stream, delimiter="\t")}
        records = list(SeqIO.parse(ROOT / "results/alignments" / f"{gene}.fasta", "fasta"))
        if len(records) != len(tips) or {r.id for r in records} != set(tips):
            raise ValueError("Alignment IDs differ")
        trees = [
            Phylo.read(ROOT / "results/trees" / f"{gene}.{suffix}", "newick")
            for suffix in ["treefile", "contree"]
        ]
        for tree in trees:
            labels = [t.name for t in tree.get_terminals()]
            if len(labels) != len(set(labels)) or set(labels) != set(tips):
                raise ValueError("Tree IDs differ")
        ml, consensus = trees
        outgroups = [
            t for t in consensus.get_terminals() if "rhizophysa" in tips[t.name]["species"].lower()
        ]
        if not outgroups or not consensus.is_monophyletic(outgroups):
            raise ValueError("Rhizophysa outgroup split absent; review rooting")
        consensus.root_with_outgroup(consensus.common_ancestor(outgroups))
        consensus.ladderize()
        (ROOT / "results/figures").mkdir(exist_ok=True)
        Phylo.write(consensus, ROOT / "results/figures" / f"{gene}.rooted.nwk", "newick")
        query = next(t for t in ml.get_terminals() if t.name == "sanger__CH887")
        closest = sorted(
            (ml.distance(query, t), t.name) for t in ml.get_terminals() if t is not query
        )
        best = closest[0][0]
        row = dict(
            gene=gene,
            tips=len(tips),
            alignment_sites=len(records[0]),
            rooting="Rhizophysa",
            nearest_tip_distance=best,
            nearest_tips=[n for d, n in closest if abs(d - best) < 1e-8],
            interpretation="Marker-tree placement does not assign genetic E/W ancestry or prove specimen identity.",
        )
        rows.append(row)
    (ROOT / "results/sanger_tree_summary.json").write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
