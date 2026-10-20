"""Validate the selected tree, root its consensus and report query/reference distances."""

import csv
import json
import re
import sys
from pathlib import Path

from Bio import Phylo, SeqIO


def main(gene: str) -> None:
    root = Path(__file__).resolve().parents[1]
    results = root / "results"
    prefix = results / "trees" / gene
    tree = Phylo.read(str(prefix) + ".treefile", "newick")
    consensus_path = Path(str(prefix) + ".contree")
    if not consensus_path.is_file():
        raise ValueError("Missing bootstrap consensus")
    consensus = Phylo.read(consensus_path, "newick")
    alignment = {r.id for r in SeqIO.parse(results / "alignments" / f"{gene}.fasta", "fasta")}
    panel = root / "inputs" / "combined" if gene in {"16S", "CO1"} else root / "inputs"
    with (panel / f"{gene}.tips.tsv").open() as stream:
        metadata = {r["tip"]: r for r in csv.DictReader(stream, delimiter="\t")}
    for candidate in (tree, consensus):
        names = [t.name for t in candidate.get_terminals()]
        if len(names) != len(set(names)) or set(names) != alignment or set(names) != set(metadata):
            raise ValueError(f"Tree/alignment/metadata tip mismatch: {gene}")
    outgroups = [
        t for t in consensus.get_terminals() if "rhizophysa" in metadata[t.name]["species"].lower()
    ]
    rooting = "midpoint"
    if outgroups and consensus.is_monophyletic(outgroups):
        consensus.root_with_outgroup(consensus.common_ancestor(outgroups))
        rooting = "Rhizophysa_outgroup"
    else:
        consensus.root_at_midpoint()
    consensus.ladderize()
    (results / "figures").mkdir(exist_ok=True)
    Phylo.write(
        consensus,
        results / "figures" / f"{gene}.rooted.nwk",
        "newick",
        format_branch_length="%.17g",
    )
    references = [t for t in tree.get_terminals() if metadata[t.name]["source"] == "NCBI_GenBank"]
    comparison = []
    for tip in tree.get_terminals():
        row = metadata[tip.name]
        if not (row["source"].startswith("October_QC") or row["source"] == "Sanger_AB1"):
            continue
        distances = [(tree.distance(tip, ref), ref) for ref in references]
        best = min(d for d, _ in distances)
        nearest = [ref for distance, ref in distances if abs(distance - best) < 1e-8]
        same = [
            (distance, ref)
            for distance, ref in distances
            if metadata[ref.name]["sample"] == row["sample"]
        ]
        comparison.append(
            {
                "sample": row["sample"],
                "local_tip": tip.name,
                "nearest_reference_distance": best,
                "nearest_reference_tips": ";".join(t.name for t in nearest),
                "nearest_reference_species": ";".join(
                    sorted({metadata[t.name]["species"] for t in nearest})
                ),
                "same_specimen_reference_distance": min(d for d, _ in same)
                if same
                else "not_available",
            }
        )
    report_text = Path(str(prefix) + ".iqtree").read_text()
    models = re.findall(
        r"(?:Best-fit model according to BIC|Model of substitution):\s*(\S+)", report_text
    )
    report = {
        "analysis": gene,
        "tips": len(alignment),
        "local_tips": len(comparison),
        "rooting": rooting,
        "model": models[0] if models else "see_iqtree_report",
        "interpretation": "Nearest-reference distances describe placement, not definitive species assignment.",
        "specimens": comparison,
    }
    (results / "tree_reports").mkdir(exist_ok=True)
    (results / "tree_reports" / f"{gene}.json").write_text(json.dumps(report, indent=2) + "\n")
    with (results / "tree_reports" / f"{gene}.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparison[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(comparison)


if __name__ == "__main__":
    main(sys.argv[1])
