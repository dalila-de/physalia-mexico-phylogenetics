"""Check tree scope and plot reference species with highlighted October specimens."""

import csv
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
    with (root / "inputs" / f"{gene}.tips.tsv").open() as stream:
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
    palette = {
        "Physalia physalis": "#2166ac",
        "Physalia utriculus": "#1b9e77",
        "Physalia megalista": "#984ea3",
        "Physalia minuta": "#e6ab02",
        "Physalia glaucus": "#e7298a",
    }
    fig, ax = plt.subplots(figsize=(16, max(12, len(alignment) * 0.16)))
    labels = {}
    colors = {}
    for tip, row in metadata.items():
        local = row["source"].startswith("October_QC")
        mexico = row["sample"] in {f"YPM-IZ-{i}" for i in range(116693, 116698)}
        label = (
            f"{row['sample']} [OCTOBER{' MEXICO' if mexico else ''}]"
            if local
            else f"{row['sample']} | {row['species'] or 'unassigned'} [reference]"
        )
        # Preserve unique display labels for alternative mitochondrial paths.
        if gene.startswith("mitogenome"):
            label += " " + tip.removeprefix("ref__").removeprefix("local__")
        labels[tip] = label
        colors[label] = (
            "#d73027" if mexico else "#000000" if local else palette.get(row["species"], "#777777")
        )
    Phylo.draw(
        consensus,
        axes=ax,
        do_show=False,
        label_func=lambda c: labels.get(c.name),
        label_colors=colors,
        show_confidence=True,
    )
    ax.set_title(
        f"{gene}: October assemblies and published reference sequences\n{rooting}; ultrafast bootstrap 1000; black = October, red = Mexico"
    )
    (results / "figures").mkdir(exist_ok=True)
    fig.savefig(results / "figures" / f"{gene}.pdf", bbox_inches="tight")
    plt.close(fig)
    Phylo.write(consensus, results / "figures" / f"{gene}.rooted.nwk", "newick")
    references = [t for t in tree.get_terminals() if metadata[t.name]["source"] == "repository"]
    comparison = []
    for tip in tree.get_terminals():
        row = metadata[tip.name]
        if not row["source"].startswith("October_QC"):
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
