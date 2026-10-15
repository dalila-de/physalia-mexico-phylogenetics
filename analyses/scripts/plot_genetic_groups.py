"""Plot existing trees with prior genome-wide genetic E/W groups."""

import argparse
import csv
import math
from copy import deepcopy
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
PALETTE = {
    "E": "#0072B2",
    "W": "#D55E00",
    "Mexico": "#CC79A7",
    "utriculus": "#009E73",
    "outgroup": "#555555",
    "Sanger": "#7851A9",
}
LABELS = {
    "E": "E — genetic group",
    "W": "W — genetic group",
    "Mexico": "Gulf (Progreso) — separate category",
    "utriculus": "P. utriculus — E/W not applicable",
    "outgroup": "Rhizophysa outgroup",
}


def read_table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main(gene, root=ROOT):
    annotations = {r["sample"]: r for r in read_table(root / "inputs/genetic_annotations.tsv")}
    panel = root / "inputs" / "combined" if gene in {"16S", "CO1"} else root / "inputs"
    tips = {r["tip"]: r for r in read_table(panel / f"{gene}.tips.tsv")}
    tree = Phylo.read(root / f"results/figures/{gene}.rooted.nwk", "newick")
    if {t.name for t in tree.get_terminals()} != set(tips):
        raise ValueError("Tree and tip metadata do not match")
    labels, colors = {}, {}
    joined = []
    for tip, row in tips.items():
        a = annotations[row["sample"]]
        group = a["plot_group"]
        label = f"{row['sample']} | {a['location'] or row['species']} [{group}]"
        if row["source"].startswith("October_QC"):
            label += " *new*"
        if a["vagrant"] == "True":
            label += " [VAGRANT]"
        labels[tip] = label
        colors[label] = PALETTE[group]
        joined.append(dict(tip=tip, **a))
    out = root / "results/genetic_groups"
    out.mkdir(exist_ok=True)
    with (out / f"{gene}.annotations.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(joined[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(joined)
    focal = deepcopy(tree)
    for tip in list(focal.get_terminals()):
        if annotations[tips[tip.name]["sample"]]["plot_group"] in {"utriculus", "outgroup"}:
            focal.prune(tip)
    # Remove the stem to excluded taxa; focal pairwise distances are unchanged.
    focal.root.branch_length = 0
    root_note = "Midpoint root" if gene == "ITS" else "Original root: Rhizophysa"
    for suffix, subtitle, view in [
        ("", " — all tips", tree),
        (".focal", " — focal view (references pruned for display)", focal),
    ]:
        leaves = view.get_terminals()
        y = {tip: len(leaves) - index for index, tip in enumerate(leaves)}
        depths = view.depths()
        topology_only = max(depths.values()) == 0
        if topology_only:
            depths = view.depths(unit_branch_lengths=True)
        for node in view.get_nonterminals(order="postorder"):
            y[node] = (y[node.clades[0]] + y[node.clades[-1]]) / 2
        xmax = max(depths.values()) or 1
        fig, ax = plt.subplots(figsize=(12, max(8, len(leaves) * 0.22 + 1.7)))
        fig.subplots_adjust(left=0.04, right=0.98, top=0.92, bottom=0.10)
        for node in view.find_clades():
            if node.clades:
                ax.plot(
                    [depths[node]] * 2,
                    [y[node.clades[0]], y[node.clades[-1]]],
                    color="#505861",
                    linewidth=0.8,
                )
                for child in node.clades:
                    ax.plot(
                        [depths[node], depths[child]],
                        [y[child]] * 2,
                        color="#505861",
                        linewidth=0.8,
                    )
            if not node.is_terminal() and node is not view.root and node.confidence is not None:
                if node.confidence >= 70:
                    ax.scatter(
                        depths[node],
                        y[node],
                        s=15,
                        zorder=4,
                        facecolor="#505861" if node.confidence >= 95 else "white",
                        edgecolor="#505861",
                        linewidth=0.7,
                    )
        for tip in leaves:
            row = tips[tip.name]
            annotation = annotations[row["sample"]]
            group = annotation["plot_group"]
            color = PALETTE[group]
            vagrant = annotation["vagrant"] == "True"
            location = annotation["location"] or row["species"].replace("Rhizophysa ", "R. ")
            name = row["sample"] + (" †" if vagrant else "")
            label = f"{name}  ·  {location}"
            ax.annotate(
                label,
                (depths[tip], y[tip]),
                xytext=(5, 0),
                textcoords="offset points",
                color=color,
                fontsize=9,
                va="center",
                weight="bold" if group in {"Mexico", "Sanger"} or vagrant else "normal",
            )
        # Reserve space for labels while keeping their position tied to each tip.
        ax.set_xlim(-xmax * 0.02, xmax * 1.65)
        ax.set_ylim(-0.9, len(leaves) + 0.8)
        ax.axis("off")
        raw_scale = xmax / 5
        unit = 10 ** math.floor(math.log10(raw_scale))
        scale = max(v for v in (1, 2, 5, 10) if v <= raw_scale / unit) * unit
        ax.plot([0, scale], [-0.15, -0.15], color="#505861", linewidth=1.4)
        scale_label = f"{scale:g} " + (
            "topological steps" if topology_only else "substitutions / site"
        )
        ax.text(scale / 2, -0.6, scale_label, ha="center", va="top", fontsize=8, color="#64707d")
        title = "Mitogenome" if gene == "mitogenome_submitted" else "COI" if gene == "CO1" else gene
        fig.text(0.045, 0.968, title, fontsize=21, weight="bold", color="#23313f")
        scope = (
            "All specimens and outgroups"
            if not suffix
            else "Focal specimens · outgroups pruned for display"
        )
        fig.text(
            0.045,
            0.941,
            f"{scope}  |  {len(leaves)} tips  |  prior genetic E/W assignments",
            fontsize=10,
            color="#64707d",
        )
        present = {annotations[tips[t.name]["sample"]]["plot_group"] for t in leaves}
        handles = [
            Line2D(
                [],
                [],
                marker="o",
                linestyle="",
                markersize=6,
                color=PALETTE[g],
                label={
                    "E": "Genetic E",
                    "W": "Genetic W",
                    "Mexico": "Gulf (Progreso)",
                    "utriculus": "P. utriculus",
                    "outgroup": "Rhizophysa",
                    "Sanger": "Sanger — E/W unassigned",
                }[g],
            )
            for g in PALETTE
            if g in present
        ]
        fig.legend(
            handles=handles,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.033),
            ncol=len(handles),
            frameon=False,
            fontsize=9,
        )
        fig.text(
            0.5,
            0.013,
            f"† Known vagrant    • Bootstrap ≥95    ○ Bootstrap 70–94    |    {root_note}",
            ha="center",
            fontsize=8,
            color="#64707d",
        )
        for ext in ("pdf", "png", "svg"):
            fig.savefig(out / f"{gene}{suffix}.{ext}", dpi=180, facecolor="white")
        plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene")
    parser.add_argument("--analysis-root", type=Path, default=ROOT)
    args = parser.parse_args()
    main(args.gene, args.analysis_root)
