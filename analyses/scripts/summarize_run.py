"""Summarize retained NCBI-based inference, search agreement and branch support."""

import json
import re

from Bio import Phylo
from common import ANALYSIS, table


def main():
    results = ANALYSIS / "results"
    summaries = {}
    for gene in ("CO1", "16S", "18S", "ITS", "mitogenome_submitted"):
        prefix = results / "trees" / gene
        text = prefix.with_suffix(".iqtree").read_text()
        block = text.split("MULTIPLE RUNS", 1)[1].split("MAXIMUM LIKELIHOOD", 1)[0]
        likelihoods = [float(v) for v in re.findall(r"^\d+\s+(-[0-9.]+)", block, re.M)]
        if len(likelihoods) < 5:
            raise ValueError(f"Incomplete replicate searches: {gene}")
        panel = ANALYSIS / "inputs" / ("combined" if gene in {"CO1", "16S"} else "")
        metadata = {r["tip"]: r for r in table(panel / f"{gene}.tips.tsv")}
        tree = Phylo.read(results / f"figures/{gene}.rooted.nwk", "newick")
        refs = [
            r["tip"]
            for r in metadata.values()
            if r["source"] == "NCBI_GenBank" and r["species"] == "Physalia physalis"
        ]
        focal = tree.common_ancestor(refs)
        new = [
            r["tip"]
            for r in metadata.values()
            if r["source"].startswith("October_QC") or r["source"] == "Sanger_AB1"
        ]
        tips = {t.name for t in focal.get_terminals()}
        local = []
        for tip in new:
            node = next(n for n in tree.get_terminals() if n.name == tip)
            parent = tree.get_path(node)[-2]
            local.append(
                dict(
                    tip=tip,
                    within_reference_clade=tip in tips,
                    immediate_group_support=parent.confidence,
                    immediate_group_tips=[n.name for n in parent.get_terminals()],
                )
            )
        summaries[gene] = dict(
            **json.loads((results / f"alignments/{gene}.json").read_text()),
            **{
                k: v
                for k, v in json.loads((results / f"tree_reports/{gene}.json").read_text()).items()
                if k in {"model", "rooting"}
            },
            searches=len(likelihoods),
            run_log_likelihoods=likelihoods,
            best_likelihood_repeats_within_0_001=sum(
                abs(v - max(likelihoods)) <= 0.001 for v in likelihoods
            ),
            physalis_reference_clade_support=focal.confidence,
            all_new_tips_within_reference_clade=set(new) <= tips,
            nonfocal_reference_species=sorted(
                {
                    metadata[t]["species"]
                    for t in tips
                    if metadata[t]["source"] == "NCBI_GenBank"
                    and metadata[t]["species"] != "Physalia physalis"
                }
            ),
            local_placements=local,
        )
    (results / "ncbi_run_summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
    print(
        json.dumps(
            {
                g: {k: v for k, v in s.items() if k != "local_placements"}
                for g, s in summaries.items()
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
