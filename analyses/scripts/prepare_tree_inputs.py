"""Reproduce upstream longest-marker selection with October assemblies as new tips."""

import argparse
import hashlib
import json
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from common import (
    ANALYSIS,
    CONFIG,
    GENES,
    IDS,
    QC,
    RESULTS,
    SAMPLES,
    UPSTREAM,
    records,
    specimen,
    table,
    ungap,
    write_table,
)


def clean_tip(identifier: str) -> str:
    # Upstream retains alternative paths distinguished only by its _R_ tag.
    # Preserve that original tag inside our ref__ namespace; strip only tags
    # introduced by the new MAFFT run during alignment validation.
    return re.sub(r"[^A-Za-z0-9_.-]", "_", identifier)


def prepare(gene: str) -> None:
    metadata = {r["seq"]: r for r in table(UPSTREAM / "tip_metadata.tsv")}
    rows, candidates, seqs = [], [], []
    restricted = CONFIG.get("reference_scope") == "published_cohort"
    published = {r["sample"] for r in SAMPLES if r["cohort"] == "published"}
    additional = set(CONFIG.get("additional_reference_samples", []))
    local_ids = [sid for sid in IDS if not restricted or sid not in published]

    def keep_reference(record):
        info = metadata.get(record.id, metadata.get(specimen(record.id), {}))
        return (
            not restricted
            or specimen(record.id) in published | additional
            or "rhizophysa" in info.get("species", "").lower()
        )

    if gene.startswith("mitogenome"):
        filename = (
            "submitted_rooted.aln.fa" if gene == "mitogenome_submitted" else "identification.aln.fa"
        )
        source = UPSTREAM / "mitogenome_tree" / filename
        reference = records(source)
        if gene == "mitogenome_all":
            reference += records(UPSTREAM / "mitogenome_tree/rhizophysa_outgroups.fasta")
        for r in reference:
            if not keep_reference(r):
                continue
            name = "ref__" + clean_tip(r.id)
            sid = specimen(r.id)
            info = metadata.get(r.id, metadata.get(sid, {}))
            seqs.append(SeqRecord(Seq(ungap(r.seq)), id=name, description=""))
            rows.append(
                {
                    "tip": name,
                    "sample": sid,
                    "source": "repository",
                    "species": info.get("species", ""),
                    "original_header": r.description,
                    "sequence_sha256": hashlib.sha256(ungap(r.seq).encode()).hexdigest(),
                }
            )
        for sid in local_ids:
            for path in sorted((QC / sid / "mitogenomes").glob("*.fasta")):
                for index, r in enumerate(records(path), 1):
                    name = f"local__{sid}__path{index}"
                    seqs.append(SeqRecord(Seq(ungap(r.seq)), id=name, description=""))
                    rows.append(
                        {
                            "tip": name,
                            "sample": sid,
                            "source": "October_QC",
                            "species": "",
                            "original_header": r.description,
                            "sequence_sha256": hashlib.sha256(ungap(r.seq).encode()).hexdigest(),
                        }
                    )
    else:
        for r in records(UPSTREAM / "gene_trees" / f"{gene}.all.fasta"):
            if not keep_reference(r):
                continue
            name = "ref__" + clean_tip(r.id)
            info = metadata.get(r.id, {})
            seqs.append(SeqRecord(Seq(ungap(r.seq)), id=name, description=""))
            rows.append(
                {
                    "tip": name,
                    "sample": r.id,
                    "source": "repository",
                    "species": info.get("species", ""),
                    "original_header": r.description,
                    "sequence_sha256": hashlib.sha256(ungap(r.seq).encode()).hexdigest(),
                }
            )
        for sid in local_ids:
            options = []
            if gene in {"16S", "CO1"}:
                mitos_gene = "rrnL" if gene == "16S" else "cox1"
                path = RESULTS / "annotations" / sid / "result.fas"
                for r in records(path):
                    if re.search(rf"(?:^|[;\s]){mitos_gene}(?:_\d+)?(?:[;\s]|$)", r.description):
                        options.append(("MITOS2", path, r))
            path = QC / sid / "markers" / f"{sid}_{GENES[gene]}.fasta"
            if path.exists():
                options += [("sharkmer", path, r) for r in records(path)]
            if not options:
                raise ValueError(f"No candidates for {sid} {gene}")
            selected = max(range(len(options)), key=lambda i: len(options[i][2]))
            for i, (origin, path, r) in enumerate(options):
                candidates.append(
                    {
                        "sample": sid,
                        "marker": gene,
                        "origin": origin,
                        "path": str(path),
                        "header": r.description,
                        "length": len(r),
                        "selected": i == selected,
                    }
                )
            origin, path, r = options[selected]
            name = "local__" + sid
            seqs.append(SeqRecord(Seq(ungap(r.seq)), id=name, description=""))
            rows.append(
                {
                    "tip": name,
                    "sample": sid,
                    "source": "October_QC_" + origin,
                    "species": "",
                    "original_header": r.description,
                    "sequence_sha256": hashlib.sha256(ungap(r.seq).encode()).hexdigest(),
                }
            )
    if restricted:
        expected = set(IDS) | additional
        if gene == "mitogenome_submitted":
            expected.remove("YPM-IZ-110972")
        ingroup = [r["sample"] for r in rows if "rhizophysa" not in r["species"].lower()]
        if set(ingroup) != expected or len(ingroup) != len(expected):
            raise ValueError(f"Missing, extra or duplicated cohort specimens: {gene}: {ingroup}")
    names = [r.id for r in seqs]
    if len(names) != len(set(names)):
        raise ValueError(f"Tip collision after normalization: {gene}")
    output = ANALYSIS / "inputs"
    output.mkdir(exist_ok=True)
    SeqIO.write(seqs, output / f"{gene}.fasta", "fasta")
    write_table(output / f"{gene}.tips.tsv", rows)
    if candidates:
        write_table(output / f"{gene}.candidates.tsv", candidates)
    if gene.startswith("mitogenome"):
        # Preserve the upstream alignment and add only new assemblies/outgroups.
        base = [r for r in records(source) if keep_reference(r)]
        for r in base:
            r.id = "ref__" + clean_tip(r.id)
            r.description = ""
        SeqIO.write(base, output / f"{gene}.base.fasta", "fasta")
        existing = {r.id for r in base}
        SeqIO.write(
            [r for r in seqs if r.id not in existing], output / f"{gene}.add.fasta", "fasta"
        )
    else:
        # These unaligned partitions document reference/query provenance; marker
        # MAFFT consumes the combined FASTA and aligns it de novo.
        SeqIO.write(
            [r for r in seqs if r.id.startswith("ref__")], output / f"{gene}.base.fasta", "fasta"
        )
        SeqIO.write(
            [r for r in seqs if r.id.startswith("local__")], output / f"{gene}.add.fasta", "fasta"
        )
    (output / f"{gene}.json").write_text(
        json.dumps(
            {
                "analysis": gene,
                "tips": len(seqs),
                "local_tips": sum(r["source"].startswith("October_QC") for r in rows),
                "reference_tips": sum(r["source"] == "repository" for r in rows),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene", choices=[*GENES, "mitogenome_submitted", "mitogenome_all"])
    prepare(parser.parse_args().gene)
