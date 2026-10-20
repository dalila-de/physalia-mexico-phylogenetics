"""Prepare accession-derived NCBI references and unchanged Mexican marker candidates."""

import argparse
import hashlib
import json
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from common import ANALYSIS, GENES, QC, QC_IDS, RESULTS, ROOT, records, table, ungap, write_table
from ncbi_references import extract_sequence, load_records, verify_identity

REFERENCE = ROOT / "references/ncbi"


def clean_tip(identifier: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", identifier)


def prepare(gene: str) -> None:
    manifest = [r for r in table(REFERENCE / "accessions.tsv") if r["marker"] == gene]
    downloaded = load_records(REFERENCE / "records.gb.gz")
    rows, candidates, seqs = [], [], []
    local_ids = QC_IDS
    for entry in manifest:
        record = downloaded[entry["accession"]]
        verify_identity(record, entry)
        sequence, feature = extract_sequence(record, entry["extraction"])
        sequence = ungap(sequence)
        if not sequence or set(sequence) - set("ACGTRYSWKMBDHVN"):
            raise ValueError(f"Invalid reference sequence: {entry}")
        name = "ref__" + clean_tip(entry["sample"])
        seqs.append(SeqRecord(Seq(sequence), id=name, description=""))
        rows.append(
            dict(
                tip=name,
                sample=entry["sample"],
                source="NCBI_GenBank",
                species=entry["species"],
                original_header=record.description,
                sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                accession=record.id,
                extraction=entry["extraction"],
                feature=feature,
            )
        )
    if gene == "mitogenome_submitted":
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
    expected = {r["sample"] for r in manifest} | set(local_ids)
    if len(rows) != len(expected) or {r["sample"] for r in rows} != expected:
        raise ValueError(f"Missing or duplicate specimens: {gene}")
    if len({r.id for r in seqs}) != len(seqs):
        raise ValueError(f"Duplicate tip IDs: {gene}")
    output = ANALYSIS / "inputs"
    output.mkdir(exist_ok=True)
    SeqIO.write(seqs, output / f"{gene}.fasta", "fasta")
    write_table(output / f"{gene}.tips.tsv", rows)
    if candidates:
        write_table(output / f"{gene}.candidates.tsv", candidates)
    (output / f"{gene}.json").write_text(
        json.dumps(
            dict(
                analysis=gene,
                tips=len(seqs),
                local_tips=len(local_ids),
                reference_tips=len(manifest),
                reference_source="NCBI GenBank accession manifest",
                alignment="de_novo",
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene", choices=[*GENES, "mitogenome_submitted"])
    prepare(parser.parse_args().gene)
