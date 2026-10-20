"""Add the already validated CH887 sequences to the new NCBI/Mexico panels."""

from Bio import SeqIO
from common import ANALYSIS, ROOT, table, write_table


def main():
    for gene in ("CO1", "16S"):
        panel = list(SeqIO.parse(ANALYSIS / "inputs" / f"{gene}.fasta", "fasta"))
        sanger = SeqIO.read(ROOT / "sanger-2026-10-14/inputs" / f"{gene}.sanger.fasta", "fasta")
        if sanger.id in {r.id for r in panel}:
            raise ValueError("Duplicate CH887 tip")
        rows = table(ANALYSIS / "inputs" / f"{gene}.tips.tsv")
        import hashlib

        rows.append(
            dict(
                tip=sanger.id,
                sample="CH887",
                source="Sanger_AB1",
                species="Physalia physalis",
                original_header="CH887-D01" if gene == "CO1" else "CH887A-D01",
                sequence_sha256=hashlib.sha256(str(sanger.seq).encode()).hexdigest(),
                accession="",
                extraction="Sanger_AB1",
                feature="",
            )
        )
        output = ANALYSIS / "inputs/combined"
        output.mkdir(exist_ok=True)
        SeqIO.write(panel + [sanger], output / f"{gene}.fasta", "fasta")
        write_table(output / f"{gene}.tips.tsv", rows)


if __name__ == "__main__":
    main()
