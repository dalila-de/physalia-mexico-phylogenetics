"""Write a QC report atomically, with explicit read-count denominators."""

import re
import sys
import zipfile
from collections.abc import Iterator
from pathlib import Path


def fasta_lengths(path: Path) -> Iterator[tuple[str, int]]:
    name, length = None, 0
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, length
            name, length = line[1:].split()[0], 0
        else:
            length += len(line.strip())
    if name is not None:
        yield name, length


def write_report(sample: str, output: Path, requested: int) -> None:
    lines = [
        "## ILLUMINA QC SUMMARY ##",
        f"ID: {sample}",
        f"Requested maximum read pairs: {requested}",
        "",
        "### FastQC Results ###",
    ]
    for mate in ("R1", "R2"):
        path = Path("fastqc") / sample / f"{sample}.{mate}.trimmed_fastqc.zip"
        with zipfile.ZipFile(path) as archive:
            members = [n for n in archive.namelist() if n.endswith("/summary.txt")]
            if len(members) != 1:
                raise ValueError(f"Missing or ambiguous FastQC summary: {path}")
            lines += [f"**{mate} Reads:**", archive.read(members[0]).decode().rstrip()]
    kraken = (Path("kraken") / f"{sample}.report.txt").read_text().splitlines()
    lines += ["", "### Kraken2 Results ###"]
    for label, taxon in [("Unclassified", "unclassified"), ("Bacterial", "Bacteria")]:
        matches = [
            line.split("\t")[0].strip() for line in kraken if line.split("\t")[-1].strip() == taxon
        ]
        if not matches and taxon == "unclassified":
            raise ValueError("Missing Kraken unclassified row")
        lines.append(f"{label} read pairs: {matches[0] if matches else '0.00'}%")
    trim = (Path("logs/trimmomatic") / f"{sample}.log").read_text()
    pairs = int(re.findall(r"Both Surviving: (\d+)", trim)[-1])
    mapped = int((Path("human") / f"{sample}.mapped_paired").read_text())
    if pairs <= 0 or not 0 <= mapped <= 2 * pairs:
        raise ValueError("Invalid human mapping numerator or paired-read denominator")
    lines += [
        "",
        "### Human Contamination Results ###",
        f"Read pairs surviving trimming: {pairs}",
        f"Read records passing human-mapping filters: {mapped}",
        f"Percent reads human: {100 * mapped / (2 * pairs):.5f}%",
        "",
        "### GetOrganelle Results ###",
    ]
    organelle = Path("getorganelle") / f"{sample}_mtgenome"
    log = (organelle / "get_org.log.txt").read_text()
    assemblies = sorted(organelle.glob("*.path_sequence.fasta"))
    if not assemblies:
        lines.append("No mitochondrial assembly path produced.")
    for assembly in assemblies:
        lines.append(f"Assembly: {assembly}")
        lines.extend(f"    {name}: {length} bp" for name, length in fasta_lengths(assembly))
    errors = [line for line in log.splitlines() if "ERROR" in line]
    lines += ["Errors from GetOrganelle:"] + (errors or ["    None recorded"])
    reads = re.findall(r"Reads used\s*=\s*(\d+)\+(\d+)", log)
    extended = organelle / "extended_1_paired.fq"
    if reads and extended.exists():
        r1_used = int(reads[-1][0])
        with extended.open() as stream:
            nlines = sum(1 for _ in stream)
        if nlines % 4 or r1_used <= 0:
            raise ValueError("Invalid GetOrganelle read counts")
        pairs_extended = nlines // 4
        lines += [
            f"Input R1 reads used by GetOrganelle: {r1_used}",
            f"Extended paired R1 reads: {pairs_extended}",
            f"Paired R1 recruitment: {100 * pairs_extended / r1_used:.5f}%",
            "Read recruitment is not a direct estimate of mitochondrial DNA abundance.",
        ]
    lines += ["", "### Sharkmer Results ###"]
    stats = Path("shark_results") / sample / f"{sample}.stats"
    if not stats.is_file() or not stats.stat().st_size:
        raise ValueError(f"Missing Sharkmer statistics: {stats}")
    products = sorted(stats.parent.glob("*.fasta"))
    if not products:
        lines.append("No PCR products assembled.")
    for product in products:
        lines.append(str(product))
        lines.extend(
            f"    {name}: {length} bp" for name, length in list(fasta_lengths(product))[:5]
        )
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_name(output.name + ".tmp")
    tmp.write_text("\n".join(lines) + "\n")
    tmp.replace(output)


if __name__ == "__main__":
    write_report(sys.argv[1], sys.argv[2], int(sys.argv[3]))
