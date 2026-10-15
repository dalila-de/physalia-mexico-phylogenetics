"""Quality-trim Sanger traces, build conservative marker sequences, and retain evidence."""

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Align, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
ALIGNER = Align.PairwiseAligner(
    mode="local", match_score=2, mismatch_score=-3, open_gap_score=-6, extend_gap_score=-1
)
IUPAC = {
    frozenset(k): v
    for k, v in [("AG", "R"), ("CT", "Y"), ("AC", "M"), ("GT", "K"), ("AT", "W"), ("CG", "S")]
}
HCO = "TAAACTTCAGGGTGACCAAAAAATCA"


def trim_bounds(qualities):
    """Maximum-scoring contiguous interval, score=0.01-error_probability (Q20)."""
    current = best = 0.0
    start = lo = hi = 0
    for i, q in enumerate(qualities):
        if current <= 0:
            current = 0.0
            start = i
        current += 0.01 - 10 ** (-q / 10)
        if current > best:
            best, lo, hi = current, start, i + 1
    return lo, hi


def consensus_pair(fwd, rev):
    """Use only bidirectional overlap; low-quality insertion artifacts may be dropped."""
    aln = ALIGNER.align(fwd.seq, rev.seq)[0]
    if aln.counts().identities < 300:
        raise ValueError("Insufficient forward/reverse agreement")
    out, evidence = [], []
    indices = aln.indices
    for col, (i, j) in enumerate(indices.T):
        b1, b2 = (str(fwd.seq[i]) if i >= 0 else "-"), (str(rev.seq[j]) if j >= 0 else "-")
        q1 = fwd.letter_annotations["phred_quality"][i] if i >= 0 else 0
        q2 = rev.letter_annotations["phred_quality"][j] if j >= 0 else 0
        if i < 0 or j < 0:
            gap_row = 0 if i < 0 else 1
            before = indices[gap_row, :col]
            after = indices[gap_row, col + 1 :]
            before = before[before >= 0]
            after = after[after >= 0]
            other = fwd if gap_row == 0 else rev
            flank_q = (
                min(
                    other.letter_annotations["phred_quality"][before[-1]],
                    other.letter_annotations["phred_quality"][after[0]],
                )
                if len(before) and len(after)
                else 0
            )
            if max(q1, q2) >= 20 or flank_q < 20:
                raise ValueError(
                    f"Unresolved indel at aligned column {col}; manual review required"
                )
            call, decision = "", "drop_low_quality_insertion_opposed_by_Q20_flanks"
        elif b1 == b2 and max(q1, q2) >= 20 and b1 in "ACGT":
            call, decision = b1, "agreement"
        elif q1 >= 20 and q2 >= 20 and b1 in "ACGT" and b2 in "ACGT":
            call, decision = IUPAC[frozenset([b1, b2])], "high_quality_conflict_ambiguity"
        elif q1 >= 20 and b1 in "ACGT":
            call, decision = b1, "forward_only_quality"
        elif q2 >= 20 and b2 in "ACGT":
            call, decision = b2, "reverse_only_quality"
        else:
            call, decision = "N", "mask_low_quality"
        out.append(call)
        evidence.append(
            dict(
                column=col + 1,
                forward_index=i,
                reverse_complement_index=j,
                forward_base=b1,
                reverse_base=b2,
                forward_quality=q1,
                reverse_quality=q2,
                consensus=call or "-",
                decision=decision,
            )
        )
    return "".join(out), evidence, str(aln)


def write_table(path, rows):
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    (REPO / "analyses/inputs/combined").mkdir(parents=True, exist_ok=True)
    out = ROOT / "results/qc"
    out.mkdir(parents=True, exist_ok=True)
    with (ROOT / "config/samples.tsv").open() as f:
        samples = list(csv.DictReader(f, delimiter="\t"))
    reports = []
    manifest = []
    overview, axes = plt.subplots(4, 1, figsize=(12, 10), layout="constrained")
    for sample_index, sample in enumerate(samples):
        traces = []
        for direction_index, direction in enumerate(["forward", "reverse"]):
            path = ROOT / "raw" / sample[direction]
            r = SeqIO.read(path, "abi")
            text = "".join(path.with_suffix(".seq").read_text().split()).upper()
            if text != str(r.seq):
                raise ValueError("SEQ and AB1 base calls differ")
            q = r.letter_annotations["phred_quality"]
            lo, hi = trim_bounds(q)
            primer = None
            motif = str(Seq(HCO).reverse_complement())
            pos = str(r.seq).find(motif)
            if sample["gene"] == "CO1" and pos >= 0 and pos > len(r) * 0.7:
                hi = min(hi, pos)
                primer = dict(
                    name="HCO2198_reverse_complement",
                    raw_start_1based=pos + 1,
                    note="Exact terminal motif; experimental primer identity unconfirmed",
                )
            t = r[lo:hi]
            traces.append(t)
            stats = dict(
                sample=sample["sample"],
                gene=sample["gene"],
                file=path.name,
                direction=direction,
                raw_length=len(r),
                raw_median_Q=float(np.median(q)),
                trim_start_1based=lo + 1,
                trim_end_1based_inclusive=hi,
                retained_length=len(t),
                retained_below_Q20=sum(x < 20 for x in t.letter_annotations["phred_quality"]),
                primer_clip=primer,
            )
            reports.append(stats)
            SeqIO.write(t, out / (path.stem + ".trimmed.fastq"), "fastq")
            ax = axes[sample_index * 2 + direction_index]
            ax.plot(np.arange(1, len(r) + 1), q, lw=0.65)
            ax.axhline(20, color="gray", ls="--", lw=0.6)
            ax.axvspan(lo + 1, hi, color="#009E73", alpha=0.2)
            ax.set(
                title=f"{path.stem}: retained {len(t)} bases",
                ylabel="Phred Q",
                xlabel="Original base-call position",
            )
            raw = r.annotations["abif_raw"]
            order = raw["FWO_1"].decode()
            peaks = np.array(raw["PLOC2"])
            colors = {"A": "#009E73", "C": "#0072B2", "G": "#333333", "T": "#D55E00"}
            with PdfPages(out / (path.stem + ".chromatogram.pdf")) as pdf:
                for page in range(0, len(r), 400):
                    fig, axs = plt.subplots(4, 1, figsize=(14, 10), layout="constrained")
                    for panel, axis in enumerate(axs):
                        begin = page + panel * 100
                        end = min(begin + 100, len(r))
                        if begin >= end:
                            axis.axis("off")
                            continue
                        left = max(0, int(peaks[begin]) - 8)
                        right = int(peaks[end - 1]) + 9
                        for channel, base in enumerate(order):
                            signal = np.array(raw[f"DATA{9 + channel}"])
                            axis.plot(
                                np.arange(left, min(right, len(signal))),
                                signal[left:right],
                                color=colors[base],
                                lw=0.7,
                            )
                        height = axis.get_ylim()[1]
                        for k in range(begin, end):
                            axis.text(
                                peaks[k],
                                height * 0.92,
                                str(r.seq[k]),
                                fontsize=5,
                                ha="center",
                                color=colors.get(str(r.seq[k]), "gray"),
                            )
                        axis.set_title(f"{path.stem}: bases {begin + 1}–{end}", fontsize=9)
                        axis.set_xticks(
                            peaks[begin:end:10], np.arange(begin + 1, end + 1, 10), fontsize=7
                        )
                    pdf.savefig(fig)
                    plt.close(fig)
            for source in (path, path.with_suffix(".seq")):
                manifest.append(
                    dict(
                        file=str(source.relative_to(ROOT)),
                        sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    )
                )
        fwd, reverse = traces
        rev = reverse.reverse_complement()
        if len(fwd) < 100:
            raise ValueError("Forward read unusable")
        if len(rev) < 100:
            sequence = "".join(
                b if q >= 20 and b in "ACGT" else "N"
                for b, q in zip(str(fwd.seq), fwd.letter_annotations["phred_quality"])
            )
            evidence = [
                dict(position=i + 1, base=b, quality=q, consensus=sequence[i])
                for i, (b, q) in enumerate(
                    zip(str(fwd.seq), fwd.letter_annotations["phred_quality"])
                )
            ]
            basis = "forward_only; reverse_failed_quality"
        else:
            sequence, evidence, alignment = consensus_pair(fwd, rev)
            (out / f"{sample['gene']}.forward_reverse.txt").write_text(alignment + "\n")
            basis = "bidirectional_overlap_only"
        write_table(out / f"{sample['gene']}.consensus_evidence.tsv", evidence)
        rec = SeqRecord(Seq(sequence), id="sanger__" + sample["sample"], description="")
        SeqIO.write(rec, ROOT / "inputs" / f"{sample['gene']}.sanger.fasta", "fasta")
        refs = list(SeqIO.parse(REPO / "analyses/inputs" / f"{sample['gene']}.fasta", "fasta"))
        if rec.id in {x.id for x in refs}:
            raise ValueError("Duplicate tip")
        SeqIO.write(
            refs + [rec], REPO / "analyses/inputs/combined" / f"{sample['gene']}.fasta", "fasta"
        )
        with (REPO / "analyses/inputs" / f"{sample['gene']}.tips.tsv").open() as stream:
            tips = list(csv.DictReader(stream, delimiter="\t"))
        tips.append(
            dict(
                tip=rec.id,
                sample=sample["sample"],
                source="Sanger_AB1",
                species="Physalia physalis",
                original_header=sample["amplicon"],
                sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
            )
        )
        write_table(REPO / "analyses/inputs/combined" / f"{sample['gene']}.tips.tsv", tips)
        comparisons = []
        for ref in refs:
            a = ALIGNER.align(ref.seq, rec.seq)[0]
            c = a.counts()
            comparisons.append(
                dict(
                    reference=ref.id,
                    score=a.score,
                    identities=c.identities,
                    mismatches=c.mismatches,
                    gap_bases=c.gaps,
                    reference_start=int(a.coordinates[0, 0]) + 1,
                    reference_end=int(a.coordinates[0, -1]),
                    query_start=int(a.coordinates[1, 0]) + 1,
                    query_end=int(a.coordinates[1, -1]),
                )
            )
        comparisons.sort(key=lambda x: -x["score"])
        write_table(out / f"{sample['gene']}.reference_matches.tsv", comparisons)
        marker_screen = []
        for marker in ["CO1", "16S", "18S", "ITS"]:
            panel = [
                r
                for r in SeqIO.parse(REPO / "analyses/inputs" / f"{marker}.fasta", "fasta")
                if r.id.startswith("ref__YPM")
            ][:3]
            for ref in panel:
                for reverse_orientation in [False, True]:
                    query = rec.seq.reverse_complement() if reverse_orientation else rec.seq
                    hit = ALIGNER.align(ref.seq, query)[0]
                    c = hit.counts()
                    marker_screen.append(
                        dict(
                            marker=marker,
                            reference=ref.id,
                            reverse_complement=reverse_orientation,
                            score=hit.score,
                            identities=c.identities,
                            mismatches=c.mismatches,
                            gap_bases=c.gaps,
                        )
                    )
        marker_screen.sort(key=lambda r: -r["score"])
        write_table(out / f"{sample['gene']}.marker_screen.tsv", marker_screen)
        if (
            marker_screen[0]["marker"] != sample["gene"]
            or marker_screen[0]["identities"] < len(sequence) * 0.85
        ):
            raise ValueError("Marker identity screen failed")
        coding = []
        if sample["gene"] == "CO1":
            for frame in range(3):
                coding_seq = rec.seq[frame : len(rec) - (len(rec) - frame) % 3]
                protein = str(coding_seq.translate(table=4))
                coding.append(
                    dict(
                        frame_offset=frame,
                        stops=protein.count("*"),
                        ambiguous_amino_acids=protein.count("X"),
                    )
                )
        (out / f"{sample['gene']}.consensus.json").write_text(
            json.dumps(
                dict(
                    sample=sample["sample"],
                    gene=sample["gene"],
                    length=len(sequence),
                    ambiguous_bases=sum(b not in "ACGT" for b in sequence),
                    basis=basis,
                    translation_table4=coding,
                ),
                indent=2,
            )
            + "\n"
        )
    with (REPO / "analyses/inputs/wgs_genetic_annotations.tsv").open() as stream:
        annotations = list(csv.DictReader(stream, delimiter="\t"))
    a = dict.fromkeys(annotations[0], "")
    a.update(
        sample="CH887",
        location=samples[0]["location"],
        plot_group="Sanger",
        basis="User identifies both amplicons as one Caribbean specimen; no genome-wide E/W assignment",
        genetic_group="unassigned",
        vagrant="False",
    )
    annotations.append(a)
    with (ROOT / "config/localities.tsv").open() as stream:
        localities = {r["sample"]: r["location"] for r in csv.DictReader(stream, delimiter="\t")}
    for row in annotations:
        if row["sample"] in localities:
            row["location"] = localities[row["sample"]]
    write_table(REPO / "analyses/inputs/genetic_annotations.tsv", annotations)
    overview.savefig(out / "quality_overview.png", dpi=160)
    plt.close(overview)
    (out / "read_qc.json").write_text(json.dumps(reports, indent=2) + "\n")
    write_table(ROOT / "provenance/raw_checksums.tsv", manifest)


if __name__ == "__main__":
    main()
