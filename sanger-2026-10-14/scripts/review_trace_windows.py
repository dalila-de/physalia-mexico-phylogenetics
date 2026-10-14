"""Render representative original trace windows for visual QC review."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
windows = [
    ("CH887-D01_M13F_B03.ab1", 75, 110),
    ("CH887-D01_M13R_B03.ab1", 100, 135),
    ("CH887A-D01_M13F_D05.ab1", 95, 130),
    ("CH887A-D01_M13R_D05.ab1", 100, 135),
]
fig, axes = plt.subplots(4, 1, figsize=(14, 11), layout="constrained")
colors = {"A": "#009E73", "C": "#0072B2", "G": "#333333", "T": "#D55E00"}
for ax, (filename, start, end) in zip(axes, windows):
    r = SeqIO.read(ROOT / "raw" / filename, "abi")
    raw = r.annotations["abif_raw"]
    peaks = np.array(raw["PLOC2"])
    order = raw["FWO_1"].decode()
    left = int(peaks[start - 1]) - 8
    right = int(peaks[end - 1]) + 9
    for k, base in enumerate(order):
        signal = np.array(raw[f"DATA{9 + k}"])
        ax.plot(np.arange(left, right), signal[left:right], color=colors[base], lw=1)
    height = ax.get_ylim()[1]
    for i in range(start - 1, end):
        ax.text(
            peaks[i],
            height * 0.95,
            str(r.seq[i]),
            ha="center",
            fontsize=9,
            color=colors.get(str(r.seq[i]), "gray"),
        )
        ax.text(
            peaks[i],
            height * 0.83,
            str(r.letter_annotations["phred_quality"][i]),
            ha="center",
            fontsize=6,
        )
    ax.set(title=f"{filename}: base call and Phred quality", xlabel="Original base-call position")
    ax.set_xticks(peaks[start - 1 : end : 5], range(start, end + 1, 5))
fig.savefig(ROOT / "results/qc/trace_review.png", dpi=150)
