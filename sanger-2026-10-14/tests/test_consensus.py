import importlib.util
import random
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

spec = importlib.util.spec_from_file_location(
    "process_sanger", Path(__file__).resolve().parents[1] / "scripts/process_sanger.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def record(seq, qualities):
    r = SeqRecord(Seq(seq))
    r.letter_annotations["phred_quality"] = qualities
    return r


def test_trim_rejects_bad_ends():
    assert module.trim_bounds([2] * 10 + [40] * 400 + [2] * 10) == (10, 410)


def test_consensus_preserves_high_quality_disagreement():
    rng = random.Random(1)
    s = "".join(rng.choices("ACGT", k=400))
    other = s[:200] + ("C" if s[200] != "C" else "A") + s[201:]
    seq, evidence, _ = module.consensus_pair(record(s, [40] * 400), record(other, [40] * 400))
    assert seq[200] not in "ACGTN"
    assert evidence[200]["decision"] == "high_quality_conflict_ambiguity"


def test_low_quality_insertion_requires_high_quality_opposing_flanks():
    rng = random.Random(2)
    s = "".join(rng.choices("ACGT", k=400))
    pos = 200
    while s[pos - 1] == "C" or s[pos] == "C":
        pos += 1
    inserted = s[:pos] + "C" + s[pos:]
    q = [40] * 401
    q[pos] = 10
    seq, evidence, _ = module.consensus_pair(record(inserted, q), record(s, [40] * 400))
    assert seq == s
    assert sum(r["decision"].startswith("drop_") for r in evidence) == 1
    with pytest.raises(ValueError, match="Unresolved indel"):
        module.consensus_pair(record(inserted, [40] * 401), record(s, [40] * 400))
