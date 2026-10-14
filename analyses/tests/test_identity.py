"""Protect specimen identity, reverse complements, and longest-record selection."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from common import sequence_match, specimen
from prepare_tree_inputs import clean_tip


def test_sequence_containment_is_not_full_identity():
    assert sequence_match("ACGT", "TTACGTCC") == "query_contained_in_reference"
    assert sequence_match("AACG", "CGTT") == "identical_reverse_complement"
    assert sequence_match("AAAC", "AAAT") == "not_exact_or_contained"


def test_mito_sample_id_retains_specimen_number():
    assert specimen("_R_YPM-IZ-110972_mtgenome_scaffold_1") == "YPM-IZ-110972"


def test_rotation_is_a_distinct_diagnostic():
    assert sequence_match("ACGTA", "GTAAC", genome=True) == "rotation_equivalent"
    assert sequence_match("ACGTA", "GTAAC") == "not_exact_or_contained"


def test_upstream_orientation_tags_can_distinguish_assembly_variants():
    assert clean_tip("_R_YPM-IZ-111760_mito_path") != clean_tip("YPM-IZ-111760_mito_path")
