"""Known-answer checks for GenBank extraction and voucher identity."""

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from ncbi_references import extract_sequence, verify_identity


def test_extracts_one_longest_cds_with_strand_not_concatenation():
    r = SeqRecord(Seq("AACCGGTTAA"), id="TEST.1")
    r.features = [
        SeqFeature(FeatureLocation(0, 3), type="CDS", qualifiers={"gene": ["cox1"]}),
        SeqFeature(FeatureLocation(0, 6, strand=-1), type="CDS", qualifiers={"gene": ["COX1"]}),
        SeqFeature(
            FeatureLocation(0, 10), type="gene", qualifiers={"gene": ["cox1"], "pseudo": [""]}
        ),
    ]
    sequence, feature = extract_sequence(r, "cox1_CDS")
    assert str(sequence) == "CCGGTT"
    assert feature == "CDS:[0:6](-)"


def test_16s_product_without_gene_name_and_missing_feature():
    r = SeqRecord(Seq("AACCGGTT"), id="TEST.1")
    r.features = [
        SeqFeature(
            FeatureLocation(0, 4, strand=-1),
            type="rRNA",
            qualifiers={"product": ["16S ribosomal RNA"]},
        )
    ]
    assert str(extract_sequence(r, "16S_rRNA")[0]) == "GGTT"
    with pytest.raises(ValueError, match="Missing annotated"):
        extract_sequence(r, "cox1_CDS")


def test_voucher_normalization_does_not_accept_different_specimen():
    r = SeqRecord(Seq("ACGT"), id="TEST.1")
    r.annotations["organism"] = "Physalia physalis"
    r.features = [
        SeqFeature(
            FeatureLocation(0, 4), type="source", qualifiers={"specimen_voucher": ["YPM:IZ:110972"]}
        )
    ]
    row = dict(accession="TEST.1", species="Physalia physalis", sample="YPM-IZ-110972")
    verify_identity(r, row)
    with pytest.raises(ValueError, match="Voucher/isolate mismatch"):
        verify_identity(r, dict(row, sample="YPM-IZ-110973"))
