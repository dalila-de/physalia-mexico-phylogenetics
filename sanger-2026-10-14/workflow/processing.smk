"""Trace QC builds combined COI/16S panels for the single shared tree workflow."""
SANGER = ROOT / "sanger-2026-10-14"
SANGER_PYTHON = config.get("analysis_python", config.get("plot_python", "python"))
rule process_sanger:
    input:
        raw=sorted(str(p) for p in (SANGER/'raw').iterdir()),
        samples=str(SANGER/'config/samples.tsv'),
        localities=str(SANGER/'config/localities.tsv'),
        script=str(SANGER/'scripts/process_sanger.py'),
        references=expand(str(ROOT/'analyses/inputs/{gene}.fasta'),gene=['CO1','16S','18S','ITS']),
        tips=expand(str(ROOT/'analyses/inputs/{gene}.tips.tsv'),gene=['CO1','16S']),
        annotations=str(ROOT/'analyses/inputs/wgs_genetic_annotations.tsv')
    output:
        fasta=expand(str(ROOT/'analyses/inputs/combined/{gene}.fasta'),gene=['CO1','16S']),
        tips=expand(str(ROOT/'analyses/inputs/combined/{gene}.tips.tsv'),gene=['CO1','16S']),
        sequences=expand(str(SANGER/'inputs/{gene}.sanger.fasta'),gene=['CO1','16S']),
        annotations=str(ROOT/'analyses/inputs/genetic_annotations.tsv'),
        qc=str(SANGER/'results/qc/read_qc.json')
    shell: '{SANGER_PYTHON:q} {input.script:q}'
