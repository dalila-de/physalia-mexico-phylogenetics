"""Download the versioned NCBI accession manifest, or verify its bundled snapshot."""

import csv
import gzip
import hashlib
import io
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "references/ncbi"


def main() -> None:
    with (REFERENCE / "accessions.tsv").open() as stream:
        accessions = sorted({r["accession"] for r in csv.DictReader(stream, delimiter="\t")})
    destination = REFERENCE / "records.gb.gz"
    provenance = REFERENCE / "download.json"
    if destination.exists():
        meta = json.loads(provenance.read_text())
        if hashlib.sha256(destination.read_bytes()).hexdigest() != meta["sha256"]:
            raise ValueError("Bundled NCBI snapshot checksum differs")
        with gzip.open(destination, "rt") as stream:
            observed = [r.id for r in SeqIO.parse(stream, "genbank")]
        if sorted(observed) != accessions:
            raise ValueError("Snapshot membership differs from accession manifest")
        print(f"Verified {len(observed)} versioned NCBI records")
        return
    chunks, batches = [], []
    for start in range(0, len(accessions), 70):
        selected = accessions[start : start + 70]
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(
            dict(db="nuccore", id=",".join(selected), rettype="gbwithparts", retmode="text")
        )
        with urllib.request.urlopen(url, timeout=60) as response:
            payload = response.read()
        observed = [r.id for r in SeqIO.parse(io.StringIO(payload.decode()), "genbank")]
        if sorted(observed) != selected:
            raise ValueError("NCBI response is incomplete or accession versions differ")
        chunks.append(payload)
        batches.append(
            dict(url=url, accessions=selected, sha256=hashlib.sha256(payload).hexdigest())
        )
        print(f"Downloaded {len(selected)} records", flush=True)
        time.sleep(0.4)
    compressed = gzip.compress(b"\n".join(chunks), mtime=0)
    digest = hashlib.sha256(compressed).hexdigest()
    if provenance.exists() and json.loads(provenance.read_text())["sha256"] != digest:
        raise ValueError("Live metadata differs from pinned snapshot; review before replacing")
    tmp = destination.with_suffix(".tmp")
    tmp.write_bytes(compressed)
    tmp.replace(destination)
    provenance.write_text(
        json.dumps(
            dict(
                retrieved_utc=datetime.now(timezone.utc).isoformat(),
                count=len(accessions),
                sha256=digest,
                batches=batches,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
