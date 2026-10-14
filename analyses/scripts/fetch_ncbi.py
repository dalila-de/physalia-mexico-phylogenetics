import argparse
import hashlib
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--out", type=Path, required=True)
p = parser.parse_args().out
p.mkdir(parents=True, exist_ok=True)
base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def get(name, args):
    url = base + name + ".fcgi?" + urllib.parse.urlencode(args)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                data = response.read()
            time.sleep(0.4)
            return data, url
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 + attempt)


query = "Physalia[Organism] AND 100:100000[SLEN]"
data, url = get("esearch", {"db": "nuccore", "term": query, "retmode": "json", "retmax": 10000})
(p / "esearch.json").write_bytes(data)
d = json.loads(data)["esearchresult"]
ids = d["idlist"]
assert len(ids) == int(d["count"])
records = []
for i in range(0, len(ids), 100):
    data, url = get(
        "efetch",
        {
            "db": "nuccore",
            "id": ",".join(ids[i : i + 100]),
            "rettype": "gbwithparts",
            "retmode": "text",
        },
    )
    if not data.startswith(b"LOCUS"):
        raise ValueError(data[:500])
    out = p / f"records_{i:04d}.gb"
    out.write_bytes(data)
    records.append(
        {
            "file": out.name,
            "url": url,
            "sha256": hashlib.sha256(data).hexdigest(),
            "ids": ids[i : i + 100],
        }
    )
    print(out.name, len(data), flush=True)
(p / "download.json").write_text(
    json.dumps(
        {
            "utc": datetime.now(timezone.utc).isoformat(),
            "query": query,
            "count": len(ids),
            "batches": records,
        },
        indent=2,
    )
    + "\n"
)
