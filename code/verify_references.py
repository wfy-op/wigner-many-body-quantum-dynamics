"""Verify DOI metadata in the book bibliography against the Crossref API.

The script uses only the Python standard library.  It records network failures
instead of silently treating an unresolved DOI as valid.
"""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from difflib import SequenceMatcher
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIB_PATH = ROOT / "bibliography" / "references.bib"
OUTPUT_PATH = ROOT / "data" / "reference_verification.json"


def parse_entries(text: str) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    for raw_line in text.splitlines():
        line = raw_line.strip()
        match = re.match(r"@\w+\{([^,]+),", line)
        if match:
            current = {"key": match.group(1)}
            entries.append(current)
            continue
        if current is None:
            continue
        field = re.match(r"(title|doi)\s*=\s*\{(.*)\},?\s*$", line)
        if field:
            current[field.group(1)] = field.group(2).rstrip("},")
    return [entry for entry in entries if entry.get("doi")]


def normalize_title(title: str) -> str:
    title = re.sub(r"<[^>]+>", " ", title)
    title = re.sub(r"\\[a-zA-Z]+", " ", title)
    title = title.replace("--", " ")
    title = re.sub(r"[{}$\\]", "", title)
    title = re.sub(r"[^0-9a-zA-Z]+", " ", title).lower()
    return " ".join(title.split())


def verify(entry: dict[str, str]) -> dict[str, object]:
    doi = entry["doi"]
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "WignerManyBodyBook-reference-audit/1.0"},
    )
    result: dict[str, object] = {"key": entry["key"], "doi": doi}
    try:
        payload = None
        for attempt in range(4):
            try:
                with urllib.request.urlopen(request, timeout=20) as response:
                    payload = json.load(response)["message"]
                break
            except urllib.error.HTTPError as error:
                if error.code != 429 and error.code < 500:
                    raise
                if attempt == 3:
                    raise
                retry_after = error.headers.get("Retry-After")
                delay = float(retry_after) if retry_after else 2.0**attempt
                time.sleep(min(max(delay, 1.0), 8.0))
        if payload is None:
            raise ValueError("Crossref returned no metadata payload")
        crossref_title = payload.get("title", [""])[0]
        similarity = SequenceMatcher(
            None,
            normalize_title(entry.get("title", "")),
            normalize_title(crossref_title),
        ).ratio()
        result.update(
            {
                "status": "verified" if similarity >= 0.75 else "metadata_warning",
                "title_similarity": round(similarity, 4),
                "crossref_title": crossref_title,
                "publisher": payload.get("publisher"),
                "type": payload.get("type"),
            }
        )
    except (urllib.error.URLError, TimeoutError, KeyError, ValueError) as error:
        result.update({"status": "network_error", "error": str(error)})
    return result


def main() -> int:
    entries = parse_entries(BIB_PATH.read_text(encoding="utf-8"))
    # A small pool plus bounded retry avoids turning a transient Crossref 429
    # into an apparently unresolved bibliography entry.
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(verify, entries))
    counts: dict[str, int] = {}
    for result in results:
        status = str(result["status"])
        counts[status] = counts.get(status, 0) + 1
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": "Crossref REST API",
        "bibliography": str(BIB_PATH.relative_to(ROOT)).replace("\\", "/"),
        "entry_count_with_doi": len(entries),
        "status_counts": counts,
        "results": results,
    }
    OUTPUT_PATH.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report["status_counts"], ensure_ascii=False))
    # A network error is unresolved evidence, not a successful verification.
    # The command succeeds only when every DOI was independently resolved and
    # its title passed the declared similarity threshold.
    return 0 if counts == {"verified": len(entries)} else 1


if __name__ == "__main__":
    sys.exit(main())
