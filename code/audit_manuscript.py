"""Static integrity checks for the LaTeX manuscript and its evidence files."""

from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "manuscript_audit.json"


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def main() -> None:
    tex_files = sorted(
        list((ROOT / "chapters").rglob("*.tex"))
        + list((ROOT / "appendices").glob("*.tex"))
        + list((ROOT / "frontmatter").glob("*.tex"))
        + list((ROOT / "config").glob("*.tex"))
        + [ROOT / "main.tex"]
    )
    chapter_files = sorted((ROOT / "chapters").rglob("ch*.tex"))
    appendix_files = sorted((ROOT / "appendices").glob("app*.tex"))
    texts = {path: read(path) for path in tex_files}

    labels = []
    citations = set()
    graphics = []
    includes = []
    placeholder_hits = []
    for path, content in texts.items():
        labels.extend(re.findall(r"\\label\{([^}]+)\}", content))
        for group in re.findall(r"\\cite\w*\{([^}]+)\}", content):
            citations.update(key.strip() for key in group.split(","))
        graphics.extend((path, item) for item in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", content))
        includes.extend((path, item) for item in re.findall(r"\\(?:input|include)\{([^}]+)\}", content))
        for pattern in ("本章初稿将在", "本附录将在", "TODO", "TBD"):
            if pattern in content:
                placeholder_hits.append({"file": str(path.relative_to(ROOT)), "pattern": pattern})

    label_counts = Counter(labels)
    duplicate_labels = sorted(key for key, count in label_counts.items() if count > 1)

    bib_text = read(ROOT / "bibliography" / "references.bib")
    bib_keys = set(re.findall(r"@\w+\{([^,]+),", bib_text))
    missing_bib_keys = sorted(citations - bib_keys)
    unused_bib_keys = sorted(bib_keys - citations)

    missing_graphics = []
    for source, item in graphics:
        candidate = ROOT / item
        options = [candidate] if candidate.suffix else [candidate.with_suffix(ext) for ext in (".pdf", ".png", ".jpg", ".jpeg")]
        if not any(option.exists() for option in options):
            missing_graphics.append({"source": str(source.relative_to(ROOT)), "target": item})

    missing_includes = []
    for source, item in includes:
        candidate = ROOT / item
        if not candidate.suffix:
            candidate = candidate.with_suffix(".tex")
        if not candidate.exists():
            missing_includes.append({"source": str(source.relative_to(ROOT)), "target": item})

    missing_reviews = [
        f"reviews/ch{number:02d}_review.md"
        for number in range(1, 19)
        if not (ROOT / f"reviews/ch{number:02d}_review.md").exists()
    ]
    test_summary = json.loads(read(ROOT / "data" / "test_summary.json"))

    audit = {
        "chapter_file_count": len(chapter_files),
        "appendix_file_count": len(appendix_files),
        "review_file_count": 18 - len(missing_reviews),
        "tex_file_count": len(tex_files),
        "label_count": len(labels),
        "duplicate_labels": duplicate_labels,
        "citation_key_count": len(citations),
        "bibliography_entry_count": len(bib_keys),
        "missing_bibliography_keys": missing_bib_keys,
        "unused_bibliography_keys": unused_bib_keys,
        "figure_reference_count": len(graphics),
        "missing_graphics": missing_graphics,
        "missing_includes": missing_includes,
        "placeholder_hits": placeholder_hits,
        "missing_chapter_reviews": missing_reviews,
        "numerical_benchmarks_all_passed": test_summary.get("all_passed", False),
        "numerical_benchmarks_passed": test_summary.get("passed_count", 0),
        "passed": not any(
            (
                len(chapter_files) != 18,
                len(appendix_files) != 4,
                duplicate_labels,
                missing_bib_keys,
                missing_graphics,
                missing_includes,
                placeholder_hits,
                missing_reviews,
                not test_summary.get("all_passed", False),
            )
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    if not audit["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
