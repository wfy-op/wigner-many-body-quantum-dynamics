"""Static integrity checks for the LaTeX manuscript and its evidence files."""

from __future__ import annotations

import hashlib
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
    control_character_hits = []
    malformed_formula_hits = []
    for path, content in texts.items():
        labels.extend(re.findall(r"\\label\{([^}]+)\}", content))
        for group in re.findall(r"\\cite\w*\{([^}]+)\}", content):
            citations.update(key.strip() for key in group.split(","))
        graphics.extend((path, item) for item in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", content))
        includes.extend((path, item) for item in re.findall(r"\\(?:input|include)\{([^}]+)\}", content))
        for pattern in ("本章初稿将在", "本附录将在", "TODO", "TBD"):
            if pattern in content:
                placeholder_hits.append({"file": str(path.relative_to(ROOT)), "pattern": pattern})
        for match in re.finditer(r"[\x00-\x08\x0b\x0c\x0e-\x1f]", content):
            control_character_hits.append(
                {
                    "file": str(path.relative_to(ROOT)),
                    "codepoint": f"U+{ord(match.group()):04X}",
                }
            )
        for pattern in (r"(?<!\\)int\\dd", r"\^Lrac", r"\\frac12\s*rac"):
            if re.search(pattern, content):
                malformed_formula_hits.append(
                    {"file": str(path.relative_to(ROOT)), "pattern": pattern}
                )

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
    reference_verification = json.loads(
        read(ROOT / "data" / "reference_verification.json")
    )
    chapter_citation_counts = {
        str(path.relative_to(ROOT)): len(
            re.findall(r"\\cite\w*\{[^}]+\}", texts[path])
        )
        for path in chapter_files
    }
    chapters_without_citations = sorted(
        path for path, count in chapter_citation_counts.items() if count == 0
    )
    chapter_exercise_block_counts = {
        str(path.relative_to(ROOT)): texts[path].count(r"\begin{exercises}")
        for path in chapter_files
    }
    chapters_without_one_exercise_block = sorted(
        path
        for path, count in chapter_exercise_block_counts.items()
        if count != 1
    )
    benchmark_count = test_summary.get("benchmark_count", 0)
    execution_passed_count = test_summary.get("execution_passed_count", 0)
    structured_validation_count = test_summary.get(
        "structured_validation_count", 0
    )
    structured_validation_passed_count = test_summary.get(
        "structured_validation_passed_count", 0
    )
    unstructured_self_check_count = test_summary.get(
        "unstructured_self_check_count", 0
    )
    benchmark_evidence_issues = []
    for result in test_summary.get("results", []):
        script = ROOT / result.get("script", "")
        metrics = ROOT / result.get("metrics_file", "")
        if not script.is_file():
            benchmark_evidence_issues.append(
                {"script": result.get("script"), "issue": "missing script"}
            )
        if not metrics.is_file():
            benchmark_evidence_issues.append(
                {"metrics": result.get("metrics_file"), "issue": "missing metrics"}
            )
            continue
        current_hash = hashlib.sha256(metrics.read_bytes()).hexdigest()
        if current_hash != result.get("metrics_sha256"):
            benchmark_evidence_issues.append(
                {
                    "metrics": result.get("metrics_file"),
                    "issue": "SHA-256 differs from data/test_summary.json",
                    "recorded": result.get("metrics_sha256"),
                    "current": current_hash,
                }
            )
        if not (
            result.get("gate_passed") is True
            and result.get("validation_level") == "structured_thresholds"
            and result.get("structured_validation_passed") is True
        ):
            benchmark_evidence_issues.append(
                {
                    "script": result.get("script"),
                    "issue": "benchmark is not a passed structured validation",
                }
            )
    references_verified = reference_verification.get("status_counts", {}).get(
        "verified", 0
    )
    references_with_doi = reference_verification.get("entry_count_with_doi", 0)

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
        "control_character_hits": control_character_hits,
        "malformed_formula_hits": malformed_formula_hits,
        "missing_chapter_reviews": missing_reviews,
        "chapter_citation_counts": chapter_citation_counts,
        "chapters_without_citations": chapters_without_citations,
        "chapter_exercise_block_counts": chapter_exercise_block_counts,
        "chapters_without_one_exercise_block": chapters_without_one_exercise_block,
        "reference_verification": {
            "entries_with_doi": references_with_doi,
            "verified": references_verified,
            "all_verified": references_with_doi == len(bib_keys) == references_verified,
        },
        "numerical_benchmarks": {
            "benchmark_count": benchmark_count,
            "execution_passed_count": execution_passed_count,
            "structured_validation_count": structured_validation_count,
            "structured_validation_passed_count": structured_validation_passed_count,
            "unstructured_self_check_count": unstructured_self_check_count,
            "overall_gate_passed": test_summary.get("overall_gate_passed", False),
            "evidence_issues": benchmark_evidence_issues,
        },
        "second_pass_review_exists": (
            ROOT / "reviews" / "revision_second_pass.md"
        ).is_file(),
        "passed": not any(
            (
                len(chapter_files) != 18,
                len(appendix_files) != 6,
                duplicate_labels,
                missing_bib_keys,
                missing_graphics,
                missing_includes,
                placeholder_hits,
                control_character_hits,
                malformed_formula_hits,
                missing_reviews,
                chapters_without_citations,
                chapters_without_one_exercise_block,
                references_with_doi != len(bib_keys),
                references_verified != len(bib_keys),
                benchmark_count != 10,
                execution_passed_count != benchmark_count,
                structured_validation_count != benchmark_count,
                structured_validation_passed_count != benchmark_count,
                unstructured_self_check_count != 0,
                benchmark_evidence_issues,
                not test_summary.get("overall_gate_passed", False),
                not (ROOT / "reviews" / "revision_second_pass.md").is_file(),
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
