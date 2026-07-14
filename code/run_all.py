"""Run all numerical benchmarks with freshness and validation-contract gates."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Mapping

from benchmark_contract import (
    METRICS_SCHEMA_NAME,
    METRICS_SCHEMA_VERSION,
    RUN_ID_ENVIRONMENT_VARIABLE,
)


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "data" / "test_summary.json"
BENCHMARKS = [
    ("code/ch05_quadratic/harmonic_wigner.py", "data/ch05/harmonic_wigner_metrics.json"),
    ("code/ch06_kerr/kerr_exact_vs_twa.py", "data/ch06/kerr_metrics.json"),
    ("code/ch08_dimer/dimer_exact_vs_twa.py", "data/ch08/dimer_metrics.json"),
    ("code/ch11_sampling/gaussian_sampling.py", "data/ch11/squeezed_sampling_metrics.json"),
    ("code/ch12_propagation/split_step_benchmark.py", "data/ch12/split_step_metrics.json"),
    (
        "code/ch13_observables/ordering_estimator_benchmark.py",
        "data/ch13/ordering_estimator_metrics.json",
    ),
    (
        "code/ch14_17_capstone/two_component_field_quench.py",
        "data/ch14_17_capstone/two_component_field_quench_metrics.json",
    ),
    ("code/ch15_quench/quench_covariance.py", "data/ch15/quench_covariance_metrics.json"),
    ("code/ch16_phase/phase_diffusion.py", "data/ch16/phase_diffusion_metrics.json"),
    ("code/ch17_order/stripe_scaling_toy.py", "data/ch17/stripe_scaling_toy_metrics.json"),
]


def utc_isoformat(timestamp_ns: int) -> str:
    return datetime.fromtimestamp(timestamp_ns / 1_000_000_000, timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def file_state(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {"exists": False, "mtime_ns": None, "size_bytes": None, "sha256": None}
    stat = path.stat()
    return {
        "exists": True,
        "mtime_ns": stat.st_mtime_ns,
        "size_bytes": stat.st_size,
        "sha256": sha256_file(path),
    }


def package_version(distribution: str) -> str | None:
    try:
        return version(distribution)
    except PackageNotFoundError:
        return None


def git_text(root: Path, *arguments: str) -> str | None:
    completed = subprocess.run(
        ["git", *arguments],
        cwd=root,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if completed.returncode != 0:
        return None
    return completed.stdout.strip()


def environment_metadata(root: Path) -> dict[str, object]:
    status = git_text(root, "status", "--porcelain=v1")
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_executable": sys.executable,
        "python_version": sys.version,
        "platform": platform.platform(),
        "packages": {
            name: package_version(name) for name in ("numpy", "scipy", "matplotlib")
        },
        "git": {
            "commit": git_text(root, "rev-parse", "HEAD"),
            "branch": git_text(root, "branch", "--show-current"),
            "dirty": None if status is None else bool(status),
            "porcelain_status": None if status is None else status.splitlines(),
        },
    }


def read_metrics(path: Path) -> tuple[dict[str, object] | None, str | None]:
    if not path.is_file():
        return None, "expected metrics file was not produced"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        return None, f"could not parse metrics JSON: {error}"
    if not isinstance(payload, dict):
        return None, "metrics JSON root is not an object"
    return payload, None


def extract_check_results(
    node: object, path: str = "validation.checks"
) -> tuple[dict[str, bool], list[str]]:
    """Recompute simple checks and every nested ``passed`` result.

    Historical benchmarks store ``validation.checks`` as a mapping from names to
    booleans.  Richer benchmarks store an object for each check and put the result
    in its ``passed`` member; such objects may themselves contain lists or nested
    check groups.  Boolean metadata inside a rich check is deliberately ignored
    unless it is named ``passed`` or lives in an explicit ``checks``/``subchecks``
    container.
    """

    results: dict[str, bool] = {}
    errors: list[str] = []

    def walk(value: object, current_path: str, collect_plain_boolean: bool) -> None:
        if type(value) is bool:
            if collect_plain_boolean:
                results[current_path] = value
            return
        if isinstance(value, dict):
            has_passed = "passed" in value
            if has_passed:
                passed = value["passed"]
                if type(passed) is bool:
                    results[f"{current_path}.passed"] = passed
                else:
                    errors.append(f"{current_path}.passed is not boolean")
            for key, child in value.items():
                if key == "passed":
                    continue
                child_path = f"{current_path}.{key}"
                if isinstance(child, (dict, list)):
                    walk(
                        child,
                        child_path,
                        collect_plain_boolean=(
                            key in {"checks", "subchecks"}
                            or (collect_plain_boolean and not has_passed)
                        ),
                    )
                elif type(child) is bool:
                    if collect_plain_boolean and not has_passed:
                        results[child_path] = child
                elif collect_plain_boolean and not has_passed:
                    errors.append(f"{child_path} is not a supported check node")
            return
        if isinstance(value, list):
            for index, child in enumerate(value):
                walk(child, f"{current_path}[{index}]", collect_plain_boolean)
            return
        if collect_plain_boolean:
            errors.append(f"{current_path} is not a supported check node")

    walk(node, path, True)
    if not results:
        errors.append("validation.checks contains no boolean check results")
    return results, errors


def validate_metrics_contract(
    metrics: dict[str, object],
    *,
    expected_run_id: str,
    expected_script: str,
    expected_script_sha256: str | None,
    orchestrator_started_ns: int,
) -> dict[str, object]:
    schema = metrics.get("schema")
    schema_valid = (
        isinstance(schema, dict)
        and schema.get("name") == METRICS_SCHEMA_NAME
        and schema.get("version") == METRICS_SCHEMA_VERSION
    )

    run_errors: list[str] = []
    run = metrics.get("run")
    if not isinstance(run, dict):
        run_errors.append("run metadata is missing or not an object")
    else:
        if run.get("run_id") != expected_run_id:
            run_errors.append("run.run_id does not match this orchestrator invocation")
        if run.get("script") != expected_script:
            run_errors.append("run.script does not match the benchmark entry")
        if run.get("script_sha256") != expected_script_sha256:
            run_errors.append("run.script_sha256 does not match the executed script")
        if run.get("script_sha256_at_start") != expected_script_sha256:
            run_errors.append("run.script_sha256_at_start does not match the executed script")
        started_ns = run.get("started_at_unix_ns")
        finalized_ns = run.get("metrics_finalized_at_unix_ns")
        if type(started_ns) is not int:
            run_errors.append("run.started_at_unix_ns is not an integer")
        elif started_ns < orchestrator_started_ns:
            run_errors.append("child run started before the orchestrator wall-clock start")
        if type(finalized_ns) is not int:
            run_errors.append("run.metrics_finalized_at_unix_ns is not an integer")
        elif type(started_ns) is int and finalized_ns < started_ns:
            run_errors.append("metrics finalization precedes child start")

    validation = metrics.get("validation")
    validation_errors: list[str] = []
    check_results: dict[str, bool] = {}
    recomputed_all_passed: bool | None = None
    reported_all_passed: bool | None = None
    if not isinstance(validation, dict):
        validation_errors.append("validation is missing or not an object")
    else:
        reported = validation.get("all_passed")
        if type(reported) is bool:
            reported_all_passed = reported
        else:
            validation_errors.append("validation.all_passed is not boolean")
        if "checks" not in validation:
            validation_errors.append("validation.checks is missing")
        else:
            check_results, check_errors = extract_check_results(validation["checks"])
            validation_errors.extend(check_errors)
            if check_results:
                recomputed_all_passed = all(check_results.values())

    if (
        reported_all_passed is not None
        and recomputed_all_passed is not None
        and reported_all_passed != recomputed_all_passed
    ):
        validation_errors.append(
            "validation.all_passed disagrees with recursively recomputed checks"
        )
    validation_consistent = (
        not validation_errors
        and reported_all_passed is not None
        and recomputed_all_passed is not None
        and reported_all_passed == recomputed_all_passed
    )

    return {
        "metrics_schema": schema,
        "metrics_schema_valid": schema_valid,
        "run_metadata_valid": not run_errors,
        "run_metadata_errors": run_errors,
        "check_results": check_results,
        "check_count": len(check_results),
        "reported_validation_all_passed": reported_all_passed,
        "recomputed_validation_all_passed": recomputed_all_passed,
        "validation_consistent": validation_consistent,
        "validation_errors": validation_errors,
    }


def run_benchmark(
    *,
    root: Path,
    relative_script: str,
    relative_metrics: str,
    python_executable: str = sys.executable,
    extra_environment: Mapping[str, str] | None = None,
) -> dict[str, object]:
    root = root.resolve()
    script_path = root / relative_script
    metrics_path = root / relative_metrics
    script_state_before = file_state(script_path)
    old_metrics_state = file_state(metrics_path)
    started_wall_ns = time.time_ns()
    run_id = str(uuid.uuid4())
    environment = os.environ.copy()
    if extra_environment:
        environment.update(extra_environment)
    environment[RUN_ID_ENVIRONMENT_VARIABLE] = run_id

    started_perf = time.perf_counter()
    completed = subprocess.run(
        [python_executable, str(script_path)],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    elapsed = time.perf_counter() - started_perf
    finished_wall_ns = time.time_ns()
    new_metrics_state = file_state(metrics_path)
    script_state_after = file_state(script_path)

    freshness_errors: list[str] = []
    if not new_metrics_state["exists"]:
        freshness_errors.append("metrics file is missing after the child process")
    else:
        new_mtime = new_metrics_state["mtime_ns"]
        if type(new_mtime) is not int or new_mtime < started_wall_ns:
            freshness_errors.append("metrics mtime_ns predates this run's wall-clock start")
        if old_metrics_state["exists"]:
            old_mtime = old_metrics_state["mtime_ns"]
            if type(old_mtime) is int and type(new_mtime) is int and new_mtime <= old_mtime:
                freshness_errors.append("metrics mtime_ns did not advance beyond the old file")
            if new_metrics_state["sha256"] == old_metrics_state["sha256"]:
                freshness_errors.append("metrics SHA-256 is unchanged from the old file")
    fresh_metrics_written = not freshness_errors

    script_unchanged = (
        script_state_before["exists"]
        and script_state_after["exists"]
        and script_state_before["sha256"] == script_state_after["sha256"]
    )
    metrics, metrics_error = read_metrics(metrics_path)
    contract: dict[str, object] = {
        "metrics_schema": None,
        "metrics_schema_valid": False,
        "run_metadata_valid": False,
        "run_metadata_errors": ["metrics could not be read"],
        "check_results": {},
        "check_count": 0,
        "reported_validation_all_passed": None,
        "recomputed_validation_all_passed": None,
        "validation_consistent": False,
        "validation_errors": ["metrics could not be read"],
    }
    if metrics is not None:
        contract = validate_metrics_contract(
            metrics,
            expected_run_id=run_id,
            expected_script=relative_script,
            expected_script_sha256=script_state_before["sha256"],
            orchestrator_started_ns=started_wall_ns,
        )

    execution_passed = completed.returncode == 0
    structured_validation_passed = (
        contract["validation_consistent"] is True
        and contract["reported_validation_all_passed"] is True
    )
    validation_level = (
        "structured_thresholds"
        if (
            contract["metrics_schema_valid"] is True
            and contract["check_count"] > 0
            and contract["reported_validation_all_passed"] is not None
        )
        else "invalid_or_missing"
    )
    gate_passed = (
        execution_passed
        and fresh_metrics_written
        and metrics_error is None
        and script_unchanged
        and contract["metrics_schema_valid"] is True
        and contract["run_metadata_valid"] is True
        and structured_validation_passed
    )

    return {
        "script": relative_script,
        "script_sha256": script_state_before["sha256"],
        "script_sha256_after": script_state_after["sha256"],
        "script_unchanged_during_run": script_unchanged,
        "metrics_file": relative_metrics,
        "run_id": run_id,
        "started_wall_clock_unix_ns": started_wall_ns,
        "started_wall_clock_utc": utc_isoformat(started_wall_ns),
        "process_finished_wall_clock_unix_ns": finished_wall_ns,
        "process_finished_wall_clock_utc": utc_isoformat(finished_wall_ns),
        "old_metrics_state": old_metrics_state,
        "new_metrics_state": new_metrics_state,
        "fresh_metrics_written": fresh_metrics_written,
        "freshness_errors": freshness_errors,
        "execution_passed": execution_passed,
        "returncode": completed.returncode,
        "metrics_read_error": metrics_error,
        "metrics_sha256": new_metrics_state["sha256"],
        **contract,
        "validation_level": validation_level,
        "structured_validation_passed": structured_validation_passed,
        "gate_passed": gate_passed,
        "elapsed_seconds": elapsed,
        "stdout_tail": completed.stdout.splitlines()[-12:],
        "stderr_tail": completed.stderr.splitlines()[-12:],
    }


def main() -> None:
    # Capture repository provenance before any benchmark rewrites tracked metrics or
    # figures.  Recording it after the loop would make every clean-source run look
    # dirty merely because the run produced its evidence artifacts.
    environment_at_start = environment_metadata(ROOT)
    results = []
    for relative_script, relative_metrics in BENCHMARKS:
        result = run_benchmark(
            root=ROOT,
            relative_script=relative_script,
            relative_metrics=relative_metrics,
        )
        results.append(result)
        print(
            f"[{'PASS' if result['gate_passed'] else 'FAIL'}] {relative_script} "
            f"({result['elapsed_seconds']:.2f} s; fresh={result['fresh_metrics_written']}; "
            f"checks={result['check_count']})"
        )

    summary = {
        "schema_version": 4,
        "environment_at_start": environment_at_start,
        "benchmark_count": len(BENCHMARKS),
        "execution_passed_count": sum(item["execution_passed"] is True for item in results),
        "fresh_metrics_count": sum(item["fresh_metrics_written"] is True for item in results),
        "metrics_schema_valid_count": sum(
            item["metrics_schema_valid"] is True for item in results
        ),
        "validation_consistent_count": sum(
            item["validation_consistent"] is True for item in results
        ),
        "structured_validation_count": sum(
            item["validation_level"] == "structured_thresholds" for item in results
        ),
        "structured_validation_passed_count": sum(
            item["structured_validation_passed"] is True for item in results
        ),
        "unstructured_self_check_count": sum(
            item["validation_level"] != "structured_thresholds" for item in results
        ),
        "overall_gate_passed": all(item["gate_passed"] is True for item in results),
        "interpretation": (
            "Every benchmark must exit successfully, freshly rewrite its metrics after this "
            "orchestrator start, match the common metrics schema and run_id, preserve the "
            "executed script SHA-256, and report validation.all_passed consistent with a "
            "recursive recomputation of validation.checks."
        ),
        "results": results,
    }
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    if not summary["overall_gate_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
