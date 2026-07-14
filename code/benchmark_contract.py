"""Shared, machine-readable provenance contract for numerical benchmark metrics."""

from __future__ import annotations

import hashlib
import os
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


METRICS_SCHEMA_NAME = "wigner-manybody-benchmark-metrics"
METRICS_SCHEMA_VERSION = 1
RUN_ID_ENVIRONMENT_VARIABLE = "BOOK_BENCHMARK_RUN_ID"


def _utc_isoformat(timestamp_ns: int) -> str:
    return datetime.fromtimestamp(timestamp_ns / 1_000_000_000, timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def begin_metrics_run(script_file: str | Path) -> dict[str, Any]:
    """Capture provenance before any scientific computation starts."""

    script_path = Path(script_file).resolve()
    started_ns = time.time_ns()
    return {
        "run_id": os.environ.get(
            RUN_ID_ENVIRONMENT_VARIABLE, f"standalone-{uuid.uuid4()}"
        ),
        "started_at_unix_ns": started_ns,
        "started_at_utc": _utc_isoformat(started_ns),
        "script_path": script_path,
        "script_sha256_at_start": _sha256(script_path),
    }


def finalize_metrics(
    metrics: dict[str, Any], run_context: dict[str, Any], script_file: str | Path
) -> dict[str, Any]:
    """Attach the common schema and run record without changing scientific fields."""

    script_path = Path(script_file).resolve()
    project_root = script_path.parents[2]
    relative_script = script_path.relative_to(project_root).as_posix()
    finalized_ns = time.time_ns()
    payload = dict(metrics)
    payload["schema"] = {
        "name": METRICS_SCHEMA_NAME,
        "version": METRICS_SCHEMA_VERSION,
    }
    payload["run"] = {
        "run_id": run_context["run_id"],
        "script": relative_script,
        "script_sha256": _sha256(script_path),
        "script_sha256_at_start": run_context["script_sha256_at_start"],
        "started_at_unix_ns": run_context["started_at_unix_ns"],
        "started_at_utc": run_context["started_at_utc"],
        "metrics_finalized_at_unix_ns": finalized_ns,
        "metrics_finalized_at_utc": _utc_isoformat(finalized_ns),
        "python_executable": sys.executable,
    }
    return payload
