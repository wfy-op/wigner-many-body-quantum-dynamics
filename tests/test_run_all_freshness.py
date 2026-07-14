from __future__ import annotations

import json
import os
import sys
import tempfile
import textwrap
import time
import unittest
from pathlib import Path
from unittest.mock import patch


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CODE_ROOT = PROJECT_ROOT / "code"
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

from benchmark_contract import (
    METRICS_SCHEMA_NAME,
    METRICS_SCHEMA_VERSION,
    RUN_ID_ENVIRONMENT_VARIABLE,
    begin_metrics_run,
    finalize_metrics,
)
from run_all import extract_check_results, run_benchmark, validate_metrics_contract


class ValidationContractTests(unittest.TestCase):
    def test_shared_helper_attaches_schema_and_script_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            script = Path(directory) / "code" / "chapter" / "example.py"
            script.parent.mkdir(parents=True)
            script.write_text("# fixture\n", encoding="utf-8")
            with patch.dict(
                os.environ, {RUN_ID_ENVIRONMENT_VARIABLE: "unit-test-run"}
            ):
                context = begin_metrics_run(script)
                payload = finalize_metrics({"answer": 42}, context, script)

            self.assertEqual(
                payload["schema"],
                {"name": METRICS_SCHEMA_NAME, "version": METRICS_SCHEMA_VERSION},
            )
            self.assertEqual(payload["run"]["run_id"], "unit-test-run")
            self.assertEqual(payload["run"]["script"], "code/chapter/example.py")
            self.assertEqual(
                payload["run"]["script_sha256"],
                payload["run"]["script_sha256_at_start"],
            )

    def test_recursive_check_recomputation_supports_both_schemas(self) -> None:
        checks = {
            "plain_boolean": True,
            "rich_check": {
                "value": 0.1,
                "maximum": 0.2,
                "diagnostic_flag": True,
                "per_case": [{"passed": True}, {"passed": False}],
                "passed": False,
            },
            "nested_group": {"checks": {"child": {"passed": True}}},
        }
        results, errors = extract_check_results(checks)
        self.assertEqual(errors, [])
        self.assertEqual(len(results), 5)
        self.assertNotIn("validation.checks.rich_check.diagnostic_flag", results)
        self.assertIn("validation.checks.rich_check.per_case[0].passed", results)
        self.assertFalse(all(results.values()))

    def test_reported_all_passed_must_match_recursive_result(self) -> None:
        metrics = {
            "schema": {"name": METRICS_SCHEMA_NAME, "version": METRICS_SCHEMA_VERSION},
            "run": {
                "run_id": "expected",
                "script": "code/example.py",
                "script_sha256": "abc",
                "script_sha256_at_start": "abc",
                "started_at_unix_ns": 20,
                "metrics_finalized_at_unix_ns": 30,
            },
            "validation": {
                "checks": {"first": True, "second": {"passed": False}},
                "all_passed": True,
            },
        }
        result = validate_metrics_contract(
            metrics,
            expected_run_id="expected",
            expected_script="code/example.py",
            expected_script_sha256="abc",
            orchestrator_started_ns=10,
        )
        self.assertFalse(result["validation_consistent"])
        self.assertFalse(result["recomputed_validation_all_passed"])
        self.assertIn("disagrees", " ".join(result["validation_errors"]))


class FreshMetricsRegressionTests(unittest.TestCase):
    def test_zero_exit_without_fresh_metrics_cannot_reuse_old_file(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "code").mkdir()
            (root / "data").mkdir()
            script = root / "code" / "does_not_write.py"
            script.write_text("# Intentionally exits zero without writing metrics.\n", encoding="utf-8")
            metrics_path = root / "data" / "metrics.json"
            metrics_path.write_text(
                json.dumps(
                    {
                        "schema": {
                            "name": METRICS_SCHEMA_NAME,
                            "version": METRICS_SCHEMA_VERSION,
                        },
                        "run": {"run_id": "old"},
                        "validation": {"checks": {"old": True}, "all_passed": True},
                    }
                ),
                encoding="utf-8",
            )
            old_ns = time.time_ns() - 2_000_000_000
            os.utime(metrics_path, ns=(old_ns, old_ns))

            result = run_benchmark(
                root=root,
                relative_script="code/does_not_write.py",
                relative_metrics="data/metrics.json",
                python_executable=sys.executable,
            )

            self.assertTrue(result["execution_passed"])
            self.assertFalse(result["fresh_metrics_written"])
            self.assertFalse(result["gate_passed"])
            self.assertEqual(
                result["old_metrics_state"]["sha256"],
                result["new_metrics_state"]["sha256"],
            )
            self.assertIn("unchanged", " ".join(result["freshness_errors"]))

    def test_fresh_contract_metrics_pass_the_gate(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "code").mkdir()
            (root / "data").mkdir()
            script = root / "code" / "writes_metrics.py"
            script.write_text(
                textwrap.dedent(
                    """
                    import hashlib
                    import json
                    import os
                    import time
                    from pathlib import Path

                    script = Path(__file__).resolve()
                    started = time.time_ns()
                    digest = hashlib.sha256(script.read_bytes()).hexdigest()
                    payload = {
                        "schema": {
                            "name": "wigner-manybody-benchmark-metrics",
                            "version": 1,
                        },
                        "run": {
                            "run_id": os.environ["BOOK_BENCHMARK_RUN_ID"],
                            "script": "code/writes_metrics.py",
                            "script_sha256": digest,
                            "script_sha256_at_start": digest,
                            "started_at_unix_ns": started,
                            "metrics_finalized_at_unix_ns": time.time_ns(),
                        },
                        "validation": {
                            "checks": {
                                "plain": True,
                                "rich": {"value": 1.0, "passed": True},
                            },
                            "all_passed": True,
                        },
                    }
                    target = script.parents[1] / "data" / "metrics.json"
                    target.write_text(json.dumps(payload), encoding="utf-8")
                    """
                ).lstrip(),
                encoding="utf-8",
            )

            result = run_benchmark(
                root=root,
                relative_script="code/writes_metrics.py",
                relative_metrics="data/metrics.json",
                python_executable=sys.executable,
            )

            self.assertTrue(result["fresh_metrics_written"])
            self.assertTrue(result["validation_consistent"])
            self.assertTrue(result["run_metadata_valid"])
            self.assertTrue(result["gate_passed"])


if __name__ == "__main__":
    unittest.main()
