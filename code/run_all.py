"""Run every numerical benchmark used by the book and record a summary."""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "data" / "test_summary.json"
SCRIPTS = [
    "code/ch05_quadratic/harmonic_wigner.py",
    "code/ch06_kerr/kerr_exact_vs_twa.py",
    "code/ch08_dimer/dimer_exact_vs_twa.py",
    "code/ch11_sampling/gaussian_sampling.py",
    "code/ch12_propagation/split_step_benchmark.py",
    "code/ch13_observables/ordering_estimator_benchmark.py",
    "code/ch15_quench/quench_covariance.py",
    "code/ch16_phase/phase_diffusion.py",
    "code/ch17_order/stripe_scaling_toy.py",
]


def main() -> None:
    results = []
    failed = False
    for relative in SCRIPTS:
        started = time.perf_counter()
        completed = subprocess.run(
            [sys.executable, str(ROOT / relative)],
            cwd=ROOT,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
        )
        elapsed = time.perf_counter() - started
        passed = completed.returncode == 0
        failed = failed or not passed
        results.append(
            {
                "script": relative,
                "passed": passed,
                "returncode": completed.returncode,
                "elapsed_seconds": elapsed,
                "stdout_tail": completed.stdout.splitlines()[-12:],
                "stderr_tail": completed.stderr.splitlines()[-12:],
            }
        )
        print(f"[{'PASS' if passed else 'FAIL'}] {relative} ({elapsed:.2f} s)")

    summary = {
        "python_executable": sys.executable,
        "benchmark_count": len(SCRIPTS),
        "passed_count": sum(item["passed"] for item in results),
        "all_passed": not failed,
        "results": results,
    }
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
