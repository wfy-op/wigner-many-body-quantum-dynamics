"""Monte Carlo check of single-mode squeezed-vacuum Wigner sampling."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


CODE_ROOT = Path(__file__).resolve().parents[1]
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))
from benchmark_contract import begin_metrics_run, finalize_metrics


ROOT = Path(__file__).resolve().parents[2]
FIGURE = ROOT / "figures" / "ch11" / "squeezed_wigner_sampling.pdf"
METRICS = ROOT / "data" / "ch11" / "squeezed_sampling_metrics.json"


def main() -> None:
    run_context = begin_metrics_run(__file__)
    seed = 20260714
    trajectories = 100_000
    squeeze_r = 0.8
    rng = np.random.default_rng(seed)

    # Vacuum Wigner samples: E|beta|^2 = 1/2.
    beta = (rng.standard_normal(trajectories) + 1j * rng.standard_normal(trajectories)) / 2
    # Book convention: r > 0 squeezes Q and gives E[alpha**2] < 0.
    alpha = np.cosh(squeeze_r) * beta - np.sinh(squeeze_r) * beta.conjugate()
    Q = np.sqrt(2) * alpha.real
    P = np.sqrt(2) * alpha.imag

    sample_covariance = np.cov(np.vstack((Q, P)), ddof=1)
    exact_covariance = np.diag(
        [0.5 * np.exp(-2 * squeeze_r), 0.5 * np.exp(2 * squeeze_r)]
    )
    max_abs_error = float(np.max(np.abs(sample_covariance - exact_covariance)))
    Q_centered = Q - np.mean(Q)
    P_centered = P - np.mean(P)
    covariance_standard_error = np.array(
        [
            [np.std(Q_centered**2, ddof=1), np.std(Q_centered * P_centered, ddof=1)],
            [np.std(Q_centered * P_centered, ddof=1), np.std(P_centered**2, ddof=1)],
        ]
    ) / np.sqrt(trajectories)
    covariance_z = np.abs(sample_covariance - exact_covariance) / covariance_standard_error
    anomalous_samples = alpha * alpha
    anomalous_mean = float(np.mean(anomalous_samples).real)
    anomalous_se = float(np.std(anomalous_samples.real, ddof=1) / np.sqrt(trajectories))
    anomalous_exact = float(-0.5 * np.sinh(2 * squeeze_r))
    anomalous_z = abs(anomalous_mean - anomalous_exact) / anomalous_se
    thresholds = {
        "max_abs_covariance_error": 0.03,
        "max_covariance_z": 5.0,
        "max_anomalous_moment_z": 5.0,
    }
    checks = {
        "covariance_absolute_error_below_cap": bool(
            max_abs_error <= thresholds["max_abs_covariance_error"]
        ),
        "covariance_within_5se": bool(np.max(covariance_z) <= thresholds["max_covariance_z"]),
        "anomalous_moment_within_5se": bool(
            anomalous_z <= thresholds["max_anomalous_moment_z"]
        ),
    }

    metrics = {
        "seed": seed,
        "trajectories": trajectories,
        "squeeze_r": squeeze_r,
        "squeezing_convention": (
            "alpha=cosh(r) beta-sinh(r) beta*; r>0 squeezes Q and E[alpha^2]<0"
        ),
        "sample_covariance": sample_covariance.tolist(),
        "exact_covariance": exact_covariance.tolist(),
        "covariance_standard_error": covariance_standard_error.tolist(),
        "max_covariance_z": float(np.max(covariance_z)),
        "max_abs_covariance_error": max_abs_error,
        "sample_anomalous_moment_real": anomalous_mean,
        "anomalous_moment_standard_error": anomalous_se,
        "exact_anomalous_moment_real": anomalous_exact,
        "anomalous_moment_z": float(anomalous_z),
        "validation": {
            "thresholds": thresholds,
            "checks": checks,
            "all_passed": bool(all(checks.values())),
        },
    }

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    METRICS.parent.mkdir(parents=True, exist_ok=True)
    metrics = finalize_metrics(metrics, run_context, __file__)
    METRICS.write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")

    keep = 5000
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    ax.scatter(Q[:keep], P[:keep], s=4, alpha=0.18, edgecolors="none", color="#1f77b4")
    ax.set_xlabel(r"$Q$")
    ax.set_ylabel(r"$P$")
    ax.set_title(rf"Squeezed-vacuum Wigner samples ($r={squeeze_r}$)")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIGURE)
    plt.close(fig)

    if not metrics["validation"]["all_passed"]:
        failed = [name for name, passed in checks.items() if not passed]
        raise RuntimeError(f"Gaussian sampling validation failed: {failed}")
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
