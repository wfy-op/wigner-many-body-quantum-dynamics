"""Reproduce the coherent-state Wigner rotation used in Chapter 5."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CODE_ROOT = Path(__file__).resolve().parents[1]
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))
from benchmark_contract import begin_metrics_run, finalize_metrics


ROOT = Path(__file__).resolve().parents[2]
FIGURE = ROOT / "figures" / "ch05" / "harmonic_wigner_rotation.pdf"
METRICS = ROOT / "data" / "ch05" / "harmonic_wigner_metrics.json"


def run(seed: int = 20260714, trajectories: int = 20_000) -> dict[str, object]:
    run_context = begin_metrics_run(__file__)
    rng = np.random.default_rng(seed)
    alpha0 = 2.0 + 0.5j
    omega = 1.0
    times = np.array([0.0, 0.5 * np.pi, np.pi])

    noise = (rng.standard_normal(trajectories) + 1j * rng.standard_normal(trajectories)) / 2.0
    initial = alpha0 + noise
    evolved = np.exp(-1j * omega * times[:, None]) * initial[None, :]

    sample_means = evolved.mean(axis=1)
    exact_means = alpha0 * np.exp(-1j * omega * times)
    number_samples = np.abs(evolved) ** 2 - 0.5
    number_estimates = np.mean(number_samples, axis=1)
    number_standard_errors = np.std(number_samples, axis=1, ddof=1) / np.sqrt(trajectories)
    mean_standard_errors_real = np.std(evolved.real, axis=1, ddof=1) / np.sqrt(trajectories)
    mean_standard_errors_imag = np.std(evolved.imag, axis=1, ddof=1) / np.sqrt(trajectories)

    quadrature_variances = []
    quadrature_variance_standard_errors = []
    for samples in evolved:
        Q = np.sqrt(2.0) * samples.real
        P = np.sqrt(2.0) * samples.imag
        variances = np.array([np.var(Q, ddof=1), np.var(P, ddof=1)])
        quadrature_variances.append(variances.tolist())
        # Exact for independent Gaussian samples up to replacing the population
        # variance by the sample variance.
        quadrature_variance_standard_errors.append(
            (variances * np.sqrt(2.0 / (trajectories - 1))).tolist()
        )

    mean_z_real = np.abs(sample_means.real - exact_means.real) / mean_standard_errors_real
    mean_z_imag = np.abs(sample_means.imag - exact_means.imag) / mean_standard_errors_imag
    number_z = np.abs(number_estimates - abs(alpha0) ** 2) / number_standard_errors
    variance_values = np.asarray(quadrature_variances)
    variance_errors = np.asarray(quadrature_variance_standard_errors)
    variance_z = np.abs(variance_values - 0.5) / variance_errors
    thresholds = {
        "max_mean_component_z": 5.0,
        "max_number_z": 5.0,
        "max_quadrature_variance_z": 5.0,
        "max_absolute_mean_error": 0.03,
    }
    checks = {
        "mean_components_within_5se": bool(
            max(np.max(mean_z_real), np.max(mean_z_imag)) <= thresholds["max_mean_component_z"]
        ),
        "number_within_5se": bool(np.max(number_z) <= thresholds["max_number_z"]),
        "quadrature_variances_within_5se": bool(
            np.max(variance_z) <= thresholds["max_quadrature_variance_z"]
        ),
        "absolute_mean_error_below_cap": bool(
            np.max(np.abs(sample_means - exact_means)) <= thresholds["max_absolute_mean_error"]
        ),
    }

    metrics: dict[str, object] = {
        "seed": seed,
        "trajectories": trajectories,
        "alpha0_real": alpha0.real,
        "alpha0_imag": alpha0.imag,
        "times": times.tolist(),
        "max_mean_error": float(np.max(np.abs(sample_means - exact_means))),
        "sample_mean_real": sample_means.real.tolist(),
        "sample_mean_imag": sample_means.imag.tolist(),
        "mean_standard_error_real": mean_standard_errors_real.tolist(),
        "mean_standard_error_imag": mean_standard_errors_imag.tolist(),
        "max_mean_component_z": float(max(np.max(mean_z_real), np.max(mean_z_imag))),
        "number_estimates": number_estimates.tolist(),
        "number_standard_errors": number_standard_errors.tolist(),
        "max_number_z": float(np.max(number_z)),
        "exact_number": float(abs(alpha0) ** 2),
        "quadrature_variances": quadrature_variances,
        "quadrature_variance_standard_errors": quadrature_variance_standard_errors,
        "max_quadrature_variance_z": float(np.max(variance_z)),
        "exact_quadrature_variance": 0.5,
        "validation": {
            "thresholds": thresholds,
            "checks": checks,
            "all_passed": bool(all(checks.values())),
        },
    }

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.0), constrained_layout=True)
    shown = min(2_000, trajectories)
    for ax, t, samples, exact_mean in zip(axes, times, evolved, exact_means, strict=True):
        Q = np.sqrt(2.0) * samples[:shown].real
        P = np.sqrt(2.0) * samples[:shown].imag
        ax.scatter(Q, P, s=3, alpha=0.18, linewidths=0, color="#2878B5")
        index = int(np.argmin(np.abs(times - t)))
        ax.errorbar(
            np.sqrt(2.0) * sample_means[index].real,
            np.sqrt(2.0) * sample_means[index].imag,
            xerr=np.sqrt(2.0) * mean_standard_errors_real[index],
            yerr=np.sqrt(2.0) * mean_standard_errors_imag[index],
            marker="o",
            markersize=4,
            capsize=2,
            color="#2878B5",
            label=r"sample mean $\pm1$ SE" if index == 0 else None,
        )
        ax.plot(np.sqrt(2.0) * exact_mean.real, np.sqrt(2.0) * exact_mean.imag,
                marker="x", markersize=8, markeredgewidth=2, color="#C82423",
                label="exact mean" if index == 0 else None)
        ax.set_title(rf"$\omega t={t:.2f}$")
        ax.set_xlabel(r"$Q$")
        ax.set_ylabel(r"$P$")
        ax.set_aspect("equal")
        ax.set_xlim(-4.2, 4.2)
        ax.set_ylim(-4.2, 4.2)
        ax.grid(alpha=0.2)
    axes[0].legend(frameon=False, fontsize=7, loc="lower left")
    fig.savefig(FIGURE)
    plt.close(fig)

    METRICS.parent.mkdir(parents=True, exist_ok=True)
    metrics = finalize_metrics(metrics, run_context, __file__)
    METRICS.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    if not metrics["validation"]["all_passed"]:
        failed = [name for name, passed in checks.items() if not passed]
        raise RuntimeError(f"Harmonic Wigner validation failed: {failed}")
    return metrics


if __name__ == "__main__":
    result = run()
    print(json.dumps(result, ensure_ascii=False, indent=2))
