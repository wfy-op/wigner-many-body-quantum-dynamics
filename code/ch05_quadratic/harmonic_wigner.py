"""Reproduce the coherent-state Wigner rotation used in Chapter 5."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FIGURE = ROOT / "figures" / "ch05" / "harmonic_wigner_rotation.pdf"
METRICS = ROOT / "data" / "ch05" / "harmonic_wigner_metrics.json"


def run(seed: int = 20260714, trajectories: int = 20_000) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    alpha0 = 2.0 + 0.5j
    omega = 1.0
    times = np.array([0.0, 0.5 * np.pi, np.pi])

    noise = (rng.standard_normal(trajectories) + 1j * rng.standard_normal(trajectories)) / 2.0
    initial = alpha0 + noise
    evolved = np.exp(-1j * omega * times[:, None]) * initial[None, :]

    sample_means = evolved.mean(axis=1)
    exact_means = alpha0 * np.exp(-1j * omega * times)
    number_estimates = np.mean(np.abs(evolved) ** 2, axis=1) - 0.5

    quadrature_variances = []
    for samples in evolved:
        q = np.sqrt(2.0) * samples.real
        p = np.sqrt(2.0) * samples.imag
        quadrature_variances.append([np.var(q, ddof=1), np.var(p, ddof=1)])

    metrics: dict[str, object] = {
        "seed": seed,
        "trajectories": trajectories,
        "alpha0_real": alpha0.real,
        "alpha0_imag": alpha0.imag,
        "times": times.tolist(),
        "max_mean_error": float(np.max(np.abs(sample_means - exact_means))),
        "number_estimates": number_estimates.tolist(),
        "exact_number": float(abs(alpha0) ** 2),
        "quadrature_variances": quadrature_variances,
        "exact_quadrature_variance": 0.5,
    }

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.0), constrained_layout=True)
    shown = min(2_000, trajectories)
    for ax, t, samples, exact_mean in zip(axes, times, evolved, exact_means, strict=True):
        q = np.sqrt(2.0) * samples[:shown].real
        p = np.sqrt(2.0) * samples[:shown].imag
        ax.scatter(q, p, s=3, alpha=0.18, linewidths=0, color="#2878B5")
        ax.plot(np.sqrt(2.0) * exact_mean.real, np.sqrt(2.0) * exact_mean.imag,
                marker="x", markersize=8, markeredgewidth=2, color="#C82423")
        ax.set_title(rf"$\omega t={t:.2f}$")
        ax.set_xlabel(r"$q$")
        ax.set_ylabel(r"$p$")
        ax.set_aspect("equal")
        ax.set_xlim(-4.2, 4.2)
        ax.set_ylim(-4.2, 4.2)
        ax.grid(alpha=0.2)
    fig.savefig(FIGURE)
    plt.close(fig)

    METRICS.parent.mkdir(parents=True, exist_ok=True)
    METRICS.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    return metrics


if __name__ == "__main__":
    result = run()
    print(json.dumps(result, ensure_ascii=False, indent=2))

