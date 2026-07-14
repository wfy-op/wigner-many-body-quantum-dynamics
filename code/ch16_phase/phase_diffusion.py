"""Gaussian phase-diffusion trajectories for the collective relative mode."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FIGURE = ROOT / "figures" / "ch16" / "phase_diffusion.pdf"
METRICS = ROOT / "data" / "ch16" / "phase_diffusion_metrics.json"


def main() -> None:
    seed = 20260714
    trajectories = 200_000
    charging_energy = 0.1
    sigma_number = 3.0
    sigma_phase = 0.15
    covariance = 0.0
    times = np.linspace(0.0, 8.0, 161)
    rng = np.random.default_rng(seed)

    covariance_matrix = np.array(
        [[sigma_phase**2, covariance], [covariance, sigma_number**2]]
    )
    initial = rng.multivariate_normal([0.0, 0.0], covariance_matrix, size=trajectories)
    phase_initial = initial[:, 0]
    number_initial = initial[:, 1]

    sampled_variance = np.empty_like(times)
    sampled_coherence = np.empty_like(times)
    for index, time in enumerate(times):
        phase = phase_initial + charging_energy * time * number_initial
        sampled_variance[index] = np.var(phase, ddof=1)
        sampled_coherence[index] = abs(np.mean(np.exp(1j * phase)))

    exact_variance = (
        sigma_phase**2
        + 2.0 * charging_energy * times * covariance
        + (charging_energy * times) ** 2 * sigma_number**2
    )
    exact_coherence = np.exp(-0.5 * exact_variance)
    max_variance_error = float(np.max(np.abs(sampled_variance - exact_variance)))
    max_coherence_error = float(np.max(np.abs(sampled_coherence - exact_coherence)))

    metrics = {
        "seed": seed,
        "trajectories": trajectories,
        "hbar": 1.0,
        "charging_energy": charging_energy,
        "sigma_number": sigma_number,
        "sigma_phase": sigma_phase,
        "initial_phase_number_covariance": covariance,
        "phase_diffusion_time": 1.0 / (charging_energy * sigma_number),
        "max_phase_variance_absolute_error": max_variance_error,
        "max_coherence_absolute_error": max_coherence_error,
    }

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    METRICS.parent.mkdir(parents=True, exist_ok=True)
    METRICS.write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4))
    axes[0].plot(times, exact_variance, color="#1f77b4", label="exact Gaussian")
    axes[0].plot(times[::8], sampled_variance[::8], "o", ms=3, mfc="none", color="#d62728", label="MC")
    axes[0].set_xlabel("time")
    axes[0].set_ylabel(r"$\mathrm{Var}(\phi)$")
    axes[0].grid(alpha=0.2)
    axes[0].legend(frameon=False)
    axes[1].plot(times, exact_coherence, color="#1f77b4", label="exact Gaussian")
    axes[1].plot(times[::8], sampled_coherence[::8], "o", ms=3, mfc="none", color="#d62728", label="MC")
    axes[1].set_xlabel("time")
    axes[1].set_ylabel(r"$|\langle e^{i\phi}\rangle|$")
    axes[1].grid(alpha=0.2)
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURE)
    plt.close(fig)

    if max_coherence_error > 0.004 or max_variance_error > 0.04:
        raise RuntimeError("Phase-diffusion Monte Carlo check failed")
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
