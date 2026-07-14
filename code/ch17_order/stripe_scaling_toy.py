"""Actual Monte Carlo toy ensembles for symmetry-restored stripe diagnostics."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FIGURE = ROOT / "figures" / "ch17" / "stripe_scaling_toy.pdf"
METRICS = ROOT / "data" / "ch17" / "stripe_scaling_toy_metrics.json"


def stripe_amplitude(density: np.ndarray, x: np.ndarray, wavevector: float, dx: float) -> complex:
    particles = dx * np.sum(density)
    return complex(dx * np.sum(density * np.exp(-1j * wavevector * x)) / particles)


def main() -> None:
    seed = 20260714
    trajectories = 2000
    lengths = np.array([32.0, 64.0, 128.0, 256.0])
    dx = 0.25
    period = 4.0
    domain_length = 4.0
    amplitude = 0.3
    density_mean = 1.0
    wavevector = 2.0 * np.pi / period
    rng = np.random.default_rng(seed)

    ordered_second = []
    ordered_mean_abs = []
    domain_second = []

    for length in lengths:
        points = round(length / dx)
        x = np.arange(points) * dx
        ordered_m = np.empty(trajectories, dtype=complex)
        domain_m = np.empty(trajectories, dtype=complex)
        points_per_domain = round(domain_length / dx)
        domains = points // points_per_domain

        for sample in range(trajectories):
            global_phase = rng.uniform(0.0, 2.0 * np.pi)
            ordered_density = density_mean * (
                1.0 + amplitude * np.cos(wavevector * x + global_phase)
            )
            ordered_m[sample] = stripe_amplitude(ordered_density, x, wavevector, dx)

            phases = rng.uniform(0.0, 2.0 * np.pi, size=domains)
            phase_field = np.repeat(phases, points_per_domain)[:points]
            domain_density = density_mean * (
                1.0 + amplitude * np.cos(wavevector * x + phase_field)
            )
            domain_m[sample] = stripe_amplitude(domain_density, x, wavevector, dx)

        ordered_second.append(float(np.mean(np.abs(ordered_m) ** 2)))
        ordered_mean_abs.append(float(abs(np.mean(ordered_m))))
        domain_second.append(float(np.mean(np.abs(domain_m) ** 2)))

    ordered_second_array = np.asarray(ordered_second)
    domain_second_array = np.asarray(domain_second)
    ordered_slope = float(np.polyfit(np.log(lengths), np.log(ordered_second_array), 1)[0])
    domain_slope = float(np.polyfit(np.log(lengths), np.log(domain_second_array), 1)[0])

    metrics = {
        "seed": seed,
        "trajectories_per_length": trajectories,
        "lengths": lengths.tolist(),
        "dx": dx,
        "stripe_period": period,
        "domain_length": domain_length,
        "stripe_amplitude": amplitude,
        "ordered_abs_mean_mq": ordered_mean_abs,
        "ordered_mean_abs_mq_squared": ordered_second,
        "finite_domain_mean_abs_mq_squared": domain_second,
        "ordered_loglog_slope": ordered_slope,
        "finite_domain_loglog_slope": domain_slope,
        "exact_ordered_plateau": amplitude**2 / 4.0,
    }

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    METRICS.parent.mkdir(parents=True, exist_ok=True)
    METRICS.write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.5))
    axes[0].plot(lengths, ordered_mean_abs, "o-", label=r"$|\mathbb{E}m_Q|$")
    axes[0].plot(lengths, np.sqrt(ordered_second_array), "s-", label=r"$\sqrt{\mathbb{E}|m_Q|^2}$")
    axes[0].axhline(amplitude / 2.0, color="black", ls="--", lw=1, label="exact amplitude")
    axes[0].set_xlabel(r"$L$")
    axes[0].set_ylabel("stripe amplitude")
    axes[0].grid(alpha=0.2)
    axes[0].legend(frameon=False)
    axes[1].loglog(lengths, ordered_second_array, "o-", label="global random phase")
    axes[1].loglog(lengths, domain_second_array, "s-", label="finite phase domains")
    axes[1].set_xlabel(r"$L$")
    axes[1].set_ylabel(r"$\mathbb{E}|m_Q|^2$")
    axes[1].grid(alpha=0.2, which="both")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURE)
    plt.close(fig)

    if abs(ordered_slope) > 0.08:
        raise RuntimeError("Ordered ensemble did not approach a plateau")
    if not (-1.12 < domain_slope < -0.88):
        raise RuntimeError("Finite-domain ensemble did not show inverse-length scaling")
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
