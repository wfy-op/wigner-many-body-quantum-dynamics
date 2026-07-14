"""Monte Carlo validation of a stable two-component Bogoliubov quench."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FIGURE = ROOT / "figures" / "ch15" / "quench_covariance.pdf"
METRICS = ROOT / "data" / "ch15" / "quench_covariance_metrics.json"


def uv_spin(k: np.ndarray, coupling: float, density: float, g_minus: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    kinetic = 0.5 * k**2
    pair = 0.5 * density * g_minus
    diagonal = kinetic + 2.0 * coupling + pair
    energy = np.sqrt((kinetic + 2.0 * coupling) * (kinetic + 2.0 * coupling + density * g_minus))
    u = np.sqrt(0.5 * (diagonal / energy + 1.0))
    v = np.sign(pair) * np.sqrt(0.5 * (diagonal / energy - 1.0))
    return u, v, energy


def main() -> None:
    seed = 20260714
    trajectories = 120_000
    density = 1.0
    g = 1.0
    g12 = 0.6
    coupling_initial = 1.0
    coupling_final = 0.1
    momenta = np.linspace(0.15, 3.0, 24)
    rng = np.random.default_rng(seed)

    ui, vi, _ = uv_spin(momenta, coupling_initial, density, g - g12)
    uf, vf, final_energy = uv_spin(momenta, coupling_final, density, g - g12)
    transform_a = uf * ui - vf * vi
    transform_b = vf * ui - uf * vi
    symplectic_error = float(np.max(np.abs(transform_a**2 - transform_b**2 - 1.0)))

    exact_occupation = transform_b**2
    exact_anomalous = transform_a * transform_b
    sampled_occupation = np.empty_like(momenta)
    sampled_anomalous = np.empty_like(momenta)

    for index in range(len(momenta)):
        b_plus = (rng.standard_normal(trajectories) + 1j * rng.standard_normal(trajectories)) / 2
        b_minus = (rng.standard_normal(trajectories) + 1j * rng.standard_normal(trajectories)) / 2
        c_plus = transform_a[index] * b_plus + transform_b[index] * b_minus.conjugate()
        c_minus = transform_a[index] * b_minus + transform_b[index] * b_plus.conjugate()
        sampled_occupation[index] = np.mean(np.abs(c_plus) ** 2) - 0.5
        sampled_anomalous[index] = np.mean(c_plus * c_minus).real

    max_occupation_error = float(np.max(np.abs(sampled_occupation - exact_occupation)))
    max_anomalous_error = float(np.max(np.abs(sampled_anomalous - exact_anomalous)))
    metrics = {
        "seed": seed,
        "trajectories_per_momentum_pair": trajectories,
        "density": density,
        "g": g,
        "g12": g12,
        "coupling_initial": coupling_initial,
        "coupling_final": coupling_final,
        "momentum_min": float(momenta.min()),
        "momentum_max": float(momenta.max()),
        "momentum_count": len(momenta),
        "minimum_final_energy": float(final_energy.min()),
        "max_symplectic_constraint_error": symplectic_error,
        "max_occupation_absolute_error": max_occupation_error,
        "max_anomalous_absolute_error": max_anomalous_error,
    }

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    METRICS.parent.mkdir(parents=True, exist_ok=True)
    METRICS.write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4))
    axes[0].plot(momenta, exact_occupation, color="#1f77b4", label="exact")
    axes[0].scatter(momenta, sampled_occupation, s=18, facecolors="none", edgecolors="#d62728", label="MC")
    axes[0].set_xlabel(r"$k$")
    axes[0].set_ylabel(r"$\langle c_k^\dagger c_k\rangle$")
    axes[0].grid(alpha=0.2)
    axes[0].legend(frameon=False)
    axes[1].plot(momenta, exact_anomalous, color="#1f77b4", label="exact")
    axes[1].scatter(momenta, sampled_anomalous, s=18, facecolors="none", edgecolors="#d62728", label="MC")
    axes[1].set_xlabel(r"$k$")
    axes[1].set_ylabel(r"$\mathrm{Re}\langle c_k c_{-k}\rangle$")
    axes[1].grid(alpha=0.2)
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURE)
    plt.close(fig)

    if symplectic_error > 1e-12:
        raise RuntimeError("Bogoliubov map failed the symplectic constraint")
    if max(max_occupation_error, max_anomalous_error) > 0.006:
        raise RuntimeError("Monte Carlo covariance check failed")
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
