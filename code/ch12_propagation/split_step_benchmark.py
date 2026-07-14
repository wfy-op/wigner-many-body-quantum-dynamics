"""Actual Strang-splitting benchmark for a 1D Wigner field trajectory."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FIGURE = ROOT / "figures" / "ch12" / "split_step_conservation.pdf"
METRICS = ROOT / "data" / "ch12" / "split_step_metrics.json"


def norm(psi: np.ndarray, dx: float) -> float:
    return float(dx * np.sum(np.abs(psi) ** 2))


def energy(
    psi: np.ndarray,
    k: np.ndarray,
    x: np.ndarray,
    dx: float,
    g: float,
    delta_c: float,
    trap_strength: float,
) -> float:
    derivative = np.fft.ifft(1j * k * np.fft.fft(psi))
    density = np.abs(psi) ** 2
    potential = trap_strength * x**2
    integrand = (
        0.5 * np.abs(derivative) ** 2
        + potential * density
        + 0.5 * g * (density**2 - 2 * delta_c * density)
    )
    return float(dx * np.sum(integrand))


def advance(
    psi: np.ndarray,
    dt: float,
    steps: int,
    k: np.ndarray,
    x: np.ndarray,
    g: float,
    delta_c: float,
    trap_strength: float,
    record: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    kinetic_half = np.exp(-0.25j * k**2 * dt)
    potential = trap_strength * x**2
    dx = float(x[1] - x[0])
    states_norm = []
    states_energy = []
    current = psi.copy()
    if record:
        states_norm.append(norm(current, dx))
        states_energy.append(energy(current, k, x, dx, g, delta_c, trap_strength))
    for _ in range(steps):
        current = np.fft.ifft(kinetic_half * np.fft.fft(current))
        local = potential + g * (np.abs(current) ** 2 - delta_c)
        current *= np.exp(-1j * local * dt)
        current = np.fft.ifft(kinetic_half * np.fft.fft(current))
        if record:
            states_norm.append(norm(current, dx))
            states_energy.append(energy(current, k, x, dx, g, delta_c, trap_strength))
    return current, np.asarray(states_norm), np.asarray(states_energy)


def relative_l2(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def main() -> None:
    seed = 20260714
    points = 64
    length = 20.0
    dx = length / points
    x = (np.arange(points) - points // 2) * dx
    k = 2 * np.pi * np.fft.fftfreq(points, d=dx)
    g = 0.05
    trap_strength = 0.02
    delta_c = 1.0 / dx
    final_time = 0.5
    rng = np.random.default_rng(seed)

    coherent = np.exp(-(x / 2.0) ** 2 / 2).astype(complex)
    coherent *= np.sqrt(100.0 / norm(coherent, dx))
    vacuum_noise = (rng.standard_normal(points) + 1j * rng.standard_normal(points)) / (
        2 * np.sqrt(dx)
    )
    initial = coherent + vacuum_noise

    step_sizes = [0.004, 0.002, 0.001, 0.0005]
    finals = {}
    for dt in step_sizes:
        steps = round(final_time / dt)
        result, _, _ = advance(
            initial, dt, steps, k, x, g, delta_c, trap_strength, record=False
        )
        finals[dt] = result

    coarse_error = relative_l2(finals[0.004], finals[0.002])
    medium_error = relative_l2(finals[0.002], finals[0.001])
    fine_error = relative_l2(finals[0.001], finals[0.0005])
    observed_order = float(np.log2(medium_error / fine_error))

    dt_record = 0.001
    steps_record = round(final_time / dt_record)
    final, norms, energies = advance(
        initial,
        dt_record,
        steps_record,
        k,
        x,
        g,
        delta_c,
        trap_strength,
        record=True,
    )
    back, _, _ = advance(
        final,
        -dt_record,
        steps_record,
        k,
        x,
        g,
        delta_c,
        trap_strength,
        record=False,
    )
    max_norm_drift = float(np.max(np.abs(norms - norms[0])) / abs(norms[0]))
    max_energy_drift = float(np.max(np.abs(energies - energies[0])) / abs(energies[0]))
    reversibility_error = relative_l2(back, initial)

    metrics = {
        "seed": seed,
        "points": points,
        "length": length,
        "dx": dx,
        "g": g,
        "trap_strength": trap_strength,
        "delta_c": delta_c,
        "final_time": final_time,
        "step_sizes": step_sizes,
        "coarse_pair_error": coarse_error,
        "medium_pair_error": medium_error,
        "fine_pair_error": fine_error,
        "observed_order_from_two_finest_pairs": observed_order,
        "max_relative_norm_drift": max_norm_drift,
        "max_relative_energy_drift": max_energy_drift,
        "reversibility_relative_l2_error": reversibility_error,
    }

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    METRICS.parent.mkdir(parents=True, exist_ok=True)
    METRICS.write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")

    times = np.linspace(0.0, final_time, len(norms))
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4))
    axes[0].plot(times, (norms - norms[0]) / norms[0], color="#1f77b4")
    axes[0].set_xlabel("time")
    axes[0].set_ylabel("relative norm drift")
    axes[0].grid(alpha=0.2)
    axes[1].plot(times, (energies - energies[0]) / abs(energies[0]), color="#d62728")
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("relative Weyl-energy drift")
    axes[1].grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIGURE)
    plt.close(fig)

    if observed_order < 1.8:
        raise RuntimeError(f"Second-order convergence check failed: {observed_order:.3f}")
    if max_norm_drift > 1e-11 or reversibility_error > 1e-10:
        raise RuntimeError("Unitary/reversibility check failed")
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
