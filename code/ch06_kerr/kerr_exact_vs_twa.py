"""Compare exact coherent-state Kerr dynamics with the truncated Wigner result."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FIGURE = ROOT / "figures" / "ch06" / "kerr_exact_vs_twa.pdf"
METRICS = ROOT / "data" / "ch06" / "kerr_metrics.json"


def exact_mean(alpha0: complex, tau: np.ndarray) -> np.ndarray:
    occupation = abs(alpha0) ** 2
    return alpha0 * np.exp(occupation * (np.exp(-1j * tau) - 1.0))


def run(seed: int = 20260714, trajectories: int = 100_000) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    occupation = 9.0
    alpha0 = np.sqrt(occupation) + 0.0j
    tau = np.linspace(0.0, 2.0 * np.pi, 401)

    initial = alpha0 + (rng.standard_normal(trajectories) + 1j * rng.standard_normal(trajectories)) / 2.0
    angular_frequency = np.abs(initial) ** 2 - 1.0

    twa = np.empty(tau.size, dtype=np.complex128)
    for index, time in enumerate(tau):
        twa[index] = np.mean(initial * np.exp(-1j * angular_frequency * time))

    exact = exact_mean(alpha0, tau)
    exact_norm = np.abs(exact) / abs(alpha0)
    twa_norm = np.abs(twa) / abs(alpha0)
    early = tau <= 0.25
    revival_index = int(np.argmin(np.abs(tau - 2.0 * np.pi)))

    metrics: dict[str, object] = {
        "seed": seed,
        "trajectories": trajectories,
        "initial_occupation": occupation,
        "max_early_coherence_error_tau_le_0.25": float(np.max(np.abs(exact_norm[early] - twa_norm[early]))),
        "exact_revival_coherence": float(exact_norm[revival_index]),
        "twa_revival_coherence": float(twa_norm[revival_index]),
        "twa_number_estimate_initial": float(np.mean(np.abs(initial) ** 2) - 0.5),
        "exact_number": occupation,
    }

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.4, 3.7), constrained_layout=True)
    ax.plot(tau, exact_norm, color="#C82423", linewidth=2.0, label="Exact")
    ax.plot(tau, twa_norm, color="#2878B5", linewidth=1.7, linestyle="--", label="TWA")
    ax.axvline(2.0 * np.pi, color="0.5", linewidth=1.0, linestyle=":")
    ax.set_xlabel(r"$\tau=\chi t$")
    ax.set_ylabel(r"$|\langle \hat a(t)\rangle|/|\alpha_0|$")
    ax.set_xlim(tau[0], tau[-1])
    ax.set_ylim(-0.02, 1.05)
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.savefig(FIGURE)
    plt.close(fig)

    METRICS.parent.mkdir(parents=True, exist_ok=True)
    METRICS.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    return metrics


if __name__ == "__main__":
    print(json.dumps(run(), ensure_ascii=False, indent=2))

