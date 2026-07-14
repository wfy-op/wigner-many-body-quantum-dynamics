"""Check normally ordered number and g2 estimators on exact Gaussian Wigner states."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FIGURE = ROOT / "figures" / "ch13" / "ordering_estimators.pdf"
METRICS = ROOT / "data" / "ch13" / "ordering_estimator_metrics.json"


def estimators(alpha: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    intensity = np.abs(alpha) ** 2
    number = intensity - 0.5
    factorial_second = intensity**2 - 2.0 * intensity + 0.5
    return number, factorial_second


def grouped_jackknife_ratio(
    number: np.ndarray, factorial_second: np.ndarray, groups: int
) -> tuple[float, float]:
    size = len(number)
    if size % groups:
        raise ValueError("Sample count must be divisible by the number of groups")
    block = size // groups
    n_sums = number.reshape(groups, block).sum(axis=1)
    f_sums = factorial_second.reshape(groups, block).sum(axis=1)
    total_n = float(n_sums.sum())
    total_f = float(f_sums.sum())
    retained = size - block
    leave_n = (total_n - n_sums) / retained
    leave_f = (total_f - f_sums) / retained
    leave_ratio = leave_f / leave_n**2
    mean_leave = float(leave_ratio.mean())
    standard_error = float(
        np.sqrt((groups - 1) / groups * np.sum((leave_ratio - mean_leave) ** 2))
    )
    estimate = float((total_f / size) / (total_n / size) ** 2)
    return estimate, standard_error


def summarize(alpha: np.ndarray, exact_n: float, exact_g2: float) -> dict[str, float]:
    number, factorial_second = estimators(alpha)
    n_mean = float(number.mean())
    g2_estimate, g2_se = grouped_jackknife_ratio(number, factorial_second, groups=200)
    f_mean = float(factorial_second.mean())
    n_se = float(number.std(ddof=1) / np.sqrt(len(number)))
    f_se = float(factorial_second.std(ddof=1) / np.sqrt(len(number)))
    return {
        "number_estimate": n_mean,
        "number_standard_error": n_se,
        "number_exact": exact_n,
        "factorial_second_estimate": f_mean,
        "factorial_second_standard_error": f_se,
        "factorial_second_exact": exact_g2 * exact_n**2,
        "g2_estimate": g2_estimate,
        "g2_jackknife_standard_error": g2_se,
        "g2_exact": exact_g2,
    }


def main() -> None:
    seed = 20260714
    trajectories = 400_000
    rng = np.random.default_rng(seed)

    coherent_n = 4.0
    coherent = np.sqrt(coherent_n) + (
        rng.standard_normal(trajectories) + 1j * rng.standard_normal(trajectories)
    ) / 2

    thermal_n = 2.0
    thermal = np.sqrt((thermal_n + 0.5) / 2) * (
        rng.standard_normal(trajectories) + 1j * rng.standard_normal(trajectories)
    )

    coherent_metrics = summarize(coherent, exact_n=coherent_n, exact_g2=1.0)
    thermal_metrics = summarize(thermal, exact_n=thermal_n, exact_g2=2.0)
    metrics = {
        "seed": seed,
        "trajectories_per_state": trajectories,
        "jackknife_groups": 200,
        "coherent": coherent_metrics,
        "thermal": thermal_metrics,
    }

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    METRICS.parent.mkdir(parents=True, exist_ok=True)
    METRICS.write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")

    labels = ["coherent", "thermal"]
    estimates = [coherent_metrics["g2_estimate"], thermal_metrics["g2_estimate"]]
    errors = [
        coherent_metrics["g2_jackknife_standard_error"],
        thermal_metrics["g2_jackknife_standard_error"],
    ]
    exact = [1.0, 2.0]
    positions = np.arange(2)
    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    ax.errorbar(
        positions,
        estimates,
        yerr=errors,
        fmt="o",
        capsize=4,
        color="#1f77b4",
        label="Wigner estimator",
    )
    ax.scatter(positions, exact, marker="x", s=70, color="#d62728", label="exact")
    ax.set_xticks(positions, labels)
    ax.set_ylabel(r"$g^{(2)}$")
    ax.set_ylim(0.8, 2.2)
    ax.grid(axis="y", alpha=0.2)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURE)
    plt.close(fig)

    for state in (coherent_metrics, thermal_metrics):
        z = abs(state["g2_estimate"] - state["g2_exact"]) / state[
            "g2_jackknife_standard_error"
        ]
        if z > 4.0:
            raise RuntimeError(f"g2 estimator differs from exact value by {z:.2f} standard errors")
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
