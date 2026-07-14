"""Monte Carlo check of single-mode squeezed-vacuum Wigner sampling."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FIGURE = ROOT / "figures" / "ch11" / "squeezed_wigner_sampling.pdf"
METRICS = ROOT / "data" / "ch11" / "squeezed_sampling_metrics.json"


def main() -> None:
    seed = 20260714
    trajectories = 100_000
    squeeze_r = 0.8
    rng = np.random.default_rng(seed)

    # Vacuum Wigner samples: E|beta|^2 = 1/2.
    beta = (rng.standard_normal(trajectories) + 1j * rng.standard_normal(trajectories)) / 2
    alpha = np.cosh(squeeze_r) * beta + np.sinh(squeeze_r) * beta.conjugate()
    q = np.sqrt(2) * alpha.real
    p = np.sqrt(2) * alpha.imag

    sample_covariance = np.cov(np.vstack((q, p)), ddof=1)
    exact_covariance = np.diag(
        [0.5 * np.exp(2 * squeeze_r), 0.5 * np.exp(-2 * squeeze_r)]
    )
    max_abs_error = float(np.max(np.abs(sample_covariance - exact_covariance)))

    metrics = {
        "seed": seed,
        "trajectories": trajectories,
        "squeeze_r": squeeze_r,
        "sample_covariance": sample_covariance.tolist(),
        "exact_covariance": exact_covariance.tolist(),
        "max_abs_covariance_error": max_abs_error,
        "sample_anomalous_moment_real": float(np.mean(alpha * alpha).real),
        "exact_anomalous_moment_real": float(0.5 * np.sinh(2 * squeeze_r)),
    }

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    METRICS.parent.mkdir(parents=True, exist_ok=True)
    METRICS.write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")

    keep = 5000
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    ax.scatter(q[:keep], p[:keep], s=4, alpha=0.18, edgecolors="none", color="#1f77b4")
    ax.set_xlabel(r"$q$")
    ax.set_ylabel(r"$p$")
    ax.set_title(rf"Squeezed-vacuum Wigner samples ($r={squeeze_r}$)")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIGURE)
    plt.close(fig)

    if max_abs_error > 0.03:
        raise RuntimeError(f"Covariance check failed: max error {max_abs_error:.4g}")
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
