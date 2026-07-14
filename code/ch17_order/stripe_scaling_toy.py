"""Actual Monte Carlo toy ensembles for symmetry-restored stripe diagnostics."""

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
FIGURE = ROOT / "figures" / "ch17" / "stripe_scaling_toy.pdf"
METRICS = ROOT / "data" / "ch17" / "stripe_scaling_toy_metrics.json"


def stripe_amplitude(density: np.ndarray, x: np.ndarray, wavevector: float, dx: float) -> complex:
    particles = dx * np.sum(density)
    return complex(dx * np.sum(density * np.exp(-1j * wavevector * x)) / particles)


def loglog_slope_and_standard_error(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    log_x = np.log(x)
    log_y = np.log(y)
    slope, intercept = np.polyfit(log_x, log_y, 1)
    residual = log_y - (slope * log_x + intercept)
    slope_se = np.sqrt(np.sum(residual**2) / (len(x) - 2) / np.sum((log_x - log_x.mean()) ** 2))
    return float(slope), float(slope_se)


def main() -> None:
    run_context = begin_metrics_run(__file__)
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
    ordered_second_se = []
    ordered_mean_abs = []
    domain_second = []
    domain_second_se = []

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

        ordered_power = np.abs(ordered_m) ** 2
        domain_power = np.abs(domain_m) ** 2
        ordered_second.append(float(np.mean(ordered_power)))
        ordered_second_se.append(float(np.std(ordered_power, ddof=1) / np.sqrt(trajectories)))
        ordered_mean_abs.append(float(abs(np.mean(ordered_m))))
        domain_second.append(float(np.mean(domain_power)))
        domain_second_se.append(float(np.std(domain_power, ddof=1) / np.sqrt(trajectories)))

    ordered_second_array = np.asarray(ordered_second)
    domain_second_array = np.asarray(domain_second)
    ordered_second_se_array = np.asarray(ordered_second_se)
    domain_second_se_array = np.asarray(domain_second_se)
    ordered_slope, ordered_slope_se = loglog_slope_and_standard_error(
        lengths, ordered_second_array
    )
    domain_slope, domain_slope_se = loglog_slope_and_standard_error(
        lengths, domain_second_array
    )
    exact_plateau = amplitude**2 / 4.0
    plateau_relative_error = float(
        np.max(np.abs(ordered_second_array - exact_plateau)) / exact_plateau
    )
    residual_mean_fraction = float(np.max(ordered_mean_abs) / (amplitude / 2.0))
    thresholds = {
        "max_abs_ordered_loglog_slope": 0.08,
        "finite_domain_slope_minimum": -1.12,
        "finite_domain_slope_maximum": -0.88,
        "max_ordered_plateau_relative_error": 0.02,
        "max_symmetry_restored_mean_fraction": 0.06,
    }
    checks = {
        "ordered_ensemble_approaches_plateau": bool(
            abs(ordered_slope) <= thresholds["max_abs_ordered_loglog_slope"]
        ),
        "finite_domains_show_inverse_length_scaling": bool(
            thresholds["finite_domain_slope_minimum"]
            <= domain_slope
            <= thresholds["finite_domain_slope_maximum"]
        ),
        "ordered_plateau_matches_exact_amplitude": bool(
            plateau_relative_error <= thresholds["max_ordered_plateau_relative_error"]
        ),
        "random_global_phase_restores_mean_order_parameter": bool(
            residual_mean_fraction <= thresholds["max_symmetry_restored_mean_fraction"]
        ),
    }

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
        "ordered_mean_abs_mq_squared_standard_error": ordered_second_se,
        "finite_domain_mean_abs_mq_squared": domain_second,
        "finite_domain_mean_abs_mq_squared_standard_error": domain_second_se,
        "ordered_loglog_slope": ordered_slope,
        "ordered_loglog_slope_standard_error": ordered_slope_se,
        "finite_domain_loglog_slope": domain_slope,
        "finite_domain_loglog_slope_standard_error": domain_slope_se,
        "exact_ordered_plateau": exact_plateau,
        "ordered_plateau_max_relative_error": plateau_relative_error,
        "symmetry_restored_mean_fraction": residual_mean_fraction,
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

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.5))
    axes[0].plot(lengths, ordered_mean_abs, "o-", label=r"$|\mathbb{E}m_Q|$")
    axes[0].plot(lengths, np.sqrt(ordered_second_array), "s-", label=r"$\sqrt{\mathbb{E}|m_Q|^2}$")
    axes[0].axhline(amplitude / 2.0, color="black", ls="--", lw=1, label="exact amplitude")
    axes[0].set_xlabel(r"$L$")
    axes[0].set_ylabel("stripe amplitude")
    axes[0].grid(alpha=0.2)
    axes[0].legend(frameon=False)
    axes[1].errorbar(
        lengths,
        ordered_second_array,
        yerr=ordered_second_se_array,
        fmt="o-",
        capsize=2,
        label="global random phase",
    )
    axes[1].errorbar(
        lengths,
        domain_second_array,
        yerr=domain_second_se_array,
        fmt="s-",
        capsize=2,
        label="finite phase domains",
    )
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"$L$")
    axes[1].set_ylabel(r"$\mathbb{E}|m_Q|^2$")
    axes[1].grid(alpha=0.2, which="both")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURE)
    plt.close(fig)

    if not metrics["validation"]["all_passed"]:
        failed = [name for name, passed in checks.items() if not passed]
        raise RuntimeError(f"Stripe-scaling validation failed: {failed}")
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
