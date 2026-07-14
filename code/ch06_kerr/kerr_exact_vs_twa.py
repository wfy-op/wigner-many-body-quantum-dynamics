"""Compare exact coherent-state Kerr dynamics with the truncated Wigner result."""

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
FIGURE = ROOT / "figures" / "ch06" / "kerr_exact_vs_twa.pdf"
METRICS = ROOT / "data" / "ch06" / "kerr_metrics.json"
TIMESERIES = ROOT / "data" / "ch06" / "kerr_timeseries.json"


def exact_mean(alpha0: complex, tau: np.ndarray) -> np.ndarray:
    occupation = abs(alpha0) ** 2
    return alpha0 * np.exp(occupation * (np.exp(-1j * tau) - 1.0))


def run(seed: int = 20260714, trajectories: int = 100_000) -> dict[str, object]:
    run_context = begin_metrics_run(__file__)
    child_count = 4
    if trajectories % child_count:
        raise ValueError("Trajectory count must be divisible by the child-stream count")
    seed_sequence = np.random.SeedSequence(seed)
    child_sequences = seed_sequence.spawn(child_count)
    occupation = 9.0
    alpha0 = np.sqrt(occupation) + 0.0j
    tau = np.linspace(0.0, 2.0 * np.pi, 401)

    group_size = trajectories // child_count
    initial_groups = []
    for child in child_sequences:
        rng = np.random.default_rng(child)
        initial_groups.append(
            alpha0
            + (rng.standard_normal(group_size) + 1j * rng.standard_normal(group_size)) / 2.0
        )
    initial = np.concatenate(initial_groups)
    angular_frequency = np.abs(initial) ** 2 - 1.0

    twa = np.empty(tau.size, dtype=np.complex128)
    twa_amplitude_se = np.empty(tau.size)
    child_twa = np.empty((child_count, tau.size), dtype=np.complex128)
    for index, time in enumerate(tau):
        values = initial * np.exp(-1j * angular_frequency * time)
        twa[index] = np.mean(values)
        rms_complex_se = np.sqrt(np.var(values.real, ddof=1) + np.var(values.imag, ddof=1))
        rms_complex_se /= np.sqrt(trajectories) * abs(alpha0)
        if abs(twa[index]) > 5.0 * rms_complex_se * abs(alpha0):
            phase = np.conj(twa[index]) / abs(twa[index])
            projected = np.real(values * phase)
            twa_amplitude_se[index] = np.std(projected, ddof=1) / (
                np.sqrt(trajectories) * abs(alpha0)
            )
        else:
            # Near zero the amplitude is non-Gaussian and its radial direction is
            # ill-defined.  The complex-mean RMS error is a conservative display scale.
            twa_amplitude_se[index] = rms_complex_se
        for group in range(child_count):
            start = group * group_size
            child_twa[group, index] = np.mean(values[start : start + group_size])

    exact = exact_mean(alpha0, tau)
    exact_norm = np.abs(exact) / abs(alpha0)
    twa_norm = np.abs(twa) / abs(alpha0)
    early = tau <= 0.25
    revival_index = int(np.argmin(np.abs(tau - 2.0 * np.pi)))
    number_samples = np.abs(initial) ** 2 - 0.5
    number_estimate = float(np.mean(number_samples))
    number_se = float(np.std(number_samples, ddof=1) / np.sqrt(trajectories))
    child_norm = np.abs(child_twa) / abs(alpha0)
    child_early_errors = np.max(np.abs(child_norm[:, early] - exact_norm[None, early]), axis=1)
    child_revivals = child_norm[:, revival_index]
    thresholds = {
        "max_early_coherence_error_tau_le_0.25": 0.02,
        "revival_upper_4se": 0.05,
        "initial_number_z": 5.0,
        "max_child_stream_early_error": 0.04,
        "exact_revival_minimum": 0.999,
    }
    number_z = abs(number_estimate - occupation) / number_se
    checks = {
        "early_window_matches_exact": bool(
            np.max(np.abs(exact_norm[early] - twa_norm[early]))
            <= thresholds["max_early_coherence_error_tau_le_0.25"]
        ),
        "twa_misses_revival_beyond_sampling_noise": bool(
            twa_norm[revival_index] + 4.0 * twa_amplitude_se[revival_index]
            <= thresholds["revival_upper_4se"]
        ),
        "initial_number_within_5se": bool(number_z <= thresholds["initial_number_z"]),
        "independent_child_streams_stable_early": bool(
            np.max(child_early_errors) <= thresholds["max_child_stream_early_error"]
        ),
        "exact_solution_revival_present": bool(
            exact_norm[revival_index] >= thresholds["exact_revival_minimum"]
        ),
    }

    metrics: dict[str, object] = {
        "seed": seed,
        "child_stream_spawn_keys": [list(child.spawn_key) for child in child_sequences],
        "trajectories": trajectories,
        "initial_occupation": occupation,
        "max_early_coherence_error_tau_le_0.25": float(np.max(np.abs(exact_norm[early] - twa_norm[early]))),
        "exact_revival_coherence": float(exact_norm[revival_index]),
        "twa_revival_coherence": float(twa_norm[revival_index]),
        "twa_revival_amplitude_standard_error": float(twa_amplitude_se[revival_index]),
        "child_stream_max_early_errors": child_early_errors.tolist(),
        "child_stream_revival_coherences": child_revivals.tolist(),
        "twa_number_estimate_initial": number_estimate,
        "twa_number_standard_error_initial": number_se,
        "twa_number_z": float(number_z),
        "exact_number": occupation,
        "uncertainty_note": (
            "Amplitude error uses a delta-method projected SE away from zero and the "
            "complex-mean RMS SE near zero; the plotted band is plus/minus 2 SE."
        ),
        "validation": {
            "thresholds": thresholds,
            "checks": checks,
            "all_passed": bool(all(checks.values())),
        },
    }

    timeseries = {
        "tau": tau.tolist(),
        "exact_coherence": exact_norm.tolist(),
        "twa_coherence": twa_norm.tolist(),
        "twa_amplitude_standard_error": twa_amplitude_se.tolist(),
    }

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.4, 3.7), constrained_layout=True)
    ax.plot(tau, exact_norm, color="#C82423", linewidth=2.0, label="Exact")
    ax.plot(tau, twa_norm, color="#2878B5", linewidth=1.7, linestyle="--", label="TWA")
    ax.fill_between(
        tau,
        np.maximum(0.0, twa_norm - 2.0 * twa_amplitude_se),
        twa_norm + 2.0 * twa_amplitude_se,
        color="#2878B5",
        alpha=0.18,
        linewidth=0,
        label=r"TWA $\pm2$ SE",
    )
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
    metrics = finalize_metrics(metrics, run_context, __file__)
    METRICS.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    TIMESERIES.write_text(json.dumps(timeseries, indent=2), encoding="utf-8")
    if not metrics["validation"]["all_passed"]:
        failed = [name for name, passed in checks.items() if not passed]
        raise RuntimeError(f"Kerr validation failed: {failed}")
    return metrics


if __name__ == "__main__":
    print(json.dumps(run(), ensure_ascii=False, indent=2))
