"""Exact two-mode Bose-Hubbard dynamics compared with TWA."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.special import gammaln


CODE_ROOT = Path(__file__).resolve().parents[1]
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))
from benchmark_contract import begin_metrics_run, finalize_metrics


ROOT = Path(__file__).resolve().parents[2]
FIGURE = ROOT / "figures" / "ch08" / "dimer_exact_vs_twa.pdf"
METRICS = ROOT / "data" / "ch08" / "dimer_metrics.json"
TIMESERIES = ROOT / "data" / "ch08" / "dimer_timeseries.json"


def annihilation(cutoff: int) -> sparse.csr_matrix:
    values = np.sqrt(np.arange(1, cutoff, dtype=float))
    return sparse.diags(values, offsets=1, shape=(cutoff, cutoff), format="csr")


def coherent_vector(alpha: complex, cutoff: int) -> np.ndarray:
    if alpha == 0:
        coefficients = np.zeros(cutoff, dtype=np.complex128)
        coefficients[0] = 1.0
        return coefficients
    n = np.arange(cutoff)
    coefficients = np.exp(-0.5 * abs(alpha) ** 2 + n * np.log(alpha + 0j) - 0.5 * gammaln(n + 1.0))
    return coefficients


def exact_dynamics(alpha_l: complex, alpha_r: complex, j_hop: float, interaction: float,
                   times: np.ndarray, cutoff: int = 24) -> tuple[np.ndarray, float]:
    a = annihilation(cutoff)
    eye = sparse.identity(cutoff, format="csr")
    a_l = sparse.kron(a, eye, format="csr")
    a_r = sparse.kron(eye, a, format="csr")
    n_l = a_l.getH() @ a_l
    n_r = a_r.getH() @ a_r
    hamiltonian = (
        -j_hop * (a_l.getH() @ a_r + a_r.getH() @ a_l)
        + 0.5 * interaction * (n_l @ (n_l - sparse.identity(cutoff**2))
                               + n_r @ (n_r - sparse.identity(cutoff**2)))
    )
    psi0 = np.kron(coherent_vector(alpha_l, cutoff), coherent_vector(alpha_r, cutoff))
    retained_norm = float(np.vdot(psi0, psi0).real)
    psi0 = psi0 / np.sqrt(retained_norm)
    states = expm_multiply(-1j * hamiltonian, psi0, start=times[0], stop=times[-1], num=times.size)
    left = np.einsum("ti,ij,tj->t", states.conj(), n_l.toarray(), states, optimize=True).real
    right = np.einsum("ti,ij,tj->t", states.conj(), n_r.toarray(), states, optimize=True).real
    return (left - right) / (left + right), retained_norm


def twa_dynamics(alpha_l: complex, alpha_r: complex, j_hop: float, interaction: float,
                 times: np.ndarray, trajectories: int, seed: int,
                 child_count: int = 4) -> tuple[np.ndarray, np.ndarray, float, np.ndarray, list[list[int]]]:
    if trajectories % child_count:
        raise ValueError("Trajectory count must be divisible by the child-stream count")
    seed_sequence = np.random.SeedSequence(seed)
    child_sequences = seed_sequence.spawn(child_count)
    group_size = trajectories // child_count
    left_groups = []
    right_groups = []
    for child in child_sequences:
        rng = np.random.default_rng(child)
        left_groups.append(
            alpha_l
            + (rng.standard_normal(group_size) + 1j * rng.standard_normal(group_size)) / 2.0
        )
        right_groups.append(
            alpha_r
            + (rng.standard_normal(group_size) + 1j * rng.standard_normal(group_size)) / 2.0
        )
    left = np.concatenate(left_groups)
    right = np.concatenate(right_groups)

    def drift(l_value: np.ndarray, r_value: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        dl = -1j * (-j_hop * r_value + interaction * (np.abs(l_value) ** 2 - 1.0) * l_value)
        dr = -1j * (-j_hop * l_value + interaction * (np.abs(r_value) ** 2 - 1.0) * r_value)
        return dl, dr

    imbalance = np.empty(times.size)
    imbalance_standard_error = np.empty(times.size)
    child_imbalance = np.empty((child_count, times.size))
    total_initial = np.mean(np.abs(left) ** 2 + np.abs(right) ** 2 - 1.0)
    for index, _ in enumerate(times):
        numerator_samples = np.abs(left) ** 2 - np.abs(right) ** 2
        denominator_samples = np.abs(left) ** 2 + np.abs(right) ** 2 - 1.0
        numerator = np.mean(numerator_samples)
        denominator = np.mean(denominator_samples)
        imbalance[index] = numerator / denominator
        influence = (numerator_samples - imbalance[index] * denominator_samples) / denominator
        imbalance_standard_error[index] = np.std(influence, ddof=1) / np.sqrt(trajectories)
        for group in range(child_count):
            start = group * group_size
            stop = start + group_size
            child_imbalance[group, index] = np.mean(numerator_samples[start:stop]) / np.mean(
                denominator_samples[start:stop]
            )
        if index == times.size - 1:
            break
        dt = times[index + 1] - times[index]
        k1l, k1r = drift(left, right)
        k2l, k2r = drift(left + 0.5 * dt * k1l, right + 0.5 * dt * k1r)
        k3l, k3r = drift(left + 0.5 * dt * k2l, right + 0.5 * dt * k2r)
        k4l, k4r = drift(left + dt * k3l, right + dt * k3r)
        left += dt * (k1l + 2.0 * k2l + 2.0 * k3l + k4l) / 6.0
        right += dt * (k1r + 2.0 * k2r + 2.0 * k3r + k4r) / 6.0
    total_final = np.mean(np.abs(left) ** 2 + np.abs(right) ** 2 - 1.0)
    return (
        imbalance,
        imbalance_standard_error,
        float(abs(total_final - total_initial)),
        child_imbalance,
        [list(child.spawn_key) for child in child_sequences],
    )


def run(seed: int = 20260714, trajectories: int = 20_000) -> dict[str, object]:
    run_context = begin_metrics_run(__file__)
    mean_particles = 8.0
    alpha_l = np.sqrt(mean_particles) + 0j
    alpha_r = 0j
    j_hop = 1.0
    interaction = 0.15
    times = np.linspace(0.0, 8.0, 401)

    exact, retained_norm = exact_dynamics(alpha_l, alpha_r, j_hop, interaction, times)
    twa, twa_se, total_drift, child_twa, child_spawn_keys = twa_dynamics(
        alpha_l, alpha_r, j_hop, interaction, times, trajectories, seed
    )
    early = times <= 2.0
    child_early_errors = np.max(np.abs(child_twa[:, early] - exact[None, early]), axis=1)
    child_full_rms_errors = np.sqrt(np.mean((child_twa - exact[None, :]) ** 2, axis=1))
    thresholds = {
        "minimum_retained_initial_norm": 0.99999,
        "max_early_imbalance_error_t_le_2": 0.02,
        "max_full_window_rms_error": 0.05,
        "max_twa_number_absolute_drift": 1.0e-5,
        "max_child_stream_early_error": 0.04,
    }
    checks = {
        "exact_fock_cutoff_retains_initial_state": bool(
            retained_norm >= thresholds["minimum_retained_initial_norm"]
        ),
        "early_window_matches_exact": bool(
            np.max(np.abs(exact[early] - twa[early]))
            <= thresholds["max_early_imbalance_error_t_le_2"]
        ),
        "full_window_rms_below_cap": bool(
            np.sqrt(np.mean((exact - twa) ** 2)) <= thresholds["max_full_window_rms_error"]
        ),
        "twa_number_conserved": bool(
            total_drift <= thresholds["max_twa_number_absolute_drift"]
        ),
        "independent_child_streams_stable_early": bool(
            np.max(child_early_errors) <= thresholds["max_child_stream_early_error"]
        ),
    }
    metrics: dict[str, object] = {
        "seed": seed,
        "child_stream_spawn_keys": child_spawn_keys,
        "trajectories": trajectories,
        "mean_particles": mean_particles,
        "J": j_hop,
        "U": interaction,
        "fock_cutoff_per_mode": 24,
        "retained_initial_norm_before_renormalization": retained_norm,
        "max_early_imbalance_error_t_le_2": float(np.max(np.abs(exact[early] - twa[early]))),
        "rms_full_window_error": float(np.sqrt(np.mean((exact - twa) ** 2))),
        "max_twa_imbalance_standard_error": float(np.max(twa_se)),
        "child_stream_max_early_errors": child_early_errors.tolist(),
        "child_stream_full_window_rms_errors": child_full_rms_errors.tolist(),
        "twa_total_number_absolute_drift": total_drift,
        "uncertainty_note": (
            "TWA imbalance is a ratio of ensemble means.  Its SE uses the paired "
            "delta-method influence function; the plotted band is plus/minus 2 SE."
        ),
        "validation": {
            "thresholds": thresholds,
            "checks": checks,
            "all_passed": bool(all(checks.values())),
        },
    }

    timeseries = {
        "Jt": times.tolist(),
        "exact_imbalance": exact.tolist(),
        "twa_imbalance": twa.tolist(),
        "twa_imbalance_standard_error": twa_se.tolist(),
    }

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.4, 3.7), constrained_layout=True)
    ax.plot(times, exact, color="#C82423", linewidth=2.0, label="Exact")
    ax.plot(times, twa, color="#2878B5", linewidth=1.6, linestyle="--", label="TWA")
    ax.fill_between(
        times,
        twa - 2.0 * twa_se,
        twa + 2.0 * twa_se,
        color="#2878B5",
        alpha=0.18,
        linewidth=0,
        label=r"TWA $\pm2$ SE",
    )
    ax.set_xlabel(r"$Jt$")
    ax.set_ylabel(r"$(\langle n_L\rangle-\langle n_R\rangle)/\langle n_L+n_R\rangle$")
    ax.set_xlim(times[0], times[-1])
    ax.set_ylim(-1.05, 1.05)
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
        raise RuntimeError(f"Dimer validation failed: {failed}")
    return metrics


if __name__ == "__main__":
    print(json.dumps(run(), ensure_ascii=False, indent=2))
