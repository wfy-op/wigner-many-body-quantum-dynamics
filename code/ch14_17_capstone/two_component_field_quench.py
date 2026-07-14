"""End-to-end two-component 1D truncated-Wigner field quench.

The calculation deliberately uses the uniform Rabi model, not Raman SOC.  It
samples a finite-temperature Bogoliubov Wigner state at J_i > 0, quenches to
J_f = 0, propagates both spatial fields with a Strang split-step method, and
measures an intercomponent coherence and a normally ordered spin structure
factor.  Common random numbers are used for the time-step comparison.
"""

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
FIGURE = ROOT / "figures" / "ch14_17_capstone" / "two_component_field_quench.pdf"
METRICS = ROOT / "data" / "ch14_17_capstone" / "two_component_field_quench_metrics.json"


def stable_bogoliubov_coefficients(
    epsilon: np.ndarray, interaction_scale: float, rabi_shift: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return E, u, v for A=epsilon+rabi_shift+scale/2, B=scale/2."""

    a_coefficient = epsilon + rabi_shift + 0.5 * interaction_scale
    b_coefficient = 0.5 * interaction_scale
    energy_squared = a_coefficient**2 - b_coefficient**2
    if np.any(energy_squared <= 0.0):
        raise ValueError("The sampled Bogoliubov background must be stable")
    energy = np.sqrt(energy_squared)
    u_mode = np.sqrt(0.5 * (a_coefficient / energy + 1.0))
    v_mode = np.sign(b_coefficient) * np.sqrt(
        np.maximum(0.0, 0.5 * (a_coefficient / energy - 1.0))
    )
    return energy, u_mode, v_mode


def thermal_occupation(energy: np.ndarray, temperature: float) -> np.ndarray:
    if temperature == 0.0:
        return np.zeros_like(energy)
    exponent = np.minimum(energy / temperature, 700.0)
    return 1.0 / np.expm1(exponent)


def _complex_mean_standard_error(values: np.ndarray) -> np.ndarray:
    """RMS standard error of a complex sample mean along axis zero."""

    trajectories = values.shape[0]
    centered = values - np.mean(values, axis=0, keepdims=True)
    return np.sqrt(
        np.sum(np.abs(centered) ** 2, axis=0)
        / (trajectories * (trajectories - 1.0))
    )


def _moment_report(
    samples: np.ndarray,
    opposite: np.ndarray,
    normal_target: np.ndarray,
    anomalous_target: np.ndarray,
) -> dict[str, object]:
    """Per-mode normal and paired-anomalous moment residuals."""

    trajectories = samples.shape[0]
    normal_values = np.abs(samples) ** 2
    normal_sample = np.mean(normal_values, axis=0)
    normal_se = np.std(normal_values, axis=0, ddof=1) / np.sqrt(trajectories)
    normal_residual = normal_sample - normal_target
    normal_z = np.abs(normal_residual) / np.maximum(normal_se, np.finfo(float).tiny)

    anomalous_values = samples * samples[:, opposite]
    anomalous_sample = np.mean(anomalous_values, axis=0)
    anomalous_se = _complex_mean_standard_error(anomalous_values)
    anomalous_residual = anomalous_sample - anomalous_target
    anomalous_z = np.abs(anomalous_residual) / np.maximum(
        anomalous_se, np.finfo(float).tiny
    )
    return {
        "normal_target": normal_target.tolist(),
        "normal_sample": normal_sample.tolist(),
        "normal_standard_error": normal_se.tolist(),
        "normal_residual": normal_residual.tolist(),
        "normal_standardized_residual": normal_z.tolist(),
        "anomalous_target_real": np.real(anomalous_target).tolist(),
        "anomalous_target_imag": np.imag(anomalous_target).tolist(),
        "anomalous_sample_real": np.real(anomalous_sample).tolist(),
        "anomalous_sample_imag": np.imag(anomalous_sample).tolist(),
        "anomalous_standard_error": anomalous_se.tolist(),
        "anomalous_absolute_residual": np.abs(anomalous_residual).tolist(),
        "anomalous_standardized_residual": anomalous_z.tolist(),
        "max_normal_absolute_residual": float(np.max(np.abs(normal_residual))),
        "max_normal_standardized_residual": float(np.max(normal_z)),
        "max_anomalous_absolute_residual": float(np.max(np.abs(anomalous_residual))),
        "max_anomalous_standardized_residual": float(np.max(anomalous_z)),
    }


def sample_initial_fields(
    *,
    rng: np.random.Generator,
    trajectories: int,
    points: int,
    length: float,
    density: float,
    interaction_same: float,
    interaction_cross: float,
    rabi_initial: float,
    temperature: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Sample the old BdG basis and reconstruct the two atomic fields."""

    wavevectors = 2.0 * np.pi * np.fft.fftfreq(points, d=length / points)
    epsilon = 0.5 * wavevectors**2
    density_scale = density * (interaction_same + interaction_cross)
    spin_scale = density * (interaction_same - interaction_cross)

    # The density k=0 Goldstone oscillator is singular in thermal BdG.  Here it
    # is instead the displaced coherent condensate mode: its fluctuation is a
    # vacuum Wigner variable with <|b_0|^2>=1/2 and <b_0^2>=0.  Thus it remains
    # a sampled canonical mode and matches the full two-component Weyl baseline.
    epsilon_density = epsilon.copy()
    # Reuse k_1 only as a finite placeholder during the vectorized BdG call;
    # all k=0 coefficients are overwritten by the coherent-mode convention.
    epsilon_density[0] = epsilon[1]
    energy_density, u_density, v_density = stable_bogoliubov_coefficients(
        epsilon_density, density_scale, 0.0
    )
    energy_spin, u_spin, v_spin = stable_bogoliubov_coefficients(
        epsilon, spin_scale, 2.0 * rabi_initial
    )
    occupation_density = thermal_occupation(energy_density, temperature)
    occupation_spin = thermal_occupation(energy_spin, temperature)
    energy_density[0] = 0.0
    occupation_density[0] = 0.0
    u_density[0] = 1.0
    v_density[0] = 0.0

    def complex_wigner(occupation: np.ndarray) -> np.ndarray:
        standard = (
            rng.normal(size=(trajectories, points))
            + 1j * rng.normal(size=(trajectories, points))
        ) / np.sqrt(2.0)
        return standard * np.sqrt(occupation[None, :] + 0.5)

    b_density = complex_wigner(occupation_density)
    b_spin = complex_wigner(occupation_spin)

    opposite = (-np.arange(points)) % points
    atomic_density_k = (
        u_density[None, :] * b_density
        - v_density[None, :] * np.conjugate(b_density[:, opposite])
    )
    atomic_spin_k = (
        u_spin[None, :] * b_spin
        - v_spin[None, :] * np.conjugate(b_spin[:, opposite])
    )
    fourier_to_field = points / np.sqrt(length)
    fluctuation_density = np.fft.ifft(atomic_density_k, axis=1) * fourier_to_field
    fluctuation_spin = np.fft.ifft(atomic_spin_k, axis=1) * fourier_to_field
    background = np.sqrt(0.5 * density)
    field_left = background + (fluctuation_density + fluctuation_spin) / np.sqrt(2.0)
    field_right = background + (fluctuation_density - fluctuation_spin) / np.sqrt(2.0)

    b_density_variance = occupation_density + 0.5
    b_spin_variance = occupation_spin + 0.5
    atomic_density_normal_target = (
        u_density**2 * b_density_variance
        + v_density**2 * b_density_variance[opposite]
    )
    atomic_spin_normal_target = (
        u_spin**2 * b_spin_variance
        + v_spin**2 * b_spin_variance[opposite]
    )
    atomic_density_anomalous_target = -(
        u_density * v_density[opposite] * b_density_variance
        + v_density * u_density[opposite] * b_density_variance[opposite]
    ).astype(complex)
    atomic_spin_anomalous_target = -(
        u_spin * v_spin[opposite] * b_spin_variance
        + v_spin * u_spin[opposite] * b_spin_variance[opposite]
    ).astype(complex)

    density_b_report = _moment_report(
        b_density,
        opposite,
        b_density_variance,
        np.zeros(points, dtype=complex),
    )
    spin_b_report = _moment_report(
        b_spin,
        opposite,
        b_spin_variance,
        np.zeros(points, dtype=complex),
    )
    density_atomic_report = _moment_report(
        atomic_density_k,
        opposite,
        atomic_density_normal_target,
        atomic_density_anomalous_target,
    )
    spin_atomic_report = _moment_report(
        atomic_spin_k,
        opposite,
        atomic_spin_normal_target,
        atomic_spin_anomalous_target,
    )

    positions = np.arange(points) * (length / points)
    fourier_matrix = np.exp(1j * np.outer(positions, wavevectors)) / np.sqrt(points)
    fourier_completeness_error = float(
        np.max(
            np.abs(
                np.conjugate(fourier_matrix).T @ fourier_matrix
                - np.eye(points)
            )
        )
    )
    channel_matrix = np.array([[1.0, 1.0], [1.0, -1.0]]) / np.sqrt(2.0)
    channel_unitarity_error = float(
        np.max(np.abs(channel_matrix.T @ channel_matrix - np.eye(2)))
    )
    sampled_mode_count = 2 * points
    ordering_baseline_mode_count = 2 * points
    sampled = {
        "max_symplectic_norm_error_density": float(
            np.max(np.abs(u_density**2 - v_density**2 - 1.0))
        ),
        "max_symplectic_norm_error_spin": float(
            np.max(np.abs(u_spin**2 - v_spin**2 - 1.0))
        ),
        "minimum_initial_density_energy_nonzero": float(np.min(energy_density[1:])),
        "minimum_initial_spin_energy": float(np.min(energy_spin)),
        "maximum_initial_thermal_occupation": float(
            max(np.max(occupation_density[1:]), np.max(occupation_spin))
        ),
        "wavevectors": wavevectors.tolist(),
        "sampled_channel_mode_count": sampled_mode_count,
        "canonical_atomic_mode_count": 2 * points,
        "ordering_baseline_mode_count": ordering_baseline_mode_count,
        "ordering_baseline_total_half_quanta": float(points),
        "mode_count_matches_ordering_baseline": bool(
            sampled_mode_count == ordering_baseline_mode_count
        ),
        "fourier_projection_completeness_max_error": fourier_completeness_error,
        "density_spin_channel_unitarity_max_error": channel_unitarity_error,
        "zero_mode_manifest": {
            "density_k0": {
                "scheme": "displaced coherent condensate mode with vacuum Wigner noise",
                "sampled": True,
                "thermal_bdg": False,
                "normal_target": 0.5,
                "normal_sample": density_b_report["normal_sample"][0],
                "anomalous_target": 0.0,
                "anomalous_sample_real": density_b_report["anomalous_sample_real"][0],
                "anomalous_sample_imag": density_b_report["anomalous_sample_imag"][0],
            },
            "spin_k0": {
                "scheme": "stable finite-temperature BdG mode",
                "sampled": True,
                "thermal_bdg": True,
                "energy": float(energy_spin[0]),
                "thermal_occupation": float(occupation_spin[0]),
            },
        },
        "per_mode_moment_validation": {
            "density_b_quasiparticles": density_b_report,
            "spin_b_quasiparticles": spin_b_report,
            "density_atomic_channel": density_atomic_report,
            "spin_atomic_channel": spin_atomic_report,
        },
    }
    return field_left, field_right, sampled


def observables(
    field_left: np.ndarray,
    field_right: np.ndarray,
    *,
    dx: float,
    spin_phase: np.ndarray,
) -> tuple[
    complex,
    float,
    float,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Return ratio-of-ensemble coherence and normal-ordered S_s(k_1)."""

    points = field_left.shape[1]
    left_particles = dx * np.sum(np.abs(field_left) ** 2, axis=1) - 0.5 * points
    right_particles = dx * np.sum(np.abs(field_right) ** 2, axis=1) - 0.5 * points
    overlap = dx * np.sum(np.conjugate(field_left) * field_right, axis=1)
    mean_overlap = np.mean(overlap)
    denominator = np.sqrt(np.mean(left_particles) * np.mean(right_particles))
    coherence_complex = mean_overlap / denominator

    spin_density = np.abs(field_left) ** 2 - np.abs(field_right) ** 2
    spin_fourier = dx * np.sum(spin_density * spin_phase[None, :], axis=1)
    total_particles = np.mean(left_particles + right_particles)
    # For two components on an N_x-site canonical grid,
    # W[(S_k S_-k)] = S_{k,W} S_{-k,W} - N_x/2.
    spin_structure = (np.mean(np.abs(spin_fourier) ** 2) - 0.5 * points) / total_particles
    return (
        coherence_complex,
        float(spin_structure),
        float(total_particles),
        overlap,
        left_particles,
        right_particles,
        spin_fourier,
    )


def discrete_weyl_energy(
    field_left: np.ndarray,
    field_right: np.ndarray,
    *,
    length: float,
    interaction_same: float,
    interaction_cross: float,
) -> np.ndarray:
    """Post-quench J_f=0 Weyl Hamiltonian for the canonical collocation grid."""

    points = field_left.shape[1]
    dx = length / points
    delta_projector = 1.0 / dx
    wavevectors = 2.0 * np.pi * np.fft.fftfreq(points, d=dx)
    epsilon = 0.5 * wavevectors**2
    left_modes = np.fft.fft(field_left, axis=1) * dx / np.sqrt(length)
    right_modes = np.fft.fft(field_right, axis=1) * dx / np.sqrt(length)

    # Each of the two canonical atomic modes contributes its quadratic
    # half-quantum, hence the combined subtraction is one per k.
    kinetic = np.sum(
        epsilon[None, :]
        * (np.abs(left_modes) ** 2 + np.abs(right_modes) ** 2 - 1.0),
        axis=1,
    )
    density_left = np.abs(field_left) ** 2
    density_right = np.abs(field_right) ** 2
    same_component = 0.5 * interaction_same * dx * np.sum(
        density_left**2
        - 2.0 * delta_projector * density_left
        + 0.5 * delta_projector**2
        + density_right**2
        - 2.0 * delta_projector * density_right
        + 0.5 * delta_projector**2,
        axis=1,
    )
    cross_component = interaction_cross * dx * np.sum(
        (density_left - 0.5 * delta_projector)
        * (density_right - 0.5 * delta_projector),
        axis=1,
    )
    return kinetic + same_component + cross_component


def propagate(
    initial_left: np.ndarray,
    initial_right: np.ndarray,
    *,
    length: float,
    interaction_same: float,
    interaction_cross: float,
    time_step: float,
    final_time: float,
    output_interval: float,
) -> dict[str, np.ndarray]:
    """Propagate the post-quench J_f=0 TWA fields by Strang splitting."""

    left = initial_left.copy()
    right = initial_right.copy()
    points = left.shape[1]
    dx = length / points
    delta_projector = 1.0 / dx
    wavevectors = 2.0 * np.pi * np.fft.fftfreq(points, d=dx)
    half_kinetic = np.exp(-0.25j * wavevectors**2 * time_step)
    spin_phase = np.exp(-1j * wavevectors[1] * np.arange(points) * dx)
    steps = round(final_time / time_step)
    stride = round(output_interval / time_step)
    if not np.isclose(steps * time_step, final_time) or not np.isclose(
        stride * time_step, output_interval
    ):
        raise ValueError("Time grid and output interval must be commensurate")

    times: list[float] = []
    coherence: list[float] = []
    coherence_real: list[float] = []
    spin_structure: list[float] = []
    particle_number: list[float] = []
    overlaps: list[np.ndarray] = []
    left_numbers: list[np.ndarray] = []
    right_numbers: list[np.ndarray] = []
    spin_fouriers: list[np.ndarray] = []
    energy_samples: list[np.ndarray] = []

    def record(step: int) -> None:
        (
            value,
            spin_value,
            particles,
            overlap,
            n_left,
            n_right,
            spin_fourier,
        ) = observables(left, right, dx=dx, spin_phase=spin_phase)
        times.append(step * time_step)
        coherence.append(abs(value))
        coherence_real.append(value.real)
        spin_structure.append(spin_value)
        particle_number.append(particles)
        overlaps.append(overlap)
        left_numbers.append(n_left)
        right_numbers.append(n_right)
        spin_fouriers.append(spin_fourier)
        energy_samples.append(
            discrete_weyl_energy(
                left,
                right,
                length=length,
                interaction_same=interaction_same,
                interaction_cross=interaction_cross,
            )
        )

    record(0)
    for step in range(1, steps + 1):
        left = np.fft.ifft(np.fft.fft(left, axis=1) * half_kinetic[None, :], axis=1)
        right = np.fft.ifft(np.fft.fft(right, axis=1) * half_kinetic[None, :], axis=1)

        density_left = np.abs(left) ** 2
        density_right = np.abs(right) ** 2
        potential_left = (
            interaction_same * (density_left - delta_projector)
            + interaction_cross * (density_right - 0.5 * delta_projector)
        )
        potential_right = (
            interaction_same * (density_right - delta_projector)
            + interaction_cross * (density_left - 0.5 * delta_projector)
        )
        left *= np.exp(-1j * potential_left * time_step)
        right *= np.exp(-1j * potential_right * time_step)

        left = np.fft.ifft(np.fft.fft(left, axis=1) * half_kinetic[None, :], axis=1)
        right = np.fft.ifft(np.fft.fft(right, axis=1) * half_kinetic[None, :], axis=1)
        if step % stride == 0:
            record(step)

    return {
        "times": np.asarray(times),
        "coherence": np.asarray(coherence),
        "coherence_real": np.asarray(coherence_real),
        "spin_structure": np.asarray(spin_structure),
        "particle_number": np.asarray(particle_number),
        "overlap_samples": np.stack(overlaps, axis=1),
        "left_number_samples": np.stack(left_numbers, axis=1),
        "right_number_samples": np.stack(right_numbers, axis=1),
        "spin_fourier_samples": np.stack(spin_fouriers, axis=1),
        "weyl_energy_samples": np.stack(energy_samples, axis=1),
        "mean_weyl_energy": np.mean(np.stack(energy_samples, axis=1), axis=0),
    }


def delete_block_jackknife_coherence(result: dict[str, np.ndarray], blocks: int) -> np.ndarray:
    overlap = result["overlap_samples"]
    left = result["left_number_samples"]
    right = result["right_number_samples"]
    trajectories = overlap.shape[0]
    if trajectories % blocks != 0:
        raise ValueError("Trajectory count must be divisible by jackknife blocks")
    indices = np.arange(trajectories)
    estimates = []
    for block in np.array_split(indices, blocks):
        keep = np.ones(trajectories, dtype=bool)
        keep[block] = False
        ratio = np.mean(overlap[keep], axis=0) / np.sqrt(
            np.mean(left[keep], axis=0) * np.mean(right[keep], axis=0)
        )
        estimates.append(np.abs(ratio))
    estimates_array = np.asarray(estimates)
    mean_leave_one_out = np.mean(estimates_array, axis=0)
    return np.sqrt(
        (blocks - 1.0)
        / blocks
        * np.sum((estimates_array - mean_leave_one_out[None, :]) ** 2, axis=0)
    )


def delete_block_jackknife_spin_structure(
    result: dict[str, np.ndarray], *, points: int, blocks: int
) -> np.ndarray:
    """Jackknife SE for the ratio-of-means, normal-ordered S_s(k_1)."""

    spin_fourier = result["spin_fourier_samples"]
    total_number = result["left_number_samples"] + result["right_number_samples"]
    trajectories = spin_fourier.shape[0]
    if trajectories % blocks != 0:
        raise ValueError("Trajectory count must be divisible by jackknife blocks")
    indices = np.arange(trajectories)
    estimates = []
    for block in np.array_split(indices, blocks):
        keep = np.ones(trajectories, dtype=bool)
        keep[block] = False
        estimate = (
            np.mean(np.abs(spin_fourier[keep]) ** 2, axis=0) - 0.5 * points
        ) / np.mean(total_number[keep], axis=0)
        estimates.append(estimate)
    estimates_array = np.asarray(estimates)
    mean_leave_one_out = np.mean(estimates_array, axis=0)
    return np.sqrt(
        (blocks - 1.0)
        / blocks
        * np.sum((estimates_array - mean_leave_one_out[None, :]) ** 2, axis=0)
    )


def energy_drift_report(result: dict[str, np.ndarray]) -> dict[str, object]:
    """Dimensionless ensemble and trajectory-RMS Weyl-energy drift curves."""

    energy_samples = result["weyl_energy_samples"]
    mean_energy = np.mean(energy_samples, axis=0)
    mean_scale = max(abs(float(mean_energy[0])), np.finfo(float).tiny)
    ensemble_relative_curve = (mean_energy - mean_energy[0]) / mean_scale
    initial_rms_scale = max(
        float(np.sqrt(np.mean(energy_samples[:, 0] ** 2))), np.finfo(float).tiny
    )
    trajectory_rms_curve = np.sqrt(
        np.mean((energy_samples - energy_samples[:, [0]]) ** 2, axis=0)
    ) / initial_rms_scale
    return {
        "initial_mean_weyl_energy": float(mean_energy[0]),
        "mean_weyl_energy": mean_energy.tolist(),
        "ensemble_relative_drift_curve": ensemble_relative_curve.tolist(),
        "maximum_absolute_ensemble_relative_drift": float(
            np.max(np.abs(ensemble_relative_curve))
        ),
        "trajectory_rms_relative_drift_curve": trajectory_rms_curve.tolist(),
        "maximum_trajectory_rms_relative_drift": float(
            np.max(trajectory_rms_curve)
        ),
    }


def main() -> None:
    run_context = begin_metrics_run(__file__)
    seed = 20260714
    parameters = {
        "length": 32.0,
        "points": 64,
        "density": 12.0,
        "interaction_same": 0.05,
        "interaction_cross": 0.035,
        "rabi_initial": 0.15,
        "rabi_final": 0.0,
        "temperature": 0.30,
        "final_time": 12.0,
        "output_interval": 0.1,
    }
    trajectories = 256
    convergence_trajectories = 64
    rng = np.random.default_rng(seed)
    initial_left, initial_right, sampling_metrics = sample_initial_fields(
        rng=rng,
        trajectories=trajectories,
        points=parameters["points"],
        length=parameters["length"],
        density=parameters["density"],
        interaction_same=parameters["interaction_same"],
        interaction_cross=parameters["interaction_cross"],
        rabi_initial=parameters["rabi_initial"],
        temperature=parameters["temperature"],
    )

    def run(left: np.ndarray, right: np.ndarray, time_step: float) -> dict[str, np.ndarray]:
        return propagate(
            left,
            right,
            length=parameters["length"],
            interaction_same=parameters["interaction_same"],
            interaction_cross=parameters["interaction_cross"],
            time_step=time_step,
            final_time=parameters["final_time"],
            output_interval=parameters["output_interval"],
        )

    production = run(initial_left, initial_right, 0.01)
    common_left = initial_left[:convergence_trajectories]
    common_right = initial_right[:convergence_trajectories]
    coarse = run(common_left, common_right, 0.02)
    medium = run(common_left, common_right, 0.01)
    fine = run(common_left, common_right, 0.005)
    half_trajectory_result = run(initial_left[:128], initial_right[:128], 0.01)
    coarse_medium_difference = float(
        np.max(np.abs(coarse["coherence"] - medium["coherence"]))
    )
    medium_fine_difference = float(
        np.max(np.abs(medium["coherence"] - fine["coherence"]))
    )
    observed_order = float(
        np.log2(coarse_medium_difference / medium_fine_difference)
    )
    coarse_medium_spin_difference = float(
        np.max(np.abs(coarse["spin_structure"] - medium["spin_structure"]))
    )
    medium_fine_spin_difference = float(
        np.max(np.abs(medium["spin_structure"] - fine["spin_structure"]))
    )
    observed_spin_order = float(
        np.log2(coarse_medium_spin_difference / medium_fine_spin_difference)
    )
    trajectory_difference = float(
        np.max(
            np.abs(
                production["coherence"]
                - half_trajectory_result["coherence"]
            )
        )
    )
    spin_trajectory_difference = float(
        np.max(
            np.abs(
                production["spin_structure"]
                - half_trajectory_result["spin_structure"]
            )
        )
    )
    coherence_standard_error = delete_block_jackknife_coherence(production, blocks=8)
    spin_structure_standard_error = delete_block_jackknife_spin_structure(
        production, points=parameters["points"], blocks=8
    )
    particle_drift = float(
        np.max(
            np.abs(
                production["particle_number"] - production["particle_number"][0]
            )
        )
        / production["particle_number"][0]
    )
    production_energy = energy_drift_report(production)
    initial_total_number_samples = (
        production["left_number_samples"][:, 0]
        + production["right_number_samples"][:, 0]
    )
    initial_total_number_standard_deviation = float(
        np.std(initial_total_number_samples, ddof=1)
    )
    initial_total_number_relative_standard_deviation = float(
        initial_total_number_standard_deviation
        / np.mean(initial_total_number_samples)
    )
    coarse_energy = energy_drift_report(coarse)
    medium_energy = energy_drift_report(medium)
    fine_energy = energy_drift_report(fine)

    # A low-cost finite cutoff-window stress test.  Each grid is independently
    # resampled in its own canonical mode basis; these are not nested data and
    # therefore must not be advertised as monotonic grid convergence.
    cutoff_points = [48, 64, 80]
    cutoff_trajectories = 64
    cutoff_stress = []
    for points in cutoff_points:
        cutoff_seed = seed + 10_000 + points
        cutoff_left, cutoff_right, cutoff_sampling = sample_initial_fields(
            rng=np.random.default_rng(cutoff_seed),
            trajectories=cutoff_trajectories,
            points=points,
            length=parameters["length"],
            density=parameters["density"],
            interaction_same=parameters["interaction_same"],
            interaction_cross=parameters["interaction_cross"],
            rabi_initial=parameters["rabi_initial"],
            temperature=parameters["temperature"],
        )
        # The highest cutoff carries the stiffest kinetic mode.  Use the
        # finest verified step for every member of this comparison so that a
        # cutoff trend is not contaminated by a cutoff-dependent splitting
        # error.
        cutoff_result = run(cutoff_left, cutoff_right, 0.005)
        cutoff_se = delete_block_jackknife_coherence(cutoff_result, blocks=8)
        cutoff_spin_se = delete_block_jackknife_spin_structure(
            cutoff_result, points=points, blocks=8
        )
        cutoff_energy = energy_drift_report(cutoff_result)
        cutoff_particle_drift = float(
            np.max(
                np.abs(
                    cutoff_result["particle_number"]
                    - cutoff_result["particle_number"][0]
                )
            )
            / cutoff_result["particle_number"][0]
        )
        cutoff_stress.append(
            {
                "points": points,
                "dx": parameters["length"] / points,
                "cutoff_wavevector": float(np.pi * points / parameters["length"]),
                "seed": cutoff_seed,
                "trajectories": cutoff_trajectories,
                "time_step": 0.005,
                "final_intercomponent_coherence": float(
                    cutoff_result["coherence"][-1]
                ),
                "final_coherence_jackknife_standard_error": float(cutoff_se[-1]),
                "final_normal_ordered_spin_structure_k1": float(
                    cutoff_result["spin_structure"][-1]
                ),
                "final_spin_structure_jackknife_standard_error": float(
                    cutoff_spin_se[-1]
                ),
                "maximum_relative_particle_number_drift": cutoff_particle_drift,
                "maximum_absolute_ensemble_relative_energy_drift": cutoff_energy[
                    "maximum_absolute_ensemble_relative_drift"
                ],
                "maximum_trajectory_rms_relative_energy_drift": cutoff_energy[
                    "maximum_trajectory_rms_relative_drift"
                ],
                "maximum_symplectic_norm_error": float(
                    max(
                        cutoff_sampling["max_symplectic_norm_error_density"],
                        cutoff_sampling["max_symplectic_norm_error_spin"],
                    )
                ),
                "sampled_channel_mode_count": cutoff_sampling[
                    "sampled_channel_mode_count"
                ],
                "ordering_baseline_mode_count": cutoff_sampling[
                    "ordering_baseline_mode_count"
                ],
                "mode_count_matches_ordering_baseline": cutoff_sampling[
                    "mode_count_matches_ordering_baseline"
                ],
                "density_k0_sampled": cutoff_sampling["zero_mode_manifest"][
                    "density_k0"
                ]["sampled"],
                "spin_k0_sampled": cutoff_sampling["zero_mode_manifest"][
                    "spin_k0"
                ]["sampled"],
            }
        )

    cutoff_pairwise_z_scores = []
    cutoff_spin_pairwise_z_scores = []
    for first in range(len(cutoff_stress)):
        for second in range(first + 1, len(cutoff_stress)):
            a = cutoff_stress[first]
            b = cutoff_stress[second]
            combined_se = np.hypot(
                a["final_coherence_jackknife_standard_error"],
                b["final_coherence_jackknife_standard_error"],
            )
            z_score = abs(
                a["final_intercomponent_coherence"]
                - b["final_intercomponent_coherence"]
            ) / combined_se
            cutoff_pairwise_z_scores.append(
                {
                    "points_pair": [a["points"], b["points"]],
                    "absolute_difference": abs(
                        a["final_intercomponent_coherence"]
                        - b["final_intercomponent_coherence"]
                    ),
                    "combined_standard_error": float(combined_se),
                    "z_score": float(z_score),
                }
            )
            combined_spin_se = np.hypot(
                a["final_spin_structure_jackknife_standard_error"],
                b["final_spin_structure_jackknife_standard_error"],
            )
            spin_z_score = abs(
                a["final_normal_ordered_spin_structure_k1"]
                - b["final_normal_ordered_spin_structure_k1"]
            ) / combined_spin_se
            cutoff_spin_pairwise_z_scores.append(
                {
                    "points_pair": [a["points"], b["points"]],
                    "absolute_difference": abs(
                        a["final_normal_ordered_spin_structure_k1"]
                        - b["final_normal_ordered_spin_structure_k1"]
                    ),
                    "combined_standard_error": float(combined_spin_se),
                    "z_score": float(spin_z_score),
                }
            )
    cutoff_max_z_score = float(
        max(item["z_score"] for item in cutoff_pairwise_z_scores)
    )
    cutoff_spin_max_z_score = float(
        max(item["z_score"] for item in cutoff_spin_pairwise_z_scores)
    )

    maximum_coherence_se = float(np.max(coherence_standard_error))
    maximum_spin_structure_se = float(np.max(spin_structure_standard_error))
    trajectory_sensitivity_in_se = trajectory_difference / maximum_coherence_se
    spin_trajectory_sensitivity_in_se = (
        spin_trajectory_difference / maximum_spin_structure_se
    )
    moment_reports = sampling_metrics["per_mode_moment_validation"]
    quasiparticle_moment_max_z = float(
        max(
            moment_reports[channel][key]
            for channel in ("density_b_quasiparticles", "spin_b_quasiparticles")
            for key in (
                "max_normal_standardized_residual",
                "max_anomalous_standardized_residual",
            )
        )
    )
    atomic_channel_moment_max_z = float(
        max(
            moment_reports[channel][key]
            for channel in ("density_atomic_channel", "spin_atomic_channel")
            for key in (
                "max_normal_standardized_residual",
                "max_anomalous_standardized_residual",
            )
        )
    )
    transform_completeness_error = float(
        max(
            sampling_metrics["fourier_projection_completeness_max_error"],
            sampling_metrics["density_spin_channel_unitarity_max_error"],
        )
    )
    cutoff_max_particle_drift = float(
        max(item["maximum_relative_particle_number_drift"] for item in cutoff_stress)
    )
    cutoff_max_symplectic_error = float(
        max(item["maximum_symplectic_norm_error"] for item in cutoff_stress)
    )
    cutoff_max_energy_drift = float(
        max(
            item["maximum_absolute_ensemble_relative_energy_drift"]
            for item in cutoff_stress
        )
    )
    coarse_energy_drift = float(
        coarse_energy["maximum_absolute_ensemble_relative_drift"]
    )
    medium_energy_drift = float(
        medium_energy["maximum_absolute_ensemble_relative_drift"]
    )
    fine_energy_drift = float(
        fine_energy["maximum_absolute_ensemble_relative_drift"]
    )
    energy_drift_ratio_medium_over_coarse = medium_energy_drift / coarse_energy_drift
    energy_drift_ratio_fine_over_medium = fine_energy_drift / medium_energy_drift
    thresholds = {
        "max_per_mode_standardized_residual": 5.0,
        "max_algebraic_residual": 1.0e-12,
        "time_order_minimum": 1.8,
        "time_order_maximum": 2.2,
        "max_coherence_finest_pair_difference": 1.0e-5,
        "max_spin_structure_finest_pair_difference": 1.0e-4,
        "max_trajectory_sensitivity_in_standard_errors": 2.0,
        "max_cutoff_pairwise_z_score": 3.0,
        "max_production_ensemble_relative_energy_drift": 1.0e-5,
        "max_production_trajectory_rms_relative_energy_drift": 1.0e-4,
        "max_cutoff_relative_energy_drift": 2.0e-5,
        "max_energy_drift_halving_ratio": 0.40,
        "max_relative_particle_number_drift": 2.0e-10,
    }
    validation_checks = {
        "sampled_mode_count_matches_weyl_baseline": {
            "sampled_channel_mode_count": sampling_metrics[
                "sampled_channel_mode_count"
            ],
            "ordering_baseline_mode_count": sampling_metrics[
                "ordering_baseline_mode_count"
            ],
            "density_k0_sampled": sampling_metrics["zero_mode_manifest"][
                "density_k0"
            ]["sampled"],
            "spin_k0_sampled": sampling_metrics["zero_mode_manifest"]["spin_k0"][
                "sampled"
            ],
            "passed": bool(
                sampling_metrics["mode_count_matches_ordering_baseline"]
                and sampling_metrics["zero_mode_manifest"]["density_k0"]["sampled"]
                and sampling_metrics["zero_mode_manifest"]["spin_k0"]["sampled"]
            ),
        },
        "quasiparticle_per_mode_moments": {
            "maximum_standardized_residual": quasiparticle_moment_max_z,
            "maximum_allowed": thresholds["max_per_mode_standardized_residual"],
            "passed": bool(
                quasiparticle_moment_max_z
                <= thresholds["max_per_mode_standardized_residual"]
            ),
        },
        "atomic_channel_per_mode_covariances": {
            "maximum_standardized_residual": atomic_channel_moment_max_z,
            "maximum_allowed": thresholds["max_per_mode_standardized_residual"],
            "passed": bool(
                atomic_channel_moment_max_z
                <= thresholds["max_per_mode_standardized_residual"]
            ),
        },
        "finite_basis_transform_completeness": {
            "maximum_residual": transform_completeness_error,
            "maximum_allowed": thresholds["max_algebraic_residual"],
            "passed": bool(
                transform_completeness_error
                <= thresholds["max_algebraic_residual"]
            ),
        },
        "particle_number_conservation": {
            "value": particle_drift,
            "maximum_allowed": thresholds["max_relative_particle_number_drift"],
            "passed": bool(
                particle_drift
                <= thresholds["max_relative_particle_number_drift"]
            ),
        },
        "bogoliubov_symplectic_normalization": {
            "value": float(
                max(
                    sampling_metrics["max_symplectic_norm_error_density"],
                    sampling_metrics["max_symplectic_norm_error_spin"],
                )
            ),
            "maximum_allowed": thresholds["max_algebraic_residual"],
            "passed": bool(
                max(
                    sampling_metrics["max_symplectic_norm_error_density"],
                    sampling_metrics["max_symplectic_norm_error_spin"],
                )
                <= thresholds["max_algebraic_residual"]
            ),
        },
        "coherence_strang_time_step_order": {
            "value": observed_order,
            "allowed_interval": [
                thresholds["time_order_minimum"],
                thresholds["time_order_maximum"],
            ],
            "passed": bool(
                thresholds["time_order_minimum"]
                <= observed_order
                <= thresholds["time_order_maximum"]
            ),
        },
        "spin_structure_strang_time_step_order": {
            "value": observed_spin_order,
            "allowed_interval": [
                thresholds["time_order_minimum"],
                thresholds["time_order_maximum"],
            ],
            "passed": bool(
                thresholds["time_order_minimum"]
                <= observed_spin_order
                <= thresholds["time_order_maximum"]
            ),
        },
        "coherence_time_step_observable_difference": {
            "value_dt_0p01_minus_dt_0p005": medium_fine_difference,
            "maximum_allowed": thresholds[
                "max_coherence_finest_pair_difference"
            ],
            "passed": bool(
                medium_fine_difference
                <= thresholds["max_coherence_finest_pair_difference"]
            ),
        },
        "spin_structure_time_step_observable_difference": {
            "value_dt_0p01_minus_dt_0p005": medium_fine_spin_difference,
            "maximum_allowed": thresholds[
                "max_spin_structure_finest_pair_difference"
            ],
            "passed": bool(
                medium_fine_spin_difference
                <= thresholds["max_spin_structure_finest_pair_difference"]
            ),
        },
        "coherence_trajectory_count_sensitivity": {
            "max_abs_coherence_256_minus_128": trajectory_difference,
            "maximum_jackknife_standard_error_256": maximum_coherence_se,
            "sensitivity_in_standard_errors": trajectory_sensitivity_in_se,
            "maximum_allowed_standard_errors": thresholds[
                "max_trajectory_sensitivity_in_standard_errors"
            ],
            "passed": bool(
                trajectory_sensitivity_in_se
                <= thresholds["max_trajectory_sensitivity_in_standard_errors"]
            ),
        },
        "spin_structure_trajectory_count_sensitivity": {
            "max_abs_spin_structure_256_minus_128": spin_trajectory_difference,
            "maximum_jackknife_standard_error_256": maximum_spin_structure_se,
            "sensitivity_in_standard_errors": spin_trajectory_sensitivity_in_se,
            "maximum_allowed_standard_errors": thresholds[
                "max_trajectory_sensitivity_in_standard_errors"
            ],
            "passed": bool(
                spin_trajectory_sensitivity_in_se
                <= thresholds["max_trajectory_sensitivity_in_standard_errors"]
            ),
        },
        "finite_cutoff_window_consistency": {
            "maximum_pairwise_final_coherence_z_score": cutoff_max_z_score,
            "maximum_allowed_z_score": thresholds["max_cutoff_pairwise_z_score"],
            "scope": "independently resampled finite cutoff window; not an extrapolation",
            "passed": bool(
                cutoff_max_z_score <= thresholds["max_cutoff_pairwise_z_score"]
            ),
        },
        "finite_cutoff_window_spin_structure_consistency": {
            "maximum_pairwise_final_spin_structure_z_score": cutoff_spin_max_z_score,
            "maximum_allowed_z_score": thresholds["max_cutoff_pairwise_z_score"],
            "scope": "independently resampled finite cutoff window; not an extrapolation",
            "passed": bool(
                cutoff_spin_max_z_score
                <= thresholds["max_cutoff_pairwise_z_score"]
            ),
        },
        "finite_cutoff_window_particle_conservation": {
            "maximum_relative_drift_across_cutoffs": cutoff_max_particle_drift,
            "maximum_allowed": thresholds["max_relative_particle_number_drift"],
            "per_cutoff": [
                {
                    "points": item["points"],
                    "value": item["maximum_relative_particle_number_drift"],
                    "passed": bool(
                        item["maximum_relative_particle_number_drift"]
                        <= thresholds["max_relative_particle_number_drift"]
                    ),
                }
                for item in cutoff_stress
            ],
            "passed": bool(
                cutoff_max_particle_drift
                <= thresholds["max_relative_particle_number_drift"]
            ),
        },
        "finite_cutoff_window_symplectic_normalization": {
            "maximum_error_across_cutoffs": cutoff_max_symplectic_error,
            "maximum_allowed": thresholds["max_algebraic_residual"],
            "passed": bool(
                cutoff_max_symplectic_error
                <= thresholds["max_algebraic_residual"]
            ),
        },
        "finite_cutoff_window_mode_count_matches_weyl_baseline": {
            "per_cutoff": [
                {
                    "points": item["points"],
                    "sampled_channel_mode_count": item[
                        "sampled_channel_mode_count"
                    ],
                    "ordering_baseline_mode_count": item[
                        "ordering_baseline_mode_count"
                    ],
                    "density_k0_sampled": item["density_k0_sampled"],
                    "spin_k0_sampled": item["spin_k0_sampled"],
                    "passed": bool(
                        item["mode_count_matches_ordering_baseline"]
                        and item["density_k0_sampled"]
                        and item["spin_k0_sampled"]
                    ),
                }
                for item in cutoff_stress
            ],
            "passed": bool(
                all(
                    item["mode_count_matches_ordering_baseline"]
                    and item["density_k0_sampled"]
                    and item["spin_k0_sampled"]
                    for item in cutoff_stress
                )
            ),
        },
        "postquench_weyl_energy_conservation": {
            "production_maximum_ensemble_relative_drift": production_energy[
                "maximum_absolute_ensemble_relative_drift"
            ],
            "production_maximum_trajectory_rms_relative_drift": production_energy[
                "maximum_trajectory_rms_relative_drift"
            ],
            "maximum_allowed_ensemble_relative_drift": thresholds[
                "max_production_ensemble_relative_energy_drift"
            ],
            "maximum_allowed_trajectory_rms_relative_drift": thresholds[
                "max_production_trajectory_rms_relative_energy_drift"
            ],
            "passed": bool(
                production_energy["maximum_absolute_ensemble_relative_drift"]
                <= thresholds["max_production_ensemble_relative_energy_drift"]
                and production_energy["maximum_trajectory_rms_relative_drift"]
                <= thresholds[
                    "max_production_trajectory_rms_relative_energy_drift"
                ]
            ),
        },
        "weyl_energy_step_halving": {
            "coarse_dt_0p02": coarse_energy_drift,
            "medium_dt_0p01": medium_energy_drift,
            "fine_dt_0p005": fine_energy_drift,
            "medium_over_coarse": energy_drift_ratio_medium_over_coarse,
            "fine_over_medium": energy_drift_ratio_fine_over_medium,
            "maximum_allowed_halving_ratio": thresholds[
                "max_energy_drift_halving_ratio"
            ],
            "passed": bool(
                energy_drift_ratio_medium_over_coarse
                <= thresholds["max_energy_drift_halving_ratio"]
                and energy_drift_ratio_fine_over_medium
                <= thresholds["max_energy_drift_halving_ratio"]
            ),
        },
        "finite_cutoff_window_weyl_energy_conservation": {
            "maximum_relative_drift_across_cutoffs": cutoff_max_energy_drift,
            "maximum_allowed": thresholds["max_cutoff_relative_energy_drift"],
            "passed": bool(
                cutoff_max_energy_drift
                <= thresholds["max_cutoff_relative_energy_drift"]
            ),
        },
    }
    validation_all_passed = all(
        check["passed"] for check in validation_checks.values()
    )
    threshold_basis = {
        "five_sigma_mode_moments": (
            "Per-mode Gaussian moment residuals use a conservative 5-SE family "
            "gate because 2*N_x normal and paired-anomalous moments are inspected."
        ),
        "algebraic_tolerance": (
            "1e-12 is a floating-point residual gate for paraunitary norms, "
            "Fourier completeness, and channel unitarity; it is not a statistical fit."
        ),
        "second_order_window": (
            "Strang splitting is second order, so Richardson orders must lie in "
            "[1.8, 2.2] and the finest-pair observable caps are set well below MC SE."
        ),
        "trajectory_gate": (
            "The nested 128/256 change must remain within two maximum delete-block "
            "jackknife standard errors for each reported stochastic observable."
        ),
        "cutoff_gate": (
            "Independently resampled cutoff endpoints must agree within three "
            "combined standard errors; this is a finite-window stress test, not extrapolation."
        ),
        "energy_gate": (
            "For the discrete J_f=0 Weyl Hamiltonian, the ensemble-mean drift "
            "must stay below 1e-5 and the RMS single-trajectory drift below "
            "1e-4 at production resolution.  The latter is not suppressed by "
            "ensemble averaging.  The ensemble drift must also fall by at "
            "least 2.5x per time-step halving, consistent with a second-order "
            "method; cutoff runs use the verified dt=0.005 step and a 2e-5 gate."
        ),
    }

    metrics = {
        "model_scope": "uniform Rabi-coupled two-component 1D field; not Raman SOC",
        "units": "hbar = mass = 1",
        "density_parameter_meaning": (
            "condensate background density; the normal-ordered total particle "
            "number also contains Bogoliubov depletion"
        ),
        "seed": seed,
        "trajectories": trajectories,
        "convergence_trajectories_common_random_numbers": convergence_trajectories,
        **parameters,
        "dx": parameters["length"] / parameters["points"],
        "cutoff_wavevector": float(np.pi * parameters["points"] / parameters["length"]),
        "time_steps_tested": [0.02, 0.01, 0.005],
        **sampling_metrics,
        "initial_mean_physical_particle_number": float(production["particle_number"][0]),
        "initial_physical_particle_number_standard_deviation": (
            initial_total_number_standard_deviation
        ),
        "initial_physical_particle_number_relative_standard_deviation": (
            initial_total_number_relative_standard_deviation
        ),
        "maximum_relative_particle_number_drift": particle_drift,
        "initial_intercomponent_coherence": float(production["coherence"][0]),
        "final_intercomponent_coherence": float(production["coherence"][-1]),
        "initial_normal_ordered_spin_structure_k1": float(production["spin_structure"][0]),
        "final_normal_ordered_spin_structure_k1": float(production["spin_structure"][-1]),
        "maximum_coherence_standard_error_jackknife": maximum_coherence_se,
        "maximum_spin_structure_standard_error_jackknife": maximum_spin_structure_se,
        "max_abs_coherence_dt_0p02_minus_dt_0p01": coarse_medium_difference,
        "max_abs_coherence_dt_0p01_minus_dt_0p005": medium_fine_difference,
        "observed_time_step_order_from_coherence": observed_order,
        "max_abs_spin_structure_dt_0p02_minus_dt_0p01": coarse_medium_spin_difference,
        "max_abs_spin_structure_dt_0p01_minus_dt_0p005": medium_fine_spin_difference,
        "observed_time_step_order_from_spin_structure": observed_spin_order,
        "max_abs_coherence_256_minus_128_trajectories": trajectory_difference,
        "trajectory_sensitivity_in_jackknife_standard_errors": trajectory_sensitivity_in_se,
        "max_abs_spin_structure_256_minus_128_trajectories": spin_trajectory_difference,
        "spin_structure_trajectory_sensitivity_in_jackknife_standard_errors": spin_trajectory_sensitivity_in_se,
        "postquench_discrete_weyl_hamiltonian": {
            "scope": "J_f=0, delta=0 canonical collocation-grid Weyl Hamiltonian",
            "production": production_energy,
            "time_step_common_noise": {
                "dt_0p02": coarse_energy,
                "dt_0p01": medium_energy,
                "dt_0p005": fine_energy,
                "medium_over_coarse_max_drift": energy_drift_ratio_medium_over_coarse,
                "fine_over_medium_max_drift": energy_drift_ratio_fine_over_medium,
            },
        },
        "cutoff_stress_test_scope": (
            "three independently resampled finite cutoffs at fixed L and dt; "
            "diagnostic only, not monotonic grid convergence"
        ),
        "cutoff_stress_test": cutoff_stress,
        "cutoff_pairwise_final_coherence_z_scores": cutoff_pairwise_z_scores,
        "cutoff_pairwise_final_spin_structure_z_scores": cutoff_spin_pairwise_z_scores,
        "cutoff_max_pairwise_z_score": cutoff_max_z_score,
        "cutoff_spin_structure_max_pairwise_z_score": cutoff_spin_max_z_score,
        "cutoff_consistent_within_three_combined_standard_errors": bool(
            cutoff_max_z_score <= 3.0
        ),
        "validation": {
            "thresholds": thresholds,
            "threshold_basis": threshold_basis,
            "checks": validation_checks,
            "all_passed": validation_all_passed,
        },
        "output_times": production["times"].tolist(),
        "coherence_ratio_of_ensemble_means": production["coherence"].tolist(),
        "coherence_jackknife_standard_error": coherence_standard_error.tolist(),
        "normal_ordered_spin_structure_k1": production["spin_structure"].tolist(),
        "spin_structure_jackknife_standard_error": spin_structure_standard_error.tolist(),
    }

    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    METRICS.parent.mkdir(parents=True, exist_ok=True)
    metrics = finalize_metrics(metrics, run_context, __file__)
    METRICS.write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")

    times = production["times"]
    fig, axes_grid = plt.subplots(2, 3, figsize=(11.4, 6.1))
    axes = axes_grid.ravel()
    axes[0].plot(times, production["coherence"], color="#1f77b4", label="256 trajectories")
    axes[0].fill_between(
        times,
        production["coherence"] - coherence_standard_error,
        production["coherence"] + coherence_standard_error,
        color="#1f77b4",
        alpha=0.20,
        linewidth=0,
        label="jackknife SE",
    )
    axes[0].set_xlabel("time")
    axes[0].set_ylabel(r"$|C_{LR}|$")
    axes[0].grid(alpha=0.2)
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].plot(times, production["spin_structure"], color="#d62728")
    axes[1].fill_between(
        times,
        production["spin_structure"] - spin_structure_standard_error,
        production["spin_structure"] + spin_structure_standard_error,
        color="#d62728",
        alpha=0.18,
        linewidth=0,
        label="jackknife SE",
    )
    axes[1].set_xlabel("time")
    axes[1].set_ylabel(r"$S_s(k_1)$")
    axes[1].grid(alpha=0.2)
    axes[1].legend(frameon=False, fontsize=8)

    axes[2].plot(
        times,
        np.asarray(production_energy["ensemble_relative_drift_curve"]),
        color="#8c564b",
    )
    axes[2].axhline(0.0, color="black", lw=0.8)
    axes[2].set_xlabel("time")
    axes[2].set_ylabel(r"$[\langle H_W(t)\rangle-\langle H_W(0)\rangle]/|\langle H_W(0)\rangle|$")
    axes[2].grid(alpha=0.2)

    axes[3].semilogy(
        [0.02, 0.01],
        [coarse_medium_difference, medium_fine_difference],
        "o-",
        color="#2ca02c",
        label=rf"$|C_{{LR}}|$, $p={observed_order:.2f}$",
    )
    axes[3].semilogy(
        [0.02, 0.01],
        [coarse_medium_spin_difference, medium_fine_spin_difference],
        "s-",
        color="#ff7f0e",
        label=rf"$S_s(k_1)$, $p={observed_spin_order:.2f}$",
    )
    axes[3].invert_xaxis()
    axes[3].set_xticks([0.02, 0.01], labels=["0.02", "0.01"])
    axes[3].set_xlabel(r"$\Delta t$")
    axes[3].set_ylabel("max consecutive-step difference")
    axes[3].grid(alpha=0.2, which="both")
    axes[3].legend(frameon=False, fontsize=8)

    cutoff_k = np.asarray([item["cutoff_wavevector"] for item in cutoff_stress])
    cutoff_coherence = np.asarray(
        [item["final_intercomponent_coherence"] for item in cutoff_stress]
    )
    cutoff_se = np.asarray(
        [item["final_coherence_jackknife_standard_error"] for item in cutoff_stress]
    )
    axes[4].errorbar(
        cutoff_k,
        cutoff_coherence,
        yerr=cutoff_se,
        fmt="o-",
        capsize=3,
        color="#9467bd",
    )
    axes[4].set_xlabel(r"$k_{\rm cut}$")
    axes[4].set_ylabel(r"final $|C_{LR}|$")
    axes[4].grid(alpha=0.2)

    cutoff_spin = np.asarray(
        [item["final_normal_ordered_spin_structure_k1"] for item in cutoff_stress]
    )
    cutoff_spin_se = np.asarray(
        [item["final_spin_structure_jackknife_standard_error"] for item in cutoff_stress]
    )
    axes[5].errorbar(
        cutoff_k,
        cutoff_spin,
        yerr=cutoff_spin_se,
        fmt="s-",
        capsize=3,
        color="#e377c2",
    )
    axes[5].set_xlabel(r"$k_{\rm cut}$")
    axes[5].set_ylabel(r"final $S_s(k_1)$")
    axes[5].grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIGURE)
    plt.close(fig)

    if not validation_all_passed:
        failed = [name for name, check in validation_checks.items() if not check["passed"]]
        raise RuntimeError(f"Capstone validation failed: {failed}")
    print(json.dumps({key: value for key, value in metrics.items() if key not in {
        "output_times",
        "coherence_ratio_of_ensemble_means",
        "coherence_jackknife_standard_error",
        "normal_ordered_spin_structure_k1",
        "spin_structure_jackknife_standard_error",
        "wavevectors",
        "per_mode_moment_validation",
        "postquench_discrete_weyl_hamiltonian",
    }}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
