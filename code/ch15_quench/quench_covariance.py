"""Monte Carlo validation of a stable two-component Bogoliubov quench."""

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


def complete_pair_paraunitary(u: float, v: float) -> np.ndarray:
    """Return S for (a_k,a_-k,a_k^dagger,a_-k^dagger)^T = S beta.

    Positive-energy columns use w=(u,-v)^T.  The off-diagonal swap accounts
    for pairing between k and -k, so S=[[U,-V],[-V,U]] for real u,v.
    """

    identity_pair = np.eye(2)
    opposite_momentum = np.array([[0.0, 1.0], [1.0, 0.0]])
    upper = np.hstack((u * identity_pair, -v * opposite_momentum))
    lower = np.hstack((-v * opposite_momentum, u * identity_pair))
    return np.vstack((upper, lower))


def symplectic_inverse(matrix: np.ndarray, metric: np.ndarray) -> np.ndarray:
    """Paraunitary inverse Sigma_z S^dagger Sigma_z."""

    return metric @ matrix.conjugate().T @ metric


def main() -> None:
    run_context = begin_metrics_run(__file__)
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
    occupation_standard_error = np.empty_like(momenta)
    anomalous_standard_error = np.empty_like(momenta)
    metric = np.diag([1.0, 1.0, -1.0, -1.0])
    identity_nambu = np.eye(4)
    identity_pair = np.eye(2)
    opposite_momentum = np.array([[0.0, 1.0], [1.0, 0.0]])
    max_complete_paraunitary_error = 0.0
    max_inverse_formula_error = 0.0
    max_inverse_identity_error = 0.0
    max_basis_transform_matrix_error = 0.0
    max_atomic_sample_reconstruction_error = 0.0
    max_atomic_covariance_basis_reconstruction_error = 0.0
    max_atomic_covariance_sample_absolute_error = 0.0
    max_atomic_covariance_sample_standardized_error = 0.0

    for index in range(len(momenta)):
        initial_basis = complete_pair_paraunitary(ui[index], vi[index])
        final_basis = complete_pair_paraunitary(uf[index], vf[index])
        for basis in (initial_basis, final_basis):
            paraunitary_residual = basis.conjugate().T @ metric @ basis - metric
            max_complete_paraunitary_error = max(
                max_complete_paraunitary_error,
                float(np.max(np.abs(paraunitary_residual))),
            )
            inverse_from_metric = symplectic_inverse(basis, metric)
            inverse_direct = np.linalg.inv(basis)
            max_inverse_formula_error = max(
                max_inverse_formula_error,
                float(np.max(np.abs(inverse_direct - inverse_from_metric))),
            )
            max_inverse_identity_error = max(
                max_inverse_identity_error,
                float(np.max(np.abs(basis @ inverse_from_metric - identity_nambu))),
                float(np.max(np.abs(inverse_from_metric @ basis - identity_nambu))),
            )

        final_inverse = symplectic_inverse(final_basis, metric)
        old_to_new = final_inverse @ initial_basis
        expected_old_to_new = np.block(
            [
                [transform_a[index] * identity_pair, transform_b[index] * opposite_momentum],
                [transform_b[index] * opposite_momentum, transform_a[index] * identity_pair],
            ]
        )
        max_basis_transform_matrix_error = max(
            max_basis_transform_matrix_error,
            float(np.max(np.abs(old_to_new - expected_old_to_new))),
        )

        b_plus = (rng.standard_normal(trajectories) + 1j * rng.standard_normal(trajectories)) / 2
        b_minus = (rng.standard_normal(trajectories) + 1j * rng.standard_normal(trajectories)) / 2
        old_nambu_samples = np.column_stack(
            (b_plus, b_minus, b_plus.conjugate(), b_minus.conjugate())
        )
        atomic_samples = old_nambu_samples @ initial_basis.T
        new_nambu_samples = old_nambu_samples @ old_to_new.T
        reconstructed_atomic_samples = new_nambu_samples @ final_basis.T
        max_atomic_sample_reconstruction_error = max(
            max_atomic_sample_reconstruction_error,
            float(np.max(np.abs(reconstructed_atomic_samples - atomic_samples))),
        )

        initial_mode_covariance = 0.5 * identity_nambu
        exact_atomic_covariance = (
            initial_basis
            @ initial_mode_covariance
            @ initial_basis.conjugate().T
        )
        new_mode_covariance = (
            old_to_new
            @ initial_mode_covariance
            @ old_to_new.conjugate().T
        )
        reconstructed_atomic_covariance = (
            final_basis @ new_mode_covariance @ final_basis.conjugate().T
        )
        max_atomic_covariance_basis_reconstruction_error = max(
            max_atomic_covariance_basis_reconstruction_error,
            float(
                np.max(
                    np.abs(reconstructed_atomic_covariance - exact_atomic_covariance)
                )
            ),
        )

        atomic_outer_samples = (
            atomic_samples[:, :, None] * atomic_samples[:, None, :].conjugate()
        )
        sampled_atomic_covariance = np.mean(atomic_outer_samples, axis=0)
        atomic_covariance_residual = sampled_atomic_covariance - exact_atomic_covariance
        max_atomic_covariance_sample_absolute_error = max(
            max_atomic_covariance_sample_absolute_error,
            float(np.max(np.abs(atomic_covariance_residual))),
        )
        covariance_real_se = np.std(
            atomic_outer_samples.real, axis=0, ddof=1
        ) / np.sqrt(trajectories)
        covariance_imag_se = np.std(
            atomic_outer_samples.imag, axis=0, ddof=1
        ) / np.sqrt(trajectories)
        real_z = np.divide(
            np.abs(atomic_covariance_residual.real),
            covariance_real_se,
            out=np.zeros_like(covariance_real_se),
            where=covariance_real_se > 1.0e-15,
        )
        imag_z = np.divide(
            np.abs(atomic_covariance_residual.imag),
            covariance_imag_se,
            out=np.zeros_like(covariance_imag_se),
            where=covariance_imag_se > 1.0e-15,
        )
        max_atomic_covariance_sample_standardized_error = max(
            max_atomic_covariance_sample_standardized_error,
            float(np.max(real_z)),
            float(np.max(imag_z)),
        )

        c_plus = new_nambu_samples[:, 0]
        c_minus = new_nambu_samples[:, 1]
        occupation_samples = np.abs(c_plus) ** 2 - 0.5
        anomalous_samples = np.real(c_plus * c_minus)
        sampled_occupation[index] = np.mean(occupation_samples)
        sampled_anomalous[index] = np.mean(anomalous_samples)
        occupation_standard_error[index] = np.std(occupation_samples, ddof=1) / np.sqrt(
            trajectories
        )
        anomalous_standard_error[index] = np.std(anomalous_samples, ddof=1) / np.sqrt(
            trajectories
        )

    max_occupation_error = float(np.max(np.abs(sampled_occupation - exact_occupation)))
    max_anomalous_error = float(np.max(np.abs(sampled_anomalous - exact_anomalous)))
    max_occupation_z = float(
        np.max(np.abs(sampled_occupation - exact_occupation) / occupation_standard_error)
    )
    max_anomalous_z = float(
        np.max(np.abs(sampled_anomalous - exact_anomalous) / anomalous_standard_error)
    )
    thresholds = {
        "minimum_final_energy": 1.0e-8,
        "max_symplectic_constraint_error": 1.0e-12,
        "max_occupation_absolute_error": 0.006,
        "max_anomalous_absolute_error": 0.006,
        "max_standardized_error": 5.0,
        "max_complete_matrix_error": 1.0e-12,
        "max_atomic_covariance_absolute_error": 0.008,
    }
    checks = {
        "final_bogoliubov_spectrum_stable": bool(
            np.min(final_energy) >= thresholds["minimum_final_energy"]
        ),
        "bogoliubov_map_symplectic": bool(
            symplectic_error <= thresholds["max_symplectic_constraint_error"]
        ),
        "occupation_absolute_error_below_cap": bool(
            max_occupation_error <= thresholds["max_occupation_absolute_error"]
        ),
        "anomalous_absolute_error_below_cap": bool(
            max_anomalous_error <= thresholds["max_anomalous_absolute_error"]
        ),
        "monte_carlo_residuals_within_5se": bool(
            max(max_occupation_z, max_anomalous_z) <= thresholds["max_standardized_error"]
        ),
        "complete_bdg_basis_paraunitary": bool(
            max_complete_paraunitary_error <= thresholds["max_complete_matrix_error"]
        ),
        "symplectic_inverse_matches_direct_inverse": bool(
            max_inverse_formula_error <= thresholds["max_complete_matrix_error"]
        ),
        "complete_bdg_basis_inverts_to_identity": bool(
            max_inverse_identity_error <= thresholds["max_complete_matrix_error"]
        ),
        "complete_basis_map_matches_pair_coefficients": bool(
            max_basis_transform_matrix_error <= thresholds["max_complete_matrix_error"]
        ),
        "atomic_samples_reconstruct_across_bases": bool(
            max_atomic_sample_reconstruction_error
            <= thresholds["max_complete_matrix_error"]
        ),
        "atomic_covariance_reconstructs_across_bases": bool(
            max_atomic_covariance_basis_reconstruction_error
            <= thresholds["max_complete_matrix_error"]
        ),
        "atomic_covariance_sampling_absolute_error_below_cap": bool(
            max_atomic_covariance_sample_absolute_error
            <= thresholds["max_atomic_covariance_absolute_error"]
        ),
        "atomic_covariance_sampling_residuals_within_5se": bool(
            max_atomic_covariance_sample_standardized_error
            <= thresholds["max_standardized_error"]
        ),
    }
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
        "nambu_convention": "positive column w=(u,-v); S=[[U,-V*],[-V,U*]]",
        "minimum_final_energy": float(final_energy.min()),
        "max_symplectic_constraint_error": symplectic_error,
        "max_complete_paraunitary_error": max_complete_paraunitary_error,
        "max_symplectic_inverse_formula_error": max_inverse_formula_error,
        "max_complete_inverse_identity_error": max_inverse_identity_error,
        "max_complete_basis_transform_matrix_error": max_basis_transform_matrix_error,
        "max_atomic_sample_reconstruction_error": max_atomic_sample_reconstruction_error,
        "max_atomic_covariance_basis_reconstruction_error": (
            max_atomic_covariance_basis_reconstruction_error
        ),
        "max_atomic_covariance_sample_absolute_error": (
            max_atomic_covariance_sample_absolute_error
        ),
        "max_atomic_covariance_sample_standardized_error": (
            max_atomic_covariance_sample_standardized_error
        ),
        "max_occupation_absolute_error": max_occupation_error,
        "max_anomalous_absolute_error": max_anomalous_error,
        "max_occupation_standardized_error": max_occupation_z,
        "max_anomalous_standardized_error": max_anomalous_z,
        "occupation_standard_error": occupation_standard_error.tolist(),
        "anomalous_standard_error": anomalous_standard_error.tolist(),
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

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4))
    axes[0].plot(momenta, exact_occupation, color="#1f77b4", label="exact")
    axes[0].errorbar(
        momenta,
        sampled_occupation,
        yerr=occupation_standard_error,
        fmt="o",
        ms=4,
        mfc="none",
        color="#d62728",
        capsize=2,
        label=r"MC $\pm1$ SE",
    )
    axes[0].set_xlabel(r"$k$")
    axes[0].set_ylabel(r"$\langle c_k^\dagger c_k\rangle$")
    axes[0].grid(alpha=0.2)
    axes[0].legend(frameon=False)
    axes[1].plot(momenta, exact_anomalous, color="#1f77b4", label="exact")
    axes[1].errorbar(
        momenta,
        sampled_anomalous,
        yerr=anomalous_standard_error,
        fmt="o",
        ms=4,
        mfc="none",
        color="#d62728",
        capsize=2,
        label=r"MC $\pm1$ SE",
    )
    axes[1].set_xlabel(r"$k$")
    axes[1].set_ylabel(r"$\mathrm{Re}\langle c_k c_{-k}\rangle$")
    axes[1].grid(alpha=0.2)
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURE)
    plt.close(fig)

    if not metrics["validation"]["all_passed"]:
        failed = [name for name, passed in checks.items() if not passed]
        raise RuntimeError(f"Bogoliubov-quench validation failed: {failed}")
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
