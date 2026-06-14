//! **Mode-superposition transient** dynamics — the modal route to a forced
//! time-history response, the companion to the direct Newmark integrator in
//! [`crate::dynamics`].
//!
//! ## What this is
//!
//! Given a linear structural system `M·ü + C·u̇ + K·u = F(t)`, instead of
//! stepping the full coupled system in time (direct integration), this:
//!
//! 1. solves the generalised eigenproblem `K·φ = ω²·M·φ` for the lowest `m`
//!    **mass-normalised** modes (`φᵀ·M·φ = 1`);
//! 2. projects the equations of motion onto those modes — with mass-normalised
//!    shapes and **modal (classical) damping** the system *decouples* into `m`
//!    independent single-DOF oscillators
//!
//!    ```text
//!      q̈ᵢ + 2·ζᵢ·ωᵢ·q̇ᵢ + ωᵢ²·qᵢ = φᵢᵀ·F(t)
//!    ```
//!
//! 3. integrates each modal oscillator (unconditionally-stable Newmark
//!    average-acceleration), and
//! 4. recombines `u(t) = Σ φᵢ·qᵢ(t)`.
//!
//! Keeping only the lowest few modes (the ones a load actually excites) makes
//! this dramatically cheaper than direct integration for large systems — the
//! standard "mode-superposition transient/harmonic" analysis of a production
//! FE code.
//!
//! ## Validated
//!
//! Against closed-form and cross-method results: a single-DOF undamped step
//! response reproduces the analytic `(F/k)(1 − cos ωt)` (dynamic amplification
//! → 2); the modal frequencies of a 2-DOF spring chain match the analytic
//! `ω² = (3 ∓ √5)/2`; a release from a pure mode shape oscillates at exactly
//! that mode's frequency with no other-mode content; and the full undamped
//! forced response agrees with an independent direct-Newmark integration of the
//! coupled system.
//!
//! ## Honest scope
//!
//! Dense `(M, K)` linear systems with classical (modal) damping, Newmark modal
//! integration — research / preliminary-design grade. It operates on assembled
//! matrices (e.g. a reduced model or a small lumped system), not directly on a
//! mesh; pair it with [`crate::assembly`] to build `M`/`K`. Non-classical
//! damping, geometric/material non-linearity, and very large sparse systems are
//! out of scope; a step toward, not an equal of, a production modal solver.

use nalgebra::{DMatrix, DVector};
use thiserror::Error;

/// Errors from the mode-superposition transient solver.
#[derive(Debug, Error)]
pub enum ModalTransientError {
    /// `M` and `K` are not square and the same size.
    #[error("mass and stiffness matrices must be square and the same size")]
    ShapeMismatch,
    /// The mass matrix is not positive-definite (its Cholesky factor failed).
    #[error("mass matrix is not positive-definite (Cholesky failed)")]
    MassNotPositiveDefinite,
    /// The symmetric eigensolve did not converge.
    #[error("symmetric eigensolve failed to converge")]
    EigenFailed,
    /// More modes were requested than the system has degrees of freedom.
    #[error("requested {requested} modes but the system has only {available} DOFs")]
    TooManyModes {
        /// Modes requested.
        requested: usize,
        /// DOFs available.
        available: usize,
    },
    /// A supplied vector had the wrong length for the system.
    #[error("vector length {got} does not match the {expected}-DOF system")]
    VectorLengthMismatch {
        /// Length supplied.
        got: usize,
        /// Length required.
        expected: usize,
    },
    /// The time step was not finite and positive.
    #[error("time step dt must be finite and positive")]
    BadTimeStep,
}

/// The mass-normalised modal basis of a linear system: the lowest angular
/// frequencies `ωᵢ` (ascending) and their mode shapes `φᵢ` solving
/// `K·φ = ω²·M·φ` with `φᵀ·M·φ = 1`.
#[derive(Clone, Debug)]
pub struct ModalBasis {
    /// Angular frequencies `ωᵢ` (rad/s), ascending.
    pub omegas: Vec<f64>,
    /// Mass-normalised mode shapes `φᵢ`, each an `n`-DOF vector.
    pub shapes: Vec<DVector<f64>>,
}

impl ModalBasis {
    /// Number of retained modes.
    pub fn n_modes(&self) -> usize {
        self.omegas.len()
    }
}

/// Compute the lowest `n_modes` mass-normalised eigenpairs of `K·φ = ω²·M·φ`.
///
/// `m` (SPD) and `k` (symmetric PSD) must be the same square size.
///
/// # Errors
///
/// See [`ModalTransientError`].
pub fn modal_basis(
    m: &DMatrix<f64>,
    k: &DMatrix<f64>,
    n_modes: usize,
) -> Result<ModalBasis, ModalTransientError> {
    let n = m.nrows();
    if m.ncols() != n || k.nrows() != n || k.ncols() != n {
        return Err(ModalTransientError::ShapeMismatch);
    }
    if n_modes == 0 || n_modes > n {
        return Err(ModalTransientError::TooManyModes {
            requested: n_modes,
            available: n,
        });
    }
    // Generalised → standard via the Cholesky factor of M: M = L·Lᵀ, the
    // substitution φ = L⁻ᵀ·ψ turns K·φ = λ·M·φ into C·ψ = λ·ψ with
    // C = L⁻¹·K·L⁻ᵀ (symmetric).
    let chol = m
        .clone()
        .cholesky()
        .ok_or(ModalTransientError::MassNotPositiveDefinite)?;
    let l = chol.l();
    let l_inv = l
        .try_inverse()
        .ok_or(ModalTransientError::MassNotPositiveDefinite)?;
    let l_inv_t = l_inv.transpose();
    let mut c = &l_inv * k * &l_inv_t;
    let c_t = c.transpose();
    c = (&c + &c_t) * 0.5; // symmetrise away float noise

    let eigen =
        nalgebra::SymmetricEigen::try_new(c, 1.0e-12, 0).ok_or(ModalTransientError::EigenFailed)?;
    let evals = &eigen.eigenvalues;
    let evecs = &eigen.eigenvectors;

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        evals[a]
            .partial_cmp(&evals[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut omegas = Vec::with_capacity(n_modes);
    let mut shapes = Vec::with_capacity(n_modes);
    for &idx in order.iter().take(n_modes) {
        let lambda = evals[idx].max(0.0); // clamp float-noise negatives
        omegas.push(lambda.sqrt());
        let psi = evecs.column(idx).into_owned();
        // φ = L⁻ᵀ·ψ is mass-normalised because ψ is unit-norm.
        shapes.push(&l_inv_t * psi);
    }
    Ok(ModalBasis { omegas, shapes })
}

/// A mode-superposition transient response: the physical displacement vector at
/// each recorded time level (`displacement[step]` is the `n`-DOF field).
#[derive(Clone, Debug)]
pub struct ModalTransientResponse {
    /// Recorded times (s); `times[0] = 0`.
    pub times: Vec<f64>,
    /// `displacement[step]` — the recombined physical displacement at each time.
    pub displacement: Vec<DVector<f64>>,
}

impl ModalTransientResponse {
    /// The physical displacement at the final recorded time.
    pub fn final_displacement(&self) -> &DVector<f64> {
        self.displacement
            .last()
            .expect("at least the initial frame is always present")
    }

    /// The time history of one degree of freedom.
    pub fn dof_history(&self, dof: usize) -> Vec<f64> {
        self.displacement
            .iter()
            .map(|u| u.get(dof).copied().unwrap_or(f64::NAN))
            .collect()
    }
}

/// Integrate the forced response `M·ü + C·u̇ + K·u = F(t)` by mode
/// superposition over `basis`, recombining to the physical displacement.
///
/// `damping_ratios[i]` is the modal damping ratio `ζᵢ` for mode `i` (the slice
/// may be shorter than the mode count — missing entries default to `0`, i.e.
/// undamped). `force(t)` returns the `n`-DOF physical load at time `t`. `u0` /
/// `v0` are the physical initial displacement / velocity. The history has
/// `n_steps + 1` frames.
///
/// # Errors
///
/// See [`ModalTransientError`].
#[allow(clippy::too_many_arguments)]
pub fn modal_transient_response<F>(
    m: &DMatrix<f64>,
    basis: &ModalBasis,
    damping_ratios: &[f64],
    force: F,
    u0: &DVector<f64>,
    v0: &DVector<f64>,
    dt: f64,
    n_steps: usize,
) -> Result<ModalTransientResponse, ModalTransientError>
where
    F: Fn(f64) -> DVector<f64>,
{
    let n = m.nrows();
    if u0.len() != n {
        return Err(ModalTransientError::VectorLengthMismatch {
            got: u0.len(),
            expected: n,
        });
    }
    if v0.len() != n {
        return Err(ModalTransientError::VectorLengthMismatch {
            got: v0.len(),
            expected: n,
        });
    }
    if !dt.is_finite() || dt <= 0.0 {
        return Err(ModalTransientError::BadTimeStep);
    }
    let nm = basis.n_modes();
    let zeta = |i: usize| damping_ratios.get(i).copied().unwrap_or(0.0);

    // Modal initial conditions: qᵢ(0) = φᵢᵀ·M·u0, q̇ᵢ(0) = φᵢᵀ·M·v0.
    let mu0 = m * u0;
    let mv0 = m * v0;
    let mut q: Vec<f64> = (0..nm).map(|i| basis.shapes[i].dot(&mu0)).collect();
    let mut qd: Vec<f64> = (0..nm).map(|i| basis.shapes[i].dot(&mv0)).collect();

    let modal_force = |f: &DVector<f64>, i: usize| basis.shapes[i].dot(f);
    let recombine = |q: &[f64]| {
        let mut u = DVector::zeros(n);
        for (phi, &qi) in basis.shapes.iter().zip(q) {
            u += phi * qi;
        }
        u
    };

    // Initial modal accelerations from the EOM (unit modal mass).
    let f0 = force(0.0);
    let mut qdd: Vec<f64> = (0..nm)
        .map(|i| {
            let w = basis.omegas[i];
            modal_force(&f0, i) - 2.0 * zeta(i) * w * qd[i] - w * w * q[i]
        })
        .collect();

    let mut times = Vec::with_capacity(n_steps + 1);
    let mut disp = Vec::with_capacity(n_steps + 1);
    times.push(0.0);
    disp.push(recombine(&q));

    // Newmark average-acceleration (β = 1/4, γ = 1/2), unit modal mass.
    let (beta, gamma) = (0.25_f64, 0.5_f64);
    let a0 = 1.0 / (beta * dt * dt);
    let a1 = gamma / (beta * dt);
    let a2 = 1.0 / (beta * dt);
    let a3 = 1.0 / (2.0 * beta) - 1.0;
    let a4 = gamma / beta - 1.0;
    let a5 = 0.5 * dt * (gamma / beta - 2.0);
    let a6 = dt * (1.0 - gamma);
    let a7 = gamma * dt;

    for s in 0..n_steps {
        let t1 = (s as f64 + 1.0) * dt;
        let f1 = force(t1);
        for i in 0..nm {
            let w = basis.omegas[i];
            let kk = w * w; // modal stiffness
            let c = 2.0 * zeta(i) * w; // modal damping (unit modal mass)
            let keff = kk + a0 + a1 * c;
            let rhs = modal_force(&f1, i)
                + (a0 * q[i] + a2 * qd[i] + a3 * qdd[i])
                + c * (a1 * q[i] + a4 * qd[i] + a5 * qdd[i]);
            let q1 = rhs / keff;
            let qdd1 = a0 * (q1 - q[i]) - a2 * qd[i] - a3 * qdd[i];
            let qd1 = qd[i] + a6 * qdd[i] + a7 * qdd1;
            q[i] = q1;
            qd[i] = qd1;
            qdd[i] = qdd1;
        }
        times.push(t1);
        disp.push(recombine(&q));
    }

    Ok(ModalTransientResponse {
        times,
        displacement: disp,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn diag2(a: f64, b: f64) -> DMatrix<f64> {
        DMatrix::from_row_slice(2, 2, &[a, 0.0, 0.0, b])
    }

    #[test]
    fn sdof_undamped_step_matches_analytic() {
        // M ü + K u = F0 (step) from rest → u(t) = (F0/K)(1 − cos ωt),
        // ω = √(K/M); dynamic amplification factor 2 (peak = 2·F0/K).
        let (mass, k, f0) = (2.0, 8.0, 4.0);
        let m = DMatrix::from_row_slice(1, 1, &[mass]);
        let kmat = DMatrix::from_row_slice(1, 1, &[k]);
        let basis = modal_basis(&m, &kmat, 1).unwrap();
        let omega = (k / mass).sqrt();
        assert!((basis.omegas[0] - omega).abs() < 1e-9);
        let dt = 2.0 * PI / omega / 400.0; // 400 steps/period
        let n = 800;
        let resp = modal_transient_response(
            &m,
            &basis,
            &[],
            |_t| DVector::from_row_slice(&[f0]),
            &DVector::from_row_slice(&[0.0]),
            &DVector::from_row_slice(&[0.0]),
            dt,
            n,
        )
        .unwrap();
        let static_disp = f0 / k;
        let mut max_u = 0.0_f64;
        for (step, u) in resp.displacement.iter().enumerate() {
            let t = step as f64 * dt;
            let analytic = static_disp * (1.0 - (omega * t).cos());
            assert!(
                (u[0] - analytic).abs() < 1e-3 * static_disp.max(1.0),
                "t={t}: u {} vs analytic {analytic}",
                u[0]
            );
            max_u = max_u.max(u[0]);
        }
        assert!(
            (max_u - 2.0 * static_disp).abs() < 0.02 * static_disp,
            "peak {max_u} vs 2·static {}",
            2.0 * static_disp
        );
    }

    #[test]
    fn two_dof_chain_modal_frequencies_match_analytic() {
        // Fixed-free 2-mass chain, unit masses, unit springs:
        // K = [[2,−1],[−1,1]], M = I → ω² = (3 ∓ √5)/2.
        let m = diag2(1.0, 1.0);
        let k = DMatrix::from_row_slice(2, 2, &[2.0, -1.0, -1.0, 1.0]);
        let basis = modal_basis(&m, &k, 2).unwrap();
        let w1_sq = (3.0 - 5.0_f64.sqrt()) / 2.0;
        let w2_sq = (3.0 + 5.0_f64.sqrt()) / 2.0;
        assert!(
            (basis.omegas[0] - w1_sq.sqrt()).abs() < 1e-9,
            "ω1 {} vs {}",
            basis.omegas[0],
            w1_sq.sqrt()
        );
        assert!(
            (basis.omegas[1] - w2_sq.sqrt()).abs() < 1e-9,
            "ω2 {} vs {}",
            basis.omegas[1],
            w2_sq.sqrt()
        );
        // Mass-normalised: φᵀ M φ = 1.
        for phi in &basis.shapes {
            let norm = (phi.transpose() * &m * phi)[(0, 0)];
            assert!((norm - 1.0).abs() < 1e-9, "φᵀMφ = {norm} ≠ 1");
        }
    }

    #[test]
    fn release_from_pure_mode_oscillates_at_that_frequency() {
        // Released from exactly mode-1's shape, the system stays in mode 1:
        // u(t) = φ1·cos(ω1 t), with no mode-2 content ever.
        let m = diag2(1.0, 1.0);
        let k = DMatrix::from_row_slice(2, 2, &[2.0, -1.0, -1.0, 1.0]);
        let basis = modal_basis(&m, &k, 2).unwrap();
        let (w1, phi1, phi2) = (
            basis.omegas[0],
            basis.shapes[0].clone(),
            basis.shapes[1].clone(),
        );
        let dt = 2.0 * PI / w1 / 400.0;
        let n = 800;
        let resp = modal_transient_response(
            &m,
            &basis,
            &[],
            |_t| DVector::zeros(2),
            &phi1, // release from mode-1 shape
            &DVector::zeros(2),
            dt,
            n,
        )
        .unwrap();
        for (step, u) in resp.displacement.iter().enumerate() {
            let t = step as f64 * dt;
            // Mode-2 modal coordinate q2 = φ2ᵀ M u must stay ~0.
            let q2 = (phi2.transpose() * &m * u)[(0, 0)];
            assert!(q2.abs() < 1e-3, "t={t}: mode-2 leak q2={q2}");
            // Mode-1 coordinate tracks cos(ω1 t) (amplitude 1, since u0 = φ1).
            let q1 = (phi1.transpose() * &m * u)[(0, 0)];
            assert!(
                (q1 - (w1 * t).cos()).abs() < 2e-3,
                "t={t}: q1 {q1} vs cos(ω1 t) {}",
                (w1 * t).cos()
            );
        }
    }

    #[test]
    fn full_modal_response_agrees_with_direct_newmark() {
        // Cross-method check: the undamped forced response by mode
        // superposition (all modes) must match an independent direct-Newmark
        // integration of the coupled 2-DOF system.
        let m = diag2(1.5, 1.0);
        let k = DMatrix::from_row_slice(2, 2, &[3.0, -1.0, -1.0, 2.0]);
        let basis = modal_basis(&m, &k, 2).unwrap();
        let force = |t: f64| DVector::from_row_slice(&[(0.7 * t).sin(), 0.0]);
        let u0 = DVector::zeros(2);
        let v0 = DVector::zeros(2);
        let dt = 1.0e-3;
        let n = 2000;
        let direct = direct_newmark(&m, &k, &force, dt, n);
        let modal = modal_transient_response(&m, &basis, &[], force, &u0, &v0, dt, n).unwrap();
        for step in (0..=n).step_by(100) {
            let diff = (&modal.displacement[step] - &direct[step]).norm();
            assert!(
                diff < 1e-3,
                "step {step}: modal {:?} vs direct {:?} (Δ={diff})",
                modal.displacement[step].as_slice(),
                direct[step].as_slice()
            );
        }
    }

    #[test]
    fn modal_damping_decays_to_the_static_solution() {
        // A damped SDOF under a step load overshoots then settles to F0/K.
        let (mass, k, f0) = (1.0, 100.0, 5.0);
        let m = DMatrix::from_row_slice(1, 1, &[mass]);
        let kmat = DMatrix::from_row_slice(1, 1, &[k]);
        let basis = modal_basis(&m, &kmat, 1).unwrap();
        let omega = (k / mass).sqrt();
        let dt = 2.0 * PI / omega / 200.0;
        let n = 4000; // many periods so damping settles it
        let resp = modal_transient_response(
            &m,
            &basis,
            &[0.1], // ζ = 0.1
            |_t| DVector::from_row_slice(&[f0]),
            &DVector::from_row_slice(&[0.0]),
            &DVector::from_row_slice(&[0.0]),
            dt,
            n,
        )
        .unwrap();
        let static_disp = f0 / k;
        let peak = resp
            .displacement
            .iter()
            .map(|u| u[0])
            .fold(f64::MIN, f64::max);
        assert!(
            peak > static_disp * 1.2,
            "underdamped step should overshoot"
        );
        let final_u = resp.final_displacement()[0];
        assert!(
            (final_u - static_disp).abs() < 0.02 * static_disp,
            "damped response should settle to F0/K = {static_disp}, got {final_u}"
        );
    }

    /// Reference direct Newmark (average-acceleration) on the coupled,
    /// undamped system `M ü + K u = F(t)` — independent of the modal path.
    fn direct_newmark<F: Fn(f64) -> DVector<f64>>(
        m: &DMatrix<f64>,
        k: &DMatrix<f64>,
        force: &F,
        dt: f64,
        n_steps: usize,
    ) -> Vec<DVector<f64>> {
        let n = m.nrows();
        let (beta, gamma) = (0.25_f64, 0.5_f64);
        let a0 = 1.0 / (beta * dt * dt);
        let a2 = 1.0 / (beta * dt);
        let a3 = 1.0 / (2.0 * beta) - 1.0;
        let a6 = dt * (1.0 - gamma);
        let a7 = gamma * dt;
        let keff = k + m * a0;
        let keff_inv = keff.clone().try_inverse().expect("Keff invertible");
        let mut u = DVector::zeros(n);
        let mut v = DVector::zeros(n);
        // a = M⁻¹ (F − K u)
        let m_inv = m.clone().try_inverse().expect("M invertible");
        let mut a = &m_inv * (force(0.0) - k * &u);
        let mut out = vec![u.clone()];
        for s in 0..n_steps {
            let t1 = (s as f64 + 1.0) * dt;
            let rhs = force(t1) + m * (a0 * &u + a2 * &v + a3 * &a);
            let u1 = &keff_inv * rhs;
            let a1v = a0 * (&u1 - &u) - a2 * &v - a3 * &a;
            let v1 = &v + a6 * &a + a7 * &a1v;
            u = u1;
            v = v1;
            a = a1v;
            out.push(u.clone());
        }
        out
    }
}
