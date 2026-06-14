//! Principal stresses + stress invariants from a 3x3 symmetric Cauchy
//! stress tensor.
//!
//! ## What this is
//!
//! A small, self-contained post-processing helper. The native solvers in
//! this crate ([`crate::native_solver`], [`crate::plasticity`],
//! [`crate::contact`], ...) recover a per-point Cauchy stress in Voigt
//! order `[σxx σyy σzz σxy σyz σzx]` and reduce it to a single scalar
//! von Mises equivalent stress. This module performs the **full principal
//! decomposition** of that same stress state — the three principal
//! stresses, their directions, and the derived invariants engineers read
//! off a stress report:
//!
//! - the three **principal stresses** `σ1 ≥ σ2 ≥ σ3` (the real
//!   eigenvalues of the symmetric stress tensor, sorted descending) and
//!   the orthonormal **principal directions** (the eigenvectors);
//! - the **maximum shear stress** `τ_max = (σ1 − σ3) / 2` — the radius
//!   of the largest Mohr's circle, the quantity a Tresca yield check
//!   compares against;
//! - the **hydrostatic (mean normal) stress** `p = (σ1 + σ2 + σ3) / 3 =
//!   tr(σ) / 3` — the pressure-like part that drives no shear;
//! - the **von Mises** equivalent stress
//!   `σ_vm = √(½·((σ1−σ2)² + (σ2−σ3)² + (σ3−σ1)²))`, computed here from
//!   the principal stresses. For a shear-free (already-principal) state
//!   this is identical to
//!   `native_solver::von_mises_from_voigt`; the two agree by
//!   construction (the off-diagonal terms in the Voigt form are exactly
//!   what the eigen-decomposition rotates away).
//!
//! ## Method
//!
//! The Cauchy stress tensor is real and symmetric, so its eigenvalues
//! are real and its eigenvectors orthogonal. The decomposition uses
//! [`nalgebra::SymmetricEigen`] — the same symmetric eigensolver the
//! modal ([`crate::modal_solver`]) and buckling ([`crate::buckling`])
//! solvers use — then sorts the eigenpairs by eigenvalue, descending, so
//! `σ1` is the algebraically-largest (most tensile) principal stress and
//! `σ3` the smallest (most compressive). No iteration, no tolerance: a
//! single symmetric eigensolve of a fixed 3x3.
//!
//! ## Honest scope
//!
//! This is a **research / preliminary-design-grade** stress-analysis
//! helper: a numerically-robust eigen-decomposition and the classical
//! invariants derived from it, validated below against the textbook
//! closed forms (uniaxial tension, pure shear). It is a pure algebraic
//! post-process of one stress tensor — it does **not** add a failure
//! theory, a fatigue model, stress averaging / extrapolation between
//! integration points, or any solver. It is not, and does not claim to
//! be, the post-processing pipeline of Ansys, Abaqus, or CATIA — those
//! couple principal stresses to calibrated failure criteria, mesh-aware
//! recovery, and full result databases this helper has no part of.

use nalgebra::{Matrix3, SymmetricEigen, Vector3};

/// The principal-stress decomposition of a single symmetric Cauchy
/// stress tensor, together with the derived scalar invariants.
///
/// Construct one with [`PrincipalStress::from_tensor`] (a
/// [`nalgebra::Matrix3`]) or [`PrincipalStress::from_voigt`] (the
/// crate's Voigt order `[σxx σyy σzz σxy σyz σzx]`). All stresses are in
/// the same units as the input (the solvers use pascals).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrincipalStress {
    /// The three principal stresses sorted **descending**: `[σ1, σ2, σ3]`
    /// with `σ1 ≥ σ2 ≥ σ3`. `σ1` is the most tensile (algebraically
    /// largest), `σ3` the most compressive.
    pub principals: [f64; 3],
    /// The unit principal directions, one per principal stress and in the
    /// same order: `directions[i]` is the eigenvector belonging to
    /// `principals[i]`. The three columns form an orthonormal triad.
    pub directions: [Vector3<f64>; 3],
}

impl PrincipalStress {
    /// Decompose a symmetric 3x3 Cauchy stress tensor into its principal
    /// stresses and directions.
    ///
    /// Only the symmetric part of `tensor` is physically meaningful for a
    /// Cauchy stress; the tensor is symmetrised (`½·(σ + σᵀ)`) before the
    /// eigensolve so a caller's round-off asymmetry cannot produce
    /// complex eigenvalues. The eigenpairs are sorted so the principal
    /// stresses come out descending (`σ1 ≥ σ2 ≥ σ3`).
    pub fn from_tensor(tensor: &Matrix3<f64>) -> Self {
        // Symmetrise to guarantee real eigenvalues regardless of any
        // floating-point asymmetry in the caller's tensor.
        let sym = (tensor + tensor.transpose()) * 0.5;
        let eigen = SymmetricEigen::new(sym);

        // Sort the (eigenvalue, eigenvector) pairs by eigenvalue,
        // descending, so index 0 is σ1 (most tensile) and 2 is σ3.
        let mut idx = [0usize, 1, 2];
        idx.sort_by(|&a, &b| {
            eigen.eigenvalues[b]
                .partial_cmp(&eigen.eigenvalues[a])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let principals = [
            eigen.eigenvalues[idx[0]],
            eigen.eigenvalues[idx[1]],
            eigen.eigenvalues[idx[2]],
        ];
        let directions = [
            eigen.eigenvectors.column(idx[0]).into_owned(),
            eigen.eigenvectors.column(idx[1]).into_owned(),
            eigen.eigenvectors.column(idx[2]).into_owned(),
        ];
        PrincipalStress {
            principals,
            directions,
        }
    }

    /// Decompose a Cauchy stress given in the crate's Voigt order
    /// `[σxx σyy σzz σxy σyz σzx]`.
    ///
    /// This is the layout the native solvers store stress in (see
    /// [`crate::plasticity::PlasticState::stress`]); it assembles the
    /// symmetric tensor
    ///
    /// ```text
    ///   | σxx σxy σzx |
    ///   | σxy σyy σyz |
    ///   | σzx σyz σzz |
    /// ```
    ///
    /// and forwards to [`PrincipalStress::from_tensor`].
    pub fn from_voigt(s: &[f64; 6]) -> Self {
        let (sxx, syy, szz, sxy, syz, szx) = (s[0], s[1], s[2], s[3], s[4], s[5]);
        let tensor = Matrix3::new(
            sxx, sxy, szx, //
            sxy, syy, syz, //
            szx, syz, szz,
        );
        Self::from_tensor(&tensor)
    }

    /// The largest (most tensile) principal stress `σ1`.
    #[inline]
    pub fn sigma1(&self) -> f64 {
        self.principals[0]
    }

    /// The intermediate principal stress `σ2`.
    #[inline]
    pub fn sigma2(&self) -> f64 {
        self.principals[1]
    }

    /// The smallest (most compressive) principal stress `σ3`.
    #[inline]
    pub fn sigma3(&self) -> f64 {
        self.principals[2]
    }

    /// The maximum shear stress `τ_max = (σ1 − σ3) / 2` — the radius of
    /// the largest Mohr's circle. This is the quantity a Tresca
    /// (maximum-shear) yield criterion compares against the shear yield
    /// strength.
    #[inline]
    pub fn max_shear(&self) -> f64 {
        (self.principals[0] - self.principals[2]) / 2.0
    }

    /// The hydrostatic (mean normal) stress
    /// `p = (σ1 + σ2 + σ3) / 3 = tr(σ) / 3` — the pressure-like part of
    /// the stress that produces no shear. The principal-sum form used
    /// here equals one third of the tensor trace, since the trace is
    /// invariant under the eigen-rotation.
    #[inline]
    pub fn hydrostatic(&self) -> f64 {
        (self.principals[0] + self.principals[1] + self.principals[2]) / 3.0
    }

    /// The von Mises equivalent stress, computed from the principal
    /// stresses:
    ///
    /// ```text
    ///   σ_vm = √( ½·((σ1−σ2)² + (σ2−σ3)² + (σ3−σ1)²) )
    /// ```
    ///
    /// For a shear-free (already-principal) input this is identical to
    /// `native_solver::von_mises_from_voigt`; for a general
    /// state both forms agree, because the von Mises stress is an
    /// invariant of the deviatoric stress and the principal stresses are
    /// just the stress tensor in its eigenbasis.
    #[inline]
    pub fn von_mises(&self) -> f64 {
        let (s1, s2, s3) = (self.principals[0], self.principals[1], self.principals[2]);
        (0.5 * ((s1 - s2).powi(2) + (s2 - s3).powi(2) + (s3 - s1).powi(2))).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// All three principal directions are mutually orthonormal — they
    /// form a proper triad regardless of the stress state.
    fn assert_orthonormal(ps: &PrincipalStress) {
        for d in &ps.directions {
            assert!(
                (d.norm() - 1.0).abs() < 1e-9,
                "principal direction must be a unit vector, |d| = {}",
                d.norm()
            );
        }
        assert!(ps.directions[0].dot(&ps.directions[1]).abs() < 1e-9);
        assert!(ps.directions[1].dot(&ps.directions[2]).abs() < 1e-9);
        assert!(ps.directions[2].dot(&ps.directions[0]).abs() < 1e-9);
    }

    #[test]
    fn uniaxial_tension_recovers_the_applied_stress() {
        // Closed form: a uniaxial stress diag(σ, 0, 0) has principal
        // stresses (σ, 0, 0), von Mises = σ, and max shear = σ/2.
        let sigma = 250.0e6;
        let ps = PrincipalStress::from_voigt(&[sigma, 0.0, 0.0, 0.0, 0.0, 0.0]);

        assert!((ps.sigma1() - sigma).abs() < 1e-3, "σ1 must be σ");
        assert!(ps.sigma2().abs() < 1e-3, "σ2 must be 0");
        assert!(ps.sigma3().abs() < 1e-3, "σ3 must be 0");
        // Sorted descending.
        assert!(ps.sigma1() >= ps.sigma2() && ps.sigma2() >= ps.sigma3());

        // von Mises of a uniaxial stress is the stress itself.
        assert!(
            (ps.von_mises() - sigma).abs() < 1e-3,
            "uniaxial von Mises must equal σ, got {}",
            ps.von_mises()
        );
        // Max shear is σ/2 (largest Mohr's circle from σ to 0).
        assert!(
            (ps.max_shear() - sigma / 2.0).abs() < 1e-3,
            "uniaxial max shear must be σ/2, got {}",
            ps.max_shear()
        );
        // Hydrostatic stress is σ/3.
        assert!((ps.hydrostatic() - sigma / 3.0).abs() < 1e-3);
        assert_orthonormal(&ps);
    }

    #[test]
    fn pure_shear_has_equal_and_opposite_principals() {
        // Closed form: a pure shear τ on the xy off-diagonal,
        // σ = [[0, τ, 0], [τ, 0, 0], [0, 0, 0]], has principal stresses
        // (τ, 0, −τ) and von Mises = τ·√3.
        let tau = 100.0e6;
        let ps = PrincipalStress::from_voigt(&[0.0, 0.0, 0.0, tau, 0.0, 0.0]);

        assert!((ps.sigma1() - tau).abs() < 1e-3, "σ1 must be +τ");
        assert!(ps.sigma2().abs() < 1e-3, "σ2 must be 0");
        assert!((ps.sigma3() + tau).abs() < 1e-3, "σ3 must be −τ");
        // Sorted descending.
        assert!(ps.sigma1() >= ps.sigma2() && ps.sigma2() >= ps.sigma3());

        // von Mises of pure shear is τ·√3.
        assert!(
            (ps.von_mises() - tau * 3.0_f64.sqrt()).abs() < 1e-2,
            "pure-shear von Mises must be τ·√3, got {}",
            ps.von_mises()
        );
        // Max shear equals the applied shear τ = (τ − (−τ))/2.
        assert!(
            (ps.max_shear() - tau).abs() < 1e-3,
            "pure-shear max shear must be τ, got {}",
            ps.max_shear()
        );
        // Pure shear is deviatoric — zero hydrostatic stress.
        assert!(
            ps.hydrostatic().abs() < 1e-3,
            "pure shear must have zero mean stress, got {}",
            ps.hydrostatic()
        );
        assert_orthonormal(&ps);
    }

    #[test]
    fn von_mises_matches_the_voigt_helper_for_a_general_state() {
        // The principal-stress von Mises must agree with the crate's
        // direct Voigt-form von Mises for an arbitrary stress, since both
        // are the same invariant.
        let s = [3.0e6, -1.0e6, 7.0e6, 2.0e6, -1.5e6, 0.5e6];
        let ps = PrincipalStress::from_voigt(&s);
        let vm_voigt = crate::native_solver::von_mises_from_voigt(&s);
        assert!(
            (ps.von_mises() - vm_voigt).abs() < 1e-2 * vm_voigt.max(1.0),
            "principal von Mises {} must match Voigt von Mises {}",
            ps.von_mises(),
            vm_voigt
        );
    }

    #[test]
    fn hydrostatic_equals_one_third_of_the_trace() {
        // p = tr(σ)/3 is invariant under the eigen-rotation, so the
        // principal-sum form must equal one third of the raw tensor
        // trace for any stress state.
        let s = [5.0e6, -2.0e6, 8.0e6, 1.0e6, 3.0e6, -2.0e6];
        let ps = PrincipalStress::from_voigt(&s);
        let trace_third = (s[0] + s[1] + s[2]) / 3.0;
        assert!(
            (ps.hydrostatic() - trace_third).abs() < 1e-3,
            "hydrostatic {} must equal tr(σ)/3 {}",
            ps.hydrostatic(),
            trace_third
        );
        // The principals must always come out sorted descending.
        assert!(ps.principals[0] >= ps.principals[1]);
        assert!(ps.principals[1] >= ps.principals[2]);
    }

    #[test]
    fn already_diagonal_tensor_is_sorted_descending() {
        // An arbitrary diagonal stress must reorder to σ1 ≥ σ2 ≥ σ3 with
        // the diagonal values as the principals.
        let ps = PrincipalStress::from_voigt(&[-4.0e6, 9.0e6, 2.0e6, 0.0, 0.0, 0.0]);
        assert!((ps.sigma1() - 9.0e6).abs() < 1e-3);
        assert!((ps.sigma2() - 2.0e6).abs() < 1e-3);
        assert!((ps.sigma3() + 4.0e6).abs() < 1e-3);
        assert_orthonormal(&ps);
    }

    #[test]
    fn from_tensor_symmetrises_a_slightly_asymmetric_input() {
        // A tensor with a small asymmetry (round-off) must still give a
        // real, sensible decomposition — the symmetric part is used.
        let tensor = Matrix3::new(
            10.0e6,
            2.0e6,
            0.0, //
            2.0e6 + 1.0,
            5.0e6,
            0.0, // tiny off-diagonal mismatch
            0.0,
            0.0,
            1.0e6,
        );
        let ps = PrincipalStress::from_tensor(&tensor);
        // Eigenvalues are finite and sorted; the sum equals the trace.
        assert!(ps.principals.iter().all(|v| v.is_finite()));
        assert!(ps.principals[0] >= ps.principals[1]);
        assert!(ps.principals[1] >= ps.principals[2]);
        let trace = 10.0e6 + 5.0e6 + 1.0e6;
        let sum = ps.principals[0] + ps.principals[1] + ps.principals[2];
        assert!((sum - trace).abs() < 1.0, "principal sum must equal trace");
    }
}
