//! Invariants of a 3×3 symmetric **stress** (or strain) tensor.
//!
//! ## What this is
//!
//! A small, self-contained post-processing helper that turns a Cauchy
//! stress tensor `σ` into the scalar quantities yield criteria,
//! pressure-dependent constitutive models, and visualisation colour-maps
//! are written in terms of — the quantities that are *independent of the
//! coordinate frame* the stress was expressed in. Given the symmetric
//! tensor it computes:
//!
//! - the three **principal invariants**
//!   `I₁ = tr σ`, `I₂`, `I₃ = det σ` — the coefficients of the
//!   characteristic polynomial `λ³ − I₁λ² + I₂λ − I₃ = 0` whose roots
//!   are the principal stresses;
//! - the **deviatoric invariants** `J₂` and `J₃` — the invariants of
//!   the stress deviator `s = σ − p·I` (the part of the stress left
//!   after the mean normal stress is removed);
//! - the **hydrostatic / mean stress** `p = I₁ / 3` (the pressure is
//!   `−p`);
//! - the **von Mises equivalent stress** `σ̄ = √(3 J₂)`, the scalar the
//!   `J2` (von Mises) yield surface in [`crate::plasticity`] is built
//!   on.
//!
//! ## Conventions
//!
//! The Voigt constructor [`StressInvariants::from_voigt`] takes the same
//! order the rest of the crate uses — `[σxx σyy σzz σxy σyz σzx]` — so
//! it drops in directly on the stress vectors
//! [`crate::native_solver`] and [`crate::elements`] produce. The von
//! Mises value it returns is identical to the crate's existing
//! `von_mises_from_voigt` / [`crate::elements::von_mises`] helpers; this
//! module adds the *rest* of the invariant set those helpers do not
//! expose (`I₁..I₃`, `J₃`, `p`, principal stresses).
//!
//! All quantities are mechanics sign-convention: tension positive, so a
//! positive `mean_stress` is net tension and the physical pressure is
//! its negation.
//!
//! ## Honest scope
//!
//! This is exact closed-form tensor algebra on a single 3×3 symmetric
//! tensor — there is no approximation and nothing to converge. It is a
//! research / preliminary-design post-processing convenience, not a
//! constitutive library: it does **not** evaluate any yield/failure
//! envelope (Drucker-Prager, Mohr-Coulomb, Hill, …), does not handle
//! anisotropy, and is not tied to any particular material. It is in no
//! way a replacement for the stress post-processing of a commercial
//! suite such as Ansys or Abaqus; those carry full failure-criterion
//! and material libraries this single-tensor helper deliberately omits.

use nalgebra::{Matrix3, Vector3};

/// The frame-independent invariants of a 3×3 symmetric stress tensor.
///
/// Build one with [`StressInvariants::from_tensor`] (from a symmetric
/// `Matrix3<f64>`) or [`StressInvariants::from_voigt`] (from a
/// `[σxx σyy σzz σxy σyz σzx]` Voigt vector). Every field is in the
/// same stress units as the input tensor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StressInvariants {
    /// First principal invariant `I₁ = tr σ = σxx + σyy + σzz`.
    pub i1: f64,
    /// Second principal invariant
    /// `I₂ = σxxσyy + σyyσzz + σzzσxx − σxy² − σyz² − σzx²`
    /// (the sum of the principal 2×2 minors of `σ`).
    pub i2: f64,
    /// Third principal invariant `I₃ = det σ`.
    pub i3: f64,
    /// Second deviatoric invariant `J₂ = ½ s:s`, where `s` is the
    /// stress deviator. Equivalently `J₂ = I₁²/3 − I₂`. Always `≥ 0`.
    pub j2: f64,
    /// Third deviatoric invariant `J₃ = det s` (the determinant of the
    /// stress deviator). Equivalently
    /// `J₃ = 2I₁³/27 − I₁I₂/3 + I₃`.
    pub j3: f64,
    /// Hydrostatic / mean normal stress `p = I₁ / 3`. The physical
    /// pressure is `−p` (tension is positive here).
    pub mean_stress: f64,
    /// von Mises equivalent stress `σ̄ = √(3 J₂)`. Always `≥ 0`.
    pub von_mises: f64,
}

impl StressInvariants {
    /// Compute the invariants of a symmetric stress tensor.
    ///
    /// Only the upper triangle of `sigma` is read, so a tensor that is
    /// symmetric apart from round-off is handled exactly as if it were
    /// perfectly symmetric (the lower triangle is ignored, not
    /// averaged).
    pub fn from_tensor(sigma: &Matrix3<f64>) -> Self {
        let sxx = sigma[(0, 0)];
        let syy = sigma[(1, 1)];
        let szz = sigma[(2, 2)];
        let sxy = sigma[(0, 1)];
        let syz = sigma[(1, 2)];
        let szx = sigma[(0, 2)];

        let i1 = sxx + syy + szz;
        let i2 = sxx * syy + syy * szz + szz * sxx - sxy * sxy - syz * syz - szx * szx;
        // det of the symmetric tensor by cofactor expansion.
        let i3 = sxx * (syy * szz - syz * syz) - sxy * (sxy * szz - syz * szx)
            + szx * (sxy * syz - syy * szx);

        // Deviatoric invariants. J2 is the second invariant of the
        // deviator; the I1²/3 − I2 identity is used to avoid forming the
        // deviator explicitly. It is non-negative by construction; clamp
        // away a tiny negative round-off so √ is always real.
        let j2 = (i1 * i1 / 3.0 - i2).max(0.0);
        // J3 = det(s) via the principal-invariant identity.
        let j3 = 2.0 * i1 * i1 * i1 / 27.0 - i1 * i2 / 3.0 + i3;

        let mean_stress = i1 / 3.0;
        let von_mises = (3.0 * j2).sqrt();

        Self {
            i1,
            i2,
            i3,
            j2,
            j3,
            mean_stress,
            von_mises,
        }
    }

    /// Compute the invariants from a Voigt stress vector in the crate's
    /// order `[σxx σyy σzz σxy σyz σzx]`.
    ///
    /// Convenience wrapper over [`StressInvariants::from_tensor`] for the
    /// stress vectors the native solvers return.
    pub fn from_voigt(s: &[f64; 6]) -> Self {
        let sigma = Matrix3::new(
            s[0], s[3], s[5], // σxx σxy σzx
            s[3], s[1], s[4], // σxy σyy σyz
            s[5], s[4], s[2], // σzx σyz σzz
        );
        Self::from_tensor(&sigma)
    }

    /// The three **principal stresses**, sorted descending
    /// (`σ₁ ≥ σ₂ ≥ σ₃`).
    ///
    /// These are the eigenvalues of the (symmetric) stress tensor — the
    /// normal stresses on the planes where the shear stress vanishes.
    /// They are reconstructed from the stored invariants by the
    /// closed-form trigonometric solution of the characteristic cubic
    /// `λ³ − I₁λ² + I₂λ − I₃ = 0`, which for a symmetric tensor always
    /// has three real roots.
    pub fn principal_stresses(&self) -> Vector3<f64> {
        // Work with the deviatoric cubic to use the stable
        // three-real-roots trigonometric form. With p_m the mean stress
        // and J2, J3 the deviatoric invariants, the deviatoric principal
        // values are 2√(J2/3)·cos(θ + 2πk/3).
        let p_m = self.mean_stress;
        let r = (self.j2 / 3.0).sqrt(); // = √(J2/3)

        if r < 1e-300 {
            // Hydrostatic: a triple root at the mean stress.
            return Vector3::new(p_m, p_m, p_m);
        }

        // cos(3θ) = J3 / (2 r³); clamp to [−1, 1] against round-off.
        let arg = (self.j3 / (2.0 * r * r * r)).clamp(-1.0, 1.0);
        let theta = arg.acos() / 3.0;

        let two_r = 2.0 * r;
        use std::f64::consts::PI;
        let s1 = two_r * theta.cos();
        let s2 = two_r * (theta - 2.0 * PI / 3.0).cos();
        let s3 = two_r * (theta + 2.0 * PI / 3.0).cos();

        // s1 ≥ s3 ≥ s2 from this branch ordering; add the mean back and
        // return sorted descending.
        let mut vals = [p_m + s1, p_m + s2, p_m + s3];
        vals.sort_by(|a, b| b.partial_cmp(a).unwrap());
        Vector3::new(vals[0], vals[1], vals[2])
    }

    /// The **maximum shear stress** `τ_max = (σ₁ − σ₃) / 2` (the Tresca
    /// shear), with `σ₁` the largest and `σ₃` the smallest principal
    /// stress.
    pub fn max_shear(&self) -> f64 {
        let p = self.principal_stresses();
        (p[0] - p[2]) / 2.0
    }
}

/// The von Mises equivalent stress of a Voigt stress vector
/// `[σxx σyy σzz σxy σyz σzx]` — `σ̄ = √(3 J₂)`.
///
/// A free-function shortcut for callers that only need the von Mises
/// value and not the whole [`StressInvariants`] set. Numerically
/// identical to [`StressInvariants::from_voigt`]'s `von_mises` field.
pub fn von_mises(s: &[f64; 6]) -> f64 {
    StressInvariants::from_voigt(s).von_mises
}

#[cfg(test)]
mod tests {
    use super::*;

    // Stresses are written at MPa scale to match the rest of the crate.

    #[test]
    fn uniaxial_i1_equals_sigma() {
        // diag(σ, 0, 0): I1 = tr σ = σ.
        let sigma = 250.0e6;
        let inv = StressInvariants::from_voigt(&[sigma, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(
            (inv.i1 - sigma).abs() < 1.0,
            "uniaxial I1 should equal σ, got {}",
            inv.i1
        );
        // For uniaxial the other principal invariants vanish.
        assert!(inv.i2.abs() < 1.0, "uniaxial I2 should be 0");
        assert!(inv.i3.abs() < 1.0, "uniaxial I3 should be 0");
        // Mean stress is σ/3.
        assert!((inv.mean_stress - sigma / 3.0).abs() < 1.0);
    }

    #[test]
    fn uniaxial_von_mises_equals_sigma() {
        // For uniaxial diag(σ,0,0): J2 = σ²/3, so √(3 J2) = |σ|.
        let sigma = 175.0e6;
        let inv = StressInvariants::from_voigt(&[sigma, 0.0, 0.0, 0.0, 0.0, 0.0]);
        // J2 closed form.
        assert!(
            (inv.j2 - sigma * sigma / 3.0).abs() / (sigma * sigma) < 1e-12,
            "uniaxial J2 should be σ²/3"
        );
        // von Mises from J2 equals σ.
        assert!(
            (inv.von_mises - sigma).abs() < 1.0,
            "uniaxial von Mises should equal σ, got {}",
            inv.von_mises
        );
    }

    #[test]
    fn pure_shear_von_mises_equals_tau_sqrt3() {
        // Pure shear σxy = τ (all else 0): J2 = τ², so √(3 J2) = τ√3.
        let tau = 100.0e6;
        let inv = StressInvariants::from_voigt(&[0.0, 0.0, 0.0, tau, 0.0, 0.0]);
        assert!((inv.i1).abs() < 1.0, "pure shear I1 should be 0");
        // I2 = −τ² for pure shear.
        assert!(
            (inv.i2 + tau * tau).abs() / (tau * tau) < 1e-12,
            "pure-shear I2 should be −τ²"
        );
        // J2 = τ².
        assert!(
            (inv.j2 - tau * tau).abs() / (tau * tau) < 1e-12,
            "pure-shear J2 should be τ²"
        );
        // von Mises = τ·√3.
        let expected = tau * 3.0_f64.sqrt();
        assert!(
            (inv.von_mises - expected).abs() < 1.0,
            "pure-shear von Mises should be τ√3, got {} vs {}",
            inv.von_mises,
            expected
        );
    }

    #[test]
    fn hydrostatic_has_zero_deviatoric_invariants() {
        // Pure hydrostatic stress p·I: deviator is zero, so J2 = J3 =
        // von Mises = 0, while I1 = 3p and mean stress = p.
        let p = 60.0e6;
        let inv = StressInvariants::from_voigt(&[p, p, p, 0.0, 0.0, 0.0]);
        assert!((inv.i1 - 3.0 * p).abs() < 1.0);
        assert!((inv.mean_stress - p).abs() < 1.0);
        assert!(inv.j2.abs() < 1.0, "hydrostatic J2 should be 0");
        assert!(inv.j3.abs() < 1e3, "hydrostatic J3 should be 0");
        assert!(inv.von_mises < 1.0, "hydrostatic von Mises should be 0");
    }

    #[test]
    fn invariants_are_rotation_invariant() {
        // The defining property: I1..I3, J2, J3, von Mises must be
        // unchanged by an orthogonal change of frame σ' = Rᵀ σ R.
        let sigma = Matrix3::new(
            120.0e6, 30.0e6, -20.0e6, //
            30.0e6, -40.0e6, 15.0e6, //
            -20.0e6, 15.0e6, 70.0e6,
        );
        let base = StressInvariants::from_tensor(&sigma);

        // A proper rotation about an arbitrary axis.
        let axis = Vector3::new(1.0, 2.0, -0.5).normalize();
        let rot = nalgebra::Rotation3::from_axis_angle(&nalgebra::Unit::new_normalize(axis), 0.9);
        let r = rot.matrix();
        let rotated = r.transpose() * sigma * r;
        let after = StressInvariants::from_tensor(&rotated);

        let scale = 1.0e6; // absolute tolerance ~1 Pa-equivalent at MPa scale
        assert!((after.i1 - base.i1).abs() < scale * 1e-6);
        assert!((after.i2 - base.i2).abs() / (scale * scale) < 1e-9);
        assert!((after.i3 - base.i3).abs() / (scale * scale * scale) < 1e-9);
        assert!((after.j2 - base.j2).abs() / (scale * scale) < 1e-9);
        assert!((after.von_mises - base.von_mises).abs() < scale * 1e-3);
    }

    #[test]
    fn principal_stresses_recover_a_diagonal_tensor() {
        // For an already-diagonal tensor the principal stresses are just
        // the sorted diagonal entries.
        let inv = StressInvariants::from_voigt(&[30.0e6, 100.0e6, -50.0e6, 0.0, 0.0, 0.0]);
        let p = inv.principal_stresses();
        assert!((p[0] - 100.0e6).abs() < 10.0, "σ1 should be 100 MPa");
        assert!((p[1] - 30.0e6).abs() < 10.0, "σ2 should be 30 MPa");
        assert!((p[2] + 50.0e6).abs() < 10.0, "σ3 should be −50 MPa");
        // I1 must equal the sum of principal stresses.
        assert!((inv.i1 - (p[0] + p[1] + p[2])).abs() < 10.0);
    }

    #[test]
    fn principal_stresses_match_nalgebra_eigenvalues() {
        // The trigonometric cubic root must agree with nalgebra's
        // symmetric eigenvalue solver on a general tensor.
        let sigma = Matrix3::new(
            120.0e6, 30.0e6, -20.0e6, //
            30.0e6, -40.0e6, 15.0e6, //
            -20.0e6, 15.0e6, 70.0e6,
        );
        let inv = StressInvariants::from_tensor(&sigma);
        let ours = inv.principal_stresses();

        let mut eig = sigma.symmetric_eigenvalues().as_slice().to_vec();
        eig.sort_by(|a, b| b.partial_cmp(a).unwrap());

        for k in 0..3 {
            assert!(
                (ours[k] - eig[k]).abs() < 10.0,
                "principal stress {k}: ours {} vs nalgebra {}",
                ours[k],
                eig[k]
            );
        }
    }

    #[test]
    fn pure_shear_max_shear_equals_tau() {
        // Pure shear τ has principal stresses (τ, 0, −τ), so
        // τ_max = (σ1 − σ3)/2 = τ.
        let tau = 80.0e6;
        let inv = StressInvariants::from_voigt(&[0.0, 0.0, 0.0, tau, 0.0, 0.0]);
        assert!(
            (inv.max_shear() - tau).abs() < 10.0,
            "max shear of pure shear should equal τ, got {}",
            inv.max_shear()
        );
    }

    #[test]
    fn j2_identity_matches_explicit_deviator() {
        // Cross-check the I1²/3 − I2 identity against ½ s:s formed
        // explicitly from the deviator.
        let s = [120.0e6, -40.0e6, 70.0e6, 30.0e6, 15.0e6, -20.0e6];
        let inv = StressInvariants::from_voigt(&s);
        let mean = (s[0] + s[1] + s[2]) / 3.0;
        let dev = [s[0] - mean, s[1] - mean, s[2] - mean, s[3], s[4], s[5]];
        // ½ s:s with engineering-shear doubling (each off-diagonal shear
        // appears twice in the tensor double-contraction).
        let j2_explicit = 0.5
            * (dev[0] * dev[0]
                + dev[1] * dev[1]
                + dev[2] * dev[2]
                + 2.0 * (dev[3] * dev[3] + dev[4] * dev[4] + dev[5] * dev[5]));
        assert!(
            (inv.j2 - j2_explicit).abs() / j2_explicit < 1e-12,
            "J2 identity should match explicit ½ s:s"
        );
    }

    #[test]
    fn free_function_matches_struct() {
        let s = [50.0e6, -10.0e6, 20.0e6, 5.0e6, -3.0e6, 8.0e6];
        assert_eq!(von_mises(&s), StressInvariants::from_voigt(&s).von_mises);
    }
}
