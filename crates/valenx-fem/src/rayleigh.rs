//! Rayleigh proportional damping coefficients for transient dynamics.
//!
//! ## What this is
//!
//! A small, self-contained helper that fits the two **Rayleigh
//! proportional damping** coefficients `alpha` and `beta` to a pair of
//! measured (or specified) modal damping targets. Implicit transient
//! structural solvers — including [`crate::dynamics`]'s Newmark-β
//! integrator — need a damping matrix `C`, but a full physical damping
//! matrix is rarely known. The standard engineering shortcut is to make
//! `C` a linear combination of the mass and stiffness matrices the
//! solver already has:
//!
//! ```text
//!   C = alpha * M + beta * K
//! ```
//!
//! This *proportional* (Rayleigh / Caughey) form is convenient because
//! it is diagonalised by the same mode shapes as `M` and `K`, so each
//! natural mode `i` (angular frequency `omega_i`) ends up with a modal
//! damping ratio that depends only on `omega_i`:
//!
//! ```text
//!   zeta(omega) = 0.5 * (alpha / omega + beta * omega)
//! ```
//!
//! The mass term (`alpha`) dominates at low frequency and the stiffness
//! term (`beta`) dominates at high frequency — the classic Rayleigh
//! "valley" shape, with a minimum near `omega = sqrt(alpha / beta)`.
//!
//! ## Fitting the coefficients
//!
//! With one frequency you can only pin one coefficient; the usual choice
//! is to match a target damping ratio at **two** frequencies, which
//! determines both. Writing `omega_1 = 2*pi*f_1`, `omega_2 = 2*pi*f_2`,
//! the two conditions `zeta(omega_1) = z_1` and `zeta(omega_2) = z_2`
//! form the linear `2x2` system
//!
//! ```text
//!   (1 / omega_1) * alpha + omega_1 * beta = 2 * z_1
//!   (1 / omega_2) * alpha + omega_2 * beta = 2 * z_2
//! ```
//!
//! whose closed-form solution (used by [`RayleighDamping::from_two_modes`])
//! is
//!
//! ```text
//!   beta  = 2 * (z_2 * omega_2 - z_1 * omega_1) / (omega_2^2 - omega_1^2)
//!   alpha = 2 * omega_1 * omega_2 * (z_1 * omega_2 - z_2 * omega_1)
//!                                  / (omega_2^2 - omega_1^2)
//! ```
//!
//! Between the two fit frequencies the actual modal damping dips *below*
//! the target ratio; outside the pair it rises above it. That is the
//! intrinsic behaviour of two-point Rayleigh damping, not an error.
//!
//! ## Honest scope — what it is and is not
//!
//! This is a **preliminary-design / research-grade** coefficient fit, not
//! a damping model in its own right. It is the exact, textbook two-point
//! Rayleigh fit and nothing more:
//!
//! - It only produces the scalar pair `(alpha, beta)`; assembling and
//!   applying `C = alpha*M + beta*K` inside a time integrator is the
//!   caller's job (this crate keeps the coefficients separate so the
//!   same numbers can drive any integrator).
//! - Proportional damping is a *modelling convenience*. Real structures
//!   rarely have damping that is genuinely proportional to mass and
//!   stiffness, and matching only two frequencies leaves every other
//!   mode under- or over-damped by the curve above. For wideband or
//!   strongly non-proportional damping a Caughey series or direct modal
//!   damping is more faithful.
//! - It makes no claim of parity with the damping models in commercial
//!   transient/MBD solvers (Ansys, Adams, Nastran). It is the same
//!   formula those tools expose for Rayleigh damping, surfaced here for
//!   the native dynamics workbench.

use std::f64::consts::PI;

use thiserror::Error;

/// Errors from fitting Rayleigh damping coefficients.
#[derive(Debug, Error, PartialEq)]
pub enum RayleighError {
    /// A supplied frequency was not finite and strictly positive. The
    /// fit divides by `omega = 2*pi*f`, so a zero, negative, or
    /// non-finite frequency is meaningless.
    #[error("frequency must be finite and positive, got {0} Hz")]
    BadFrequency(f64),
    /// A supplied damping ratio was not finite. (Negative ratios are
    /// permitted — they correspond to negative coefficients — but `NaN`
    /// or infinite ratios are rejected.)
    #[error("damping ratio must be finite, got {0}")]
    BadRatio(f64),
    /// The two target frequencies were equal (to within a tiny
    /// tolerance). Two distinct frequencies are required to determine
    /// both coefficients; a single frequency leaves the system rank-1.
    #[error("the two target frequencies must be distinct, got {0} Hz twice")]
    DuplicateFrequency(f64),
}

/// Rayleigh proportional damping coefficients for `C = alpha*M + beta*K`.
///
/// `alpha` multiplies the mass matrix (it damps low-frequency / rigid-ish
/// motion) and `beta` multiplies the stiffness matrix (it damps
/// high-frequency motion). Build one with
/// [`RayleighDamping::from_two_modes`] from two `(frequency, damping
/// ratio)` targets, or construct it directly from known coefficients.
///
/// The resulting modal damping ratio at any angular frequency is
/// [`RayleighDamping::ratio_at_omega`] (or
/// [`RayleighDamping::ratio_at_hz`] for cyclic frequency in hertz).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayleighDamping {
    /// Mass-proportional coefficient `alpha` in `C = alpha*M + beta*K`.
    /// Units: 1/second. Dominates the damping ratio at low frequency.
    pub alpha: f64,
    /// Stiffness-proportional coefficient `beta` in `C = alpha*M + beta*K`.
    /// Units: second. Dominates the damping ratio at high frequency.
    pub beta: f64,
}

impl RayleighDamping {
    /// Construct directly from known coefficients `alpha` (mass term)
    /// and `beta` (stiffness term).
    ///
    /// No fitting is performed; this is the plain constructor for when
    /// the coefficients are already known (for example, copied from a
    /// reference model). Use [`RayleighDamping::from_two_modes`] to fit
    /// them to damping-ratio targets instead.
    pub fn new(alpha: f64, beta: f64) -> Self {
        Self { alpha, beta }
    }

    /// Fit `alpha` and `beta` so the modal damping ratio equals `z1` at
    /// frequency `f1` (hertz) and `z2` at frequency `f2` (hertz).
    ///
    /// Solves the `2x2` linear system described in the [module
    /// docs](self) in closed form. The two frequencies must be finite,
    /// strictly positive, and distinct; the ratios must be finite.
    ///
    /// # Errors
    ///
    /// Returns [`RayleighError::BadFrequency`] if either frequency is not
    /// finite and positive, [`RayleighError::BadRatio`] if either ratio
    /// is not finite, or [`RayleighError::DuplicateFrequency`] if the two
    /// frequencies coincide.
    pub fn from_two_modes(f1: f64, z1: f64, f2: f64, z2: f64) -> Result<Self, RayleighError> {
        if !f1.is_finite() || f1 <= 0.0 {
            return Err(RayleighError::BadFrequency(f1));
        }
        if !f2.is_finite() || f2 <= 0.0 {
            return Err(RayleighError::BadFrequency(f2));
        }
        if !z1.is_finite() {
            return Err(RayleighError::BadRatio(z1));
        }
        if !z2.is_finite() {
            return Err(RayleighError::BadRatio(z2));
        }

        let w1 = 2.0 * PI * f1;
        let w2 = 2.0 * PI * f2;

        // Distinctness is judged on the angular frequencies, scaled by
        // their magnitude so the check is meaningful across all ranges.
        let denom = w2 * w2 - w1 * w1;
        if denom.abs() <= f64::EPSILON * w2 * w2.max(w1 * w1) {
            return Err(RayleighError::DuplicateFrequency(f1));
        }

        // From zeta(w) = 0.5*(alpha/w + beta*w) matched at the two modes:
        //   beta  = 2 (z2 w2 - z1 w1) / (w2^2 - w1^2)
        //   alpha = 2 w1 w2 (z1 w2 - z2 w1) / (w2^2 - w1^2)
        let beta = 2.0 * (z2 * w2 - z1 * w1) / denom;
        let alpha = 2.0 * w1 * w2 * (z1 * w2 - z2 * w1) / denom;

        Ok(Self { alpha, beta })
    }

    /// Modal damping ratio at angular frequency `omega` (radians per
    /// second): `zeta = 0.5 * (alpha / omega + beta * omega)`.
    ///
    /// `omega` must be strictly positive (the mass term divides by it);
    /// a non-positive `omega` returns `f64::NAN`.
    pub fn ratio_at_omega(&self, omega: f64) -> f64 {
        if omega <= 0.0 {
            return f64::NAN;
        }
        0.5 * (self.alpha / omega + self.beta * omega)
    }

    /// Modal damping ratio at cyclic frequency `f` (hertz).
    ///
    /// Converts to angular frequency `omega = 2*pi*f` and defers to
    /// [`RayleighDamping::ratio_at_omega`]. A non-positive `f` returns
    /// `f64::NAN`.
    pub fn ratio_at_hz(&self, f: f64) -> f64 {
        self.ratio_at_omega(2.0 * PI * f)
    }

    /// The angular frequency `omega_min = sqrt(alpha / beta)` at which the
    /// Rayleigh damping curve reaches its minimum, when one exists.
    ///
    /// The two-term curve `0.5*(alpha/omega + beta*omega)` has an interior
    /// minimum only when `alpha` and `beta` are both strictly positive;
    /// otherwise (a mass-only, stiffness-only, or sign-mixed fit) there is
    /// no interior minimum and this returns `None`.
    pub fn min_ratio_omega(&self) -> Option<f64> {
        if self.alpha > 0.0 && self.beta > 0.0 {
            Some((self.alpha / self.beta).sqrt())
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The fit reproduces both target damping ratios at the two target
    /// frequencies to within ~1e-9 — the defining closed-form property.
    #[test]
    fn fit_recovers_both_targets() {
        let f1 = 2.0;
        let z1 = 0.02;
        let f2 = 50.0;
        let z2 = 0.05;
        let r = RayleighDamping::from_two_modes(f1, z1, f2, z2).unwrap();

        assert!(
            (r.ratio_at_hz(f1) - z1).abs() < 1e-9,
            "zeta(f1) = {} expected {z1}",
            r.ratio_at_hz(f1)
        );
        assert!(
            (r.ratio_at_hz(f2) - z2).abs() < 1e-9,
            "zeta(f2) = {} expected {z2}",
            r.ratio_at_hz(f2)
        );
    }

    /// Equal target ratios at two frequencies give the same recovered
    /// ratio at both — a sanity check on the symmetric case.
    #[test]
    fn equal_targets_recover_equally() {
        let z = 0.03;
        let r = RayleighDamping::from_two_modes(5.0, z, 80.0, z).unwrap();
        assert!((r.ratio_at_hz(5.0) - z).abs() < 1e-9);
        assert!((r.ratio_at_hz(80.0) - z).abs() < 1e-9);
        // Both coefficients are positive when both targets are positive,
        // so the curve has an interior minimum strictly between the
        // fit frequencies (where the ratio dips below z).
        let w_min = r.min_ratio_omega().expect("interior minimum exists");
        let w1 = 2.0 * PI * 5.0;
        let w2 = 2.0 * PI * 80.0;
        assert!(w_min > w1 && w_min < w2);
        assert!(r.ratio_at_omega(w_min) < z);
    }

    /// Mass-only damping (beta = 0): the damping ratio is exactly
    /// proportional to 1/omega, i.e. `zeta(omega)*omega` is constant.
    #[test]
    fn mass_only_scales_as_inverse_omega() {
        let r = RayleighDamping::new(0.6, 0.0);
        assert_eq!(r.beta, 0.0);
        assert!(r.min_ratio_omega().is_none());

        // zeta(omega) = 0.5*alpha/omega  =>  zeta*omega = 0.5*alpha (const).
        let expect = 0.5 * r.alpha;
        for &w in &[1.0_f64, 7.5, 42.0, 1234.0] {
            let prod = r.ratio_at_omega(w) * w;
            assert!(
                (prod - expect).abs() < 1e-12,
                "zeta*omega = {prod} expected {expect} at omega = {w}"
            );
        }
        // Doubling omega halves the ratio.
        let a = r.ratio_at_omega(10.0);
        let b = r.ratio_at_omega(20.0);
        assert!((a / b - 2.0).abs() < 1e-12);
    }

    /// Stiffness-only damping (alpha = 0): the damping ratio is exactly
    /// proportional to omega, i.e. `zeta(omega)/omega` is constant.
    #[test]
    fn stiffness_only_scales_as_omega() {
        let r = RayleighDamping::new(0.0, 1.0e-4);
        assert_eq!(r.alpha, 0.0);
        assert!(r.min_ratio_omega().is_none());

        // zeta(omega) = 0.5*beta*omega  =>  zeta/omega = 0.5*beta (const).
        let expect = 0.5 * r.beta;
        for &w in &[1.0_f64, 7.5, 42.0, 1234.0] {
            let ratio = r.ratio_at_omega(w) / w;
            assert!(
                (ratio - expect).abs() < 1e-15,
                "zeta/omega = {ratio} expected {expect} at omega = {w}"
            );
        }
        // Doubling omega doubles the ratio.
        let a = r.ratio_at_omega(10.0);
        let b = r.ratio_at_omega(20.0);
        assert!((b / a - 2.0).abs() < 1e-12);
    }

    /// A pure mass-only target pair (z proportional to 1/f) fits back to
    /// beta = 0, confirming the solver recovers the degenerate case.
    #[test]
    fn fit_recovers_mass_only_case() {
        // Choose z so that zeta*omega is constant: z_i = c / omega_i.
        let c = 0.4; // = 0.5 * alpha
        let f1 = 3.0;
        let f2 = 60.0;
        let w1 = 2.0 * PI * f1;
        let w2 = 2.0 * PI * f2;
        let r = RayleighDamping::from_two_modes(f1, c / w1, f2, c / w2).unwrap();
        assert!(r.beta.abs() < 1e-12, "beta = {} expected ~0", r.beta);
        assert!(
            (r.alpha - 2.0 * c).abs() < 1e-9,
            "alpha = {} expected {}",
            r.alpha,
            2.0 * c
        );
    }

    /// Bad inputs are rejected with the right error variants.
    #[test]
    fn rejects_invalid_inputs() {
        assert_eq!(
            RayleighDamping::from_two_modes(0.0, 0.02, 50.0, 0.05),
            Err(RayleighError::BadFrequency(0.0))
        );
        assert_eq!(
            RayleighDamping::from_two_modes(-1.0, 0.02, 50.0, 0.05),
            Err(RayleighError::BadFrequency(-1.0))
        );
        // NaN != NaN, so match the variant rather than compare equality.
        assert!(matches!(
            RayleighDamping::from_two_modes(2.0, f64::NAN, 50.0, 0.05),
            Err(RayleighError::BadRatio(r)) if r.is_nan()
        ));
        assert_eq!(
            RayleighDamping::from_two_modes(10.0, 0.02, 10.0, 0.05),
            Err(RayleighError::DuplicateFrequency(10.0))
        );
    }

    /// Non-positive frequencies return NaN from the evaluators rather
    /// than dividing by zero.
    #[test]
    fn non_positive_frequency_is_nan() {
        let r = RayleighDamping::new(0.5, 1e-4);
        assert!(r.ratio_at_omega(0.0).is_nan());
        assert!(r.ratio_at_omega(-3.0).is_nan());
        assert!(r.ratio_at_hz(0.0).is_nan());
    }
}
