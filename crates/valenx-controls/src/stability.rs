//! Second-order characteristic-equation stability.
//!
//! A continuous linear time-invariant system is **asymptotically stable**
//! iff every root of its characteristic polynomial lies strictly in the
//! open left half of the complex plane (all roots have negative real
//! part). For a quadratic
//!
//! ```text
//! a*s^2 + b*s + c = 0   (a != 0)
//! ```
//!
//! the Routh-Hurwitz criterion collapses to a famously simple statement:
//! **the system is stable iff all three coefficients share the same
//! (non-zero) sign** — equivalently, after normalising so `a > 0`, iff
//! `b > 0` and `c > 0`.
//!
//! ## Model
//!
//! This module provides:
//!
//! - [`QuadraticChar`] — a validated quadratic characteristic polynomial.
//! - [`QuadraticChar::is_stable`] — the Routh-Hurwitz test above.
//! - [`QuadraticChar::roots`] — the two (possibly complex) roots in
//!   closed form, so a caller can inspect the real parts directly and
//!   cross-check the criterion.
//!
//! A [`SecondOrder`](crate::second_order::SecondOrder) system
//! `s^2 + 2*zeta*wn*s + wn^2` has `a = 1`, `b = 2*zeta*wn`, `c = wn^2`;
//! with `wn > 0` it is stable exactly when `zeta > 0`, which this module
//! and [`SecondOrder`](crate::second_order::SecondOrder) agree on.

use serde::{Deserialize, Serialize};

use crate::error::{ControlsError, Result};

/// A root of the quadratic characteristic polynomial.
///
/// Real roots are represented with `imaginary == 0.0`; a complex-
/// conjugate pair shares the same `real` with `imaginary` of opposite
/// sign.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Root {
    /// Real part of the root.
    pub real: f64,
    /// Imaginary part of the root.
    pub imaginary: f64,
}

/// A validated quadratic characteristic polynomial `a*s^2 + b*s + c`.
///
/// Construct with [`QuadraticChar::new`], which rejects a non-finite
/// coefficient or a zero leading coefficient (`a == 0`, which would make
/// the polynomial first-order, not second).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct QuadraticChar {
    /// Leading coefficient `a` (coefficient of `s^2`), non-zero.
    pub a: f64,
    /// Coefficient `b` (coefficient of `s`).
    pub b: f64,
    /// Constant term `c`.
    pub c: f64,
}

impl QuadraticChar {
    /// Construct a quadratic characteristic polynomial.
    ///
    /// # Errors
    ///
    /// Returns [`ControlsError::InvalidParameter`] if any coefficient is
    /// non-finite, or if the leading coefficient `a` is zero (the
    /// polynomial would not be second-order).
    pub fn new(a: f64, b: f64, c: f64) -> Result<Self> {
        for (name, v) in [("a", a), ("b", b), ("c", c)] {
            if !v.is_finite() {
                return Err(ControlsError::invalid(name, "coefficient must be finite"));
            }
        }
        if a == 0.0 {
            return Err(ControlsError::invalid(
                "a",
                "leading coefficient must be non-zero for a second-order polynomial",
            ));
        }
        Ok(Self { a, b, c })
    }

    /// Build the characteristic polynomial of a standard second-order
    /// system `s^2 + 2*zeta*wn*s + wn^2`.
    ///
    /// # Errors
    ///
    /// Returns [`ControlsError::InvalidParameter`] if `wn` is non-finite
    /// or `zeta` is non-finite. (Unlike
    /// [`SecondOrder`](crate::second_order::SecondOrder), this does not
    /// require `wn > 0` or `zeta >= 0`, so callers can probe the
    /// stability boundary from either side.)
    pub fn from_wn_zeta(wn: f64, zeta: f64) -> Result<Self> {
        if !wn.is_finite() {
            return Err(ControlsError::invalid("wn", "must be finite"));
        }
        if !zeta.is_finite() {
            return Err(ControlsError::invalid("zeta", "must be finite"));
        }
        Self::new(1.0, 2.0 * zeta * wn, wn * wn)
    }

    /// Routh-Hurwitz stability test for the quadratic.
    ///
    /// Returns `true` iff both roots have strictly negative real part,
    /// which for a quadratic is exactly: all three coefficients share the
    /// same non-zero sign. Implemented by normalising so the leading
    /// coefficient is positive and checking `b > 0 && c > 0`.
    pub fn is_stable(&self) -> bool {
        // Normalise so the s^2 coefficient is positive; multiplying the
        // whole polynomial by -1 does not move its roots.
        let sign = if self.a > 0.0 { 1.0 } else { -1.0 };
        let b = sign * self.b;
        let c = sign * self.c;
        b > 0.0 && c > 0.0
    }

    /// Discriminant `b^2 - 4*a*c` of the quadratic.
    ///
    /// Positive -> two distinct real roots; zero -> a repeated real root;
    /// negative -> a complex-conjugate pair.
    pub fn discriminant(&self) -> f64 {
        self.b * self.b - 4.0 * self.a * self.c
    }

    /// The two roots of the quadratic in closed form.
    ///
    /// For a non-negative discriminant both roots are real
    /// (`imaginary == 0.0`); for a negative discriminant they form a
    /// complex-conjugate pair `real ± i*imaginary`.
    pub fn roots(&self) -> [Root; 2] {
        let disc = self.discriminant();
        let two_a = 2.0 * self.a;
        if disc >= 0.0 {
            let sqrt_disc = disc.sqrt();
            [
                Root {
                    real: (-self.b + sqrt_disc) / two_a,
                    imaginary: 0.0,
                },
                Root {
                    real: (-self.b - sqrt_disc) / two_a,
                    imaginary: 0.0,
                },
            ]
        } else {
            let real = -self.b / two_a;
            let imag = (-disc).sqrt() / two_a;
            [
                Root {
                    real,
                    imaginary: imag,
                },
                Root {
                    real,
                    imaginary: -imag,
                },
            ]
        }
    }

    /// The maximum real part across both roots.
    ///
    /// Strictly negative iff the system is stable; this is an alternative
    /// to [`is_stable`](Self::is_stable) that returns the stability margin
    /// rather than a boolean.
    pub fn max_real_part(&self) -> f64 {
        let [r0, r1] = self.roots();
        r0.real.max(r1.real)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn rejects_bad_coefficients() {
        assert!(QuadraticChar::new(0.0, 1.0, 1.0).is_err()); // a = 0
        assert!(QuadraticChar::new(f64::NAN, 1.0, 1.0).is_err());
        assert!(QuadraticChar::new(1.0, f64::INFINITY, 1.0).is_err());
        assert!(QuadraticChar::new(1.0, 2.0, 1.0).is_ok());
    }

    #[test]
    fn stable_iff_all_coefficients_positive() {
        // VALIDATE: char eqn stable iff coeffs share a sign.
        assert!(QuadraticChar::new(1.0, 3.0, 2.0).unwrap().is_stable());
        // A zero or negative middle coefficient -> unstable / marginal.
        assert!(!QuadraticChar::new(1.0, -3.0, 2.0).unwrap().is_stable());
        assert!(!QuadraticChar::new(1.0, 0.0, 2.0).unwrap().is_stable());
        // Negative constant term -> a positive real root -> unstable.
        assert!(!QuadraticChar::new(1.0, 3.0, -2.0).unwrap().is_stable());
        // All-negative coefficients are the same polynomial * -1: stable.
        assert!(QuadraticChar::new(-1.0, -3.0, -2.0).unwrap().is_stable());
    }

    #[test]
    fn is_stable_agrees_with_root_real_parts() {
        // Stable example: roots of s^2 + 3s + 2 are -1 and -2.
        let stable = QuadraticChar::new(1.0, 3.0, 2.0).unwrap();
        assert!(stable.is_stable());
        let roots = stable.roots();
        assert!(roots[0].real < 0.0 && roots[1].real < 0.0);
        assert!(stable.max_real_part() < 0.0);
        // Roots are exactly -1 and -2 (order: +sqrt first).
        assert!((roots[0].real - (-1.0)).abs() < EPS, "{:?}", roots[0]);
        assert!((roots[1].real - (-2.0)).abs() < EPS, "{:?}", roots[1]);

        // Unstable example: s^2 - 3s + 2 has roots +1 and +2.
        let unstable = QuadraticChar::new(1.0, -3.0, 2.0).unwrap();
        assert!(!unstable.is_stable());
        assert!(unstable.max_real_part() > 0.0);
    }

    #[test]
    fn complex_roots_are_conjugate_pair() {
        // s^2 + 2s + 5 -> roots -1 ± 2i (discriminant 4 - 20 = -16 < 0).
        let q = QuadraticChar::new(1.0, 2.0, 5.0).unwrap();
        assert!(q.discriminant() < 0.0);
        let [r0, r1] = q.roots();
        assert!((r0.real - (-1.0)).abs() < EPS, "{r0:?}");
        assert!((r1.real - (-1.0)).abs() < EPS, "{r1:?}");
        assert!((r0.imaginary - 2.0).abs() < EPS, "{r0:?}");
        assert!((r1.imaginary - (-2.0)).abs() < EPS, "{r1:?}");
        // Negative real part -> stable, and the criterion agrees.
        assert!(q.is_stable());
    }

    #[test]
    fn repeated_real_root_at_critical_damping() {
        // s^2 + 2s + 1 = (s+1)^2 -> double root at -1, discriminant 0.
        let q = QuadraticChar::new(1.0, 2.0, 1.0).unwrap();
        assert!(q.discriminant().abs() < EPS);
        let [r0, r1] = q.roots();
        assert!((r0.real - (-1.0)).abs() < EPS);
        assert!((r1.real - (-1.0)).abs() < EPS);
        assert!(r0.imaginary.abs() < EPS && r1.imaginary.abs() < EPS);
        assert!(q.is_stable());
    }

    #[test]
    fn second_order_helper_matches_wn_zeta() {
        // s^2 + 2*zeta*wn*s + wn^2 with wn = 4, zeta = 0.5:
        // -> a=1, b=4, c=16.
        let q = QuadraticChar::from_wn_zeta(4.0, 0.5).unwrap();
        assert!((q.a - 1.0).abs() < EPS);
        assert!((q.b - 4.0).abs() < EPS);
        assert!((q.c - 16.0).abs() < EPS);
        // Positive damping -> stable; zero damping -> marginal (b = 0).
        assert!(q.is_stable());
        assert!(!QuadraticChar::from_wn_zeta(4.0, 0.0).unwrap().is_stable());
        // Negative damping -> unstable (right-half-plane poles).
        assert!(!QuadraticChar::from_wn_zeta(4.0, -0.3).unwrap().is_stable());
    }

    #[test]
    fn roots_round_trip_through_json() {
        let q = QuadraticChar::new(1.0, 2.0, 5.0).unwrap();
        let json = serde_json::to_string(&q).unwrap();
        let back: QuadraticChar = serde_json::from_str(&json).unwrap();
        assert_eq!(q, back);
    }
}
