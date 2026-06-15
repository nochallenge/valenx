//! Steinhart-Hart three-coefficient thermistor model.
//!
//! ## Model
//!
//! The Steinhart-Hart equation models the resistance-temperature curve
//! of a thermistor far more accurately than the single-`beta`
//! approximation, using three coefficients `A`, `B`, `C`:
//!
//! ```text
//! 1/T = A + B*ln(R) + C*ln(R)^3
//! ```
//!
//! with `T` in kelvin and `R` in ohms. Temperature-from-resistance is
//! the direct evaluation above. Resistance-from-temperature has a
//! closed form via Cardano's solution of the depressed cubic in
//! `x = ln(R)`:
//!
//! ```text
//! x^3 + (B/C) x + (A - 1/T)/C = 0
//! ```
//!
//! Writing `y = (A - 1/T)/C` and `p = B/C`, the single real root is
//!
//! ```text
//! x = cbrt(-y/2 + sqrt(y^2/4 + p^3/27))
//!   + cbrt(-y/2 - sqrt(y^2/4 + p^3/27))
//! ```
//!
//! and `R = exp(x)`. For physical NTC coefficients (`C > 0`,
//! `p = B/C > 0`) the discriminant is positive and exactly one real
//! root exists, so the inverse is unambiguous.
//!
//! ## Honest scope
//!
//! Textbook closed-form model. The three coefficients are usually fit
//! from three calibration points (see
//! [`SteinhartHart::fit_three_point`]); accuracy degrades outside the
//! span those points bracket. Kelvin and ohms throughout. No
//! self-heating, lead-resistance, or tolerance modelling, and no
//! datasheet lookup.

use crate::error::{check_resistance, check_temperature, ThermistorError};
use serde::{Deserialize, Serialize};

/// Steinhart-Hart coefficients `(A, B, C)` for one thermistor.
///
/// Build directly with [`SteinhartHart::new`] when the coefficients are
/// known, or fit them from three measured points with
/// [`SteinhartHart::fit_three_point`]. Convert in either direction with
/// [`SteinhartHart::temperature_at`] and
/// [`SteinhartHart::resistance_at`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SteinhartHart {
    /// The `A` coefficient (units: 1/K).
    a: f64,
    /// The `B` coefficient (units: 1/K).
    b: f64,
    /// The `C` coefficient (units: 1/K).
    c: f64,
}

impl SteinhartHart {
    /// Build a model from explicit coefficients.
    ///
    /// All three must be finite. `C` must be non-zero, since the
    /// resistance-from-temperature inverse divides by `C` and a zero
    /// `C` collapses the model to the (different) beta-style form.
    ///
    /// # Errors
    ///
    /// Returns [`ThermistorError::BadParameter`] if any coefficient is
    /// non-finite, or if `c` is zero.
    pub fn new(a: f64, b: f64, c: f64) -> Result<Self, ThermistorError> {
        for (name, v) in [("a", a), ("b", b), ("c", c)] {
            if !v.is_finite() {
                return Err(ThermistorError::BadParameter {
                    name,
                    value: v,
                    reason: "coefficient must be finite",
                });
            }
        }
        if c == 0.0 {
            return Err(ThermistorError::BadParameter {
                name: "c",
                value: c,
                reason: "C must be non-zero for the cubic inverse",
            });
        }
        Ok(SteinhartHart { a, b, c })
    }

    /// Fit the three coefficients from three measured resistance/
    /// temperature pairs.
    ///
    /// With `L_i = ln(R_i)` and `y_i = 1/T_i`, the model is linear in
    /// `(A, B, C)`:
    ///
    /// ```text
    /// y_i = A + B*L_i + C*L_i^3   for i = 1, 2, 3
    /// ```
    ///
    /// This 3x3 system is solved by Cramer's rule. The three sample
    /// resistances must be pairwise distinct, otherwise the coefficient
    /// matrix is singular and the fit is undetermined.
    ///
    /// # Errors
    ///
    /// Returns [`ThermistorError::BadParameter`] for any out-of-domain
    /// input, and [`ThermistorError::Degenerate`] if two resistances
    /// coincide (singular system). Returns [`ThermistorError::NonFinite`]
    /// if the solved coefficients are not all finite.
    pub fn fit_three_point(points: [(f64, f64); 3]) -> Result<Self, ThermistorError> {
        // Validate and pre-compute ln(R), 1/T for each point.
        let mut l = [0.0_f64; 3];
        let mut y = [0.0_f64; 3];
        for (i, (r, t)) in points.iter().enumerate() {
            let name = R_NAMES[i];
            let r = check_resistance(name, *r)?;
            let t = check_temperature(T_NAMES[i], *t)?;
            l[i] = r.ln();
            y[i] = 1.0 / t;
        }

        // The design matrix [[1, L, L^3]; ...] is a Vandermonde-style
        // matrix in `L = ln(R)`; it is singular exactly when two of the
        // `L` values coincide. Test that directly (scale-relative) so a
        // duplicate resistance is rejected cleanly, rather than slipping
        // through as a tiny-but-nonzero floating-point determinant that
        // then yields garbage coefficients.
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            let scale = l[i].abs().max(l[j].abs()).max(1.0);
            if (l[i] - l[j]).abs() <= scale * LN_R_DISTINCT_EPS {
                return Err(ThermistorError::Degenerate(
                    "the three calibration resistances must be pairwise distinct",
                ));
            }
        }

        // Rows of the design matrix M = [[1, L, L^3]; ...].
        let rows = [
            [1.0, l[0], l[0].powi(3)],
            [1.0, l[1], l[1].powi(3)],
            [1.0, l[2], l[2].powi(3)],
        ];
        let det = det3(&rows);
        if det == 0.0 {
            return Err(ThermistorError::Degenerate(
                "the three calibration resistances must be pairwise distinct",
            ));
        }

        // Cramer's rule: replace column j with the RHS vector y.
        let a = det3(&replace_col(&rows, 0, &y)) / det;
        let b = det3(&replace_col(&rows, 1, &y)) / det;
        let c = det3(&replace_col(&rows, 2, &y)) / det;
        if !(a.is_finite() && b.is_finite() && c.is_finite()) {
            return Err(ThermistorError::NonFinite(
                "fitted Steinhart-Hart coefficients are not all finite",
            ));
        }
        SteinhartHart::new(a, b, c)
    }

    /// The `A` coefficient (units: 1/K).
    pub fn a(&self) -> f64 {
        self.a
    }

    /// The `B` coefficient (units: 1/K).
    pub fn b(&self) -> f64 {
        self.b
    }

    /// The `C` coefficient (units: 1/K).
    pub fn c(&self) -> f64 {
        self.c
    }

    /// Absolute temperature at resistance `r_ohms`, in kelvin.
    ///
    /// Direct evaluation of `1/T = A + B*ln(R) + C*ln(R)^3` followed by
    /// reciprocation.
    ///
    /// # Errors
    ///
    /// Returns [`ThermistorError::BadParameter`] if `r_ohms` is out of
    /// domain, or [`ThermistorError::NonFinite`] if the polynomial
    /// yields a non-positive or non-finite `1/T` (a non-physical
    /// temperature).
    pub fn temperature_at(&self, r_ohms: f64) -> Result<f64, ThermistorError> {
        let r = check_resistance("r_ohms", r_ohms)?;
        let ln_r = r.ln();
        let inv_t = self.a + self.b * ln_r + self.c * ln_r.powi(3);
        if !inv_t.is_finite() || inv_t <= 0.0 {
            return Err(ThermistorError::NonFinite(
                "Steinhart-Hart polynomial implies a non-physical temperature",
            ));
        }
        Ok(1.0 / inv_t)
    }

    /// Resistance at absolute temperature `t_kelvin`, in ohms.
    ///
    /// Inverts the equation through Cardano's formula for the depressed
    /// cubic in `x = ln(R)` described in the module docs, then returns
    /// `exp(x)`. This is the analytic inverse of
    /// [`temperature_at`](SteinhartHart::temperature_at).
    ///
    /// # Errors
    ///
    /// Returns [`ThermistorError::BadParameter`] if `t_kelvin` is out of
    /// domain, or [`ThermistorError::NonFinite`] if the cubic solution
    /// or the final exponential leaves the representable range.
    pub fn resistance_at(&self, t_kelvin: f64) -> Result<f64, ThermistorError> {
        let t = check_temperature("t_kelvin", t_kelvin)?;
        let target = 1.0 / t; // the value of A + B*x + C*x^3 we must hit

        // Near-linear limit. When the cubic coefficient `C` is
        // negligible relative to the linear coefficient `B` over the
        // working range of `x = ln(R)`, the equation degenerates to the
        // beta-style line `1/T = A + B*x`, whose inverse is exact and
        // well conditioned. This is precisely the case where a
        // Steinhart-Hart model has been fit to a pure beta curve (then
        // `C` rounds to ~0), where the cubic formula below would divide
        // by an almost-zero `C` and overflow. We detect it by checking
        // the cubic term's contribution at the linear-estimate root.
        if self.b != 0.0 {
            let x_lin = (target - self.a) / self.b;
            let cubic_term = (self.c * x_lin.powi(3)).abs();
            let linear_term = (self.b * x_lin).abs();
            if cubic_term <= linear_term * CUBIC_NEGLIGIBLE {
                return finite_exp(x_lin);
            }
        }

        // General depressed cubic x^3 + p x + q = 0 with
        // p = B/C, q = (A - 1/T)/C.
        let p = self.b / self.c;
        let q = (self.a - target) / self.c;
        // Cardano: discriminant term under the square root.
        let disc = q * q / 4.0 + p * p * p / 27.0;
        let x = if disc >= 0.0 {
            // One real root.
            let s = disc.sqrt();
            cbrt(-q / 2.0 + s) + cbrt(-q / 2.0 - s)
        } else {
            // Three real roots (trigonometric form), which requires
            // p < 0. Enumerate all three and pick the one that best
            // satisfies the original equation, so the inverse stays a
            // true inverse regardless of which branch is physical.
            let m = 2.0 * (-p / 3.0).sqrt();
            let phi = (3.0 * q / (p * m)).clamp(-1.0, 1.0).acos();
            let candidates = [
                m * (phi / 3.0).cos(),
                m * ((phi + 2.0 * std::f64::consts::PI) / 3.0).cos(),
                m * ((phi + 4.0 * std::f64::consts::PI) / 3.0).cos(),
            ];
            let residual = |x: f64| (x * x * x + p * x + q).abs();
            candidates
                .into_iter()
                .min_by(|a, b| residual(*a).total_cmp(&residual(*b)))
                .expect("three candidate roots")
        };
        if !x.is_finite() {
            return Err(ThermistorError::NonFinite(
                "Steinhart-Hart cubic root is non-finite",
            ));
        }
        finite_exp(x)
    }
}

/// Relative threshold below which the Steinhart-Hart cubic term is
/// treated as numerically absent and the model is inverted as the
/// linear beta-style law instead. Chosen well above `f64` epsilon so a
/// genuinely cubic curve is never mistaken for a linear one, yet small
/// enough that a `C` rounded to ~1e-20 (a beta curve fit to
/// Steinhart-Hart) is caught.
const CUBIC_NEGLIGIBLE: f64 = 1e-9;

/// `exp(x)` guarded against overflow to a non-finite resistance.
fn finite_exp(x: f64) -> Result<f64, ThermistorError> {
    let r = x.exp();
    if !r.is_finite() {
        return Err(ThermistorError::NonFinite(
            "Steinhart-Hart resistance overflowed to a non-finite value",
        ));
    }
    Ok(r)
}

/// Relative tolerance below which two `ln(R)` values are treated as the
/// same point, making the three-point fit singular. Far above `f64`
/// epsilon so distinct-but-close calibration resistances still fit,
/// while exact (or floating-point-identical) duplicates are rejected.
const LN_R_DISTINCT_EPS: f64 = 1e-12;

/// Parameter names used when validating the three fit resistances, so
/// error messages point at the offending point.
const R_NAMES: [&str; 3] = ["r1_ohms", "r2_ohms", "r3_ohms"];
/// Parameter names used when validating the three fit temperatures.
const T_NAMES: [&str; 3] = ["t1_kelvin", "t2_kelvin", "t3_kelvin"];

/// Real cube root, defined for negative arguments (unlike `powf(1/3)`).
fn cbrt(x: f64) -> f64 {
    x.cbrt()
}

/// Determinant of a 3x3 matrix given as three row arrays.
fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Return a copy of `m` with column `col` replaced by the vector `v`.
fn replace_col(m: &[[f64; 3]; 3], col: usize, v: &[f64; 3]) -> [[f64; 3]; 3] {
    let mut out = *m;
    for (row, vi) in out.iter_mut().zip(v.iter()) {
        row[col] = *vi;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Published Steinhart-Hart coefficients for a Vishay-style 10 kohm
    /// NTC (widely tabulated reference set). Used as a fixed model whose
    /// curve we sample and re-fit.
    fn reference() -> SteinhartHart {
        SteinhartHart::new(1.009_249_522e-3, 2.378_405_444e-4, 2.019_202_697e-7)
            .expect("valid coefficients")
    }

    #[test]
    fn temperature_inverts_resistance_round_trip() {
        let m = reference();
        for r in [100.0_f64, 1_000.0, 10_000.0, 32_650.0, 100_000.0] {
            let t = m.temperature_at(r).unwrap();
            let back = m.resistance_at(t).unwrap();
            // The cubic inverse is analytic; tolerate only float noise.
            let rel = (back - r).abs() / r;
            assert!(rel < 1e-9, "round trip failed at R={r}: got {back}");
        }
    }

    #[test]
    fn resistance_inverts_temperature_round_trip() {
        let m = reference();
        for t in [273.15_f64, 298.15, 310.15, 333.15, 373.15] {
            let r = m.resistance_at(t).unwrap();
            let back = m.temperature_at(r).unwrap();
            assert!(
                (back - t).abs() < 1e-7,
                "round trip failed at T={t}: got {back}"
            );
        }
    }

    #[test]
    fn ntc_resistance_falls_as_temperature_rises() {
        let m = reference();
        let r0 = m.resistance_at(273.15).unwrap();
        let r25 = m.resistance_at(298.15).unwrap();
        let r50 = m.resistance_at(323.15).unwrap();
        assert!(
            r0 > r25 && r25 > r50,
            "expected monotone NTC: {r0} {r25} {r50}"
        );
    }

    #[test]
    fn known_room_temperature_resistance_is_about_10k() {
        // This published coefficient set is a ~10 kohm part: near room
        // temperature its resistance is in the low-10k band. The exact
        // 10 kohm crossing for these particular coefficients is at
        // ~24.68 C (297.83 K), not a round 25 C — assert both facts so
        // the test pins the curve, not a wished-for value.
        let m = reference();
        let r25 = m.resistance_at(298.15).unwrap(); // 25 C
        assert!(
            (9_000.0..=11_000.0).contains(&r25),
            "expected ~10k near room temp, got {r25}"
        );
        // The temperature at exactly 10 kohm, recovered by the model.
        let t_at_10k = m.temperature_at(10_000.0).unwrap();
        assert!(
            (t_at_10k - 297.8313).abs() < 1e-3,
            "10k crossing at {t_at_10k} K"
        );
    }

    #[test]
    fn three_point_fit_recovers_coefficients() {
        // Sample the reference curve at three temperatures, then re-fit.
        let truth = reference();
        let temps = [273.15_f64, 298.15, 323.15];
        let mut pts = [(0.0, 0.0); 3];
        for (i, t) in temps.iter().enumerate() {
            pts[i] = (truth.resistance_at(*t).unwrap(), *t);
        }
        let fitted = SteinhartHart::fit_three_point(pts).unwrap();
        assert!(
            (fitted.a() - truth.a()).abs() < 1e-12,
            "A {} vs {}",
            fitted.a(),
            truth.a()
        );
        assert!(
            (fitted.b() - truth.b()).abs() < 1e-12,
            "B {} vs {}",
            fitted.b(),
            truth.b()
        );
        assert!(
            (fitted.c() - truth.c()).abs() < 1e-15,
            "C {} vs {}",
            fitted.c(),
            truth.c()
        );
    }

    #[test]
    fn near_linear_model_inverts_via_fallback() {
        // A Steinhart-Hart model whose cubic coefficient is essentially
        // zero is just the beta-style line 1/T = A + B*ln(R). The
        // resistance-from-temperature inverse must take the linear
        // fallback and still round-trip exactly, rather than dividing by
        // the near-zero C and overflowing.
        let a = 1.022_284_695e-3_f64;
        let b = 2.531_645_570e-4_f64;
        let c = 1e-20; // non-zero (constructor requires it) but negligible
        let m = SteinhartHart::new(a, b, c).unwrap();
        for t in [273.15_f64, 298.15, 313.15, 333.15] {
            let r = m.resistance_at(t).unwrap();
            assert!(r.is_finite(), "fallback produced non-finite R at {t}");
            let back = m.temperature_at(r).unwrap();
            assert!(
                (back - t).abs() < 1e-6,
                "near-linear round trip failed at {t}: got {back}"
            );
        }
    }

    #[test]
    fn fit_then_predicts_fourth_point() {
        // Fit on three points, then check a held-out fourth temperature
        // reproduces the reference curve.
        let truth = reference();
        let temps = [273.15_f64, 298.15, 323.15];
        let mut pts = [(0.0, 0.0); 3];
        for (i, t) in temps.iter().enumerate() {
            pts[i] = (truth.resistance_at(*t).unwrap(), *t);
        }
        let fitted = SteinhartHart::fit_three_point(pts).unwrap();
        let held_out_t = 310.15;
        let r_truth = truth.resistance_at(held_out_t).unwrap();
        let r_fit = fitted.resistance_at(held_out_t).unwrap();
        assert!((r_truth - r_fit).abs() / r_truth < 1e-9);
    }

    #[test]
    fn fit_rejects_duplicate_resistance() {
        let pts = [(10_000.0, 298.15), (10_000.0, 310.15), (5_000.0, 323.15)];
        let err = SteinhartHart::fit_three_point(pts).unwrap_err();
        assert_eq!(err.code(), "thermistor.degenerate");
    }

    #[test]
    fn constructor_rejects_zero_c_and_non_finite() {
        assert!(SteinhartHart::new(1e-3, 2e-4, 0.0).is_err());
        assert!(SteinhartHart::new(f64::NAN, 2e-4, 2e-7).is_err());
    }

    #[test]
    fn det3_matches_hand_value() {
        // Identity determinant is 1; a simple known matrix checks signs.
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!((det3(&id) - 1.0).abs() < 1e-15);
        let m = [[2.0, 0.0, 1.0], [3.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        // Expand: det = 3 (computed by hand).
        assert!((det3(&m) - 3.0).abs() < 1e-12, "got {}", det3(&m));
    }
}
