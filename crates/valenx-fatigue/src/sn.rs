//! Basquin stress-life (S-N) curve.
//!
//! ## Model
//!
//! The high-cycle fatigue life of a material under fully-reversed
//! loading is described by **Basquin's power law**, a straight line on a
//! log-log plot of stress amplitude `S` against cycles-to-failure `N`:
//!
//! ```text
//! S = a * N^b
//! ```
//!
//! - `a` is the **fatigue strength coefficient** (units of stress); it
//!   is the value of `S` extrapolated to `N = 1` cycle, so `a > 0`.
//! - `b` is the **fatigue strength exponent** (dimensionless, Basquin's
//!   exponent); because life *grows* as stress *falls*, `b < 0`.
//!   Typical metals fall in `-0.05 .. -0.12`.
//!
//! Inverting for the life at a given stress amplitude:
//!
//! ```text
//! N = (S / a)^(1/b)
//! ```
//!
//! Many materials (notably ferrous alloys and titanium) show a fatigue
//! or **endurance limit** `Se`: a stress below which the part survives
//! indefinitely. When an endurance limit is configured, the curve is
//! capped horizontally — stress amplitudes at or below `Se` return an
//! "infinite life" verdict rather than a finite cycle count.
//!
//! ## Honest scope
//!
//! This is the textbook closed-form Basquin relation with an optional
//! flat endurance cutoff. It is **not** a full strain-life
//! (Coffin-Manson) low-cycle model, carries no scatter band / reliability
//! factor, and applies no Marin surface / size / loading / temperature
//! knock-down factors — the caller supplies an already-corrected curve.
//! Research/educational grade, not a production design tool.

use crate::error::{FatigueError, Result};
use serde::{Deserialize, Serialize};

/// A Basquin stress-life curve `S = a * N^b`, optionally capped by a
/// horizontal endurance limit.
///
/// Build one with [`SnCurve::new`] (which validates the coefficients) or
/// with [`SnCurve::from_two_points`] (which fits `a` and `b` to two
/// measured `(N, S)` points). Use [`SnCurve::with_endurance_limit`] to
/// attach a fatigue limit afterwards.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct SnCurve {
    /// Fatigue strength coefficient `a` (stress at `N = 1`), `a > 0`.
    pub a: f64,
    /// Fatigue strength exponent `b` (Basquin's exponent), `b < 0`.
    pub b: f64,
    /// Optional endurance limit `Se`: stress amplitudes at or below this
    /// value are treated as infinite-life. `None` means the power law
    /// extrapolates forever (no fatigue limit, e.g. many aluminium
    /// alloys).
    pub endurance_limit: Option<f64>,
}

/// The outcome of evaluating [`SnCurve::cycles_to_failure`].
///
/// A stress amplitude at or below the configured endurance limit yields
/// [`Life::Infinite`]; anything above yields a finite [`Life::Finite`]
/// cycle count.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Life {
    /// A finite predicted number of cycles to failure.
    Finite(f64),
    /// The stress amplitude is at or below the endurance limit — the
    /// part is predicted to survive indefinitely.
    Infinite,
}

impl Life {
    /// The finite cycle count, or `None` for [`Life::Infinite`].
    pub fn finite(self) -> Option<f64> {
        match self {
            Life::Finite(n) => Some(n),
            Life::Infinite => None,
        }
    }

    /// `true` for [`Life::Infinite`].
    pub fn is_infinite(self) -> bool {
        matches!(self, Life::Infinite)
    }
}

impl SnCurve {
    /// Build a curve from raw Basquin coefficients.
    ///
    /// # Errors
    ///
    /// Returns [`FatigueError::Invalid`] if `a` is not strictly positive
    /// or not finite, or if `b` is not strictly negative or not finite.
    pub fn new(a: f64, b: f64) -> Result<Self> {
        if !a.is_finite() || a <= 0.0 {
            return Err(FatigueError::invalid(
                "a",
                format!("strength coefficient must be finite and > 0, got {a}"),
            ));
        }
        if !b.is_finite() || b >= 0.0 {
            return Err(FatigueError::invalid(
                "b",
                format!("Basquin exponent must be finite and < 0, got {b}"),
            ));
        }
        Ok(SnCurve {
            a,
            b,
            endurance_limit: None,
        })
    }

    /// Attach an endurance limit `se` to the curve (consuming `self`).
    ///
    /// # Errors
    ///
    /// Returns [`FatigueError::Invalid`] if `se` is not strictly positive
    /// or not finite.
    pub fn with_endurance_limit(mut self, se: f64) -> Result<Self> {
        if !se.is_finite() || se <= 0.0 {
            return Err(FatigueError::invalid(
                "endurance_limit",
                format!("endurance limit must be finite and > 0, got {se}"),
            ));
        }
        self.endurance_limit = Some(se);
        Ok(self)
    }

    /// Fit `a` and `b` to two measured `(N, S)` points on the curve.
    ///
    /// Solves the log-log straight line through `(n1, s1)` and
    /// `(n2, s2)`: `b = ln(s2/s1) / ln(n2/n1)` and `a = s1 / n1^b`.
    ///
    /// A common construction takes the point `(1e3, 0.9 Su)` and the
    /// fatigue-limit point `(1e6, Se)` for a steel.
    ///
    /// # Errors
    ///
    /// Returns [`FatigueError::Invalid`] if any stress or cycle value is
    /// not strictly positive / finite, or if the two cycle counts are
    /// equal (a vertical line has no slope). Returns
    /// [`FatigueError::Domain`] if the fitted exponent is not negative
    /// (the second point does not have a lower stress at a higher life).
    pub fn from_two_points(n1: f64, s1: f64, n2: f64, s2: f64) -> Result<Self> {
        for (name, v) in [("n1", n1), ("s1", s1), ("n2", n2), ("s2", s2)] {
            if !v.is_finite() || v <= 0.0 {
                return Err(FatigueError::invalid(
                    name,
                    format!("must be finite and > 0, got {v}"),
                ));
            }
        }
        if (n1 - n2).abs() <= f64::EPSILON * n1.max(n2) {
            return Err(FatigueError::invalid(
                "n2",
                "the two cycle counts must differ to define a slope".to_string(),
            ));
        }
        let b = (s2 / s1).ln() / (n2 / n1).ln();
        if !b.is_finite() || b >= 0.0 {
            return Err(FatigueError::domain(format!(
                "fitted Basquin exponent is not negative (b = {b}); the \
                 higher-life point must have the lower stress"
            )));
        }
        let a = s1 / n1.powf(b);
        SnCurve::new(a, b)
    }

    /// Stress amplitude `S = a * N^b` predicted at a given finite life
    /// `n` (cycles).
    ///
    /// # Errors
    ///
    /// Returns [`FatigueError::Invalid`] if `n` is not strictly positive
    /// or not finite.
    pub fn stress_at_cycles(&self, n: f64) -> Result<f64> {
        if !n.is_finite() || n <= 0.0 {
            return Err(FatigueError::invalid(
                "n",
                format!("cycle count must be finite and > 0, got {n}"),
            ));
        }
        Ok(self.a * n.powf(self.b))
    }

    /// Predicted life at a stress amplitude `s`.
    ///
    /// If an endurance limit is configured and `s <= Se`, returns
    /// [`Life::Infinite`]. Otherwise returns the finite Basquin life
    /// `N = (s / a)^(1/b)` wrapped in [`Life::Finite`].
    ///
    /// # Errors
    ///
    /// Returns [`FatigueError::Invalid`] if `s` is not strictly positive
    /// or not finite.
    pub fn cycles_to_failure(&self, s: f64) -> Result<Life> {
        if !s.is_finite() || s <= 0.0 {
            return Err(FatigueError::invalid(
                "s",
                format!("stress amplitude must be finite and > 0, got {s}"),
            ));
        }
        if let Some(se) = self.endurance_limit {
            if s <= se {
                return Ok(Life::Infinite);
            }
        }
        // N = (s / a)^(1/b). With a > 0, s > 0 and b < 0 this is always a
        // finite positive number.
        let n = (s / self.a).powf(1.0 / self.b);
        Ok(Life::Finite(n))
    }

    /// The finite life at a stress amplitude `s`, ignoring any endurance
    /// limit (the bare Basquin inverse).
    ///
    /// Useful when the caller wants the extrapolated power-law life even
    /// below the fatigue limit.
    ///
    /// # Errors
    ///
    /// Returns [`FatigueError::Invalid`] if `s` is not strictly positive
    /// or not finite.
    pub fn cycles_to_failure_unbounded(&self, s: f64) -> Result<f64> {
        if !s.is_finite() || s <= 0.0 {
            return Err(FatigueError::invalid(
                "s",
                format!("stress amplitude must be finite and > 0, got {s}"),
            ));
        }
        Ok((s / self.a).powf(1.0 / self.b))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tight epsilon for analytic round-trips, scaled relative.
    fn close(x: f64, y: f64) {
        let tol = 1e-9 * x.abs().max(y.abs()).max(1.0);
        assert!((x - y).abs() < tol, "expected {x} ~= {y}");
    }

    #[test]
    fn constructor_rejects_bad_coefficients() {
        assert!(SnCurve::new(0.0, -0.1).is_err());
        assert!(SnCurve::new(-1.0, -0.1).is_err());
        assert!(SnCurve::new(f64::NAN, -0.1).is_err());
        assert!(SnCurve::new(1000.0, 0.0).is_err());
        assert!(SnCurve::new(1000.0, 0.1).is_err());
        assert!(SnCurve::new(1000.0, f64::INFINITY).is_err());
        assert!(SnCurve::new(1000.0, -0.1).is_ok());
    }

    #[test]
    fn endurance_limit_rejects_bad_values() {
        let c = SnCurve::new(1000.0, -0.1).unwrap();
        assert!(c.with_endurance_limit(0.0).is_err());
        assert!(c.with_endurance_limit(-5.0).is_err());
        assert!(c.with_endurance_limit(f64::NAN).is_err());
        assert!(c.with_endurance_limit(200.0).is_ok());
    }

    /// S-N round-trip: stress -> cycles -> stress recovers the input.
    #[test]
    fn stress_cycles_round_trip() {
        let c = SnCurve::new(1200.0, -0.085).unwrap();
        for &s in &[1000.0, 600.0, 350.0, 150.0] {
            let n = c.cycles_to_failure_unbounded(s).unwrap();
            let s_back = c.stress_at_cycles(n).unwrap();
            close(s_back, s);
        }
    }

    /// The inverse-direction round-trip: cycles -> stress -> cycles.
    #[test]
    fn cycles_stress_round_trip() {
        let c = SnCurve::new(1200.0, -0.085).unwrap();
        for &n in &[1.0e3, 1.0e4, 1.0e5, 1.0e6, 5.0e6] {
            let s = c.stress_at_cycles(n).unwrap();
            let n_back = c.cycles_to_failure_unbounded(s).unwrap();
            close(n_back, n);
        }
    }

    /// At N = 1 cycle, S = a exactly (the strength coefficient).
    #[test]
    fn stress_at_one_cycle_equals_a() {
        let c = SnCurve::new(1837.0, -0.0977).unwrap();
        close(c.stress_at_cycles(1.0).unwrap(), 1837.0);
    }

    /// Fitting two points and reading them back is exact.
    #[test]
    fn from_two_points_passes_through_both() {
        // Classic steel construction: (1e3, 0.9*Su) and (1e6, Se).
        let su = 1000.0;
        let se = 0.5 * su;
        let s1 = 0.9 * su;
        let c = SnCurve::from_two_points(1.0e3, s1, 1.0e6, se).unwrap();
        close(c.stress_at_cycles(1.0e3).unwrap(), s1);
        close(c.stress_at_cycles(1.0e6).unwrap(), se);
        // And the analytic slope b = log10(s2/s1)/3 over 3 decades.
        let b_expected = (se / s1).ln() / (1.0e6_f64 / 1.0e3).ln();
        close(c.b, b_expected);
    }

    #[test]
    fn from_two_points_rejects_non_decreasing() {
        // Higher life with a HIGHER stress -> positive slope -> rejected.
        assert!(SnCurve::from_two_points(1.0e3, 400.0, 1.0e6, 500.0).is_err());
        // Equal cycle counts -> no slope -> rejected.
        assert!(SnCurve::from_two_points(1.0e3, 500.0, 1.0e3, 400.0).is_err());
        // Non-positive inputs -> rejected.
        assert!(SnCurve::from_two_points(0.0, 500.0, 1.0e6, 250.0).is_err());
    }

    /// Endurance limit caps life: below Se is infinite, above is finite.
    #[test]
    fn endurance_limit_gives_infinite_life() {
        let c = SnCurve::new(1000.0, -0.1)
            .unwrap()
            .with_endurance_limit(200.0)
            .unwrap();
        assert!(c.cycles_to_failure(150.0).unwrap().is_infinite());
        assert!(c.cycles_to_failure(200.0).unwrap().is_infinite()); // at the limit
        let life = c.cycles_to_failure(400.0).unwrap();
        assert!(!life.is_infinite());
        assert!(life.finite().unwrap() > 0.0);
    }

    /// Without an endurance limit the curve extrapolates to a finite life
    /// at any positive stress.
    #[test]
    fn no_endurance_limit_is_always_finite() {
        let c = SnCurve::new(1000.0, -0.1).unwrap();
        assert!(!c.cycles_to_failure(10.0).unwrap().is_infinite());
    }

    /// Lower stress amplitude always predicts a longer life (monotonic).
    #[test]
    fn lower_stress_means_longer_life() {
        let c = SnCurve::new(1200.0, -0.09).unwrap();
        let n_high = c.cycles_to_failure_unbounded(600.0).unwrap();
        let n_low = c.cycles_to_failure_unbounded(300.0).unwrap();
        assert!(
            n_low > n_high,
            "n_low {n_low} should exceed n_high {n_high}"
        );
    }

    #[test]
    fn evaluation_rejects_bad_inputs() {
        let c = SnCurve::new(1000.0, -0.1).unwrap();
        assert!(c.stress_at_cycles(0.0).is_err());
        assert!(c.stress_at_cycles(-5.0).is_err());
        assert!(c.cycles_to_failure(0.0).is_err());
        assert!(c.cycles_to_failure_unbounded(f64::NAN).is_err());
    }
}
