//! One-compartment intravenous (IV) bolus pharmacokinetics.
//!
//! The body is idealised as a single well-stirred compartment of apparent
//! volume `V` (volume of distribution) into which the whole `dose` is
//! placed instantaneously at `t = 0`, and from which the drug is cleared
//! by a first-order process of clearance `CL`. The amount `A(t)` then
//! obeys `dA/dt = -k·A` with `k = CL/V`, so the concentration is the
//! decaying exponential
//!
//! `C(t) = (dose / V) · exp(-k·t)`.
//!
//! From that single expression every standard descriptor follows in
//! closed form:
//!
//! - elimination rate constant `k = CL / V`
//! - peak concentration `Cmax = C(0) = dose / V`
//! - terminal half-life `t½ = ln(2) / k`
//! - total exposure `AUC(0→∞) = dose / CL`
//! - the time to fall to a threshold concentration `c`,
//!   `t = (1/k)·ln(Cmax / c)`.
//!
//! All quantities are in consistent units chosen by the caller (e.g.
//! dose mg, volume L, clearance L/h ⇒ concentration mg/L, time h); the
//! crate is unit-agnostic and never converts.

use crate::error::{require_non_negative, require_positive, Result};
use serde::{Deserialize, Serialize};

/// Parameters of a one-compartment IV-bolus model.
///
/// Construct with [`OneCompartmentIv::new`], which validates that the
/// volume and clearance are strictly positive and the dose is
/// non-negative; the fields are then guaranteed sane for every method.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct OneCompartmentIv {
    dose: f64,
    volume: f64,
    clearance: f64,
}

impl OneCompartmentIv {
    /// Build a model from a `dose`, an apparent volume of distribution
    /// `volume` (`V`), and a `clearance` (`CL`).
    ///
    /// # Errors
    ///
    /// Returns [`PkError::NotPositive`](crate::PkError::NotPositive) if
    /// `volume` or `clearance` is not strictly positive (or non-finite),
    /// and [`PkError::Negative`](crate::PkError::Negative) if `dose` is
    /// negative (or non-finite).
    pub fn new(dose: f64, volume: f64, clearance: f64) -> Result<Self> {
        let dose = require_non_negative("dose", dose)?;
        let volume = require_positive("volume", volume)?;
        let clearance = require_positive("clearance", clearance)?;
        Ok(Self {
            dose,
            volume,
            clearance,
        })
    }

    /// The administered dose.
    #[inline]
    pub fn dose(&self) -> f64 {
        self.dose
    }

    /// The apparent volume of distribution `V`.
    #[inline]
    pub fn volume(&self) -> f64 {
        self.volume
    }

    /// The clearance `CL`.
    #[inline]
    pub fn clearance(&self) -> f64 {
        self.clearance
    }

    /// First-order elimination rate constant `k = CL / V`.
    #[inline]
    pub fn elimination_rate(&self) -> f64 {
        self.clearance / self.volume
    }

    /// Peak concentration `Cmax = C(0) = dose / V`, reached immediately at
    /// the moment of the bolus.
    #[inline]
    pub fn cmax(&self) -> f64 {
        self.dose / self.volume
    }

    /// Concentration `C(t) = (dose / V)·exp(-k·t)` at time `t >= 0`.
    ///
    /// A negative or non-finite `t` would extrapolate the model outside
    /// its domain, so it is rejected.
    ///
    /// # Errors
    ///
    /// Returns [`PkError::Negative`](crate::PkError::Negative) if `t` is
    /// negative or non-finite.
    pub fn concentration(&self, t: f64) -> Result<f64> {
        let t = require_non_negative("t", t)?;
        let k = self.elimination_rate();
        Ok(self.cmax() * (-k * t).exp())
    }

    /// Terminal elimination half-life `t½ = ln(2) / k`, the time for the
    /// concentration to halve.
    #[inline]
    pub fn half_life(&self) -> f64 {
        std::f64::consts::LN_2 / self.elimination_rate()
    }

    /// Total area under the concentration-time curve
    /// `AUC(0→∞) = dose / CL`.
    ///
    /// This is the integral of [`concentration`](Self::concentration)
    /// over `[0, ∞)`; for the mono-exponential it equals `Cmax / k`,
    /// which simplifies to `dose / CL`.
    #[inline]
    pub fn auc(&self) -> f64 {
        self.dose / self.clearance
    }

    /// Time at which the falling concentration first reaches the
    /// threshold `c`, i.e. the `t` solving `C(t) = c`:
    /// `t = (1/k)·ln(Cmax / c)`.
    ///
    /// The threshold must be strictly positive (the curve never reaches
    /// zero in finite time) and must not exceed `Cmax` (the concentration
    /// only decreases, so a threshold above the peak is never reached for
    /// `t >= 0`).
    ///
    /// # Errors
    ///
    /// Returns [`PkError::NotPositive`](crate::PkError::NotPositive) if
    /// `c` is not strictly positive, and
    /// [`PkError::NotPositive`](crate::PkError::NotPositive) with name
    /// `"threshold-above-cmax"` if `c > Cmax` (unreachable while
    /// declining).
    pub fn time_to_threshold(&self, c: f64) -> Result<f64> {
        let c = require_positive("threshold", c)?;
        let cmax = self.cmax();
        if c > cmax {
            // A descending curve never climbs back to a level above its
            // starting peak; reuse the "not positive" variant to flag an
            // out-of-domain threshold with a descriptive name.
            return Err(crate::PkError::NotPositive {
                name: "threshold-above-cmax",
                value: c,
            });
        }
        let k = self.elimination_rate();
        Ok((cmax / c).ln() / k)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-12;

    fn model() -> OneCompartmentIv {
        // dose 100 mg, V 10 L, CL 5 L/h → k = 0.5 /h, Cmax 10 mg/L.
        OneCompartmentIv::new(100.0, 10.0, 5.0).unwrap()
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(OneCompartmentIv::new(100.0, 0.0, 5.0).is_err());
        assert!(OneCompartmentIv::new(100.0, -1.0, 5.0).is_err());
        assert!(OneCompartmentIv::new(100.0, 10.0, 0.0).is_err());
        assert!(OneCompartmentIv::new(-1.0, 10.0, 5.0).is_err());
        assert!(OneCompartmentIv::new(f64::NAN, 10.0, 5.0).is_err());
        assert!(OneCompartmentIv::new(100.0, f64::INFINITY, 5.0).is_err());
        // A zero dose is a degenerate but valid model.
        assert!(OneCompartmentIv::new(0.0, 10.0, 5.0).is_ok());
    }

    #[test]
    fn elimination_rate_is_cl_over_v() {
        let k = model().elimination_rate();
        assert!((k - 0.5).abs() < TOL, "k = {k}");
    }

    #[test]
    fn concentration_at_zero_is_dose_over_volume() {
        // ANALYTIC: C(0) = dose / V.
        let m = model();
        let c0 = m.concentration(0.0).unwrap();
        assert!((c0 - 10.0).abs() < TOL, "C(0) = {c0}");
        assert!((c0 - m.dose() / m.volume()).abs() < TOL);
        // And it equals Cmax.
        assert!((c0 - m.cmax()).abs() < TOL);
    }

    #[test]
    fn concentration_at_half_life_is_half_of_c0() {
        // ANALYTIC: C(t½) = C(0) / 2.
        let m = model();
        let c0 = m.concentration(0.0).unwrap();
        let c = m.concentration(m.half_life()).unwrap();
        assert!(
            (c - c0 / 2.0).abs() < TOL,
            "C(t½) = {c}, C(0)/2 = {}",
            c0 / 2.0
        );
    }

    #[test]
    fn half_life_matches_ln2_over_k() {
        let m = model();
        let expected = std::f64::consts::LN_2 / 0.5;
        assert!(
            (m.half_life() - expected).abs() < TOL,
            "t½ = {}",
            m.half_life()
        );
    }

    #[test]
    fn auc_is_dose_over_clearance() {
        // ANALYTIC: AUC = dose / CL = 100 / 5 = 20.
        let m = model();
        assert!((m.auc() - 20.0).abs() < TOL, "AUC = {}", m.auc());
        assert!((m.auc() - m.dose() / m.clearance()).abs() < TOL);
    }

    #[test]
    fn auc_equals_cmax_over_k() {
        // Cross-check the two closed forms agree: dose/CL == Cmax/k.
        let m = model();
        let via_k = m.cmax() / m.elimination_rate();
        assert!((m.auc() - via_k).abs() < TOL, "{} vs {via_k}", m.auc());
    }

    #[test]
    fn concentration_decays_monotonically() {
        let m = model();
        let mut prev = f64::INFINITY;
        for i in 0..20 {
            let t = i as f64 * 0.5;
            let c = m.concentration(t).unwrap();
            assert!(c < prev, "C should decrease: C({t}) = {c}, prev = {prev}");
            prev = c;
        }
    }

    #[test]
    fn time_to_threshold_round_trips() {
        // ANALYTIC: solving C(t) = c and then evaluating C(t) recovers c.
        let m = model();
        let target = 2.5; // mg/L (= Cmax / 4 → t = 2·t½).
        let t = m.time_to_threshold(target).unwrap();
        // Cmax/4 → ln(4)/k = 2·ln2/0.5 = 4·ln2 ≈ 2.7726.
        let expected_t = (10.0_f64 / 2.5).ln() / 0.5;
        assert!((t - expected_t).abs() < TOL, "t = {t}");
        let c = m.concentration(t).unwrap();
        assert!((c - target).abs() < 1e-10, "round-trip C = {c}");
    }

    #[test]
    fn time_to_threshold_rejects_unreachable_levels() {
        let m = model();
        // Above the peak: never reached while declining.
        assert!(m.time_to_threshold(11.0).is_err());
        // Non-positive: the curve never reaches zero in finite time.
        assert!(m.time_to_threshold(0.0).is_err());
        assert!(m.time_to_threshold(-1.0).is_err());
        // Exactly Cmax → t = 0.
        let t = m.time_to_threshold(m.cmax()).unwrap();
        assert!(t.abs() < TOL, "t at Cmax = {t}");
    }

    #[test]
    fn concentration_rejects_negative_time() {
        assert!(model().concentration(-1.0).is_err());
        assert!(model().concentration(f64::NAN).is_err());
    }
}
