//! Darcy friction factor correlations.
//!
//! The (Darcy-Weisbach) friction factor `f` is the dimensionless
//! coefficient that closes the head-loss equation in
//! [`crate::headloss`]. Two regimes are covered exactly:
//!
//! - **Laminar** (`Re < 2300`): the exact analytic result for fully
//!   developed circular-pipe flow,
//!   ```text
//!     f = 64 / Re
//!   ```
//! - **Turbulent** (`Re >= 4000`): the explicit **Haaland** correlation,
//!   ```text
//!     1 / sqrt(f) = -1.8 * log10[ (eps_rel / 3.7)^1.11 + 6.9 / Re ]
//!   ```
//!   where `eps_rel = epsilon / D` is the relative pipe roughness. The
//!   Haaland equation (S. E. Haaland, 1983) is an explicit algebraic
//!   approximation to the implicit Colebrook-White equation
//!   ```text
//!     1 / sqrt(f) = -2.0 * log10[ eps_rel / 3.7 + 2.51 / (Re * sqrt(f)) ]
//!   ```
//!   and reproduces it to within a couple of percent over the whole
//!   turbulent range, without iteration.
//!
//! In the **transitional** band `2300 <= Re < 4000` no correlation is
//! reliable; the [`friction_factor`] dispatcher returns the Haaland
//! value there but flags the regime so the caller can decide how much to
//! trust it.

use crate::error::{require_in_closed, require_positive, PipeFlowError};
use crate::reynolds::FlowRegime;

/// Darcy friction factor for fully developed **laminar** circular-pipe
/// flow: the exact `f = 64 / Re`.
///
/// `re` must be finite and strictly positive. This is the analytic
/// Hagen-Poiseuille result and carries no roughness dependence — laminar
/// friction is independent of wall roughness.
///
/// # Examples
///
/// ```
/// use valenx_pipeflow::friction::laminar_friction_factor;
///
/// let f = laminar_friction_factor(2000.0).unwrap();
/// assert!((f - 0.032).abs() < 1e-12); // 64/2000 = 0.032
/// ```
pub fn laminar_friction_factor(re: f64) -> Result<f64, PipeFlowError> {
    let re = require_positive("re", re)?;
    Ok(64.0 / re)
}

/// Darcy friction factor for **turbulent** flow by the explicit Haaland
/// correlation.
///
/// `re` must be finite and strictly positive; `relative_roughness`
/// (`epsilon / D`) must be finite and in `[0, 0.05]` (the range over
/// which the Moody chart and these correlations are defined — `0` is a
/// hydraulically smooth pipe).
///
/// Returns the Darcy (not Fanning) friction factor. To recover the
/// Fanning factor divide by four.
///
/// # Examples
///
/// ```
/// use valenx_pipeflow::friction::haaland_friction_factor;
///
/// // Smooth pipe at Re = 1e5.
/// let f = haaland_friction_factor(1.0e5, 0.0).unwrap();
/// assert!((f - 0.0178).abs() < 1e-3);
/// ```
pub fn haaland_friction_factor(re: f64, relative_roughness: f64) -> Result<f64, PipeFlowError> {
    let re = require_positive("re", re)?;
    let eps_rel = require_in_closed(
        "relative_roughness",
        relative_roughness,
        0.0,
        0.05,
        "[0, 0.05]",
    )?;
    let inv_sqrt_f = -1.8 * ((eps_rel / 3.7).powf(1.11) + 6.9 / re).log10();
    Ok(1.0 / (inv_sqrt_f * inv_sqrt_f))
}

/// Darcy friction factor dispatched on the flow regime.
///
/// - `Re < 2300` -> exact laminar `64/Re` (roughness ignored).
/// - `Re >= 2300` -> Haaland correlation with the supplied relative
///   roughness (this includes the transitional band, where the value is
///   the best available smooth estimate but should be treated with
///   caution — see [`FrictionResult::regime`]).
///
/// The returned [`FrictionResult`] bundles the friction factor with the
/// classified [`FlowRegime`] so callers do not silently apply a
/// turbulent correlation to laminar flow (or vice versa).
pub fn friction_factor(re: f64, relative_roughness: f64) -> Result<FrictionResult, PipeFlowError> {
    let re_checked = require_positive("re", re)?;
    let regime = FlowRegime::classify(re_checked);
    let f = if regime.is_laminar() {
        laminar_friction_factor(re_checked)?
    } else {
        haaland_friction_factor(re_checked, relative_roughness)?
    };
    Ok(FrictionResult {
        reynolds: re_checked,
        friction_factor: f,
        regime,
    })
}

/// Result of a [`friction_factor`] evaluation: the friction factor
/// alongside the Reynolds number and the regime it was computed for.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct FrictionResult {
    /// The Reynolds number the factor was evaluated at.
    pub reynolds: f64,
    /// The Darcy friction factor (dimensionless).
    pub friction_factor: f64,
    /// The flow regime the factor was computed under.
    pub regime: FlowRegime,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference implementation: solve the implicit Colebrook-White
    /// equation by fixed-point iteration. Used purely to validate that
    /// the explicit Haaland correlation tracks it to within a few
    /// percent.
    fn colebrook_white(re: f64, eps_rel: f64) -> f64 {
        // 1/sqrt(f) = -2 log10( eps/3.7 + 2.51 / (Re sqrt(f)) ).
        let mut f: f64 = 0.02; // Reasonable seed.
        for _ in 0..200 {
            let rhs = -2.0 * (eps_rel / 3.7 + 2.51 / (re * f.sqrt())).log10();
            f = 1.0 / (rhs * rhs);
        }
        f
    }

    /// `f = 64/Re` exactly across the laminar range.
    #[test]
    fn laminar_is_exactly_64_over_re() {
        for &re in &[100.0_f64, 500.0, 1000.0, 1500.0, 2000.0, 2299.0] {
            let f = laminar_friction_factor(re).unwrap();
            assert!((f - 64.0 / re).abs() < 1e-12, "Re={re}");
        }
    }

    /// At the laminar limit Re=2300, f = 64/2300.
    #[test]
    fn laminar_at_2300() {
        let f = laminar_friction_factor(2300.0).unwrap();
        assert!((f - 64.0 / 2300.0).abs() < 1e-12);
        assert!((f - 0.027_826_086_956_521_74).abs() < 1e-12);
    }

    /// Haaland stays within a few percent of Colebrook-White at a
    /// canonical turbulent test point (Re=1e5, eps/D=1e-4).
    #[test]
    fn haaland_within_a_few_percent_of_colebrook() {
        let re = 1.0e5;
        let eps_rel = 1.0e-4;
        let f_haaland = haaland_friction_factor(re, eps_rel).unwrap();
        let f_colebrook = colebrook_white(re, eps_rel);
        let rel_err = (f_haaland - f_colebrook).abs() / f_colebrook;
        // Haaland's stated accuracy is within ~2% of Colebrook.
        assert!(
            rel_err < 0.02,
            "rel_err={rel_err}, haaland={f_haaland}, colebrook={f_colebrook}"
        );
    }

    /// The agreement holds across a sweep of Reynolds numbers and
    /// roughnesses spanning the turbulent regime.
    #[test]
    fn haaland_tracks_colebrook_across_regime() {
        let cases = [
            (4000.0, 0.0),
            (1.0e4, 0.0),
            (1.0e5, 1.0e-4),
            (1.0e6, 1.0e-3),
            (1.0e7, 1.0e-2),
            (5.0e5, 5.0e-3),
        ];
        for &(re, eps_rel) in &cases {
            let f_h = haaland_friction_factor(re, eps_rel).unwrap();
            let f_c = colebrook_white(re, eps_rel);
            let rel_err = (f_h - f_c).abs() / f_c;
            assert!(
                rel_err < 0.03,
                "Re={re}, eps_rel={eps_rel}: rel_err={rel_err}"
            );
        }
    }

    /// In the fully rough limit, friction becomes independent of Re
    /// (the von Karman rough-pipe asymptote): doubling Re barely moves f.
    #[test]
    fn fully_rough_is_reynolds_independent() {
        let eps_rel = 0.02;
        let f_lo = haaland_friction_factor(1.0e7, eps_rel).unwrap();
        let f_hi = haaland_friction_factor(1.0e8, eps_rel).unwrap();
        let rel_change = (f_hi - f_lo).abs() / f_lo;
        assert!(rel_change < 0.02, "rel_change={rel_change}");
    }

    /// Rougher pipe gives a larger friction factor at fixed Re.
    #[test]
    fn rougher_pipe_has_higher_friction() {
        let re = 1.0e5;
        let smooth = haaland_friction_factor(re, 0.0).unwrap();
        let rough = haaland_friction_factor(re, 1.0e-2).unwrap();
        assert!(rough > smooth);
    }

    /// Higher Reynolds number lowers the (smooth-pipe) friction factor.
    #[test]
    fn higher_re_lowers_smooth_friction() {
        let f_low_re = haaland_friction_factor(1.0e4, 0.0).unwrap();
        let f_high_re = haaland_friction_factor(1.0e6, 0.0).unwrap();
        assert!(f_high_re < f_low_re);
    }

    /// The dispatcher picks laminar below 2300 and matches `64/Re`.
    #[test]
    fn dispatch_uses_laminar_below_2300() {
        let r = friction_factor(1500.0, 0.001).unwrap();
        assert_eq!(r.regime, FlowRegime::Laminar);
        assert!((r.friction_factor - 64.0 / 1500.0).abs() < 1e-12);
    }

    /// The dispatcher picks Haaland in the turbulent range and matches
    /// the direct Haaland call.
    #[test]
    fn dispatch_uses_haaland_when_turbulent() {
        let re = 5.0e5;
        let eps_rel = 1.0e-4;
        let r = friction_factor(re, eps_rel).unwrap();
        assert_eq!(r.regime, FlowRegime::Turbulent);
        let direct = haaland_friction_factor(re, eps_rel).unwrap();
        assert!((r.friction_factor - direct).abs() < 1e-12);
        assert!((r.reynolds - re).abs() < 1e-6);
    }

    /// Bad inputs are rejected by every entry point.
    #[test]
    fn rejects_bad_inputs() {
        assert!(laminar_friction_factor(0.0).is_err());
        assert!(laminar_friction_factor(f64::NAN).is_err());
        assert!(haaland_friction_factor(-1.0, 0.0).is_err());
        // Relative roughness out of the [0, 0.05] range.
        assert!(haaland_friction_factor(1.0e5, 0.1).is_err());
        assert!(haaland_friction_factor(1.0e5, -0.01).is_err());
        assert!(friction_factor(0.0, 0.0).is_err());
    }

    /// Roughness exactly at the upper bound 0.05 is accepted.
    #[test]
    fn relative_roughness_at_upper_bound_ok() {
        assert!(haaland_friction_factor(1.0e6, 0.05).is_ok());
    }
}
