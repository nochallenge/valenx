//! Melting curve — temperature-dependent ensemble.
//!
//! As temperature rises, base pairs break: an RNA "melts". The
//! melting curve plots an order parameter — here the *expected
//! fraction of paired bases* — against temperature. The midpoint of
//! the transition is the melting temperature `Tm`.
//!
//! For each temperature the ensemble is summarised by a McCaskill
//! partition function ([`crate::ensemble::partition`]); the fraction
//! paired is `(Σᵢ p_paired(i)) / n` from that ensemble.
//!
//! ## v1 scope
//!
//! The Turner parameters are referenced at 37 °C and are *not*
//! re-fitted at other temperatures — only the Boltzmann factor
//! `exp(-E/RT)` carries the temperature. This is the standard
//! first-order treatment for a melting-curve preview; an absolute
//! `Tm` should be read as approximate.

use crate::ensemble::partition::partition_function_at;
use crate::error::{Result, RnaStructError};
use crate::fold::constraint::FoldConstraints;
use crate::rna::RnaSeq;

/// Conversion offset between Celsius and Kelvin.
pub const KELVIN_OFFSET: f64 = 273.15;

/// A single point on a melting curve.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MeltPoint {
    /// Temperature in degrees Celsius.
    pub temperature_c: f64,
    /// The expected fraction of bases that are paired, in `[0, 1]`.
    pub fraction_paired: f64,
}

/// A computed melting curve.
#[derive(Clone, Debug, PartialEq)]
pub struct MeltingCurve {
    /// The curve points, ordered by ascending temperature.
    pub points: Vec<MeltPoint>,
}

impl MeltingCurve {
    /// Estimates the melting temperature `Tm` — the temperature at
    /// which the paired fraction crosses the midpoint between its
    /// low-temperature and high-temperature plateaus.
    ///
    /// Returns `None` if the curve is empty or never crosses the
    /// midpoint (e.g. an RNA that does not pair at any temperature).
    pub fn tm(&self) -> Option<f64> {
        if self.points.len() < 2 {
            return None;
        }
        let first = self.points.first()?.fraction_paired;
        let last = self.points.last()?.fraction_paired;
        let mid = 0.5 * (first + last);
        // Find the first descending crossing of `mid`.
        for w in self.points.windows(2) {
            let (a, b) = (w[0], w[1]);
            let hi = a.fraction_paired.max(b.fraction_paired);
            let lo = a.fraction_paired.min(b.fraction_paired);
            if lo <= mid && mid <= hi && (a.fraction_paired - b.fraction_paired).abs() > 1e-12 {
                // linear interpolation
                let t = (mid - a.fraction_paired) / (b.fraction_paired - a.fraction_paired);
                return Some(a.temperature_c + t * (b.temperature_c - a.temperature_c));
            }
        }
        None
    }
}

/// Computes a melting curve over `[min_c, max_c]` °C in `steps`
/// equal intervals (`steps + 1` points).
///
/// # Errors
/// [`RnaStructError::Invalid`] if `steps == 0`, `min_c >= max_c`, or
/// either temperature is at or below absolute zero.
pub fn melting_curve(seq: &RnaSeq, min_c: f64, max_c: f64, steps: usize) -> Result<MeltingCurve> {
    if steps == 0 {
        return Err(RnaStructError::invalid("steps", "need at least one step"));
    }
    if min_c.is_nan() || max_c.is_nan() || min_c >= max_c {
        return Err(RnaStructError::invalid(
            "range",
            "min temperature must be below max temperature",
        ));
    }
    if min_c + KELVIN_OFFSET <= 0.0 || max_c + KELVIN_OFFSET <= 0.0 {
        return Err(RnaStructError::invalid(
            "temperature",
            "temperatures must be above absolute zero",
        ));
    }
    let n = seq.len();
    let cons = FoldConstraints::none(n);
    let mut points = Vec::with_capacity(steps + 1);
    for s in 0..=steps {
        let temperature_c = min_c + (max_c - min_c) * (s as f64 / steps as f64);
        let temperature_k = temperature_c + KELVIN_OFFSET;
        let pf = partition_function_at(seq, temperature_k, &cons)?;
        let mut paired = 0.0;
        for i in 0..n {
            paired += 1.0 - pf.unpaired_probability(i);
        }
        let fraction_paired = if n > 0 { paired / n as f64 } else { 0.0 };
        points.push(MeltPoint {
            temperature_c,
            fraction_paired: fraction_paired.clamp(0.0, 1.0),
        });
    }
    Ok(MeltingCurve { points })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn curve_has_the_right_number_of_points() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let c = melting_curve(&seq, 10.0, 90.0, 8).unwrap();
        assert_eq!(c.points.len(), 9);
    }

    #[test]
    fn pairing_decreases_with_temperature() {
        // A stable stem should be more paired when cold than when hot.
        let seq = RnaSeq::parse("GGGGGGGAAAACCCCCCC").unwrap();
        let c = melting_curve(&seq, 0.0, 100.0, 10).unwrap();
        let cold = c.points.first().unwrap().fraction_paired;
        let hot = c.points.last().unwrap().fraction_paired;
        assert!(
            cold >= hot,
            "RNA should be at least as paired when cold ({cold}) as hot ({hot})"
        );
    }

    #[test]
    fn tm_lies_within_the_scanned_range() {
        let seq = RnaSeq::parse("GGGGGGGAAAACCCCCCC").unwrap();
        let c = melting_curve(&seq, 0.0, 120.0, 24).unwrap();
        if let Some(tm) = c.tm() {
            assert!((0.0..=120.0).contains(&tm), "Tm {tm} outside scan range");
        }
    }

    #[test]
    fn unpairable_rna_has_flat_curve() {
        let seq = RnaSeq::parse("AAAAAAAAAAAA").unwrap();
        let c = melting_curve(&seq, 10.0, 90.0, 6).unwrap();
        for p in &c.points {
            assert!(p.fraction_paired < 1e-6, "should never pair");
        }
        assert!(c.tm().is_none());
    }

    #[test]
    fn rejects_bad_ranges() {
        let seq = RnaSeq::parse("GGGGCCCC").unwrap();
        assert!(melting_curve(&seq, 90.0, 10.0, 5).is_err());
        assert!(melting_curve(&seq, 10.0, 90.0, 0).is_err());
        assert!(melting_curve(&seq, -300.0, -280.0, 5).is_err());
    }
}
