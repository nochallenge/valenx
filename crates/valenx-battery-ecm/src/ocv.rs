//! Open-circuit-voltage (OCV) versus state-of-charge (SoC) lookup table.
//!
//! ## Model
//!
//! A real cell's *rest* (no-load) terminal voltage is a monotone
//! function of its state of charge: `OCV = f(SoC)`. There is no
//! closed-form expression for a given chemistry, so the relationship is
//! supplied as a small table of `(SoC, OCV)` breakpoints — typically
//! measured by slowly (de)charging a cell and letting it relax — and
//! evaluated by **piecewise-linear interpolation** between the nearest
//! two breakpoints.
//!
//! [`OcvSocTable`] enforces the two physical invariants that make the
//! interpolation well-defined and monotone:
//!
//! - the SoC breakpoints are **strictly increasing** (so every query
//!   lands in exactly one bracket), and
//! - the OCV column is **monotonically non-decreasing** with SoC (the
//!   rest voltage of a real cell never falls as it charges).
//!
//! Queries are **clamped** to the table's SoC span: a SoC below the
//! first breakpoint returns the first OCV, a SoC above the last returns
//! the last. Because the table is monotone, the interpolated OCV is
//! itself a monotone, continuous function of SoC.

use crate::error::{EcmError, Result};
use serde::{Deserialize, Serialize};

/// Absolute tolerance used when checking table monotonicity.
///
/// A tiny negative step (within this band) is treated as flat rather
/// than as a monotonicity violation, so floating-point round-off in a
/// hand-entered table does not spuriously reject an otherwise-valid
/// curve.
const MONOTONIC_EPS: f64 = 1e-12;

/// A monotone open-circuit-voltage / state-of-charge lookup table.
///
/// Build one with [`OcvSocTable::new`]; evaluate it with
/// [`OcvSocTable::ocv_at`]. The breakpoints are stored sorted and
/// validated, so evaluation is infallible.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OcvSocTable {
    /// Strictly-increasing state-of-charge breakpoints, each in `[0, 1]`.
    soc: Vec<f64>,
    /// Open-circuit voltage at each breakpoint, monotonically
    /// non-decreasing and parallel to `soc`.
    ocv: Vec<f64>,
}

impl OcvSocTable {
    /// Build a validated table from parallel `(SoC, OCV)` columns.
    ///
    /// # Errors
    ///
    /// Returns [`EcmError::Table`] if there are fewer than two points,
    /// the two columns differ in length, the SoC breakpoints are not
    /// strictly increasing, or the OCV column is not monotonically
    /// non-decreasing. Returns [`EcmError::SocOutOfRange`] if any SoC
    /// breakpoint lies outside `[0, 1]`, and [`EcmError::Invalid`] if
    /// any value is non-finite.
    pub fn new(soc: Vec<f64>, ocv: Vec<f64>) -> Result<Self> {
        if soc.len() != ocv.len() {
            return Err(EcmError::table(format!(
                "SoC column has {} points but OCV column has {}",
                soc.len(),
                ocv.len()
            )));
        }
        if soc.len() < 2 {
            return Err(EcmError::table(format!(
                "need at least 2 breakpoints, got {}",
                soc.len()
            )));
        }
        for (i, (&s, &v)) in soc.iter().zip(ocv.iter()).enumerate() {
            if !s.is_finite() {
                return Err(EcmError::invalid(
                    "ocv_table.soc",
                    format!("breakpoint {i} is not finite"),
                ));
            }
            if !v.is_finite() {
                return Err(EcmError::invalid(
                    "ocv_table.ocv",
                    format!("value {i} is not finite"),
                ));
            }
            if !(0.0..=1.0).contains(&s) {
                return Err(EcmError::SocOutOfRange {
                    soc: s,
                    context: "OCV breakpoint",
                });
            }
        }
        for i in 1..soc.len() {
            if soc[i] <= soc[i - 1] {
                return Err(EcmError::table(format!(
                    "SoC breakpoints must strictly increase: soc[{}]={} <= soc[{}]={}",
                    i,
                    soc[i],
                    i - 1,
                    soc[i - 1]
                )));
            }
            if ocv[i] < ocv[i - 1] - MONOTONIC_EPS {
                return Err(EcmError::table(format!(
                    "OCV must be non-decreasing in SoC: ocv[{}]={} < ocv[{}]={}",
                    i,
                    ocv[i],
                    i - 1,
                    ocv[i - 1]
                )));
            }
        }
        Ok(Self { soc, ocv })
    }

    /// Number of breakpoints in the table.
    pub fn len(&self) -> usize {
        self.soc.len()
    }

    /// Always `false` — a valid table has at least two breakpoints. This
    /// method exists only to satisfy the `len`/`is_empty` clippy lint.
    pub fn is_empty(&self) -> bool {
        self.soc.is_empty()
    }

    /// The lowest tabulated state of charge.
    pub fn min_soc(&self) -> f64 {
        self.soc[0]
    }

    /// The highest tabulated state of charge.
    pub fn max_soc(&self) -> f64 {
        self.soc[self.soc.len() - 1]
    }

    /// Open-circuit voltage at a given state of charge.
    ///
    /// Evaluated by piecewise-linear interpolation between the bracketing
    /// breakpoints. A `soc` below [`min_soc`](Self::min_soc) returns the
    /// first OCV and a `soc` above [`max_soc`](Self::max_soc) returns the
    /// last (the query is **clamped** to the table span). Because the
    /// table is monotone, the result is a monotone continuous function of
    /// `soc`.
    pub fn ocv_at(&self, soc: f64) -> f64 {
        let n = self.soc.len();
        if soc <= self.soc[0] {
            return self.ocv[0];
        }
        if soc >= self.soc[n - 1] {
            return self.ocv[n - 1];
        }
        // Find the bracket [i-1, i] with soc[i-1] < soc <= soc[i].
        // `partition_point` returns the first index whose breakpoint is
        // strictly greater than `soc`; it is in `1..n` here because the
        // two clamping guards above have excluded the ends.
        let i = self.soc.partition_point(|&s| s <= soc);
        let (s0, s1) = (self.soc[i - 1], self.soc[i]);
        let (v0, v1) = (self.ocv[i - 1], self.ocv[i]);
        let t = (soc - s0) / (s1 - s0);
        v0 + t * (v1 - v0)
    }

    /// Borrow the SoC breakpoint column.
    pub fn soc_points(&self) -> &[f64] {
        &self.soc
    }

    /// Borrow the OCV column.
    pub fn ocv_points(&self) -> &[f64] {
        &self.ocv
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A simple, valid 5-point curve rising from 3.0 V to 4.2 V.
    fn sample() -> OcvSocTable {
        OcvSocTable::new(
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            vec![3.0, 3.5, 3.7, 3.9, 4.2],
        )
        .unwrap()
    }

    #[test]
    fn exact_breakpoints_recovered() {
        let t = sample();
        for (&s, &v) in t.soc_points().iter().zip(t.ocv_points()) {
            assert!((t.ocv_at(s) - v).abs() < 1e-12, "at soc={s}");
        }
    }

    #[test]
    fn linear_midpoint_interpolation() {
        let t = sample();
        // Midway between (0.0, 3.0) and (0.25, 3.5) -> 3.25 V.
        let got = t.ocv_at(0.125);
        assert!((got - 3.25).abs() < 1e-12, "got {got}");
        // Midway between (0.5, 3.7) and (0.75, 3.9) -> 3.8 V.
        let got = t.ocv_at(0.625);
        assert!((got - 3.8).abs() < 1e-12, "got {got}");
    }

    #[test]
    fn queries_clamp_to_span() {
        let t = sample();
        assert!((t.ocv_at(-1.0) - 3.0).abs() < 1e-12);
        assert!((t.ocv_at(0.0) - 3.0).abs() < 1e-12);
        assert!((t.ocv_at(2.0) - 4.2).abs() < 1e-12);
        assert!((t.ocv_at(1.0) - 4.2).abs() < 1e-12);
    }

    #[test]
    fn interpolation_is_monotone_non_decreasing() {
        let t = sample();
        let mut prev = f64::NEG_INFINITY;
        // Sweep SoC across and beyond the span; OCV never decreases.
        for k in -10..=110 {
            let soc = k as f64 / 100.0;
            let v = t.ocv_at(soc);
            assert!(v >= prev - 1e-12, "non-monotone at soc={soc}: {v} < {prev}");
            prev = v;
        }
    }

    #[test]
    fn rejects_too_few_points() {
        let err = OcvSocTable::new(vec![0.5], vec![3.7]).unwrap_err();
        assert_eq!(err.code(), "battery_ecm.table");
    }

    #[test]
    fn rejects_length_mismatch() {
        let err = OcvSocTable::new(vec![0.0, 1.0], vec![3.0]).unwrap_err();
        assert_eq!(err.code(), "battery_ecm.table");
    }

    #[test]
    fn rejects_non_increasing_soc() {
        let err = OcvSocTable::new(vec![0.0, 0.5, 0.5], vec![3.0, 3.5, 3.7]).unwrap_err();
        assert_eq!(err.code(), "battery_ecm.table");
        let err = OcvSocTable::new(vec![0.0, 0.6, 0.5], vec![3.0, 3.5, 3.7]).unwrap_err();
        assert_eq!(err.code(), "battery_ecm.table");
    }

    #[test]
    fn rejects_decreasing_ocv() {
        let err = OcvSocTable::new(vec![0.0, 0.5, 1.0], vec![3.0, 3.5, 3.4]).unwrap_err();
        assert_eq!(err.code(), "battery_ecm.table");
    }

    #[test]
    fn rejects_soc_out_of_unit_interval() {
        let err = OcvSocTable::new(vec![-0.1, 0.5, 1.0], vec![3.0, 3.5, 4.2]).unwrap_err();
        assert_eq!(err.code(), "battery_ecm.soc_out_of_range");
        let err = OcvSocTable::new(vec![0.0, 0.5, 1.1], vec![3.0, 3.5, 4.2]).unwrap_err();
        assert_eq!(err.code(), "battery_ecm.soc_out_of_range");
    }

    #[test]
    fn flat_ocv_segment_allowed() {
        // A plateau (constant OCV across a SoC span) is physical and
        // must be accepted: non-decreasing, not strictly increasing.
        let t = OcvSocTable::new(vec![0.0, 0.5, 1.0], vec![3.5, 3.7, 3.7]).unwrap();
        assert!((t.ocv_at(0.75) - 3.7).abs() < 1e-12);
    }

    #[test]
    fn span_accessors() {
        let t = sample();
        assert!((t.min_soc() - 0.0).abs() < 1e-12);
        assert!((t.max_soc() - 1.0).abs() < 1e-12);
        assert_eq!(t.len(), 5);
        assert!(!t.is_empty());
    }
}
