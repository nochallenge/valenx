//! Log-mean temperature difference (LMTD).
//!
//! For a two-stream heat exchanger with terminal temperature differences
//! `dt1` and `dt2` (the hot-minus-cold gap at each end of the unit), the
//! log-mean temperature difference is
//!
//! ```text
//! LMTD = (dt1 - dt2) / ln(dt1 / dt2)
//! ```
//!
//! which is the correct mean driving temperature for the integrated
//! `Q = U A LMTD` relation under the usual textbook assumptions (constant
//! `U`, constant specific heats, no phase change, negligible losses).
//!
//! When `dt1 == dt2` the closed form is the indeterminate `0/0`; the
//! analytic limit is simply the common value, so this module returns
//! `dt1` in that case (and stays numerically stable in a neighbourhood of
//! it via the same exact limit).
//!
//! ## Sign convention
//!
//! Both terminal differences must be strictly positive. A non-positive
//! terminal difference means the streams have crossed (a temperature
//! cross / pinch), which the single-mean LMTD form cannot represent, so
//! it is reported as [`HxError::InfeasibleTemperatureProfile`].

use crate::error::HxError;

/// The two terminal temperature differences of an exchanger, in kelvin
/// (equivalently degrees Celsius, since this is a *difference*).
///
/// `dt1` and `dt2` are the hot-stream-minus-cold-stream temperature gaps
/// at the two ends of the exchanger. For a pure counter-current unit with
/// hot side `Th_in -> Th_out` and cold side `Tc_in -> Tc_out`:
///
/// ```text
/// dt1 = Th_in  - Tc_out
/// dt2 = Th_out - Tc_in
/// ```
///
/// for a co-current unit:
///
/// ```text
/// dt1 = Th_in  - Tc_in
/// dt2 = Th_out - Tc_out
/// ```
///
/// The LMTD itself is symmetric in `dt1` and `dt2`, so the labelling of
/// which end is "1" does not affect the result.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct TerminalDeltas {
    /// Temperature difference at end 1 (kelvin), strictly positive.
    pub dt1: f64,
    /// Temperature difference at end 2 (kelvin), strictly positive.
    pub dt2: f64,
}

impl TerminalDeltas {
    /// Construct and validate a pair of terminal temperature differences.
    ///
    /// # Errors
    ///
    /// Returns [`HxError::BadParameter`] if either value is non-finite,
    /// and [`HxError::InfeasibleTemperatureProfile`] if either value is
    /// non-positive (a temperature cross the LMTD form cannot model).
    pub fn new(dt1: f64, dt2: f64) -> Result<Self, HxError> {
        for (name, v) in [("dt1", dt1), ("dt2", dt2)] {
            if !v.is_finite() {
                return Err(HxError::bad(name, format!("must be finite, got {v}")));
            }
            if v <= 0.0 {
                return Err(HxError::InfeasibleTemperatureProfile(format!(
                    "terminal difference `{name}` = {v} must be > 0 (streams have crossed)"
                )));
            }
        }
        Ok(Self { dt1, dt2 })
    }

    /// Terminal differences for a **counter-current** exchanger from the
    /// four stream temperatures (kelvin).
    ///
    /// `th_in`/`th_out` are the hot-stream inlet/outlet, `tc_in`/`tc_out`
    /// the cold-stream inlet/outlet. Validation matches [`Self::new`].
    ///
    /// # Errors
    ///
    /// See [`Self::new`].
    pub fn counter_current(
        th_in: f64,
        th_out: f64,
        tc_in: f64,
        tc_out: f64,
    ) -> Result<Self, HxError> {
        Self::new(th_in - tc_out, th_out - tc_in)
    }

    /// Terminal differences for a **co-current** (parallel-flow)
    /// exchanger from the four stream temperatures (kelvin).
    ///
    /// # Errors
    ///
    /// See [`Self::new`].
    pub fn co_current(th_in: f64, th_out: f64, tc_in: f64, tc_out: f64) -> Result<Self, HxError> {
        Self::new(th_in - tc_in, th_out - tc_out)
    }

    /// The log-mean temperature difference for these terminals (kelvin).
    ///
    /// Uses the exact analytic limit `LMTD = dt1` when the two terminal
    /// differences are equal (within a tiny relative tolerance), avoiding
    /// the `0/0` of the closed form.
    pub fn lmtd(&self) -> f64 {
        lmtd(self.dt1, self.dt2)
    }
}

/// Log-mean temperature difference of two strictly-positive terminal
/// differences `dt1` and `dt2` (kelvin).
///
/// This is the free-function core used by [`TerminalDeltas::lmtd`]. It
/// assumes both arguments are already validated as finite and positive
/// (the [`TerminalDeltas`] constructors guarantee this). For equal
/// terminals it returns the analytic limit `dt1` rather than `0/0`.
///
/// ```
/// use valenx_shellandtube::lmtd::lmtd;
/// // dt1 = 30 K, dt2 = 10 K  ->  20 / ln 3 ≈ 18.2043 K.
/// let v = lmtd(30.0, 10.0);
/// assert!((v - 18.204_784_532_536_746).abs() < 1e-9);
/// ```
pub fn lmtd(dt1: f64, dt2: f64) -> f64 {
    // Relative closeness test: when the two terminal gaps coincide the
    // closed form is 0/0 and the true limit is their common value.
    let denom = dt1 - dt2;
    let scale = dt1.abs().max(dt2.abs());
    if denom.abs() <= 1e-12 * scale {
        // Both are positive and (near) equal; their mean is the limit.
        return 0.5 * (dt1 + dt2);
    }
    denom / (dt1 / dt2).ln()
}
