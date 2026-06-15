//! The LMTD correction factor `F`.
//!
//! A pure counter-current exchanger is the thermodynamic best case: it
//! extracts the full log-mean driving temperature. Real shell-and-tube
//! units mix co-current and counter-current passes, so the *effective*
//! mean driving temperature is reduced by a dimensionless correction
//! factor `F`:
//!
//! ```text
//! Q = U A F LMTD_counter-current
//! ```
//!
//! By construction `F` lies in the half-open interval `(0, 1]`:
//!
//! - `F = 1` recovers ideal counter-current (or any phase-change stream,
//!   where the configuration no longer matters).
//! - `0 < F < 1` is the realistic multi-pass case; lower `F` means a
//!   larger area is required for the same duty.
//! - `F <= 0` is unphysical and a design red flag.
//!
//! This crate does **not** read `F` off the classic `P`/`R` Bowman charts;
//! it treats `F` as a validated user input ([`CorrectionFactor`]) so the
//! sizing relations stay exact and chart-source-independent.

use crate::error::HxError;

/// A validated LMTD correction factor in the interval `(0, 1]`.
///
/// Construct it with [`CorrectionFactor::new`] (or the convenience
/// [`CorrectionFactor::ideal`]); the wrapped value is then guaranteed to
/// be a finite number strictly greater than zero and at most one, which
/// the sizing functions rely on.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CorrectionFactor(f64);

impl CorrectionFactor {
    /// Validate and wrap a correction factor.
    ///
    /// # Errors
    ///
    /// Returns [`HxError::CorrectionFactorOutOfRange`] if `f` is
    /// non-finite, `<= 0`, or `> 1`.
    pub fn new(f: f64) -> Result<Self, HxError> {
        if !f.is_finite() {
            return Err(HxError::CorrectionFactorOutOfRange {
                value: f,
                reason: "must be a finite number".to_string(),
            });
        }
        if f <= 0.0 {
            return Err(HxError::CorrectionFactorOutOfRange {
                value: f,
                reason: "must be > 0 (an exchanger with F <= 0 is unphysical)".to_string(),
            });
        }
        if f > 1.0 {
            return Err(HxError::CorrectionFactorOutOfRange {
                value: f,
                reason: "must be <= 1 (F = 1 is the counter-current ideal)".to_string(),
            });
        }
        Ok(Self(f))
    }

    /// The ideal counter-current correction factor, `F = 1`.
    ///
    /// Useful as a baseline and for the phase-change case where the
    /// multi-pass geometry no longer penalises the driving temperature.
    pub fn ideal() -> Self {
        Self(1.0)
    }

    /// The wrapped factor as a plain `f64` in `(0, 1]`.
    pub fn value(self) -> f64 {
        self.0
    }

    /// The **effective** (F-corrected) log-mean temperature difference,
    /// `F * lmtd`, in kelvin.
    ///
    /// This is the mean driving temperature the exchanger actually sees,
    /// and is always `<=` the raw `lmtd` because `F <= 1`.
    pub fn effective_lmtd(self, lmtd: f64) -> f64 {
        self.0 * lmtd
    }
}
