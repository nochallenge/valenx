//! The lossless transmission line itself: characteristic impedance
//! derived from the distributed line constants.
//!
//! ## Model
//!
//! For an ideal **lossless** line (series resistance `R = 0`, shunt
//! conductance `G = 0`) the characteristic impedance reduces to the
//! purely real
//!
//! ```text
//! Z0 = sqrt(L / C)
//! ```
//!
//! where `L` is the series inductance per unit length (H/m) and `C` is
//! the shunt capacitance per unit length (F/m). The per-unit-length
//! basis cancels, so the same ratio holds whether `L` and `C` are given
//! per metre, per foot, or as the totals for a fixed length.

use serde::{Deserialize, Serialize};

use crate::error::{ensure_positive, TlError};
use crate::reflection::{Load, Reflection};

/// A lossless transmission line, characterised by its real
/// characteristic impedance `Z0` (ohms).
///
/// Construct from distributed constants with [`Line::from_lc`], or
/// directly from a known impedance with [`Line::from_z0`] (the common
/// case of a catalogued 50 Ω or 75 Ω line).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Line {
    /// Characteristic impedance `Z0` in ohms. Always finite and `> 0`.
    z0_ohms: f64,
}

impl Line {
    /// Build a line from its distributed series inductance and shunt
    /// capacitance per unit length, via `Z0 = sqrt(L / C)`.
    ///
    /// `inductance_per_m` is `L` (H/m) and `capacitance_per_m` is `C`
    /// (F/m). Both must be finite and strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`TlError::NonPositive`] if either constant is `<= 0`, or
    /// [`TlError::NotFinite`] if either is `NaN` / `±∞`.
    ///
    /// # Examples
    ///
    /// A line with `L = 250 nH/m` and `C = 100 pF/m` has
    /// `Z0 = sqrt(250e-9 / 100e-12) = 50 Ω`:
    ///
    /// ```
    /// use valenx_transmissionline::Line;
    ///
    /// let line = Line::from_lc(250e-9, 100e-12).unwrap();
    /// assert!((line.z0_ohms() - 50.0).abs() < 1e-9);
    /// ```
    pub fn from_lc(inductance_per_m: f64, capacitance_per_m: f64) -> Result<Self, TlError> {
        let l = ensure_positive("inductance_per_m", inductance_per_m)?;
        let c = ensure_positive("capacitance_per_m", capacitance_per_m)?;
        // L, C > 0 ⇒ L / C > 0 ⇒ sqrt is finite and positive.
        Self::from_z0((l / c).sqrt())
    }

    /// Build a line directly from a known characteristic impedance
    /// `z0_ohms` (e.g. `50.0` or `75.0`).
    ///
    /// `z0_ohms` must be finite and strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`TlError::NonPositive`] if `z0_ohms <= 0`, or
    /// [`TlError::NotFinite`] if it is not finite.
    pub fn from_z0(z0_ohms: f64) -> Result<Self, TlError> {
        let z0_ohms = ensure_positive("z0_ohms", z0_ohms)?;
        Ok(Self { z0_ohms })
    }

    /// The characteristic impedance `Z0` in ohms.
    #[inline]
    pub fn z0_ohms(&self) -> f64 {
        self.z0_ohms
    }

    /// Compute the full reflection report for a purely resistive `load`
    /// terminating this line.
    ///
    /// This is the convenience entry point that ties the line to a
    /// [`Load`]: it evaluates the reflection coefficient and every
    /// derived figure of merit (VSWR, return loss, mismatch loss) in one
    /// call. See [`Reflection`] for the individual quantities.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_transmissionline::{Line, Load};
    ///
    /// let line = Line::from_z0(50.0).unwrap();
    /// let matched = line.reflection(Load::resistive(50.0).unwrap());
    /// assert!(matched.gamma().abs() < 1e-12);
    /// assert!((matched.vswr().unwrap() - 1.0).abs() < 1e-12);
    /// ```
    pub fn reflection(&self, load: Load) -> Reflection {
        Reflection::from_line_load(*self, load)
    }
}
