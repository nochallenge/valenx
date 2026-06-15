//! Area and tube-count sizing (the LMTD design method).
//!
//! Given a thermal duty `Q`, an overall heat-transfer coefficient `U`, a
//! correction factor `F` and the terminal temperatures, the required
//! heat-transfer surface area is
//!
//! ```text
//! A = Q / (U * F * LMTD)
//! ```
//!
//! Once the area is known, the number of plain tubes needed for a chosen
//! tube outside diameter `d` and effective length `L` follows from the
//! per-tube cylindrical surface `pi * d * L`:
//!
//! ```text
//! n_tubes_real = A / (pi * d * L)
//! n_tubes      = ceil(n_tubes_real)   // a whole-tube bundle
//! ```
//!
//! All quantities are SI: `Q` in watts, `U` in `W/(m^2 K)`, temperatures
//! and their differences in kelvin, lengths in metres, area in `m^2`.

use crate::correction::CorrectionFactor;
use crate::error::HxError;
use crate::lmtd::TerminalDeltas;
use std::f64::consts::PI;

/// Inputs for an LMTD-method sizing calculation.
///
/// Build one with [`SizingInput::new`], which validates every field, then
/// call [`size`] to obtain the [`SizingResult`].
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SizingInput {
    /// Thermal duty / heat-transfer rate `Q` (watts), strictly positive.
    pub duty_w: f64,
    /// Overall heat-transfer coefficient `U` (`W/(m^2 K)`), strictly
    /// positive. Referenced to the same area being solved for.
    pub u_w_per_m2k: f64,
    /// LMTD correction factor `F` in `(0, 1]`.
    pub correction: CorrectionFactor,
    /// The two terminal temperature differences.
    pub deltas: TerminalDeltas,
}

impl SizingInput {
    /// Validate and assemble a sizing input from raw scalars.
    ///
    /// `f` is validated into a [`CorrectionFactor`] and the terminal
    /// differences into a [`TerminalDeltas`]; `duty_w` and `u_w_per_m2k`
    /// must be finite and strictly positive.
    ///
    /// # Errors
    ///
    /// Propagates [`HxError`] from any of the field validations:
    /// [`HxError::BadParameter`] for a bad duty or `U`,
    /// [`HxError::CorrectionFactorOutOfRange`] for a bad `F`, and
    /// [`HxError::InfeasibleTemperatureProfile`] for a temperature cross.
    pub fn new(
        duty_w: f64,
        u_w_per_m2k: f64,
        f: f64,
        deltas: TerminalDeltas,
    ) -> Result<Self, HxError> {
        let duty_w = HxError::require_positive("duty_w", duty_w)?;
        let u_w_per_m2k = HxError::require_positive("u_w_per_m2k", u_w_per_m2k)?;
        let correction = CorrectionFactor::new(f)?;
        Ok(Self {
            duty_w,
            u_w_per_m2k,
            correction,
            deltas,
        })
    }
}

/// A chosen tube geometry: outside diameter and effective length, both in
/// metres and strictly positive.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct TubeGeometry {
    /// Tube outside diameter `d` (metres), strictly positive.
    pub outside_diameter_m: f64,
    /// Effective (heat-transfer) tube length `L` (metres), strictly
    /// positive.
    pub length_m: f64,
}

impl TubeGeometry {
    /// Validate and build a tube geometry.
    ///
    /// # Errors
    ///
    /// Returns [`HxError::BadParameter`] if either dimension is
    /// non-finite, zero or negative.
    pub fn new(outside_diameter_m: f64, length_m: f64) -> Result<Self, HxError> {
        let outside_diameter_m =
            HxError::require_positive("outside_diameter_m", outside_diameter_m)?;
        let length_m = HxError::require_positive("length_m", length_m)?;
        Ok(Self {
            outside_diameter_m,
            length_m,
        })
    }

    /// Heat-transfer surface area of a single tube, `pi * d * L` (`m^2`).
    pub fn area_per_tube_m2(&self) -> f64 {
        PI * self.outside_diameter_m * self.length_m
    }

    /// Continuous (fractional) number of tubes required to provide
    /// `total_area_m2` of surface, `A / (pi * d * L)`.
    ///
    /// This is the un-rounded value; [`SizingResult::tube_count`] rounds
    /// it up to a whole bundle.
    pub fn tubes_for_area(&self, total_area_m2: f64) -> f64 {
        total_area_m2 / self.area_per_tube_m2()
    }
}

/// Result of an LMTD-method sizing calculation.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SizingResult {
    /// Raw counter-current log-mean temperature difference (kelvin).
    pub lmtd_k: f64,
    /// Effective `F`-corrected mean driving temperature, `F * LMTD`
    /// (kelvin). Always `<=` [`SizingResult::lmtd_k`].
    pub effective_lmtd_k: f64,
    /// Required heat-transfer surface area `A` (`m^2`).
    pub area_m2: f64,
}

impl SizingResult {
    /// Continuous number of tubes for a given [`TubeGeometry`],
    /// `A / (pi * d * L)` (not yet rounded).
    pub fn tubes_real(&self, tube: &TubeGeometry) -> f64 {
        tube.tubes_for_area(self.area_m2)
    }

    /// Whole-tube count for a given [`TubeGeometry`]: the continuous tube
    /// requirement rounded **up** so the bundle delivers at least the
    /// required area.
    pub fn tube_count(&self, tube: &TubeGeometry) -> u64 {
        // area_m2 > 0 and area_per_tube > 0, so tubes_real > 0 and ceil
        // is >= 1; cast is safe for any realistic bundle size.
        self.tubes_real(tube).ceil() as u64
    }
}

/// Solve the LMTD-method sizing for the required surface area.
///
/// Returns the raw LMTD, the `F`-corrected effective LMTD, and the
/// required area `A = Q / (U F LMTD)`. The result never includes a tube
/// count — that depends on a chosen [`TubeGeometry`] and is obtained from
/// [`SizingResult::tube_count`].
///
/// ```
/// use valenx_shellandtube::{size, SizingInput, TerminalDeltas};
///
/// let deltas = TerminalDeltas::new(30.0, 10.0).unwrap();
/// let input = SizingInput::new(100_000.0, 500.0, 0.9, deltas).unwrap();
/// let r = size(&input);
/// // LMTD = 20 / ln 3 ≈ 18.2048 K, effective ≈ 16.3843 K.
/// assert!((r.lmtd_k - 18.204_784_532_536_746).abs() < 1e-9);
/// assert!(r.area_m2 > 0.0);
/// ```
pub fn size(input: &SizingInput) -> SizingResult {
    let lmtd_k = input.deltas.lmtd();
    let effective_lmtd_k = input.correction.effective_lmtd(lmtd_k);
    // All factors are strictly positive (validated on construction), so
    // the denominator is strictly positive and the area is finite.
    let area_m2 = input.duty_w / (input.u_w_per_m2k * effective_lmtd_k);
    SizingResult {
        lmtd_k,
        effective_lmtd_k,
        area_m2,
    }
}
