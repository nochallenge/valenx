//! The transistor device model and the beta current relations.
//!
//! A [`Transistor`] bundles the three numbers the DC hand-analysis
//! needs: the forward current gain `beta`, the base-emitter turn-on
//! voltage `VBE`, and the collector-emitter saturation floor
//! `Vce_sat`. From `beta` alone the terminal currents are related by
//!
//! > `Ic = beta * Ib`,  `Ie = (beta + 1) * Ib`,  `Ie = Ic + Ib`.
//!
//! See the [crate-level docs](crate) for the modelling assumptions.

use crate::error::BjtError;
use serde::{Deserialize, Serialize};

/// A bipolar junction transistor described by its DC hand-analysis
/// parameters.
///
/// All voltages are in volts. `beta` is dimensionless. Construct with
/// [`Transistor::new`] (validated) or [`Transistor::silicon`] for a
/// silicon default (`VBE = 0.7 V`, `Vce_sat = 0.2 V`).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Transistor {
    /// Forward DC current gain `beta = Ic / Ib` (also written `hFE`).
    pub beta: f64,
    /// Base-emitter turn-on voltage drop while conducting, in volts
    /// (≈ 0.7 V for silicon, ≈ 0.3 V for germanium).
    pub vbe: f64,
    /// Collector-emitter voltage at the edge of saturation, in volts
    /// (≈ 0.2 V for a small-signal silicon BJT). The device is treated
    /// as active while `Vce > vce_sat`.
    pub vce_sat: f64,
}

impl Transistor {
    /// Build a validated transistor.
    ///
    /// # Errors
    ///
    /// Returns [`BjtError::BadParameter`] if `beta <= 0`, if `vbe < 0`,
    /// if `vce_sat < 0`, or if any argument is not finite.
    pub fn new(beta: f64, vbe: f64, vce_sat: f64) -> Result<Self, BjtError> {
        if !beta.is_finite() || beta <= 0.0 {
            return Err(BjtError::bad_parameter(
                "beta",
                "current gain must be finite and strictly positive",
                beta,
            ));
        }
        if !vbe.is_finite() || vbe < 0.0 {
            return Err(BjtError::bad_parameter(
                "vbe",
                "turn-on voltage must be finite and non-negative",
                vbe,
            ));
        }
        if !vce_sat.is_finite() || vce_sat < 0.0 {
            return Err(BjtError::bad_parameter(
                "vce_sat",
                "saturation voltage must be finite and non-negative",
                vce_sat,
            ));
        }
        Ok(Self { beta, vbe, vce_sat })
    }

    /// A silicon BJT with the textbook defaults `VBE = 0.7 V` and
    /// `Vce_sat = 0.2 V`, parameterised only by its gain `beta`.
    ///
    /// # Errors
    ///
    /// Returns [`BjtError::BadParameter`] if `beta` is not finite and
    /// strictly positive.
    pub fn silicon(beta: f64) -> Result<Self, BjtError> {
        Self::new(beta, 0.7, 0.2)
    }

    /// Collector current from base current: `Ic = beta * Ib`.
    ///
    /// The sign of `ib` is preserved; for a forward-active NPN bias the
    /// caller passes a non-negative base current.
    pub fn collector_current(&self, ib: f64) -> f64 {
        self.beta * ib
    }

    /// Emitter current from base current: `Ie = (beta + 1) * Ib`.
    pub fn emitter_current(&self, ib: f64) -> f64 {
        (self.beta + 1.0) * ib
    }

    /// Base current that yields a given collector current:
    /// `Ib = Ic / beta`.
    pub fn base_current_for_collector(&self, ic: f64) -> f64 {
        ic / self.beta
    }

    /// Recover the gain from a measured collector / base current pair:
    /// `beta = Ic / Ib`.
    ///
    /// # Errors
    ///
    /// Returns [`BjtError::BadParameter`] if `ib` is zero or not finite
    /// (the ratio would be undefined).
    pub fn beta_from_currents(ic: f64, ib: f64) -> Result<f64, BjtError> {
        if !ib.is_finite() || ib == 0.0 {
            return Err(BjtError::bad_parameter(
                "ib",
                "base current must be finite and non-zero to define beta",
                ib,
            ));
        }
        Ok(ic / ib)
    }
}

/// Which DC operating region the transistor sits in.
///
/// This crate distinguishes the two regions relevant to a conducting
/// amplifier/switch bias point. Cut-off (no base drive) is surfaced as
/// the [`BjtError::CutOff`] error from the bias solvers rather than as a
/// region here, because a cut-off device has no meaningful collector
/// current to report.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum Region {
    /// Forward-active: `Vce > Vce_sat`. The collector current follows
    /// `Ic = beta * Ib` and the device acts as a (roughly) linear
    /// current source — the amplifier region.
    Active,
    /// Saturation: the active-region collector current would pull `Vce`
    /// to or below `Vce_sat`, so `Vce` is pinned at the floor and the
    /// collector current is limited by the external resistors rather
    /// than by `beta` — the closed-switch region.
    Saturation,
}
