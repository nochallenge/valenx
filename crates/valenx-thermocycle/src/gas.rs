//! The working-fluid heat-capacity ratio `γ = c_p / c_v`.
//!
//! The air-standard Otto, Diesel and Brayton efficiencies all depend on
//! the ratio of specific heats of the working fluid, conventionally
//! written `γ` (gamma) and sometimes `k`. This module wraps it in a
//! validated [`HeatCapacityRatio`] newtype so an out-of-range value is
//! rejected once, at construction, rather than producing a silent
//! non-physical efficiency downstream.
//!
//! A real (classically modelled) gas always has `γ > 1`: the constant-
//! pressure specific heat exceeds the constant-volume one by exactly the
//! specific gas constant, `c_p = c_v + R`. Monatomic gases sit at
//! `γ = 5/3 ≈ 1.667`, diatomic gases (and air near room temperature) at
//! `γ = 7/5 = 1.4`.

use crate::error::{CycleError, Result};
use serde::{Deserialize, Serialize};

/// A validated heat-capacity ratio `γ = c_p / c_v`, guaranteed to be a
/// finite value strictly greater than `1`.
///
/// Construct one with [`HeatCapacityRatio::new`] (validating) or via the
/// named presets [`HeatCapacityRatio::air`],
/// [`HeatCapacityRatio::monatomic`] and [`HeatCapacityRatio::diatomic`].
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct HeatCapacityRatio(f64);

impl HeatCapacityRatio {
    /// The air-standard value used throughout introductory thermodynamics,
    /// `γ = 1.4`. Air is predominantly diatomic (N₂, O₂) and near room
    /// temperature behaves as an ideal diatomic gas.
    pub const AIR: f64 = 1.4;

    /// Ideal monatomic-gas ratio, `γ = 5/3 ≈ 1.6667` (e.g. helium, argon).
    pub const MONATOMIC: f64 = 5.0 / 3.0;

    /// Ideal diatomic-gas ratio, `γ = 7/5 = 1.4` (e.g. N₂, O₂, H₂ near
    /// room temperature). Numerically equal to [`Self::AIR`].
    pub const DIATOMIC: f64 = 7.0 / 5.0;

    /// Validate and wrap a raw `γ`.
    ///
    /// # Errors
    ///
    /// Returns [`CycleError::NotFinite`] if `gamma` is `NaN` / infinite,
    /// or [`CycleError::GammaTooLow`] if `gamma <= 1`.
    pub fn new(gamma: f64) -> Result<Self> {
        if !gamma.is_finite() {
            return Err(CycleError::NotFinite {
                name: "gamma",
                value: gamma,
            });
        }
        if gamma <= 1.0 {
            return Err(CycleError::GammaTooLow { value: gamma });
        }
        Ok(Self(gamma))
    }

    /// The air-standard ratio, `γ = 1.4`.
    ///
    /// Infallible: the constant is in range by construction.
    pub fn air() -> Self {
        Self(Self::AIR)
    }

    /// The ideal monatomic-gas ratio, `γ = 5/3`.
    pub fn monatomic() -> Self {
        Self(Self::MONATOMIC)
    }

    /// The ideal diatomic-gas ratio, `γ = 7/5`.
    pub fn diatomic() -> Self {
        Self(Self::DIATOMIC)
    }

    /// The wrapped numeric value of `γ`.
    pub fn value(self) -> f64 {
        self.0
    }
}

impl Default for HeatCapacityRatio {
    /// Defaults to the air-standard value, `γ = 1.4`.
    fn default() -> Self {
        Self::air()
    }
}
