//! The air-standard Brayton (Joule) cycle — the ideal gas-turbine /
//! jet engine.
//!
//! The Brayton cycle idealises a gas turbine as isentropic compression,
//! constant-pressure heat addition (combustion), isentropic expansion
//! through the turbine, and constant-pressure heat rejection. Under the
//! air-standard assumptions its thermal efficiency depends only on the
//! **pressure ratio** `r_p` across the compressor and the heat-capacity
//! ratio `γ`:
//!
//! ```text
//! η_brayton = 1 - 1 / r_p^((γ - 1) / γ)
//! ```
//!
//! Efficiency rises monotonically with the pressure ratio, which is why
//! gas turbines chase ever-higher pressure ratios — again bounded above
//! by the Carnot limit for the same temperature extremes. (The *net work*
//! per unit mass, by contrast, peaks at an intermediate pressure ratio,
//! but that work optimum is outside the scope of this efficiency-only
//! crate.)

use crate::error::{CycleError, Result};
use crate::gas::HeatCapacityRatio;
use serde::{Deserialize, Serialize};

/// A validated air-standard Brayton cycle: a compressor pressure ratio
/// `r_p > 1` and a working-fluid heat-capacity ratio `γ`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Brayton {
    pressure_ratio: f64,
    gamma: HeatCapacityRatio,
}

impl Brayton {
    /// Build a Brayton cycle from its pressure ratio and working-fluid
    /// `γ`.
    ///
    /// # Errors
    ///
    /// - [`CycleError::NotFinite`] if `pressure_ratio` is non-finite.
    /// - [`CycleError::RatioTooLow`] if `pressure_ratio <= 1` (no
    ///   compression means no net work).
    pub fn new(pressure_ratio: f64, gamma: HeatCapacityRatio) -> Result<Self> {
        let rp = crate::error::finite("pressure_ratio", pressure_ratio)?;
        if rp <= 1.0 {
            return Err(CycleError::RatioTooLow {
                name: "pressure_ratio",
                value: rp,
            });
        }
        Ok(Self {
            pressure_ratio: rp,
            gamma,
        })
    }

    /// Build a Brayton cycle for the air-standard working fluid
    /// (`γ = 1.4`).
    ///
    /// # Errors
    ///
    /// Same as [`Brayton::new`].
    pub fn with_air(pressure_ratio: f64) -> Result<Self> {
        Self::new(pressure_ratio, HeatCapacityRatio::air())
    }

    /// The compressor pressure ratio `r_p`.
    pub fn pressure_ratio(self) -> f64 {
        self.pressure_ratio
    }

    /// The working-fluid heat-capacity ratio `γ`.
    pub fn gamma(self) -> HeatCapacityRatio {
        self.gamma
    }

    /// The air-standard Brayton thermal efficiency,
    /// `η = 1 - 1 / r_p^((γ - 1) / γ)`.
    ///
    /// Always strictly inside `(0, 1)`: with `r_p > 1` and `γ > 1` the
    /// exponent `(γ - 1) / γ` is in `(0, 1)`, so `r_p^((γ-1)/γ) > 1` and
    /// the subtracted reciprocal lies in `(0, 1)`.
    pub fn efficiency(self) -> f64 {
        let g = self.gamma.value();
        let exponent = (g - 1.0) / g;
        1.0 - 1.0 / self.pressure_ratio.powf(exponent)
    }
}

/// The air-standard Brayton efficiency `1 - 1/r_p^((γ-1)/γ)` from a
/// pressure ratio and the working-fluid `γ`.
///
/// Convenience wrapper around [`Brayton::new`] + [`Brayton::efficiency`].
///
/// # Errors
///
/// Propagates the validation errors of [`Brayton::new`] and
/// [`HeatCapacityRatio::new`].
pub fn brayton_efficiency(pressure_ratio: f64, gamma: f64) -> Result<f64> {
    let gamma = HeatCapacityRatio::new(gamma)?;
    Ok(Brayton::new(pressure_ratio, gamma)?.efficiency())
}
