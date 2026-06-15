//! The air-standard Otto cycle — the ideal spark-ignition (petrol)
//! engine.
//!
//! The Otto cycle idealises a four-stroke spark-ignition engine as two
//! adiabatic (isentropic) strokes bracketing two constant-volume
//! heat-transfer processes. Under the air-standard assumptions
//! (working fluid is an ideal gas with constant specific heats), its
//! thermal efficiency depends only on the **compression ratio** `r` (the
//! ratio of cylinder volume at bottom-dead-centre to that at
//! top-dead-centre) and the heat-capacity ratio `γ`:
//!
//! ```text
//! η_otto = 1 - 1 / r^(γ - 1)
//! ```
//!
//! Efficiency rises monotonically with `r` (more compression extracts
//! more work), which is exactly why higher-compression engines are more
//! efficient — bounded above by the Carnot limit for the same
//! temperature extremes.

use crate::error::{CycleError, Result};
use crate::gas::HeatCapacityRatio;
use serde::{Deserialize, Serialize};

/// A validated air-standard Otto cycle: a compression ratio `r > 1` and a
/// working-fluid heat-capacity ratio `γ`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Otto {
    compression_ratio: f64,
    gamma: HeatCapacityRatio,
}

impl Otto {
    /// Build an Otto cycle from its compression ratio and working-fluid
    /// `γ`.
    ///
    /// # Errors
    ///
    /// - [`CycleError::NotFinite`] if `compression_ratio` is non-finite.
    /// - [`CycleError::RatioTooLow`] if `compression_ratio <= 1` (an
    ///   engine that does not compress its charge does no net work).
    pub fn new(compression_ratio: f64, gamma: HeatCapacityRatio) -> Result<Self> {
        let r = crate::error::finite("compression_ratio", compression_ratio)?;
        if r <= 1.0 {
            return Err(CycleError::RatioTooLow {
                name: "compression_ratio",
                value: r,
            });
        }
        Ok(Self {
            compression_ratio: r,
            gamma,
        })
    }

    /// Build an Otto cycle for the air-standard working fluid
    /// (`γ = 1.4`).
    ///
    /// # Errors
    ///
    /// Same as [`Otto::new`].
    pub fn with_air(compression_ratio: f64) -> Result<Self> {
        Self::new(compression_ratio, HeatCapacityRatio::air())
    }

    /// The compression ratio `r`.
    pub fn compression_ratio(self) -> f64 {
        self.compression_ratio
    }

    /// The working-fluid heat-capacity ratio `γ`.
    pub fn gamma(self) -> HeatCapacityRatio {
        self.gamma
    }

    /// The air-standard Otto thermal efficiency,
    /// `η = 1 - 1 / r^(γ - 1)`.
    ///
    /// Always strictly inside `(0, 1)`: with `r > 1` and `γ > 1` the term
    /// `r^(γ - 1)` is strictly greater than one, so the subtracted
    /// reciprocal lies in `(0, 1)`.
    pub fn efficiency(self) -> f64 {
        1.0 - 1.0 / self.compression_ratio.powf(self.gamma.value() - 1.0)
    }
}

/// The air-standard Otto efficiency `1 - 1/r^(γ-1)` from a compression
/// ratio and the working-fluid `γ`.
///
/// Convenience wrapper around [`Otto::new`] + [`Otto::efficiency`].
///
/// # Errors
///
/// Propagates the validation errors of [`Otto::new`] and
/// [`HeatCapacityRatio::new`].
pub fn otto_efficiency(compression_ratio: f64, gamma: f64) -> Result<f64> {
    let gamma = HeatCapacityRatio::new(gamma)?;
    Ok(Otto::new(compression_ratio, gamma)?.efficiency())
}
