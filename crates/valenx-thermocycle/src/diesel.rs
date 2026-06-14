//! The air-standard Diesel cycle — the ideal compression-ignition engine.
//!
//! The Diesel cycle differs from the [`Otto`](crate::otto) cycle in that
//! heat is added at **constant pressure** (fuel injected into hot,
//! compressed air and burned as the piston begins its power stroke)
//! rather than at constant volume. Under the air-standard assumptions its
//! efficiency depends on the compression ratio `r`, the **cutoff ratio**
//! `r_c` (the volume ratio over which combustion occurs, `V₃ / V₂ ≥ 1`)
//! and the heat-capacity ratio `γ`:
//!
//! ```text
//!                1     r_c^γ - 1
//! η_diesel = 1 - ─── · ─────────────
//!                r^(γ-1)   γ (r_c - 1)
//! ```
//!
//! The bracketed factor `(r_c^γ - 1) / (γ (r_c - 1))` is always greater
//! than one for `r_c > 1`, so for a *given* compression ratio the Diesel
//! cycle is slightly **less** efficient than the Otto cycle — the
//! constant-pressure heat addition is the source of that penalty. As the
//! cutoff ratio approaches unity (`r_c → 1`) the factor tends to one and
//! the Diesel efficiency recovers the Otto value, since the constant-
//! pressure heat-addition shrinks to a constant-volume one.

use crate::error::{CycleError, Result};
use crate::gas::HeatCapacityRatio;
use serde::{Deserialize, Serialize};

/// A validated air-standard Diesel cycle: compression ratio `r > 1`,
/// cutoff ratio `r_c >= 1`, and working-fluid heat-capacity ratio `γ`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Diesel {
    compression_ratio: f64,
    cutoff_ratio: f64,
    gamma: HeatCapacityRatio,
}

impl Diesel {
    /// Build a Diesel cycle from its compression ratio, cutoff ratio and
    /// working-fluid `γ`.
    ///
    /// # Errors
    ///
    /// - [`CycleError::NotFinite`] if either ratio is non-finite.
    /// - [`CycleError::RatioTooLow`] if `compression_ratio <= 1`.
    /// - [`CycleError::NotPositive`] if `cutoff_ratio < 1` — the cutoff
    ///   ratio is `V₃ / V₂` with `V₃ ≥ V₂` (heat is added as the gas
    ///   expands), so it cannot be below one. A value of exactly `1`
    ///   degenerates to the Otto cycle and is accepted.
    pub fn new(
        compression_ratio: f64,
        cutoff_ratio: f64,
        gamma: HeatCapacityRatio,
    ) -> Result<Self> {
        let r = crate::error::finite("compression_ratio", compression_ratio)?;
        if r <= 1.0 {
            return Err(CycleError::RatioTooLow {
                name: "compression_ratio",
                value: r,
            });
        }
        let rc = crate::error::finite("cutoff_ratio", cutoff_ratio)?;
        if rc < 1.0 {
            return Err(CycleError::NotPositive {
                name: "cutoff_ratio",
                floor: 1.0,
                value: rc,
            });
        }
        Ok(Self {
            compression_ratio: r,
            cutoff_ratio: rc,
            gamma,
        })
    }

    /// Build a Diesel cycle for the air-standard working fluid
    /// (`γ = 1.4`).
    ///
    /// # Errors
    ///
    /// Same as [`Diesel::new`].
    pub fn with_air(compression_ratio: f64, cutoff_ratio: f64) -> Result<Self> {
        Self::new(compression_ratio, cutoff_ratio, HeatCapacityRatio::air())
    }

    /// The compression ratio `r`.
    pub fn compression_ratio(self) -> f64 {
        self.compression_ratio
    }

    /// The cutoff ratio `r_c = V₃ / V₂`.
    pub fn cutoff_ratio(self) -> f64 {
        self.cutoff_ratio
    }

    /// The working-fluid heat-capacity ratio `γ`.
    pub fn gamma(self) -> HeatCapacityRatio {
        self.gamma
    }

    /// The Diesel-relative-to-Otto **cutoff penalty factor**
    /// `(r_c^γ - 1) / (γ (r_c - 1))`.
    ///
    /// This is `>= 1` for all `r_c >= 1`, equalling one exactly in the
    /// `r_c = 1` limit (handled analytically to avoid the `0 / 0` form),
    /// and it is what makes the Diesel cycle less efficient than an Otto
    /// cycle of the same compression ratio.
    pub fn cutoff_penalty(self) -> f64 {
        let g = self.gamma.value();
        let rc = self.cutoff_ratio;
        // At r_c = 1 the expression is the indeterminate 0/0; its limit is
        // exactly 1 (L'Hopital), which is also the Otto cycle's factor.
        if (rc - 1.0).abs() <= f64::EPSILON {
            1.0
        } else {
            (rc.powf(g) - 1.0) / (g * (rc - 1.0))
        }
    }

    /// The air-standard Diesel thermal efficiency,
    /// `η = 1 - (1 / r^(γ-1)) · (r_c^γ - 1) / (γ (r_c - 1))`.
    ///
    /// Always strictly inside `(0, 1)`: the Otto base `1 - 1/r^(γ-1)` is
    /// in `(0, 1)`, and multiplying the subtracted reciprocal term by the
    /// cutoff penalty (which is `>= 1`) keeps the whole expression below
    /// one and above zero for physical inputs.
    pub fn efficiency(self) -> f64 {
        let base = 1.0 / self.compression_ratio.powf(self.gamma.value() - 1.0);
        1.0 - base * self.cutoff_penalty()
    }
}

/// The air-standard Diesel efficiency from compression ratio, cutoff ratio
/// and working-fluid `γ`.
///
/// Convenience wrapper around [`Diesel::new`] + [`Diesel::efficiency`].
///
/// # Errors
///
/// Propagates the validation errors of [`Diesel::new`] and
/// [`HeatCapacityRatio::new`].
pub fn diesel_efficiency(compression_ratio: f64, cutoff_ratio: f64, gamma: f64) -> Result<f64> {
    let gamma = HeatCapacityRatio::new(gamma)?;
    Ok(Diesel::new(compression_ratio, cutoff_ratio, gamma)?.efficiency())
}
