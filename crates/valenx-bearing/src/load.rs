//! Dynamic equivalent load `P = X·Fr + Y·Fa`.
//!
//! A bearing in service usually carries a combined radial *and* axial
//! load, but the life formula takes a single number — the **dynamic
//! equivalent load** `P`: the constant, purely radial (for radial
//! bearings) load that would give the same life as the real combined
//! load. ISO 281 / SKF write it as a linear combination
//!
//! ```text
//! P = X · Fr + Y · Fa
//! ```
//!
//! where `Fr` is the radial component, `Fa` the axial (thrust)
//! component, and `X`, `Y` are the dimensionless radial and axial load
//! factors read from the bearing manufacturer's table (they depend on
//! the bearing series and on the `Fa/Fr` ratio relative to the limit
//! `e`). This crate does **not** hide the `X`/`Y` selection: you pass
//! the factors in, and [`EquivalentLoad`] applies the formula.

use serde::{Deserialize, Serialize};

use crate::error::{require_non_negative, BearingError};

/// The dynamic equivalent radial load and the inputs it was built
/// from.
///
/// Construct it with [`EquivalentLoad::new`] (validated) and read the
/// combined value back with [`EquivalentLoad::value`]. The struct
/// keeps the constituent radial / axial loads and the `X` / `Y`
/// factors so a report can show the full derivation, not just the
/// answer.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EquivalentLoad {
    /// Radial load component `Fr` (newtons), `>= 0`.
    pub radial: f64,
    /// Axial (thrust) load component `Fa` (newtons), `>= 0`.
    pub axial: f64,
    /// Dimensionless radial load factor `X`, `>= 0`.
    pub x_factor: f64,
    /// Dimensionless axial load factor `Y`, `>= 0`.
    pub y_factor: f64,
}

impl EquivalentLoad {
    /// Build a validated dynamic-equivalent-load case from its radial
    /// and axial components and their load factors.
    ///
    /// All four arguments must be finite and non-negative: a load or a
    /// factor of zero is physically meaningful (a purely radial load
    /// has `Fa = 0`; a radial-only treatment uses `Y = 0`), but a
    /// negative value is not.
    ///
    /// # Errors
    ///
    /// Returns [`BearingError`] when any argument is `NaN`, infinite,
    /// or negative.
    ///
    /// ```
    /// use valenx_bearing::EquivalentLoad;
    /// // Fr = 8000 N, Fa = 3000 N, X = 0.56, Y = 1.6
    /// let p = EquivalentLoad::new(8000.0, 3000.0, 0.56, 1.6).unwrap();
    /// // P = 0.56*8000 + 1.6*3000 = 4480 + 4800 = 9280 N
    /// assert!((p.value() - 9280.0).abs() < 1e-9);
    /// ```
    pub fn new(
        radial: f64,
        axial: f64,
        x_factor: f64,
        y_factor: f64,
    ) -> Result<Self, BearingError> {
        let radial = require_non_negative("radial", radial)?;
        let axial = require_non_negative("axial", axial)?;
        let x_factor = require_non_negative("x_factor", x_factor)?;
        let y_factor = require_non_negative("y_factor", y_factor)?;
        Ok(Self {
            radial,
            axial,
            x_factor,
            y_factor,
        })
    }

    /// A purely radial load: the equivalent load is the radial force
    /// itself (`X = 1`, `Y = 0`, `Fa = 0`).
    ///
    /// This is the common shortcut for a deep-groove ball bearing
    /// carrying no thrust, where `P = Fr`.
    ///
    /// # Errors
    ///
    /// Returns [`BearingError`] when `radial` is `NaN`, infinite, or
    /// negative.
    ///
    /// ```
    /// use valenx_bearing::EquivalentLoad;
    /// let p = EquivalentLoad::radial_only(5000.0).unwrap();
    /// assert!((p.value() - 5000.0).abs() < 1e-9);
    /// ```
    pub fn radial_only(radial: f64) -> Result<Self, BearingError> {
        Self::new(radial, 0.0, 1.0, 0.0)
    }

    /// The dynamic equivalent load `P = X·Fr + Y·Fa` in newtons.
    ///
    /// Because every constituent was validated non-negative at
    /// construction, the result is always finite and non-negative.
    #[must_use]
    pub fn value(&self) -> f64 {
        self.x_factor * self.radial + self.y_factor * self.axial
    }
}
