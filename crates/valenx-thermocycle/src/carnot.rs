//! The Carnot cycle — the upper bound on heat-engine efficiency.
//!
//! Carnot's theorem states that no heat engine operating between two
//! thermal reservoirs at absolute temperatures `T_c < T_h` can be more
//! efficient than a reversible (Carnot) engine running between the same
//! two reservoirs, whose efficiency is
//!
//! ```text
//! η_carnot = 1 - T_c / T_h
//! ```
//!
//! with both temperatures in an absolute scale (kelvin). This is the
//! benchmark every other cycle in this crate is compared against.

use crate::error::{CycleError, Result};
use serde::{Deserialize, Serialize};

/// A pair of reservoir temperatures defining a Carnot engine, validated so
/// that `0 < T_c < T_h` with both finite.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Carnot {
    t_cold: f64,
    t_hot: f64,
}

impl Carnot {
    /// Build a Carnot engine from its cold and hot reservoir absolute
    /// temperatures, both in **kelvin**.
    ///
    /// # Errors
    ///
    /// - [`CycleError::NotFinite`] if either temperature is non-finite.
    /// - [`CycleError::NotPositive`] if either temperature is `<= 0`
    ///   (the kelvin scale starts at absolute zero).
    /// - [`CycleError::TemperatureOrder`] if `t_cold >= t_hot`, leaving
    ///   the engine no temperature drop to exploit.
    pub fn new(t_cold: f64, t_hot: f64) -> Result<Self> {
        let t_cold = crate::error::positive("t_cold", t_cold)?;
        let t_hot = crate::error::positive("t_hot", t_hot)?;
        if t_cold >= t_hot {
            return Err(CycleError::TemperatureOrder { t_cold, t_hot });
        }
        Ok(Self { t_cold, t_hot })
    }

    /// The cold-reservoir absolute temperature in kelvin.
    pub fn t_cold(self) -> f64 {
        self.t_cold
    }

    /// The hot-reservoir absolute temperature in kelvin.
    pub fn t_hot(self) -> f64 {
        self.t_hot
    }

    /// The Carnot (maximum reversible) thermal efficiency,
    /// `η = 1 - T_c / T_h`.
    ///
    /// Because the constructor guarantees `0 < T_c < T_h`, the result is
    /// always strictly inside the open interval `(0, 1)`.
    pub fn efficiency(self) -> f64 {
        1.0 - self.t_cold / self.t_hot
    }

    /// The coefficient of performance of the *reversed* Carnot cycle run
    /// as a **refrigerator** (heat pumped out of the cold reservoir per
    /// unit work input), `COP_ref = T_c / (T_h - T_c)`.
    ///
    /// Unlike the efficiency, this is unbounded above: as `T_h → T_c⁺` an
    /// ideal refrigerator approaches infinite COP.
    pub fn cop_refrigerator(self) -> f64 {
        self.t_cold / (self.t_hot - self.t_cold)
    }

    /// The coefficient of performance of the reversed Carnot cycle run as
    /// a **heat pump** (heat delivered to the hot reservoir per unit work
    /// input), `COP_hp = T_h / (T_h - T_c)`.
    ///
    /// Exactly one greater than the refrigerator COP, reflecting energy
    /// conservation `Q_h = Q_c + W`.
    pub fn cop_heat_pump(self) -> f64 {
        self.t_hot / (self.t_hot - self.t_cold)
    }
}

/// The Carnot efficiency `1 - T_c / T_h` from two reservoir temperatures
/// in kelvin, validating them in passing.
///
/// Convenience wrapper around [`Carnot::new`] + [`Carnot::efficiency`].
///
/// # Errors
///
/// Propagates the validation errors of [`Carnot::new`].
pub fn carnot_efficiency(t_cold: f64, t_hot: f64) -> Result<f64> {
    Ok(Carnot::new(t_cold, t_hot)?.efficiency())
}
