//! The Wheatstone bridge that turns a gauge's resistance change into a
//! measurable voltage.
//!
//! ## Model
//!
//! A strain gauge is wired into one (or more) legs of a Wheatstone
//! bridge excited by `Vin`. To first order in the strain `ε`, the
//! normalised bridge output is
//!
//! ```text
//! Vout / Vin = (N / 4) · GF · ε
//! ```
//!
//! where `N` is the number of *active* gauges arranged so their outputs
//! add:
//!
//! 1. **Quarter bridge** — one active gauge, three fixed resistors.
//!    `Vout/Vin = GF · ε / 4`.
//! 2. **Half bridge** — two active gauges (`N = 2`); twice the quarter
//!    output.
//! 3. **Full bridge** — four active gauges (`N = 4`); four times the
//!    quarter output.
//!
//! At zero strain every leg is balanced and the output is exactly zero
//! regardless of `GF`, `Vin` or configuration. This linearised relation
//! is the standard textbook bridge equation; see the crate-level
//! "Honest scope" note for the second-order and parasitic effects it
//! deliberately omits.

use serde::{Deserialize, Serialize};

use crate::error::{Result, StrainGaugeError};
use crate::gauge::Gauge;

/// Wheatstone-bridge wiring configuration, i.e. how many gauge legs are
/// active.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum BridgeConfig {
    /// One active gauge (`N = 1`): `Vout/Vin = GF · ε / 4`.
    Quarter,
    /// Two active gauges (`N = 2`): twice the quarter-bridge output.
    Half,
    /// Four active gauges (`N = 4`): four times the quarter-bridge
    /// output.
    Full,
}

impl BridgeConfig {
    /// Number of active gauges `N` for this configuration.
    ///
    /// ```
    /// use valenx_straingauge::BridgeConfig;
    ///
    /// assert_eq!(BridgeConfig::Quarter.active_arms(), 1);
    /// assert_eq!(BridgeConfig::Half.active_arms(), 2);
    /// assert_eq!(BridgeConfig::Full.active_arms(), 4);
    /// ```
    pub fn active_arms(self) -> u32 {
        match self {
            BridgeConfig::Quarter => 1,
            BridgeConfig::Half => 2,
            BridgeConfig::Full => 4,
        }
    }

    /// Output gain factor `N / 4` relative to `GF · ε`.
    ///
    /// The full bridge equation is `Vout/Vin = gain · GF · ε`, so the
    /// quarter bridge has gain `0.25`, the half bridge `0.5`, and the
    /// full bridge `1.0`.
    ///
    /// ```
    /// use valenx_straingauge::BridgeConfig;
    ///
    /// assert!((BridgeConfig::Quarter.gain() - 0.25).abs() < 1e-15);
    /// assert!((BridgeConfig::Full.gain() - 1.0).abs() < 1e-15);
    /// ```
    pub fn gain(self) -> f64 {
        self.active_arms() as f64 / 4.0
    }

    /// Human-readable label.
    pub fn label(self) -> &'static str {
        match self {
            BridgeConfig::Quarter => "Quarter bridge",
            BridgeConfig::Half => "Half bridge",
            BridgeConfig::Full => "Full bridge",
        }
    }
}

/// A strain-gauge Wheatstone bridge: a [`Gauge`] plus a wiring
/// [`BridgeConfig`].
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Bridge {
    /// The active strain gauge.
    pub gauge: Gauge,
    /// How many of its legs are active.
    pub config: BridgeConfig,
}

impl Bridge {
    /// Build a bridge from a gauge and a configuration.
    pub fn new(gauge: Gauge, config: BridgeConfig) -> Self {
        Self { gauge, config }
    }

    /// Normalised bridge output `Vout/Vin = (N/4) · GF · ε` for the
    /// supplied strain.
    ///
    /// Strain is signed and must be finite; the sign carries through to
    /// the output ratio. At `ε = 0` the result is exactly `0.0` (the
    /// balanced bridge) for every configuration.
    ///
    /// ```
    /// use valenx_straingauge::{Bridge, BridgeConfig, Gauge};
    ///
    /// let q = Bridge::new(Gauge::constantan(), BridgeConfig::Quarter);
    /// // GF = 2, ε = 1e-3  ⇒  Vout/Vin = 2 · 1e-3 / 4 = 5e-4.
    /// let ratio = q.output_ratio(1.0e-3).unwrap();
    /// assert!((ratio - 5.0e-4).abs() < 1e-15);
    /// ```
    pub fn output_ratio(&self, strain: f64) -> Result<f64> {
        let strain = StrainGaugeError::finite("strain", strain)?;
        Ok(self.config.gain() * self.gauge.gauge_factor * strain)
    }

    /// Absolute bridge output voltage `Vout = Vin · (N/4) · GF · ε` in
    /// the same unit as `excitation_voltage`.
    ///
    /// The excitation voltage `Vin` must be finite and strictly
    /// positive; strain must be finite.
    ///
    /// ```
    /// use valenx_straingauge::{Bridge, BridgeConfig, Gauge};
    ///
    /// let full = Bridge::new(Gauge::constantan(), BridgeConfig::Full);
    /// // Vin = 5 V, GF = 2, ε = 1e-3, full bridge ⇒
    /// // Vout = 5 · 1.0 · 2 · 1e-3 = 0.01 V = 10 mV.
    /// let v = full.output_voltage(5.0, 1.0e-3).unwrap();
    /// assert!((v - 0.01).abs() < 1e-12);
    /// ```
    pub fn output_voltage(&self, excitation_voltage: f64, strain: f64) -> Result<f64> {
        let vin = StrainGaugeError::positive("excitation_voltage", excitation_voltage)?;
        let ratio = self.output_ratio(strain)?;
        Ok(vin * ratio)
    }

    /// Recover the strain from a measured normalised bridge output,
    /// `ε = (Vout/Vin) / ((N/4) · GF)` — the inverse of
    /// [`Bridge::output_ratio`].
    ///
    /// ```
    /// use valenx_straingauge::{Bridge, BridgeConfig, Gauge};
    ///
    /// let half = Bridge::new(Gauge::constantan(), BridgeConfig::Half);
    /// let eps = half.strain_from_output_ratio(1.0e-3).unwrap();
    /// // Vout/Vin = 0.5 · 2 · ε  ⇒  ε = ratio / 1.0 = 1e-3.
    /// assert!((eps - 1.0e-3).abs() < 1e-15);
    /// ```
    pub fn strain_from_output_ratio(&self, output_ratio: f64) -> Result<f64> {
        let ratio = StrainGaugeError::finite("output_ratio", output_ratio)?;
        // gain() is 0.25/0.5/1.0 (non-zero) and gauge_factor is
        // guaranteed positive, so the denominator is non-zero and the
        // division is finite.
        let denom = self.config.gain() * self.gauge.gauge_factor;
        Ok(ratio / denom)
    }
}
