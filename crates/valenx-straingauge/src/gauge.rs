//! The strain gauge itself: gauge factor, fractional resistance change,
//! and the uniaxial Hooke's-law stress that the measured strain implies.
//!
//! ## Model
//!
//! A bonded metallic-foil strain gauge changes its resistance in
//! proportion to the mechanical strain `ε` it experiences. The constant
//! of proportionality is the dimensionless **gauge factor**
//!
//! ```text
//! GF = (ΔR / R) / ε
//! ```
//!
//! so the fractional resistance change is `ΔR/R = GF · ε`. Typical
//! constantan foil gauges have `GF ≈ 2.0`.
//!
//! Strain is reported here as a pure ratio (m/m). One *microstrain*
//! (`µε`) is `1e-6`; the helper [`microstrain`] converts.
//!
//! Given the material's Young's modulus `E`, the uniaxial stress that
//! produced the strain follows from Hooke's law `σ = E · ε`. The stress
//! is returned in the same pressure unit as `E` (pass `E` in pascals to
//! get stress in pascals, or in MPa to get MPa — the relation is
//! unit-agnostic as long as you are consistent).

use serde::{Deserialize, Serialize};

use crate::error::{Result, StrainGaugeError};

/// A bonded strain gauge characterised by its gauge factor.
///
/// The gauge factor is the only intrinsic property the resistance model
/// needs; nominal resistance (commonly 120 Ω or 350 Ω) only matters when
/// you want an *absolute* resistance change rather than the fractional
/// one, so it is supplied per-call to [`Gauge::resistance_change_ohm`]
/// rather than stored.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Gauge {
    /// Dimensionless gauge factor `GF = (ΔR/R)/ε`. Must be `> 0`.
    pub gauge_factor: f64,
}

impl Gauge {
    /// Build a gauge from its gauge factor, validating that the factor
    /// is finite and strictly positive.
    ///
    /// ```
    /// use valenx_straingauge::Gauge;
    ///
    /// let g = Gauge::new(2.0).unwrap();
    /// assert!((g.gauge_factor - 2.0).abs() < 1e-12);
    /// assert!(Gauge::new(0.0).is_err());
    /// ```
    pub fn new(gauge_factor: f64) -> Result<Self> {
        let gauge_factor = StrainGaugeError::positive("gauge_factor", gauge_factor)?;
        Ok(Self { gauge_factor })
    }

    /// A standard constantan foil gauge (`GF = 2.0`).
    pub fn constantan() -> Self {
        // Safe: 2.0 is finite and positive.
        Self { gauge_factor: 2.0 }
    }

    /// Fractional resistance change `ΔR/R = GF · ε` for the supplied
    /// strain.
    ///
    /// Strain is signed: positive in tension, negative in compression,
    /// and the returned ratio carries the same sign. The strain must be
    /// finite (zero is allowed and yields `0.0`).
    ///
    /// ```
    /// use valenx_straingauge::Gauge;
    ///
    /// let g = Gauge::constantan();
    /// // GF = 2, ε = 1000 µε = 1e-3  ⇒  ΔR/R = 2e-3.
    /// let dr = g.fractional_resistance_change(1.0e-3).unwrap();
    /// assert!((dr - 2.0e-3).abs() < 1e-15);
    /// ```
    pub fn fractional_resistance_change(&self, strain: f64) -> Result<f64> {
        let strain = StrainGaugeError::finite("strain", strain)?;
        Ok(self.gauge_factor * strain)
    }

    /// Absolute resistance change `ΔR = R · GF · ε` (in the same unit as
    /// `nominal_resistance_ohm`).
    ///
    /// The nominal resistance must be finite and strictly positive
    /// (e.g. `120.0` or `350.0` ohms).
    ///
    /// ```
    /// use valenx_straingauge::Gauge;
    ///
    /// let g = Gauge::constantan();
    /// // 350 Ω gauge, ε = 1e-3  ⇒  ΔR = 350 · 2 · 1e-3 = 0.7 Ω.
    /// let dr = g.resistance_change_ohm(350.0, 1.0e-3).unwrap();
    /// assert!((dr - 0.7).abs() < 1e-12);
    /// ```
    pub fn resistance_change_ohm(&self, nominal_resistance_ohm: f64, strain: f64) -> Result<f64> {
        let r = StrainGaugeError::positive("nominal_resistance_ohm", nominal_resistance_ohm)?;
        let frac = self.fractional_resistance_change(strain)?;
        Ok(r * frac)
    }

    /// Recover the strain implied by a measured fractional resistance
    /// change, `ε = (ΔR/R) / GF` — the inverse of
    /// [`Gauge::fractional_resistance_change`].
    ///
    /// ```
    /// use valenx_straingauge::Gauge;
    ///
    /// let g = Gauge::constantan();
    /// let eps = g.strain_from_fractional_resistance_change(2.0e-3).unwrap();
    /// assert!((eps - 1.0e-3).abs() < 1e-15);
    /// ```
    pub fn strain_from_fractional_resistance_change(
        &self,
        fractional_resistance_change: f64,
    ) -> Result<f64> {
        let frac =
            StrainGaugeError::finite("fractional_resistance_change", fractional_resistance_change)?;
        // gauge_factor is guaranteed positive (and hence non-zero) by the
        // constructor, so this division never produces NaN/Inf.
        Ok(frac / self.gauge_factor)
    }
}

/// Convert a strain expressed in *microstrain* (`µε`, units of `1e-6`)
/// to the pure ratio (m/m) used throughout this crate.
///
/// ```
/// use valenx_straingauge::microstrain;
///
/// assert!((microstrain(1000.0) - 1.0e-3).abs() < 1e-15);
/// ```
pub fn microstrain(micro: f64) -> f64 {
    micro * 1.0e-6
}

/// Uniaxial Hooke's-law stress `σ = E · ε`.
///
/// `youngs_modulus` is Young's modulus `E` and must be finite and
/// strictly positive. The strain `ε` is signed and must be finite. The
/// returned stress is in the same pressure unit as `E` (pass `E` in
/// pascals for stress in pascals).
///
/// ```
/// use valenx_straingauge::stress;
///
/// // Steel: E = 200 GPa, ε = 1e-3  ⇒  σ = 200e6 Pa = 200 MPa.
/// let sigma = stress(200.0e9, 1.0e-3).unwrap();
/// assert!((sigma - 200.0e6).abs() < 1.0e-3);
/// ```
pub fn stress(youngs_modulus: f64, strain: f64) -> Result<f64> {
    let e = StrainGaugeError::positive("youngs_modulus", youngs_modulus)?;
    let eps = StrainGaugeError::finite("strain", strain)?;
    Ok(e * eps)
}

/// Recover the strain from a measured uniaxial stress, `ε = σ / E` — the
/// inverse of [`stress`].
///
/// `youngs_modulus` must be finite and strictly positive; `stress` must
/// be finite.
///
/// ```
/// use valenx_straingauge::strain_from_stress;
///
/// let eps = strain_from_stress(200.0e9, 200.0e6).unwrap();
/// assert!((eps - 1.0e-3).abs() < 1e-15);
/// ```
pub fn strain_from_stress(youngs_modulus: f64, stress: f64) -> Result<f64> {
    let e = StrainGaugeError::positive("youngs_modulus", youngs_modulus)?;
    let sigma = StrainGaugeError::finite("stress", stress)?;
    // e is guaranteed positive (non-zero), so this is finite.
    Ok(sigma / e)
}
