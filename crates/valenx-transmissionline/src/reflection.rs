//! Reflection at a resistive load and the standing-wave figures of
//! merit derived from it.
//!
//! ## Model
//!
//! When a line of real characteristic impedance `Z0` is terminated in a
//! purely resistive load `ZL`, the voltage **reflection coefficient** is
//! the real number
//!
//! ```text
//! gamma = (ZL - Z0) / (ZL + Z0)
//! ```
//!
//! bounded by `-1 <= gamma <= 1` for any passive resistive load:
//!
//! - `ZL = Z0` (matched): `gamma = 0` — no reflection.
//! - `ZL = 0` (short): `gamma = -1` — full reflection, inverted.
//! - `ZL -> infinity` (open): `gamma = +1` — full reflection, in phase.
//!
//! From `|gamma|` follow the standard derived quantities:
//!
//! ```text
//! VSWR        = (1 + |gamma|) / (1 - |gamma|)
//! return loss = -20 * log10(|gamma|)            (dB, >= 0)
//! mismatch    = -10 * log10(1 - |gamma|^2)      (dB, >= 0)
//! ```
//!
//! `VSWR >= 1` always; it diverges to `+infinity` as `|gamma| -> 1`.
//! Return loss is `+infinity` for a perfect match (`|gamma| = 0`) and
//! `0 dB` for total reflection (`|gamma| = 1`).

use serde::{Deserialize, Serialize};

use crate::error::{ensure_non_negative, ensure_positive, TlError};
use crate::line::Line;

/// A purely resistive termination presented to a transmission line.
///
/// The model here is the classic real-load case: `ZL` is a
/// non-negative resistance in ohms. `0 Ω` denotes a short circuit; an
/// idealised open circuit is represented separately by [`Load::Open`]
/// so that the infinite-impedance limit is exact rather than a large
/// finite stand-in.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Load {
    /// A finite resistive load of the given magnitude in ohms (`>= 0`).
    ///
    /// Use [`Load::resistive`] to construct one with validation.
    Resistive {
        /// Load resistance in ohms. Non-negative and finite.
        ohms: f64,
    },
    /// An idealised open circuit (`ZL -> infinity`), giving
    /// `gamma = +1` exactly.
    Open,
}

impl Load {
    /// Construct a finite resistive load of `ohms` ohms.
    ///
    /// `ohms` must be finite and non-negative. A value of `0.0` is a
    /// valid short circuit.
    ///
    /// # Errors
    ///
    /// Returns [`TlError::Negative`] if `ohms < 0`, or
    /// [`TlError::NotFinite`] if `ohms` is `NaN` / `±∞`. (For the
    /// infinite-impedance limit use [`Load::Open`] rather than passing
    /// `f64::INFINITY`.)
    pub fn resistive(ohms: f64) -> Result<Self, TlError> {
        let ohms = ensure_non_negative("load_ohms", ohms)?;
        Ok(Load::Resistive { ohms })
    }

    /// A short circuit (`ZL = 0`), giving `gamma = -1`.
    ///
    /// Convenience alias for `Load::resistive(0.0)`.
    pub fn short() -> Self {
        Load::Resistive { ohms: 0.0 }
    }
}

/// The result of reflecting an incident wave off a load on a line.
///
/// Holds the signed reflection coefficient together with the line and
/// load it was computed from, and exposes every derived figure of merit
/// as a method. Construct via [`Line::reflection`] or
/// [`Reflection::from_line_load`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Reflection {
    z0_ohms: f64,
    /// `None` for an open-circuit load (`ZL -> infinity`); otherwise the
    /// finite load resistance in ohms.
    load_ohms: Option<f64>,
    /// Signed voltage reflection coefficient in `[-1, 1]`.
    gamma: f64,
}

impl Reflection {
    /// Evaluate the reflection of `load` on `line`.
    ///
    /// For a finite resistive load this applies
    /// `gamma = (ZL - Z0) / (ZL + Z0)`; for [`Load::Open`] it returns
    /// the exact open-circuit limit `gamma = +1`. The result is always
    /// finite with `|gamma| <= 1`.
    pub fn from_line_load(line: Line, load: Load) -> Self {
        let z0 = line.z0_ohms();
        match load {
            Load::Open => Self {
                z0_ohms: z0,
                load_ohms: None,
                gamma: 1.0,
            },
            Load::Resistive { ohms } => {
                // z0 > 0 and ohms >= 0 ⇒ (ohms + z0) > 0, so the divide
                // is well defined and gamma ∈ [-1, 1].
                let gamma = (ohms - z0) / (ohms + z0);
                Self {
                    z0_ohms: z0,
                    load_ohms: Some(ohms),
                    gamma,
                }
            }
        }
    }

    /// Build a reflection report directly from a known signed
    /// reflection coefficient `gamma`, bypassing the line/load model.
    ///
    /// Useful when `gamma` is measured (e.g. read off a vector network
    /// analyser) rather than computed from impedances. The associated
    /// line and load impedances are then unknown and reported as `None`.
    ///
    /// `gamma` must be finite with `|gamma| <= 1`.
    ///
    /// # Errors
    ///
    /// Returns [`TlError::GammaOutOfRange`] if `|gamma| > 1`, or
    /// [`TlError::NotFinite`] if `gamma` is not finite.
    pub fn from_gamma(gamma: f64) -> Result<Self, TlError> {
        if !gamma.is_finite() {
            return Err(TlError::NotFinite {
                name: "gamma",
                value: gamma,
            });
        }
        if gamma.abs() > 1.0 {
            return Err(TlError::GammaOutOfRange { value: gamma });
        }
        Ok(Self {
            z0_ohms: f64::NAN,
            load_ohms: None,
            gamma,
        })
    }

    /// The signed voltage reflection coefficient `gamma` in `[-1, 1]`.
    ///
    /// Negative for `ZL < Z0` (phase-inverted reflection, e.g. a short),
    /// positive for `ZL > Z0` (in-phase reflection, e.g. an open), and
    /// `0` for a perfect match.
    #[inline]
    pub fn gamma(&self) -> f64 {
        self.gamma
    }

    /// The reflection magnitude `|gamma|` in `[0, 1]`.
    #[inline]
    pub fn gamma_magnitude(&self) -> f64 {
        self.gamma.abs()
    }

    /// The fraction of incident **power** that is reflected,
    /// `|gamma|^2` in `[0, 1]`.
    #[inline]
    pub fn power_reflected_fraction(&self) -> f64 {
        self.gamma * self.gamma
    }

    /// The fraction of incident **power** delivered to the load,
    /// `1 - |gamma|^2` in `[0, 1]`.
    #[inline]
    pub fn power_transmitted_fraction(&self) -> f64 {
        1.0 - self.power_reflected_fraction()
    }

    /// `true` if the load is a perfect match (`gamma == 0`), within an
    /// exact floating-point comparison.
    #[inline]
    pub fn is_matched(&self) -> bool {
        self.gamma == 0.0
    }

    /// `true` if the reflection is total (`|gamma| == 1`), i.e. an open
    /// or short circuit.
    #[inline]
    pub fn is_total_reflection(&self) -> bool {
        self.gamma_magnitude() == 1.0
    }

    /// The voltage standing-wave ratio,
    /// `VSWR = (1 + |gamma|) / (1 - |gamma|)`.
    ///
    /// Returns `None` for total reflection (`|gamma| = 1`), where the
    /// ratio diverges to `+infinity` (open / short circuit). Otherwise
    /// the value is finite and `>= 1`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_transmissionline::{Line, Load};
    ///
    /// let line = Line::from_z0(50.0).unwrap();
    /// // 100 Ω on a 50 Ω line ⇒ |gamma| = 1/3 ⇒ VSWR = 2.
    /// let r = line.reflection(Load::resistive(100.0).unwrap());
    /// assert!((r.vswr().unwrap() - 2.0).abs() < 1e-12);
    /// // Open circuit ⇒ VSWR diverges.
    /// assert!(line.reflection(Load::Open).vswr().is_none());
    /// ```
    pub fn vswr(&self) -> Option<f64> {
        let mag = self.gamma_magnitude();
        if mag >= 1.0 {
            None
        } else {
            Some((1.0 + mag) / (1.0 - mag))
        }
    }

    /// The return loss in decibels, `-20 * log10(|gamma|)`.
    ///
    /// Return loss is non-negative and *larger is better* (less
    /// reflected power). Returns `None` for a perfect match
    /// (`|gamma| = 0`), where it is `+infinity`. For total reflection
    /// (`|gamma| = 1`) it is exactly `0 dB`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_transmissionline::{Line, Load};
    ///
    /// let line = Line::from_z0(50.0).unwrap();
    /// // Matched ⇒ infinite return loss ⇒ None.
    /// assert!(line.reflection(Load::resistive(50.0).unwrap()).return_loss_db().is_none());
    /// // Short ⇒ |gamma| = 1 ⇒ 0 dB.
    /// let rl = line.reflection(Load::short()).return_loss_db().unwrap();
    /// assert!(rl.abs() < 1e-12);
    /// ```
    pub fn return_loss_db(&self) -> Option<f64> {
        let mag = self.gamma_magnitude();
        if mag == 0.0 {
            None
        } else {
            Some(-20.0 * mag.log10())
        }
    }

    /// The mismatch loss in decibels, `-10 * log10(1 - |gamma|^2)`.
    ///
    /// This is the power lost relative to a matched load purely due to
    /// reflection. It is `0 dB` for a perfect match and diverges to
    /// `+infinity` for total reflection (returned as `None`), where no
    /// power reaches the load.
    pub fn mismatch_loss_db(&self) -> Option<f64> {
        let transmitted = self.power_transmitted_fraction();
        if transmitted <= 0.0 {
            None
        } else {
            Some(-10.0 * transmitted.log10())
        }
    }

    /// The characteristic impedance `Z0` (ohms) this reflection was
    /// computed against, or `None` if the reflection was built directly
    /// from a measured `gamma` via [`Reflection::from_gamma`].
    pub fn z0_ohms(&self) -> Option<f64> {
        if self.z0_ohms.is_finite() {
            Some(self.z0_ohms)
        } else {
            None
        }
    }

    /// The load resistance `ZL` (ohms), or `None` for an open circuit or
    /// for a reflection built from a measured `gamma`.
    pub fn load_ohms(&self) -> Option<f64> {
        self.load_ohms
    }
}

/// Compute a load resistance from a known characteristic impedance and a
/// signed reflection coefficient — the algebraic inverse of
/// `gamma = (ZL - Z0) / (ZL + Z0)`.
///
/// Solving for `ZL` gives `ZL = Z0 * (1 + gamma) / (1 - gamma)`.
///
/// `z0_ohms` must be finite and strictly positive; `gamma` must be
/// finite with `-1 <= gamma < 1` (the open-circuit limit `gamma = 1`
/// would require infinite resistance and is rejected).
///
/// # Errors
///
/// Returns [`TlError::NonPositive`] / [`TlError::NotFinite`] for a bad
/// `z0_ohms`, or [`TlError::GammaOutOfRange`] if `gamma` is outside
/// `[-1, 1)`.
///
/// # Examples
///
/// ```
/// use valenx_transmissionline::load_from_gamma;
///
/// // gamma = 1/3 on a 50 Ω line ⇒ ZL = 100 Ω.
/// let zl = load_from_gamma(50.0, 1.0 / 3.0).unwrap();
/// assert!((zl - 100.0).abs() < 1e-9);
/// ```
pub fn load_from_gamma(z0_ohms: f64, gamma: f64) -> Result<f64, TlError> {
    let z0 = ensure_positive("z0_ohms", z0_ohms)?;
    if !gamma.is_finite() {
        return Err(TlError::NotFinite {
            name: "gamma",
            value: gamma,
        });
    }
    // gamma == 1 ⇒ open circuit ⇒ infinite ZL; reject it explicitly so
    // we never divide by zero or return a non-finite resistance.
    if !(-1.0..1.0).contains(&gamma) {
        return Err(TlError::GammaOutOfRange { value: gamma });
    }
    Ok(z0 * (1.0 + gamma) / (1.0 - gamma))
}
