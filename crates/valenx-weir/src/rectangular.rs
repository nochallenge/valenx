//! Sharp-crested **rectangular** weir discharge.
//!
//! A rectangular weir is a horizontal crest of width (crest length) `L`
//! over which water spills with an upstream head `H` measured from the
//! crest to the undisturbed free surface. Integrating the ideal
//! velocity `v(z) = √(2 g z)` over the rectangular opening and applying
//! a lumped discharge coefficient `Cd` gives the standard weir
//! equation
//!
//! ```text
//!   Q = Cd · (2/3) · √(2 g) · L · H^(3/2)
//! ```
//!
//! so discharge scales with the **3/2 power of head** and **linearly**
//! with both the crest length and the discharge coefficient.

use crate::error::{require_positive, WeirError};
use crate::G_STANDARD;
use serde::{Deserialize, Serialize};

/// The dimensionless `2/3` prefactor in the rectangular weir equation.
///
/// It arises from integrating `√(2 g z)` over `z ∈ [0, H]`:
/// `∫₀ᴴ √z dz = (2/3) H^(3/2)`.
pub const RECT_COEFFICIENT: f64 = 2.0 / 3.0;

/// A sharp-crested rectangular weir, validated on construction.
///
/// The struct stores only the three independent quantities that define
/// the discharge: the crest length `L`, the discharge coefficient `Cd`,
/// and the gravitational acceleration `g`. The head `H` is supplied per
/// evaluation to [`discharge`](RectangularWeir::discharge) because a
/// single weir is rated across a range of heads.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RectangularWeir {
    /// Crest length (weir width) `L`, in metres. Strictly positive.
    crest_length_m: f64,
    /// Dimensionless discharge coefficient `Cd`. Strictly positive;
    /// physically `≈ 0.6`–`0.65` for a ventilated sharp-crested weir.
    discharge_coefficient: f64,
    /// Gravitational acceleration `g`, in m·s⁻². Strictly positive.
    gravity: f64,
}

impl RectangularWeir {
    /// Construct a rectangular weir at standard gravity
    /// (`g = `[`G_STANDARD`]).
    ///
    /// # Errors
    ///
    /// Returns [`WeirError::NonPositive`] / [`WeirError::NotFinite`] if
    /// `crest_length_m` or `discharge_coefficient` is not a finite,
    /// strictly-positive number.
    pub fn new(crest_length_m: f64, discharge_coefficient: f64) -> Result<Self, WeirError> {
        Self::with_gravity(crest_length_m, discharge_coefficient, G_STANDARD)
    }

    /// Construct a rectangular weir with an explicit gravitational
    /// acceleration `gravity` (m·s⁻²).
    ///
    /// # Errors
    ///
    /// Returns a [`WeirError`] if any of `crest_length_m`,
    /// `discharge_coefficient` or `gravity` is not a finite,
    /// strictly-positive number.
    pub fn with_gravity(
        crest_length_m: f64,
        discharge_coefficient: f64,
        gravity: f64,
    ) -> Result<Self, WeirError> {
        Ok(Self {
            crest_length_m: require_positive("crest_length", crest_length_m)?,
            discharge_coefficient: require_positive(
                "discharge_coefficient",
                discharge_coefficient,
            )?,
            gravity: require_positive("gravity", gravity)?,
        })
    }

    /// Crest length (weir width) `L`, in metres.
    pub fn crest_length_m(&self) -> f64 {
        self.crest_length_m
    }

    /// Dimensionless discharge coefficient `Cd`.
    pub fn discharge_coefficient(&self) -> f64 {
        self.discharge_coefficient
    }

    /// Gravitational acceleration `g`, in m·s⁻².
    pub fn gravity(&self) -> f64 {
        self.gravity
    }

    /// Volumetric discharge `Q` (m³·s⁻¹) at upstream head
    /// `head_m` (metres):
    ///
    /// ```text
    ///   Q = Cd · (2/3) · √(2 g) · L · H^(3/2)
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`WeirError`] if `head_m` is not a finite,
    /// strictly-positive number. A zero or negative head is rejected
    /// rather than returning `Q = 0`, because it indicates the weir is
    /// not actually flowing and the caller almost certainly has a bug.
    pub fn discharge(&self, head_m: f64) -> Result<f64, WeirError> {
        let h = require_positive("head", head_m)?;
        Ok(self.discharge_coefficient
            * RECT_COEFFICIENT
            * (2.0 * self.gravity).sqrt()
            * self.crest_length_m
            * h.powf(1.5))
    }
}
