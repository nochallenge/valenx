//! Footing geometry: width and founding depth.
//!
//! [`Footing`] carries the two geometric inputs the strip-footing
//! Terzaghi equation needs — the footing width `B` and the depth `Df`
//! of the founding plane below the adjacent ground surface. Build it
//! through [`Footing::new`], which validates both.

use serde::{Deserialize, Serialize};

use crate::error::SoilBearingError;

/// Geometry of a shallow strip footing.
///
/// - `width` is the footing width `B` (length units, e.g. m).
/// - `depth` is the founding depth `Df` below grade (length units),
///   used to form the surcharge `q = gamma * Df`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Footing {
    width: f64,
    depth: f64,
}

impl Footing {
    /// Build a validated footing geometry.
    ///
    /// # Parameters
    ///
    /// - `width` (`B`) must be finite and strictly positive.
    /// - `depth` (`Df`) must be finite and non-negative (`0` models a
    ///   surface footing, which contributes no surcharge term).
    ///
    /// # Errors
    ///
    /// Returns [`SoilBearingError::NotFinite`] for any NaN/infinite
    /// input, or [`SoilBearingError::InvalidParameter`] for an
    /// out-of-domain value.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_soilbearing::Footing;
    ///
    /// let f = Footing::new(2.0, 1.5).unwrap();
    /// assert!((f.width() - 2.0).abs() < 1e-12);
    /// assert!((f.depth() - 1.5).abs() < 1e-12);
    ///
    /// assert!(Footing::new(0.0, 1.0).is_err()); // zero width
    /// assert!(Footing::new(2.0, -0.5).is_err()); // negative depth
    /// ```
    pub fn new(width: f64, depth: f64) -> Result<Self, SoilBearingError> {
        let width = SoilBearingError::require_finite("width", width)?;
        let depth = SoilBearingError::require_finite("depth", depth)?;

        if width <= 0.0 {
            return Err(SoilBearingError::invalid(
                "width",
                width,
                "footing width B must be strictly positive",
            ));
        }
        if depth < 0.0 {
            return Err(SoilBearingError::invalid(
                "depth",
                depth,
                "founding depth Df must be non-negative",
            ));
        }

        Ok(Footing { width, depth })
    }

    /// Footing width `B` (length units).
    pub fn width(&self) -> f64 {
        self.width
    }

    /// Founding depth `Df` below grade (length units).
    pub fn depth(&self) -> f64 {
        self.depth
    }

    /// Effective surcharge `q = gamma * Df` at the founding plane.
    ///
    /// `gamma` is the unit weight of the soil above the founding level;
    /// here we use the bearing soil's own `unit_weight`, the common
    /// single-stratum simplification.
    pub fn surcharge(&self, unit_weight: f64) -> f64 {
        unit_weight * self.depth
    }
}
