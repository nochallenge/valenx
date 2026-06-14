//! Newton's-law-of-cooling (surface convection) resistance.
//!
//! ## Model
//!
//! Convective heat transfer between a surface of area `A` and a fluid,
//! characterised by a convection (film) coefficient `h`, obeys Newton's
//! law of cooling `Q = h * A * ΔT`. Cast as a thermal resistance this
//! is
//!
//! ```text
//! R_conv = 1 / (h * A)        [K/W]
//! ```
//!
//! A larger surface area or a more vigorous film coefficient lowers the
//! convective resistance. This is the standard surface-convection
//! resistance used in series/parallel thermal circuits (Incropera
//! §3.1).

use serde::{Deserialize, Serialize};

use crate::error::{require_finite, require_positive, Result};

/// A convecting surface: area `A` exchanging heat with a fluid through
/// a film coefficient `h`.
///
/// SI units: `area_m2` in m², `h_w_per_m2k` in W/(m²·K).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ConvectiveSurface {
    /// Surface area `A` in contact with the fluid (m²).
    pub area_m2: f64,
    /// Convection (film) coefficient `h` (W/(m²·K)).
    pub h_w_per_m2k: f64,
}

impl ConvectiveSurface {
    /// Build a validated convecting surface.
    ///
    /// # Errors
    ///
    /// Returns [`HeatTransferError::BadParameter`](crate::HeatTransferError::BadParameter)
    /// if `area` or `h` is not finite and strictly positive.
    pub fn new(area_m2: f64, h_w_per_m2k: f64) -> Result<Self> {
        Ok(Self {
            area_m2: require_positive("area_m2", area_m2)?,
            h_w_per_m2k: require_positive("h_w_per_m2k", h_w_per_m2k)?,
        })
    }

    /// Convective thermal resistance `R = 1 / (h * A)` in K/W.
    pub fn resistance(&self) -> f64 {
        1.0 / (self.h_w_per_m2k * self.area_m2)
    }

    /// Convective heat rate `Q = h * A * ΔT` (W) for a surface/fluid
    /// temperature pair.
    ///
    /// Equivalent to `ΔT / R`. A positive result means heat leaves the
    /// (hotter) surface into the fluid.
    ///
    /// # Errors
    ///
    /// Returns an error if either temperature is non-finite.
    pub fn heat_rate(&self, t_surface: f64, t_fluid: f64) -> Result<f64> {
        let t_surface = require_finite("t_surface", t_surface)?;
        let t_fluid = require_finite("t_fluid", t_fluid)?;
        Ok(self.h_w_per_m2k * self.area_m2 * (t_surface - t_fluid))
    }
}

/// Free-function convective resistance `R = 1 / (h * A)` (K/W).
///
/// # Errors
///
/// Returns an error if `h` or `area` is not finite and strictly
/// positive.
pub fn convection_resistance(h_w_per_m2k: f64, area_m2: f64) -> Result<f64> {
    Ok(ConvectiveSurface::new(area_m2, h_w_per_m2k)?.resistance())
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn resistance_matches_closed_form() {
        // h = 25 W/m2K, A = 2 m^2 -> R = 1/(25*2) = 0.02 K/W.
        let surf = ConvectiveSurface::new(2.0, 25.0).unwrap();
        assert!((surf.resistance() - 0.02).abs() < EPS);
    }

    #[test]
    fn heat_rate_equals_delta_t_over_r() {
        let surf = ConvectiveSurface::new(2.0, 25.0).unwrap();
        let q_direct = surf.heat_rate(60.0, 20.0).unwrap();
        let q_via_r = (60.0 - 20.0) / surf.resistance();
        assert!((q_direct - q_via_r).abs() < 1e-9);
        // h*A*ΔT = 25*2*40 = 2000 W.
        assert!((q_direct - 2000.0).abs() < 1e-9);
    }

    #[test]
    fn larger_h_lowers_resistance() {
        let calm = ConvectiveSurface::new(1.0, 5.0).unwrap();
        let windy = ConvectiveSurface::new(1.0, 50.0).unwrap();
        assert!(windy.resistance() < calm.resistance());
    }

    #[test]
    fn rejects_non_positive_inputs() {
        assert!(ConvectiveSurface::new(0.0, 10.0).is_err());
        assert!(ConvectiveSurface::new(1.0, 0.0).is_err());
        assert!(convection_resistance(f64::INFINITY, 1.0).is_err());
    }
}
