//! Plane-wall (1D) steady conduction.
//!
//! ## Model
//!
//! For one-dimensional, steady conduction through a homogeneous plane
//! wall of thickness `L`, cross-sectional area `A` and thermal
//! conductivity `k` (Fourier's law, no heat generation), the
//! conductive **thermal resistance** is
//!
//! ```text
//! R_cond = L / (k * A)        [K/W]
//! ```
//!
//! The heat rate driven by a temperature difference `ΔT` across that
//! resistance is `Q = ΔT / R_cond`, and the temperature varies
//! **linearly** with position through the wall:
//!
//! ```text
//! T(x) = T_hot - (T_hot - T_cold) * x / L
//! ```
//!
//! with `x` measured from the hot face (`x = 0`) to the cold face
//! (`x = L`).
//!
//! These are the standard closed-form results in any heat-transfer
//! text (Incropera §3.1, Cengel §3-1); they hold exactly only for the
//! idealised 1D, steady, constant-`k`, no-generation case.

use serde::{Deserialize, Serialize};

use crate::error::{require_finite, require_positive, Result};

/// A homogeneous plane wall (slab) carrying 1D steady conduction.
///
/// All quantities are SI: `thickness_m` in metres, `area_m2` in m²,
/// `conductivity_w_per_mk` in W/(m·K).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PlaneWall {
    /// Wall thickness `L` along the conduction direction (m).
    pub thickness_m: f64,
    /// Cross-sectional area `A` normal to the heat flow (m²).
    pub area_m2: f64,
    /// Thermal conductivity `k` (W/(m·K)).
    pub conductivity_w_per_mk: f64,
}

impl PlaneWall {
    /// Build a validated plane wall.
    ///
    /// # Errors
    ///
    /// Returns [`HeatTransferError::BadParameter`](crate::HeatTransferError::BadParameter)
    /// if any of `thickness`, `area` or `conductivity` is not finite and
    /// strictly positive.
    pub fn new(thickness_m: f64, area_m2: f64, conductivity_w_per_mk: f64) -> Result<Self> {
        Ok(Self {
            thickness_m: require_positive("thickness_m", thickness_m)?,
            area_m2: require_positive("area_m2", area_m2)?,
            conductivity_w_per_mk: require_positive(
                "conductivity_w_per_mk",
                conductivity_w_per_mk,
            )?,
        })
    }

    /// Conductive thermal resistance `R = L / (k * A)` in K/W.
    ///
    /// Increasing thickness (more insulation) raises the resistance;
    /// increasing conductivity or area lowers it.
    pub fn resistance(&self) -> f64 {
        self.thickness_m / (self.conductivity_w_per_mk * self.area_m2)
    }

    /// Steady conductive heat rate `Q = ΔT / R` (W) driven by a hot/cold
    /// face temperature pair.
    ///
    /// A positive result means heat flows from `t_hot` to `t_cold`.
    ///
    /// # Errors
    ///
    /// Returns an error if either temperature is non-finite.
    pub fn heat_rate(&self, t_hot: f64, t_cold: f64) -> Result<f64> {
        let t_hot = require_finite("t_hot", t_hot)?;
        let t_cold = require_finite("t_cold", t_cold)?;
        Ok((t_hot - t_cold) / self.resistance())
    }

    /// Temperature `T(x)` at position `x` measured from the hot face.
    ///
    /// The profile is linear: `T(0) = t_hot` and `T(L) = t_cold`.
    ///
    /// # Errors
    ///
    /// Returns an error if `x` is non-finite or outside `[0, L]`, or if
    /// either face temperature is non-finite.
    pub fn temperature_at(&self, x_m: f64, t_hot: f64, t_cold: f64) -> Result<f64> {
        let x = require_finite("x_m", x_m)?;
        let t_hot = require_finite("t_hot", t_hot)?;
        let t_cold = require_finite("t_cold", t_cold)?;
        if x < 0.0 || x > self.thickness_m {
            return Err(crate::HeatTransferError::BadParameter {
                name: "x_m",
                reason: format!("must lie in [0, {L}], got {x}", L = self.thickness_m),
            });
        }
        Ok(t_hot - (t_hot - t_cold) * x / self.thickness_m)
    }

    /// Sample the linear temperature profile at `n` equally spaced
    /// stations from the hot face (`x = 0`) to the cold face (`x = L`).
    ///
    /// Returns `n` pairs `(x, T(x))`. With `n = 1` only the hot face is
    /// returned.
    ///
    /// # Errors
    ///
    /// Returns an error if `n == 0` or if either face temperature is
    /// non-finite.
    pub fn profile(&self, n: usize, t_hot: f64, t_cold: f64) -> Result<Vec<(f64, f64)>> {
        if n == 0 {
            return Err(crate::HeatTransferError::BadParameter {
                name: "n",
                reason: "must be >= 1".to_string(),
            });
        }
        let t_hot = require_finite("t_hot", t_hot)?;
        let t_cold = require_finite("t_cold", t_cold)?;
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            // Guard the n == 1 case (denominator would be zero).
            let frac = if n == 1 {
                0.0
            } else {
                i as f64 / (n as f64 - 1.0)
            };
            let x = frac * self.thickness_m;
            let t = t_hot - (t_hot - t_cold) * frac;
            out.push((x, t));
        }
        Ok(out)
    }
}

/// Free-function conductive resistance `R = L / (k * A)` (K/W).
///
/// Thin wrapper over [`PlaneWall::resistance`] for callers that only
/// need the scalar.
///
/// # Errors
///
/// Returns an error if any argument is not finite and strictly
/// positive.
pub fn conduction_resistance(
    thickness_m: f64,
    conductivity_w_per_mk: f64,
    area_m2: f64,
) -> Result<f64> {
    Ok(PlaneWall::new(thickness_m, area_m2, conductivity_w_per_mk)?.resistance())
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn resistance_matches_closed_form() {
        // L = 0.1 m, k = 2 W/mK, A = 4 m^2 -> R = 0.1 / (2*4) = 0.0125 K/W.
        let wall = PlaneWall::new(0.1, 4.0, 2.0).unwrap();
        assert!((wall.resistance() - 0.0125).abs() < EPS);
    }

    #[test]
    fn heat_rate_is_delta_t_over_r() {
        let wall = PlaneWall::new(0.1, 4.0, 2.0).unwrap();
        // Q = (100-20)/0.0125 = 6400 W.
        let q = wall.heat_rate(100.0, 20.0).unwrap();
        assert!((q - 6400.0).abs() < 1e-9);
    }

    #[test]
    fn thicker_wall_lowers_heat_loss() {
        // Same k, A, ΔT: a thicker (more-insulating) wall conducts less.
        let thin = PlaneWall::new(0.05, 1.0, 0.04).unwrap();
        let thick = PlaneWall::new(0.20, 1.0, 0.04).unwrap();
        let q_thin = thin.heat_rate(30.0, 0.0).unwrap();
        let q_thick = thick.heat_rate(30.0, 0.0).unwrap();
        assert!(q_thick < q_thin);
        // Quadrupling thickness quarters the heat loss.
        assert!((q_thin / q_thick - 4.0).abs() < 1e-9);
    }

    #[test]
    fn profile_is_linear_and_hits_both_faces() {
        let wall = PlaneWall::new(0.2, 1.0, 1.0).unwrap();
        let mid = wall.temperature_at(0.1, 100.0, 0.0).unwrap();
        // Midplane of a linear profile is the mean of the faces.
        assert!((mid - 50.0).abs() < 1e-9);
        assert!((wall.temperature_at(0.0, 100.0, 0.0).unwrap() - 100.0).abs() < 1e-9);
        assert!((wall.temperature_at(0.2, 100.0, 0.0).unwrap() - 0.0).abs() < 1e-9);
    }

    #[test]
    fn profile_samples_endpoints() {
        let wall = PlaneWall::new(1.0, 1.0, 1.0).unwrap();
        let pts = wall.profile(5, 80.0, 20.0).unwrap();
        assert_eq!(pts.len(), 5);
        assert!((pts[0].0 - 0.0).abs() < EPS);
        assert!((pts[0].1 - 80.0).abs() < EPS);
        assert!((pts[4].0 - 1.0).abs() < EPS);
        assert!((pts[4].1 - 20.0).abs() < EPS);
        // Even spacing in x.
        for i in 0..4 {
            let dx = pts[i + 1].0 - pts[i].0;
            assert!((dx - 0.25).abs() < 1e-9);
        }
    }

    #[test]
    fn rejects_non_positive_inputs() {
        assert!(PlaneWall::new(0.0, 1.0, 1.0).is_err());
        assert!(PlaneWall::new(1.0, -1.0, 1.0).is_err());
        assert!(PlaneWall::new(1.0, 1.0, f64::NAN).is_err());
    }

    #[test]
    fn temperature_outside_wall_is_rejected() {
        let wall = PlaneWall::new(0.2, 1.0, 1.0).unwrap();
        assert!(wall.temperature_at(-0.01, 100.0, 0.0).is_err());
        assert!(wall.temperature_at(0.3, 100.0, 0.0).is_err());
    }
}
