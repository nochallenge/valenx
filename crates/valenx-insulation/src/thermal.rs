//! Single-layer and surface-film thermal resistances.
//!
//! The fundamental unit of every building heat-loss calculation is the
//! *area-specific thermal resistance* (the "R-value"), measured in
//! `m^2.K/W`. For a solid layer of thickness `L` (metres) and thermal
//! conductivity `k` (`W/(m.K)`) the one-dimensional steady-state
//! conduction resistance is
//!
//! `R = L / k`
//!
//! For a surface exposed to a fluid (air) the convective/radiative
//! boundary is captured by a *surface film coefficient* `h`
//! (`W/(m^2.K)`); the film's resistance is its reciprocal
//!
//! `R_film = 1 / h`
//!
//! Both kinds of resistance live in the same `m^2.K/W` units and add in
//! series (see [`crate::wall`]).

use serde::{Deserialize, Serialize};

use crate::error::InsulationError;

/// A single solid layer of a wall assembly: a homogeneous material of
/// known thickness and thermal conductivity.
///
/// Construct with [`Layer::new`], which validates that both inputs are
/// finite and strictly positive. The area-specific conduction
/// resistance is then [`Layer::resistance`].
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Layer {
    /// Layer thickness `L` in metres.
    pub thickness_m: f64,
    /// Thermal conductivity `k` in `W/(m.K)`.
    pub conductivity_w_per_m_k: f64,
}

impl Layer {
    /// Build a [`Layer`] from a thickness (m) and a thermal
    /// conductivity (`W/(m.K)`).
    ///
    /// # Errors
    ///
    /// Returns [`InsulationError::NonPositive`] if either argument is
    /// not a finite, strictly positive number.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_insulation::Layer;
    ///
    /// // 100 mm of EPS foam, k = 0.035 W/(m.K).
    /// let l = Layer::new(0.10, 0.035).unwrap();
    /// assert!((l.resistance() - (0.10 / 0.035)).abs() < 1e-12);
    /// ```
    pub fn new(thickness_m: f64, conductivity_w_per_m_k: f64) -> Result<Self, InsulationError> {
        let thickness_m = InsulationError::require_positive("thickness_m", thickness_m)?;
        let conductivity_w_per_m_k =
            InsulationError::require_positive("conductivity_w_per_m_k", conductivity_w_per_m_k)?;
        Ok(Self {
            thickness_m,
            conductivity_w_per_m_k,
        })
    }

    /// Area-specific conduction resistance `R = L / k`, in `m^2.K/W`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_insulation::Layer;
    ///
    /// let l = Layer::new(0.2, 0.5).unwrap();
    /// assert!((l.resistance() - 0.4).abs() < 1e-12);
    /// ```
    pub fn resistance(&self) -> f64 {
        self.thickness_m / self.conductivity_w_per_m_k
    }
}

/// A convective/radiative surface film on the inside or outside face of
/// a wall, described by its film (surface) coefficient `h`.
///
/// Construct with [`SurfaceFilm::new`] from `h`, or with the
/// convenience constructors [`SurfaceFilm::interior_default`] /
/// [`SurfaceFilm::exterior_default`] which use the common ISO 6946
/// horizontal-heat-flow reference coefficients. The film's
/// area-specific resistance is [`SurfaceFilm::resistance`].
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct SurfaceFilm {
    /// Surface film coefficient `h` in `W/(m^2.K)`.
    pub coefficient_w_per_m2_k: f64,
}

impl SurfaceFilm {
    /// Build a [`SurfaceFilm`] from a film coefficient `h`
    /// (`W/(m^2.K)`).
    ///
    /// # Errors
    ///
    /// Returns [`InsulationError::NonPositive`] if `h` is not a finite,
    /// strictly positive number.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_insulation::SurfaceFilm;
    ///
    /// let f = SurfaceFilm::new(8.0).unwrap();
    /// assert!((f.resistance() - 0.125).abs() < 1e-12);
    /// ```
    pub fn new(coefficient_w_per_m2_k: f64) -> Result<Self, InsulationError> {
        let coefficient_w_per_m2_k =
            InsulationError::require_positive("coefficient_w_per_m2_k", coefficient_w_per_m2_k)?;
        Ok(Self {
            coefficient_w_per_m2_k,
        })
    }

    /// Interior surface film with the ISO 6946 reference coefficient for
    /// horizontal heat flow, `h = 7.69 W/(m^2.K)` (`R_si = 0.13
    /// m^2.K/W`).
    ///
    /// This is a fixed textbook reference value and is always valid, so
    /// it is infallible.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_insulation::SurfaceFilm;
    ///
    /// let si = SurfaceFilm::interior_default();
    /// assert!((si.resistance() - 0.13).abs() < 1e-3);
    /// ```
    pub fn interior_default() -> Self {
        // h chosen so that 1/h = 0.13 m^2.K/W (ISO 6946 R_si).
        Self {
            coefficient_w_per_m2_k: 1.0 / 0.13,
        }
    }

    /// Exterior surface film with the ISO 6946 reference coefficient for
    /// horizontal heat flow, `h = 25 W/(m^2.K)` (`R_se = 0.04
    /// m^2.K/W`).
    ///
    /// This is a fixed textbook reference value and is always valid, so
    /// it is infallible.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_insulation::SurfaceFilm;
    ///
    /// let se = SurfaceFilm::exterior_default();
    /// assert!((se.resistance() - 0.04).abs() < 1e-12);
    /// ```
    pub fn exterior_default() -> Self {
        Self {
            coefficient_w_per_m2_k: 25.0,
        }
    }

    /// Area-specific surface resistance `R = 1 / h`, in `m^2.K/W`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_insulation::SurfaceFilm;
    ///
    /// let f = SurfaceFilm::new(4.0).unwrap();
    /// assert!((f.resistance() - 0.25).abs() < 1e-12);
    /// ```
    pub fn resistance(&self) -> f64 {
        1.0 / self.coefficient_w_per_m2_k
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layer_resistance_is_l_over_k() {
        // Ground truth: R = L/k. 0.05 m / 0.04 W/(m.K) = 1.25 m^2.K/W.
        let l = Layer::new(0.05, 0.04).unwrap();
        assert!((l.resistance() - 1.25).abs() < 1e-12);
    }

    #[test]
    fn thicker_layer_has_more_resistance() {
        // Same material; double the thickness must double the R-value.
        let thin = Layer::new(0.05, 0.04).unwrap();
        let thick = Layer::new(0.10, 0.04).unwrap();
        assert!(thick.resistance() > thin.resistance());
        assert!((thick.resistance() - 2.0 * thin.resistance()).abs() < 1e-12);
    }

    #[test]
    fn lower_conductivity_has_more_resistance() {
        // Same thickness; a better insulator (lower k) gives a higher R.
        let conductor = Layer::new(0.05, 0.50).unwrap();
        let insulator = Layer::new(0.05, 0.04).unwrap();
        assert!(insulator.resistance() > conductor.resistance());
        // R scales as 1/k: (0.50/0.04) = 12.5x.
        assert!((insulator.resistance() / conductor.resistance() - 12.5).abs() < 1e-9);
    }

    #[test]
    fn film_resistance_is_reciprocal_of_h() {
        // Ground truth: R_film = 1/h. h = 8 -> R = 0.125 m^2.K/W.
        let f = SurfaceFilm::new(8.0).unwrap();
        assert!((f.resistance() - 0.125).abs() < 1e-12);
    }

    #[test]
    fn iso_default_films_have_reference_resistances() {
        // ISO 6946 horizontal-flow reference: R_si = 0.13, R_se = 0.04.
        let si = SurfaceFilm::interior_default();
        let se = SurfaceFilm::exterior_default();
        assert!((si.resistance() - 0.13).abs() < 1e-3);
        assert!((se.resistance() - 0.04).abs() < 1e-12);
    }

    #[test]
    fn constructors_reject_non_positive_inputs() {
        assert!(Layer::new(0.0, 0.04).is_err());
        assert!(Layer::new(0.05, 0.0).is_err());
        assert!(Layer::new(-0.05, 0.04).is_err());
        assert!(SurfaceFilm::new(0.0).is_err());
        assert!(SurfaceFilm::new(-1.0).is_err());
        assert!(SurfaceFilm::new(f64::INFINITY).is_err());
    }
}
