//! Composite walls: series resistance, U-value, and heat-loss rate.
//!
//! A composite wall is an ordered stack of solid [`Layer`]s sandwiched
//! between an inner and outer [`SurfaceFilm`]. Because the heat flow is
//! modelled as one-dimensional and steady, the thermal resistances of
//! all elements add *in series*:
//!
//! `R_total = R_si + sum(L_i / k_i) + R_se`
//!
//! The area-specific thermal transmittance ("U-value") is the
//! reciprocal of the total resistance,
//!
//! `U = 1 / R_total`     (`W/(m^2.K)`)
//!
//! and the steady-state heat-loss rate through a wall of area `A`
//! (m^2) under a temperature difference `dT` (K) is
//!
//! `Q = U * A * dT`      (W)
//!
//! ## Honest scope
//!
//! This is a textbook one-dimensional steady-state model. It ignores
//! two-dimensional thermal bridging, moisture transport, air leakage,
//! and the temperature dependence of `k`. It is a teaching/estimation
//! aid, not a code-compliance energy model.

use serde::{Deserialize, Serialize};

use crate::error::InsulationError;
use crate::thermal::{Layer, SurfaceFilm};

/// An ordered, multi-layer wall assembly with optional inner and outer
/// surface films.
///
/// Build one with [`CompositeWall::builder`], push layers and films,
/// then call [`CompositeWall::total_resistance`],
/// [`CompositeWall::u_value`], or [`CompositeWall::heat_loss`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CompositeWall {
    /// Solid layers, ordered from the inner face to the outer face.
    layers: Vec<Layer>,
    /// Interior (inner-face) surface film, if modelled.
    interior_film: Option<SurfaceFilm>,
    /// Exterior (outer-face) surface film, if modelled.
    exterior_film: Option<SurfaceFilm>,
}

impl CompositeWall {
    /// Start building a wall. See [`CompositeWallBuilder`].
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_insulation::{CompositeWall, Layer};
    ///
    /// let wall = CompositeWall::builder()
    ///     .layer(Layer::new(0.2, 1.0).unwrap())
    ///     .build()
    ///     .unwrap();
    /// assert!((wall.total_resistance() - 0.2).abs() < 1e-12);
    /// ```
    pub fn builder() -> CompositeWallBuilder {
        CompositeWallBuilder::default()
    }

    /// The solid layers, inner face first.
    pub fn layers(&self) -> &[Layer] {
        &self.layers
    }

    /// The interior surface film, if one was added.
    pub fn interior_film(&self) -> Option<SurfaceFilm> {
        self.interior_film
    }

    /// The exterior surface film, if one was added.
    pub fn exterior_film(&self) -> Option<SurfaceFilm> {
        self.exterior_film
    }

    /// Total area-specific resistance, summed in series over the inner
    /// film, every solid layer, and the outer film:
    ///
    /// `R_total = R_si + sum(L_i / k_i) + R_se`   (`m^2.K/W`)
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_insulation::{CompositeWall, Layer, SurfaceFilm};
    ///
    /// // Two layers (R = 0.2 + 2.5) between R_si = 0.13 and R_se = 0.04.
    /// let wall = CompositeWall::builder()
    ///     .interior_film(SurfaceFilm::interior_default())
    ///     .layer(Layer::new(0.2, 1.0).unwrap())   // R = 0.2
    ///     .layer(Layer::new(0.1, 0.04).unwrap())  // R = 2.5
    ///     .exterior_film(SurfaceFilm::exterior_default())
    ///     .build()
    ///     .unwrap();
    /// let expected = 0.13 + 0.2 + 2.5 + 0.04;
    /// assert!((wall.total_resistance() - expected).abs() < 1e-3);
    /// ```
    pub fn total_resistance(&self) -> f64 {
        let mut r = 0.0;
        if let Some(f) = self.interior_film {
            r += f.resistance();
        }
        for layer in &self.layers {
            r += layer.resistance();
        }
        if let Some(f) = self.exterior_film {
            r += f.resistance();
        }
        r
    }

    /// Area-specific thermal transmittance `U = 1 / R_total`, in
    /// `W/(m^2.K)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_insulation::{CompositeWall, Layer};
    ///
    /// // A single R = 4.0 layer -> U = 0.25 W/(m^2.K).
    /// let wall = CompositeWall::builder()
    ///     .layer(Layer::new(0.16, 0.04).unwrap())
    ///     .build()
    ///     .unwrap();
    /// assert!((wall.u_value() - 0.25).abs() < 1e-9);
    /// ```
    pub fn u_value(&self) -> f64 {
        1.0 / self.total_resistance()
    }

    /// Steady-state heat-loss rate `Q = U * A * dT`, in watts.
    ///
    /// `area_m2` is the wall area (m^2) and `delta_t_k` is the
    /// inside-minus-outside temperature difference (K). Both must be
    /// finite and strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`InsulationError::NonPositive`] if `area_m2` or
    /// `delta_t_k` is not a finite, strictly positive number.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_insulation::{CompositeWall, Layer};
    ///
    /// // U = 0.25, A = 10 m^2, dT = 20 K -> Q = 50 W.
    /// let wall = CompositeWall::builder()
    ///     .layer(Layer::new(0.16, 0.04).unwrap())
    ///     .build()
    ///     .unwrap();
    /// let q = wall.heat_loss(10.0, 20.0).unwrap();
    /// assert!((q - 50.0).abs() < 1e-9);
    /// ```
    pub fn heat_loss(&self, area_m2: f64, delta_t_k: f64) -> Result<f64, InsulationError> {
        let area_m2 = InsulationError::require_positive("area_m2", area_m2)?;
        let delta_t_k = InsulationError::require_positive("delta_t_k", delta_t_k)?;
        Ok(self.u_value() * area_m2 * delta_t_k)
    }
}

/// Builder for a [`CompositeWall`].
///
/// Layers are appended inner-face first via [`CompositeWallBuilder::layer`].
/// Surface films are optional; set them with
/// [`CompositeWallBuilder::interior_film`] and
/// [`CompositeWallBuilder::exterior_film`]. Finalise with
/// [`CompositeWallBuilder::build`], which rejects an assembly that has
/// neither a layer nor a film.
#[derive(Clone, Debug, Default)]
pub struct CompositeWallBuilder {
    layers: Vec<Layer>,
    interior_film: Option<SurfaceFilm>,
    exterior_film: Option<SurfaceFilm>,
}

impl CompositeWallBuilder {
    /// Append one solid [`Layer`] to the stack (inner face first).
    pub fn layer(mut self, layer: Layer) -> Self {
        self.layers.push(layer);
        self
    }

    /// Set the interior (inner-face) surface film.
    pub fn interior_film(mut self, film: SurfaceFilm) -> Self {
        self.interior_film = Some(film);
        self
    }

    /// Set the exterior (outer-face) surface film.
    pub fn exterior_film(mut self, film: SurfaceFilm) -> Self {
        self.exterior_film = Some(film);
        self
    }

    /// Finalise the [`CompositeWall`].
    ///
    /// # Errors
    ///
    /// Returns [`InsulationError::EmptyAssembly`] if no layer and no
    /// surface film were added, since the wall would then have zero
    /// resistance and an infinite U-value.
    pub fn build(self) -> Result<CompositeWall, InsulationError> {
        if self.layers.is_empty() && self.interior_film.is_none() && self.exterior_film.is_none() {
            return Err(InsulationError::EmptyAssembly);
        }
        Ok(CompositeWall {
            layers: self.layers,
            interior_film: self.interior_film,
            exterior_film: self.exterior_film,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a single-layer wall with no films, returning resistance R.
    fn single_layer_wall(thickness: f64, k: f64) -> CompositeWall {
        CompositeWall::builder()
            .layer(Layer::new(thickness, k).unwrap())
            .build()
            .unwrap()
    }

    #[test]
    fn total_resistance_adds_layers_in_series() {
        // R = 0.2 + 2.5 = 2.7 m^2.K/W, no films.
        let wall = CompositeWall::builder()
            .layer(Layer::new(0.2, 1.0).unwrap()) // R = 0.2
            .layer(Layer::new(0.1, 0.04).unwrap()) // R = 2.5
            .build()
            .unwrap();
        assert!((wall.total_resistance() - 2.7).abs() < 1e-12);
    }

    #[test]
    fn total_resistance_includes_both_films() {
        // R = R_si + R_layer + R_se = 0.13 + 0.2 + 0.04 = 0.37.
        let wall = CompositeWall::builder()
            .interior_film(SurfaceFilm::interior_default())
            .layer(Layer::new(0.2, 1.0).unwrap())
            .exterior_film(SurfaceFilm::exterior_default())
            .build()
            .unwrap();
        assert!((wall.total_resistance() - 0.37).abs() < 1e-3);
    }

    #[test]
    fn series_resistance_equals_sum_of_parts() {
        // The series rule: R_total must equal the arithmetic sum of the
        // individual element resistances, exactly.
        let l1 = Layer::new(0.012, 0.25).unwrap(); // plasterboard, R ~ 0.048
        let l2 = Layer::new(0.10, 0.035).unwrap(); // insulation,  R ~ 2.857
        let l3 = Layer::new(0.10, 0.77).unwrap(); // brick,        R ~ 0.130
        let si = SurfaceFilm::interior_default();
        let se = SurfaceFilm::exterior_default();
        let wall = CompositeWall::builder()
            .interior_film(si)
            .layer(l1)
            .layer(l2)
            .layer(l3)
            .exterior_film(se)
            .build()
            .unwrap();
        let parts =
            si.resistance() + l1.resistance() + l2.resistance() + l3.resistance() + se.resistance();
        assert!((wall.total_resistance() - parts).abs() < 1e-12);
    }

    #[test]
    fn u_value_is_reciprocal_of_total_resistance() {
        // R = 4.0 -> U = 0.25 W/(m^2.K).
        let wall = single_layer_wall(0.16, 0.04); // R = 4.0
        assert!((wall.total_resistance() - 4.0).abs() < 1e-12);
        assert!((wall.u_value() - 0.25).abs() < 1e-12);
        // And U * R_total == 1 by construction.
        assert!((wall.u_value() * wall.total_resistance() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn heat_loss_is_u_times_area_times_delta_t() {
        // U = 0.25, A = 10, dT = 20 -> Q = 50 W.
        let wall = single_layer_wall(0.16, 0.04); // R = 4.0, U = 0.25
        let q = wall.heat_loss(10.0, 20.0).unwrap();
        assert!((q - 50.0).abs() < 1e-9);
    }

    #[test]
    fn heat_loss_scales_linearly_with_area_and_delta_t() {
        let wall = single_layer_wall(0.16, 0.04);
        let base = wall.heat_loss(10.0, 20.0).unwrap();
        let double_area = wall.heat_loss(20.0, 20.0).unwrap();
        let double_dt = wall.heat_loss(10.0, 40.0).unwrap();
        assert!((double_area - 2.0 * base).abs() < 1e-9);
        assert!((double_dt - 2.0 * base).abs() < 1e-9);
    }

    #[test]
    fn thicker_insulation_lowers_heat_loss() {
        // More insulation -> more R -> lower U -> less heat loss for the
        // same area and temperature difference.
        let thin = single_layer_wall(0.05, 0.04); // R = 1.25
        let thick = single_layer_wall(0.20, 0.04); // R = 5.0
        let q_thin = thin.heat_loss(10.0, 20.0).unwrap();
        let q_thick = thick.heat_loss(10.0, 20.0).unwrap();
        assert!(q_thick < q_thin);
    }

    #[test]
    fn lower_conductivity_lowers_heat_loss() {
        // Same geometry; a better insulator (lower k) loses less heat.
        let conductive = single_layer_wall(0.10, 0.50); // R = 0.2
        let insulating = single_layer_wall(0.10, 0.04); // R = 2.5
        let q_cond = conductive.heat_loss(10.0, 20.0).unwrap();
        let q_ins = insulating.heat_loss(10.0, 20.0).unwrap();
        assert!(q_ins < q_cond);
    }

    #[test]
    fn empty_assembly_is_rejected() {
        let err = CompositeWall::builder().build().unwrap_err();
        assert!(matches!(err, InsulationError::EmptyAssembly));
    }

    #[test]
    fn film_only_assembly_is_allowed() {
        // A bare interior film alone is a valid (if trivial) assembly.
        let wall = CompositeWall::builder()
            .interior_film(SurfaceFilm::interior_default())
            .build()
            .unwrap();
        assert!((wall.total_resistance() - 0.13).abs() < 1e-3);
    }

    #[test]
    fn heat_loss_rejects_non_positive_inputs() {
        let wall = single_layer_wall(0.16, 0.04);
        assert!(wall.heat_loss(0.0, 20.0).is_err());
        assert!(wall.heat_loss(10.0, 0.0).is_err());
        assert!(wall.heat_loss(-1.0, 20.0).is_err());
        assert!(wall.heat_loss(10.0, f64::NAN).is_err());
    }

    #[test]
    fn worked_example_typical_cavity_wall() {
        // A textbook-style filled-cavity wall, checked end to end.
        //   R_si  = 0.13
        //   plaster  12 mm, k 0.57  -> 0.02105...
        //   block   100 mm, k 0.15  -> 0.66667...
        //   EPS      80 mm, k 0.035 -> 2.28571...
        //   brick   100 mm, k 0.77  -> 0.12987...
        //   R_se  = 0.04
        let si = SurfaceFilm::interior_default();
        let se = SurfaceFilm::exterior_default();
        let plaster = Layer::new(0.012, 0.57).unwrap();
        let block = Layer::new(0.10, 0.15).unwrap();
        let eps = Layer::new(0.08, 0.035).unwrap();
        let brick = Layer::new(0.10, 0.77).unwrap();
        let wall = CompositeWall::builder()
            .interior_film(si)
            .layer(plaster)
            .layer(block)
            .layer(eps)
            .layer(brick)
            .exterior_film(se)
            .build()
            .unwrap();

        let expected_r = 0.13 + 0.012 / 0.57 + 0.10 / 0.15 + 0.08 / 0.035 + 0.10 / 0.77 + 0.04;
        assert!((wall.total_resistance() - expected_r).abs() < 1e-9);

        let expected_u = 1.0 / expected_r;
        assert!((wall.u_value() - expected_u).abs() < 1e-9);

        // 12 m^2 wall, 21 C inside / 1 C outside -> dT = 20 K.
        let expected_q = expected_u * 12.0 * 20.0;
        let q = wall.heat_loss(12.0, 20.0).unwrap();
        assert!((q - expected_q).abs() < 1e-9);
        // Sanity band: such a wall is roughly U ~ 0.33, Q ~ 78 W.
        assert!(wall.u_value() > 0.30 && wall.u_value() < 0.36);
        assert!(q > 70.0 && q < 85.0);
    }
}
