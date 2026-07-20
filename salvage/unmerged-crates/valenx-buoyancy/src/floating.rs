//! Flotation of a rectangular box (a simple barge / pontoon) on a free
//! surface.
//!
//! ## Governing equations
//!
//! A right rectangular prism of waterplane `length * beam` floating
//! upright sinks until the displaced volume supports its weight. With
//! mass `m`, fluid density `rho`, and waterplane area `A_wp = L * B`,
//! Archimedes gives the displaced volume `V = m / rho`, so the upright
//! draft (submerged depth) is
//!
//! ```text
//! V = m / rho
//! T = V / (L * B)
//! ```
//!
//! Equivalently the draft follows from the floating-body density ratio
//! `T = depth * (rho_body / rho_fluid)` when the box has a known
//! overall `depth` and mean density — the classic "fraction submerged
//! equals the density ratio" result. The box floats only if the
//! required draft does not exceed its depth; otherwise it is swamped.

use crate::error::{BuoyancyError, Result};

/// Geometry and loading of an upright floating rectangular box.
///
/// All lengths are in metres, `mass_kg` in kilograms. `depth_m` is the
/// box's full vertical extent (used only to decide whether it stays
/// afloat); the computed [`draft`](FloatingBox::draft) is how deep it
/// actually sits.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FloatingBox {
    /// Waterplane length (m), along the x axis.
    pub length_m: f64,
    /// Waterplane beam / width (m), along the y axis.
    pub beam_m: f64,
    /// Full vertical depth of the box (m).
    pub depth_m: f64,
    /// Total mass of the box and its contents (kg).
    pub mass_kg: f64,
}

impl FloatingBox {
    /// Construct and validate a [`FloatingBox`]. Every dimension and the
    /// mass must be finite and strictly positive.
    pub fn new(length_m: f64, beam_m: f64, depth_m: f64, mass_kg: f64) -> Result<Self> {
        Ok(FloatingBox {
            length_m: BuoyancyError::positive("length_m", length_m)?,
            beam_m: BuoyancyError::positive("beam_m", beam_m)?,
            depth_m: BuoyancyError::positive("depth_m", depth_m)?,
            mass_kg: BuoyancyError::positive("mass_kg", mass_kg)?,
        })
    }

    /// Waterplane area `A_wp = L * B` (m^2).
    pub fn waterplane_area(&self) -> f64 {
        self.length_m * self.beam_m
    }

    /// Volume of fluid the box must displace to float, `V = m / rho`
    /// (m^3), in a fluid of density `rho_fluid` (kg/m^3).
    ///
    /// `rho_fluid` must be finite and strictly positive.
    pub fn displaced_volume(&self, rho_fluid: f64) -> Result<f64> {
        let rho = BuoyancyError::positive("rho_fluid", rho_fluid)?;
        Ok(self.mass_kg / rho)
    }

    /// Upright draft `T = V / A_wp` (m) — how deep the box sits in a
    /// fluid of density `rho_fluid` (kg/m^3).
    ///
    /// `rho_fluid` must be finite and strictly positive.
    ///
    /// ```
    /// use valenx_buoyancy::floating::FloatingBox;
    /// // 10 m x 4 m x 2 m box, 40 000 kg, in fresh water:
    /// // V = 40 m^3, A = 40 m^2 -> draft exactly 1 m.
    /// let b = FloatingBox::new(10.0, 4.0, 2.0, 40_000.0).unwrap();
    /// let t = b.draft(1000.0).unwrap();
    /// assert!((t - 1.0).abs() < 1e-9);
    /// ```
    pub fn draft(&self, rho_fluid: f64) -> Result<f64> {
        let v = self.displaced_volume(rho_fluid)?;
        Ok(v / self.waterplane_area())
    }

    /// Mean density of the box, `rho_body = m / (L * B * depth)`
    /// (kg/m^3). For an upright box the fraction of its depth that is
    /// submerged equals `rho_body / rho_fluid`.
    pub fn mean_density(&self) -> f64 {
        self.mass_kg / (self.waterplane_area() * self.depth_m)
    }

    /// Whether the box floats upright without being swamped — i.e. the
    /// required [`draft`](FloatingBox::draft) is at most its
    /// `depth_m`. Equivalently, whether `rho_body <= rho_fluid`.
    ///
    /// `rho_fluid` must be finite and strictly positive.
    pub fn floats(&self, rho_fluid: f64) -> Result<bool> {
        let t = self.draft(rho_fluid)?;
        Ok(t <= self.depth_m)
    }

    /// Freeboard `depth - draft` (m): the unsubmerged height above the
    /// waterline. Negative if the box is swamped.
    ///
    /// `rho_fluid` must be finite and strictly positive.
    pub fn freeboard(&self, rho_fluid: f64) -> Result<f64> {
        let t = self.draft(rho_fluid)?;
        Ok(self.depth_m - t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn barge() -> FloatingBox {
        // 10 x 4 x 2 m, 40 t.
        FloatingBox::new(10.0, 4.0, 2.0, 40_000.0).unwrap()
    }

    #[test]
    fn waterplane_area_is_length_times_beam() {
        assert!((barge().waterplane_area() - 40.0).abs() < EPS);
    }

    #[test]
    fn displaced_volume_is_mass_over_density() {
        // 40 000 kg / 1000 = 40 m^3.
        assert!((barge().displaced_volume(1000.0).unwrap() - 40.0).abs() < EPS);
    }

    #[test]
    fn draft_known_value_fresh_water() {
        // V = 40 m^3, A = 40 m^2 -> T = 1 m exactly.
        assert!((barge().draft(1000.0).unwrap() - 1.0).abs() < EPS);
    }

    #[test]
    fn draft_equals_depth_times_density_ratio() {
        // Independent ground truth: submerged fraction = rho_body / rho_fluid.
        let b = barge();
        let rho_fluid = 1025.0;
        let predicted = b.depth_m * (b.mean_density() / rho_fluid);
        let actual = b.draft(rho_fluid).unwrap();
        assert!((actual - predicted).abs() < 1e-9, "{actual} vs {predicted}");
    }

    #[test]
    fn denser_fluid_gives_smaller_draft() {
        let b = barge();
        let t_fresh = b.draft(1000.0).unwrap();
        let t_sea = b.draft(1025.0).unwrap();
        assert!(t_sea < t_fresh);
        // Ratio of drafts is inverse ratio of densities.
        assert!((t_sea / t_fresh - 1000.0 / 1025.0).abs() < 1e-9);
    }

    #[test]
    fn mean_density_known_value() {
        // 40 000 / (40 * 2) = 500 kg/m^3.
        assert!((barge().mean_density() - 500.0).abs() < EPS);
    }

    #[test]
    fn floats_when_lighter_than_fluid_swamps_when_heavier() {
        // 500 kg/m^3 < 1000 -> floats with 1 m freeboard.
        let b = barge();
        assert!(b.floats(1000.0).unwrap());
        assert!((b.freeboard(1000.0).unwrap() - 1.0).abs() < EPS);

        // A box denser than the fluid cannot float: 10 x 4 x 2 m,
        // 100 t -> rho_body = 1250 > 1000.
        let heavy = FloatingBox::new(10.0, 4.0, 2.0, 100_000.0).unwrap();
        assert!(!heavy.floats(1000.0).unwrap());
        assert!(heavy.freeboard(1000.0).unwrap() < 0.0);
    }

    #[test]
    fn neutral_at_equal_density_sits_flush() {
        // rho_body == rho_fluid -> draft == depth, freeboard == 0
        // (the limiting flotation case).
        let b = FloatingBox::new(2.0, 2.0, 1.0, 4000.0).unwrap(); // rho_body = 1000
        assert!((b.mean_density() - 1000.0).abs() < EPS);
        assert!((b.draft(1000.0).unwrap() - b.depth_m).abs() < EPS);
        assert!(b.freeboard(1000.0).unwrap().abs() < EPS);
        assert!(b.floats(1000.0).unwrap()); // T <= depth holds with equality
    }

    #[test]
    fn constructor_and_methods_reject_bad_input() {
        assert!(FloatingBox::new(0.0, 4.0, 2.0, 1.0).is_err());
        assert!(FloatingBox::new(10.0, -4.0, 2.0, 1.0).is_err());
        assert!(FloatingBox::new(10.0, 4.0, 2.0, f64::INFINITY).is_err());
        assert!(barge().draft(0.0).is_err());
        assert!(barge().displaced_volume(-1.0).is_err());
    }
}
