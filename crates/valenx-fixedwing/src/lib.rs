//! # valenx-fixedwing
//!
//! Fixed-wing aircraft **preliminary point-performance** from a parabolic
//! drag polar.
//!
//! ## What
//!
//! Given a wing area, mass, maximum lift coefficient, aspect ratio,
//! zero-lift drag coefficient, Oswald efficiency and air density, this
//! crate computes:
//!
//! - lift `L = 0.5*rho*V^2*S*CL` ([`Aircraft::lift`]) and the wing loading
//!   `W/S` ([`Aircraft::wing_loading`]);
//! - the stall speed `Vs = sqrt(2W / (rho*S*CLmax))`
//!   ([`Aircraft::stall_speed`]) and the lift coefficient required for
//!   level flight at a speed ([`Aircraft::level_flight_cl`]);
//! - the parabolic drag polar `CD = CD0 + CL^2/(pi*AR*e)`
//!   ([`Aircraft::drag_coefficient`]);
//! - the maximum lift-to-drag ratio and the lift coefficient that achieves
//!   it ([`Aircraft::max_lift_to_drag`], [`Aircraft::cl_at_max_ld`]); and
//! - the best (still-air) glide ratio and glide range from an altitude
//!   ([`Aircraft::best_glide_ratio`], [`Aircraft::glide_range`]).
//!
//! ## Model
//!
//! Steady level flight balances lift against weight, `L = W = m*g`, with
//! `L = 0.5*rho*V^2*S*CL`. Setting `CL = CLmax` gives the minimum (stall)
//! speed `Vs = sqrt(2W / (rho*S*CLmax))`. Drag uses the textbook parabolic
//! polar
//!
//! `CD = CD0 + k*CL^2`,  `k = 1/(pi*AR*e)`,
//!
//! with `CD0` the zero-lift drag coefficient, `AR` the wing aspect ratio
//! and `e` the Oswald span efficiency. The lift-to-drag ratio `CL/CD` is
//! maximised at `CL_opt = sqrt(CD0/k) = sqrt(pi*AR*e*CD0)`, where
//! `CD = 2*CD0`, giving
//!
//! `(L/D)max = CL_opt / (2*CD0) = 0.5*sqrt(pi*AR*e / CD0)`.
//!
//! In an unpowered glide the best glide ratio equals `(L/D)max`, so from
//! altitude `h` the still-air glide range is `h * (L/D)max`.
//!
//! Units are SI: m^2 for `S`, kg for mass, kg/m^3 for `rho`, m/s for speed;
//! results are newtons, N/m^2, m/s and metres.
//!
//! ## Honest scope
//!
//! Research/educational grade. This is first-order steady-level / glide
//! point-performance from a single parabolic polar with constant
//! coefficients. It omits compressibility and Reynolds-number effects, the
//! CL-dependence of `CD0` (drag rise / separation near `CLmax`), propulsion
//! and thrust-required matching, weight and balance, manoeuvre and gust
//! loads, and all stability and control. It is not flight-grade and is not
//! a substitute for wind-tunnel / CFD data (see `valenx-aero`) or certified
//! flight-performance analysis.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

/// Standard gravity (m/s^2).
pub const GRAVITY: f64 = 9.80665;
/// Nominal sea-level air density (kg/m^3, ISA 15 C).
pub const SEA_LEVEL_AIR_DENSITY: f64 = 1.225;

/// An out-of-domain aircraft input.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum FixedWingError {
    /// A quantity that must be finite and strictly positive was not.
    #[error("{quantity} must be finite and positive, got {value}")]
    NonPositive {
        /// The offending quantity's name.
        quantity: &'static str,
        /// The offending value.
        value: f64,
    },
    /// The Oswald efficiency was outside the physical range `(0, 1]`.
    #[error("Oswald efficiency must be in (0, 1], got {value}")]
    OswaldEfficiency {
        /// The offending value.
        value: f64,
    },
}

/// Return `value` when finite and strictly positive, else an error.
fn require_positive(quantity: &'static str, value: f64) -> Result<f64, FixedWingError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(FixedWingError::NonPositive { quantity, value })
    }
}

/// A fixed-wing aircraft, validated on construction.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Aircraft {
    /// Reference wing area `S` (m^2).
    pub wing_area_m2: f64,
    /// All-up mass `m` (kg).
    pub mass_kg: f64,
    /// Maximum (clean) lift coefficient `CLmax`.
    pub cl_max: f64,
    /// Wing aspect ratio `AR` (span^2 / area).
    pub aspect_ratio: f64,
    /// Zero-lift drag coefficient `CD0`.
    pub cd0: f64,
    /// Oswald span efficiency `e` in `(0, 1]`.
    pub oswald_efficiency: f64,
    /// Air density `rho` (kg/m^3).
    pub air_density: f64,
}

impl Aircraft {
    /// Build a validated aircraft. Wing area, mass, `CLmax`, aspect ratio,
    /// `CD0` and air density must be finite and positive; the Oswald
    /// efficiency must lie in `(0, 1]`.
    ///
    /// # Errors
    ///
    /// Returns [`FixedWingError`] when any input is out of its physical
    /// domain.
    pub fn new(
        wing_area_m2: f64,
        mass_kg: f64,
        cl_max: f64,
        aspect_ratio: f64,
        cd0: f64,
        oswald_efficiency: f64,
        air_density: f64,
    ) -> Result<Self, FixedWingError> {
        let wing_area_m2 = require_positive("wing area", wing_area_m2)?;
        let mass_kg = require_positive("mass", mass_kg)?;
        let cl_max = require_positive("CLmax", cl_max)?;
        let aspect_ratio = require_positive("aspect ratio", aspect_ratio)?;
        let cd0 = require_positive("CD0", cd0)?;
        let air_density = require_positive("air density", air_density)?;
        if !(oswald_efficiency.is_finite() && oswald_efficiency > 0.0 && oswald_efficiency <= 1.0) {
            return Err(FixedWingError::OswaldEfficiency {
                value: oswald_efficiency,
            });
        }
        Ok(Self {
            wing_area_m2,
            mass_kg,
            cl_max,
            aspect_ratio,
            cd0,
            oswald_efficiency,
            air_density,
        })
    }

    /// Weight `W = m*g` (N).
    pub fn weight(&self) -> f64 {
        self.mass_kg * GRAVITY
    }

    /// Wing loading `W/S` (N/m^2).
    pub fn wing_loading(&self) -> f64 {
        self.weight() / self.wing_area_m2
    }

    /// Induced-drag factor `k = 1/(pi*AR*e)`.
    pub fn induced_drag_factor(&self) -> f64 {
        1.0 / (PI * self.aspect_ratio * self.oswald_efficiency)
    }

    /// Lift `L = 0.5*rho*V^2*S*CL` (N) at airspeed `v` (m/s) and lift
    /// coefficient `cl`.
    pub fn lift(&self, speed_m_s: f64, cl: f64) -> f64 {
        0.5 * self.air_density * speed_m_s * speed_m_s * self.wing_area_m2 * cl
    }

    /// Stall speed `Vs = sqrt(2W / (rho*S*CLmax))` (m/s).
    pub fn stall_speed(&self) -> f64 {
        (2.0 * self.weight() / (self.air_density * self.wing_area_m2 * self.cl_max)).sqrt()
    }

    /// Lift coefficient required to hold level flight at airspeed `v`:
    /// `CL = W / (0.5*rho*V^2*S)`.
    pub fn level_flight_cl(&self, speed_m_s: f64) -> f64 {
        self.weight() / (0.5 * self.air_density * speed_m_s * speed_m_s * self.wing_area_m2)
    }

    /// Drag coefficient from the parabolic polar `CD = CD0 + k*CL^2`.
    pub fn drag_coefficient(&self, cl: f64) -> f64 {
        self.cd0 + self.induced_drag_factor() * cl * cl
    }

    /// Lift coefficient at the best lift-to-drag ratio,
    /// `CL_opt = sqrt(CD0/k) = sqrt(pi*AR*e*CD0)`.
    pub fn cl_at_max_ld(&self) -> f64 {
        (self.cd0 / self.induced_drag_factor()).sqrt()
    }

    /// Maximum lift-to-drag ratio `(L/D)max = 0.5*sqrt(pi*AR*e / CD0)`.
    pub fn max_lift_to_drag(&self) -> f64 {
        0.5 * (PI * self.aspect_ratio * self.oswald_efficiency / self.cd0).sqrt()
    }

    /// Best (still-air) glide ratio, equal to `(L/D)max`.
    pub fn best_glide_ratio(&self) -> f64 {
        self.max_lift_to_drag()
    }

    /// Still-air glide range (m) from altitude `h`: `h * (L/D)max`.
    pub fn glide_range(&self, altitude_m: f64) -> f64 {
        altitude_m.max(0.0) * self.best_glide_ratio()
    }

    /// Compute the full performance report in one call.
    pub fn performance(&self) -> Performance {
        Performance {
            weight_n: self.weight(),
            wing_loading_pa: self.wing_loading(),
            stall_speed_m_s: self.stall_speed(),
            cl_at_max_ld: self.cl_at_max_ld(),
            max_lift_to_drag: self.max_lift_to_drag(),
        }
    }
}

/// An aircraft's computed point-performance.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Performance {
    /// Weight `W` (N).
    pub weight_n: f64,
    /// Wing loading `W/S` (N/m^2).
    pub wing_loading_pa: f64,
    /// Stall speed `Vs` (m/s).
    pub stall_speed_m_s: f64,
    /// Lift coefficient at the best L/D.
    pub cl_at_max_ld: f64,
    /// Maximum lift-to-drag ratio.
    pub max_lift_to_drag: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-6 * b.abs().max(1.0)
    }

    /// A light single-engine GA aircraft.
    fn ga() -> Aircraft {
        Aircraft::new(16.0, 1100.0, 1.6, 7.5, 0.027, 0.8, SEA_LEVEL_AIR_DENSITY).unwrap()
    }

    #[test]
    fn lift_at_stall_speed_with_clmax_equals_weight() {
        let a = ga();
        let vs = a.stall_speed();
        assert!(close(a.lift(vs, a.cl_max), a.weight()));
    }

    #[test]
    fn level_flight_cl_reproduces_weight() {
        let a = ga();
        let v = 40.0;
        let cl = a.level_flight_cl(v);
        assert!(close(a.lift(v, cl), a.weight()));
    }

    #[test]
    fn max_ld_matches_cl_over_cd_at_the_optimum() {
        let a = ga();
        let cl_opt = a.cl_at_max_ld();
        // At CL_opt the polar gives CD = 2*CD0, so L/D = CL_opt/(2*CD0).
        assert!(close(a.drag_coefficient(cl_opt), 2.0 * a.cd0));
        let ld = cl_opt / a.drag_coefficient(cl_opt);
        assert!(close(a.max_lift_to_drag(), ld));
    }

    #[test]
    fn ld_is_maximised_at_cl_opt() {
        let a = ga();
        let cl_opt = a.cl_at_max_ld();
        let ld_opt = cl_opt / a.drag_coefficient(cl_opt);
        for &cl in &[0.3_f64, 0.5, 0.9, 1.2] {
            let ld = cl / a.drag_coefficient(cl);
            assert!(ld <= ld_opt + 1e-9, "L/D at CL={cl} exceeded the optimum");
        }
    }

    #[test]
    fn glide_range_is_altitude_times_glide_ratio() {
        let a = ga();
        assert!(close(a.glide_range(3000.0), 3000.0 * a.best_glide_ratio()));
        // Higher aspect ratio -> better glide ratio.
        let glider =
            Aircraft::new(16.0, 1100.0, 1.6, 20.0, 0.027, 0.85, SEA_LEVEL_AIR_DENSITY).unwrap();
        assert!(glider.best_glide_ratio() > a.best_glide_ratio());
    }

    #[test]
    fn heavier_aircraft_stalls_faster() {
        let light = ga();
        let heavy =
            Aircraft::new(16.0, 1400.0, 1.6, 7.5, 0.027, 0.8, SEA_LEVEL_AIR_DENSITY).unwrap();
        assert!(heavy.stall_speed() > light.stall_speed());
    }

    #[test]
    fn rejects_out_of_domain_inputs() {
        assert!(Aircraft::new(0.0, 1100.0, 1.6, 7.5, 0.027, 0.8, SEA_LEVEL_AIR_DENSITY).is_err());
        assert!(Aircraft::new(16.0, 1100.0, 1.6, 7.5, 0.0, 0.8, SEA_LEVEL_AIR_DENSITY).is_err());
        assert!(Aircraft::new(16.0, 1100.0, 1.6, 7.5, 0.027, 1.5, SEA_LEVEL_AIR_DENSITY).is_err());
        assert!(Aircraft::new(16.0, 1100.0, 1.6, 7.5, 0.027, 0.0, SEA_LEVEL_AIR_DENSITY).is_err());
    }

    #[test]
    fn performance_is_serde_round_trippable() {
        let p = ga().performance();
        let json = serde_json::to_string(&p).unwrap();
        let back: Performance = serde_json::from_str(&json).unwrap();
        // Tolerance compare (serde_json float parse is ~1 ULP approximate).
        assert!(close(back.weight_n, p.weight_n));
        assert!(close(back.stall_speed_m_s, p.stall_speed_m_s));
        assert!(close(back.cl_at_max_ld, p.cl_at_max_ld));
        assert!(close(back.max_lift_to_drag, p.max_lift_to_drag));
    }
}
