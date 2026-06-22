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
//!   ([`Aircraft::best_glide_ratio`], [`Aircraft::glide_range`]); and
//! - the lift-curve slope `dCL/dα` — the thin-airfoil 2D value `2π`
//!   ([`THIN_AIRFOIL_LIFT_SLOPE_PER_RAD`]) and its finite-wing
//!   correction `a = a0/(1 + a0/(π·e·AR))`
//!   ([`finite_wing_lift_slope_per_rad`], [`Aircraft::lift_curve_slope_per_rad`]).
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

/// Two-dimensional (infinite-span) lift-curve slope `a0 = dCl/dα` of a
/// thin airfoil, **per radian**, from thin-airfoil theory: exactly `2π`
/// (≈ 6.283 /rad, ≈ 0.10966 /deg). This is the inviscid small-angle
/// result for a thin section, independent of camber (camber only shifts
/// the zero-lift angle, not the slope).
///
/// Reference: Anderson, *Fundamentals of Aerodynamics*, thin-airfoil
/// theory (`cl = 2π(α − α_{L=0})`).
pub const THIN_AIRFOIL_LIFT_SLOPE_PER_RAD: f64 = 2.0 * PI;

/// The thin-airfoil 2D lift-curve slope expressed **per degree**:
/// `2π / 180·π = π²/90 ≈ 0.10966 /deg`.
///
/// Reference: Anderson, *Fundamentals of Aerodynamics*.
pub fn thin_airfoil_lift_slope_per_deg() -> f64 {
    THIN_AIRFOIL_LIFT_SLOPE_PER_RAD.to_radians()
}

/// Finite-wing lift-curve slope `a = dCL/dα` (per radian) from Prandtl
/// lifting-line theory:
///
/// ```text
/// a = a0 / (1 + a0 / (π · e · AR))
/// ```
///
/// where `a0` is the 2D section lift slope (per radian, e.g.
/// [`THIN_AIRFOIL_LIFT_SLOPE_PER_RAD`]), `aspect_ratio` is the wing
/// `AR = b²/S`, and `span_efficiency` is the span-efficiency factor `e`
/// (≈ 1 for an elliptic loading; here the Oswald-type factor in `(0, 1]`).
/// The finite wing always has a **shallower** slope than its 2D section
/// (`a < a0`), approaching `a0` only as `AR → ∞`.
///
/// Reference: Anderson, *Fundamentals of Aerodynamics* (Prandtl
/// lifting-line theory). Returns `0.0` for any non-positive / non-finite
/// argument.
pub fn finite_wing_lift_slope_per_rad(a0: f64, aspect_ratio: f64, span_efficiency: f64) -> f64 {
    if ![a0, aspect_ratio, span_efficiency]
        .iter()
        .all(|x| x.is_finite() && *x > 0.0)
    {
        return 0.0;
    }
    a0 / (1.0 + a0 / (PI * span_efficiency * aspect_ratio))
}

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

    /// This wing's finite-span lift-curve slope `dCL/dα` (per radian),
    /// from Prandtl lifting-line theory applied to a thin-airfoil 2D
    /// section (`a0 = 2π`):
    ///
    /// ```text
    /// a = 2π / (1 + 2π / (π · e · AR))
    /// ```
    ///
    /// using this aircraft's aspect ratio and span efficiency. Always
    /// less than the 2D `2π`, tending to it as `AR → ∞`. See
    /// [`finite_wing_lift_slope_per_rad`].
    ///
    /// Reference: Anderson, *Fundamentals of Aerodynamics*.
    pub fn lift_curve_slope_per_rad(&self) -> f64 {
        finite_wing_lift_slope_per_rad(
            THIN_AIRFOIL_LIFT_SLOPE_PER_RAD,
            self.aspect_ratio,
            self.oswald_efficiency,
        )
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

    /// Thin-airfoil theory: the 2D lift-curve slope is exactly `dCl/dα =
    /// 2π` per radian (≈ 0.10966 per degree). Reference: Anderson,
    /// *Fundamentals of Aerodynamics*, thin-airfoil theory.
    #[test]
    fn thin_airfoil_lift_slope_is_two_pi_per_radian() {
        // Per radian: exactly 2π (≈ 6.28319 /rad).
        assert_eq!(THIN_AIRFOIL_LIFT_SLOPE_PER_RAD, 2.0 * PI);
        // Per degree: 2π·(π/180) = π²/90 ≈ 0.109662.
        let per_deg = thin_airfoil_lift_slope_per_deg();
        assert!(
            (per_deg - 0.109_662).abs() < 1e-5,
            "got {per_deg} /deg, expected ≈0.10966"
        );

        // Cross-check by finite difference of cl = 2π·α over a small
        // angle window (the defining relation of thin-airfoil lift).
        let cl = |alpha_rad: f64| THIN_AIRFOIL_LIFT_SLOPE_PER_RAD * alpha_rad;
        let a1 = 1.0_f64.to_radians();
        let a2 = 3.0_f64.to_radians();
        let slope_fd = (cl(a2) - cl(a1)) / (a2 - a1);
        assert!((slope_fd - 2.0 * PI).abs() < 1e-9, "fd slope {slope_fd}");
    }

    /// Finite-wing (Prandtl lifting-line) lift slope
    /// `a = a0/(1 + a0/(π·e·AR))` must be below the 2D value and rise
    /// monotonically toward `2π` as aspect ratio grows. Reference:
    /// Anderson, *Fundamentals of Aerodynamics* (lifting-line theory).
    #[test]
    fn finite_wing_lift_slope_trends_toward_two_pi() {
        let a0 = THIN_AIRFOIL_LIFT_SLOPE_PER_RAD;

        // GA wing AR=7.5, e=0.8: a = 2π/(1 + 2π/(π·0.8·7.5)) = 4.71239 /rad.
        let a_ga = finite_wing_lift_slope_per_rad(a0, 7.5, 0.8);
        assert!((a_ga - 4.712_389).abs() < 1e-4, "GA slope {a_ga} /rad");

        // Strictly shallower than the 2D section.
        assert!(a_ga < a0);

        // Monotone increasing with aspect ratio, all below 2π.
        let mut prev = 0.0;
        for &ar in &[4.0_f64, 8.0, 16.0, 32.0, 1000.0] {
            let a = finite_wing_lift_slope_per_rad(a0, ar, 1.0);
            assert!(a > prev, "slope not increasing at AR={ar}: {a} <= {prev}");
            assert!(a < a0 + 1e-9, "slope {a} exceeded 2π at AR={ar}");
            prev = a;
        }
        // Approaches 2π in the high-AR limit.
        let a_huge = finite_wing_lift_slope_per_rad(a0, 1.0e6, 1.0);
        assert!(
            (a_huge - a0).abs() < 1e-3,
            "high-AR limit {a_huge} vs 2π {a0}"
        );

        // The Aircraft method agrees with the free function for that wing.
        let a = ga(); // AR=7.5, e=0.8
        assert!((a.lift_curve_slope_per_rad() - a_ga).abs() < 1e-12);

        // Guards: non-positive / non-finite arguments → 0.0.
        assert_eq!(finite_wing_lift_slope_per_rad(a0, 0.0, 1.0), 0.0);
        assert_eq!(finite_wing_lift_slope_per_rad(a0, 8.0, 0.0), 0.0);
        assert_eq!(finite_wing_lift_slope_per_rad(f64::NAN, 8.0, 1.0), 0.0);
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
