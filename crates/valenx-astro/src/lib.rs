//! # valenx-astro
//!
//! A native **launch-vehicle ascent + trajectory simulator** — point a
//! rocket at the sky, fly it to orbit, and read back the engineering
//! answer: the orbit it reached, its `Δv` budget, peak dynamic
//! pressure, peak g-load and the full flight profile.
//!
//! ## What this is
//!
//! Describe a vehicle as a stack of [`vehicle::Stage`]s (dry mass,
//! propellant, thrust and `Isp` at sea level and vacuum) plus a payload
//! and a drag model, choose a [`guidance::GuidanceProgram`], and
//! [`simulate_ascent`] integrates the flight from the pad to orbital
//! insertion.
//!
//! The physics is a **planar 3-DOF point-mass** model in an
//! Earth-centred inertial frame:
//!
//! - **Inverse-square gravity** about a spherical Earth (WGS-84 `μ`).
//! - **US Standard Atmosphere 1976** ([`atmosphere`]) for density,
//!   pressure and the speed of sound (drives Mach-dependent drag and
//!   pressure-dependent thrust).
//! - **Mach-dependent aerodynamic drag** against the **co-rotating
//!   atmosphere**, so the Earth-rotation launch boost and the airflow
//!   the vehicle actually feels are both modelled.
//! - **Pressure-corrected thrust** interpolated between each stage's
//!   sea-level and vacuum ratings, with a constant propellant
//!   mass-flow rate and automatic **staging** when a stage runs dry.
//! - A classic **vertical-rise → pitch-kick → gravity-turn** steering
//!   law ([`guidance`]), or **closed-loop orbital insertion**
//!   ([`mission`]) that flies ascent → coast-to-apoapsis →
//!   circularisation burn to reach a near-circular target orbit.
//!
//! Integration is a fixed-step **RK4** ([`sim`]); the resulting state
//! is converted to Keplerian [`orbit`] elements (apoapsis / periapsis,
//! eccentricity, period, specific energy) so you can see exactly what
//! orbit — or re-entry arc — the vehicle is on at cutoff.
//!
//! ```
//! use valenx_astro::{presets, simulate_ascent};
//!
//! let vehicle = presets::two_stage_medium_lift();
//! let config = presets::leo_ascent_config();
//! let result = simulate_ascent(&vehicle, &config).expect("valid case");
//! println!(
//!     "apoapsis {:.0} km, periapsis {:.0} km, Δv budget {:.0} m/s",
//!     result.apoapsis_km(),
//!     result.periapsis_km(),
//!     result.ideal_delta_v,
//! );
//! ```
//!
//! ## Honest scope — a real v1, not a flight GNC stack
//!
//! Every model here is the genuine article — the atmosphere recovers
//! the standard tables, the rocket equation and orbital elements are
//! exact, the RK4 integrator conserves orbital energy, and the example
//! vehicle reaches a bound orbit. It is deliberately a **v1**:
//!
//! - **Planar (2-D) point-mass ascent**, embedded into a 3-D orbit by
//!   launch geometry ([`flight3d`]); the *powered* flight is not yet a
//!   native 3-D integrator, and there is no rigid-body attitude / 6-DOF
//!   or thrust-vector / aero-moment modelling.
//! - **Two guidance modes**: an open-loop gravity turn, and a
//!   closed-loop ascent → coast → circularise insertion that targets a
//!   circular altitude. Neither is a full powered-explicit-guidance
//!   (PEG) flight computer with continuous on-line targeting.
//! - **On-orbit** mechanics include 3-D classical elements and **J2**
//!   ([`orbit3d`]); winds ([`wind`]) perturb the ascent drag. Still no
//!   higher-order geopotential, drag/third-body decay, or geodetic
//!   (non-spherical) Earth.
//! - **First-order propulsion** ([`propulsion`]): ideal-rocket nozzle
//!   thrust/Isp and pressure interpolation; no finite-rate combustion,
//!   engine transients, or throttling schedule.
//!
//! None of those omissions makes the result meaningless — the orbit,
//! the `Δv` budget, max-Q and the staging timeline are all real
//! engineering numbers. Each is a documented, well-understood
//! extension on the way toward a fuller flight-mechanics suite.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod aero;
pub mod atmosphere;
pub mod budget;
pub mod config;
pub mod constants;
pub mod dynamics;
pub mod error;
pub mod flight3d;
pub mod flight6dof;
pub mod groundtrack;
pub mod guidance;
pub mod influence;
pub mod lambert;
pub mod landing;
pub mod launch;
pub mod maneuver;
pub mod mass;
pub mod mission;
pub mod orbit;
pub mod orbit3d;
pub mod presets;
pub mod propulsion;
pub mod recovery;
pub mod reentry;
pub mod rendezvous;
pub mod result;
pub mod rigidbody;
pub mod sim;
pub mod vehicle;
pub mod wind;
pub mod windows;

pub use config::{AscentConfig, GuidanceMode};
pub use error::AstroError;
pub use flight3d::{ascent_to_orbit, Ascent3d};
pub use flight6dof::{ControlGains, State6dof};
pub use guidance::GuidanceProgram;
pub use influence::{hill_sphere_radius, sphere_of_influence_radius};
pub use lambert::lambert;
pub use maneuver::{bielliptic_transfer, hohmann_transfer, Transfer};
pub use orbit::{elements, OrbitElements};
pub use orbit3d::{ClassicalElements, StateVector};
pub use propulsion::{EngineDesign, EnginePerformance};
pub use result::{AscentResult, FlightEvent, Outcome, TrajectorySample};
pub use rigidbody::{AttitudeState, Inertia};
pub use sim::{propagate_two_body, simulate_ascent};
pub use vehicle::{DragModel, Stage, Vehicle};
pub use wind::WindModel;

#[cfg(test)]
mod integration_tests {
    use super::*;

    #[test]
    fn medium_lift_reaches_orbit() {
        let vehicle = presets::two_stage_medium_lift();
        let config = presets::leo_ascent_config();
        let r = simulate_ascent(&vehicle, &config).expect("run");

        // The vehicle should finish on main-engine cutoff, not impact
        // or time out.
        assert_eq!(r.outcome, Outcome::Meco, "outcome {:?}", r.outcome);

        // Plenty of ideal Δv for orbit.
        assert!(r.ideal_delta_v > 9_400.0, "Δv budget {}", r.ideal_delta_v);

        // It should clear the Kármán line and be moving near orbital
        // speed at cutoff.
        assert!(r.reached_space, "apoapsis only {:.1} km", r.apoapsis_km());
        assert!(
            r.final_speed_inertial > 7_000.0,
            "MECO speed {:.0} m/s",
            r.final_speed_inertial
        );

        // The trajectory must be a genuine bound orbit with its
        // periapsis above the atmosphere — not a re-entry arc.
        assert!(r.reached_orbit, "periapsis only {:.1} km", r.periapsis_km());
        assert!(r.orbit.is_bound && r.orbit.eccentricity < 1.0);
        assert!(
            r.periapsis_km() > 100.0,
            "periapsis {:.1} km",
            r.periapsis_km()
        );

        // Max-Q must be positive and occur down in the lower
        // atmosphere (below ~30 km), as it does for real launchers.
        assert!(r.max_dynamic_pressure > 0.0);
        assert!(
            r.max_q_altitude_m < 30_000.0,
            "max-Q at {:.0} m",
            r.max_q_altitude_m
        );

        // Staging + MECO must both be recorded.
        assert!(r.events.iter().any(|e| e.kind.contains("Staging")));
        assert!(r.events.iter().any(|e| e.kind.contains("MECO")));

        // A non-trivial flight produces a trajectory series.
        assert!(r.samples.len() > 10);
    }

    #[test]
    fn underpowered_vehicle_does_not_reach_orbit() {
        // Strip the upper stage to a stub: nowhere near orbital Δv.
        let mut vehicle = presets::two_stage_medium_lift();
        vehicle.stages[1].propellant_mass = 1_000.0;
        let config = presets::leo_ascent_config();
        let r = simulate_ascent(&vehicle, &config).expect("run");
        assert!(!r.reached_orbit);
    }

    #[test]
    fn wind_perturbs_trajectory_and_raises_max_q() {
        // Still air vs. a strong steady wind across the whole lower
        // atmosphere. During the early near-vertical climb the vehicle
        // co-rotates with the air, so a steady wind shows up almost
        // entirely as extra air-relative speed -> higher max-Q and a
        // measurably different trajectory.
        let vehicle = presets::two_stage_medium_lift();
        let calm = presets::leo_ascent_config();
        let mut windy = calm;
        windy.wind = WindModel::Constant(180.0);

        let r_calm = simulate_ascent(&vehicle, &calm).expect("calm");
        let r_windy = simulate_ascent(&vehicle, &windy).expect("windy");

        assert!(
            r_windy.max_dynamic_pressure > r_calm.max_dynamic_pressure,
            "windy max-Q {} should exceed calm {}",
            r_windy.max_dynamic_pressure,
            r_calm.max_dynamic_pressure
        );
        // The trajectory is measurably perturbed.
        let dx = r_windy.final_position_m[0] - r_calm.final_position_m[0];
        let dy = r_windy.final_position_m[1] - r_calm.final_position_m[1];
        assert!((dx * dx + dy * dy).sqrt() > 1.0, "trajectory should differ");
        assert!(r_windy.reached_space);

        // Still air must reproduce Phase 0 exactly (no accidental drift).
        let r_calm2 = simulate_ascent(&vehicle, &presets::leo_ascent_config()).expect("calm2");
        assert_eq!(r_calm.max_dynamic_pressure, r_calm2.max_dynamic_pressure);
    }

    #[test]
    fn closed_loop_inserts_into_a_near_circular_orbit() {
        // The closed-loop insertion config flies ascent -> coast to
        // apoapsis -> circularisation burn, and should reach an almost
        // circular ~300 km orbit (a real LEO, not the eccentric orbit
        // the open-loop gravity turn settles into).
        let vehicle = presets::two_stage_medium_lift();
        let config = presets::leo_insertion_config();
        let r = simulate_ascent(&vehicle, &config).expect("run");

        assert_eq!(r.outcome, Outcome::Meco, "outcome {:?}", r.outcome);
        assert!(r.reached_orbit, "periapsis only {:.1} km", r.periapsis_km());

        // Near-circular: low eccentricity, apoapsis and periapsis both
        // close to the 300 km target.
        assert!(
            r.orbit.eccentricity < 0.02,
            "ecc {:.4}",
            r.orbit.eccentricity
        );
        assert!(
            (r.apoapsis_km() - 300.0).abs() < 60.0,
            "apoapsis {:.0} km",
            r.apoapsis_km()
        );
        assert!(
            (r.periapsis_km() - 300.0).abs() < 60.0,
            "periapsis {:.0} km",
            r.periapsis_km()
        );

        // The flight must show the coast (engine-off) then the
        // circularisation completion event.
        assert!(r
            .events
            .iter()
            .any(|e| e.kind.contains("Orbit insertion complete")));
    }
}
