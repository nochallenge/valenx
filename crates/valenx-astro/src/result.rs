//! Output types for an ascent simulation.

use serde::{Deserialize, Serialize};

use crate::orbit::OrbitElements;

/// A discrete moment in the flight, recorded for plotting / inspection.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TrajectorySample {
    /// Mission-elapsed time (s).
    pub time: f64,
    /// Geometric altitude above the equatorial radius (m).
    pub altitude_m: f64,
    /// Great-circle downrange distance from the launch site (m).
    pub downrange_m: f64,
    /// Inertial speed (m/s).
    pub speed_inertial: f64,
    /// Speed relative to the co-rotating atmosphere (m/s).
    pub speed_relative: f64,
    /// Mach number relative to the local air.
    pub mach: f64,
    /// Current total mass (kg).
    pub mass: f64,
    /// Dynamic pressure (Pa).
    pub dynamic_pressure: f64,
    /// Sensed (non-gravitational) acceleration in g.
    pub acceleration_g: f64,
}

/// A discrete flight milestone (liftoff, staging, burnout, …).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FlightEvent {
    /// Mission-elapsed time (s).
    pub time: f64,
    /// Altitude at the event (m).
    pub altitude_m: f64,
    /// Inertial speed at the event (m/s).
    pub speed: f64,
    /// What happened.
    pub kind: String,
}

/// How a run ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Outcome {
    /// All stages burned out; the vehicle is on the reported orbit.
    Meco,
    /// The vehicle returned to the surface before all stages burned out.
    Impact,
    /// The simulated-time cap was reached first.
    TimedOut,
}

/// The full result of an ascent run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AscentResult {
    /// How the run terminated.
    pub outcome: Outcome,
    /// Liftoff gross mass (kg).
    pub liftoff_mass: f64,
    /// Ideal (Tsiolkovsky) `Δv` budget of the vehicle (m/s).
    pub ideal_delta_v: f64,
    /// Mission-elapsed time at main-engine cutoff / termination (s).
    pub final_time: f64,
    /// Altitude at termination (m).
    pub final_altitude_m: f64,
    /// Inertial speed at termination (m/s).
    pub final_speed_inertial: f64,
    /// Orbital elements of the trajectory at termination.
    pub orbit: OrbitElements,
    /// Peak dynamic pressure encountered (Pa).
    pub max_dynamic_pressure: f64,
    /// Altitude at which peak dynamic pressure occurred (m).
    pub max_q_altitude_m: f64,
    /// Peak sensed (non-gravitational) acceleration (g).
    pub max_acceleration_g: f64,
    /// True if the apoapsis cleared the Kármán line (100 km).
    pub reached_space: bool,
    /// True if the trajectory is bound with a periapsis above 100 km —
    /// a stable low orbit, not a re-entry arc.
    pub reached_orbit: bool,
    /// Flight milestones in chronological order.
    pub events: Vec<FlightEvent>,
    /// Down-sampled trajectory series.
    pub samples: Vec<TrajectorySample>,
    /// Insertion position at termination in the planar launch-plane
    /// frame (m): `[radial, downrange]`. Used to embed the in-plane
    /// trajectory into a 3-D orbital plane ([`crate::flight3d`]).
    pub final_position_m: [f64; 2],
    /// Insertion velocity at termination in the planar launch-plane
    /// frame (m/s): `[radial, downrange]`.
    pub final_velocity_ms: [f64; 2],
}

impl AscentResult {
    /// Convenience: apoapsis altitude in km (`f64::INFINITY` if unbound).
    pub fn apoapsis_km(&self) -> f64 {
        self.orbit.apoapsis_altitude / 1000.0
    }

    /// Convenience: periapsis altitude in km.
    pub fn periapsis_km(&self) -> f64 {
        self.orbit.periapsis_altitude / 1000.0
    }
}
