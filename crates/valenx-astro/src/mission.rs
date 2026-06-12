//! The per-step flight controller that turns a [`GuidanceMode`] into a
//! concrete engine command.
//!
//! For the open-loop mode this is trivial — engines always on, fly the
//! gravity turn, let the sim's propellant bookkeeping cut off at
//! depletion. For closed-loop orbital insertion it is a small state
//! machine: **ascend** under the gravity turn until the osculating
//! apoapsis reaches the target, **coast** (engines off) to apoapsis,
//! then **circularise** with a prograde-horizontal burn until the
//! periapsis is raised to the target.

use nalgebra::Vector2;

use crate::config::GuidanceMode;
use crate::guidance::GuidanceProgram;
use crate::orbit;

/// Which direction the engine points when it is firing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThrustPolicy {
    /// Engine off (coast); direction is irrelevant.
    Off,
    /// Fly the open-loop gravity-turn pitch program.
    GravityTurn,
    /// Thrust prograde and horizontal (perpendicular to the radius, in
    /// the direction of motion) — the circularising direction.
    ProgradeHorizontal,
}

/// The command issued for one integration step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Command {
    /// Whether the engine should fire this step (subject to remaining
    /// propellant, which the sim enforces).
    pub engine_on: bool,
    /// The direction policy used to aim the thrust.
    pub policy: ThrustPolicy,
    /// When true the mission is complete (orbit achieved) and the run
    /// should terminate with a controlled cutoff.
    pub terminate: bool,
}

/// Phase of a closed-loop insertion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Phase {
    Ascent,
    Coast,
    Circularize,
}

/// Per-step flight controller.
#[derive(Debug, Clone)]
pub struct Controller {
    mode: GuidanceMode,
    program: GuidanceProgram,
    phase: Phase,
}

/// Minimum altitude (m) before a closed-loop ascent is allowed to cut
/// off for the coast — guards against cutting off while still in the
/// dense atmosphere if the apoapsis target is set very low.
const MIN_CUTOFF_ALTITUDE_M: f64 = 30_000.0;

/// How close to the target the periapsis must be raised before the
/// circularisation burn is considered done (m). A few km of slop keeps
/// the final eccentricity tiny without splitting integration steps.
const CIRCULARIZE_TOLERANCE_M: f64 = 5_000.0;

impl Controller {
    /// Create a controller for the given mode and ascent pitch program.
    pub fn new(mode: GuidanceMode, program: GuidanceProgram) -> Self {
        Self {
            mode,
            program,
            phase: Phase::Ascent,
        }
    }

    /// The unit thrust direction for a policy at a given state. Pure —
    /// safe to call repeatedly within an RK4 sub-step evaluation.
    pub fn direction(
        &self,
        policy: ThrustPolicy,
        position: Vector2<f64>,
        velocity: Vector2<f64>,
        t: f64,
    ) -> Vector2<f64> {
        match policy {
            ThrustPolicy::Off => safe_radial(position),
            ThrustPolicy::GravityTurn => self.program.thrust_direction(position, velocity, t),
            ThrustPolicy::ProgradeHorizontal => prograde_horizontal(position, velocity),
        }
    }

    /// Advance the controller one step and return the engine command.
    /// The current state drives the closed-loop phase transitions (the
    /// gravity-turn pitch program's time dependence lives in
    /// [`Controller::direction`]).
    pub fn command(&mut self, position: Vector2<f64>, velocity: Vector2<f64>) -> Command {
        let target_alt = match self.mode {
            GuidanceMode::OpenLoopGravityTurn => {
                return Command {
                    engine_on: true,
                    policy: ThrustPolicy::GravityTurn,
                    terminate: false,
                };
            }
            GuidanceMode::ClosedLoopInsertion { target_altitude_m } => target_altitude_m,
        };

        // Per-step closed-loop call on the integrated state, which is
        // finite and above the surface; use the validation-free core so
        // `command` stays infallible (matches the hot-path convention in
        // `rigidbody::propagate_unchecked`).
        let o = orbit::elements_with_mu_unchecked(position, velocity, crate::constants::MU_EARTH);
        let radial_hat = safe_radial(position);
        let altitude = position.norm() - crate::constants::R_EARTH;
        let radial_speed = velocity.dot(&radial_hat);

        match self.phase {
            Phase::Ascent => {
                // Cut off and coast once the apoapsis reaches the target
                // (and we are safely out of the dense atmosphere).
                if o.is_bound
                    && o.apoapsis_altitude >= target_alt
                    && altitude >= MIN_CUTOFF_ALTITUDE_M
                {
                    self.phase = Phase::Coast;
                    Command {
                        engine_on: false,
                        policy: ThrustPolicy::Off,
                        terminate: false,
                    }
                } else {
                    Command {
                        engine_on: true,
                        policy: ThrustPolicy::GravityTurn,
                        terminate: false,
                    }
                }
            }
            Phase::Coast => {
                // Coast (engines off) until apoapsis, i.e. until the
                // vehicle stops climbing (radial speed crosses zero).
                if radial_speed <= 0.0 {
                    self.phase = Phase::Circularize;
                    Command {
                        engine_on: true,
                        policy: ThrustPolicy::ProgradeHorizontal,
                        terminate: false,
                    }
                } else {
                    Command {
                        engine_on: false,
                        policy: ThrustPolicy::Off,
                        terminate: false,
                    }
                }
            }
            Phase::Circularize => {
                // Burn prograde-horizontal until the periapsis is raised
                // to the target — then the orbit is circular and we stop.
                if o.is_bound && o.periapsis_altitude >= target_alt - CIRCULARIZE_TOLERANCE_M {
                    Command {
                        engine_on: false,
                        policy: ThrustPolicy::Off,
                        terminate: true,
                    }
                } else {
                    Command {
                        engine_on: true,
                        policy: ThrustPolicy::ProgradeHorizontal,
                        terminate: false,
                    }
                }
            }
        }
    }
}

/// Radial ("up") unit vector, with a safe fallback at the origin.
fn safe_radial(position: Vector2<f64>) -> Vector2<f64> {
    let n = position.norm();
    if n > 1e-9 {
        position / n
    } else {
        Vector2::new(1.0, 0.0)
    }
}

/// Prograde-horizontal unit vector: the component of velocity
/// perpendicular to the radius, normalised. Falls back to the local
/// "east" (radial rotated +90°) if the velocity is purely radial.
pub fn prograde_horizontal(position: Vector2<f64>, velocity: Vector2<f64>) -> Vector2<f64> {
    let radial_hat = safe_radial(position);
    let v_tangential = velocity - velocity.dot(&radial_hat) * radial_hat;
    let n = v_tangential.norm();
    if n > 1e-6 {
        v_tangential / n
    } else {
        Vector2::new(-radial_hat.y, radial_hat.x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::R_EARTH;

    #[test]
    fn open_loop_always_burns_gravity_turn() {
        let mut c = Controller::new(
            GuidanceMode::OpenLoopGravityTurn,
            GuidanceProgram::default(),
        );
        let cmd = c.command(
            Vector2::new(R_EARTH + 1000.0, 0.0),
            Vector2::new(0.0, 100.0),
        );
        assert!(cmd.engine_on);
        assert_eq!(cmd.policy, ThrustPolicy::GravityTurn);
        assert!(!cmd.terminate);
    }

    #[test]
    fn prograde_horizontal_is_perpendicular_to_radius() {
        let pos = Vector2::new(R_EARTH + 300_000.0, 0.0);
        let vel = Vector2::new(50.0, 7_000.0); // mostly tangential + a little radial
        let dir = prograde_horizontal(pos, vel);
        // Perpendicular to the radius (here +x), so x-component ~ 0.
        assert!(dir.x.abs() < 1e-9, "dir {dir:?}");
        // Same tangential sense as the motion (+y).
        assert!(dir.y > 0.0);
        assert!((dir.norm() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn closed_loop_walks_ascent_to_coast_to_circularize() {
        let mut c = Controller::new(
            GuidanceMode::ClosedLoopInsertion {
                target_altitude_m: 300_000.0,
            },
            GuidanceProgram::default(),
        );

        // Low + slow -> still ascending under power.
        let lo = c.command(
            Vector2::new(R_EARTH + 10_000.0, 0.0),
            Vector2::new(0.0, 200.0),
        );
        assert!(lo.engine_on && lo.policy == ThrustPolicy::GravityTurn);

        // High altitude, apoapsis already past target, still climbing
        // -> cut off and coast.
        let r = R_EARTH + 150_000.0;
        // Choose a tangential speed giving an apoapsis well above 300 km
        // while radial speed is still positive (climbing).
        let pos = Vector2::new(r, 0.0);
        let vel = Vector2::new(500.0, 8_200.0);
        let coast = c.command(pos, vel);
        assert!(!coast.engine_on && coast.policy == ThrustPolicy::Off);

        // Now at apoapsis (radial speed <= 0) -> circularisation burn.
        let at_apo = c.command(Vector2::new(r, 0.0), Vector2::new(-1.0, 7_500.0));
        assert!(at_apo.engine_on && at_apo.policy == ThrustPolicy::ProgradeHorizontal);
    }
}
