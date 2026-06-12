//! The ascent simulation loop: a fixed-step RK4 integrator over the
//! planar point-mass equations of motion, with discrete staging events
//! and pressure-dependent thrust.

use nalgebra::Vector2;

use crate::config::AscentConfig;
use crate::constants::{G0, R_EARTH};
use crate::dynamics;
use crate::error::AstroError;
use crate::mission::Controller;
use crate::orbit;
use crate::result::{AscentResult, FlightEvent, Outcome, TrajectorySample};
use crate::vehicle::Vehicle;

/// Absolute ceiling on the number of integration steps any single
/// simulation or propagation call may take. A fixed-step run derives its
/// step count from caller-supplied `max_time / dt`; this caps that so an
/// absurd configuration (tiny `dt`, huge `max_time`, or an enormous
/// explicit step count) is rejected up front instead of hanging.
///
/// 100 million steps is far beyond any physical ascent or orbit
/// propagation (a 4-month flight at a 0.1 ms step is ~3.5e10 steps, well
/// past this) yet finite, so the loop is always bounded.
pub const MAX_SIM_STEPS: u64 = 100_000_000;

/// Absolute ceiling on the number of recorded trajectory samples a run
/// may retain, bounding the output `Vec`'s memory.
pub const MAX_SAMPLES: u64 = 10_000_000;

/// Reject a propagation step count that exceeds [`MAX_SIM_STEPS`].
///
/// Shared by every fixed-step propagator so an absurd `steps` argument
/// (up to `u64::MAX`) is refused up front rather than silently truncated
/// or run to a hang.
pub(crate) fn check_step_count(steps: u64) -> Result<(), AstroError> {
    if steps > MAX_SIM_STEPS {
        return Err(AstroError::OutOfRange {
            what: "steps",
            value: steps,
            max: MAX_SIM_STEPS,
        });
    }
    Ok(())
}

/// The continuously-integrated part of the state (position + velocity);
/// mass and stage bookkeeping are stepped discretely alongside it.
#[derive(Clone, Copy, Debug)]
struct State {
    pos: Vector2<f64>,
    vel: Vector2<f64>,
}

/// Simulate a launch ascent and return the full flight record.
///
/// The vehicle lifts off due-east from the equator, flies the
/// configured open-loop pitch program, stages automatically as each
/// stage's propellant is exhausted, and the run ends at main-engine
/// cutoff (all stages burned), at surface impact, or at the simulated-
/// time cap — whichever comes first. The reported [`crate::orbit`]
/// elements describe the conic the vehicle is on at termination.
pub fn simulate_ascent(
    vehicle: &Vehicle,
    config: &AscentConfig,
) -> Result<AscentResult, AstroError> {
    vehicle.validate()?;
    config.validate()?;

    let dt = config.time_step;
    let liftoff_mass = vehicle.initial_mass();
    // `vehicle.validate()` above already guarantees every stage's
    // `isp_vac` and masses are physical, so this cannot actually error
    // here; propagate it anyway to keep the budget non-NaN by construction.
    let ideal_delta_v = vehicle.ideal_delta_v()?;

    // Launch on the +x axis; local "up" is +x, downrange ("east") is +y.
    let r0 = R_EARTH + config.launch_altitude_m;
    let pos0 = Vector2::new(r0, 0.0);
    // At rest on the ground means co-rotating with the atmosphere.
    let vel0 = dynamics::atmosphere_velocity(pos0);
    let launch_angle = pos0.y.atan2(pos0.x);

    let mut state = State {
        pos: pos0,
        vel: vel0,
    };
    let mut mass = liftoff_mass;
    let mut stage_idx = 0usize;
    let mut prop_remaining = vehicle.stages[0].propellant_mass;
    let mut t = 0.0f64;

    let mut events: Vec<FlightEvent> = Vec::new();
    let mut samples: Vec<TrajectorySample> = Vec::new();
    let mut max_q = 0.0f64;
    let mut max_q_alt = dynamics::altitude(pos0);
    let mut max_accel_g = 0.0f64;
    let mut last_sample_t = f64::NEG_INFINITY;

    events.push(FlightEvent {
        time: 0.0,
        altitude_m: dynamics::altitude(pos0),
        speed: vel0.norm(),
        kind: "Liftoff".into(),
    });

    let max_steps = (config.max_time / dt).ceil() as u64 + 2;
    let mut outcome = Outcome::TimedOut;
    let mut steps = 0u64;
    let mut controller = Controller::new(config.mode, config.guidance);

    loop {
        if steps > max_steps {
            return Err(AstroError::StepBudgetExceeded(max_steps));
        }
        steps += 1;

        let alt = dynamics::altitude(state.pos);
        let atmos = dynamics::atmosphere_at(state.pos);

        // Horizontal wind at this altitude, along the local east
        // (downrange) direction. Zero for the default still-air model.
        let east_hat = {
            let n = state.pos.norm();
            if n > 0.0 {
                Vector2::new(-state.pos.y / n, state.pos.x / n)
            } else {
                Vector2::new(0.0, 1.0)
            }
        };
        let wind_vec = config.wind.speed_at(alt) * east_hat;

        // Ask the flight controller what the engine should do this step.
        let cmd = controller.command(state.pos, state.vel);
        let burning = cmd.engine_on && stage_idx < vehicle.stages.len() && prop_remaining > 0.0;
        let (thrust_mag, mdot) = if burning {
            let st = &vehicle.stages[stage_idx];
            // Validated vehicle => mass_flow is finite & positive; `?`
            // never fires here but keeps the divisor non-NaN by contract.
            (st.thrust(atmos.pressure), st.mass_flow()?)
        } else {
            (0.0, 0.0)
        };
        let dir = controller.direction(cmd.policy, state.pos, state.vel, t);

        // Sensed (non-gravitational) acceleration, for the g-load peak.
        let drag_a = dynamics::drag_accel(
            state.pos,
            state.vel,
            mass,
            vehicle.reference_area,
            &vehicle.drag,
            &atmos,
            wind_vec,
        );
        let thrust_a = dynamics::thrust_accel(thrust_mag, dir, mass);
        let sensed_g = (drag_a + thrust_a).norm() / G0;
        max_accel_g = max_accel_g.max(sensed_g);

        let q = dynamics::dynamic_pressure(state.pos, state.vel, &atmos, wind_vec);
        if q > max_q {
            max_q = q;
            max_q_alt = alt;
        }

        // Record a down-sampled trajectory point.
        if samples.is_empty() || t - last_sample_t >= config.sample_interval - 1e-9 {
            let v_rel = state.vel - dynamics::atmosphere_velocity(state.pos);
            let speed_rel = v_rel.norm();
            let mach = if atmos.speed_of_sound > 0.0 {
                speed_rel / atmos.speed_of_sound
            } else {
                0.0
            };
            let angle = state.pos.y.atan2(state.pos.x) - launch_angle;
            samples.push(TrajectorySample {
                time: t,
                altitude_m: alt,
                downrange_m: R_EARTH * angle,
                speed_inertial: state.vel.norm(),
                speed_relative: speed_rel,
                mach,
                mass,
                dynamic_pressure: q,
                acceleration_g: sensed_g,
            });
            last_sample_t = t;
        }

        // Closed-loop insertion complete: the orbit is circularised.
        if cmd.terminate {
            events.push(FlightEvent {
                time: t,
                altitude_m: alt,
                speed: state.vel.norm(),
                kind: "Orbit insertion complete".into(),
            });
            outcome = Outcome::Meco;
            break;
        }

        // Surface impact (only fires once airborne work has begun).
        if alt < 0.0 {
            outcome = Outcome::Impact;
            events.push(FlightEvent {
                time: t,
                altitude_m: alt,
                speed: state.vel.norm(),
                kind: "Impact".into(),
            });
            break;
        }
        if t >= config.max_time {
            // `outcome` is already `TimedOut` from initialisation.
            break;
        }

        // RK4 step. Mass and burn state are frozen across the step; the
        // thrust direction + ambient thrust are re-evaluated at each
        // sub-state so a steep pressure/velocity gradient is captured.
        let policy = cmd.policy;
        let accel = |s: &State| -> Vector2<f64> {
            let atmos_s = dynamics::atmosphere_at(s.pos);
            let thrust_s = if burning {
                vehicle.stages[stage_idx].thrust(atmos_s.pressure)
            } else {
                0.0
            };
            let dir_s = controller.direction(policy, s.pos, s.vel, t);
            dynamics::total_accel(
                s.pos, s.vel, mass, vehicle, &atmos_s, thrust_s, dir_s, wind_vec,
            )
        };
        state = rk4_step(state, dt, &accel);

        // Discrete propellant burn + staging.
        if burning {
            let burned = mdot * dt;
            if burned >= prop_remaining {
                // Stage exhausted this step: consume the remaining
                // propellant, then either stage or shut down.
                mass -= prop_remaining;
                let dry = vehicle.stages[stage_idx].dry_mass;
                let alt_now = dynamics::altitude(state.pos);
                let spd_now = state.vel.norm();
                if stage_idx + 1 < vehicle.stages.len() {
                    mass -= dry; // jettison spent stage
                    events.push(FlightEvent {
                        time: t + dt,
                        altitude_m: alt_now,
                        speed: spd_now,
                        kind: format!("Staging: jettison {}", vehicle.stages[stage_idx].name),
                    });
                    stage_idx += 1;
                    prop_remaining = vehicle.stages[stage_idx].propellant_mass;
                } else {
                    events.push(FlightEvent {
                        time: t + dt,
                        altitude_m: alt_now,
                        speed: spd_now,
                        kind: "MECO (final stage burnout)".into(),
                    });
                    outcome = Outcome::Meco;
                    t += dt;
                    break;
                }
            } else {
                mass -= burned;
                prop_remaining -= burned;
            }
        }

        t += dt;
    }

    // The integrated ascent state is finite and well above the surface,
    // so use the validation-free core (the public `elements` would only
    // reject zero/non-finite/parabolic-singular states, none of which a
    // physical burned-out trajectory produces).
    let orbit = orbit::elements_with_mu_unchecked(state.pos, state.vel, crate::constants::MU_EARTH);
    let final_alt = dynamics::altitude(state.pos);
    let reached_space = orbit.apoapsis_altitude >= 100_000.0;
    let reached_orbit = orbit.is_bound && orbit.periapsis_altitude >= 100_000.0;

    Ok(AscentResult {
        outcome,
        liftoff_mass,
        ideal_delta_v,
        final_time: t,
        final_altitude_m: final_alt,
        final_speed_inertial: state.vel.norm(),
        orbit,
        max_dynamic_pressure: max_q,
        max_q_altitude_m: max_q_alt,
        max_acceleration_g: max_accel_g,
        reached_space,
        reached_orbit,
        events,
        samples,
        final_position_m: [state.pos.x, state.pos.y],
        final_velocity_ms: [state.vel.x, state.vel.y],
    })
}

/// One classical RK4 step of the second-order system
/// `ẋ = v, v̇ = a(x, v)`.
fn rk4_step<F>(s: State, dt: f64, accel: &F) -> State
where
    F: Fn(&State) -> Vector2<f64>,
{
    let k1x = s.vel;
    let k1v = accel(&s);

    let s2 = State {
        pos: s.pos + 0.5 * dt * k1x,
        vel: s.vel + 0.5 * dt * k1v,
    };
    let k2x = s2.vel;
    let k2v = accel(&s2);

    let s3 = State {
        pos: s.pos + 0.5 * dt * k2x,
        vel: s.vel + 0.5 * dt * k2v,
    };
    let k3x = s3.vel;
    let k3v = accel(&s3);

    let s4 = State {
        pos: s.pos + dt * k3x,
        vel: s.vel + dt * k3v,
    };
    let k4x = s4.vel;
    let k4v = accel(&s4);

    State {
        pos: s.pos + dt / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x),
        vel: s.vel + dt / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v),
    }
}

/// Propagate a ballistic (gravity-only, no thrust, no drag) two-body
/// state forward by `steps` RK4 steps of size `dt`. Useful for coasting
/// an insertion state and for validating the integrator.
///
/// # Errors
///
/// Returns [`AstroError::OutOfRange`] if `steps` exceeds
/// [`MAX_SIM_STEPS`], so an absurd step count is refused rather than run
/// to a hang.
pub fn propagate_two_body(
    position: Vector2<f64>,
    velocity: Vector2<f64>,
    dt: f64,
    steps: u64,
) -> Result<(Vector2<f64>, Vector2<f64>), AstroError> {
    check_step_count(steps)?;
    let mut s = State {
        pos: position,
        vel: velocity,
    };
    let accel = |st: &State| dynamics::gravity_accel(st.pos);
    for _ in 0..steps {
        s = rk4_step(s, dt, &accel);
    }
    Ok((s.pos, s.vel))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::MU_EARTH;

    #[test]
    fn rk4_conserves_energy_on_circular_orbit() {
        // A circular LEO propagated for one full period should return
        // close to its start, and specific energy must be conserved.
        let radius = R_EARTH + 400_000.0;
        let v = (MU_EARTH / radius).sqrt();
        let pos = Vector2::new(radius, 0.0);
        let vel = Vector2::new(0.0, v);
        let energy0 = vel.norm_squared() / 2.0 - MU_EARTH / pos.norm();

        let period = 2.0 * std::f64::consts::PI * (radius.powi(3) / MU_EARTH).sqrt();
        let dt = 1.0;
        let steps = (period / dt).round() as u64;
        let (pf, vf) = propagate_two_body(pos, vel, dt, steps).expect("valid step count");

        let energy1 = vf.norm_squared() / 2.0 - MU_EARTH / pf.norm();
        assert!(
            (energy1 - energy0).abs() / energy0.abs() < 1e-6,
            "energy drift: {energy0} -> {energy1}"
        );
        // Position should return to roughly the start (< 5 km of a
        // ~6778 km radius after a full revolution).
        assert!(
            (pf - pos).norm() < 5_000.0,
            "closure error {}",
            (pf - pos).norm()
        );
    }

    #[test]
    fn rk4_matches_analytical_kepler_on_eccentric_orbit() {
        // GROUND-TRUTH VALIDATION. Propagate an eccentric (e = 0.6) orbit
        // with the RK4 two-body integrator and compare, at several orbital
        // phases, against the EXACT analytical Kepler solution: solve
        // `M = E - e·sinE` for the eccentric anomaly, then the closed-form
        // perifocal position `(a(cosE - e), a√(1-e²)·sinE)`. The gap is the
        // integrator's true error against the analytic answer — the real
        // measure of "is the orbital solver correct?", far stronger than a
        // circular orbit (constant speed) or mere energy conservation.
        let mu = MU_EARTH;
        let a = R_EARTH + 20_000_000.0; // ~26 378 km semi-major axis (GTO-like)
        let e = 0.6;

        // Initial state at periapsis (on +x, moving +y → prograde).
        let r_p = a * (1.0 - e);
        let v_p = (mu / a).sqrt() * ((1.0 + e) / (1.0 - e)).sqrt();
        let pos0 = Vector2::new(r_p, 0.0);
        let vel0 = Vector2::new(0.0, v_p);

        let n = (mu / a.powi(3)).sqrt(); // mean motion
        let period = 2.0 * std::f64::consts::PI / n;

        // Exact analytic position at elapsed time `t` (periapsis at t = 0).
        let kepler = |t: f64| -> Vector2<f64> {
            let m = n * t;
            let mut ea = m; // Newton on Kepler's equation
            for _ in 0..60 {
                ea -= (ea - e * ea.sin() - m) / (1.0 - e * ea.cos());
            }
            Vector2::new(a * (ea.cos() - e), a * (1.0 - e * e).sqrt() * ea.sin())
        };

        let dt = 1.0;
        let mut max_err = 0.0_f64;
        for frac in [0.1_f64, 0.25, 0.5, 0.75, 0.9, 1.0] {
            let steps = (frac * period / dt).round() as u64;
            let (pf, _vf) = propagate_two_body(pos0, vel0, dt, steps).expect("valid steps");
            let err = (pf - kepler(steps as f64 * dt)).norm();
            max_err = max_err.max(err);
        }
        println!(
            "VALIDATION rk4-vs-Kepler: a={:.0} km e={e} period={:.0} s  max position error = {max_err:.4e} m",
            a / 1000.0,
            period
        );
        // Measured: ~1.7e-5 m (≈17 µm, relative ≈6e-13) — floating-point-
        // limited. RK4 at a 1 s step tracks the exact two-body solution to
        // sub-millimetre over a full ~11.8 h eccentric revolution. The 1 m
        // bound is a generous regression guard against future degradation.
        assert!(
            max_err < 1.0,
            "RK4 deviates from analytic Kepler by {max_err:.3e} m (> 1 m)"
        );
    }

    #[test]
    fn simulate_ascent_rejects_unbounded_config_without_hanging() {
        // The H1 hang/OOM repro: a 1 ns step with a 1e15 s cap would be
        // ~u64::MAX steps and an unbounded sample Vec. It must return an
        // Err immediately, not loop.
        let vehicle = crate::presets::two_stage_medium_lift();
        let config = AscentConfig {
            time_step: 1e-9,
            max_time: 1e15,
            ..crate::presets::leo_ascent_config()
        };
        assert!(
            simulate_ascent(&vehicle, &config).is_err(),
            "unbounded config must be rejected, not run"
        );
    }

    #[test]
    fn radius_stays_constant_on_circular_orbit() {
        let radius = R_EARTH + 500_000.0;
        let v = (MU_EARTH / radius).sqrt();
        let mut pos = Vector2::new(radius, 0.0);
        let mut vel = Vector2::new(0.0, v);
        for _ in 0..600 {
            let (p, w) = propagate_two_body(pos, vel, 1.0, 1).expect("valid step count");
            pos = p;
            vel = w;
            assert!((pos.norm() - radius).abs() < 2_000.0);
        }
    }

    #[test]
    fn propagate_two_body_rejects_absurd_step_count() {
        // u64::MAX steps would hang; must be a clean Err instead.
        let pos = Vector2::new(R_EARTH + 400_000.0, 0.0);
        let vel = Vector2::new(0.0, 7_700.0);
        let r = propagate_two_body(pos, vel, 1.0, u64::MAX);
        assert!(
            matches!(r, Err(AstroError::OutOfRange { what: "steps", .. })),
            "u64::MAX steps must be rejected, got {r:?}"
        );
        // A valid step count still works.
        assert!(propagate_two_body(pos, vel, 1.0, 10).is_ok());
    }
}
