//! # valenx-projectile — point-mass projectile ballistics
//!
//! Elementary exterior ballistics for a point mass launched over flat,
//! level ground: the exact drag-free (vacuum) trajectory in closed form,
//! and a drag-affected trajectory integrated numerically. Pure
//! algorithms, no external processes, no platform dependencies.
//!
//! ## What
//!
//! - **Vacuum motion** ([`vacuum`]) — a [`VacuumShot`] yields the
//!   closed-form time of flight, horizontal [`range`](VacuumShot::range),
//!   apex height and time-to-apex. Helpers
//!   [`optimal_vacuum_angle_rad`] and [`max_vacuum_range`] give the
//!   range-maximising 45° launch and the resulting `v0²/g` peak range,
//!   and [`vacuum_angles_for_range`] inverts the range relation to the
//!   two (complementary) launch angles that reach a target range.
//! - **Quadratic-drag motion** ([`drag`]) — a [`DragShot`] integrates the
//!   gravity-plus-quadratic-drag equations of motion with classical
//!   fourth-order Runge–Kutta ([`DragShot::integrate`]) and reports the
//!   landing [`Trajectory`] (range, flight time, apex). The
//!   range-maximising launch angle under drag is found by golden-section
//!   search ([`optimal_drag_angle`]); because drag penalises the longer,
//!   higher 45° arc, that optimum falls **below** 45°.
//!
//! ## Model
//!
//! Flat-Earth, uniform gravity `g`, single point mass, launch and impact
//! at the same height. The drag force is the high-Reynolds-number
//! quadratic law `F = -½ρC_dA|v|v`, whose constants are folded into one
//! lumped coefficient `k = ρC_dA/(2m)` (units 1/m). Setting `k = 0`
//! recovers exact vacuum motion, which the test-suite uses to validate
//! the integrator against the closed forms.
//!
//! Ground-truth relations the tests pin down:
//!
//! - `R = v0²·sin(2θ)/g`, maximised at `θ = 45°` with `R_max = v0²/g`;
//!   inverting it, a sub-maximal target range is reached by two
//!   complementary angles `½·arcsin(Rg/v0²)` and its complement.
//! - `H = v0²·sin²(θ)/(2g)` for the apex.
//! - With drag, `R < R_vacuum` at the same angle, and the optimal angle
//!   is `< 45°`.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are textbook closed-form
//! point-mass equations and a well-established explicit numerical
//! integrator (RK4), **not** a clinical/medical tool and **not** a
//! production engineering, weapons, fire-control, or
//! ballistic-certification tool. The model holds `g`, air density,
//! drag coefficient and reference area constant over the flight and
//! ignores altitude-varying density, Mach/Reynolds-dependent drag, wind,
//! spin/Magnus lift, the Coriolis effect, and any terrain other than a
//! level plane. Use it for teaching and order-of-magnitude estimates,
//! not for anything where a real trajectory has to be relied upon.
//!
//! ## Errors
//!
//! Every fallible entry point returns
//! [`Result<_, ProjectileError>`](error::ProjectileError); inputs that
//! would otherwise produce a silent `NaN`/`Inf` (non-physical speed,
//! gravity, drag coefficient, or launch angle) are rejected up front by
//! the validated constructors.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod drag;
pub mod error;
pub mod vacuum;

// --- Convenience re-exports of the most-used items --------------------

pub use drag::{optimal_drag_angle, DragShot, OptimalAngle, OptimizeConfig, Trajectory};
pub use error::{ProjectileError, Result};
pub use vacuum::{
    max_vacuum_range, optimal_vacuum_angle_rad, vacuum_angles_for_range, VacuumShot,
    STANDARD_GRAVITY,
};

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::FRAC_PI_4;

    /// End-to-end: the same launch evaluated in vacuum and with a small
    /// drag must agree on the qualitative facts — vacuum range peaks at
    /// 45°, drag shortens it, and the drag-optimal angle is below 45°.
    #[test]
    fn end_to_end_vacuum_vs_drag() {
        let v0 = 75.0;
        let g = STANDARD_GRAVITY;

        // Vacuum: 45° is optimal and gives v0²/g.
        let vac45 = VacuumShot::new(v0, optimal_vacuum_angle_rad(), g)
            .unwrap()
            .range();
        assert!((vac45 - max_vacuum_range(v0, g).unwrap()).abs() < 1e-6);

        // Drag at 45° is shorter than vacuum at 45°.
        let drag45 = DragShot::new(v0, FRAC_PI_4, g, 0.008)
            .unwrap()
            .integrate(1e-3, 5_000_000)
            .unwrap()
            .range;
        assert!(drag45 < vac45);

        // The drag-optimal angle is below 45° and does at least as well
        // as the (suboptimal) 45° drag shot.
        let cfg = OptimizeConfig::new(v0, g, 0.008, 1e-3, 5_000_000).unwrap();
        let opt = optimal_drag_angle(&cfg, 1.0_f64.to_radians(), 89.0_f64.to_radians(), 1e-5, 200)
            .unwrap();
        assert!(opt.angle_rad < FRAC_PI_4);
        assert!(opt.range >= drag45 - 1e-6);
    }

    /// The public config / result structs round-trip through JSON.
    #[test]
    fn structs_serialize_round_trip() {
        let shot = DragShot::at_standard_gravity(30.0, FRAC_PI_4, 0.01).unwrap();
        let json = serde_json::to_string(&shot).unwrap();
        let back: DragShot = serde_json::from_str(&json).unwrap();
        assert_eq!(shot, back);

        let traj = shot.integrate(1e-3, 5_000_000).unwrap();
        let tj = serde_json::to_string(&traj).unwrap();
        let traj_back: Trajectory = serde_json::from_str(&tj).unwrap();
        assert_eq!(traj, traj_back);
    }
}
