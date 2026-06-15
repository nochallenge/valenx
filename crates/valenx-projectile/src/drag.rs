//! Projectile motion **with quadratic aerodynamic drag**, integrated
//! numerically by classical fourth-order Runge–Kutta (RK4).
//!
//! ## Model
//!
//! A point mass moves in the vertical plane under uniform gravity `g`
//! and a drag force opposing its velocity with magnitude proportional to
//! the *square* of the speed (the high-Reynolds-number regime that holds
//! for everyday projectiles: balls, shells, stones):
//!
//! ```text
//! F_drag = -½ ρ C_d A |v| v        (vector, opposes motion)
//! ```
//!
//! Dividing by the mass `m` and folding the constants into a single
//! lumped **drag coefficient**
//!
//! ```text
//! k = ρ C_d A / (2 m)              (units: 1/m)
//! ```
//!
//! the equations of motion become the first-order system
//!
//! ```text
//! dx/dt = vx
//! dy/dt = vy
//! dvx/dt = -k·|v|·vx
//! dvy/dt = -g - k·|v|·vy        with |v| = sqrt(vx² + vy²)
//! ```
//!
//! Setting `k = 0` recovers exact vacuum motion, which the tests use as
//! a cross-check against the closed forms in [`crate::vacuum`].
//!
//! ## Honest scope
//!
//! `g`, `ρ`, `C_d` and `A` are all treated as constant over the flight
//! (no altitude-varying density, no Mach- or Reynolds-dependent drag
//! coefficient, no wind, no Magnus/spin, no Coriolis). This is the
//! standard textbook quadratic-drag trajectory, suitable for teaching
//! and rough estimates — not a ballistic-table or fire-control tool.

use serde::{Deserialize, Serialize};

use crate::error::{ProjectileError, Result};
use crate::vacuum::STANDARD_GRAVITY;

/// A launch condition for a drag-affected trajectory on flat, level
/// ground. SI units throughout: speed m/s, angle radians, gravity m/s²,
/// drag coefficient `k` in 1/m (see the module docs for its definition).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DragShot {
    /// Launch speed `v0` (m/s), strictly positive.
    pub speed: f64,
    /// Launch elevation angle `θ` (radians) in `[0, π/2]`.
    pub angle_rad: f64,
    /// Gravitational acceleration `g` (m/s²), strictly positive.
    pub gravity: f64,
    /// Lumped quadratic-drag coefficient `k = ρ·C_d·A/(2m)` (1/m),
    /// non-negative. `k = 0` is the vacuum limit.
    pub drag_coeff: f64,
}

/// The integrated outcome of a drag trajectory: where and when it came
/// back down to the launch height, and how high it climbed.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Trajectory {
    /// Horizontal distance (m) at which `y` returned to the launch
    /// height (`y = 0`), found by linear interpolation across the step
    /// in which the sign of `y` changed.
    pub range: f64,
    /// Total time of flight (s) to that range.
    pub time_of_flight: f64,
    /// Maximum height (m) reached during the flight (apex).
    pub apex_height: f64,
    /// Time (s) at which the apex was reached.
    pub apex_time: f64,
    /// Number of RK4 steps actually integrated.
    pub steps: u64,
}

/// One phase-space state `[x, y, vx, vy]` used by the integrator.
type State = [f64; 4];

impl DragShot {
    /// Build a validated drag shot.
    ///
    /// # Errors
    ///
    /// Returns [`ProjectileError::InvalidParameter`] if `speed` or
    /// `gravity` is not finite-and-positive or `drag_coeff` is not
    /// finite-and-non-negative, or [`ProjectileError::AngleOutOfRange`]
    /// if `angle_rad` is outside `[0, π/2]`.
    pub fn new(speed: f64, angle_rad: f64, gravity: f64, drag_coeff: f64) -> Result<Self> {
        let speed = ProjectileError::require_positive("speed", speed)?;
        let gravity = ProjectileError::require_positive("gravity", gravity)?;
        let drag_coeff = ProjectileError::require_non_negative("drag_coeff", drag_coeff)?;
        let angle_rad = ProjectileError::require_angle(angle_rad)?;
        Ok(Self {
            speed,
            angle_rad,
            gravity,
            drag_coeff,
        })
    }

    /// Build a drag shot using [`STANDARD_GRAVITY`] for `g`.
    ///
    /// # Errors
    ///
    /// As [`DragShot::new`], minus the gravity check.
    pub fn at_standard_gravity(speed: f64, angle_rad: f64, drag_coeff: f64) -> Result<Self> {
        Self::new(speed, angle_rad, STANDARD_GRAVITY, drag_coeff)
    }

    /// The time-derivative of the state `[x, y, vx, vy]` under gravity +
    /// quadratic drag.
    #[inline]
    fn deriv(&self, s: &State) -> State {
        let (vx, vy) = (s[2], s[3]);
        let speed = (vx * vx + vy * vy).sqrt();
        let k = self.drag_coeff;
        [vx, vy, -k * speed * vx, -self.gravity - k * speed * vy]
    }

    /// Integrate the trajectory with a fixed RK4 step `dt` until the
    /// projectile falls back to the launch height (`y ≤ 0` after having
    /// risen), then linearly interpolate the landing point.
    ///
    /// `dt` is the integration step in seconds; `max_steps` caps the run
    /// so a pathological input cannot loop forever.
    ///
    /// # Errors
    ///
    /// - [`ProjectileError::InvalidParameter`] if `dt` is not
    ///   finite-and-positive.
    /// - [`ProjectileError::StepBudgetExceeded`] if the projectile has
    ///   not landed within `max_steps` steps.
    pub fn integrate(&self, dt: f64, max_steps: u64) -> Result<Trajectory> {
        let dt = ProjectileError::require_positive("dt", dt)?;

        // Initial state at the origin with the launch velocity.
        let (sin, cos) = self.angle_rad.sin_cos();
        let mut s: State = [0.0, 0.0, self.speed * cos, self.speed * sin];
        let mut t = 0.0_f64;

        let mut apex_height = 0.0_f64;
        let mut apex_time = 0.0_f64;

        // A purely horizontal (θ = 0) launch never rises, so it "lands"
        // immediately: report a zero-length flight rather than looping.
        if sin == 0.0 {
            return Ok(Trajectory {
                range: 0.0,
                time_of_flight: 0.0,
                apex_height: 0.0,
                apex_time: 0.0,
                steps: 0,
            });
        }

        for step in 1..=max_steps {
            let prev = s;
            let prev_t = t;

            // Classical RK4 over one step of size dt.
            let k1 = self.deriv(&s);
            let s2 = add_scaled(&s, &k1, dt * 0.5);
            let k2 = self.deriv(&s2);
            let s3 = add_scaled(&s, &k2, dt * 0.5);
            let k3 = self.deriv(&s3);
            let s4 = add_scaled(&s, &k3, dt);
            let k4 = self.deriv(&s4);
            for i in 0..4 {
                s[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
            t += dt;

            // Track the apex (peak y). With monotone-up-then-down motion
            // this captures the true maximum to step resolution.
            if s[1] > apex_height {
                apex_height = s[1];
                apex_time = t;
            }

            // Landing: y crossed from >= 0 to < 0. Interpolate linearly
            // in y across the step to recover the ground intercept.
            if s[1] < 0.0 && prev[1] >= 0.0 {
                let y0 = prev[1];
                let y1 = s[1];
                let frac = y0 / (y0 - y1); // in [0, 1)
                let range = prev[0] + frac * (s[0] - prev[0]);
                let tof = prev_t + frac * (t - prev_t);
                return Ok(Trajectory {
                    range,
                    time_of_flight: tof,
                    apex_height,
                    apex_time,
                    steps: step,
                });
            }
        }

        Err(ProjectileError::StepBudgetExceeded(max_steps))
    }
}

/// `out = a + scale·b`, componentwise over a 4-vector.
#[inline]
fn add_scaled(a: &State, b: &State, scale: f64) -> State {
    [
        a[0] + scale * b[0],
        a[1] + scale * b[1],
        a[2] + scale * b[2],
        a[3] + scale * b[3],
    ]
}

/// Settings shared by every trajectory evaluated during an
/// optimal-angle search.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct OptimizeConfig {
    /// Launch speed `v0` (m/s).
    pub speed: f64,
    /// Gravity `g` (m/s²).
    pub gravity: f64,
    /// Lumped drag coefficient `k` (1/m).
    pub drag_coeff: f64,
    /// RK4 step (s) used for each candidate trajectory.
    pub dt: f64,
    /// Per-trajectory step budget.
    pub max_steps: u64,
}

impl OptimizeConfig {
    /// Build a validated optimisation configuration.
    ///
    /// # Errors
    ///
    /// [`ProjectileError::InvalidParameter`] for any non-physical scalar
    /// (`speed`/`gravity`/`dt` not positive, `drag_coeff` negative).
    pub fn new(speed: f64, gravity: f64, drag_coeff: f64, dt: f64, max_steps: u64) -> Result<Self> {
        let speed = ProjectileError::require_positive("speed", speed)?;
        let gravity = ProjectileError::require_positive("gravity", gravity)?;
        let drag_coeff = ProjectileError::require_non_negative("drag_coeff", drag_coeff)?;
        let dt = ProjectileError::require_positive("dt", dt)?;
        Ok(Self {
            speed,
            gravity,
            drag_coeff,
            dt,
            max_steps,
        })
    }

    /// Evaluate the level-ground range for a single launch angle.
    fn range_at(&self, angle_rad: f64) -> Result<f64> {
        let shot = DragShot::new(self.speed, angle_rad, self.gravity, self.drag_coeff)?;
        Ok(shot.integrate(self.dt, self.max_steps)?.range)
    }
}

/// The angle (radians) that maximises range under drag, plus the range
/// achieved there.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct OptimalAngle {
    /// Range-maximising launch angle (radians) within the search bracket.
    pub angle_rad: f64,
    /// Range (m) achieved at [`OptimalAngle::angle_rad`].
    pub range: f64,
}

impl OptimalAngle {
    /// The optimal launch angle in degrees.
    #[must_use]
    pub fn angle_deg(&self) -> f64 {
        self.angle_rad.to_degrees()
    }
}

/// Find the range-maximising launch angle by **golden-section search**
/// over `[lo_rad, hi_rad]`.
///
/// The drag-trajectory range is a smooth, single-peaked (unimodal)
/// function of launch angle on `(0, π/2)`, so golden-section search — a
/// bracket-shrinking method that needs only function evaluations, no
/// derivatives — converges robustly to the interior maximum. The bracket
/// width is reduced by the golden ratio each iteration until it is below
/// `tol_rad` (or `max_iter` is hit).
///
/// # Errors
///
/// - [`ProjectileError::InvalidBracket`] if the bracket is degenerate or
///   leaves `[0, π/2]`, or if `tol_rad` is not positive.
/// - Any error propagated from evaluating a candidate trajectory.
pub fn optimal_drag_angle(
    cfg: &OptimizeConfig,
    lo_rad: f64,
    hi_rad: f64,
    tol_rad: f64,
    max_iter: u32,
) -> Result<OptimalAngle> {
    if !lo_rad.is_finite() || !hi_rad.is_finite() || lo_rad >= hi_rad {
        return Err(ProjectileError::InvalidBracket(
            "lower bound must be finite and strictly below the finite upper bound",
        ));
    }
    if lo_rad < 0.0 || hi_rad > core::f64::consts::FRAC_PI_2 {
        return Err(ProjectileError::InvalidBracket(
            "bracket must lie within [0, π/2]",
        ));
    }
    if !tol_rad.is_finite() || tol_rad <= 0.0 {
        return Err(ProjectileError::InvalidBracket(
            "tolerance must be finite and > 0",
        ));
    }

    // Golden ratio constants for the interior probe points.
    let inv_phi: f64 = (5.0_f64.sqrt() - 1.0) / 2.0; // 1/φ  ≈ 0.6180339887
    let inv_phi2: f64 = (3.0 - 5.0_f64.sqrt()) / 2.0; // 1/φ² ≈ 0.3819660113

    let mut a = lo_rad;
    let mut b = hi_rad;
    let mut c = a + inv_phi2 * (b - a);
    let mut d = a + inv_phi * (b - a);
    let mut fc = cfg.range_at(c)?;
    let mut fd = cfg.range_at(d)?;

    let mut iters = 0_u32;
    while (b - a) > tol_rad && iters < max_iter {
        if fc > fd {
            // Maximum is in [a, d]; discard the upper part.
            b = d;
            d = c;
            fd = fc;
            c = a + inv_phi2 * (b - a);
            fc = cfg.range_at(c)?;
        } else {
            // Maximum is in [c, b]; discard the lower part.
            a = c;
            c = d;
            fc = fd;
            d = a + inv_phi * (b - a);
            fd = cfg.range_at(d)?;
        }
        iters += 1;
    }

    let angle_rad = 0.5 * (a + b);
    let range = cfg.range_at(angle_rad)?;
    Ok(OptimalAngle { angle_rad, range })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vacuum::VacuumShot;
    use core::f64::consts::{FRAC_PI_2, FRAC_PI_4};

    /// With zero drag the RK4 integrator must reproduce the closed-form
    /// vacuum range/apex/time to a few parts in ten thousand — this is
    /// the integrator's ground-truth correctness check.
    #[test]
    fn zero_drag_matches_vacuum_closed_form() {
        let v0 = 40.0;
        let g = 9.81;
        for deg in [15.0_f64, 30.0, 45.0, 60.0, 75.0] {
            let theta = deg.to_radians();
            let exact = VacuumShot::new(v0, theta, g).unwrap();
            let num = DragShot::new(v0, theta, g, 0.0)
                .unwrap()
                .integrate(1e-3, 5_000_000)
                .unwrap();

            let rel_range = (num.range - exact.range()).abs() / exact.range();
            assert!(rel_range < 1e-3, "range rel err {rel_range} at {deg}°");

            let rel_tof =
                (num.time_of_flight - exact.time_of_flight()).abs() / exact.time_of_flight();
            assert!(rel_tof < 1e-3, "tof rel err {rel_tof} at {deg}°");

            // Apex is captured to step resolution (vy·dt ~ a few cm here).
            let abs_apex = (num.apex_height - exact.apex_height()).abs();
            assert!(abs_apex < 0.5, "apex abs err {abs_apex} m at {deg}°");
        }
    }

    /// Drag must shorten the flight: at a fixed angle the drag range is
    /// strictly less than the vacuum range.
    #[test]
    fn drag_range_is_less_than_vacuum() {
        let v0 = 60.0;
        let g = 9.81;
        let theta = FRAC_PI_4;
        let vac = VacuumShot::new(v0, theta, g).unwrap().range();
        let drag = DragShot::new(v0, theta, g, 0.01)
            .unwrap()
            .integrate(1e-3, 5_000_000)
            .unwrap()
            .range;
        assert!(
            drag < vac,
            "expected drag range {drag} < vacuum range {vac}"
        );
        assert!(drag > 0.0);
    }

    /// Heavier drag shortens the range monotonically.
    #[test]
    fn more_drag_means_less_range() {
        let v0 = 50.0;
        let g = 9.81;
        let theta = FRAC_PI_4;
        let mut last = f64::INFINITY;
        for k in [0.0, 0.005, 0.01, 0.02, 0.04] {
            let r = DragShot::new(v0, theta, g, k)
                .unwrap()
                .integrate(1e-3, 5_000_000)
                .unwrap()
                .range;
            assert!(r < last, "range did not decrease at k = {k}");
            last = r;
        }
    }

    /// The headline validation: under quadratic drag the optimal launch
    /// angle is strictly **below** 45°.
    #[test]
    fn optimal_drag_angle_is_below_45_degrees() {
        let cfg = OptimizeConfig::new(80.0, 9.81, 0.01, 1e-3, 5_000_000).unwrap();
        let opt = optimal_drag_angle(&cfg, 1.0_f64.to_radians(), 89.0_f64.to_radians(), 1e-5, 200)
            .unwrap();
        assert!(
            opt.angle_rad < FRAC_PI_4,
            "drag optimum {} rad ({}°) was not below 45°",
            opt.angle_rad,
            opt.angle_deg()
        );
        // Sanity: a sensible launch-angle band, not pinned to a bound.
        assert!(opt.angle_deg() > 20.0 && opt.angle_deg() < 44.9);
    }

    /// The golden-section optimum is genuinely the best: it beats a dense
    /// brute-force angle sweep (to within the search tolerance).
    #[test]
    fn golden_section_beats_brute_force_sweep() {
        let cfg = OptimizeConfig::new(70.0, 9.81, 0.012, 1e-3, 5_000_000).unwrap();
        let opt = optimal_drag_angle(&cfg, 1.0_f64.to_radians(), 89.0_f64.to_radians(), 1e-5, 200)
            .unwrap();

        // Brute-force 1-degree sweep.
        let mut best = 0.0_f64;
        for deg in 1..=89 {
            let r = cfg.range_at((deg as f64).to_radians()).unwrap();
            best = best.max(r);
        }
        // The golden-section result is at least as good as the best grid
        // point, minus a small slack for the grid being coarse.
        assert!(
            opt.range >= best - 1e-3,
            "golden {} < brute {}",
            opt.range,
            best
        );
    }

    /// In the vacuum limit (`k = 0`) the optimiser must recover 45°,
    /// tying the numerical search back to the analytic result.
    #[test]
    fn optimal_angle_recovers_45_when_drag_free() {
        let cfg = OptimizeConfig::new(45.0, 9.81, 0.0, 1e-3, 5_000_000).unwrap();
        let opt = optimal_drag_angle(&cfg, 1.0_f64.to_radians(), 89.0_f64.to_radians(), 1e-6, 300)
            .unwrap();
        assert!(
            (opt.angle_rad - FRAC_PI_4).abs() < 1e-3,
            "drag-free optimum {} rad not ~45°",
            opt.angle_rad
        );
    }

    /// A horizontal launch returns a degenerate zero-flight trajectory
    /// rather than spinning to the step budget.
    #[test]
    fn horizontal_launch_is_degenerate() {
        let traj = DragShot::new(30.0, 0.0, 9.81, 0.01)
            .unwrap()
            .integrate(1e-3, 1_000)
            .unwrap();
        assert!(traj.range.abs() < 1e-12);
        assert_eq!(traj.steps, 0);
    }

    /// Too small a step budget surfaces as a clean error, not a hang.
    #[test]
    fn exhausting_step_budget_errors() {
        let err = DragShot::new(50.0, FRAC_PI_4, 9.81, 0.0)
            .unwrap()
            .integrate(1e-3, 5) // nowhere near enough to land
            .unwrap_err();
        assert!(matches!(err, ProjectileError::StepBudgetExceeded(5)));
    }

    /// Bad optimisation brackets are rejected.
    #[test]
    fn rejects_bad_bracket() {
        let cfg = OptimizeConfig::new(50.0, 9.81, 0.01, 1e-3, 100_000).unwrap();
        assert!(matches!(
            optimal_drag_angle(&cfg, 1.0, 1.0, 1e-5, 100),
            Err(ProjectileError::InvalidBracket(_))
        ));
        assert!(matches!(
            optimal_drag_angle(&cfg, -0.1, 1.0, 1e-5, 100),
            Err(ProjectileError::InvalidBracket(_))
        ));
        assert!(matches!(
            optimal_drag_angle(&cfg, 0.1, FRAC_PI_2 + 0.1, 1e-5, 100),
            Err(ProjectileError::InvalidBracket(_))
        ));
        assert!(matches!(
            optimal_drag_angle(&cfg, 0.1, 1.0, 0.0, 100),
            Err(ProjectileError::InvalidBracket(_))
        ));
    }

    /// Invalid scalars are rejected by the config constructor.
    #[test]
    fn rejects_invalid_config() {
        assert!(matches!(
            OptimizeConfig::new(0.0, 9.81, 0.01, 1e-3, 100),
            Err(ProjectileError::InvalidParameter { field: "speed", .. })
        ));
        assert!(matches!(
            OptimizeConfig::new(10.0, 9.81, -1.0, 1e-3, 100),
            Err(ProjectileError::InvalidParameter {
                field: "drag_coeff",
                ..
            })
        ));
    }
}
