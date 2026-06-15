//! Fixed-step classical fourth-order Runge-Kutta (RK4) integrator.
//!
//! ## Model
//!
//! Every dynamical system in this crate is an autonomous first-order
//! ODE on a fixed-length state vector,
//!
//! ```text
//! dy/dt = f(t, y),   y(t0) = y0,   y in R^D.
//! ```
//!
//! The classical RK4 step advances `y_n -> y_{n+1}` over a step `h`
//! with four slope evaluations:
//!
//! ```text
//! k1 = f(t,        y)
//! k2 = f(t + h/2,  y + (h/2) k1)
//! k3 = f(t + h/2,  y + (h/2) k2)
//! k4 = f(t + h,    y + h     k3)
//! y_{n+1} = y_n + (h/6)(k1 + 2 k2 + 2 k3 + k4).
//! ```
//!
//! This is a fourth-order method: the local truncation error per step
//! is `O(h^5)` and the global error over a fixed interval is `O(h^4)`,
//! which the [`tests`](self) verify against the analytic solution of
//! `y' = y` (whose error must fall by ~16x when `h` is halved).
//!
//! The integrator is fixed-step and explicit — appropriate for the
//! smooth, non-stiff population and epidemic models here. It is *not*
//! adaptive and carries no stiff-system safeguards.

use crate::error::{PopError, Result};

/// A fixed-length state vector of `D` real components.
///
/// Population models map their compartments onto the slots: logistic
/// uses `D = 1` (`[N]`), SIR uses `D = 3` (`[S, I, R]`), Lotka-Volterra
/// uses `D = 2` (`[prey, predator]`).
pub type State<const D: usize> = [f64; D];

/// One sample of an integrated trajectory: the time and the state at
/// that time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sample<const D: usize> {
    /// Time of this sample.
    pub t: f64,
    /// State vector at time [`t`](Self::t).
    pub y: State<D>,
}

/// Hard ceiling on the number of RK4 steps a single
/// [`integrate`] call may take. Guards against a tiny `dt` over a huge
/// horizon silently allocating an enormous trajectory.
pub const MAX_STEPS: u64 = 100_000_000;

/// Advance one classical RK4 step of size `h` from `(t, y)` under the
/// right-hand side `f`.
///
/// `f(t, y)` returns `dy/dt`. The returned value is `y(t + h)` to
/// fourth order. This is the building block of [`integrate`]; it is
/// exposed for callers that want to drive their own stepping loop.
///
/// # Example
///
/// ```
/// use valenx_popdynamics::rk4::rk4_step;
/// // Exponential growth y' = y, exact y(h) = y0 * e^h.
/// let y = rk4_step(0.0, [1.0], 0.1, |_t, s| [s[0]]);
/// let exact = 0.1_f64.exp();
/// assert!((y[0] - exact).abs() < 1e-6);
/// ```
pub fn rk4_step<const D: usize, F>(t: f64, y: State<D>, h: f64, f: F) -> State<D>
where
    F: Fn(f64, &State<D>) -> State<D>,
{
    let k1 = f(t, &y);

    let mut y2 = y;
    for i in 0..D {
        y2[i] = y[i] + 0.5 * h * k1[i];
    }
    let k2 = f(t + 0.5 * h, &y2);

    let mut y3 = y;
    for i in 0..D {
        y3[i] = y[i] + 0.5 * h * k2[i];
    }
    let k3 = f(t + 0.5 * h, &y3);

    let mut y4 = y;
    for i in 0..D {
        y4[i] = y[i] + h * k3[i];
    }
    let k4 = f(t + h, &y4);

    let mut out = y;
    for i in 0..D {
        out[i] = y[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    out
}

/// Integrate `dy/dt = f(t, y)` from `t_start` to `t_end` with fixed
/// step `dt`, returning the trajectory including both endpoints.
///
/// The number of whole steps is `n = floor((t_end - t_start) / dt)`;
/// the returned vector has `n + 1` samples, the first at `t_start` and
/// the last at `t_start + n * dt`. If `(t_end - t_start)` is not an
/// exact multiple of `dt`, the final fractional remainder is *not*
/// integrated — pick a `dt` that divides the window for an exact
/// landing on `t_end`.
///
/// # Errors
///
/// - [`PopError::Invalid`] if `dt <= 0` or `t_end <= t_start`.
/// - [`PopError::TooManySteps`] if the implied step count exceeds
///   [`MAX_STEPS`].
pub fn integrate<const D: usize, F>(
    f: F,
    y0: State<D>,
    t_start: f64,
    t_end: f64,
    dt: f64,
) -> Result<Vec<Sample<D>>>
where
    F: Fn(f64, &State<D>) -> State<D>,
{
    if dt.is_nan() || dt <= 0.0 {
        return Err(PopError::invalid("dt", "step size must be positive"));
    }
    if t_start.is_nan() || t_end.is_nan() || t_end <= t_start {
        return Err(PopError::invalid(
            "t_end",
            "end time must be strictly greater than start time",
        ));
    }

    let span = t_end - t_start;
    // floor of span/dt, guarded so the cast cannot wrap.
    let n_float = (span / dt).floor();
    if !n_float.is_finite() || n_float < 0.0 || n_float > MAX_STEPS as f64 {
        return Err(PopError::TooManySteps {
            requested: if n_float.is_finite() && n_float >= 0.0 {
                n_float as u64
            } else {
                u64::MAX
            },
            ceiling: MAX_STEPS,
        });
    }
    let n = n_float as u64;

    let mut out = Vec::with_capacity(n as usize + 1);
    let mut t = t_start;
    let mut y = y0;
    out.push(Sample { t, y });

    for _ in 0..n {
        y = rk4_step(t, y, dt, &f);
        t += dt;
        out.push(Sample { t, y });
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// RHS for the scalar exponential test problem y' = y.
    fn exp_rhs(_t: f64, y: &State<1>) -> State<1> {
        [y[0]]
    }

    #[test]
    fn single_step_matches_exponential() {
        // y' = y, y(0) = 1 => y(h) = e^h. RK4 single-step error ~ O(h^5).
        let h = 0.05;
        let y = rk4_step(0.0, [1.0], h, exp_rhs);
        let exact = h.exp();
        assert!((y[0] - exact).abs() < 1e-7, "y={y:?} exact={exact}");
    }

    #[test]
    fn integrate_exponential_global_accuracy() {
        // Integrate y' = y from 0 to 1; exact y(1) = e.
        let traj = integrate(exp_rhs, [1.0], 0.0, 1.0, 0.01).unwrap();
        let last = traj.last().unwrap();
        let e = std::f64::consts::E;
        assert!(
            (last.y[0] - e).abs() < 1e-6,
            "got {got} want {e}",
            got = last.y[0]
        );
        // 0..1 in steps of 0.01 => 100 steps => 101 samples.
        assert_eq!(traj.len(), 101);
        assert!((traj[0].t - 0.0).abs() < 1e-12);
        assert!((last.t - 1.0).abs() < 1e-9);
    }

    #[test]
    fn fourth_order_convergence() {
        // Halving dt must cut the global error by roughly 2^4 = 16.
        let e = std::f64::consts::E;
        let err = |dt: f64| {
            let traj = integrate(exp_rhs, [1.0], 0.0, 1.0, dt).unwrap();
            (traj.last().unwrap().y[0] - e).abs()
        };
        let e_coarse = err(0.1);
        let e_fine = err(0.05);
        let ratio = e_coarse / e_fine;
        // Fourth-order => ratio ~ 16; allow a generous band.
        assert!(
            (12.0..20.0).contains(&ratio),
            "convergence ratio {ratio} not ~16 (coarse={e_coarse:e}, fine={e_fine:e})"
        );
    }

    #[test]
    fn rejects_nonpositive_dt() {
        let err = integrate(exp_rhs, [1.0], 0.0, 1.0, 0.0).unwrap_err();
        assert_eq!(err.code(), "popdynamics.invalid");
        let err = integrate(exp_rhs, [1.0], 0.0, 1.0, -0.1).unwrap_err();
        assert_eq!(err.code(), "popdynamics.invalid");
    }

    #[test]
    fn rejects_bad_window() {
        let err = integrate(exp_rhs, [1.0], 1.0, 1.0, 0.1).unwrap_err();
        assert_eq!(err.code(), "popdynamics.invalid");
        let err = integrate(exp_rhs, [1.0], 2.0, 1.0, 0.1).unwrap_err();
        assert_eq!(err.code(), "popdynamics.invalid");
    }

    #[test]
    fn rejects_too_many_steps() {
        // 1e6 span / 1e-6 dt = 1e12 steps >> MAX_STEPS.
        let err = integrate(exp_rhs, [1.0], 0.0, 1.0e6, 1.0e-6).unwrap_err();
        assert_eq!(err.code(), "popdynamics.too_many_steps");
        assert_eq!(err.category(), crate::error::ErrorCategory::Limit);
    }
}
