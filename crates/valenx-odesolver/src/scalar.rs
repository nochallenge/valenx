//! Scalar (single-state) fixed-step integration of `dy/dt = f(t, y)`.
//!
//! Use [`integrate`] to drive a chosen [`Method`] forward by a fixed number
//! of steps, or the individual [`euler_step`] / [`heun_step`] / [`rk4_step`]
//! helpers to advance a single step. The derivative is any callable
//! `Fn(f64, f64) -> f64`, taking the current time and state and returning
//! `dy/dt`.

use crate::error::OdeError;
use crate::method::Method;

/// One explicit-Euler step of size `h` from `(t, y)`.
///
/// Returns `y + h * f(t, y)`. See the [method module](crate::method) for the
/// full scheme.
///
/// # Examples
///
/// ```
/// use valenx_odesolver::scalar::euler_step;
/// // dy/dt = y, one step of h = 0.5 from y = 1.0 gives 1 + 0.5*1 = 1.5.
/// let y1 = euler_step(|_t, y| y, 0.0, 1.0, 0.5);
/// assert!((y1 - 1.5).abs() < 1e-12);
/// ```
pub fn euler_step<F>(f: F, t: f64, y: f64, h: f64) -> f64
where
    F: Fn(f64, f64) -> f64,
{
    y + h * f(t, y)
}

/// One Heun (improved-Euler / RK2) step of size `h` from `(t, y)`.
///
/// # Examples
///
/// ```
/// use valenx_odesolver::scalar::heun_step;
/// // dy/dt = y from y = 1.0, h = 1.0:
/// // k1 = 1, k2 = f(1, 1 + 1) = 2, y1 = 1 + 0.5*(1 + 2) = 2.5.
/// let y1 = heun_step(|_t, y| y, 0.0, 1.0, 1.0);
/// assert!((y1 - 2.5).abs() < 1e-12);
/// ```
pub fn heun_step<F>(f: F, t: f64, y: f64, h: f64) -> f64
where
    F: Fn(f64, f64) -> f64,
{
    let k1 = f(t, y);
    let k2 = f(t + h, y + h * k1);
    y + 0.5 * h * (k1 + k2)
}

/// One classical RK4 step of size `h` from `(t, y)`.
///
/// # Examples
///
/// ```
/// use valenx_odesolver::scalar::rk4_step;
/// // dy/dt = y from y = 1.0, h = 1.0 approximates e ~= 2.71828; RK4's
/// // single-step estimate is 1 + 1 + 1/2 + 1/6 + 1/24 = 2.708333...
/// let y1 = rk4_step(|_t, y| y, 0.0, 1.0, 1.0);
/// assert!((y1 - 2.708_333_333_333_333).abs() < 1e-12);
/// ```
pub fn rk4_step<F>(f: F, t: f64, y: f64, h: f64) -> f64
where
    F: Fn(f64, f64) -> f64,
{
    let k1 = f(t, y);
    let k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
    let k3 = f(t + 0.5 * h, y + 0.5 * h * k2);
    let k4 = f(t + h, y + h * k3);
    y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
}

/// Advance a single step with the scheme selected by `method`.
fn step<F>(method: Method, f: &F, t: f64, y: f64, h: f64) -> f64
where
    F: Fn(f64, f64) -> f64,
{
    match method {
        Method::Euler => euler_step(f, t, y, h),
        Method::Heun => heun_step(f, t, y, h),
        Method::Rk4 => rk4_step(f, t, y, h),
    }
}

/// Integrate `dy/dt = f(t, y)` from `t0` for `steps` fixed steps of size `dt`.
///
/// Returns the full trajectory of states `[y0, y1, ..., y_steps]` (length
/// `steps + 1`); element `i` is the approximation at time `t0 + i*dt`. The
/// derivative is sampled with the scheme chosen by `method`.
///
/// # Errors
///
/// Returns [`OdeError::BadStep`] if `dt` is not strictly positive and finite,
/// [`OdeError::BadStepCount`] if `steps` is zero, or [`OdeError::NonFinite`]
/// if `y0` is not finite.
///
/// Note that the integration itself is *not* re-validated for finiteness at
/// every step: a stiff problem run with too large a step can legitimately
/// overflow to `±∞`, and surfacing that as data (rather than an error) lets
/// the caller observe the blow-up. Inputs, however, are always checked.
///
/// # Examples
///
/// ```
/// use valenx_odesolver::{scalar::integrate, Method};
/// // dy/dt = y, y(0) = 1 -> y(t) = e^t. After 100 RK4 steps to t = 1,
/// // the final state should match e to high accuracy.
/// let traj = integrate(Method::Rk4, |_t, y| y, 0.0, 1.0, 0.01, 100).unwrap();
/// assert_eq!(traj.len(), 101);
/// let e = std::f64::consts::E;
/// assert!((traj[100] - e).abs() < 1e-8);
/// ```
pub fn integrate<F>(
    method: Method,
    f: F,
    t0: f64,
    y0: f64,
    dt: f64,
    steps: usize,
) -> Result<Vec<f64>, OdeError>
where
    F: Fn(f64, f64) -> f64,
{
    if let Some(e) = OdeError::bad_step(dt) {
        return Err(e);
    }
    if let Some(e) = OdeError::bad_step_count(steps) {
        return Err(e);
    }
    if let Some(e) = OdeError::non_finite("y0", y0) {
        return Err(e);
    }
    if let Some(e) = OdeError::non_finite("t0", t0) {
        return Err(e);
    }

    let mut out = Vec::with_capacity(steps + 1);
    out.push(y0);
    let mut t = t0;
    let mut y = y0;
    for _ in 0..steps {
        y = step(method, &f, t, y, dt);
        t += dt;
        out.push(y);
    }
    Ok(out)
}

/// Convenience wrapper around [`integrate`] returning only the final state at
/// `t0 + steps*dt`.
///
/// # Errors
///
/// Propagates every error from [`integrate`].
///
/// # Examples
///
/// ```
/// use valenx_odesolver::{scalar::integrate_final, Method};
/// // dy/dt = -y decays: y(1) = e^{-1} ~= 0.3679.
/// let yf = integrate_final(Method::Rk4, |_t, y| -y, 0.0, 1.0, 0.001, 1000).unwrap();
/// assert!((yf - (-1.0_f64).exp()).abs() < 1e-8);
/// ```
pub fn integrate_final<F>(
    method: Method,
    f: F,
    t0: f64,
    y0: f64,
    dt: f64,
    steps: usize,
) -> Result<f64, OdeError>
where
    F: Fn(f64, f64) -> f64,
{
    let traj = integrate(method, f, t0, y0, dt, steps)?;
    // `integrate` guarantees a non-empty trajectory (steps >= 1), so the last
    // element always exists.
    Ok(*traj.last().expect("trajectory is non-empty"))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference assertion helper; never uses `==` on floats.
    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn rk4_matches_exp_to_high_accuracy() {
        // dy/dt = y, y(0) = 1 -> exp(t). RK4 with a modest step nails it.
        let yf = integrate_final(Method::Rk4, |_t, y| y, 0.0, 1.0, 0.01, 100).unwrap();
        let exact = std::f64::consts::E;
        assert!(
            close(yf, exact, 1e-8),
            "rk4 exp endpoint {yf} vs exact {exact}"
        );
    }

    #[test]
    fn euler_is_less_accurate_than_rk4_at_same_step() {
        // Same problem, same step, same horizon: RK4's endpoint error must be
        // strictly smaller than Euler's.
        let exact = std::f64::consts::E;
        let dt = 0.05;
        let steps = 20; // t: 0 -> 1
        let y_euler = integrate_final(Method::Euler, |_t, y| y, 0.0, 1.0, dt, steps).unwrap();
        let y_rk4 = integrate_final(Method::Rk4, |_t, y| y, 0.0, 1.0, dt, steps).unwrap();
        let err_euler = (y_euler - exact).abs();
        let err_rk4 = (y_rk4 - exact).abs();
        assert!(
            err_rk4 < err_euler,
            "rk4 err {err_rk4} should beat euler err {err_euler}"
        );
        // Heun should sit strictly between the two in accuracy.
        let y_heun = integrate_final(Method::Heun, |_t, y| y, 0.0, 1.0, dt, steps).unwrap();
        let err_heun = (y_heun - exact).abs();
        assert!(
            err_rk4 < err_heun && err_heun < err_euler,
            "expected rk4 {err_rk4} < heun {err_heun} < euler {err_euler}"
        );
    }

    #[test]
    fn decay_solution_decreases_monotonically() {
        // dy/dt = -y, y(0) = 1 -> exp(-t): strictly decreasing, stays positive.
        let traj = integrate(Method::Rk4, |_t, y| -y, 0.0, 1.0, 0.01, 200).unwrap();
        for w in traj.windows(2) {
            assert!(w[1] < w[0], "decay not monotone: {} -> {}", w[0], w[1]);
        }
        let last = *traj.last().unwrap();
        assert!(last > 0.0, "decay should stay positive, got {last}");
        let exact = (-2.0_f64).exp(); // t = 2.0
        assert!(
            close(last, exact, 1e-8),
            "decay endpoint {last} vs exact {exact}"
        );
    }

    #[test]
    fn rk4_step_halving_cuts_error_by_about_sixteen() {
        // Fourth-order global error ~ O(h^4): halving h cuts error ~16x.
        let exact = std::f64::consts::E;
        // Coarse: dt = 0.1, 10 steps -> t = 1.
        let coarse = integrate_final(Method::Rk4, |_t, y| y, 0.0, 1.0, 0.1, 10).unwrap();
        // Fine: dt = 0.05, 20 steps -> t = 1.
        let fine = integrate_final(Method::Rk4, |_t, y| y, 0.0, 1.0, 0.05, 20).unwrap();
        let err_coarse = (coarse - exact).abs();
        let err_fine = (fine - exact).abs();
        let ratio = err_coarse / err_fine;
        // Theoretical 16; allow a generous band for finite-precision /
        // higher-order remainder terms.
        assert!(
            ratio > 12.0 && ratio < 20.0,
            "rk4 halving ratio {ratio} (coarse {err_coarse}, fine {err_fine}) not ~16"
        );
    }

    #[test]
    fn euler_step_halving_cuts_error_by_about_two() {
        // First-order global error ~ O(h): halving h cuts error ~2x.
        let exact = std::f64::consts::E;
        let coarse = integrate_final(Method::Euler, |_t, y| y, 0.0, 1.0, 0.001, 1000).unwrap();
        let fine = integrate_final(Method::Euler, |_t, y| y, 0.0, 1.0, 0.0005, 2000).unwrap();
        let err_coarse = (coarse - exact).abs();
        let err_fine = (fine - exact).abs();
        let ratio = err_coarse / err_fine;
        assert!(
            ratio > 1.8 && ratio < 2.2,
            "euler halving ratio {ratio} (coarse {err_coarse}, fine {err_fine}) not ~2"
        );
    }

    #[test]
    fn single_steps_match_hand_computed_values() {
        // dy/dt = y from y = 1.0, h = 1.0.
        let euler = euler_step(|_t, y| y, 0.0, 1.0, 1.0); // 2.0
        let heun = heun_step(|_t, y| y, 0.0, 1.0, 1.0); // 2.5
        let rk4 = rk4_step(|_t, y| y, 0.0, 1.0, 1.0); // 2.7083333...
        assert!(close(euler, 2.0, 1e-12), "euler one-step {euler}");
        assert!(close(heun, 2.5, 1e-12), "heun one-step {heun}");
        assert!(
            close(rk4, 2.708_333_333_333_333_3, 1e-12),
            "rk4 one-step {rk4}"
        );
    }

    #[test]
    fn time_dependent_rhs_integrates_polynomial_exactly_with_rk4() {
        // dy/dt = 2t, y(0) = 0 -> y = t^2. RK4 is exact for cubics, so a
        // quadratic solution is reproduced to round-off.
        let yf = integrate_final(Method::Rk4, |t, _y| 2.0 * t, 0.0, 0.0, 0.1, 30).unwrap();
        let t = 3.0;
        assert!(
            close(yf, t * t, 1e-10),
            "polynomial endpoint {yf} vs {}",
            t * t
        );
    }

    #[test]
    fn bad_inputs_are_rejected() {
        assert!(matches!(
            integrate(Method::Rk4, |_t, y| y, 0.0, 1.0, 0.0, 10),
            Err(OdeError::BadStep { .. })
        ));
        assert!(matches!(
            integrate(Method::Rk4, |_t, y| y, 0.0, 1.0, -0.1, 10),
            Err(OdeError::BadStep { .. })
        ));
        assert!(matches!(
            integrate(Method::Rk4, |_t, y| y, 0.0, 1.0, 0.1, 0),
            Err(OdeError::BadStepCount)
        ));
        assert!(matches!(
            integrate(Method::Rk4, |_t, y| y, 0.0, f64::NAN, 0.1, 10),
            Err(OdeError::NonFinite { .. })
        ));
    }

    #[test]
    fn trajectory_length_and_first_element() {
        let traj = integrate(Method::Heun, |_t, y| y, 0.0, 2.0, 0.1, 7).unwrap();
        assert_eq!(traj.len(), 8);
        assert!(close(traj[0], 2.0, 1e-15), "first element {}", traj[0]);
    }
}
