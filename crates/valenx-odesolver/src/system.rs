//! Vec-of-state (system) fixed-step integration of `dy/dt = f(t, y)`.
//!
//! Here the state `y` is a `Vec<f64>` of arbitrary (fixed) length `n`, and the
//! derivative closure `Fn(f64, &[f64]) -> Vec<f64>` returns the per-component
//! rates. This is the form used for higher-order ODEs reduced to first order
//! (e.g. the harmonic oscillator `[x, v]` with `dx/dt = v`,
//! `dv/dt = -ω² x`).
//!
//! Use [`integrate`] for the whole trajectory or [`integrate_final`] for just
//! the endpoint. The single-step helpers ([`euler_step`], [`heun_step`],
//! [`rk4_step`]) advance one fixed step.

use crate::error::OdeError;
use crate::method::Method;

/// Componentwise `out = a + scale * b`. Both slices must share length `n`.
#[inline]
fn axpy(a: &[f64], scale: f64, b: &[f64]) -> Vec<f64> {
    a.iter().zip(b).map(|(ai, bi)| ai + scale * bi).collect()
}

/// Validate that the derivative output length matches the state length.
fn check_dim(expected: usize, got: &[f64]) -> Result<(), OdeError> {
    match OdeError::dimension_mismatch(expected, got.len()) {
        Some(e) => Err(e),
        None => Ok(()),
    }
}

/// One explicit-Euler step of size `h` from `(t, y)` for a system.
///
/// # Errors
///
/// Returns [`OdeError::DimensionMismatch`] if the derivative closure returns a
/// vector whose length differs from `y`.
///
/// # Examples
///
/// ```
/// use valenx_odesolver::system::euler_step;
/// // dy/dt = y (componentwise), one step h = 0.5 from [1, 2] -> [1.5, 3.0].
/// let y1 = euler_step(|_t, y| y.to_vec(), 0.0, &[1.0, 2.0], 0.5).unwrap();
/// assert!((y1[0] - 1.5).abs() < 1e-12 && (y1[1] - 3.0).abs() < 1e-12);
/// ```
pub fn euler_step<F>(f: F, t: f64, y: &[f64], h: f64) -> Result<Vec<f64>, OdeError>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let k1 = f(t, y);
    check_dim(y.len(), &k1)?;
    Ok(axpy(y, h, &k1))
}

/// One Heun (improved-Euler / RK2) step of size `h` from `(t, y)` for a system.
///
/// # Errors
///
/// Returns [`OdeError::DimensionMismatch`] if any derivative evaluation
/// returns a vector whose length differs from `y`.
///
/// # Examples
///
/// ```
/// use valenx_odesolver::system::heun_step;
/// // dy/dt = y from [1.0], h = 1.0: k1 = 1, k2 = 2, y1 = 1 + 0.5*(1+2) = 2.5.
/// let y1 = heun_step(|_t, y| y.to_vec(), 0.0, &[1.0], 1.0).unwrap();
/// assert!((y1[0] - 2.5).abs() < 1e-12);
/// ```
pub fn heun_step<F>(f: F, t: f64, y: &[f64], h: f64) -> Result<Vec<f64>, OdeError>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let n = y.len();
    let k1 = f(t, y);
    check_dim(n, &k1)?;
    let y_pred = axpy(y, h, &k1);
    let k2 = f(t + h, &y_pred);
    check_dim(n, &k2)?;
    Ok((0..n).map(|i| y[i] + 0.5 * h * (k1[i] + k2[i])).collect())
}

/// One classical RK4 step of size `h` from `(t, y)` for a system.
///
/// # Errors
///
/// Returns [`OdeError::DimensionMismatch`] if any derivative evaluation
/// returns a vector whose length differs from `y`.
///
/// # Examples
///
/// ```
/// use valenx_odesolver::system::rk4_step;
/// // dy/dt = y from [1.0], h = 1.0 -> 1 + 1 + 1/2 + 1/6 + 1/24 = 2.708333...
/// let y1 = rk4_step(|_t, y| y.to_vec(), 0.0, &[1.0], 1.0).unwrap();
/// assert!((y1[0] - 2.708_333_333_333_333).abs() < 1e-12);
/// ```
pub fn rk4_step<F>(f: F, t: f64, y: &[f64], h: f64) -> Result<Vec<f64>, OdeError>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let n = y.len();
    let k1 = f(t, y);
    check_dim(n, &k1)?;
    let y2 = axpy(y, 0.5 * h, &k1);
    let k2 = f(t + 0.5 * h, &y2);
    check_dim(n, &k2)?;
    let y3 = axpy(y, 0.5 * h, &k2);
    let k3 = f(t + 0.5 * h, &y3);
    check_dim(n, &k3)?;
    let y4 = axpy(y, h, &k3);
    let k4 = f(t + h, &y4);
    check_dim(n, &k4)?;
    Ok((0..n)
        .map(|i| y[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
        .collect())
}

/// Advance a single system step with the scheme selected by `method`.
fn step<F>(method: Method, f: &F, t: f64, y: &[f64], h: f64) -> Result<Vec<f64>, OdeError>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    match method {
        Method::Euler => euler_step(f, t, y, h),
        Method::Heun => heun_step(f, t, y, h),
        Method::Rk4 => rk4_step(f, t, y, h),
    }
}

/// Integrate the system `dy/dt = f(t, y)` from `t0` for `steps` fixed steps of
/// size `dt`.
///
/// Returns the trajectory `[y0, y1, ..., y_steps]` as a `Vec` of state vectors
/// (outer length `steps + 1`); each inner vector has the same length as `y0`.
///
/// # Errors
///
/// Returns [`OdeError::BadStep`] if `dt` is not strictly positive and finite,
/// [`OdeError::BadStepCount`] if `steps` is zero, [`OdeError::NonFinite`] if any
/// component of `y0` (or `t0`) is not finite, and
/// [`OdeError::DimensionMismatch`] if the derivative closure ever returns a
/// vector of the wrong length.
///
/// # Examples
///
/// ```
/// use valenx_odesolver::{system::integrate, Method};
/// // Harmonic oscillator with omega = 1: state [x, v], x(0)=1, v(0)=0.
/// // dx/dt = v, dv/dt = -x. Energy E = 0.5*(v^2 + x^2) is conserved.
/// let f = |_t: f64, y: &[f64]| vec![y[1], -y[0]];
/// let traj = integrate(Method::Rk4, f, 0.0, &[1.0, 0.0], 0.01, 100).unwrap();
/// assert_eq!(traj.len(), 101);
/// let e0 = 0.5 * (traj[0][1].powi(2) + traj[0][0].powi(2));
/// let ef = 0.5 * (traj[100][1].powi(2) + traj[100][0].powi(2));
/// assert!((ef - e0).abs() < 1e-6);
/// ```
pub fn integrate<F>(
    method: Method,
    f: F,
    t0: f64,
    y0: &[f64],
    dt: f64,
    steps: usize,
) -> Result<Vec<Vec<f64>>, OdeError>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    if let Some(e) = OdeError::bad_step(dt) {
        return Err(e);
    }
    if let Some(e) = OdeError::bad_step_count(steps) {
        return Err(e);
    }
    if let Some(e) = OdeError::non_finite("t0", t0) {
        return Err(e);
    }
    for (i, &c) in y0.iter().enumerate() {
        // Tag the offending component; `"y0[i]"` would need allocation, so the
        // static name plus the value is enough to localise in practice.
        if !c.is_finite() {
            return Err(OdeError::NonFinite {
                name: "y0_component",
                value: c,
            });
        }
        let _ = i;
    }

    let mut out: Vec<Vec<f64>> = Vec::with_capacity(steps + 1);
    out.push(y0.to_vec());
    let mut t = t0;
    let mut y = y0.to_vec();
    for _ in 0..steps {
        y = step(method, &f, t, &y, dt)?;
        t += dt;
        out.push(y.clone());
    }
    Ok(out)
}

/// Convenience wrapper around [`integrate`] returning only the final state.
///
/// # Errors
///
/// Propagates every error from [`integrate`].
///
/// # Examples
///
/// ```
/// use valenx_odesolver::{system::integrate_final, Method};
/// let f = |_t: f64, y: &[f64]| vec![y[1], -y[0]];
/// let yf = integrate_final(Method::Rk4, f, 0.0, &[1.0, 0.0], 0.001, 1000).unwrap();
/// // At t = 1, x = cos(1), v = -sin(1).
/// assert!((yf[0] - 1.0_f64.cos()).abs() < 1e-6);
/// assert!((yf[1] + 1.0_f64.sin()).abs() < 1e-6);
/// ```
pub fn integrate_final<F>(
    method: Method,
    f: F,
    t0: f64,
    y0: &[f64],
    dt: f64,
    steps: usize,
) -> Result<Vec<f64>, OdeError>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let traj = integrate(method, f, t0, y0, dt, steps)?;
    Ok(traj.into_iter().last().expect("trajectory is non-empty"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    /// Undamped unit-frequency oscillator: dx/dt = v, dv/dt = -x.
    fn oscillator(_t: f64, y: &[f64]) -> Vec<f64> {
        vec![y[1], -y[0]]
    }

    /// Mechanical energy E = 0.5 (v^2 + x^2) for the unit oscillator.
    fn energy(state: &[f64]) -> f64 {
        0.5 * (state[1] * state[1] + state[0] * state[0])
    }

    #[test]
    fn harmonic_oscillator_is_roughly_energy_conserving_rk4() {
        // Integrate many periods; RK4 keeps energy almost flat (tiny secular
        // drift). x(0) = 1, v(0) = 0 -> E0 = 0.5.
        let traj = integrate(Method::Rk4, oscillator, 0.0, &[1.0, 0.0], 0.01, 6283).unwrap();
        let e0 = energy(&traj[0]);
        assert!(close(e0, 0.5, 1e-15), "initial energy {e0}");
        // Worst-case relative energy drift over ~10 periods stays small.
        let mut max_rel = 0.0_f64;
        for s in &traj {
            let rel = ((energy(s) - e0) / e0).abs();
            max_rel = max_rel.max(rel);
        }
        assert!(
            max_rel < 1e-3,
            "rk4 oscillator energy drift {max_rel} too large"
        );
    }

    #[test]
    fn harmonic_oscillator_tracks_analytic_solution() {
        // x(t) = cos(t), v(t) = -sin(t) for x(0)=1, v(0)=0.
        let yf = integrate_final(Method::Rk4, oscillator, 0.0, &[1.0, 0.0], 0.001, 1571).unwrap();
        let t: f64 = 1.571; // ~ pi/2: x ~ 0, v ~ -1.
        assert!(
            close(yf[0], t.cos(), 1e-6),
            "x endpoint {} vs cos {}",
            yf[0],
            t.cos()
        );
        assert!(
            close(yf[1], -t.sin(), 1e-6),
            "v endpoint {} vs -sin {}",
            yf[1],
            -t.sin()
        );
    }

    #[test]
    fn euler_drifts_more_than_rk4_on_oscillator_energy() {
        // Explicit Euler injects energy into the undamped oscillator
        // (the classic Euler instability); RK4 does not. Compare energy
        // growth at a horizon where the difference is unambiguous.
        let dt = 0.05;
        let steps = 200; // ~ 1.6 periods
        let euler = integrate(Method::Euler, oscillator, 0.0, &[1.0, 0.0], dt, steps).unwrap();
        let rk4 = integrate(Method::Rk4, oscillator, 0.0, &[1.0, 0.0], dt, steps).unwrap();
        let e0 = 0.5;
        let euler_growth = (energy(euler.last().unwrap()) - e0).abs();
        let rk4_growth = (energy(rk4.last().unwrap()) - e0).abs();
        assert!(
            rk4_growth < euler_growth,
            "rk4 energy growth {rk4_growth} should be below euler {euler_growth}"
        );
        // And Euler should visibly gain energy.
        assert!(
            energy(euler.last().unwrap()) > e0,
            "explicit euler should gain energy on the oscillator"
        );
    }

    #[test]
    fn system_decay_matches_exp_per_component() {
        // dy/dt = -y on a 3-vector with different initial amplitudes; each
        // component decays like its own exp(-t).
        let f = |_t: f64, y: &[f64]| y.iter().map(|v| -v).collect::<Vec<_>>();
        let y0 = [1.0, 2.0, -3.0];
        let yf = integrate_final(Method::Rk4, f, 0.0, &y0, 0.001, 1000).unwrap();
        let factor = (-1.0_f64).exp();
        for (i, &start) in y0.iter().enumerate() {
            let exact = start * factor;
            assert!(
                close(yf[i], exact, 1e-7),
                "component {i}: {} vs {exact}",
                yf[i]
            );
        }
    }

    #[test]
    fn system_rk4_step_halving_cuts_error_by_about_sixteen() {
        // Same O(h^4) law in the vector case, measured on the oscillator's x.
        let exact_x = 1.0_f64.cos(); // x(1) = cos(1)
        let coarse =
            integrate_final(Method::Rk4, oscillator, 0.0, &[1.0, 0.0], 0.1, 10).unwrap()[0];
        let fine = integrate_final(Method::Rk4, oscillator, 0.0, &[1.0, 0.0], 0.05, 20).unwrap()[0];
        let err_coarse = (coarse - exact_x).abs();
        let err_fine = (fine - exact_x).abs();
        let ratio = err_coarse / err_fine;
        assert!(
            ratio > 12.0 && ratio < 20.0,
            "system rk4 halving ratio {ratio} (coarse {err_coarse}, fine {err_fine}) not ~16"
        );
    }

    #[test]
    fn dimension_mismatch_is_reported() {
        // Closure returns a wrong-length derivative -> structured error.
        let f = |_t: f64, _y: &[f64]| vec![0.0]; // length 1, state length 2
        let err = integrate(Method::Euler, f, 0.0, &[1.0, 0.0], 0.1, 1).unwrap_err();
        assert!(
            matches!(
                err,
                OdeError::DimensionMismatch {
                    expected: 2,
                    actual: 1
                }
            ),
            "unexpected error {err:?}"
        );
    }

    #[test]
    fn single_system_steps_match_hand_values() {
        // dy/dt = y on [1.0], h = 1.0: euler 2.0, heun 2.5, rk4 2.708333...
        let g = |_t: f64, y: &[f64]| y.to_vec();
        let e = euler_step(g, 0.0, &[1.0], 1.0).unwrap();
        let h = heun_step(g, 0.0, &[1.0], 1.0).unwrap();
        let r = rk4_step(g, 0.0, &[1.0], 1.0).unwrap();
        assert!(close(e[0], 2.0, 1e-12), "euler {}", e[0]);
        assert!(close(h[0], 2.5, 1e-12), "heun {}", h[0]);
        assert!(close(r[0], 2.708_333_333_333_333_3, 1e-12), "rk4 {}", r[0]);
    }

    #[test]
    fn bad_system_inputs_are_rejected() {
        let f = |_t: f64, y: &[f64]| y.to_vec();
        assert!(matches!(
            integrate(Method::Rk4, f, 0.0, &[1.0], 0.0, 5),
            Err(OdeError::BadStep { .. })
        ));
        assert!(matches!(
            integrate(Method::Rk4, f, 0.0, &[1.0], 0.1, 0),
            Err(OdeError::BadStepCount)
        ));
        assert!(matches!(
            integrate(Method::Rk4, f, 0.0, &[f64::NAN], 0.1, 5),
            Err(OdeError::NonFinite { .. })
        ));
    }

    #[test]
    fn trajectory_shape_is_consistent() {
        let traj = integrate(Method::Heun, oscillator, 0.0, &[1.0, 0.0], 0.1, 5).unwrap();
        assert_eq!(traj.len(), 6);
        for s in &traj {
            assert_eq!(s.len(), 2);
        }
    }
}
