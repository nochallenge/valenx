//! Geodesic integration: orbits, light bending, perihelion precession and
//! photon capture.
//!
//! A geodesic obeys `d²xᵘ/dλ² = −Γᵘ_αβ (dxᵃ/dλ)(dxᵝ/dλ)`. We integrate the
//! first-order system `(xᵘ, uᵘ)` with an adaptive Runge–Kutta-4 step (error
//! controlled by step-doubling), pulling the Christoffel symbols from the
//! curvature engine at each stage. Everything is in the equatorial plane
//! `θ = π/2` (where the spacetimes are reflection-symmetric, so a geodesic
//! started there stays there).
//!
//! Geometrized units `G = c = 1`.

// The acceleration is a Christoffel double-contraction over the 4 indices;
// explicit index loops mirror the formula and index several tensors at once.
#![allow(clippy::needless_range_loop)]

use std::f64::consts::{FRAC_PI_2, PI};

use crate::curvature::christoffel_at;
use crate::metric::Spacetime;
use crate::spacetimes::KerrNewman;
use crate::{RelativityError, Result};

/// Whether a geodesic is that of light or of a massive particle.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Kind {
    /// Null geodesic (photon): `g_μν uᵘ uᵛ = 0`.
    Null,
    /// Timelike geodesic (massive particle), affine parameter = proper time:
    /// `g_μν uᵘ uᵛ = −1`.
    Timelike,
}

impl Kind {
    /// The conserved norm `g_μν uᵘ uᵛ` for this kind (`0` null, `−1` timelike).
    pub fn norm(self) -> f64 {
        match self {
            Kind::Null => 0.0,
            Kind::Timelike => -1.0,
        }
    }
}

/// Position `xᵘ = (t, r, θ, φ)` and 4-velocity `uᵘ = dxᵘ/dλ` of a geodesic.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct GeodesicState {
    /// Coordinates `(t, r, θ, φ)`.
    pub x: [f64; 4],
    /// 4-velocity `dxᵘ/dλ`.
    pub u: [f64; 4],
}

/// Why [`integrate_geodesic`] stopped.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum StopReason {
    /// Fell to/below the capture radius (plunging toward the hole).
    Captured,
    /// Reached the escape radius while moving outward.
    Escaped,
    /// Hit the affine-parameter or step-count limit.
    Limit,
}

/// Tuning for [`integrate_geodesic`].
#[derive(Clone, Copy, Debug)]
pub struct GeodesicOptions {
    /// Initial affine step.
    pub step: f64,
    /// Per-step local error tolerance for the adaptive controller.
    pub tol: f64,
    /// Stop (`Captured`) once `r ≤ r_capture`.
    pub r_capture: f64,
    /// Stop (`Escaped`) once `r ≥ r_escape` while moving outward.
    pub r_escape: f64,
    /// Hard cap on total affine length.
    pub max_lambda: f64,
    /// Hard cap on number of accepted steps.
    pub max_steps: usize,
}

impl Default for GeodesicOptions {
    fn default() -> Self {
        Self {
            step: 1.0,
            tol: 1e-11,
            r_capture: 0.0,
            r_escape: 1e9,
            max_lambda: 1e12,
            max_steps: 2_000_000,
        }
    }
}

/// An integrated geodesic path.
#[derive(Clone, Debug)]
pub struct Trajectory {
    /// Sampled states, one per accepted step (including the initial state).
    pub states: Vec<GeodesicState>,
    /// Affine parameter at each sample.
    pub lambda: Vec<f64>,
    /// Why integration stopped.
    pub stop: StopReason,
}

impl Trajectory {
    /// The final sampled state.
    pub fn last(&self) -> GeodesicState {
        *self
            .states
            .last()
            .expect("trajectory always has the initial state")
    }
}

/// Conserved energy `E = −g_{tμ} uᵘ` (from the time-translation Killing vector).
pub fn energy(bh: &KerrNewman, st: &GeodesicState) -> f64 {
    let g = bh.metric::<f64>(st.x);
    -(0..4).map(|m| g[0][m] * st.u[m]).sum::<f64>()
}

/// Conserved axial angular momentum `L = g_{φμ} uᵘ` (from the axial Killing
/// vector).
pub fn angular_momentum(bh: &KerrNewman, st: &GeodesicState) -> f64 {
    let g = bh.metric::<f64>(st.x);
    (0..4).map(|m| g[3][m] * st.u[m]).sum::<f64>()
}

/// The norm `g_μν uᵘ uᵛ` (should stay `0` for null, `−1` for timelike).
pub fn norm(bh: &KerrNewman, st: &GeodesicState) -> f64 {
    let g = bh.metric::<f64>(st.x);
    let mut s = 0.0;
    for a in 0..4 {
        for b in 0..4 {
            s += g[a][b] * st.u[a] * st.u[b];
        }
    }
    s
}

/// Build an equatorial geodesic state from conserved `(E, L)` at radius `r`.
///
/// Solves the metric for `(uᵗ, uᵠ)` given `E` and `L`, sets `uᶿ = 0`, and gets
/// `uʳ` from the norm condition (sign chosen by `ingoing`). At a radial turning
/// point `uʳ = 0` regardless of `ingoing`.
///
/// # Errors
/// [`RelativityError::OutsideDomain`] if no real `uʳ` exists there (the radius
/// is forbidden for that `(E, L)`), or the metric is degenerate.
pub fn equatorial_state(
    bh: &KerrNewman,
    r: f64,
    e: f64,
    l: f64,
    kind: Kind,
    ingoing: bool,
) -> Result<GeodesicState> {
    let x = [0.0, r, FRAC_PI_2, 0.0];
    let g = bh.metric::<f64>(x);
    let (g_tt, g_tp, g_pp, g_rr) = (g[0][0], g[0][3], g[3][3], g[1][1]);
    // [ -g_tt  -g_tp ] [u^t]   [E]
    // [  g_tp   g_pp ] [u^φ] = [L]
    let det = -g_tt * g_pp + g_tp * g_tp;
    if det.abs() < 1e-300 {
        return Err(RelativityError::OutsideDomain(format!(
            "degenerate t-φ metric block at r={r}"
        )));
    }
    let ut = (e * g_pp + g_tp * l) / det;
    let up = (-g_tt * l - g_tp * e) / det;
    let rhs = kind.norm() - (g_tt * ut * ut + 2.0 * g_tp * ut * up + g_pp * up * up);
    let ur2 = rhs / g_rr;
    // A genuinely forbidden radius gives a clearly negative uʳ²; tiny negative
    // values are floating-point noise at a turning point (uʳ = 0), so clamp.
    if !ur2.is_finite() || ur2 < -1e-9 {
        return Err(RelativityError::OutsideDomain(format!(
            "no real radial velocity at r={r} (uʳ² = {ur2})"
        )));
    }
    let ur_mag = ur2.max(0.0).sqrt();
    let ur = if ingoing { -ur_mag } else { ur_mag };
    Ok(GeodesicState {
        x,
        u: [ut, ur, 0.0, up],
    })
}

/// Right-hand side of the geodesic system as an 8-vector `(dx, du)`.
fn deriv(bh: &KerrNewman, y: &[f64; 8]) -> Result<[f64; 8]> {
    let x = [y[0], y[1], y[2], y[3]];
    let u = [y[4], y[5], y[6], y[7]];
    let gamma = christoffel_at(bh, x)?;
    let mut a = [0.0; 4];
    for mu in 0..4 {
        let mut acc = 0.0;
        for al in 0..4 {
            for be in 0..4 {
                acc -= gamma[mu][al][be] * u[al] * u[be];
            }
        }
        a[mu] = acc;
    }
    Ok([u[0], u[1], u[2], u[3], a[0], a[1], a[2], a[3]])
}

/// One classical RK4 step of size `h`.
fn rk4(bh: &KerrNewman, y: &[f64; 8], h: f64) -> Result<[f64; 8]> {
    let k1 = deriv(bh, y)?;
    let y2: [f64; 8] = std::array::from_fn(|i| y[i] + 0.5 * h * k1[i]);
    let k2 = deriv(bh, &y2)?;
    let y3: [f64; 8] = std::array::from_fn(|i| y[i] + 0.5 * h * k2[i]);
    let k3 = deriv(bh, &y3)?;
    let y4: [f64; 8] = std::array::from_fn(|i| y[i] + h * k3[i]);
    let k4 = deriv(bh, &y4)?;
    Ok(std::array::from_fn(|i| {
        y[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
    }))
}

/// One accepted adaptive step via step-doubling. Returns `(new_y, step_used,
/// next_step)`.
fn adaptive_step(bh: &KerrNewman, y: &[f64; 8], h0: f64, tol: f64) -> Result<([f64; 8], f64, f64)> {
    let mut h = h0;
    loop {
        let big = rk4(bh, y, h)?;
        let mid = rk4(bh, y, 0.5 * h)?;
        let two = rk4(bh, &mid, 0.5 * h)?;
        let err = (0..8).map(|i| (two[i] - big[i]).abs()).fold(0.0, f64::max);
        if err <= tol || h <= 1e-9 {
            // Richardson-extrapolated (5th-order) accepted value.
            let y_new: [f64; 8] = std::array::from_fn(|i| two[i] + (two[i] - big[i]) / 15.0);
            let grow = if err > 0.0 {
                (0.9 * (tol / err).powf(0.2)).clamp(0.2, 4.0)
            } else {
                4.0
            };
            return Ok((y_new, h, h * grow));
        }
        h *= (0.9 * (tol / err).powf(0.2)).clamp(0.1, 0.9);
    }
}

/// Integrate a geodesic from `init` until a stop condition in `opts` fires.
///
/// # Errors
/// [`RelativityError::CoordinateSingularity`] if the path reaches a point where
/// the metric breaks down (e.g. a horizon) before a clean capture/escape.
pub fn integrate_geodesic(
    bh: &KerrNewman,
    init: GeodesicState,
    opts: GeodesicOptions,
) -> Result<Trajectory> {
    let mut y = [
        init.x[0], init.x[1], init.x[2], init.x[3], init.u[0], init.u[1], init.u[2], init.u[3],
    ];
    let mut lambda = 0.0;
    let mut h = opts.step;
    let mut states = vec![init];
    let mut lambdas = vec![0.0];

    for _ in 0..opts.max_steps {
        let r = y[1];
        let ur = y[5];
        if r <= opts.r_capture {
            return Ok(Trajectory {
                states,
                lambda: lambdas,
                stop: StopReason::Captured,
            });
        }
        if r >= opts.r_escape && ur > 0.0 {
            return Ok(Trajectory {
                states,
                lambda: lambdas,
                stop: StopReason::Escaped,
            });
        }
        if lambda >= opts.max_lambda {
            return Ok(Trajectory {
                states,
                lambda: lambdas,
                stop: StopReason::Limit,
            });
        }
        let (yn, used, next) = adaptive_step(bh, &y, h, opts.tol)?;
        y = yn;
        lambda += used;
        h = next;
        states.push(GeodesicState {
            x: [y[0], y[1], y[2], y[3]],
            u: [y[4], y[5], y[6], y[7]],
        });
        lambdas.push(lambda);
    }
    Ok(Trajectory {
        states,
        lambda: lambdas,
        stop: StopReason::Limit,
    })
}

/// Weak-field light-deflection angle `4M/b` (radians) — the analytic reference.
pub fn light_deflection_weak_field(mass: f64, impact_b: f64) -> f64 {
    4.0 * mass / impact_b
}

/// Schwarzschild perihelion advance per orbit `6πM/p` (radians), with
/// semi-latus rectum `p = a(1 − e²)` — the analytic (weak-field) reference.
pub fn perihelion_advance_per_orbit(mass: f64, semilatus_p: f64) -> f64 {
    6.0 * PI * mass / semilatus_p
}

/// Numerically integrate a photon of impact parameter `impact_b` past a
/// non-rotating black hole and return its total deflection angle (radians),
/// with the flat-space baseline at the finite start radius removed. In the weak
/// field this approaches [`light_deflection_weak_field`] (`4M/b`).
///
/// # Errors
/// [`RelativityError::Unsupported`] for a rotating hole; [`RelativityError::
/// GeodesicNonConvergence`] if the photon is captured instead of escaping.
pub fn deflection_angle(bh: &KerrNewman, impact_b: f64) -> Result<f64> {
    if bh.spin != 0.0 {
        return Err(RelativityError::Unsupported(
            "deflection_angle is for non-rotating holes (use the general integrator otherwise)"
                .into(),
        ));
    }
    if bh.mass <= 0.0 || impact_b <= 0.0 {
        return Err(RelativityError::InvalidParameter(
            "mass and impact parameter must be positive".into(),
        ));
    }
    let r0 = (1.0e5_f64).max(100.0 * impact_b);
    let init = equatorial_state(bh, r0, 1.0, impact_b, Kind::Null, true)?;
    let r_capture = horizon_capture_radius(bh);
    let opts = GeodesicOptions {
        step: r0 / 100.0,
        tol: 1e-12,
        r_capture,
        r_escape: r0,
        max_lambda: 1e9,
        max_steps: 2_000_000,
    };
    let traj = integrate_geodesic(bh, init, opts)?;
    match traj.stop {
        StopReason::Escaped => {
            let phi = traj.last().x[3] - init.x[3];
            let flat_baseline = 2.0 * (impact_b / r0).acos();
            Ok(phi - flat_baseline)
        }
        StopReason::Captured => Err(RelativityError::GeodesicNonConvergence(
            "photon was captured (impact parameter below the critical value)".into(),
        )),
        StopReason::Limit => Err(RelativityError::GeodesicNonConvergence(
            "photon did not escape within the integration limits".into(),
        )),
    }
}

/// Numerically integrate one radial period of a bound, equatorial, *Schwarzschild*
/// orbit with the given perihelion and aphelion radii, and return the perihelion
/// advance per orbit (radians). Validates against [`perihelion_advance_per_orbit`]
/// (`6πM/p`) in the weak field.
///
/// # Errors
/// [`RelativityError::Unsupported`] if the hole spins or is charged;
/// [`RelativityError::InvalidParameter`] for bad radii;
/// [`RelativityError::GeodesicNonConvergence`] if no perihelion return is found.
pub fn orbit_precession(bh: &KerrNewman, r_peri: f64, r_apo: f64) -> Result<f64> {
    if bh.spin != 0.0 || bh.charge != 0.0 {
        return Err(RelativityError::Unsupported(
            "orbit_precession uses the Schwarzschild closed-form orbit setup".into(),
        ));
    }
    if bh.mass <= 0.0 || r_peri <= 0.0 || r_apo <= r_peri {
        return Err(RelativityError::InvalidParameter(
            "need 0 < r_peri < r_apo and positive mass".into(),
        ));
    }
    let m = bh.mass;
    let f = |r: f64| 1.0 - 2.0 * m / r;
    let (r1, r2) = (r_peri, r_apo);
    // Solve the two turning-point conditions E² = f(r)(1 + L²/r²) for E², L².
    let l2 = (f(r2) - f(r1)) / (f(r1) / (r1 * r1) - f(r2) / (r2 * r2));
    let e2 = f(r1) * (1.0 + l2 / (r1 * r1));
    if !l2.is_finite() || !e2.is_finite() || l2 <= 0.0 || e2 <= 0.0 {
        return Err(RelativityError::InvalidParameter(
            "no bound orbit for these radii".into(),
        ));
    }
    let init = equatorial_state(bh, r1, e2.sqrt(), l2.sqrt(), Kind::Timelike, false)?;
    let mut y = [
        init.x[0], init.x[1], init.x[2], init.x[3], init.u[0], init.u[1], init.u[2], init.u[3],
    ];
    let mut h = (r2 - r1) / 50.0;
    let tol = 1e-12;
    let mut lambda = 0.0;
    let max_lambda = 1e9;
    // First push it off the turning point (u^r becomes positive), then detect
    // the next perihelion: u^r crossing from negative back to positive.
    for _ in 0..5_000_000 {
        let prev = y;
        let (yn, used, next) = adaptive_step(bh, &y, h, tol)?;
        y = yn;
        h = next;
        lambda += used;
        if prev[5] < 0.0 && y[5] >= 0.0 {
            // Interpolate φ at the u^r = 0 crossing.
            let t = -prev[5] / (y[5] - prev[5]);
            let phi_cross = prev[3] + t * (y[3] - prev[3]);
            return Ok(phi_cross - 2.0 * PI);
        }
        if lambda > max_lambda {
            break;
        }
    }
    Err(RelativityError::GeodesicNonConvergence(
        "no perihelion return found within limits".into(),
    ))
}

/// A capture radius just outside the outer horizon (or `2M` if no horizon),
/// used as the plunge cutoff so integration stops before the coordinate
/// singularity.
fn horizon_capture_radius(bh: &KerrNewman) -> f64 {
    let disc = bh.horizon_discriminant();
    let rplus = if disc >= 0.0 {
        bh.mass + disc.sqrt()
    } else {
        2.0 * bh.mass
    };
    rplus * 1.0001 + 1e-6
}
