//! Black-hole shadow ray-tracer.
//!
//! Backward ray-tracing: for each pixel of a distant observer's image plane we
//! launch a null geodesic *inward* and follow it. If it plunges through the
//! horizon the pixel is part of the **shadow** (dark); if it turns around and
//! escapes, the pixel sees the sky behind the hole. The boundary between the
//! two is the photon ring.
//!
//! The photon's initial conserved quantities are set from the image-plane
//! ("celestial") coordinates `(α, β)` using the standard distant-observer
//! relations (Bardeen 1973):
//!
//! ```text
//! λ = L/E = −α sin θ_obs
//! η = Q/E² = β² + cos²θ_obs (α² − a²)
//! ```
//!
//! For Schwarzschild (`a = 0`) the shadow is a centred disk of radius
//! `b = 3√3 M`; for Kerr seen edge-on it is shifted and flattened on the side
//! co-rotating with the spin — the asymmetry the Event Horizon Telescope images
//! of M87* and Sgr A* probe.
//!
//! Geometrized units `G = c = 1`.

// Raising an index is a small contraction over the 4 spacetime indices;
// explicit index loops mirror the formula and touch several tensors at once.
#![allow(clippy::needless_range_loop)]

use crate::geodesics::{integrate_geodesic, GeodesicOptions, GeodesicState, StopReason};
use crate::metric::Spacetime;
use crate::spacetimes::KerrNewman;
use crate::tensor;
use crate::{RelativityError, Result};

/// A rendered shadow image: per-pixel "is this pixel inside the shadow?" plus
/// the radius at which each escaping ray was last seen (useful for shading the
/// sky). Row-major, `height × width`, with `+β` (image "up") in the first row.
#[derive(Clone, Debug)]
pub struct ShadowImage {
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
    /// Half-width of the image plane in units of `M`; the plane spans
    /// `[−half_extent, +half_extent]` on both `α` and `β` axes.
    pub half_extent: f64,
    /// `true` where the ray was captured by the hole (shadow), row-major.
    pub shadow: Vec<bool>,
}

impl ShadowImage {
    /// Whether pixel `(col, row)` lies in the shadow.
    pub fn is_shadow(&self, col: usize, row: usize) -> bool {
        self.shadow[row * self.width + col]
    }

    /// Fraction of pixels that fall inside the shadow.
    pub fn shadow_fraction(&self) -> f64 {
        let n = self.shadow.iter().filter(|&&s| s).count();
        n as f64 / (self.width * self.height) as f64
    }
}

/// Trace a single image-plane pixel `(alpha, beta)` for an observer at radius
/// `r_obs` and polar angle `theta_obs`, returning how the ray ended.
///
/// # Errors
/// [`RelativityError`] propagated from the integrator (e.g. a coordinate
/// singularity), or [`RelativityError::InvalidParameter`] for a bad observer.
pub fn trace_pixel(
    bh: &KerrNewman,
    r_obs: f64,
    theta_obs: f64,
    alpha: f64,
    beta: f64,
) -> Result<StopReason> {
    let init = observer_photon(bh, r_obs, theta_obs, alpha, beta)?;
    let r_capture = {
        let disc = bh.horizon_discriminant();
        let rplus = if disc >= 0.0 {
            bh.mass + disc.sqrt()
        } else {
            2.0 * bh.mass
        };
        rplus * 1.01 + 1e-6
    };
    let opts = GeodesicOptions {
        step: r_obs / 200.0,
        tol: 1e-9,
        r_capture,
        r_escape: r_obs * 1.001,
        max_lambda: 50.0 * r_obs,
        max_steps: 200_000,
    };
    Ok(integrate_geodesic(bh, init, opts)?.stop)
}

/// Build the inward-going photon state at the observer for image coordinates
/// `(alpha, beta)`.
fn observer_photon(
    bh: &KerrNewman,
    r_obs: f64,
    theta_obs: f64,
    alpha: f64,
    beta: f64,
) -> Result<GeodesicState> {
    if !r_obs.is_finite() || r_obs <= 3.0 * bh.mass {
        return Err(RelativityError::InvalidParameter(
            "observer radius must be well outside the hole".into(),
        ));
    }
    let a = bh.spin;
    let x = [0.0, r_obs, theta_obs, 0.0];
    let g = bh.metric::<f64>(x);
    let ginv = tensor::inverse(&g).ok_or_else(|| {
        RelativityError::CoordinateSingularity("degenerate metric at observer".into())
    })?;

    let sin0 = theta_obs.sin();
    let cos0 = theta_obs.cos();
    if sin0.abs() < 1e-9 {
        return Err(RelativityError::InvalidParameter(
            "observer on the polar axis (sin θ = 0) is unsupported".into(),
        ));
    }

    // Conserved quantities from celestial coordinates (E = 1).
    let lambda = -alpha * sin0; // L
    let eta = beta * beta + cos0 * cos0 * (alpha * alpha - a * a); // Carter Q

    // Covariant momentum components.
    let p_t = -1.0;
    let p_phi = lambda;
    let ptheta2 = eta - (lambda * lambda / (sin0 * sin0) - a * a) * cos0 * cos0;
    let theta_sign = if beta >= 0.0 { 1.0 } else { -1.0 };
    let p_theta = theta_sign * ptheta2.max(0.0).sqrt();

    // Null condition g^{μν} p_μ p_ν = 0 → solve for p_r (ingoing).
    let s = ginv[0][0] * p_t * p_t
        + 2.0 * ginv[0][3] * p_t * p_phi
        + ginv[3][3] * p_phi * p_phi
        + ginv[2][2] * p_theta * p_theta;
    let pr2 = -s / ginv[1][1];
    if !pr2.is_finite() || pr2 < 0.0 {
        return Err(RelativityError::InvalidParameter(format!(
            "no real inward photon for (α={alpha}, β={beta})"
        )));
    }
    let p_r = -pr2.sqrt(); // ingoing

    let p = [p_t, p_r, p_theta, p_phi];
    // Raise: u^μ = g^{μν} p_ν.
    let mut u = [0.0; 4];
    for mu in 0..4 {
        u[mu] = (0..4).map(|nu| ginv[mu][nu] * p[nu]).sum();
    }
    Ok(GeodesicState { x, u })
}

/// Render a full shadow image: trace `width × height` pixels spanning
/// `[−half_extent, half_extent]` (in `M`) on both axes, for an observer at
/// `(r_obs, theta_obs)`.
///
/// # Errors
/// Propagates integrator/observer errors from [`trace_pixel`].
pub fn render_shadow(
    bh: &KerrNewman,
    r_obs: f64,
    theta_obs: f64,
    half_extent: f64,
    width: usize,
    height: usize,
) -> Result<ShadowImage> {
    let mut shadow = vec![false; width * height];
    for row in 0..height {
        // +β at the top row.
        let beta = half_extent * (1.0 - 2.0 * (row as f64 + 0.5) / height as f64);
        for col in 0..width {
            let alpha = half_extent * (2.0 * (col as f64 + 0.5) / width as f64 - 1.0);
            let captured = matches!(
                trace_pixel(bh, r_obs, theta_obs, alpha, beta)?,
                StopReason::Captured
            );
            shadow[row * width + col] = captured;
        }
    }
    Ok(ShadowImage {
        width,
        height,
        half_extent,
        shadow,
    })
}

/// The shadow's left and right edges along the `β = 0` line of an edge-on
/// (`θ_obs = π/2`) image, found by bisection on `α`. For Schwarzschild these are
/// `∓3√3 M` (symmetric); for Kerr they are asymmetric (frame dragging shifts the
/// shadow toward the co-rotating side).
///
/// Returns `(alpha_left, alpha_right)` with `alpha_left < 0 < alpha_right`.
///
/// # Errors
/// Propagates ray-tracing errors.
pub fn equatorial_shadow_edges(bh: &KerrNewman, search_max: f64) -> Result<(f64, f64)> {
    let r_obs = 1000.0_f64.max(50.0 * bh.mass);
    let captured = |alpha: f64| -> Result<bool> {
        Ok(matches!(
            trace_pixel(bh, r_obs, std::f64::consts::FRAC_PI_2, alpha, 0.0)?,
            StopReason::Captured
        ))
    };
    // The shadow is a connected interval in α containing 0; bisect each edge
    // between a captured inner point and an escaping outer point.
    let edge = |sign: f64| -> Result<f64> {
        let (mut inside, mut outside) = (0.0_f64, sign * search_max);
        // Require the inner point captured and the outer point escaping.
        if !captured(inside)? {
            return Err(RelativityError::GeodesicNonConvergence(
                "image centre is not in shadow".into(),
            ));
        }
        for _ in 0..60 {
            let mid = 0.5 * (inside + outside);
            if captured(mid)? {
                inside = mid;
            } else {
                outside = mid;
            }
        }
        Ok(0.5 * (inside + outside))
    };
    Ok((edge(-1.0)?, edge(1.0)?))
}
