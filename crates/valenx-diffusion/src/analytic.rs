//! Closed-form 1-D diffusion: Fick's first law and the Gaussian
//! point-source (instantaneous plane-source) Green's function.
//!
//! ## Fick's first law
//!
//! The diffusive flux is proportional to the negative concentration
//! gradient,
//!
//! ```text
//!   J = -D dC/dx
//! ```
//!
//! so matter moves *down* the gradient: where `dC/dx > 0`, `J < 0`.
//! [`first_law_flux`] evaluates this for an explicit gradient;
//! [`flux_central`] estimates the gradient by a centred difference from
//! two samples a distance `dx` apart.
//!
//! ## Instantaneous point source
//!
//! Releasing a mass per unit area `M` at the origin at `t = 0` into an
//! infinite 1-D medium gives the heat-kernel / fundamental solution of
//! Fick's second law,
//!
//! ```text
//!   C(x, t) = M / sqrt(4 pi D t) * exp( -x^2 / (4 D t) )
//! ```
//!
//! a Gaussian in `x` whose integral over all `x` is the conserved `M`
//! for every `t > 0`, and whose spatial variance grows linearly,
//!
//! ```text
//!   var(t) = 2 D t        (so the RMS spread is sqrt(2 D t)).
//! ```
//!
//! [`gaussian_point_source`] evaluates `C(x, t)`,
//! [`gaussian_variance`] returns `2 D t`, and [`gaussian_std`] returns
//! its square root. The inverses [`time_to_reach_variance`] /
//! [`time_to_reach_std`] answer the design question — how long until the
//! spread reaches a prescribed width.

use crate::error::{DiffusionError, Result};
use std::f64::consts::PI;

/// Fick's first law: the diffusive flux for an explicit concentration
/// gradient.
///
/// Returns `J = -d * grad`, where `grad = dC/dx`. The flux points down
/// the gradient: a positive gradient yields a negative (leftward) flux.
///
/// # Errors
///
/// [`DiffusionError::BadParameter`] if `d` is not strictly positive, or
/// if either argument is non-finite.
///
/// # Examples
///
/// ```
/// use valenx_diffusion::first_law_flux;
///
/// // A positive gradient drives flux in the -x direction.
/// let j = first_law_flux(2.0, 3.0).unwrap();
/// assert!((j - (-6.0)).abs() < 1e-12);
/// ```
pub fn first_law_flux(d: f64, grad: f64) -> Result<f64> {
    check_diffusivity(d)?;
    if !grad.is_finite() {
        return Err(DiffusionError::bad_parameter("grad", "must be finite"));
    }
    Ok(-d * grad)
}

/// Fick's first law with the gradient estimated by a centred
/// difference between two concentration samples `dx` apart.
///
/// `c_left` is the concentration at `x - dx/2` (more precisely, the
/// lower-`x` sample) and `c_right` the concentration at the higher-`x`
/// sample; the points are separated by `dx`. The estimated gradient is
/// `(c_right - c_left) / dx` and the returned flux is `-d` times that.
///
/// # Errors
///
/// [`DiffusionError::BadParameter`] if `d` or `dx` is not strictly
/// positive, or if any argument is non-finite.
///
/// # Examples
///
/// ```
/// use valenx_diffusion::flux_central;
///
/// // C rises from 0 to 1 over dx = 0.5, so dC/dx = 2 and J = -D*2.
/// let j = flux_central(1.0, 0.0, 1.0, 0.5).unwrap();
/// assert!((j - (-2.0)).abs() < 1e-12);
/// ```
pub fn flux_central(d: f64, c_left: f64, c_right: f64, dx: f64) -> Result<f64> {
    check_diffusivity(d)?;
    check_spacing(dx)?;
    if !c_left.is_finite() || !c_right.is_finite() {
        return Err(DiffusionError::bad_parameter(
            "concentration",
            "samples must be finite",
        ));
    }
    let grad = (c_right - c_left) / dx;
    Ok(-d * grad)
}

/// The instantaneous point-source (plane-source) concentration at
/// position `x` and time `t`.
///
/// Evaluates
///
/// ```text
///   C(x, t) = mass / sqrt(4 pi D t) * exp( -x^2 / (4 D t) )
/// ```
///
/// the Green's function of Fick's second law for a mass per unit area
/// `mass` released at the origin at `t = 0`.
///
/// # Errors
///
/// [`DiffusionError::BadParameter`] if `d` is not strictly positive, if
/// `t` is not strictly positive (the kernel is singular at `t = 0`), or
/// if any argument is non-finite. `mass` may be any finite value
/// (including zero or negative, for superposing sources).
///
/// # Examples
///
/// ```
/// use valenx_diffusion::gaussian_point_source;
///
/// // Peak at the origin equals mass / sqrt(4 pi D t).
/// let peak = gaussian_point_source(1.0, 1.0, 0.0, 1.0).unwrap();
/// let expected = 1.0 / (4.0 * std::f64::consts::PI).sqrt();
/// assert!((peak - expected).abs() < 1e-12);
/// ```
pub fn gaussian_point_source(mass: f64, d: f64, x: f64, t: f64) -> Result<f64> {
    check_diffusivity(d)?;
    check_time_positive(t)?;
    if !mass.is_finite() {
        return Err(DiffusionError::bad_parameter("mass", "must be finite"));
    }
    if !x.is_finite() {
        return Err(DiffusionError::bad_parameter("x", "must be finite"));
    }
    let four_dt = 4.0 * d * t;
    let norm = mass / (PI * four_dt).sqrt();
    Ok(norm * (-x * x / four_dt).exp())
}

/// The spatial variance of the Gaussian point-source profile at time
/// `t`: `var = 2 D t`.
///
/// This is the exact second central moment of
/// [`gaussian_point_source`]; the spreading is purely a function of the
/// product `D t`.
///
/// # Errors
///
/// [`DiffusionError::BadParameter`] if `d` is not strictly positive, if
/// `t` is negative, or if either argument is non-finite. `t = 0` is
/// allowed and returns `0`.
///
/// # Examples
///
/// ```
/// use valenx_diffusion::gaussian_variance;
///
/// assert!((gaussian_variance(0.5, 4.0).unwrap() - 4.0).abs() < 1e-12);
/// ```
pub fn gaussian_variance(d: f64, t: f64) -> Result<f64> {
    check_diffusivity(d)?;
    check_time_nonneg(t)?;
    Ok(2.0 * d * t)
}

/// The standard deviation (RMS spread) of the Gaussian point-source
/// profile at time `t`: `sqrt(2 D t)`.
///
/// # Errors
///
/// Same conditions as [`gaussian_variance`].
///
/// # Examples
///
/// ```
/// use valenx_diffusion::gaussian_std;
///
/// // var = 2*0.5*4 = 4, so std = 2.
/// assert!((gaussian_std(0.5, 4.0).unwrap() - 2.0).abs() < 1e-12);
/// ```
pub fn gaussian_std(d: f64, t: f64) -> Result<f64> {
    Ok(gaussian_variance(d, t)?.sqrt())
}

/// The time for the point-source spreading to reach a given spatial
/// **variance**, inverting `var = 2 D t`:
///
/// ```text
///   t = target_variance / (2 D).
/// ```
///
/// The practical design question — how long until a tracer has spread to
/// a prescribed width — answered in closed form.
///
/// # Errors
///
/// [`DiffusionError::BadParameter`] if `d` is not strictly positive, or if
/// `target_variance` is negative or non-finite (`0` is allowed and returns
/// `0`).
///
/// # Examples
///
/// ```
/// use valenx_diffusion::{gaussian_variance, time_to_reach_variance};
///
/// // Reaching variance 4 with D = 0.5 takes t = 4 / (2*0.5) = 4.
/// let t = time_to_reach_variance(0.5, 4.0).unwrap();
/// assert!((t - 4.0).abs() < 1e-12);
/// // ... which the forward law reproduces.
/// assert!((gaussian_variance(0.5, t).unwrap() - 4.0).abs() < 1e-12);
/// ```
pub fn time_to_reach_variance(d: f64, target_variance: f64) -> Result<f64> {
    check_diffusivity(d)?;
    if !target_variance.is_finite() || target_variance < 0.0 {
        return Err(DiffusionError::bad_parameter(
            "target_variance",
            "must be finite and non-negative",
        ));
    }
    Ok(target_variance / (2.0 * d))
}

/// The time for the RMS spread (standard deviation) to reach `target_std`,
/// inverting `std = sqrt(2 D t)`:
///
/// ```text
///   t = target_std^2 / (2 D).
/// ```
///
/// Because the spread grows as the square root of time, the time scales as
/// the *square* of the target distance — diffusing twice as far takes four
/// times as long.
///
/// # Errors
///
/// [`DiffusionError::BadParameter`] if `d` is not strictly positive, or if
/// `target_std` is negative or non-finite (`0` is allowed and returns `0`).
///
/// # Examples
///
/// ```
/// use valenx_diffusion::{gaussian_std, time_to_reach_std};
///
/// // Reaching an RMS spread of 2 with D = 0.5 takes t = 4 / (2*0.5) = 4.
/// let t = time_to_reach_std(0.5, 2.0).unwrap();
/// assert!((t - 4.0).abs() < 1e-12);
/// assert!((gaussian_std(0.5, t).unwrap() - 2.0).abs() < 1e-12);
/// ```
pub fn time_to_reach_std(d: f64, target_std: f64) -> Result<f64> {
    check_diffusivity(d)?;
    if !target_std.is_finite() || target_std < 0.0 {
        return Err(DiffusionError::bad_parameter(
            "target_std",
            "must be finite and non-negative",
        ));
    }
    Ok(target_std * target_std / (2.0 * d))
}

/// The distance from the source at which the instantaneous point-source
/// profile has fallen to a given **fraction of its peak** — the spatial
/// inverse of [`gaussian_point_source`].
///
/// Since `C(x, t) / C(0, t) = exp(-x^2 / (4 D t))`, setting that ratio to
/// `fraction` and solving for the non-negative distance gives
///
/// ```text
///   x = sqrt( 4 D t * ln(1 / fraction) ).
/// ```
///
/// It answers "how far out has the tracer reached the `fraction` contour?"
/// At `fraction = 1` it is `0` (the peak sits at the origin); at the
/// special value `fraction = exp(-1/2)` it returns exactly the RMS spread
/// [`gaussian_std`] (`sqrt(2 D t)`); and as `fraction -> 0` it grows
/// without bound. By symmetry the same concentration also occurs at `-x`.
/// This complements the temporal inverses [`time_to_reach_variance`] /
/// [`time_to_reach_std`] with a spatial one.
///
/// # Errors
///
/// [`DiffusionError::BadParameter`] if `d` is not strictly positive, if
/// `t` is not strictly positive, or if `fraction` is not finite and in
/// `(0, 1]`.
///
/// # Examples
///
/// ```
/// use valenx_diffusion::{distance_for_concentration_fraction, gaussian_std};
///
/// // The exp(-1/2) contour sits exactly one standard deviation out.
/// let (d, t) = (1.5, 2.0);
/// let x = distance_for_concentration_fraction(d, t, (-0.5_f64).exp()).unwrap();
/// assert!((x - gaussian_std(d, t).unwrap()).abs() < 1e-12);
/// ```
pub fn distance_for_concentration_fraction(d: f64, t: f64, fraction: f64) -> Result<f64> {
    check_diffusivity(d)?;
    check_time_positive(t)?;
    if !fraction.is_finite() || fraction <= 0.0 || fraction > 1.0 {
        return Err(DiffusionError::bad_parameter(
            "fraction",
            "must be finite and in (0, 1]",
        ));
    }
    let ln_inv = (1.0 / fraction).ln();
    Ok((4.0 * d * t * ln_inv).sqrt())
}

// --- internal validators -------------------------------------------------

/// Reject a non-positive or non-finite diffusion coefficient.
fn check_diffusivity(d: f64) -> Result<()> {
    if !d.is_finite() || d <= 0.0 {
        return Err(DiffusionError::bad_parameter(
            "D",
            "diffusion coefficient must be finite and strictly positive",
        ));
    }
    Ok(())
}

/// Reject a non-positive or non-finite cell spacing.
fn check_spacing(dx: f64) -> Result<()> {
    if !dx.is_finite() || dx <= 0.0 {
        return Err(DiffusionError::bad_parameter(
            "dx",
            "cell spacing must be finite and strictly positive",
        ));
    }
    Ok(())
}

/// Reject a non-positive or non-finite time (for the singular kernel).
fn check_time_positive(t: f64) -> Result<()> {
    if !t.is_finite() || t <= 0.0 {
        return Err(DiffusionError::bad_parameter(
            "t",
            "time must be finite and strictly positive (the point-source kernel is singular at t = 0)",
        ));
    }
    Ok(())
}

/// Reject a negative or non-finite time (where `t = 0` is admissible).
fn check_time_nonneg(t: f64) -> Result<()> {
    if !t.is_finite() || t < 0.0 {
        return Err(DiffusionError::bad_parameter(
            "t",
            "time must be finite and non-negative",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn flux_opposes_gradient() {
        // Positive gradient -> negative (down-gradient) flux.
        let j = first_law_flux(2.5, 4.0).unwrap();
        assert!(j < 0.0);
        assert!((j - (-10.0)).abs() < EPS, "got {j}");

        // Negative gradient -> positive flux.
        let j = first_law_flux(2.5, -4.0).unwrap();
        assert!(j > 0.0);
        assert!((j - 10.0).abs() < EPS, "got {j}");

        // Zero gradient -> zero flux.
        let j = first_law_flux(2.5, 0.0).unwrap();
        assert!(j.abs() < EPS, "got {j}");
    }

    #[test]
    fn central_flux_matches_explicit_gradient() {
        // C goes 1 -> 3 over dx = 0.5: gradient = 4, flux = -D*4.
        let d = 0.75;
        let j = flux_central(d, 1.0, 3.0, 0.5).unwrap();
        assert!((j - (-d * 4.0)).abs() < EPS, "got {j}");
    }

    #[test]
    fn point_source_peak_value() {
        // At x = 0 the exponential is 1, so C = M / sqrt(4 pi D t).
        let mass = 5.0;
        let d = 2.0;
        let t = 3.0;
        let c = gaussian_point_source(mass, d, 0.0, t).unwrap();
        let expected = mass / (4.0 * PI * d * t).sqrt();
        assert!((c - expected).abs() < EPS, "got {c}, want {expected}");
    }

    #[test]
    fn point_source_is_symmetric_and_decays() {
        let (mass, d, t) = (1.0, 1.0, 0.5);
        let c_pos = gaussian_point_source(mass, d, 0.7, t).unwrap();
        let c_neg = gaussian_point_source(mass, d, -0.7, t).unwrap();
        // Even in x.
        assert!((c_pos - c_neg).abs() < EPS, "{c_pos} vs {c_neg}");
        // Strictly below the peak.
        let peak = gaussian_point_source(mass, d, 0.0, t).unwrap();
        assert!(c_pos < peak);
    }

    #[test]
    fn point_source_value_at_one_sigma() {
        // At |x| = sigma = sqrt(2 D t) the profile is exp(-1/2) of peak.
        let (mass, d, t): (f64, f64, f64) = (1.0, 1.5, 2.0);
        let sigma = (2.0 * d * t).sqrt();
        let peak = gaussian_point_source(mass, d, 0.0, t).unwrap();
        let at_sigma = gaussian_point_source(mass, d, sigma, t).unwrap();
        let ratio = at_sigma / peak;
        assert!((ratio - (-0.5_f64).exp()).abs() < 1e-10, "ratio {ratio}");
    }

    #[test]
    fn variance_grows_as_two_d_t() {
        let d = 0.3;
        for &t in &[0.0, 1.0, 5.0, 12.5] {
            let v = gaussian_variance(d, t).unwrap();
            assert!((v - 2.0 * d * t).abs() < EPS, "t={t} got {v}");
        }
        // Doubling t doubles the variance.
        let v1 = gaussian_variance(d, 4.0).unwrap();
        let v2 = gaussian_variance(d, 8.0).unwrap();
        assert!((v2 - 2.0 * v1).abs() < EPS, "{v1} {v2}");
    }

    #[test]
    fn variance_recovered_by_numeric_moment() {
        // Numerically integrate x^2 * C(x) / integral(C) and check it
        // matches the analytic 2 D t. Trapezoid over a wide window.
        let (mass, d, t): (f64, f64, f64) = (1.0, 1.0, 1.0);
        let sigma = (2.0 * d * t).sqrt();
        let half = 8.0 * sigma;
        let n = 20_001;
        let dx = 2.0 * half / (n as f64 - 1.0);
        let mut m0 = 0.0;
        let mut m2 = 0.0;
        for i in 0..n {
            let x = -half + dx * i as f64;
            let c = gaussian_point_source(mass, d, x, t).unwrap();
            let w = if i == 0 || i == n - 1 { 0.5 } else { 1.0 };
            m0 += w * c * dx;
            m2 += w * x * x * c * dx;
        }
        // Zeroth moment recovers the conserved mass.
        assert!((m0 - mass).abs() < 1e-4, "m0 = {m0}");
        let var = m2 / m0;
        assert!((var - 2.0 * d * t).abs() < 1e-3, "var = {var}");
    }

    #[test]
    fn std_is_sqrt_of_variance() {
        let (d, t) = (2.0, 3.0);
        let s = gaussian_std(d, t).unwrap();
        let v = gaussian_variance(d, t).unwrap();
        assert!((s * s - v).abs() < EPS, "s={s} v={v}");
    }

    #[test]
    fn time_to_reach_variance_inverts_the_forward_law() {
        let d = 0.3;
        // Forward then inverse recovers the time.
        for &t in &[0.0, 1.0, 5.0, 12.5] {
            let v = gaussian_variance(d, t).unwrap();
            let back = time_to_reach_variance(d, v).unwrap();
            assert!((back - t).abs() < EPS, "t={t} -> v={v} -> {back}");
        }
        // Inverse then forward recovers the variance.
        for &v in &[0.0, 2.0, 9.0] {
            let t = time_to_reach_variance(d, v).unwrap();
            assert!((gaussian_variance(d, t).unwrap() - v).abs() < EPS, "v={v}");
        }
    }

    #[test]
    fn time_to_reach_std_inverts_and_scales_with_distance_squared() {
        let d = 0.5;
        // Inverts the std law: D=0.5, sigma=2 -> t = 4/(2*0.5) = 4.
        let t = time_to_reach_std(d, 2.0).unwrap();
        assert!((t - 4.0).abs() < EPS, "t = {t}");
        assert!((gaussian_std(d, t).unwrap() - 2.0).abs() < EPS);
        // Diffusing twice as far takes four times as long.
        let t1 = time_to_reach_std(d, 1.0).unwrap();
        let t2 = time_to_reach_std(d, 2.0).unwrap();
        assert!((t2 - 4.0 * t1).abs() < EPS, "{t1} {t2}");
        // Consistent with the variance form (var = std^2).
        let sigma = 1.7;
        assert!(
            (time_to_reach_std(d, sigma).unwrap()
                - time_to_reach_variance(d, sigma * sigma).unwrap())
            .abs()
                < EPS
        );
    }

    #[test]
    fn inverse_time_rejects_bad_parameters() {
        assert!(time_to_reach_variance(0.0, 1.0).is_err());
        assert!(time_to_reach_variance(1.0, -1.0).is_err());
        assert!(time_to_reach_std(-1.0, 1.0).is_err());
        assert!(time_to_reach_std(1.0, f64::NAN).is_err());
        // Zero target -> zero time.
        assert!(time_to_reach_variance(1.0, 0.0).unwrap().abs() < EPS);
        assert!(time_to_reach_std(1.0, 0.0).unwrap().abs() < EPS);
    }

    #[test]
    fn distance_for_fraction_round_trips_through_the_profile() {
        // The distance to a fraction f, fed back into the profile, recovers
        // exactly that fraction of the peak.
        let (mass, d, t) = (3.0, 1.2, 2.5);
        let peak = gaussian_point_source(mass, d, 0.0, t).unwrap();
        for &f in &[1.0_f64, 0.8, 0.5, 0.1, 0.01] {
            let x = distance_for_concentration_fraction(d, t, f).unwrap();
            let ratio = gaussian_point_source(mass, d, x, t).unwrap() / peak;
            assert!((ratio - f).abs() < 1e-12, "f={f} -> x={x} -> ratio={ratio}");
        }
    }

    #[test]
    fn distance_for_one_over_sqrt_e_fraction_is_one_sigma() {
        // The exp(-1/2) contour is exactly one standard deviation out.
        let (d, t) = (1.5, 2.0);
        let x = distance_for_concentration_fraction(d, t, (-0.5_f64).exp()).unwrap();
        assert!((x - gaussian_std(d, t).unwrap()).abs() < EPS, "x={x}");
    }

    #[test]
    fn distance_for_full_fraction_is_zero_and_grows_as_fraction_falls() {
        let (d, t) = (0.7, 3.0);
        // The peak (fraction 1) sits at the origin.
        assert!(
            distance_for_concentration_fraction(d, t, 1.0)
                .unwrap()
                .abs()
                < EPS
        );
        // A lower target fraction lies farther out.
        let near = distance_for_concentration_fraction(d, t, 0.5).unwrap();
        let far = distance_for_concentration_fraction(d, t, 0.1).unwrap();
        assert!(far > near && near > 0.0, "near={near} far={far}");
    }

    #[test]
    fn distance_for_fraction_rejects_bad_parameters() {
        assert!(distance_for_concentration_fraction(0.0, 1.0, 0.5).is_err()); // D
        assert!(distance_for_concentration_fraction(1.0, 0.0, 0.5).is_err()); // t = 0
        assert!(distance_for_concentration_fraction(1.0, 1.0, 0.0).is_err()); // fraction 0
        assert!(distance_for_concentration_fraction(1.0, 1.0, 1.5).is_err()); // > 1
        assert!(distance_for_concentration_fraction(1.0, 1.0, f64::NAN).is_err());
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(first_law_flux(0.0, 1.0).is_err());
        assert!(first_law_flux(-1.0, 1.0).is_err());
        assert!(first_law_flux(1.0, f64::NAN).is_err());
        assert!(gaussian_point_source(1.0, 1.0, 0.0, 0.0).is_err());
        assert!(gaussian_point_source(1.0, 1.0, 0.0, -1.0).is_err());
        assert!(gaussian_variance(1.0, -1.0).is_err());
        assert!(flux_central(1.0, 0.0, 1.0, 0.0).is_err());
        // t = 0 is allowed for the variance and returns 0.
        assert!((gaussian_variance(1.0, 0.0).unwrap()).abs() < EPS);
    }
}
