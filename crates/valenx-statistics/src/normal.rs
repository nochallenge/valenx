//! The standard normal distribution `N(0, 1)`.
//!
//! Provides the probability density [`pdf`], the cumulative distribution
//! [`cdf`], and the survival function [`sf`]. The density is the exact
//! closed form; the cumulative distribution is expressed through the error
//! function [`erf`], which is evaluated with the **Abramowitz & Stegun
//! formula 7.1.26** rational-times-exponential approximation (maximum
//! absolute error `≈ 1.5 × 10⁻⁷`).
//!
//! ## Model
//!
//! For the standard normal:
//!
//! - density `φ(z) = (1 / √(2π)) · exp(−z² / 2)`,
//! - distribution `Φ(z) = ½ · [1 + erf(z / √2)]`,
//! - survival `1 − Φ(z) = Φ(−z)` by symmetry.
//!
//! ## Honest scope
//!
//! The A&S 7.1.26 erf is a *seven-significant-figure* approximation, not a
//! correctly-rounded special function. It is ample for teaching, plots and
//! the inferential helpers in this crate, but it is not a substitute for a
//! certified statistics library in a clinical or production setting.

use std::f64::consts::{FRAC_1_SQRT_2, PI};

/// The Gauss error function `erf(x)` via Abramowitz & Stegun 7.1.26.
///
/// `erf(x) = (2 / √π) · ∫₀ˣ exp(−t²) dt`. The approximation evaluates a
/// degree-5 polynomial in `t = 1 / (1 + p·|x|)` times `exp(−x²)` for
/// `x >= 0`, then uses the odd symmetry `erf(−x) = −erf(x)`. The maximum
/// absolute error is about `1.5 × 10⁻⁷` across the whole real line.
///
/// ## Examples
///
/// ```
/// use valenx_statistics::normal::erf;
/// // erf is odd and pins the origin.
/// assert!(erf(0.0).abs() < 1e-12);
/// assert!((erf(-1.0) + erf(1.0)).abs() < 1e-12);
/// // Saturates to ±1 in the tails.
/// assert!((erf(5.0) - 1.0).abs() < 1e-6);
/// ```
pub fn erf(x: f64) -> f64 {
    // `erf` is exactly zero at the origin (it is an odd function). The
    // A&S 7.1.26 polynomial only sums to ~1 − 1e-9 at `t = 1`, so pin the
    // analytic value here rather than carry that approximation error — this
    // also makes `cdf(0) = 0.5` exact.
    if x == 0.0 {
        return 0.0;
    }

    // Coefficients of A&S 7.1.26.
    const A1: f64 = 0.254_829_592;
    const A2: f64 = -0.284_496_736;
    const A3: f64 = 1.421_413_741;
    const A4: f64 = -1.453_152_027;
    const A5: f64 = 1.061_405_429;
    const P: f64 = 0.327_591_1;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let ax = x.abs();
    let t = 1.0 / (1.0 + P * ax);
    // Horner evaluation of the degree-5 polynomial in `t`.
    let poly = ((((A5 * t + A4) * t + A3) * t + A2) * t + A1) * t;
    let y = 1.0 - poly * (-ax * ax).exp();
    sign * y
}

/// The complementary error function `erfc(x) = 1 − erf(x)`.
///
/// Provided for symmetry with [`erf`]; shares its `≈ 1.5 × 10⁻⁷` accuracy
/// and is what backs the survival function [`sf`].
///
/// ## Examples
///
/// ```
/// use valenx_statistics::normal::erfc;
/// assert!((erfc(0.0) - 1.0).abs() < 1e-12);
/// ```
pub fn erfc(x: f64) -> f64 {
    1.0 - erf(x)
}

/// The standard-normal probability density function `φ(z)`.
///
/// `φ(z) = (1 / √(2π)) · exp(−z² / 2)`. It is strictly positive, integrates
/// to one over the real line, peaks at `z = 0` with value `1 / √(2π)`, and
/// is **even**: `φ(−z) = φ(z)`.
///
/// ## Examples
///
/// ```
/// use valenx_statistics::normal::pdf;
/// // Peak height at the origin is 1/sqrt(2*pi).
/// let peak = 1.0 / (2.0 * std::f64::consts::PI).sqrt();
/// assert!((pdf(0.0) - peak).abs() < 1e-12);
/// // Symmetric about zero.
/// assert!((pdf(-1.3) - pdf(1.3)).abs() < 1e-12);
/// ```
pub fn pdf(z: f64) -> f64 {
    let inv_sqrt_2pi = 1.0 / (2.0 * PI).sqrt();
    inv_sqrt_2pi * (-0.5 * z * z).exp()
}

/// The standard-normal cumulative distribution function `Φ(z)`.
///
/// `Φ(z) = P(Z ≤ z) = ½ · [1 + erf(z / √2)]`. It is monotone increasing
/// from `0` to `1`, passes through `Φ(0) = ½`, and satisfies the reflection
/// identity `Φ(−z) = 1 − Φ(z)`. Accuracy follows the A&S [`erf`] used
/// underneath (`≈ 1.5 × 10⁻⁷`).
///
/// ## Examples
///
/// ```
/// use valenx_statistics::normal::cdf;
/// // The median of the standard normal is exactly one half.
/// assert!((cdf(0.0) - 0.5).abs() < 1e-12);
/// // ~68% of mass lies within one standard deviation.
/// let within_1sd = cdf(1.0) - cdf(-1.0);
/// assert!((within_1sd - 0.682_689_5).abs() < 1e-6);
/// ```
pub fn cdf(z: f64) -> f64 {
    0.5 * (1.0 + erf(z * FRAC_1_SQRT_2))
}

/// The standard-normal survival function `1 − Φ(z) = P(Z > z)`.
///
/// Computed as `½ · erfc(z / √2)`, equal by symmetry to `Φ(−z)`. Useful for
/// upper-tail probabilities and one-sided p-values.
///
/// ## Examples
///
/// ```
/// use valenx_statistics::normal::{cdf, sf};
/// // Survival and CDF partition the unit mass.
/// let z = 0.8;
/// assert!((sf(z) + cdf(z) - 1.0).abs() < 1e-12);
/// ```
pub fn sf(z: f64) -> f64 {
    0.5 * erfc(z * FRAC_1_SQRT_2)
}
