//! Exterior free-field radiation: the pulsating sphere (monopole) and the
//! point dipole.
//!
//! Where the rest of this crate is *interior* / statistical (room modes,
//! Sabine reverberation) or single-point (SPL, Doppler), this module covers
//! the complementary case of **sound radiated into free space** by a compact
//! source — the everyday acoustic-radiation textbook results.
//!
//! ## Model
//!
//! ### Pulsating sphere (monopole)
//!
//! A sphere of radius `a` whose surface oscillates radially with velocity
//! amplitude `U` (peak), immersed in a fluid of density `rho` and sound
//! speed `c`, radiates a spherically-spreading wave. With wavenumber
//! `k = omega / c = 2*pi*f / c`, the **peak** acoustic pressure amplitude at
//! a radius `r >= a` is
//!
//! ```text
//! |p(r)| = rho * c * k * a^2 * U / (r * sqrt(1 + (k a)^2))      [Pa, peak]
//! ```
//!
//! (Kinsler & Frey, *Fundamentals of Acoustics*; Pierce, *Acoustics*.)
//! The crucial structural fact is the explicit `1 / r`: outside the sphere
//! the pressure falls off as the inverse of the distance, so **doubling the
//! distance drops the level by exactly `20*log10(2) ≈ 6.02` dB** — the
//! defining signature of a point (monopole) source in free field. The
//! `(k a)^2 / (1 + (k a)^2)` factor is the radiation efficiency: a sphere
//! that is small compared with the wavelength (`k a << 1`) is a poor
//! radiator, while `k a >> 1` approaches the planar limit.
//!
//! The time-averaged **radiated acoustic power**, obtained by integrating
//! the far-field intensity over a great sphere, is
//!
//! ```text
//! W = 2*pi * rho * c * a^2 * U^2 * (k a)^2 / (1 + (k a)^2)       [W]
//! ```
//!
//! with `U` the **peak** surface velocity (equivalently
//! `4*pi * rho * c * a^2 * U_rms^2 * (k a)^2 / (1 + (k a)^2)` with
//! `U_rms = U / sqrt(2)`). The corresponding **far-field intensity**
//! `I(r) = |p_rms(r)|^2 / (rho c) = W / (4*pi*r^2)` then recovers the
//! inverse-square spreading of power through concentric spheres, and the two
//! formulas above are mutually consistent through that identity (verified in
//! the module tests).
//!
//! ### Point dipole (directivity)
//!
//! Two equal out-of-phase monopoles separated by a small distance form a
//! dipole: a figure-of-eight source whose far-field pressure magnitude
//! carries a `cos(theta)` directivity, `theta` measured from the dipole
//! axis. Only the **directivity factor** `|cos(theta)|` is modelled here
//! (normalised to `1` on axis): broadside (`theta = 90 degrees`) is a null,
//! and the level relative to the on-axis maximum is
//! `20*log10(|cos theta|)` dB.
//!
//! ## Honest scope
//!
//! Research / educational, closed-form, single-frequency (time-harmonic)
//! magnitudes only — no phase, no transient build-up, no multi-source
//! interference field, no baffle / half-space correction, and no
//! near-field reactive detail beyond the `1 + (k a)^2` factor already in the
//! sphere formula. The dipole helper returns a *relative* directivity, not
//! an absolute pressure. The medium is a lossless homogeneous fluid (no
//! atmospheric absorption). Use it to reason about spreading, radiation
//! efficiency and directivity — not to certify a measurement.

use crate::error::{AcousticsError, Result};

/// Reference density of air at `20` degrees Celsius and one atmosphere, in
/// kilograms per cubic metre — a convenient default for the `rho` argument
/// of [`monopole_pressure`] and [`monopole_radiated_power`].
pub const RHO_AIR: f64 = 1.204;

fn check_positive_finite(name: &'static str, value: f64) -> Result<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(AcousticsError::InvalidDimension { name, value })
    }
}

fn check_frequency(name: &'static str, value: f64) -> Result<()> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(AcousticsError::InvalidFrequency { name, value })
    }
}

fn check_speed_of_sound(name: &'static str, value: f64) -> Result<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(AcousticsError::InvalidSpeedOfSound { name, value })
    }
}

/// Acoustic wavenumber `k = 2*pi*f / c` (radians per metre) for frequency
/// `frequency_hz` in a medium of sound speed `speed_of_sound`.
///
/// # Errors
///
/// - [`AcousticsError::InvalidFrequency`] if `frequency_hz` is negative or
///   non-finite.
/// - [`AcousticsError::InvalidSpeedOfSound`] if `speed_of_sound` is not
///   finite and strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_acoustics::radiation::wavenumber;
/// // At 343 m/s a 343 Hz tone has a 1 m wavelength, so k = 2*pi.
/// let k = wavenumber(343.0, 343.0).unwrap();
/// assert!((k - 2.0 * std::f64::consts::PI).abs() < 1e-9);
/// ```
pub fn wavenumber(frequency_hz: f64, speed_of_sound: f64) -> Result<f64> {
    check_frequency("frequency_hz", frequency_hz)?;
    check_speed_of_sound("speed_of_sound", speed_of_sound)?;
    Ok(2.0 * core::f64::consts::PI * frequency_hz / speed_of_sound)
}

/// Peak acoustic pressure amplitude (pascals) radiated by a pulsating
/// sphere of radius `radius_a` with surface velocity amplitude
/// `surface_velocity` (peak), evaluated at a field radius `field_r >= a`.
///
/// `|p(r)| = rho * c * k * a^2 * U / (r * sqrt(1 + (k a)^2))`, with
/// `k = 2*pi*f / c` the wavenumber.
///
/// The `1 / r` dependence means the level falls by `≈ 6.02` dB per distance
/// doubling — the analytic free-field monopole law validated in this
/// module's tests.
///
/// # Errors
///
/// - [`AcousticsError::InvalidDimension`] if `radius_a`, `field_r`,
///   `density`, or `surface_velocity` is non-positive or non-finite, **or**
///   if `field_r < radius_a` (the field point would be inside the sphere,
///   where this exterior solution does not apply).
/// - [`AcousticsError::InvalidFrequency`] if `frequency_hz` is negative or
///   non-finite.
/// - [`AcousticsError::InvalidSpeedOfSound`] if `speed_of_sound` is not
///   finite and strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_acoustics::radiation::monopole_pressure;
/// // Pressure at 2 m is half that at 1 m (1/r spreading).
/// let p1 = monopole_pressure(0.1, 0.01, 1000.0, 343.0, 1.204, 1.0).unwrap();
/// let p2 = monopole_pressure(0.1, 0.01, 1000.0, 343.0, 1.204, 2.0).unwrap();
/// assert!((p2 - 0.5 * p1).abs() < 1e-9 * p1);
/// ```
#[allow(clippy::too_many_arguments)]
pub fn monopole_pressure(
    radius_a: f64,
    surface_velocity: f64,
    frequency_hz: f64,
    speed_of_sound: f64,
    density: f64,
    field_r: f64,
) -> Result<f64> {
    check_positive_finite("radius_a", radius_a)?;
    check_positive_finite("surface_velocity", surface_velocity)?;
    check_positive_finite("density", density)?;
    check_positive_finite("field_r", field_r)?;
    check_frequency("frequency_hz", frequency_hz)?;
    check_speed_of_sound("speed_of_sound", speed_of_sound)?;
    if field_r < radius_a {
        return Err(AcousticsError::InvalidDimension {
            name: "field_r",
            value: field_r,
        });
    }
    let k = 2.0 * core::f64::consts::PI * frequency_hz / speed_of_sound;
    let ka = k * radius_a;
    Ok(
        density * speed_of_sound * k * radius_a * radius_a * surface_velocity
            / (field_r * (1.0 + ka * ka).sqrt()),
    )
}

/// Time-averaged acoustic power (watts) radiated by a pulsating sphere of
/// radius `radius_a` with surface velocity amplitude `surface_velocity`
/// (peak).
///
/// `W = 2*pi * rho * c * a^2 * U^2 * (k a)^2 / (1 + (k a)^2)`, with `U` the
/// peak surface velocity and `k = 2*pi*f / c` (equivalently the
/// `U_rms = U / sqrt(2)` form `4*pi * rho * c * a^2 * U_rms^2 * ...`).
///
/// The trailing `(k a)^2 / (1 + (k a)^2)` is the sphere's radiation
/// efficiency: it rises monotonically from `≈ (k a)^2` for a compact source
/// (`k a << 1`) toward `1` in the short-wavelength limit (`k a >> 1`).
///
/// # Errors
///
/// Same domain checks as [`monopole_pressure`] (minus the field-radius
/// argument, which this power integral does not take): non-positive or
/// non-finite `radius_a` / `surface_velocity` / `density` give
/// [`AcousticsError::InvalidDimension`]; a bad `frequency_hz` gives
/// [`AcousticsError::InvalidFrequency`]; a bad `speed_of_sound` gives
/// [`AcousticsError::InvalidSpeedOfSound`].
///
/// # Examples
///
/// ```
/// use valenx_acoustics::radiation::monopole_radiated_power;
/// let w = monopole_radiated_power(0.05, 0.01, 500.0, 343.0, 1.204).unwrap();
/// assert!(w > 0.0);
/// ```
pub fn monopole_radiated_power(
    radius_a: f64,
    surface_velocity: f64,
    frequency_hz: f64,
    speed_of_sound: f64,
    density: f64,
) -> Result<f64> {
    check_positive_finite("radius_a", radius_a)?;
    check_positive_finite("surface_velocity", surface_velocity)?;
    check_positive_finite("density", density)?;
    check_frequency("frequency_hz", frequency_hz)?;
    check_speed_of_sound("speed_of_sound", speed_of_sound)?;
    let k = 2.0 * core::f64::consts::PI * frequency_hz / speed_of_sound;
    let ka = k * radius_a;
    // W = 2*pi*rho*c*a^2 * U^2 * (ka)^2 / (1 + (ka)^2), U peak.
    Ok(2.0
        * core::f64::consts::PI
        * density
        * speed_of_sound
        * radius_a
        * radius_a
        * surface_velocity
        * surface_velocity
        * (ka * ka)
        / (1.0 + ka * ka))
}

/// Far-field relative directivity `|cos(theta)|` of a point dipole, `theta`
/// (radians) measured from the dipole axis.
///
/// Normalised to `1.0` on axis (`theta = 0`). Broadside (`theta = pi/2`) is
/// the dipole null, returning `0.0`. Combine with [`directivity_db`] to get
/// the on-axis-relative level in decibels.
///
/// # Errors
///
/// [`AcousticsError::InvalidVelocity`] (reused as the "non-finite angle"
/// signal) if `theta_rad` is not finite.
///
/// # Examples
///
/// ```
/// use valenx_acoustics::radiation::dipole_directivity;
/// use std::f64::consts::PI;
/// assert!((dipole_directivity(0.0).unwrap() - 1.0).abs() < 1e-12); // on axis
/// assert!(dipole_directivity(PI / 2.0).unwrap() < 1e-12);          // null
/// ```
pub fn dipole_directivity(theta_rad: f64) -> Result<f64> {
    if !theta_rad.is_finite() {
        return Err(AcousticsError::InvalidVelocity {
            name: "theta_rad",
            value: theta_rad,
        });
    }
    Ok(theta_rad.cos().abs())
}

/// Convert a linear far-field directivity ratio `d` (e.g. the value from
/// [`dipole_directivity`]) into a level relative to the on-axis maximum, in
/// decibels: `20*log10(d)`.
///
/// A ratio of `1.0` is `0` dB (on axis); a ratio of `0.5` is `≈ -6.02` dB.
///
/// # Errors
///
/// [`AcousticsError::InvalidPressure`] if `ratio` is non-positive or
/// non-finite (the logarithm of a non-positive ratio is undefined). The
/// exact broadside null (`ratio = 0`) therefore reports an error rather than
/// returning `-inf`.
///
/// # Examples
///
/// ```
/// use valenx_acoustics::radiation::directivity_db;
/// assert!((directivity_db(1.0).unwrap()).abs() < 1e-12);
/// assert!((directivity_db(0.5).unwrap() + 6.0206).abs() < 1e-3);
/// ```
pub fn directivity_db(ratio: f64) -> Result<f64> {
    if !ratio.is_finite() || ratio <= 0.0 {
        return Err(AcousticsError::InvalidPressure {
            name: "ratio",
            value: ratio,
        });
    }
    Ok(20.0 * ratio.log10())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const C: f64 = 343.0;

    /// Wavenumber matches `k = 2*pi/lambda` for a chosen wavelength.
    #[test]
    fn wavenumber_matches_wavelength() {
        // 1 m wavelength at 343 m/s -> 343 Hz -> k = 2*pi.
        let k = wavenumber(343.0, C).unwrap();
        assert!((k - 2.0 * PI).abs() < 1e-9, "got {k}");
        // Zero frequency -> zero wavenumber.
        assert!((wavenumber(0.0, C).unwrap()).abs() < 1e-15);
    }

    /// NAMED ANALYTIC VALIDATION — free-field monopole `1/r` spreading:
    /// the pressure magnitude halves per distance doubling, i.e. exactly
    /// `-6.02` dB per doubling, independent of frequency and sphere size.
    #[test]
    fn monopole_obeys_inverse_r_minus_six_db_per_doubling() {
        let a = 0.1;
        let u = 0.01;
        for &f in &[100.0f64, 1000.0, 5000.0] {
            for &r in &[1.0f64, 2.5, 10.0] {
                let p_r = monopole_pressure(a, u, f, C, RHO_AIR, r).unwrap();
                let p_2r = monopole_pressure(a, u, f, C, RHO_AIR, 2.0 * r).unwrap();
                // Pressure ratio is exactly 1/2.
                assert!(
                    (p_2r / p_r - 0.5).abs() < 1e-12,
                    "f={f} r={r}: ratio {}",
                    p_2r / p_r
                );
                // ... which is -6.0206 dB.
                let drop_db = 20.0 * (p_r / p_2r).log10();
                assert!(
                    (drop_db - 6.020_599_913).abs() < 1e-6,
                    "f={f} r={r}: drop {drop_db} dB"
                );
            }
        }
    }

    /// The exact pressure formula reproduces an independent hand evaluation.
    #[test]
    fn monopole_pressure_matches_closed_form() {
        let (a, u, f, rho, r) = (0.05, 0.02, 800.0, RHO_AIR, 3.0);
        let k = 2.0 * PI * f / C;
        let ka = k * a;
        let expected = rho * C * k * a * a * u / (r * (1.0 + ka * ka).sqrt());
        let got = monopole_pressure(a, u, f, C, rho, r).unwrap();
        assert!((got - expected).abs() < 1e-12 * expected, "got {got}");
    }

    /// Far-field intensity from the radiated power matches `|p_rms|^2/(rho c)`:
    /// `W/(4 pi r^2) == p_rms^2 / (rho c)` in the `r >> a`, `k a` regime
    /// where the simple `1/r` field dominates. We check the power and the
    /// pressure are mutually consistent through the intensity identity.
    #[test]
    fn power_and_pressure_are_intensity_consistent() {
        let (a, u, f, rho) = (0.02, 0.03, 2000.0, RHO_AIR);
        let r = 50.0; // deep far field
        let w = monopole_radiated_power(a, u, f, C, rho).unwrap();
        let p_peak = monopole_pressure(a, u, f, C, rho, r).unwrap();
        let p_rms_sq = p_peak * p_peak / 2.0;
        let i_from_pressure = p_rms_sq / (rho * C);
        let i_from_power = w / (4.0 * PI * r * r);
        let rel = (i_from_pressure - i_from_power).abs() / i_from_power;
        assert!(rel < 1e-9, "intensity mismatch rel={rel}");
    }

    /// Radiation efficiency rises with `k a`: a larger / higher-frequency
    /// sphere radiates more power for the same surface velocity.
    #[test]
    fn radiated_power_increases_with_ka() {
        let u = 0.01;
        let w_low = monopole_radiated_power(0.01, u, 100.0, C, RHO_AIR).unwrap();
        let w_high = monopole_radiated_power(0.10, u, 4000.0, C, RHO_AIR).unwrap();
        assert!(w_high > w_low, "{w_high} !> {w_low}");
        assert!(w_low > 0.0);
    }

    /// Dipole directivity: unity on axis, zero broadside, and `|cos|`
    /// symmetric front/back.
    #[test]
    fn dipole_directivity_pattern() {
        assert!((dipole_directivity(0.0).unwrap() - 1.0).abs() < 1e-12);
        assert!(dipole_directivity(PI / 2.0).unwrap() < 1e-12);
        assert!((dipole_directivity(PI).unwrap() - 1.0).abs() < 1e-12);
        // 60 degrees from axis -> cos 60 = 0.5.
        assert!((dipole_directivity(PI / 3.0).unwrap() - 0.5).abs() < 1e-12);
    }

    /// Directivity in dB: 0 dB on axis, -6.02 dB at the half-amplitude angle.
    #[test]
    fn directivity_db_values() {
        assert!(directivity_db(1.0).unwrap().abs() < 1e-12);
        assert!((directivity_db(0.5).unwrap() + 6.020_599_913).abs() < 1e-6);
    }

    // ---- bad-input: every domain guard must fail loud -------------------

    #[test]
    fn negative_frequency_fails_loud() {
        assert!(wavenumber(-1.0, C).is_err());
        assert!(monopole_pressure(0.1, 0.01, -10.0, C, RHO_AIR, 1.0).is_err());
        assert!(monopole_radiated_power(0.1, 0.01, -10.0, C, RHO_AIR).is_err());
    }

    #[test]
    fn nonpositive_radius_fails_loud() {
        assert!(monopole_pressure(0.0, 0.01, 1000.0, C, RHO_AIR, 1.0).is_err());
        assert!(monopole_pressure(-0.1, 0.01, 1000.0, C, RHO_AIR, 1.0).is_err());
        assert!(monopole_radiated_power(0.0, 0.01, 1000.0, C, RHO_AIR).is_err());
    }

    #[test]
    fn field_inside_sphere_fails_loud() {
        // r < a is non-physical for the exterior solution.
        let err = monopole_pressure(1.0, 0.01, 1000.0, C, RHO_AIR, 0.5).unwrap_err();
        assert_eq!(err.code(), "acoustics.invalid_dimension");
        // r == a is allowed (on the surface).
        assert!(monopole_pressure(1.0, 0.01, 1000.0, C, RHO_AIR, 1.0).is_ok());
    }

    #[test]
    fn bad_medium_and_velocity_fail_loud() {
        assert!(monopole_pressure(0.1, 0.0, 1000.0, C, RHO_AIR, 1.0).is_err()); // U=0
        assert!(monopole_pressure(0.1, 0.01, 1000.0, 0.0, RHO_AIR, 1.0).is_err()); // c=0
        assert!(monopole_pressure(0.1, 0.01, 1000.0, C, 0.0, 1.0).is_err()); // rho=0
        assert!(monopole_pressure(0.1, 0.01, 1000.0, C, RHO_AIR, f64::NAN).is_err());
        assert!(dipole_directivity(f64::NAN).is_err());
        assert!(dipole_directivity(f64::INFINITY).is_err());
    }

    #[test]
    fn directivity_db_rejects_nonpositive() {
        assert!(directivity_db(0.0).is_err()); // exact null -> error, not -inf
        assert!(directivity_db(-0.5).is_err());
        assert!(directivity_db(f64::NAN).is_err());
    }
}
