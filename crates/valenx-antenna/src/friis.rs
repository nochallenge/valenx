//! The Friis transmission equation and free-space path loss.
//!
//! ## Model
//!
//! For two antennas separated by a line-of-sight distance `d` in free
//! space, with transmit gain `Gt` and receive gain `Gr` (both linear
//! ratios) and wavelength `lambda`, the received-to-transmitted power
//! ratio is the **Friis transmission equation**
//!
//! ```text
//! Pr / Pt = Gt * Gr * (lambda / (4 * pi * d))^2
//! ```
//!
//! The factor `(lambda / (4*pi*d))^2` is the free-space spreading term;
//! its reciprocal (for isotropic antennas) is the **free-space path
//! loss** `FSPL = (4*pi*d / lambda)^2`. In decibels,
//!
//! ```text
//! FSPL_dB = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
//! ```
//!
//! The equation is valid in the **far field** of both antennas; this
//! crate does not enforce a far-field distance — that is the caller's
//! responsibility — but the relations themselves are the canonical
//! closed forms.

use crate::error::{require_non_negative_gain, require_positive, AntennaError};

/// The free-space spreading factor `(lambda / (4*pi*d))^2`
/// (dimensionless), i.e. the Friis ratio for isotropic
/// (`Gt = Gr = 1`) antennas.
///
/// # Errors
///
/// Returns an error if `wavelength_m` or `distance_m` is not finite and
/// strictly positive.
pub fn free_space_factor(wavelength_m: f64, distance_m: f64) -> Result<f64, AntennaError> {
    let lambda = require_positive("wavelength_m", wavelength_m)?;
    let d = require_positive("distance_m", distance_m)?;
    let s = lambda / (4.0 * core::f64::consts::PI * d);
    Ok(s * s)
}

/// The Friis power ratio `Pr/Pt = Gt * Gr * (lambda/(4*pi*d))^2`
/// (dimensionless).
///
/// `gt` and `gr` are **linear** gain ratios (use
/// [`crate::gain::from_dbi`] to convert from dBi first).
///
/// # Errors
///
/// Returns an error if `gt` or `gr` is not finite and non-negative, or
/// if `wavelength_m` / `distance_m` is not finite and strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_antenna::friis::power_ratio;
/// // Isotropic antennas, lambda = 1 m, d = 1 m.
/// let r = power_ratio(1.0, 1.0, 1.0, 1.0).unwrap();
/// let expected = (1.0_f64 / (4.0 * std::f64::consts::PI)).powi(2);
/// assert!((r - expected).abs() < 1e-12);
/// ```
pub fn power_ratio(
    gt: f64,
    gr: f64,
    wavelength_m: f64,
    distance_m: f64,
) -> Result<f64, AntennaError> {
    let gt = require_non_negative_gain("gt", gt)?;
    let gr = require_non_negative_gain("gr", gr)?;
    let factor = free_space_factor(wavelength_m, distance_m)?;
    Ok(gt * gr * factor)
}

/// Received power (watts) for a transmitted power `pt_w` over a Friis
/// link.
///
/// `Pr = Pt * Gt * Gr * (lambda/(4*pi*d))^2`.
///
/// # Errors
///
/// Returns an error if `pt_w` is not finite and strictly positive, the
/// gains are not finite and non-negative, or `wavelength_m` /
/// `distance_m` is not finite and strictly positive.
pub fn received_power(
    pt_w: f64,
    gt: f64,
    gr: f64,
    wavelength_m: f64,
    distance_m: f64,
) -> Result<f64, AntennaError> {
    let pt = require_positive("pt_w", pt_w)?;
    let ratio = power_ratio(gt, gr, wavelength_m, distance_m)?;
    Ok(pt * ratio)
}

/// Free-space path loss as a linear ratio `(4*pi*d/lambda)^2`
/// (dimensionless, `>= 1` for `d >> lambda`). This is the reciprocal of
/// [`free_space_factor`].
///
/// # Errors
///
/// Returns an error if `wavelength_m` or `distance_m` is not finite and
/// strictly positive.
pub fn free_space_path_loss(wavelength_m: f64, distance_m: f64) -> Result<f64, AntennaError> {
    let lambda = require_positive("wavelength_m", wavelength_m)?;
    let d = require_positive("distance_m", distance_m)?;
    let l = 4.0 * core::f64::consts::PI * d / lambda;
    Ok(l * l)
}

/// Free-space path loss in decibels:
/// `FSPL_dB = 10*log10((4*pi*d/lambda)^2) = 20*log10(4*pi*d/lambda)`.
///
/// # Errors
///
/// Returns an error if `wavelength_m` or `distance_m` is not finite and
/// strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_antenna::friis::free_space_path_loss_db;
/// use valenx_antenna::wave::wavelength_from_frequency;
/// // ~92.45 dB at 1 km and 1 GHz (the textbook reference value).
/// let lambda = wavelength_from_frequency(1.0e9).unwrap();
/// let fspl = free_space_path_loss_db(lambda, 1_000.0).unwrap();
/// assert!((fspl - 92.448).abs() < 0.01);
/// ```
pub fn free_space_path_loss_db(wavelength_m: f64, distance_m: f64) -> Result<f64, AntennaError> {
    let l = free_space_path_loss(wavelength_m, distance_m)?;
    // l > 0 by construction, so the log is well-defined.
    Ok(10.0 * l.log10())
}

/// The maximum line-of-sight range (metres) at which a Friis link still
/// delivers at least a target received power, inverting
/// [`received_power`] for the distance `d`:
///
/// ```text
/// d_max = (lambda / (4*pi)) * sqrt(Pt * Gt * Gr / Pr_min)
/// ```
///
/// Given a transmit power `pt_w`, the two **linear** gains `gt` / `gr`,
/// the wavelength, and a receiver sensitivity `pr_min_w` (the minimum
/// usable received power), this is the distance at which the received
/// power has fallen to exactly `pr_min_w`; beyond it the link drops below
/// sensitivity. Range grows as the square root of transmit power and
/// gain, linearly with wavelength, and as the inverse square root of the
/// required sensitivity.
///
/// # Errors
///
/// Returns an error if `pt_w` or `pr_min_w` is not finite and strictly
/// positive, the gains are not finite and non-negative, or
/// `wavelength_m` is not finite and strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_antenna::friis::{max_range_m, received_power};
/// // The received power at d_max equals the sensitivity exactly.
/// let d = max_range_m(5.0, 10.0, 4.0, 0.125, 1.0e-9).unwrap();
/// let pr = received_power(5.0, 10.0, 4.0, 0.125, d).unwrap();
/// assert!((pr - 1.0e-9).abs() / 1.0e-9 < 1e-9);
/// ```
pub fn max_range_m(
    pt_w: f64,
    gt: f64,
    gr: f64,
    wavelength_m: f64,
    pr_min_w: f64,
) -> Result<f64, AntennaError> {
    let pt = require_positive("pt_w", pt_w)?;
    let gt = require_non_negative_gain("gt", gt)?;
    let gr = require_non_negative_gain("gr", gr)?;
    let lambda = require_positive("wavelength_m", wavelength_m)?;
    let pr_min = require_positive("pr_min_w", pr_min_w)?;
    let d = (lambda / (4.0 * core::f64::consts::PI)) * (pt * gt * gr / pr_min).sqrt();
    Ok(d)
}

/// The transmit power (watts) needed to deliver a target received power
/// across a Friis link of known length — the link-budget sizing inverse
/// of [`received_power`] solved for `Pt`:
///
/// ```text
/// Pt = Pr_min * (4*pi*d / lambda)^2 / (Gt * Gr) = Pr_min / power_ratio
/// ```
///
/// Given the receiver sensitivity `pr_min_w`, the two **linear** gains,
/// the wavelength and the link distance `distance_m`, this is the transmit
/// power at which the received power equals exactly `pr_min_w`. It scales
/// with the free-space path loss (as `d^2` and `1/lambda^2`) and inversely
/// with the gain product, mirroring the way [`max_range_m`] inverts the
/// same link for distance instead of power.
///
/// # Errors
///
/// Returns an error if `pr_min_w`, `wavelength_m` or `distance_m` is not
/// finite and strictly positive, if `gt` or `gr` is not finite and
/// non-negative, or if the gain product is zero (a link with `Gt*Gr = 0`
/// radiates nothing toward the receiver, so no finite transmit power can
/// close it).
///
/// # Examples
///
/// ```
/// use valenx_antenna::friis::{received_power, required_transmit_power};
/// // Size the transmit power, then confirm the link delivers exactly the
/// // target sensitivity at that power.
/// let pt = required_transmit_power(1.0e-9, 10.0, 4.0, 0.125, 1_500.0).unwrap();
/// let pr = received_power(pt, 10.0, 4.0, 0.125, 1_500.0).unwrap();
/// assert!((pr - 1.0e-9).abs() / 1.0e-9 < 1e-9);
/// ```
pub fn required_transmit_power(
    pr_min_w: f64,
    gt: f64,
    gr: f64,
    wavelength_m: f64,
    distance_m: f64,
) -> Result<f64, AntennaError> {
    let pr_min = require_positive("pr_min_w", pr_min_w)?;
    let ratio = power_ratio(gt, gr, wavelength_m, distance_m)?;
    // The ratio is zero only when a gain is zero; such a link cannot close.
    let ratio = require_positive("power_ratio", ratio)?;
    Ok(pr_min / ratio)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::wave::wavelength_from_frequency;
    use core::f64::consts::PI;

    const EPS: f64 = 1e-12;

    #[test]
    fn isotropic_unit_case() {
        let r = power_ratio(1.0, 1.0, 1.0, 1.0).unwrap();
        let expected = (1.0 / (4.0 * PI)).powi(2);
        assert!((r - expected).abs() < EPS);
    }

    #[test]
    fn ratio_scales_as_inverse_distance_squared() {
        // Doubling distance quarters the received-power ratio.
        let near = power_ratio(2.0, 3.0, 0.3, 10.0).unwrap();
        let far = power_ratio(2.0, 3.0, 0.3, 20.0).unwrap();
        assert!((far - near / 4.0).abs() < 1e-15);
        // Tripling distance -> 1/9.
        let far3 = power_ratio(2.0, 3.0, 0.3, 30.0).unwrap();
        assert!((far3 - near / 9.0).abs() < 1e-15);
    }

    #[test]
    fn ratio_scales_as_inverse_frequency_squared() {
        // Pr/Pt ~ lambda^2 ~ (1/f)^2 at fixed gain & distance, so
        // doubling the frequency quarters the ratio.
        let d = 100.0;
        let lambda_lo = wavelength_from_frequency(1.0e9).unwrap();
        let lambda_hi = wavelength_from_frequency(2.0e9).unwrap();
        let r_lo = power_ratio(1.0, 1.0, lambda_lo, d).unwrap();
        let r_hi = power_ratio(1.0, 1.0, lambda_hi, d).unwrap();
        assert!(
            (r_hi - r_lo / 4.0).abs() / r_lo < 1e-12,
            "1/f^2 scaling failed: r_lo={r_lo}, r_hi={r_hi}"
        );
    }

    #[test]
    fn ratio_linear_in_each_gain() {
        let base = power_ratio(1.0, 1.0, 0.5, 50.0).unwrap();
        let gt5 = power_ratio(5.0, 1.0, 0.5, 50.0).unwrap();
        let gr7 = power_ratio(1.0, 7.0, 0.5, 50.0).unwrap();
        let both = power_ratio(5.0, 7.0, 0.5, 50.0).unwrap();
        assert!((gt5 - 5.0 * base).abs() < 1e-15);
        assert!((gr7 - 7.0 * base).abs() < 1e-15);
        assert!((both - 35.0 * base).abs() < 1e-15);
    }

    #[test]
    fn received_power_is_pt_times_ratio() {
        let pt = 10.0;
        let ratio = power_ratio(2.0, 4.0, 0.12, 250.0).unwrap();
        let pr = received_power(pt, 2.0, 4.0, 0.12, 250.0).unwrap();
        assert!((pr - pt * ratio).abs() < 1e-15);
    }

    #[test]
    fn path_loss_is_reciprocal_of_factor() {
        let lambda = 0.3;
        let d = 1234.0;
        let factor = free_space_factor(lambda, d).unwrap();
        let loss = free_space_path_loss(lambda, d).unwrap();
        assert!((factor * loss - 1.0).abs() < 1e-9);
    }

    #[test]
    fn fspl_db_reference_1km_1ghz() {
        // Canonical textbook value: FSPL(1 km, 1 GHz) ~ 92.45 dB.
        let lambda = wavelength_from_frequency(1.0e9).unwrap();
        let fspl = free_space_path_loss_db(lambda, 1_000.0).unwrap();
        assert!((fspl - 92.448).abs() < 0.01, "got {fspl} dB");
    }

    #[test]
    fn fspl_db_reference_1km_2400mhz() {
        // FSPL(1 km, 2.4 GHz) ~ 100.05 dB.
        let lambda = wavelength_from_frequency(2.4e9).unwrap();
        let fspl = free_space_path_loss_db(lambda, 1_000.0).unwrap();
        assert!((fspl - 100.05).abs() < 0.02, "got {fspl} dB");
    }

    #[test]
    fn fspl_db_adds_six_db_per_distance_doubling() {
        let lambda = 0.1;
        let l1 = free_space_path_loss_db(lambda, 100.0).unwrap();
        let l2 = free_space_path_loss_db(lambda, 200.0).unwrap();
        // 20*log10(2) ~ 6.0206 dB.
        assert!((l2 - l1 - 6.020_599_913).abs() < 1e-6);
    }

    #[test]
    fn fspl_db_adds_six_db_per_frequency_doubling() {
        let l1 =
            free_space_path_loss_db(wavelength_from_frequency(1.0e9).unwrap(), 1_000.0).unwrap();
        let l2 =
            free_space_path_loss_db(wavelength_from_frequency(2.0e9).unwrap(), 1_000.0).unwrap();
        assert!((l2 - l1 - 6.020_599_913).abs() < 1e-6);
    }

    #[test]
    fn friis_consistent_with_path_loss() {
        // For isotropic antennas, Pr/Pt should equal 1 / FSPL_linear.
        let lambda = 0.21;
        let d = 4321.0;
        let ratio = power_ratio(1.0, 1.0, lambda, d).unwrap();
        let loss = free_space_path_loss(lambda, d).unwrap();
        assert!((ratio - 1.0 / loss).abs() / ratio < 1e-12);
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(power_ratio(-1.0, 1.0, 0.1, 10.0).is_err());
        assert!(power_ratio(1.0, 1.0, 0.0, 10.0).is_err());
        assert!(power_ratio(1.0, 1.0, 0.1, 0.0).is_err());
        assert!(received_power(0.0, 1.0, 1.0, 0.1, 10.0).is_err());
        assert!(free_space_path_loss_db(0.1, -5.0).is_err());
    }

    #[test]
    fn max_range_round_trips_through_received_power() {
        // The defining property: at d_max the received power is exactly
        // the sensitivity, so received_power inverts max_range_m.
        let (pt, gt, gr) = (5.0, 10.0, 4.0);
        let lambda = wavelength_from_frequency(2.4e9).unwrap();
        let pr_min = 1.0e-9; // 1 nW.
        let d = max_range_m(pt, gt, gr, lambda, pr_min).unwrap();
        let pr = received_power(pt, gt, gr, lambda, d).unwrap();
        assert!((pr - pr_min).abs() / pr_min < 1e-9, "got {pr}");
    }

    #[test]
    fn max_range_matches_closed_form() {
        // Isotropic, lambda = 1, Pt = 1, Pr_min = (1/4pi)^2 -> d_max = 1
        // (the received power one metre away for unit isotropic link).
        let pr_min = (1.0 / (4.0 * PI)).powi(2);
        let d = max_range_m(1.0, 1.0, 1.0, 1.0, pr_min).unwrap();
        assert!((d - 1.0).abs() < 1e-9, "got {d}");
    }

    #[test]
    fn max_range_scaling_laws() {
        let (gt, gr, lambda, pr) = (3.0, 6.0, 0.125, 2.0e-10);
        let base = max_range_m(1.0, gt, gr, lambda, pr).unwrap();
        // 4x transmit power -> 2x range (square-root law).
        let quad_pt = max_range_m(4.0, gt, gr, lambda, pr).unwrap();
        assert!((quad_pt / base - 2.0).abs() < 1e-9, "power scaling");
        // 2x wavelength -> 2x range (linear).
        let dbl_lambda = max_range_m(1.0, gt, gr, 2.0 * lambda, pr).unwrap();
        assert!((dbl_lambda / base - 2.0).abs() < 1e-9, "wavelength scaling");
        // 4x stricter sensitivity -> half the range.
        let worse_rx = max_range_m(1.0, gt, gr, lambda, 4.0 * pr).unwrap();
        assert!((worse_rx / base - 0.5).abs() < 1e-9, "sensitivity scaling");
    }

    #[test]
    fn max_range_rejects_bad_inputs() {
        let lambda = 0.1;
        assert!(max_range_m(0.0, 1.0, 1.0, lambda, 1e-9).is_err()); // pt <= 0
        assert!(max_range_m(1.0, -1.0, 1.0, lambda, 1e-9).is_err()); // negative gain
        assert!(max_range_m(1.0, 1.0, 1.0, 0.0, 1e-9).is_err()); // lambda <= 0
        assert!(max_range_m(1.0, 1.0, 1.0, lambda, 0.0).is_err()); // pr_min <= 0
        assert!(max_range_m(f64::NAN, 1.0, 1.0, lambda, 1e-9).is_err()); // non-finite
    }

    #[test]
    fn required_transmit_power_round_trips_received_power() {
        // forward∘inverse: size Pt for a target Pr, confirm received_power
        // reproduces it; then the reverse direction.
        let (gt, gr) = (10.0, 4.0);
        let lambda = wavelength_from_frequency(2.4e9).unwrap();
        let d = 1_500.0;
        let pr_target = 2.0e-9;
        let pt = required_transmit_power(pr_target, gt, gr, lambda, d).unwrap();
        let pr = received_power(pt, gt, gr, lambda, d).unwrap();
        assert!((pr - pr_target).abs() / pr_target < 1e-9, "got {pr}");
        // Pt -> Pr -> Pt.
        let pt0 = 7.5;
        let pr0 = received_power(pt0, gt, gr, lambda, d).unwrap();
        let pt_back = required_transmit_power(pr0, gt, gr, lambda, d).unwrap();
        assert!((pt_back - pt0).abs() / pt0 < 1e-9, "got {pt_back}");
    }

    #[test]
    fn required_transmit_power_matches_closed_form() {
        // Pt = Pr_min * FSPL / (Gt*Gr).
        let (gt, gr) = (3.0, 6.0);
        let lambda = 0.125;
        let d = 800.0;
        let pr_min = 1.0e-10;
        let pt = required_transmit_power(pr_min, gt, gr, lambda, d).unwrap();
        let fspl = free_space_path_loss(lambda, d).unwrap();
        let expected = pr_min * fspl / (gt * gr);
        assert!(
            (pt - expected).abs() / expected < 1e-12,
            "got {pt} vs {expected}"
        );
    }

    #[test]
    fn required_transmit_power_consistent_with_max_range() {
        // The two inverses agree: at d = max_range_m(pt, ...), the transmit
        // power required to reach pr_min is exactly pt again.
        let (gt, gr) = (8.0, 5.0);
        let lambda = wavelength_from_frequency(900.0e6).unwrap();
        let pr_min = 5.0e-10;
        let pt = 12.0;
        let d_max = max_range_m(pt, gt, gr, lambda, pr_min).unwrap();
        let pt_back = required_transmit_power(pr_min, gt, gr, lambda, d_max).unwrap();
        assert!((pt_back - pt).abs() / pt < 1e-9, "got {pt_back}");
    }

    #[test]
    fn required_transmit_power_scaling_laws() {
        let (gt, gr, lambda, pr) = (4.0, 2.0, 0.1, 1.0e-9);
        let base = required_transmit_power(pr, gt, gr, lambda, 100.0).unwrap();
        // 2x distance -> 4x power (d^2).
        let dbl_d = required_transmit_power(pr, gt, gr, lambda, 200.0).unwrap();
        assert!((dbl_d / base - 4.0).abs() < 1e-9, "distance scaling");
        // 2x wavelength -> 1/4 power (1/lambda^2).
        let dbl_l = required_transmit_power(pr, gt, gr, 2.0 * lambda, 100.0).unwrap();
        assert!((dbl_l / base - 0.25).abs() < 1e-9, "wavelength scaling");
        // 2x sensitivity floor -> 2x power (linear in Pr_min).
        let dbl_pr = required_transmit_power(2.0 * pr, gt, gr, lambda, 100.0).unwrap();
        assert!((dbl_pr / base - 2.0).abs() < 1e-9, "sensitivity scaling");
        // 2x gain product -> half power.
        let dbl_g = required_transmit_power(pr, 8.0, gr, lambda, 100.0).unwrap();
        assert!((dbl_g / base - 0.5).abs() < 1e-9, "gain scaling");
    }

    #[test]
    fn required_transmit_power_rejects_bad_inputs() {
        let lambda = 0.1;
        assert!(required_transmit_power(0.0, 1.0, 1.0, lambda, 10.0).is_err()); // pr_min <= 0
        assert!(required_transmit_power(1e-9, -1.0, 1.0, lambda, 10.0).is_err()); // negative gain
        assert!(required_transmit_power(1e-9, 1.0, 1.0, 0.0, 10.0).is_err()); // lambda <= 0
        assert!(required_transmit_power(1e-9, 1.0, 1.0, lambda, 0.0).is_err()); // d <= 0
        assert!(required_transmit_power(1e-9, 0.0, 1.0, lambda, 10.0).is_err()); // zero gain product
        assert!(required_transmit_power(f64::NAN, 1.0, 1.0, lambda, 10.0).is_err());
        // non-finite
    }
}
