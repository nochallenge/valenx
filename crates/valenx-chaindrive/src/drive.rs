//! Single-stage roller-chain drive kinematics and torque.
//!
//! Given a validated [`SprocketPair`] this module computes the four
//! quantities a chain-drive designer reaches for first:
//!
//! - **Chain (link) velocity** from a sprocket's rotational speed —
//!   [`chain_velocity_m_per_s`].
//! - **Driven-shaft speed** from the input speed and the ratio —
//!   [`driven_speed_rpm`].
//! - **Required chain length** in whole pitches (links) for a given
//!   centre distance — [`chain_length_pitches`].
//! - **Output torque** delivered to the driven shaft —
//!   [`output_torque_n_m`].
//!
//! [`analyze`] bundles all of them into a single [`DriveResult`].

use crate::error::ChainDriveError;
use crate::spec::SprocketPair;
use serde::{Deserialize, Serialize};

/// Linear speed of the chain (the speed at which the links travel along
/// the chain line), in metres per second.
///
/// The chain advances one pitch for every tooth that passes the
/// engagement point, so in one shaft revolution the chain travels
/// `pitch × teeth`. At `rpm` revolutions per minute that is
///
/// ```text
/// v = p · z · n / 60_000   [m/s]
/// ```
///
/// with `p` the pitch in **millimetres**, `z` the tooth count and `n`
/// the speed in rev/min (the `/ 60_000` converts mm/min to m/s).
///
/// Because the same chain runs over both sprockets, evaluating this for
/// the driver `(z1, n1)` and the driven `(z2, n2)` gives the **same**
/// chain velocity — that identity is exactly the gear-ratio relation
/// `n2 = n1 · z1 / z2`, and the crate tests check it directly.
///
/// # Errors
///
/// Returns [`ChainDriveError::BadParameter`] if `teeth` is zero or `rpm`
/// is non-finite or negative.
pub fn chain_velocity_m_per_s(pitch_mm: f64, teeth: u32, rpm: f64) -> Result<f64, ChainDriveError> {
    if !pitch_mm.is_finite() || pitch_mm <= 0.0 {
        return Err(ChainDriveError::bad_parameter(
            "pitch_mm",
            format!("must be a finite value > 0, got {pitch_mm}"),
        ));
    }
    if teeth == 0 {
        return Err(ChainDriveError::bad_parameter("teeth", "must be >= 1"));
    }
    if !rpm.is_finite() || rpm < 0.0 {
        return Err(ChainDriveError::bad_parameter(
            "rpm",
            format!("must be a finite value >= 0, got {rpm}"),
        ));
    }
    Ok(pitch_mm * teeth as f64 * rpm / 60_000.0)
}

/// Rotational speed of the **driven** shaft, in rev/min, for a given
/// input (driver) speed.
///
/// From the chain-velocity identity above, `n2 = n1 / ratio = n1 ·
/// N1 / N2`. A reduction (`ratio > 1`) therefore makes the output turn
/// **slower** than the input; a step-up (`ratio < 1`, the bigger
/// sprocket driving) makes it turn **faster**.
///
/// # Errors
///
/// Returns [`ChainDriveError::BadParameter`] if `input_rpm` is
/// non-finite or negative.
pub fn driven_speed_rpm(pair: &SprocketPair, input_rpm: f64) -> Result<f64, ChainDriveError> {
    if !input_rpm.is_finite() || input_rpm < 0.0 {
        return Err(ChainDriveError::bad_parameter(
            "input_rpm",
            format!("must be a finite value >= 0, got {input_rpm}"),
        ));
    }
    Ok(input_rpm / pair.ratio())
}

/// Number of chain pitches (links) required to wrap both sprockets at a
/// given shaft **centre distance**, rounded up to a whole even link
/// count.
///
/// The standard textbook chain-length approximation, with the centre
/// distance `C` expressed in pitches, is
///
/// ```text
///                z1 + z2     (z2 − z1)²
/// L_pitches ≈ 2C + ─────── + ──────────
///                     2        4 π² C
/// ```
///
/// where `z1`, `z2` are the tooth counts. The first term is the two
/// straight runs, the second the average sprocket wrap, and the third a
/// small correction for the wrap difference between unequal sprockets.
///
/// Roller chain is supplied in whole links and an **even** link count
/// avoids an offset (cranked) link, so the real-valued result is rounded
/// **up** to the next even integer; the raw real value is available via
/// [`chain_length_pitches_exact`].
///
/// # Errors
///
/// Returns [`ChainDriveError::BadParameter`] if `center_distance_mm` is
/// non-finite or not strictly positive, and
/// [`ChainDriveError::Degenerate`] if the centre distance is too small
/// for the two sprockets' pitch circles to clear one another (their
/// pitch radii would overlap).
pub fn chain_length_pitches(
    pair: &SprocketPair,
    center_distance_mm: f64,
) -> Result<u32, ChainDriveError> {
    let exact = chain_length_pitches_exact(pair, center_distance_mm)?;
    // Round up to the next whole link, then to the next even link.
    let whole = exact.ceil() as u64;
    let even = whole + (whole % 2);
    Ok(even as u32)
}

/// The raw, real-valued chain length in pitches — the
/// [`chain_length_pitches`] approximation **before** rounding up to a
/// whole even link count.
///
/// Useful when you want the continuous value (e.g. to compare two centre
/// distances) rather than a buildable link count.
///
/// # Errors
///
/// Same validation as [`chain_length_pitches`]: rejects a non-finite or
/// non-positive centre distance, and a centre distance so small the
/// pitch circles overlap.
pub fn chain_length_pitches_exact(
    pair: &SprocketPair,
    center_distance_mm: f64,
) -> Result<f64, ChainDriveError> {
    if !center_distance_mm.is_finite() || center_distance_mm <= 0.0 {
        return Err(ChainDriveError::bad_parameter(
            "center_distance_mm",
            format!("must be a finite value > 0, got {center_distance_mm}"),
        ));
    }

    // Pitch circles must clear: centre distance > sum of pitch radii.
    let min_center = 0.5 * (pair.driver_pitch_diameter_mm() + pair.driven_pitch_diameter_mm());
    if center_distance_mm <= min_center {
        return Err(ChainDriveError::degenerate(format!(
            "centre distance {center_distance_mm} mm <= sum of pitch radii \
             {min_center} mm; sprockets overlap"
        )));
    }

    let c = center_distance_mm / pair.pitch_mm; // centre distance in pitches
    let z1 = pair.driver_teeth as f64;
    let z2 = pair.driven_teeth as f64;
    let dz = z2 - z1;

    let length = 2.0 * c
        + (z1 + z2) / 2.0
        + (dz * dz) / (4.0 * std::f64::consts::PI * std::f64::consts::PI * c);
    Ok(length)
}

/// Output torque delivered to the **driven** shaft, in newton-metres.
///
/// An ideal (loss-free) chain drive conserves power, so torque scales by
/// the same ratio the speed is divided by:
///
/// ```text
/// T_out = T_in · ratio = T_in · N2 / N1
/// ```
///
/// A reduction (`ratio > 1`) therefore **multiplies** torque while
/// reducing speed; a step-up reduces torque. This is the loss-free
/// upper bound — a real drive's output is lower by the chain's
/// mechanical efficiency, which this textbook model does not apply.
///
/// # Errors
///
/// Returns [`ChainDriveError::BadParameter`] if `input_torque_n_m` is
/// non-finite or negative.
pub fn output_torque_n_m(
    pair: &SprocketPair,
    input_torque_n_m: f64,
) -> Result<f64, ChainDriveError> {
    if !input_torque_n_m.is_finite() || input_torque_n_m < 0.0 {
        return Err(ChainDriveError::bad_parameter(
            "input_torque_n_m",
            format!("must be a finite value >= 0, got {input_torque_n_m}"),
        ));
    }
    Ok(input_torque_n_m * pair.ratio())
}

/// The complete set of operating-point quantities for a single-stage
/// chain drive, produced by [`analyze`].
///
/// All fields are derived analytically from the inputs; there is no
/// hidden state or iteration.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct DriveResult {
    /// Speed-reduction ratio `N2 / N1` (driven ÷ driver teeth).
    pub ratio: f64,
    /// Linear chain velocity, m/s (identical whether evaluated on the
    /// driver or the driven sprocket).
    pub chain_velocity_m_per_s: f64,
    /// Driven-shaft rotational speed, rev/min.
    pub driven_speed_rpm: f64,
    /// Output torque at the driven shaft, N·m (loss-free).
    pub output_torque_n_m: f64,
    /// Buildable chain length, whole even number of pitches (links).
    pub chain_length_pitches: u32,
}

/// Analyse a chain drive at one operating point.
///
/// Combines [`SprocketPair::ratio`], [`chain_velocity_m_per_s`]
/// (evaluated on the driver sprocket), [`driven_speed_rpm`],
/// [`output_torque_n_m`] and [`chain_length_pitches`] into a single
/// [`DriveResult`].
///
/// # Errors
///
/// Propagates any [`ChainDriveError`] from the underlying calculations:
/// invalid input speed / torque / centre distance, or a degenerate
/// geometry.
pub fn analyze(
    pair: &SprocketPair,
    input_rpm: f64,
    input_torque_n_m: f64,
    center_distance_mm: f64,
) -> Result<DriveResult, ChainDriveError> {
    Ok(DriveResult {
        ratio: pair.ratio(),
        chain_velocity_m_per_s: chain_velocity_m_per_s(
            pair.pitch_mm,
            pair.driver_teeth,
            input_rpm,
        )?,
        driven_speed_rpm: driven_speed_rpm(pair, input_rpm)?,
        output_torque_n_m: output_torque_n_m(pair, input_torque_n_m)?,
        chain_length_pitches: chain_length_pitches(pair, center_distance_mm)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::SprocketPair;

    const EPS: f64 = 1e-9;

    fn pair(n1: u32, n2: u32, pitch: f64) -> SprocketPair {
        SprocketPair::new(n1, n2, pitch).expect("valid pair")
    }

    #[test]
    fn chain_velocity_matches_hand_calc() {
        // p = 12.7 mm, z = 17, n = 1000 rpm.
        // v = 12.7 * 17 * 1000 / 60_000 = 215_900 / 60_000
        //   = 3.598333... m/s.
        let v = chain_velocity_m_per_s(12.7, 17, 1000.0).unwrap();
        assert!((v - 3.598_333_333_333).abs() < 1e-9, "v = {v}");
    }

    #[test]
    fn chain_velocity_scales_linearly_with_speed() {
        let v1 = chain_velocity_m_per_s(12.7, 17, 500.0).unwrap();
        let v2 = chain_velocity_m_per_s(12.7, 17, 1000.0).unwrap();
        assert!((v2 - 2.0 * v1).abs() < EPS);
    }

    #[test]
    fn chain_velocity_zero_speed_is_zero() {
        let v = chain_velocity_m_per_s(12.7, 17, 0.0).unwrap();
        assert!(v.abs() < EPS);
    }

    #[test]
    fn chain_velocity_identical_on_both_sprockets() {
        // The defining property of a chain: one chain, one velocity.
        // Driver 17T @ 1000 rpm; the matched driven speed is
        // n2 = n1 * z1 / z2. Evaluating the chain velocity on each
        // sprocket must give the same number.
        let p = pair(17, 34, 12.7);
        let n1 = 1000.0;
        let n2 = n1 * p.driver_teeth as f64 / p.driven_teeth as f64;
        let v1 = chain_velocity_m_per_s(p.pitch_mm, p.driver_teeth, n1).unwrap();
        let v2 = chain_velocity_m_per_s(p.pitch_mm, p.driven_teeth, n2).unwrap();
        assert!((v1 - v2).abs() < 1e-9, "v1 = {v1}, v2 = {v2}");
    }

    #[test]
    fn driven_speed_is_input_over_ratio() {
        // 2:1 reduction halves the speed.
        let p = pair(17, 34, 12.7);
        let n2 = driven_speed_rpm(&p, 1000.0).unwrap();
        assert!((n2 - 500.0).abs() < EPS, "n2 = {n2}");
    }

    #[test]
    fn larger_driven_sprocket_runs_slower() {
        // Same input speed; the pair with the bigger driven sprocket
        // (larger ratio) must spin its output more slowly.
        let small = pair(17, 25, 12.7);
        let large = pair(17, 51, 12.7);
        let n_small = driven_speed_rpm(&small, 1000.0).unwrap();
        let n_large = driven_speed_rpm(&large, 1000.0).unwrap();
        assert!(n_large < n_small);
    }

    #[test]
    fn step_up_drive_runs_output_faster_than_input() {
        // Driver bigger than driven -> overdrive -> output faster.
        let p = pair(40, 20, 12.7);
        let n2 = driven_speed_rpm(&p, 1000.0).unwrap();
        assert!(n2 > 1000.0);
        assert!((n2 - 2000.0).abs() < EPS, "n2 = {n2}");
    }

    #[test]
    fn output_torque_is_input_times_ratio() {
        // 2:1 reduction doubles the torque.
        let p = pair(17, 34, 12.7);
        let t = output_torque_n_m(&p, 50.0).unwrap();
        assert!((t - 100.0).abs() < EPS, "t = {t}");
    }

    #[test]
    fn output_torque_drops_in_a_step_up() {
        let p = pair(40, 20, 12.7); // ratio 0.5
        let t = output_torque_n_m(&p, 100.0).unwrap();
        assert!((t - 50.0).abs() < EPS, "t = {t}");
    }

    #[test]
    fn power_is_conserved_speed_down_torque_up() {
        // Loss-free: P_in = P_out. With P = T * omega and
        // omega = 2*pi*n/60, the product T*n must be invariant.
        let p = pair(17, 51, 12.7); // 3:1
        let t_in = 30.0;
        let n_in = 1500.0;
        let t_out = output_torque_n_m(&p, t_in).unwrap();
        let n_out = driven_speed_rpm(&p, n_in).unwrap();
        assert!((t_in * n_in - t_out * n_out).abs() < 1e-6);
    }

    #[test]
    fn chain_length_equal_sprockets_matches_simple_loop() {
        // For equal sprockets the correction term vanishes
        // (dz = 0) and the exact length is 2C + z.
        // C = 30 pitches, z = 17 each -> 2*30 + 17 = 77 pitches.
        let p = pair(17, 17, 10.0);
        let c_mm = 30.0 * 10.0; // 30 pitches
        let exact = chain_length_pitches_exact(&p, c_mm).unwrap();
        assert!((exact - 77.0).abs() < 1e-9, "exact = {exact}");
    }

    #[test]
    fn chain_length_unequal_matches_full_formula() {
        // Independent reference evaluation of the three-term formula.
        // z1 = 15, z2 = 45, pitch = 12.7, C = 500 mm.
        // C_pitches = 500 / 12.7 = 39.37008 pitches.
        // L = 2C + (z1+z2)/2 + (z2-z1)^2 / (4*pi^2*C)
        //   = 78.74016 + 30 + 900 / (39.4784 * 39.37008)
        //   = 78.74016 + 30 + 900 / 1554.366
        //   = 78.74016 + 30 + 0.579020
        //   = 109.319185 pitches.
        let p = pair(15, 45, 12.7);
        let exact = chain_length_pitches_exact(&p, 500.0).unwrap();
        let c = 500.0 / 12.7;
        let expected =
            2.0 * c + 30.0 + 900.0 / (4.0 * std::f64::consts::PI * std::f64::consts::PI * c);
        assert!((exact - expected).abs() < 1e-9, "exact = {exact}");
        assert!((exact - 109.319_185).abs() < 1e-3, "exact = {exact}");
    }

    #[test]
    fn chain_length_rounds_up_to_even_links() {
        // exact = 77.0 (already whole & odd) -> rounds to 78.
        let p = pair(17, 17, 10.0);
        let n = chain_length_pitches(&p, 300.0).unwrap();
        assert_eq!(n, 78);
        assert_eq!(n % 2, 0);
    }

    #[test]
    fn chain_length_grows_with_center_distance() {
        let p = pair(15, 45, 12.7);
        let short = chain_length_pitches_exact(&p, 400.0).unwrap();
        let long = chain_length_pitches_exact(&p, 800.0).unwrap();
        assert!(long > short);
    }

    #[test]
    fn chain_length_rejects_overlapping_sprockets() {
        // Centre distance below the sum of pitch radii is degenerate.
        let p = pair(40, 60, 12.7);
        let too_close = 1.0; // mm — far smaller than the pitch radii
        let err = chain_length_pitches(&p, too_close).unwrap_err();
        assert_eq!(err.code(), "chaindrive.degenerate");
    }

    #[test]
    fn chain_length_rejects_bad_center_distance() {
        let p = pair(17, 34, 12.7);
        let err = chain_length_pitches(&p, 0.0).unwrap_err();
        assert_eq!(err.code(), "chaindrive.bad_parameter");
        let err = chain_length_pitches(&p, f64::INFINITY).unwrap_err();
        assert_eq!(err.code(), "chaindrive.bad_parameter");
    }

    #[test]
    fn analyze_bundles_consistent_values() {
        let p = pair(17, 34, 12.7);
        let r = analyze(&p, 1000.0, 50.0, 500.0).unwrap();
        assert!((r.ratio - 2.0).abs() < EPS);
        assert!((r.driven_speed_rpm - 500.0).abs() < EPS);
        assert!((r.output_torque_n_m - 100.0).abs() < EPS);
        // chain velocity matches the standalone helper.
        let v = chain_velocity_m_per_s(12.7, 17, 1000.0).unwrap();
        assert!((r.chain_velocity_m_per_s - v).abs() < EPS);
        assert!(r.chain_length_pitches > 0);
        assert_eq!(r.chain_length_pitches % 2, 0);
    }

    #[test]
    fn drive_result_serde_round_trip() {
        let p = pair(17, 34, 12.7);
        let r = analyze(&p, 1000.0, 50.0, 500.0).unwrap();
        let json = serde_json::to_string(&r).expect("serialize");
        let back: DriveResult = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(r, back);
    }

    #[test]
    fn rejects_negative_input_speed_and_torque() {
        let p = pair(17, 34, 12.7);
        assert_eq!(
            driven_speed_rpm(&p, -1.0).unwrap_err().code(),
            "chaindrive.bad_parameter"
        );
        assert_eq!(
            output_torque_n_m(&p, -1.0).unwrap_err().code(),
            "chaindrive.bad_parameter"
        );
    }
}
