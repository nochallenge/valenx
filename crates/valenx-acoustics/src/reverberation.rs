//! Statistical-acoustics **reverberation time** (Sabine & Eyring).
//!
//! Reverberation time `RT60` is the time for a room's sound-energy density
//! to decay by 60 dB (to one-millionth of its steady value) after the
//! source stops. Unlike the [modal](crate::room) view, which resolves
//! individual standing waves, this is the *statistical* (diffuse-field)
//! description that governs how "live" or "dead" a space sounds.
//!
//! ## Model
//!
//! Sabine's law assumes a perfectly diffuse field and writes the decay in
//! terms of the room volume `V`, the speed of sound `c`, and the **total
//! absorption** `A = Σ Sᵢ αᵢ` (in sabins, m²) — surface area times
//! absorption coefficient, summed over every surface:
//!
//! ```text
//! RT60(Sabine) = 24 ln(10) · V / (c · A)
//! ```
//!
//! The constant `24 ln(10) ≈ 55.26`, divided by `c = 343 m/s`, recovers
//! the familiar textbook coefficient `RT60 ≈ 0.161 · V / A`.
//!
//! Eyring–Norris corrects Sabine for highly absorptive rooms by replacing
//! `A` with `−S·ln(1 − ᾱ)`, where `S` is the total surface area and `ᾱ`
//! the mean absorption coefficient:
//!
//! ```text
//! RT60(Eyring) = 24 ln(10) · V / ( −c · S · ln(1 − ᾱ) )
//! ```
//!
//! For small `ᾱ`, `−ln(1 − ᾱ) ≈ ᾱ` so Eyring collapses onto Sabine (with
//! `A = S ᾱ`); as `ᾱ → 1` (a fully absorptive room) Eyring correctly tends
//! to zero while Sabine over-predicts.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are the classical diffuse-field
//! statistical formulas: they assume a well-mixed sound field, frequency-
//! independent absorption (no per-band machinery here — pass the
//! band-specific `A`/`ᾱ` yourself), and ignore air absorption, source
//! directivity, room shape beyond `V`/`S`, and the Schroeder transition
//! below which the modal picture dominates. Use them to learn and
//! estimate, not to certify an acoustic design.

use crate::error::{AcousticsError, Result};

/// `24 ln(10) ≈ 55.262`, the numerator constant of the reverberation-time
/// law (so `24 ln10 · V / (c·A)` reduces to the textbook `0.161 · V/A` at
/// `c = 343 m/s`).
const REVERB_CONST: f64 = 24.0 * std::f64::consts::LN_10;

/// Validate a strictly-positive, finite dimensional quantity, surfacing it
/// through [`AcousticsError::InvalidDimension`] under `name` on failure.
fn check_positive(name: &'static str, value: f64) -> Result<f64> {
    if !value.is_finite() || value <= 0.0 {
        return Err(AcousticsError::InvalidDimension { name, value });
    }
    Ok(value)
}

/// Total absorption `A = Σ Sᵢ αᵢ` (sabins, m²) of a set of `(area,
/// coefficient)` surfaces — the quantity the [`sabine_reverberation_time`]
/// needs.
///
/// # Errors
///
/// [`AcousticsError::InvalidDimension`] if any area is negative / non-finite
/// or any absorption coefficient lies outside `[0, 1]`.
pub fn total_absorption(surfaces: &[(f64, f64)]) -> Result<f64> {
    let mut total = 0.0;
    for &(area, alpha) in surfaces {
        if !area.is_finite() || area < 0.0 {
            return Err(AcousticsError::InvalidDimension {
                name: "surface_area",
                value: area,
            });
        }
        if !(0.0..=1.0).contains(&alpha) {
            return Err(AcousticsError::InvalidDimension {
                name: "absorption_coefficient",
                value: alpha,
            });
        }
        total += area * alpha;
    }
    Ok(total)
}

/// **Sabine** reverberation time `RT60 = 24 ln(10) · V / (c · A)` (seconds)
/// for a room of volume `volume` (m³) with total absorption
/// `total_absorption` (sabins, m²) at speed of sound `speed_of_sound`
/// (m/s).
///
/// # Errors
///
/// [`AcousticsError::InvalidDimension`] if `volume` or `total_absorption`
/// is not finite and strictly positive; [`AcousticsError::InvalidSpeedOfSound`]
/// if `speed_of_sound` is not finite and strictly positive.
pub fn sabine_reverberation_time(
    volume: f64,
    total_absorption: f64,
    speed_of_sound: f64,
) -> Result<f64> {
    let v = check_positive("volume", volume)?;
    let a = check_positive("total_absorption", total_absorption)?;
    if !speed_of_sound.is_finite() || speed_of_sound <= 0.0 {
        return Err(AcousticsError::InvalidSpeedOfSound {
            name: "speed_of_sound",
            value: speed_of_sound,
        });
    }
    Ok(REVERB_CONST * v / (speed_of_sound * a))
}

/// **Required total absorption** to hit a target Sabine reverberation time —
/// the exact inverse of [`sabine_reverberation_time`]:
///
/// ```text
/// A = 24 ln(10) · V / (c · RT60)
/// ```
///
/// in sabins (m²). This is the room-treatment sizing question: choose the
/// reverberation time you want (`target_rt60`, e.g. ~0.8 s for a classroom,
/// ~2 s for a concert hall) and read off how many sabins of absorption to
/// install in a room of volume `volume` (m³) at speed of sound
/// `speed_of_sound` (m/s). Feeding the result back into
/// [`sabine_reverberation_time`] reproduces `target_rt60`.
///
/// # Errors
///
/// [`AcousticsError::InvalidDimension`] if `volume` or `target_rt60` is not
/// finite and strictly positive; [`AcousticsError::InvalidSpeedOfSound`] if
/// `speed_of_sound` is not finite and strictly positive.
pub fn absorption_for_reverberation_time(
    volume: f64,
    target_rt60: f64,
    speed_of_sound: f64,
) -> Result<f64> {
    let v = check_positive("volume", volume)?;
    let rt = check_positive("target_rt60", target_rt60)?;
    if !speed_of_sound.is_finite() || speed_of_sound <= 0.0 {
        return Err(AcousticsError::InvalidSpeedOfSound {
            name: "speed_of_sound",
            value: speed_of_sound,
        });
    }
    Ok(REVERB_CONST * v / (speed_of_sound * rt))
}

/// **Eyring–Norris** reverberation time
/// `RT60 = 24 ln(10) · V / (−c · S · ln(1 − ᾱ))` (seconds) for a room of
/// volume `volume` (m³), total surface area `surface_area` (m²), mean
/// absorption coefficient `mean_absorption` (`0 < ᾱ < 1`) at speed of sound
/// `speed_of_sound` (m/s).
///
/// More accurate than Sabine for absorptive rooms; it collapses onto Sabine
/// for small `ᾱ`.
///
/// # Errors
///
/// [`AcousticsError::InvalidDimension`] if `volume` / `surface_area` is not
/// finite and strictly positive, or `mean_absorption` is not strictly inside
/// `(0, 1)`; [`AcousticsError::InvalidSpeedOfSound`] for a bad
/// `speed_of_sound`.
pub fn eyring_reverberation_time(
    volume: f64,
    surface_area: f64,
    mean_absorption: f64,
    speed_of_sound: f64,
) -> Result<f64> {
    let v = check_positive("volume", volume)?;
    let s = check_positive("surface_area", surface_area)?;
    if !mean_absorption.is_finite() || mean_absorption <= 0.0 || mean_absorption >= 1.0 {
        return Err(AcousticsError::InvalidDimension {
            name: "mean_absorption",
            value: mean_absorption,
        });
    }
    if !speed_of_sound.is_finite() || speed_of_sound <= 0.0 {
        return Err(AcousticsError::InvalidSpeedOfSound {
            name: "speed_of_sound",
            value: speed_of_sound,
        });
    }
    Ok(REVERB_CONST * v / (-speed_of_sound * s * (1.0 - mean_absorption).ln()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::doppler::speed_of_sound;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn sabine_matches_textbook_value_and_constant() {
        let c = speed_of_sound(20.0).unwrap(); // ~343 m/s
        let rt = sabine_reverberation_time(200.0, 40.0, c).unwrap();
        // Exact closed form.
        assert!(close(rt, REVERB_CONST * 200.0 / (c * 40.0), 1e-12));
        // Textbook coefficient RT60 ~ 0.161 V/A: 0.161*200/40 = 0.805 s.
        assert!(close(rt, 0.805, 5e-3), "RT60 = {rt}");
        // The lumped constant 24 ln10 / c is ~0.161 at room temperature.
        assert!(close(REVERB_CONST / c, 0.161, 1e-3));
    }

    #[test]
    fn sabine_scales_inversely_with_absorption_and_with_volume() {
        let c = 343.0;
        let base = sabine_reverberation_time(150.0, 30.0, c).unwrap();
        // Double the absorption -> half the time.
        let more_abs = sabine_reverberation_time(150.0, 60.0, c).unwrap();
        assert!(close(more_abs, base / 2.0, 1e-9));
        // Double the volume -> double the time.
        let bigger = sabine_reverberation_time(300.0, 30.0, c).unwrap();
        assert!(close(bigger, base * 2.0, 1e-9));
    }

    #[test]
    fn absorption_inverse_round_trips_with_sabine() {
        let c = 343.0;
        // inverse -> forward: size A for a target RT60, Sabine reproduces it.
        let target = 0.8;
        let a = absorption_for_reverberation_time(200.0, target, c).unwrap();
        let rt = sabine_reverberation_time(200.0, a, c).unwrap();
        assert!(close(rt, target, 1e-9), "rt {rt} vs target {target}");
        // forward -> inverse: A recovered from the RT60 it produces.
        let a0 = 35.0;
        let rt0 = sabine_reverberation_time(180.0, a0, c).unwrap();
        let a_back = absorption_for_reverberation_time(180.0, rt0, c).unwrap();
        assert!(close(a_back, a0, 1e-9), "A {a_back} vs {a0}");
    }

    #[test]
    fn absorption_for_target_matches_textbook_coefficient() {
        // Inverse of the textbook forward case: V=200, RT60=0.805 s -> A ~ 40
        // sabins (0.161 V / RT60 = 0.161*200/0.805 = 40.0).
        let c = speed_of_sound(20.0).unwrap();
        let a = absorption_for_reverberation_time(200.0, 0.805, c).unwrap();
        assert!(close(a, REVERB_CONST * 200.0 / (c * 0.805), 1e-12));
        assert!(close(a, 40.0, 0.2), "A = {a}");
    }

    #[test]
    fn absorption_inverse_rejects_bad_inputs() {
        let c = 343.0;
        assert!(absorption_for_reverberation_time(0.0, 0.8, c).is_err()); // volume
        assert!(absorption_for_reverberation_time(200.0, 0.0, c).is_err()); // rt60
        assert!(absorption_for_reverberation_time(200.0, 0.8, 0.0).is_err()); // c
    }

    #[test]
    fn total_absorption_sums_area_times_coefficient() {
        // 10 m² at 0.1 + 20 m² at 0.3 = 1 + 6 = 7 sabins.
        let a = total_absorption(&[(10.0, 0.1), (20.0, 0.3)]).unwrap();
        assert!(close(a, 7.0, 1e-12));
        // Feeding it into Sabine is consistent with the bare-A call.
        let c = 343.0;
        let rt = sabine_reverberation_time(120.0, a, c).unwrap();
        assert!(close(rt, REVERB_CONST * 120.0 / (c * 7.0), 1e-12));
    }

    #[test]
    fn eyring_collapses_to_sabine_at_low_absorption_and_is_lower_when_high() {
        let c = 343.0;
        let s = 240.0;
        // Low mean absorption: Eyring ~ Sabine with A = S*alpha.
        let alpha_lo = 0.02;
        let eyr_lo = eyring_reverberation_time(200.0, s, alpha_lo, c).unwrap();
        let sab_lo = sabine_reverberation_time(200.0, s * alpha_lo, c).unwrap();
        assert!(
            (eyr_lo - sab_lo).abs() / sab_lo < 0.02,
            "eyr {eyr_lo} vs sab {sab_lo}"
        );
        // High mean absorption: Eyring predicts a shorter (more correct) tail.
        let alpha_hi = 0.5;
        let eyr_hi = eyring_reverberation_time(200.0, s, alpha_hi, c).unwrap();
        let sab_hi = sabine_reverberation_time(200.0, s * alpha_hi, c).unwrap();
        assert!(
            eyr_hi < sab_hi,
            "eyring {eyr_hi} should be below sabine {sab_hi}"
        );
    }

    #[test]
    fn rejects_bad_inputs() {
        let c = 343.0;
        assert!(sabine_reverberation_time(0.0, 40.0, c).is_err());
        assert!(sabine_reverberation_time(200.0, 0.0, c).is_err());
        assert!(sabine_reverberation_time(200.0, 40.0, 0.0).is_err());
        assert!(eyring_reverberation_time(200.0, 240.0, 0.0, c).is_err()); // alpha <= 0
        assert!(eyring_reverberation_time(200.0, 240.0, 1.0, c).is_err()); // alpha >= 1
        assert!(total_absorption(&[(10.0, 1.5)]).is_err()); // coeff > 1
        assert!(total_absorption(&[(-1.0, 0.5)]).is_err()); // negative area
    }
}
