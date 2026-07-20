//! Initial (small-angle) transverse stability: the waterplane second
//! moment, the metacentric radius `BM`, and the metacentric height
//! `GM`.
//!
//! ## Governing equations
//!
//! For small heel angles the transverse metacentre `M` sits a height
//! `BM` above the centre of buoyancy `B`, where
//!
//! ```text
//! BM = I_wp / V_disp
//! ```
//!
//! with `I_wp` the second moment of the waterplane area about the
//! longitudinal (fore-aft) centreline and `V_disp` the displaced
//! volume. The metacentric height — the lever between the centre of
//! gravity `G` and the metacentre `M` — is
//!
//! ```text
//! GM = KM - KG = (KB + BM) - KG = BM - BG
//! ```
//!
//! where `BG = KG - KB` is the (signed) height of `G` above `B`. A
//! body is in stable initial equilibrium when `GM > 0`, neutral when
//! `GM == 0`, and unstable when `GM < 0`.
//!
//! ## Rectangular-barge ground truth
//!
//! For a box of beam `B` floating at draft `T`, the waterplane is a
//! rectangle `L * B`, so its second moment about the centreline is
//! `I_wp = L * B^3 / 12` and the displaced volume is `V = L * B * T`.
//! Their ratio collapses to the classic closed form
//!
//! ```text
//! BM = (L * B^3 / 12) / (L * B * T) = B^2 / (12 * T)
//! ```
//!
//! independent of length `L`.

use crate::error::{BuoyancyError, Result};

/// Second moment of a rectangular waterplane about its own longitudinal
/// centreline, `I = L * B^3 / 12` (m^4), for a waterplane of length
/// `length_m` and beam `beam_m`.
///
/// This is the area moment of inertia that governs transverse (roll)
/// stability — note the beam is cubed, so it dominates.
///
/// Both arguments must be finite and strictly positive.
///
/// ```
/// use valenx_buoyancy::stability::rectangular_waterplane_inertia;
/// // L = 10, B = 4 -> 10 * 64 / 12 = 53.333...
/// let i = rectangular_waterplane_inertia(10.0, 4.0).unwrap();
/// assert!((i - 640.0 / 12.0).abs() < 1e-9);
/// ```
pub fn rectangular_waterplane_inertia(length_m: f64, beam_m: f64) -> Result<f64> {
    let l = BuoyancyError::positive("length_m", length_m)?;
    let b = BuoyancyError::positive("beam_m", beam_m)?;
    Ok(l * b.powi(3) / 12.0)
}

/// Metacentric radius `BM = I_wp / V_disp` (m): the height of the
/// transverse metacentre above the centre of buoyancy, from the
/// waterplane second moment `i_waterplane_m4` and displaced volume
/// `displaced_volume_m3`.
///
/// Both arguments must be finite; `i_waterplane_m4` must be
/// non-negative and `displaced_volume_m3` strictly positive.
///
/// ```
/// use valenx_buoyancy::stability::metacentric_radius;
/// // I = 53.333..., V = 40 -> BM = 1.333...
/// let bm = metacentric_radius(640.0 / 12.0, 40.0).unwrap();
/// assert!((bm - (640.0 / 12.0) / 40.0).abs() < 1e-9);
/// ```
pub fn metacentric_radius(i_waterplane_m4: f64, displaced_volume_m3: f64) -> Result<f64> {
    let i = BuoyancyError::non_negative("i_waterplane_m4", i_waterplane_m4)?;
    let v = BuoyancyError::positive("displaced_volume_m3", displaced_volume_m3)?;
    Ok(i / v)
}

/// Metacentric radius of an upright rectangular box (barge) floating at
/// draft `draft_m` with beam `beam_m`, via the closed form
/// `BM = B^2 / (12 * T)` (m).
///
/// This is the analytic ground truth that
/// [`metacentric_radius`]`(rectangular_waterplane_inertia(L, B), L*B*T)`
/// must reproduce for any length `L`.
///
/// Both arguments must be finite and strictly positive.
///
/// ```
/// use valenx_buoyancy::stability::barge_metacentric_radius;
/// // B = 4, T = 1 -> 16 / 12 = 1.333...
/// let bm = barge_metacentric_radius(4.0, 1.0).unwrap();
/// assert!((bm - 16.0 / 12.0).abs() < 1e-9);
/// ```
pub fn barge_metacentric_radius(beam_m: f64, draft_m: f64) -> Result<f64> {
    let b = BuoyancyError::positive("beam_m", beam_m)?;
    let t = BuoyancyError::positive("draft_m", draft_m)?;
    Ok(b.powi(2) / (12.0 * t))
}

/// Metacentric height `GM = BM - BG` (m), where `bm_m` is the
/// metacentric radius and `bg_m` is the height of the centre of gravity
/// `G` above the centre of buoyancy `B` (`BG = KG - KB`).
///
/// A positive `GM` is the small-angle condition for stable floating
/// equilibrium; see [`stability_verdict`]. `bg_m` may be negative when
/// `G` lies below `B` (an inherently stable arrangement), so it is only
/// required to be finite.
///
/// Both arguments must be finite.
///
/// ```
/// use valenx_buoyancy::stability::metacentric_height;
/// // BM = 1.333..., BG = 0.5 -> GM = 0.833...
/// let gm = metacentric_height(16.0 / 12.0, 0.5).unwrap();
/// assert!((gm - (16.0 / 12.0 - 0.5)).abs() < 1e-9);
/// ```
pub fn metacentric_height(bm_m: f64, bg_m: f64) -> Result<f64> {
    let bm = BuoyancyError::finite("bm_m", bm_m)?;
    let bg = BuoyancyError::finite("bg_m", bg_m)?;
    Ok(bm - bg)
}

/// Small-angle stability verdict from the sign of `GM`, returned by
/// [`stability_verdict`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum StabilityVerdict {
    /// `GM > tol`: the body returns toward upright when heeled.
    Stable,
    /// `|GM| <= tol`: indifferent equilibrium.
    Neutral,
    /// `GM < -tol`: the body capsizes away from upright.
    Unstable,
}

/// Classify initial stability from `gm_m` (the metacentric height, m).
/// `GM` greater than `tol_m` is [`StabilityVerdict::Stable`], within
/// `+/- tol_m` is [`StabilityVerdict::Neutral`], and below `-tol_m` is
/// [`StabilityVerdict::Unstable`].
///
/// `gm_m` must be finite; `tol_m` must be finite and non-negative.
///
/// ```
/// use valenx_buoyancy::stability::{stability_verdict, StabilityVerdict};
/// assert_eq!(stability_verdict(0.8, 1e-6).unwrap(), StabilityVerdict::Stable);
/// assert_eq!(stability_verdict(0.0, 1e-6).unwrap(), StabilityVerdict::Neutral);
/// assert_eq!(stability_verdict(-0.3, 1e-6).unwrap(), StabilityVerdict::Unstable);
/// ```
pub fn stability_verdict(gm_m: f64, tol_m: f64) -> Result<StabilityVerdict> {
    let gm = BuoyancyError::finite("gm_m", gm_m)?;
    let tol = BuoyancyError::non_negative("tol_m", tol_m)?;
    if gm > tol {
        Ok(StabilityVerdict::Stable)
    } else if gm < -tol {
        Ok(StabilityVerdict::Unstable)
    } else {
        Ok(StabilityVerdict::Neutral)
    }
}

/// One-shot small-angle stability of an upright rectangular barge.
///
/// Given beam `beam_m`, draft `draft_m`, and the height of `G` above
/// `B` `bg_m`, computes `BM = B^2 / (12 T)` and returns
/// `GM = BM - BG`. The sign of the result feeds
/// [`stability_verdict`].
///
/// `beam_m` and `draft_m` must be finite and strictly positive; `bg_m`
/// must be finite.
///
/// ```
/// use valenx_buoyancy::stability::barge_metacentric_height;
/// // B = 4, T = 1, BG = 0.5 -> GM = 16/12 - 0.5 = 0.8333...
/// let gm = barge_metacentric_height(4.0, 1.0, 0.5).unwrap();
/// assert!((gm - (16.0 / 12.0 - 0.5)).abs() < 1e-9);
/// assert!(gm > 0.0); // stable
/// ```
pub fn barge_metacentric_height(beam_m: f64, draft_m: f64, bg_m: f64) -> Result<f64> {
    let bm = barge_metacentric_radius(beam_m, draft_m)?;
    metacentric_height(bm, bg_m)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn rectangular_inertia_known_value() {
        // L = 12, B = 6 -> 12 * 216 / 12 = 216.
        let i = rectangular_waterplane_inertia(12.0, 6.0).unwrap();
        assert!((i - 216.0).abs() < EPS);
    }

    #[test]
    fn metacentric_radius_is_inertia_over_volume() {
        let i = 216.0;
        let v = 72.0;
        let bm = metacentric_radius(i, v).unwrap();
        assert!((bm - 3.0).abs() < EPS);
    }

    #[test]
    fn barge_bm_matches_closed_form() {
        // B = 6, T = 2 -> 36 / 24 = 1.5.
        let bm = barge_metacentric_radius(6.0, 2.0).unwrap();
        assert!((bm - 1.5).abs() < EPS);
    }

    #[test]
    fn barge_closed_form_equals_general_inertia_ratio() {
        // GROUND TRUTH: BM = B^2/(12T) must equal I_wp / V for the same
        // box, for any length L. Check across several L.
        let beam = 5.0;
        let draft = 1.5;
        let closed = barge_metacentric_radius(beam, draft).unwrap();
        for &l in &[3.0_f64, 7.5, 20.0, 100.0] {
            let i = rectangular_waterplane_inertia(l, beam).unwrap();
            let v = l * beam * draft; // box displaced volume
            let general = metacentric_radius(i, v).unwrap();
            assert!(
                (general - closed).abs() < 1e-9,
                "L = {l}: general {general} vs closed {closed}"
            );
        }
    }

    #[test]
    fn metacentric_height_is_bm_minus_bg() {
        let gm = metacentric_height(2.5, 1.0).unwrap();
        assert!((gm - 1.5).abs() < EPS);
        // G below B (negative BG) raises GM.
        let gm2 = metacentric_height(2.5, -0.5).unwrap();
        assert!((gm2 - 3.0).abs() < EPS);
    }

    #[test]
    fn barge_metacentric_height_known_value() {
        // B = 4, T = 1, BG = 0.5 -> 16/12 - 0.5 = 0.83333...
        let gm = barge_metacentric_height(4.0, 1.0, 0.5).unwrap();
        assert!((gm - (4.0_f64.powi(2) / 12.0 - 0.5)).abs() < EPS);
    }

    #[test]
    fn verdict_stable_neutral_unstable_by_sign() {
        assert_eq!(
            stability_verdict(1.0, 1e-6).unwrap(),
            StabilityVerdict::Stable
        );
        assert_eq!(
            stability_verdict(0.0, 1e-6).unwrap(),
            StabilityVerdict::Neutral
        );
        assert_eq!(
            stability_verdict(-1.0, 1e-6).unwrap(),
            StabilityVerdict::Unstable
        );
    }

    #[test]
    fn verdict_exact_zero_gm_is_neutral() {
        // GM == 0 (M coincides with G) is the neutral boundary case.
        assert_eq!(
            stability_verdict(0.0, 0.0).unwrap(),
            StabilityVerdict::Neutral
        );
    }

    #[test]
    fn wide_shallow_barge_is_stable_tall_narrow_is_unstable() {
        // Wide & shallow: big BM, easily exceeds a modest BG -> stable.
        let gm_wide = barge_metacentric_height(8.0, 1.0, 1.0).unwrap(); // 64/12 - 1 = 4.33
        assert!(gm_wide > 0.0);
        assert_eq!(
            stability_verdict(gm_wide, 1e-6).unwrap(),
            StabilityVerdict::Stable
        );

        // Narrow & deep with a high G: small BM, large BG -> unstable.
        let gm_narrow = barge_metacentric_height(2.0, 3.0, 2.0).unwrap(); // 4/36 - 2 < 0
        assert!(gm_narrow < 0.0);
        assert_eq!(
            stability_verdict(gm_narrow, 1e-6).unwrap(),
            StabilityVerdict::Unstable
        );
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(rectangular_waterplane_inertia(0.0, 4.0).is_err());
        assert!(metacentric_radius(1.0, 0.0).is_err()); // V must be > 0
        assert!(metacentric_radius(-1.0, 5.0).is_err()); // I must be >= 0
        assert!(barge_metacentric_radius(4.0, 0.0).is_err());
        assert!(metacentric_height(f64::NAN, 1.0).is_err());
        assert!(stability_verdict(1.0, -1e-6).is_err()); // tol must be >= 0
        assert!(stability_verdict(f64::INFINITY, 1e-6).is_err());
    }
}
