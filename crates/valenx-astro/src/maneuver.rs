//! Impulsive orbital-maneuver planning: Δv budgets for the standard
//! transfers used in mission design.
//!
//! Built on the [`crate::orbit3d`] layer, this computes the classic
//! two-body impulsive maneuvers — **Hohmann** and **bi-elliptic**
//! transfers between circular orbits, and **plane changes** — returning
//! the burn magnitudes, total Δv, and transfer time. These are the
//! numbers that size an upper stage or a satellite's propellant for a
//! given mission (e.g. "LEO → GEO costs ~3.9 km/s").
//!
//! The maneuvers are modelled as ideal impulses (instantaneous Δv); the
//! formulas are exact two-body results, so the unit tests validate
//! against both an independent recomputation and known textbook values.

use serde::{Deserialize, Serialize};

use crate::constants::MU_EARTH;
use crate::error::AstroError;

/// Reject a non-physical orbital radius (non-finite or `<= 0`), which
/// would otherwise feed a NaN/Inf into `circular_speed` / `vis_viva`.
fn check_radius(field: &'static str, r: f64) -> Result<(), AstroError> {
    if !r.is_finite() || r <= 0.0 {
        return Err(AstroError::InvalidGuidance(field));
    }
    Ok(())
}

/// The Δv breakdown and timing of a transfer.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Transfer {
    /// First burn magnitude (m/s).
    pub delta_v1: f64,
    /// Second burn magnitude (m/s).
    pub delta_v2: f64,
    /// Third burn magnitude (m/s); zero for a two-burn Hohmann.
    pub delta_v3: f64,
    /// Total Δv (m/s).
    pub total_delta_v: f64,
    /// Time of flight along the transfer (s).
    pub transfer_time: f64,
}

impl Transfer {
    fn two_burn(dv1: f64, dv2: f64, time: f64) -> Self {
        Self {
            delta_v1: dv1,
            delta_v2: dv2,
            delta_v3: 0.0,
            total_delta_v: dv1 + dv2,
            transfer_time: time,
        }
    }
}

/// Circular orbital speed at radius `r` (m).
fn circular_speed(r: f64) -> f64 {
    (MU_EARTH / r).sqrt()
}

/// Speed at radius `r` on an ellipse of semi-major axis `a` (vis-viva).
fn vis_viva(r: f64, a: f64) -> f64 {
    (MU_EARTH * (2.0 / r - 1.0 / a)).sqrt()
}

/// Hohmann transfer between two **circular, coplanar** orbits of radii
/// `r1` and `r2` (m). The minimum-Δv two-impulse transfer for most
/// radius ratios.
///
/// `delta_v1` raises (or lowers) the orbit onto the transfer ellipse at
/// `r1`; `delta_v2` circularises at `r2`. The transfer time is half the
/// period of the transfer ellipse.
///
/// # Errors
///
/// Returns [`AstroError::InvalidGuidance`] if either radius is
/// non-finite or non-positive (which would otherwise produce a NaN/Inf
/// Δv).
pub fn hohmann_transfer(r1: f64, r2: f64) -> Result<Transfer, AstroError> {
    check_radius("hohmann r1 must be > 0", r1)?;
    check_radius("hohmann r2 must be > 0", r2)?;
    let v1 = circular_speed(r1);
    let v2 = circular_speed(r2);
    let a_t = 0.5 * (r1 + r2);
    let v_peri = vis_viva(r1, a_t); // speed on the ellipse at r1
    let v_apo = vis_viva(r2, a_t); // speed on the ellipse at r2

    let dv1 = (v_peri - v1).abs();
    let dv2 = (v2 - v_apo).abs();
    let time = std::f64::consts::PI * (a_t.powi(3) / MU_EARTH).sqrt();
    Ok(Transfer::two_burn(dv1, dv2, time))
}

/// Bi-elliptic transfer between two circular orbits via an intermediate
/// apoapsis radius `r_b` (m). For very large radius ratios (roughly
/// `r2/r1 > 11.94`) a bi-elliptic transfer with a high enough `r_b` can
/// beat the Hohmann total Δv, at the cost of much longer transfer time.
///
/// # Errors
///
/// Returns [`AstroError::InvalidGuidance`] if any of `r1`, `r2`, `r_b`
/// is non-finite or non-positive.
pub fn bielliptic_transfer(r1: f64, r2: f64, r_b: f64) -> Result<Transfer, AstroError> {
    check_radius("bielliptic r1 must be > 0", r1)?;
    check_radius("bielliptic r2 must be > 0", r2)?;
    check_radius("bielliptic r_b must be > 0", r_b)?;
    let v1 = circular_speed(r1);
    let v2 = circular_speed(r2);

    // First ellipse: r1 -> r_b.
    let a1 = 0.5 * (r1 + r_b);
    let dv1 = (vis_viva(r1, a1) - v1).abs();

    // Second ellipse: r_b -> r2 (burn at r_b).
    let a2 = 0.5 * (r2 + r_b);
    let dv2 = (vis_viva(r_b, a2) - vis_viva(r_b, a1)).abs();

    // Circularise at r2.
    let dv3 = (v2 - vis_viva(r2, a2)).abs();

    let time =
        std::f64::consts::PI * ((a1.powi(3) / MU_EARTH).sqrt() + (a2.powi(3) / MU_EARTH).sqrt());

    Ok(Transfer {
        delta_v1: dv1,
        delta_v2: dv2,
        delta_v3: dv3,
        total_delta_v: dv1 + dv2 + dv3,
        transfer_time: time,
    })
}

/// Δv (m/s) for a pure inclination change of `delta_inclination` (rad)
/// performed at orbital speed `v`: `Δv = 2 v sin(Δi/2)`.
pub fn plane_change_dv(orbital_speed: f64, delta_inclination: f64) -> f64 {
    2.0 * orbital_speed * (delta_inclination / 2.0).sin()
}

/// Δv (m/s) for a circular orbit of radius `r` (m) changing inclination
/// by `delta_inclination` (rad) — `plane_change_dv` at circular speed.
///
/// # Errors
///
/// Returns [`AstroError::InvalidGuidance`] if `r` is non-finite or
/// non-positive, which would otherwise make the circular speed
/// `√(μ/r)` — and hence the returned Δv — a silent `NaN`/`Inf`.
pub fn circular_plane_change_dv(r: f64, delta_inclination: f64) -> Result<f64, AstroError> {
    check_radius("circular_plane_change_dv r must be > 0", r)?;
    Ok(plane_change_dv(circular_speed(r), delta_inclination))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::R_EARTH;

    #[test]
    fn hohmann_leo_to_geo_matches_known_budget() {
        let r1 = R_EARTH + 300_000.0; // ~300 km LEO
        let r2 = 42_164_000.0; // GEO radius
        let t = hohmann_transfer(r1, r2).expect("valid radii");

        // Textbook LEO->GEO Hohmann is ~3.9 km/s total, split ~2.4 + ~1.5.
        assert!((t.delta_v1 - 2_425.0).abs() < 40.0, "dv1 {}", t.delta_v1);
        assert!((t.delta_v2 - 1_467.0).abs() < 40.0, "dv2 {}", t.delta_v2);
        assert!(
            (t.total_delta_v - 3_892.0).abs() < 60.0,
            "total {}",
            t.total_delta_v
        );
        // Transfer time ~5.27 h.
        assert!(
            (t.transfer_time - 18_990.0).abs() < 200.0,
            "t {}",
            t.transfer_time
        );
    }

    #[test]
    fn hohmann_is_symmetric_in_total_dv() {
        // Going up then back down costs the same total Δv.
        let r1 = R_EARTH + 400_000.0;
        let r2 = R_EARTH + 2_000_000.0;
        let up = hohmann_transfer(r1, r2).expect("valid radii");
        let down = hohmann_transfer(r2, r1).expect("valid radii");
        assert!((up.total_delta_v - down.total_delta_v).abs() < 1e-6);
        // Same-radius transfer costs nothing.
        let none = hohmann_transfer(r1, r1).expect("valid radii");
        assert!(none.total_delta_v < 1e-9);
    }

    #[test]
    fn bielliptic_beats_hohmann_for_large_ratio() {
        // For a radius ratio well above ~11.94, a high-apoapsis
        // bi-elliptic transfer has a lower total Δv than the Hohmann.
        let r1 = R_EARTH + 300_000.0;
        let r2 = r1 * 15.0;
        let h = hohmann_transfer(r1, r2).expect("valid radii");
        let b = bielliptic_transfer(r1, r2, r2 * 5.0).expect("valid radii");
        assert!(
            b.total_delta_v < h.total_delta_v,
            "bielliptic {} should beat hohmann {}",
            b.total_delta_v,
            h.total_delta_v
        );
        // ...but it takes much longer.
        assert!(b.transfer_time > h.transfer_time);
    }

    #[test]
    fn plane_change_is_expensive_at_leo_speed() {
        // A 28.5° plane change in LEO costs a large Δv (~3.8 km/s).
        let r = R_EARTH + 300_000.0;
        let dv = circular_plane_change_dv(r, 28.5_f64.to_radians()).expect("valid radius");
        let v = circular_speed(r);
        let expected = 2.0 * v * (28.5_f64.to_radians() / 2.0).sin();
        assert!((dv - expected).abs() < 1e-9);
        assert!((3_700.0..3_900.0).contains(&dv), "dv {dv}");
    }

    #[test]
    fn zero_plane_change_is_free() {
        let dv = circular_plane_change_dv(R_EARTH + 500_000.0, 0.0).expect("valid radius");
        assert!(dv.abs() < 1e-12);
    }

    #[test]
    fn transfers_reject_nonphysical_radii_instead_of_nan() {
        // r <= 0 / non-finite used to feed straight into circular_speed /
        // vis_viva and emit a NaN/Inf Δv. They must now error.
        let r_ok = R_EARTH + 300_000.0;
        assert!(hohmann_transfer(0.0, r_ok).is_err());
        assert!(hohmann_transfer(r_ok, -1.0).is_err());
        assert!(hohmann_transfer(f64::NAN, r_ok).is_err());
        assert!(bielliptic_transfer(r_ok, r_ok * 15.0, 0.0).is_err());
        assert!(bielliptic_transfer(-1.0, r_ok * 15.0, r_ok * 50.0).is_err());
        // The circular plane-change sibling has the same √(μ/r) hazard.
        assert!(matches!(
            circular_plane_change_dv(0.0, 0.5),
            Err(AstroError::InvalidGuidance(_))
        ));
        assert!(circular_plane_change_dv(-1.0, 0.5).is_err());
        assert!(circular_plane_change_dv(f64::NAN, 0.5).is_err());
        // A valid transfer still produces finite numbers.
        let t = hohmann_transfer(r_ok, r_ok * 2.0).expect("valid");
        assert!(t.total_delta_v.is_finite() && t.transfer_time.is_finite());
    }
}
