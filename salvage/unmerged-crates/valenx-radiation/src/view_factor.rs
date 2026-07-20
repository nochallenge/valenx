//! View factors (configuration factors) and two-surface enclosure
//! exchange.
//!
//! The view factor `F12` is the fraction of radiation leaving surface 1
//! that is intercepted by surface 2. Two algebraic relations let one
//! deduce unknown factors from known ones:
//!
//! Reciprocity:
//!
//! ```text
//! A1 * F12 = A2 * F21
//! ```
//!
//! Summation (an enclosure of `N` surfaces, each surface `i` seeing the
//! whole enclosure including possibly itself):
//!
//! ```text
//! sum_j F_ij = 1
//! ```
//!
//! For a two-surface enclosure of *diffuse gray* surfaces with areas
//! `A1, A2`, emissivities `eps1, eps2`, and view factor `F12`, the net
//! radiant exchange from 1 to 2 is the textbook series-resistance
//! result
//!
//! ```text
//!            sigma * (T1^4 - T2^4)
//! Q12 = -----------------------------------------------
//!        (1-eps1)/(eps1 A1) + 1/(A1 F12) + (1-eps2)/(eps2 A2)
//! ```
//!
//! Ground-truth checks: reciprocity is exact by construction;
//! two black bodies (`eps1 = eps2 = 1`) reduce the resistance network
//! to `Q12 = A1 F12 sigma (T1^4 - T2^4)`; two infinite parallel gray
//! plates (`A1 = A2`, `F12 = 1`) reduce to the classic
//! `q'' = sigma (T1^4 - T2^4) / (1/eps1 + 1/eps2 - 1)`.

use crate::error::{
    check_area, check_temperature, check_unit_interval, finite, RadiationError, Result,
};
use crate::stefan_boltzmann::SIGMA;

/// Reciprocity: given `A1`, `F12`, `A2`, return `F21 = A1 F12 / A2`.
///
/// # Arguments
///
/// `area1` — area of surface 1 (`> 0`);
/// `f12` — view factor from 1 to 2, in `[0, 1]`;
/// `area2` — area of surface 2 (`> 0`).
///
/// # Returns
///
/// `F21`, the view factor from 2 to 1.
///
/// # Errors
///
/// Domain-checks all three inputs. Returns
/// [`RadiationError::Inconsistent`] if the implied `F21` would exceed 1
/// (which signals physically incompatible areas / `F12`).
///
/// # Example
///
/// ```
/// use valenx_radiation::view_factor::reciprocal_view_factor;
/// // A small surface (A1 = 1) fully seeing a large one (A2 = 4):
/// // F12 = 1  =>  F21 = A1/A2 = 0.25.
/// let f21 = reciprocal_view_factor(1.0, 1.0, 4.0).unwrap();
/// assert!((f21 - 0.25).abs() < 1e-12);
/// ```
pub fn reciprocal_view_factor(area1: f64, f12: f64, area2: f64) -> Result<f64> {
    let a1 = check_area("area1", area1)?;
    let f = check_unit_interval("f12", f12)?;
    let a2 = check_area("area2", area2)?;
    let f21 = a1 * f / a2;
    if f21 > 1.0 + 1e-12 {
        return Err(RadiationError::Inconsistent {
            reason: format!("reciprocity yields F21 = {f21} > 1 for A1={a1}, F12={f}, A2={a2}"),
        });
    }
    // Clamp a tiny floating overshoot of the [0,1] domain back to 1.
    Ok(f21.min(1.0))
}

/// Summation rule for a two-surface enclosure: the self-view factor
/// `F11 = 1 - F12`.
///
/// Returns the fraction of radiation leaving surface 1 that returns to
/// surface 1 (nonzero only for a concave surface). For a flat or convex
/// surface `F11 = 0`, so `F12 = 1`.
///
/// # Errors
///
/// Returns [`RadiationError::OutOfDomain`] if `f12` is outside `[0, 1]`
/// or non-finite.
pub fn self_view_factor_two_surface(f12: f64) -> Result<f64> {
    let f = check_unit_interval("f12", f12)?;
    Ok(1.0 - f)
}

/// Closing view factor of an `N`-surface enclosure: from
/// `sum_j F_ij = 1`, the last unknown factor is
/// `F_last = 1 - sum(known_factors)`.
///
/// # Arguments
///
/// `known_factors` — the already-known `F_ij` for surface `i`, each in
/// `[0, 1]`; their sum must not exceed 1.
///
/// # Errors
///
/// Returns [`RadiationError::OutOfDomain`] if any entry is outside
/// `[0, 1]` or non-finite, and [`RadiationError::Inconsistent`] if the
/// known factors already sum to more than 1.
pub fn closing_view_factor(known_factors: &[f64]) -> Result<f64> {
    let mut sum = 0.0;
    for (i, &f) in known_factors.iter().enumerate() {
        let v = finite("known_factor", f)?;
        if !(0.0..=1.0).contains(&v) {
            return Err(RadiationError::OutOfDomain {
                what: "known_factor",
                value: v,
                reason: "each view factor must lie in [0, 1]",
            });
        }
        let _ = i;
        sum += v;
    }
    if sum > 1.0 + 1e-12 {
        return Err(RadiationError::Inconsistent {
            reason: format!("known view factors sum to {sum} > 1"),
        });
    }
    Ok((1.0 - sum).max(0.0))
}

/// Net radiant exchange `Q12` (watts) between two diffuse gray surfaces
/// forming an enclosure, via the series-resistance network.
///
/// ```text
///            sigma (T1^4 - T2^4)
/// Q12 = ------------------------------------------------
///        R_surf1 + R_space + R_surf2
/// ```
///
/// with surface resistances `R_surf = (1-eps)/(eps A)` and space
/// resistance `R_space = 1/(A1 F12)`.
///
/// # Arguments
///
/// `area1`, `area2` — surface areas in `m^2` (`> 0`);
/// `eps1`, `eps2` — emissivities in `(0, 1]` (a zero emissivity makes
/// the surface resistance infinite, so it is rejected here);
/// `f12` — view factor from 1 to 2 in `(0, 1]` (a zero factor means the
/// surfaces do not see each other, giving infinite space resistance);
/// `t1_k`, `t2_k` — temperatures in kelvin (`>= 0`).
///
/// # Returns
///
/// Net power from surface 1 to surface 2 in watts; positive when
/// `T1 > T2`.
///
/// # Errors
///
/// Domain-checks every argument; emissivities and `f12` must be
/// strictly positive (returns [`RadiationError::OutOfDomain`]
/// otherwise) so the resistances stay finite.
#[allow(clippy::too_many_arguments)]
pub fn two_surface_exchange(
    area1: f64,
    eps1: f64,
    area2: f64,
    eps2: f64,
    f12: f64,
    t1_k: f64,
    t2_k: f64,
) -> Result<f64> {
    let a1 = check_area("area1", area1)?;
    let a2 = check_area("area2", area2)?;
    let e1 = positive_unit("eps1", eps1)?;
    let e2 = positive_unit("eps2", eps2)?;
    let f = positive_unit("f12", f12)?;
    let t1 = check_temperature("t1_k", t1_k)?;
    let t2 = check_temperature("t2_k", t2_k)?;

    let r_surf1 = (1.0 - e1) / (e1 * a1);
    let r_space = 1.0 / (a1 * f);
    let r_surf2 = (1.0 - e2) / (e2 * a2);
    let total_r = r_surf1 + r_space + r_surf2;

    Ok(SIGMA * (t1.powi(4) - t2.powi(4)) / total_r)
}

/// Net radiant *flux* between two infinite parallel gray plates,
/// `q'' = sigma (T1^4 - T2^4) / (1/eps1 + 1/eps2 - 1)` in `W / m^2`.
///
/// This is the `A1 = A2`, `F12 = 1` limit of [`two_surface_exchange`]
/// divided by the (common) area, and is the canonical textbook result
/// for radiation between large parallel surfaces.
///
/// # Errors
///
/// Emissivities must be strictly positive and `<= 1`; temperatures
/// `>= 0`.
pub fn parallel_plates_flux(eps1: f64, eps2: f64, t1_k: f64, t2_k: f64) -> Result<f64> {
    let e1 = positive_unit("eps1", eps1)?;
    let e2 = positive_unit("eps2", eps2)?;
    let t1 = check_temperature("t1_k", t1_k)?;
    let t2 = check_temperature("t2_k", t2_k)?;
    let denom = 1.0 / e1 + 1.0 / e2 - 1.0;
    Ok(SIGMA * (t1.powi(4) - t2.powi(4)) / denom)
}

/// Like [`check_unit_interval`](crate::error::check_unit_interval) but
/// rejects exactly zero — used where a zero value would make a
/// resistance infinite.
fn positive_unit(what: &'static str, x: f64) -> Result<f64> {
    let v = check_unit_interval(what, x)?;
    if v <= 0.0 {
        return Err(RadiationError::OutOfDomain {
            what,
            value: v,
            reason: "must lie in (0, 1] (zero makes a resistance infinite)",
        });
    }
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn reciprocity_is_exact() {
        // A1 F12 == A2 F21 by construction.
        let (a1, f12, a2) = (2.0, 0.3, 5.0);
        let f21 = reciprocal_view_factor(a1, f12, a2).unwrap();
        assert!((a1 * f12 - a2 * f21).abs() < EPS, "f21={f21}");
    }

    #[test]
    fn small_surface_in_large_enclosure() {
        // Convex A1 fully enclosed by A2: F12 = 1, so F21 = A1/A2.
        let f21 = reciprocal_view_factor(1.0, 1.0, 4.0).unwrap();
        assert!((f21 - 0.25).abs() < 1e-12, "got {f21}");
    }

    #[test]
    fn reciprocity_rejects_impossible_geometry() {
        // A1 F12 / A2 = 4*1/1 = 4 > 1 is non-physical.
        let e = reciprocal_view_factor(4.0, 1.0, 1.0).unwrap_err();
        assert_eq!(e.code(), "radiation.inconsistent");
    }

    #[test]
    fn summation_two_surface() {
        assert!((self_view_factor_two_surface(1.0).unwrap()).abs() < EPS);
        assert!((self_view_factor_two_surface(0.3).unwrap() - 0.7).abs() < EPS);
    }

    #[test]
    fn closing_factor_sums_to_one() {
        // Enclosure surface seeing three others 0.2, 0.3, 0.1 => closes at 0.4.
        let last = closing_view_factor(&[0.2, 0.3, 0.1]).unwrap();
        assert!((last - 0.4).abs() < EPS, "got {last}");
        // Already-complete set closes at 0.
        let z = closing_view_factor(&[0.6, 0.4]).unwrap();
        assert!(z.abs() < EPS, "got {z}");
    }

    #[test]
    fn closing_factor_rejects_oversum() {
        assert!(closing_view_factor(&[0.8, 0.5]).is_err());
        assert!(closing_view_factor(&[1.2]).is_err());
        assert!(closing_view_factor(&[f64::NAN]).is_err());
    }

    #[test]
    fn two_black_bodies_reduce_to_a1_f12_eb() {
        // eps1 = eps2 = 1 kills the surface resistances:
        // Q12 = A1 F12 sigma (T1^4 - T2^4).
        let q = two_surface_exchange(2.0, 1.0, 3.0, 1.0, 0.5, 600.0, 300.0).unwrap();
        let truth = 2.0 * 0.5 * SIGMA * (600.0f64.powi(4) - 300.0f64.powi(4));
        assert!((q - truth).abs() < 1e-6, "got {q}");
    }

    #[test]
    fn two_surface_isothermal_zero() {
        let q = two_surface_exchange(2.0, 0.8, 2.0, 0.9, 1.0, 500.0, 500.0).unwrap();
        assert!(q.abs() < EPS, "got {q}");
    }

    #[test]
    fn enclosure_limit_matches_parallel_plates_formula() {
        // A1 = A2 = A, F12 = 1  =>  Q12 / A must equal the
        // parallel-plate flux sigma(T1^4-T2^4)/(1/e1 + 1/e2 - 1).
        let (a, e1, e2, t1, t2) = (3.0, 0.6, 0.8, 800.0, 500.0);
        let q = two_surface_exchange(a, e1, a, e2, 1.0, t1, t2).unwrap();
        let flux = parallel_plates_flux(e1, e2, t1, t2).unwrap();
        assert!((q / a - flux).abs() < 1e-9, "q/A={} flux={flux}", q / a);
    }

    #[test]
    fn parallel_plates_black_recovers_stefan_boltzmann() {
        // eps1 = eps2 = 1 => denom = 1 => q'' = sigma (T1^4 - T2^4).
        let flux = parallel_plates_flux(1.0, 1.0, 1000.0, 0.0).unwrap();
        let truth = SIGMA * 1000.0f64.powi(4);
        assert!((flux - truth).abs() < 1e-6, "got {flux}");
    }

    #[test]
    fn parallel_plates_hand_worked() {
        // e1=e2=0.5 => denom = 2 + 2 - 1 = 3.
        // T1=600, T2=400: 600^4-400^4 = 1.296e11 - 2.56e10 = 1.04e11.
        // q'' = sigma * 1.04e11 / 3 = 1965.73... W/m^2.
        let flux = parallel_plates_flux(0.5, 0.5, 600.0, 400.0).unwrap();
        let truth = SIGMA * (600.0f64.powi(4) - 400.0f64.powi(4)) / 3.0;
        assert!((flux - truth).abs() < 1e-9, "got {flux}");
        assert!((flux - 1965.73).abs() < 1e-1, "got {flux}");
    }

    #[test]
    fn two_surface_rejects_zero_emissivity_or_factor() {
        assert!(two_surface_exchange(1.0, 0.0, 1.0, 0.5, 1.0, 400.0, 300.0).is_err());
        assert!(two_surface_exchange(1.0, 0.5, 1.0, 0.5, 0.0, 400.0, 300.0).is_err());
        assert!(parallel_plates_flux(0.0, 0.5, 400.0, 300.0).is_err());
        assert!(two_surface_exchange(0.0, 0.5, 1.0, 0.5, 1.0, 400.0, 300.0).is_err());
    }
}
