//! Winkler-Bach bending stress.
//!
//! With the neutral-axis radius `R_n` and eccentricity `e = R - R_n` known
//! (see [`crate::section`]), the bending stress at a fibre of radius `r` is
//!
//! ```text
//! sigma(r) = M * y / (A * e * (R_n - y)) = M * y / (A * e * r)
//! ```
//!
//! where `y = R_n - r` is the distance from the neutral axis measured
//! *toward* the centre of curvature, so `R_n - y = r`. The inner fibre has
//! `y > 0` and the outer fibre has `y < 0`.
//!
//! ## Sign convention
//!
//! `M` is taken positive when the moment tends to *decrease* the curvature
//! (straighten the beam). With that convention a positive `M` puts the inner
//! fibre in tension (`sigma_inner > 0`) and the outer fibre in compression
//! (`sigma_outer < 0`), and the magnitude is always larger on the inside.
//! This is the crane-hook / chain-link stress concentration on the inner
//! edge.

use crate::error::{require_finite, require_positive, CurvedBeamError};
use crate::section::{Section, SectionProps};

/// Bending stress at an arbitrary radius `r` within the section.
///
/// # Errors
///
/// Returns [`CurvedBeamError::NotFinite`] if `moment` is not finite, and
/// [`CurvedBeamError::NotPositive`] if `r` is not strictly positive.
pub fn stress_at_radius(section: &Section, moment: f64, r: f64) -> Result<f64, CurvedBeamError> {
    let m = require_finite("moment", moment)?;
    let r = require_positive("r", r)?;
    let p = section.props();
    let y = p.r_neutral - r;
    Ok(m * y / (p.area * p.eccentricity * r))
}

/// Stress at the inner fibre (`r = r_i`).
///
/// With the crate sign convention a positive `moment` yields a positive
/// (tensile) result here.
///
/// # Errors
///
/// Returns [`CurvedBeamError::NotFinite`] if `moment` is not finite.
pub fn stress_inner(section: &Section, moment: f64) -> Result<f64, CurvedBeamError> {
    stress_at_radius(section, moment, section.r_inner())
}

/// Stress at the outer fibre (`r = r_o`).
///
/// With the crate sign convention a positive `moment` yields a negative
/// (compressive) result here.
///
/// # Errors
///
/// Returns [`CurvedBeamError::NotFinite`] if `moment` is not finite.
pub fn stress_outer(section: &Section, moment: f64) -> Result<f64, CurvedBeamError> {
    stress_at_radius(section, moment, section.r_outer())
}

/// Inner- and outer-fibre stresses bundled together with the section
/// properties used to compute them.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FibreStresses {
    /// Inner-fibre stress (tensile-positive for positive moment).
    pub inner: f64,
    /// Outer-fibre stress (compressive-negative for positive moment).
    pub outer: f64,
    /// The section properties that drove the calculation.
    pub props: SectionProps,
}

/// Compute both extreme-fibre stresses in one shot.
///
/// # Errors
///
/// Returns [`CurvedBeamError::NotFinite`] if `moment` is not finite.
pub fn fibre_stresses(section: &Section, moment: f64) -> Result<FibreStresses, CurvedBeamError> {
    Ok(FibreStresses {
        inner: stress_inner(section, moment)?,
        outer: stress_outer(section, moment)?,
        props: section.props(),
    })
}

/// Straight (Euler-Bernoulli) bending stress `sigma = M y / I` at distance
/// `y` from the centroidal axis, for the same section.
///
/// This is the large-radius limit the curved-beam formula must converge to;
/// it exists so tests (and callers) can quantify how far a given beam departs
/// from straight-beam theory. Here `y` is measured from the centroid and is
/// positive toward the centre of curvature, matching the inner fibre being
/// positive under a positive moment.
///
/// # Errors
///
/// Returns [`CurvedBeamError::NotFinite`] if `moment` or `y` is not finite.
pub fn straight_beam_stress(
    section: &Section,
    moment: f64,
    y: f64,
) -> Result<f64, CurvedBeamError> {
    let m = require_finite("moment", moment)?;
    let y = require_finite("y", y)?;
    Ok(m * y / section.inertia())
}

/// Inner-fibre stress predicted by straight-beam theory, `M (R - r_i) / I`.
///
/// Distance from the centroid to the inner fibre is `R - r_i`, taken
/// positive (toward the centre of curvature).
///
/// # Errors
///
/// Returns [`CurvedBeamError::NotFinite`] if `moment` is not finite.
pub fn straight_beam_inner(section: &Section, moment: f64) -> Result<f64, CurvedBeamError> {
    let y = section.r_centroid() - section.r_inner();
    straight_beam_stress(section, moment, y)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a mildly curved rectangle for the common cases.
    fn rect(r_i: f64, r_o: f64) -> Section {
        Section::rectangle(1.0, r_i, r_o).unwrap()
    }

    #[test]
    fn inner_fibre_sign_follows_moment() {
        // Ground truth: positive (straightening) moment -> inner tension.
        let s = rect(4.0, 6.0);
        let si = stress_inner(&s, 100.0).unwrap();
        assert!(si > 0.0, "inner fibre should be tensile for +M, got {si}");
        let si_neg = stress_inner(&s, -100.0).unwrap();
        assert!(si_neg < 0.0);
        // Exactly antisymmetric in M.
        assert!((si + si_neg).abs() < 1e-9);
    }

    #[test]
    fn outer_fibre_opposes_inner() {
        let s = rect(4.0, 6.0);
        let si = stress_inner(&s, 100.0).unwrap();
        let so = stress_outer(&s, 100.0).unwrap();
        assert!(so < 0.0, "outer fibre should be compressive for +M");
        // Inner magnitude exceeds outer for a curved beam.
        assert!(si.abs() > so.abs());
    }

    #[test]
    fn hand_worked_winkler_value() {
        // b=1, r_i=4, r_o=6, M=100.
        // A = 2, R = 5, R_n = 2/ln(1.5) = 4.9326177914..., e = R - R_n.
        // sigma_inner = M (R_n - r_i) / (A e r_i).
        let s = rect(4.0, 6.0);
        let r_n = 2.0 / 1.5_f64.ln();
        let e = 5.0 - r_n;
        let expected_inner = 100.0 * (r_n - 4.0) / (2.0 * e * 4.0);
        let expected_outer = 100.0 * (r_n - 6.0) / (2.0 * e * 6.0);
        assert!((stress_inner(&s, 100.0).unwrap() - expected_inner).abs() < 1e-9);
        assert!((stress_outer(&s, 100.0).unwrap() - expected_outer).abs() < 1e-9);
        // Independent decimal pins (full-precision Winkler-Bach evaluation).
        assert!((stress_inner(&s, 100.0).unwrap() - 172.97899697648992).abs() < 1e-9);
        assert!((stress_outer(&s, 100.0).unwrap() - (-131.98599798432662)).abs() < 1e-9);
    }

    #[test]
    fn stress_at_radius_matches_fibre_helpers() {
        let s = rect(4.0, 6.0);
        let a = stress_at_radius(&s, 50.0, s.r_inner()).unwrap();
        let b = stress_inner(&s, 50.0).unwrap();
        assert!((a - b).abs() < 1e-12);
    }

    #[test]
    fn neutral_axis_has_zero_stress() {
        // Ground truth: at r = R_n the fibre stress vanishes (y = 0).
        let s = rect(4.0, 6.0);
        let r_n = s.r_neutral();
        let sigma = stress_at_radius(&s, 1234.0, r_n).unwrap();
        assert!(
            sigma.abs() < 1e-9,
            "stress at neutral axis must be ~0, got {sigma}"
        );
    }

    #[test]
    fn straight_beam_limit_large_radius_inner() {
        // Ground truth: as R grows, Winkler inner stress -> M (R - r_i) / I.
        // Use a large radius with fixed depth h = 2, so curvature is gentle.
        let r_i = 1000.0;
        let r_o = 1002.0;
        let s = rect(r_i, r_o);
        let curved = stress_inner(&s, 100.0).unwrap();
        let straight = straight_beam_inner(&s, 100.0).unwrap();
        // Relative error should be tiny (order (h/R)^2 ~ 1e-6).
        let rel = (curved - straight).abs() / straight.abs();
        assert!(rel < 1e-3, "rel err {rel} too large for large-radius limit");
    }

    #[test]
    fn straight_beam_limit_circle() {
        // Same convergence for a solid-circular section.
        let s = Section::circle(1.0, 1000.0).unwrap();
        let curved = stress_inner(&s, 100.0).unwrap();
        let straight = straight_beam_inner(&s, 100.0).unwrap();
        let rel = (curved - straight).abs() / straight.abs();
        assert!(rel < 1e-3, "circle large-radius rel err {rel}");
    }

    #[test]
    fn sharply_curved_beam_departs_from_straight() {
        // For a tightly curved beam the two theories disagree substantially;
        // Winkler predicts a markedly higher inner-fibre stress.
        let s = rect(4.0, 6.0);
        let curved = stress_inner(&s, 100.0).unwrap();
        let straight = straight_beam_inner(&s, 100.0).unwrap();
        // Straight-beam: M (R - r_i)/I = 100 * 1 / (1*8/12) = 150.
        assert!((straight - 150.0).abs() < 1e-9);
        // Winkler inner (~173 here) is markedly higher; differs by > 10%.
        assert!(
            curved > straight,
            "curved inner stress should exceed straight"
        );
        let rel = (curved - straight).abs() / straight.abs();
        assert!(rel > 0.1, "expected large departure, got rel {rel}");
    }

    #[test]
    fn fibre_stresses_bundle_consistent() {
        let s = rect(4.0, 6.0);
        let f = fibre_stresses(&s, 100.0).unwrap();
        assert!((f.inner - stress_inner(&s, 100.0).unwrap()).abs() < 1e-12);
        assert!((f.outer - stress_outer(&s, 100.0).unwrap()).abs() < 1e-12);
        assert!((f.props.r_neutral - s.r_neutral()).abs() < 1e-12);
    }

    #[test]
    fn zero_moment_gives_zero_stress() {
        let s = rect(4.0, 6.0);
        assert!(stress_inner(&s, 0.0).unwrap().abs() < 1e-12);
        assert!(stress_outer(&s, 0.0).unwrap().abs() < 1e-12);
    }

    #[test]
    fn rejects_non_finite_moment() {
        let s = rect(4.0, 6.0);
        let err = stress_inner(&s, f64::INFINITY).unwrap_err();
        assert_eq!(err.code(), "curvedbeam.not_finite");
    }

    #[test]
    fn rejects_non_positive_radius() {
        let s = rect(4.0, 6.0);
        let err = stress_at_radius(&s, 10.0, -1.0).unwrap_err();
        assert_eq!(err.code(), "curvedbeam.not_positive");
    }
}
