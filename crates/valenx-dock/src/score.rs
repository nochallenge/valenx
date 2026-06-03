//! Vina scoring function — five terms + weighted sum.
//!
//! Weights from Trott & Olson, J. Comput. Chem. 31 (2010) 455-461,
//! supplementary Table S1. All terms are pure functions of the
//! surface distance `d = ||r_i - r_j|| - vdw_i - vdw_j`.

use nalgebra::Vector3;

use crate::atom_type::Ad4AtomType;
use crate::PAIR_CUTOFF;

/// Vina inter-atomic weights, frozen at the published values.
pub mod weights {
    /// Steep gaussian centered at 0 Å surface separation.
    pub const GAUSS1: f64 = -0.035579;
    /// Wider gaussian centered at 3 Å surface separation.
    pub const GAUSS2: f64 = -0.005156;
    /// Soft-wall repulsion below 0 Å surface separation.
    pub const REPULSION: f64 = 0.840245;
    /// Hydrophobic contact bonus (carbons + halogens).
    pub const HYDROPHOBIC: f64 = -0.035069;
    /// Hydrogen-bond bonus (donor-acceptor).
    pub const HBOND: f64 = -0.587439;
    /// Torsion-entropy denominator weight.
    pub const N_ROT: f64 = 0.05846;
}

/// Surface distance: cartesian distance minus the sum of VDW radii.
pub fn surface_distance(r_i: Vector3<f64>, r_j: Vector3<f64>, vdw_i: f64, vdw_j: f64) -> f64 {
    (r_i - r_j).norm() - vdw_i - vdw_j
}

/// True when the centre-to-centre distance is within the Vina pair cutoff.
pub fn within_cutoff(r_i: Vector3<f64>, r_j: Vector3<f64>) -> bool {
    (r_i - r_j).norm_squared() <= PAIR_CUTOFF * PAIR_CUTOFF
}

/// Vina's narrow attractive gaussian: exp(-(d-0)^2 / 0.5^2).
/// Surface distance `d` in Å.
pub fn gauss1(d: f64) -> f64 {
    let sigma = 0.5;
    (-(d * d) / (sigma * sigma)).exp()
}

/// Vina's broad attractive gaussian: exp(-(d-3)^2 / 2.0^2).
pub fn gauss2(d: f64) -> f64 {
    let sigma = 2.0;
    let dd = d - 3.0;
    (-(dd * dd) / (sigma * sigma)).exp()
}

/// Soft-wall repulsion: d^2 when d < 0, else 0. Models steric clash.
pub fn repulsion(d: f64) -> f64 {
    if d < 0.0 {
        d * d
    } else {
        0.0
    }
}

/// Hydrophobic distance factor: piecewise linear from 1 at d ≤ 0.5
/// to 0 at d ≥ 1.5, linearly interpolated between.
pub fn hydrophobic_factor(d: f64) -> f64 {
    if d <= 0.5 {
        1.0
    } else if d >= 1.5 {
        0.0
    } else {
        1.0 - (d - 0.5)
    }
}

/// Hydrophobic pair contribution: factor * (both atoms hydrophobic).
pub fn hydrophobic_pair(a: Ad4AtomType, b: Ad4AtomType, d: f64) -> f64 {
    if a.props().is_hydrophobic && b.props().is_hydrophobic {
        hydrophobic_factor(d)
    } else {
        0.0
    }
}

/// Hydrogen-bond distance factor: piecewise linear from 1 at d ≤ -0.7
/// to 0 at d ≥ 0, linearly interpolated between.
pub fn hbond_factor(d: f64) -> f64 {
    if d <= -0.7 {
        1.0
    } else if d >= 0.0 {
        0.0
    } else {
        d / -0.7
    }
}

/// Hydrogen-bond pair contribution. Requires one donor + one acceptor.
pub fn hbond_pair(a: Ad4AtomType, b: Ad4AtomType, d: f64) -> f64 {
    let pa = a.props();
    let pb = b.props();
    let donor_acceptor = (pa.is_donor && pb.is_acceptor) || (pa.is_acceptor && pb.is_donor);
    if donor_acceptor {
        hbond_factor(d)
    } else {
        0.0
    }
}

/// Total Vina pair score: weighted sum of all five terms.
pub fn pair_score(a: Ad4AtomType, b: Ad4AtomType, d: f64) -> f64 {
    weights::GAUSS1 * gauss1(d)
        + weights::GAUSS2 * gauss2(d)
        + weights::REPULSION * repulsion(d)
        + weights::HYDROPHOBIC * hydrophobic_pair(a, b, d)
        + weights::HBOND * hbond_pair(a, b, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn surface_distance_centre_to_centre_minus_radii() {
        let a = Vector3::new(0.0, 0.0, 0.0);
        let b = Vector3::new(5.0, 0.0, 0.0);
        let d = surface_distance(a, b, 1.9, 1.9);
        assert!((d - (5.0 - 3.8)).abs() < 1e-12);
    }

    #[test]
    fn surface_distance_can_be_negative_inside_radii() {
        let a = Vector3::new(0.0, 0.0, 0.0);
        let b = Vector3::new(1.0, 0.0, 0.0);
        let d = surface_distance(a, b, 1.9, 1.9);
        assert!(d < 0.0);
        assert!((d - (1.0 - 3.8)).abs() < 1e-12);
    }

    #[test]
    fn cutoff_excludes_far_pairs() {
        let a = Vector3::new(0.0, 0.0, 0.0);
        let close = Vector3::new(7.9, 0.0, 0.0);
        let far = Vector3::new(8.1, 0.0, 0.0);
        assert!(within_cutoff(a, close));
        assert!(!within_cutoff(a, far));
    }

    #[test]
    fn unused_atom_type_import_silenced() {
        // Make sure the Ad4AtomType import path stays exercised even
        // before subsequent tasks reference it directly.
        let _ = Ad4AtomType::C;
    }

    #[test]
    fn gauss1_peaks_at_zero_surface_distance() {
        // exp(-(d-0)^2 / 0.5^2). At d=0 the gaussian is 1.0; at
        // d = 0.5 (one sigma) it drops to 1/e.
        assert!((gauss1(0.0) - 1.0).abs() < 1e-12);
        let expected_at_sigma = (-1.0_f64).exp();
        assert!((gauss1(0.5) - expected_at_sigma).abs() < 1e-12);
    }

    #[test]
    fn gauss2_peaks_at_three_angstrom() {
        assert!((gauss2(3.0) - 1.0).abs() < 1e-12);
        // At d = 3.0 + 2.0 sigma the gaussian is 1/e.
        assert!((gauss2(5.0) - (-1.0_f64).exp()).abs() < 1e-12);
    }

    #[test]
    fn repulsion_zero_at_contact_quadratic_inside() {
        assert_eq!(repulsion(0.0), 0.0);
        assert_eq!(repulsion(1.0), 0.0); // outside the wall, no penalty
        assert_eq!(repulsion(-0.5), 0.25);
        assert_eq!(repulsion(-2.0), 4.0);
    }

    #[test]
    fn hydrophobic_ramp_endpoints_and_midpoint() {
        assert_eq!(hydrophobic_factor(0.0), 1.0);
        assert_eq!(hydrophobic_factor(0.5), 1.0);
        assert_eq!(hydrophobic_factor(1.5), 0.0);
        assert_eq!(hydrophobic_factor(2.0), 0.0);
        // d = 1.0 is the midpoint of [0.5, 1.5] — linear ramp at 0.5.
        assert!((hydrophobic_factor(1.0) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn hydrophobic_pair_filter() {
        // Both hydrophobic — full factor.
        assert_eq!(hydrophobic_pair(Ad4AtomType::C, Ad4AtomType::C, 0.0), 1.0);
        // Mixed — zero.
        assert_eq!(hydrophobic_pair(Ad4AtomType::C, Ad4AtomType::OA, 0.0), 0.0);
    }

    #[test]
    fn hbond_ramp_endpoints_and_midpoint() {
        assert_eq!(hbond_factor(-1.0), 1.0);
        assert_eq!(hbond_factor(-0.7), 1.0);
        assert_eq!(hbond_factor(0.0), 0.0);
        assert_eq!(hbond_factor(0.5), 0.0);
        // d = -0.35 is the midpoint of [-0.7, 0.0] — linear ramp at 0.5.
        assert!((hbond_factor(-0.35) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn hbond_pair_requires_donor_acceptor() {
        // HD donor + OA acceptor → full factor at d = -1.
        assert!((hbond_pair(Ad4AtomType::HD, Ad4AtomType::OA, -1.0) - 1.0).abs() < 1e-12);
        // OA + OA: two acceptors, no donor → zero.
        assert_eq!(hbond_pair(Ad4AtomType::OA, Ad4AtomType::OA, -1.0), 0.0);
        // C + C: neither donor nor acceptor → zero.
        assert_eq!(hbond_pair(Ad4AtomType::C, Ad4AtomType::C, -1.0), 0.0);
    }

    #[test]
    fn pair_score_two_carbons_at_contact() {
        // Two aliphatic carbons (C/C) at surface distance 0:
        //   gauss1(0) = 1, gauss2(0) = exp(-9/4) ≈ 0.105399
        //   repulsion(0) = 0
        //   hydrophobic(0) = 1
        //   hbond = 0 (no donor/acceptor)
        // Score = w_g1*1 + w_g2*0.105399 + w_rep*0 + w_hyd*1 + w_hb*0
        let a = Ad4AtomType::C;
        let b = Ad4AtomType::C;
        let expected = weights::GAUSS1 * 1.0
            + weights::GAUSS2 * (-9.0_f64 / 4.0).exp()
            + weights::HYDROPHOBIC * 1.0;
        let got = pair_score(a, b, 0.0);
        assert!(
            (got - expected).abs() < 1e-9,
            "got {got}, expected {expected}"
        );
    }

    #[test]
    fn pair_score_zero_far_apart() {
        // At d = 100 Å, both gaussians are ~0, ramps are 0, repulsion 0.
        // Total should be effectively zero (just gaussian numerical leak).
        let s = pair_score(Ad4AtomType::C, Ad4AtomType::C, 100.0);
        assert!(s.abs() < 1e-30);
    }
}
