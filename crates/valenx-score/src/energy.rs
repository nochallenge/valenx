//! Validated molecular-energy primitives (kcal/mol; charges in `e`, distances
//! in Å).

use crate::error::{require_finite, require_positive, ScoreError};

/// Coulomb constant in kcal·Å·mol⁻¹·e⁻² (so `coulomb` is in kcal/mol for
/// charges in elementary charge and distances in Å).
pub const K_COULOMB: f64 = 332.0637;

/// Coulomb electrostatic energy `k·q1·q2 / (eps·r)`.
///
/// `r` must be `> 0` and the dielectric `eps >= 1`. Like charges give a
/// positive (repulsive) energy; opposite charges a negative one.
pub fn coulomb(q1: f64, q2: f64, r: f64, dielectric: f64) -> Result<f64, ScoreError> {
    require_finite("charge", q1)?;
    require_finite("charge", q2)?;
    require_positive("distance", r)?;
    if !dielectric.is_finite() || dielectric < 1.0 {
        return Err(ScoreError::DielectricTooSmall { value: dielectric });
    }
    Ok(K_COULOMB * q1 * q2 / (dielectric * r))
}

/// The 12-6 Lennard-Jones potential `4·eps·((sigma/r)¹² − (sigma/r)⁶)`.
///
/// `eps`, `sigma`, `r` must be `> 0`. Zero at `r = sigma`, minimum `−eps` at
/// `r = 2^(1/6)·sigma`, strongly repulsive as `r → 0`.
pub fn lennard_jones(eps: f64, sigma: f64, r: f64) -> Result<f64, ScoreError> {
    require_positive("epsilon", eps)?;
    require_positive("sigma", sigma)?;
    require_positive("distance", r)?;
    let sr6 = (sigma / r).powi(6);
    Ok(4.0 * eps * (sr6 * sr6 - sr6))
}

/// The distance `2^(1/6)·sigma` at which [`lennard_jones`] reaches its minimum.
pub fn lj_min_distance(sigma: f64) -> f64 {
    2.0_f64.powf(1.0 / 6.0) * sigma
}

/// The Born self-solvation free energy of a charge `q` in a sphere of
/// `radius`: `−(k/2)·(1/eps_in − 1/eps_out)·q² / radius`.
///
/// Negative (favorable) when moving a charge from a low-dielectric interior
/// (`eps_in`, e.g. 1) into a high-dielectric solvent (`eps_out`, e.g. 78.5);
/// zero when the two dielectrics are equal.
pub fn born_solvation(q: f64, radius: f64, eps_in: f64, eps_out: f64) -> Result<f64, ScoreError> {
    require_finite("charge", q)?;
    require_positive("born radius", radius)?;
    check_dielectric(eps_in)?;
    check_dielectric(eps_out)?;
    Ok(-0.5 * K_COULOMB * (1.0 / eps_in - 1.0 / eps_out) * q * q / radius)
}

/// The Still generalized-Born effective distance
/// `f_GB = sqrt(r² + a1·a2·exp(−r²/(4·a1·a2)))`.
///
/// Reduces to the Born radius `a` when `r = 0` and `a1 = a2 = a`, and to `r`
/// for large `r`. `a1`, `a2` must be `> 0`; `r >= 0`.
pub fn gb_function(r: f64, a1: f64, a2: f64) -> Result<f64, ScoreError> {
    require_positive("born radius", a1)?;
    require_positive("born radius", a2)?;
    require_finite("distance", r)?;
    if r < 0.0 {
        return Err(ScoreError::NonPositive {
            what: "distance",
            value: r,
        });
    }
    let prod = a1 * a2;
    Ok((r * r + prod * (-r * r / (4.0 * prod)).exp()).sqrt())
}

/// The generalized-Born pairwise polar solvation energy
/// `−(k/2)·(1/eps_in − 1/eps_out)·q1·q2 / f_GB`.
///
/// With `q1 = q2`, `r = 0`, `a1 = a2` this equals [`born_solvation`].
pub fn gb_pair_polar(
    q1: f64,
    q2: f64,
    r: f64,
    a1: f64,
    a2: f64,
    eps_in: f64,
    eps_out: f64,
) -> Result<f64, ScoreError> {
    require_finite("charge", q1)?;
    require_finite("charge", q2)?;
    check_dielectric(eps_in)?;
    check_dielectric(eps_out)?;
    let f = gb_function(r, a1, a2)?;
    Ok(-0.5 * K_COULOMB * (1.0 / eps_in - 1.0 / eps_out) * q1 * q2 / f)
}

/// A SASA-linear nonpolar solvation energy `gamma·sasa + beta` (`sasa` in Å²,
/// `gamma` in kcal·mol⁻¹·Å⁻²).
pub fn sasa_nonpolar(sasa: f64, gamma: f64, beta: f64) -> Result<f64, ScoreError> {
    if !sasa.is_finite() || sasa < 0.0 {
        return Err(ScoreError::NonPositive {
            what: "sasa",
            value: sasa,
        });
    }
    require_finite("gamma", gamma)?;
    require_finite("beta", beta)?;
    Ok(gamma * sasa + beta)
}

fn check_dielectric(eps: f64) -> Result<(), ScoreError> {
    if !eps.is_finite() || eps < 1.0 {
        return Err(ScoreError::DielectricTooSmall { value: eps });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coulomb_sign_and_scaling() {
        // opposite charges attract (negative)
        assert!(coulomb(1.0, -1.0, 1.0, 1.0).unwrap() < 0.0);
        // like charges repel (positive), magnitude k at unit separation
        assert!((coulomb(1.0, 1.0, 1.0, 1.0).unwrap() - K_COULOMB).abs() < 1e-9);
        // 1/r scaling: doubling r halves the magnitude
        let e1 = coulomb(1.0, 1.0, 1.0, 1.0).unwrap();
        let e2 = coulomb(1.0, 1.0, 2.0, 1.0).unwrap();
        assert!((e2 - e1 / 2.0).abs() < 1e-9);
        // dielectric screens
        assert!((coulomb(1.0, 1.0, 1.0, 80.0).unwrap() - K_COULOMB / 80.0).abs() < 1e-9);
    }

    #[test]
    fn lj_well_bottom_and_zero_crossing() {
        let (eps, sigma) = (0.5, 3.4);
        assert!((lennard_jones(eps, sigma, sigma).unwrap()).abs() < 1e-9); // zero at sigma
        let rmin = lj_min_distance(sigma);
        assert!((lennard_jones(eps, sigma, rmin).unwrap() + eps).abs() < 1e-9); // -eps at min
                                                                                // repulsive (large positive) well inside sigma
        assert!(lennard_jones(eps, sigma, 0.5 * sigma).unwrap() > 0.0);
    }

    #[test]
    fn born_is_favorable_and_zero_when_dielectrics_match() {
        assert!(born_solvation(1.0, 2.0, 1.0, 78.5).unwrap() < 0.0);
        assert!(born_solvation(1.0, 2.0, 78.5, 78.5).unwrap().abs() < 1e-12);
        // hand value: -(332.0637/2)*(1 - 1/78.5)*1 / 2
        let expect = -0.5 * K_COULOMB * (1.0 - 1.0 / 78.5) / 2.0;
        assert!((born_solvation(1.0, 2.0, 1.0, 78.5).unwrap() - expect).abs() < 1e-9);
    }

    #[test]
    fn gb_function_limits() {
        // r = 0, equal radii -> the radius itself
        assert!((gb_function(0.0, 2.5, 2.5).unwrap() - 2.5).abs() < 1e-12);
        // r = 0, unequal -> sqrt(a1*a2)
        assert!((gb_function(0.0, 2.0, 8.0).unwrap() - 4.0).abs() < 1e-12);
        // large r -> approaches r
        let big = gb_function(100.0, 2.0, 2.0).unwrap();
        assert!((big - 100.0).abs() < 1e-3);
    }

    #[test]
    fn gb_pair_reduces_to_born_self() {
        let born = born_solvation(1.5, 2.0, 1.0, 78.5).unwrap();
        let gb = gb_pair_polar(1.5, 1.5, 0.0, 2.0, 2.0, 1.0, 78.5).unwrap();
        assert!((born - gb).abs() < 1e-9);
    }

    #[test]
    fn sasa_is_linear() {
        assert!((sasa_nonpolar(100.0, 0.0072, 0.0).unwrap() - 0.72).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(
            coulomb(1.0, 1.0, 0.0, 1.0).unwrap_err().code(),
            "non_positive"
        );
        assert_eq!(
            coulomb(1.0, 1.0, 1.0, 0.5).unwrap_err().code(),
            "dielectric_too_small"
        );
        assert_eq!(
            lennard_jones(-1.0, 1.0, 1.0).unwrap_err().code(),
            "non_positive"
        );
        assert_eq!(
            sasa_nonpolar(-1.0, 0.1, 0.0).unwrap_err().code(),
            "non_positive"
        );
    }
}
