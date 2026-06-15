//! End-to-end rosette analysis: gauge readings to strain, principal
//! strains, and (optionally) plane stress in one struct.

use serde::{Deserialize, Serialize};

use crate::material::{ElasticMaterial, PlaneStress};
use crate::rosette::{
    principal_strains, reduce, CartesianStrain, PrincipalStrain, RosetteReadings,
};

/// Complete reduction of a rectangular rosette together with the stress
/// state implied by an [`ElasticMaterial`].
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RosetteAnalysis {
    /// Recovered Cartesian strain state.
    pub strain: CartesianStrain,
    /// Principal strains and the major-axis orientation.
    pub principal: PrincipalStrain,
    /// In-plane stress from 2D Hooke's law.
    pub stress: PlaneStress,
}

/// Run the full pipeline: reduce the gauges, solve the principal
/// strains, and map to plane stress with `material`.
///
/// This is the convenience entry point most callers want. Each stage is
/// also available standalone ([`reduce`], [`principal_strains`],
/// [`ElasticMaterial::plane_stress`]).
///
/// # Examples
///
/// ```
/// use valenx_strainrosette::{analyze, ElasticMaterial, RosetteReadings};
///
/// let mat = ElasticMaterial::new(200_000.0, 0.3).unwrap();
/// let a = analyze(RosetteReadings::new(0.0010, 0.0006, 0.0002), &mat);
/// // eps_x is the 0-degree gauge.
/// assert!((a.strain.eps_x - 0.0010).abs() < 1e-12);
/// // Major principal strain is at least the minor one.
/// assert!(a.principal.eps_1 >= a.principal.eps_2);
/// ```
pub fn analyze(readings: RosetteReadings, material: &ElasticMaterial) -> RosetteAnalysis {
    let strain = reduce(readings);
    let principal = principal_strains(strain);
    let stress = material.plane_stress(strain);
    RosetteAnalysis {
        strain,
        principal,
        stress,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn pipeline_stages_agree_with_standalone_calls() {
        let mat = ElasticMaterial::new(210_000.0, 0.3).unwrap();
        let r = RosetteReadings::new(9.0e-4, 4.0e-4, 1.0e-4);
        let a = analyze(r, &mat);

        let s = reduce(r);
        assert_eq!(a.strain, s);
        assert_eq!(a.principal, principal_strains(s));
        assert_eq!(a.stress, mat.plane_stress(s));
    }

    #[test]
    fn uniaxial_end_to_end_recovers_applied_stress() {
        // Drive a known uniaxial-stress field. Apply eps_x, with the
        // physically consistent lateral strain eps_y = -nu * eps_x and a
        // 45-degree gauge reading of (eps_x + eps_y)/2 (gamma = 0). The
        // pipeline must recover sigma_x = E * eps_x, sigma_y = 0.
        let e = 200_000.0;
        let nu = 0.3;
        let mat = ElasticMaterial::new(e, nu).unwrap();
        let eps_x = 1.0e-3;
        let eps_y = -nu * eps_x;
        let eps_45 = 0.5 * (eps_x + eps_y); // gamma_xy = 0 for aligned principal axes
        let a = analyze(RosetteReadings::new(eps_x, eps_45, eps_y), &mat);

        assert!((a.strain.eps_x - eps_x).abs() < 1e-12);
        assert!((a.strain.eps_y - eps_y).abs() < 1e-12);
        assert!(
            a.strain.gamma_xy.abs() < 1e-12,
            "gamma = {g}",
            g = a.strain.gamma_xy
        );

        // Principal strains are the applied axial and lateral strains.
        assert!((a.principal.eps_1 - eps_x).abs() < 1e-12);
        assert!((a.principal.eps_2 - eps_y).abs() < 1e-12);
        assert!(a.principal.theta_p.abs() < 1e-12);

        // Stress: uniaxial.
        assert!(
            (a.stress.sigma_x - e * eps_x).abs() < 1e-6,
            "sigma_x = {sx}",
            sx = a.stress.sigma_x
        );
        assert!(
            a.stress.sigma_y.abs() < EPS,
            "sigma_y = {sy}",
            sy = a.stress.sigma_y
        );
        assert!(
            a.stress.tau_xy.abs() < EPS,
            "tau_xy = {t}",
            t = a.stress.tau_xy
        );
    }
}
