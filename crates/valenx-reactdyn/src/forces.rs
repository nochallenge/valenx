//! Nuclear forces for AIMD — computed **numerically** by central finite
//! differences of `valenx-qchem`'s single-point energy.
//!
//! qchem ships no analytic nuclear gradient (`GeometryOptRequest::run` is
//! a documented stub), so the force on every atomic coordinate is
//! `F_iα = −[E(r + δ·e_iα) − E(r − δ·e_iα)] / (2δ)`. This costs `6N`
//! single-point energies per force evaluation, which is why AIMD here is
//! scoped to small systems — but it is a legitimate, well-defined method
//! that uses only the existing, validated single-point core.
//!
//! Positions are in **bohr** and energies in **hartree** (qchem-native),
//! so the returned forces are in hartree/bohr with no conversion.

use valenx_qchem::dft::{Functional, GridQuality};
use valenx_qchem::driver::{run_dft, run_rhf, run_rhf_embedded, run_uhf};
use valenx_qchem::element::Element;
use valenx_qchem::geometry::{Atom, MolecularGeometry};
use valenx_qchem::scf::rhf::ScfSettings;

use crate::error::{ReactDynError, Result};

/// The electronic-structure method used for each single-point energy.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum Method {
    /// Restricted Hartree-Fock (closed-shell).
    #[default]
    Rhf,
    /// Unrestricted Hartree-Fock (open-shell).
    Uhf,
    /// Restricted Kohn-Sham DFT with the B3LYP functional.
    Dft,
}

/// Single-point electronic energy (hartree) at a geometry given in bohr.
///
/// # Errors
/// Returns [`ReactDynError::Qchem`] if the SCF fails, a basis is missing,
/// or the geometry is rejected by qchem.
pub fn single_point_energy(
    elements: &[Element],
    pos_bohr: &[[f64; 3]],
    charge: i32,
    multiplicity: u32,
    method: Method,
    basis: &str,
) -> Result<f64> {
    let atoms: Vec<Atom> = elements
        .iter()
        .zip(pos_bohr)
        .map(|(e, p)| Atom::new(*e, *p))
        .collect();
    let geom = MolecularGeometry::with_charge_multiplicity(atoms, charge, multiplicity);
    let settings = ScfSettings::default();
    let report = match method {
        Method::Rhf => run_rhf(&geom, basis, settings),
        Method::Uhf => run_uhf(&geom, basis, settings),
        Method::Dft => run_dft(
            &geom,
            basis,
            Functional::B3lyp,
            GridQuality::default(),
            settings,
        ),
    }
    .map_err(|e| ReactDynError::Qchem(e.to_string()))?;
    Ok(report.total_energy)
}

/// Single-point **RHF** energy of the QM region in the field of external
/// point charges `(q, position_bohr)` — electrostatic QM/MM embedding.
/// The MM charges enter the SCF and polarize the density; the returned
/// energy includes the electron–charge interaction (the nuclei–charge
/// term is the caller's). v1 is RHF (closed-shell) only.
pub fn single_point_energy_embedded(
    elements: &[Element],
    pos_bohr: &[[f64; 3]],
    charge: i32,
    multiplicity: u32,
    basis: &str,
    external_charges: &[(f64, [f64; 3])],
) -> Result<f64> {
    let atoms: Vec<Atom> = elements
        .iter()
        .zip(pos_bohr)
        .map(|(e, p)| Atom::new(*e, *p))
        .collect();
    let geom = MolecularGeometry::with_charge_multiplicity(atoms, charge, multiplicity);
    let report = run_rhf_embedded(&geom, basis, ScfSettings::default(), external_charges)
        .map_err(|e| ReactDynError::Qchem(e.to_string()))?;
    Ok(report.total_energy)
}

/// Numerical nuclear forces (hartree/bohr) by central finite difference
/// of the single-point energy. `delta_bohr` is the displacement; ≈ 0.01
/// bohr is a good default (small enough for accuracy, large enough to
/// stay above SCF noise).
///
/// # Errors
/// Propagates any [`ReactDynError::Qchem`] from a perturbed single-point
/// energy.
pub fn numerical_forces(
    elements: &[Element],
    pos_bohr: &[[f64; 3]],
    charge: i32,
    multiplicity: u32,
    method: Method,
    basis: &str,
    delta_bohr: f64,
) -> Result<Vec<[f64; 3]>> {
    if delta_bohr <= 0.0 || !delta_bohr.is_finite() {
        return Err(ReactDynError::Invalid {
            reason: format!("finite-difference delta must be positive (got {delta_bohr})"),
        });
    }
    let n = elements.len();
    let mut forces = vec![[0.0_f64; 3]; n];
    let mut p = pos_bohr.to_vec();
    for i in 0..n {
        for d in 0..3 {
            let orig = p[i][d];
            p[i][d] = orig + delta_bohr;
            let e_plus = single_point_energy(elements, &p, charge, multiplicity, method, basis)?;
            p[i][d] = orig - delta_bohr;
            let e_minus = single_point_energy(elements, &p, charge, multiplicity, method, basis)?;
            p[i][d] = orig;
            forces[i][d] = -(e_plus - e_minus) / (2.0 * delta_bohr);
        }
    }
    Ok(forces)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn h2_force_sign_is_physical() {
        let h = Element::from_symbol("H").unwrap();
        let elems = [h, h];
        let delta = 0.01; // bohr
        let (charge, mult) = (0, 1);

        // Compressed bond (0.5 bohr; equilibrium ≈ 1.4 bohr) → repulsion
        // pushes the atoms apart along z.
        let compressed = [[0.0, 0.0, 0.25], [0.0, 0.0, -0.25]];
        let f = numerical_forces(
            &elems,
            &compressed,
            charge,
            mult,
            Method::Rhf,
            "STO-3G",
            delta,
        )
        .unwrap();
        assert!(
            f[0][2] > 0.0,
            "compressed: top atom should be pushed +z, got {}",
            f[0][2]
        );
        assert!(
            f[1][2] < 0.0,
            "compressed: bottom atom should be pushed -z, got {}",
            f[1][2]
        );

        // Stretched bond (2.0 bohr) → attraction pulls the atoms together.
        let stretched = [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];
        let f = numerical_forces(
            &elems,
            &stretched,
            charge,
            mult,
            Method::Rhf,
            "STO-3G",
            delta,
        )
        .unwrap();
        assert!(
            f[0][2] < 0.0,
            "stretched: top atom should be pulled -z, got {}",
            f[0][2]
        );
        assert!(
            f[1][2] > 0.0,
            "stretched: bottom atom should be pulled +z, got {}",
            f[1][2]
        );
    }

    #[test]
    fn non_positive_delta_fails_loud() {
        let h = Element::from_symbol("H").unwrap();
        let r = numerical_forces(
            &[h, h],
            &[[0.0; 3], [0.0, 0.0, 1.4]],
            0,
            1,
            Method::Rhf,
            "STO-3G",
            0.0,
        );
        assert!(matches!(r, Err(ReactDynError::Invalid { .. })));
    }
}
