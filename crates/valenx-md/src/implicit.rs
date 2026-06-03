//! Implicit solvent — generalized Born (**roadmap feature 30b**).
//!
//! Explicit water is expensive — most atoms in a solvated simulation
//! are solvent. An **implicit solvent** replaces the water with a
//! continuum dielectric and adds an analytic *solvation free energy*
//! term to the force field. The **generalized Born** (GB) model is the
//! standard fast implicit electrostatics.
//!
//! GB works in two stages:
//!
//! 1. **Born radii.** Each atom `i` gets an effective Born radius
//!    `Rᵢ` — loosely, how deeply buried it is. A surface atom has a
//!    radius close to its intrinsic (van der Waals) radius; a buried
//!    atom has a much larger one. This v1 uses the **Hawkins-Cramer-
//!    Truhlar (HCT)** pairwise-descreening approximation: every other
//!    atom "descreens" atom `i`, reducing its inverse Born radius.
//!
//! 2. **GB energy.** The electrostatic solvation free energy is the
//!    **Still** GB pair kernel:
//!
//!    ```text
//!    ΔG_pol = −½·f·(1/ε_in − 1/ε_out)·Σᵢⱼ qᵢqⱼ / f_GB(rᵢⱼ, Rᵢ, Rⱼ)
//!    f_GB = √( rᵢⱼ² + RᵢRⱼ·exp(−rᵢⱼ²/(4RᵢRⱼ)) )
//!    ```
//!
//!    `ε_in` is the solute dielectric (≈ 1), `ε_out` the solvent
//!    (≈ 78.5 for water). The `i = j` self-term is the Born
//!    self-energy of each charge.
//!
//! ## v1 caveat — honest scope
//!
//! This is a real GB/HCT implementation: HCT pairwise descreening for
//! the Born radii and the Still kernel for the energy. What it does
//! *not* include — and a production GBSA model would — is the
//! **nonpolar / surface-area** (SASA-scaled) cavitation term, and the
//! later **GBn / GBn2** neck corrections that improve buried-atom
//! radii. The polar electrostatic solvation energy is the real v1; the
//! nonpolar term is the documented next step.

use nalgebra::Vector3;

use crate::error::{MdError, Result};
use crate::units::COULOMB;

/// Dielectric constant of liquid water at room temperature — the usual
/// solvent dielectric for a GB model.
pub const WATER_DIELECTRIC: f64 = 78.5;

/// Parameters of the generalized-Born model.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct GbParams {
    /// Solute (interior) dielectric constant ε_in.
    pub solute_dielectric: f64,
    /// Solvent (exterior) dielectric constant ε_out.
    pub solvent_dielectric: f64,
    /// HCT overlap scale factor (the descreening strength). The
    /// canonical HCT value is ~0.8.
    pub hct_scale: f64,
    /// A small offset (nm) subtracted from intrinsic radii when
    /// forming the descreening integral — the "dielectric offset" of
    /// the GB-HCT model.
    pub dielectric_offset: f64,
}

impl Default for GbParams {
    fn default() -> Self {
        GbParams {
            solute_dielectric: 1.0,
            solvent_dielectric: WATER_DIELECTRIC,
            hct_scale: 0.8,
            dielectric_offset: 0.009,
        }
    }
}

impl GbParams {
    /// Validates the parameters.
    fn check(&self) -> Result<()> {
        if !(self.solute_dielectric.is_finite() && self.solute_dielectric >= 1.0) {
            return Err(MdError::invalid(
                "solute_dielectric",
                "must be finite and at least 1",
            ));
        }
        if !(self.solvent_dielectric.is_finite() && self.solvent_dielectric > 1.0) {
            return Err(MdError::invalid(
                "solvent_dielectric",
                "must be finite and greater than 1",
            ));
        }
        if !(self.hct_scale.is_finite() && self.hct_scale > 0.0) {
            return Err(MdError::invalid("hct_scale", "must be finite and positive"));
        }
        Ok(())
    }
}

/// One atom as seen by the GB model: a position, a charge, and an
/// intrinsic (van der Waals) radius.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct GbAtom {
    /// Position (nm).
    pub position: Vector3<f64>,
    /// Partial charge (e).
    pub charge: f64,
    /// Intrinsic / van der Waals radius (nm).
    pub radius: f64,
}

impl GbAtom {
    /// Builds a GB atom.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if the radius is not finite and positive
    /// or the charge is not finite.
    pub fn new(position: Vector3<f64>, charge: f64, radius: f64) -> Result<Self> {
        if !(radius.is_finite() && radius > 0.0) {
            return Err(MdError::invalid("radius", "must be finite and positive"));
        }
        if !charge.is_finite() {
            return Err(MdError::invalid("charge", "must be finite"));
        }
        Ok(GbAtom {
            position,
            charge,
            radius,
        })
    }
}

/// Computes the effective Born radii of a set of atoms with the
/// Hawkins-Cramer-Truhlar pairwise-descreening approximation.
///
/// Returns one radius per atom (nm). A more buried atom gets a larger
/// radius.
///
/// # Errors
/// [`MdError::Invalid`] for bad parameters.
pub fn born_radii(atoms: &[GbAtom], params: &GbParams) -> Result<Vec<f64>> {
    params.check()?;
    let n = atoms.len();
    let mut radii = Vec::with_capacity(n);
    for i in 0..n {
        // Intrinsic radius minus the dielectric offset.
        let rho_i = (atoms[i].radius - params.dielectric_offset).max(1e-3);
        // Start the inverse Born radius at the isolated-atom value.
        let mut inv_born = 1.0 / rho_i;
        // Subtract every other atom's HCT descreening contribution.
        for j in 0..n {
            if i == j {
                continue;
            }
            let r = (atoms[i].position - atoms[j].position).norm();
            if r < 1e-9 {
                continue;
            }
            let rho_j =
                ((atoms[j].radius - params.dielectric_offset) * params.hct_scale).max(1e-3);
            inv_born -= hct_descreen(r, rho_i, rho_j);
        }
        // Clamp: the Born radius cannot be smaller than the intrinsic
        // radius and stays finite.
        let born = (1.0 / inv_born.max(1.0 / 50.0)).max(rho_i);
        radii.push(born);
    }
    Ok(radii)
}

/// One atom's HCT descreening contribution to another's inverse Born
/// radius. This is the standard HCT integral piece (returned as a
/// non-negative amount to subtract).
fn hct_descreen(r: f64, rho_i: f64, rho_j: f64) -> f64 {
    // The descreener only matters when its sphere reaches atom i.
    if r >= rho_i + rho_j {
        // Fully separated: the closed-form HCT term for the
        // non-overlapping case.
        let l = r - rho_j;
        let u = r + rho_j;
        if l <= rho_i {
            // Partial: clamp the lower limit to rho_i.
            return descreen_integral(r, rho_i, rho_j, rho_i.max(l), u);
        }
        return descreen_integral(r, rho_i, rho_j, l, u);
    }
    // Overlapping spheres: integrate from rho_i outward.
    let u = r + rho_j;
    descreen_integral(r, rho_i, rho_j, rho_i, u.max(rho_i))
}

/// The HCT descreening integral evaluated between limits `l` and `u`.
///
/// This is the standard Hawkins-Cramer-Truhlar closed form for the
/// volume integral of a descreening sphere of radius `rho_j` whose
/// centre is a distance `r` away.
fn descreen_integral(r: f64, _rho_i: f64, rho_j: f64, l: f64, u: f64) -> f64 {
    if u <= l {
        return 0.0;
    }
    let value = 1.0 / l - 1.0 / u
        + 0.25 * (r - rho_j * rho_j / r) * (1.0 / (u * u) - 1.0 / (l * l))
        + 0.5 / r * (l / u).ln();
    (0.5 * value).max(0.0)
}

/// Computes the generalized-Born polar solvation free energy (kJ/mol)
/// of a set of atoms.
///
/// Negative — solvation stabilises a charged solute. Uses the Still GB
/// kernel with the [`born_radii`] effective radii.
///
/// # Errors
/// [`MdError::Invalid`] for bad parameters.
pub fn gb_solvation_energy(atoms: &[GbAtom], params: &GbParams) -> Result<f64> {
    params.check()?;
    if atoms.len() < 2 {
        // A single charge still has a self-energy; handle n = 0/1.
        if atoms.len() == 1 {
            let radii = born_radii(atoms, params)?;
            let prefactor =
                -0.5 * COULOMB * (1.0 / params.solute_dielectric - 1.0 / params.solvent_dielectric);
            let q = atoms[0].charge;
            return Ok(prefactor * q * q / radii[0]);
        }
        return Ok(0.0);
    }
    let radii = born_radii(atoms, params)?;
    let prefactor =
        -0.5 * COULOMB * (1.0 / params.solute_dielectric - 1.0 / params.solvent_dielectric);
    let n = atoms.len();
    let mut energy = 0.0;
    for i in 0..n {
        for j in 0..n {
            let qi = atoms[i].charge;
            let qj = atoms[j].charge;
            if qi == 0.0 || qj == 0.0 {
                continue;
            }
            let r2 = if i == j {
                0.0
            } else {
                (atoms[i].position - atoms[j].position).norm_squared()
            };
            let f_gb = still_kernel(r2, radii[i], radii[j]);
            energy += prefactor * qi * qj / f_gb;
        }
    }
    Ok(energy)
}

/// The Still generalized-Born interaction kernel
/// `f_GB = √(r² + Rᵢ·Rⱼ·exp(−r²/(4·Rᵢ·Rⱼ)))`.
fn still_kernel(r2: f64, born_i: f64, born_j: f64) -> f64 {
    let rr = born_i * born_j;
    (r2 + rr * (-r2 / (4.0 * rr)).exp()).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gb_atom_validates_input() {
        assert!(GbAtom::new(Vector3::zeros(), 1.0, 0.0).is_err());
        assert!(GbAtom::new(Vector3::zeros(), f64::NAN, 0.15).is_err());
        assert!(GbAtom::new(Vector3::zeros(), 1.0, 0.15).is_ok());
    }

    #[test]
    fn params_validate() {
        let mut p = GbParams::default();
        assert!(p.check().is_ok());
        p.solvent_dielectric = 0.5;
        assert!(p.check().is_err());
        let p2 = GbParams {
            hct_scale: -1.0,
            ..GbParams::default()
        };
        assert!(p2.check().is_err());
    }

    #[test]
    fn buried_atom_has_larger_born_radius() {
        // A central atom surrounded by neighbours vs an isolated atom.
        let r = 0.17;
        let mut crowded = vec![GbAtom::new(Vector3::zeros(), 0.0, r).unwrap()];
        for (dx, dy, dz) in [
            (0.35, 0.0, 0.0),
            (-0.35, 0.0, 0.0),
            (0.0, 0.35, 0.0),
            (0.0, -0.35, 0.0),
            (0.0, 0.0, 0.35),
            (0.0, 0.0, -0.35),
        ] {
            crowded.push(
                GbAtom::new(Vector3::new(dx, dy, dz), 0.0, r).unwrap(),
            );
        }
        let radii_crowded = born_radii(&crowded, &GbParams::default()).unwrap();

        let isolated = vec![GbAtom::new(Vector3::zeros(), 0.0, r).unwrap()];
        let radii_isolated = born_radii(&isolated, &GbParams::default()).unwrap();

        // The buried central atom should have a larger effective Born
        // radius than the same atom in isolation.
        assert!(
            radii_crowded[0] > radii_isolated[0],
            "buried {} not > isolated {}",
            radii_crowded[0],
            radii_isolated[0]
        );
    }

    #[test]
    fn born_radius_is_at_least_intrinsic() {
        let atoms = vec![
            GbAtom::new(Vector3::zeros(), 1.0, 0.15).unwrap(),
            GbAtom::new(Vector3::new(0.4, 0.0, 0.0), -1.0, 0.15).unwrap(),
        ];
        let radii = born_radii(&atoms, &GbParams::default()).unwrap();
        for r in radii {
            assert!(r >= 0.15 - 0.01, "Born radius {r} below intrinsic");
            assert!(r.is_finite());
        }
    }

    #[test]
    fn solvation_energy_is_negative_for_charges() {
        // A pair of opposite charges in water: GB solvation stabilises.
        let atoms = vec![
            GbAtom::new(Vector3::zeros(), 1.0, 0.17).unwrap(),
            GbAtom::new(Vector3::new(0.5, 0.0, 0.0), -1.0, 0.17).unwrap(),
        ];
        let e = gb_solvation_energy(&atoms, &GbParams::default()).unwrap();
        assert!(e < 0.0, "GB solvation energy = {e}");
        assert!(e.is_finite());
    }

    #[test]
    fn single_charge_has_a_self_energy() {
        let atoms = vec![GbAtom::new(Vector3::zeros(), 1.0, 0.17).unwrap()];
        let e = gb_solvation_energy(&atoms, &GbParams::default()).unwrap();
        // The Born self-energy of a charge is negative (favourable).
        assert!(e < 0.0);
    }

    #[test]
    fn neutral_system_has_zero_solvation_energy() {
        let atoms = vec![
            GbAtom::new(Vector3::zeros(), 0.0, 0.17).unwrap(),
            GbAtom::new(Vector3::new(0.5, 0.0, 0.0), 0.0, 0.17).unwrap(),
        ];
        let e = gb_solvation_energy(&atoms, &GbParams::default()).unwrap();
        assert!(e.abs() < 1e-12);
    }

    #[test]
    fn still_kernel_endpoints() {
        // At r = 0 the kernel reduces to √(RᵢRⱼ).
        let k0 = still_kernel(0.0, 0.2, 0.3);
        assert!((k0 - (0.2_f64 * 0.3).sqrt()).abs() < 1e-9);
        // At large r the kernel approaches r.
        let r = 5.0;
        let klarge = still_kernel(r * r, 0.2, 0.3);
        assert!((klarge - r).abs() < 1e-3);
    }
}
