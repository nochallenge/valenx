//! A transparent MM-GBSA-style endpoint binding-energy decomposition.

use serde::{Deserialize, Serialize};

use crate::energy::{coulomb, gb_pair_polar, lennard_jones};
use crate::error::{require_finite, require_positive, ScoreError};

/// A single interaction site ("bead"): a point charge with Lennard-Jones and
/// Born parameters. A deliberately coarse stand-in for an atom or a group.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Bead {
    /// Partial charge in elementary charge `e`.
    pub charge: f64,
    /// Lennard-Jones well depth (kcal/mol), `> 0`.
    pub lj_epsilon: f64,
    /// Lennard-Jones diameter (Å), `> 0`.
    pub lj_sigma: f64,
    /// Effective Born radius (Å), `> 0`.
    pub born_radius: f64,
}

impl Bead {
    /// A validated bead.
    pub fn new(
        charge: f64,
        lj_epsilon: f64,
        lj_sigma: f64,
        born_radius: f64,
    ) -> Result<Self, ScoreError> {
        require_finite("charge", charge)?;
        require_positive("lj_epsilon", lj_epsilon)?;
        require_positive("lj_sigma", lj_sigma)?;
        require_positive("born_radius", born_radius)?;
        Ok(Self {
            charge,
            lj_epsilon,
            lj_sigma,
            born_radius,
        })
    }
}

/// An endpoint binding-energy decomposition (kcal/mol). The headline number is
/// [`BindingEnergy::total`]; the components are kept visible so a reviewer can
/// see *why* a candidate ranks where it does.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BindingEnergy {
    /// Van der Waals (Lennard-Jones) interaction.
    pub vdw: f64,
    /// Gas-phase electrostatic (Coulomb) interaction.
    pub electrostatic: f64,
    /// Polar solvation / screening (generalized Born).
    pub polar_solvation: f64,
    /// Nonpolar solvation from buried surface area.
    pub nonpolar_solvation: f64,
}

impl BindingEnergy {
    /// The total binding energy: the sum of the four components.
    pub fn total(&self) -> f64 {
        self.vdw + self.electrostatic + self.polar_solvation + self.nonpolar_solvation
    }

    /// Assemble a decomposition for a minimal two-bead "ligand + receptor"
    /// contact at center-to-center `separation` (Å) in a solvent of dielectric
    /// `solvent_dielectric`, burying `buried_sasa` Å² of surface
    /// (`gamma` = nonpolar surface tension, kcal·mol⁻¹·Å⁻²).
    ///
    /// Lorentz-Berthelot combining rules are used for the LJ parameters
    /// (`sigma = (σ_l+σ_r)/2`, `epsilon = sqrt(ε_l·ε_r)`). This is the
    /// **illustrative endpoint** described in the crate docs — a single-contact
    /// estimate, not a validated ΔG.
    pub fn two_bead_estimate(
        ligand: &Bead,
        receptor: &Bead,
        separation: f64,
        solvent_dielectric: f64,
        buried_sasa: f64,
        gamma: f64,
    ) -> Result<Self, ScoreError> {
        require_positive("separation", separation)?;
        require_finite("gamma", gamma)?;
        if !buried_sasa.is_finite() || buried_sasa < 0.0 {
            return Err(ScoreError::NonPositive {
                what: "buried_sasa",
                value: buried_sasa,
            });
        }
        let sigma_ij = 0.5 * (ligand.lj_sigma + receptor.lj_sigma);
        let eps_ij = (ligand.lj_epsilon * receptor.lj_epsilon).sqrt();
        let vdw = lennard_jones(eps_ij, sigma_ij, separation)?;
        let electrostatic = coulomb(ligand.charge, receptor.charge, separation, 1.0)?;
        let polar_solvation = gb_pair_polar(
            ligand.charge,
            receptor.charge,
            separation,
            ligand.born_radius,
            receptor.born_radius,
            1.0,
            solvent_dielectric,
        )?;
        // Burial removes exposed surface, a favorable (negative) nonpolar term.
        let nonpolar_solvation = -gamma * buried_sasa;
        Ok(Self {
            vdw,
            electrostatic,
            polar_solvation,
            nonpolar_solvation,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::energy::lj_min_distance;

    fn neutral_bead() -> Bead {
        Bead::new(0.0, 0.5, 3.4, 2.0).unwrap()
    }

    #[test]
    fn total_is_component_sum() {
        let be = BindingEnergy {
            vdw: -1.0,
            electrostatic: -2.0,
            polar_solvation: 0.5,
            nonpolar_solvation: -0.3,
        };
        assert!((be.total() - (-2.8)).abs() < 1e-12);
    }

    #[test]
    fn neutral_beads_at_lj_min_give_minus_epsilon() {
        let b = neutral_bead();
        let sigma_ij = b.lj_sigma; // identical beads
        let eps_ij = b.lj_epsilon;
        let r = lj_min_distance(sigma_ij);
        let be = BindingEnergy::two_bead_estimate(&b, &b, r, 78.5, 0.0, 0.0072).unwrap();
        assert!((be.vdw + eps_ij).abs() < 1e-9); // -epsilon
        assert!(be.electrostatic.abs() < 1e-12); // neutral
        assert!(be.polar_solvation.abs() < 1e-12); // neutral
        assert!(be.nonpolar_solvation.abs() < 1e-12); // no burial
        assert!((be.total() + eps_ij).abs() < 1e-9);
    }

    #[test]
    fn opposite_charges_give_negative_electrostatic() {
        let pos = Bead::new(1.0, 0.5, 3.4, 2.0).unwrap();
        let neg = Bead::new(-1.0, 0.5, 3.4, 2.0).unwrap();
        let be = BindingEnergy::two_bead_estimate(&pos, &neg, 4.0, 78.5, 0.0, 0.0072).unwrap();
        assert!(be.electrostatic < 0.0);
    }

    #[test]
    fn burial_gives_favorable_nonpolar() {
        let b = neutral_bead();
        let be = BindingEnergy::two_bead_estimate(&b, &b, 4.0, 78.5, 50.0, 0.0072).unwrap();
        assert!((be.nonpolar_solvation + 0.0072 * 50.0).abs() < 1e-12);
        assert!(be.nonpolar_solvation < 0.0);
    }

    #[test]
    fn rejects_bad_beads_and_separation() {
        assert_eq!(
            Bead::new(0.0, -1.0, 3.4, 2.0).unwrap_err().code(),
            "non_positive"
        );
        let b = neutral_bead();
        assert_eq!(
            BindingEnergy::two_bead_estimate(&b, &b, 0.0, 78.5, 0.0, 0.0072)
                .unwrap_err()
                .code(),
            "non_positive"
        );
    }
}
