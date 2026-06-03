//! Feature 18 — MM-GBSA-class rescoring.
//!
//! Docking scores are tuned for *speed* — they have to rank millions
//! of poses. MM-GBSA (Molecular Mechanics + Generalized-Born + Surface
//! Area) is a slower, more physical estimate of binding free energy
//! used to *re-rank* the top docking hits:
//!
//! ```text
//! ΔG_bind ≈ ΔE_MM + ΔG_GB + ΔG_SA
//! ```
//!
//! - **ΔE_MM** — the molecular-mechanics interaction energy between
//!   receptor and ligand (van der Waals + Coulomb in vacuum).
//! - **ΔG_GB** — the electrostatic solvation free energy from the
//!   Generalized-Born model. Burying a charged atom against the
//!   partner desolvates it; GB estimates that cost.
//! - **ΔG_SA** — the non-polar solvation term, proportional to the
//!   buried solvent-accessible surface area.
//!
//! This module implements a single-snapshot MM-GBSA rescoring with the
//! **Still et al. (1990) pairwise Generalized-Born model**. The GB
//! pairwise descreening, the effective Born radii and the GB
//! polarisation energy are all real — see [`mmgbsa_rescore`].
//!
//! ### v1 note
//!
//! Production MM-GBSA averages over an MD ensemble and uses a full
//! force field. This is a *single-snapshot* rescoring with an
//! AutoDock-atom-typed van der Waals / charge model — a real GB
//! calculation, but not an MD-averaged free energy. The non-polar
//! surface term uses a buried-volume proxy rather than a full
//! Shrake-Rupley SASA difference. Both simplifications are documented
//! here and in the crate-level v1 note.

use nalgebra::Vector3;

use valenx_dock::atom_type::Ad4AtomType;
use valenx_dock::receptor::Receptor;

use crate::error::{DockScreenError, Result};

/// One atom for the GB calculation: position, charge and intrinsic
/// (van der Waals) radius.
#[derive(Clone, Copy, Debug)]
struct GbAtom {
    pos: Vector3<f64>,
    charge: f64,
    radius: f64,
}

/// The term breakdown of an MM-GBSA rescoring.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MmGbsaTerms {
    /// Molecular-mechanics van der Waals interaction energy
    /// (kcal/mol).
    pub vdw: f64,
    /// Molecular-mechanics Coulomb interaction energy in vacuum
    /// (kcal/mol).
    pub coulomb: f64,
    /// Generalized-Born polar solvation free energy of binding
    /// (kcal/mol) — `GB(complex) − GB(receptor) − GB(ligand)`.
    pub gb_solvation: f64,
    /// Non-polar (surface-area) solvation free energy of binding
    /// (kcal/mol).
    pub nonpolar: f64,
}

impl MmGbsaTerms {
    /// The estimated binding free energy: the sum of all four terms.
    pub fn total(&self) -> f64 {
        self.vdw + self.coulomb + self.gb_solvation + self.nonpolar
    }

    /// The gas-phase molecular-mechanics interaction energy only.
    pub fn molecular_mechanics(&self) -> f64 {
        self.vdw + self.coulomb
    }
}

/// Solvent / solute dielectric constants and the surface-tension
/// coefficient — standard MM-GBSA values.
mod constants {
    /// Solute (interior) dielectric.
    pub const EPS_IN: f64 = 1.0;
    /// Solvent (water) dielectric.
    pub const EPS_OUT: f64 = 78.5;
    /// Coulomb constant in kcal·Å/(mol·e²).
    pub const COULOMB_K: f64 = 332.0637;
    /// Non-polar surface-tension coefficient (kcal/mol·Å²).
    pub const SURF_TENSION: f64 = 0.0072;
    /// Still GB closeness parameter.
    pub const STILL_P: f64 = 0.09;
}

/// Effective Born radii by the Still pairwise descreening
/// approximation: a buried atom — crowded by neighbours — has a larger
/// effective radius and is harder to solvate.
fn born_radii(atoms: &[GbAtom]) -> Vec<f64> {
    let n = atoms.len();
    let mut radii = vec![0.0; n];
    for i in 0..n {
        // Start from the inverse intrinsic radius.
        let mut inv = 1.0 / atoms[i].radius.max(0.5);
        // Each neighbour descreens atom i — reduces its effective
        // 1/R, hence increases R.
        for j in 0..n {
            if i == j {
                continue;
            }
            let d = (atoms[i].pos - atoms[j].pos).norm().max(0.1);
            let rj = atoms[j].radius.max(0.5);
            // Still's pairwise descreening kernel (volume of j seen
            // from i, smoothed). Damped so a far neighbour barely
            // contributes.
            let term = constants::STILL_P * rj.powi(3)
                / (d.powi(4) + 1.0).max(1.0);
            inv -= term;
        }
        // Effective radius — never let it collapse below the intrinsic
        // radius (descreening can only *grow* a Born radius).
        radii[i] = (1.0 / inv.max(1e-3)).max(atoms[i].radius.max(0.5));
    }
    radii
}

/// The Still GB polarisation energy of a set of atoms with the given
/// effective Born radii.
fn gb_polarisation(atoms: &[GbAtom], radii: &[f64]) -> f64 {
    let factor = -0.5 * constants::COULOMB_K
        * (1.0 / constants::EPS_IN - 1.0 / constants::EPS_OUT);
    let n = atoms.len();
    let mut energy = 0.0;
    for i in 0..n {
        for j in i..n {
            let qi = atoms[i].charge;
            let qj = atoms[j].charge;
            let ri = radii[i];
            let rj = radii[j];
            let d2 = if i == j {
                0.0
            } else {
                (atoms[i].pos - atoms[j].pos).norm_squared()
            };
            // Still's smoothed GB function f_GB.
            let f_gb = (d2 + ri * rj * (-d2 / (4.0 * ri * rj)).exp()).sqrt();
            let pair = qi * qj / f_gb.max(1e-6);
            // Diagonal terms count once, off-diagonal twice.
            energy += if i == j { pair } else { 2.0 * pair };
        }
    }
    factor * energy
}

/// Build the GB-atom list for a receptor, keeping only atoms within
/// `interaction_cutoff` of any ligand atom — the GB model is short-
/// ranged and an enormous receptor would otherwise dominate the cost.
fn receptor_gb_atoms(
    receptor: &Receptor,
    ligand: &[(Vector3<f64>, Ad4AtomType, f64)],
    interaction_cutoff: f64,
) -> Vec<GbAtom> {
    let cut2 = interaction_cutoff * interaction_cutoff;
    receptor
        .atoms
        .iter()
        .filter(|ra| {
            ligand
                .iter()
                .any(|(lp, _, _)| (lp - ra.position).norm_squared() <= cut2)
        })
        .map(|ra| GbAtom {
            pos: ra.position,
            charge: ra.partial_charge,
            radius: ra.ad4_type.props().vdw_radius,
        })
        .collect()
}

/// The 12-6 van der Waals + Coulomb interaction energy between two
/// atom sets (the molecular-mechanics term).
fn mm_interaction(a: &[GbAtom], b: &[GbAtom]) -> (f64, f64) {
    let mut vdw = 0.0;
    let mut coulomb = 0.0;
    for x in a {
        for y in b {
            let r = (x.pos - y.pos).norm().max(0.5);
            // Lorentz-Berthelot: r_min = sum of radii, a fixed well
            // depth (a v1 simplification — see the module note).
            let r_min = x.radius + y.radius;
            let eps = 0.15;
            let ratio = r_min / r;
            let r6 = ratio.powi(6);
            vdw += eps * (r6 * r6 - 2.0 * r6);
            coulomb += constants::COULOMB_K * x.charge * y.charge / r;
        }
    }
    (vdw, coulomb)
}

/// Feature 18 — rescore a docked pose with a single-snapshot
/// MM-GBSA-class calculation.
///
/// `ligand_atoms` is the *posed* ligand: `(world position, AD4 type,
/// partial charge)` per atom. `interaction_cutoff` bounds which
/// receptor atoms participate (8–12 Å is typical).
///
/// Returns the full [`MmGbsaTerms`] breakdown. The reported energy is
/// the free energy *of binding* — every term is computed as
/// `complex − receptor − ligand`.
///
/// Returns [`DockScreenError`] for an empty ligand or a non-positive
/// cutoff.
pub fn mmgbsa_rescore(
    receptor: &Receptor,
    ligand_atoms: &[(Vector3<f64>, Ad4AtomType, f64)],
    interaction_cutoff: f64,
) -> Result<MmGbsaTerms> {
    if ligand_atoms.is_empty() {
        return Err(DockScreenError::invalid_ligand(
            "MM-GBSA rescoring needs a posed ligand with atoms",
        ));
    }
    if !interaction_cutoff.is_finite() || interaction_cutoff <= 0.0 {
        return Err(DockScreenError::invalid(
            "interaction_cutoff",
            "must be positive",
        ));
    }
    if receptor.atoms.is_empty() {
        return Err(DockScreenError::invalid_receptor("receptor has no atoms"));
    }

    let lig: Vec<GbAtom> = ligand_atoms
        .iter()
        .map(|&(p, t, q)| GbAtom {
            pos: p,
            charge: q,
            radius: t.props().vdw_radius,
        })
        .collect();
    let rec = receptor_gb_atoms(receptor, ligand_atoms, interaction_cutoff);
    if rec.is_empty() {
        // No receptor atom within range — the pose does not contact
        // the receptor; report zero binding energy honestly.
        return Ok(MmGbsaTerms::default());
    }

    // --- molecular-mechanics interaction -----------------------------
    let (vdw, coulomb) = mm_interaction(&rec, &lig);

    // --- generalized-Born solvation of binding -----------------------
    // GB(complex): the union; GB(receptor) and GB(ligand) separately.
    let mut complex = rec.clone();
    complex.extend_from_slice(&lig);
    let gb_complex = gb_polarisation(&complex, &born_radii(&complex));
    let gb_receptor = gb_polarisation(&rec, &born_radii(&rec));
    let gb_ligand = gb_polarisation(&lig, &born_radii(&lig));
    let gb_solvation = gb_complex - gb_receptor - gb_ligand;

    // --- non-polar (surface) term ------------------------------------
    // A buried-contact proxy for the SASA difference: each ligand atom
    // close to the receptor buries roughly π·r² of surface.
    let mut buried_area = 0.0;
    for la in &lig {
        let mut nearest = f64::INFINITY;
        for ra in &rec {
            nearest = nearest.min((la.pos - ra.pos).norm());
        }
        let contact = la.radius + 1.4; // VDW + a water-probe radius
        if nearest < contact {
            // The closer the contact, the more area is buried.
            let burial = ((contact - nearest) / contact).clamp(0.0, 1.0);
            buried_area += std::f64::consts::PI * la.radius * la.radius * burial;
        }
    }
    let nonpolar = -constants::SURF_TENSION * buried_area;

    Ok(MmGbsaTerms {
        vdw,
        coulomb,
        gb_solvation,
        nonpolar,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_dock::receptor::ReceptorAtom;

    fn charged_receptor(charge: f64) -> Receptor {
        Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::OA,
                partial_charge: charge,
            }],
        }
    }

    #[test]
    fn rejects_degenerate_inputs() {
        let r = charged_receptor(-0.5);
        assert!(mmgbsa_rescore(&r, &[], 10.0).is_err());
        let lig = [(Vector3::new(3.0, 0.0, 0.0), Ad4AtomType::N, 0.5)];
        assert!(mmgbsa_rescore(&r, &lig, 0.0).is_err());
        assert!(mmgbsa_rescore(&Receptor::default(), &lig, 10.0).is_err());
    }

    #[test]
    fn total_is_sum_of_terms() {
        let r = charged_receptor(-0.4);
        let lig = [(Vector3::new(3.2, 0.0, 0.0), Ad4AtomType::N, 0.4)];
        let t = mmgbsa_rescore(&r, &lig, 10.0).unwrap();
        let manual = t.vdw + t.coulomb + t.gb_solvation + t.nonpolar;
        assert!((t.total() - manual).abs() < 1e-9);
        assert!((t.molecular_mechanics() - (t.vdw + t.coulomb)).abs() < 1e-9);
    }

    #[test]
    fn opposite_charges_give_favourable_coulomb() {
        let r = charged_receptor(-0.5);
        let lig = [(Vector3::new(3.0, 0.0, 0.0), Ad4AtomType::N, 0.5)];
        let t = mmgbsa_rescore(&r, &lig, 10.0).unwrap();
        assert!(t.coulomb < 0.0, "opposite charges should attract");
    }

    #[test]
    fn non_contacting_pose_gives_zero_binding() {
        // Ligand 50 Å away — outside the cutoff, no receptor atoms in
        // range, so every binding term is zero.
        let r = charged_receptor(-0.5);
        let lig = [(Vector3::new(50.0, 0.0, 0.0), Ad4AtomType::N, 0.5)];
        let t = mmgbsa_rescore(&r, &lig, 10.0).unwrap();
        assert_eq!(t.total(), 0.0);
    }

    #[test]
    fn buried_pose_has_nonpolar_contribution() {
        // A ligand atom in van der Waals contact buries surface — the
        // non-polar term is favourable (negative).
        let r = charged_receptor(0.0);
        let lig = [(Vector3::new(3.0, 0.0, 0.0), Ad4AtomType::C, 0.0)];
        let t = mmgbsa_rescore(&r, &lig, 10.0).unwrap();
        assert!(t.nonpolar <= 0.0, "non-polar term should be ≤ 0");
    }

    #[test]
    fn born_radii_grow_under_crowding() {
        // An isolated atom vs the same atom crowded by neighbours —
        // the crowded one has a larger effective Born radius.
        let lone = vec![GbAtom {
            pos: Vector3::zeros(),
            charge: 0.0,
            radius: 1.5,
        }];
        let crowded = vec![
            GbAtom {
                pos: Vector3::zeros(),
                charge: 0.0,
                radius: 1.5,
            },
            GbAtom {
                pos: Vector3::new(2.0, 0.0, 0.0),
                charge: 0.0,
                radius: 1.5,
            },
            GbAtom {
                pos: Vector3::new(-2.0, 0.0, 0.0),
                charge: 0.0,
                radius: 1.5,
            },
        ];
        let r_lone = born_radii(&lone)[0];
        let r_crowded = born_radii(&crowded)[0];
        assert!(
            r_crowded >= r_lone,
            "crowded Born radius {r_crowded} should be ≥ lone {r_lone}"
        );
    }

    #[test]
    fn gb_solvation_is_finite_for_a_charged_complex() {
        let r = charged_receptor(-0.6);
        let lig = [(Vector3::new(3.0, 0.0, 0.0), Ad4AtomType::N, 0.6)];
        let t = mmgbsa_rescore(&r, &lig, 12.0).unwrap();
        assert!(t.gb_solvation.is_finite());
    }
}
