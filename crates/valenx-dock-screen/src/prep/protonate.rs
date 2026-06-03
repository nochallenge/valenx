//! Feature 3 — protonation / partial-charge assignment.
//!
//! Docking needs two things from a ligand's chemistry: the right
//! *protonation state* at the assay pH (a carboxylate should be
//! deprotonated at pH 7, an aliphatic amine protonated), and a set of
//! *partial charges* for the electrostatics term of the AutoDock4
//! score.
//!
//! This module assigns both:
//!
//! - [`protonate_at_ph`] applies a small set of pKa rules to set
//!   formal charges for a target pH. A group whose pKa is below the
//!   pH loses a proton (acids); a group whose pKa is above the pH
//!   gains one (bases). This is the standard "major microspecies at
//!   pH X" rule used by ligand-prep tools — *not* a
//!   Poisson-Boltzmann pKa solver.
//! - [`assign_charges`] then runs [`valenx_cheminf::charge::gasteiger_charges`]
//!   on the protonated molecule (the Gasteiger PEOE method seeds from
//!   formal charge, so protonation feeds directly into it).
//!
//! The combined entry point is [`prepare_ligand`].

use valenx_cheminf::charge::gasteiger_charges;
use valenx_cheminf::molecule::{BondOrder, Molecule};

use crate::error::{DockScreenError, Result};

/// Physiological pH — the default protonation target.
pub const PHYSIOLOGICAL_PH: f64 = 7.4;

/// Which partial-charge model to use. Only Gasteiger is implemented as
/// a real v1; the variant exists so a future MMFF / AM1-BCC model can
/// be added without an API break.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum ChargeModel {
    /// Gasteiger-Marsili PEOE charges (the AutoDock default family).
    #[default]
    Gasteiger,
}

/// The result of preparing a ligand: a protonated molecule, its
/// per-atom partial charges and a record of which groups changed
/// protonation state.
#[derive(Clone, Debug, PartialEq)]
pub struct ProtonationResult {
    /// The molecule after protonation-state assignment. Formal charges
    /// reflect the major microspecies at the target pH.
    pub molecule: Molecule,
    /// Per-atom partial charges, parallel to `molecule.atoms`.
    pub charges: Vec<f64>,
    /// The pH the protonation rules were applied at.
    pub ph: f64,
    /// Number of acidic groups deprotonated.
    pub deprotonated: usize,
    /// Number of basic groups protonated.
    pub protonated: usize,
}

impl ProtonationResult {
    /// Net formal charge of the protonated molecule.
    pub fn net_charge(&self) -> i32 {
        self.molecule
            .atoms
            .iter()
            .map(|a| i32::from(a.formal_charge))
            .sum()
    }
}

/// Apply pH-aware protonation rules to a molecule, returning a new
/// molecule with updated formal charges plus a count of the changes.
///
/// The rules (representative, not exhaustive — the common
/// drug-relevant groups):
///
/// | group | pKa | at pH 7.4 |
/// |---|---|---|
/// | carboxylic acid –COOH | ~4.5 | deprotonated → –COO⁻ |
/// | aliphatic amine –NH₂ | ~10.5 | protonated → –NH₃⁺ |
/// | phosphate / phosphonate | ~2 / ~7 | deprotonated |
/// | tetrazole / sulfonic acid | ~1–4 | deprotonated |
///
/// Returns the protonated molecule and `(deprotonated, protonated)`
/// counts.
pub fn protonate_at_ph(mol: &Molecule, ph: f64) -> Result<(Molecule, usize, usize)> {
    if !ph.is_finite() || !(0.0..=14.0).contains(&ph) {
        return Err(DockScreenError::invalid(
            "ph",
            format!("pH must be a finite value in 0..=14, got {ph}"),
        ));
    }
    let mut out = mol.clone();
    let mut deprotonated = 0usize;
    let mut protonated = 0usize;

    for i in 0..out.atoms.len() {
        let z = out.atoms[i].atomic_number;
        // --- Acids: deprotonate when group pKa < pH ----------------
        // Carboxyl / phosphate / sulfonate oxygen: an O that is part
        // of an acid group and currently neutral.
        if z == 8 && out.atoms[i].formal_charge == 0 {
            if let Some(pka) = oxyacid_pka(&out, i) {
                if pka < ph {
                    out.atoms[i].formal_charge = -1;
                    // Drop the acidic proton if present as a node.
                    drop_one_hydrogen(&mut out, i);
                    deprotonated += 1;
                }
            }
        }
        // --- Bases: protonate when group pKa > pH ------------------
        // Aliphatic amine nitrogen: an sp3 N, not amide, not aromatic.
        if z == 7 && out.atoms[i].formal_charge == 0 {
            if let Some(pka) = amine_pka(&out, i) {
                if pka > ph {
                    out.atoms[i].formal_charge = 1;
                    protonated += 1;
                }
            }
        }
    }
    Ok((out, deprotonated, protonated))
}

/// Estimate the pKa of an oxygen that may be part of an oxyacid group
/// (carboxyl, phosphate, sulfonate). Returns `None` if the oxygen is
/// not a recognised acid oxygen.
fn oxyacid_pka(mol: &Molecule, oxygen: usize) -> Option<f64> {
    // The oxygen must be a hydroxyl-like terminal O: degree 1 (bonded
    // only to the acid centre) or carrying a hydrogen.
    let heavy_nbrs: Vec<usize> = mol
        .neighbors(oxygen)
        .into_iter()
        .filter(|&n| !mol.atoms[n].is_hydrogen())
        .collect();
    if heavy_nbrs.len() != 1 {
        return None;
    }
    let center = heavy_nbrs[0];
    let cz = mol.atoms[center].atomic_number;
    // The acid centre must itself bear a double-bonded oxygen (the
    // C=O / P=O / S=O that defines the acid).
    let has_carbonyl_o = mol.bonds.iter().any(|b| {
        b.touches(center)
            && b.order == BondOrder::Double
            && b.other(center).is_some_and(|o| {
                mol.atoms[o].atomic_number == 8 && o != oxygen
            })
    });
    if !has_carbonyl_o {
        return None;
    }
    match cz {
        6 => Some(4.5),  // carboxylic acid
        15 => Some(2.0), // phosphate / phosphonate
        16 => Some(1.5), // sulfonic acid
        _ => None,
    }
}

/// Estimate the pKa of a nitrogen that may be a protonatable amine.
/// Returns `None` for amide N, aromatic N, nitro N etc.
fn amine_pka(mol: &Molecule, nitrogen: usize) -> Option<f64> {
    if mol.atoms[nitrogen].aromatic {
        // Aromatic nitrogens (pyridine pKa ~5.2) — below pH 7.4, so
        // they stay neutral. Treat as non-protonatable for v1.
        return None;
    }
    // Amide N: bonded to a carbonyl carbon → not basic.
    let bonded_to_carbonyl = mol.neighbors(nitrogen).into_iter().any(|c| {
        mol.atoms[c].atomic_number == 6
            && mol.bonds.iter().any(|b| {
                b.touches(c)
                    && b.order == BondOrder::Double
                    && b.other(c).is_some_and(|o| mol.atoms[o].atomic_number == 8)
            })
    });
    if bonded_to_carbonyl {
        return None;
    }
    // Nitrogen double / triple bonded (imine, nitrile, nitro) — not a
    // simple amine.
    let has_multiple_bond = mol
        .bonds
        .iter()
        .any(|b| b.touches(nitrogen) && matches!(b.order, BondOrder::Double | BondOrder::Triple));
    if has_multiple_bond {
        return None;
    }
    // A plain sp3 amine — pKa ~10.5 (primary/secondary/tertiary all
    // sit comfortably above pH 7.4 for v1 purposes).
    Some(10.5)
}

/// Remove one explicit-hydrogen *node* bonded to `atom`, building a
/// fresh molecule. If the hydrogen is implicit (`Atom::implicit_h`),
/// decrement that count instead. No-op if there is no hydrogen.
fn drop_one_hydrogen(mol: &mut Molecule, atom: usize) {
    // Prefer dropping an explicit-H node so the graph stays exact.
    if let Some(h) = mol
        .neighbors(atom)
        .into_iter()
        .find(|&n| mol.atoms[n].is_hydrogen())
    {
        let keep: Vec<usize> = (0..mol.atoms.len()).filter(|&i| i != h).collect();
        *mol = mol.subgraph(&keep);
        return;
    }
    // Otherwise decrement the implicit / bracket-explicit count.
    let a = &mut mol.atoms[atom];
    if a.implicit_h > 0 {
        a.implicit_h -= 1;
    } else if a.explicit_h > 0 {
        a.explicit_h -= 1;
    }
}

/// Assign partial charges to a (already-protonated) molecule.
pub fn assign_charges(mol: &Molecule, model: ChargeModel) -> Vec<f64> {
    match model {
        ChargeModel::Gasteiger => gasteiger_charges(mol),
    }
}

/// Full ligand prep: protonate at `ph`, then assign charges with
/// `model`. The combined feature-3 entry point.
pub fn prepare_ligand(mol: &Molecule, ph: f64, model: ChargeModel) -> Result<ProtonationResult> {
    let (protonated_mol, deprotonated, protonated) = protonate_at_ph(mol, ph)?;
    let charges = assign_charges(&protonated_mol, model);
    Ok(ProtonationResult {
        molecule: protonated_mol,
        charges,
        ph,
        deprotonated,
        protonated,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cheminf::mol_from_smiles;

    #[test]
    fn carboxylic_acid_deprotonates_at_physiological_ph() {
        // Acetic acid CC(=O)O — the hydroxyl O should pick up a -1
        // formal charge at pH 7.4 (pKa ~4.5 < 7.4).
        let m = mol_from_smiles("CC(=O)O").unwrap();
        let (out, deprot, _prot) = protonate_at_ph(&m, PHYSIOLOGICAL_PH).unwrap();
        assert_eq!(deprot, 1, "expected one deprotonation");
        let net: i32 = out.atoms.iter().map(|a| i32::from(a.formal_charge)).sum();
        assert_eq!(net, -1, "carboxylate should carry net -1");
    }

    #[test]
    fn carboxylic_acid_stays_neutral_at_low_ph() {
        // At pH 2 the carboxyl (pKa ~4.5) is NOT deprotonated.
        let m = mol_from_smiles("CC(=O)O").unwrap();
        let (out, deprot, _) = protonate_at_ph(&m, 2.0).unwrap();
        assert_eq!(deprot, 0);
        let net: i32 = out.atoms.iter().map(|a| i32::from(a.formal_charge)).sum();
        assert_eq!(net, 0);
    }

    #[test]
    fn aliphatic_amine_protonates_at_physiological_ph() {
        // Methylamine CN — sp3 amine, pKa ~10.5 > 7.4 → protonated.
        let m = mol_from_smiles("CN").unwrap();
        let (out, _deprot, prot) = protonate_at_ph(&m, PHYSIOLOGICAL_PH).unwrap();
        assert_eq!(prot, 1, "expected one protonation");
        let net: i32 = out.atoms.iter().map(|a| i32::from(a.formal_charge)).sum();
        assert_eq!(net, 1, "ammonium should carry net +1");
    }

    #[test]
    fn amide_nitrogen_is_not_protonated() {
        // Acetamide CC(=O)N — the N is amide, must stay neutral.
        let m = mol_from_smiles("CC(=O)N").unwrap();
        let (out, _deprot, prot) = protonate_at_ph(&m, PHYSIOLOGICAL_PH).unwrap();
        assert_eq!(prot, 0, "amide N must not protonate");
        let net: i32 = out.atoms.iter().map(|a| i32::from(a.formal_charge)).sum();
        assert_eq!(net, 0);
    }

    #[test]
    fn protonate_rejects_out_of_range_ph() {
        let m = mol_from_smiles("CCO").unwrap();
        assert!(protonate_at_ph(&m, -1.0).is_err());
        assert!(protonate_at_ph(&m, 20.0).is_err());
        assert!(protonate_at_ph(&m, f64::NAN).is_err());
    }

    #[test]
    fn prepare_ligand_returns_charges_for_every_atom() {
        let m = mol_from_smiles("CC(=O)O").unwrap();
        let result = prepare_ligand(&m, PHYSIOLOGICAL_PH, ChargeModel::Gasteiger).unwrap();
        assert_eq!(result.charges.len(), result.molecule.atom_count());
        assert_eq!(result.ph, PHYSIOLOGICAL_PH);
        // The carboxylate is net -1.
        assert_eq!(result.net_charge(), -1);
    }

    #[test]
    fn charges_sum_to_net_formal_charge() {
        // Gasteiger conserves charge — the partial charges of a
        // protonated molecule sum to its net formal charge.
        let m = mol_from_smiles("CN").unwrap();
        let result = prepare_ligand(&m, PHYSIOLOGICAL_PH, ChargeModel::Gasteiger).unwrap();
        let total: f64 = result.charges.iter().sum();
        assert!(
            (total - result.net_charge() as f64).abs() < 1e-4,
            "charge total {total} should match net {}",
            result.net_charge()
        );
    }

    #[test]
    fn neutral_alcohol_is_untouched() {
        // Ethanol CCO — neither acid nor base, stays neutral.
        let m = mol_from_smiles("CCO").unwrap();
        let result = prepare_ligand(&m, PHYSIOLOGICAL_PH, ChargeModel::Gasteiger).unwrap();
        assert_eq!(result.deprotonated, 0);
        assert_eq!(result.protonated, 0);
        assert_eq!(result.net_charge(), 0);
    }
}
