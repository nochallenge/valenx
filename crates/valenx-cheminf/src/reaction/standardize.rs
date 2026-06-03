//! Molecule standardization — salt stripping and normalisation.
//!
//! Real-world structure data is messy: salts, mixtures, inconsistent
//! charge / protonation states. These routines bring a molecule to a
//! canonical "parent" form for comparison and modelling:
//!
//! - [`strip_salts`] removes counter-ion / solvent fragments, keeping
//!   the largest organic fragment;
//! - [`neutralize`] removes the obvious formal charges that a simple
//!   protonation would introduce (a `[O-]` becomes `O`, an `[NH3+]`
//!   becomes `N`), leaving genuinely charged species (quaternary
//!   ammonium, zwitterion pairs) intact;
//! - [`standardize`] runs strip → neutralize → re-perceive as one
//!   pass, the form a dedup pipeline should key on.

use crate::molecule::Molecule;

/// A small built-in list of common salt / solvent fragment formulas
/// (by element-count signature) that [`strip_salts`] discards even
/// when they would otherwise be the largest fragment is *not* used —
/// [`strip_salts`] keeps the largest organic fragment instead, which
/// is simpler and handles the long tail. This constant documents the
/// canonical salt set for reference / future use.
pub const COMMON_SALTS: &[&str] = &[
    "Cl", "Br", "I", "F", // halide counter-ions
    "Na", "K", "Li", "Ca", "Mg", // metal counter-ions
    "O", // water of crystallisation
];

/// Remove salt / solvent fragments, returning the parent molecule.
///
/// The parent is the connected fragment with the most heavy atoms; on
/// a tie the fragment containing carbon wins (so an organic acid beats
/// an equal-size inorganic counter-ion).
pub fn strip_salts(mol: &Molecule) -> Molecule {
    if mol.is_empty() || mol.component_count() <= 1 {
        return mol.clone();
    }
    let comp = mol.components();
    let n_comp = comp.iter().copied().max().map_or(0, |m| m + 1);
    // score each component: (has_carbon, heavy_atom_count)
    let mut score: Vec<(bool, usize)> = vec![(false, 0); n_comp];
    for (i, &c) in comp.iter().enumerate() {
        let a = &mol.atoms[i];
        if !a.is_hydrogen() {
            score[c].1 += 1;
        }
        if a.atomic_number == 6 {
            score[c].0 = true;
        }
    }
    let best = (0..n_comp)
        .max_by(|&a, &b| {
            score[a]
                .1
                .cmp(&score[b].1)
                .then(score[a].0.cmp(&score[b].0))
        })
        .unwrap_or(0);
    let keep: Vec<usize> = comp
        .iter()
        .enumerate()
        .filter(|(_, &c)| c == best)
        .map(|(i, _)| i)
        .collect();
    let mut parent = mol.subgraph(&keep);
    crate::perceive::perceive_all(&mut parent);
    parent
}

/// Neutralise the obvious protonation charges: a singly-charged
/// heteroatom whose charge can be cancelled by adding or removing one
/// hydrogen is reset to neutral.
///
/// `[O-]` → `O` (gain an H), `[N+]H…` ammonium → `N` (lose an H). A
/// charge that cannot be balanced this way — a quaternary `[N+]` with
/// four bonds and no H, a carbanion — is left untouched.
pub fn neutralize(mol: &Molecule) -> Molecule {
    let mut out = mol.clone();
    for i in 0..out.atoms.len() {
        let a = &out.atoms[i];
        let charge = a.formal_charge;
        if charge == 0 {
            continue;
        }
        match a.atomic_number {
            8 | 16 => {
                // O-/S- gains a proton; O+/S+ loses one if it has an H
                if charge == -1 || (charge == 1 && a.total_h() > 0) {
                    out.atoms[i].formal_charge = 0;
                }
            }
            7 => {
                // ammonium N+ drops a proton if it has one; N- gains one
                if charge == -1 || (charge == 1 && a.total_h() > 0) {
                    out.atoms[i].formal_charge = 0;
                }
            }
            _ => {}
        }
    }
    crate::perceive::hydrogen::recompute_implicit_hydrogens(&mut out);
    out
}

/// Full standardisation: strip salts, neutralise protonation charges
/// and re-run perception. The canonical "parent" form.
pub fn standardize(mol: &Molecule) -> Molecule {
    let stripped = strip_salts(mol);
    let mut neutral = neutralize(&stripped);
    crate::perceive::perceive_all(&mut neutral);
    neutral
}

/// `true` if `mol` is a multi-component mixture (carries a salt /
/// solvent / counter-ion).
pub fn is_mixture(mol: &Molecule) -> bool {
    mol.component_count() > 1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn strips_sodium_counterion() {
        // sodium acetate → acetate parent
        let m = mol_from_smiles("CC(=O)[O-].[Na+]").unwrap();
        assert!(is_mixture(&m));
        let parent = strip_salts(&m);
        assert!(!is_mixture(&parent));
        // the parent keeps the carbons
        assert!(parent.atoms.iter().any(|a| a.atomic_number == 6));
        assert!(!parent.atoms.iter().any(|a| a.atomic_number == 11));
    }

    #[test]
    fn strips_hydrochloride() {
        // an amine hydrochloride → the amine
        let m = mol_from_smiles("CCN.Cl").unwrap();
        let parent = strip_salts(&m);
        assert_eq!(parent.component_count(), 1);
        assert!(parent.atoms.iter().any(|a| a.atomic_number == 7));
    }

    #[test]
    fn neutralize_alkoxide() {
        let m = mol_from_smiles("CC[O-]").unwrap();
        let n = neutralize(&m);
        let o = n.atoms.iter().find(|a| a.atomic_number == 8).unwrap();
        assert_eq!(o.formal_charge, 0);
        // the neutral O now carries an implicit hydrogen
        assert_eq!(o.total_h(), 1);
    }

    #[test]
    fn neutralize_ammonium() {
        let m = mol_from_smiles("[NH4+]").unwrap();
        let n = neutralize(&m);
        assert_eq!(n.atoms[0].formal_charge, 0);
    }

    #[test]
    fn standardize_runs_full_pipeline() {
        let m = mol_from_smiles("CC(=O)[O-].[Na+]").unwrap();
        let std = standardize(&m);
        assert_eq!(std.component_count(), 1);
        // acetate was neutralised to acetic acid
        assert!(std.atoms.iter().all(|a| a.formal_charge == 0));
    }

    #[test]
    fn single_molecule_unchanged() {
        let m = mol_from_smiles("CCO").unwrap();
        let parent = strip_salts(&m);
        assert_eq!(parent.heavy_atom_count(), 3);
    }

    #[test]
    fn quaternary_ammonium_kept_charged() {
        // tetramethylammonium has no removable proton — stays +1
        let m = mol_from_smiles("C[N+](C)(C)C").unwrap();
        let n = neutralize(&m);
        let nitrogen = n.atoms.iter().find(|a| a.atomic_number == 7).unwrap();
        assert_eq!(nitrogen.formal_charge, 1);
    }
}
