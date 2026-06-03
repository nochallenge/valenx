//! Aromaticity perception — the Hückel 4n+2 rule on ring systems.
//!
//! For each SSSR ring the perceiver counts the π electrons each ring
//! atom contributes:
//!
//! - an atom in a ring double bond contributes **1** (one π electron
//!   of the double bond);
//! - a heteroatom with a lone pair available for the ring (pyrrole-type
//!   N, furan O, thiophene S) contributes **2**;
//! - an `sp²` atom with an empty p-orbital (a carbocation centre, a
//!   boron) contributes **0**;
//! - an exocyclic double bond (`C=O` on the ring) contributes **0**.
//!
//! A ring is aromatic iff every atom is `sp²`-eligible and the total π
//! count is `4n + 2`. Fused systems are then re-checked: a ring all of
//! whose atoms already sit in an aromatic ring is itself flagged
//! aromatic (handles naphthalene / indole-type fusion).
//!
//! **v1 simplification:** this is per-SSSR-ring Hückel, not a full
//! "aromatic system" perception — antiaromatic exclusion and the rare
//! macrocyclic / azulene cases are handled approximately. Most drug-like
//! ring systems (benzene, pyridine, pyrrole, furan, thiophene,
//! imidazole, naphthalene, indole, quinoline, purine) are classified
//! correctly.

use super::rings::{Ring, RingInfo};
use crate::molecule::{BondOrder, Molecule};

/// Run aromaticity perception, mutating `mol`: every atom and bond that
/// the Hückel test places in an aromatic ring gets its `aromatic` flag
/// set; atoms / bonds outside aromatic rings get it cleared. Ring bonds
/// in an aromatic ring are upgraded to [`BondOrder::Aromatic`].
///
/// `rings` must be the SSSR of `mol` (see [`super::rings::sssr`]).
pub fn perceive_aromaticity(mol: &mut Molecule, rings: &RingInfo) {
    // start from a clean slate so re-perception is idempotent
    for a in &mut mol.atoms {
        a.aromatic = false;
    }
    for b in &mut mol.bonds {
        b.aromatic = false;
    }

    let mut aromatic_ring = vec![false; rings.rings.len()];
    // first pass: independent Hückel test on each ring
    for (ri, ring) in rings.rings.iter().enumerate() {
        if huckel_aromatic(mol, ring) {
            aromatic_ring[ri] = true;
        }
    }
    // Fusion pass: a ring fused to an already-aromatic ring — sharing a
    // bond with it — is itself aromatic provided every one of its atoms
    // is sp²-eligible. This is the standard perception of fused
    // polycyclic aromatics: naphthalene, anthracene, indole, quinoline,
    // purine. (A per-ring Hückel count cannot do this on its own — a
    // Kekulé fused system places the double bonds so that one SSSR ring
    // has only two ring-double-bonds within its own atoms, even though
    // the whole π system is aromatic.) Requiring every atom to be
    // sp²-eligible keeps a *saturated* fused ring — e.g. the CH₂ ring of
    // indane or tetralin — correctly non-aromatic. Iterated to a
    // fixpoint over the small SSSR set so a chain of fused rings all
    // light up.
    loop {
        let mut changed = false;
        for ri in 0..rings.rings.len() {
            if aromatic_ring[ri] {
                continue;
            }
            let ring = &rings.rings[ri];
            // every atom of this ring must be able to join a π system
            if !ring.atoms.iter().all(|&a| sp2_eligible(mol, a)) {
                continue;
            }
            // …and it must share at least one bond with an aromatic ring
            let shares_aromatic_bond = (0..rings.rings.len()).any(|rj| {
                aromatic_ring[rj]
                    && ring
                        .bonds
                        .iter()
                        .any(|&b| rings.rings[rj].bonds.contains(&b))
            });
            if shares_aromatic_bond {
                aromatic_ring[ri] = true;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }

    // apply flags
    for (ri, &is_arom) in aromatic_ring.iter().enumerate() {
        if !is_arom {
            continue;
        }
        let ring = &rings.rings[ri];
        for &a in &ring.atoms {
            mol.atoms[a].aromatic = true;
        }
        for &b in &ring.bonds {
            mol.bonds[b].aromatic = true;
            mol.bonds[b].order = BondOrder::Aromatic;
        }
    }
}

/// Hückel test for a single ring: count π electrons, require `4n+2`.
fn huckel_aromatic(mol: &Molecule, ring: &Ring) -> bool {
    if ring.size() < 3 || ring.size() > 8 {
        // restrict to common ring sizes; macrocyclic aromaticity is
        // out of v1 scope
        return false;
    }
    if !ring.atoms.iter().all(|&a| sp2_eligible(mol, a)) {
        return false;
    }
    let mut pi = 0i32;
    for &a in &ring.atoms {
        pi += pi_contribution(mol, a, ring);
        if pi < 0 {
            return false; // atom kills aromaticity
        }
    }
    pi >= 2 && (pi - 2) % 4 == 0
}

/// Is `atom` capable of being part of an aromatic π system at all?
/// Excludes `sp³` carbons (4 σ bonds + Hs), atoms with too many
/// neighbours, and atoms of elements that never aromatise.
fn sp2_eligible(mol: &Molecule, atom: usize) -> bool {
    let a = &mol.atoms[atom];
    // only the aromatic-capable elements
    if !crate::element::AROMATIC_ELEMENTS.contains(&a.atomic_number) {
        return false;
    }
    let degree = mol.degree(atom);
    // a ring atom in an aromatic ring has 2 or 3 heavy neighbours
    if !(2..=3).contains(&degree) {
        return false;
    }
    // a saturated sp3 carbon (4 single bonds incl. Hs and no double
    // bond) cannot be aromatic
    let has_multiple = mol
        .bonds_on(atom)
        .iter()
        .any(|&bi| matches!(mol.bonds[bi].order, BondOrder::Double | BondOrder::Aromatic));
    let total_bonds = mol.explicit_valence(atom) + f64::from(a.total_h());
    if a.atomic_number == 6 && !has_multiple && total_bonds >= 4.0 {
        return false;
    }
    true
}

/// π-electron contribution of `atom` to `ring`. Negative means the
/// atom prevents aromaticity (an `sp³` interruption).
fn pi_contribution(mol: &Molecule, atom: usize, ring: &Ring) -> i32 {
    let a = &mol.atoms[atom];
    // a ring double / aromatic bond → this atom donates 1 π electron
    let in_ring_pi = mol.bonds_on(atom).iter().any(|&bi| {
        ring.contains_bond(bi)
            && matches!(mol.bonds[bi].order, BondOrder::Double | BondOrder::Aromatic)
    });
    // an exocyclic double bond (C=O on the ring) → carbonyl carbon
    // donates 0 (its p-orbital is in the C=O π, not the ring)
    let exocyclic_double = mol.bonds_on(atom).iter().any(|&bi| {
        !ring.contains_bond(bi) && mol.bonds[bi].order == BondOrder::Double
    });

    match a.atomic_number {
        6 => {
            // carbon
            if in_ring_pi {
                1
            } else if exocyclic_double {
                0 // carbonyl-type sp2 carbon
            } else if a.formal_charge < 0 {
                2 // carbanion lone pair (cyclopentadienyl)
            } else if a.formal_charge > 0 {
                0 // carbocation empty p-orbital (tropylium)
            } else {
                // an sp3 CH2 in the ring breaks aromaticity
                -100
            }
        }
        7 => {
            // nitrogen
            if in_ring_pi {
                1 // pyridine-type
            } else {
                2 // pyrrole-type lone pair
            }
        }
        8 | 16 | 34 | 52 => {
            // O / S / Se / Te — furan / thiophene heteroatom lone pair
            if in_ring_pi {
                1
            } else {
                2
            }
        }
        5 => 0, // boron — empty p-orbital
        15 => {
            if in_ring_pi {
                1
            } else {
                2
            }
        }
        _ => -100,
    }
}

/// Convenience: is the molecule's atom `atom` aromatic after
/// perception? (a thin alias for `mol.atoms[atom].aromatic`).
pub fn is_aromatic_atom(mol: &Molecule, atom: usize) -> bool {
    mol.atoms.get(atom).map(|a| a.aromatic).unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::super::rings::sssr;
    use super::*;
    use crate::smiles::parse_smiles;

    fn aromatic_atom_count(smiles: &str) -> usize {
        let mut m = parse_smiles(smiles).unwrap();
        let r = sssr(&m);
        perceive_aromaticity(&mut m, &r);
        m.atoms.iter().filter(|a| a.aromatic).count()
    }

    #[test]
    fn benzene_is_aromatic() {
        assert_eq!(aromatic_atom_count("C1=CC=CC=C1"), 6);
    }

    #[test]
    fn cyclohexane_is_not() {
        assert_eq!(aromatic_atom_count("C1CCCCC1"), 0);
    }

    #[test]
    fn pyridine_kekule_aromatic() {
        assert_eq!(aromatic_atom_count("C1=CC=NC=C1"), 6);
    }

    #[test]
    fn pyrrole_aromatic() {
        // pyrrole: 4 carbons (2 double bonds → 4 pi) + N lone pair (2)
        assert_eq!(aromatic_atom_count("C1=CC=CN1"), 5);
    }

    #[test]
    fn furan_aromatic() {
        assert_eq!(aromatic_atom_count("C1=CC=CO1"), 5);
    }

    #[test]
    fn cyclobutadiene_not_aromatic() {
        // 4 pi electrons → 4n, antiaromatic, must NOT be flagged
        assert_eq!(aromatic_atom_count("C1=CC=C1"), 0);
    }

    #[test]
    fn naphthalene_all_aromatic() {
        assert_eq!(aromatic_atom_count("C1=CC=C2C=CC=CC2=C1"), 10);
    }

    #[test]
    fn aromatic_smiles_input_stays_aromatic() {
        // lowercase input — perception should keep all 6 aromatic
        assert_eq!(aromatic_atom_count("c1ccccc1"), 6);
    }
}
