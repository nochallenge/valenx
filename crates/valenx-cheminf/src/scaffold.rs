//! Molecular scaffolds and maximum common substructure.
//!
//! - [`murcko_scaffold`] extracts the Bemis-Murcko scaffold of a
//!   molecule: keep every ring atom and every atom on a *linker* path
//!   between two rings, drop the rest (the terminal side-chains). The
//!   result is the "molecular framework" used to cluster a library by
//!   core structure.
//! - [`generic_scaffold`] additionally flattens every scaffold atom to
//!   carbon and every bond to single — the Bemis-Murcko *graph*
//!   framework, which groups e.g. a pyridine and a benzene scaffold
//!   together.
//! - [`maximum_common_substructure`] finds the largest substructure
//!   shared by two molecules — a connected common subgraph by a
//!   greedy seed-and-extend over compatible atom pairs.
//!
//! **v1 simplification.** The MCS is a connected, greedy maximum
//! common *induced* subgraph keyed on element + bond order; it is not
//! the exhaustive maximum-common-edge-subgraph an NP-hard exact
//! algorithm (or RDKit's `rdFMCS`) would return. It finds a large
//! common core quickly and deterministically, which is what scaffold
//! analysis and series alignment need; it can miss the true optimum on
//! adversarial inputs.

use crate::molecule::{BondOrder, Molecule};
use crate::perceive::rings::sssr;

/// Extract the Bemis-Murcko scaffold of `mol` — ring systems plus the
/// linker atoms joining them, as a fresh molecule.
///
/// A molecule with no rings yields an empty scaffold.
pub fn murcko_scaffold(mol: &Molecule) -> Molecule {
    let keep = scaffold_atom_set(mol);
    mol.subgraph(&keep)
}

/// The Bemis-Murcko *generic* (graph) framework — the Murcko scaffold
/// with every atom set to carbon and every bond to single.
pub fn generic_scaffold(mol: &Molecule) -> Molecule {
    let mut s = murcko_scaffold(mol);
    for a in &mut s.atoms {
        a.atomic_number = 6;
        a.aromatic = false;
        a.formal_charge = 0;
        a.isotope = None;
        a.explicit_h = 0;
        a.implicit_h = 0;
    }
    for b in &mut s.bonds {
        b.order = BondOrder::Single;
        b.aromatic = false;
    }
    crate::perceive::hydrogen::recompute_implicit_hydrogens(&mut s);
    s
}

/// The set of atom indices that belong to the Murcko scaffold of
/// `mol`: every ring atom, plus every atom lying on a shortest path
/// between two ring systems (a linker).
pub fn scaffold_atom_set(mol: &Molecule) -> Vec<usize> {
    let rings = sssr(mol);
    let n = mol.atoms.len();
    let mut in_scaffold = vec![false; n];
    // all ring atoms
    for (i, flag) in in_scaffold.iter_mut().enumerate() {
        if rings.atom_in_ring(i) {
            *flag = true;
        }
    }
    if !in_scaffold.iter().any(|&x| x) {
        return Vec::new(); // acyclic → empty scaffold
    }

    // Iteratively prune: a terminal atom (degree ≤ 1 within the kept
    // set) that is not itself a ring atom is a side-chain leaf — drop
    // it. Repeat to a fixpoint. Start by keeping everything, then peel.
    let mut keep = vec![true; n];
    loop {
        let mut changed = false;
        for i in 0..n {
            if !keep[i] || in_scaffold[i] {
                continue; // ring atoms are never pruned
            }
            let live_degree = mol.neighbors(i).iter().filter(|&&v| keep[v]).count();
            if live_degree <= 1 {
                keep[i] = false;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    (0..n).filter(|&i| keep[i]).collect()
}

/// A maximum (greedy, connected) common substructure of two molecules,
/// returned as the list of `(atom_in_a, atom_in_b)` correspondences.
///
/// Atoms match when they share an element; bonds match when they share
/// a bond order. The result is a connected mapping; its length is the
/// MCS size in atoms.
pub fn maximum_common_substructure(a: &Molecule, b: &Molecule) -> Vec<(usize, usize)> {
    let mut best: Vec<(usize, usize)> = Vec::new();
    // try every compatible starting pair as a seed
    for ia in 0..a.atoms.len() {
        for ib in 0..b.atoms.len() {
            if !atoms_compatible(a, ia, b, ib) {
                continue;
            }
            let mapping = grow_mcs(a, b, (ia, ib));
            if mapping.len() > best.len() {
                best = mapping;
            }
        }
        // small early-out: if we already mapped all of `a`, stop
        if best.len() == a.atoms.len().min(b.atoms.len()) {
            break;
        }
    }
    best
}

/// The size (atom count) of the MCS — a similarity-style descriptor.
pub fn mcs_size(a: &Molecule, b: &Molecule) -> usize {
    maximum_common_substructure(a, b).len()
}

/// Greedily grow a connected common subgraph from one seed pair.
fn grow_mcs(a: &Molecule, b: &Molecule, seed: (usize, usize)) -> Vec<(usize, usize)> {
    let mut mapping = vec![seed];
    let mut used_a = vec![false; a.atoms.len()];
    let mut used_b = vec![false; b.atoms.len()];
    used_a[seed.0] = true;
    used_b[seed.1] = true;

    loop {
        // find the best extension: a bond from a mapped `a`-atom to an
        // unused `a`-atom that can be mirrored in `b`
        let mut chosen: Option<(usize, usize)> = None;
        'search: for &(ma, mb) in &mapping {
            for bi_a in a.bonds_on(ma) {
                let bond_a = &a.bonds[bi_a];
                let Some(na) = bond_a.other(ma) else {
                    continue;
                };
                if used_a[na] {
                    continue;
                }
                // need a `b`-bond from `mb` of the same order to a
                // compatible unused `b`-atom
                for bi_b in b.bonds_on(mb) {
                    let bond_b = &b.bonds[bi_b];
                    let Some(nb) = bond_b.other(mb) else {
                        continue;
                    };
                    if used_b[nb] {
                        continue;
                    }
                    if bond_a.order != bond_b.order {
                        continue;
                    }
                    if !atoms_compatible(a, na, b, nb) {
                        continue;
                    }
                    // also require that every already-mapped neighbour
                    // relationship is preserved (induced-subgraph
                    // consistency)
                    if mapping_consistent(a, b, &mapping, na, nb) {
                        chosen = Some((na, nb));
                        break 'search;
                    }
                }
            }
        }
        match chosen {
            Some((na, nb)) => {
                mapping.push((na, nb));
                used_a[na] = true;
                used_b[nb] = true;
            }
            None => break,
        }
    }
    mapping
}

/// Would adding `(na, nb)` keep the mapping a valid induced subgraph?
/// For every already-mapped pair, `na`-to-that-`a`-atom must be a bond
/// iff `nb`-to-that-`b`-atom is, and the orders must agree.
fn mapping_consistent(
    a: &Molecule,
    b: &Molecule,
    mapping: &[(usize, usize)],
    na: usize,
    nb: usize,
) -> bool {
    for &(ma, mb) in mapping {
        let ba = a.bond_between(na, ma);
        let bb = b.bond_between(nb, mb);
        match (ba, bb) {
            (Some(x), Some(y)) => {
                if a.bonds[x].order != b.bonds[y].order {
                    return false;
                }
            }
            (None, None) => {}
            _ => return false, // a bond on one side but not the other
        }
    }
    true
}

/// Atoms match for MCS / scaffold purposes: same element and same
/// aromatic state.
fn atoms_compatible(a: &Molecule, ia: usize, b: &Molecule, ib: usize) -> bool {
    let x = &a.atoms[ia];
    let y = &b.atoms[ib];
    x.atomic_number == y.atomic_number && x.aromatic == y.aromatic
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn benzene_scaffold_is_the_ring() {
        // toluene → the scaffold is just the benzene ring (drop CH3)
        let m = mol_from_smiles("Cc1ccccc1").unwrap();
        let s = murcko_scaffold(&m);
        assert_eq!(s.heavy_atom_count(), 6);
        assert!(s.atoms.iter().all(|a| a.atomic_number == 6));
    }

    #[test]
    fn acyclic_has_empty_scaffold() {
        let m = mol_from_smiles("CCCCCC").unwrap();
        let s = murcko_scaffold(&m);
        assert_eq!(s.atom_count(), 0);
    }

    #[test]
    fn biphenyl_keeps_linker() {
        // two rings joined directly — both rings survive
        let m = mol_from_smiles("c1ccccc1-c1ccccc1").unwrap();
        let s = murcko_scaffold(&m);
        assert_eq!(s.heavy_atom_count(), 12);
    }

    #[test]
    fn linker_atom_is_kept() {
        // two phenyl rings joined by a CH2 — the CH2 is a linker
        let m = mol_from_smiles("c1ccccc1Cc1ccccc1").unwrap();
        let s = murcko_scaffold(&m);
        // 12 ring carbons + 1 linker carbon
        assert_eq!(s.heavy_atom_count(), 13);
    }

    #[test]
    fn generic_scaffold_flattens() {
        let m = mol_from_smiles("c1ccncc1").unwrap(); // pyridine
        let g = generic_scaffold(&m);
        // every atom becomes carbon
        assert!(g.atoms.iter().all(|a| a.atomic_number == 6));
        assert_eq!(g.heavy_atom_count(), 6);
    }

    #[test]
    fn mcs_of_identical_is_whole_molecule() {
        let a = mol_from_smiles("CCO").unwrap();
        let b = mol_from_smiles("CCO").unwrap();
        assert_eq!(mcs_size(&a, &b), 3);
    }

    #[test]
    fn mcs_finds_shared_core() {
        // ethanol vs propanol share an ethanol-sized C-C-O core
        let a = mol_from_smiles("CCO").unwrap();
        let b = mol_from_smiles("CCCO").unwrap();
        let size = mcs_size(&a, &b);
        assert!(size >= 3, "shared C-C-O core, got {size}");
    }

    #[test]
    fn mcs_of_dissimilar_is_small() {
        let a = mol_from_smiles("CCO").unwrap();
        let b = mol_from_smiles("c1ccccc1").unwrap();
        // very little in common
        assert!(mcs_size(&a, &b) <= 2);
    }
}
