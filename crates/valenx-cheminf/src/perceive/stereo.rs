//! Stereochemistry perception — tetrahedral R/S and double-bond E/Z.
//!
//! Two jobs:
//!
//! - [`assign_cip_labels`] looks at every atom carrying a parsed
//!   [`Chirality::Cw`] / [`Chirality::Ccw`] tag, ranks its neighbours
//!   by CIP priority, and converts the as-written parity into an
//!   absolute [`Chirality::R`] / [`Chirality::S`] label.
//! - [`assign_double_bond_stereo`] looks at every double bond whose
//!   substituent single bonds carry SMILES direction marks (`/` `\` →
//!   [`BondStereo::Up`] / [`BondStereo::Down`]) and converts them to an
//!   absolute [`BondStereo::E`] / [`BondStereo::Z`] label.
//!
//! **CIP priority — v1 rule.** The substituent ranking is a real
//! recursive CIP comparison: compare atomic number; on a tie compare
//! the sorted multiset of neighbour atomic numbers; recurse breadth-
//! first to a bounded depth. Duplicate atoms for double bonds are
//! handled by phantom neighbours. This resolves the overwhelming
//! majority of real centres correctly. It does **not** implement the
//! full 2013 IUPAC CIP rules (isotope tie-breaks at rule 2, the
//! like/unlike rules 4–5, auxiliary descriptors) — those rare cases
//! are documented as a v1 limitation.

use crate::molecule::{BondOrder, BondStereo, Chirality, Molecule};

/// Per-atom CIP priority used when ranking the neighbours of a stereo
/// centre. Higher = higher priority.
fn cip_rank_pair(mol: &Molecule, root: usize, nbr: usize) -> CipKey {
    // BFS sphere expansion from `nbr`, never stepping back onto `root`
    // through the first bond. Each sphere is a sorted descending list
    // of atomic numbers; the keys compare lexicographically.
    const MAX_DEPTH: usize = 6;
    let n = mol.atoms.len();
    let mut spheres: Vec<Vec<u16>> = Vec::with_capacity(MAX_DEPTH);
    // frontier holds (atom, came_from) so we can add double-bond
    // phantom duplicates
    let mut frontier: Vec<(usize, usize)> = vec![(nbr, root)];
    let mut visited = vec![false; n];
    visited[root] = true;
    visited[nbr] = true;
    for _ in 0..MAX_DEPTH {
        if frontier.is_empty() {
            break;
        }
        let mut sphere: Vec<u16> = Vec::new();
        let mut next: Vec<(usize, usize)> = Vec::new();
        for &(atom, _from) in &frontier {
            let a = &mol.atoms[atom];
            // record the atom itself plus its hydrogens (atomic no. 1)
            sphere.push(u16::from(a.atomic_number));
            sphere.extend(std::iter::repeat_n(1u16, usize::from(a.total_h())));
            for bi in mol.bonds_on(atom) {
                let bond = &mol.bonds[bi];
                if let Some(v) = bond.other(atom) {
                    // double / triple bonds add phantom duplicate atoms
                    let multiplicity = match bond.order {
                        BondOrder::Double => 2,
                        BondOrder::Triple => 3,
                        _ => 1,
                    };
                    for k in 0..multiplicity {
                        if k == 0 {
                            if !visited[v] {
                                visited[v] = true;
                                next.push((v, atom));
                            }
                        } else {
                            // phantom duplicate of v — contributes its
                            // atomic number to the next sphere but is
                            // not traversed further
                            sphere.push(u16::from(mol.atoms[v].atomic_number));
                        }
                    }
                }
            }
        }
        sphere.sort_unstable_by(|a, b| b.cmp(a));
        spheres.push(sphere);
        frontier = next;
    }
    CipKey { spheres }
}

/// A comparable CIP key: a list of spheres, each a descending atomic-
/// number multiset. Compared sphere-by-sphere, then element-by-element.
#[derive(Clone, Debug, PartialEq, Eq)]
struct CipKey {
    spheres: Vec<Vec<u16>>,
}

impl PartialOrd for CipKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for CipKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        let depth = self.spheres.len().max(other.spheres.len());
        for d in 0..depth {
            let empty = Vec::new();
            let a = self.spheres.get(d).unwrap_or(&empty);
            let b = other.spheres.get(d).unwrap_or(&empty);
            // compare the descending multisets lexicographically
            let len = a.len().max(b.len());
            for i in 0..len {
                let av = a.get(i).copied().unwrap_or(0);
                let bv = b.get(i).copied().unwrap_or(0);
                match av.cmp(&bv) {
                    Ordering::Equal => {}
                    ord => return ord,
                }
            }
        }
        Ordering::Equal
    }
}

/// Rank the heavy + hydrogen neighbours of `center` by CIP priority,
/// returning neighbour atom indices ordered **highest priority first**.
/// Implicit hydrogens are represented by a sentinel `usize::MAX`.
pub fn cip_neighbor_order(mol: &Molecule, center: usize) -> Vec<usize> {
    let mut nbrs: Vec<usize> = mol.neighbors(center);
    // Hydrogens on the centre — explicit (bracket `[C@H]`) or implicit —
    // count as lowest-priority phantom neighbours; use the TOTAL H count.
    let n_h = mol.atoms[center].total_h();
    let mut keyed: Vec<(CipKey, usize)> = nbrs
        .drain(..)
        .map(|v| (cip_rank_pair(mol, center, v), v))
        .collect();
    for _ in 0..n_h {
        keyed.push((
            CipKey {
                spheres: vec![vec![1]],
            },
            usize::MAX,
        ));
    }
    // highest priority first
    keyed.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    keyed.into_iter().map(|(_, v)| v).collect()
}

/// Assign absolute R/S labels to every atom with a parsed `@` / `@@`
/// tetrahedral tag. Atoms already labelled R/S, or with no tag, are
/// left alone. Mutates `mol`.
pub fn assign_cip_labels(mol: &mut Molecule) {
    let n = mol.atoms.len();
    for center in 0..n {
        let parsed = mol.atoms[center].chirality;
        let parity = match parsed {
            Chirality::Cw | Chirality::Ccw => parsed,
            _ => continue,
        };
        // A real tetrahedral centre needs 3 heavy neighbours + 1 H, or
        // 4 heavy neighbours. The hydrogen on a stereocentre may be
        // written explicitly inside a bracket (`[C@H]` → `explicit_h`)
        // or be valence-implicit (`implicit_h`); both must count, so the
        // 4-coordination test uses the TOTAL hydrogen count, not just
        // `implicit_h`.
        let degree = mol.degree(center);
        let n_h = usize::from(mol.atoms[center].total_h());
        if degree + n_h != 4 {
            // not a classic 4-coordinate centre — leave the raw tag
            continue;
        }
        // neighbours in the order they appear in the bond list — that
        // is the SMILES "as written" order the parity refers to
        let written: Vec<usize> = {
            let mut v: Vec<usize> = mol
                .bonds
                .iter()
                .filter_map(|b| b.other(center))
                .collect();
            for _ in 0..n_h {
                v.push(usize::MAX);
            }
            v
        };
        let priority = cip_neighbor_order(mol, center);
        // map each written neighbour to its priority rank (0 = highest)
        let rank = |atom: usize| -> usize {
            priority.iter().position(|&p| p == atom).unwrap_or(3)
        };
        let perm: Vec<usize> = written.iter().map(|&a| rank(a)).collect();
        // parity of the permutation that sorts `perm` ascending
        let swaps = permutation_parity(&perm);
        // `@` (Ccw) with an even sort-parity → looking from lowest
        // priority, the 1-2-3 order is anticlockwise → S; flip on odd.
        let base_ccw = swaps % 2 == 0;
        let is_ccw = match parity {
            Chirality::Ccw => base_ccw,
            Chirality::Cw => !base_ccw,
            _ => unreachable!(),
        };
        mol.atoms[center].chirality = if is_ccw { Chirality::S } else { Chirality::R };
    }
}

/// Parity (number of transpositions, mod 2 not taken) needed to sort a
/// permutation of `0..k` ascending — counts inversions.
fn permutation_parity(perm: &[usize]) -> usize {
    let mut inv = 0usize;
    for i in 0..perm.len() {
        for j in i + 1..perm.len() {
            if perm[i] > perm[j] {
                inv += 1;
            }
        }
    }
    inv
}

/// Assign absolute E/Z to every double bond whose two single-bond
/// substituents carry SMILES `/` `\` direction marks. Mutates `mol`.
///
/// The rule: for `F/C=C/F` the two reference bonds point the *same*
/// way around the double bond → the high-priority substituents are
/// trans → `E`; `F/C=C\F` → cis → `Z`.
pub fn assign_double_bond_stereo(mol: &mut Molecule) {
    let double_bonds: Vec<usize> = mol
        .bonds
        .iter()
        .enumerate()
        .filter(|(_, b)| b.order == BondOrder::Double)
        .map(|(i, _)| i)
        .collect();

    for di in double_bonds {
        let (c1, c2) = (mol.bonds[di].a, mol.bonds[di].b);
        // find a directional reference bond on each end
        let ref1 = directional_ref(mol, c1, di);
        let ref2 = directional_ref(mol, c2, di);
        let (Some((b1, dir1, sub1)), Some((b2, dir2, sub2))) = (ref1, ref2) else {
            continue;
        };
        let _ = (b1, b2);
        // normalise direction so it is "as seen from c1 / c2"
        // dir is Up if the slash points up *into* the double-bond atom
        let up1 = dir1 == BondStereo::Up;
        let up2 = dir2 == BondStereo::Up;
        // same orientation of the two reference single bonds → the two
        // reference substituents are on opposite sides (trans)
        let ref_trans = up1 == up2;

        // are the reference substituents the CIP-highest on each end?
        let high1 = highest_priority_substituent(mol, c1, c2);
        let high2 = highest_priority_substituent(mol, c2, c1);
        let ref1_is_high = high1 == Some(sub1);
        let ref2_is_high = high2 == Some(sub2);
        // flip the trans/cis sense once for each end whose reference
        // substituent is NOT the high-priority one
        let mut high_trans = ref_trans;
        if !ref1_is_high {
            high_trans = !high_trans;
        }
        if !ref2_is_high {
            high_trans = !high_trans;
        }
        mol.bonds[di].stereo = if high_trans {
            BondStereo::E
        } else {
            BondStereo::Z
        };
    }
}

/// Find a directional (`/` `\`) single bond on `atom` other than the
/// double bond `db`. Returns `(bond_index, direction, other_atom)`.
fn directional_ref(
    mol: &Molecule,
    atom: usize,
    db: usize,
) -> Option<(usize, BondStereo, usize)> {
    for bi in mol.bonds_on(atom) {
        if bi == db {
            continue;
        }
        let b = &mol.bonds[bi];
        if matches!(b.stereo, BondStereo::Up | BondStereo::Down) {
            let other = b.other(atom)?;
            // direction is recorded relative to the bond's `a`->`b`;
            // normalise so "Up" means pointing toward `atom`
            let dir = if b.a == atom {
                // bond written atom -> other, slash as recorded
                b.stereo
            } else {
                // bond written other -> atom, slash sense reverses
                match b.stereo {
                    BondStereo::Up => BondStereo::Down,
                    BondStereo::Down => BondStereo::Up,
                    s => s,
                }
            };
            return Some((bi, dir, other));
        }
    }
    None
}

/// The CIP-highest-priority substituent of `atom` that is not the other
/// double-bond carbon `partner`.
fn highest_priority_substituent(
    mol: &Molecule,
    atom: usize,
    partner: usize,
) -> Option<usize> {
    let order = cip_neighbor_order(mol, atom);
    order
        .into_iter()
        .find(|&v| v != partner && v != usize::MAX)
}

/// Count the perceived tetrahedral stereo centres (atoms with an R or S
/// label after [`assign_cip_labels`]).
pub fn stereocenter_count(mol: &Molecule) -> usize {
    mol.atoms
        .iter()
        .filter(|a| matches!(a.chirality, Chirality::R | Chirality::S))
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn cip_ranks_by_atomic_number() {
        // bromochlorofluoromethane centre: Br > Cl > F > H
        let m = parse_smiles("[C@H](F)(Cl)Br").unwrap();
        let order = cip_neighbor_order(&m, 0);
        // first three are heavy neighbours, highest first
        let zs: Vec<u8> = order
            .iter()
            .filter(|&&v| v != usize::MAX)
            .map(|&v| m.atoms[v].atomic_number)
            .collect();
        assert_eq!(zs, vec![35, 17, 9], "Br,Cl,F order");
    }

    #[test]
    fn assigns_rs_label() {
        let mut m = parse_smiles("[C@H](F)(Cl)Br").unwrap();
        assign_cip_labels(&mut m);
        assert!(matches!(m.atoms[0].chirality, Chirality::R | Chirality::S));
        assert_eq!(stereocenter_count(&m), 1);

        // the mirror image must get the opposite label
        let mut m2 = parse_smiles("[C@@H](F)(Cl)Br").unwrap();
        assign_cip_labels(&mut m2);
        assert_ne!(m.atoms[0].chirality, m2.atoms[0].chirality);
    }

    #[test]
    fn non_stereocenter_left_alone() {
        // a carbon with two identical substituents is not a centre
        let mut m = parse_smiles("CC(C)O").unwrap();
        assign_cip_labels(&mut m);
        assert_eq!(stereocenter_count(&m), 0);
    }

    #[test]
    fn double_bond_ez() {
        // trans-2-butene F/C=C/F analogue with methyls
        let mut trans = parse_smiles("C/C=C/C").unwrap();
        assign_double_bond_stereo(&mut trans);
        let mut cis = parse_smiles("C/C=C\\C").unwrap();
        assign_double_bond_stereo(&mut cis);
        let t = trans.bonds.iter().find(|b| b.order == BondOrder::Double).unwrap().stereo;
        let c = cis.bonds.iter().find(|b| b.order == BondOrder::Double).unwrap().stereo;
        assert!(matches!(t, BondStereo::E | BondStereo::Z));
        assert!(matches!(c, BondStereo::E | BondStereo::Z));
        assert_ne!(t, c, "cis and trans must differ");
    }
}

