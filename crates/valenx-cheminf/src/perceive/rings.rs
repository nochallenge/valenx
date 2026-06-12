//! Ring perception — the smallest set of smallest rings (SSSR).
//!
//! The SSSR is the chemist's working ring set: a set of `B − A + C`
//! rings (B = bonds, A = atoms, C = connected components — the
//! circuit rank) chosen so each is as small as possible. It is not a
//! unique invariant for fused / spiro systems (the "SSSR ambiguity"),
//! but it is the set RDKit, Open Babel and every depiction tool use,
//! and it is exactly right for aromaticity perception and ring-count
//! descriptors.
//!
//! Algorithm: for each bond, run a BFS that finds the shortest cycle
//! through that bond; collect candidate rings sorted by size; then
//! greedily accept a ring iff it adds a new linearly-independent row
//! to the cycle space over GF(2) (Gaussian elimination on bond-set bit
//! vectors). This is the standard Figueras / "smallest-rings" approach
//! and yields a correct SSSR for the molecules this crate targets.

use crate::molecule::Molecule;
use std::collections::VecDeque;

/// One ring: the ordered list of atom indices around the cycle, and the
/// set of bond indices it uses.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Ring {
    /// Atom indices in cycle order (first atom is *not* repeated).
    pub atoms: Vec<usize>,
    /// Bond indices forming the ring, sorted ascending.
    pub bonds: Vec<usize>,
}

impl Ring {
    /// Ring size — the number of atoms (== number of bonds).
    pub fn size(&self) -> usize {
        self.atoms.len()
    }

    /// `true` if `atom` lies on this ring.
    pub fn contains_atom(&self, atom: usize) -> bool {
        self.atoms.contains(&atom)
    }

    /// `true` if `bond` lies on this ring.
    pub fn contains_bond(&self, bond: usize) -> bool {
        self.bonds.binary_search(&bond).is_ok()
    }
}

/// The full ring analysis of a molecule.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RingInfo {
    /// The SSSR rings.
    pub rings: Vec<Ring>,
    /// Per-atom: how many SSSR rings the atom belongs to.
    pub atom_ring_count: Vec<u32>,
    /// Per-bond: how many SSSR rings the bond belongs to.
    pub bond_ring_count: Vec<u32>,
}

impl RingInfo {
    /// `true` if `atom` is in any ring.
    pub fn atom_in_ring(&self, atom: usize) -> bool {
        self.atom_ring_count.get(atom).copied().unwrap_or(0) > 0
    }

    /// `true` if `bond` is in any ring.
    pub fn bond_in_ring(&self, bond: usize) -> bool {
        self.bond_ring_count.get(bond).copied().unwrap_or(0) > 0
    }

    /// Number of SSSR rings.
    pub fn ring_count(&self) -> usize {
        self.rings.len()
    }

    /// `true` if `atom` belongs to two or more rings (a fusion atom).
    pub fn is_fusion_atom(&self, atom: usize) -> bool {
        self.atom_ring_count.get(atom).copied().unwrap_or(0) >= 2
    }
}

/// Circuit rank `B − A + C` — the number of independent rings the SSSR
/// must contain.
pub fn circuit_rank(mol: &Molecule) -> usize {
    let b = mol.bond_count();
    let a = mol.atom_count();
    let c = mol.component_count();
    (b + c).saturating_sub(a)
}

/// Compute the smallest set of smallest rings for `mol`.
pub fn sssr(mol: &Molecule) -> RingInfo {
    let n_atoms = mol.atoms.len();
    let n_bonds = mol.bonds.len();
    let mut info = RingInfo {
        rings: Vec::new(),
        atom_ring_count: vec![0; n_atoms],
        bond_ring_count: vec![0; n_bonds],
    };
    let target = circuit_rank(mol);
    if target == 0 {
        return info;
    }

    // Candidate rings: for each bond, the shortest cycle through it.
    let mut candidates: Vec<Ring> = Vec::new();
    for bi in 0..n_bonds {
        if let Some(ring) = shortest_cycle_through_bond(mol, bi) {
            candidates.push(ring);
        }
    }
    // Smallest first; deterministic on ties via the bond-set.
    candidates.sort_by(|a, b| {
        a.atoms
            .len()
            .cmp(&b.atoms.len())
            .then_with(|| a.bonds.cmp(&b.bonds))
    });
    candidates.dedup_by(|a, b| a.bonds == b.bonds);

    // Greedy GF(2) independence selection.
    let mut basis: Vec<Vec<u64>> = Vec::new();
    let words = n_bonds.div_ceil(64);
    for cand in candidates {
        if info.rings.len() == target {
            break;
        }
        let mut vec = vec![0u64; words];
        for &b in &cand.bonds {
            vec[b / 64] |= 1u64 << (b % 64);
        }
        if is_independent(&basis, &vec) {
            basis.push(vec);
            info.rings.push(cand);
        }
    }

    for ring in &info.rings {
        for &a in &ring.atoms {
            info.atom_ring_count[a] += 1;
        }
        for &b in &ring.bonds {
            info.bond_ring_count[b] += 1;
        }
    }
    info
}

/// BFS shortest cycle that uses bond `bi`: drop `bi`, find the shortest
/// path between its endpoints, then re-add `bi` to close the cycle.
fn shortest_cycle_through_bond(mol: &Molecule, bi: usize) -> Option<Ring> {
    let (start, goal) = (mol.bonds[bi].a, mol.bonds[bi].b);
    // BFS from `start` to `goal` not using bond `bi`.
    let n = mol.atoms.len();
    let mut prev_atom = vec![usize::MAX; n];
    let mut prev_bond = vec![usize::MAX; n];
    let mut visited = vec![false; n];
    let mut q = VecDeque::new();
    q.push_back(start);
    visited[start] = true;
    while let Some(u) = q.pop_front() {
        if u == goal {
            break;
        }
        for bj in mol.bonds_on(u) {
            if bj == bi {
                continue;
            }
            if let Some(v) = mol.bonds[bj].other(u) {
                if !visited[v] {
                    visited[v] = true;
                    prev_atom[v] = u;
                    prev_bond[v] = bj;
                    q.push_back(v);
                }
            }
        }
    }
    if !visited[goal] {
        return None;
    }
    // reconstruct path goal -> start
    let mut atoms = vec![goal];
    let mut bonds = vec![bi];
    let mut cur = goal;
    while cur != start {
        bonds.push(prev_bond[cur]);
        cur = prev_atom[cur];
        atoms.push(cur);
    }
    // atoms currently goal..start; that is a valid cycle order
    bonds.sort_unstable();
    bonds.dedup();
    Some(Ring { atoms, bonds })
}

/// GF(2) linear-independence test: can `vec` be reduced to zero by the
/// existing `basis`? Returns `true` if it is independent (non-zero
/// after reduction).
fn is_independent(basis: &[Vec<u64>], vec: &[u64]) -> bool {
    let mut v = vec.to_vec();
    for row in basis {
        // pivot = lowest set bit of `row`
        if let Some(pivot) = lowest_set_bit(row) {
            let (w, b) = (pivot / 64, pivot % 64);
            if v[w] & (1u64 << b) != 0 {
                for i in 0..v.len() {
                    v[i] ^= row[i];
                }
            }
        }
    }
    v.iter().any(|&w| w != 0)
}

fn lowest_set_bit(words: &[u64]) -> Option<usize> {
    for (i, &w) in words.iter().enumerate() {
        if w != 0 {
            return Some(i * 64 + w.trailing_zeros() as usize);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn benzene_one_ring() {
        let m = parse_smiles("c1ccccc1").unwrap();
        let info = sssr(&m);
        assert_eq!(info.ring_count(), 1);
        assert_eq!(info.rings[0].size(), 6);
        assert!(info.atom_in_ring(0));
    }

    #[test]
    fn acyclic_no_rings() {
        let m = parse_smiles("CCCCC").unwrap();
        let info = sssr(&m);
        assert_eq!(info.ring_count(), 0);
        assert_eq!(circuit_rank(&m), 0);
    }

    #[test]
    fn naphthalene_two_rings() {
        let m = parse_smiles("c1ccc2ccccc2c1").unwrap();
        let info = sssr(&m);
        assert_eq!(info.ring_count(), 2);
        // both rings size 6
        assert!(info.rings.iter().all(|r| r.size() == 6));
        // the two fusion atoms belong to both rings
        let fused = (0..m.atom_count())
            .filter(|&a| info.is_fusion_atom(a))
            .count();
        assert_eq!(fused, 2);
    }

    #[test]
    fn cyclohexane_ring_size() {
        let m = parse_smiles("C1CCCCC1").unwrap();
        let info = sssr(&m);
        assert_eq!(info.ring_count(), 1);
        assert_eq!(info.rings[0].size(), 6);
    }

    #[test]
    fn bicyclic_circuit_rank() {
        // norbornane-ish: two fused rings sharing a bridge
        let m = parse_smiles("C1CC2CCC1C2").unwrap();
        let info = sssr(&m);
        assert_eq!(circuit_rank(&m), 2);
        assert_eq!(info.ring_count(), 2);
    }
}
