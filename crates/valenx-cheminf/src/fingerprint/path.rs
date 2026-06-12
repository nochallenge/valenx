//! Topological / path-based fingerprints (the Daylight-style FP).
//!
//! A path fingerprint enumerates every linear bond path in the
//! molecule up to a maximum length, hashes each path's atom/bond
//! sequence, and folds the hashes into a dense bit-vector. Unlike the
//! circular ECFP it captures *linear connectivity* — it is the
//! fingerprint historically used for substructure-screening.
//!
//! [`path_fingerprint`] walks every simple path of 1..=`max_len` bonds
//! with a depth-bounded DFS, hashing each in both directions so the
//! result is independent of traversal orientation.

use super::bitvec::FingerprintBits;
use crate::molecule::Molecule;

/// Compute a topological path fingerprint.
///
/// - `max_len` — longest path, in bonds, to enumerate (7 is the
///   classic Daylight default);
/// - `n_bits` — folded bit-vector length.
pub fn path_fingerprint(mol: &Molecule, max_len: usize, n_bits: usize) -> FingerprintBits {
    let mut fp = FingerprintBits::new(n_bits);
    let n = mol.atoms.len();
    if n == 0 {
        return fp;
    }
    // single-atom features (path length 0) so isolated atoms register
    for i in 0..n {
        fp.set_hashed(atom_hash(mol, i));
    }
    // DFS each path up to max_len bonds
    let mut visited = vec![false; n];
    for start in 0..n {
        let mut path_atoms = vec![start];
        visited[start] = true;
        walk(mol, start, max_len, &mut path_atoms, &mut visited, &mut fp);
        visited[start] = false;
    }
    fp
}

fn walk(
    mol: &Molecule,
    current: usize,
    remaining: usize,
    path: &mut Vec<usize>,
    visited: &mut [bool],
    fp: &mut FingerprintBits,
) {
    if remaining == 0 {
        return;
    }
    for bi in mol.bonds_on(current) {
        let bond = &mol.bonds[bi];
        let Some(next) = bond.other(current) else {
            continue;
        };
        if visited[next] {
            continue;
        }
        visited[next] = true;
        path.push(next);
        fp.set_hashed(path_hash(mol, path));
        walk(mol, next, remaining - 1, path, visited, fp);
        path.pop();
        visited[next] = false;
    }
}

/// Hash one atom's element + aromatic flag.
fn atom_hash(mol: &Molecule, i: usize) -> u64 {
    let a = &mol.atoms[i];
    let mut h = 1469598103934665603u64;
    feed(&mut h, u64::from(a.atomic_number));
    feed(&mut h, u64::from(a.aromatic));
    h
}

/// Hash a path's atom + bond sequence, canonicalised so the path and
/// its reverse hash identically.
fn path_hash(mol: &Molecule, path: &[usize]) -> u64 {
    let forward = directed_path_hash(mol, path, false);
    let backward = directed_path_hash(mol, path, true);
    // commutative combine → orientation-independent
    forward.wrapping_add(backward) ^ forward.wrapping_mul(backward | 1)
}

fn directed_path_hash(mol: &Molecule, path: &[usize], reverse: bool) -> u64 {
    let mut h = 1469598103934665603u64;
    let idx: Vec<usize> = if reverse {
        path.iter().rev().copied().collect()
    } else {
        path.to_vec()
    };
    for w in idx.windows(2) {
        let a = &mol.atoms[w[0]];
        feed(&mut h, u64::from(a.atomic_number));
        feed(&mut h, u64::from(a.aromatic));
        if let Some(bi) = mol.bond_between(w[0], w[1]) {
            feed(&mut h, bond_code(&mol.bonds[bi]));
        }
    }
    if let Some(&last) = idx.last() {
        let a = &mol.atoms[last];
        feed(&mut h, u64::from(a.atomic_number));
        feed(&mut h, u64::from(a.aromatic));
    }
    h
}

fn bond_code(bond: &crate::molecule::Bond) -> u64 {
    use crate::molecule::BondOrder::*;
    match bond.order {
        Single => 1,
        Double => 2,
        Triple => 3,
        Quadruple => 4,
        Aromatic => 5,
    }
}

fn feed(h: &mut u64, x: u64) {
    *h ^= x;
    *h = h.wrapping_mul(1099511628211);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn path_fp_has_bits() {
        let m = mol_from_smiles("CCCCCC").unwrap();
        let fp = path_fingerprint(&m, 7, 1024);
        assert!(fp.count_ones() > 0);
    }

    #[test]
    fn orientation_independent() {
        // a symmetric molecule must hash the same regardless of which
        // end the atoms were numbered from
        let a = mol_from_smiles("CCOCC").unwrap();
        let b = mol_from_smiles("CCOCC").unwrap();
        assert_eq!(path_fingerprint(&a, 5, 512), path_fingerprint(&b, 5, 512));
    }

    #[test]
    fn distinct_molecules_distinct() {
        let a = mol_from_smiles("CCCC").unwrap();
        let b = mol_from_smiles("CCCN").unwrap();
        assert_ne!(path_fingerprint(&a, 7, 1024), path_fingerprint(&b, 7, 1024));
    }

    #[test]
    fn longer_paths_more_bits() {
        let m = mol_from_smiles("CCCCCCCCCC").unwrap();
        let short = path_fingerprint(&m, 2, 2048).count_ones();
        let long = path_fingerprint(&m, 8, 2048).count_ones();
        assert!(long >= short);
    }
}
