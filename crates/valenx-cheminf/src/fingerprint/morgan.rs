//! Morgan / ECFP circular fingerprints.
//!
//! Extended-Connectivity FingerPrints encode the circular environment
//! of every atom out to a given bond radius:
//!
//! 1. each atom starts with an identifier hashed from its invariants
//!    (atomic number, degree, hydrogen count, charge, ring flag — the
//!    "Daylight" atom invariants);
//! 2. for `radius` rounds, every atom's identifier is updated to a
//!    hash of `(round, own id, sorted (bond order, neighbour id)
//!    pairs)` — exactly the Morgan relaxation;
//! 3. each identifier seen at each round is a *feature*; the set of
//!    features is folded into a dense [`FingerprintBits`].
//!
//! `ecfp(mol, 2, 2048)` is the classic ECFP4 (radius 2 → diameter 4).

use super::bitvec::FingerprintBits;
use crate::molecule::Molecule;

/// Compute a Morgan / ECFP fingerprint.
///
/// - `radius` — number of relaxation rounds (ECFP*n* uses
///   `radius = n/2`; radius 2 is ECFP4);
/// - `n_bits` — length of the folded bit-vector.
pub fn ecfp(mol: &Molecule, radius: usize, n_bits: usize) -> FingerprintBits {
    let mut fp = FingerprintBits::new(n_bits);
    for feature in ecfp_features(mol, radius) {
        fp.set_hashed(feature);
    }
    fp
}

/// The raw (un-folded) set of ECFP feature identifiers — useful for
/// sparse / count fingerprints and for similarity that wants the exact
/// feature set rather than a folded bit-vector.
pub fn ecfp_features(mol: &Molecule, radius: usize) -> Vec<u64> {
    let n = mol.atoms.len();
    if n == 0 {
        return Vec::new();
    }
    // round-0 identifiers
    let mut ids: Vec<u64> = (0..n).map(|i| atom_identifier(mol, i)).collect();
    let mut features: Vec<u64> = ids.clone();

    for round in 1..=radius {
        let mut next = vec![0u64; n];
        for i in 0..n {
            // gather (bond order code, neighbour id) pairs, sorted
            let mut pairs: Vec<(u64, u64)> = mol
                .bonds_on(i)
                .iter()
                .filter_map(|&bi| {
                    let bond = &mol.bonds[bi];
                    bond.other(i).map(|v| (bond_code(bond), ids[v]))
                })
                .collect();
            pairs.sort_unstable();
            let mut h = fnv_init();
            fnv_feed(&mut h, round as u64);
            fnv_feed(&mut h, ids[i]);
            for (bc, nid) in pairs {
                fnv_feed(&mut h, bc);
                fnv_feed(&mut h, nid);
            }
            next[i] = h;
        }
        ids = next;
        features.extend_from_slice(&ids);
    }
    features.sort_unstable();
    features.dedup();
    features
}

/// FCFP-style: a functional-class circular fingerprint that seeds atoms
/// by *pharmacophoric role* (donor / acceptor / aromatic / halogen /
/// …) instead of element, giving a coarser, scaffold-tolerant
/// fingerprint.
pub fn fcfp(mol: &Molecule, radius: usize, n_bits: usize) -> FingerprintBits {
    let n = mol.atoms.len();
    let mut fp = FingerprintBits::new(n_bits);
    if n == 0 {
        return fp;
    }
    let mut ids: Vec<u64> = (0..n).map(|i| functional_class(mol, i)).collect();
    for f in &ids {
        fp.set_hashed(*f);
    }
    for round in 1..=radius {
        let mut next = vec![0u64; n];
        for i in 0..n {
            let mut pairs: Vec<(u64, u64)> = mol
                .bonds_on(i)
                .iter()
                .filter_map(|&bi| {
                    let bond = &mol.bonds[bi];
                    bond.other(i).map(|v| (bond_code(bond), ids[v]))
                })
                .collect();
            pairs.sort_unstable();
            let mut h = fnv_init();
            fnv_feed(&mut h, 0xF000 + round as u64);
            fnv_feed(&mut h, ids[i]);
            for (bc, nid) in pairs {
                fnv_feed(&mut h, bc);
                fnv_feed(&mut h, nid);
            }
            next[i] = h;
        }
        ids = next;
        for f in &ids {
            fp.set_hashed(*f);
        }
    }
    fp
}

/// The round-0 atom identifier — a hash of the Daylight invariants.
fn atom_identifier(mol: &Molecule, i: usize) -> u64 {
    let a = &mol.atoms[i];
    let mut h = fnv_init();
    fnv_feed(&mut h, u64::from(a.atomic_number));
    fnv_feed(&mut h, mol.degree(i) as u64);
    fnv_feed(&mut h, u64::from(a.total_h()));
    fnv_feed(&mut h, (i64::from(a.formal_charge) + 16) as u64);
    fnv_feed(&mut h, u64::from(a.aromatic));
    fnv_feed(&mut h, u64::from(a.isotope.unwrap_or(0)));
    h
}

/// A coarse pharmacophoric class id for FCFP seeding.
fn functional_class(mol: &Molecule, i: usize) -> u64 {
    let a = &mol.atoms[i];
    let mut class = 0u64;
    if a.aromatic {
        class |= 1;
    }
    // H-bond donor: N/O carrying at least one hydrogen
    if matches!(a.atomic_number, 7 | 8) && a.total_h() > 0 {
        class |= 2;
    }
    // H-bond acceptor: N/O (lone-pair bearing)
    if matches!(a.atomic_number, 7 | 8) {
        class |= 4;
    }
    // halogen
    if matches!(a.atomic_number, 9 | 17 | 35 | 53) {
        class |= 8;
    }
    // charged
    if a.formal_charge != 0 {
        class |= 16;
    }
    // ring / chain
    if mol.degree(i) >= 3 {
        class |= 32;
    }
    let mut h = fnv_init();
    fnv_feed(&mut h, class);
    h
}

/// Encode a bond's order into a small integer for hashing.
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

const FNV_OFFSET: u64 = 0xcbf29ce484222325;
const FNV_PRIME: u64 = 0x100000001b3;

fn fnv_init() -> u64 {
    FNV_OFFSET
}

fn fnv_feed(h: &mut u64, x: u64) {
    for shift in (0..64).step_by(8) {
        *h ^= (x >> shift) & 0xff;
        *h = h.wrapping_mul(FNV_PRIME);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn fingerprint_has_bits_set() {
        let m = mol_from_smiles("CCO").unwrap();
        let fp = ecfp(&m, 2, 1024);
        assert!(fp.count_ones() > 0);
        assert_eq!(fp.len(), 1024);
    }

    #[test]
    fn same_molecule_same_fingerprint() {
        let a = mol_from_smiles("CCO").unwrap();
        let b = mol_from_smiles("CCO").unwrap();
        assert_eq!(ecfp(&a, 2, 512), ecfp(&b, 2, 512));
    }

    #[test]
    fn different_molecules_differ() {
        let a = mol_from_smiles("CCO").unwrap();
        let b = mol_from_smiles("c1ccccc1").unwrap();
        assert_ne!(ecfp(&a, 2, 1024), ecfp(&b, 2, 1024));
    }

    #[test]
    fn larger_radius_more_features() {
        let m = mol_from_smiles("CCCCCCCC").unwrap();
        let f0 = ecfp_features(&m, 0).len();
        let f2 = ecfp_features(&m, 2).len();
        assert!(f2 >= f0, "radius 2 has at least as many features");
    }

    #[test]
    fn fcfp_works() {
        let m = mol_from_smiles("CCN").unwrap();
        let fp = fcfp(&m, 2, 512);
        assert!(fp.count_ones() > 0);
    }
}
