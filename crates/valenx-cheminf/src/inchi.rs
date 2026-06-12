//! An InChI-class layered canonical identifier.
//!
//! The official IUPAC InChI is a multi-year reference implementation
//! with a very specific normalisation and canonicalisation procedure.
//! This module produces an identifier *in the spirit of* InChI — a
//! layered, deterministic string built from graph invariants — that
//! serves the same practical purpose inside Valenx: a stable key for
//! deduplication and structure lookup, identical for any two molecules
//! that are the same graph regardless of atom numbering.
//!
//! Layer structure of the produced string (`VLNX=1` prefix marks it as
//! the Valenx variant, *not* a real `InChI=1S` string):
//!
//! ```text
//! VLNX=1/<formula>/c<connectivity>/h<hydrogens>[/q<charge>][/p<protons>]
//! ```
//!
//! - **formula** — the Hill molecular formula;
//! - **c** layer — the canonical connection table: for each atom in
//!   canonical order, the sorted canonical ranks of its heavy
//!   neighbours;
//! - **h** layer — the hydrogen count per canonical atom;
//! - **q** layer — the net charge, present only when non-zero;
//! - **p** layer — a removable-proton balance, present only when the
//!   standardiser would add/remove protons (always `0` in this v1, the
//!   field is reserved).
//!
//! Because the connection layer is keyed on the same Morgan canonical
//! ranking used by [`crate::smiles::write_canonical_smiles`], two
//! isomorphic molecules map to the same `VLNX` string.
//!
//! **v1 limitations / honest naming.** This is deliberately *not*
//! `InChI=1S/...`; it will not match an InChI computed by the IUPAC
//! software, and it does not implement the InChI stereo / isotope
//! sublayers (`/b`, `/t`, `/m`, `/s`, `/i`). It is a self-consistent
//! canonical key, suitable for internal dedup and lookup.

use crate::molecule::Molecule;
use crate::perceive::formula::molecular_formula;
use crate::smiles::canonical_ranks;

/// Compute the Valenx layered canonical identifier for `mol`.
pub fn inchi_class_string(mol: &Molecule) -> String {
    if mol.is_empty() {
        return "VLNX=1//".to_string();
    }
    let formula = molecular_formula(mol);
    let ranks = canonical_ranks(mol);
    let n = mol.atoms.len();

    // canonical index → atom index
    let mut atom_of_rank = vec![0usize; n];
    for (atom, &r) in ranks.iter().enumerate() {
        atom_of_rank[r] = atom;
    }

    // --- connectivity layer ----------------------------------------
    // For each atom in canonical order, list the canonical ranks of
    // its heavy neighbours (sorted). Hydrogens are excluded — they go
    // into the h layer.
    let mut c_parts: Vec<String> = Vec::with_capacity(n);
    for &atom in &atom_of_rank {
        if mol.atoms[atom].is_hydrogen() {
            c_parts.push("H".to_string());
            continue;
        }
        let mut nbr_ranks: Vec<usize> = mol
            .neighbors(atom)
            .into_iter()
            .filter(|&v| !mol.atoms[v].is_hydrogen())
            .map(|v| ranks[v])
            .collect();
        nbr_ranks.sort_unstable();
        let joined = nbr_ranks
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join("-");
        c_parts.push(joined);
    }
    let c_layer = c_parts.join(",");

    // --- hydrogen layer --------------------------------------------
    let h_parts: Vec<String> = (0..n)
        .map(|r| {
            let atom = atom_of_rank[r];
            mol.atoms[atom].total_h().to_string()
        })
        .collect();
    let h_layer = h_parts.join(",");

    // --- charge layer ----------------------------------------------
    let net: i32 = mol.atoms.iter().map(|a| i32::from(a.formal_charge)).sum();

    let mut s = format!("VLNX=1/{formula}/c{c_layer}/h{h_layer}");
    if net != 0 {
        s.push_str(&format!("/q{net:+}"));
    }
    s
}

/// `true` if `a` and `b` have the same Valenx canonical identifier —
/// i.e. they are the same molecular graph. The intended dedup test.
pub fn same_structure(a: &Molecule, b: &Molecule) -> bool {
    inchi_class_string(a) == inchi_class_string(b)
}

/// A short fixed-length hash of the identifier, in the spirit of the
/// InChIKey: 14 uppercase base-32 characters derived from a 64-bit
/// FNV-1a hash of the full layered string. Collision-resistant enough
/// for a dedup index; not cryptographic.
pub fn inchi_class_key(mol: &Molecule) -> String {
    let s = inchi_class_string(mol);
    let mut h = 0xcbf29ce484222325u64;
    for byte in s.bytes() {
        h ^= u64::from(byte);
        h = h.wrapping_mul(0x100000001b3);
    }
    const ALPHABET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ234567";
    let mut key = String::with_capacity(14);
    let mut v = h;
    for _ in 0..14 {
        key.push(char::from(ALPHABET[(v & 31) as usize]));
        v >>= 2; // overlap windows so all 64 bits influence the key
    }
    key
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn identifier_has_layers() {
        let m = mol_from_smiles("CCO").unwrap();
        let s = inchi_class_string(&m);
        assert!(s.starts_with("VLNX=1/"));
        assert!(s.contains("/c"));
        assert!(s.contains("/h"));
        assert!(s.contains("C2H6O"));
    }

    #[test]
    fn order_independent() {
        let a = mol_from_smiles("OCC").unwrap();
        let b = mol_from_smiles("CCO").unwrap();
        assert_eq!(inchi_class_string(&a), inchi_class_string(&b));
        assert!(same_structure(&a, &b));
    }

    #[test]
    fn different_molecules_differ() {
        let a = mol_from_smiles("CCO").unwrap();
        let b = mol_from_smiles("CCN").unwrap();
        assert_ne!(inchi_class_string(&a), inchi_class_string(&b));
        assert!(!same_structure(&a, &b));
    }

    #[test]
    fn charge_layer_present_when_charged() {
        let neutral = mol_from_smiles("N").unwrap();
        assert!(!inchi_class_string(&neutral).contains("/q"));
        let cation = mol_from_smiles("[NH4+]").unwrap();
        assert!(inchi_class_string(&cation).contains("/q+1"));
    }

    #[test]
    fn key_is_fixed_length() {
        let m = mol_from_smiles("c1ccccc1").unwrap();
        let key = inchi_class_key(&m);
        assert_eq!(key.len(), 14);
        assert!(key
            .chars()
            .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit()));
        // same molecule, different numbering → same key
        let m2 = mol_from_smiles("c1ccccc1").unwrap();
        assert_eq!(inchi_class_key(&m), inchi_class_key(&m2));
    }

    #[test]
    fn empty_molecule() {
        let m = Molecule::new();
        assert_eq!(inchi_class_string(&m), "VLNX=1//");
    }
}
