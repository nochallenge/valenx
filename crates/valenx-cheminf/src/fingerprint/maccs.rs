//! MACCS-class structural keys.
//!
//! The MDL MACCS set is a fixed dictionary of substructure questions —
//! "is there a 6-membered ring?", "is there an `O=C`?" — each setting
//! one bit. The official set is 166 keys; this module ships a
//! ~80-question subset expressed as SMARTS patterns plus a handful of
//! compositional / counting keys. It is enough to give a meaningful
//! structural-key fingerprint for similarity and clustering.
//!
//! [`maccs_keys`] returns a [`FingerprintBits`] one bit per defined
//! key. The set is deliberately documented as *MACCS-class*, not the
//! exact MDL definition (some original keys depend on properties this
//! v1 does not perceive — those slots are simply omitted).

use super::bitvec::FingerprintBits;
use crate::molecule::Molecule;
use crate::perceive::rings::sssr;
use crate::smarts::SmartsPattern;

/// Each MACCS-class key: a label and the SMARTS that answers it.
struct KeyDef {
    smarts: &'static str,
}

/// The SMARTS-defined structural keys. Order is fixed — bit `i`
/// corresponds to `KEY_DEFS[i]`.
const KEY_DEFS: &[KeyDef] = &[
    KeyDef { smarts: "[#8]" },  // oxygen present
    KeyDef { smarts: "[#7]" },  // nitrogen present
    KeyDef { smarts: "[#16]" }, // sulfur present
    KeyDef { smarts: "[#15]" }, // phosphorus present
    KeyDef {
        smarts: "[F,Cl,Br,I]",
    }, // halogen present
    KeyDef {
        smarts: "[#6]=[#8]",
    }, // carbonyl
    KeyDef {
        smarts: "[#6]=[#7]",
    }, // imine
    KeyDef {
        smarts: "[#6]#[#7]",
    }, // nitrile
    KeyDef {
        smarts: "[#6]#[#6]",
    }, // alkyne
    KeyDef {
        smarts: "[#6]=[#6]",
    }, // alkene
    KeyDef { smarts: "[#8H1]" }, // hydroxyl-ish (O with H)
    KeyDef {
        smarts: "[#7H1,#7H2]",
    }, // N-H
    KeyDef {
        smarts: "[#6][#8H1]",
    }, // C-O-H alcohol
    KeyDef {
        smarts: "[#6](=[#8])[#8]",
    }, // carboxyl / ester O
    KeyDef {
        smarts: "[#6](=[#8])[#7]",
    }, // amide
    KeyDef {
        smarts: "[#6](=[#8])[#6]",
    }, // ketone
    KeyDef {
        smarts: "[#7](=[#8])",
    }, // nitro-ish N=O
    KeyDef {
        smarts: "[#8]=[#16]",
    }, // S=O sulfoxide/sulfone
    KeyDef { smarts: "[#6][#7]" }, // C-N
    KeyDef { smarts: "[#6][#8]" }, // C-O
    KeyDef {
        smarts: "[#6][#16]",
    }, // C-S
    KeyDef { smarts: "[#7][#7]" }, // N-N
    KeyDef { smarts: "[#8][#8]" }, // O-O peroxide
    KeyDef { smarts: "[#7][#8]" }, // N-O
    KeyDef { smarts: "c" },     // aromatic atom
    KeyDef { smarts: "c[#7]" }, // aromatic C-N
    KeyDef { smarts: "c[#8]" }, // aromatic C-O
    KeyDef {
        smarts: "c[F,Cl,Br,I]",
    }, // aromatic halide
    KeyDef {
        smarts: "[#6][#6][#6][#6]",
    }, // 4-carbon chain
    KeyDef {
        smarts: "[#6][#6][#6][#6][#6][#6]",
    }, // 6-carbon chain
    KeyDef { smarts: "[CH3]" }, // methyl
    KeyDef { smarts: "[CH2]" }, // methylene
    KeyDef { smarts: "[CH]" },  // methine
    KeyDef { smarts: "[#6X4]" }, // sp3 carbon
    KeyDef { smarts: "[#7X3]" }, // trivalent nitrogen
    KeyDef { smarts: "[#8X2]" }, // divalent oxygen (ether)
    KeyDef { smarts: "[+]" },   // a cation
    KeyDef { smarts: "[-]" },   // an anion
    KeyDef { smarts: "[R]" },   // any ring atom
    KeyDef { smarts: "[r3]" },  // 3-membered ring
    KeyDef { smarts: "[r4]" },  // 4-membered ring
    KeyDef { smarts: "[r5]" },  // 5-membered ring
    KeyDef { smarts: "[r6]" },  // 6-membered ring
    KeyDef { smarts: "[r7]" },  // 7-membered ring
    KeyDef { smarts: "[#7R]" }, // nitrogen in a ring
    KeyDef { smarts: "[#8R]" }, // oxygen in a ring
    KeyDef { smarts: "[#16R]" }, // sulfur in a ring
    KeyDef { smarts: "c1ccccc1" }, // benzene ring
    KeyDef {
        smarts: "[#6][#6](=[#8])[#8H1]",
    }, // carboxylic acid
    KeyDef {
        smarts: "[#7][#6]=[#8]",
    }, // amide N-C=O
    KeyDef {
        smarts: "[#6][#7][#6]",
    }, // C-N-C
    KeyDef {
        smarts: "[#6][#8][#6]",
    }, // ether C-O-C
    KeyDef {
        smarts: "[#6][#16][#6]",
    }, // thioether
    KeyDef {
        smarts: "[#8H1][#6]=[#8]",
    }, // (another carboxyl form)
    KeyDef { smarts: "[#7H2]" }, // primary amine NH2
    KeyDef {
        smarts: "[#6]=[#6][#6]=[#6]",
    }, // conjugated diene
    KeyDef {
        smarts: "[#6R]=[#6R]",
    }, // ring double bond
    KeyDef {
        smarts: "[#6](-[#6])(-[#6])(-[#6])-[#6]",
    }, // quaternary C
    KeyDef { smarts: "[#9]" },  // fluorine
    KeyDef { smarts: "[#17]" }, // chlorine
    KeyDef { smarts: "[#35]" }, // bromine
    KeyDef { smarts: "[#53]" }, // iodine
];

/// Compute the MACCS-class structural-key fingerprint.
///
/// The fingerprint length is `KEY_DEFS.len() + 4`: one bit per SMARTS
/// key plus four compositional count keys at the end (≥1 ring, ≥2
/// rings, ≥8 heavy atoms, ≥16 heavy atoms).
pub fn maccs_keys(mol: &Molecule) -> FingerprintBits {
    let n_struct = KEY_DEFS.len();
    let mut fp = FingerprintBits::new(n_struct + 4);
    for (i, def) in KEY_DEFS.iter().enumerate() {
        if let Ok(pat) = SmartsPattern::parse(def.smarts) {
            if pat.matches(mol) {
                fp.set(i);
            }
        }
    }
    // compositional count keys
    let rings = sssr(mol).ring_count();
    let heavy = mol.heavy_atom_count();
    if rings >= 1 {
        fp.set(n_struct);
    }
    if rings >= 2 {
        fp.set(n_struct + 1);
    }
    if heavy >= 8 {
        fp.set(n_struct + 2);
    }
    if heavy >= 16 {
        fp.set(n_struct + 3);
    }
    fp
}

/// Number of MACCS-class keys (the fingerprint length).
pub fn maccs_key_count() -> usize {
    KEY_DEFS.len() + 4
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn benzene_sets_aromatic_keys() {
        let m = mol_from_smiles("c1ccccc1").unwrap();
        let fp = maccs_keys(&m);
        assert!(fp.count_ones() > 0);
        assert_eq!(fp.len(), maccs_key_count());
        // aromatic-atom key (index 24) and benzene-ring key (index 47)
        assert!(fp.get(24), "aromatic atom key");
        assert!(fp.get(47), "benzene ring key");
    }

    #[test]
    fn carboxylic_acid_keys() {
        let m = mol_from_smiles("CC(=O)O").unwrap();
        let fp = maccs_keys(&m);
        // carbonyl key (index 5) and oxygen key (index 0)
        assert!(fp.get(0), "oxygen present");
        assert!(fp.get(5), "carbonyl");
    }

    #[test]
    fn alkane_no_heteroatom_keys() {
        let m = mol_from_smiles("CCCCCC").unwrap();
        let fp = maccs_keys(&m);
        assert!(!fp.get(0), "no oxygen");
        assert!(!fp.get(1), "no nitrogen");
        assert!(!fp.get(24), "no aromatic atom");
    }

    #[test]
    fn distinct_molecules_distinct_keys() {
        let a = mol_from_smiles("c1ccccc1").unwrap();
        let b = mol_from_smiles("CCCCCC").unwrap();
        assert_ne!(maccs_keys(&a), maccs_keys(&b));
    }
}
