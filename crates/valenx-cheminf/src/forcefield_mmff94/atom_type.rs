//! MMFF94 atom typing.
//!
//! Implements the perception of the [Merck Molecular Force Field][halgren]
//! atom types (Halgren, J. Comp. Chem. 17 (1996), part II, "MMFF94 atom
//! types"). Each heavy atom is assigned a [`MmffType`] integer code in
//! the published MMFF94 numbering (e.g. `1` = `CR` sp³ carbon, `37` =
//! `CB` aromatic carbon, `7` = `O=C` carbonyl oxygen) and a parallel
//! ASCII symbol (`"CR"`, `"CB"`, `"O=C"`) for diagnostics. Hydrogens get
//! their MMFF94 type from the heavy atom they're attached to (`5` =
//! `HC` carbon-attached H, `21` = `HOR` alcohol H, …).
//!
//! [halgren]: https://doi.org/10.1002/(SICI)1096-987X(199604)17:5/6%3C490::AID-JCC1%3E3.0.CO;2-P
//!
//! ## Coverage
//!
//! A *representative subset* of MMFF94's ~95 published types covering
//! the common organic-chemistry subspace this crate targets — C / H / N
//! / O / S / P / halogens. The full list of types this perceiver can
//! produce (numbered as in Halgren table I):
//!
//! | # | symbol | meaning |
//! |---|---|---|
//! | 1 | `CR` | sp³ alkyl carbon |
//! | 2 | `C=C` | sp² olefinic carbon |
//! | 3 | `C=O`, `C=N` carbonyl / imine carbon (also amide, ester, carboxylate) |
//! | 4 | `CSP` | sp acetylenic carbon |
//! | 6 | `OR` | divalent ether / alcohol / phenol oxygen |
//! | 7 | `O=C` | carbonyl / amide / carboxyl =O |
//! | 8 | `NR` | sp³ amine nitrogen |
//! | 9 | `N=C` | sp² imine / amidine nitrogen |
//! | 10 | `NC=O` | amide / sulfonamide nitrogen |
//! | 11 | `F` | fluorine |
//! | 12 | `CL` | chlorine |
//! | 13 | `BR` | bromine |
//! | 14 | `I` | iodine |
//! | 15 | `S` | divalent thioether / thiol sulfur |
//! | 16 | `=S` | doubly-bonded sulfur (thione) |
//! | 17 | `S=O` | sulfoxide sulfur |
//! | 18 | `SO2` | sulfone / sulfonate sulfur |
//! | 25 | `P` | phosphine / phosphate P |
//! | 26 | `-P=C` | phosphorus in `-P=` system |
//! | 32 | `O2CM` | carboxylate / nitro / sulfate terminal O |
//! | 35 | `OM` | oxide / alkoxide O⁻ |
//! | 37 | `CB` | aromatic carbon in 6-ring (benzene) |
//! | 38 | `NPYD` | pyridine-type aromatic N |
//! | 39 | `NPYL` | pyrrole-type aromatic N |
//! | 40 | `NC=C` | enamine / vinylogous amide N |
//! | 41 | `CO2M` | carboxylate / sulfonate carbon |
//! | 5 | `HC` | H attached to carbon |
//! | 21 | `HOR` | H on alcohol / phenol O |
//! | 23 | `HNR` | H on amine / amide N |
//! | 27 | `HOCO` | H on carboxylic-acid O |
//! | 28 | `HOCC` | H on enol / phenol O attached to sp² C |
//! | 29 | `HOS` | H on sulfonic-acid / phosphate O |
//! | 71 | `HS` | H on sulfur |
//!
//! Atoms outside this subset (uncommon heteroatom states, the rare
//! types for boron and metals, the inorganic anion shells) fall back to
//! [`MmffType::UNKNOWN`] — a documented gap that callers may detect.
//!
//! Atom typing is invoked by [`type_molecule`] which returns a vector
//! parallel to `mol.atoms`. It needs aromaticity + ring perception to
//! have run (see [`crate::perceive::perceive_all`]).

use crate::molecule::{BondOrder, Molecule};
use crate::perceive::rings::{sssr, RingInfo};

/// A perceived MMFF94 atom type with its published integer index and
/// short string symbol.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct MmffType {
    /// MMFF94 numeric type (Halgren table I).
    pub number: u8,
    /// Short ASCII symbol (`"CR"`, `"O=C"`).
    pub symbol: &'static str,
}

impl MmffType {
    /// Sentinel for an atom outside the implemented subset. Carries
    /// number `0`.
    pub const UNKNOWN: MmffType = MmffType {
        number: 0,
        symbol: "??",
    };

    /// `true` if this is the [`MmffType::UNKNOWN`] sentinel.
    pub fn is_unknown(self) -> bool {
        self.number == 0
    }
}

// --- the MMFF94 type table (subset implemented here) ----------------

macro_rules! mt {
    ($n:expr, $s:expr) => {
        MmffType {
            number: $n,
            symbol: $s,
        }
    };
}

// Carbons
const CR: MmffType = mt!(1, "CR"); // sp3 alkyl
const C_E_C: MmffType = mt!(2, "C=C"); // sp2 olefin
const C_E_O: MmffType = mt!(3, "C=O"); // carbonyl / amide / ester / imine sp2 C
const CSP: MmffType = mt!(4, "CSP"); // sp acetylene
const CB: MmffType = mt!(37, "CB"); // aromatic 6-ring C (benzene)
const CO2M: MmffType = mt!(41, "CO2M"); // carboxylate / sulfonate central C

// Hydrogens
const HC: MmffType = mt!(5, "HC"); // H on carbon
const HOR: MmffType = mt!(21, "HOR"); // alcohol / ether H
const HNR: MmffType = mt!(23, "HNR"); // amine / amide H
const HOCO: MmffType = mt!(27, "HOCO"); // carboxylic-acid H
const HOCC: MmffType = mt!(28, "HOCC"); // enol / phenol H
const HOS: MmffType = mt!(29, "HOS"); // sulfonic / phosphate H
const HS: MmffType = mt!(71, "HS"); // H on sulfur

// Oxygens
const OR: MmffType = mt!(6, "OR"); // ether / alcohol / phenol divalent
const O_E_C: MmffType = mt!(7, "O=C"); // carbonyl / amide / carboxyl =O
const O2CM: MmffType = mt!(32, "O2CM"); // carboxylate / nitro / sulfonate terminal
const OM: MmffType = mt!(35, "OM"); // alkoxide / oxide O-

// Nitrogens
const NR: MmffType = mt!(8, "NR"); // sp3 amine
const N_E_C: MmffType = mt!(9, "N=C"); // sp2 imine
const NC_E_O: MmffType = mt!(10, "NC=O"); // amide / sulfonamide
const NPYD: MmffType = mt!(38, "NPYD"); // pyridine-type aromatic N
const NPYL: MmffType = mt!(39, "NPYL"); // pyrrole-type aromatic N
const NC_E_C: MmffType = mt!(40, "NC=C"); // enamine / vinylogous amide N

// Sulfur
const SR: MmffType = mt!(15, "S"); // divalent thiother / thiol
const S_E_C: MmffType = mt!(16, "=S"); // thione
const SO: MmffType = mt!(17, "S=O"); // sulfoxide
const SO2: MmffType = mt!(18, "SO2"); // sulfone / sulfonate
                                      // Phosphorus
const PR: MmffType = mt!(25, "P"); // phosphine / phosphate
const PEC: MmffType = mt!(26, "-P=C"); // P in -P= system
                                       // Halogens
const F_TYPE: MmffType = mt!(11, "F");
const CL_TYPE: MmffType = mt!(12, "CL");
const BR_TYPE: MmffType = mt!(13, "BR");
const I_TYPE: MmffType = mt!(14, "I");

/// Run MMFF94 atom typing on `mol`. The result is parallel to
/// `mol.atoms`; an atom outside the implemented subset gets
/// [`MmffType::UNKNOWN`].
///
/// The molecule's aromaticity must already be perceived (call
/// [`crate::perceive::perceive_all`] first). Hydrogens are typed from
/// their heavy neighbour, so call after [`crate::perceive::hydrogen::add_explicit_hydrogens`]
/// if you want H types populated.
pub fn type_molecule(mol: &Molecule) -> Vec<MmffType> {
    let rings = sssr(mol);
    let mut types = vec![MmffType::UNKNOWN; mol.atoms.len()];
    // First pass: heavy atoms.
    for (i, ty) in types.iter_mut().enumerate() {
        if mol.atoms[i].atomic_number == 1 {
            continue;
        }
        *ty = perceive_heavy(mol, i, &rings);
    }
    // Second pass: hydrogens, by their heavy neighbour's type. Read
    // from the (now-filled) `types` slice via clone-and-overwrite.
    let snapshot = types.clone();
    for (i, ty) in types.iter_mut().enumerate() {
        if mol.atoms[i].atomic_number != 1 {
            continue;
        }
        *ty = perceive_hydrogen(mol, i, &snapshot);
    }
    types
}

fn perceive_heavy(mol: &Molecule, i: usize, rings: &RingInfo) -> MmffType {
    let a = &mol.atoms[i];
    match a.atomic_number {
        6 => type_carbon(mol, i, rings),
        7 => type_nitrogen(mol, i, rings),
        8 => type_oxygen(mol, i),
        15 => type_phosphorus(mol, i),
        16 => type_sulfur(mol, i),
        9 => F_TYPE,
        17 => CL_TYPE,
        35 => BR_TYPE,
        53 => I_TYPE,
        _ => MmffType::UNKNOWN,
    }
}

// --- Carbon -----------------------------------------------------------

fn type_carbon(mol: &Molecule, i: usize, _rings: &RingInfo) -> MmffType {
    // Aromatic six-membered ring → CB (37). Also covers fused aromatic
    // rings (naphthalene, indole etc.) — every aromatic ring atom maps
    // to CB in this subset because we don't separately enumerate the
    // five-ring carbon variants.
    if mol.atoms[i].aromatic {
        return CB;
    }
    let bonds = mol.bonds_on(i);
    // Pick the highest-order bond from this carbon to drive sp/sp²/sp³
    // classification.
    let has_triple = bonds
        .iter()
        .any(|&b| mol.bonds[b].order == BondOrder::Triple);
    if has_triple {
        return CSP;
    }
    let has_double = bonds
        .iter()
        .any(|&b| mol.bonds[b].order == BondOrder::Double);
    if has_double {
        // Carbonyl / imine / olefin discrimination.
        let mut bonded_to_o_double = false;
        let mut bonded_to_n_double = false;
        let mut bonded_to_c_double = false;
        for &bi in &bonds {
            let other = mol.bonds[bi].other(i).unwrap();
            if mol.bonds[bi].order == BondOrder::Double {
                match mol.atoms[other].atomic_number {
                    8 => bonded_to_o_double = true,
                    7 => bonded_to_n_double = true,
                    6 => bonded_to_c_double = true,
                    _ => {}
                }
            }
        }
        // Carboxylate / nitro test: C bound to two equivalent O via a
        // resonance-delocalised pair → CO2M (41). Detect as: C bonded
        // to two oxygens via a (single / double) pair where at least
        // one of the oxygens has formal charge −1 (or is a terminal
        // O2CM bound only to this carbon).
        if bonded_to_o_double && is_carboxylate_center(mol, i) {
            return CO2M;
        }
        if bonded_to_o_double || bonded_to_n_double {
            return C_E_O;
        }
        if bonded_to_c_double {
            return C_E_C;
        }
    }
    // Single-bonded only — but a carboxylate carbon in its neutral acid
    // form (one C=O, one C–OH) is *also* CO2M in MMFF94's typing for
    // the partial charges to come out right. We covered the double-bond
    // path above; here we cover the deprotonated case where one O is
    // O⁻ and the other is the C=O (matched above).
    CR
}

/// `true` when carbon `i` is the central carbon of a carboxylate /
/// sulfonate-style two-oxygen anion: it has at least two single+double
/// bonds to oxygens of which one is a charged terminal O⁻ (or a real
/// O2CM via a charged neighbour). The simple +/− carboxylate (`C(=O)[O-]`).
fn is_carboxylate_center(mol: &Molecule, i: usize) -> bool {
    let bonds = mol.bonds_on(i);
    let mut o_single_minus = false;
    let mut o_double = false;
    for &bi in &bonds {
        let other = mol.bonds[bi].other(i).unwrap();
        if mol.atoms[other].atomic_number != 8 {
            continue;
        }
        match mol.bonds[bi].order {
            BondOrder::Single if mol.atoms[other].formal_charge < 0 => o_single_minus = true,
            BondOrder::Double => o_double = true,
            _ => {}
        }
    }
    o_single_minus && o_double
}

// --- Nitrogen --------------------------------------------------------

fn type_nitrogen(mol: &Molecule, i: usize, _rings: &RingInfo) -> MmffType {
    let a = &mol.atoms[i];
    // Aromatic: pyridine vs pyrrole. Pyridine-type N (NPYD) is the lone
    // pair sticking out — it has an in-ring π bond (a real C=N in the
    // Kekulé form, or contributes 1 e to the ring). Pyrrole-type N
    // (NPYL) donates its lone pair (2 e) to the ring — it is single-
    // bonded only inside the ring; in the Kekulé form it has no in-ring
    // double bond. Heuristic: if the aromatic N has any double / triple
    // ring bond it is NPYD; otherwise NPYL.
    if a.aromatic {
        // NPYL (pyrrole-type) donates its lone pair to the ring, so it
        // carries a substituent (H or R) → degree ≥ 3 (or total_h ≥ 1
        // when the H is still implicit). NPYD (pyridine-type) has the
        // lone pair out of the ring → degree 2 (just the two ring
        // neighbours) and no H.
        //
        // `mol.neighbors` counts whatever atoms are nodes (so after
        // add_explicit_hydrogens the H counts). `a.total_h()` counts
        // implicit + bracket-explicit H — the pre-expansion form. We
        // sum both so this works in either representation.
        let total_subs = mol.neighbors(i).len() + usize::from(a.total_h());
        if total_subs >= 3 {
            return NPYL;
        }
        return NPYD;
    }
    // Amide / sulfonamide N: a single-bonded N adjacent to a C=O / S(=O).
    if is_amide_nitrogen(mol, i) {
        return NC_E_O;
    }
    // Enamine / vinylogous amide N: N single-bonded to an sp² carbon.
    if is_enamine_nitrogen(mol, i) {
        return NC_E_C;
    }
    let bonds = mol.bonds_on(i);
    let has_double_or_triple = bonds
        .iter()
        .any(|&b| matches!(mol.bonds[b].order, BondOrder::Double | BondOrder::Triple));
    if has_double_or_triple {
        return N_E_C;
    }
    NR
}

fn is_amide_nitrogen(mol: &Molecule, n: usize) -> bool {
    for &nbr in &mol.neighbors(n) {
        // Look for an adjacent sp² carbon whose double-bond partner is
        // an oxygen (amide / urea) or for an adjacent sulfur in S(=O).
        match mol.atoms[nbr].atomic_number {
            6 => {
                for &bi in &mol.bonds_on(nbr) {
                    if mol.bonds[bi].order == BondOrder::Double {
                        let other = mol.bonds[bi].other(nbr).unwrap();
                        if mol.atoms[other].atomic_number == 8 {
                            return true;
                        }
                    }
                }
            }
            16 => {
                // sulfonamide
                let oxy = mol
                    .bonds_on(nbr)
                    .iter()
                    .filter(|&&b| mol.bonds[b].order == BondOrder::Double)
                    .filter_map(|&b| {
                        let o = mol.bonds[b].other(nbr).unwrap();
                        (mol.atoms[o].atomic_number == 8).then_some(o)
                    })
                    .count();
                if oxy >= 1 {
                    return true;
                }
            }
            _ => {}
        }
    }
    false
}

fn is_enamine_nitrogen(mol: &Molecule, n: usize) -> bool {
    for &nbr in &mol.neighbors(n) {
        if mol.atoms[nbr].atomic_number != 6 {
            continue;
        }
        // sp² C: has a double bond to another C (and not amide/imine)
        for &bi in &mol.bonds_on(nbr) {
            if mol.bonds[bi].order == BondOrder::Double {
                let other = mol.bonds[bi].other(nbr).unwrap();
                if mol.atoms[other].atomic_number == 6 && other != n {
                    return true;
                }
            }
        }
    }
    false
}

// --- Oxygen ----------------------------------------------------------

fn type_oxygen(mol: &Molecule, i: usize) -> MmffType {
    let a = &mol.atoms[i];
    // O⁻ alkoxide / oxide / carboxylate terminal.
    if a.formal_charge <= -1 {
        let nbrs = mol.neighbors(i);
        if nbrs.len() == 1 {
            let nb = nbrs[0];
            // O2CM if the central atom has another oxygen bond
            // (carboxylate / sulfonate / nitro / phosphate pattern).
            // The other oxygen may be single- or double-bonded — the
            // resonance equalises them.
            let central_bonds = mol.bonds_on(nb);
            let other_o = central_bonds.iter().any(|&bi| {
                let bond = &mol.bonds[bi];
                let other = bond.other(nb).unwrap();
                other != i && mol.atoms[other].atomic_number == 8
            });
            if other_o {
                return O2CM;
            }
            return OM;
        }
        return OM;
    }
    let bonds = mol.bonds_on(i);
    let has_double = bonds
        .iter()
        .any(|&b| mol.bonds[b].order == BondOrder::Double);
    if has_double {
        // O=C on a carboxylate centre is itself O2CM (equivalent O
        // pair). Detect: there's another oxygen bonded to the same
        // carbon that carries a negative charge or is single-bonded.
        let nbrs = mol.neighbors(i);
        if nbrs.len() == 1 {
            let nb = nbrs[0];
            // Look for a sister oxygen with O⁻ or just a single-bonded O
            // alongside this =O on the same carbon — that's the
            // carboxylate pattern.
            let central_bonds = mol.bonds_on(nb);
            let charged_o_sister = central_bonds.iter().any(|&bi| {
                let bond = &mol.bonds[bi];
                let other = bond.other(nb).unwrap();
                other != i
                    && mol.atoms[other].atomic_number == 8
                    && mol.atoms[other].formal_charge < 0
            });
            if charged_o_sister {
                return O2CM;
            }
        }
        return O_E_C;
    }
    OR
}

// --- Sulfur ----------------------------------------------------------

fn type_sulfur(mol: &Molecule, i: usize) -> MmffType {
    let bonds = mol.bonds_on(i);
    let oxy_double = bonds
        .iter()
        .filter(|&&b| mol.bonds[b].order == BondOrder::Double)
        .filter(|&&b| {
            let other = mol.bonds[b].other(i).unwrap();
            mol.atoms[other].atomic_number == 8
        })
        .count();
    match oxy_double {
        0 => {
            // S=C thione vs divalent SR.
            let has_double = bonds
                .iter()
                .any(|&b| mol.bonds[b].order == BondOrder::Double);
            if has_double {
                S_E_C
            } else {
                SR
            }
        }
        1 => SO,
        _ => SO2,
    }
}

// --- Phosphorus ------------------------------------------------------

fn type_phosphorus(mol: &Molecule, i: usize) -> MmffType {
    let bonds = mol.bonds_on(i);
    let has_pec = bonds.iter().any(|&b| {
        mol.bonds[b].order == BondOrder::Double && {
            let other = mol.bonds[b].other(i).unwrap();
            mol.atoms[other].atomic_number == 6
        }
    });
    if has_pec {
        PEC
    } else {
        PR
    }
}

// --- Hydrogens --------------------------------------------------------

fn perceive_hydrogen(mol: &Molecule, i: usize, types: &[MmffType]) -> MmffType {
    let nbrs = mol.neighbors(i);
    if nbrs.is_empty() {
        return MmffType::UNKNOWN;
    }
    let heavy = nbrs[0];
    let parent_type = types[heavy];
    match mol.atoms[heavy].atomic_number {
        6 => HC,
        7 => HNR,
        8 => {
            // Discriminate alcohol / phenol / carboxyl / enol H.
            // Look at what the oxygen is attached to.
            for &on in &mol.neighbors(heavy) {
                if on == i {
                    continue;
                }
                let z = mol.atoms[on].atomic_number;
                if z == 16 || z == 15 {
                    return HOS;
                }
                if z == 6 && parent_type == OR {
                    // Carboxyl acid OH (parent C is CO2M / C=O) →
                    // HOCO. Enol / phenol OH (parent C is CB or
                    // C=C) → HOCC.
                    let other_type = types[on];
                    if other_type == C_E_O || other_type == CO2M {
                        return HOCO;
                    }
                    if other_type == CB || other_type == C_E_C {
                        return HOCC;
                    }
                    return HOR;
                }
            }
            HOR
        }
        16 => HS,
        _ => MmffType::UNKNOWN,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;
    use crate::perceive::hydrogen::add_explicit_hydrogens;

    fn typed(smiles: &str) -> (Vec<MmffType>, crate::molecule::Molecule) {
        let m = mol_from_smiles(smiles).expect("parse");
        let m = add_explicit_hydrogens(&m);
        (type_molecule(&m), m)
    }

    fn heavy_types(smiles: &str) -> Vec<&'static str> {
        let (t, m) = typed(smiles);
        m.atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.atomic_number != 1)
            .map(|(i, _)| t[i].symbol)
            .collect()
    }

    #[test]
    fn ethane_is_all_cr_hc() {
        let (t, m) = typed("CC");
        let heavy: Vec<_> = m
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.atomic_number != 1)
            .map(|(i, _)| t[i])
            .collect();
        assert!(
            heavy.iter().all(|x| *x == CR),
            "ethane heavies should be CR"
        );
        let h: Vec<_> = m
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.atomic_number == 1)
            .map(|(i, _)| t[i])
            .collect();
        assert!(h.iter().all(|x| *x == HC), "ethane H's should be HC");
    }

    #[test]
    fn water_is_or_hor() {
        let (t, m) = typed("O");
        let h: Vec<_> = m
            .atoms
            .iter()
            .enumerate()
            .map(|(i, a)| (a.atomic_number, t[i]))
            .collect();
        let o = h.iter().find(|(z, _)| *z == 8).unwrap().1;
        assert_eq!(o, OR);
        for (z, ty) in &h {
            if *z == 1 {
                assert_eq!(*ty, HOR);
            }
        }
    }

    #[test]
    fn methanol_oxygen_is_or() {
        let s = heavy_types("CO");
        assert!(s.contains(&"CR"));
        assert!(s.contains(&"OR"));
    }

    #[test]
    fn benzene_is_cb() {
        let s = heavy_types("c1ccccc1");
        assert_eq!(s.iter().filter(|&&x| x == "CB").count(), 6);
    }

    #[test]
    fn methylamine_is_cr_nr() {
        let s = heavy_types("CN");
        assert!(s.contains(&"CR"));
        assert!(s.contains(&"NR"));
    }

    #[test]
    fn acetic_acid_has_co2m_or_carboxyl_carbon() {
        // CC(=O)O — neutral acid: carbonyl C is "C=O" (3), carboxyl O
        // is OR (6). Even when not deprotonated MMFF94 does not assign
        // CO2M to a neutral COOH — that goes to the [O-] form below.
        let s = heavy_types("CC(=O)O");
        // Carbonyl carbon must be C=O (3); the OH oxygen is OR (6); the
        // double-bonded oxygen is O=C (7); the methyl carbon is CR.
        assert!(s.contains(&"C=O"));
        assert!(s.contains(&"O=C"));
        assert!(s.contains(&"OR"));
        assert!(s.contains(&"CR"));
    }

    #[test]
    fn acetate_anion_uses_co2m_o2cm() {
        // CC(=O)[O-] — deprotonated. Central carbon becomes CO2M (41),
        // both oxygens become O2CM (32).
        let s = heavy_types("CC(=O)[O-]");
        assert!(s.contains(&"CO2M"));
        assert_eq!(s.iter().filter(|&&x| x == "O2CM").count(), 2);
    }

    #[test]
    fn methanethiol_is_s_hs() {
        let (t, m) = typed("CS");
        let s = m
            .atoms
            .iter()
            .enumerate()
            .find(|(_, a)| a.atomic_number == 16)
            .unwrap()
            .0;
        assert_eq!(t[s], SR);
        // the S-bound H
        let has_hs = m
            .atoms
            .iter()
            .enumerate()
            .any(|(i, a)| a.atomic_number == 1 && t[i] == HS);
        assert!(has_hs, "methanethiol should have an HS hydrogen");
    }

    #[test]
    fn methyl_chloride_has_cl_type() {
        let s = heavy_types("CCl");
        assert!(s.contains(&"CL"));
    }

    #[test]
    fn methylphosphine_is_p() {
        let s = heavy_types("CP");
        assert!(s.contains(&"P"));
    }

    #[test]
    fn acetone_carbonyl_is_c_eq_o() {
        // (CH3)2C=O — central C is C=O, the =O is O=C.
        let s = heavy_types("CC(=O)C");
        assert_eq!(s.iter().filter(|&&x| x == "C=O").count(), 1);
        assert_eq!(s.iter().filter(|&&x| x == "O=C").count(), 1);
    }

    #[test]
    fn acetonitrile_has_sp_carbon() {
        let s = heavy_types("CC#N");
        assert!(s.contains(&"CSP"));
        // and the N
        assert!(s.contains(&"N=C"));
    }

    #[test]
    fn pyridine_n_is_npyd() {
        let s = heavy_types("c1ccncc1");
        assert!(s.contains(&"NPYD"));
    }

    #[test]
    fn pyrrole_n_is_npyl() {
        // Kekulé pyrrole — N has an N-H, so degree 3 after H expansion,
        // → NPYL.
        let s = heavy_types("C1=CC=CN1");
        assert!(s.contains(&"NPYL"), "got {s:?}");
    }

    #[test]
    fn acetamide_n_is_nc_eq_o() {
        // CH3-C(=O)-NH2
        let s = heavy_types("CC(=O)N");
        assert!(s.contains(&"NC=O"));
    }

    #[test]
    fn dimethyl_sulfone_has_so2() {
        // (CH3)2S(=O)(=O)
        let s = heavy_types("CS(=O)(=O)C");
        assert!(s.contains(&"SO2"));
    }

    #[test]
    fn dimethyl_sulfoxide_has_so() {
        let s = heavy_types("CS(=O)C");
        assert!(s.contains(&"S=O"));
    }
}
