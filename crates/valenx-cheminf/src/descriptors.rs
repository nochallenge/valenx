//! Physicochemical descriptors and drug-likeness rules.
//!
//! The descriptors a medicinal chemist computes first:
//!
//! - [`crippen_logp`] — the octanol/water partition coefficient by
//!   Wildman-Crippen atom contributions;
//! - [`tpsa`] — the topological polar surface area (Ertl fragment
//!   contributions);
//! - [`hbd`] / [`hba`] — hydrogen-bond donor / acceptor counts;
//! - [`rotatable_bonds`] — count of freely-rotating single bonds;
//! - [`ring_count`], [`aromatic_ring_count`], [`heavy_atom_count`],
//!   [`fraction_csp3`], [`formal_charge`];
//!
//! plus [`lipinski`] and [`veber`], which roll the descriptors into the
//! classic oral-bioavailability rule sets.
//!
//! **v1 fidelity.** The Crippen logP uses a representative subset of
//! the published atom-type contributions keyed on element +
//! aromaticity + neighbour pattern; TPSA uses the Ertl polar-fragment
//! table. Both are within a small offset of RDKit on drug-like
//! molecules but are not bit-identical (RDKit ships the full ~110-type
//! Crippen table).

use crate::molecule::{BondOrder, Molecule};
use crate::perceive::rings::sssr;

/// Wildman-Crippen logP — sum of per-atom octanol/water contributions.
///
/// Hydrogens contribute via their heavy-atom carrier (the contributions
/// below already fold a typical hydrogen environment in).
pub fn crippen_logp(mol: &Molecule) -> f64 {
    let mut sum = 0.0;
    for i in 0..mol.atoms.len() {
        sum += atom_logp_contribution(mol, i);
    }
    sum
}

/// Per-atom Crippen logP contribution.
fn atom_logp_contribution(mol: &Molecule, i: usize) -> f64 {
    let a = &mol.atoms[i];
    let h = f64::from(a.total_h());
    // hydrogens carry a small positive contribution each
    let h_contrib = 0.123 * h;
    let heavy = match a.atomic_number {
        6 => {
            // carbon: aromatic vs aliphatic, with a small polar-
            // neighbour penalty
            let polar_nbr = mol
                .neighbors(i)
                .iter()
                .any(|&v| matches!(mol.atoms[v].atomic_number, 7 | 8));
            if a.aromatic {
                if polar_nbr {
                    0.151
                } else {
                    0.296
                }
            } else if polar_nbr {
                -0.038
            } else {
                0.218
            }
        }
        7 => {
            // nitrogen: amines are quite hydrophilic
            if a.aromatic {
                -0.329
            } else if h > 0.0 {
                -0.806
            } else {
                -0.448
            }
        }
        8 => {
            // oxygen: hydroxyl vs ether/carbonyl
            if h > 0.0 {
                -0.467
            } else if mol
                .bonds_on(i)
                .iter()
                .any(|&b| mol.bonds[b].order == BondOrder::Double)
            {
                0.007 // carbonyl O
            } else {
                -0.295 // ether O
            }
        }
        9 => 0.412,   // fluorine — lipophilic
        17 => 0.638,  // chlorine
        35 => 0.825,  // bromine
        53 => 1.009,  // iodine
        16 => 0.640,  // sulfur
        15 => 0.230,  // phosphorus
        1 => 0.0,     // explicit H node handled by h_contrib of carrier
        _ => 0.0,
    };
    // an explicit-H node contributes on its own
    if a.atomic_number == 1 {
        return 0.123;
    }
    heavy + h_contrib
}

/// Topological polar surface area (Å²) — Ertl fragment contributions
/// for the polar N / O atoms.
pub fn tpsa(mol: &Molecule) -> f64 {
    let mut sum = 0.0;
    for i in 0..mol.atoms.len() {
        let a = &mol.atoms[i];
        let h = a.total_h();
        match a.atomic_number {
            7 => sum += nitrogen_tpsa(mol, i, h),
            8 => sum += oxygen_tpsa(mol, i, h),
            16 => sum += sulfur_tpsa(mol, i, h),
            15 => sum += 13.59, // generic trivalent P
            _ => {}
        }
    }
    sum
}

fn nitrogen_tpsa(mol: &Molecule, i: usize, h: u8) -> f64 {
    let a = &mol.atoms[i];
    let degree = mol.degree(i);
    let charged = a.formal_charge != 0;
    if a.aromatic {
        match h {
            0 => 12.89, // pyridine-type
            _ => 15.79, // pyrrole-type N-H
        }
    } else if charged {
        24.0 // protonated amine, rough
    } else {
        match (degree, h) {
            (1, 2) => 26.02, // primary amine
            (2, 1) => 12.03, // secondary amine
            (3, 0) => 3.24,  // tertiary amine
            (1, 1) => 23.85, // =N-H imine-ish
            _ => 12.36,      // generic
        }
    }
}

fn oxygen_tpsa(mol: &Molecule, i: usize, h: u8) -> f64 {
    let a = &mol.atoms[i];
    if a.formal_charge < 0 {
        return 23.06; // carboxylate / alkoxide O-
    }
    if a.aromatic {
        return 13.14; // furan-type
    }
    let has_double = mol
        .bonds_on(i)
        .iter()
        .any(|&b| mol.bonds[b].order == BondOrder::Double);
    if has_double {
        17.07 // carbonyl O
    } else if h > 0 {
        20.23 // hydroxyl
    } else {
        9.23 // ether
    }
}

fn sulfur_tpsa(mol: &Molecule, i: usize, _h: u8) -> f64 {
    let doubles = mol
        .bonds_on(i)
        .iter()
        .filter(|&&b| mol.bonds[b].order == BondOrder::Double)
        .count();
    match doubles {
        0 => 25.30, // thioether / thiol
        1 => 19.21, // sulfoxide
        _ => 8.38,  // sulfone
    }
}

/// Hydrogen-bond donor count — N/O atoms bearing at least one hydrogen
/// (the Lipinski definition: count of N-H and O-H).
pub fn hbd(mol: &Molecule) -> usize {
    mol.atoms
        .iter()
        .filter(|a| matches!(a.atomic_number, 7 | 8) && a.total_h() > 0)
        .count()
}

/// Hydrogen-bond acceptor count — every nitrogen and oxygen (the
/// Lipinski definition: count of N + O atoms).
pub fn hba(mol: &Molecule) -> usize {
    mol.atoms
        .iter()
        .filter(|a| matches!(a.atomic_number, 7 | 8))
        .count()
}

/// Count of rotatable bonds — a single, acyclic bond between two
/// non-terminal heavy atoms, excluding amide C-N bonds.
pub fn rotatable_bonds(mol: &Molecule) -> usize {
    let rings = sssr(mol);
    let mut count = 0;
    for (bi, b) in mol.bonds.iter().enumerate() {
        if b.order != BondOrder::Single {
            continue;
        }
        if rings.bond_in_ring(bi) {
            continue;
        }
        // both endpoints must be non-terminal heavy atoms
        let (a, c) = (b.a, b.b);
        if mol.atoms[a].is_hydrogen() || mol.atoms[c].is_hydrogen() {
            continue;
        }
        let deg_a = mol.degree(a);
        let deg_c = mol.degree(c);
        if deg_a < 2 || deg_c < 2 {
            continue;
        }
        // skip amide bonds: a C(=O)-N single bond
        if is_amide_bond(mol, a, c) {
            continue;
        }
        count += 1;
    }
    count
}

/// Is the single bond between `a` and `c` an amide C-N bond?
fn is_amide_bond(mol: &Molecule, a: usize, c: usize) -> bool {
    let pair = [(a, c), (c, a)];
    for &(carbon, nitrogen) in &pair {
        if mol.atoms[carbon].atomic_number == 6 && mol.atoms[nitrogen].atomic_number == 7 {
            // the carbon must have a double bond to an oxygen
            let carbonyl = mol.bonds_on(carbon).iter().any(|&bi| {
                let bd = &mol.bonds[bi];
                bd.order == BondOrder::Double
                    && bd
                        .other(carbon)
                        .map(|o| mol.atoms[o].atomic_number == 8)
                        .unwrap_or(false)
            });
            if carbonyl {
                return true;
            }
        }
    }
    false
}

/// Number of explicit double bonds (`BondOrder::Double`, e.g. C=C, C=O). Aromatic ring
/// bonds are `BondOrder::Aromatic` (not Double) so they are NOT counted: benzene = 0,
/// ethylene = 1, CO₂ = 2.
pub fn double_bond_count(mol: &Molecule) -> usize {
    mol.bonds
        .iter()
        .filter(|b| b.order == BondOrder::Double)
        .count()
}

/// Number of explicit triple bonds (`BondOrder::Triple`, e.g. C#C, C#N). Aromatic ring
/// bonds are `BondOrder::Aromatic` (not Triple) so they are NOT counted.
pub fn triple_bond_count(mol: &Molecule) -> usize {
    mol.bonds
        .iter()
        .filter(|b| b.order == BondOrder::Triple)
        .count()
}

/// Number of SSSR rings.
pub fn ring_count(mol: &Molecule) -> usize {
    sssr(mol).ring_count()
}

/// Number of SSSR rings that are aromatic (every atom aromatic).
pub fn aromatic_ring_count(mol: &Molecule) -> usize {
    let info = sssr(mol);
    info.rings
        .iter()
        .filter(|r| r.atoms.iter().all(|&a| mol.atoms[a].aromatic))
        .count()
}

/// Number of atoms flagged aromatic by Hückel (4n+2) perception — the atom-level
/// companion to [`aromatic_ring_count`] (e.g. 6 for benzene, 10 for naphthalene,
/// 0 for cyclohexane).
pub fn aromatic_atom_count(mol: &Molecule) -> usize {
    mol.atoms.iter().filter(|a| a.aromatic).count()
}

/// Number of bonds whose **both** endpoint atoms are aromatic — the bond-level
/// companion to [`aromatic_atom_count`] (atoms) and [`aromatic_ring_count`] (rings),
/// e.g. 6 for benzene, 11 for naphthalene, 0 for cyclohexane.
pub fn aromatic_bond_count(mol: &Molecule) -> usize {
    mol.bonds
        .iter()
        .filter(|b| mol.atoms[b.a].aromatic && mol.atoms[b.b].aromatic)
        .count()
}

/// Heavy (non-hydrogen) atom count.
pub fn heavy_atom_count(mol: &Molecule) -> usize {
    mol.heavy_atom_count()
}

/// Heteroatom count — atoms that are neither carbon nor hydrogen (atomic number
/// outside {1, 6}), e.g. 1 for ethanol (O), 2 for acetic acid (2 O), 0 for benzene.
pub fn heteroatom_count(mol: &Molecule) -> usize {
    mol.atoms
        .iter()
        .filter(|a| a.atomic_number != 1 && a.atomic_number != 6)
        .count()
}

/// Halogen atom count — atoms that are F, Cl, Br or I (atomic number ∈ {9, 17, 35, 53});
/// a distinct subset of [`heteroatom_count`]. E.g. 1 for chlorobenzene, 3 for chloroform,
/// 4 for CF₄, 0 for benzene.
pub fn halogen_count(mol: &Molecule) -> usize {
    mol.atoms
        .iter()
        .filter(|a| matches!(a.atomic_number, 9 | 17 | 35 | 53))
        .count()
}

/// Fraction of carbons that are `sp³` (Csp3) — a flatness / 3D-shape
/// descriptor (`Fsp3`).
pub fn fraction_csp3(mol: &Molecule) -> f64 {
    let carbons: Vec<usize> = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.atomic_number == 6)
        .map(|(i, _)| i)
        .collect();
    if carbons.is_empty() {
        return 0.0;
    }
    let sp3 = carbons
        .iter()
        .filter(|&&c| {
            !mol.atoms[c].aromatic
                && !mol
                    .bonds_on(c)
                    .iter()
                    .any(|&b| matches!(mol.bonds[b].order, BondOrder::Double | BondOrder::Triple))
        })
        .count();
    sp3 as f64 / carbons.len() as f64
}

/// Net formal charge.
pub fn formal_charge(mol: &Molecule) -> i32 {
    mol.atoms.iter().map(|a| i32::from(a.formal_charge)).sum()
}

/// The Lipinski "rule of five" verdict.
#[derive(Clone, Debug, PartialEq)]
pub struct LipinskiResult {
    /// Average molecular weight (g/mol). Rule: ≤ 500.
    pub molecular_weight: f64,
    /// Crippen logP. Rule: ≤ 5.
    pub logp: f64,
    /// H-bond donor count. Rule: ≤ 5.
    pub hbd: usize,
    /// H-bond acceptor count. Rule: ≤ 10.
    pub hba: usize,
    /// Number of the four rules violated (`0` = fully compliant).
    pub violations: u8,
}

impl LipinskiResult {
    /// `true` if the molecule passes (≤ 1 violation, the usual
    /// "drug-like" threshold).
    pub fn passes(&self) -> bool {
        self.violations <= 1
    }
}

/// Evaluate the Lipinski rule of five.
pub fn lipinski(mol: &Molecule) -> LipinskiResult {
    let mw = crate::perceive::formula::average_molecular_weight(mol);
    let logp = crippen_logp(mol);
    let d = hbd(mol);
    let a = hba(mol);
    let mut v = 0u8;
    if mw > 500.0 {
        v += 1;
    }
    if logp > 5.0 {
        v += 1;
    }
    if d > 5 {
        v += 1;
    }
    if a > 10 {
        v += 1;
    }
    LipinskiResult {
        molecular_weight: mw,
        logp,
        hbd: d,
        hba: a,
        violations: v,
    }
}

/// The Veber oral-bioavailability rule verdict.
#[derive(Clone, Debug, PartialEq)]
pub struct VeberResult {
    /// Rotatable-bond count. Rule: ≤ 10.
    pub rotatable_bonds: usize,
    /// Topological polar surface area. Rule: ≤ 140 Å².
    pub tpsa: f64,
    /// `true` if both Veber criteria are met.
    pub passes: bool,
}

/// Evaluate the Veber rules (rotatable bonds ≤ 10, TPSA ≤ 140 Å²).
pub fn veber(mol: &Molecule) -> VeberResult {
    let rb = rotatable_bonds(mol);
    let psa = tpsa(mol);
    VeberResult {
        rotatable_bonds: rb,
        tpsa: psa,
        passes: rb <= 10 && psa <= 140.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;

    #[test]
    fn hbd_hba_counts() {
        // ethanol: 1 donor (O-H), 1 acceptor (O)
        let m = mol_from_smiles("CCO").unwrap();
        assert_eq!(hbd(&m), 1);
        assert_eq!(hba(&m), 1);

        // acetic acid: 1 donor, 2 oxygens → 2 acceptors
        let m = mol_from_smiles("CC(=O)O").unwrap();
        assert_eq!(hbd(&m), 1);
        assert_eq!(hba(&m), 2);
    }

    #[test]
    fn rotatable_bond_count() {
        // n-butane: C-C-C-C, the central C-C is rotatable
        let m = mol_from_smiles("CCCC").unwrap();
        assert_eq!(rotatable_bonds(&m), 1);
        // benzene: ring bonds are not rotatable
        let m = mol_from_smiles("c1ccccc1").unwrap();
        assert_eq!(rotatable_bonds(&m), 0);
    }

    #[test]
    fn double_bond_counts() {
        // ethylene C=C: 1 double; CO₂ O=C=O: 2; acetic acid CC(=O)O: 1.
        assert_eq!(double_bond_count(&mol_from_smiles("C=C").unwrap()), 1);
        assert_eq!(double_bond_count(&mol_from_smiles("O=C=O").unwrap()), 2);
        assert_eq!(double_bond_count(&mol_from_smiles("CC(=O)O").unwrap()), 1);
        // ethane: none; benzene: aromatic bonds are NOT Double → 0.
        assert_eq!(double_bond_count(&mol_from_smiles("CC").unwrap()), 0);
        assert_eq!(double_bond_count(&mol_from_smiles("c1ccccc1").unwrap()), 0);
    }

    #[test]
    fn triple_bond_counts() {
        // acetylene C#C: 1 triple; N₂ N#N: 1; HCN C#N: 1.
        assert_eq!(triple_bond_count(&mol_from_smiles("C#C").unwrap()), 1);
        assert_eq!(triple_bond_count(&mol_from_smiles("N#N").unwrap()), 1);
        assert_eq!(triple_bond_count(&mol_from_smiles("C#N").unwrap()), 1);
        // ethane (single), ethylene (double), benzene (aromatic) → 0 triples.
        assert_eq!(triple_bond_count(&mol_from_smiles("CC").unwrap()), 0);
        assert_eq!(triple_bond_count(&mol_from_smiles("C=C").unwrap()), 0);
        assert_eq!(triple_bond_count(&mol_from_smiles("c1ccccc1").unwrap()), 0);
    }

    #[test]
    fn amide_bond_not_rotatable() {
        // N-methylacetamide CC(=O)NC has NO rotatable bonds: every
        // single bond either touches a terminal methyl (a bond to a
        // terminal atom is not rotatable — rotating it changes no
        // conformation) or is the amide C(=O)-N bond (excluded for its
        // partial double-bond character). This matches RDKit's
        // NumRotatableBonds.
        let m = mol_from_smiles("CC(=O)NC").unwrap();
        assert_eq!(rotatable_bonds(&m), 0);

        // N-ethyl propanamide CCC(=O)NCC: here two bonds join
        // non-terminal heavy atoms — the C-C(=O) and the N-C — so it has
        // 2 rotatable bonds, while the amide C(=O)-N stays excluded.
        // Without the amide exclusion the count would be 3.
        let amide = mol_from_smiles("CCC(=O)NCC").unwrap();
        assert_eq!(rotatable_bonds(&amide), 2);
    }

    #[test]
    fn tpsa_is_positive_for_polar() {
        let polar = mol_from_smiles("CCO").unwrap();
        assert!(tpsa(&polar) > 15.0);
        // pure hydrocarbon has zero TPSA
        let hydro = mol_from_smiles("CCCCCC").unwrap();
        assert_eq!(tpsa(&hydro), 0.0);
    }

    #[test]
    fn logp_hydrocarbon_positive_alcohol_lower() {
        let hexane = crippen_logp(&mol_from_smiles("CCCCCC").unwrap());
        let hexanol = crippen_logp(&mol_from_smiles("CCCCCCO").unwrap());
        assert!(hexane > 0.0);
        assert!(hexanol < hexane, "the OH lowers logP");
    }

    #[test]
    fn lipinski_passes_small_molecule() {
        let m = mol_from_smiles("CCO").unwrap();
        let r = lipinski(&m);
        assert_eq!(r.violations, 0);
        assert!(r.passes());
    }

    #[test]
    fn veber_passes_simple() {
        let m = mol_from_smiles("c1ccccc1").unwrap();
        let r = veber(&m);
        assert!(r.passes);
    }

    #[test]
    fn ring_counts() {
        let naph = mol_from_smiles("c1ccc2ccccc2c1").unwrap();
        assert_eq!(ring_count(&naph), 2);
        assert_eq!(aromatic_ring_count(&naph), 2);
    }

    #[test]
    fn aromatic_atom_counts() {
        // benzene: all 6 ring carbons aromatic.
        assert_eq!(aromatic_atom_count(&mol_from_smiles("c1ccccc1").unwrap()), 6);
        // naphthalene: two fused aromatic rings → 10 aromatic carbons.
        assert_eq!(aromatic_atom_count(&mol_from_smiles("c1ccc2ccccc2c1").unwrap()), 10);
        // pyridine: 5 C + 1 N in the aromatic ring → 6.
        assert_eq!(aromatic_atom_count(&mol_from_smiles("c1ccncc1").unwrap()), 6);
        // cyclohexane (aliphatic) and ethanol: no aromatic atoms.
        assert_eq!(aromatic_atom_count(&mol_from_smiles("C1CCCCC1").unwrap()), 0);
        assert_eq!(aromatic_atom_count(&mol_from_smiles("CCO").unwrap()), 0);
    }

    #[test]
    fn fraction_csp3_extremes() {
        // hexane: all sp3
        assert!((fraction_csp3(&mol_from_smiles("CCCCCC").unwrap()) - 1.0).abs() < 1e-9);
        // benzene: all aromatic, no sp3
        assert_eq!(fraction_csp3(&mol_from_smiles("c1ccccc1").unwrap()), 0.0);
    }

    #[test]
    fn heteroatom_count_exact() {
        // ethanol: 1 O; acetic acid: 2 O; aniline: 1 N; chlorobenzene: 1 Cl.
        assert_eq!(heteroatom_count(&mol_from_smiles("CCO").unwrap()), 1);
        assert_eq!(heteroatom_count(&mol_from_smiles("CC(=O)O").unwrap()), 2);
        assert_eq!(heteroatom_count(&mol_from_smiles("Nc1ccccc1").unwrap()), 1);
        assert_eq!(heteroatom_count(&mol_from_smiles("Clc1ccccc1").unwrap()), 1);
        // pure hydrocarbons: no heteroatoms.
        assert_eq!(heteroatom_count(&mol_from_smiles("c1ccccc1").unwrap()), 0);
        assert_eq!(heteroatom_count(&mol_from_smiles("CCCCCC").unwrap()), 0);
    }

    #[test]
    fn aromatic_bond_count_exact() {
        // benzene: 6 aromatic ring bonds; naphthalene: 11 (two fused rings share 1); pyridine: 6.
        assert_eq!(aromatic_bond_count(&mol_from_smiles("c1ccccc1").unwrap()), 6);
        assert_eq!(aromatic_bond_count(&mol_from_smiles("c1ccc2ccccc2c1").unwrap()), 11);
        assert_eq!(aromatic_bond_count(&mol_from_smiles("c1ccncc1").unwrap()), 6);
        // chlorobenzene still has 6 aromatic ring bonds (the C–Cl bond is not aromatic).
        assert_eq!(aromatic_bond_count(&mol_from_smiles("Clc1ccccc1").unwrap()), 6);
        // no aromatic atoms → no aromatic bonds.
        assert_eq!(aromatic_bond_count(&mol_from_smiles("C1CCCCC1").unwrap()), 0);
    }

    #[test]
    fn halogen_count_exact() {
        // chlorobenzene 1 Cl; chloroform 3 Cl; CF₄ 4 F; benzene + ethanol 0.
        assert_eq!(halogen_count(&mol_from_smiles("Clc1ccccc1").unwrap()), 1);
        assert_eq!(halogen_count(&mol_from_smiles("C(Cl)(Cl)Cl").unwrap()), 3);
        assert_eq!(halogen_count(&mol_from_smiles("FC(F)(F)F").unwrap()), 4);
        assert_eq!(halogen_count(&mol_from_smiles("c1ccccc1").unwrap()), 0);
        assert_eq!(halogen_count(&mol_from_smiles("CCO").unwrap()), 0);
    }
}
