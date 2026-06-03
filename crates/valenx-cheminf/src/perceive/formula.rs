//! Molecular formula and molecular weight.
//!
//! [`molecular_formula`] returns a Hill-ordered formula string
//! (carbon first, hydrogen second, then every other element
//! alphabetically — the universal convention). Hydrogens are counted
//! from explicit H *nodes*, bracket `explicit_h` and perception-added
//! `implicit_h`, so the formula is correct whether or not the molecule
//! has been hydrogen-expanded.
//!
//! [`average_molecular_weight`] sums conventional atomic weights;
//! [`monoisotopic_mass`] sums the most-abundant-isotope exact masses
//! (an explicit isotope on an atom overrides the natural value when it
//! is a tabulated common isotope).

use crate::element;
use crate::molecule::Molecule;
use std::collections::BTreeMap;

/// Element-count map for a molecule, keyed by atomic number, hydrogens
/// folded in. The building block of every formula / weight routine.
pub fn element_counts(mol: &Molecule) -> BTreeMap<u8, u32> {
    let mut counts: BTreeMap<u8, u32> = BTreeMap::new();
    let mut h_total: u32 = 0;
    for a in &mol.atoms {
        if a.is_dummy() {
            continue;
        }
        *counts.entry(a.atomic_number).or_insert(0) += 1;
        if a.atomic_number != 1 {
            h_total += u32::from(a.total_h());
        }
    }
    if h_total > 0 {
        *counts.entry(1).or_insert(0) += h_total;
    }
    counts
}

/// Hill-system molecular formula, e.g. `"C2H6O"` for ethanol.
pub fn molecular_formula(mol: &Molecule) -> String {
    let counts = element_counts(mol);
    let mut out = String::new();

    let emit = |z: u8, n: u32, out: &mut String| {
        if n == 0 {
            return;
        }
        if let Some(e) = element::by_number(z) {
            out.push_str(e.symbol);
        }
        if n > 1 {
            out.push_str(&n.to_string());
        }
    };

    // Hill order: C, then H, then the rest alphabetically by symbol.
    let has_carbon = counts.get(&6).copied().unwrap_or(0) > 0;
    if has_carbon {
        emit(6, counts.get(&6).copied().unwrap_or(0), &mut out);
        emit(1, counts.get(&1).copied().unwrap_or(0), &mut out);
    }
    let mut rest: Vec<(&'static str, u8, u32)> = counts
        .iter()
        .filter(|(&z, _)| {
            if has_carbon {
                z != 6 && z != 1
            } else {
                true
            }
        })
        .filter_map(|(&z, &n)| element::by_number(z).map(|e| (e.symbol, z, n)))
        .collect();
    rest.sort_by(|a, b| a.0.cmp(b.0));
    for (_, z, n) in rest {
        emit(z, n, &mut out);
    }
    if out.is_empty() {
        out.push_str("(empty)");
    }
    out
}

/// Average (conventional) molecular weight in g/mol.
pub fn average_molecular_weight(mol: &Molecule) -> f64 {
    element_counts(mol)
        .iter()
        .filter_map(|(&z, &n)| element::by_number(z).map(|e| e.average_mass * f64::from(n)))
        .sum()
}

/// Monoisotopic (exact) mass in u — the most-abundant isotope of each
/// element, except where an atom carries an explicit `isotope` matching
/// a tabulated common isotope.
pub fn monoisotopic_mass(mol: &Molecule) -> f64 {
    let mut mass = 0.0;
    let mut h_count: u32 = 0;
    for a in &mol.atoms {
        if a.is_dummy() {
            continue;
        }
        if a.atomic_number == 1 {
            mass += hydrogen_mass(a.isotope);
        } else {
            mass += heavy_mass(a.atomic_number, a.isotope);
            h_count += u32::from(a.total_h());
        }
    }
    mass + f64::from(h_count) * 1.007_825_032
}

fn hydrogen_mass(isotope: Option<u16>) -> f64 {
    match isotope {
        Some(2) => 2.014_101_778, // deuterium
        Some(3) => 3.016_049_28,  // tritium
        _ => 1.007_825_032,
    }
}

fn heavy_mass(z: u8, isotope: Option<u16>) -> f64 {
    let e = match element::by_number(z) {
        Some(e) => e,
        None => return 0.0,
    };
    match isotope {
        Some(iso) if iso == e.monoisotopic_number => e.monoisotopic_mass,
        // a non-default isotope: approximate by its mass number
        Some(iso) => f64::from(iso),
        None => e.monoisotopic_mass,
    }
}

/// Total formal charge of the molecule (sum over atoms).
pub fn net_charge(mol: &Molecule) -> i32 {
    mol.atoms.iter().map(|a| i32::from(a.formal_charge)).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn ethanol_formula() {
        let m = parse_smiles("CCO").unwrap();
        assert_eq!(molecular_formula(&m), "C2H6O");
    }

    #[test]
    fn benzene_formula() {
        let m = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(molecular_formula(&m), "C6H6");
    }

    #[test]
    fn water_no_carbon() {
        let m = parse_smiles("O").unwrap();
        assert_eq!(molecular_formula(&m), "H2O");
    }

    #[test]
    fn average_weight_reasonable() {
        let m = parse_smiles("CCO").unwrap();
        let mw = average_molecular_weight(&m);
        assert!((mw - 46.07).abs() < 0.05, "ethanol MW = {mw}");
    }

    #[test]
    fn monoisotopic_vs_average() {
        let m = parse_smiles("c1ccccc1").unwrap();
        let mono = monoisotopic_mass(&m);
        let avg = average_molecular_weight(&m);
        // benzene monoisotopic ~78.0470, average ~78.11
        assert!((mono - 78.0470).abs() < 0.01, "mono = {mono}");
        assert!(avg > mono);
    }

    #[test]
    fn charge_sum() {
        let m = parse_smiles("[Na+].[Cl-]").unwrap();
        assert_eq!(net_charge(&m), 0);
        let m = parse_smiles("[NH4+]").unwrap();
        assert_eq!(net_charge(&m), 1);
    }
}

/// Reference-value validation: molecular formula and average molecular
/// weight of well-known molecules parsed from their SMILES, checked
/// against published values.
#[cfg(test)]
mod validation {
    use super::{average_molecular_weight, molecular_formula};
    use crate::smiles::parser::parse_smiles;

    /// Molecules whose composition is textbook-known: Hill formula and
    /// average molecular weight must match the published values.
    #[test]
    fn known_molecules_have_the_right_formula_and_weight() {
        let cases: &[(&str, &str, f64)] = &[
            // Benzene — C6H6, 78.11 g/mol.
            ("c1ccccc1", "C6H6", 78.11),
            // Water — H2O, 18.015 g/mol.
            ("O", "H2O", 18.015),
            // Ethanol — C2H6O, 46.07 g/mol.
            ("CCO", "C2H6O", 46.07),
            // Caffeine — C8H10N4O2, 194.19 g/mol.
            ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "C8H10N4O2", 194.19),
            // Aspirin (acetylsalicylic acid) — C9H8O4, 180.16 g/mol.
            ("CC(=O)OC1=CC=CC=C1C(=O)O", "C9H8O4", 180.16),
        ];
        for &(smiles, formula, mw) in cases {
            let m = parse_smiles(smiles)
                .unwrap_or_else(|e| panic!("parse {smiles}: {e:?}"));
            assert_eq!(
                molecular_formula(&m),
                formula,
                "Hill formula of {smiles}"
            );
            let got = average_molecular_weight(&m);
            assert!(
                (got - mw).abs() < 0.05,
                "{smiles}: MW {got:.3} vs published {mw}"
            );
        }
    }
}
