//! Implicit / explicit hydrogen handling.
//!
//! The SMILES parser already fills [`Atom::implicit_h`](crate::molecule::Atom::implicit_h)
//! from valence; this module gives the two transformations every
//! cheminformatics toolkit needs:
//!
//! - [`add_explicit_hydrogens`] — turn every implicit hydrogen into a
//!   real H atom node with a single bond, so 3D-coordinate generation
//!   and force fields see them;
//! - [`remove_explicit_hydrogens`] — collapse plain terminal H nodes
//!   back into the heavy atom's implicit count (the inverse), keeping
//!   isotopically-labelled or charged hydrogens as nodes;
//! - [`recompute_implicit_hydrogens`] — re-derive the implicit count
//!   from the current bond graph (used after an edit / a reaction).

use crate::element;
use crate::molecule::{Atom, Bond, Molecule};

/// Return a copy of `mol` with every implicit hydrogen materialised as
/// an explicit H atom + single bond. The heavy atoms keep their indices
/// (new H nodes are appended); each heavy atom's `implicit_h` becomes 0.
pub fn add_explicit_hydrogens(mol: &Molecule) -> Molecule {
    let mut out = mol.clone();
    let original = out.atoms.len();
    for i in 0..original {
        let need = out.atoms[i].implicit_h;
        if need == 0 || out.atoms[i].is_hydrogen() {
            continue;
        }
        out.atoms[i].implicit_h = 0;
        for _ in 0..need {
            let h = out.add_atom(Atom::new(1));
            out.add_bond(Bond::single(i, h));
            if !out.coords.is_empty() {
                // place the new H at the heavy atom for now; a layout /
                // conformer pass will move it. keeps lengths consistent.
                out.coords.push(out.coords[i]);
            }
        }
    }
    out
}

/// Return a copy of `mol` with plain terminal hydrogens collapsed back
/// into their heavy neighbour's implicit count. A hydrogen is *kept* as
/// a node when it carries an isotope label, a formal charge, or is not
/// a simple single-bonded terminal atom (e.g. an `H2` molecule).
pub fn remove_explicit_hydrogens(mol: &Molecule) -> Molecule {
    let n = mol.atoms.len();
    let mut keep = vec![true; n];
    // count implicit additions for each surviving heavy atom
    let mut extra_h = vec![0u8; n];

    for (i, a) in mol.atoms.iter().enumerate() {
        if !a.is_hydrogen() {
            continue;
        }
        if a.isotope.is_some() || a.formal_charge != 0 {
            continue; // labelled / charged H — keep as a node
        }
        let nbrs = mol.neighbors(i);
        if nbrs.len() != 1 {
            continue; // H2, bridging H — keep
        }
        let heavy = nbrs[0];
        if mol.atoms[heavy].is_hydrogen() {
            continue; // H-H — keep both
        }
        // a plain terminal H → collapse it
        keep[i] = false;
        extra_h[heavy] = extra_h[heavy].saturating_add(1);
    }

    let survivors: Vec<usize> = (0..n).filter(|&i| keep[i]).collect();
    let mut out = mol.subgraph(&survivors);
    // re-apply the collapsed hydrogen counts on the renumbered atoms
    let mut remap = vec![usize::MAX; n];
    for (new, &old) in survivors.iter().enumerate() {
        remap[old] = new;
    }
    for (old, &add) in extra_h.iter().enumerate() {
        if add > 0 && remap[old] != usize::MAX {
            out.atoms[remap[old]].implicit_h += add;
        }
    }
    out
}

/// Re-derive every heavy atom's implicit-hydrogen count from the
/// current bond graph + formal charge. Mutates `mol` in place. Use
/// after a bond edit / a reaction transform so valences stay correct.
pub fn recompute_implicit_hydrogens(mol: &mut Molecule) {
    for i in 0..mol.atoms.len() {
        let a = &mol.atoms[i];
        if a.is_dummy() || a.is_hydrogen() {
            continue;
        }
        if a.explicit_h > 0 {
            // an atom that declared an explicit count keeps it
            mol.atoms[i].implicit_h = 0;
            continue;
        }
        let z = a.atomic_number;
        let valences = element::standard_valences(z);
        if valences.is_empty() {
            mol.atoms[i].implicit_h = 0;
            continue;
        }
        let bond_sum = mol.explicit_valence(i);
        let charge = i32::from(mol.atoms[i].formal_charge);
        let target = pick_target_valence(valences, bond_sum, charge, z);
        let h = (target - bond_sum).round();
        mol.atoms[i].implicit_h = if h > 0.0 { h as u8 } else { 0 };
    }
}

/// Smallest standard valence ≥ the current bond sum, charge-adjusted.
fn pick_target_valence(valences: &[u8], bond_sum: f64, charge: i32, z: u8) -> f64 {
    let offset = match z {
        7 | 15 | 8 | 16 => charge,
        6 => -charge.abs(),
        _ => charge,
    };
    for &v in valences {
        let t = f64::from(v) + f64::from(offset);
        if t >= bond_sum - 1e-6 {
            return t.max(0.0);
        }
    }
    (f64::from(*valences.last().unwrap()) + f64::from(offset)).max(bond_sum)
}

/// Total hydrogen count of the molecule — explicit H nodes plus every
/// heavy atom's implicit + bracket-explicit count.
pub fn total_hydrogen_count(mol: &Molecule) -> u32 {
    let nodes = mol.atoms.iter().filter(|a| a.is_hydrogen()).count() as u32;
    let attached: u32 = mol
        .atoms
        .iter()
        .filter(|a| !a.is_hydrogen() && !a.is_dummy())
        .map(|a| u32::from(a.total_h()))
        .sum();
    nodes + attached
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::BondOrder;
    use crate::smiles::parse_smiles;

    #[test]
    fn add_then_remove_round_trips() {
        let m = parse_smiles("CCO").unwrap();
        let heavy = m.heavy_atom_count();
        let h_total = total_hydrogen_count(&m);

        let expanded = add_explicit_hydrogens(&m);
        // ethanol C2H6O → 6 explicit H nodes
        assert_eq!(expanded.atoms.iter().filter(|a| a.is_hydrogen()).count(), 6);
        assert_eq!(expanded.atom_count(), heavy + 6);
        assert!(expanded
            .atoms
            .iter()
            .all(|a| a.implicit_h == 0 || a.is_hydrogen()));

        let collapsed = remove_explicit_hydrogens(&expanded);
        assert_eq!(collapsed.heavy_atom_count(), heavy);
        assert_eq!(total_hydrogen_count(&collapsed), h_total);
    }

    #[test]
    fn methane_expands_to_four_h() {
        let m = parse_smiles("C").unwrap();
        let e = add_explicit_hydrogens(&m);
        assert_eq!(e.atom_count(), 5);
        assert_eq!(e.degree(0), 4);
    }

    #[test]
    fn recompute_after_edit() {
        let mut m = parse_smiles("CC").unwrap(); // ethane, each C has 3 H
        assert_eq!(m.atoms[0].implicit_h, 3);
        // upgrade the C-C bond to a double bond → ethene, 2 H each
        m.bonds[0].order = BondOrder::Double;
        recompute_implicit_hydrogens(&mut m);
        assert_eq!(m.atoms[0].implicit_h, 2);
        assert_eq!(m.atoms[1].implicit_h, 2);
    }

    #[test]
    fn labelled_hydrogen_is_kept() {
        // a deuterium on a carbon must survive remove_explicit_hydrogens
        let mut m = parse_smiles("C").unwrap();
        let mut d = Atom::new(1);
        d.isotope = Some(2);
        let di = m.add_atom(d);
        m.add_bond(Bond::single(0, di));
        let collapsed = remove_explicit_hydrogens(&m);
        assert!(collapsed.atoms.iter().any(|a| a.isotope == Some(2)));
    }
}
