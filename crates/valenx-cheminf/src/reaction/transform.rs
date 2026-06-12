//! SMIRKS-class reaction transforms and combinatorial enumeration.
//!
//! A [`ReactionTransform`] is a reaction template `pattern >> product`
//! where both halves are SMARTS-style graphs and atom-map numbers
//! (`[C:1]`) link a reactant atom to the product atom it becomes.
//! [`ReactionTransform::apply`] runs the template on a molecule:
//!
//! 1. find every substructure match of the reactant pattern (VF2);
//! 2. for each match, build the product by copying the input, then
//!    rewiring the mapped atoms / bonds to follow the product
//!    template, dropping the bonds the template breaks and adding the
//!    bonds it forms;
//! 3. recompute implicit hydrogens so valences stay correct.
//!
//! [`enumerate_library`] takes one transform and several reactant
//! lists and produces the full combinatorial product set — a one-pot
//! virtual library.
//!
//! **v1 scope.** The transform supports single-reactant templates and
//! the common two-reactant case (the pattern is matched against a
//! `.`-joined pair). Bond-order changes and bond make/break between
//! mapped atoms are applied; stereochemistry and the creation of brand
//! new unmapped product atoms beyond a small fixed set are not modeled
//! — those need a richer SMIRKS engine and are documented as out of
//! v1 scope.

use crate::error::{CheminfError, Result};
use crate::molecule::{Atom, Bond, BondOrder, Molecule};
use crate::smarts::SmartsPattern;

/// A reaction template: a reactant SMARTS, a product SMARTS, and the
/// raw SMIRKS string it was parsed from.
#[derive(Clone, Debug)]
pub struct ReactionTransform {
    /// The reactant-side query pattern.
    pub reactant: SmartsPattern,
    /// The product-side template pattern.
    pub product: SmartsPattern,
    /// The original SMIRKS text.
    pub source: String,
}

impl ReactionTransform {
    /// Parse a SMIRKS-class string `reactantSMARTS>>productSMARTS`.
    pub fn parse(smirks: &str) -> Result<Self> {
        let s = smirks.trim();
        let parts: Vec<&str> = s.split(">>").collect();
        if parts.len() != 2 {
            return Err(CheminfError::parse(
                "smarts",
                "SMIRKS must contain exactly one '>>' separator",
            ));
        }
        let reactant = SmartsPattern::parse(parts[0])?;
        let product = SmartsPattern::parse(parts[1])?;
        // every mapped product atom must have a mapped reactant atom
        for pa in &product.atoms {
            if pa.map != 0 && !reactant.atoms.iter().any(|ra| ra.map == pa.map) {
                return Err(CheminfError::parse(
                    "smarts",
                    format!("product atom-map :{} has no reactant counterpart", pa.map),
                ));
            }
        }
        Ok(ReactionTransform {
            reactant,
            product,
            source: s.to_string(),
        })
    }

    /// Apply the transform to `mol`, returning every distinct product.
    ///
    /// One product per reactant-pattern match; duplicate products
    /// (from symmetry-equivalent matches) are de-duplicated by
    /// canonical SMILES.
    pub fn apply(&self, mol: &Molecule) -> Vec<Molecule> {
        let matches = self.reactant.find_matches(mol);
        let mut products: Vec<Molecule> = Vec::new();
        let mut seen: Vec<String> = Vec::new();
        for m in matches {
            if let Some(prod) = self.apply_match(mol, &m) {
                let key = crate::smiles::write_canonical_smiles(&prod);
                if !seen.contains(&key) {
                    seen.push(key);
                    products.push(prod);
                }
            }
        }
        products
    }

    /// Apply the transform to one specific match.
    fn apply_match(&self, mol: &Molecule, mapping: &[usize]) -> Option<Molecule> {
        let mut out = mol.clone();

        // map: atom-map number → target atom index
        let mut map_to_target: std::collections::HashMap<u32, usize> =
            std::collections::HashMap::new();
        for (qi, &ti) in mapping.iter().enumerate() {
            let m = self.reactant.atoms[qi].map;
            if m != 0 {
                map_to_target.insert(m, ti);
            }
        }

        // 1. remove bonds present between mapped reactant atoms but
        //    absent in the product template
        let reactant_bonds = mapped_bond_set(&self.reactant);
        let product_bonds = mapped_bond_set(&self.product);
        let mut bonds_to_remove: Vec<(usize, usize)> = Vec::new();
        for &(m1, m2, _) in &reactant_bonds {
            if !product_bonds
                .iter()
                .any(|&(p1, p2, _)| (p1 == m1 && p2 == m2) || (p1 == m2 && p2 == m1))
            {
                if let (Some(&t1), Some(&t2)) = (map_to_target.get(&m1), map_to_target.get(&m2)) {
                    bonds_to_remove.push((t1, t2));
                }
            }
        }

        // 2. apply product bonds: change order or add new bond
        for &(m1, m2, order) in &product_bonds {
            if let (Some(&t1), Some(&t2)) = (map_to_target.get(&m1), map_to_target.get(&m2)) {
                match out.bond_between(t1, t2) {
                    Some(bi) => {
                        out.bonds[bi].order = order;
                        out.bonds[bi].aromatic = order == BondOrder::Aromatic;
                    }
                    None => {
                        out.add_bond(Bond {
                            a: t1,
                            b: t2,
                            order,
                            aromatic: order == BondOrder::Aromatic,
                            stereo: crate::molecule::BondStereo::None,
                        });
                    }
                }
            }
        }

        // 3. apply per-atom changes from the product template
        //    (element / charge changes on a mapped atom)
        for (pi, pa) in self.product.atoms.iter().enumerate() {
            if pa.map == 0 {
                continue;
            }
            if let Some(&ti) = map_to_target.get(&pa.map) {
                apply_atom_template(&mut out.atoms[ti], &self.product, pi);
            }
        }

        // 4. remove the broken bonds (after additions, so indices we
        //    captured for removal are still valid)
        remove_bonds(&mut out, &bonds_to_remove);

        crate::perceive::hydrogen::recompute_implicit_hydrogens(&mut out);
        if out.validate().is_err() {
            return None;
        }
        Some(out)
    }
}

/// The mapped bonds of a SMARTS pattern as `(map_a, map_b, order)`.
fn mapped_bond_set(pat: &SmartsPattern) -> Vec<(u32, u32, BondOrder)> {
    let mut out = Vec::new();
    for qb in &pat.bonds {
        let ma = pat.atoms[qb.a].map;
        let mb = pat.atoms[qb.b].map;
        if ma != 0 && mb != 0 {
            out.push((ma, mb, bond_expr_order(qb.expr)));
        }
    }
    out
}

/// Map a SMARTS [`crate::smarts::BondExpr`] to a concrete bond order
/// for use on the product side.
fn bond_expr_order(expr: crate::smarts::BondExpr) -> BondOrder {
    use crate::smarts::BondExpr::*;
    match expr {
        Double => BondOrder::Double,
        Triple => BondOrder::Triple,
        Aromatic => BondOrder::Aromatic,
        _ => BondOrder::Single,
    }
}

/// Copy the element / charge a product template atom prescribes onto a
/// concrete output atom.
fn apply_atom_template(atom: &mut Atom, product: &SmartsPattern, pi: usize) {
    use crate::smarts::AtomExpr;
    // walk the (possibly And-wrapped) expression for an Element / Charge
    fn visit(expr: &AtomExpr, atom: &mut Atom) {
        match expr {
            AtomExpr::Element { z, aromatic } => {
                if *z != 0 {
                    atom.atomic_number = *z;
                }
                if let Some(arom) = aromatic {
                    atom.aromatic = *arom;
                }
            }
            AtomExpr::Charge(c) => {
                atom.formal_charge = (*c).clamp(-9, 9) as i8;
            }
            AtomExpr::And(subs) => {
                for s in subs {
                    visit(s, atom);
                }
            }
            _ => {}
        }
    }
    visit(&product.atoms[pi].expr, atom);
}

/// Rebuild `mol` without the bonds in `pairs`.
fn remove_bonds(mol: &mut Molecule, pairs: &[(usize, usize)]) {
    if pairs.is_empty() {
        return;
    }
    mol.bonds.retain(|b| {
        !pairs
            .iter()
            .any(|&(a, c)| (b.a == a && b.b == c) || (b.a == c && b.b == a))
    });
}

/// Combinatorially enumerate a virtual library: apply `transform` to
/// every combination drawn one-from-each of `reactant_sets`.
///
/// For a single-reactant transform pass one set. For a two-reactant
/// transform pass two sets; each `(r1, r2)` pair is joined into a
/// single disconnected molecule before the template is applied, so a
/// bond-forming template can stitch them.
pub fn enumerate_library(
    transform: &ReactionTransform,
    reactant_sets: &[Vec<Molecule>],
) -> Vec<Molecule> {
    if reactant_sets.is_empty() {
        return Vec::new();
    }
    // build the cartesian product of indices
    let mut combos: Vec<Vec<usize>> = vec![Vec::new()];
    for set in reactant_sets {
        let mut next = Vec::new();
        for combo in &combos {
            for i in 0..set.len() {
                let mut c = combo.clone();
                c.push(i);
                next.push(c);
            }
        }
        combos = next;
    }

    // With two or more reactant sets the template is a coupling: a
    // productive match bridges the merged fragments, so a *connected*
    // product is required. (A free substructure match of a
    // multi-component pattern could otherwise place every component
    // inside one reactant, leaving the others untouched and the product
    // disconnected — not a real library member.) A single-set
    // enumeration keeps every product, since a one-reactant transform
    // may legitimately fragment its input.
    let require_connected = reactant_sets.len() >= 2;

    let mut products: Vec<Molecule> = Vec::new();
    let mut seen: Vec<String> = Vec::new();
    for combo in combos {
        // join the chosen reactants into one (possibly disconnected)
        // molecule so a cross-reactant template can match
        let mut merged = Molecule::new();
        for (set_idx, &mol_idx) in combo.iter().enumerate() {
            merge_into(&mut merged, &reactant_sets[set_idx][mol_idx]);
        }
        for prod in transform.apply(&merged) {
            if require_connected && !prod.is_connected() {
                continue;
            }
            let key = crate::smiles::write_canonical_smiles(&prod);
            if !seen.contains(&key) {
                seen.push(key);
                products.push(prod);
            }
        }
    }
    products
}

/// Append `src`'s atoms and bonds into `dst`, offsetting bond indices.
fn merge_into(dst: &mut Molecule, src: &Molecule) {
    let offset = dst.atoms.len();
    for a in &src.atoms {
        dst.atoms.push(a.clone());
    }
    for b in &src.bonds {
        dst.bonds.push(Bond {
            a: b.a + offset,
            b: b.b + offset,
            ..b.clone()
        });
    }
    if !src.coords.is_empty() && dst.coords.len() == offset {
        dst.coords.extend_from_slice(&src.coords);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_from_smiles;
    use crate::smiles::write_canonical_smiles;

    #[test]
    fn parse_smirks() {
        let t = ReactionTransform::parse("[C:1]=[C:2]>>[C:1]-[C:2]").unwrap();
        assert_eq!(t.reactant.atom_count(), 2);
        assert_eq!(t.product.atom_count(), 2);
    }

    #[test]
    fn rejects_bad_smirks() {
        assert!(ReactionTransform::parse("CC>CO").is_err()); // no '>>'
                                                             // product map with no reactant counterpart
        assert!(ReactionTransform::parse("[C:1]>>[C:1][O:2]").is_err());
    }

    #[test]
    fn hydrogenation_transform() {
        // reduce a C=C to a C-C
        let t = ReactionTransform::parse("[C:1]=[C:2]>>[C:1]-[C:2]").unwrap();
        let ethene = mol_from_smiles("C=C").unwrap();
        let products = t.apply(&ethene);
        assert_eq!(products.len(), 1);
        // product is ethane
        let p = &products[0];
        assert!(p.bonds.iter().all(|b| b.order == BondOrder::Single));
        // each carbon now carries 3 H (ethane)
        assert!(p.atoms.iter().all(|a| a.implicit_h == 3));
    }

    #[test]
    fn oxidation_changes_bond_order() {
        // oxidise an alcohol C-O to a carbonyl C=O
        let t = ReactionTransform::parse("[C:1]-[O:2]>>[C:1]=[O:2]").unwrap();
        let methanol = mol_from_smiles("CO").unwrap();
        let products = t.apply(&methanol);
        assert!(!products.is_empty());
        let p = &products[0];
        assert!(p.bonds.iter().any(|b| b.order == BondOrder::Double));
    }

    #[test]
    fn no_match_yields_no_products() {
        let t = ReactionTransform::parse("[C:1]#[C:2]>>[C:1]=[C:2]").unwrap();
        let ethane = mol_from_smiles("CC").unwrap();
        assert!(t.apply(&ethane).is_empty());
    }

    #[test]
    fn library_enumeration() {
        // a bond-forming amide-style coupling: join a C-mapped atom of
        // one reactant to an O-mapped atom of another
        let t = ReactionTransform::parse("[C:1].[O:2]>>[C:1]-[O:2]").unwrap();
        let acids = vec![mol_from_smiles("C").unwrap()];
        let alcohols = vec![
            mol_from_smiles("O").unwrap(),
            mol_from_smiles("CO").unwrap(),
        ];
        let lib = enumerate_library(&t, &[acids, alcohols]);
        // two alcohol partners → up to two distinct products
        assert!(!lib.is_empty());
        for m in &lib {
            assert!(
                m.is_connected(),
                "product {} disconnected",
                write_canonical_smiles(m)
            );
        }
    }
}
