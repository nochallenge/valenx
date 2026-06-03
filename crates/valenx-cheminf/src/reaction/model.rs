//! Reaction data model and reaction-SMILES I/O.
//!
//! A [`Reaction`] is three lists of [`Molecule`]s — reactants, agents
//! and products — exactly the reaction-SMILES layout
//! `reactants > agents > products`. [`parse_reaction_smiles`] reads
//! that string (each role a `.`-separated set of molecule SMILES);
//! [`write_reaction_smiles`] writes it back.

use crate::error::{CheminfError, Result};
use crate::molecule::Molecule;
use crate::smiles::{parse_smiles, write_smiles};

/// A chemical reaction: reactant, agent (catalyst / solvent) and
/// product molecule sets.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Reaction {
    /// Reactant molecules (the `>`-left set).
    pub reactants: Vec<Molecule>,
    /// Agent molecules — catalysts, solvents (the middle set).
    pub agents: Vec<Molecule>,
    /// Product molecules (the `>`-right set).
    pub products: Vec<Molecule>,
}

impl Reaction {
    /// An empty reaction.
    pub fn new() -> Self {
        Reaction::default()
    }

    /// Total atom count across reactants — useful for atom-economy
    /// style checks.
    pub fn reactant_atom_count(&self) -> usize {
        self.reactants.iter().map(|m| m.heavy_atom_count()).sum()
    }

    /// Total atom count across products.
    pub fn product_atom_count(&self) -> usize {
        self.products.iter().map(|m| m.heavy_atom_count()).sum()
    }

    /// `true` if the heavy-atom count is conserved reactants → products
    /// (a necessary, not sufficient, balance check — agents excluded).
    pub fn heavy_atoms_balanced(&self) -> bool {
        self.reactant_atom_count() == self.product_atom_count()
    }

    /// Atom economy — fraction of reactant heavy atoms that end up in
    /// the products (`0.0` if there are no reactant atoms).
    pub fn atom_economy(&self) -> f64 {
        let r = self.reactant_atom_count();
        if r == 0 {
            return 0.0;
        }
        self.product_atom_count() as f64 / r as f64
    }
}

/// Parse a reaction-SMILES string `reactants>agents>products`.
///
/// Each of the three fields is a `.`-separated set of molecule SMILES;
/// any field may be empty (`C>>CO` has no agents).
pub fn parse_reaction_smiles(input: &str) -> Result<Reaction> {
    let s = input.trim();
    let parts: Vec<&str> = s.split('>').collect();
    if parts.len() != 3 {
        return Err(CheminfError::parse(
            "smiles",
            "reaction SMILES must have exactly two '>' separators",
        ));
    }
    let parse_set = |field: &str| -> Result<Vec<Molecule>> {
        if field.trim().is_empty() {
            return Ok(Vec::new());
        }
        field
            .split('.')
            .filter(|p| !p.trim().is_empty())
            .map(parse_smiles)
            .collect()
    };
    Ok(Reaction {
        reactants: parse_set(parts[0])?,
        agents: parse_set(parts[1])?,
        products: parse_set(parts[2])?,
    })
}

/// Write a [`Reaction`] back to a reaction-SMILES string.
pub fn write_reaction_smiles(rxn: &Reaction) -> String {
    let join = |set: &[Molecule]| {
        set.iter()
            .map(write_smiles)
            .collect::<Vec<_>>()
            .join(".")
    };
    format!(
        "{}>{}>{}",
        join(&rxn.reactants),
        join(&rxn.agents),
        join(&rxn.products)
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_simple_reaction() {
        // esterification: acid + alcohol > > ester (+ water)
        let rxn = parse_reaction_smiles("CC(=O)O.CCO>>CC(=O)OCC.O").unwrap();
        assert_eq!(rxn.reactants.len(), 2);
        assert_eq!(rxn.agents.len(), 0);
        assert_eq!(rxn.products.len(), 2);
    }

    #[test]
    fn parse_with_agents() {
        let rxn = parse_reaction_smiles("C=C.[H][H]>[Pd]>CC").unwrap();
        assert_eq!(rxn.reactants.len(), 2);
        assert_eq!(rxn.agents.len(), 1);
        assert_eq!(rxn.products.len(), 1);
    }

    #[test]
    fn round_trips() {
        let rxn = parse_reaction_smiles("CCO>>CC=O").unwrap();
        let s = write_reaction_smiles(&rxn);
        let back = parse_reaction_smiles(&s).unwrap();
        assert_eq!(back.reactants.len(), 1);
        assert_eq!(back.products.len(), 1);
    }

    #[test]
    fn atom_balance_and_economy() {
        // a balanced rearrangement
        let rxn = parse_reaction_smiles("CCO>>CCO").unwrap();
        assert!(rxn.heavy_atoms_balanced());
        assert!((rxn.atom_economy() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_malformed() {
        assert!(parse_reaction_smiles("CCO>CC=O").is_err()); // one '>'
        assert!(parse_reaction_smiles("CCO>>>CC").is_err()); // three '>'
    }
}
