//! Quaternary structure / biological-assembly generation.
//!
//! A crystallographic coordinate file stores one copy (the
//! *asymmetric unit*); the biologically relevant molecule is built by
//! applying the `BIOMT` / symmetry operators from `REMARK 350`
//! (PDB) or `_pdbx_struct_oper_list` (mmCIF).
//!
//! [`generate_assembly`] applies a set of [`SymmetryOperator`]s to a
//! structure, producing a new [`Structure`] with one transformed copy
//! per operator. Each copy's chains are renamed so the assembly has
//! unique chain identifiers.

use crate::error::{BiostructError, Result};
use crate::structure::{Atom, Chain, Model, Residue, Structure, SymmetryOperator};

/// Apply a single symmetry operator to a model, returning a new model
/// with every atom transformed.
pub fn transform_model(model: &Model, op: &SymmetryOperator) -> Model {
    let mut out = Model::new(model.serial);
    for chain in &model.chains {
        let mut new_chain = Chain::new(&chain.id);
        new_chain.seqres = chain.seqres.clone();
        for residue in &chain.residues {
            let mut new_res = Residue::new(&residue.name, residue.seq_num);
            new_res.ins_code = residue.ins_code;
            new_res.hetatm = residue.hetatm;
            for atom in &residue.atoms {
                new_res.atoms.push(Atom {
                    coord: op.apply(&atom.coord),
                    ..atom.clone()
                });
            }
            new_chain.residues.push(new_res);
        }
        out.chains.push(new_chain);
    }
    out
}

/// Generate a biological assembly by applying `operators` to the
/// first model of `structure`.
///
/// The result has a single model containing one transformed copy per
/// operator. Copies after the first get a numeric suffix on each
/// chain id (`A` → `A`, `A_2`, `A_3`, …) so chain ids stay unique.
/// An empty operator list is rejected.
pub fn generate_assembly(
    structure: &Structure,
    operators: &[SymmetryOperator],
) -> Result<Structure> {
    if operators.is_empty() {
        return Err(BiostructError::invalid(
            "operators",
            "assembly generation needs at least one operator",
        ));
    }
    let base_model = structure.first_model();
    let mut assembled = Model::new(1);

    for (copy, op) in operators.iter().enumerate() {
        let transformed = transform_model(base_model, op);
        for mut chain in transformed.chains {
            if copy > 0 {
                chain.id = format!("{}_{}", chain.id, copy + 1);
            }
            assembled.chains.push(chain);
        }
    }

    let mut out = Structure::new(format!("{}_assembly", structure.id));
    out.title = if structure.title.is_empty() {
        "biological assembly".to_string()
    } else {
        format!("{} (biological assembly)", structure.title)
    };
    out.models = vec![assembled];
    Ok(out)
}

/// Generate the assembly declared by the structure's own
/// `assembly_operators` (the parsed `REMARK 350` / mmCIF operators).
///
/// When the structure carries no operators this returns the
/// asymmetric unit unchanged (wrapped as a fresh structure) — the
/// asymmetric unit *is* the biological unit for many entries.
pub fn generate_declared_assembly(structure: &Structure) -> Result<Structure> {
    if structure.assembly_operators.is_empty() {
        return generate_assembly(structure, &[SymmetryOperator::identity()]);
    }
    generate_assembly(structure, &structure.assembly_operators)
}

/// Count the chains a declared assembly would produce without
/// building it — a cheap "is this a monomer / dimer / …" check.
pub fn assembly_chain_count(structure: &Structure) -> usize {
    let ops = structure.assembly_operators.len().max(1);
    structure.first_model().chains.len() * ops
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    fn one_atom_structure() -> Structure {
        let mut s = Structure::new("TEST");
        let mut chain = Chain::new("A");
        let mut r = Residue::new("ALA", 1);
        r.atoms
            .push(Atom::new("CA", "C", Point3::new(1.0, 2.0, 3.0)));
        chain.residues.push(r);
        s.first_model_mut().chains.push(chain);
        s
    }

    #[test]
    fn identity_assembly_preserves_coordinates() {
        let s = one_atom_structure();
        let asm = generate_assembly(&s, &[SymmetryOperator::identity()]).unwrap();
        let a = asm.first_model().atoms().next().unwrap();
        assert!((a.coord - Point3::new(1.0, 2.0, 3.0)).norm() < 1e-12);
        assert_eq!(asm.first_model().chains.len(), 1);
    }

    #[test]
    fn two_operator_assembly_doubles_chains() {
        let s = one_atom_structure();
        let mut shifted = SymmetryOperator::identity();
        shifted.serial = 2;
        shifted.translation = [10.0, 0.0, 0.0];
        let asm = generate_assembly(&s, &[SymmetryOperator::identity(), shifted]).unwrap();
        assert_eq!(asm.first_model().chains.len(), 2);
        // The two chains have distinct ids.
        let ids: Vec<&str> = asm
            .first_model()
            .chains
            .iter()
            .map(|c| c.id.as_str())
            .collect();
        assert_eq!(ids, vec!["A", "A_2"]);
        // The second copy is shifted by +10 in x.
        let second = asm.first_model().chains[1].residues[0].atoms[0].coord;
        assert!((second - Point3::new(11.0, 2.0, 3.0)).norm() < 1e-9);
    }

    #[test]
    fn transform_model_rotates() {
        let s = one_atom_structure();
        // 90-degree rotation about z.
        let mut op = SymmetryOperator::identity();
        op.rotation = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let t = transform_model(s.first_model(), &op);
        let p = t.chains[0].residues[0].atoms[0].coord;
        // (1,2,3) rotated 90 deg about z -> (-2,1,3).
        assert!((p - Point3::new(-2.0, 1.0, 3.0)).norm() < 1e-9);
    }

    #[test]
    fn declared_assembly_uses_structure_operators() {
        let mut s = one_atom_structure();
        let mut op2 = SymmetryOperator::identity();
        op2.serial = 2;
        op2.translation = [5.0, 0.0, 0.0];
        s.assembly_operators = vec![SymmetryOperator::identity(), op2];
        let asm = generate_declared_assembly(&s).unwrap();
        assert_eq!(asm.first_model().chains.len(), 2);
        assert_eq!(assembly_chain_count(&s), 2);
    }

    #[test]
    fn no_operators_falls_back_to_asu() {
        let s = one_atom_structure();
        let asm = generate_declared_assembly(&s).unwrap();
        assert_eq!(asm.first_model().chains.len(), 1);
    }

    #[test]
    fn empty_operator_list_errors() {
        let s = one_atom_structure();
        assert!(generate_assembly(&s, &[]).is_err());
    }
}
