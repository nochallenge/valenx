//! Feature 2 — ligand torsion-tree extraction.
//!
//! A flexible ligand docks as a *torsion tree*: one rigid root group
//! of atoms, plus a set of branches each of which is a rigid group of
//! atoms rotatable about a single bond relative to its parent. This is
//! exactly the structure AutoDockTools writes into a PDBQT
//! `ROOT` / `BRANCH` block and that [`valenx_dock::ligand::Ligand`]
//! parses.
//!
//! This module does the *perception* side: given a molecular graph
//! (atoms + bonds), it finds which bonds are rotatable and assembles
//! them into a [`TorsionTree`]. The rotatable-bond rule is the
//! standard one used by AutoDockTools, RDKit and Open Babel:
//!
//! A bond is **rotatable** when it is
//!
//! - a single (σ) bond,
//! - not in a ring,
//! - not terminal (each endpoint has a heavy-atom degree ≥ 2), and
//! - not an amide C–N bond (the partial double-bond character of a
//!   peptide bond makes it effectively rigid).

use std::collections::{BTreeSet, VecDeque};

use valenx_cheminf::molecule::{BondOrder, Molecule};
use valenx_cheminf::perceive::rings::sssr;

use crate::error::{DockScreenError, Result};

/// A single rotatable bond, identified by its two endpoint atom
/// indices into the source molecule.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RotatableBond {
    /// The bond's index in [`Molecule::bonds`].
    pub bond_index: usize,
    /// One endpoint atom index.
    pub atom_a: usize,
    /// The other endpoint atom index.
    pub atom_b: usize,
}

/// A rigid group of atoms in a torsion tree: the root, or a branch
/// rotatable about a single bond relative to its parent group.
#[derive(Clone, Debug, PartialEq)]
pub struct TorsionGroup {
    /// Atom indices (into the source molecule) that belong rigidly to
    /// this group.
    pub atoms: Vec<usize>,
    /// `(parent_atom, child_atom)` defining the rotatable bond this
    /// group turns about. `None` for the root group.
    pub axis: Option<(usize, usize)>,
    /// Child-group indices into [`TorsionTree::groups`].
    pub children: Vec<usize>,
}

/// A ligand's torsion tree: a rigid root plus rotatable branches.
#[derive(Clone, Debug, PartialEq)]
pub struct TorsionTree {
    /// Every rigid group; `root` indexes the root.
    pub groups: Vec<TorsionGroup>,
    /// Index of the root group in [`TorsionTree::groups`].
    pub root: usize,
    /// Every rotatable bond, in discovery order. `rotatable.len()`
    /// equals the ligand's torsional degrees of freedom (`TORSDOF`).
    pub rotatable: Vec<RotatableBond>,
}

impl TorsionTree {
    /// Number of rotatable bonds — the ligand's torsional DOF.
    pub fn torsion_count(&self) -> usize {
        self.rotatable.len()
    }

    /// Number of rigid groups (root + branches).
    pub fn group_count(&self) -> usize {
        self.groups.len()
    }

    /// Extract the torsion tree of a molecule.
    ///
    /// The heaviest connected fragment is taken as the ligand (so a
    /// salt counter-ion does not confuse the tree). The largest rigid
    /// group becomes the root; remaining groups hang off it through
    /// their rotatable bonds.
    ///
    /// Returns [`DockScreenError::InvalidLigand`] if the molecule has
    /// no heavy atoms.
    pub fn from_molecule(mol: &Molecule) -> Result<Self> {
        if mol.heavy_atom_count() == 0 {
            return Err(DockScreenError::invalid_ligand(
                "molecule has no heavy atoms",
            ));
        }
        let rings = sssr(mol);
        // 1. Find the rotatable bonds.
        let mut rotatable: Vec<RotatableBond> = Vec::new();
        for (bi, bond) in mol.bonds.iter().enumerate() {
            if is_rotatable_bond(mol, bi, &rings) {
                rotatable.push(RotatableBond {
                    bond_index: bi,
                    atom_a: bond.a,
                    atom_b: bond.b,
                });
            }
        }
        // 2. Cut the molecule at every rotatable bond — the resulting
        //    connected components are the rigid groups.
        let rotatable_set: BTreeSet<usize> = rotatable.iter().map(|r| r.bond_index).collect();
        let rigid_components = rigid_components(mol, &rotatable_set);
        // Component id per atom; usize::MAX for atoms in no component
        // (cannot happen — every atom is in exactly one component).
        let mut comp_of = vec![usize::MAX; mol.atoms.len()];
        for (cid, atoms) in rigid_components.iter().enumerate() {
            for &a in atoms {
                comp_of[a] = cid;
            }
        }
        // 3. The largest component (by atom count) is the root.
        let root_comp = rigid_components
            .iter()
            .enumerate()
            .max_by_key(|(_, atoms)| atoms.len())
            .map(|(cid, _)| cid)
            .unwrap_or(0);
        // 4. Build the group tree by BFS over the component graph,
        //    where two components are linked by a rotatable bond.
        let mut groups: Vec<TorsionGroup> = rigid_components
            .iter()
            .map(|atoms| TorsionGroup {
                atoms: atoms.clone(),
                axis: None,
                children: Vec::new(),
            })
            .collect();
        let mut visited = vec![false; groups.len()];
        let mut queue: VecDeque<usize> = VecDeque::new();
        visited[root_comp] = true;
        queue.push_back(root_comp);
        let mut ordered_rotatable: Vec<RotatableBond> = Vec::new();
        while let Some(parent) = queue.pop_front() {
            // Find rotatable bonds with one endpoint in `parent`.
            for rb in &rotatable {
                let ca = comp_of[rb.atom_a];
                let cb = comp_of[rb.atom_b];
                let (parent_atom, child_atom, child_comp) = if ca == parent && cb != parent {
                    (rb.atom_a, rb.atom_b, cb)
                } else if cb == parent && ca != parent {
                    (rb.atom_b, rb.atom_a, ca)
                } else {
                    continue;
                };
                if visited[child_comp] {
                    continue;
                }
                visited[child_comp] = true;
                groups[child_comp].axis = Some((parent_atom, child_atom));
                groups[parent].children.push(child_comp);
                ordered_rotatable.push(*rb);
                queue.push_back(child_comp);
            }
        }
        Ok(TorsionTree {
            groups,
            root: root_comp,
            rotatable: ordered_rotatable,
        })
    }

    /// Project the per-group atom membership of an already-parsed
    /// [`valenx_dock::ligand::Ligand`]. The dock crate has already
    /// built the tree from a PDBQT `ROOT` / `BRANCH` block; this just
    /// re-expresses it in this crate's [`TorsionTree`] shape so the
    /// rest of the pipeline has one type.
    pub fn from_dock_ligand(lig: &valenx_dock::ligand::Ligand) -> Self {
        let groups: Vec<TorsionGroup> = lig
            .groups
            .iter()
            .map(|g| {
                let axis = if g.axis.0 == usize::MAX || g.axis.1 == usize::MAX {
                    None
                } else {
                    Some(g.axis)
                };
                TorsionGroup {
                    atoms: g.atom_indices.clone(),
                    axis,
                    children: g.children.clone(),
                }
            })
            .collect();
        let rotatable: Vec<RotatableBond> = lig
            .groups
            .iter()
            .enumerate()
            .filter_map(|(gi, g)| {
                if g.axis.0 == usize::MAX || g.axis.1 == usize::MAX {
                    None
                } else {
                    Some(RotatableBond {
                        bond_index: gi, // group index doubles as a stable id here
                        atom_a: g.axis.0,
                        atom_b: g.axis.1,
                    })
                }
            })
            .collect();
        TorsionTree {
            groups,
            root: lig.root_group,
            rotatable,
        }
    }
}

/// Count the rotatable bonds of a molecule without building the full
/// tree — a cheap descriptor used by screening filters.
pub fn rotatable_bond_count(mol: &Molecule) -> usize {
    let rings = sssr(mol);
    (0..mol.bonds.len())
        .filter(|&bi| is_rotatable_bond(mol, bi, &rings))
        .count()
}

/// The rotatable-bond predicate. See the module docs for the rule.
fn is_rotatable_bond(
    mol: &Molecule,
    bond_index: usize,
    rings: &valenx_cheminf::perceive::rings::RingInfo,
) -> bool {
    let bond = &mol.bonds[bond_index];
    // Single (σ) bonds only.
    if bond.order != BondOrder::Single || bond.aromatic {
        return false;
    }
    // Ring bonds are not rotatable (the ring constrains them).
    if rings.bond_in_ring(bond_index) {
        return false;
    }
    let a = bond.a;
    let b = bond.b;
    // Skip bonds to hydrogen.
    if mol.atoms[a].is_hydrogen() || mol.atoms[b].is_hydrogen() {
        return false;
    }
    // Terminal bonds (one endpoint has only this heavy neighbour) are
    // not rotatable — rotating a terminal group changes nothing.
    if heavy_degree(mol, a) < 2 || heavy_degree(mol, b) < 2 {
        return false;
    }
    // Amide C–N bond: peptide-bond resonance makes it effectively
    // rigid. Detect a C and an N where the C also carries a =O.
    if is_amide_bond(mol, a, b) {
        return false;
    }
    true
}

/// Heavy-atom (non-hydrogen) degree of an atom.
fn heavy_degree(mol: &Molecule, atom: usize) -> usize {
    mol.neighbors(atom)
        .into_iter()
        .filter(|&n| !mol.atoms[n].is_hydrogen())
        .count()
}

/// `true` if the `a`–`b` bond is an amide C–N bond: one endpoint a
/// carbon double-bonded to an oxygen, the other a nitrogen.
fn is_amide_bond(mol: &Molecule, a: usize, b: usize) -> bool {
    let amide = |c: usize, n: usize| -> bool {
        mol.atoms[c].atomic_number == 6
            && mol.atoms[n].atomic_number == 7
            && mol.bonds.iter().any(|bd| {
                bd.touches(c)
                    && bd.order == BondOrder::Double
                    && bd.other(c).is_some_and(|o| mol.atoms[o].atomic_number == 8)
            })
    };
    amide(a, b) || amide(b, a)
}

/// Connected components of the molecule after deleting every bond in
/// `cut_bonds`. Each component is a sorted atom-index list.
fn rigid_components(mol: &Molecule, cut_bonds: &BTreeSet<usize>) -> Vec<Vec<usize>> {
    let n = mol.atoms.len();
    // Adjacency restricted to non-cut bonds.
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (bi, bond) in mol.bonds.iter().enumerate() {
        if cut_bonds.contains(&bi) {
            continue;
        }
        adj[bond.a].push(bond.b);
        adj[bond.b].push(bond.a);
    }
    let mut comp = vec![usize::MAX; n];
    let mut next = 0usize;
    for start in 0..n {
        if comp[start] != usize::MAX {
            continue;
        }
        let mut stack = vec![start];
        comp[start] = next;
        while let Some(u) = stack.pop() {
            for &v in &adj[u] {
                if comp[v] == usize::MAX {
                    comp[v] = next;
                    stack.push(v);
                }
            }
        }
        next += 1;
    }
    let mut out: Vec<Vec<usize>> = vec![Vec::new(); next];
    for (atom, &c) in comp.iter().enumerate() {
        out[c].push(atom);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cheminf::mol_from_smiles;

    #[test]
    fn ethane_has_no_rotatable_bond() {
        // CC — the single C-C bond is terminal on both ends.
        let m = mol_from_smiles("CC").unwrap();
        assert_eq!(rotatable_bond_count(&m), 0);
    }

    #[test]
    fn butane_has_one_rotatable_bond() {
        // CCCC — only the central C-C bond is non-terminal.
        let m = mol_from_smiles("CCCC").unwrap();
        assert_eq!(rotatable_bond_count(&m), 1);
    }

    #[test]
    fn benzene_ring_bonds_are_not_rotatable() {
        let m = mol_from_smiles("c1ccccc1").unwrap();
        assert_eq!(rotatable_bond_count(&m), 0);
    }

    #[test]
    fn biphenyl_inter_ring_bond_is_rotatable() {
        // Two phenyls joined by one single bond — that bond rotates.
        let m = mol_from_smiles("c1ccccc1-c1ccccc1").unwrap();
        assert_eq!(rotatable_bond_count(&m), 1);
    }

    #[test]
    fn amide_bond_is_treated_as_rigid() {
        // N-methylacetamide CC(=O)NC — the C(=O)-N bond is amide.
        // The remaining single bonds: methyl C-C(=O) is terminal,
        // N-C(methyl) is terminal — so zero rotatable bonds total.
        let m = mol_from_smiles("CC(=O)NC").unwrap();
        assert_eq!(rotatable_bond_count(&m), 0);
    }

    #[test]
    fn torsion_tree_root_plus_one_branch_for_butane() {
        let m = mol_from_smiles("CCCC").unwrap();
        let tree = TorsionTree::from_molecule(&m).unwrap();
        assert_eq!(tree.torsion_count(), 1);
        // One rotatable bond → exactly two rigid groups.
        assert_eq!(tree.group_count(), 2);
        // The root has one child branch.
        assert_eq!(tree.groups[tree.root].children.len(), 1);
        // Every heavy atom appears in exactly one group.
        let mut all: Vec<usize> = tree
            .groups
            .iter()
            .flat_map(|g| g.atoms.iter().copied())
            .collect();
        all.sort();
        all.dedup();
        assert_eq!(all.len(), m.atom_count());
    }

    #[test]
    fn torsion_tree_rejects_empty_molecule() {
        let m = Molecule::new();
        assert!(TorsionTree::from_molecule(&m).is_err());
    }

    #[test]
    fn torsion_tree_branch_axis_links_parent_to_child() {
        let m = mol_from_smiles("CCCC").unwrap();
        let tree = TorsionTree::from_molecule(&m).unwrap();
        let branch = &tree.groups[tree.groups[tree.root].children[0]];
        let (pa, ca) = branch.axis.expect("branch must carry an axis");
        // The parent atom is in the root group, the child atom in the
        // branch group.
        assert!(tree.groups[tree.root].atoms.contains(&pa));
        assert!(branch.atoms.contains(&ca));
    }

    #[test]
    fn from_dock_ligand_roundtrips_group_membership() {
        let pdbqt = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ATOM      2  C2  LIG A   1       1.500   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
BRANCH   2   3
ATOM      3  C3  LIG A   1       3.000   0.000   0.000  1.00  0.00     0.000 C
ENDBRANCH   2   3
TORSDOF 1
";
        let lig = valenx_dock::ligand::Ligand::from_pdbqt(pdbqt).unwrap();
        let tree = TorsionTree::from_dock_ligand(&lig);
        assert_eq!(tree.torsion_count(), 1);
        assert_eq!(tree.group_count(), 2);
        // Root has atoms 0,1; branch has atom 2.
        assert_eq!(tree.groups[tree.root].atoms, vec![0, 1]);
    }
}
