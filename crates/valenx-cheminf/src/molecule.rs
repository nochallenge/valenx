//! The [`Molecule`] graph model — the central data structure every
//! other module reads or builds.
//!
//! A molecule is an undirected graph: an arena of [`Atom`]s and an
//! arena of [`Bond`]s, each addressed by a `usize` index. Indices are
//! stable for the lifetime of the molecule (this crate never deletes an
//! atom in place — operations that drop atoms build a fresh molecule).
//!
//! The model is deliberately RDKit-shaped: each atom carries an
//! element, a formal charge, an optional isotope, an aromatic flag, an
//! explicit-hydrogen count and (after perception) an implicit-hydrogen
//! count and an optional tetrahedral parity; each bond carries a bond
//! order, an aromatic flag and an optional double-bond / wedge stereo
//! marker.

use crate::element::{self, Element};
use serde::{Deserialize, Serialize};

/// A bond order. Aromatic bonds carry [`BondOrder::Aromatic`] *and*
/// have [`Bond::aromatic`] set; a Kekulé form uses `Single` / `Double`
/// with the aromatic flag still set on ring bonds.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondOrder {
    /// Single bond (σ only).
    Single,
    /// Double bond.
    Double,
    /// Triple bond.
    Triple,
    /// Quadruple bond (rare; metal-metal).
    Quadruple,
    /// Delocalised aromatic bond (order ≈ 1.5).
    Aromatic,
}

impl BondOrder {
    /// Numeric bond order used by valence sums: aromatic counts as
    /// `1.5`, the rest are their integer order.
    pub fn as_float(self) -> f64 {
        match self {
            BondOrder::Single => 1.0,
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::Quadruple => 4.0,
            BondOrder::Aromatic => 1.5,
        }
    }

    /// Integer bond order for valence counting where aromatic bonds
    /// have been Kekulé-assigned; aromatic falls back to `1` (the σ
    /// contribution). Use [`as_float`](Self::as_float) for fractional
    /// sums.
    pub fn as_int(self) -> u8 {
        match self {
            BondOrder::Single | BondOrder::Aromatic => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
        }
    }
}

/// Double-bond geometric / wedge stereo marker on a [`Bond`].
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondStereo {
    /// No stereo recorded.
    None,
    /// `E` (trans / opposite) double-bond configuration.
    E,
    /// `Z` (cis / same-side) double-bond configuration.
    Z,
    /// SMILES `/` directional single bond (used to derive E/Z).
    Up,
    /// SMILES `\` directional single bond.
    Down,
}

/// Tetrahedral chirality at an atom — the parity of the neighbour
/// ordering as written, plus (after CIP perception) an R/S label.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Chirality {
    /// No stereo centre / not perceived.
    None,
    /// SMILES `@` — anticlockwise neighbour order.
    Ccw,
    /// SMILES `@@` — clockwise neighbour order.
    Cw,
    /// CIP-assigned `R` configuration.
    R,
    /// CIP-assigned `S` configuration.
    S,
}

/// A single atom in a [`Molecule`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Atom {
    /// Atomic number (1 = H, 6 = C, …). `0` is reserved for a SMARTS
    /// wildcard / dummy atom.
    pub atomic_number: u8,
    /// Formal charge (signed).
    pub formal_charge: i8,
    /// Isotope mass number, or `None` for natural-abundance.
    pub isotope: Option<u16>,
    /// `true` if the atom is part of an aromatic system.
    pub aromatic: bool,
    /// Hydrogens written explicitly in the source notation (the `H3`
    /// in `[CH3]`). Implicit hydrogens go in [`Atom::implicit_h`].
    pub explicit_h: u8,
    /// Hydrogens added by valence-based perception. Zero until
    /// [`crate::perceive`] runs.
    pub implicit_h: u8,
    /// Tetrahedral chirality, if any.
    pub chirality: Chirality,
    /// Optional atom-map number (`[C:1]`), used by reaction SMARTS.
    pub map_number: u32,
}

impl Atom {
    /// A neutral, non-aromatic, natural-abundance atom of element `z`.
    pub fn new(z: u8) -> Self {
        Atom {
            atomic_number: z,
            formal_charge: 0,
            isotope: None,
            aromatic: false,
            explicit_h: 0,
            implicit_h: 0,
            chirality: Chirality::None,
            map_number: 0,
        }
    }

    /// The [`Element`] record, or `None` for the dummy atom (`z == 0`).
    pub fn element(&self) -> Option<&'static Element> {
        element::by_number(self.atomic_number)
    }

    /// IUPAC element symbol; `"*"` for the dummy atom.
    pub fn symbol(&self) -> &'static str {
        self.element().map(|e| e.symbol).unwrap_or("*")
    }

    /// Total hydrogen count — explicit plus perception-added implicit.
    pub fn total_h(&self) -> u8 {
        self.explicit_h + self.implicit_h
    }

    /// `true` for a SMARTS / dummy wildcard atom.
    pub fn is_dummy(&self) -> bool {
        self.atomic_number == 0
    }

    /// `true` if this is a hydrogen atom (an *explicit* H node).
    pub fn is_hydrogen(&self) -> bool {
        self.atomic_number == 1
    }
}

/// A single bond between two atoms in a [`Molecule`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Bond {
    /// Index of the first endpoint atom.
    pub a: usize,
    /// Index of the second endpoint atom.
    pub b: usize,
    /// Bond order.
    pub order: BondOrder,
    /// `true` if the bond is part of an aromatic ring system.
    pub aromatic: bool,
    /// Double-bond / directional stereo marker.
    pub stereo: BondStereo,
}

impl Bond {
    /// A plain single bond between `a` and `b`.
    pub fn single(a: usize, b: usize) -> Self {
        Bond {
            a,
            b,
            order: BondOrder::Single,
            aromatic: false,
            stereo: BondStereo::None,
        }
    }

    /// The endpoint of this bond that is *not* `atom`, or `None` if
    /// `atom` is not an endpoint.
    pub fn other(&self, atom: usize) -> Option<usize> {
        if self.a == atom {
            Some(self.b)
        } else if self.b == atom {
            Some(self.a)
        } else {
            None
        }
    }

    /// `true` if `atom` is one of this bond's endpoints.
    pub fn touches(&self, atom: usize) -> bool {
        self.a == atom || self.b == atom
    }
}

/// A molecular graph: atoms, bonds, an adjacency index and an optional
/// name / set of free-form properties (carried through SDF I/O).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Molecule {
    /// Atom arena, indexed by atom id.
    pub atoms: Vec<Atom>,
    /// Bond arena, indexed by bond id.
    pub bonds: Vec<Bond>,
    /// Optional molecule title (SDF record name).
    pub name: String,
    /// Free-form SDF `<property>` key/value pairs, insertion-ordered.
    pub properties: Vec<(String, String)>,
    /// Per-atom 2D or 3D coordinates, parallel to [`Molecule::atoms`].
    /// Empty until a layout / conformer routine fills it; the `z`
    /// component is `0.0` for a 2D depiction.
    pub coords: Vec<[f64; 3]>,
    /// `true` once the coordinates above are genuine 3D (a conformer),
    /// `false` for a 2D depiction or no coordinates.
    pub coords_3d: bool,
}

impl Molecule {
    /// An empty molecule.
    pub fn new() -> Self {
        Molecule::default()
    }

    /// Number of atoms (including any explicit hydrogens present as
    /// nodes).
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Number of bonds.
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// `true` if the molecule has no atoms.
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Append an atom, returning its index.
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        self.atoms.push(atom);
        self.atoms.len() - 1
    }

    /// Append a bond, returning its index. Does **not** validate the
    /// endpoint indices — call [`Molecule::validate`] afterwards.
    pub fn add_bond(&mut self, bond: Bond) -> usize {
        self.bonds.push(bond);
        self.bonds.len() - 1
    }

    /// Convenience: add a bond of the given order between two atoms.
    pub fn bond(&mut self, a: usize, b: usize, order: BondOrder) -> usize {
        self.add_bond(Bond {
            a,
            b,
            order,
            aromatic: order == BondOrder::Aromatic,
            stereo: BondStereo::None,
        })
    }

    /// Indices of every bond incident on `atom`.
    pub fn bonds_on(&self, atom: usize) -> Vec<usize> {
        self.bonds
            .iter()
            .enumerate()
            .filter(|(_, b)| b.touches(atom))
            .map(|(i, _)| i)
            .collect()
    }

    /// Indices of every atom directly bonded to `atom`.
    pub fn neighbors(&self, atom: usize) -> Vec<usize> {
        self.bonds
            .iter()
            .filter_map(|b| b.other(atom))
            .collect()
    }

    /// Degree of `atom` — the number of incident bonds (heavy + the
    /// explicit-H nodes; *not* counting implicit hydrogens).
    pub fn degree(&self, atom: usize) -> usize {
        self.bonds.iter().filter(|b| b.touches(atom)).count()
    }

    /// The bond joining `a` and `b`, if one exists.
    pub fn bond_between(&self, a: usize, b: usize) -> Option<usize> {
        self.bonds.iter().position(|bd| {
            (bd.a == a && bd.b == b) || (bd.a == b && bd.b == a)
        })
    }

    /// Sum of bond orders around `atom` (aromatic bonds count `1.5`),
    /// **not** including hydrogens.
    pub fn explicit_valence(&self, atom: usize) -> f64 {
        self.bonds
            .iter()
            .filter(|b| b.touches(atom))
            .map(|b| b.order.as_float())
            .sum()
    }

    /// Total valence of `atom` = explicit bond-order sum + every
    /// hydrogen (explicit-H nodes already count in the bond sum, so we
    /// add only [`Atom::implicit_h`] and the bracket-`explicit_h`).
    pub fn total_valence(&self, atom: usize) -> f64 {
        let a = &self.atoms[atom];
        self.explicit_valence(atom) + f64::from(a.implicit_h) + f64::from(a.explicit_h)
    }

    /// Indices of every heavy (non-hydrogen) atom.
    pub fn heavy_atoms(&self) -> Vec<usize> {
        self.atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| !a.is_hydrogen())
            .map(|(i, _)| i)
            .collect()
    }

    /// Number of heavy (non-hydrogen) atoms.
    pub fn heavy_atom_count(&self) -> usize {
        self.atoms.iter().filter(|a| !a.is_hydrogen()).count()
    }

    /// Connected components, returned as a `component_id` per atom.
    ///
    /// Two atoms share a component iff a bond path joins them. The ids
    /// are dense `0..n_components` in atom order.
    pub fn components(&self) -> Vec<usize> {
        let n = self.atoms.len();
        let mut comp = vec![usize::MAX; n];
        let mut next = 0usize;
        for start in 0..n {
            if comp[start] != usize::MAX {
                continue;
            }
            let mut stack = vec![start];
            comp[start] = next;
            while let Some(u) = stack.pop() {
                for v in self.neighbors(u) {
                    if comp[v] == usize::MAX {
                        comp[v] = next;
                        stack.push(v);
                    }
                }
            }
            next += 1;
        }
        comp
    }

    /// Number of connected components (`0` for an empty molecule).
    pub fn component_count(&self) -> usize {
        self.components().into_iter().max().map_or(0, |m| m + 1)
    }

    /// `true` if every atom is reachable from atom 0 — i.e. the
    /// molecule is a single connected fragment.
    pub fn is_connected(&self) -> bool {
        self.component_count() <= 1
    }

    /// Validate structural invariants: every bond endpoint in range and
    /// no atom bonded to itself. Returns the molecule's own error type.
    pub fn validate(&self) -> crate::error::Result<()> {
        let n = self.atoms.len();
        for (i, b) in self.bonds.iter().enumerate() {
            if b.a >= n || b.b >= n {
                return Err(crate::error::CheminfError::invalid_molecule(format!(
                    "bond {i} references atom out of range (n_atoms = {n})"
                )));
            }
            if b.a == b.b {
                return Err(crate::error::CheminfError::invalid_molecule(format!(
                    "bond {i} is a self-loop on atom {}",
                    b.a
                )));
            }
        }
        if !self.coords.is_empty() && self.coords.len() != n {
            return Err(crate::error::CheminfError::invalid_molecule(format!(
                "coords length {} != atom count {n}",
                self.coords.len()
            )));
        }
        Ok(())
    }

    /// Build a new molecule keeping only the atoms whose index is in
    /// `keep`, renumbering atoms and bonds densely. Bonds with an
    /// endpoint outside `keep` are dropped. Coordinates and per-atom
    /// data follow their atoms. Used by salt-stripping, scaffold
    /// extraction and MCS.
    pub fn subgraph(&self, keep: &[usize]) -> Molecule {
        let mut remap = vec![usize::MAX; self.atoms.len()];
        let mut out = Molecule {
            name: self.name.clone(),
            properties: self.properties.clone(),
            coords_3d: self.coords_3d,
            ..Molecule::default()
        };
        for &old in keep {
            if old < self.atoms.len() && remap[old] == usize::MAX {
                remap[old] = out.atoms.len();
                out.atoms.push(self.atoms[old].clone());
                if !self.coords.is_empty() {
                    out.coords.push(self.coords[old]);
                }
            }
        }
        for b in &self.bonds {
            if b.a < remap.len() && b.b < remap.len() {
                let (na, nb) = (remap[b.a], remap[b.b]);
                if na != usize::MAX && nb != usize::MAX {
                    out.bonds.push(Bond {
                        a: na,
                        b: nb,
                        ..b.clone()
                    });
                }
            }
        }
        out
    }

    /// The largest connected fragment (by heavy-atom count), as a fresh
    /// molecule. Returns an empty molecule when `self` is empty. The
    /// workhorse behind salt stripping.
    pub fn largest_fragment(&self) -> Molecule {
        if self.is_empty() {
            return Molecule::new();
        }
        let comp = self.components();
        let n_comp = comp.iter().copied().max().map_or(0, |m| m + 1);
        let mut heavy = vec![0usize; n_comp];
        for (i, &c) in comp.iter().enumerate() {
            if !self.atoms[i].is_hydrogen() {
                heavy[c] += 1;
            }
        }
        let best = heavy
            .iter()
            .enumerate()
            .max_by_key(|(_, &h)| h)
            .map(|(c, _)| c)
            .unwrap_or(0);
        let keep: Vec<usize> = comp
            .iter()
            .enumerate()
            .filter(|(_, &c)| c == best)
            .map(|(i, _)| i)
            .collect();
        self.subgraph(&keep)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Ethanol skeleton C-C-O (heavy atoms only).
    fn ethanol() -> Molecule {
        let mut m = Molecule::new();
        let c0 = m.add_atom(Atom::new(6));
        let c1 = m.add_atom(Atom::new(6));
        let o = m.add_atom(Atom::new(8));
        m.bond(c0, c1, BondOrder::Single);
        m.bond(c1, o, BondOrder::Single);
        m
    }

    #[test]
    fn build_and_query() {
        let m = ethanol();
        assert_eq!(m.atom_count(), 3);
        assert_eq!(m.bond_count(), 2);
        assert_eq!(m.degree(1), 2);
        assert_eq!(m.neighbors(1), vec![0, 2]);
        assert_eq!(m.bond_between(0, 1), Some(0));
        assert_eq!(m.bond_between(0, 2), None);
        assert!(m.validate().is_ok());
        assert!(m.is_connected());
    }

    #[test]
    fn valence_sums() {
        let m = ethanol();
        assert_eq!(m.explicit_valence(1), 2.0);
        assert_eq!(m.explicit_valence(0), 1.0);
    }

    #[test]
    fn bond_order_arithmetic() {
        assert_eq!(BondOrder::Aromatic.as_float(), 1.5);
        assert_eq!(BondOrder::Triple.as_int(), 3);
        assert_eq!(BondOrder::Aromatic.as_int(), 1);
    }

    #[test]
    fn validate_rejects_bad_bonds() {
        let mut m = Molecule::new();
        m.add_atom(Atom::new(6));
        m.add_bond(Bond::single(0, 9));
        assert!(m.validate().is_err());

        let mut m2 = Molecule::new();
        m2.add_atom(Atom::new(6));
        m2.add_bond(Bond::single(0, 0));
        assert!(m2.validate().is_err());
    }

    #[test]
    fn components_and_largest_fragment() {
        // Two disconnected fragments: ethanol + a lone Na atom.
        let mut m = ethanol();
        m.add_atom(Atom::new(11));
        assert_eq!(m.component_count(), 2);
        assert!(!m.is_connected());

        let big = m.largest_fragment();
        assert_eq!(big.atom_count(), 3);
        assert_eq!(big.heavy_atom_count(), 3);
    }

    #[test]
    fn subgraph_renumbers() {
        let m = ethanol();
        let sub = m.subgraph(&[1, 2]);
        assert_eq!(sub.atom_count(), 2);
        assert_eq!(sub.bond_count(), 1);
        assert_eq!(sub.atoms[0].atomic_number, 6);
        assert_eq!(sub.atoms[1].atomic_number, 8);
    }
}
