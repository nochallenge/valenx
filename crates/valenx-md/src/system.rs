//! The molecular-system data model.
//!
//! **Roadmap feature 1.** Three layered types:
//!
//! - [`Atom`] — one particle: an element / type label, mass, partial
//!   charge, residue / chain bookkeeping. Static properties, no
//!   coordinates.
//! - [`Topology`] — the bonded graph: which atom indices are joined by
//!   [`Bond`]s, [`Angle`]s, proper [`Dihedral`]s and [`Improper`]s,
//!   plus the per-atom [`Atom`] list.
//! - [`System`] — a `Topology` made dynamical: per-atom position and
//!   velocity vectors plus a [`SimBox`].
//!
//! Coordinates live in nm and velocities in nm/ps (see [`crate::units`]).
//! The topology stores bonded *connectivity* only; the numeric force-
//! field parameters that decorate it live in [`crate::forcefield`].

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::{MdError, Result};
use crate::pbc::SimBox;

/// One particle in the system. Carries static properties only —
/// position and velocity live in the owning [`System`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Atom {
    /// Force-field atom-type name (e.g. `"CT"`, `"OW"`, `"HW"`).
    /// Used to look up Lennard-Jones parameters in a force field.
    pub type_name: String,
    /// Element symbol (e.g. `"C"`, `"O"`, `"H"`). Best-effort; may be
    /// empty if a reader could not infer it.
    pub element: String,
    /// Mass in atomic mass units (u = g/mol). Must be positive.
    pub mass: f64,
    /// Partial charge in elementary-charge units (e).
    pub charge: f64,
    /// Residue name this atom belongs to (e.g. `"ALA"`, `"SOL"`).
    pub residue: String,
    /// Residue sequence number.
    pub residue_id: i32,
    /// Chain / segment identifier.
    pub chain: String,
    /// Atom name within its residue (e.g. `"CA"`, `"OW"`).
    pub name: String,
}

impl Atom {
    /// Builds an atom with the given type, mass and charge; residue /
    /// chain fields default to empty / zero and the element is left
    /// blank.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if the mass is not finite and positive.
    pub fn new(type_name: impl Into<String>, mass: f64, charge: f64) -> Result<Self> {
        if !(mass.is_finite() && mass > 0.0) {
            return Err(MdError::invalid(
                "mass",
                format!("must be finite and positive, got {mass}"),
            ));
        }
        if !charge.is_finite() {
            return Err(MdError::invalid("charge", "must be finite"));
        }
        Ok(Atom {
            type_name: type_name.into(),
            element: String::new(),
            mass,
            charge,
            residue: String::new(),
            residue_id: 0,
            chain: String::new(),
            name: String::new(),
        })
    }

    /// Builder-style setter for the element symbol.
    pub fn with_element(mut self, element: impl Into<String>) -> Self {
        self.element = element.into();
        self
    }

    /// Builder-style setter for the atom name within its residue.
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = name.into();
        self
    }

    /// Builder-style setter for the residue name + id.
    pub fn with_residue(mut self, residue: impl Into<String>, residue_id: i32) -> Self {
        self.residue = residue.into();
        self.residue_id = residue_id;
        self
    }

    /// Builder-style setter for the chain identifier.
    pub fn with_chain(mut self, chain: impl Into<String>) -> Self {
        self.chain = chain.into();
        self
    }
}

/// A harmonic bond between two atoms (connectivity only — the spring
/// constant lives in the force field).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Bond {
    /// First atom index.
    pub i: usize,
    /// Second atom index.
    pub j: usize,
}

impl Bond {
    /// A bond between atoms `i` and `j`.
    pub fn new(i: usize, j: usize) -> Self {
        Bond { i, j }
    }
}

/// An angle defined by three atoms `i`-`j`-`k`, vertex at `j`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Angle {
    /// First leg atom.
    pub i: usize,
    /// Vertex atom.
    pub j: usize,
    /// Second leg atom.
    pub k: usize,
}

impl Angle {
    /// An angle `i`-`j`-`k` with the vertex at `j`.
    pub fn new(i: usize, j: usize, k: usize) -> Self {
        Angle { i, j, k }
    }
}

/// A proper dihedral / torsion defined by four atoms `i`-`j`-`k`-`l`,
/// the angle measured about the central `j`-`k` bond.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Dihedral {
    /// First atom.
    pub i: usize,
    /// Second (central-bond) atom.
    pub j: usize,
    /// Third (central-bond) atom.
    pub k: usize,
    /// Fourth atom.
    pub l: usize,
}

impl Dihedral {
    /// A dihedral `i`-`j`-`k`-`l`.
    pub fn new(i: usize, j: usize, k: usize, l: usize) -> Self {
        Dihedral { i, j, k, l }
    }
}

/// An improper dihedral — four atoms holding a planarity / chirality
/// constraint. The convention here: the angle is measured between the
/// `i`-`j`-`k` and `j`-`k`-`l` planes, identical geometry to a proper
/// dihedral but typically with a harmonic (not periodic) potential.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Improper {
    /// Central atom (usually the one whose planarity is constrained).
    pub i: usize,
    /// Neighbour atom.
    pub j: usize,
    /// Neighbour atom.
    pub k: usize,
    /// Neighbour atom.
    pub l: usize,
}

impl Improper {
    /// An improper dihedral over atoms `i`-`j`-`k`-`l`.
    pub fn new(i: usize, j: usize, k: usize, l: usize) -> Self {
        Improper { i, j, k, l }
    }
}

/// The bonded graph: an atom list plus the bond / angle / dihedral /
/// improper index lists.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Topology {
    /// Per-atom static properties.
    pub atoms: Vec<Atom>,
    /// Two-body bonds.
    pub bonds: Vec<Bond>,
    /// Three-body angles.
    pub angles: Vec<Angle>,
    /// Four-body proper dihedrals.
    pub dihedrals: Vec<Dihedral>,
    /// Four-body improper dihedrals.
    pub impropers: Vec<Improper>,
}

impl Topology {
    /// An empty topology.
    pub fn new() -> Self {
        Topology::default()
    }

    /// Number of atoms.
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// Whether the topology has no atoms.
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Appends an atom, returning its index.
    pub fn push_atom(&mut self, atom: Atom) -> usize {
        self.atoms.push(atom);
        self.atoms.len() - 1
    }

    /// Adds a bond after bounds-checking both indices.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if either index is out of range or the
    /// two indices are equal.
    pub fn add_bond(&mut self, i: usize, j: usize) -> Result<()> {
        self.check_index("bond", i)?;
        self.check_index("bond", j)?;
        if i == j {
            return Err(MdError::invalid("bond", "an atom cannot bond to itself"));
        }
        self.bonds.push(Bond::new(i, j));
        Ok(())
    }

    /// Adds an angle after bounds-checking all three indices.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if any index is out of range or two of
    /// them coincide.
    pub fn add_angle(&mut self, i: usize, j: usize, k: usize) -> Result<()> {
        for idx in [i, j, k] {
            self.check_index("angle", idx)?;
        }
        if i == j || j == k || i == k {
            return Err(MdError::invalid("angle", "indices must be distinct"));
        }
        self.angles.push(Angle::new(i, j, k));
        Ok(())
    }

    /// Adds a proper dihedral after bounds-checking all four indices.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if any index is out of range.
    pub fn add_dihedral(&mut self, i: usize, j: usize, k: usize, l: usize) -> Result<()> {
        for idx in [i, j, k, l] {
            self.check_index("dihedral", idx)?;
        }
        self.dihedrals.push(Dihedral::new(i, j, k, l));
        Ok(())
    }

    /// Adds an improper dihedral after bounds-checking all four
    /// indices.
    ///
    /// # Errors
    /// [`MdError::Invalid`] if any index is out of range.
    pub fn add_improper(&mut self, i: usize, j: usize, k: usize, l: usize) -> Result<()> {
        for idx in [i, j, k, l] {
            self.check_index("improper", idx)?;
        }
        self.impropers.push(Improper::new(i, j, k, l));
        Ok(())
    }

    /// Total mass of the system (u).
    pub fn total_mass(&self) -> f64 {
        self.atoms.iter().map(|a| a.mass).sum()
    }

    /// Net charge of the system (e).
    pub fn net_charge(&self) -> f64 {
        self.atoms.iter().map(|a| a.charge).sum()
    }

    /// The set of atom-index pairs that are 1-2 (bonded) or 1-3
    /// (angle-related) neighbours — the pairs conventionally
    /// **excluded** from the nonbonded sum.
    ///
    /// Returns each pair as `(min, max)` so the result can seed a
    /// lookup set directly.
    pub fn nonbonded_exclusions(&self) -> Vec<(usize, usize)> {
        let mut set = std::collections::BTreeSet::new();
        let mut add = |a: usize, b: usize| {
            if a != b {
                set.insert((a.min(b), a.max(b)));
            }
        };
        for b in &self.bonds {
            add(b.i, b.j);
        }
        for a in &self.angles {
            add(a.i, a.j);
            add(a.j, a.k);
            add(a.i, a.k);
        }
        set.into_iter().collect()
    }

    /// The **1-4 pairs** — the two end atoms `(i, l)` of every proper
    /// dihedral that are not already 1-2 or 1-3 neighbours. In AMBER/OPLS
    /// force fields these take a *scaled* nonbonded interaction (see
    /// [`crate::nonbonded::ScaledPairs14`]), not the full one. Returned as
    /// sorted, deduplicated `(min, max)` pairs.
    pub fn one_four_pairs(&self) -> Vec<(usize, usize)> {
        let excluded: std::collections::BTreeSet<(usize, usize)> =
            self.nonbonded_exclusions().into_iter().collect();
        let mut set = std::collections::BTreeSet::new();
        for d in &self.dihedrals {
            let (a, b) = (d.i.min(d.l), d.i.max(d.l));
            if a != b && !excluded.contains(&(a, b)) {
                set.insert((a, b));
            }
        }
        set.into_iter().collect()
    }

    fn check_index(&self, what: &'static str, idx: usize) -> Result<()> {
        if idx >= self.atoms.len() {
            Err(MdError::invalid(
                what,
                format!(
                    "atom index {idx} out of range (system has {} atoms)",
                    self.atoms.len()
                ),
            ))
        } else {
            Ok(())
        }
    }

    /// Validates internal consistency: every bonded-term index is in
    /// range and every mass is positive.
    ///
    /// # Errors
    /// The first inconsistency found.
    pub fn validate(&self) -> Result<()> {
        for (n, a) in self.atoms.iter().enumerate() {
            if !(a.mass.is_finite() && a.mass > 0.0) {
                return Err(MdError::invalid(
                    "mass",
                    format!("atom {n} has non-positive mass {}", a.mass),
                ));
            }
        }
        let n = self.atoms.len();
        let check = |what: &'static str, idx: usize| -> Result<()> {
            if idx >= n {
                Err(MdError::invalid(
                    what,
                    format!("index {idx} out of range ({n} atoms)"),
                ))
            } else {
                Ok(())
            }
        };
        for b in &self.bonds {
            check("bond", b.i)?;
            check("bond", b.j)?;
        }
        for a in &self.angles {
            check("angle", a.i)?;
            check("angle", a.j)?;
            check("angle", a.k)?;
        }
        for d in &self.dihedrals {
            for idx in [d.i, d.j, d.k, d.l] {
                check("dihedral", idx)?;
            }
        }
        for im in &self.impropers {
            for idx in [im.i, im.j, im.k, im.l] {
                check("improper", idx)?;
            }
        }
        Ok(())
    }
}

/// A dynamical system: a [`Topology`] plus per-atom positions and
/// velocities and a periodic [`SimBox`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct System {
    /// The bonded graph and atom properties.
    pub topology: Topology,
    /// Per-atom positions (nm). Length equals `topology.len()`.
    pub positions: Vec<Vector3<f64>>,
    /// Per-atom velocities (nm/ps). Length equals `topology.len()`.
    pub velocities: Vec<Vector3<f64>>,
    /// The periodic simulation cell.
    pub cell: SimBox,
}

impl System {
    /// Builds a system from a topology and a position list, with zero
    /// initial velocities and a non-periodic box.
    ///
    /// # Errors
    /// [`MdError::DimensionMismatch`] if the position count differs
    /// from the atom count; propagates [`Topology::validate`] errors.
    pub fn new(topology: Topology, positions: Vec<Vector3<f64>>) -> Result<Self> {
        topology.validate()?;
        if positions.len() != topology.len() {
            return Err(MdError::dimension(format!(
                "{} positions for {} atoms",
                positions.len(),
                topology.len()
            )));
        }
        let n = topology.len();
        Ok(System {
            topology,
            positions,
            velocities: vec![Vector3::zeros(); n],
            cell: SimBox::none(),
        })
    }

    /// Number of atoms.
    pub fn len(&self) -> usize {
        self.topology.len()
    }

    /// Whether the system has no atoms.
    pub fn is_empty(&self) -> bool {
        self.topology.is_empty()
    }

    /// Builder-style setter for the periodic box.
    pub fn with_cell(mut self, cell: SimBox) -> Self {
        self.cell = cell;
        self
    }

    /// Sets the velocity array.
    ///
    /// # Errors
    /// [`MdError::DimensionMismatch`] if the length is wrong.
    pub fn set_velocities(&mut self, velocities: Vec<Vector3<f64>>) -> Result<()> {
        if velocities.len() != self.len() {
            return Err(MdError::dimension(format!(
                "{} velocities for {} atoms",
                velocities.len(),
                self.len()
            )));
        }
        self.velocities = velocities;
        Ok(())
    }

    /// Number of mechanical degrees of freedom — `3N` minus the 3
    /// removed when the centre-of-mass motion is constrained.
    ///
    /// `constraints` is the number of additional holonomic
    /// constraints (e.g. SHAKE-fixed bonds) the caller is applying;
    /// pass 0 if none.
    pub fn degrees_of_freedom(&self, constraints: usize) -> usize {
        let raw = 3 * self.len();
        raw.saturating_sub(3).saturating_sub(constraints)
    }

    /// Total kinetic energy `Σ ½ m v²` (kJ/mol).
    pub fn kinetic_energy(&self) -> f64 {
        self.topology
            .atoms
            .iter()
            .zip(&self.velocities)
            .map(|(a, v)| 0.5 * a.mass * v.norm_squared())
            .sum()
    }

    /// Instantaneous temperature (K) from the equipartition theorem,
    /// `T = 2·KE / (dof·k_B)`.
    ///
    /// `constraints` is forwarded to [`degrees_of_freedom`](Self::degrees_of_freedom).
    pub fn temperature(&self, constraints: usize) -> f64 {
        let dof = self.degrees_of_freedom(constraints);
        if dof == 0 {
            return 0.0;
        }
        2.0 * self.kinetic_energy() / (dof as f64 * crate::units::BOLTZMANN)
    }

    /// Centre of mass of the system (nm).
    pub fn center_of_mass(&self) -> Vector3<f64> {
        let total = self.topology.total_mass();
        if total <= 0.0 {
            return Vector3::zeros();
        }
        let mut com = Vector3::zeros();
        for (a, p) in self.topology.atoms.iter().zip(&self.positions) {
            com += a.mass * p;
        }
        com / total
    }

    /// Total linear momentum `Σ m v` (u·nm/ps).
    pub fn linear_momentum(&self) -> Vector3<f64> {
        let mut p = Vector3::zeros();
        for (a, v) in self.topology.atoms.iter().zip(&self.velocities) {
            p += a.mass * v;
        }
        p
    }

    /// Removes net centre-of-mass translation by subtracting the
    /// mass-weighted mean velocity from every atom. Called before a
    /// run so the conserved total momentum is zero.
    pub fn remove_com_motion(&mut self) {
        let total = self.topology.total_mass();
        if total <= 0.0 {
            return;
        }
        let v_com = self.linear_momentum() / total;
        for v in &mut self.velocities {
            *v -= v_com;
        }
    }

    /// Wraps every atom position into the primary periodic cell.
    pub fn wrap_into_cell(&mut self) {
        if !self.cell.is_periodic() {
            return;
        }
        for p in &mut self.positions {
            *p = self.cell.wrap(*p);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn diatomic() -> System {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 12.0, -0.2).unwrap());
        top.push_atom(Atom::new("B", 1.0, 0.2).unwrap());
        top.add_bond(0, 1).unwrap();
        let pos = vec![Vector3::zeros(), Vector3::new(0.1, 0.0, 0.0)];
        System::new(top, pos).unwrap()
    }

    #[test]
    fn atom_rejects_bad_mass_and_charge() {
        assert!(Atom::new("X", 0.0, 0.0).is_err());
        assert!(Atom::new("X", -1.0, 0.0).is_err());
        assert!(Atom::new("X", 1.0, f64::NAN).is_err());
        assert!(Atom::new("X", 1.0, 0.5).is_ok());
    }

    #[test]
    fn topology_bounds_checks_bonded_terms() {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 1.0, 0.0).unwrap());
        assert!(top.add_bond(0, 5).is_err());
        assert!(top.add_bond(0, 0).is_err());
        top.push_atom(Atom::new("B", 1.0, 0.0).unwrap());
        assert!(top.add_bond(0, 1).is_ok());
        assert!(top.add_angle(0, 1, 0).is_err());
    }

    #[test]
    fn system_rejects_size_mismatch() {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 1.0, 0.0).unwrap());
        let err = System::new(top, vec![]).unwrap_err();
        assert_eq!(err.category(), "input");
    }

    #[test]
    fn net_charge_and_total_mass() {
        let s = diatomic();
        assert!((s.topology.net_charge()).abs() < 1e-12);
        assert!((s.topology.total_mass() - 13.0).abs() < 1e-12);
    }

    #[test]
    fn kinetic_energy_and_temperature() {
        let mut s = diatomic();
        s.set_velocities(vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 2.0, 0.0),
        ])
        .unwrap();
        // KE = 0.5*12*1 + 0.5*1*4 = 6 + 2 = 8 kJ/mol.
        assert!((s.kinetic_energy() - 8.0).abs() < 1e-9);
        // dof = 3*2 - 3 = 3; T = 2*8/(3*kB) > 0.
        assert!(s.temperature(0) > 0.0);
    }

    #[test]
    fn com_motion_removal_zeros_momentum() {
        let mut s = diatomic();
        s.set_velocities(vec![
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(0.0, 5.0, 0.0),
        ])
        .unwrap();
        s.remove_com_motion();
        assert!(s.linear_momentum().norm() < 1e-9);
    }

    #[test]
    fn exclusions_cover_1_2_and_1_3() {
        let mut top = Topology::new();
        for _ in 0..3 {
            top.push_atom(Atom::new("A", 1.0, 0.0).unwrap());
        }
        top.add_bond(0, 1).unwrap();
        top.add_bond(1, 2).unwrap();
        top.add_angle(0, 1, 2).unwrap();
        let ex = top.nonbonded_exclusions();
        // (0,1) and (1,2) bonded, (0,2) is the 1-3 pair.
        assert!(ex.contains(&(0, 1)));
        assert!(ex.contains(&(1, 2)));
        assert!(ex.contains(&(0, 2)));
    }

    #[test]
    fn degrees_of_freedom_accounts_for_constraints() {
        let s = diatomic();
        assert_eq!(s.degrees_of_freedom(0), 3);
        assert_eq!(s.degrees_of_freedom(1), 2);
    }
}
