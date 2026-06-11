//! Nonbonded interactions — the intermolecular force terms.
//!
//! **Roadmap features 7–11.** The pairwise interactions between atoms
//! that are *not* directly bonded:
//!
//! - [`lj`] — Lennard-Jones van der Waals, with a hard cutoff and an
//!   energy shift so the potential is continuous at the cutoff
//!   (feature 7).
//! - [`coulomb`] — electrostatics: a direct Coulomb sum and a
//!   reaction-field correction for the cut-off long-range tail
//!   (feature 8).
//! - [`pme`] — an Ewald / particle-mesh-Ewald v1, splitting the
//!   Coulomb sum into a screened real-space part and a reciprocal-
//!   space structure-factor sum (feature 9).
//! - [`neighbor`] — a cell list and a Verlet neighbour list so the
//!   pair loop is `O(N)` rather than `O(N²)` (feature 10).
//! - The minimum-image convention (feature 11) lives in
//!   [`crate::pbc`] and is used by every term here.
//!
//! Every nonbonded term skips the **excluded pairs** — the 1-2 and
//! 1-3 neighbours from [`crate::system::Topology::nonbonded_exclusions`]
//! — because those interactions are already covered by the bonded
//! terms. The exclusion set is passed in explicitly so the caller
//! controls the policy.

pub mod coulomb;
pub mod lj;
pub mod neighbor;
pub mod pme;
pub mod scaled14;

pub use scaled14::ScaledPairs14;

use std::collections::HashSet;

use crate::system::Topology;

/// A set of atom-index pairs excluded from the nonbonded sum.
///
/// Wraps a `HashSet<(usize, usize)>` keyed by the ordered pair
/// `(min, max)` so membership tests are `O(1)`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ExclusionSet {
    pairs: HashSet<(usize, usize)>,
}

impl ExclusionSet {
    /// An empty exclusion set — every pair interacts.
    pub fn none() -> Self {
        ExclusionSet::default()
    }

    /// The standard 1-2 / 1-3 exclusion set for a topology.
    pub fn from_topology(topology: &Topology) -> Self {
        ExclusionSet {
            pairs: topology.nonbonded_exclusions().into_iter().collect(),
        }
    }

    /// Adds an excluded pair (order-independent).
    pub fn insert(&mut self, a: usize, b: usize) {
        if a != b {
            self.pairs.insert((a.min(b), a.max(b)));
        }
    }

    /// Whether the pair `(a, b)` is excluded.
    pub fn contains(&self, a: usize, b: usize) -> bool {
        self.pairs.contains(&(a.min(b), a.max(b)))
    }

    /// Number of excluded pairs.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    /// Whether the set is empty.
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::{Atom, Topology};

    #[test]
    fn exclusion_set_is_order_independent() {
        let mut ex = ExclusionSet::none();
        ex.insert(3, 1);
        assert!(ex.contains(1, 3));
        assert!(ex.contains(3, 1));
        assert!(!ex.contains(1, 2));
    }

    #[test]
    fn from_topology_picks_up_bonds_and_angles() {
        let mut top = Topology::new();
        for _ in 0..3 {
            top.push_atom(Atom::new("A", 1.0, 0.0).unwrap());
        }
        top.add_bond(0, 1).unwrap();
        top.add_bond(1, 2).unwrap();
        top.add_angle(0, 1, 2).unwrap();
        let ex = ExclusionSet::from_topology(&top);
        assert!(ex.contains(0, 1));
        assert!(ex.contains(0, 2)); // 1-3 pair
    }
}
