//! Mayer / Wiberg bond orders.
//!
//! A *bond order* is a quantum-mechanical count of the electron pairs
//! shared between two atoms — roughly 1 for a single bond, 2 for a
//! double bond. The Mayer bond order generalises the Wiberg index to a
//! non-orthogonal basis:
//!
//! ```text
//! B_AB = Σ_{μ∈A} Σ_{ν∈B} (D S)_{μν} (D S)_{νμ}
//! ```
//!
//! where `D` is the total density and `S` the overlap. (In an
//! orthonormal basis `S = 1` and this is exactly the Wiberg index.)
//! The *valence* of an atom is the sum of its bond orders to all other
//! atoms.
//!
//! For a closed-shell RHF density the spin-resolved Mayer formula
//! collapses to the simple expression above, which is what
//! [`mayer_bond_orders`] evaluates.

use crate::basis::BasisSet;
use crate::geometry::MolecularGeometry;
use nalgebra::DMatrix;

/// A table of Mayer bond orders for a molecule.
#[derive(Clone, Debug)]
pub struct BondOrderMatrix {
    /// Number of atoms.
    pub n_atoms: usize,
    /// Symmetric `n_atoms × n_atoms` bond-order matrix, row-major.
    orders: Vec<f64>,
}

impl BondOrderMatrix {
    /// The bond order between atoms `a` and `b`.
    pub fn order(&self, a: usize, b: usize) -> f64 {
        self.orders[a * self.n_atoms + b]
    }

    /// The Mayer valence of atom `a` — the sum of its bond orders to
    /// every other atom.
    pub fn valence(&self, a: usize) -> f64 {
        (0..self.n_atoms)
            .filter(|&b| b != a)
            .map(|b| self.order(a, b))
            .sum()
    }

    /// Every atom pair `(a, b, order)` with `a < b` whose bond order
    /// exceeds `threshold` — a quick "what is bonded to what" list.
    pub fn significant_bonds(&self, threshold: f64) -> Vec<(usize, usize, f64)> {
        let mut out = Vec::new();
        for a in 0..self.n_atoms {
            for b in (a + 1)..self.n_atoms {
                let o = self.order(a, b);
                if o > threshold {
                    out.push((a, b, o));
                }
            }
        }
        out
    }
}

/// Compute the Mayer bond-order matrix from a (total) density matrix.
///
/// `density` is the full electron density (RHF `D`, or UHF `Dᵅ + Dᵝ`).
pub fn mayer_bond_orders(
    geometry: &MolecularGeometry,
    basis: &BasisSet,
    density: &DMatrix<f64>,
    overlap: &DMatrix<f64>,
) -> BondOrderMatrix {
    let ds = density * overlap;
    let owners: Vec<usize> = basis.functions.iter().map(|f| f.atom_index).collect();
    let n_atoms = geometry.n_atoms();
    let mut orders = vec![0.0; n_atoms * n_atoms];

    let nbf = ds.nrows();
    for mu in 0..nbf {
        let a = owners[mu];
        for nu in 0..nbf {
            let b = owners[nu];
            // B_AB accumulates (DS)_{μν} (DS)_{νμ}.
            orders[a * n_atoms + b] += ds[(mu, nu)] * ds[(nu, mu)];
        }
    }
    // The diagonal A==A is not a bond order — zero it.
    for a in 0..n_atoms {
        orders[a * n_atoms + a] = 0.0;
    }
    BondOrderMatrix { n_atoms, orders }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Atom;
    use crate::integrals::IntegralSet;
    use crate::scf::rhf::{run_rhf_scf, ScfSettings};

    #[test]
    fn h2_has_a_single_bond() {
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7414]).unwrap(),
        ]);
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 2, ScfSettings::default()).unwrap();
        let bo = mayer_bond_orders(&geom, &basis, &res.density, &ints.overlap);
        // H2: bond order should be very close to 1.
        assert!(
            (bo.order(0, 1) - 1.0).abs() < 0.05,
            "H-H bond order = {}",
            bo.order(0, 1)
        );
        // Each hydrogen has valence ≈ 1.
        assert!((bo.valence(0) - 1.0).abs() < 0.05);
    }

    #[test]
    fn water_oxygen_has_valence_near_two() {
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("O", [0.0, 0.0, 0.1173]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.7572, -0.4692]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, -0.7572, -0.4692]).unwrap(),
        ]);
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 10, ScfSettings::default()).unwrap();
        let bo = mayer_bond_orders(&geom, &basis, &res.density, &ints.overlap);
        // Two O-H single bonds → oxygen valence near 2.
        assert!(bo.valence(0) > 1.6 && bo.valence(0) < 2.3, "O valence {}", bo.valence(0));
        // The two O-H bonds are the significant ones.
        let bonds = bo.significant_bonds(0.5);
        assert_eq!(bonds.len(), 2);
        // The H-H "bond" is negligible.
        assert!(bo.order(1, 2) < 0.3);
    }

    #[test]
    fn bond_order_matrix_is_symmetric() {
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("O", [0.0, 0.0, 0.1173]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.7572, -0.4692]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, -0.7572, -0.4692]).unwrap(),
        ]);
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 10, ScfSettings::default()).unwrap();
        let bo = mayer_bond_orders(&geom, &basis, &res.density, &ints.overlap);
        for a in 0..3 {
            for b in 0..3 {
                assert!((bo.order(a, b) - bo.order(b, a)).abs() < 1.0e-10);
            }
        }
    }
}
