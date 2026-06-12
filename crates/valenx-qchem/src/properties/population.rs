//! Population analysis — Mulliken and Löwdin atomic partial charges.
//!
//! Population analysis partitions the molecule's electrons among its
//! atoms, which gives a set of atomic partial charges. Two schemes:
//!
//! - **Mulliken** — assigns the overlap population `D S` to atoms by
//!   the basis function each row belongs to. Cheap and standard, but
//!   notoriously basis-set-sensitive.
//! - **Löwdin** — first symmetrically orthogonalises the basis,
//!   `D' = S^{1/2} D S^{1/2}`, then assigns the diagonal of `D'`. Less
//!   basis-set-sensitive because the orthogonalised functions are
//!   atom-localised.
//!
//! For each scheme the atomic population is summed over the atom's
//! basis functions and the partial charge is `Z_A − population_A`.

use crate::basis::BasisSet;
use crate::geometry::MolecularGeometry;
use crate::scf::linalg::symmetric_eigh;
use nalgebra::DMatrix;

/// Per-atom population-analysis output.
#[derive(Clone, Debug)]
pub struct PopulationAnalysis {
    /// Gross electronic population on each atom (electrons).
    pub atomic_populations: Vec<f64>,
    /// Net partial charge on each atom (`Z − population`), in `e`.
    pub partial_charges: Vec<f64>,
}

impl PopulationAnalysis {
    /// The total charge — the partial charges summed. Equals the
    /// molecular charge to numerical precision.
    pub fn total_charge(&self) -> f64 {
        self.partial_charges.iter().sum()
    }
}

/// Map each basis function to its owning atom.
fn function_atoms(basis: &BasisSet) -> Vec<usize> {
    basis.functions.iter().map(|f| f.atom_index).collect()
}

/// Sum per-basis-function gross populations onto atoms and form the
/// partial charges.
fn assemble(geometry: &MolecularGeometry, basis: &BasisSet, gross: &[f64]) -> PopulationAnalysis {
    let owners = function_atoms(basis);
    let n_atoms = geometry.n_atoms();
    let mut atomic = vec![0.0; n_atoms];
    for (f, &pop) in gross.iter().enumerate() {
        atomic[owners[f]] += pop;
    }
    let charges = geometry
        .atoms
        .iter()
        .zip(&atomic)
        .map(|(a, &pop)| a.element.nuclear_charge() - pop)
        .collect();
    PopulationAnalysis {
        atomic_populations: atomic,
        partial_charges: charges,
    }
}

/// Mulliken population analysis from a (total) density matrix.
///
/// The gross population of basis function `μ` is `(D S)_{μμ}`.
pub fn mulliken(
    geometry: &MolecularGeometry,
    basis: &BasisSet,
    density: &DMatrix<f64>,
    overlap: &DMatrix<f64>,
) -> PopulationAnalysis {
    let ds = density * overlap;
    let gross: Vec<f64> = (0..ds.nrows()).map(|i| ds[(i, i)]).collect();
    assemble(geometry, basis, &gross)
}

/// Löwdin population analysis from a (total) density matrix.
///
/// The gross population of basis function `μ` is the diagonal of the
/// symmetrically-orthogonalised density `S^{1/2} D S^{1/2}`.
pub fn lowdin(
    geometry: &MolecularGeometry,
    basis: &BasisSet,
    density: &DMatrix<f64>,
    overlap: &DMatrix<f64>,
) -> PopulationAnalysis {
    let s_half = sqrt_spd(overlap);
    let d_ortho = &s_half * density * &s_half;
    let gross: Vec<f64> = (0..d_ortho.nrows()).map(|i| d_ortho[(i, i)]).collect();
    assemble(geometry, basis, &gross)
}

/// The symmetric positive-definite square root `S^{1/2}` via the
/// eigendecomposition of `S`.
pub(crate) fn sqrt_spd(s: &DMatrix<f64>) -> DMatrix<f64> {
    let (vals, vecs) = symmetric_eigh(s);
    let n = s.nrows();
    let mut sqrt = DMatrix::<f64>::zeros(n, n);
    for i in 0..n {
        sqrt[(i, i)] = vals[i].max(0.0).sqrt();
    }
    &vecs * sqrt * vecs.transpose()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Atom;
    use crate::integrals::IntegralSet;
    use crate::scf::rhf::{run_rhf_scf, ScfSettings};

    fn h2() -> MolecularGeometry {
        MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7414]).unwrap(),
        ])
    }

    fn water() -> MolecularGeometry {
        MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("O", [0.0, 0.0, 0.1173]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.7572, -0.4692]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, -0.7572, -0.4692]).unwrap(),
        ])
    }

    #[test]
    fn h2_is_symmetric_and_neutral() {
        let geom = h2();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 2, ScfSettings::default()).unwrap();
        let pop = mulliken(&geom, &basis, &res.density, &ints.overlap);
        // Homonuclear: each H carries one electron, zero charge.
        assert!((pop.partial_charges[0]).abs() < 1.0e-8);
        assert!((pop.partial_charges[1]).abs() < 1.0e-8);
        assert!(pop.total_charge().abs() < 1.0e-7);
    }

    #[test]
    fn water_oxygen_is_negative() {
        let geom = water();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 10, ScfSettings::default()).unwrap();
        let pop = mulliken(&geom, &basis, &res.density, &ints.overlap);
        // Oxygen pulls electron density: it should be negative,
        // the hydrogens positive.
        assert!(
            pop.partial_charges[0] < 0.0,
            "O charge {}",
            pop.partial_charges[0]
        );
        assert!(pop.partial_charges[1] > 0.0);
        assert!(pop.total_charge().abs() < 1.0e-6);
    }

    #[test]
    fn lowdin_also_conserves_charge() {
        let geom = water();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 10, ScfSettings::default()).unwrap();
        let pop = lowdin(&geom, &basis, &res.density, &ints.overlap);
        assert!(pop.total_charge().abs() < 1.0e-6);
        assert!(pop.partial_charges[0] < 0.0);
    }

    #[test]
    fn sqrt_spd_squares_back() {
        let s = DMatrix::from_row_slice(2, 2, &[1.0, 0.25, 0.25, 1.0]);
        let root = sqrt_spd(&s);
        let back = &root * &root;
        for i in 0..2 {
            for j in 0..2 {
                assert!((back[(i, j)] - s[(i, j)]).abs() < 1.0e-12);
            }
        }
    }
}
