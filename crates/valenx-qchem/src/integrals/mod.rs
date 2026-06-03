//! Molecular integrals over Gaussian basis functions.
//!
//! This module is the numerical engine beneath the SCF: it turns a
//! [`MolecularGeometry`] and a [`BasisSet`] into the matrices the
//! Hartree-Fock equations need.
//!
//! - [`boys`] — the Boys function `F_n(x)`.
//! - [`mcmurchie`] — the McMurchie-Davidson Hermite recursions shared
//!   by every integral.
//! - [`one_electron`] — overlap, kinetic-energy, nuclear-attraction and
//!   multipole integrals.
//! - [`two_electron`] — the four-index electron-repulsion tensor.
//!
//! The top-level [`IntegralSet`] bundles the assembled matrices for one
//! molecule + basis: the overlap `S`, kinetic `T`, nuclear-attraction
//! `V` and the dipole-moment integral matrices, plus the
//! [`EriTensor`]. [`nuclear_repulsion`] supplies the classical
//! nucleus-nucleus energy.

pub mod boys;
pub mod mcmurchie;
pub mod one_electron;
pub mod two_electron;

use crate::basis::BasisSet;
use crate::geometry::MolecularGeometry;
use nalgebra::DMatrix;
use two_electron::EriTensor;

/// The classical nuclear-repulsion energy
/// `E_nuc = Σ_{A<B} Z_A Z_B / R_AB` in Hartree.
pub fn nuclear_repulsion(geometry: &MolecularGeometry) -> f64 {
    let mut e = 0.0;
    let atoms = &geometry.atoms;
    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            let zi = atoms[i].element.nuclear_charge();
            let zj = atoms[j].element.nuclear_charge();
            let d2: f64 = (0..3)
                .map(|k| {
                    let d = atoms[i].position[k] - atoms[j].position[k];
                    d * d
                })
                .sum();
            e += zi * zj / d2.sqrt();
        }
    }
    e
}

/// The one-electron potential matrix from external point charges
/// `(q, position_bohr)` — the same operator as the nuclear attraction,
/// with the external charges as the centres. Add it to a core
/// Hamiltonian for **electrostatic QM/MM embedding**: the MM charges then
/// enter the SCF and polarize the electron density. (The classical
/// nuclei–charge term is the caller's, as for nuclear repulsion.)
pub fn external_charge_potential(basis: &BasisSet, charges: &[(f64, [f64; 3])]) -> DMatrix<f64> {
    let n = basis.n_functions();
    let f = &basis.functions;
    let mut v = DMatrix::<f64>::zeros(n, n);
    for mu in 0..n {
        for nu in 0..=mu {
            let val = one_electron::nuclear_attraction(&f[mu], &f[nu], charges);
            v[(mu, nu)] = val;
            v[(nu, mu)] = val;
        }
    }
    v
}

/// The assembled one- and two-electron integrals for a molecule + basis.
#[derive(Clone, Debug)]
pub struct IntegralSet {
    /// Overlap matrix `S` (`n × n`).
    pub overlap: DMatrix<f64>,
    /// Kinetic-energy matrix `T` (`n × n`).
    pub kinetic: DMatrix<f64>,
    /// Nuclear-attraction matrix `V` (`n × n`).
    pub nuclear: DMatrix<f64>,
    /// Dipole-moment integral matrices `[μx, μy, μz]` about the origin.
    pub dipole: [DMatrix<f64>; 3],
    /// The four-index electron-repulsion tensor.
    pub eri: EriTensor,
    /// Nuclear-repulsion energy (Hartree).
    pub e_nuclear: f64,
}

impl IntegralSet {
    /// The basis dimension `n`.
    #[inline]
    pub fn n(&self) -> usize {
        self.overlap.nrows()
    }

    /// The core Hamiltonian `H = T + V`.
    pub fn core_hamiltonian(&self) -> DMatrix<f64> {
        &self.kinetic + &self.nuclear
    }

    /// Compute every integral for `geometry` in `basis`.
    ///
    /// The nuclear-attraction integrals use every atom's nuclear charge
    /// as a point charge; the dipole integrals are taken about the
    /// Cartesian origin.
    pub fn compute(geometry: &MolecularGeometry, basis: &BasisSet) -> IntegralSet {
        let n = basis.n_functions();
        let f = &basis.functions;

        // Point charges for the nuclear-attraction operator.
        let charges: Vec<(f64, [f64; 3])> = geometry
            .atoms
            .iter()
            .map(|a| (a.element.nuclear_charge(), a.position))
            .collect();

        let mut s = DMatrix::<f64>::zeros(n, n);
        let mut t = DMatrix::<f64>::zeros(n, n);
        let mut v = DMatrix::<f64>::zeros(n, n);
        let mut dx = DMatrix::<f64>::zeros(n, n);
        let mut dy = DMatrix::<f64>::zeros(n, n);
        let mut dz = DMatrix::<f64>::zeros(n, n);

        for mu in 0..n {
            for nu in 0..=mu {
                let s_v = one_electron::overlap(&f[mu], &f[nu]);
                let t_v = one_electron::kinetic(&f[mu], &f[nu]);
                let v_v = one_electron::nuclear_attraction(&f[mu], &f[nu], &charges);
                let dx_v =
                    one_electron::multipole(&f[mu], &f[nu], [0.0; 3], (1, 0, 0));
                let dy_v =
                    one_electron::multipole(&f[mu], &f[nu], [0.0; 3], (0, 1, 0));
                let dz_v =
                    one_electron::multipole(&f[mu], &f[nu], [0.0; 3], (0, 0, 1));
                for &(i, j) in &[(mu, nu), (nu, mu)] {
                    s[(i, j)] = s_v;
                    t[(i, j)] = t_v;
                    v[(i, j)] = v_v;
                    dx[(i, j)] = dx_v;
                    dy[(i, j)] = dy_v;
                    dz[(i, j)] = dz_v;
                }
            }
        }

        IntegralSet {
            overlap: s,
            kinetic: t,
            nuclear: v,
            dipole: [dx, dy, dz],
            eri: EriTensor::build(basis),
            e_nuclear: nuclear_repulsion(geometry),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Atom;

    fn h2_geometry() -> MolecularGeometry {
        MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7414]).unwrap(),
        ])
    }

    #[test]
    fn nuclear_repulsion_of_h2() {
        // Z=1 each; E_nuc = 1/R with R the bond length in bohr.
        let geom = h2_geometry();
        let e = nuclear_repulsion(&geom);
        let r_bohr = 0.7414 * crate::geometry::BOHR_PER_ANGSTROM;
        assert!((e - 1.0 / r_bohr).abs() < 1.0e-10, "E_nuc = {e}");
    }

    #[test]
    fn nuclear_repulsion_single_atom_is_zero() {
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("He", [0.0, 0.0, 0.0]).unwrap(),
        ]);
        assert_eq!(nuclear_repulsion(&geom), 0.0);
    }

    #[test]
    fn integral_set_dimensions_and_symmetry() {
        let geom = h2_geometry();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        assert_eq!(ints.n(), 2);
        // S, T, V are symmetric.
        for m in [&ints.overlap, &ints.kinetic, &ints.nuclear] {
            assert!((m[(0, 1)] - m[(1, 0)]).abs() < 1.0e-13);
        }
        // S has unit diagonal (normalised basis).
        assert!((ints.overlap[(0, 0)] - 1.0).abs() < 1.0e-10);
        assert!((ints.overlap[(1, 1)] - 1.0).abs() < 1.0e-10);
    }

    #[test]
    fn core_hamiltonian_is_t_plus_v() {
        let geom = h2_geometry();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let h = ints.core_hamiltonian();
        assert!((h[(0, 0)] - (ints.kinetic[(0, 0)] + ints.nuclear[(0, 0)])).abs() < 1.0e-13);
    }

    #[test]
    fn h2_overlap_off_diagonal_is_physical() {
        // STO-3G H2 overlap S_12 is a known ~0.659 at this geometry.
        let geom = h2_geometry();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let s12 = ints.overlap[(0, 1)];
        assert!(s12 > 0.6 && s12 < 0.7, "S_12 = {s12}");
    }
}
