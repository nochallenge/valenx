//! The molecular dipole moment.
//!
//! The dipole moment has an electronic part and a nuclear part:
//!
//! ```text
//! μ = μ_nuclear + μ_electronic
//! μ_nuclear   = Σ_A Z_A R_A
//! μ_electronic = − Σ_{μν} D_{μν} ⟨μ| r |ν⟩
//! ```
//!
//! the electronic part being negative because electrons carry negative
//! charge. The dipole-integral matrices `⟨μ|x|ν⟩` etc. come from the
//! [`IntegralSet`]; both parts are taken about the same origin (the
//! Cartesian origin), so for a neutral molecule the total is
//! origin-independent.
//!
//! The result is reported in atomic units (`e·a₀`) and in debye —
//! [`AU_TO_DEBYE`] is the conversion (1 a.u. = 2.541746 D).

use crate::geometry::MolecularGeometry;
use crate::integrals::IntegralSet;
use nalgebra::DMatrix;

/// Atomic-unit dipole (`e·a₀`) per debye. Multiply a.u. by this to get
/// debye.
pub const AU_TO_DEBYE: f64 = 2.541_746_473;

/// A molecular dipole moment.
#[derive(Copy, Clone, Debug)]
pub struct DipoleMoment {
    /// Cartesian components `[x, y, z]` in atomic units (`e·a₀`).
    pub vector_au: [f64; 3],
}

impl DipoleMoment {
    /// The dipole magnitude in atomic units.
    pub fn magnitude_au(&self) -> f64 {
        (self.vector_au[0].powi(2)
            + self.vector_au[1].powi(2)
            + self.vector_au[2].powi(2))
        .sqrt()
    }

    /// The dipole magnitude in debye.
    pub fn magnitude_debye(&self) -> f64 {
        self.magnitude_au() * AU_TO_DEBYE
    }

    /// The Cartesian components in debye.
    pub fn vector_debye(&self) -> [f64; 3] {
        [
            self.vector_au[0] * AU_TO_DEBYE,
            self.vector_au[1] * AU_TO_DEBYE,
            self.vector_au[2] * AU_TO_DEBYE,
        ]
    }
}

/// Compute the dipole moment from a (total) density matrix.
///
/// `density` is the full electron density — for an RHF calculation that
/// is the closed-shell `D`; for UHF it is `Dᵅ + Dᵝ`.
pub fn dipole_moment(
    geometry: &MolecularGeometry,
    integrals: &IntegralSet,
    density: &DMatrix<f64>,
) -> DipoleMoment {
    // Nuclear part: Σ_A Z_A R_A.
    let mut mu = [0.0f64; 3];
    for atom in &geometry.atoms {
        let z = atom.element.nuclear_charge();
        for (k, m) in mu.iter_mut().enumerate() {
            *m += z * atom.position[k];
        }
    }
    // Electronic part: − Σ_{μν} D_{μν} ⟨μ|r_k|ν⟩.
    for (k, m) in mu.iter_mut().enumerate() {
        let dip_k = &integrals.dipole[k];
        let mut elec = 0.0;
        for i in 0..density.nrows() {
            for j in 0..density.ncols() {
                elec += density[(i, j)] * dip_k[(i, j)];
            }
        }
        *m -= elec;
    }
    DipoleMoment { vector_au: mu }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::BasisSet;
    use crate::geometry::Atom;
    use crate::scf::rhf::{run_rhf_scf, ScfSettings};

    #[test]
    fn homonuclear_h2_has_no_dipole() {
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7414]).unwrap(),
        ]);
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 2, ScfSettings::default()).unwrap();
        let mu = dipole_moment(&geom, &ints, &res.density);
        assert!(mu.magnitude_au() < 1.0e-7, "|μ| = {}", mu.magnitude_au());
    }

    #[test]
    fn water_has_a_dipole() {
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("O", [0.0, 0.0, 0.1173]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.7572, -0.4692]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, -0.7572, -0.4692]).unwrap(),
        ]);
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 10, ScfSettings::default()).unwrap();
        let mu = dipole_moment(&geom, &ints, &res.density);
        // STO-3G water dipole is roughly 1.7 D — at least clearly
        // nonzero and pointing along z by C2v symmetry.
        assert!(mu.magnitude_debye() > 1.0, "|μ| = {} D", mu.magnitude_debye());
        assert!(mu.vector_au[0].abs() < 1.0e-6);
        assert!(mu.vector_au[1].abs() < 1.0e-6);
    }

    #[test]
    fn debye_conversion_is_consistent() {
        let d = DipoleMoment {
            vector_au: [1.0, 0.0, 0.0],
        };
        assert!((d.magnitude_debye() - AU_TO_DEBYE).abs() < 1.0e-12);
        assert!((d.vector_debye()[0] - AU_TO_DEBYE).abs() < 1.0e-12);
    }
}
