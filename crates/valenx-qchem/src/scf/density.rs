//! Density-matrix sanity checks — electron count and idempotency.
//!
//! A converged Hartree-Fock density matrix obeys two exact algebraic
//! identities, and checking them is the cheapest way to catch a broken
//! SCF:
//!
//! 1. **Electron count** — `tr(D S)` equals the number of electrons.
//!    For a closed-shell density `D` carries the factor of 2, so
//!    `tr(D S) = N`.
//! 2. **Idempotency** — the spin density is a projector onto the
//!    occupied space: `Dˢ S Dˢ = Dˢ`. For the closed-shell `D = 2 Dˢ`
//!    this reads `½ D S D = D`.
//!
//! [`DensityCheck::closed_shell`] evaluates both and bundles the
//! residuals into a [`DensityCheckReport`].

use nalgebra::DMatrix;

/// Tolerance below which a density residual counts as "satisfied".
pub const DENSITY_CHECK_TOL: f64 = 1.0e-6;

/// The numerical residuals of the density-matrix identities.
#[derive(Copy, Clone, Debug)]
pub struct DensityCheckReport {
    /// `tr(D S)` — should equal the electron count.
    pub electron_count: f64,
    /// The electron count the density *should* integrate to.
    pub expected_electrons: f64,
    /// Largest absolute element of the idempotency residual.
    pub idempotency_residual: f64,
}

impl DensityCheckReport {
    /// The absolute error in the integrated electron count.
    pub fn electron_count_error(&self) -> f64 {
        (self.electron_count - self.expected_electrons).abs()
    }

    /// `true` when both identities hold to [`DENSITY_CHECK_TOL`].
    pub fn is_valid(&self) -> bool {
        self.electron_count_error() < DENSITY_CHECK_TOL
            && self.idempotency_residual < DENSITY_CHECK_TOL
    }
}

/// Density-matrix checks.
pub struct DensityCheck;

impl DensityCheck {
    /// Check a closed-shell (RHF) density `D` against the overlap `S`
    /// for an `n_electrons`-electron molecule.
    ///
    /// The closed-shell density satisfies `tr(D S) = N` and
    /// `½ D S D = D`.
    pub fn closed_shell(
        density: &DMatrix<f64>,
        overlap: &DMatrix<f64>,
        n_electrons: u32,
    ) -> DensityCheckReport {
        let ds = density * overlap;
        let electron_count = ds.trace();
        // Idempotency: ½ D S D − D.
        let dsd = 0.5 * &ds * density;
        let residual = (&dsd - density).abs().max();
        DensityCheckReport {
            electron_count,
            expected_electrons: f64::from(n_electrons),
            idempotency_residual: residual,
        }
    }

    /// Check a single-spin (UHF) density `Dˢ`, which satisfies
    /// `tr(Dˢ S) = nˢ` and `Dˢ S Dˢ = Dˢ`.
    pub fn single_spin(
        density: &DMatrix<f64>,
        overlap: &DMatrix<f64>,
        n_spin: u32,
    ) -> DensityCheckReport {
        let ds = density * overlap;
        let electron_count = ds.trace();
        let dsd = &ds * density;
        let residual = (&dsd - density).abs().max();
        DensityCheckReport {
            electron_count,
            expected_electrons: f64::from(n_spin),
            idempotency_residual: residual,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::BasisSet;
    use crate::geometry::{Atom, MolecularGeometry};
    use crate::integrals::IntegralSet;
    use crate::scf::rhf::{run_rhf_scf, ScfSettings};

    #[test]
    fn converged_rhf_density_passes_checks() {
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7414]).unwrap(),
        ]);
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 2, ScfSettings::default()).unwrap();
        let report = DensityCheck::closed_shell(&res.density, &ints.overlap, 2);
        assert!(report.is_valid(), "report = {report:?}");
        assert!(report.electron_count_error() < 1.0e-8);
    }

    #[test]
    fn water_density_integrates_to_ten_electrons() {
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("O", [0.0, 0.0, 0.1173]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.7572, -0.4692]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, -0.7572, -0.4692]).unwrap(),
        ]);
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 10, ScfSettings::default()).unwrap();
        let report = DensityCheck::closed_shell(&res.density, &ints.overlap, 10);
        assert!((report.electron_count - 10.0).abs() < 1.0e-7);
    }

    #[test]
    fn broken_density_fails_idempotency() {
        // An arbitrary non-projector matrix should fail the check.
        let s = DMatrix::<f64>::identity(2, 2);
        let bad = DMatrix::from_row_slice(2, 2, &[1.0, 0.5, 0.5, 1.0]);
        let report = DensityCheck::closed_shell(&bad, &s, 2);
        assert!(!report.is_valid());
    }
}
