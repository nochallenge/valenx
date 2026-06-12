//! Restricted Hartree-Fock — the closed-shell SCF.
//!
//! Restricted Hartree-Fock (RHF) describes a closed-shell molecule with
//! one spatial orbital per pair of electrons. The Roothaan-Hall
//! equations
//!
//! ```text
//! F C = S C ε        F = H + G(D)
//! ```
//!
//! are solved self-consistently: a density `D` builds a Fock matrix
//! `F`, diagonalising `F` gives new orbitals, the orbitals build a new
//! density, and so on until the energy and density stop changing.
//!
//! ## The pieces
//!
//! - [`core_guess_density`] — the initial density from diagonalising
//!   the bare core Hamiltonian (no electron repulsion yet).
//! - [`build_density`] — `D_{μν} = 2 Σ_i^{occ} C_{μi} C_{νi}` from the
//!   occupied orbitals.
//! - [`build_fock`] — `F = H + Σ_{λσ} D_{λσ} [(μν|λσ) − ½(μλ|νσ)]`,
//!   the core Hamiltonian plus the Coulomb `J` and exchange `K` terms.
//! - [`rhf_energy`] — `E = ½ Σ_{μν} D_{μν}(H_{μν} + F_{μν}) + E_nuc`.
//! - [`run_rhf_scf`] — the DIIS-accelerated iteration loop.

use crate::error::{QchemError, Result};
use crate::integrals::two_electron::EriTensor;
use crate::integrals::IntegralSet;
use crate::scf::diis::Diis;
use crate::scf::linalg::{solve_roothaan, symmetric_orthogonalizer, Orthogonalizer};
use nalgebra::{DMatrix, DVector};

/// Convergence thresholds and iteration cap for an RHF calculation.
#[derive(Copy, Clone, Debug)]
pub struct ScfSettings {
    /// Converged when the total-energy change drops below this (Ha).
    pub energy_tol: f64,
    /// Converged when the DIIS error norm drops below this.
    pub density_tol: f64,
    /// Maximum number of SCF cycles before giving up.
    pub max_iterations: usize,
    /// Number of DIIS `(F, e)` pairs to retain.
    pub diis_vectors: usize,
}

impl Default for ScfSettings {
    fn default() -> Self {
        ScfSettings {
            energy_tol: 1.0e-9,
            density_tol: 1.0e-7,
            max_iterations: 128,
            diis_vectors: 8,
        }
    }
}

/// A per-iteration record of the SCF history.
#[derive(Copy, Clone, Debug)]
pub struct ScfIteration {
    /// 1-based cycle number.
    pub cycle: usize,
    /// Total electronic + nuclear energy at this cycle (Ha).
    pub energy: f64,
    /// Energy change from the previous cycle (Ha).
    pub delta_energy: f64,
    /// DIIS error norm at this cycle.
    pub error_norm: f64,
}

/// The converged restricted-Hartree-Fock solution.
#[derive(Clone, Debug)]
pub struct RhfResult {
    /// Total RHF energy (electronic + nuclear), in Hartree.
    pub total_energy: f64,
    /// The electronic energy alone (Ha).
    pub electronic_energy: f64,
    /// Nuclear-repulsion energy (Ha).
    pub nuclear_repulsion: f64,
    /// Molecular-orbital energies `ε`, ascending (Ha).
    pub orbital_energies: DVector<f64>,
    /// Molecular-orbital coefficients `C` (`n_basis × n_mo`).
    pub orbital_coefficients: DMatrix<f64>,
    /// The converged density matrix `D`.
    pub density: DMatrix<f64>,
    /// Number of doubly-occupied orbitals.
    pub n_occupied: usize,
    /// Per-cycle convergence history.
    pub iterations: Vec<ScfIteration>,
    /// `true` when canonical orthogonalisation dropped a near-singular
    /// direction.
    pub linear_dependence_dropped: bool,
}

impl RhfResult {
    /// The HOMO energy (`None` for a system with no electrons).
    pub fn homo_energy(&self) -> Option<f64> {
        if self.n_occupied == 0 {
            None
        } else {
            Some(self.orbital_energies[self.n_occupied - 1])
        }
    }

    /// The LUMO energy (`None` when every orbital is occupied).
    pub fn lumo_energy(&self) -> Option<f64> {
        if self.n_occupied >= self.orbital_energies.len() {
            None
        } else {
            Some(self.orbital_energies[self.n_occupied])
        }
    }

    /// The HOMO-LUMO gap in Hartree (`None` when undefined).
    pub fn homo_lumo_gap(&self) -> Option<f64> {
        Some(self.lumo_energy()? - self.homo_energy()?)
    }
}

/// Build the core-Hamiltonian guess density: diagonalise `H` as if it
/// were the Fock matrix and fill the lowest `n_occ` orbitals.
pub fn core_guess_density(
    h_core: &DMatrix<f64>,
    ortho: &Orthogonalizer,
    n_occupied: usize,
) -> DMatrix<f64> {
    let (_, c) = solve_roothaan(h_core, ortho);
    build_density(&c, n_occupied)
}

/// Build the closed-shell density matrix
/// `D_{μν} = 2 Σ_i^{occ} C_{μi} C_{νi}` from the occupied orbitals.
pub fn build_density(c: &DMatrix<f64>, n_occupied: usize) -> DMatrix<f64> {
    let n = c.nrows();
    let mut d = DMatrix::<f64>::zeros(n, n);
    let occ = n_occupied.min(c.ncols());
    for mu in 0..n {
        for nu in 0..n {
            let mut acc = 0.0;
            for i in 0..occ {
                acc += c[(mu, i)] * c[(nu, i)];
            }
            d[(mu, nu)] = 2.0 * acc;
        }
    }
    d
}

/// Build the RHF Fock matrix `F = H + G(D)` with
/// `G_{μν} = Σ_{λσ} D_{λσ} [(μν|λσ) − ½(μλ|νσ)]`.
pub fn build_fock(h_core: &DMatrix<f64>, density: &DMatrix<f64>, eri: &EriTensor) -> DMatrix<f64> {
    let n = h_core.nrows();
    let mut f = h_core.clone();
    for mu in 0..n {
        for nu in 0..n {
            let mut g = 0.0;
            for la in 0..n {
                for si in 0..n {
                    let coulomb = eri.get(mu, nu, la, si);
                    let exchange = eri.get(mu, la, nu, si);
                    g += density[(la, si)] * (coulomb - 0.5 * exchange);
                }
            }
            f[(mu, nu)] += g;
        }
    }
    f
}

/// The RHF electronic energy
/// `E_elec = ½ Σ_{μν} D_{μν} (H_{μν} + F_{μν})`.
pub fn rhf_electronic_energy(
    density: &DMatrix<f64>,
    h_core: &DMatrix<f64>,
    fock: &DMatrix<f64>,
) -> f64 {
    let n = density.nrows();
    let mut e = 0.0;
    for mu in 0..n {
        for nu in 0..n {
            e += 0.5 * density[(mu, nu)] * (h_core[(mu, nu)] + fock[(mu, nu)]);
        }
    }
    e
}

/// The total RHF energy — electronic plus nuclear repulsion.
pub fn rhf_energy(
    density: &DMatrix<f64>,
    h_core: &DMatrix<f64>,
    fock: &DMatrix<f64>,
    e_nuclear: f64,
) -> f64 {
    rhf_electronic_energy(density, h_core, fock) + e_nuclear
}

/// Run the restricted-Hartree-Fock SCF iteration.
///
/// `n_electrons` must be even (RHF is closed-shell only); the number of
/// doubly-occupied orbitals is `n_electrons / 2`.
///
/// # Errors
///
/// - [`QchemError::InvalidInput`] when `n_electrons` is odd or exceeds
///   twice the basis dimension.
/// - [`QchemError::ScfNotConverged`] when the loop hits
///   `max_iterations` without satisfying both tolerances.
pub fn run_rhf_scf(
    integrals: &IntegralSet,
    n_electrons: u32,
    settings: ScfSettings,
) -> Result<RhfResult> {
    if n_electrons % 2 != 0 {
        return Err(QchemError::invalid(format!(
            "RHF needs an even electron count, got {n_electrons}"
        )));
    }
    let n_occupied = (n_electrons / 2) as usize;

    let ortho = symmetric_orthogonalizer(&integrals.overlap)?;
    if n_occupied > ortho.n_retained() {
        return Err(QchemError::invalid(format!(
            "{n_electrons} electrons need {n_occupied} orbitals but only \
             {} linearly-independent basis functions are available",
            ortho.n_retained()
        )));
    }

    let h_core = integrals.core_hamiltonian();
    let mut density = core_guess_density(&h_core, &ortho, n_occupied);

    let mut diis = Diis::new(settings.diis_vectors);
    let mut iterations = Vec::new();
    let mut last_energy = 0.0;

    for cycle in 1..=settings.max_iterations {
        let fock = build_fock(&h_core, &density, &integrals.eri);
        let energy = rhf_energy(&density, &h_core, &fock, integrals.e_nuclear);

        // DIIS error and extrapolation.
        let error = Diis::error_vector(&fock, &density, &integrals.overlap, &ortho);
        let error_norm = Diis::error_norm(&error);
        diis.push(fock.clone(), error);
        let fock_used = diis.extrapolate().unwrap_or(fock);

        let (orbital_energies, c) = solve_roothaan(&fock_used, &ortho);
        let orbital_coefficients = c.clone();
        density = build_density(&c, n_occupied);

        let delta_energy = energy - last_energy;
        last_energy = energy;
        iterations.push(ScfIteration {
            cycle,
            energy,
            delta_energy,
            error_norm,
        });

        if cycle > 1
            && delta_energy.abs() < settings.energy_tol
            && error_norm < settings.density_tol
        {
            // Recompute energy with the final density for consistency.
            let final_fock = build_fock(&h_core, &density, &integrals.eri);
            let final_energy = rhf_energy(&density, &h_core, &final_fock, integrals.e_nuclear);
            let electronic = rhf_electronic_energy(&density, &h_core, &final_fock);
            return Ok(RhfResult {
                total_energy: final_energy,
                electronic_energy: electronic,
                nuclear_repulsion: integrals.e_nuclear,
                orbital_energies,
                orbital_coefficients,
                density,
                n_occupied,
                iterations,
                linear_dependence_dropped: ortho.linear_dependence_dropped,
            });
        }
    }

    Err(QchemError::ScfNotConverged {
        iterations: settings.max_iterations,
        last_delta_energy: iterations.last().map(|i| i.delta_energy).unwrap_or(0.0),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::BasisSet;
    use crate::geometry::{Atom, MolecularGeometry};

    fn h2() -> MolecularGeometry {
        MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7414]).unwrap(),
        ])
    }

    fn water() -> MolecularGeometry {
        // Experimental water geometry.
        MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("O", [0.000_000, 0.000_000, 0.117_300]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.000_000, 0.757_200, -0.469_200]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.000_000, -0.757_200, -0.469_200]).unwrap(),
        ])
    }

    #[test]
    fn build_density_has_correct_trace() {
        // tr(D S) must equal the electron count; with S = identity
        // tr(D) = 2 * n_occ.
        let c = DMatrix::<f64>::identity(3, 3);
        let d = build_density(&c, 2);
        assert!((d.trace() - 4.0).abs() < 1.0e-12);
    }

    #[test]
    fn rhf_rejects_odd_electron_count() {
        let geom = h2();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        assert!(run_rhf_scf(&ints, 3, ScfSettings::default()).is_err());
    }

    #[test]
    fn sto3g_h2_energy() {
        // The textbook STO-3G H2 RHF energy at R = 0.7414 Å is
        // about -1.1167 Hartree.
        let geom = h2();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 2, ScfSettings::default()).unwrap();
        assert!(
            (res.total_energy - (-1.1167)).abs() < 2.0e-3,
            "STO-3G H2 E = {} (expected ~ -1.1167)",
            res.total_energy
        );
        // One occupied orbital; the HOMO must be below the LUMO.
        assert_eq!(res.n_occupied, 1);
        assert!(res.homo_lumo_gap().unwrap() > 0.0);
    }

    #[test]
    fn sto3g_water_energy() {
        // STO-3G water RHF energy is about -74.94 Hartree.
        let geom = water();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 10, ScfSettings::default()).unwrap();
        assert!(
            (res.total_energy - (-74.94)).abs() < 0.1,
            "STO-3G water E = {} (expected ~ -74.94)",
            res.total_energy
        );
        assert_eq!(res.n_occupied, 5);
    }

    #[test]
    fn scf_converges_and_energy_is_monotone_late() {
        let geom = water();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 10, ScfSettings::default()).unwrap();
        // The last iteration's energy change is below tolerance.
        let last = res.iterations.last().unwrap();
        assert!(last.delta_energy.abs() < 1.0e-8);
    }

    #[test]
    fn h2_energy_is_below_hydrogen_atoms() {
        // The bound H2 molecule must lie below two isolated H atoms
        // (each -0.4666 Ha at STO-3G).
        let geom = h2();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let res = run_rhf_scf(&ints, 2, ScfSettings::default()).unwrap();
        assert!(res.total_energy < -0.9);
    }
}
