//! The Kohn-Sham self-consistent-field loop.
//!
//! Kohn-Sham DFT replaces the Hartree-Fock exchange operator with the
//! exchange-correlation potential `V_xc`. The Kohn-Sham matrix is
//!
//! ```text
//! F^KS = H_core + J(D) + V_xc(D)   (+ a·K(D) for a hybrid)
//! ```
//!
//! where `H_core` is the one-electron core Hamiltonian, `J` the Coulomb
//! matrix, `V_xc` the exchange-correlation matrix integrated on the
//! molecular grid, and — for a hybrid functional — `a·K` a fraction of
//! the Hartree-Fock exchange. The matrix is diagonalised, a new density
//! built, and the cycle repeated to self-consistency, accelerated by
//! the crate's existing Pulay DIIS.
//!
//! ## The `V_xc` matrix
//!
//! For an LDA the XC matrix element is the grid quadrature
//!
//! ```text
//! V_xc[μν] = Σ_g w_g · v_xc(ρ_g) · φ_μ(r_g) φ_ν(r_g).
//! ```
//!
//! For a GGA the energy also depends on `∇ρ`, and the functional
//! derivative carries a divergence term `−∇·(2 ∂f/∂σ ∇ρ)`. Rather than
//! differentiate the potential, this module uses the standard
//! integration-by-parts form: the GGA contribution to the matrix
//! element is
//!
//! ```text
//! Σ_g w_g · 2 (∂f/∂σ) ∇ρ_g · ( φ_μ ∇φ_ν + φ_ν ∇φ_μ )
//! ```
//!
//! (with `∂f/∂σ = (∂f/∂|∇ρ|) / (2|∇ρ|)`), which moves the derivative
//! onto the basis functions — exactly what the grid already supplies
//! through [`GridDensity`]. The matrix is symmetrised, so it is a valid
//! Fock-like operator.
//!
//! ## The XC energy
//!
//! The total XC energy is the direct quadrature
//! `E_xc = Σ_g w_g ρ_g ε_xc(ρ_g, |∇ρ_g|)`. For a hybrid the
//! exact-exchange energy `−(a/2)·tr(D·K)`-style term is added through
//! the standard `½ tr(D F)` energy expression with the hybrid Fock
//! matrix.

use crate::basis::BasisSet;
use crate::dft::functional::Functional;
use crate::dft::grid::{GridDensity, GridQuality, MolecularGrid};
use crate::error::{QchemError, Result};
use crate::integrals::two_electron::EriTensor;
use crate::integrals::IntegralSet;
use crate::scf::diis::Diis;
use crate::scf::linalg::{solve_roothaan, symmetric_orthogonalizer};
use crate::scf::rhf::{build_density, ScfIteration, ScfSettings};
use nalgebra::{DMatrix, DVector};

/// The converged Kohn-Sham DFT solution.
#[derive(Clone, Debug)]
pub struct KsResult {
    /// Total KS-DFT energy (electronic + nuclear), in Hartree.
    pub total_energy: f64,
    /// The electronic energy alone (Ha).
    pub electronic_energy: f64,
    /// Nuclear-repulsion energy (Ha).
    pub nuclear_repulsion: f64,
    /// The exchange-correlation energy `E_xc` (Ha).
    pub xc_energy: f64,
    /// The functional used.
    pub functional: Functional,
    /// Kohn-Sham orbital energies `ε`, ascending (Ha).
    pub orbital_energies: DVector<f64>,
    /// Kohn-Sham orbital coefficients `C`.
    pub orbital_coefficients: DMatrix<f64>,
    /// The converged density matrix `D`.
    pub density: DMatrix<f64>,
    /// Number of doubly-occupied Kohn-Sham orbitals.
    pub n_occupied: usize,
    /// The number of electrons recovered by integrating `ρ` on the
    /// grid — a grid-quality diagnostic; should equal the true count.
    pub grid_electron_count: f64,
    /// Per-cycle convergence history.
    pub iterations: Vec<ScfIteration>,
}

impl KsResult {
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

/// Build the closed-shell Coulomb matrix `J_{μν} = Σ_{λσ} D_{λσ}(μν|λσ)`.
fn build_coulomb(density: &DMatrix<f64>, eri: &EriTensor) -> DMatrix<f64> {
    let n = density.nrows();
    let mut j = DMatrix::<f64>::zeros(n, n);
    for mu in 0..n {
        for nu in 0..n {
            let mut acc = 0.0;
            for la in 0..n {
                for si in 0..n {
                    acc += density[(la, si)] * eri.get(mu, nu, la, si);
                }
            }
            j[(mu, nu)] = acc;
        }
    }
    j
}

/// Build the closed-shell Hartree-Fock exchange matrix
/// `K_{μν} = Σ_{λσ} D_{λσ}(μλ|νσ)`.
///
/// Used only for a hybrid functional. The Kohn-Sham exchange term is
/// `−(a/2) K` — the `½` is the closed-shell exchange factor and `a` the
/// hybrid mixing fraction.
fn build_exchange(density: &DMatrix<f64>, eri: &EriTensor) -> DMatrix<f64> {
    let n = density.nrows();
    let mut k = DMatrix::<f64>::zeros(n, n);
    for mu in 0..n {
        for nu in 0..n {
            let mut acc = 0.0;
            for la in 0..n {
                for si in 0..n {
                    acc += density[(la, si)] * eri.get(mu, la, nu, si);
                }
            }
            k[(mu, nu)] = acc;
        }
    }
    k
}

/// The exchange-correlation matrix and energy from a grid quadrature.
struct XcResult {
    /// The XC matrix `V_xc` (`n × n`, symmetric).
    matrix: DMatrix<f64>,
    /// The XC energy `E_xc`.
    energy: f64,
}

/// Build the exchange-correlation matrix `V_xc` and energy `E_xc` by
/// integrating the functional on the molecular grid.
///
/// `gd` carries the density, gradient and cached basis values/gradients
/// for the *current* density matrix.
fn build_xc(
    grid: &MolecularGrid,
    gd: &GridDensity,
    n_basis: usize,
    functional: Functional,
) -> XcResult {
    let mut v_xc = DMatrix::<f64>::zeros(n_basis, n_basis);
    let mut e_xc = 0.0;

    for (pi, gp) in grid.points.iter().enumerate() {
        let rho = gd.rho[pi];
        if rho <= 1.0e-12 {
            continue;
        }
        let grad_norm = gd.grad_norm(pi);
        let xc = functional.evaluate(rho, grad_norm);
        let w = gp.weight;

        // E_xc = Σ w ρ ε_xc.
        e_xc += w * rho * xc.energy_density;

        let phi = &gd.phi[pi];
        let dphi = &gd.dphi[pi];
        let v_rho = xc.potential; // ∂f/∂ρ

        // GGA divergence term: 2 ∂f/∂σ ∇ρ, with ∂f/∂σ = (∂f/∂|∇ρ|)/(2|∇ρ|).
        // The vector that multiplies (φ_μ ∇φ_ν + φ_ν ∇φ_μ) is
        //   gga_vec = 2 ∂f/∂σ ∇ρ = (∂f/∂|∇ρ|/|∇ρ|) ∇ρ.
        let gga_vec = if grad_norm > 1.0e-10 && xc.gradient_potential != 0.0 {
            let factor = xc.gradient_potential / grad_norm;
            let g = gd.grad[pi];
            [factor * g[0], factor * g[1], factor * g[2]]
        } else {
            [0.0; 3]
        };

        // Accumulate the matrix element contribution.
        //   V_rho part:  w · v_rho · φ_μ φ_ν
        //   V_gga part:  w · gga_vec · (φ_μ ∇φ_ν + φ_ν ∇φ_μ)
        for mu in 0..n_basis {
            let phi_mu = phi[mu];
            let dphi_mu = dphi[mu];
            if phi_mu == 0.0 && dphi_mu == [0.0; 3] {
                continue;
            }
            for nu in 0..n_basis {
                let phi_nu = phi[nu];
                let dphi_nu = dphi[nu];
                let mut contrib = v_rho * phi_mu * phi_nu;
                // GGA: gga_vec · (φ_μ ∇φ_ν + φ_ν ∇φ_μ).
                contrib += gga_vec[0] * (phi_mu * dphi_nu[0] + phi_nu * dphi_mu[0]);
                contrib += gga_vec[1] * (phi_mu * dphi_nu[1] + phi_nu * dphi_mu[1]);
                contrib += gga_vec[2] * (phi_mu * dphi_nu[2] + phi_nu * dphi_mu[2]);
                v_xc[(mu, nu)] += w * contrib;
            }
        }
    }

    // Symmetrise to guarantee a valid Fock-like operator.
    let sym = 0.5 * (&v_xc + v_xc.transpose());
    XcResult {
        matrix: sym,
        energy: e_xc,
    }
}

/// Run the restricted Kohn-Sham SCF for a closed-shell molecule.
///
/// `n_electrons` must be even. `functional` selects the
/// exchange-correlation functional; `quality` the integration grid.
///
/// # Errors
///
/// - [`QchemError::InvalidInput`] when `n_electrons` is odd or exceeds
///   twice the basis dimension.
/// - [`QchemError::ScfNotConverged`] when the loop hits
///   `max_iterations` without satisfying both tolerances.
pub fn run_ks_scf(
    integrals: &IntegralSet,
    basis: &BasisSet,
    grid: &MolecularGrid,
    n_electrons: u32,
    functional: Functional,
    settings: ScfSettings,
) -> Result<KsResult> {
    if n_electrons % 2 != 0 {
        return Err(QchemError::invalid(format!(
            "restricted Kohn-Sham needs an even electron count, got \
             {n_electrons}"
        )));
    }
    let n_occupied = (n_electrons / 2) as usize;
    let n_basis = basis.n_functions();

    let ortho = symmetric_orthogonalizer(&integrals.overlap)?;
    if n_occupied > ortho.n_retained() {
        return Err(QchemError::invalid(format!(
            "{n_electrons} electrons need {n_occupied} orbitals but only \
             {} linearly-independent basis functions are available",
            ortho.n_retained()
        )));
    }

    let h_core = integrals.core_hamiltonian();
    let a_exact = functional.exact_exchange_fraction();

    // Core-Hamiltonian density guess.
    let (_, c0) = solve_roothaan(&h_core, &ortho);
    let mut density = build_density(&c0, n_occupied);

    let mut diis = Diis::new(settings.diis_vectors);
    let mut iterations = Vec::new();
    let mut last_energy = 0.0;
    let mut final_xc_energy = 0.0;
    let mut final_grid_n = 0.0;

    for cycle in 1..=settings.max_iterations {
        // Coulomb matrix.
        let j = build_coulomb(&density, &integrals.eri);
        // XC matrix + energy on the grid (uses the current density).
        let gd = GridDensity::evaluate(grid, basis, &density);
        let xc = build_xc(grid, &gd, n_basis, functional);
        let grid_n = gd.integrate_electrons(grid);

        // Kohn-Sham matrix F = H + J + V_xc (+ a·exchange for a hybrid).
        let mut fock = &h_core + &j + &xc.matrix;
        let k = if a_exact > 0.0 {
            let k = build_exchange(&density, &integrals.eri);
            // Hybrid: subtract a·(½K) — the closed-shell exchange.
            fock -= (a_exact * 0.5) * &k;
            Some(k)
        } else {
            None
        };

        // KS-DFT energy. The electronic energy is
        //   E = ½ Σ D (H + H + J)  +  E_xc  (− ¼ a Σ D K for a hybrid).
        // i.e. the one-electron + Coulomb pieces from the standard
        // ½tr[D(H+F_J)] form, plus E_xc, plus the hybrid exact-exchange
        // energy −(a/4) tr(D K).
        let mut e_elec = 0.0;
        for mu in 0..n_basis {
            for nu in 0..n_basis {
                e_elec += 0.5 * density[(mu, nu)] * (2.0 * h_core[(mu, nu)] + j[(mu, nu)]);
            }
        }
        e_elec += xc.energy;
        if let Some(k) = &k {
            let mut e_x = 0.0;
            for mu in 0..n_basis {
                for nu in 0..n_basis {
                    e_x += density[(mu, nu)] * k[(mu, nu)];
                }
            }
            // Closed-shell exact exchange energy is −¼ tr(D K); scale
            // by the hybrid fraction a.
            e_elec -= 0.25 * a_exact * e_x;
        }
        let energy = e_elec + integrals.e_nuclear;

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
        final_xc_energy = xc.energy;
        final_grid_n = grid_n;
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
            // Recompute everything with the final density.
            let j = build_coulomb(&density, &integrals.eri);
            let gd = GridDensity::evaluate(grid, basis, &density);
            let xc = build_xc(grid, &gd, n_basis, functional);
            let grid_n = gd.integrate_electrons(grid);
            let mut e_elec = 0.0;
            for mu in 0..n_basis {
                for nu in 0..n_basis {
                    e_elec += 0.5 * density[(mu, nu)] * (2.0 * h_core[(mu, nu)] + j[(mu, nu)]);
                }
            }
            e_elec += xc.energy;
            if a_exact > 0.0 {
                let k = build_exchange(&density, &integrals.eri);
                let mut e_x = 0.0;
                for mu in 0..n_basis {
                    for nu in 0..n_basis {
                        e_x += density[(mu, nu)] * k[(mu, nu)];
                    }
                }
                e_elec -= 0.25 * a_exact * e_x;
            }
            return Ok(KsResult {
                total_energy: e_elec + integrals.e_nuclear,
                electronic_energy: e_elec,
                nuclear_repulsion: integrals.e_nuclear,
                xc_energy: xc.energy,
                functional,
                orbital_energies,
                orbital_coefficients,
                density,
                n_occupied,
                grid_electron_count: grid_n,
                iterations,
            });
        }
    }

    let _ = (final_xc_energy, final_grid_n);
    Err(QchemError::ScfNotConverged {
        iterations: settings.max_iterations,
        last_delta_energy: iterations.last().map(|i| i.delta_energy).unwrap_or(0.0),
    })
}

/// Convenience: build the grid and run the restricted Kohn-Sham SCF in
/// one call.
///
/// # Errors
///
/// Propagates every error from [`run_ks_scf`].
pub fn run_ks(
    integrals: &IntegralSet,
    basis: &BasisSet,
    geometry: &crate::geometry::MolecularGeometry,
    n_electrons: u32,
    functional: Functional,
    quality: GridQuality,
    settings: ScfSettings,
) -> Result<KsResult> {
    let grid = MolecularGrid::build(geometry, quality);
    run_ks_scf(integrals, basis, &grid, n_electrons, functional, settings)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Atom, MolecularGeometry};

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
    fn lda_scf_converges_for_h2() {
        let geom = h2();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let r = run_ks(
            &ints,
            &basis,
            &geom,
            2,
            Functional::Lda,
            GridQuality::Medium,
            ScfSettings::default(),
        )
        .unwrap();
        // The SCF converged.
        let last = r.iterations.last().unwrap();
        assert!(last.delta_energy.abs() < 1.0e-7);
        // The grid recovers 2 electrons.
        assert!((r.grid_electron_count - 2.0).abs() < 1.0e-2);
        // One occupied orbital, positive HOMO-LUMO gap.
        assert_eq!(r.n_occupied, 1);
        assert!(r.homo_lumo_gap().unwrap() > 0.0);
    }

    #[test]
    fn pbe_scf_converges_for_h2() {
        let geom = h2();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let r = run_ks(
            &ints,
            &basis,
            &geom,
            2,
            Functional::Pbe,
            GridQuality::Medium,
            ScfSettings::default(),
        )
        .unwrap();
        assert!(r.iterations.last().unwrap().delta_energy.abs() < 1.0e-7);
        assert!((r.grid_electron_count - 2.0).abs() < 1.0e-2);
        // XC energy is negative.
        assert!(r.xc_energy < 0.0);
    }

    #[test]
    fn b3lyp_scf_converges_for_h2() {
        let geom = h2();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let r = run_ks(
            &ints,
            &basis,
            &geom,
            2,
            Functional::B3lyp,
            GridQuality::Medium,
            ScfSettings::default(),
        )
        .unwrap();
        assert!(r.iterations.last().unwrap().delta_energy.abs() < 1.0e-7);
        assert!((r.grid_electron_count - 2.0).abs() < 1.0e-2);
        assert_eq!(r.functional, Functional::B3lyp);
    }

    #[test]
    fn ks_rejects_odd_electron_count() {
        let geom = h2();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let grid = MolecularGrid::build(&geom, GridQuality::Coarse);
        assert!(run_ks_scf(
            &ints,
            &basis,
            &grid,
            3,
            Functional::Lda,
            ScfSettings::default()
        )
        .is_err());
    }

    #[test]
    fn lda_water_converges_and_grid_integrates() {
        let geom = water();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let r = run_ks(
            &ints,
            &basis,
            &geom,
            10,
            Functional::Lda,
            GridQuality::Medium,
            ScfSettings::default(),
        )
        .unwrap();
        // 10 electrons recovered by the grid.
        assert!(
            (r.grid_electron_count - 10.0).abs() < 5.0e-2,
            "grid N = {}",
            r.grid_electron_count
        );
        assert_eq!(r.n_occupied, 5);
    }

    /// The KS-DFT total energy is the sum of its documented pieces —
    /// one-electron + Coulomb + XC + nuclear repulsion — and the XC
    /// energy is a substantial negative fraction of the total. (DFT
    /// total energy is *not* universally below Hartree-Fock at a fixed
    /// finite basis: at the minimal STO-3G basis, LDA exchange is less
    /// accurate than exact exchange and the LDA water energy sits
    /// *above* HF/STO-3G — the over-binding only shows at a complete
    /// basis. So this test checks internal consistency, not an
    /// HF-vs-DFT ordering.)
    #[test]
    fn dft_energy_decomposition_is_consistent() {
        let geom = water();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let dft = run_ks(
            &ints,
            &basis,
            &geom,
            10,
            Functional::Lda,
            GridQuality::Medium,
            ScfSettings::default(),
        )
        .unwrap();
        // total = electronic + nuclear repulsion.
        assert!(
            (dft.total_energy - (dft.electronic_energy + dft.nuclear_repulsion)).abs() < 1.0e-10
        );
        // The XC energy is negative and a chemically sensible
        // magnitude for water (exchange-correlation is roughly
        // −8 to −9 Ha at this level).
        assert!(dft.xc_energy < 0.0);
        assert!(
            dft.xc_energy > -12.0 && dft.xc_energy < -5.0,
            "E_xc = {} out of physical band",
            dft.xc_energy
        );
    }
}
