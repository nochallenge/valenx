//! Unrestricted Hartree-Fock — the open-shell SCF.
//!
//! Unrestricted Hartree-Fock (UHF) gives the α and β electrons their
//! own spatial orbitals, so it describes radicals, cations and any
//! other open-shell system that restricted Hartree-Fock cannot.
//!
//! ## The Pople-Nesbet equations
//!
//! Two coupled Roothaan-style problems, one per spin:
//!
//! ```text
//! Fᵅ Cᵅ = S Cᵅ εᵅ      Fᵅ = H + J(Dᵅ + Dᵝ) − Kᵅ(Dᵅ)
//! Fᵝ Cᵝ = S Cᵝ εᵝ      Fᵝ = H + J(Dᵅ + Dᵝ) − Kᵝ(Dᵝ)
//! ```
//!
//! The Coulomb term sees the *total* density `Dᵅ + Dᵝ`; the exchange
//! term sees only the same-spin density. The two problems are coupled
//! through the shared Coulomb term and are iterated together to
//! self-consistency, accelerated by an independent DIIS history per
//! spin.
//!
//! For a closed-shell singlet UHF reproduces the RHF energy exactly
//! (`Dᵅ = Dᵝ`); the open-shell case is the reason it exists.

use crate::error::{QchemError, Result};
use crate::integrals::two_electron::EriTensor;
use crate::integrals::IntegralSet;
use crate::scf::diis::Diis;
use crate::scf::linalg::{solve_roothaan, symmetric_orthogonalizer};
use crate::scf::rhf::{ScfIteration, ScfSettings};
use nalgebra::{DMatrix, DVector};

/// The converged unrestricted-Hartree-Fock solution.
#[derive(Clone, Debug)]
pub struct UhfResult {
    /// Total UHF energy (electronic + nuclear), in Hartree.
    pub total_energy: f64,
    /// The electronic energy alone (Ha).
    pub electronic_energy: f64,
    /// Nuclear-repulsion energy (Ha).
    pub nuclear_repulsion: f64,
    /// α molecular-orbital energies, ascending (Ha).
    pub alpha_orbital_energies: DVector<f64>,
    /// β molecular-orbital energies, ascending (Ha).
    pub beta_orbital_energies: DVector<f64>,
    /// α molecular-orbital coefficients.
    pub alpha_coefficients: DMatrix<f64>,
    /// β molecular-orbital coefficients.
    pub beta_coefficients: DMatrix<f64>,
    /// The converged α density matrix.
    pub alpha_density: DMatrix<f64>,
    /// The converged β density matrix.
    pub beta_density: DMatrix<f64>,
    /// Number of α (then β) occupied orbitals.
    pub n_alpha: usize,
    /// Number of β occupied orbitals.
    pub n_beta: usize,
    /// `⟨S²⟩` expectation value — a spin-contamination diagnostic.
    pub s_squared: f64,
    /// Per-cycle convergence history.
    pub iterations: Vec<ScfIteration>,
}

impl UhfResult {
    /// The total spin-density matrix `Dᵅ − Dᵝ`.
    pub fn spin_density(&self) -> DMatrix<f64> {
        &self.alpha_density - &self.beta_density
    }

    /// The exact `⟨S²⟩` for the assigned spin state — `S(S+1)` with
    /// `S = (nᵅ − nᵝ)/2`. The difference from
    /// [`s_squared`](Self::s_squared) measures spin contamination.
    pub fn exact_s_squared(&self) -> f64 {
        let s = (self.n_alpha as f64 - self.n_beta as f64) / 2.0;
        s * (s + 1.0)
    }
}

/// Build a single-spin density `Dˢ_{μν} = Σ_i^{occ} C_{μi} C_{νi}` (no
/// factor of 2 — one spin only).
fn spin_density(c: &DMatrix<f64>, n_occupied: usize) -> DMatrix<f64> {
    let n = c.nrows();
    let mut d = DMatrix::<f64>::zeros(n, n);
    let occ = n_occupied.min(c.ncols());
    for mu in 0..n {
        for nu in 0..n {
            let mut acc = 0.0;
            for i in 0..occ {
                acc += c[(mu, i)] * c[(nu, i)];
            }
            d[(mu, nu)] = acc;
        }
    }
    d
}

/// Build a single-spin Fock matrix
/// `Fˢ = H + J(D_total) − Kˢ(Dˢ)`.
fn build_spin_fock(
    h_core: &DMatrix<f64>,
    d_total: &DMatrix<f64>,
    d_spin: &DMatrix<f64>,
    eri: &EriTensor,
) -> DMatrix<f64> {
    let n = h_core.nrows();
    let mut f = h_core.clone();
    for mu in 0..n {
        for nu in 0..n {
            let mut g = 0.0;
            for la in 0..n {
                for si in 0..n {
                    let coulomb = eri.get(mu, nu, la, si);
                    let exchange = eri.get(mu, la, nu, si);
                    // J sees total density, K sees same-spin density.
                    g += d_total[(la, si)] * coulomb - d_spin[(la, si)] * exchange;
                }
            }
            f[(mu, nu)] += g;
        }
    }
    f
}

/// The UHF electronic energy
/// `E = ½ Σ [D_total·H + Dᵅ·Fᵅ + Dᵝ·Fᵝ]`.
fn uhf_electronic_energy(
    d_alpha: &DMatrix<f64>,
    d_beta: &DMatrix<f64>,
    h_core: &DMatrix<f64>,
    f_alpha: &DMatrix<f64>,
    f_beta: &DMatrix<f64>,
) -> f64 {
    let n = h_core.nrows();
    let mut e = 0.0;
    for mu in 0..n {
        for nu in 0..n {
            let d_tot = d_alpha[(mu, nu)] + d_beta[(mu, nu)];
            e += 0.5
                * (d_tot * h_core[(mu, nu)]
                    + d_alpha[(mu, nu)] * f_alpha[(mu, nu)]
                    + d_beta[(mu, nu)] * f_beta[(mu, nu)]);
        }
    }
    e
}

/// Approximate `⟨S²⟩` for a UHF determinant (Szabo-Ostlund eq. 3.296):
/// `⟨S²⟩ = S_z(S_z+1) + nᵝ − Σ_ij |⟨αᵢ|βⱼ⟩|²`, the sum over occupied
/// α / β orbital pairs through the overlap metric `S`.
fn compute_s_squared(
    c_alpha: &DMatrix<f64>,
    c_beta: &DMatrix<f64>,
    overlap: &DMatrix<f64>,
    n_alpha: usize,
    n_beta: usize,
) -> f64 {
    let sz = (n_alpha as f64 - n_beta as f64) / 2.0;
    // Overlap of occupied α and β orbitals: Sᵅᵝ = Cᵅᵀ S Cᵝ.
    let s_ab = c_alpha.transpose() * overlap * c_beta;
    let mut contamination = 0.0;
    for i in 0..n_alpha.min(s_ab.nrows()) {
        for j in 0..n_beta.min(s_ab.ncols()) {
            contamination += s_ab[(i, j)] * s_ab[(i, j)];
        }
    }
    sz * (sz + 1.0) + n_beta as f64 - contamination
}

/// Run the unrestricted-Hartree-Fock SCF iteration.
///
/// `n_alpha` and `n_beta` are the α and β electron counts; for a
/// closed-shell system pass `n_alpha == n_beta` and UHF reproduces RHF.
///
/// # Errors
///
/// - [`QchemError::InvalidInput`] when either spin needs more orbitals
///   than the basis provides.
/// - [`QchemError::ScfNotConverged`] when the loop hits
///   `max_iterations` without satisfying both tolerances.
pub fn run_uhf_scf(
    integrals: &IntegralSet,
    n_alpha: u32,
    n_beta: u32,
    settings: ScfSettings,
) -> Result<UhfResult> {
    let na = n_alpha as usize;
    let nb = n_beta as usize;

    let ortho = symmetric_orthogonalizer(&integrals.overlap)?;
    if na.max(nb) > ortho.n_retained() {
        return Err(QchemError::invalid(format!(
            "UHF needs {} orbitals but only {} linearly-independent basis \
             functions are available",
            na.max(nb),
            ortho.n_retained()
        )));
    }

    let h_core = integrals.core_hamiltonian();

    // Core-Hamiltonian guess — same orbitals for both spins initially.
    let (_, c0) = solve_roothaan(&h_core, &ortho);
    let mut d_alpha = spin_density(&c0, na);
    let mut d_beta = spin_density(&c0, nb);

    let mut diis_a = Diis::new(settings.diis_vectors);
    let mut diis_b = Diis::new(settings.diis_vectors);
    let mut iterations = Vec::new();
    let mut last_energy = 0.0;

    for cycle in 1..=settings.max_iterations {
        let d_total = &d_alpha + &d_beta;
        let f_alpha = build_spin_fock(&h_core, &d_total, &d_alpha, &integrals.eri);
        let f_beta = build_spin_fock(&h_core, &d_total, &d_beta, &integrals.eri);

        let energy = uhf_electronic_energy(
            &d_alpha, &d_beta, &h_core, &f_alpha, &f_beta,
        ) + integrals.e_nuclear;

        // Independent DIIS per spin; the larger error drives convergence.
        let err_a =
            Diis::error_vector(&f_alpha, &d_alpha, &integrals.overlap, &ortho);
        let err_b =
            Diis::error_vector(&f_beta, &d_beta, &integrals.overlap, &ortho);
        let error_norm =
            Diis::error_norm(&err_a).max(Diis::error_norm(&err_b));
        diis_a.push(f_alpha.clone(), err_a);
        diis_b.push(f_beta.clone(), err_b);
        let fa_used = diis_a.extrapolate().unwrap_or(f_alpha);
        let fb_used = diis_b.extrapolate().unwrap_or(f_beta);

        let (ea, c_alpha) = solve_roothaan(&fa_used, &ortho);
        let (eb, c_beta) = solve_roothaan(&fb_used, &ortho);
        let ca = c_alpha.clone();
        let cb = c_beta.clone();
        d_alpha = spin_density(&c_alpha, na);
        d_beta = spin_density(&c_beta, nb);

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
            let d_total = &d_alpha + &d_beta;
            let f_alpha =
                build_spin_fock(&h_core, &d_total, &d_alpha, &integrals.eri);
            let f_beta =
                build_spin_fock(&h_core, &d_total, &d_beta, &integrals.eri);
            let electronic = uhf_electronic_energy(
                &d_alpha, &d_beta, &h_core, &f_alpha, &f_beta,
            );
            let s2 = compute_s_squared(&ca, &cb, &integrals.overlap, na, nb);
            return Ok(UhfResult {
                total_energy: electronic + integrals.e_nuclear,
                electronic_energy: electronic,
                nuclear_repulsion: integrals.e_nuclear,
                alpha_orbital_energies: ea,
                beta_orbital_energies: eb,
                alpha_coefficients: ca,
                beta_coefficients: cb,
                alpha_density: d_alpha,
                beta_density: d_beta,
                n_alpha: na,
                n_beta: nb,
                s_squared: s2,
                iterations,
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
    use crate::scf::rhf::run_rhf_scf;

    fn h2() -> MolecularGeometry {
        MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7414]).unwrap(),
        ])
    }

    #[test]
    fn uhf_matches_rhf_for_closed_shell() {
        // For a closed-shell singlet, UHF must reproduce the RHF energy.
        let geom = h2();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let rhf = run_rhf_scf(&ints, 2, ScfSettings::default()).unwrap();
        let uhf = run_uhf_scf(&ints, 1, 1, ScfSettings::default()).unwrap();
        assert!(
            (uhf.total_energy - rhf.total_energy).abs() < 1.0e-6,
            "UHF {} vs RHF {}",
            uhf.total_energy,
            rhf.total_energy
        );
    }

    #[test]
    fn closed_shell_uhf_has_zero_spin_contamination() {
        let geom = h2();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let uhf = run_uhf_scf(&ints, 1, 1, ScfSettings::default()).unwrap();
        // Singlet: ⟨S²⟩ should be ~0.
        assert!(uhf.s_squared.abs() < 1.0e-6, "S² = {}", uhf.s_squared);
        assert_eq!(uhf.exact_s_squared(), 0.0);
    }

    #[test]
    fn hydrogen_atom_doublet() {
        // A single H atom is a doublet: 1 alpha electron, 0 beta.
        let geom = MolecularGeometry::with_charge_multiplicity(
            vec![Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap()],
            0,
            2,
        );
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let uhf = run_uhf_scf(&ints, 1, 0, ScfSettings::default()).unwrap();
        // STO-3G hydrogen atom energy is about -0.4666 Hartree.
        assert!(
            (uhf.total_energy - (-0.4666)).abs() < 1.0e-3,
            "H atom E = {}",
            uhf.total_energy
        );
        // Doublet ⟨S²⟩ should be 0.75.
        assert!((uhf.exact_s_squared() - 0.75).abs() < 1.0e-12);
    }

    #[test]
    fn h2_cation_open_shell() {
        // H2+ : 1 electron. UHF should converge and be bound vs H + H+.
        let geom = MolecularGeometry::with_charge_multiplicity(
            h2().atoms.clone(),
            1,
            2,
        );
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let uhf = run_uhf_scf(&ints, 1, 0, ScfSettings::default()).unwrap();
        // H2+ STO-3G energy is roughly -0.58 Hartree.
        assert!(uhf.total_energy < -0.4, "H2+ E = {}", uhf.total_energy);
        assert_eq!(uhf.n_alpha, 1);
        assert_eq!(uhf.n_beta, 0);
    }
}
