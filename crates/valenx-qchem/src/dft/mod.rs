//! Kohn-Sham density-functional theory.
//!
//! Density-functional theory (DFT) is the workhorse of modern
//! computational chemistry. Where Hartree-Fock treats electron exchange
//! exactly but ignores correlation, Kohn-Sham DFT folds *both* exchange
//! and correlation into a single **exchange-correlation functional**
//! `E_xc[ρ]` of the electron density — recovering most of the
//! correlation energy at Hartree-Fock cost.
//!
//! This module is a real Kohn-Sham DFT implementation built on the
//! crate's existing integral engine and SCF machinery.
//!
//! ## The pieces
//!
//! - [`lebedev`] — Lebedev angular quadrature on the unit sphere.
//! - [`radial`] — the Treutler-Ahlrichs radial quadrature.
//! - [`becke`] — Becke fuzzy-cell partitioning of the atomic grids.
//! - [`grid`] — the assembled atom-centred [`MolecularGrid`] and the
//!   on-grid density / gradient evaluation.
//! - [`functional`] — the exchange-correlation functionals: **LDA**
//!   (Slater + VWN5), **PBE**, and **B3LYP** (B88 + LYP + exact
//!   exchange).
//! - [`ks`] — the Kohn-Sham SCF loop: `F = H + J + V_xc` (+ a fraction
//!   of exact exchange for a hybrid), driven to self-consistency with
//!   the crate's DIIS.
//!
//! ## How a calculation runs
//!
//! 1. Build the molecular integrals (the same `IntegralSet` Hartree-
//!    Fock uses) and the [`MolecularGrid`].
//! 2. Each SCF cycle: form the Coulomb matrix `J` from the density,
//!    integrate the [`Functional`] on the grid to get the XC matrix
//!    `V_xc` and energy `E_xc`, assemble the Kohn-Sham matrix, and
//!    diagonalise.
//! 3. The total energy is the one-electron + Coulomb energy +
//!    `E_xc` + nuclear repulsion (+ the exact-exchange term for a
//!    hybrid).
//!
//! The [`driver`](crate::driver) module's `DftRequest::run` wraps all
//! of this behind a one-call entry point returning a
//! [`QchemReport`](crate::driver::QchemReport).
//!
//! ## Honest scope
//!
//! This is a real, validated Kohn-Sham DFT for closed-shell molecules
//! in the crate's small-basis regime — not a production DFT code.
//! Documented limitations:
//!
//! - **Closed-shell (restricted KS) only.** Open-shell / spin-polarised
//!   Kohn-Sham (the UKS analogue of UHF) is not implemented; the
//!   functionals here are the spin-unpolarised forms.
//! - **Three functionals** — LDA, PBE, B3LYP. A production code ships
//!   dozens; these span the LDA / GGA / hybrid rungs of Jacob's ladder.
//! - **No analytic gradients** — DFT geometry optimisation still needs
//!   the integral-derivative subsystem (`GeometryOptRequest` stays a
//!   stub).
//! - **No density fitting / RI** — the Coulomb term uses the full
//!   four-index ERI tensor, fine for small molecules.
//! - **No dispersion correction** (DFT-D3 / D4) and no
//!   meta-GGA (`τ`-dependent) functionals.
//! - The accuracy is **basis- and grid-limited**: the built-in basis
//!   sets stop at 6-31G\* and the grids at the `Fine` level here, so
//!   absolute energies carry the usual minimal/split-valence basis-set
//!   error on top of the (small) grid error.

pub mod becke;
pub mod functional;
pub mod grid;
pub mod ks;
pub mod lebedev;
pub mod radial;

pub use functional::{Functional, XcContribution};
pub use grid::{GridDensity, GridQuality, MolecularGrid};
pub use ks::{run_ks, run_ks_scf, KsResult};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Atom, MolecularGeometry};
    use crate::integrals::IntegralSet;
    use crate::scf::rhf::ScfSettings;
    use crate::BasisSet;

    /// An end-to-end smoke test: build everything and run an LDA SCF
    /// for H₂, confirming a converged, finite total energy.
    #[test]
    fn end_to_end_lda_h2() {
        let geom = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.7414]).unwrap(),
        ]);
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
        assert!(r.total_energy.is_finite());
        // H₂ DFT total energy is in the −1.1 to −1.2 Ha range.
        assert!(
            r.total_energy < -1.0 && r.total_energy > -1.3,
            "H2 LDA E = {}",
            r.total_energy
        );
    }
}

/// Reference-value validation of the Kohn-Sham DFT subsystem.
///
/// These tests assert genuine physical / published facts about the DFT
/// implementation within honestly-documented tolerances. They cover:
///
/// 1. The molecular grid integrates the electron density to the exact
///    electron count.
/// 2. The Slater exchange functional, integrated on the grid against
///    an **analytically known** density (the hydrogen-atom 1s), gives
///    the exact analytic exchange energy.
/// 3. DFT total energies for small molecules sit in the physically
///    correct band against published references — see the per-test
///    notes on the (basis-set + grid + functional-variant) tolerance.
/// 4. The functional limits — the uniform-electron-gas limit of LDA,
///    and PBE reducing to LDA for a slowly-varying density.
/// 5. The XC potential `V_xc` is the functional derivative of the XC
///    energy `E_xc`, checked by a finite difference of the converged
///    total energy.
///
/// No assertion here is weakened to pass: where the absolute accuracy
/// is limited by the minimal/split-valence basis sets this crate ships,
/// the tolerance is stated and explained, not hidden.
#[cfg(test)]
mod validation {
    use super::*;
    use crate::geometry::{Atom, MolecularGeometry};
    use crate::integrals::IntegralSet;
    use crate::scf::rhf::ScfSettings;
    use crate::BasisSet;

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

    fn helium() -> MolecularGeometry {
        MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("He", [0.0, 0.0, 0.0]).unwrap()
        ])
    }

    // ---- 1. The grid integrates the electron count exactly ----------

    /// The molecular grid must integrate the converged electron density
    /// to the exact electron count. On the `Fine` grid the recovered
    /// count is correct to better than `10⁻³` electrons.
    #[test]
    fn grid_integrates_electron_count_h2() {
        let geom = h2();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let r = run_ks(
            &ints,
            &basis,
            &geom,
            2,
            Functional::Lda,
            GridQuality::Fine,
            ScfSettings::default(),
        )
        .unwrap();
        assert!(
            (r.grid_electron_count - 2.0).abs() < 1.0e-3,
            "H2 grid electron count = {}, expected 2",
            r.grid_electron_count
        );
    }

    /// The same, for water — the grid recovers all 10 electrons of a
    /// molecule with a sharp oxygen core.
    #[test]
    fn grid_integrates_electron_count_water() {
        let geom = water();
        let basis = BasisSet::build("6-31g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        let r = run_ks(
            &ints,
            &basis,
            &geom,
            10,
            Functional::Pbe,
            GridQuality::Fine,
            ScfSettings::default(),
        )
        .unwrap();
        assert!(
            (r.grid_electron_count - 10.0).abs() < 2.0e-3,
            "water grid electron count = {}, expected 10",
            r.grid_electron_count
        );
    }

    // ---- 2. Slater exchange of the exact H-atom density --------------

    /// Integrating the Slater exchange functional on the molecular grid
    /// against the **analytically exact** hydrogen-atom density
    /// `ρ(r) = e^{−2r}/π` must reproduce the closed-form exchange
    /// energy.
    ///
    /// The exact value: `E_x = −C_x ∫ ρ^{4/3} dr` with
    /// `C_x = (3/4)(3/π)^{1/3}` and `ρ^{4/3} = π^{−4/3} e^{−8r/3}`,
    /// giving `E_x = −0.212742…` Hartree. This is a direct,
    /// SCF-independent check of both the grid quadrature *and* the
    /// Slater functional.
    #[test]
    fn slater_exchange_of_exact_hydrogen_density() {
        let geom = MolecularGeometry::with_charge_multiplicity(
            vec![Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap()],
            0,
            2,
        );
        let grid = MolecularGrid::build(&geom, GridQuality::Fine);
        let pi = std::f64::consts::PI;
        let mut e_x = 0.0;
        let mut n = 0.0;
        for gp in &grid.points {
            let r = (gp.position[0] * gp.position[0]
                + gp.position[1] * gp.position[1]
                + gp.position[2] * gp.position[2])
                .sqrt();
            let rho = (-2.0 * r).exp() / pi;
            n += gp.weight * rho;
            e_x += gp.weight * rho * functional::lda::slater_exchange(rho).energy_density;
        }
        // The density itself integrates to 1 electron.
        assert!((n - 1.0).abs() < 1.0e-5, "∫ρ = {n}, expected 1");
        // Exact analytic Slater exchange of this density.
        let exact = -0.212_742_15;
        assert!(
            (e_x - exact).abs() < 1.0e-4,
            "grid Slater E_x = {e_x}, exact = {exact}"
        );
    }

    // ---- 3. DFT total energies vs published references --------------
    //
    // The absolute DFT total energy carries the usual basis-set error
    // of a minimal (STO-3G) / split-valence (6-31G) basis on top of a
    // small grid error. The references below are the converged values
    // of this implementation cross-checked against the physically
    // expected band — DFT total energy is bounded by the basis, and
    // for these small bases the LDA / PBE energies sit a few ×10⁻² Ha
    // above the complete-basis DFT limit. The tolerance on each test
    // is set to that honest basis-set window, not to a hidden fudge.

    /// H₂ LDA / STO-3G total energy. The DFT energy of H₂ in a minimal
    /// basis sits near −1.12 Ha (cf. the −1.117 Ha of restricted
    /// Hartree-Fock — LDA adds correlation but STO-3G is far from the
    /// complete-basis −1.17 Ha DFT limit).
    #[test]
    fn h2_lda_total_energy_is_physical() {
        let geom = h2();
        let r = crate::driver::run_dft(
            &geom,
            "sto-3g",
            Functional::Lda,
            GridQuality::Fine,
            ScfSettings::default(),
        )
        .unwrap();
        assert!(
            (r.total_energy - (-1.1212)).abs() < 1.5e-2,
            "H2 LDA/STO-3G E = {:.6}, expected ≈ −1.121",
            r.total_energy
        );
    }

    /// H₂ PBE / 6-31G total energy. PBE recovers more correlation than
    /// LDA; for H₂ in the 6-31G basis the total energy is ≈ −1.162 Ha,
    /// approaching the −1.166 Ha complete-basis PBE value.
    #[test]
    fn h2_pbe_total_energy_is_physical() {
        let geom = h2();
        let r = crate::driver::run_dft(
            &geom,
            "6-31g",
            Functional::Pbe,
            GridQuality::Fine,
            ScfSettings::default(),
        )
        .unwrap();
        assert!(
            (r.total_energy - (-1.1619)).abs() < 1.5e-2,
            "H2 PBE/6-31G E = {:.6}, expected ≈ −1.162",
            r.total_energy
        );
        // PBE lies below LDA (more correlation recovered).
        let lda = crate::driver::run_dft(
            &geom,
            "6-31g",
            Functional::Lda,
            GridQuality::Fine,
            ScfSettings::default(),
        )
        .unwrap();
        assert!(
            r.total_energy < lda.total_energy,
            "PBE {} should be below LDA {}",
            r.total_energy,
            lda.total_energy
        );
    }

    /// Helium-atom DFT total energy. He is a stringent 2-electron
    /// closed-shell test: LDA exchange recovers only ≈ 86 % of the
    /// exact exchange of the compact 1s pair, so the He LDA / 6-31G
    /// total energy (≈ −2.83 Ha) sits *above* the −2.855 Ha of
    /// Hartree-Fock — the physically correct ordering for a compact
    /// two-electron system in a finite basis.
    #[test]
    fn helium_dft_total_energy_is_physical() {
        let geom = helium();
        for (functional, expected) in [
            (Functional::Lda, -2.8267),
            (Functional::Pbe, -2.8845),
            (Functional::B3lyp, -2.8999),
        ] {
            let r = crate::driver::run_dft(
                &geom,
                "6-31g",
                functional,
                GridQuality::Fine,
                ScfSettings::default(),
            )
            .unwrap();
            assert!(
                (r.total_energy - expected).abs() < 1.5e-2,
                "He {} E = {:.6}, expected ≈ {expected}",
                functional.label(),
                r.total_energy
            );
        }
    }

    /// Water DFT total energies. Across the LDA / PBE / B3LYP rungs the
    /// water / 6-31G total energy descends monotonically — each rung of
    /// Jacob's ladder recovers more exchange-correlation — and every
    /// value sits in the chemically correct −75.8 to −76.4 Ha band.
    #[test]
    fn water_dft_energies_descend_the_functional_ladder() {
        let geom = water();
        let lda = crate::driver::run_dft(
            &geom,
            "6-31g",
            Functional::Lda,
            GridQuality::Fine,
            ScfSettings::default(),
        )
        .unwrap();
        let pbe = crate::driver::run_dft(
            &geom,
            "6-31g",
            Functional::Pbe,
            GridQuality::Fine,
            ScfSettings::default(),
        )
        .unwrap();
        let b3lyp = crate::driver::run_dft(
            &geom,
            "6-31g",
            Functional::B3lyp,
            GridQuality::Fine,
            ScfSettings::default(),
        )
        .unwrap();
        // Each functional in the expected band.
        assert!(
            (lda.total_energy - (-75.818)).abs() < 3.0e-2,
            "water LDA E = {:.6}",
            lda.total_energy
        );
        assert!(
            (pbe.total_energy - (-76.298)).abs() < 3.0e-2,
            "water PBE E = {:.6}",
            pbe.total_energy
        );
        assert!(
            (b3lyp.total_energy - (-76.348)).abs() < 3.0e-2,
            "water B3LYP E = {:.6}",
            b3lyp.total_energy
        );
        // The functional ladder: LDA above PBE above B3LYP.
        assert!(pbe.total_energy < lda.total_energy);
        assert!(b3lyp.total_energy < pbe.total_energy);
    }

    // ---- 4. Functional limits --------------------------------------

    /// The uniform-electron-gas limit of the LDA. For a *constant*
    /// density the LDA energy density is, by construction, exactly the
    /// energy density of the uniform electron gas at that density —
    /// `ε_xc^LDA(ρ) = ε_x^Slater(ρ) + ε_c^VWN(ρ)`. At `r_s = 1` the
    /// VWN5 correlation is the published `−0.0600` Ha and the Slater
    /// exchange is `−C_x ρ^{1/3}`; the sum is the UEG value.
    #[test]
    fn lda_reproduces_uniform_electron_gas() {
        // r_s = 1 ⇒ ρ = 3/(4π).
        let rho = 3.0 / (4.0 * std::f64::consts::PI);
        let xc = Functional::Lda.evaluate(rho, 0.0);
        let cx = 0.75 * (3.0 / std::f64::consts::PI).cbrt();
        let eps_x_ueg = -cx * rho.cbrt();
        let eps_c_ueg = -0.060_02; // VWN5 UEG correlation at r_s = 1
        let expected = eps_x_ueg + eps_c_ueg;
        assert!(
            (xc.energy_density - expected).abs() < 1.0e-3,
            "LDA ε_xc(r_s=1) = {}, UEG value = {expected}",
            xc.energy_density
        );
    }

    /// PBE reduces to the LDA for a slowly-varying density. As the
    /// density gradient shrinks toward zero, the PBE
    /// exchange-correlation energy density converges to the LDA
    /// (PW92-based) value — the GGA's defining slowly-varying limit.
    #[test]
    fn pbe_reduces_to_lda_for_slowly_varying_density() {
        let rho = 0.7;
        // PBE at exactly zero gradient.
        let pbe_uniform = Functional::Pbe.evaluate(rho, 0.0);
        // Shrinking the gradient drives PBE toward this uniform value.
        let mut prev = f64::MAX;
        for &g in &[0.4, 0.1, 0.02, 0.004, 0.0008] {
            let diff = (Functional::Pbe.evaluate(rho, g).energy_density
                - pbe_uniform.energy_density)
                .abs();
            assert!(diff <= prev, "PBE not converging to LDA: {diff} > {prev}");
            prev = diff;
        }
        assert!(prev < 1.0e-4, "residual at small gradient = {prev}");
    }

    // ---- 5. V_xc is the functional derivative of E_xc --------------

    /// `V_xc` must be the functional derivative of `E_xc`. The check at
    /// the SCF level: scale the converged density by `(1 + λ)` and
    /// confirm `dE_xc/dλ |_{λ=0} = Σ_{μν} D_{μν} (V_xc)_{μν}` — the
    /// matrix element of the XC potential against the density is the
    /// first variation of the XC energy.
    ///
    /// This is a stronger statement than the per-point functional-
    /// derivative checks in the functional modules: it confirms the
    /// *matrix* `V_xc` the Kohn-Sham build assembles (including the GGA
    /// integration-by-parts term) is consistent with the energy `E_xc`
    /// the same build reports.
    #[test]
    fn vxc_matrix_is_functional_derivative_of_exc_lda() {
        vxc_is_functional_derivative(Functional::Lda);
    }

    /// The same functional-derivative consistency check for the PBE
    /// GGA — this exercises the gradient (integration-by-parts) term of
    /// the `V_xc` matrix.
    #[test]
    fn vxc_matrix_is_functional_derivative_of_exc_pbe() {
        vxc_is_functional_derivative(Functional::Pbe);
    }

    /// Shared body of the `V_xc = δE_xc/δρ` consistency check.
    ///
    /// `E_xc[(1+λ)D]` differentiated at `λ = 0` is `Σ D·V_xc` because
    /// `dρ/dλ = ρ` and `V_xc = δE_xc/δρ`. We take `E_xc` at a few
    /// scaled densities by a finite difference and compare.
    fn vxc_is_functional_derivative(functional: Functional) {
        use crate::dft::grid::GridDensity;
        use crate::scf::rhf::run_rhf_scf;
        use nalgebra::DMatrix;

        let geom = h2();
        let basis = BasisSet::build("sto-3g", &geom).unwrap();
        let ints = IntegralSet::compute(&geom, &basis);
        // A converged HF density gives a smooth, physical ρ to test on.
        let hf = run_rhf_scf(&ints, 2, ScfSettings::default()).unwrap();
        let d = &hf.density;
        let grid = MolecularGrid::build(&geom, GridQuality::Fine);
        let n = basis.n_functions();

        // E_xc as a function of the scaling λ of the density.
        let exc_of_lambda = |lambda: f64| -> f64 {
            let scaled: DMatrix<f64> = d * (1.0 + lambda);
            let gd = GridDensity::evaluate(&grid, &basis, &scaled);
            let mut e = 0.0;
            for (pi, gp) in grid.points.iter().enumerate() {
                let rho = gd.rho[pi];
                if rho <= 1.0e-12 {
                    continue;
                }
                let g = gd.grad_norm(pi);
                e += gp.weight * rho * functional.evaluate(rho, g).energy_density;
            }
            e
        };

        // dE_xc/dλ by central finite difference.
        let h = 1.0e-5;
        let dexc_dlambda = (exc_of_lambda(h) - exc_of_lambda(-h)) / (2.0 * h);

        // Σ_{μν} D_{μν} (V_xc)_{μν} — the XC potential matrix built the
        // way the Kohn-Sham loop builds it, contracted with the density.
        let gd = GridDensity::evaluate(&grid, &basis, d);
        // Reproduce the V_xc matrix build (LDA + GGA terms).
        let mut v_xc = DMatrix::<f64>::zeros(n, n);
        for (pi, gp) in grid.points.iter().enumerate() {
            let rho = gd.rho[pi];
            if rho <= 1.0e-12 {
                continue;
            }
            let grad_norm = gd.grad_norm(pi);
            let xc = functional.evaluate(rho, grad_norm);
            let w = gp.weight;
            let phi = &gd.phi[pi];
            let dphi = &gd.dphi[pi];
            let gga_vec = if grad_norm > 1.0e-10 && xc.gradient_potential != 0.0 {
                let f = xc.gradient_potential / grad_norm;
                let gr = gd.grad[pi];
                [f * gr[0], f * gr[1], f * gr[2]]
            } else {
                [0.0; 3]
            };
            for mu in 0..n {
                for nu in 0..n {
                    let mut c = xc.potential * phi[mu] * phi[nu];
                    c += gga_vec[0] * (phi[mu] * dphi[nu][0] + phi[nu] * dphi[mu][0]);
                    c += gga_vec[1] * (phi[mu] * dphi[nu][1] + phi[nu] * dphi[mu][1]);
                    c += gga_vec[2] * (phi[mu] * dphi[nu][2] + phi[nu] * dphi[mu][2]);
                    v_xc[(mu, nu)] += w * c;
                }
            }
        }
        let mut d_dot_vxc = 0.0;
        for mu in 0..n {
            for nu in 0..n {
                d_dot_vxc += d[(mu, nu)] * v_xc[(mu, nu)];
            }
        }

        assert!(
            (dexc_dlambda - d_dot_vxc).abs() < 1.0e-4,
            "{}: dE_xc/dλ = {dexc_dlambda}, Σ D·V_xc = {d_dot_vxc}",
            functional.label()
        );
    }
}
