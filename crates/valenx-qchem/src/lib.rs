//! # valenx-qchem — quantum-chemistry core
//!
//! A native-Rust restricted /
//! unrestricted Hartree-Fock program for small molecules — the
//! ab-initio core that Psi4, NWChem, PySCF and GAMESS-US provide, in
//! the small-system regime. Pure algorithms, no external processes.
//!
//! ## What it does
//!
//! - **Geometry & basis** — a [`geometry::MolecularGeometry`] (atoms,
//!   Cartesian coordinates, charge, multiplicity, xyz I/O) and a
//!   [`basis`] Gaussian-basis data model with built-in STO-3G, 3-21G,
//!   6-31G and 6-31G\* sets for H–Ne.
//! - **Integrals** ([`integrals`]) — overlap, kinetic-energy,
//!   nuclear-attraction and dipole one-electron integrals plus the full
//!   electron-repulsion tensor, all via the McMurchie-Davidson
//!   recursion and a series / asymptotic Boys function; the
//!   nuclear-repulsion energy.
//! - **SCF** ([`scf`]) — the core Hamiltonian, Löwdin symmetric
//!   orthogonalisation, a core-Hamiltonian density guess, the
//!   Roothaan-Hall restricted-Hartree-Fock loop with Pulay DIIS
//!   acceleration, and an unrestricted-Hartree-Fock loop for open-shell
//!   systems.
//! - **Properties** ([`properties`]) — Mulliken and Löwdin population
//!   analysis, the molecular dipole moment, molecular-orbital energies
//!   with the HOMO-LUMO gap, orbital evaluation on a 3D grid, and
//!   Mayer / Wiberg bond orders.
//! - **DFT** ([`dft`]) — real Kohn-Sham density-functional theory: an
//!   atom-centred molecular integration grid (Treutler-Ahlrichs radial
//!   × Lebedev angular, Becke fuzzy-cell partitioning), the **LDA**
//!   (Slater + VWN5), **PBE** and **B3LYP** exchange-correlation
//!   functionals, and a Kohn-Sham SCF loop.
//! - **Correlation & semi-empirical** — the MP2 correlation energy in
//!   [`post`] and the extended-Hückel method in [`semiempirical`].
//! - **Drivers** ([`driver`]) — top-level `run_rhf` / `run_uhf` /
//!   `run_mp2` / `run_dft` and a [`QchemReport`].
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, QchemError>`](error::QchemError). The error type carries
//! stable [`code`](error::QchemError::code) and
//! [`category`](error::QchemError::category) accessors for telemetry.
//!
//! ## v1 scope — honest caveats
//!
//! This is a real working program, not production parity with Psi4. The
//! Hartree-Fock core (geometry, basis, integrals, RHF, UHF, DIIS,
//! properties, MP2, extended Hückel) and the Kohn-Sham DFT subsystem
//! (the grid, the LDA / PBE / B3LYP functionals, the KS-SCF) are real
//! and validated. One genuine subsystem remains an honest typed stub
//! that returns [`QchemError::NotYetImplemented`]: **geometry
//! optimisation** (needs analytic energy gradients) — see
//! [`driver::GeometryOptRequest`].
//!
//! Other simplifications, documented inline in each module: the basis
//! library covers H–Ne only; angular momentum runs s, p, d (Cartesian
//! d functions, no spherical-harmonic transform); the ERI build is a
//! dense `O(n⁴)` loop with no Schwarz screening; MP2 is the canonical
//! closed-shell formula from the RHF reference only; DFT is
//! closed-shell (restricted Kohn-Sham) with three functionals and no
//! analytic gradients or dispersion correction.

#![forbid(unsafe_code)]

pub mod basis;
pub mod dft;
pub mod driver;
pub mod element;
pub mod error;
pub mod geometry;
pub mod integrals;
pub mod post;
pub mod properties;
pub mod scf;
pub mod semiempirical;

// --- Convenience re-exports of the most-used types --------------------

pub use basis::BasisSet;
pub use dft::{Functional, GridQuality, KsResult};
pub use driver::{run_dft, run_mp2, run_rhf, run_uhf, QchemReport};
pub use element::Element;
pub use error::{ErrorCategory, QchemError, Result};
pub use geometry::{Atom, MolecularGeometry};
pub use scf::{RhfResult, UhfResult};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crate_error_reexported() {
        let e = QchemError::not_yet("dft", "needs an XC library");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn geometry_reexport_round_trips() {
        let m = MolecularGeometry::new(vec![
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.0]).unwrap(),
            Atom::from_symbol_angstrom("H", [0.0, 0.0, 0.74]).unwrap(),
        ]);
        assert_eq!(m.n_atoms(), 2);
    }
}
