//! The self-consistent-field engine.
//!
//! This module turns an [`IntegralSet`](crate::integrals::IntegralSet)
//! into a converged Hartree-Fock wavefunction.
//!
//! - [`linalg`] — Löwdin and canonical orthogonalisation and the
//!   generalized symmetric eigenproblem `FC = SCε`.
//! - [`diis`] — Pulay DIIS convergence acceleration.
//! - [`rhf`] — the restricted (closed-shell) Hartree-Fock loop.
//! - [`uhf`] — the unrestricted (open-shell) Hartree-Fock loop.
//! - [`density`] — density-matrix sanity checks (electron count,
//!   idempotency).
//!
//! The two entry points are
//! [`run_rhf_scf`](rhf::run_rhf_scf) and
//! [`run_uhf_scf`](uhf::run_uhf_scf); the [`driver`](crate::driver)
//! module wraps them with geometry / basis setup and reporting.

pub mod density;
pub mod diis;
pub mod linalg;
pub mod rhf;
pub mod uhf;

pub use density::{DensityCheck, DensityCheckReport};
pub use rhf::{RhfResult, ScfIteration, ScfSettings};
pub use uhf::UhfResult;
