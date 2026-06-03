//! Trajectory and ensemble analysis.
//!
//! **Roadmap features 26–29.** The observables an MD run is actually
//! run *for*:
//!
//! - [`reporters`] — instantaneous thermodynamic observables: kinetic,
//!   potential and total energy; the equipartition temperature; the
//!   virial pressure; and an [`reporters::ObservableLog`] that records
//!   them over a run (feature 26).
//! - [`rdf`] — the radial distribution function `g(r)`, the
//!   pair-correlation structure of a fluid (feature 27).
//! - [`msd`] — the mean-squared displacement and the Einstein
//!   diffusion coefficient `D = MSD/(6t)` (feature 28).
//! - [`rmsd`] — RMSD and RMSF after an optimal Kabsch superposition,
//!   the standard structural-deviation measures (feature 29).
//!
//! These work on a stored [`crate::io::trajectory::Trajectory`] or
//! directly on [`crate::system::System`] snapshots.

pub mod msd;
pub mod rdf;
pub mod reporters;
pub mod rmsd;
