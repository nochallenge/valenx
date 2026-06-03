//! Molecular properties derived from a converged wavefunction.
//!
//! Once the SCF has produced a density matrix and a set of orbitals,
//! this module extracts the chemically interesting numbers:
//!
//! - [`population`] — Mulliken and Löwdin population analysis and the
//!   resulting atomic partial charges.
//! - [`dipole`] — the electronic-plus-nuclear dipole moment.
//! - [`orbitals`] — molecular-orbital energies, occupations and the
//!   HOMO-LUMO gap.
//! - [`grid`] — orbital and density evaluation on a 3-D grid for
//!   visualisation (Gaussian-cube-style volumetric data).
//! - [`bondorder`] — Mayer / Wiberg bond orders and atomic valences.

pub mod bondorder;
pub mod dipole;
pub mod grid;
pub mod orbitals;
pub mod population;

pub use bondorder::{mayer_bond_orders, BondOrderMatrix};
pub use dipole::{dipole_moment, DipoleMoment};
pub use grid::{density_grid, orbital_grid, VolumetricGrid};
pub use orbitals::{restricted_summary, unrestricted_spin_summary, OrbitalSummary};
pub use population::{lowdin, mulliken, PopulationAnalysis};
