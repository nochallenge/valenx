//! Post-Hartree-Fock correlation methods.
//!
//! Hartree-Fock is a mean-field theory — each electron feels only the
//! average field of the others. *Post-HF* methods add the missing
//! electron correlation on top of a converged HF reference.
//!
//! - [`mp2`] — second-order Møller-Plesset perturbation theory, the
//!   cheapest systematic correlation correction, built from the RHF
//!   reference via an AO→MO integral transformation.

pub mod mp2;

pub use mp2::{ao_to_mo_transform, mp2_energy, Mp2Result};
