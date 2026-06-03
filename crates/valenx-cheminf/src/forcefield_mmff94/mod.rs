//! Production-depth MMFF94 (Merck Molecular Force Field) implementation.
//!
//! This module replaces the v1 reduced [`crate::forcefield`] for the
//! conformer-cleanup hot path. It implements the published MMFF94
//! parameterisation (Halgren, J. Comp. Chem. 17, 490-641, 1996) — an
//! atom-typed force field with bond, angle, stretch-bend, torsion,
//! buffered-14-7 van-der-Waals and Coulomb-electrostatic energy terms
//! and an analytic gradient suitable for fast minimisation.
//!
//! ## Submodules
//!
//! - [`atom_type`] — perceive each atom's MMFF94 type from the
//!   molecular graph (Halgren table I).
//! - [`params`] — the tabulated bond / angle / stretch-bend / torsion
//!   / vdW parameters (transcribed from the published tables for the
//!   common organic atom-type subset).
//! - [`mod@energy`] — the energy expression, gradient and minimiser. The
//!   [`energy::Mmff94Setup`] caches all per-bond / per-angle
//!   parameters and partial charges so the inner loop of minimisation
//!   never re-types or re-charges.
//!
//! ## Coverage and v1 honesty
//!
//! - **Atom types.** A representative subset of MMFF94's ~95 published
//!   atom types covering common organic chemistry (C / H / N / O / S
//!   / P / halogens). The full subset is listed in
//!   [`atom_type`]. Atoms outside the subset get
//!   [`MmffType::UNKNOWN`](atom_type::MmffType::UNKNOWN) and fall
//!   back to rule-based parameters.
//! - **Parameter tables** are transcribed from Halgren 1996 parts
//!   II-V for the typed atoms. Missing combinations fall back to the
//!   covalent-radius / hybridisation-rule fallbacks MMFF94 documents.
//! - **Charges** use the Gasteiger-PEOE model in [`crate::charge`] as
//!   a substitute for MMFF94's full bond-charge-increment table —
//!   that's the largest single data file the parameterisation needs
//!   and the main remaining gap vs full MMFF94.
//! - **Out-of-plane bending** is not yet included (a small term for
//!   the typed subset; documented gap).
//!
//! ## Public API
//!
//! Most callers want [`clean_up_geometry`] — set up parameters,
//! minimise to a local energy minimum in place. The energy /
//! gradient / minimiser primitives are re-exported for users who want
//! finer-grained control (e.g. multi-conformer search).

pub mod atom_type;
pub mod energy;
pub mod params;

pub use atom_type::{type_molecule, MmffType};
pub use energy::{energy, gradient, minimize, setup, Mmff94Energy, Mmff94Setup};

use crate::molecule::Molecule;

/// MMFF94 cleanup: type the atoms, compute partial charges, build all
/// the parameter lists once, and minimise the energy. Mutates the
/// coordinates in place; returns the final [`Mmff94Energy`].
pub fn clean_up_geometry(mol: &mut Molecule, max_steps: usize) -> Mmff94Energy {
    if mol.coords.len() != mol.atoms.len() || mol.atoms.len() < 2 {
        let s = setup(mol);
        return energy(mol, &s);
    }
    let s = setup(mol);
    minimize(mol, &s, max_steps)
}
