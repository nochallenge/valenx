//! Structural geometry: distances, angles, Ramachandran analysis,
//! contact maps, clash detection, bulk shape descriptors and SASA.
//!
//! The submodules are independent — pull in only what a given
//! analysis needs:
//!
//! - [`distance`](mod@distance) — pairwise distance + a
//!   grid-accelerated neighbour search ([`distance::NeighborGrid`]).
//! - [`angles`] — bond angles, dihedrals, backbone φ/ψ/ω and
//!   sidechain χ torsions.
//! - [`ramachandran`] — φ/ψ region classification.
//! - [`contacts`] — residue contact maps and steric-clash detection.
//! - [`shape`] — centre of mass, radius of gyration, principal axes.
//! - [`sasa`] — Shrake-Rupley solvent-accessible surface area.
//! - [`vdw`] — van der Waals / covalent radii tables.

pub mod angles;
pub mod contacts;
pub mod distance;
pub mod ramachandran;
pub mod sasa;
pub mod shape;
pub mod vdw;

pub use angles::{backbone_torsions, bond_angle, dihedral, sidechain_chi, BackboneTorsions};
pub use contacts::{contact_map, detect_clashes, Clash, ContactMap, ContactMode};
pub use distance::{distance, NeighborGrid};
pub use ramachandran::{classify, RamachandranRegion};
pub use sasa::{shrake_rupley, SasaResult, WATER_PROBE};
pub use shape::{center_of_mass, principal_axes, radius_of_gyration, PrincipalAxes};
