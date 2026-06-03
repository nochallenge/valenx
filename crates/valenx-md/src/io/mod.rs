//! Structure and trajectory I/O.
//!
//! **Roadmap features 3–6.** Readers and writers for the file formats
//! an MD workflow exchanges with the rest of the world:
//!
//! - [`pdb`] — the Protein Data Bank format: `ATOM` / `HETATM`
//!   records, the canonical structure interchange format (feature 3).
//! - [`xyz`] — the minimal XYZ format: an atom count, a comment, then
//!   `element x y z` lines (feature 4).
//! - [`gro`] — the GROMACS `.gro` format: fixed-column residue / atom
//!   records with positions, optional velocities and a box line
//!   (feature 5).
//! - [`trajectory`] — a [`trajectory::Trajectory`] container plus a
//!   compact binary DCD-class writer/reader and a human-readable
//!   framed-text format (feature 6).
//!
//! Coordinates in PDB files are Ångström (the format's unit) and are
//! converted to/from the crate's nm on read/write; XYZ and GRO are
//! handled in their native units (XYZ Ångström, GRO nm).

pub mod gro;
pub mod pdb;
pub mod trajectory;
pub mod xyz;

/// Ångström per nanometre — the PDB / XYZ length unit is Å, the
/// crate works in nm.
pub const ANGSTROM_PER_NM: f64 = 10.0;
