//! Format-specific readers + writers for canonical [`crate::Mesh`].
//!
//! Each submodule covers one file format. Pull only the modules your
//! adapter needs to keep build times tight.
//!
//! - [`obj`] — Wavefront OBJ. Minimal subset (`v`, `f`), good enough
//!   for downstream tooling that just needs triangles + positions.
//! - [`ply`] — Stanford PLY (ASCII + binary little/big-endian).
//!   Header-driven, parses any `element vertex N` / `element face M`
//!   declaration.
//! - [`three_mf`] — 3MF stub. Returns an error documenting the
//!   missing `zip` dependency. v1.5 lands the implementation once
//!   `zip` lands in workspace deps.

pub mod obj;
pub mod ply;
pub mod three_mf;
