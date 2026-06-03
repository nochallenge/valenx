//! # valenx-fillet-brep
//!
//! True-BRep fillet and chamfer for Valenx (Phase 14).
//!
//! Where `valenx-fillet` (Phase 3) tessellates the input solid and
//! patches strip triangles onto sharp edges, this crate operates one
//! level up — on the BRep faces and edges themselves — and produces a
//! [`valenx_cad::Solid::Brep`] whose tangent surfaces round-trip
//! through STEP/IGES and survive further BRep ops (booleans, sweeps,
//! more fillets).
//!
//! # Honest scope
//!
//! A production-grade BRep fillet (rolling-ball trajectory + offset
//! surface intersection + N-way corner blends) is multi-month-to-year
//! research work that ships as 50–100 kLOC inside Parasolid, ACIS, and
//! OCCT. What ships here is the **genuinely-achievable slice**: the
//! single convex planar-edge fillet, built as real BRep constructive
//! solid geometry (see [`brep_build`]).
//!
//! - **Edge:** must be a *convex straight edge* shared by exactly two
//!   faces — fillets between curved faces, fillets at concave edges,
//!   and "face-fillet" (between non-adjacent faces) are explicit
//!   non-goals.
//! - **Faces:** both adjacent faces must be planar. The bisector
//!   normal between two planes is well-defined, the tangent points
//!   are exact, and the fillet surface is a single circular-arc
//!   sweep. Curved adjacent faces would require offset-surface
//!   intersection (deferred — Tier 3).
//! - **Radius:** constant *or* **linearly variable** along the edge.
//!   [`fillet::fillet_planar_edge`] is the constant-radius fillet;
//!   [`fillet::fillet_variable_radius_edge`] (Phase 14.6) tapers the
//!   radius linearly from a `radius_start` at one endpoint to a
//!   `radius_end` at the other, lofting the cutter / fillet-bar between
//!   the two end cross-sections. A *general* radius law (a spline of
//!   radius-vs-arc-length) would loft through intermediate stations
//!   and is a bounded follow-up.
//! - **Corner blending:** when 3+ filleted edges meet at a vertex (the
//!   usual case for "fillet every edge of a cube"), the
//!   **orthogonal convex 3-edge corner** — the corner of a box — is
//!   blended with a real **rolling-ball corner** (Phase 14.7): the
//!   radius-`r` ball seated tangent to all three faces, added back as
//!   a spherical octant via `(solid − corner_cutter) ∪ corner_ball`
//!   (see [`corner_build`]). [`fillet::fillet_solid_edges`] composes
//!   the per-edge fillets and then blends every box-style corner.
//!   General N-edge, non-orthogonal, and concave corners remain a
//!   Tier-3 problem (ChFi3d in OCCT is 50k+ LOC) — those corners fall
//!   back to the independent per-edge fillets.
//! - **Boolean robustness:** the fillet is `(solid − cutter) ∪ bar`
//!   via the real `truck_shapeops` booleans; the cutter trims flush
//!   with the adjacent faces, and coincident faces are the fragile
//!   input for any boolean kernel. On geometry the kernel cannot
//!   resolve, the fillet surfaces a soft error and the caller falls
//!   through to the mesh-domain pipeline.
//!
//! # Architecture
//!
//! - [`error::FilletBrepError`] — typed error enum mirroring
//!   `valenx-fillet`'s `FilletError` but tailored to BRep failure
//!   modes.
//! - [`topology`] — borrow-helpers for the underlying truck Solid (pull
//!   the inner type out of [`valenx_cad::Solid`], walk faces/edges,
//!   extract endpoints).
//! - [`edge_classify`] — adjacent-face lookup, planar-face detection,
//!   convex-edge test. Mirrors the mesh-domain classification in
//!   `valenx-fillet` but at the BRep level.
//! - [`fillet`] — the per-edge fillet itself ([`fillet::fillet_planar_edge`])
//!   and the batch helper ([`fillet::fillet_solid_edges`]).
//! - [`brep_build`] — the **real BRep fillet surgery** (Phase 14.5):
//!   builds the cutter + fillet-bar prisms and evaluates
//!   `(solid − cutter) ∪ bar` with the `truck_shapeops` booleans.
//! - [`chamfer`] — flat-bevel analog of [`fillet`].
//! - [`corner`] — corner-blend **detection**: classifies a vertex as
//!   the supported orthogonal-convex 3-edge corner and finds every
//!   blendable corner of a filleted-edge set.
//! - [`corner_build`] — the **real BRep corner-blend surgery** (Phase
//!   14.7): builds the corner cutter + seated ball and evaluates
//!   `(solid − corner_cutter) ∪ corner_ball` with the `truck_shapeops`
//!   booleans, producing the rolling-ball spherical corner.
//! - [`bridge`] — convenience selectors used by the
//!   `valenx-feature-tree` dispatcher: pick which edges to fillet
//!   given an angle threshold (parallels Phase 3's mesh-domain
//!   selector but on BRep faces).
//!
//! # Fallback policy
//!
//! Callers (the feature-tree `ops::fillet` / `ops::chamfer` evaluators)
//! should treat `FilletBrepError::NotPlanarFaces` and
//! `FilletBrepError::NonConvexEdge` as **soft errors** — the input is
//! geometrically valid, the BRep path just does not support it; fall
//! through to the Phase 3 mesh-domain pipeline.
//! `FilletBrepError::BadParameter` / `RadiusTooLarge` are hard errors
//! (bad input) and should be surfaced to the user verbatim.
//!
//! `FilletBrepError::TruckOp` raised by the fillet's boolean stage
//! means `truck_shapeops` could not resolve the flush cutter faces
//! for this particular geometry. The BRep fillet genuinely failed, so
//! a robustness-minded dispatcher may also fall through to the
//! mesh-domain pipeline on it — that produces the same fillet shape
//! by strip insertion. See [`error::FilletBrepError::is_soft`].

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod brep_build;
pub mod bridge;
pub mod chamfer;
pub mod corner;
pub mod corner_build;
pub mod edge_classify;
pub mod error;
pub mod fillet;
pub mod topology;

pub use error::{ErrorCategory, FilletBrepError};
