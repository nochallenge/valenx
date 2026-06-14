//! # valenx-surface
//!
//! NURBS curves + surfaces for Valenx (Phase 9 of the FreeCAD-parity
//! roadmap, extended in Phase 19).
//!
//! ## What's in here
//!
//! - [`NurbsCurve`] / [`NurbsSurface`] — the canonical data
//!   structures, with validated [`NurbsCurve::new`] /
//!   [`NurbsSurface::new`] constructors.
//! - Cox-de Boor basis-function evaluation
//!   ([`nurbs_curve::basis_functions`]) + knot-span lookup
//!   ([`nurbs_curve::find_knot_span`]) — building blocks every
//!   downstream module consumes.
//! - [`NurbsCurve::evaluate`] / [`NurbsSurface::evaluate`] —
//!   rational evaluation by tensor product.
//! - [`coons::fill`] — Coons patch from four boundary curves.
//! - [`sew::stitch`] (G0) / [`sew::g2_stitch`] (Phase 19C, G2
//!   continuity adjusting 3 rows of CPs per side) /
//!   [`extend::extrapolate`] — boundary manipulation.
//! - [`intersect::surface_surface`] (Phase 9 tessellation v1) /
//!   [`intersect::true_ssi`] (Phase 19B hybrid: tessellation seed
//!   then Newton midpoint-snap refinement + cubic NURBS LSQ fit) /
//!   [`march_ssi::marching_ssi`] (Phase 19F production: continuous
//!   Bajaj-style trace in parametric `(u, v)` of both surfaces with
//!   adaptive step + Newton closest-foot correction + cubic LSQ fit;
//!   handles smooth curved intersections, boundary termination, and
//!   loop closure) /
//!   [`trim::by_curve`] (legacy world-xy v1) /
//!   [`trim::by_curve_in_uv`] (Phase 9.5 parametric (u, v) domain
//!   trim — works on any surface, projects trim curve into (u, v)
//!   via Gauss-Newton foot-point) — boundary curve operations.
//! - [`knot_ops`] (Phase 19A) — Boehm knot insertion, Tiller knot
//!   removal, Bezier-decomposition degree elevation; curve + surface
//!   variants.
//! - [`fit`] (Phase 19D) — NURBS curve LSQ fit through points;
//!   tensor-product surface fit from a structured grid; scattered
//!   point cloud fit via plane projection.
//! - [`scatter_fit`] (Phase 19F) — production scattered point-cloud
//!   NURBS fitting with PCA principal-plane parameterisation,
//!   alternating parameter-vs-surface refinement (Newton closest-foot
//!   reprojection between iterations), and feature-edge detection
//!   that inserts a C0 knot line where the cloud bends sharply.
//! - [`blend`] (Phase 19F) — rolling-ball blend surface: traces the
//!   ball-center curve on the bisector of two surfaces, emits the
//!   contact curves on each surface + the blend strip swept between
//!   them as a tensor-product NURBS surface.
//! - [`ruled`] (Phase 19E) — ruled surface between two curves,
//!   linear extrusion along a vector, cone from a curve to an apex.
//! - [`tessellate::surface`] — sample a NURBS surface into a
//!   [`valenx_mesh::Mesh`] for the viewport.
//! - [`persist::SurfaceFile`] — RON envelope for round-tripping
//!   surfaces to disk.
//!
//! ## Output is mesh-backed
//!
//! For v1, every NURBS surface that leaves this crate to participate
//! in downstream Part Design / Mesh ops is tessellated and wrapped
//! via [`valenx_cad::Solid::from_mesh`]. True NURBS-to-BRep
//! conversion (so STEP export can serialise the rational basis) is
//! tracked under Phase 9.5.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod blend;
pub mod continuity;
pub mod coons;
pub mod curvature;
pub mod curvature_comb;
pub mod draft_angle;
pub mod error;
pub mod extend;
pub mod fit;
pub mod intersect;
pub mod knot_ops;
pub mod march_ssi;
pub mod nurbs_curve;
pub mod nurbs_surface;
pub mod persist;
pub mod ruled;
pub mod scatter_fit;
pub mod sew;
pub mod tessellate;
pub mod trim;

pub use continuity::{measure_edge_continuity, ContinuityReport};
pub use curvature::{local_shape, principal_curvatures, LocalShape};
pub use curvature_comb::{curvature_comb, CurvatureComb};
pub use draft_angle::{draft_angle, draft_report, DraftReport};
pub use error::SurfaceError;
pub use nurbs_curve::NurbsCurve;
pub use nurbs_surface::NurbsSurface;
pub use sew::Edge;
pub use trim::{TrimSide, UvTrimParams};
