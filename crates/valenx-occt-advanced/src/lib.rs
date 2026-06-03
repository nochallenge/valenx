//! # valenx-occt-advanced
//!
//! Round 5 Block 3 (Phases 131-160) — OpenCASCADE (OCCT) **advanced
//! operations** feature parity for Valenx.
//!
//! OCCT's advanced-operations layer is the second-largest subsystem
//! after the surface-modeling kernel. It includes:
//!
//! - **`BRepOffsetAPI_*`** sketch-driven sweep/loft/draft builders
//!   that go beyond the basic ThruSections / DraftAngle in Block 1.
//! - **`BRepFeat_*`** feature-based variants (Pad / Pocket /
//!   sketch-driven extrude) that integrate features into an existing
//!   solid rather than creating standalone bodies.
//! - **`ShapeAnalysis_*`** read-only inspectors that report
//!   topological / geometric defects (open shells, reversed faces,
//!   degenerate edges, non-monotonic parameterisations, …).
//! - **`ShapeUpgrade_*`** writers that consume an analyzer's findings
//!   and emit a repaired shape (unify same-domain faces, close open
//!   wires, split at C1 discontinuities, …).
//! - **`GeomLib_*`** geometric-query primitives that lift queries off
//!   raw NURBS evaluation (tangent / normal / curvature / evolute,
//!   point-on-curve intersection).
//! - **`BRepFeat_*` local operations** that surgically replace,
//!   remove, or split faces / solids without re-running the full
//!   feature tree.
//!
//! This crate provides a Rust-native API surface for 30 of those, one
//! module per OCCT-equivalent feature so each can advance from
//! "scaffold" → "v1 implementation" → "production parity"
//! independently.
//!
//! ## v1 strategy
//!
//! Each of the 30 functions either:
//!
//! 1. Maps cleanly onto an existing `valenx-cad` / `valenx-surface` /
//!    `truck-modeling` capability, in which case it implements
//!    honestly (parameter validation + delegation + typed errors); or
//! 2. Returns [`OcctAdvancedError::NotYetImplemented`] with rustdoc
//!    describing what the real OCCT API does and which Phase `N.5`
//!    follow-up will deliver the deep implementation.
//!
//! Both kinds carry a public function signature so downstream crates
//! (the toolbox UI, the feature tree, integration tests) can be
//! written against the final API today and have the stubs fill in
//! later without churn.
//!
//! ## Feature catalogue (Phases 131-160)
//!
//! ### Sweep + Loft advanced (Phases 131-138)
//! - [`offset_api_thru_sections_with_guides()`],
//!   [`offset_api_draft_angle_with_neutral_plane()`],
//!   [`feat_make_prism_with_sketch()`],
//!   [`feat_make_revol_with_sketch()`],
//!   [`feat_make_pipe_with_path_constraint()`],
//!   [`feat_make_dgreater_pad()`],
//!   [`feat_make_dsubtract_pocket()`],
//!   [`feat_make_loft_with_rails()`].
//!
//! ### Shape analysis (Phases 139-146)
//! - [`shape_analysis_freebounds()`],
//!   [`shape_analysis_orientedclosedsolid()`],
//!   [`shape_analysis_curve_validity()`],
//!   [`shape_analysis_surface_validity()`],
//!   [`shape_analysis_wireorder()`],
//!   [`shape_analysis_check_dist_ratio()`],
//!   [`shape_analysis_topology()`],
//!   [`shape_analysis_fix_shape()`].
//!
//! ### Shape upgrade (Phases 147-152)
//! - [`shape_upgrade_unifysamedomain()`],
//!   [`shape_upgrade_shapeconvert_revolution_to_bspline()`],
//!   [`shape_upgrade_remove_internal_wires()`],
//!   [`shape_upgrade_split_continuity()`],
//!   [`shape_upgrade_close_open_wires()`],
//!   [`shape_upgrade_face_division()`].
//!
//! ### Geometric library (Phases 153-157)
//! - [`geom_lib_intersect_point_curve()`],
//!   [`geom_lib_normal_at_point()`],
//!   [`geom_lib_tangent_at_point()`],
//!   [`geom_lib_curvature_at_point()`],
//!   [`geom_lib_evolute_curve()`].
//!
//! ### Local operations (Phases 158-160)
//! - [`local_op_replace_face()`],
//!   [`local_op_remove_face()`],
//!   [`local_op_split_solid_with_plane()`].
//!
//! ## Error model
//!
//! All public APIs return [`Result<_, OcctAdvancedError>`]. See
//! [`OcctAdvancedError::code`] / [`OcctAdvancedError::category`] for
//! the stable taxonomy.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;

// Shared guide-curve loft machinery for Phases 131 + 138.
pub mod guide_loft;

// Sweep + Loft advanced (Phases 131-138)
pub mod offset_api_thru_sections_with_guides;
pub mod offset_api_draft_angle_with_neutral_plane;
pub mod feat_make_prism_with_sketch;
pub mod feat_make_revol_with_sketch;
pub mod feat_make_pipe_with_path_constraint;
pub mod feat_make_dgreater_pad;
pub mod feat_make_dsubtract_pocket;
pub mod feat_make_loft_with_rails;

// Shape analysis (Phases 139-146)
pub mod shape_analysis_freebounds;
pub mod shape_analysis_orientedclosedsolid;
pub mod shape_analysis_curve_validity;
pub mod shape_analysis_surface_validity;
pub mod shape_analysis_wireorder;
pub mod shape_analysis_check_dist_ratio;
pub mod shape_analysis_topology;
pub mod shape_analysis_fix_shape;

// Shape upgrade (Phases 147-152)
pub mod shape_upgrade_unifysamedomain;
pub mod shape_upgrade_shapeconvert_revolution_to_bspline;
pub mod shape_upgrade_remove_internal_wires;
pub mod shape_upgrade_split_continuity;
pub mod shape_upgrade_close_open_wires;
pub mod shape_upgrade_face_division;

// Geometric library (Phases 153-157)
pub mod geom_lib_intersect_point_curve;
pub mod geom_lib_normal_at_point;
pub mod geom_lib_tangent_at_point;
pub mod geom_lib_curvature_at_point;
pub mod geom_lib_evolute_curve;

// Local operations (Phases 158-160)
pub mod local_op_replace_face;
pub mod local_op_remove_face;
pub mod local_op_split_solid_with_plane;

pub use error::{ErrorCategory, OcctAdvancedError};

// Re-export the entry points so callers can
// `use valenx_occt_advanced::shape_analysis_freebounds;`
// instead of the full module-path mouthful.
pub use feat_make_dgreater_pad::feat_make_dgreater_pad;
pub use feat_make_dsubtract_pocket::feat_make_dsubtract_pocket;
pub use feat_make_loft_with_rails::feat_make_loft_with_rails;
pub use feat_make_pipe_with_path_constraint::{
    feat_make_pipe_with_path_constraint, FrameLaw,
};
pub use feat_make_prism_with_sketch::feat_make_prism_with_sketch;
pub use feat_make_revol_with_sketch::feat_make_revol_with_sketch;
pub use geom_lib_curvature_at_point::geom_lib_curvature_at_point;
pub use geom_lib_evolute_curve::geom_lib_evolute_curve;
pub use geom_lib_intersect_point_curve::geom_lib_intersect_point_curve;
pub use geom_lib_normal_at_point::geom_lib_normal_at_point;
pub use geom_lib_tangent_at_point::geom_lib_tangent_at_point;
pub use local_op_remove_face::local_op_remove_face;
pub use local_op_replace_face::local_op_replace_face;
pub use local_op_split_solid_with_plane::local_op_split_solid_with_plane;
pub use offset_api_draft_angle_with_neutral_plane::offset_api_draft_angle_with_neutral_plane;
pub use offset_api_thru_sections_with_guides::{
    offset_api_thru_sections_with_guides, offset_api_thru_sections_with_guides_ex,
};
pub use shape_analysis_check_dist_ratio::shape_analysis_check_dist_ratio;
pub use shape_analysis_curve_validity::shape_analysis_curve_validity;
pub use shape_analysis_fix_shape::shape_analysis_fix_shape;
pub use shape_analysis_freebounds::shape_analysis_freebounds;
pub use shape_analysis_orientedclosedsolid::shape_analysis_orientedclosedsolid;
pub use shape_analysis_surface_validity::shape_analysis_surface_validity;
pub use shape_analysis_topology::shape_analysis_topology;
pub use shape_analysis_wireorder::shape_analysis_wireorder;
pub use shape_upgrade_close_open_wires::{
    shape_upgrade_close_open_wires, shape_upgrade_close_open_wires_arc,
};
pub use shape_upgrade_face_division::shape_upgrade_face_division;
pub use shape_upgrade_remove_internal_wires::shape_upgrade_remove_internal_wires;
pub use shape_upgrade_shapeconvert_revolution_to_bspline::shape_upgrade_shapeconvert_revolution_to_bspline;
pub use shape_upgrade_split_continuity::shape_upgrade_split_continuity;
pub use shape_upgrade_unifysamedomain::shape_upgrade_unifysamedomain;
