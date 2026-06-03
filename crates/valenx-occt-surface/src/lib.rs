//! # valenx-occt-surface
//!
//! Phases 70-100 — OpenCASCADE (OCCT) **surface
//! modeling** feature parity for Valenx.
//!
//! OCCT ships ~1M lines of C++ across 200+ packages and roughly 30
//! years of accumulated geometry-kernel work. This crate provides a
//! Rust-native API surface for 31 of the most commonly invoked OCCT
//! surface-modeling APIs, organised as one module per feature so each
//! can advance from "scaffold" → "v1 implementation" → "production
//! parity" independently.
//!
//! ## v1 strategy
//!
//! Each of the 31 functions either:
//!
//! 1. Maps cleanly onto an existing `valenx-cad` / `valenx-surface` /
//!    `truck-modeling` capability, in which case it implements
//!    honestly (parameter validation + delegation + typed errors); or
//! 2. Returns [`OcctSurfaceError::NotYetImplemented`] with rustdoc
//!    describing what the real OCCT API does and which Phase `N.5`
//!    follow-up will deliver the deep implementation.
//!
//! Both kinds carry a public function signature so downstream crates
//! (the toolbox UI, the feature tree, integration tests) can be
//! written against the final API today and have the stubs fill in
//! later without churn.
//!
//! ## Feature catalogue (Phases 70-100)
//!
//! - **Sectioning / intersection:** [`algo_section()`], [`cut_api_section()`].
//! - **Sweeps:** [`pipe_shell()`], [`sweep_api_pipe()`],
//!   [`sweep_api_pipe_shell()`], [`sweep_api_thru_sections()`],
//!   [`sweep_api_evolved()`].
//! - **Offsets / shells:** [`offset_surface()`],
//!   [`offset_api_make_offset()`], [`offset_api_draft_angle()`],
//!   [`offset_api_filling()`].
//! - **Sewing / stitching:** [`builder_sewing()`].
//! - **Blends / fillets:** [`chfi3d_corner_blends()`].
//! - **Surface filling / fitting:** [`geom_fill_bsurf_filling()`],
//!   [`geom_fill_section_law()`], [`approx_surface_fit()`],
//!   [`approx_curve_fit()`].
//! - **Primitives:** [`prim_api_box()`], [`prim_api_cylinder()`],
//!   [`prim_api_cone()`], [`prim_api_sphere()`], [`prim_api_torus()`],
//!   [`prim_api_prism()`], [`prim_api_revol()`], [`prim_api_wedge()`],
//!   [`prim_api_half_space()`].
//! - **Boolean cuts:** [`cut_api_cut_with_warp()`].
//! - **Feature-based modeling:** [`feat_make_prism()`],
//!   [`feat_make_revol()`], [`feat_make_pipe()`], [`feat_make_draft()`].
//!
//! ## Error model
//!
//! All public APIs return [`Result<_, OcctSurfaceError>`]. See
//! [`OcctSurfaceError::code`] / [`OcctSurfaceError::category`] for
//! the stable taxonomy.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;

// Shared machinery for the feature-based ops (Phases 97-100).
pub mod feat_support;

// Shared frame-transport machinery for the sweep family (Phases 71, 89, 91).
pub mod sweep_support;

// Phase 70-79
pub mod algo_section;
pub mod pipe_shell;
pub mod offset_surface;
pub mod builder_sewing;
pub mod chfi3d_corner_blends;
pub mod geom_fill_bsurf_filling;
pub mod geom_fill_section_law;
pub mod approx_surface_fit;
pub mod approx_curve_fit;
pub mod prim_api_cylinder;

// Phase 80-89
pub mod prim_api_cone;
pub mod prim_api_sphere;
pub mod prim_api_torus;
pub mod prim_api_box;
pub mod prim_api_prism;
pub mod prim_api_revol;
pub mod prim_api_wedge;
pub mod prim_api_half_space;
pub mod sweep_api_pipe;
pub mod sweep_api_pipe_shell;

// Phase 90-100
pub mod sweep_api_thru_sections;
pub mod sweep_api_evolved;
pub mod offset_api_make_offset;
pub mod offset_api_draft_angle;
pub mod offset_api_filling;
pub mod cut_api_section;
pub mod cut_api_cut_with_warp;
pub mod feat_make_prism;
pub mod feat_make_revol;
pub mod feat_make_pipe;
pub mod feat_make_draft;

pub use error::{ErrorCategory, OcctSurfaceError};

// Re-export the entry points so callers can `use valenx_occt_surface::pipe_shell;`
// instead of the full module-path mouthful.
pub use algo_section::algo_section;
pub use approx_curve_fit::approx_curve_fit;
pub use approx_surface_fit::approx_surface_fit;
pub use builder_sewing::builder_sewing;
pub use chfi3d_corner_blends::chfi3d_corner_blends;
pub use cut_api_cut_with_warp::cut_api_cut_with_warp;
pub use cut_api_section::cut_api_section;
pub use feat_make_draft::feat_make_draft;
pub use feat_make_pipe::feat_make_pipe;
pub use feat_make_prism::feat_make_prism;
pub use feat_make_revol::feat_make_revol;
pub use geom_fill_bsurf_filling::geom_fill_bsurf_filling;
pub use geom_fill_section_law::geom_fill_section_law;
pub use offset_api_draft_angle::offset_api_draft_angle;
pub use offset_api_filling::offset_api_filling;
pub use offset_api_make_offset::offset_api_make_offset;
pub use offset_surface::offset_surface;
pub use pipe_shell::pipe_shell;
pub use prim_api_box::prim_api_box;
pub use prim_api_cone::prim_api_cone;
pub use prim_api_cylinder::{prim_api_cylinder, prim_api_cylinder_on_axis};
pub use prim_api_half_space::{prim_api_half_space, HALF_SPACE_EXTENT};
pub use prim_api_prism::prim_api_prism;
pub use prim_api_revol::prim_api_revol;
pub use prim_api_sphere::prim_api_sphere;
pub use prim_api_torus::prim_api_torus;
pub use prim_api_wedge::prim_api_wedge;
pub use sweep_api_evolved::sweep_api_evolved;
pub use sweep_api_pipe::sweep_api_pipe;
pub use sweep_api_pipe_shell::sweep_api_pipe_shell;
pub use sweep_api_thru_sections::sweep_api_thru_sections;
