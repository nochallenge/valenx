//! Real Rust-native CAD kernel for Valenx, built on the `truck` BRep
//! library. Replaces the "edit in FreeCAD" subprocess escape hatch
//! for the most common solid-modeling operations.
//!
//! Capabilities
//! ============
//!
//! - **Primitives** — [`box_solid`], [`cylinder`], [`sphere`],
//!   [`cone`], [`torus`], [`prism`]. All built from
//!   [`truck_modeling::builder`] sweeps so the resulting [`Solid`] is
//!   a proper closed BRep (faces / edges / vertices form a closed
//!   2-manifold), not just a triangle soup.
//! - **Boolean ops** — [`union`], [`difference`], [`intersection`]
//!   between two solids, routed through `truck-shapeops`'s `or` / `and`
//!   primitives. Difference is implemented as `A AND (NOT B)` after
//!   inverting B's face orientations.
//! - **Fillets** — [`fillet_edges`] is the documented stub: truck 0.6
//!   does NOT ship an edge-fillet algorithm, so the function returns a
//!   typed [`CadError::NotImplemented`] rather than silently leaking
//!   the unmodified solid. Callers must surface this back to the user
//!   so they can fall back to "Open in FreeCAD" for now.
//! - **Tessellation** — [`solid_to_mesh`] runs truck-meshalgo's
//!   constrained-Delaunay tessellator and converts the resulting
//!   `truck_polymesh::PolygonMesh` into a [`valenx_mesh::Mesh`] of
//!   `Tri3` elements suitable for the egui viewport renderer.
//! - **Measurement** — [`solid_volume`], [`solid_area`],
//!   [`is_closed_solid`], [`euler_characteristic`] compute mass
//!   properties + structural validity. They are how the validation
//!   suite proves a constructed solid against its analytic ground
//!   truth — see [`measure`].
//!
//! Tolerance notes
//! ===============
//!
//! Both tessellation and the boolean operators take an explicit
//! tolerance in model units. The toolbox UI hardcodes `0.5` for
//! tessellation and `0.05` for booleans for v1; we can expose those
//! through the form when users complain.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod boolean;
pub mod fillet;
pub mod measure;
pub mod primitives;
pub mod solid;
pub mod tessellate;

pub use boolean::{difference, intersection, union};
pub use fillet::fillet_edges;
pub use measure::{
    euler_characteristic, is_closed_solid, is_closed_solid_tol, solid_area, solid_area_tol,
    solid_bounding_box, solid_bounding_box_aspect_ratio, solid_bounding_box_aspect_ratio_tol,
    solid_bounding_box_center, solid_bounding_box_center_tol,
    solid_bounding_box_corey_shape_factor, solid_bounding_box_corey_shape_factor_tol,
    solid_bounding_box_diagonal, solid_bounding_box_diagonal_tol,
    solid_bounding_box_elongation, solid_bounding_box_elongation_tol,
    solid_bounding_box_extents, solid_bounding_box_extents_tol,
    solid_bounding_box_fill_ratio, solid_bounding_box_fill_ratio_tol,
    solid_bounding_box_flatness, solid_bounding_box_flatness_tol,
    solid_bounding_box_intercept_sphericity, solid_bounding_box_intercept_sphericity_tol,
    solid_bounding_box_intermediate_extent, solid_bounding_box_intermediate_extent_tol,
    solid_bounding_box_longest_extent, solid_bounding_box_longest_extent_tol,
    solid_bounding_box_max_projection_sphericity, solid_bounding_box_max_projection_sphericity_tol,
    solid_bounding_box_shortest_extent, solid_bounding_box_shortest_extent_tol,
    solid_bounding_box_surface_area, solid_bounding_box_surface_area_tol, solid_bounding_box_tol,
    solid_bounding_box_volume, solid_bounding_box_volume_tol,
    solid_centroid, solid_centroid_tol,
    solid_equivalent_cube_edge, solid_equivalent_cube_edge_tol,
    solid_equivalent_ellipsoid_semi_axes, solid_equivalent_ellipsoid_semi_axes_tol,
    solid_equivalent_sphere_area, solid_equivalent_sphere_area_tol,
    solid_equivalent_sphere_diameter, solid_equivalent_sphere_diameter_tol,
    solid_equivalent_sphere_radius, solid_equivalent_sphere_radius_tol,
    solid_genus, solid_inertia_anisotropy, solid_inertia_anisotropy_tol,
    solid_inertia_trace, solid_inertia_trace_tol,
    solid_mean_chord_length, solid_mean_chord_length_tol,
    solid_principal_moments, solid_principal_moments_tol,
    solid_principal_radii_of_gyration, solid_principal_radii_of_gyration_tol,
    solid_radius_of_gyration,
    solid_radius_of_gyration_tol, solid_sauter_mean_diameter, solid_sauter_mean_diameter_tol,
    solid_specific_surface_area, solid_specific_surface_area_tol,
    solid_sphericity, solid_sphericity_tol, solid_total_gaussian_curvature, solid_volume,
    solid_volume_tol,
};
pub use primitives::{box_solid, cone, cylinder, prism, sphere, torus};
pub use solid::{CadError, Solid};
pub use tessellate::solid_to_mesh;

/// Default tolerance used for tessellation when the caller doesn't
/// supply one. Matches the value the desktop toolbox passes.
pub const DEFAULT_TESS_TOLERANCE: f64 = 0.5;

/// Default tolerance used for boolean operations.
pub const DEFAULT_BOOL_TOLERANCE: f64 = 0.05;
