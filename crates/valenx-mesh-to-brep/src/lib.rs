//! # valenx-mesh-to-brep
//!
//! Mesh → BRep reverse engineering for Valenx (Phase 23).
//!
//! Workflow:
//!
//! 1. **Region detection.** [`feature_detect::detect_planes`] runs a
//!    region-growing fit on the mesh's triangle plane equations.
//!    [`feature_detect::detect_cylinders`] and
//!    [`feature_detect::detect_spheres`] are **real RANSAC / algebraic
//!    fits** (Phase 23 v2): the cylinder detector seeds axes from
//!    cross-products of non-parallel normals and confirms with an
//!    algebraic circle fit; the sphere detector runs an algebraic
//!    least-squares sphere fit. Both reject a candidate whose fit
//!    residual is too large, so a cube never mis-classifies as a
//!    curved primitive. [`feature_detect::detect_cones`] is deferred.
//! 2. **Surface fitting.** Each region becomes a NURBS surface.
//!    [`reconstruct::nurbs_from_region`] emits a flat quad patch for a
//!    quick visualisation; the **production path** —
//!    [`reconstruct::fit_planar_region_nurbs`] and
//!    [`reconstruct::fit_cylindrical_region_nurbs`] — fits a genuine
//!    tensor-product NURBS surface through the region's vertices via
//!    `valenx-surface` and returns a [`reconstruct::RegionFit`]
//!    **tolerance report** (RMS + worst-case deviation + sample
//!    count).
//! 3. **BRep reconstruction.** [`reconstruct::brep_from_mesh`] glues
//!    every fitted surface into one [`valenx_cad::Solid::Mesh`]
//!    container.
//! 4. **Closed-BRep sewing.** [`sew::sew_regions`] stitches the fitted
//!    regions into a watertight solid — it recognises a fitted **box**
//!    (six planar regions → a real closed `Solid::Brep`) and a fitted
//!    **cylinder** (lateral region + caps → a real closed
//!    `Solid::Brep`), and welds the patches into a watertight
//!    mesh-backed shell otherwise.
//!
//! ## Honest scope
//!
//! - Cone detection is deferred (a robust cone fit is genuinely harder
//!   than the plane / cylinder / sphere cases); un-classified
//!   triangles route through the plane fallback.
//! - Closed-BRep sewing ([`sew::sew_regions`]) reconstructs a genuine
//!   `Solid::Brep` for the recognised box / cylinder cases; every
//!   other shape sews in the mesh domain (a watertight mesh-backed
//!   `Solid` when the welded patches close, an open patch set
//!   otherwise — reported via [`sew::SewOutcome`]). A general
//!   parametric trim-and-stitch of arbitrary fitted NURBS faces into a
//!   `Solid::Brep` shell stays a Tier-3 follow-up gated on the
//!   parametric-BRep substrate.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod feature_detect;
pub mod reconstruct;
pub mod sew;

pub use feature_detect::{
    detect_cones, detect_cylinders, detect_planes, detect_spheres, ConicalRegion,
    CylindricalRegion, PlanarRegion, SphericalRegion,
};
pub use reconstruct::{
    brep_from_mesh, fit_cylindrical_region_nurbs, fit_planar_region_nurbs, nurbs_from_region,
    ReconstructError, RegionFit,
};
pub use sew::{sew_regions, SewOutcome, SewReport, SewResult};
