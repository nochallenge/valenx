//! # valenx-curvedbeam
//!
//! Winkler-Bach bending of sharply curved beams (crane hooks, chain links,
//! C-frames, press frames), with the straight-beam limit recovered as a
//! cross-check.
//!
//! ## What
//!
//! When a beam's radius of curvature is comparable to its depth, plane
//! sections still stay plane but the strain distribution becomes
//! hyperbolic, not linear: the neutral axis shifts *inward* from the
//! centroid toward the centre of curvature, and the inner-fibre stress is
//! substantially higher than ordinary `My/I` straight-beam theory predicts.
//! This crate computes, in closed form, the neutral-axis shift and the
//! inner/outer fibre stresses for the two classic textbook cross-sections
//! (rectangle and solid circle), and exposes the straight-beam stress for
//! the same section so the two can be compared directly.
//!
//! Modules:
//!
//! [`section`] builds a cross-[`Section`] and its derived
//! [`SectionProps`] (area, centroidal radius `R`, neutral-axis radius
//! `R_n`, eccentricity `e`).
//!
//! [`stress`] evaluates [`stress_at_radius`], [`stress_inner`],
//! [`stress_outer`], the [`fibre_stresses`] bundle, and the limiting
//! [`straight_beam_stress`].
//!
//! ## Model
//!
//! The neutral axis sits at radius
//!
//! ```text
//! R_n = A / integral(dA / r)
//! ```
//!
//! (`A` the area, `r` the radial distance of an area element from the
//! centre of curvature), giving the closed forms
//!
//! ```text
//! rectangle:      R_n = (r_o - r_i) / ln(r_o / r_i)
//! solid circle:   R_n = (R + sqrt(R^2 - c^2)) / 2
//! ```
//!
//! The eccentricity is `e = R - R_n >= 0`, and the bending stress at a fibre
//! of radius `r`, with `y = R_n - r` the distance from the neutral axis
//! (positive toward the centre of curvature), is the Winkler-Bach formula
//!
//! ```text
//! sigma = M * y / (A * e * (R_n - y)) = M * y / (A * e * r).
//! ```
//!
//! `M` is positive when it tends to decrease the curvature (straighten the
//! beam); a positive `M` then puts the inner fibre in tension and the outer
//! fibre in compression. As `R / depth -> infinity`, `e -> I / (A R)` and
//! the formula collapses to the straight-beam result `sigma = M y / I`.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are the standard textbook closed-form
//! Winkler-Bach equations (Boresi; Shigley), validated in the test-suite
//! against analytic ground truth: the hand-worked neutral-axis value, the
//! sign of the inner-fibre stress, the zero-stress neutral axis, and the
//! large-radius limit converging to `My/I`. It is NOT a clinical/medical or
//! production-certified engineering tool: it models a single in-plane
//! bending moment on a prismatic, linear-elastic, homogeneous section and
//! accounts for no component tolerances, safety factors, fatigue or thermal
//! limits, stress-concentration details at holes/keyways, transverse-shear
//! or direct-axial superposition, or any code-compliance check. Do not size
//! load-bearing hardware from its output.
//!
//! ## Example
//!
//! ```rust
//! use valenx_curvedbeam::{Section, fibre_stresses, straight_beam_inner};
//!
//! // Rectangular curved beam: width 1, inner radius 4, outer radius 6,
//! // under a straightening moment M = 100 (consistent units).
//! let section = Section::rectangle(1.0, 4.0, 6.0).unwrap();
//! let f = fibre_stresses(&section, 100.0).unwrap();
//!
//! // Neutral axis lies inside the centroid (R = 5).
//! assert!(f.props.r_neutral < f.props.r_centroid);
//!
//! // Inner fibre is in tension, outer fibre in compression.
//! assert!(f.inner > 0.0);
//! assert!(f.outer < 0.0);
//!
//! // The inner-fibre stress beats the straight-beam My/I estimate — the
//! // hallmark of a sharply curved beam.
//! let straight = straight_beam_inner(&section, 100.0).unwrap();
//! assert!(f.inner > straight);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod section;
pub mod stress;

pub use error::{CurvedBeamError, ErrorCategory};
pub use section::{Section, SectionProps};
pub use stress::{
    fibre_stresses, straight_beam_inner, straight_beam_stress, stress_at_radius, stress_inner,
    stress_outer, FibreStresses,
};
