//! # valenx-relativity
//!
//! An in-house general-relativity engine, built from scratch in Rust — not a
//! wrapper around SageMath, Mathematica or the Einstein Toolkit.
//!
//! The core is a **metric-agnostic tensor-calculus engine**: give it any
//! metric `g_μν(x)` (as a function generic over the scalar type) and it returns
//! the Christoffel symbols, the Riemann/Ricci/Einstein tensors, and the Ricci
//! and Kretschmann curvature scalars — using forward-mode automatic
//! differentiation, so derivatives are exact to machine precision with no
//! finite-difference step to tune. See [`curvature_at`].
//!
//! It ships the full **black-hole family** as built-in, validated spacetimes —
//! [`schwarzschild`], [`kerr`], [`reissner_nordstrom`], and the general
//! [`KerrNewman`] master metric they specialise — plus flat [`Minkowski`] as a
//! zero-curvature baseline.
//!
//! ## Units
//!
//! Geometrized units `G = c = 1` throughout; mass, spin and charge all have
//! dimensions of length, and the natural scale is `M = 1` (so the Schwarzschild
//! horizon sits at `r = 2`). An SI conversion layer is provided separately.
//!
//! ## Honesty
//!
//! Every fallible path returns [`Result`]; the engine never emits a silent
//! `NaN`/`∞`. Each quantity is cross-checked against closed-form ground truth
//! in the test suite (Schwarzschild is a vacuum solution, so `R_{ab} = 0` and
//! `K = 48 M²/r⁶`; the same vacuum check holds for Kerr).
//!
//! ## Example
//!
//! ```
//! use valenx_relativity::{schwarzschild, curvature_at};
//!
//! // A unit-mass Schwarzschild black hole, evaluated at r = 5M on the equator.
//! let bh = schwarzschild(1.0);
//! let c = curvature_at(&bh, [0.0, 5.0, std::f64::consts::FRAC_PI_2, 0.0]).unwrap();
//!
//! // Vacuum solution: Ricci tensor vanishes.
//! assert!(c.ricci_scalar.abs() < 1e-9);
//! // Kretschmann scalar K = 48 M² / r⁶.
//! let expected = 48.0 / 5.0_f64.powi(6);
//! assert!((c.kretschmann - expected).abs() / expected < 1e-9);
//! ```

#![forbid(unsafe_code)]

mod autodiff;
pub mod curvature;
pub mod metric;
pub mod observables;
pub mod spacetimes;
pub mod tensor;
pub mod thermo;
pub mod units;

pub use autodiff::{Dual, HyperDual, Scalar};
pub use curvature::{curvature_at, Curvature};
pub use metric::{CoordSystem, Spacetime};
pub use observables::{
    ergosphere_radius, gravitational_redshift, horizons, isco, photon_sphere, shadow_radius,
    Horizons, OrbitSense,
};
pub use spacetimes::{kerr, reissner_nordstrom, schwarzschild, KerrNewman, Minkowski};
pub use thermo::{thermodynamics, Thermodynamics};

use thiserror::Error;

/// Errors returned by the relativity engine.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum RelativityError {
    /// A physical parameter was invalid — e.g. a non-positive mass, or a
    /// super-extremal spin/charge (`a² + Q² > M²`) that would describe a naked
    /// singularity rather than a black hole.
    #[error("invalid parameter: {0}")]
    InvalidParameter(String),

    /// The metric is singular or non-finite at the requested point
    /// (`det g → 0`), as happens on a horizon or coordinate axis where
    /// Boyer–Lindquist coordinates break down.
    #[error("coordinate singularity: {0}")]
    CoordinateSingularity(String),

    /// A geodesic integration failed to converge within its limits.
    #[error("geodesic did not converge: {0}")]
    GeodesicNonConvergence(String),

    /// The requested point lies outside the metric's valid domain.
    #[error("point outside domain: {0}")]
    OutsideDomain(String),

    /// No closed form is implemented for this quantity at these parameters
    /// (e.g. an ISCO with both spin and charge); the numerical geodesic path
    /// covers those cases instead.
    #[error("unsupported: {0}")]
    Unsupported(String),
}

/// Result alias for fallible relativity computations.
pub type Result<T> = core::result::Result<T, RelativityError>;
