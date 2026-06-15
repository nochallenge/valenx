//! # valenx-kinematics
//!
//! Planar (2-D) linkage and cam kinematics: the closed-form pieces of
//! classical machine theory that need no solver and no mesh.
//!
//! ## What
//!
//! Three self-contained models, depending only on `serde` (for the
//! data-type derives) and `thiserror` (for the error enum):
//!
//! - **Four-bar linkage** ([`FourBar`]) — Grashof mobility
//!   classification from the link lengths ([`FourBar::grashof_class`])
//!   and full loop-closure: given the crank angle, solve the coupler
//!   and rocker angles and the moving-pin positions on either assembly
//!   branch ([`FourBar::solve`]).
//! - **Cam follower** ([`CamRise`]) — the simple-harmonic and
//!   cycloidal *rise* displacement laws, with their first two cam-angle
//!   derivatives ([`CamRise::evaluate`]).
//! - **Errors** ([`KinematicsError`]) — validated constructors reject
//!   non-physical inputs and non-assemblable configurations up front.
//!
//! ## Model
//!
//! - *Grashof criterion.* With the four link lengths sorted so that
//!   `s` is shortest and `l` longest and `p`, `q` the middle pair, the
//!   linkage is a Grashof (crank) mechanism when `s + l < p + q`, a
//!   change-point mechanism when `s + l == p + q`, and a non-Grashof
//!   double-rocker when `s + l > p + q`.
//! - *Loop closure.* The vector loop `r2 + r3 = r1 + r4` (crank +
//!   coupler = ground + rocker) is solved geometrically: the crank pin
//!   is placed from the crank angle, and the coupler pin is found as a
//!   circle–circle intersection (the two solutions are the *open* and
//!   *crossed* assemblies). See [`fourbar`] for the full sign
//!   convention.
//! - *Cam laws.* For a rise of lift `L` over a cam rotation `β`, the
//!   SHM law is `s = (L/2)[1 − cos(πθ/β)]` and the cycloidal law is
//!   `s = L[θ/β − sin(2πθ/β)/2π]`. Cycloidal motion has zero velocity
//!   *and* zero acceleration at both ends; SHM has zero end velocity
//!   but a finite end acceleration. See [`cam`].
//!
//! ## Honest scope
//!
//! Research/educational grade. These are textbook closed-form planar
//! models from the standard machine-theory literature (Norton,
//! *Design of Machinery*; Shigley; Erdman & Sandor). They cover
//! idealised, rigid, frictionless, single-degree-of-freedom planar
//! mechanisms only. There is no spatial (3-D) linkage support, no
//! dynamics (forces, inertia, balancing), no contact/pressure-angle or
//! undercutting analysis for the cam, and no tolerance or
//! manufacturing modelling. This crate is NOT a clinical/medical or
//! production engineering tool: do not use it as the sole basis for a
//! load-bearing or safety-critical mechanism design.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod cam;
pub mod error;
pub mod fourbar;

pub use cam::{CamMotion, CamRise, FollowerState};
pub use error::{ErrorCategory, KinematicsError};
pub use fourbar::{Assembly, FourBar, FourBarPose, GrashofClass};
