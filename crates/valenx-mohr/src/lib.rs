//! # valenx-mohr
//!
//! 2D plane-stress transformation and **Mohr's circle**: turn the three
//! in-plane stress components into the principal stresses, the maximum
//! in-plane shear, the principal orientation, and the stress acting on
//! any inclined plane.
//!
//! ## What
//!
//! Build a [`StressState`] from the normal stresses `sx`, `sy` and the
//! in-plane shear `txy`, then read off the standard results:
//!
//! - [`StressState::principal_stresses`] — the principal stresses
//!   `s1 >= s2`.
//! - [`StressState::max_shear`] — the maximum in-plane shear stress.
//! - [`StressState::absolute_max_shear`] — the absolute (out-of-plane)
//!   maximum shear stress, accounting for the zero third principal
//!   stress of a plane-stress state (the Tresca-governing value).
//! - [`StressState::principal_angle`] — the orientation (radians) of
//!   the plane carrying `s1`.
//! - [`StressState::stress_on_plane`] — the normal and shear stress on
//!   a plane inclined at an arbitrary angle.
//! - [`StressState::mean_normal`] / [`StressState::radius`] — the centre
//!   and radius of Mohr's circle.
//!
//! ```
//! use valenx_mohr::StressState;
//!
//! let s = StressState::new(-20.0, 90.0, 60.0).expect("finite inputs");
//! let p = s.principal_stresses();
//! println!("s1 = {:.2}, s2 = {:.2}", p.s1, p.s2);
//! println!("max shear = {:.2}", s.max_shear());
//! // The shear vanishes on the principal plane:
//! let on = s.stress_on_plane(s.principal_angle()).expect("finite angle");
//! assert!(on.shear.abs() < 1e-6);
//! ```
//!
//! ## Model
//!
//! For the symmetric in-plane stress tensor
//! `[[sx, txy], [txy, sy]]`, the principal stresses are
//!
//! ```text
//! s1,2 = (sx + sy)/2 +/- sqrt(((sx - sy)/2)^2 + txy^2)
//! ```
//!
//! which are exactly the eigenvalues of that tensor. The mean normal
//! stress `(sx + sy)/2` is the circle centre, the square-root term is
//! the circle radius `R` (equal to the maximum in-plane shear), and the
//! principal direction satisfies `tan(2 theta_p) = 2 txy / (sx - sy)`.
//!
//! The stress on a plane whose outward normal is rotated `theta`
//! counter-clockwise from the `x` axis follows the double-angle
//! transformation
//!
//! ```text
//! sx'  = (sx + sy)/2 + (sx - sy)/2 * cos(2 theta) + txy * sin(2 theta)
//! txy' = -(sx - sy)/2 * sin(2 theta) + txy * cos(2 theta)
//! ```
//!
//! so the sum `sx' + sy'` is invariant (`= sx + sy = s1 + s2`) and every
//! transformed `(sx', txy')` pair lies on Mohr's circle.
//!
//! Treating the plane-stress state as 3D with a zero out-of-plane
//! principal stress (`s3 = 0`), the absolute maximum shear is
//!
//! ```text
//! tau_abs = (max(s1, s2, 0) - min(s1, s2, 0)) / 2,
//! ```
//!
//! which equals the in-plane `(s1 - s2)/2` only when `s1` and `s2`
//! straddle zero, and otherwise exceeds it.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are the textbook closed-form,
//! well-established stress-transformation formulae of plane-stress
//! mechanics of materials — the same relations taught in any
//! strength-of-materials course. The crate covers the **2D in-plane**
//! case plus the plane-stress absolute maximum shear (which uses the
//! zero out-of-plane principal stress); it does not compute the full 3D
//! principal-stress problem for a general triaxial state, failure-
//! criterion checks, plasticity, or stress-concentration factors. It is
//! NOT a clinical, medical, or production engineering tool and must not
//! be used as the sole basis for safety-critical design.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod stress;

pub use error::{ErrorCategory, MohrError};
pub use stress::{PlaneStress, PrincipalStresses, StressState};
