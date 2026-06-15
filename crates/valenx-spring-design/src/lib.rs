//! # valenx-spring-design
//!
//! Closed-form design calculator for round-wire **helical compression
//! springs**.
//!
//! ## What
//!
//! Given a spring's geometry (wire diameter `d`, mean coil diameter
//! `D`, active-coil count `N`) and material shear modulus `G`, this
//! crate computes the quantities a mechanical designer reaches for
//! first:
//!
//! - the spring **index** `C = D / d`,
//! - the **Wahl** curvature-correction factor `K_w`,
//! - the spring **rate** (stiffness) `k`,
//! - the axial **deflection** `delta` under a given force (and its
//!   inverse, the force for a given deflection),
//! - the Wahl-corrected torsional **shear stress** `tau` (and the
//!   uncorrected value for comparison).
//!
//! The whole surface is the single validated type
//! [`HelicalSpring`] plus the [`SpringError`] taxonomy.
//!
//! ## Model
//!
//! Standard linear-elastic helical-spring mechanics (Shigley / Wahl),
//! all closed form:
//!
//! ```text
//!   C   = D / d                                   (spring index)
//!   K_w = (4C - 1)/(4C - 4) + 0.615/C             (Wahl factor)
//!   k   = G d^4 / (8 D^3 N)                        (spring rate)
//!   delta = F / k                                  (deflection)
//!   tau = K_w * 8 F D / (pi d^3)                   (max shear stress)
//! ```
//!
//! Consequences worth noting and exercised in the tests: rate scales
//! as `d^4`, as `1/D^3`, and as `1/N`; the Wahl factor is always
//! greater than `1` and tends to `1` as the index grows; the shear
//! stress is linear in force and in the slenderness group `D / d^3`.
//!
//! The formulas are dimensionally consistent — use any single coherent
//! unit system. The examples use mm / N / MPa, giving a rate in `N/mm`,
//! deflection in `mm`, and stress in `MPa`.
//!
//! ## Honest scope
//!
//! Research/educational grade. This implements textbook closed-form
//! analytic models only. It is **NOT** a clinical, medical, or
//! production engineering tool. It deliberately omits everything a real
//! spring qualification requires: fatigue / endurance life, solid /
//! free / installed lengths and set, buckling and lateral stability,
//! end-coil effects on the *effective* active-coil count, pitch and
//! helix-angle corrections, temperature and relaxation, surge /
//! resonance, and manufacturing tolerances. Do not use it for
//! life-safety, load-bearing, or any production design without
//! independent verification by a qualified engineer.
//!
//! ## Example
//!
//! ```
//! use valenx_spring_design::HelicalSpring;
//!
//! // d = 2 mm wire, D = 16 mm mean coil, 10 active coils, spring steel.
//! let spring = HelicalSpring::new(2.0, 16.0, 10.0, 79_300.0).unwrap();
//!
//! assert!((spring.spring_index() - 8.0).abs() < 1e-12);
//! assert!(spring.wahl_factor() > 1.0);
//!
//! let k = spring.rate(); // N/mm
//! let delta = spring.deflection(40.0).unwrap(); // mm under 40 N
//! assert!((delta - 40.0 / k).abs() < 1e-12);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod spring;

pub use error::{ErrorCategory, SpringError};
pub use spring::HelicalSpring;
