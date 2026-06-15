//! # valenx-clutch
//!
//! Torque-capacity and transmissible-power calculator for **dry plate /
//! disc friction clutches**, built on the two classic annular-friction
//! theories.
//!
//! ## What
//!
//! Describe a clutch as a friction-face [`clutch::ClutchGeometry`] (an
//! inner and outer radius), a coefficient of friction `mu`, and the
//! number of friction surfaces in contact `N`, then ask the resulting
//! [`clutch::FrictionClutch`] two questions:
//!
//! - how much **torque** can it carry before it slips, for a given
//!   axial clamp force `F`, under either contact-pressure idealisation?
//! - how much **power** can it transmit once engaged at angular velocity
//!   `omega`?
//!
//! ```
//! use valenx_clutch::{ClutchGeometry, FrictionClutch, PressureModel, rpm_to_rad_per_s};
//!
//! // Single-plate clutch: 100/200 mm face, mu = 0.3, both faces grip (N = 2).
//! let geom = ClutchGeometry::new(100.0, 200.0).unwrap();
//! let clutch = FrictionClutch::new(geom, 0.3, 2).unwrap();
//!
//! let clamp_n = 5_000.0; // axial clamp force, newtons
//! let t_wear = clutch.torque_uniform_wear(clamp_n).unwrap();
//! let t_pres = clutch.torque_uniform_pressure(clamp_n).unwrap();
//! assert!(t_pres > t_wear); // new-plate theory is the optimistic bound
//!
//! let omega = rpm_to_rad_per_s(3_000.0).unwrap();
//! let power_w = clutch
//!     .power(PressureModel::UniformWear, clamp_n, omega)
//!     .unwrap();
//! assert!(power_w > 0.0);
//! ```
//!
//! ## Model
//!
//! For an annular friction face of inner radius `ri` and outer radius
//! `ro`, axial clamp force `F`, coefficient of friction `mu`, and `N`
//! friction surfaces in contact:
//!
//! ```text
//! uniform wear:     T = mu * F * N * (ro + ri) / 2
//! uniform pressure: T = (2/3) * mu * F * N * (ro^3 - ri^3) / (ro^2 - ri^2)
//! transmitted power: P = T * omega
//! ```
//!
//! Both torque laws are `mu * F * N * r_eff`, differing only in the
//! effective lever arm `r_eff`: the **arithmetic** mean radius for a
//! worn-in face (uniform wear, conservative) versus the larger
//! **area-weighted** mean radius for a new flat face (uniform pressure,
//! optimistic). The full derivation lives in the [`clutch`] module
//! documentation.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are the **textbook closed-form**
//! rigid-body, steady-state capacity equations (Shigley, Juvinall,
//! Norton) — the limiting torque a dry friction clutch can carry without
//! slipping. They do **not** model engagement transients, slip heating
//! or thermal fade, the friction-coefficient drop with temperature and
//! sliding speed, plate-flatness / pressure-distribution departures, wet
//! (lubricated) clutches, wear life, or fatigue. This crate is for
//! first-order sizing and teaching; it is **not** a clinical, medical,
//! or production engineering tool and is no substitute for a validated
//! clutch test rig or a certified design code.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod clutch;
pub mod error;

pub use clutch::{rpm_to_rad_per_s, ClutchGeometry, FrictionClutch, PressureModel};
pub use error::ClutchError;
