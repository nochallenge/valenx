//! # valenx-flywheel
//!
//! Closed-form sizing math for flywheel (rotational kinetic) energy
//! storage and engine-governor speed-smoothing.
//!
//! ## What
//!
//! Given a rotor geometry and an operating speed, this crate computes
//! four things. (1) The **stored kinetic energy** `E = 1/2 I omega^2`.
//! (2) The rotor **mass moment of inertia** `I` for the canonical
//! axisymmetric shapes (solid disk, thin ring/rim, annular disk, point
//! mass, and a thin solid rod about its centre). (3) The **usable /
//! extractable energy** between two operating speeds, and the
//! governor-sizing form `dE = I omega_avg^2 Cs` driven by a target
//! coefficient of fluctuation `Cs`. (4) A first-order **rim (hoop)
//! stress** estimate `sigma = rho (omega r)^2 = rho v^2` for a thin
//! rotating ring.
//!
//! All quantities use SI base units: kilograms, metres, radians per
//! second, joules, pascals. Helpers convert revolutions per minute to
//! radians per second where convenient.
//!
//! ## Model
//!
//! The rotor is treated as a **rigid body of uniform density** spinning
//! about its geometric axis of symmetry. Inertias come straight from the
//! standard rigid-body integrals. A solid disk / cylinder about its axis
//! has `I = 1/2 m r^2`; a thin ring / rim with all mass at radius `r`
//! has `I = m r^2`; an annular (hollow) disk between radii `r_in` and
//! `r_out` has `I = 1/2 m (r_in^2 + r_out^2)`; a thin uniform rod of
//! length `L` about its centre has `I = 1/12 m L^2`; and a point mass at
//! radius `r` has `I = m r^2`.
//!
//! Energy follows the rotational analogue of `1/2 m v^2`, namely
//! `E = 1/2 I omega^2`. Because energy scales with `omega^2`, **doubling
//! the angular speed quadruples the stored energy** — a property the
//! test suite pins explicitly.
//!
//! Governor sizing uses the coefficient of fluctuation
//! `Cs = (omega_max - omega_min) / omega_avg`. The energy released as the
//! rotor slows from `omega_max` to `omega_min` is exactly
//! `dE = 1/2 I (omega_max^2 - omega_min^2)`, which factorises to
//! `dE = I omega_avg^2 Cs` using `omega_avg = (omega_max + omega_min)/2`.
//!
//! The rim-stress model is the elementary thin-rotating-ring result: a
//! ring of material density `rho` whose particles travel at rim speed
//! `v = omega r` sustains a hoop (tangential) tensile stress
//! `sigma = rho v^2`. This is the radius-independent leading-order term
//! that sets the burst speed of rim-type flywheels.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are **textbook closed-form models**
//! with strong idealising assumptions: a perfectly rigid rotor of
//! uniform density, inertias evaluated for ideal axisymmetric shapes,
//! and a thin-rim hoop-stress formula that ignores radial stress,
//! Poisson coupling, disk-thickness gradients, the press-fit / hub
//! interface, fibre-composite anisotropy, fatigue, creep, temperature,
//! windage, and bearing losses. Real flywheel rotors are sized with
//! detailed elasticity / FEA, validated safety factors, and burst
//! testing.
//!
//! This crate is **NOT a clinical, medical, or production engineering
//! tool** and must not be used to certify, qualify, or operate a
//! physical rotor where failure could cause injury. Use it for learning,
//! first-order feasibility estimates, and as a sanity-check oracle for
//! higher-fidelity solvers.
//!
//! ## Surface
//!
//! [`Rotor`] is the rotor geometry / mass enum, with validated
//! constructors and [`Rotor::moment_of_inertia`]. [`Flywheel`] pairs a
//! rotor with a material density and exposes energy / stress /
//! governor-sizing methods. The free functions [`kinetic_energy`],
//! [`usable_energy`], and [`flywheel_inertia_for_energy`] give
//! geometry-independent forms of the core relations;
//! [`coefficient_of_fluctuation`] and [`energy_fluctuation`] cover
//! governor sizing; and [`rim_stress`], [`rim_speed`], [`rpm_to_rad_s`],
//! [`rad_s_to_rpm`] are the stress and unit helpers. [`FlywheelError`]
//! and [`ErrorCategory`] form the error taxonomy.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod energy;
pub mod error;
pub mod flywheel;
pub mod rotor;
pub mod stress;

pub use energy::{
    coefficient_of_fluctuation, energy_fluctuation, flywheel_inertia_for_energy, kinetic_energy,
    rad_s_to_rpm, rpm_to_rad_s, usable_energy,
};
pub use error::{ErrorCategory, FlywheelError};
pub use flywheel::Flywheel;
pub use rotor::Rotor;
pub use stress::{rim_speed, rim_stress};
