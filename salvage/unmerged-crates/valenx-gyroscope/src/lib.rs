//! # valenx-gyroscope
//!
//! ## What
//!
//! A tiny, self-contained library for the gyroscopic couple of a
//! symmetric rigid rotor and the two textbook applications most often
//! worked in a *Theory of Machines* course: an aeroplane making a turn
//! and a ship steering or pitching. It computes the couple magnitude, the
//! precession rate, the spin angular momentum, and the active-versus-
//! reactive couple sense, all from validated inputs.
//!
//! ## Model
//!
//! A rotor of polar mass moment of inertia `I` (kg*m^2) spinning at
//! angular velocity `w` (rad/s) carries spin angular momentum `L = I*w`
//! along its axis. Forcing that axis to precess at angular velocity `wp`
//! (rad/s) about a perpendicular axis changes the momentum vector at the
//! rate
//!
//! `C = dL/dt = I * w * wp`,
//!
//! the magnitude of the **active** couple that must be applied. The rotor
//! reacts with an equal and opposite **reactive** couple on its frame.
//! Rearranging gives the precession rate `wp = C / (I*w)` and the spin
//! `w = C / (I*wp)`. In the application cases the precession rate is set
//! by the vehicle motion: a level turn gives `wp = v / R`, steady steering
//! gives `wp` directly, and simple-harmonic pitching of amplitude `phi0`
//! and period `T` gives a peak `wp_max = phi0 * (2*pi / T)`.
//!
//! ## Honest scope
//!
//! This is research/educational grade: standard textbook closed-form
//! models, validated against analytic ground truth (hand-worked
//! Khurmi/Rattan *Theory of Machines* couple numbers and known limiting
//! cases). It is NOT a clinical/medical or production-certified
//! engineering tool — there are no component tolerances, safety factors,
//! fatigue/thermal limits, bearing-load ratings, or code-compliance
//! checks. Do not size real hardware with it.
//!
//! ## Example
//!
//! ```rust
//! use valenx_gyroscope::{
//!     angular_momentum, gyroscopic_couple, precession_rate, rpm_to_rad_s,
//! };
//!
//! // A propeller: I = 20 kg*m^2 spinning at 1000 rpm.
//! let inertia = 20.0_f64;
//! let spin = rpm_to_rad_s(1000.0).unwrap(); // ~104.72 rad/s
//!
//! // The aeroplane turns: v = 50 m/s on a 100 m radius => wp = 0.5 rad/s.
//! let precession = 50.0 / 100.0;
//! let couple = gyroscopic_couple(inertia, spin, precession).unwrap();
//! assert!((couple - 1047.197551).abs() < 1e-4);
//!
//! // The relation is invertible: recover wp from the couple.
//! let wp = precession_rate(couple, inertia, spin).unwrap();
//! assert!((wp - precession).abs() < 1e-9);
//!
//! // Spin angular momentum L = I*w.
//! let l = angular_momentum(inertia, spin).unwrap();
//! assert!((l - inertia * spin).abs() < 1e-9);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod applications;
pub mod couple;
pub mod error;

pub use applications::{
    aeroplane_turn_couple, ship_pitching_peak_couple, ship_steering_couple, CoupleSense,
};
pub use couple::{
    angular_momentum, gyroscopic_couple, precession_rate, rad_s_to_rpm, rpm_to_rad_s,
    spin_for_couple,
};
pub use error::{GyroError, Result};
