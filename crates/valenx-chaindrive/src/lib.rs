//! # valenx-chaindrive
//!
//! Closed-form **roller-chain drive** calculator: from a pair of
//! sprockets and a chain pitch, work out the speed ratio, the chain
//! velocity, the chain length you need to cut for a given shaft spacing,
//! and the torque the drive delivers to the output shaft.
//!
//! ## What
//!
//! A single-stage chain drive is a driver sprocket on the input shaft, a
//! driven sprocket on the output shaft, and one roller chain wrapped
//! around both. Describe it as a [`SprocketPair`] (driver teeth, driven
//! teeth, chain pitch) and this crate answers the four questions a chain
//! designer asks first:
//!
//! - the **gear ratio** [`SprocketPair::ratio`] — driven ÷ driver teeth;
//! - the **chain velocity** [`chain_velocity_m_per_s`] — how fast the
//!   links travel;
//! - the **driven-shaft speed** [`driven_speed_rpm`] and **output
//!   torque** [`output_torque_n_m`];
//! - the **chain length** [`chain_length_pitches`] — how many links to
//!   buy for a chosen centre distance.
//!
//! [`analyze`] rolls all of them into one [`DriveResult`].
//!
//! ```
//! use valenx_chaindrive::{analyze, SprocketPair};
//!
//! // 17-tooth driver, 34-tooth driven, ANSI 40 (12.7 mm) chain.
//! let pair = SprocketPair::new(17, 34, 12.7).expect("valid sprockets");
//! let result = analyze(&pair, 1000.0, 50.0, 500.0).expect("valid drive");
//!
//! assert!((result.ratio - 2.0).abs() < 1e-9); // 2:1 reduction
//! assert!((result.driven_speed_rpm - 500.0).abs() < 1e-9); // half speed
//! assert!((result.output_torque_n_m - 100.0).abs() < 1e-9); // double torque
//! ```
//!
//! ## Model
//!
//! Every formula is the standard textbook closed form for a roller-chain
//! drive; nothing is iterated or fitted.
//!
//! - **Ratio.** Because the same chain links engage both sprockets, the
//!   ratio is purely `i = N2 / N1` (driven over driver tooth count) and
//!   does not depend on the pitch.
//! - **Chain velocity.** The chain advances one pitch per tooth that
//!   passes engagement, so `v = p · z · n / 60_000` m/s for pitch `p`
//!   (mm), tooth count `z` and speed `n` (rev/min). Evaluated on either
//!   sprocket the result is identical — that identity *is* the ratio
//!   relation `n2 = n1 · N1 / N2`.
//! - **Chain length.** The classic approximation in pitches,
//!   `L ≈ 2C + (z1 + z2)/2 + (z2 − z1)² / (4π²C)`, with the centre
//!   distance `C` expressed in pitches; the buildable result is rounded
//!   up to a whole **even** link count.
//! - **Output torque.** Loss-free power conservation gives `T_out =
//!   T_in · i`: a reduction multiplies torque while dividing speed.
//! - **Pitch diameter.** `d = p / sin(π / z)`, the sprocket pitch circle
//!   used for the geometric-overlap (degenerate-centre-distance) check.
//!
//! ## Honest scope
//!
//! Research/educational grade: textbook closed-form / numerical models,
//! **NOT** a clinical/medical/production engineering tool. This crate
//! computes ideal single-stage kinematics and the loss-free torque upper
//! bound only. It deliberately does **not** model:
//!
//! - chain efficiency, friction, or transmission losses — the real
//!   output torque is lower than [`output_torque_n_m`] reports;
//! - chordal (polygonal) speed variation, so the chain velocity is the
//!   mean, not the instantaneous, value;
//! - chain tension, catenary sag, slack-side dynamics, or required idler
//!   take-up;
//! - power / fatigue ratings, lubrication, wear elongation, or service
//!   factors;
//! - multi-strand, multi-sprocket, or serpentine layouts (single driver
//!   and single driven only).
//!
//! Do not size a real power-transmission chain from these numbers
//! without applying the manufacturer ratings and service factors.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod drive;
pub mod error;
pub mod spec;

pub use drive::{
    analyze, chain_length_pitches, chain_length_pitches_exact, chain_velocity_m_per_s,
    driven_speed_rpm, output_torque_n_m, DriveResult,
};
pub use error::{ChainDriveError, ErrorCategory};
pub use spec::{SprocketPair, MIN_SPROCKET_TEETH};
