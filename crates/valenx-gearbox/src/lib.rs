//! # valenx-gearbox
//!
//! Closed-form kinematics for gear trains.
//!
//! ## What
//!
//! Three composable building blocks for analysing power transmission
//! through gears:
//!
//! - [`GearStage`] — one meshing gear pair. Gear ratio
//!   `teeth_out / teeth_in`, with the inverse speed / torque transform
//!   and a mechanical efficiency.
//! - [`CompoundTrain`] — several stages in series. Overall ratio is the
//!   product of the stage ratios; overall efficiency the product of the
//!   stage efficiencies.
//! - [`PlanetarySet`] — a fixed-ring planetary (epicyclic) set with the
//!   sun as input and the carrier as output. Reduction ratio
//!   `1 + ring / sun` from the Willis fundamental train equation.
//!
//! All inputs flow through validated constructors that return
//! [`GearboxError`] rather than panicking.
//!
//! ## Model
//!
//! Each member is a rigid ideal gear: speed scales by `1 / ratio`,
//! torque scales by `ratio` and is then reduced by an efficiency factor
//! in `(0, 1]`. For a lossless member (`efficiency = 1`) power is
//! conserved, because the speed reduction and torque multiplication are
//! exact reciprocals. Tooth counts are integers; ratios, speeds, and
//! torques are dimensionless `f64` magnitudes in whatever consistent
//! unit system the caller chooses (rpm and N·m, rad/s and lbf·ft, ...).
//!
//! Governing relations:
//!
//! ```text
//! stage:      ratio = teeth_out / teeth_in
//! compound:   ratio = product of stage ratios
//! planetary:  ratio = 1 + ring / sun          (fixed ring, Willis)
//! speed:      out   = in / ratio
//! torque:     out   = in * ratio * efficiency
//! ```
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form
//! rigid-body gear models — exactly the formulas in an introductory
//! machine-design text. This crate is **not** a clinical/medical or
//! production engineering tool. It does **not** model tooth bending or
//! Hertzian contact stress, backlash or transmission error, gear
//! efficiency from first principles (efficiency is a user-supplied
//! knob), lubrication, thermal effects, dynamic / resonance behaviour,
//! load sharing among planets, or fatigue life. Do not use it to size,
//! certify, or sign off a real drivetrain.
//!
//! ## Example
//!
//! ```
//! use valenx_gearbox::{CompoundTrain, GearStage, PlanetarySet};
//!
//! // A single 4:1 reduction stage.
//! let stage = GearStage::new(10, 40)?;
//! assert!((stage.ratio() - 4.0).abs() < 1e-12);
//! assert!((stage.output_speed(1000.0) - 250.0).abs() < 1e-12);
//!
//! // Two stages in series: 4:1 then 3:1 -> 12:1 overall.
//! let train = CompoundTrain::new(vec![GearStage::new(10, 40)?, GearStage::new(10, 30)?])?;
//! assert!((train.ratio() - 12.0).abs() < 1e-12);
//!
//! // Fixed-ring planetary: sun 24, ring 72 -> 1 + 72/24 = 4:1.
//! let planetary = PlanetarySet::new(24, 72)?;
//! assert!((planetary.ratio() - 4.0).abs() < 1e-12);
//! # Ok::<(), valenx_gearbox::GearboxError>(())
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod compound;
pub mod error;
pub mod planetary;
pub mod stage;

pub use compound::CompoundTrain;
pub use error::{ErrorCategory, GearboxError};
pub use planetary::PlanetarySet;
pub use stage::GearStage;
