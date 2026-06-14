//! # valenx-muscle
//!
//! Native **skeletal-muscle mechanics** — a Hill-type muscle model plus a
//! **myostatin loss-of-function** parameterization that answers the question
//! "if you remove the muscle-growth brake, do you get *bigger* muscle or
//! *stronger* muscle?" (Spoiler, and the whole point: those are not the same.)
//!
//! ## The model
//!
//! A classic **Hill-type** contractile model ([`hill`]). A muscle's force is
//! its peak isometric force `F₀ = σ · PCSA` (specific tension ×
//! physiological cross-sectional area) scaled by three dimensionless factors:
//!
//! - **active force-length** ([`active_force_length`]) — force peaks at the
//!   optimal fibre length and falls off as the muscle shortens or stretches
//!   (Gordon-Huxley-Julian 1966; modelled here as a Gaussian),
//! - **force-velocity** ([`force_velocity`]) — A. V. Hill's 1938 hyperbola:
//!   force drops as a muscle shortens faster, and rises above `F₀` when it is
//!   forcibly lengthened,
//! - **passive elasticity** ([`passive_force_length`]) — the connective-tissue
//!   spring resisting stretch beyond the optimal length,
//!
//! all projected onto the tendon through the **pennation** angle. Joint
//! strength is then `torque = force · moment arm` ([`joint_torque`]).
//!
//! ## Bigger ≠ stronger — the myostatin lesson
//!
//! [`MyostatinKnockout`] applies what loss of myostatin actually does to a
//! muscle, parameterized from the animal literature: PCSA (size) goes **up**,
//! but specific tension (force *per unit area*) goes **down**, fibres shift
//! toward fast / fatigable types, and oxidative capacity drops (Amthor et al.,
//! *PNAS* 2007). [`compare_wild_vs_knockout`] then shows the headline result:
//! the muscle can gain a lot of **size** while gaining far less **strength** —
//! exactly why myostatin-blocking drugs increased lean mass but disappointed
//! on strength in human trials. To get genuinely stronger you have to improve
//! the *whole chain* — specific tension, neural drive (voluntary activation),
//! tendon — not just lift the size cap.
//!
//! ```
//! use valenx_muscle::{compare_wild_vs_knockout, Muscle, MyostatinKnockout};
//!
//! let quad = Muscle::human_quadriceps();
//! let cmp = compare_wild_vs_knockout(&quad, &MyostatinKnockout::null());
//! // Big size jump, much smaller strength jump.
//! assert!(cmp.size_gain_pct > cmp.strength_gain_pct);
//! ```
//!
//! ## Honest scope
//!
//! Research / **educational** grade. The model classes are the textbook ones
//! (the same lineage as the Thelen-2003 / OpenSim Hill muscles), validated
//! here against their defining behaviours — the force-length peak, the Hill
//! force-velocity curve, specific tension ≈ 25 N/cm², and the
//! size-vs-strength dissociation. It is **not** a clinical predictor of any
//! specific person's outcome, and **not** a production musculoskeletal suite
//! (OpenSim). It models *principles*, not patients.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod hill;
pub mod myostatin;
pub mod strength;

pub use error::MuscleError;
pub use hill::{
    active_force_length, force_velocity, passive_force_length, Muscle,
    DEFAULT_SPECIFIC_TENSION_N_CM2,
};
pub use myostatin::MyostatinKnockout;
pub use strength::{compare_wild_vs_knockout, joint_torque, StrengthComparison};
