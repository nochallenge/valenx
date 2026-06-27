//! # valenx-survivability
//!
//! **Defensive survivability & protection** analysis — the M7 track of the
//! Valenx defense modeling-&-simulation roadmap
//! (`docs/design/2026-06-23-valenx-defense-ms-roadmap.md`). It answers one
//! question, and only that question:
//!
//! > **Given a threat (a blast at a stand-off, or an impacting projectile),
//! > what is the minimum protection that lets a structure and its occupants
//! > *survive*?**
//!
//! This is the same physics as **civil blast-resistant building design**
//! (UFC 3-340-02, ASCE 59) and **automotive crash safety** — the *protective*
//! side of the shared impact/blast physics. Every output is framed as
//! "minimum protection to survive threat X".
//!
//! ## Dual-use boundary (hard line)
//!
//! This crate is **protective / defensive only**. It models how to *protect*
//! against a threat, **never** how to penetrate, defeat, or optimize a weapon
//! against protection. Concretely, there are **no** warhead, lethality,
//! fragmentation, shaped-charge, or penetration-optimization models. The blast
//! *source* appears only as the standard idealized far-field overpressure curve
//! used as a design *load*; the armor models size the *minimum plate that keeps
//! a threat out*, never a means of getting through one. If a question is "how do
//! I defeat / penetrate X", it is out of scope for this crate by construction.
//!
//! ## What it computes
//!
//! - **Blast loading** ([`blast`]) — free-field air-blast peak side-on
//!   overpressure and positive-phase impulse versus the Hopkinson–Cranz scaled
//!   distance `Z = R / W^(1/3)`, via the **Brode (1955)** overpressure fit and
//!   the **Newmark–Hansen (1961)** impulse relation, assembled into a
//!   Friedlander pressure pulse. See [`BlastLoad`].
//! - **Protective structural response** ([`response`]) — the single-degree-of-
//!   freedom dynamic response of a protective element (wall / plate / panel) to
//!   the blast pulse: peak deflection, ductility ratio, dynamic amplification
//!   ([`sdof_response`]), and the canonical **pressure–impulse (P–I)
//!   iso-damage diagram** ([`PressureImpulseDiagram`]) with its impulsive and
//!   quasi-static asymptotes. The time integration **reuses the validated
//!   Friedlander / Newmark-β SDOF transient solver in [`valenx_fem`]** rather
//!   than reimplementing it.
//! - **Armor / impact protection sizing** ([`armor`]) — the **minimum plate
//!   thickness and areal density** that just defeats a projectile of a given
//!   kinetic energy, via a plug-shear ballistic-limit energy balance. See
//!   [`minimum_protection`].
//! - **Occupant survivability** ([`occupant`]) — a simple peak-acceleration
//!   injury screen: peak g versus a supplied whole-body tolerance limit. See
//!   [`assess_occupant`].
//!
//! ## Reuse
//!
//! The structural-response core (the [`valenx_fem::FriedlanderPulse`],
//! [`valenx_fem::solve_sdof_blast`], and the Newmark-β
//! [`valenx_fem::solve_transient_response`]) already exists in `valenx-fem`'s
//! defensively-scoped `transient` module; this crate consumes it. The *new*
//! contributions here are the blast-loading correlations, the P–I diagram, the
//! armor sizing, and the occupant screen.
//!
//! ## Status
//!
//! **Research / educational grade, validation-pending.** Every model is a
//! well-established open-literature protective-design equation, implemented with
//! `#![forbid(unsafe_code)]` and fail-loud input validation (degenerate inputs —
//! zero charge, stand-off, mass, area, or tolerance — return `Err`, never a
//! `NaN` or a panic). It is **not** a certified blast or ballistic engineering
//! tool; absolute numbers depend on uncertain dynamic material data and the
//! assumed failure mechanism, so treat the outputs as fast trade-study
//! *screens*.
//!
//! ```
//! use valenx_survivability::{BlastLoad, PressureImpulseDiagram, Threat, minimum_protection};
//!
//! // Free-field load of 100 kg TNT-equivalent at 20 m stand-off.
//! let load = BlastLoad::tnt_free_air(100.0, 20.0).expect("physical");
//! assert!(load.peak_overpressure_pa > 0.0);
//! assert!(load.impulse_pa_s > 0.0);
//!
//! // P–I iso-damage diagram for an elastic protective panel.
//! let pi = PressureImpulseDiagram::elastic(500.0, 5.0e6, 0.02).expect("physical");
//! assert!(pi.pressure_asymptote_pa > 0.0); // quasi-static asymptote
//! assert!(pi.impulse_asymptote_pa_s > 0.0); // impulsive asymptote
//!
//! // Minimum protective plate to stop a 10 g projectile at 800 m/s.
//! let threat = Threat::new(0.01, 800.0, 0.008).expect("physical");
//! let sizing = minimum_protection(&threat, 7850.0, 5.0e8).expect("physical");
//! assert!(sizing.min_thickness_m > 0.0);
//! ```

#![forbid(unsafe_code)]

pub mod armor;
pub mod blast;
pub mod error;
pub mod occupant;
pub mod response;

pub use armor::{minimum_protection, ArmorSizing, Threat};
pub use blast::{
    peak_overpressure_pa, positive_impulse_pa_s, scaled_distance, BlastLoad, Z_MAX, Z_MIN,
};
pub use error::SurvivabilityError;
pub use occupant::{assess_occupant, OccupantAssessment, G0};
pub use response::{sdof_response, PiPoint, PressureImpulseDiagram, SdofResponse};

// Re-export the reused FEM transient types so a survivability caller has the
// whole pipeline (pulse + controls) in one place without depending on the names
// of the FEM crate's internal module layout.
pub use valenx_fem::{FriedlanderPulse, NewmarkBeta, TransientControls, TransientResponse};
