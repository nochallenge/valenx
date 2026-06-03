//! Prime editing — Group C (features 11–15).
//!
//! Prime editing writes a programmed substitution, insertion or
//! deletion with a Cas9 nickase–reverse-transcriptase fusion guided
//! by a pegRNA — no double-strand break, no donor DNA. This block
//! designs prime-editing experiments:
//!
//! - [`editor`] — feature 11: the prime-editor database
//!   (PE2, PE3, PE3b, PEmax).
//! - [`pegrna`] — features 12–14: pegRNA design (spacer + PBS + RT
//!   template) for a substitution / insertion / deletion
//!   ([`pegrna::design_pegrna`]), a PBS / RT-template length-optimisation
//!   scan ([`pegrna::scan_pbs_rt`]) and PE3 / PE3b nicking-guide design
//!   ([`pegrna::design_nicking_guide`]).
//! - [`strategy`] — feature 15: prime-editing strategy modelling and a
//!   transparent efficiency heuristic ([`strategy::prime_efficiency`],
//!   [`strategy::model_strategy`]).
//!
//! ## v1 scope
//!
//! The pegRNA designer requires the spacer / PAM on the forward strand
//! of the supplied reference window. Both the length-scan score and
//! the whole-strategy efficiency score are transparent
//! feature-weighted heuristics (PBS Tm optimum, RT-template length
//! preference, homopolymer penalties, editor architecture, edit
//! size), documented as such — there is no trained model anywhere in
//! this block.

pub mod editor;
pub mod pegrna;
pub mod strategy;
