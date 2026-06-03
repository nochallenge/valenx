//! Reaction-network modelling — features 1 through 7.
//!
//! This module is the data-model layer of the crate. Everything else
//! (the ODE integrators, the stochastic simulators, the FBA solver,
//! the synthetic-biology designers) consumes the types defined here.
//!
//! - [`network`] — the central [`Model`]: species, reactions,
//!   compartments and global parameters, with structural validation
//!   and on-demand stoichiometry-matrix materialisation (feature 1).
//! - [`sbml`] — an SBML Level-3-core *subset* reader and writer that
//!   round-trips a kinetic reaction network (features 2 and 3).
//! - [`kinetics`] — the [`RateLaw`] enum: mass-action,
//!   Michaelis-Menten, Hill and constant-flux laws (features 4, 5, 6).
//! - [`rulebased`] — a BioNetGen-class rule expander that floods a
//!   site-state rule set into a flat [`Model`] (feature 7).

pub mod events;
pub mod expr;
pub mod kinetics;
pub mod network;
pub mod rulebased;
pub mod sbml;

pub use events::{AssignmentRule, EventAssignment, RateRule, SbmlEvent, SbmlRules, VarRef};
pub use expr::Expr;
pub use kinetics::RateLaw;
pub use network::{Compartment, Model, Parameter, Reaction, Species};
pub use rulebased::{Microstate, MoleculeTemplate, Rule, RuleModel, Site};
pub use sbml::{read_sbml, write_sbml, SbmlReadReport};
