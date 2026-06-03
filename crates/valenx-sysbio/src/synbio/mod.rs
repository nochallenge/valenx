//! Synthetic biology — features 24 through 29.
//!
//! The design-and-build half of the crate: data models and planners
//! for engineered genetic systems, as opposed to the analysis of
//! natural ones.
//!
//! - [`sbol`] — an SBOL-class genetic-design data model: parts,
//!   components, devices and sequence annotations (feature 24).
//! - [`circuit`] — a genetic logic-circuit model with Cello-class
//!   Boolean and analog truth-table simulation (feature 25).
//! - [`grn`] — a Hill-kinetics gene-regulatory-network ODE model that
//!   compiles to a reaction network and simulates (feature 26).
//! - [`assembly`] — Gibson, Golden-Gate and BioBrick DNA-assembly
//!   planners (features 27 and 28).
//! - [`library`] — a standard part library, sequence annotation and
//!   codon-context design helpers (feature 29).

pub mod assembly;
pub mod circuit;
pub mod grn;
pub mod library;
pub mod sbol;

pub use assembly::{
    design_gibson_fragments, plan_biobrick, plan_gibson, plan_golden_gate,
    plan_golden_gate_checked, wallace_tm, BioBrickPlan, GibsonPlan, GoldenGatePlan,
};
pub use circuit::{Circuit, Gate, InputSignal, LogicNode, LogicOp, Wire};
pub use grn::{Gene, GeneNetwork, RegEdge, Regulation};
pub use library::{
    annotate_sequence, codon_adaptation_index, relative_codon_frequency,
    standard_part_library, Annotation, PartLibrary,
};
pub use sbol::{Component, Device, Module, Part, PartRole, SequenceAnnotation};
