//! Subprocess adapters for the neural-network tools.
//!
//! # The "no neural networks" rule
//!
//! **This crate ships zero reimplemented neural networks.** AlphaFold,
//! ColabFold, ESMFold, RoseTTAFold, OmegaFold, Boltz, Chai,
//! RFdiffusion, ProteinMPNN, ESM-IF, Chroma, DiffDock, GNINA, RELION,
//! cryoSPARC and EMAN2 are all *external programs* — they need trained
//! weights, GPUs, and decades of one team's engineering. Valenx never
//! reimplements them.
//!
//! What this module provides instead is honest *subprocess adapter
//! scaffolding*: for each tool, a typed request struct, a typed result
//! struct, and a `run_*` function that
//!
//! 1. locates the external binary on `PATH` (with a Windows `.exe` /
//!    `PATHEXT` fallback);
//! 2. if it is missing, returns
//!    [`DockScreenError::ToolNotAvailable`](crate::error::DockScreenError::ToolNotAvailable)
//!    with an install hint — the honest failure mode;
//! 3. if it is present, constructs the command line that would invoke
//!    it ([`AdapterCommand`]) so the caller (or the Valenx job runner)
//!    can execute it.
//!
//! The `run_*` functions deliberately *return the command* rather than
//! spawning it — spawning an external docking / prediction tool is a
//! long-running, GPU-bound job that belongs in the application's job
//! runner, not in a library call. This also keeps the whole module
//! unit-testable (no subprocess is ever launched from a unit test).
//!
//! The modules:
//!
//! - [`structure_prediction`] (feature 26) — AlphaFold / ColabFold /
//!   ESMFold / RoseTTAFold / OmegaFold / Boltz / Chai.
//! - [`generative_design`] (feature 27) — RFdiffusion / ProteinMPNN /
//!   ESM-IF / Chroma.
//! - [`nn_docking`] (feature 28) — DiffDock / GNINA.
//! - [`cryo_em`] (feature 29) — RELION / cryoSPARC / EMAN2.

pub mod common;
pub mod cryo_em;
pub mod generative_design;
pub mod nn_docking;
pub mod structure_prediction;

pub use common::{find_executable, AdapterCommand, ToolStatus};
pub use cryo_em::{run_cryo_em, CryoEmRequest, CryoEmResult, CryoEmTool};
pub use generative_design::{
    run_generative_design, GenerativeDesignRequest, GenerativeDesignResult, GenerativeTool,
};
pub use nn_docking::{run_nn_docking, NnDockingRequest, NnDockingResult, NnDockingTool};
pub use structure_prediction::{
    run_structure_prediction, StructurePredictionRequest, StructurePredictionResult,
    StructurePredictionTool,
};
