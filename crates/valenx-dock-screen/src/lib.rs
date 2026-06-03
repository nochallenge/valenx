//! # valenx-dock-screen — docking & virtual screening
//!
//! A native-Rust extension of the [`valenx_dock`] AutoDock Vina port
//! into a full classical docking and structure-based virtual-screening
//! toolkit: a replacement for the classical core of AutoDock 4,
//! AutoDock Vina, smina, rDock and DOCK 6, plus honest subprocess
//! adapters for the neural-network tools (AlphaFold, ProteinMPNN,
//! DiffDock, RELION and friends).
//!
//! ## What it does
//!
//! - **Preparation** ([`prep`]) — a receptor grid-box definition, a
//!   ligand torsion-tree extractor (rotatable-bond detection), a
//!   pH-aware Gasteiger-style protonation / partial-charge assignment
//!   (reusing [`valenx_cheminf`]), and flexible-sidechain selection.
//! - **Scoring** ([`score`]) — a Vina-class empirical scoring function
//!   (gauss / repulsion / hydrophobic / H-bond + the rotatable-bond
//!   penalty, built on [`valenx_dock`]'s frozen Vina weights), an
//!   AutoDock4-class 12-6 force-field scoring function with directional
//!   H-bonds, electrostatics and desolvation, per-receptor-atom-type
//!   affinity-map precomputation and trilinear-interpolated fast grid
//!   scoring.
//! - **Search** ([`search`]) — an AutoDock4-class Lamarckian genetic
//!   algorithm, a Monte-Carlo / simulated-annealing pose search, a
//!   Vina-class iterated local search (reusing [`valenx_dock`]'s BFGS),
//!   and rigid + flexible docking drivers.
//! - **Screening** ([`mod@screen`]) — batch virtual screening over a
//!   ligand library, RMSD-based pose clustering and ranking, ensemble
//!   docking against multiple receptor conformations, consensus
//!   scoring with rank aggregation, and an MM-GBSA-class
//!   generalized-Born rescoring.
//! - **Analysis** ([`analyze`]) — a geometric fpocket-class
//!   binding-pocket detector, a protein-ligand interaction
//!   fingerprint, pose-RMSD redocking validation, pharmacophore-based
//!   and ADMET-lite library filtering (reusing [`valenx_cheminf`]),
//!   PDBQT-style docking result I/O, and a covalent-docking driver.
//! - **Adapters** ([`adapters`]) — typed request / result structs and
//!   `run_*` functions for the neural-network tools. **These never
//!   reimplement a model.** Each `run_*` builds a subprocess command
//!   and returns [`error::DockScreenError::ToolNotAvailable`]
//!   when the external binary is absent from `PATH`.
//! - **Drivers** ([`driver`]) — top-level [`driver::dock`] /
//!   [`driver::screen`] entry points, [`driver::DockingReport`]
//!   / [`driver::ScreeningReport`], and an
//!   [`driver::AdapterRegistry`] that probes which
//!   external tools are on `PATH`.
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, DockScreenError>`](error::DockScreenError). The error
//! type carries stable [`code`](error::DockScreenError::code) and
//! [`category`](error::DockScreenError::category) accessors for
//! telemetry.
//!
//! ## The "no neural networks" rule
//!
//! This crate ships **zero
//! reimplemented neural networks**. AlphaFold, ColabFold, ESMFold,
//! RoseTTAFold, RFdiffusion, ProteinMPNN, DiffDock, GNINA, RELION and
//! cryoSPARC are all *adapter-only*: [`adapters`] provides a typed
//! request / result struct and a `run_*` function that constructs the
//! subprocess command and shells out to the real tool when it is
//! installed. The classical AutoDock-class docking — scoring, the
//! genetic algorithm, Monte Carlo, iterated local search, grid maps —
//! *is* implemented as a real working v1.
//!
//! ## v1 scope
//!
//! This is a real working v1, not production parity with the 30-year
//! reference implementations. Each module documents its own
//! simplifications inline; the notable ones: the AutoDock4 force field
//! uses the published 12-6 / 12-10 functional forms with a built-in
//! parameter table, not the full AD4.2 `.dat` parameterisation;
//! protonation is a pH-rule + Gasteiger-style charge model, not a
//! Poisson-Boltzmann pKa solver; the genetic algorithm is a real
//! population / crossover / mutation / local-search loop but is not
//! tuned to AutoDock's exact operator schedule; the MM-GBSA rescoring
//! uses a single-snapshot generalized-Born / SASA model, not an MD
//! ensemble; the pocket detector is a geometric grid / sphere method
//! in the spirit of fpocket, not the alpha-sphere Voronoi
//! construction; covalent docking anchors the ligand at a reactive
//! receptor atom without bond-formation chemistry.

#![forbid(unsafe_code)]

pub mod adapters;
pub mod analyze;
pub mod driver;
pub mod error;
pub mod prep;
pub mod score;
pub mod screen;
pub mod search;

// --- Convenience re-exports of the most-used types --------------------

pub use driver::{dock, screen, AdapterRegistry, DockingReport, ScreeningReport};
pub use error::{DockScreenError, ErrorCategory, Result};
pub use prep::{GridBox, TorsionTree};
pub use score::{ScoringFunction, VinaTerms};

#[cfg(test)]
mod tests {
    use super::*;

    /// A minimal one-atom PDBQT ligand reused across crate-level tests.
    const MINI_LIGAND: &str = "\
ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";

    /// A minimal two-atom PDBQT receptor.
    const MINI_RECEPTOR: &str = "\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ATOM      2  CB  ALA A   1       1.500   0.000   0.000  1.00  0.00     0.000 C
";

    #[test]
    fn crate_error_reexported() {
        let e = DockScreenError::not_yet("x");
        assert_eq!(e.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn end_to_end_dock_smoke() {
        // The top-level driver parses receptor + ligand, builds grids,
        // searches, clusters and returns a report — all in-process,
        // no subprocess, no file dialogs.
        let report =
            dock(MINI_RECEPTOR, MINI_LIGAND, &driver::DockParams::fast()).expect("dock runs");
        assert!(!report.poses.is_empty(), "expected at least one pose");
    }

    #[test]
    fn grid_box_reexported() {
        let gb = GridBox::new([0.0, 0.0, 0.0], [20.0, 20.0, 20.0]);
        assert!(gb.is_ok());
    }
}
