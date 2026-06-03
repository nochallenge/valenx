//! The fourteen Genetics-workbench panels.
//!
//! One module per Round-6 computational-biology crate. Each module
//! exposes a `*Panel` state struct (a field on
//! [`crate::genetics_workbench::GeneticsWorkbenchState`]) and a `draw`
//! function that renders the panel into the workbench side panel.
//!
//! Every panel follows the same shape — a tool palette, input forms,
//! a Run action calling the crate's real native API, and a results
//! area — using the shared helpers in [`common`]. The exception is
//! [`rna_designer`], an end-to-end RNA / mRNA design workbench whose
//! five sections span `valenx-rnastruct`, `valenx-rnadesign` and
//! `valenx-genediting`.

pub mod common;
pub mod molecule_view;

pub mod alignment;
pub mod biostruct;
pub mod cheminf;
pub mod docking;
pub mod genediting;
pub mod genomics;
pub mod md;
pub mod phylogenetics;
pub mod popgen;
pub mod qchem;
pub mod rna_designer;
pub mod rnastruct;
pub mod sequence;
pub mod structpredict;
pub mod sysbio;

use crate::genetics_workbench::GeneticsPanel;
use crate::ValenxApp;

/// Run the primary action of the currently-active genetics panel.
/// Wired into the global Ctrl+R shortcut by `update::dispatch_shortcut`
/// so a keyboard run-press routes to the right panel.
///
/// Panels without a clear "primary action" (read-only result displays)
/// are no-ops — the shortcut handler ignores the no-op silently.
pub fn run_active_panel(app: &mut ValenxApp) {
    match app.genetics.active {
        GeneticsPanel::Sequence => sequence::run_primary_shortcut(app),
        GeneticsPanel::RnaStructure => rnastruct::run_primary_shortcut(app),
        GeneticsPanel::Cheminformatics => cheminf::run_primary_shortcut(app),
        GeneticsPanel::Alignment => alignment::run_primary_shortcut(app),
        // The rest of the panels' Run buttons are still mouse-only —
        // wiring them is a follow-up; the keyboard-cheat-sheet entry
        // notes "Ctrl+R: active panel's main button (where wired)".
        _ => {}
    }
}
