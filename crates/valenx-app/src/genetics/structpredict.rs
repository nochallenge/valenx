//! Panel 15 — **Structure Prediction** (`valenx-structpredict`).
//!
//! Classical, weights-free protein structure prediction + design — the
//! Modeller / Rosetta-class algorithms that need no trained network.
//! Two modes:
//!
//! - **Ab-initio fold** — fold a sequence from scratch by Monte-Carlo
//!   fragment assembly ([`valenx_structpredict::driver::predict_abinitio`]).
//! - **Fixed-backbone design** — redesign a sequence onto a supplied
//!   backbone ([`valenx_structpredict::driver::design_sequence`]).
//!
//! Honest framing, surfaced in the UI: these are real, useful classical
//! methods but they are **not** AlphaFold-accuracy — that accuracy comes
//! from a trained network's learnt co-evolutionary signal, shipped here
//! as subprocess adapters only.

use eframe::egui;

use valenx_structpredict::driver::{design_sequence, predict_abinitio, StructPredictReport};

use super::common;
use crate::ValenxApp;

/// Which structure-prediction job the panel runs.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(crate) enum SpMode {
    /// Ab-initio fragment-assembly fold from sequence.
    #[default]
    AbInitio,
    /// Fixed-backbone (inverse-folding) sequence design.
    Design,
}

/// Snapshot of every editable input the panel owns (for undo / redo).
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub(crate) struct StructPredictSnapshot {
    mode: SpMode,
    sequence: String,
    structure_text: String,
    chain: String,
    centroid_moves: usize,
    repack_moves: usize,
    design_moves: usize,
    seed: u64,
}

/// Form + result state for the Structure Prediction panel.
pub struct StructPredictPanel {
    mode: SpMode,
    /// Target sequence (ab-initio).
    sequence: String,
    /// Backbone PDB / mmCIF text (design).
    structure_text: String,
    /// Chain id to redesign (design).
    chain: String,
    /// Monte-Carlo move counts (ab-initio centroid + all-atom repack).
    centroid_moves: usize,
    repack_moves: usize,
    /// Monte-Carlo move count for the design search.
    design_moves: usize,
    /// RNG seed — fixes every stochastic step so runs reproduce.
    seed: u64,
    error: Option<String>,
    result: String,
    /// Undo / redo over the form inputs, reversing the last Run.
    history: crate::undo::History<StructPredictSnapshot>,
}

impl StructPredictPanel {
    fn snapshot(&self) -> StructPredictSnapshot {
        StructPredictSnapshot {
            mode: self.mode,
            sequence: self.sequence.clone(),
            structure_text: self.structure_text.clone(),
            chain: self.chain.clone(),
            centroid_moves: self.centroid_moves,
            repack_moves: self.repack_moves,
            design_moves: self.design_moves,
            seed: self.seed,
        }
    }
    fn restore(&mut self, s: StructPredictSnapshot) {
        self.mode = s.mode;
        self.sequence = s.sequence;
        self.structure_text = s.structure_text;
        self.chain = s.chain;
        self.centroid_moves = s.centroid_moves;
        self.repack_moves = s.repack_moves;
        self.design_moves = s.design_moves;
        self.seed = s.seed;
    }
    pub fn undo_edit(&mut self) -> bool {
        let current = self.snapshot();
        if let Some(prev) = self.history.undo(current) {
            self.restore(prev);
            self.error = None;
            true
        } else {
            false
        }
    }
    pub fn redo_edit(&mut self) -> bool {
        let current = self.snapshot();
        if let Some(next) = self.history.redo(current) {
            self.restore(next);
            self.error = None;
            true
        } else {
            false
        }
    }
    pub fn can_undo(&self) -> bool {
        self.history.can_undo()
    }
    pub fn can_redo(&self) -> bool {
        self.history.can_redo()
    }
}

impl Default for StructPredictPanel {
    fn default() -> Self {
        StructPredictPanel {
            mode: SpMode::AbInitio,
            // Trp-cage TC5b (20 aa) — the classic ab-initio mini-protein
            // benchmark; small enough to fold responsively.
            sequence: "NLYIQWLKDGGPSSGRPPPS".to_string(),
            structure_text: String::new(),
            chain: "A".to_string(),
            centroid_moves: 800,
            repack_moves: 150,
            design_moves: 500,
            seed: 42,
            error: None,
            result: String::new(),
            history: crate::undo::History::new(),
        }
    }
}

/// Render the Structure Prediction panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.structpredict;

    common::section(ui, "Method");
    ui.horizontal_wrapped(|ui| {
        ui.radio_value(&mut p.mode, SpMode::AbInitio, "Ab-initio fold")
            .on_hover_text("Monte-Carlo fragment assembly from sequence — no template.");
        ui.radio_value(&mut p.mode, SpMode::Design, "Fixed-backbone design")
            .on_hover_text("Redesign a sequence onto a supplied backbone (inverse folding).");
    });
    ui.label(
        egui::RichText::new(
            "Classical, weights-free methods (Modeller / Rosetta-class). Real and useful, \
             but NOT AlphaFold-accuracy — use the AlphaFold / ESMFold adapters for that.",
        )
        .weak()
        .small(),
    );
    ui.separator();

    match p.mode {
        SpMode::AbInitio => {
            common::section(ui, "Target sequence (one-letter amino acids)");
            common::seq_input(ui, "sp_seq", "sequence", &mut p.sequence, 4);
            ui.horizontal(|ui| {
                ui.label("centroid moves");
                ui.add(egui::DragValue::new(&mut p.centroid_moves).speed(10.0));
                ui.label("repack moves");
                ui.add(egui::DragValue::new(&mut p.repack_moves).speed(5.0));
            });
            ui.horizontal(|ui| {
                ui.label("seed");
                ui.add(egui::DragValue::new(&mut p.seed).speed(1.0));
            });
            ui.label(
                egui::RichText::new("Folding runs synchronously — a larger protein or more moves takes longer.")
                    .weak()
                    .small(),
            );
            if common::run_button(ui, "Predict (ab-initio)") {
                let snap = p.snapshot();
                p.history.record(snap);
                run_abinitio(p);
            }
        }
        SpMode::Design => {
            common::section(ui, "Backbone structure (PDB / mmCIF)");
            ui.horizontal(|ui| {
                if ui.small_button("Load structure…").clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("structure", &["pdb", "cif", "mmcif", "ent", "txt"])
                        .pick_file()
                    {
                        match valenx_core::io_caps::read_capped_to_string(
                            &path,
                            valenx_core::io_caps::MAX_GENETICS_FILE_BYTES as usize,
                        ) {
                            Ok(t) => p.structure_text = t,
                            Err(e) => p.error = Some(format!("read: {e}")),
                        }
                    }
                }
            });
            ui.add(
                egui::TextEdit::multiline(&mut p.structure_text)
                    .id_source("sp_structure")
                    .font(egui::TextStyle::Monospace)
                    .desired_width(f32::INFINITY)
                    .desired_rows(6)
                    .hint_text("paste a PDB / mmCIF backbone here"),
            );
            ui.horizontal(|ui| {
                ui.label("chain");
                ui.add(egui::TextEdit::singleline(&mut p.chain).desired_width(40.0));
                ui.label("design moves");
                ui.add(egui::DragValue::new(&mut p.design_moves).speed(10.0));
                ui.label("seed");
                ui.add(egui::DragValue::new(&mut p.seed).speed(1.0));
            });
            if common::run_button(ui, "Design sequence") {
                let snap = p.snapshot();
                p.history.record(snap);
                run_design(p);
            }
        }
    }

    ui.horizontal(|ui| {
        let (u, r) = common::undo_redo_inline(ui, p.can_undo(), p.can_redo());
        if u {
            p.undo_edit();
        }
        if r {
            p.redo_edit();
        }
        ui.label(
            egui::RichText::new("Ctrl+Z / Ctrl+Y reverses last Run")
                .weak()
                .small(),
        );
    });

    common::error_line(ui, &p.error);
    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Result");
        common::mono_output(ui, "sp_result", &p.result, 16);
    }
}

/// Run ab-initio prediction — extracted from the button closure so it is
/// callable from the headless UI tests.
fn run_abinitio(p: &mut StructPredictPanel) {
    p.error = None;
    let seq = common::clean_sequence(&p.sequence);
    if seq.is_empty() {
        p.error = Some("enter a protein sequence (one-letter amino acids)".into());
        return;
    }
    match predict_abinitio(&seq, p.centroid_moves, p.repack_moves, p.seed) {
        Ok(report) => p.result = format_report(&report),
        Err(e) => p.error = Some(format!("predict: {e}")),
    }
}

/// Run fixed-backbone design — extracted for the headless UI tests.
fn run_design(p: &mut StructPredictPanel) {
    p.error = None;
    if p.structure_text.trim().is_empty() {
        p.error = Some("paste or load a backbone structure (PDB / mmCIF)".into());
        return;
    }
    let chain = if p.chain.trim().is_empty() {
        "A"
    } else {
        p.chain.trim()
    };
    match design_sequence(&p.structure_text, chain, p.design_moves, p.seed) {
        Ok(report) => p.result = format_report(&report),
        Err(e) => p.error = Some(format!("design: {e}")),
    }
}

/// Format a [`StructPredictReport`] for the result block. The report's
/// own `notes` field carries the honest provenance + caveat summary.
fn format_report(report: &StructPredictReport) -> String {
    let seq = report.sequence.clone().unwrap_or_default();
    let mut out = format!("job      : {:?}\nresidues : {}\n", report.kind, seq.len());
    if !seq.is_empty() {
        out.push_str(&format!("sequence : {seq}\n"));
    }
    out.push('\n');
    out.push_str(&report.notes);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn abinitio_folds_the_default_peptide() {
        // Low move counts keep the headless test fast; the 20-aa default
        // (Trp-cage) is a valid ab-initio target.
        let mut p = StructPredictPanel {
            centroid_moves: 50,
            repack_moves: 10,
            ..Default::default()
        };
        run_abinitio(&mut p);
        assert!(p.error.is_none(), "unexpected error: {:?}", p.error);
        assert!(p.result.contains("residues"));
    }

    #[test]
    fn empty_sequence_fails_loud() {
        let mut p = StructPredictPanel {
            sequence: String::new(),
            ..Default::default()
        };
        run_abinitio(&mut p);
        assert!(p.error.is_some(), "empty sequence must surface an error");
    }

    #[test]
    fn design_without_structure_fails_loud() {
        let mut p = StructPredictPanel {
            mode: SpMode::Design,
            structure_text: String::new(),
            ..Default::default()
        };
        run_design(&mut p);
        assert!(p.error.is_some(), "design with no backbone must surface an error");
    }
}
