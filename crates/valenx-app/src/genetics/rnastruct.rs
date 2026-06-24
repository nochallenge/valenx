//! Panel 5 — **RNA Structure** (`valenx-rnastruct`).
//!
//! Fold an RNA sequence — Zuker minimum-free-energy folding (dot-bracket
//! plus MFE) and Nussinov maximum-base-pairing — and compute the
//! McCaskill partition function with base-pair probabilities. Every
//! computation is a native `valenx-rnastruct` call.

use eframe::egui;

use valenx_rnastruct::ensemble::partition::partition_function;
use valenx_rnastruct::fold::nussinov::fold as nussinov_fold;
use valenx_rnastruct::fold::zuker::mfe;
use valenx_rnastruct::rna::RnaSeq;

use super::common;
use crate::ValenxApp;

/// Which folding method is showing.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum Tool {
    #[default]
    Zuker,
    Nussinov,
    Partition,
}

/// Form + result state for the RNA Structure panel.
pub struct RnaStructPanel {
    tool: Tool,
    /// RNA / DNA sequence (DNA is transcribed automatically).
    seq_text: String,
    /// Base-pair-probability significance threshold for the partition
    /// function tool.
    bpp_threshold: f64,
    error: Option<String>,
    result: String,
    /// Undo / redo history over the editable RNA sequence text.
    history: crate::undo::History<String>,
}

impl Default for RnaStructPanel {
    fn default() -> Self {
        RnaStructPanel {
            tool: Tool::Zuker,
            seq_text: "GGGAAAUCCUCUUUACCCGGAAGAGGGAAACCC".to_string(),
            bpp_threshold: 0.10,
            error: None,
            result: String::new(),
            history: crate::undo::History::new(),
        }
    }
}

impl RnaStructPanel {
    /// Undo the last sequence edit. Returns `true` if a snapshot was
    /// popped.
    pub fn undo_edit(&mut self) -> bool {
        if let Some(prev) = self.history.undo(self.seq_text.clone()) {
            self.seq_text = prev;
            self.error = None;
            true
        } else {
            false
        }
    }
    /// Redo the last undone sequence edit.
    pub fn redo_edit(&mut self) -> bool {
        if let Some(next) = self.history.redo(self.seq_text.clone()) {
            self.seq_text = next;
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

/// Keyboard-shortcut entry: run the active folding tool.
pub fn run_primary_shortcut(app: &mut crate::ValenxApp) {
    run_fold(&mut app.genetics.rnastruct);
}

/// Render the RNA Structure panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.rnastruct;

    common::section(ui, "Sequence");
    common::seq_input(
        ui,
        "rna_seq_input",
        "RNA (or DNA — auto-transcribed):",
        &mut p.seq_text,
        4,
    );
    ui.label(
        egui::RichText::new(format!(
            "length: {}",
            common::clean_sequence(&p.seq_text).len()
        ))
        .weak()
        .small(),
    );
    ui.separator();

    common::section(ui, "Folding method");
    ui.horizontal_wrapped(|ui| {
        ui.selectable_value(&mut p.tool, Tool::Zuker, "Zuker MFE");
        ui.selectable_value(&mut p.tool, Tool::Nussinov, "Nussinov");
        ui.selectable_value(&mut p.tool, Tool::Partition, "Partition function");
    });
    ui.add_space(4.0);

    if p.tool == Tool::Partition {
        ui.horizontal(|ui| {
            let lbl = ui.label("BPP threshold:");
            ui.add(
                egui::DragValue::new(&mut p.bpp_threshold)
                    .speed(0.01)
                    .range(0.0..=1.0),
            )
            .labelled_by(lbl.id);
        });
    }

    let run_label = match p.tool {
        Tool::Zuker => "Fold (Zuker MFE)",
        Tool::Nussinov => "Fold (Nussinov)",
        Tool::Partition => "Compute partition function",
    };
    if common::run_button(ui, run_label) {
        run_fold(p);
    }

    common::error_line(ui, &p.error);
    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Result");
        common::mono_output(ui, "rna_result", &p.result, 14);
    }
}

/// Run the selected RNA folding tool — extracted from the button
/// closure so it is callable from the headless UI tests.
fn run_fold(p: &mut RnaStructPanel) {
    p.error = None;
    let cleaned = common::clean_sequence(&p.seq_text);
    if cleaned.is_empty() {
        p.error = Some("sequence is empty".into());
        return;
    }
    match RnaSeq::parse(cleaned.as_bytes()) {
        Ok(rna) => match p.tool {
            Tool::Zuker => match mfe(&rna) {
                Ok(res) => {
                    let db = res.structure.to_dot_bracket();
                    p.result = format!(
                        "Zuker minimum-free-energy fold\n\
                             MFE          : {:.2} kcal/mol\n\
                             base pairs   : {}\n\n\
                             {}\n{}",
                        res.energy,
                        res.structure.pairs().len(),
                        rna.as_str(),
                        db,
                    );
                }
                Err(e) => p.error = Some(e.to_string()),
            },
            Tool::Nussinov => match nussinov_fold(&rna) {
                Ok(res) => {
                    let db = res.structure.to_dot_bracket();
                    p.result = format!(
                        "Nussinov maximum-base-pairing fold\n\
                             base pairs   : {}\n\n{}\n{}",
                        res.structure.pairs().len(),
                        rna.as_str(),
                        db,
                    );
                }
                Err(e) => p.error = Some(e.to_string()),
            },
            Tool::Partition => match partition_function(&rna) {
                Ok(pf) => {
                    let mut sig = pf.significant_pairs(p.bpp_threshold);
                    sig.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
                    let mut out = format!(
                        "McCaskill partition function\n\
                             ensemble ΔG  : {:.3} kcal/mol\n\
                             {} base pairs with p ≥ {:.2}:\n",
                        pf.ensemble_free_energy(),
                        sig.len(),
                        p.bpp_threshold,
                    );
                    for (i, j, prob) in sig.iter().take(40) {
                        let bar = "#".repeat((prob * 20.0).round() as usize);
                        out.push_str(&format!(
                            "  {:>3}-{:<3}  p={:.3}  {}\n",
                            i + 1,
                            j + 1,
                            prob,
                            bar,
                        ));
                    }
                    p.result = out;
                }
                Err(e) => p.error = Some(e.to_string()),
            },
        },
        Err(e) => p.error = Some(e.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_sequence_folds_with_zuker() {
        let p = RnaStructPanel::default();
        let cleaned = common::clean_sequence(&p.seq_text);
        let rna = RnaSeq::parse(cleaned.as_bytes()).unwrap();
        let res = mfe(&rna).unwrap();
        // The dot-bracket string is always the sequence's length.
        assert_eq!(res.structure.to_dot_bracket().len(), rna.len());
    }

    #[test]
    fn partition_function_runs() {
        let p = RnaStructPanel::default();
        let rna = RnaSeq::parse(common::clean_sequence(&p.seq_text).as_bytes()).unwrap();
        let pf = partition_function(&rna).unwrap();
        assert!(pf.q() > 0.0);
    }
}

/// Headless egui UI-logic tests for the RNA Structure panel.
#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use crate::genetics_workbench::GeneticsPanel;
    use crate::ValenxApp;

    fn draw_headless(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                super::draw(app, ui);
            });
        });
    }

    fn app_with_panel() -> ValenxApp {
        let mut app = ValenxApp::default();
        app.genetics.active = GeneticsPanel::RnaStructure;
        app
    }

    #[test]
    fn draws_every_tool_without_panic() {
        for tool in [Tool::Zuker, Tool::Nussinov, Tool::Partition] {
            let mut app = app_with_panel();
            app.genetics.rnastruct.tool = tool;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_post_run_and_error_states_without_panic() {
        let mut app = app_with_panel();
        app.genetics.rnastruct.result = "Zuker fold\nMFE : -8.20 kcal/mol\n".to_string();
        draw_headless(&mut app);
        let mut app = app_with_panel();
        app.genetics.rnastruct.error = Some("sequence is empty".to_string());
        draw_headless(&mut app);
        let mut app = app_with_panel();
        app.genetics.rnastruct.seq_text.clear();
        draw_headless(&mut app);
    }

    #[test]
    fn run_fold_with_each_method_produces_a_structure() {
        // Each folding method calls the real valenx-rnastruct API.
        for tool in [Tool::Zuker, Tool::Nussinov, Tool::Partition] {
            let mut p = RnaStructPanel {
                tool,
                ..RnaStructPanel::default()
            };
            run_fold(&mut p);
            assert!(p.error.is_none(), "{tool:?} errored: {:?}", p.error);
            assert!(!p.result.is_empty());
        }
    }

    #[test]
    fn run_fold_zuker_reports_an_mfe() {
        let mut p = RnaStructPanel {
            tool: Tool::Zuker,
            ..RnaStructPanel::default()
        };
        run_fold(&mut p);
        assert!(p.result.contains("MFE"));
        assert!(p.result.contains("kcal/mol"));
    }

    #[test]
    fn run_fold_surfaces_error_on_empty_input() {
        let mut p = RnaStructPanel {
            seq_text: String::new(),
            ..RnaStructPanel::default()
        };
        run_fold(&mut p);
        assert!(p.error.is_some(), "fold should error on empty input");
        assert!(p.result.is_empty());
    }

    #[test]
    fn run_fold_handles_malformed_input() {
        // A sequence of non-nucleotide letters is malformed RNA — the
        // panel must surface an error rather than panicking.
        let mut p = RnaStructPanel {
            seq_text: "ZZZXXXQQQ".to_string(),
            ..RnaStructPanel::default()
        };
        run_fold(&mut p);
        assert!(p.error.is_some(), "fold should error on malformed input");
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        use egui::accesskit::Role;
        // The BPP threshold DragValue is only visible in Partition mode.
        let mut app = app_with_panel();
        app.genetics.rnastruct.tool = Tool::Partition;
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                super::draw(&mut app, ui);
            });
        });
        let nodes = out
            .platform_output
            .accesskit_update
            .expect("accesskit tree produced")
            .nodes;
        let spin_buttons: Vec<_> = nodes
            .iter()
            .filter(|(_, n)| n.role() == Role::SpinButton)
            .collect();
        assert!(
            !spin_buttons.is_empty(),
            "rnastruct Partition mode should expose at least one SpinButton"
        );
        assert!(
            spin_buttons
                .iter()
                .all(|(_, n)| !n.labelled_by().is_empty()),
            "every rnastruct DragValue must be labelled_by its caption (AI-drivable name)"
        );
    }
}
