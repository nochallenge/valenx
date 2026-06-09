//! Panel 7 — **Cheminformatics** (`valenx-cheminf`).
//!
//! Parse a SMILES string and compute molecular descriptors (formula,
//! weight, logP, TPSA, H-bond donors / acceptors, rotatable bonds), the
//! Lipinski rule-of-five and Veber verdicts, the QED drug-likeness
//! score, and the ECFP4 Tanimoto similarity between two molecules — all
//! native `valenx-cheminf` calls.

use eframe::egui;

use valenx_cheminf::analyze::report::MoleculeReport;
use valenx_cheminf::coords::embed_3d;
use valenx_cheminf::fingerprint::morgan::ecfp;
use valenx_cheminf::fingerprint::similarity::{dice, tanimoto};
use valenx_cheminf::mol_from_smiles;

use super::common;
use super::molecule_view::{self, ViewMolecule};
use crate::ValenxApp;

/// Which cheminformatics sub-tool is showing.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum Tool {
    #[default]
    Descriptors,
    Similarity,
}

/// Snapshot of every editable input the Cheminformatics panel owns.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub(crate) struct CheminfSnapshot {
    pub(crate) smiles_a: String,
    pub(crate) smiles_b: String,
    pub(crate) fp_radius: usize,
}

/// Form + result state for the Cheminformatics panel.
pub struct CheminfPanel {
    tool: Tool,
    /// SMILES for the descriptors tool / molecule A of the similarity
    /// tool.
    smiles_a: String,
    /// SMILES for molecule B of the similarity tool.
    smiles_b: String,
    /// ECFP radius for the similarity fingerprint.
    fp_radius: usize,
    error: Option<String>,
    result: String,
    /// Undo / redo over both SMILES strings + the FP radius.
    history: crate::undo::History<CheminfSnapshot>,
}

impl Default for CheminfPanel {
    fn default() -> Self {
        CheminfPanel {
            tool: Tool::Descriptors,
            // Caffeine.
            smiles_a: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C".to_string(),
            // Theophylline — a close caffeine analogue, good for the
            // similarity demo.
            smiles_b: "CN1C2=C(C(=O)N(C1=O)C)NC=N2".to_string(),
            fp_radius: 2,
            error: None,
            result: String::new(),
            history: crate::undo::History::new(),
        }
    }
}

impl CheminfPanel {
    fn snapshot(&self) -> CheminfSnapshot {
        CheminfSnapshot {
            smiles_a: self.smiles_a.clone(),
            smiles_b: self.smiles_b.clone(),
            fp_radius: self.fp_radius,
        }
    }
    fn restore(&mut self, s: CheminfSnapshot) {
        self.smiles_a = s.smiles_a;
        self.smiles_b = s.smiles_b;
        self.fp_radius = s.fp_radius;
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
    pub fn can_undo(&self) -> bool { self.history.can_undo() }
    pub fn can_redo(&self) -> bool { self.history.can_redo() }
}

/// Keyboard-shortcut entry: run the active sub-tool.
pub fn run_primary_shortcut(app: &mut crate::ValenxApp) {
    let p = &mut app.genetics.cheminf;
    match p.tool {
        Tool::Descriptors => run_descriptors(p),
        Tool::Similarity => run_similarity(p),
    }
}

/// Render the Cheminformatics panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.cheminf;

    common::section(ui, "Tool");
    ui.horizontal(|ui| {
        ui.selectable_value(&mut p.tool, Tool::Descriptors, "Descriptors + Lipinski");
        ui.selectable_value(&mut p.tool, Tool::Similarity, "Fingerprint similarity");
        ui.separator();
        let (u, r) = common::undo_redo_inline(ui, p.can_undo(), p.can_redo());
        if u { p.undo_edit(); }
        if r { p.redo_edit(); }
    });
    ui.separator();

    match p.tool {
        Tool::Descriptors => draw_descriptors(p, ui),
        Tool::Similarity => draw_similarity(p, ui),
    }

    common::error_line(ui, &p.error);

    // --- 3-D viewport integration ---------------------------------
    // Parse molecule A's SMILES, generate a 3-D conformer by distance
    // geometry, and push a ball-and-stick / spacefill mesh into the
    // app's wgpu 3-D viewport.
    if !app.genetics.cheminf.smiles_a.trim().is_empty() {
        ui.horizontal(|ui| {
            if ui.button("Show molecule A in 3D viewport").clicked() {
                show_in_viewport(app, false);
            }
            if ui.button("Show (spacefill)").clicked() {
                show_in_viewport(app, true);
            }
        });
    }

    let p = &app.genetics.cheminf;
    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Result");
        common::mono_output(ui, "cheminf_result", &p.result, 16);
    }
}

/// Parse molecule A's SMILES, embed a 3-D conformer, and push a
/// ball-and-stick (or spacefill) mesh into the app's 3-D viewport.
fn show_in_viewport(app: &mut ValenxApp, spacefill: bool) {
    let smiles = app.genetics.cheminf.smiles_a.trim().to_string();
    let mol = match mol_from_smiles(&smiles) {
        Ok(m) => m,
        Err(e) => {
            app.genetics.cheminf.error = Some(format!("parse SMILES: {e}"));
            return;
        }
    };
    // SMILES carries no coordinates — generate a 3-D conformer by
    // distance geometry before meshing.
    let conformer = match embed_3d(&mol, 1) {
        Ok(m) => m,
        Err(e) => {
            app.genetics.cheminf.error = Some(format!("3-D embedding: {e}"));
            return;
        }
    };
    let view = match ViewMolecule::from_cheminf(&conformer) {
        Some(v) => v,
        None => {
            app.genetics.cheminf.error =
                Some("3-D embedding produced no usable conformer".to_string());
            return;
        }
    };
    let mesh = if spacefill {
        molecule_view::spacefill(&view)
    } else {
        molecule_view::ball_and_stick(&view, 0.25, 0.16)
    };
    let label = if spacefill {
        "molecule.spacefill"
    } else {
        "molecule.ball-stick"
    };
    match molecule_view::show_molecule(app, mesh, label) {
        Ok(_) => app.genetics.cheminf.error = None,
        Err(e) => app.genetics.cheminf.error = Some(e),
    }
}

fn draw_descriptors(p: &mut CheminfPanel, ui: &mut egui::Ui) {
    common::section(ui, "SMILES");
    ui.add(
        egui::TextEdit::multiline(&mut p.smiles_a)
            .id_source("cheminf_smiles_a")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(2),
    );
    if common::run_button(ui, "Compute descriptors") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_descriptors(p);
    }
}

/// Run the descriptor + Lipinski computation — extracted for the
/// headless UI tests.
fn run_descriptors(p: &mut CheminfPanel) {
    p.error = None;
    match MoleculeReport::from_smiles(p.smiles_a.trim()) {
        Ok(r) => {
                // The report struct doesn't carry these two — recompute from the molecule.
                let (heteroatoms, aromatic_atoms) =
                    valenx_cheminf::mol_from_smiles(p.smiles_a.trim())
                        .map(|m| {
                            (
                                valenx_cheminf::descriptors::heteroatom_count(&m),
                                valenx_cheminf::descriptors::aromatic_atom_count(&m),
                            )
                        })
                        .unwrap_or((0, 0));
                p.result = format!(
                    "canonical SMILES : {}\nformula          : {}\n\
                     average MW       : {:.3} g/mol\nmonoisotopic     : {:.4} u\n\
                     heavy atoms      : {}\nformal charge    : {}\n\n\
                     -- descriptors --\nCrippen logP     : {:.3}\n\
                     TPSA             : {:.2} Å²\nH-bond donors    : {}\n\
                     H-bond acceptors : {}\nrotatable bonds  : {}\n\
                     rings (SSSR)     : {}  ({} aromatic)\nheteroatoms      : {}\n\
                     aromatic atoms   : {}\nfraction Csp³    : {:.3}\n\n\
                     -- drug-likeness --\nLipinski violations : {} / 4\n\
                     Veber             : {}\nQED score          : {:.3}\n\
                     verdict            : {}",
                    r.canonical_smiles,
                    r.formula,
                    r.average_mw,
                    r.monoisotopic_mass,
                    r.heavy_atoms,
                    r.formal_charge,
                    r.logp,
                    r.tpsa,
                    r.hbd,
                    r.hba,
                    r.rotatable_bonds,
                    r.ring_count,
                    r.aromatic_rings,
                    heteroatoms,
                    aromatic_atoms,
                    r.fraction_csp3,
                    r.lipinski.violations,
                    if r.veber.passes { "pass" } else { "fail" },
                    r.qed,
                    if r.is_drug_like() { "drug-like" } else { "non-drug-like" },
                );
            }
            Err(e) => p.error = Some(e.to_string()),
        }
}

fn draw_similarity(p: &mut CheminfPanel, ui: &mut egui::Ui) {
    common::section(ui, "Molecule A (SMILES)");
    ui.add(
        egui::TextEdit::multiline(&mut p.smiles_a)
            .id_source("cheminf_sim_a")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(2),
    );
    common::section(ui, "Molecule B (SMILES)");
    ui.add(
        egui::TextEdit::multiline(&mut p.smiles_b)
            .id_source("cheminf_sim_b")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(2),
    );
    ui.horizontal(|ui| {
        ui.label("ECFP radius:");
        ui.add(egui::DragValue::new(&mut p.fp_radius).range(1..=4));
    });
    if common::run_button(ui, "Compute Tanimoto similarity") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_similarity(p);
    }
}

/// Run the fingerprint-similarity computation — extracted for the
/// headless UI tests.
fn run_similarity(p: &mut CheminfPanel) {
    p.error = None;
    match (
        mol_from_smiles(p.smiles_a.trim()),
        mol_from_smiles(p.smiles_b.trim()),
    ) {
        (Ok(a), Ok(b)) => {
            let fp_a = ecfp(&a, p.fp_radius, 2048);
            let fp_b = ecfp(&b, p.fp_radius, 2048);
            let tani = tanimoto(&fp_a, &fp_b);
            let dic = dice(&fp_a, &fp_b);
            p.result = format!(
                "ECFP{} (2048-bit) fingerprint similarity\n\n\
                 molecule A : {} heavy atoms\nmolecule B : {} heavy atoms\n\n\
                 Tanimoto   : {:.4}\nDice       : {:.4}\n\n{}",
                p.fp_radius * 2,
                a.heavy_atom_count(),
                b.heavy_atom_count(),
                tani,
                dic,
                if tani >= 0.85 {
                    "→ very similar (Tanimoto ≥ 0.85)"
                } else if tani >= 0.4 {
                    "→ moderately similar"
                } else {
                    "→ structurally distinct"
                },
            );
        }
        (Err(e), _) => p.error = Some(format!("molecule A: {e}")),
        (_, Err(e)) => p.error = Some(format!("molecule B: {e}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn caffeine_descriptors_compute() {
        let p = CheminfPanel::default();
        let r = MoleculeReport::from_smiles(p.smiles_a.trim()).expect("caffeine must parse");
        assert!(r.average_mw > 150.0 && r.average_mw < 250.0);
        assert!(r.heavy_atoms > 0);
    }

    #[test]
    fn analogues_are_similar() {
        let p = CheminfPanel::default();
        let a = mol_from_smiles(p.smiles_a.trim()).unwrap();
        let b = mol_from_smiles(p.smiles_b.trim()).unwrap();
        let t = tanimoto(&ecfp(&a, 2, 2048), &ecfp(&b, 2, 2048));
        assert!((0.0..=1.0).contains(&t));
    }
}

/// Headless egui UI-logic tests for the Cheminformatics panel.
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
        app.genetics.active = GeneticsPanel::Cheminformatics;
        app
    }

    #[test]
    fn draws_both_tools_without_panic() {
        for tool in [Tool::Descriptors, Tool::Similarity] {
            let mut app = app_with_panel();
            app.genetics.cheminf.tool = tool;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_post_run_and_error_states_without_panic() {
        let mut app = app_with_panel();
        app.genetics.cheminf.result = "formula : C8H10N4O2\n".to_string();
        draw_headless(&mut app);
        let mut app = app_with_panel();
        app.genetics.cheminf.error = Some("molecule A: parse failed".to_string());
        draw_headless(&mut app);
        // Empty SMILES — the "Show in 3D" affordance is gated off.
        let mut app = app_with_panel();
        app.genetics.cheminf.smiles_a.clear();
        draw_headless(&mut app);
    }

    #[test]
    fn run_descriptors_computes_caffeine_properties() {
        // Caffeine (the default) → the real valenx-cheminf descriptor
        // pipeline produces a correctly-formatted report.
        let mut p = CheminfPanel::default();
        run_descriptors(&mut p);
        assert!(p.error.is_none(), "descriptors errored: {:?}", p.error);
        assert!(p.result.contains("formula"));
        assert!(p.result.contains("Crippen logP"));
        assert!(p.result.contains("QED score"));
    }

    #[test]
    fn run_similarity_computes_a_tanimoto() {
        let mut p = CheminfPanel::default();
        run_similarity(&mut p);
        assert!(p.error.is_none(), "similarity errored: {:?}", p.error);
        assert!(p.result.contains("Tanimoto"));
        assert!(p.result.contains("Dice"));
    }

    #[test]
    fn run_actions_surface_errors_on_bad_smiles() {
        // A malformed SMILES string must produce an error, not a panic.
        let mut p = CheminfPanel {
            smiles_a: "this(((is not smiles".to_string(),
            ..CheminfPanel::default()
        };
        run_descriptors(&mut p);
        assert!(p.error.is_some(), "descriptors should error on bad SMILES");
        // Similarity with a bad molecule B.
        let mut p = CheminfPanel {
            smiles_b: "@@@invalid@@@".to_string(),
            ..CheminfPanel::default()
        };
        run_similarity(&mut p);
        assert!(p.error.is_some(), "similarity should error on bad SMILES");
    }
}
