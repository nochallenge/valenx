//! Panel 3 — **Phylogenetics** (`valenx-phylo`).
//!
//! Build a phylogenetic tree (neighbor-joining / UPGMA / WPGMA / BIONJ)
//! from aligned sequences via a substitution-corrected distance matrix,
//! read / write Newick, and render the tree as an ASCII cladogram —
//! all native `valenx-phylo` (+ `valenx-align` for the MSA) calls.

use eframe::egui;

use valenx_align::matrix::{GapCost, ScoringScheme, SubstitutionMatrix};
use valenx_align::msa::progressive::align as msa_align;
use valenx_align::msa::Msa;
use valenx_phylo::distance::cluster::{bionj, neighbor_joining, upgma, wpgma};
use valenx_phylo::distance::matrix::{distance_matrix, DistanceModel};
use valenx_phylo::io::newick::{read_newick, write_newick};
use valenx_phylo::render::render_ascii;
use valenx_phylo::tree::Tree;

use super::common;
use crate::ValenxApp;

/// Tree-building method.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum Method {
    #[default]
    NeighborJoining,
    Upgma,
    Wpgma,
    Bionj,
}

/// Panel sub-tool: build from sequences, or import a Newick string.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum Tool {
    #[default]
    Build,
    NewickIo,
}

/// Snapshot of every editable input the Phylogenetics panel owns.
/// One snapshot per `History` entry — `Ctrl+Z` rewinds every input
/// to its prior state atomically, mirroring the Alignment panel's
/// snapshot pattern.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct PhylogeneticsSnapshot {
    pub seq_text: String,
    pub newick: String,
}

/// Form + result state for the Phylogenetics panel.
pub struct PhylogeneticsPanel {
    tool: Tool,
    method: Method,
    model: DistanceModel,
    /// Aligned / unaligned sequences (multi-FASTA).
    seq_text: String,
    /// Newick text for the import tool / the build output.
    newick: String,
    error: Option<String>,
    result: String,
    /// Undo / redo history over both editable text inputs. See
    /// [`PhylogeneticsSnapshot`].
    history: crate::undo::History<PhylogeneticsSnapshot>,
}

impl Default for PhylogeneticsPanel {
    fn default() -> Self {
        PhylogeneticsPanel {
            tool: Tool::Build,
            method: Method::NeighborJoining,
            model: DistanceModel::JukesCantor,
            seq_text: ">human\nACGTACGTAAGGTTCCAACGTACGT\n\
                       >chimp\nACGTACGTAAGCTTCCAACGTACGT\n\
                       >mouse\nACGTACGAAAGGTTCGAACGTACGT\n\
                       >chicken\nACGTTCGTAAGGATCCAACGAACGT\n"
                .to_string(),
            newick: "((A:0.1,B:0.2):0.05,(C:0.15,D:0.1):0.05);".to_string(),
            error: None,
            result: String::new(),
            history: crate::undo::History::new(),
        }
    }
}

impl PhylogeneticsPanel {
    fn snapshot(&self) -> PhylogeneticsSnapshot {
        PhylogeneticsSnapshot {
            seq_text: self.seq_text.clone(),
            newick: self.newick.clone(),
        }
    }
    fn restore(&mut self, s: PhylogeneticsSnapshot) {
        self.seq_text = s.seq_text;
        self.newick = s.newick;
    }
    /// Undo the most recent input edit.
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
    /// Redo the most recently undone input edit.
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

/// Render the Phylogenetics panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.phylogenetics;

    common::section(ui, "Tool");
    ui.horizontal(|ui| {
        ui.selectable_value(&mut p.tool, Tool::Build, "Build tree")
            .on_hover_text(
                "Build a phylogenetic tree from a multi-FASTA via NJ / BIONJ / UPGMA / WPGMA.",
            );
        ui.selectable_value(&mut p.tool, Tool::NewickIo, "Newick in/out")
            .on_hover_text("Round-trip a Newick string + render the tree as ASCII / SVG.");
        // Inline undo / redo so Ctrl+Z reverses the last input edit.
        ui.separator();
        let (undo, redo) = common::undo_redo_inline(ui, p.can_undo(), p.can_redo());
        if undo {
            p.undo_edit();
        }
        if redo {
            p.redo_edit();
        }
    });
    ui.separator();

    match p.tool {
        Tool::Build => draw_build(p, ui),
        Tool::NewickIo => draw_newick(p, ui),
    }

    common::error_line(ui, &p.error);
    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Result");
        common::mono_output(ui, "phylo_result", &p.result, 16);
    }
}

fn draw_build(p: &mut PhylogeneticsPanel, ui: &mut egui::Ui) {
    common::section(ui, "Sequences (multi-FASTA)");
    ui.add(
        egui::TextEdit::multiline(&mut p.seq_text)
            .id_source("phylo_seq_input")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(7),
    );

    common::section(ui, "Method + model");
    ui.horizontal_wrapped(|ui| {
        ui.radio_value(&mut p.method, Method::NeighborJoining, "NJ")
            .on_hover_text("Neighbor-Joining (Saitou & Nei, 1987). Unrooted, distance-based.");
        ui.radio_value(&mut p.method, Method::Bionj, "BIONJ")
            .on_hover_text("BIONJ (Gascuel, 1997) — variance-weighted NJ improvement.");
        ui.radio_value(&mut p.method, Method::Upgma, "UPGMA")
            .on_hover_text("UPGMA — rooted ultrametric tree assuming a constant molecular clock.");
        ui.radio_value(&mut p.method, Method::Wpgma, "WPGMA")
            .on_hover_text("WPGMA — UPGMA with equal weighting at each merge step.");
    });
    ui.horizontal(|ui| {
        ui.label("Distance model:")
            .on_hover_text("How pairwise distances are computed from the alignment.");
        egui::ComboBox::from_id_source("phylo_model")
            .selected_text(format!("{:?}", p.model))
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut p.model, DistanceModel::PDistance, "p-distance")
                    .on_hover_text("Fraction of differing positions — uncorrected, fast.");
                ui.selectable_value(&mut p.model, DistanceModel::JukesCantor, "Jukes-Cantor")
                    .on_hover_text(
                        "Jukes-Cantor 1969 — single substitution rate, equal base frequencies.",
                    );
                ui.selectable_value(&mut p.model, DistanceModel::Kimura2P, "Kimura-2P")
                    .on_hover_text("Kimura 1980 — distinct transition / transversion rates.");
                ui.selectable_value(&mut p.model, DistanceModel::TamuraNei, "Tamura-Nei")
                    .on_hover_text(
                        "Tamura-Nei 1993 — separate purine / pyrimidine transition rates.",
                    );
            });
    });

    if common::run_button(ui, "Build phylogeny") {
        // Snapshot the inputs the user just executed so a later Ctrl+Z
        // can recover them after subsequent edits.
        let snap = p.snapshot();
        p.history.record(snap);
        run_build(p);
    }
}

/// Run the phylogeny build — extracted for the headless UI tests.
fn run_build(p: &mut PhylogeneticsPanel) {
    p.error = None;
    let recs = common::parse_fasta(&p.seq_text);
    if recs.len() < 2 {
        p.error = Some("need at least two named sequences (>name FASTA)".into());
        return;
    }
    // The distance matrix needs an equal-width MSA. If the input
    // rows already share a length use them directly; otherwise run
    // a quick progressive alignment first.
    let labels: Vec<String> = recs
        .iter()
        .enumerate()
        .map(|(i, (l, _))| {
            if l.is_empty() {
                format!("seq{}", i + 1)
            } else {
                l.clone()
            }
        })
        .collect();
    let owned: Vec<Vec<u8>> = recs
        .iter()
        .map(|(_, s)| s.to_ascii_uppercase().into_bytes())
        .collect();
    let same_width = owned
        .first()
        .map(|f| owned.iter().all(|r| r.len() == f.len()))
        .unwrap_or(false);
    let msa: Msa = if same_width {
        match Msa::new(owned.clone()) {
            Ok(m) => m,
            Err(e) => {
                p.error = Some(e.to_string());
                return;
            }
        }
    } else {
        let refs: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
        let scheme = ScoringScheme::new(SubstitutionMatrix::nuc44(), GapCost::new(10, 1));
        match msa_align(&refs, &scheme) {
            Ok(m) => m,
            Err(e) => {
                p.error = Some(format!("alignment: {e}"));
                return;
            }
        }
    };
    let dm = match distance_matrix(&msa, &labels, p.model) {
        Ok(d) => d,
        Err(e) => {
            p.error = Some(format!("distance matrix: {e}"));
            return;
        }
    };
    let tree: Result<Tree, _> = match p.method {
        Method::NeighborJoining => neighbor_joining(&dm),
        Method::Bionj => bionj(&dm),
        Method::Upgma => upgma(&dm),
        Method::Wpgma => wpgma(&dm),
    };
    match tree {
        Ok(t) => {
            p.newick = write_newick(&t);
            p.result = format!(
                    "{:?} tree · {} leaves · {} nodes · {} cherries · depth {} · sackin {} · tree length {:.3} · binary {}\n\nCladogram:\n{}\n\nNewick:\n{}",
                    p.method,
                    t.leaf_count(),
                    t.node_count(),
                    t.cherry_count(),
                    t.max_depth(),
                    t.sackin_index(),
                    t.total_length(),
                    t.is_binary(),
                    render_ascii(&t, 56),
                    p.newick,
                );
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

fn draw_newick(p: &mut PhylogeneticsPanel, ui: &mut egui::Ui) {
    common::section(ui, "Newick string");
    ui.add(
        egui::TextEdit::multiline(&mut p.newick)
            .id_source("phylo_newick_input")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(4),
    );
    ui.horizontal(|ui| {
        if ui.small_button("Load Newick…").clicked() {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("Newick", &["nwk", "newick", "tree", "txt"])
                .pick_file()
            {
                // Round-21 H1: see biostruct loader.
                match valenx_core::io_caps::read_capped_to_string(
                    &path,
                    valenx_core::io_caps::MAX_GENETICS_FILE_BYTES as usize,
                ) {
                    Ok(t) => p.newick = t,
                    Err(e) => p.error = Some(format!("read: {e}")),
                }
            }
        }
        if ui.small_button("Save Newick…").clicked() {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("Newick", &["nwk"])
                .set_file_name("tree.nwk")
                .save_file()
            {
                if let Err(e) = valenx_core::io_caps::atomic_write_str(&path, &p.newick) {
                    p.error = Some(format!("write: {e}"));
                }
            }
        }
    });
    if common::run_button(ui, "Parse + render") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_newick(p);
    }
}

/// Run the Newick parse + render — extracted for the headless UI tests.
fn run_newick(p: &mut PhylogeneticsPanel) {
    p.error = None;
    match read_newick(&p.newick) {
        Ok(t) => {
            p.result = format!(
                "parsed · {} leaves · {} nodes · {} cherries · depth {} · sackin {} · tree length {:.3} · binary {}\nleaves: {}\n\nCladogram:\n{}",
                t.leaf_count(),
                t.node_count(),
                t.cherry_count(),
                t.max_depth(),
                t.sackin_index(),
                t.total_length(),
                t.is_binary(),
                t.leaf_labels().join(", "),
                render_ascii(&t, 56),
            );
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_newick_parses() {
        let p = PhylogeneticsPanel::default();
        let t = read_newick(&p.newick).expect("default Newick must parse");
        assert_eq!(t.leaf_count(), 4);
    }

    #[test]
    fn default_sequences_are_equal_width() {
        // The four demo sequences are pre-aligned (25 nt each) so the
        // build path can skip the progressive-alignment step.
        let p = PhylogeneticsPanel::default();
        let recs = common::parse_fasta(&p.seq_text);
        assert_eq!(recs.len(), 4);
        let w = recs[0].1.len();
        assert!(recs.iter().all(|(_, s)| s.len() == w));
    }
}

/// Headless egui UI-logic tests for the Phylogenetics panel.
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
        app.genetics.active = GeneticsPanel::Phylogenetics;
        app
    }

    #[test]
    fn draws_every_tool_without_panic() {
        for tool in [Tool::Build, Tool::NewickIo] {
            let mut app = app_with_panel();
            app.genetics.phylogenetics.tool = tool;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_post_run_and_error_states_without_panic() {
        let mut app = app_with_panel();
        app.genetics.phylogenetics.result = "NJ tree · 4 leaves\n".to_string();
        draw_headless(&mut app);
        let mut app = app_with_panel();
        app.genetics.phylogenetics.error = Some("need at least two named sequences".to_string());
        draw_headless(&mut app);
    }

    #[test]
    fn run_build_constructs_a_tree_with_each_method() {
        // Each tree-building method calls the real valenx-phylo API.
        for method in [
            Method::NeighborJoining,
            Method::Bionj,
            Method::Upgma,
            Method::Wpgma,
        ] {
            let mut p = PhylogeneticsPanel {
                method,
                ..PhylogeneticsPanel::default()
            };
            run_build(&mut p);
            assert!(p.error.is_none(), "{method:?} errored: {:?}", p.error);
            assert!(p.result.contains("leaves"));
            assert!(p.result.contains("Newick"));
        }
    }

    #[test]
    fn run_newick_parses_and_renders() {
        let mut p = PhylogeneticsPanel::default();
        run_newick(&mut p);
        assert!(p.error.is_none(), "Newick parse errored: {:?}", p.error);
        assert!(p.result.contains("parsed"));
        assert!(p.result.contains("Cladogram"));
    }

    #[test]
    fn run_actions_surface_errors_on_bad_input() {
        // Build with a single sequence (needs ≥ 2).
        let mut p = PhylogeneticsPanel {
            seq_text: ">only\nACGTACGT\n".to_string(),
            ..PhylogeneticsPanel::default()
        };
        run_build(&mut p);
        assert!(p.error.is_some(), "build should error with one sequence");
        // Newick parse with a malformed string.
        let mut p = PhylogeneticsPanel {
            newick: "((((not a tree".to_string(),
            ..PhylogeneticsPanel::default()
        };
        run_newick(&mut p);
        assert!(p.error.is_some(), "parse should error on malformed Newick");
    }
}
