//! Panel 2 — **Alignment** (`valenx-align`).
//!
//! Pairwise alignment (Needleman-Wunsch / Smith-Waterman / Gotoh
//! affine) with a rendered alignment block, progressive multiple
//! sequence alignment, and a k-mer seed search — all native
//! `valenx-align` calls.

use eframe::egui;

use valenx_align::matrix::{GapCost, ScoringScheme, SubstitutionMatrix};
use valenx_align::msa::progressive::align as msa_align;
use valenx_align::pairwise::global::{gotoh, needleman_wunsch};
use valenx_align::pairwise::local::smith_waterman;
use valenx_align::search::kmer::KmerIndex;

use super::common;
use crate::ValenxApp;

/// Which alignment sub-tool is showing.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum Tool {
    #[default]
    Pairwise,
    Multiple,
    KmerSearch,
}

/// Pairwise algorithm choice.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum PairAlgo {
    #[default]
    NeedlemanWunsch,
    SmithWaterman,
    GotohAffine,
}

/// Substitution-matrix choice.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum Matrix {
    #[default]
    Nuc44,
    Blosum62,
    Blosum80,
    Identity,
}

impl Matrix {
    fn build(self) -> SubstitutionMatrix {
        match self {
            Matrix::Nuc44 => SubstitutionMatrix::nuc44(),
            Matrix::Blosum62 => SubstitutionMatrix::blosum62(),
            Matrix::Blosum80 => SubstitutionMatrix::blosum80(),
            Matrix::Identity => SubstitutionMatrix::identity(2, -1),
        }
    }
}

/// Form + result state for the Alignment panel.
pub struct AlignmentPanel {
    tool: Tool,
    algo: PairAlgo,
    matrix: Matrix,
    gap_open: i32,
    gap_extend: i32,
    /// Pairwise sequence A / B.
    seq_a: String,
    seq_b: String,
    /// MSA input — multi-FASTA or one sequence per line.
    msa_text: String,
    /// k-mer search — index sequences + query.
    kmer_db: String,
    kmer_query: String,
    kmer_k: usize,
    error: Option<String>,
    result: String,
    /// Undo / redo history for the currently-active input field. The
    /// stack captures whichever input the user was editing most
    /// recently (sequence A / B for pairwise, msa_text for MSA, the
    /// pair (kmer_db, kmer_query) bundled as a single snapshot for
    /// the k-mer tool).
    history: crate::undo::History<AlignmentSnapshot>,
}

/// Snapshot of every editable input the alignment panel owns. The
/// undo / redo stack stores these so a Ctrl+Z reverses the most
/// recent mutation regardless of which sub-tool was active.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct AlignmentSnapshot {
    pub seq_a: String,
    pub seq_b: String,
    pub msa_text: String,
    pub kmer_db: String,
    pub kmer_query: String,
}

impl Default for AlignmentPanel {
    fn default() -> Self {
        AlignmentPanel {
            tool: Tool::Pairwise,
            algo: PairAlgo::NeedlemanWunsch,
            matrix: Matrix::Nuc44,
            gap_open: 10,
            gap_extend: 1,
            seq_a: "ACGTACGTACGTAAGGTTCCAACGT".to_string(),
            seq_b: "ACGTAAGTACGTAACGTTCCAACGT".to_string(),
            msa_text: ">a\nACGTACGTAAGGTT\n>b\nACGTAAGTAAGGTT\n>c\nACGTACGTAACGTT\n"
                .to_string(),
            kmer_db: ">ref\nACGTACGTACGTAAGGTTCCAACGTACGTACGT\n".to_string(),
            kmer_query: "AAGGTTCCAACGT".to_string(),
            kmer_k: 8,
            error: None,
            result: String::new(),
            history: crate::undo::History::new(),
        }
    }
}

impl AlignmentPanel {
    /// Capture the current state of every editable input. Used by
    /// undo_edit / redo_edit to swap state in/out atomically.
    fn snapshot(&self) -> AlignmentSnapshot {
        AlignmentSnapshot {
            seq_a: self.seq_a.clone(),
            seq_b: self.seq_b.clone(),
            msa_text: self.msa_text.clone(),
            kmer_db: self.kmer_db.clone(),
            kmer_query: self.kmer_query.clone(),
        }
    }

    /// Apply a snapshot — restores every input from the captured state.
    fn restore(&mut self, s: AlignmentSnapshot) {
        self.seq_a = s.seq_a;
        self.seq_b = s.seq_b;
        self.msa_text = s.msa_text;
        self.kmer_db = s.kmer_db;
        self.kmer_query = s.kmer_query;
    }

    /// Try to undo the most recent input edit. Returns `true` if a
    /// snapshot was restored.
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

    /// Try to redo the most recently undone input edit.
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

    /// Whether an undo / redo would do something.
    pub fn can_undo(&self) -> bool { self.history.can_undo() }
    pub fn can_redo(&self) -> bool { self.history.can_redo() }
}

/// Keyboard-shortcut entry: run the active sub-tool.
pub fn run_primary_shortcut(app: &mut crate::ValenxApp) {
    let p = &mut app.genetics.alignment;
    match p.tool {
        Tool::Pairwise => run_pairwise(p),
        Tool::Multiple => run_msa(p),
        Tool::KmerSearch => run_kmer(p),
    }
}

/// Render the Alignment panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.alignment;

    common::section(ui, "Tool");
    ui.horizontal_wrapped(|ui| {
        ui.selectable_value(&mut p.tool, Tool::Pairwise, "Pairwise");
        ui.selectable_value(&mut p.tool, Tool::Multiple, "Multiple (MSA)");
        ui.selectable_value(&mut p.tool, Tool::KmerSearch, "k-mer search");
    });
    ui.separator();

    match p.tool {
        Tool::Pairwise => draw_pairwise(p, ui),
        Tool::Multiple => draw_msa(p, ui),
        Tool::KmerSearch => draw_kmer(p, ui),
    }

    common::error_line(ui, &p.error);
    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Result");
        common::mono_output(ui, "align_result", &p.result, 14);
    }
}

fn scoring_scheme(p: &AlignmentPanel) -> ScoringScheme {
    ScoringScheme::new(p.matrix.build(), GapCost::new(p.gap_open, p.gap_extend))
}

fn draw_pairwise(p: &mut AlignmentPanel, ui: &mut egui::Ui) {
    common::section(ui, "Sequences");
    common::seq_input(ui, "align_a", "Sequence A:", &mut p.seq_a, 3);
    common::seq_input(ui, "align_b", "Sequence B:", &mut p.seq_b, 3);

    common::section(ui, "Algorithm + scoring");
    ui.horizontal(|ui| {
        ui.radio_value(&mut p.algo, PairAlgo::NeedlemanWunsch, "NW (global)");
        ui.radio_value(&mut p.algo, PairAlgo::SmithWaterman, "SW (local)");
        ui.radio_value(&mut p.algo, PairAlgo::GotohAffine, "Gotoh (affine)");
    });
    ui.horizontal(|ui| {
        ui.label("Matrix:");
        egui::ComboBox::from_id_source("align_matrix")
            .selected_text(format!("{:?}", p.matrix))
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut p.matrix, Matrix::Nuc44, "NUC.4.4 (DNA)");
                ui.selectable_value(&mut p.matrix, Matrix::Blosum62, "BLOSUM62");
                ui.selectable_value(&mut p.matrix, Matrix::Blosum80, "BLOSUM80");
                ui.selectable_value(&mut p.matrix, Matrix::Identity, "identity ±");
            });
    });
    ui.horizontal(|ui| {
        ui.label("Gap open:");
        ui.add(egui::DragValue::new(&mut p.gap_open).range(0..=50));
        ui.label("extend:");
        ui.add(egui::DragValue::new(&mut p.gap_extend).range(0..=50));
    });

    if common::run_button(ui, "Align") {
        run_pairwise(p);
    }
}

/// Run the pairwise alignment — extracted for the headless UI tests.
fn run_pairwise(p: &mut AlignmentPanel) {
    p.error = None;
    let a = common::clean_sequence(&p.seq_a);
    let b = common::clean_sequence(&p.seq_b);
    if a.is_empty() || b.is_empty() {
        p.error = Some("both sequences must be non-empty".into());
        return;
    }
    let scheme = scoring_scheme(p);
    let result = match p.algo {
        PairAlgo::NeedlemanWunsch => needleman_wunsch(a.as_bytes(), b.as_bytes(), &scheme),
        PairAlgo::SmithWaterman => smith_waterman(a.as_bytes(), b.as_bytes(), &scheme),
        PairAlgo::GotohAffine => gotoh(a.as_bytes(), b.as_bytes(), &scheme),
    };
    match result {
        Ok(aln) => {
            let stats = aln.stats(&scheme.matrix);
            p.result = format!(
                "score        : {}\ncolumns      : {}\nidentity     : {:.1} %\n\
                 similarity   : {:.1} %\ngaps         : {} ({} opens)\n\n{}",
                aln.score,
                stats.columns,
                stats.percent_identity(),
                stats.percent_similarity(),
                stats.gaps,
                stats.gap_opens,
                aln.pretty(60, &scheme.matrix),
            );
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

fn draw_msa(p: &mut AlignmentPanel, ui: &mut egui::Ui) {
    common::section(ui, "Sequences (multi-FASTA)");
    ui.add(
        egui::TextEdit::multiline(&mut p.msa_text)
            .id_source("align_msa_input")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(6),
    );
    ui.horizontal(|ui| {
        ui.label("Matrix:");
        egui::ComboBox::from_id_source("align_msa_matrix")
            .selected_text(format!("{:?}", p.matrix))
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut p.matrix, Matrix::Nuc44, "NUC.4.4 (DNA)");
                ui.selectable_value(&mut p.matrix, Matrix::Blosum62, "BLOSUM62");
                ui.selectable_value(&mut p.matrix, Matrix::Identity, "identity ±");
            });
    });
    if common::run_button(ui, "Build MSA (progressive)") {
        run_msa(p);
    }
}

/// Run the progressive MSA — extracted for the headless UI tests.
fn run_msa(p: &mut AlignmentPanel) {
    p.error = None;
    let recs = common::parse_fasta(&p.msa_text);
    if recs.len() < 2 {
        p.error = Some("need at least two sequences (use >name FASTA headers)".into());
        return;
    }
    let owned: Vec<Vec<u8>> = recs
        .iter()
        .map(|(_, s)| s.to_ascii_uppercase().into_bytes())
        .collect();
    let refs: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
    let scheme = scoring_scheme(p);
    match msa_align(&refs, &scheme) {
        Ok(msa) => {
            let mut out = format!(
                "{} sequences · {} columns · SP score {}\n\n",
                msa.depth(),
                msa.width(),
                msa.sum_of_pairs(&scheme),
            );
            for (i, (label, _)) in recs.iter().enumerate() {
                let name = if label.is_empty() {
                    format!("seq{}", i + 1)
                } else {
                    label.clone()
                };
                let row = msa.row_str(i).unwrap_or("");
                out.push_str(&format!("{name:<12} {row}\n"));
            }
            p.result = out;
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

fn draw_kmer(p: &mut AlignmentPanel, ui: &mut egui::Ui) {
    common::section(ui, "Reference (multi-FASTA)");
    ui.add(
        egui::TextEdit::multiline(&mut p.kmer_db)
            .id_source("align_kmer_db")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(4),
    );
    common::seq_input(ui, "align_kmer_query", "Query:", &mut p.kmer_query, 2);
    ui.horizontal(|ui| {
        ui.label("k:");
        ui.add(egui::DragValue::new(&mut p.kmer_k).range(2..=31));
    });
    if common::run_button(ui, "Seed search") {
        run_kmer(p);
    }
}

/// Run the k-mer seed search — extracted for the headless UI tests.
fn run_kmer(p: &mut AlignmentPanel) {
    p.error = None;
    let recs = common::parse_fasta(&p.kmer_db);
    let query = common::clean_sequence(&p.kmer_query);
    if recs.is_empty() || query.is_empty() {
        p.error = Some("need a reference and a non-empty query".into());
        return;
    }
    let owned: Vec<Vec<u8>> = recs
        .iter()
        .map(|(_, s)| s.to_ascii_uppercase().into_bytes())
        .collect();
    let refs: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
    match KmerIndex::build_many(&refs, p.kmer_k) {
        Ok(index) => {
            let hits = index.seed_query(query.as_bytes());
            let mut out = format!(
                "index: {} sequences · {} distinct {}-mers\n{} seed hits:\n",
                index.sequence_count(),
                index.distinct_kmers(),
                p.kmer_k,
                hits.len(),
            );
            for (q_off, hit) in hits.iter().take(60) {
                let name = recs
                    .get(hit.seq_id)
                    .map(|(l, _)| l.as_str())
                    .filter(|l| !l.is_empty())
                    .map(|l| l.to_string())
                    .unwrap_or_else(|| format!("seq{}", hit.seq_id));
                out.push_str(&format!(
                    "  query@{q_off:<5} → {name} @{}\n",
                    hit.offset,
                ));
            }
            p.result = out;
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_pairwise_aligns() {
        let p = AlignmentPanel::default();
        let a = common::clean_sequence(&p.seq_a);
        let b = common::clean_sequence(&p.seq_b);
        let scheme = scoring_scheme(&p);
        let aln = needleman_wunsch(a.as_bytes(), b.as_bytes(), &scheme).unwrap();
        assert!(aln.len() >= a.len().max(b.len()));
    }

    #[test]
    fn default_msa_has_three_sequences() {
        let p = AlignmentPanel::default();
        let recs = common::parse_fasta(&p.msa_text);
        assert_eq!(recs.len(), 3);
    }
}

/// Headless egui UI-logic tests for the Alignment panel.
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
        app.genetics.active = GeneticsPanel::Alignment;
        app
    }

    #[test]
    fn draws_every_tool_without_panic() {
        for tool in [Tool::Pairwise, Tool::Multiple, Tool::KmerSearch] {
            let mut app = app_with_panel();
            app.genetics.alignment.tool = tool;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_post_run_and_error_states_without_panic() {
        let mut app = app_with_panel();
        app.genetics.alignment.result = "score : 42\n".to_string();
        draw_headless(&mut app);
        let mut app = app_with_panel();
        app.genetics.alignment.error = Some("both sequences must be non-empty".to_string());
        draw_headless(&mut app);
    }

    #[test]
    fn run_pairwise_aligns_with_each_algorithm() {
        // Each pairwise algorithm calls the real valenx-align API and
        // produces a correctly-formatted result.
        for algo in [
            PairAlgo::NeedlemanWunsch,
            PairAlgo::SmithWaterman,
            PairAlgo::GotohAffine,
        ] {
            let mut p = AlignmentPanel {
                algo,
                ..AlignmentPanel::default()
            };
            run_pairwise(&mut p);
            assert!(p.error.is_none(), "{algo:?} errored: {:?}", p.error);
            assert!(p.result.contains("score"));
            assert!(p.result.contains("identity"));
        }
    }

    #[test]
    fn run_msa_builds_a_multiple_alignment() {
        let mut p = AlignmentPanel::default();
        run_msa(&mut p);
        assert!(p.error.is_none(), "MSA errored: {:?}", p.error);
        assert!(p.result.contains("sequences"));
        assert!(p.result.contains("columns"));
    }

    #[test]
    fn run_kmer_finds_seed_hits() {
        let mut p = AlignmentPanel::default();
        run_kmer(&mut p);
        assert!(p.error.is_none(), "k-mer search errored: {:?}", p.error);
        assert!(p.result.contains("seed hits"));
    }

    #[test]
    fn run_actions_surface_errors_on_empty_input() {
        // Pairwise with empty sequences.
        let mut p = AlignmentPanel {
            seq_a: String::new(),
            seq_b: String::new(),
            ..AlignmentPanel::default()
        };
        run_pairwise(&mut p);
        assert!(p.error.is_some(), "pairwise should error on empty input");
        // MSA with a single sequence (needs ≥ 2).
        let mut p = AlignmentPanel {
            msa_text: ">only\nACGTACGT\n".to_string(),
            ..AlignmentPanel::default()
        };
        run_msa(&mut p);
        assert!(p.error.is_some(), "MSA should error with one sequence");
        // k-mer search with an empty query.
        let mut p = AlignmentPanel {
            kmer_query: String::new(),
            ..AlignmentPanel::default()
        };
        run_kmer(&mut p);
        assert!(p.error.is_some(), "k-mer search should error on empty query");
    }
}
