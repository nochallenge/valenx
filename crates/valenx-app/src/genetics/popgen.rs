//! Panel 4 — **Population Genetics** (`valenx-popgen`).
//!
//! Run a forward-time Wright-Fisher simulation or a backward Kingman
//! coalescent, then compute summary statistics — the site-frequency
//! spectrum (rendered as a bar chart), nucleotide diversity π,
//! Watterson's θ, Tajima's D, Fu & Li's D and Hudson's Fst between
//! two sub-samples — all native `valenx-popgen` calls.

use eframe::egui;
use egui_plot::{Bar, BarChart, Plot};

use valenx_popgen::coalescent::kingman::{coalescent, PopHistory};
use valenx_popgen::forward::wright_fisher::{SimulationConfig, WrightFisher};
use valenx_popgen::infer::GenotypeMatrix;
use valenx_popgen::stats::diversity::{
    expected_heterozygosity, fay_wu_h, fu_li_d, minor_allele_frequency, nucleotide_diversity,
    tajimas_d, wattersons_theta,
};
use valenx_popgen::stats::fst::fst_hudson;
use valenx_popgen::stats::sfs::site_frequency_spectrum;

use super::common;
use crate::ValenxApp;

/// Simulation engine choice.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(crate) enum Engine {
    #[default]
    WrightFisher,
    Coalescent,
}

/// Snapshot of every editable input. Numeric-only — Engine choice
/// + every parameter slider — so `Ctrl+Z` reverses a parameter
///   sweep the user wants to back out of.
#[derive(Clone, Debug, PartialEq, Default)]
pub(crate) struct PopgenSnapshot {
    pub(crate) engine: Engine,
    pub(crate) pop_size: usize,
    pub(crate) generations: usize,
    pub(crate) mutation_rate: f64,
    pub(crate) recombination_rate: f64,
    pub(crate) seq_length: f64,
    pub(crate) seed: u64,
}

/// Form + result state for the Population Genetics panel.
pub struct PopgenPanel {
    engine: Engine,
    /// Population / sample size.
    pop_size: usize,
    /// Generations to simulate (Wright-Fisher).
    generations: usize,
    /// Per-base mutation rate.
    mutation_rate: f64,
    /// Per-base recombination rate.
    recombination_rate: f64,
    /// Simulated segment length.
    seq_length: f64,
    /// RNG seed.
    seed: u64,
    error: Option<String>,
    result: String,
    /// SFS bar heights for the plot (`None` until a run produces one).
    sfs: Option<Vec<usize>>,
    /// Undo / redo over the simulation-parameter inputs. Snapshot
    /// recorded on Run so the user can roll back a parameter sweep.
    history: crate::undo::History<PopgenSnapshot>,
}

impl Default for PopgenPanel {
    fn default() -> Self {
        PopgenPanel {
            engine: Engine::WrightFisher,
            pop_size: 50,
            generations: 100,
            mutation_rate: 1.0e-3,
            recombination_rate: 1.0e-4,
            seq_length: 1000.0,
            seed: 42,
            error: None,
            result: String::new(),
            sfs: None,
            history: crate::undo::History::new(),
        }
    }
}

impl PopgenPanel {
    fn snapshot(&self) -> PopgenSnapshot {
        PopgenSnapshot {
            engine: self.engine,
            pop_size: self.pop_size,
            generations: self.generations,
            mutation_rate: self.mutation_rate,
            recombination_rate: self.recombination_rate,
            seq_length: self.seq_length,
            seed: self.seed,
        }
    }
    fn restore(&mut self, s: PopgenSnapshot) {
        self.engine = s.engine;
        self.pop_size = s.pop_size;
        self.generations = s.generations;
        self.mutation_rate = s.mutation_rate;
        self.recombination_rate = s.recombination_rate;
        self.seq_length = s.seq_length;
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

/// Compute + format the summary stats for a genotype matrix.
fn summarize(gm: &GenotypeMatrix) -> (String, Vec<usize>) {
    let mut out = String::new();
    out.push_str(&format!(
        "samples      : {}\nsites        : {} total\nsegregating  : {} sites\n",
        gm.n_samples(),
        gm.n_sites(),
        gm.segregating_sites(),
    ));
    let fmt = |label: &str, r: valenx_popgen::Result<f64>| match r {
        Ok(v) => format!("{label:<13}: {v:.5}\n"),
        Err(e) => format!("{label:<13}: (n/a — {e})\n"),
    };
    out.push_str(&fmt("π (diversity)", nucleotide_diversity(gm)));
    out.push_str(&fmt("Watterson θ", wattersons_theta(gm)));
    out.push_str(&fmt("Tajima's D", tajimas_d(gm)));
    out.push_str(&fmt("Fu & Li's D", fu_li_d(gm)));
    out.push_str(&fmt("Fay & Wu's H", fay_wu_h(gm)));
    out.push_str(&fmt("He (gene div)", expected_heterozygosity(gm)));
    out.push_str(&fmt("MAF (mean)", minor_allele_frequency(gm)));

    // Hudson's Fst between the first and second halves of the sample.
    let rows = gm.rows();
    if rows.len() >= 4 {
        let mid = rows.len() / 2;
        let pos = gm.positions().to_vec();
        let p1 = GenotypeMatrix::from_rows(rows[..mid].to_vec(), pos.clone());
        let p2 = GenotypeMatrix::from_rows(rows[mid..].to_vec(), pos);
        if let (Ok(a), Ok(b)) = (p1, p2) {
            match fst_hudson(&a, &b) {
                Ok(v) => out.push_str(&format!("Fst (halves) : {v:.5}\n")),
                Err(e) => out.push_str(&format!("Fst (halves) : (n/a — {e})\n")),
            }
        }
    }

    let sfs = match site_frequency_spectrum(gm) {
        Ok(s) => {
            out.push_str(&format!(
                "\nSite-frequency spectrum (n={}):\n  singletons: {}\n",
                s.sample_size(),
                s.singletons(),
            ));
            s.counts().to_vec()
        }
        Err(e) => {
            out.push_str(&format!("\nSFS: (n/a — {e})\n"));
            Vec::new()
        }
    };
    (out, sfs)
}

/// Render the Population Genetics panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.popgen;

    common::section(ui, "Engine");
    ui.horizontal(|ui| {
        ui.radio_value(
            &mut p.engine,
            Engine::WrightFisher,
            "Wright-Fisher (forward)",
        )
        .on_hover_text("Forward-time simulator — Wright-Fisher with mutation + recombination.");
        ui.radio_value(&mut p.engine, Engine::Coalescent, "Coalescent (backward)")
            .on_hover_text("Backward-time simulator — Kingman coalescent (Hudson-style).");
    });
    ui.separator();

    common::section(ui, "Parameters");
    egui::Grid::new("popgen_params")
        .num_columns(2)
        .spacing([8.0, 4.0])
        .show(ui, |ui| {
            let lbl = ui
                .label(if p.engine == Engine::WrightFisher {
                    "Population size N"
                } else {
                    "Sample size"
                })
                .on_hover_text(if p.engine == Engine::WrightFisher {
                    "Diploid effective population size."
                } else {
                    "Number of haplotypes sampled from the coalescent."
                });
            ui.add(egui::DragValue::new(&mut p.pop_size).range(2..=2000))
                .labelled_by(lbl.id)
                .on_hover_text("Integer population / sample size.");
            ui.end_row();
            if p.engine == Engine::WrightFisher {
                let lbl = ui
                    .label("Generations")
                    .on_hover_text("Number of forward-time generations to simulate.");
                ui.add(egui::DragValue::new(&mut p.generations).range(1..=5000))
                    .labelled_by(lbl.id)
                    .on_hover_text("Generations of Wright-Fisher iteration.");
                ui.end_row();
                let lbl = ui.label("Mutation rate").on_hover_text(
                    "Per-site, per-generation mutation rate μ (dimensionless 0..=1).",
                );
                ui.add(
                    egui::DragValue::new(&mut p.mutation_rate)
                        .speed(1.0e-4)
                        .range(0.0..=1.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Per-site μ.");
                ui.end_row();
                let lbl = ui
                    .label("Recombination rate")
                    .on_hover_text("Per-pair-of-sites, per-generation recombination probability.");
                ui.add(
                    egui::DragValue::new(&mut p.recombination_rate)
                        .speed(1.0e-5)
                        .range(0.0..=1.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Per-pair-of-sites r.");
                ui.end_row();
                let lbl = ui.label("Segment length").on_hover_text(
                    "Genome segment length in base pairs (used to scale recombination).",
                );
                ui.add(egui::DragValue::new(&mut p.seq_length).range(1.0..=1.0e7))
                    .labelled_by(lbl.id)
                    .on_hover_text("Segment length (bp).");
                ui.end_row();
            }
            let lbl = ui.label("RNG seed");
            ui.add(egui::DragValue::new(&mut p.seed))
                .labelled_by(lbl.id);
            ui.end_row();
        });

    let run_label = match p.engine {
        Engine::WrightFisher => "Run Wright-Fisher",
        Engine::Coalescent => "Run coalescent",
    };
    if common::run_button(ui, run_label) {
        let snap = p.snapshot();
        p.history.record(snap);
        run_engine(p);
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

    if let Some(sfs) = &p.sfs {
        ui.separator();
        common::section(ui, "Site-frequency spectrum");
        let bars: Vec<Bar> = sfs
            .iter()
            .enumerate()
            .skip(1) // class 0 is the monomorphic bin
            .map(|(i, &count)| Bar::new(i as f64, count as f64))
            .collect();
        Plot::new("popgen_sfs_plot")
            .height(140.0)
            .show_axes([true, true])
            .allow_zoom(false)
            .allow_drag(false)
            .show(ui, |plot_ui| {
                plot_ui.bar_chart(BarChart::new(bars).name("derived-allele count"));
            });
    }

    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Summary statistics");
        common::mono_output(ui, "popgen_result", &p.result, 14);
    }
}

/// Run the selected population-genetics engine — extracted from the
/// button closure so it is callable from the headless UI tests.
fn run_engine(p: &mut PopgenPanel) {
    p.error = None;
    p.sfs = None;
    match p.engine {
        Engine::WrightFisher => {
            let cfg = SimulationConfig::neutral(
                p.pop_size,
                p.generations,
                p.mutation_rate,
                p.recombination_rate,
                p.seq_length,
                p.seed,
            );
            match cfg.and_then(WrightFisher::new).and_then(|wf| wf.run()) {
                Ok(pop) => {
                    let gm = pop.genotype_matrix();
                    let (text, sfs) = summarize(&gm);
                    p.result = format!(
                        "Wright-Fisher · N={} · {} generations\n\n{}",
                        pop.size(),
                        p.generations,
                        text,
                    );
                    if !sfs.is_empty() {
                        p.sfs = Some(sfs);
                    }
                }
                Err(e) => p.error = Some(e.to_string()),
            }
        }
        Engine::Coalescent => {
            let labels: Vec<String> = (0..p.pop_size).map(|i| format!("n{i}")).collect();
            // A constant-size diploid history (Ne = N).
            let hist = PopHistory::Constant(p.pop_size as f64);
            match coalescent(&labels, &hist, p.seed) {
                Ok(tree) => {
                    p.result = format!(
                        "Kingman coalescent · {} lineages\n\
                         tree: {} leaves · {} nodes\n\n\
                         A coalescent genealogy is a valenx-phylo Tree — \
                         switch to the Phylogenetics panel to render it, \
                         or compute statistics from a Wright-Fisher run.",
                        p.pop_size,
                        tree.leaf_count(),
                        tree.node_count(),
                    );
                }
                Err(e) => p.error = Some(e.to_string()),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_wright_fisher_runs() {
        let p = PopgenPanel::default();
        let cfg = SimulationConfig::neutral(
            p.pop_size,
            p.generations,
            p.mutation_rate,
            p.recombination_rate,
            p.seq_length,
            p.seed,
        )
        .unwrap();
        let pop = WrightFisher::new(cfg).unwrap().run().unwrap();
        assert!(pop.size() > 0);
    }

    #[test]
    fn coalescent_history_is_valid() {
        let hist = PopHistory::Constant(100.0);
        assert!(hist.validate().is_ok());
    }
}

/// Headless egui UI-logic tests for the Population Genetics panel.
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
        app.genetics.active = GeneticsPanel::PopulationGenetics;
        app
    }

    /// A small Wright-Fisher case so the run stays fast in tests.
    fn small_wf() -> PopgenPanel {
        PopgenPanel {
            engine: Engine::WrightFisher,
            pop_size: 20,
            generations: 20,
            ..PopgenPanel::default()
        }
    }

    #[test]
    fn draws_both_engines_without_panic() {
        for engine in [Engine::WrightFisher, Engine::Coalescent] {
            let mut app = app_with_panel();
            app.genetics.popgen.engine = engine;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_post_run_state_with_sfs_without_panic() {
        // Post-run: a result string + a populated SFS must both render
        // (the SFS drives an egui_plot bar chart).
        let mut app = app_with_panel();
        app.genetics.popgen.result = "Wright-Fisher · N=20\n".to_string();
        app.genetics.popgen.sfs = Some(vec![0, 3, 1, 2]);
        draw_headless(&mut app);
        // Error state.
        let mut app = app_with_panel();
        app.genetics.popgen.error = Some("simulation failed".to_string());
        draw_headless(&mut app);
    }

    #[test]
    fn run_wright_fisher_produces_summary_statistics() {
        // A Wright-Fisher run calls the real valenx-popgen API and
        // produces a correctly-formatted summary.
        let mut p = small_wf();
        run_engine(&mut p);
        assert!(p.error.is_none(), "Wright-Fisher errored: {:?}", p.error);
        assert!(p.result.contains("Wright-Fisher"));
        assert!(p.result.contains("samples"));
    }

    #[test]
    fn run_coalescent_produces_a_genealogy() {
        let mut p = PopgenPanel {
            engine: Engine::Coalescent,
            pop_size: 12,
            ..PopgenPanel::default()
        };
        run_engine(&mut p);
        assert!(p.error.is_none(), "coalescent errored: {:?}", p.error);
        assert!(p.result.contains("Kingman coalescent"));
        assert!(p.result.contains("leaves"));
    }

    #[test]
    fn run_coalescent_handles_a_degenerate_sample() {
        // A single-lineage coalescent is a degenerate / malformed
        // request — the panel must surface an error or an empty-ish
        // result rather than panicking.
        let mut p = PopgenPanel {
            engine: Engine::Coalescent,
            pop_size: 1,
            ..PopgenPanel::default()
        };
        run_engine(&mut p);
        // Either an error is surfaced or a result is produced — never
        // a panic, and never a half-built state.
        assert!(p.error.is_some() || !p.result.is_empty());
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        use egui::accesskit::Role;
        // WrightFisher renders more DragValues (generations, mutation rate, etc.).
        let mut app = app_with_panel();
        app.genetics.popgen.engine = Engine::WrightFisher;
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
            "popgen WrightFisher should expose at least one SpinButton"
        );
        assert!(
            spin_buttons
                .iter()
                .all(|(_, n)| !n.labelled_by().is_empty()),
            "every popgen DragValue must be labelled_by its caption (AI-drivable name)"
        );
    }
}
