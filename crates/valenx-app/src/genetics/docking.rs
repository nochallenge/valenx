//! Panel 12 — **Docking** (`valenx-dock-screen`).
//!
//! Set up a receptor and ligand (PDBQT), run a single docking job, and
//! run structure-based virtual screening over a small ligand library —
//! all native `valenx-dock-screen` calls. The search box defaults to
//! one enclosing the receptor.

use eframe::egui;

use valenx_dock_screen::driver::{dock, screen, DockParams};

use super::common;
use super::molecule_view::{self, ViewAtom, ViewMolecule};
use crate::ValenxApp;

/// Which docking sub-tool is showing.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum Tool {
    #[default]
    Dock,
    Screen,
}

/// Snapshot of every editable input the Docking panel owns.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub(crate) struct DockingSnapshot {
    pub(crate) receptor: String,
    pub(crate) ligand: String,
    pub(crate) library: String,
    pub(crate) n_runs: usize,
    pub(crate) seed: u64,
}

/// Form + result state for the Docking panel.
pub struct DockingPanel {
    tool: Tool,
    /// Receptor PDBQT text.
    receptor: String,
    /// Ligand PDBQT text (single-dock tool).
    ligand: String,
    /// Ligand library for the screening tool — multiple `>name`
    /// FASTA-style records, each body a PDBQT block.
    library: String,
    /// Number of independent search runs.
    n_runs: usize,
    /// RNG seed.
    seed: u64,
    error: Option<String>,
    result: String,
    /// Undo / redo over every editable input.
    history: crate::undo::History<DockingSnapshot>,
}

impl DockingPanel {
    fn snapshot(&self) -> DockingSnapshot {
        DockingSnapshot {
            receptor: self.receptor.clone(),
            ligand: self.ligand.clone(),
            library: self.library.clone(),
            n_runs: self.n_runs,
            seed: self.seed,
        }
    }
    fn restore(&mut self, s: DockingSnapshot) {
        self.receptor = s.receptor;
        self.ligand = s.ligand;
        self.library = s.library;
        self.n_runs = s.n_runs;
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

/// A minimal valid receptor PDBQT — a few well-spaced atoms.
const DEMO_RECEPTOR: &str = "\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ATOM      2  OG1 SER A   2       2.500   0.500   0.000  1.00  0.00    -0.300 OA
ATOM      3  ND1 HIS A   3       1.200   2.500   0.500  1.00  0.00    -0.250 NA
ATOM      4  CA  PHE A   4      -1.500   1.000   1.000  1.00  0.00     0.000 A
ATOM      5  CA  LEU A   5       0.500  -2.000   1.500  1.00  0.00     0.000 C
";

/// A minimal valid ligand PDBQT — two heavy atoms, one rotatable bond.
const DEMO_LIGAND: &str = "\
ROOT
ATOM      1  C1  LIG A   1       0.500   0.500   0.500  1.00  0.00     0.000 C
ENDROOT
BRANCH   1   2
ATOM      2  O1  LIG A   1       1.800   0.500   0.500  1.00  0.00    -0.300 OA
ENDBRANCH   1   2
TORSDOF 1
";

impl Default for DockingPanel {
    fn default() -> Self {
        DockingPanel {
            tool: Tool::Dock,
            receptor: DEMO_RECEPTOR.to_string(),
            ligand: DEMO_LIGAND.to_string(),
            library: format!(">ligand_A\n{DEMO_LIGAND}\n>ligand_B\n{DEMO_LIGAND}"),
            n_runs: 4,
            seed: 1,
            error: None,
            result: String::new(),
            history: crate::undo::History::new(),
        }
    }
}

/// Render the Docking panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.docking;

    common::section(ui, "Tool");
    ui.horizontal(|ui| {
        ui.selectable_value(&mut p.tool, Tool::Dock, "Dock one ligand")
            .on_hover_text("Dock a single ligand into a receptor and rank poses.");
        ui.selectable_value(&mut p.tool, Tool::Screen, "Virtual screening")
            .on_hover_text("Dock a library of ligands and rank them by best-pose affinity.");
        ui.separator();
        let (u, r) = common::undo_redo_inline(ui, p.can_undo(), p.can_redo());
        if u {
            p.undo_edit();
        }
        if r {
            p.redo_edit();
        }
    });
    ui.separator();

    common::section(ui, "Receptor (PDBQT)");
    ui.horizontal(|ui| {
        if ui.small_button("Load receptor…").clicked() {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("PDBQT", &["pdbqt"])
                .pick_file()
            {
                // Round-21 H1: bounded read so a multi-GB PDBQT
                // path can't OOM the renderer (sister to the dock
                // panel's MCP-side cap on the same format).
                match valenx_core::io_caps::read_capped_to_string(
                    &path,
                    valenx_core::io_caps::MAX_GENETICS_FILE_BYTES as usize,
                ) {
                    Ok(t) => p.receptor = t,
                    Err(e) => p.error = Some(format!("read: {e}")),
                }
            }
        }
        ui.label(format!("{} lines", p.receptor.lines().count()));
    });
    ui.add(
        egui::TextEdit::multiline(&mut p.receptor)
            .id_source("docking_receptor")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(5),
    );

    match p.tool {
        Tool::Dock => draw_dock(p, ui),
        Tool::Screen => draw_screen(p, ui),
    }

    common::error_line(ui, &p.error);

    // --- 3-D viewport integration ---------------------------------
    // Push the receptor (single-dock tool also adds the ligand atoms)
    // into the app's wgpu 3-D viewport as a ball-and-stick model.
    if !app.genetics.docking.receptor.trim().is_empty() {
        ui.horizontal(|ui| {
            let with_ligand = app.genetics.docking.tool == Tool::Dock;
            let btn = if with_ligand {
                "Show receptor + ligand in 3D viewport"
            } else {
                "Show receptor in 3D viewport"
            };
            if ui.button(btn).clicked() {
                show_in_viewport(app);
            }
        });
    }

    let p = &app.genetics.docking;
    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Result");
        common::mono_output(ui, "docking_result", &p.result, 14);
    }
}

/// Parse `ATOM` / `HETATM` records out of a PDBQT block into view
/// atoms.
///
/// PDBQT shares the PDB fixed-column coordinate layout (columns 31-54,
/// 1-based); the trailing AutoDock atom-type token is ignored and the
/// element is guessed from the atom name (column 13-16). Bonds are
/// inferred downstream by the covalent-radius rule.
fn parse_pdbqt_atoms(text: &str) -> Vec<ViewAtom> {
    let mut atoms = Vec::new();
    for line in text.lines() {
        let record = line.get(0..6).unwrap_or("").trim();
        if record != "ATOM" && record != "HETATM" {
            continue;
        }
        let col =
            |a: usize, b: usize| -> &str { line.get(a..b.min(line.len())).unwrap_or("").trim() };
        let name = col(12, 16);
        let x = col(30, 38).parse::<f32>();
        let y = col(38, 46).parse::<f32>();
        let z = col(46, 54).parse::<f32>();
        if let (Ok(x), Ok(y), Ok(z)) = (x, y, z) {
            atoms.push(ViewAtom::new([x, y, z], guess_element(name)));
        }
    }
    atoms
}

/// Best-effort element symbol from a PDB(QT) atom name — the two-letter
/// biomolecular elements first, otherwise the first alphabetic char.
fn guess_element(atom_name: &str) -> String {
    let upper = atom_name.trim().to_ascii_uppercase();
    for two in ["CL", "NA", "MG", "FE", "ZN", "CA", "BR"] {
        if upper.starts_with(two) {
            return two.to_string();
        }
    }
    upper
        .chars()
        .find(|c| c.is_ascii_alphabetic())
        .map(|c| c.to_string())
        .unwrap_or_default()
}

/// Build the receptor (plus the ligand for the single-dock tool) into
/// a ball-and-stick mesh and push it into the app's 3-D viewport.
fn show_in_viewport(app: &mut ValenxApp) {
    let p = &app.genetics.docking;
    let mut atoms = parse_pdbqt_atoms(&p.receptor);
    if p.tool == Tool::Dock {
        atoms.extend(parse_pdbqt_atoms(&p.ligand));
    }
    if atoms.is_empty() {
        app.genetics.docking.error = Some("no ATOM records found in the PDBQT text".to_string());
        return;
    }
    let bonds = molecule_view::detect_bonds(&atoms);
    let view = ViewMolecule { atoms, bonds };
    let mesh = molecule_view::ball_and_stick(&view, 0.26, 0.17);
    match molecule_view::show_molecule(app, mesh, "docking-complex.ball-stick") {
        Ok(_) => app.genetics.docking.error = None,
        Err(e) => app.genetics.docking.error = Some(e),
    }
}

fn draw_dock(p: &mut DockingPanel, ui: &mut egui::Ui) {
    common::section(ui, "Ligand (PDBQT)");
    ui.horizontal(|ui| {
        if ui.small_button("Load ligand…").clicked() {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("PDBQT", &["pdbqt"])
                .pick_file()
            {
                // Round-21 H1: see receptor loader above.
                match valenx_core::io_caps::read_capped_to_string(
                    &path,
                    valenx_core::io_caps::MAX_GENETICS_FILE_BYTES as usize,
                ) {
                    Ok(t) => p.ligand = t,
                    Err(e) => p.error = Some(format!("read: {e}")),
                }
            }
        }
    });
    ui.add(
        egui::TextEdit::multiline(&mut p.ligand)
            .id_source("docking_ligand")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(5),
    );
    ui.horizontal(|ui| {
        let lbl_runs = ui.label("Search runs:").on_hover_text(
            "Number of independent docking restarts. More = better pose coverage but slower.",
        );
        ui.add(egui::DragValue::new(&mut p.n_runs).range(1..=64))
            .labelled_by(lbl_runs.id)
            .on_hover_text("Independent search restarts (typical 8–16).");
        let lbl_seed = ui
            .label("seed:")
            .on_hover_text("RNG seed — same seed reproduces the same poses.");
        ui.add(egui::DragValue::new(&mut p.seed))
            .labelled_by(lbl_seed.id)
            .on_hover_text("Reproducibility seed.");
    });

    if common::run_button(ui, "Dock") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_dock(p);
    }
}

/// Run the single-ligand docking — extracted for the headless UI tests.
fn run_dock(p: &mut DockingPanel) {
    p.error = None;
    let params = DockParams {
        n_runs: p.n_runs,
        seed: p.seed,
        ..DockParams::fast()
    };
    match dock(&p.receptor, &p.ligand, &params) {
        Ok(report) => {
            let mut out = format!(
                "docking complete\nbinding modes : {}\nposes sampled : {}\n\n",
                report.n_binding_modes(),
                report.poses.len(),
            );
            match report.best_score() {
                Some(s) => out.push_str(&format!("best score    : {s:.3} kcal/mol\n\n")),
                None => out.push_str("best score    : (none)\n\n"),
            }
            out.push_str("-- pose clusters (best first) --\n");
            for (i, cluster) in report.clusters.iter().take(20).enumerate() {
                out.push_str(&format!(
                    "  #{:<3} score {:>9.3} kcal/mol · mean {:>9.3} kcal/mol · {} member(s)\n",
                    i + 1,
                    cluster.best_score,
                    cluster.mean_score(),
                    cluster.members.len(),
                ));
            }
            p.result = out;
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

fn draw_screen(p: &mut DockingPanel, ui: &mut egui::Ui) {
    common::section(ui, "Ligand library");
    ui.label(
        egui::RichText::new("multiple >name records, each body a PDBQT block")
            .weak()
            .small(),
    );
    ui.add(
        egui::TextEdit::multiline(&mut p.library)
            .id_source("docking_library")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(7),
    );
    ui.horizontal(|ui| {
        let lbl_runs = ui.label("Search runs:").on_hover_text(
            "Number of independent docking restarts. More = better pose coverage but slower.",
        );
        ui.add(egui::DragValue::new(&mut p.n_runs).range(1..=64))
            .labelled_by(lbl_runs.id)
            .on_hover_text("Independent search restarts (typical 8–16).");
        let lbl_seed = ui
            .label("seed:")
            .on_hover_text("RNG seed — same seed reproduces the same poses.");
        ui.add(egui::DragValue::new(&mut p.seed))
            .labelled_by(lbl_seed.id)
            .on_hover_text("Reproducibility seed.");
    });

    if common::run_button(ui, "Run virtual screen") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_screen(p);
    }
}

/// Run the virtual screen — extracted for the headless UI tests.
fn run_screen(p: &mut DockingPanel) {
    p.error = None;
    let library = parse_pdbqt_library(&p.library);
    if library.is_empty() {
        p.error = Some("library is empty — add >name records with PDBQT bodies".into());
        return;
    }
    let params = DockParams {
        n_runs: p.n_runs,
        seed: p.seed,
        ..DockParams::fast()
    };
    match screen(&p.receptor, &library, &params) {
        Ok(report) => {
            let mut out = format!(
                "virtual screen complete\nscreened      : {}\nfailed        : {}\n\n\
                     -- ranked hits (best score first) --\n",
                report.n_screened, report.n_failed,
            );
            let mut hits = report.hits.clone();
            hits.sort_by(|a, b| {
                let av = a.best_score.unwrap_or(f64::INFINITY);
                let bv = b.best_score.unwrap_or(f64::INFINITY);
                av.partial_cmp(&bv).unwrap_or(std::cmp::Ordering::Equal)
            });
            for (i, hit) in hits.iter().enumerate() {
                match (&hit.best_score, &hit.failure) {
                    (Some(s), _) => out.push_str(&format!(
                        "  #{:<3} {:<16} {:>9.3} kcal/mol\n",
                        i + 1,
                        hit.name,
                        s,
                    )),
                    (None, Some(f)) => out.push_str(&format!(
                        "  #{:<3} {:<16} FAILED — {}\n",
                        i + 1,
                        hit.name,
                        f,
                    )),
                    (None, None) => {
                        out.push_str(&format!("  #{:<3} {:<16} (no score)\n", i + 1, hit.name,))
                    }
                }
            }
            p.result = out;
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

/// Split a `>name`-delimited blob into `(name, pdbqt)` library entries.
/// Unlike a FASTA parse the PDBQT bodies keep their internal newlines
/// and spacing intact.
fn parse_pdbqt_library(raw: &str) -> Vec<(String, String)> {
    let mut out: Vec<(String, String)> = Vec::new();
    let mut name = String::new();
    let mut body = String::new();
    let mut started = false;
    for line in raw.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if started {
                out.push((std::mem::take(&mut name), std::mem::take(&mut body)));
            }
            name = rest.trim().to_string();
            started = true;
        } else if started {
            body.push_str(line);
            body.push('\n');
        }
    }
    if started {
        out.push((name, body));
    }
    out.into_iter()
        .filter(|(_, b)| b.contains("ATOM"))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_dock_runs() {
        let p = DockingPanel::default();
        let params = DockParams {
            n_runs: 2,
            seed: 1,
            ..DockParams::fast()
        };
        let report = dock(&p.receptor, &p.ligand, &params).expect("demo dock must run");
        assert!(!report.poses.is_empty());
    }

    #[test]
    fn library_parses_two_entries() {
        let p = DockingPanel::default();
        let lib = parse_pdbqt_library(&p.library);
        assert_eq!(lib.len(), 2);
        assert!(lib.iter().all(|(_, b)| b.contains("ROOT")));
    }
}

/// Headless egui UI-logic tests for the Docking panel.
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
        app.genetics.active = GeneticsPanel::Docking;
        app
    }

    #[test]
    fn draws_both_tools_without_panic() {
        for tool in [Tool::Dock, Tool::Screen] {
            let mut app = app_with_panel();
            app.genetics.docking.tool = tool;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_post_run_and_error_states_without_panic() {
        let mut app = app_with_panel();
        app.genetics.docking.result = "docking complete\nbinding modes : 2\n".to_string();
        draw_headless(&mut app);
        let mut app = app_with_panel();
        app.genetics.docking.error = Some("library is empty".to_string());
        draw_headless(&mut app);
        // Empty receptor — the "Show in 3D" affordance is gated off.
        let mut app = app_with_panel();
        app.genetics.docking.receptor.clear();
        draw_headless(&mut app);
    }

    #[test]
    fn run_dock_docks_the_demo_complex() {
        // The demo receptor + ligand → the real valenx-dock-screen
        // docking driver produces a correctly-formatted report.
        let mut p = DockingPanel {
            tool: Tool::Dock,
            n_runs: 2,
            ..DockingPanel::default()
        };
        run_dock(&mut p);
        assert!(p.error.is_none(), "dock errored: {:?}", p.error);
        assert!(p.result.contains("docking complete"));
        assert!(p.result.contains("poses sampled"));
    }

    #[test]
    fn run_screen_screens_the_demo_library() {
        let mut p = DockingPanel {
            tool: Tool::Screen,
            n_runs: 2,
            ..DockingPanel::default()
        };
        run_screen(&mut p);
        assert!(p.error.is_none(), "screen errored: {:?}", p.error);
        assert!(p.result.contains("virtual screen complete"));
    }

    #[test]
    fn run_screen_surfaces_error_on_empty_library() {
        // An empty library is malformed input — the panel must
        // surface an error rather than panicking.
        let mut p = DockingPanel {
            tool: Tool::Screen,
            library: String::new(),
            ..DockingPanel::default()
        };
        run_screen(&mut p);
        assert!(p.error.is_some(), "screen should error on empty library");
        assert!(p.result.is_empty());
    }

    #[test]
    fn run_dock_surfaces_error_on_malformed_receptor() {
        // A receptor with no ATOM records cannot dock — the panel
        // surfaces an error rather than panicking.
        let mut p = DockingPanel {
            tool: Tool::Dock,
            receptor: "not a pdbqt file".to_string(),
            n_runs: 2,
            ..DockingPanel::default()
        };
        run_dock(&mut p);
        assert!(
            p.error.is_some(),
            "dock should error on a malformed receptor"
        );
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        use egui::accesskit::Role;
        for tool in [Tool::Dock, Tool::Screen] {
            let mut app = app_with_panel();
            app.genetics.docking.tool = tool;
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
                "docking tool {tool:?} should expose at least one SpinButton"
            );
            assert!(
                spin_buttons
                    .iter()
                    .all(|(_, n)| !n.labelled_by().is_empty()),
                "every docking DragValue ({tool:?}) must be labelled_by its caption"
            );
        }
    }
}
