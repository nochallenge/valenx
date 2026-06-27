//! Panel 13 — **Gene Editing** (`valenx-genediting`).
//!
//! Design CRISPR nuclease guide RNAs, base-editing guides, prime-editing
//! pegRNAs and mRNA therapeutic constructs for a target, and run the
//! edit-strategy advisor that ranks editing modalities for a desired
//! change — all native `valenx-genediting` calls. Every score is a
//! transparent rule-based heuristic (no trained weights).

use eframe::egui;

use valenx_genediting::base_edit::design::{design_base_edit, BaseEditRequest};
use valenx_genediting::base_edit::editor::BaseEditorId;
use valenx_genediting::crispr::guide_design::{design_guides, GuideDesignRequest};
use valenx_genediting::crispr::nuclease::NucleaseId;
use valenx_genediting::mrna::design::{design_mrna, MrnaDesignRequest};
use valenx_genediting::mrna::tailcap::MrnaUseCase;
use valenx_genediting::prime_edit::editor::PrimeEditorId;
use valenx_genediting::prime_edit::pegrna::{design_pegrna, PegRnaRequest, PrimeEdit};
use valenx_genediting::workflow::advisor::{advise_strategy, DesiredChange};

use super::common;
use crate::ValenxApp;

/// Which gene-editing sub-tool is showing.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(crate) enum Tool {
    #[default]
    Crispr,
    BaseEdit,
    PrimeEdit,
    Mrna,
    Advisor,
}

/// Snapshot of every editable input the Gene Editing panel owns.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct GeneEditingSnapshot {
    pub(crate) target: String,
    pub(crate) nuclease: NucleaseId,
    pub(crate) base_editor: BaseEditorId,
    pub(crate) prime_editor: PrimeEditorId,
    pub(crate) edit_pos: usize,
    pub(crate) from_base: char,
    pub(crate) to_base: char,
    pub(crate) cds: String,
    pub(crate) mrna_use: MrnaUseCase,
    pub(crate) advisor_change: AdvisorChange,
}

/// Form + result state for the Gene Editing panel.
pub struct GeneEditingPanel {
    tool: Tool,
    /// Target DNA for guide / base / prime editing.
    target: String,
    /// Nuclease for CRISPR guide design.
    nuclease: NucleaseId,
    /// Base editor for the base-edit tool.
    base_editor: BaseEditorId,
    /// Prime editor for the prime-edit tool.
    prime_editor: PrimeEditorId,
    /// Edit position (0-based) for base / prime editing.
    edit_pos: usize,
    /// `from` / `to` bases for a point edit.
    from_base: char,
    to_base: char,
    /// CDS for the mRNA design tool.
    cds: String,
    /// mRNA use-case.
    mrna_use: MrnaUseCase,
    /// Advisor: desired-change kind (0..5 enum-ish selector).
    advisor_change: AdvisorChange,
    error: Option<String>,
    result: String,
    /// Undo / redo over every editable input.
    history: crate::undo::History<GeneEditingSnapshot>,
}

impl GeneEditingPanel {
    fn snapshot(&self) -> GeneEditingSnapshot {
        GeneEditingSnapshot {
            target: self.target.clone(),
            nuclease: self.nuclease,
            base_editor: self.base_editor,
            prime_editor: self.prime_editor,
            edit_pos: self.edit_pos,
            from_base: self.from_base,
            to_base: self.to_base,
            cds: self.cds.clone(),
            mrna_use: self.mrna_use,
            advisor_change: self.advisor_change,
        }
    }
    fn restore(&mut self, s: GeneEditingSnapshot) {
        self.target = s.target;
        self.nuclease = s.nuclease;
        self.base_editor = s.base_editor;
        self.prime_editor = s.prime_editor;
        self.edit_pos = s.edit_pos;
        self.from_base = s.from_base;
        self.to_base = s.to_base;
        self.cds = s.cds;
        self.mrna_use = s.mrna_use;
        self.advisor_change = s.advisor_change;
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

/// The advisor's desired-change menu choice.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(crate) enum AdvisorChange {
    #[default]
    PointTransition,
    PointTransversion,
    SmallInsertion,
    SmallDeletion,
    LargeKnockIn,
    Knockout,
}

impl Default for GeneEditingPanel {
    fn default() -> Self {
        GeneEditingPanel {
            tool: Tool::Crispr,
            target: "ACGTGGAATTCCGGAAGGTTCCAACGGGATGCATGCATGCTAGCTAGCTAGGCATGCATGC".to_string(),
            nuclease: NucleaseId::SpCas9,
            base_editor: BaseEditorId::Be4Max,
            prime_editor: PrimeEditorId::Pe3,
            // Position 35 of the default `target` is a `C` that a
            // guide can place inside the BE4Max activity window, so
            // the base-edit tool (with `from_base` = `C`) opens with a
            // runnable example — and the position is valid for prime
            // editing too.
            edit_pos: 35,
            from_base: 'C',
            to_base: 'T',
            cds: "ATGGCTAGCAAAGGCGAAGAACTGTTCACCGGCGTGGTGCCAATTCTGGTGGAACTGGAT\
                  GGCGATGTGAACGGCCATAAATTCAGCGTGTAA"
                .to_string(),
            mrna_use: MrnaUseCase::Vaccine,
            advisor_change: AdvisorChange::PointTransition,
            error: None,
            result: String::new(),
            history: crate::undo::History::new(),
        }
    }
}

/// Render the Gene Editing panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.genediting;

    common::section(ui, "Tool");
    ui.horizontal_wrapped(|ui| {
        ui.selectable_value(&mut p.tool, Tool::Crispr, "CRISPR guides")
            .on_hover_text(
                "Design guide RNAs against a target with PAM filtering + off-target scoring.",
            );
        ui.selectable_value(&mut p.tool, Tool::BaseEdit, "Base editing")
            .on_hover_text("Pick a base-editor and identify edit windows + bystander positions.");
        ui.selectable_value(&mut p.tool, Tool::PrimeEdit, "Prime editing")
            .on_hover_text("Design pegRNA + nick site for a precise sequence edit.");
        ui.selectable_value(&mut p.tool, Tool::Mrna, "mRNA design");
        ui.selectable_value(&mut p.tool, Tool::Advisor, "Strategy advisor");
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

    match p.tool {
        Tool::Crispr => draw_crispr(p, ui),
        Tool::BaseEdit => draw_base_edit(p, ui),
        Tool::PrimeEdit => draw_prime_edit(p, ui),
        Tool::Mrna => draw_mrna(p, ui),
        Tool::Advisor => draw_advisor(p, ui),
    }

    common::error_line(ui, &p.error);
    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Result");
        common::mono_output(ui, "genediting_result", &p.result, 16);
    }
}

fn draw_crispr(p: &mut GeneEditingPanel, ui: &mut egui::Ui) {
    common::seq_input(ui, "ge_crispr_target", "Target DNA:", &mut p.target, 3);
    ui.horizontal(|ui| {
        ui.label("Nuclease:");
        egui::ComboBox::from_id_source("ge_nuclease")
            .selected_text(format!("{:?}", p.nuclease))
            .show_ui(ui, |ui| {
                for id in NucleaseId::all() {
                    ui.selectable_value(&mut p.nuclease, id, format!("{id:?}"));
                }
            });
    });
    if common::run_button(ui, "Design guide RNAs") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_crispr_design(p);
    }
}

/// Run the CRISPR guide design — extracted for the headless UI tests.
fn run_crispr_design(p: &mut GeneEditingPanel) {
    p.error = None;
    let target = common::clean_sequence(&p.target);
    if target.is_empty() {
        p.error = Some("target sequence is empty".into());
        return;
    }
    let req = GuideDesignRequest::new(target.into_bytes(), p.nuclease);
    match design_guides(&req) {
        Ok(guides) => {
            let mut out = format!(
                "{:?} · {} guide candidate(s) (best design score first):\n",
                p.nuclease,
                guides.len(),
            );
            for (i, g) in guides.iter().take(30).enumerate() {
                out.push_str(&format!(
                    "#{:<3} {} {} @{:<5} {}  on-target {:.2}  spec {:.2}  score {:.3}\n",
                    i + 1,
                    g.protospacer,
                    g.pam,
                    g.start,
                    if g.reverse { "rev" } else { "fwd" },
                    g.on_target,
                    g.specificity,
                    g.design_score,
                ));
            }
            p.result = out;
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

fn draw_base_edit(p: &mut GeneEditingPanel, ui: &mut egui::Ui) {
    common::seq_input(ui, "ge_be_target", "Reference DNA:", &mut p.target, 3);
    egui::Grid::new("ge_be_params")
        .num_columns(2)
        .spacing([8.0, 4.0])
        .show(ui, |ui| {
            ui.label("Editor");
            egui::ComboBox::from_id_source("ge_base_editor")
                .selected_text(format!("{:?}", p.base_editor))
                .show_ui(ui, |ui| {
                    for id in BaseEditorId::all() {
                        ui.selectable_value(&mut p.base_editor, id, format!("{id:?}"));
                    }
                });
            ui.end_row();
            let lbl = ui.label("Edit position");
            ui.add(egui::DragValue::new(&mut p.edit_pos).range(0..=100_000))
                .labelled_by(lbl.id);
            ui.end_row();
            ui.label("From base");
            edit_base_combo(ui, "ge_be_from", &mut p.from_base);
            ui.end_row();
            ui.label("To base");
            edit_base_combo(ui, "ge_be_to", &mut p.to_base);
            ui.end_row();
        });
    if common::run_button(ui, "Design base-editing guide") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_base_edit(p);
    }
}

/// Run the base-editing guide design — extracted for the headless UI
/// tests.
fn run_base_edit(p: &mut GeneEditingPanel) {
    p.error = None;
    let reference = common::clean_sequence(&p.target);
    if p.edit_pos >= reference.len() {
        p.error = Some("edit position is past the end of the reference".into());
        return;
    }
    let req = BaseEditRequest {
        reference: reference.into_bytes(),
        edit_pos: p.edit_pos,
        from_base: p.from_base as u8,
        to_base: p.to_base as u8,
        editor: p.base_editor,
    };
    match design_base_edit(&req) {
        Ok(d) => {
            p.result = format!(
                "{:?} base-editing design\n\n\
                 protospacer    : {}\nPAM            : {}\n\
                 guide start    : {} ({})\n\
                 edit window pos: {}\n\
                 bystanders     : {} other editable base(s)\n\
                 clean edit     : {}\nproduct purity : {:.3}",
                p.base_editor,
                d.protospacer,
                d.pam,
                d.guide_start,
                if d.guide_reverse {
                    "reverse"
                } else {
                    "forward"
                },
                d.window_position,
                d.window.bystander_positions.len(),
                if d.window.clean {
                    "yes"
                } else {
                    "no — bystanders present"
                },
                d.product_purity,
            );
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

fn draw_prime_edit(p: &mut GeneEditingPanel, ui: &mut egui::Ui) {
    common::seq_input(ui, "ge_pe_target", "Reference DNA:", &mut p.target, 3);
    egui::Grid::new("ge_pe_params")
        .num_columns(2)
        .spacing([8.0, 4.0])
        .show(ui, |ui| {
            ui.label("Editor");
            egui::ComboBox::from_id_source("ge_prime_editor")
                .selected_text(format!("{:?}", p.prime_editor))
                .show_ui(ui, |ui| {
                    for id in PrimeEditorId::all() {
                        ui.selectable_value(&mut p.prime_editor, id, format!("{id:?}"));
                    }
                });
            ui.end_row();
            let lbl = ui.label("Edit position");
            ui.add(egui::DragValue::new(&mut p.edit_pos).range(0..=100_000))
                .labelled_by(lbl.id);
            ui.end_row();
            ui.label("Substitute to");
            edit_base_combo(ui, "ge_pe_to", &mut p.to_base);
            ui.end_row();
        });
    if common::run_button(ui, "Design pegRNA") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_prime_edit(p);
    }
}

/// Run the pegRNA design — extracted for the headless UI tests.
fn run_prime_edit(p: &mut GeneEditingPanel) {
    p.error = None;
    let reference = common::clean_sequence(&p.target);
    if p.edit_pos >= reference.len() {
        p.error = Some("edit position is past the end of the reference".into());
        return;
    }
    // A single-base substitution edit.
    let edit = PrimeEdit::Substitution {
        len: 1,
        to: vec![p.to_base as u8],
    };
    let req = PegRnaRequest::new(reference.into_bytes(), p.edit_pos, edit, p.prime_editor);
    match design_pegrna(&req) {
        Ok(peg) => {
            let s = |v: &[u8]| String::from_utf8_lossy(v).into_owned();
            p.result = format!(
                "{:?} pegRNA design\n\n\
                 spacer         : {}\nPAM            : {}\n\
                 spacer start   : {}\nnick position  : {}\n\
                 PBS            : {} ({} nt)\n\
                 RT template    : {} ({} nt)\n\
                 3' extension   : {} ({} nt)\n\n\
                 full pegRNA    :\n{}",
                p.prime_editor,
                s(&peg.spacer),
                s(&peg.pam),
                peg.spacer_start,
                peg.nick_pos,
                s(&peg.pbs),
                peg.pbs.len(),
                s(&peg.rt_template),
                peg.rt_template.len(),
                s(&peg.three_prime_extension),
                peg.three_prime_extension.len(),
                s(&peg.full_pegrna),
            );
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

fn draw_mrna(p: &mut GeneEditingPanel, ui: &mut egui::Ui) {
    common::seq_input(ui, "ge_mrna_cds", "Coding sequence (CDS):", &mut p.cds, 4);
    ui.horizontal(|ui| {
        ui.label("Use case:");
        egui::ComboBox::from_id_source("ge_mrna_use")
            .selected_text(format!("{:?}", p.mrna_use))
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut p.mrna_use, MrnaUseCase::Vaccine, "Vaccine");
                ui.selectable_value(
                    &mut p.mrna_use,
                    MrnaUseCase::ProteinReplacement,
                    "Protein replacement",
                );
                ui.selectable_value(
                    &mut p.mrna_use,
                    MrnaUseCase::TransientEditing,
                    "Transient editing",
                );
            });
    });
    if common::run_button(ui, "Design mRNA construct") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_mrna(p);
    }
}

/// Run the mRNA construct design — extracted for the headless UI tests.
fn run_mrna(p: &mut GeneEditingPanel) {
    p.error = None;
    let cds = common::clean_sequence(&p.cds);
    if cds.is_empty() {
        p.error = Some("CDS is empty".into());
        return;
    }
    let req = MrnaDesignRequest::new(cds.into_bytes(), p.mrna_use);
    match design_mrna(&req) {
        Ok(report) => {
            let mut out = format!(
                "mRNA design — {:?}\n\n\
                 construct length : {} nt\n\
                 CAI before       : {:.3}\nCAI after        : {:.3}\n\
                 uridine fraction : {:.3}\nstart openness   : {:.3}\n\
                 Kozak score      : {:.3}\npoly-A length    : {} nt\n\
                 cap              : {:?}\noverall score    : {:.3}\n",
                p.mrna_use,
                report.construct.len(),
                report.cai_before,
                report.cai_after,
                report.uridine_fraction,
                report.start_openness,
                report.kozak_score,
                report.poly_a_len,
                report.cap,
                report.overall_score,
            );
            if !report.notes.is_empty() {
                out.push_str("\nnotes:\n");
                for n in &report.notes {
                    out.push_str(&format!("  • {n}\n"));
                }
            }
            p.result = out;
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

fn draw_advisor(p: &mut GeneEditingPanel, ui: &mut egui::Ui) {
    common::section(ui, "Desired change");
    egui::ComboBox::from_id_source("ge_advisor_change")
        .selected_text(advisor_change_label(p.advisor_change))
        .width(ui.available_width())
        .show_ui(ui, |ui| {
            for c in [
                AdvisorChange::PointTransition,
                AdvisorChange::PointTransversion,
                AdvisorChange::SmallInsertion,
                AdvisorChange::SmallDeletion,
                AdvisorChange::LargeKnockIn,
                AdvisorChange::Knockout,
            ] {
                ui.selectable_value(&mut p.advisor_change, c, advisor_change_label(c));
            }
        });
    ui.label(
        egui::RichText::new("ranks base / prime / nuclease editing for the change")
            .weak()
            .small(),
    );
    if common::run_button(ui, "Advise editing strategy") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_advisor(p);
    }
}

/// Run the edit-strategy advisor — extracted for the headless UI tests.
fn run_advisor(p: &mut GeneEditingPanel) {
    p.error = None;
    let change = match p.advisor_change {
        // A C->T transition is base-editor reachable.
        AdvisorChange::PointTransition => DesiredChange::PointMutation {
            from: b'C',
            to: b'T',
        },
        // A C->A transversion is not base-editable.
        AdvisorChange::PointTransversion => DesiredChange::PointMutation {
            from: b'C',
            to: b'A',
        },
        AdvisorChange::SmallInsertion => DesiredChange::SmallInsertion { size: 3 },
        AdvisorChange::SmallDeletion => DesiredChange::SmallDeletion { size: 3 },
        AdvisorChange::LargeKnockIn => DesiredChange::LargeKnockIn { size: 800 },
        AdvisorChange::Knockout => DesiredChange::GeneKnockout,
    };
    let advice = advise_strategy(&change, true);
    let mut out = format!(
        "change   : {}\nsummary  : {}\n\n-- ranked approaches --\n",
        advice.change.label(),
        advice.summary,
    );
    for (i, sa) in advice.ranked.iter().enumerate() {
        out.push_str(&format!(
            "#{}  {:<26} suitability {:.2}\n     {}\n",
            i + 1,
            sa.approach.name(),
            sa.suitability,
            sa.rationale,
        ));
    }
    p.result = out;
}

/// A small A/C/G/T combo for an edit-base field.
fn edit_base_combo(ui: &mut egui::Ui, id: &str, base: &mut char) {
    egui::ComboBox::from_id_source(id)
        .selected_text(base.to_string())
        .width(56.0)
        .show_ui(ui, |ui| {
            for b in ['A', 'C', 'G', 'T'] {
                ui.selectable_value(base, b, b.to_string());
            }
        });
}

fn advisor_change_label(c: AdvisorChange) -> &'static str {
    match c {
        AdvisorChange::PointTransition => "Point mutation — transition (C→T)",
        AdvisorChange::PointTransversion => "Point mutation — transversion (C→A)",
        AdvisorChange::SmallInsertion => "Small insertion (3 bp)",
        AdvisorChange::SmallDeletion => "Small deletion (3 bp)",
        AdvisorChange::LargeKnockIn => "Large knock-in (800 bp)",
        AdvisorChange::Knockout => "Gene knockout",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_crispr_designs_guides() {
        let p = GeneEditingPanel::default();
        let target = common::clean_sequence(&p.target);
        let req = GuideDesignRequest::new(target.into_bytes(), NucleaseId::SpCas9);
        // The demo target contains NGG sites — guide design must
        // return at least one candidate (or an empty list, never an
        // error for a valid target).
        assert!(design_guides(&req).is_ok());
    }

    #[test]
    fn advisor_ranks_all_approaches() {
        let advice = advise_strategy(
            &DesiredChange::PointMutation {
                from: b'C',
                to: b'T',
            },
            true,
        );
        // Four modalities are always scored.
        assert_eq!(advice.ranked.len(), 4);
    }

    #[test]
    fn default_mrna_designs() {
        let p = GeneEditingPanel::default();
        let cds = common::clean_sequence(&p.cds);
        let req = MrnaDesignRequest::new(cds.into_bytes(), MrnaUseCase::Vaccine);
        let report = design_mrna(&req).expect("demo CDS must design");
        assert!(!report.construct.is_empty());
    }
}

/// Headless egui UI-logic tests for the Gene Editing panel.
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
        app.genetics.active = GeneticsPanel::GeneEditing;
        app
    }

    #[test]
    fn draws_every_tool_without_panic() {
        for tool in [
            Tool::Crispr,
            Tool::BaseEdit,
            Tool::PrimeEdit,
            Tool::Mrna,
            Tool::Advisor,
        ] {
            let mut app = app_with_panel();
            app.genetics.genediting.tool = tool;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_post_run_and_error_states_without_panic() {
        let mut app = app_with_panel();
        app.genetics.genediting.result = "SpCas9 · 4 guide candidate(s)\n".to_string();
        draw_headless(&mut app);
        let mut app = app_with_panel();
        app.genetics.genediting.error = Some("target sequence is empty".to_string());
        draw_headless(&mut app);
    }

    #[test]
    fn run_crispr_design_designs_guides() {
        // The demo target → the real valenx-genediting guide-design
        // API runs without error.
        let mut p = GeneEditingPanel {
            tool: Tool::Crispr,
            ..GeneEditingPanel::default()
        };
        run_crispr_design(&mut p);
        assert!(p.error.is_none(), "guide design errored: {:?}", p.error);
        assert!(p.result.contains("guide candidate"));
    }

    #[test]
    fn run_base_edit_designs_a_guide() {
        let mut p = GeneEditingPanel {
            tool: Tool::BaseEdit,
            ..GeneEditingPanel::default()
        };
        run_base_edit(&mut p);
        assert!(p.error.is_none(), "base edit errored: {:?}", p.error);
        assert!(p.result.contains("base-editing design"));
    }

    #[test]
    fn run_prime_edit_designs_a_pegrna() {
        let mut p = GeneEditingPanel {
            tool: Tool::PrimeEdit,
            ..GeneEditingPanel::default()
        };
        run_prime_edit(&mut p);
        assert!(p.error.is_none(), "prime edit errored: {:?}", p.error);
        assert!(p.result.contains("pegRNA design"));
    }

    #[test]
    fn run_mrna_designs_a_construct() {
        let mut p = GeneEditingPanel {
            tool: Tool::Mrna,
            ..GeneEditingPanel::default()
        };
        run_mrna(&mut p);
        assert!(p.error.is_none(), "mRNA design errored: {:?}", p.error);
        assert!(p.result.contains("mRNA design"));
        assert!(p.result.contains("overall score"));
    }

    #[test]
    fn run_advisor_ranks_strategies() {
        // The advisor ranks every editing modality for the change.
        for change in [
            AdvisorChange::PointTransition,
            AdvisorChange::PointTransversion,
            AdvisorChange::SmallInsertion,
            AdvisorChange::SmallDeletion,
            AdvisorChange::LargeKnockIn,
            AdvisorChange::Knockout,
        ] {
            let mut p = GeneEditingPanel {
                tool: Tool::Advisor,
                advisor_change: change,
                ..GeneEditingPanel::default()
            };
            run_advisor(&mut p);
            assert!(p.error.is_none(), "advisor errored: {:?}", p.error);
            assert!(p.result.contains("ranked approaches"));
        }
    }

    #[test]
    fn run_actions_surface_errors_on_bad_input() {
        // CRISPR with an empty target.
        let mut p = GeneEditingPanel {
            tool: Tool::Crispr,
            target: String::new(),
            ..GeneEditingPanel::default()
        };
        run_crispr_design(&mut p);
        assert!(
            p.error.is_some(),
            "guide design should error on empty target"
        );
        // Base edit with the edit position past the reference end.
        let mut p = GeneEditingPanel {
            tool: Tool::BaseEdit,
            edit_pos: 100_000,
            ..GeneEditingPanel::default()
        };
        run_base_edit(&mut p);
        assert!(
            p.error.is_some(),
            "base edit should error on an out-of-range position"
        );
        // mRNA design with an empty CDS.
        let mut p = GeneEditingPanel {
            tool: Tool::Mrna,
            cds: String::new(),
            ..GeneEditingPanel::default()
        };
        run_mrna(&mut p);
        assert!(
            p.error.is_some(),
            "mRNA design should error on an empty CDS"
        );
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        use egui::accesskit::Role;
        // BaseEdit and PrimeEdit both expose the edit_pos DragValue.
        for tool in [Tool::BaseEdit, Tool::PrimeEdit] {
            let mut app = app_with_panel();
            app.genetics.genediting.tool = tool;
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
                "genediting tool {tool:?} should expose at least one SpinButton"
            );
            assert!(
                spin_buttons
                    .iter()
                    .all(|(_, n)| !n.labelled_by().is_empty()),
                "every genediting DragValue ({tool:?}) must be labelled_by its caption"
            );
        }
    }
}
