//! The right-side **Genetics Workbench** panel.
//!
//! Valenx ships fifteen native computational-biology library crates
//! (`valenx-bioseq`, `valenx-align`, `valenx-phylo`, `valenx-popgen`,
//! `valenx-rnastruct`, `valenx-rnadesign`, `valenx-md`, `valenx-cheminf`,
//! `valenx-biostruct`, `valenx-qchem`, `valenx-genomics`,
//! `valenx-sysbio`, `valenx-dock-screen`, `valenx-genediting`,
//! `valenx-structpredict`) — but before this module they were libraries
//! + MCP APIs only, with no desktop UI.
//!
//! This module surfaces them as polished egui panels, mirroring the CAD
//! side's [`crate::mesh_toolbox`] idiom: a resizable right-hand
//! [`egui::SidePanel`], a panel selector at the top, and exactly one
//! crate-focused panel shown at a time. Each panel is an organized tool
//! palette + input forms + a Run action that calls the crate's real
//! API + a results area.
//!
//! The **RNA Designer** panel surfaces `valenx-rnadesign` (the unified
//! synthetic-RNA design workflow) as a guided start-to-finish wizard;
//! unlike the other panels it steps through several Back/Next stages
//! rather than a single tool palette.
//!
//! The workbench owns one [`GeneticsWorkbenchState`] (a field on
//! [`ValenxApp`]). Each panel owns a small sub-state struct in the
//! [`crate::genetics`] sub-module; the panels never run external
//! processes — every computation is a pure native call into a Round-6
//! crate.

use eframe::egui;

use crate::genetics;
use crate::ValenxApp;

/// Which of the fifteen bio-crate panels the workbench is showing.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum GeneticsPanel {
    /// `valenx-bioseq` — sequence editing, translation, ORFs, cloning.
    #[default]
    Sequence,
    /// `valenx-align` — pairwise / multiple alignment, k-mer search.
    Alignment,
    /// `valenx-phylo` — distance-method phylogenetic trees.
    Phylogenetics,
    /// `valenx-popgen` — Wright-Fisher / coalescent + summary stats.
    PopulationGenetics,
    /// `valenx-rnastruct` — RNA secondary-structure folding.
    RnaStructure,
    /// `valenx-rnadesign` — guided synthetic-RNA design wizard.
    RnaDesigner,
    /// `valenx-md` — classical molecular-dynamics engine.
    MolecularDynamics,
    /// `valenx-cheminf` — cheminformatics descriptors + fingerprints.
    Cheminformatics,
    /// `valenx-biostruct` — macromolecular structure analysis.
    MacromolecularStructure,
    /// `valenx-qchem` — Hartree-Fock quantum chemistry.
    QuantumChemistry,
    /// `valenx-genomics` — NGS / variant / CRISPR tooling.
    Genomics,
    /// `valenx-sysbio` — systems-biology reaction networks.
    SystemsBiology,
    /// `valenx-dock-screen` — docking + virtual screening.
    Docking,
    /// `valenx-genediting` — CRISPR / base / prime / mRNA design.
    GeneEditing,
    /// `valenx-structpredict` — classical protein structure prediction
    /// + fixed-backbone design (no trained weights).
    StructurePrediction,
}

impl GeneticsPanel {
    /// Every panel, in display order.
    pub const ALL: [GeneticsPanel; 15] = [
        GeneticsPanel::Sequence,
        GeneticsPanel::Alignment,
        GeneticsPanel::Phylogenetics,
        GeneticsPanel::PopulationGenetics,
        GeneticsPanel::RnaStructure,
        GeneticsPanel::RnaDesigner,
        GeneticsPanel::MolecularDynamics,
        GeneticsPanel::Cheminformatics,
        GeneticsPanel::MacromolecularStructure,
        GeneticsPanel::QuantumChemistry,
        GeneticsPanel::Genomics,
        GeneticsPanel::SystemsBiology,
        GeneticsPanel::Docking,
        GeneticsPanel::GeneEditing,
        GeneticsPanel::StructurePrediction,
    ];

    /// Short tab label for the panel selector.
    pub fn label(self) -> &'static str {
        match self {
            GeneticsPanel::Sequence => "Sequence",
            GeneticsPanel::Alignment => "Alignment",
            GeneticsPanel::Phylogenetics => "Phylogenetics",
            GeneticsPanel::PopulationGenetics => "Population Genetics",
            GeneticsPanel::RnaStructure => "RNA Structure",
            GeneticsPanel::RnaDesigner => "RNA Designer",
            GeneticsPanel::MolecularDynamics => "Molecular Dynamics",
            GeneticsPanel::Cheminformatics => "Cheminformatics",
            GeneticsPanel::MacromolecularStructure => "Macromolecular Structure",
            GeneticsPanel::QuantumChemistry => "Quantum Chemistry",
            GeneticsPanel::Genomics => "Genomics",
            GeneticsPanel::SystemsBiology => "Systems Biology",
            GeneticsPanel::Docking => "Docking",
            GeneticsPanel::GeneEditing => "Gene Editing",
            GeneticsPanel::StructurePrediction => "Structure Prediction",
        }
    }

    /// The backing crate name — shown under the panel header.
    pub fn crate_name(self) -> &'static str {
        match self {
            GeneticsPanel::Sequence => "valenx-bioseq",
            GeneticsPanel::Alignment => "valenx-align",
            GeneticsPanel::Phylogenetics => "valenx-phylo",
            GeneticsPanel::PopulationGenetics => "valenx-popgen",
            GeneticsPanel::RnaStructure => "valenx-rnastruct",
            GeneticsPanel::RnaDesigner => "valenx-rnadesign",
            GeneticsPanel::MolecularDynamics => "valenx-md",
            GeneticsPanel::Cheminformatics => "valenx-cheminf",
            GeneticsPanel::MacromolecularStructure => "valenx-biostruct",
            GeneticsPanel::QuantumChemistry => "valenx-qchem",
            GeneticsPanel::Genomics => "valenx-genomics",
            GeneticsPanel::SystemsBiology => "valenx-sysbio",
            GeneticsPanel::Docking => "valenx-dock-screen",
            GeneticsPanel::GeneEditing => "valenx-genediting",
            GeneticsPanel::StructurePrediction => "valenx-structpredict",
        }
    }
}

/// All fifteen panels' form + result state, plus the active-panel
/// selector. One instance lives on the `ValenxApp` (the private
/// `genetics` field), exactly as the CAD-side `MeshToolboxState` does.
#[derive(Default)]
pub struct GeneticsWorkbenchState {
    /// Which panel is currently shown.
    pub active: GeneticsPanel,

    /// Per-panel state. Each is a small struct in the `genetics`
    /// sub-module holding that panel's form inputs + last results.
    pub sequence: genetics::sequence::SequencePanel,
    pub alignment: genetics::alignment::AlignmentPanel,
    pub phylogenetics: genetics::phylogenetics::PhylogeneticsPanel,
    pub popgen: genetics::popgen::PopgenPanel,
    pub rnastruct: genetics::rnastruct::RnaStructPanel,
    pub rna_designer: genetics::rna_designer::RnaDesignerPanel,
    pub md: genetics::md::MdPanel,
    pub cheminf: genetics::cheminf::CheminfPanel,
    pub biostruct: genetics::biostruct::BiostructPanel,
    pub qchem: genetics::qchem::QchemPanel,
    pub genomics: genetics::genomics::GenomicsPanel,
    pub sysbio: genetics::sysbio::SysbioPanel,
    pub docking: genetics::docking::DockingPanel,
    pub genediting: genetics::genediting::GeneEditingPanel,
    pub structpredict: genetics::structpredict::StructPredictPanel,
}

/// Draw the Genetics Workbench right-side panel.
///
/// Mirrors [`crate::mesh_toolbox::draw_mesh_toolbox`]: a no-op when the
/// `show_genetics_workbench` toggle is off, otherwise a resizable
/// [`egui::SidePanel`] mounted before the central viewport so egui
/// docks it to the right.
pub fn draw_genetics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_genetics_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_genetics_workbench",
        "Genetics Workbench",
        genetics_workbench_body,
    );
    if close {
        app.show_genetics_workbench = false;
    }
}

/// The Genetics workbench body — the panel selector plus the dispatch to
/// each of the 15 native computational-biology toolkits. Extracted from
/// [`draw_genetics_workbench`] so it can be hosted by the classic
/// [`crate::workbench_chrome::workbench_shell`] *or* the opt-in dockable
/// tile layout ([`crate::dock_layout`]) without duplicating logic.
pub(crate) fn genetics_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new("15 native computational-biology toolkits")
            .weak()
            .small(),
    );
    ui.separator();

    // --- Panel selector --------------------------------------
    // A combo box for compact switching plus a wrapped grid of
    // selectable chips so the 15 panels are all one click away
    // — the same "tool palette" feel as the CAD workbenches.
    //
    // Each chip carries a hover tooltip sourced from
    // crate::panel_help so the user can see what a panel does
    // before switching to it. F1 once inside a panel opens
    // the full help popup.
    let panel_lbl = ui.label("Bio tool panel");
    egui::ComboBox::from_id_source("genetics_panel_combo")
        .selected_text(app.genetics.active.label())
        .width(ui.available_width())
        .show_ui(ui, |ui| {
            for panel in GeneticsPanel::ALL {
                ui.selectable_value(&mut app.genetics.active, panel, panel.label())
                    .on_hover_text(crate::panel_help::short_summary(panel.label()));
            }
        })
        .response
        .labelled_by(panel_lbl.id);
    ui.add_space(2.0);
    ui.horizontal_wrapped(|ui| {
        for panel in GeneticsPanel::ALL {
            let selected = app.genetics.active == panel;
            if ui
                .selectable_label(selected, panel.label())
                .on_hover_text(crate::panel_help::short_summary(panel.label()))
                .clicked()
            {
                app.genetics.active = panel;
            }
        }
    });
    ui.separator();

    // --- Active panel header ---------------------------------
    let active = app.genetics.active;
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new(active.label()).heading());
    });
    ui.label(
        egui::RichText::new(format!("backed by `{}`", active.crate_name()))
            .weak()
            .small(),
    );
    ui.separator();

    // --- Active panel body -----------------------------------
    //
    // Fade-in animation on panel switch — egui's
    // `animate_bool_with_time` interpolates 0..1 on a state flip
    // (the `false → true` transition). We key the animation id
    // on the active-panel label so switching panels makes a
    // *new* id whose stored value is `false`, then jumps to
    // `true` on this frame → the value crosses 0 to 1 over
    // 0.15 s. Stale ids are GC'd by egui automatically.
    let anim_id = egui::Id::new(("valenx_genetics_panel_switch", active.label()));
    let t = ui.ctx().animate_bool_with_time(anim_id, true, 0.15);
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.scope(|ui| {
                ui.set_opacity(t.clamp(0.0, 1.0));
                match active {
                    GeneticsPanel::Sequence => genetics::sequence::draw(app, ui),
                    GeneticsPanel::Alignment => genetics::alignment::draw(app, ui),
                    GeneticsPanel::Phylogenetics => genetics::phylogenetics::draw(app, ui),
                    GeneticsPanel::PopulationGenetics => genetics::popgen::draw(app, ui),
                    GeneticsPanel::RnaStructure => genetics::rnastruct::draw(app, ui),
                    GeneticsPanel::RnaDesigner => genetics::rna_designer::draw(app, ui),
                    GeneticsPanel::MolecularDynamics => genetics::md::draw(app, ui),
                    GeneticsPanel::Cheminformatics => genetics::cheminf::draw(app, ui),
                    GeneticsPanel::MacromolecularStructure => genetics::biostruct::draw(app, ui),
                    GeneticsPanel::QuantumChemistry => genetics::qchem::draw(app, ui),
                    GeneticsPanel::Genomics => genetics::genomics::draw(app, ui),
                    GeneticsPanel::SystemsBiology => genetics::sysbio::draw(app, ui),
                    GeneticsPanel::Docking => genetics::docking::draw(app, ui),
                    GeneticsPanel::GeneEditing => genetics::genediting::draw(app, ui),
                    GeneticsPanel::StructurePrediction => genetics::structpredict::draw(app, ui),
                }
            });
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_panels_listed_once() {
        // The ALL array must contain every variant exactly once — the
        // selector iterates it, so a missing variant would be
        // unreachable in the UI.
        assert_eq!(GeneticsPanel::ALL.len(), 15);
        for panel in GeneticsPanel::ALL {
            let n = GeneticsPanel::ALL.iter().filter(|p| **p == panel).count();
            assert_eq!(n, 1, "{panel:?} appears {n} times");
        }
    }

    #[test]
    fn every_panel_has_labels() {
        // Labels + crate names feed the UI; none may be empty.
        for panel in GeneticsPanel::ALL {
            assert!(!panel.label().is_empty());
            assert!(panel.crate_name().starts_with("valenx-"));
        }
    }

    #[test]
    fn default_panel_is_sequence() {
        assert_eq!(GeneticsPanel::default(), GeneticsPanel::Sequence);
        assert_eq!(
            GeneticsWorkbenchState::default().active,
            GeneticsPanel::Sequence
        );
    }
}

/// Headless egui UI-logic tests for the Genetics Workbench host panel.
///
/// These mount the whole right-side workbench panel in a windowless
/// [`egui::Context`] and switch through every one of the 15 panels —
/// the workbench-level draw never opens a window and never reaches
/// `rfd::FileDialog`.
#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use crate::ValenxApp;

    /// Run the whole workbench panel once in a headless context.
    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_genetics_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        // With the toggle off the workbench draws nothing and never
        // panics — the default state.
        let mut app = ValenxApp::default();
        assert!(!app.show_genetics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_every_panel_without_panic() {
        // With the workbench shown, mount it once per panel — the panel
        // selector + the active panel body must all render headlessly.
        for panel in GeneticsPanel::ALL {
            let mut app = ValenxApp::default();
            app.show_genetics_workbench = true;
            app.genetics.active = panel;
            draw_workbench(&mut app);
        }
    }
}
