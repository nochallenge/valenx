//! Panel 11 — **Systems Biology** (`valenx-sysbio`).
//!
//! A free-form reaction-network editor. The user adds and removes
//! species rows (name + initial concentration) and reaction rows
//! (reactant + product species with stoichiometry, plus a rate law —
//! mass-action, Michaelis-Menten or Hill — with its parameters). The
//! panel assembles that into a real [`valenx_sysbio::Model`] and runs
//! it three ways:
//!
//! - a deterministic ODE time course (RK45) via [`TimeCourse`],
//! - a stochastic Gillespie SSA via [`StochasticModel`],
//! - flux-balance analysis via [`FbaProblem`] when the network has the
//!   structure FBA needs.
//!
//! The earlier v1 hard-wired a fixed `A → B → C` chain; this editor
//! replaces that with a genuine model the user builds from nothing.

use eframe::egui;
use egui_plot::{Legend, Line, PlotPoints};

use valenx_sysbio::fba::FbaProblem;
use valenx_sysbio::model::{Model, RateLaw, Reaction, Species};
use valenx_sysbio::ode::TimeCourse;
use valenx_sysbio::stochastic::StochasticModel;

use super::common;
use crate::plot_ui::managed_plot_mem_cfg;
use crate::ValenxApp;

/// Which simulation engine is showing.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(crate) enum Engine {
    #[default]
    Ode,
    Gillespie,
    Fba,
}

/// The kinetic-law family a [`ReactionRow`] uses. Stored separately
/// from the parameters so the UI can switch laws without losing them.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(crate) enum LawKind {
    /// Zeroth-order constant flux.
    Constant,
    /// Mass action: `k · ∏ [reactant]^stoich`.
    #[default]
    MassAction,
    /// Michaelis-Menten: `Vmax · [S] / (Km + [S])`.
    MichaelisMenten,
    /// Hill: `Vmax · [S]^n / (Kd^n + [S]^n)` (activation).
    Hill,
}

impl LawKind {
    fn label(self) -> &'static str {
        match self {
            LawKind::Constant => "Constant flux",
            LawKind::MassAction => "Mass action",
            LawKind::MichaelisMenten => "Michaelis-Menten",
            LawKind::Hill => "Hill",
        }
    }
}

/// One editable species row.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct SpeciesRow {
    /// Species identifier (must be unique + non-empty to assemble).
    name: String,
    /// Initial amount / concentration.
    initial: f64,
}

/// One editable reaction row.
///
/// `reactants` / `products` are `(species_index, stoichiometry)`
/// pairs into the species list. The rate law is described by `law`
/// plus the four generic parameter slots `k`, `param2`, `param3` and
/// `regulator` — which of them is read depends on `law`:
///
/// - `Constant`: `k` is the flux.
/// - `MassAction`: `k` is the rate constant; the kinetic orders are
///   the reactant stoichiometries.
/// - `MichaelisMenten`: `k` = Vmax, `param2` = Km, the substrate is
///   the first reactant.
/// - `Hill`: `k` = Vmax, `param2` = Kd, `param3` = n, `regulator` is
///   the regulator species index.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ReactionRow {
    /// Reaction identifier.
    name: String,
    /// Reactant `(species index, stoichiometry)` pairs.
    reactants: Vec<(usize, f64)>,
    /// Product `(species index, stoichiometry)` pairs.
    products: Vec<(usize, f64)>,
    /// Rate-law family.
    law: LawKind,
    /// First parameter — `k` / `Vmax` / flux depending on `law`.
    k: f64,
    /// Second parameter — `Km` / `Kd`.
    param2: f64,
    /// Third parameter — Hill coefficient `n`.
    param3: f64,
    /// Regulator species index (Hill law only).
    regulator: usize,
}

impl ReactionRow {
    /// A blank mass-action reaction.
    fn blank(name: impl Into<String>) -> Self {
        ReactionRow {
            name: name.into(),
            reactants: Vec::new(),
            products: Vec::new(),
            law: LawKind::MassAction,
            k: 0.1,
            param2: 1.0,
            param3: 2.0,
            regulator: 0,
        }
    }
}

/// Snapshot of every editable input — the species + reactions lists,
/// engine choice, and the numeric controls.
#[derive(Clone, Debug, PartialEq, Default)]
pub(crate) struct SysbioSnapshot {
    pub(crate) engine: Engine,
    pub(crate) species: Vec<SpeciesRow>,
    pub(crate) reactions: Vec<ReactionRow>,
    pub(crate) t_end: f64,
    pub(crate) seed: u64,
    pub(crate) fba_objective: usize,
}

/// Form + result state for the Systems Biology panel.
pub struct SysbioPanel {
    engine: Engine,
    /// Editable species list.
    species: Vec<SpeciesRow>,
    /// Editable reaction list.
    reactions: Vec<ReactionRow>,
    /// Simulation end time (ODE / Gillespie).
    t_end: f64,
    /// RNG seed for the stochastic run.
    seed: u64,
    /// Index of the reaction maximised by FBA.
    fba_objective: usize,
    error: Option<String>,
    result: String,
    /// Time-course series: `(species_name, points)`.
    series: Vec<(String, Vec<[f64; 2]>)>,
    /// Undo / redo over the whole model definition.
    history: crate::undo::History<SysbioSnapshot>,
}

impl SysbioPanel {
    fn snapshot(&self) -> SysbioSnapshot {
        SysbioSnapshot {
            engine: self.engine,
            species: self.species.clone(),
            reactions: self.reactions.clone(),
            t_end: self.t_end,
            seed: self.seed,
            fba_objective: self.fba_objective,
        }
    }
    fn restore(&mut self, s: SysbioSnapshot) {
        self.engine = s.engine;
        self.species = s.species;
        self.reactions = s.reactions;
        self.t_end = s.t_end;
        self.seed = s.seed;
        self.fba_objective = s.fba_objective;
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

impl Default for SysbioPanel {
    fn default() -> Self {
        // Seed with the classic A → B → C chain so the panel opens
        // with a runnable example — but now every part is editable.
        SysbioPanel {
            engine: Engine::Ode,
            species: vec![
                SpeciesRow {
                    name: "A".to_string(),
                    initial: 100.0,
                },
                SpeciesRow {
                    name: "B".to_string(),
                    initial: 0.0,
                },
                SpeciesRow {
                    name: "C".to_string(),
                    initial: 0.0,
                },
            ],
            reactions: vec![
                ReactionRow {
                    name: "influx".to_string(),
                    reactants: vec![],
                    products: vec![(0, 1.0)],
                    law: LawKind::Constant,
                    k: 5.0,
                    param2: 1.0,
                    param3: 2.0,
                    regulator: 0,
                },
                ReactionRow {
                    name: "a_to_b".to_string(),
                    reactants: vec![(0, 1.0)],
                    products: vec![(1, 1.0)],
                    law: LawKind::MassAction,
                    k: 0.30,
                    param2: 1.0,
                    param3: 2.0,
                    regulator: 0,
                },
                ReactionRow {
                    name: "b_to_c".to_string(),
                    reactants: vec![(1, 1.0)],
                    products: vec![(2, 1.0)],
                    law: LawKind::MassAction,
                    k: 0.15,
                    param2: 1.0,
                    param3: 2.0,
                    regulator: 0,
                },
            ],
            t_end: 30.0,
            seed: 1,
            fba_objective: 2,
            error: None,
            result: String::new(),
            series: Vec::new(),
            history: crate::undo::History::new(),
        }
    }
}

/// Assemble a [`valenx_sysbio::Model`] from the editor rows.
///
/// Validates that the network is well-formed enough to simulate:
/// at least one species, every species name unique and non-empty,
/// every reaction name non-empty. Returns the assembled model or a
/// human-readable reason it could not be built.
fn assemble_model(p: &SysbioPanel) -> Result<Model, String> {
    if p.species.is_empty() {
        return Err("add at least one species".to_string());
    }
    for (i, s) in p.species.iter().enumerate() {
        if s.name.trim().is_empty() {
            return Err(format!("species {} has no name", i + 1));
        }
        if p.species.iter().skip(i + 1).any(|o| o.name == s.name) {
            return Err(format!("duplicate species name `{}`", s.name));
        }
    }
    if p.reactions.is_empty() {
        return Err("add at least one reaction".to_string());
    }

    let mut model = Model::new("genetics_network");
    for s in &p.species {
        model.add_species(Species::new(s.name.trim(), s.initial));
    }
    let ns = p.species.len();
    for (ri, r) in p.reactions.iter().enumerate() {
        if r.name.trim().is_empty() {
            return Err(format!("reaction {} has no name", ri + 1));
        }
        // Guard every species index — a removed species could leave a
        // stale reference behind.
        let in_range = |pairs: &[(usize, f64)]| pairs.iter().all(|&(i, _)| i < ns);
        if !in_range(&r.reactants) || !in_range(&r.products) {
            return Err(format!(
                "reaction `{}` references a deleted species — re-pick its species",
                r.name
            ));
        }
        let rate_law = build_rate_law(r)?;
        model.add_reaction(Reaction {
            id: r.name.trim().to_string(),
            reactants: r.reactants.clone(),
            products: r.products.clone(),
            rate_law,
            reversible: false,
        });
    }
    model
        .validate()
        .map_err(|e| format!("invalid model: {e}"))?;
    Ok(model)
}

/// Build the [`RateLaw`] for one reaction row.
fn build_rate_law(r: &ReactionRow) -> Result<RateLaw, String> {
    match r.law {
        LawKind::Constant => Ok(RateLaw::Constant { rate: r.k }),
        LawKind::MassAction => Ok(RateLaw::MassAction {
            k: r.k,
            // Kinetic orders are the reactant stoichiometries.
            reactants: r.reactants.clone(),
        }),
        LawKind::MichaelisMenten => {
            let substrate = r.reactants.first().map(|&(i, _)| i).ok_or_else(|| {
                format!(
                    "reaction `{}` is Michaelis-Menten but has no reactant \
                         to use as the substrate",
                    r.name
                )
            })?;
            Ok(RateLaw::MichaelisMenten {
                vmax: r.k,
                km: r.param2,
                substrate,
            })
        }
        LawKind::Hill => Ok(RateLaw::Hill {
            vmax: r.k,
            kd: r.param2,
            n: r.param3,
            regulator: r.regulator,
            repress: false,
        }),
    }
}

/// Render the Systems Biology panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.sysbio;

    draw_species_editor(p, ui);
    ui.separator();
    draw_reaction_editor(p, ui);
    ui.separator();

    common::section(ui, "Engine");
    ui.horizontal_wrapped(|ui| {
        ui.radio_value(&mut p.engine, Engine::Ode, "ODE (RK45)");
        ui.radio_value(&mut p.engine, Engine::Gillespie, "Gillespie SSA");
        ui.radio_value(&mut p.engine, Engine::Fba, "Flux-balance analysis");
    });
    if p.engine != Engine::Fba {
        ui.horizontal(|ui| {
            let lbl = ui.label("End time:");
            ui.add(
                egui::DragValue::new(&mut p.t_end)
                    .speed(1.0)
                    .range(1.0..=1000.0),
            )
            .labelled_by(lbl.id);
            if p.engine == Engine::Gillespie {
                let lbl = ui.label("seed:");
                ui.add(egui::DragValue::new(&mut p.seed))
                    .labelled_by(lbl.id);
            }
        });
    } else if !p.reactions.is_empty() {
        ui.horizontal(|ui| {
            ui.label("Maximise flux of:");
            let sel = p
                .reactions
                .get(p.fba_objective)
                .map(|r| r.name.clone())
                .unwrap_or_else(|| "(pick)".to_string());
            egui::ComboBox::from_id_source("sysbio_fba_obj")
                .selected_text(sel)
                .show_ui(ui, |ui| {
                    for (i, r) in p.reactions.iter().enumerate() {
                        ui.selectable_value(&mut p.fba_objective, i, &r.name);
                    }
                });
        });
    }

    let run_label = match p.engine {
        Engine::Ode => "Run ODE time course",
        Engine::Gillespie => "Run Gillespie SSA",
        Engine::Fba => "Run flux-balance analysis",
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

    if !p.series.is_empty() {
        ui.separator();
        common::section(ui, "Time course");
        managed_plot_mem_cfg(
            ui,
            "sysbio_plot",
            170.0,
            |plot| plot.legend(Legend::default()),
            |plot_ui| {
                for (name, pts) in &p.series {
                    plot_ui.line(Line::new(PlotPoints::from(pts.clone())).name(name));
                }
            },
        );
    }

    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Result");
        common::mono_output(ui, "sysbio_result", &p.result, 12);
    }
}

/// Draw the add/remove species editor.
fn draw_species_editor(p: &mut SysbioPanel, ui: &mut egui::Ui) {
    common::section(ui, "Species");
    let mut remove: Option<usize> = None;
    egui::Grid::new("sysbio_species_grid")
        .num_columns(3)
        .spacing([8.0, 4.0])
        .show(ui, |ui| {
            ui.label(egui::RichText::new("name").weak().small());
            let lbl_initial_hdr = ui.label(egui::RichText::new("initial amount").weak().small());
            ui.label("");
            ui.end_row();
            for (i, s) in p.species.iter_mut().enumerate() {
                ui.add(
                    egui::TextEdit::singleline(&mut s.name)
                        .id_source(format!("sysbio_sp_name_{i}"))
                        .desired_width(80.0),
                );
                ui.add(
                    egui::DragValue::new(&mut s.initial)
                        .speed(1.0)
                        .range(0.0..=1.0e9),
                )
                .labelled_by(lbl_initial_hdr.id);
                if ui.small_button("✖").clicked() {
                    remove = Some(i);
                }
                ui.end_row();
            }
        });
    if let Some(i) = remove {
        p.species.remove(i);
        prune_species_refs(p, i);
    }
    if ui.small_button("➕ Add species").clicked() {
        let name = next_name(
            "S",
            &p.species.iter().map(|s| s.name.clone()).collect::<Vec<_>>(),
        );
        p.species.push(SpeciesRow { name, initial: 0.0 });
    }
}

/// Draw the add/remove reaction editor.
fn draw_reaction_editor(p: &mut SysbioPanel, ui: &mut egui::Ui) {
    common::section(ui, "Reactions");
    let species_names: Vec<String> = p.species.iter().map(|s| s.name.clone()).collect();
    let mut remove: Option<usize> = None;
    for ri in 0..p.reactions.len() {
        let header = {
            let r = &p.reactions[ri];
            format!("{}  ({})", r.name, r.law.label())
        };
        egui::CollapsingHeader::new(header)
            .id_source(format!("sysbio_rxn_{ri}"))
            .default_open(ri == 0)
            .show(ui, |ui| {
                let r = &mut p.reactions[ri];
                ui.horizontal(|ui| {
                    ui.label("name:");
                    ui.add(
                        egui::TextEdit::singleline(&mut r.name)
                            .id_source(format!("sysbio_rxn_name_{ri}"))
                            .desired_width(110.0),
                    );
                    if ui.small_button("Delete reaction").clicked() {
                        remove = Some(ri);
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("rate law:");
                    egui::ComboBox::from_id_source(format!("sysbio_rxn_law_{ri}"))
                        .selected_text(r.law.label())
                        .show_ui(ui, |ui| {
                            for k in [
                                LawKind::Constant,
                                LawKind::MassAction,
                                LawKind::MichaelisMenten,
                                LawKind::Hill,
                            ] {
                                ui.selectable_value(&mut r.law, k, k.label());
                            }
                        });
                });
                draw_stoich_editor(ui, "Reactants", &mut r.reactants, &species_names, ri, 0);
                draw_stoich_editor(ui, "Products", &mut r.products, &species_names, ri, 1);
                draw_law_params(ui, r, &species_names, ri);
            });
    }
    if let Some(ri) = remove {
        p.reactions.remove(ri);
    }
    if ui.small_button("➕ Add reaction").clicked() {
        let name = next_name(
            "r",
            &p.reactions
                .iter()
                .map(|r| r.name.clone())
                .collect::<Vec<_>>(),
        );
        p.reactions.push(ReactionRow::blank(name));
    }
}

/// Draw the `(species, stoichiometry)` list editor for one side of a
/// reaction. `side` distinguishes reactants (0) from products (1) for
/// stable widget ids.
fn draw_stoich_editor(
    ui: &mut egui::Ui,
    label: &str,
    pairs: &mut Vec<(usize, f64)>,
    species_names: &[String],
    rxn: usize,
    side: usize,
) {
    ui.label(egui::RichText::new(label).strong().small());
    if species_names.is_empty() {
        ui.label(egui::RichText::new("add a species first").weak().small());
        return;
    }
    let mut remove: Option<usize> = None;
    for (idx, (sp, stoich)) in pairs.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            // Clamp a stale index left by a deleted species.
            if *sp >= species_names.len() {
                *sp = 0;
            }
            egui::ComboBox::from_id_source(format!("sysbio_st_{rxn}_{side}_{idx}"))
                .selected_text(&species_names[*sp])
                .width(80.0)
                .show_ui(ui, |ui| {
                    for (si, name) in species_names.iter().enumerate() {
                        ui.selectable_value(sp, si, name);
                    }
                });
            let lbl = ui.label("×");
            ui.add(egui::DragValue::new(stoich).speed(0.1).range(0.0..=100.0))
                .labelled_by(lbl.id);
            if ui.small_button("✖").clicked() {
                remove = Some(idx);
            }
        });
    }
    if let Some(idx) = remove {
        pairs.remove(idx);
    }
    if ui
        .small_button(format!("➕ add {}", label.to_ascii_lowercase()))
        .clicked()
    {
        pairs.push((0, 1.0));
    }
}

/// Draw the rate-law parameter widgets for one reaction. Which fields
/// show depends on the law family.
fn draw_law_params(ui: &mut egui::Ui, r: &mut ReactionRow, species_names: &[String], rxn: usize) {
    ui.label(egui::RichText::new("Parameters").strong().small());
    egui::Grid::new(format!("sysbio_law_params_{rxn}"))
        .num_columns(2)
        .spacing([8.0, 4.0])
        .show(ui, |ui| match r.law {
            LawKind::Constant => {
                let lbl = ui.label("flux");
                ui.add(egui::DragValue::new(&mut r.k).speed(0.1).range(0.0..=1.0e6))
                    .labelled_by(lbl.id);
                ui.end_row();
            }
            LawKind::MassAction => {
                let lbl = ui.label("k (rate constant)");
                ui.add(
                    egui::DragValue::new(&mut r.k)
                        .speed(0.01)
                        .range(0.0..=1.0e6),
                )
                .labelled_by(lbl.id);
                ui.end_row();
            }
            LawKind::MichaelisMenten => {
                let lbl = ui.label("Vmax");
                ui.add(egui::DragValue::new(&mut r.k).speed(0.1).range(0.0..=1.0e6))
                    .labelled_by(lbl.id);
                ui.end_row();
                let lbl = ui.label("Km");
                ui.add(
                    egui::DragValue::new(&mut r.param2)
                        .speed(0.1)
                        .range(1.0e-6..=1.0e6),
                )
                .labelled_by(lbl.id);
                ui.end_row();
            }
            LawKind::Hill => {
                let lbl = ui.label("Vmax");
                ui.add(egui::DragValue::new(&mut r.k).speed(0.1).range(0.0..=1.0e6))
                    .labelled_by(lbl.id);
                ui.end_row();
                let lbl = ui.label("Kd");
                ui.add(
                    egui::DragValue::new(&mut r.param2)
                        .speed(0.1)
                        .range(1.0e-6..=1.0e6),
                )
                .labelled_by(lbl.id);
                ui.end_row();
                let lbl = ui.label("Hill coefficient n");
                ui.add(
                    egui::DragValue::new(&mut r.param3)
                        .speed(0.1)
                        .range(0.1..=12.0),
                )
                .labelled_by(lbl.id);
                ui.end_row();
                ui.label("regulator");
                if !species_names.is_empty() {
                    if r.regulator >= species_names.len() {
                        r.regulator = 0;
                    }
                    egui::ComboBox::from_id_source(format!("sysbio_hill_reg_{rxn}"))
                        .selected_text(&species_names[r.regulator])
                        .show_ui(ui, |ui| {
                            for (si, name) in species_names.iter().enumerate() {
                                ui.selectable_value(&mut r.regulator, si, name);
                            }
                        });
                } else {
                    ui.label("(add a species)");
                }
                ui.end_row();
            }
        });
}

/// Run the selected engine against the assembled model.
fn run_engine(p: &mut SysbioPanel) {
    p.error = None;
    p.series.clear();
    p.result.clear();

    let model = match assemble_model(p) {
        Ok(m) => m,
        Err(e) => {
            p.error = Some(e);
            return;
        }
    };
    let names: Vec<String> = p.species.iter().map(|s| s.name.clone()).collect();

    match p.engine {
        Engine::Ode => run_ode(p, &model, &names),
        Engine::Gillespie => run_gillespie(p, &model, &names),
        Engine::Fba => run_fba(p, &model),
    }
}

/// Run the deterministic ODE time course and fill the plot + summary.
fn run_ode(p: &mut SysbioPanel, model: &Model, names: &[String]) {
    let tc = TimeCourse::new(p.t_end);
    match tc.run(model) {
        Ok(traj) => {
            for (i, name) in names.iter().enumerate() {
                let pts: Vec<[f64; 2]> = traj
                    .times
                    .iter()
                    .zip(traj.states.iter())
                    .map(|(t, s)| [*t, s.get(i).copied().unwrap_or(0.0)])
                    .collect();
                p.series.push((name.clone(), pts));
            }
            let fin = traj.final_state().unwrap_or(&[]);
            let mut out = format!(
                "ODE time course (RK45)\nsamples : {}\nt range : [0, {}]\n\n\
                 -- final state --\n",
                traj.len(),
                p.t_end,
            );
            for (i, name) in names.iter().enumerate() {
                out.push_str(&format!(
                    "  {name:<12} {:>12.4}\n",
                    fin.get(i).copied().unwrap_or(0.0)
                ));
            }
            p.result = out;
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

/// Run the Gillespie SSA and fill the plot + summary.
fn run_gillespie(p: &mut SysbioPanel, model: &Model, names: &[String]) {
    match StochasticModel::from_model(model) {
        Ok(sm) => match sm.gillespie(p.t_end, p.seed, 1_000_000) {
            Ok(trace) => {
                for (i, name) in names.iter().enumerate() {
                    let pts: Vec<[f64; 2]> = trace
                        .times
                        .iter()
                        .zip(trace.states.iter())
                        .map(|(t, s)| [*t, s.get(i).copied().unwrap_or(0) as f64])
                        .collect();
                    p.series.push((name.clone(), pts));
                }
                let fin = trace.final_state().unwrap_or(&[]);
                let mut out = format!(
                    "Gillespie SSA\nreaction events : {}\nt range : [0, {}]\n\n\
                     -- final counts --\n",
                    trace.len(),
                    p.t_end,
                );
                for (i, name) in names.iter().enumerate() {
                    out.push_str(&format!(
                        "  {name:<12} {:>10}\n",
                        fin.get(i).copied().unwrap_or(0)
                    ));
                }
                p.result = out;
            }
            Err(e) => p.error = Some(e.to_string()),
        },
        Err(e) => p.error = Some(e.to_string()),
    }
}

/// Run flux-balance analysis maximising the chosen reaction's flux.
fn run_fba(p: &mut SysbioPanel, model: &Model) {
    let objective = match p.reactions.get(p.fba_objective) {
        Some(r) => r.name.trim().to_string(),
        None => {
            p.error = Some("pick a reaction to maximise".to_string());
            return;
        }
    };
    match FbaProblem::from_model(model) {
        Ok(mut fba) => {
            if let Err(e) = fba.set_objective(&objective) {
                p.error = Some(format!("objective: {e}"));
                return;
            }
            match fba.optimize() {
                Ok(sol) => {
                    let mut out = format!(
                        "Flux-balance analysis\nobjective   : maximise flux({objective})\n\
                         feasible    : {}\nobjective Z : {:.5}\n\n-- reaction fluxes --\n",
                        sol.feasible, sol.objective_value,
                    );
                    for (r, flux) in p.reactions.iter().zip(sol.fluxes.iter()) {
                        out.push_str(&format!("  {:<14} {flux:>12.5}\n", r.name));
                    }
                    p.result = out;
                }
                Err(e) => p.error = Some(e.to_string()),
            }
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

/// After a species is removed, drop reaction stoichiometry entries
/// referencing it and shift higher indices down by one.
fn prune_species_refs(p: &mut SysbioPanel, removed: usize) {
    let fix = |pairs: &mut Vec<(usize, f64)>| {
        pairs.retain(|&(i, _)| i != removed);
        for (i, _) in pairs.iter_mut() {
            if *i > removed {
                *i -= 1;
            }
        }
    };
    for r in &mut p.reactions {
        fix(&mut r.reactants);
        fix(&mut r.products);
        if r.regulator == removed {
            r.regulator = 0;
        } else if r.regulator > removed {
            r.regulator -= 1;
        }
    }
}

/// Generate the next free `<prefix><n>` name not already in `existing`.
fn next_name(prefix: &str, existing: &[String]) -> String {
    for n in 1..10_000 {
        let candidate = format!("{prefix}{n}");
        if !existing.iter().any(|e| e == &candidate) {
            return candidate;
        }
    }
    format!("{prefix}_new")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_network_assembles_and_validates() {
        let p = SysbioPanel::default();
        let model = assemble_model(&p).expect("default network must assemble");
        assert_eq!(model.species.len(), 3);
        assert_eq!(model.reactions.len(), 3);
        assert!(model.validate().is_ok());
    }

    #[test]
    fn default_ode_runs_on_the_assembled_model() {
        let p = SysbioPanel::default();
        let model = assemble_model(&p).unwrap();
        let traj = TimeCourse::new(p.t_end).run(&model).unwrap();
        assert!(traj.len() > 1);
        assert_eq!(traj.states[0].len(), 3);
    }

    #[test]
    fn default_gillespie_runs_on_the_assembled_model() {
        let p = SysbioPanel::default();
        let model = assemble_model(&p).unwrap();
        let sm = StochasticModel::from_model(&model).unwrap();
        let trace = sm.gillespie(p.t_end, p.seed, 1_000_000).unwrap();
        assert!(!trace.times.is_empty());
    }

    #[test]
    fn editing_adds_a_species_and_reaction() {
        // Build a fresh two-species network from nothing.
        let mut p = SysbioPanel {
            species: vec![],
            reactions: vec![],
            ..SysbioPanel::default()
        };
        p.species.push(SpeciesRow {
            name: "X".to_string(),
            initial: 50.0,
        });
        p.species.push(SpeciesRow {
            name: "Y".to_string(),
            initial: 0.0,
        });
        let mut rxn = ReactionRow::blank("x_to_y");
        rxn.reactants = vec![(0, 1.0)];
        rxn.products = vec![(1, 1.0)];
        p.reactions.push(rxn);
        let model = assemble_model(&p).expect("user-built network must assemble");
        assert_eq!(model.species.len(), 2);
        assert_eq!(model.reactions.len(), 1);
    }

    #[test]
    fn assemble_rejects_duplicate_species_names() {
        let mut p = SysbioPanel::default();
        p.species[1].name = "A".to_string(); // collide with species 0
        assert!(assemble_model(&p).is_err());
    }

    #[test]
    fn assemble_rejects_empty_species_name() {
        let mut p = SysbioPanel::default();
        p.species[0].name = "  ".to_string();
        assert!(assemble_model(&p).is_err());
    }

    #[test]
    fn michaelis_menten_law_needs_a_substrate() {
        // An MM reaction with no reactant cannot pick a substrate.
        let mut r = ReactionRow::blank("mm");
        r.law = LawKind::MichaelisMenten;
        assert!(build_rate_law(&r).is_err());
        r.reactants = vec![(0, 1.0)];
        assert!(build_rate_law(&r).is_ok());
    }

    #[test]
    fn hill_law_builds_with_its_parameters() {
        let mut r = ReactionRow::blank("hill");
        r.law = LawKind::Hill;
        r.param3 = 4.0;
        match build_rate_law(&r).unwrap() {
            RateLaw::Hill { n, repress, .. } => {
                assert_eq!(n, 4.0);
                assert!(!repress);
            }
            other => panic!("expected Hill, got {other:?}"),
        }
    }

    #[test]
    fn prune_species_refs_drops_and_shifts_indices() {
        let mut p = SysbioPanel::default();
        // Reaction b_to_c references species 1 (B) and 2 (C).
        // Removing species 0 (A) must shift those down to 0 and 1.
        p.species.remove(0);
        prune_species_refs(&mut p, 0);
        let b_to_c = &p.reactions[2];
        assert_eq!(b_to_c.reactants, vec![(0, 1.0)]); // was (1, ..)
        assert_eq!(b_to_c.products, vec![(1, 1.0)]); // was (2, ..)
                                                     // The influx reaction's product (was species 0) is dropped.
        assert!(p.reactions[0].products.is_empty());
    }

    #[test]
    fn next_name_skips_taken_names() {
        let taken = vec!["S1".to_string(), "S2".to_string()];
        assert_eq!(next_name("S", &taken), "S3");
        assert_eq!(next_name("r", &[]), "r1");
    }

    #[test]
    fn fba_runs_on_the_default_network() {
        let p = SysbioPanel {
            fba_objective: 2, // b_to_c
            ..SysbioPanel::default()
        };
        let model = assemble_model(&p).unwrap();
        let mut fba = FbaProblem::from_model(&model).unwrap();
        assert!(fba.set_objective("b_to_c").is_ok());
        // optimize() must not panic; feasibility is model-dependent.
        let _ = fba.optimize();
    }
}

/// Headless egui UI-logic tests for the Systems Biology panel.
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
        app.genetics.active = GeneticsPanel::SystemsBiology;
        app
    }

    #[test]
    fn draws_every_engine_without_panic() {
        // The species + reaction editors and the per-engine controls
        // must all render without panicking.
        for engine in [Engine::Ode, Engine::Gillespie, Engine::Fba] {
            let mut app = app_with_panel();
            app.genetics.sysbio.engine = engine;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_an_empty_network_without_panic() {
        // A network with no species / reactions exercises the
        // "add a species first" empty-state branches.
        let mut app = app_with_panel();
        app.genetics.sysbio.species.clear();
        app.genetics.sysbio.reactions.clear();
        draw_headless(&mut app);
    }

    #[test]
    fn draws_post_run_state_with_series_without_panic() {
        // Post-run: a time-course series drives an egui_plot line chart.
        let mut app = app_with_panel();
        app.genetics.sysbio.series = vec![
            ("A".to_string(), vec![[0.0, 100.0], [1.0, 70.0]]),
            ("B".to_string(), vec![[0.0, 0.0], [1.0, 25.0]]),
        ];
        app.genetics.sysbio.result = "ODE time course (RK45)\n".to_string();
        draw_headless(&mut app);
        // Error state.
        let mut app = app_with_panel();
        app.genetics.sysbio.error = Some("duplicate species name".to_string());
        draw_headless(&mut app);
    }

    #[test]
    fn run_engine_ode_produces_a_time_course() {
        // The ODE engine calls the real valenx-sysbio TimeCourse API
        // and fills the plot series + summary.
        let mut p = SysbioPanel {
            engine: Engine::Ode,
            ..SysbioPanel::default()
        };
        run_engine(&mut p);
        assert!(p.error.is_none(), "ODE run errored: {:?}", p.error);
        assert!(!p.series.is_empty(), "ODE produced no series");
        assert!(p.result.contains("ODE time course"));
    }

    #[test]
    fn run_engine_gillespie_produces_a_trace() {
        let mut p = SysbioPanel {
            engine: Engine::Gillespie,
            ..SysbioPanel::default()
        };
        run_engine(&mut p);
        assert!(p.error.is_none(), "Gillespie run errored: {:?}", p.error);
        assert!(p.result.contains("Gillespie SSA"));
    }

    #[test]
    fn run_engine_fba_runs_without_panic() {
        // FBA runs end-to-end; feasibility is model-dependent, but the
        // run must never panic and must leave a well-formed state.
        let mut p = SysbioPanel {
            engine: Engine::Fba,
            fba_objective: 2,
            ..SysbioPanel::default()
        };
        run_engine(&mut p);
        assert!(p.error.is_some() || !p.result.is_empty());
    }

    #[test]
    fn run_engine_surfaces_error_on_a_malformed_network() {
        // A network with no species cannot be assembled — the panel
        // must surface an error rather than panicking.
        let mut p = SysbioPanel {
            species: vec![],
            reactions: vec![],
            ..SysbioPanel::default()
        };
        run_engine(&mut p);
        assert!(p.error.is_some(), "run should error on an empty network");
        assert!(p.series.is_empty());
        // A duplicate species name is also malformed.
        let mut p = SysbioPanel::default();
        p.species[1].name = "A".to_string(); // collide with species 0
        run_engine(&mut p);
        assert!(p.error.is_some(), "run should error on duplicate species");
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        use egui::accesskit::Role;
        // ODE mode exposes t_end; Gillespie additionally exposes seed.
        for engine in [Engine::Ode, Engine::Gillespie] {
            let mut app = app_with_panel();
            app.genetics.sysbio.engine = engine;
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
                "sysbio {engine:?} should expose at least one SpinButton"
            );
            assert!(
                spin_buttons
                    .iter()
                    .all(|(_, n)| !n.labelled_by().is_empty()),
                "every sysbio DragValue ({engine:?}) must be labelled_by its caption"
            );
        }
    }
}
