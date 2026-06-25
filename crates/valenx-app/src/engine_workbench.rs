//! The right-side **Engine** workbench panel — a reactive
//! **design → analyze → optimize → export** loop over a liquid-rocket engine.
//!
//! Mirrors the other single-file workbenches (`rocket_workbench`,
//! `fem_workbench`, …): a resizable [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_engine_workbench`, toggled from the View menu.
//!
//! Editing any design knob (chamber pressure, temperature, gas γ / molar
//! mass, throat area, expansion ratio) re-runs the ideal-nozzle performance
//! ([`valenx_astro::EngineDesign`]) and the first-order Bartz regen-cooling
//! balance live, with a SAFE / OVER-FLUX verdict against a target cooling
//! margin. A one-click optimizer searches chamber-pressure × expansion-ratio
//! for the highest sea-level Isp that still clears the cooling margin, the
//! 3-D nozzle contour renders in the central viewport, and the contour is
//! exportable to STL.
//!
//! Honest scope: first-order **preliminary-design** physics — ideal-nozzle
//! thermodynamics + a Bartz heat-flux estimate + a Rao parabolic-approximation
//! nozzle contour. Not a combustion-CFD / conjugate-heat-transfer /
//! generative-geometry tool.

use std::path::PathBuf;

use eframe::egui;

use crate::types::LoadedMesh;
use crate::ValenxApp;
use valenx_astro::{
    combust, optimize_engine, solve_cycle, CoolingInputs, CoolingPerformance, CycleInputs,
    EngineDesign, EngineOptimum, EnginePerformance, Propellant,
};
use valenx_rocket_demo::nozzle::nozzle_mesh;

/// Ideal-nozzle performance + cooling balance for the current design.
struct EngineReport {
    vacuum: EnginePerformance,
    sea_level: EnginePerformance,
    cooling: CoolingPerformance,
}

/// Persistent form + result state for the Engine workbench.
pub struct EngineWorkbenchState {
    /// The engine design knobs (a kerolox-class default).
    design: EngineDesign,
    /// Cooling-circuit inputs (RP-1 regen default).
    cooling: CoolingInputs,
    /// Required cooling margin for the verdict + the optimizer constraint.
    target_margin: f64,
    /// Last optimizer result, if any.
    opt: Option<EngineOptimum>,
    /// Deferred request to load the 3-D nozzle into the central viewport.
    show_3d_request: bool,
    /// Deferred request to load the detailed 3-D engine (chamber + regen
    /// nozzle + powerhead) into the central viewport.
    show_engine_3d_request: bool,
    /// Last STL-export outcome message (path on success, error otherwise).
    last_export: Option<String>,
    /// Propellant for the equilibrium-combustion chamber prediction.
    propellant: Propellant,
    /// Oxidizer/fuel mass mixture ratio for the combustion prediction.
    mixture_ratio: f64,
}

impl Default for EngineWorkbenchState {
    fn default() -> Self {
        Self {
            // A kerolox gas-generator-class design point.
            design: EngineDesign {
                chamber_pressure: 9.7e6,
                chamber_temperature: 3_500.0,
                gamma: 1.2,
                molar_mass: 22.0,
                throat_area: 0.05,
                expansion_ratio: 16.0,
            },
            cooling: CoolingInputs::default(),
            target_margin: 1.5,
            opt: None,
            show_3d_request: false,
            show_engine_3d_request: false,
            last_export: None,
            propellant: Propellant::Ch4Lox,
            mixture_ratio: 3.6,
        }
    }
}

impl EngineWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`
    /// so an agent can discover the name space.
    ///
    /// These are the **settable design / chemistry knobs** — they match exactly
    /// what the workbench body draws (and what each control is `labelled_by` /
    /// captioned via a `Slider`/`ComboBox`). Output-only readouts (performance,
    /// Bartz cooling, the staged-combustion cycle), the optimizer, and the
    /// 3-D-geometry / STL-export buttons are actions, not settable values, and
    /// are intentionally **not** listed here.
    ///
    /// Two captions carry display-unit conversions handled by [`agent_set`]:
    /// `"chamber pressure"` is set in **bar** (stored internally in Pa) and
    /// `"throat area"` is set in **cm²** (stored internally in m²), matching the
    /// `DragValue` units a user edits.
    ///
    /// [`agent_set`]: Self::agent_set
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "chamber pressure",
            "chamber temp",
            "gas γ",
            "molar mass",
            "throat area",
            "expansion ratio ε",
            "target cooling margin",
            "propellant",
            "mixture ratio (O/F)",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. The caption strings match exactly what the workbench
    /// body draws, so an agent can set a parameter by the same name a user reads.
    ///
    /// Fail-loud: an unknown caption or a value of the wrong type returns
    /// `Err(String)` (the bridge turns it into a `warn` feed note) — never a
    /// panic, and no field is written on error. Numeric captions read
    /// [`AgentValue::as_f64`]; the one enum caption (`propellant`) reads
    /// [`AgentValue::as_str`] and matches a small set of combo names.
    ///
    /// Unit handling mirrors the UI: `"chamber pressure"` is given in **bar**
    /// (the field stores Pa, so it is multiplied by `1e5`) and `"throat area"`
    /// is given in **cm²** (the field stores m², so it is multiplied by `1e-4`).
    /// Output-only quantities (thrust / Isp / heat flux / cycle), the optimizer,
    /// and the 3-D / export buttons are not settable and are rejected as unknown.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        match name {
            // -- Design knobs (numeric DragValues) --
            // Caption is in bar; the field stores Pa (UI multiplies by 1e5).
            "chamber pressure" => self.design.chamber_pressure = value.as_f64()? * 1.0e5,
            "chamber temp" => self.design.chamber_temperature = value.as_f64()?,
            "gas γ" => self.design.gamma = value.as_f64()?,
            "molar mass" => self.design.molar_mass = value.as_f64()?,
            // Caption is in cm²; the field stores m² (UI multiplies by 1e-4).
            "throat area" => self.design.throat_area = value.as_f64()? * 1.0e-4,
            "expansion ratio ε" => self.design.expansion_ratio = value.as_f64()?,
            // -- Verdict / optimizer constraint (Slider) --
            "target cooling margin" => self.target_margin = value.as_f64()?,
            // -- Combustion chemistry --
            "propellant" => self.propellant = parse_propellant(value.as_str()?)?,
            "mixture ratio (O/F)" => self.mixture_ratio = value.as_f64()?,
            other => return Err(format!("unknown Engine control: {other:?}")),
        }
        Ok(())
    }
}

/// Human-readable label for a propellant in the picker.
fn propellant_label(p: Propellant) -> &'static str {
    match p {
        Propellant::H2Lox => "H₂ / LOX (hydrolox)",
        Propellant::Rp1Lox => "RP-1 / LOX (kerolox)",
        Propellant::Ch4Lox => "CH₄ / LOX (methalox)",
    }
}

/// Parse a propellant-combo name (for the agent `SetControl` bridge) into a
/// [`Propellant`]. Case-insensitive and tolerant of the common shorthand
/// (the chemistry, the slang, and the picker words) so an agent can pick the
/// combo by the same name a user reads. Fail-loud on an unrecognised name so
/// a typo becomes a `warn` note, not a silent no-op.
fn parse_propellant(s: &str) -> Result<Propellant, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "h2lox" | "h2/lox" | "h2 / lox" | "hydrolox" | "lh2/lox" | "h2" => Ok(Propellant::H2Lox),
        "rp1lox" | "rp1/lox" | "rp-1/lox" | "rp-1 / lox" | "kerolox" | "rp1" | "rp-1" => {
            Ok(Propellant::Rp1Lox)
        }
        "ch4lox" | "ch4/lox" | "ch4 / lox" | "methalox" | "ch4" | "methane" => {
            Ok(Propellant::Ch4Lox)
        }
        other => Err(format!(
            "unknown propellant '{other}' (expected 'hydrolox', 'kerolox', or 'methalox')"
        )),
    }
}

/// Ideal-nozzle performance + Bartz cooling for `design`; `Err` (as a
/// human-readable string) for a non-physical design point.
fn analyze(design: &EngineDesign, cooling: &CoolingInputs) -> Result<EngineReport, String> {
    Ok(EngineReport {
        vacuum: design.vacuum().map_err(|e| e.to_string())?,
        sea_level: design.sea_level().map_err(|e| e.to_string())?,
        cooling: design.cooling(cooling).map_err(|e| e.to_string())?,
    })
}

/// Build the 3-D nozzle mesh for the current design and load it into the
/// central viewport (replacing any current STL / mesh) so it can be orbited.
fn load_nozzle_3d(app: &mut ValenxApp) {
    let d = app.engine.design;
    let mesh = nozzle_mesh(d.throat_area, d.expansion_ratio, 48);
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<engine>/nozzle"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Build the procedurally-detailed 3-D engine — combustion chamber, fluted
/// regen-cooled nozzle, injector dome and the full powerhead (twin
/// turbopumps, twin preburners, hot-gas manifold + plumbing) — and load it
/// into the central viewport so it can be orbited in 3-D.
fn load_engine_3d(app: &mut ValenxApp) {
    let mesh = crate::rocket_mesh::detailed_engine_mesh();
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<engine>/detailed"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Draw the Engine workbench right-side panel. A no-op when the
/// `show_engine_workbench` toggle is off.
pub fn draw_engine_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_engine_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_engine_workbench",
        "Engine — design → analyze",
        engine_workbench_body,
    );
    if close {
        app.show_engine_workbench = false;
    }

    drain_deferred(app);
}

/// Drain the Engine workbench's deferred requests (outside any panel
/// borrow): load the 3-D nozzle / detailed powerhead into the central
/// viewport when the body set the request flag. Called by both
/// [`draw_engine_workbench`] and the dockable layout
/// ([`crate::dock_layout`]) so the 3-D buttons work in either host.
pub(crate) fn drain_deferred(app: &mut ValenxApp) {
    if app.engine.show_3d_request {
        app.engine.show_3d_request = false;
        load_nozzle_3d(app);
    }
    if app.engine.show_engine_3d_request {
        app.engine.show_engine_3d_request = false;
        load_engine_3d(app);
    }
}

/// The Engine workbench body — design knobs, live analysis, the staged-
/// combustion cycle readout, the optimizer, and the 3-D / export controls.
/// Extracted from [`draw_engine_workbench`] so it can be hosted *either* by
/// the classic right-side [`crate::workbench_chrome::workbench_shell`] *or*
/// by the opt-in dockable tile layout ([`crate::dock_layout`]) — no
/// duplicated logic.
pub(crate) fn engine_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new("ideal-nozzle performance + Bartz regen-cooling · valenx-astro")
            .weak()
            .small(),
    );
    ui.separator();

    let s = &mut app.engine;
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            // ── Design knobs ─────────────────────────────────────
            ui.label(egui::RichText::new("Design point").strong());
            ui.horizontal(|ui| {
                let lbl = ui.label("chamber pressure");
                let mut bar = s.design.chamber_pressure / 1.0e5;
                if ui
                    .add(
                        egui::DragValue::new(&mut bar)
                            .speed(1.0)
                            .range(1.0..=400.0)
                            .suffix(" bar"),
                    )
                    .labelled_by(lbl.id)
                    .changed()
                {
                    s.design.chamber_pressure = bar * 1.0e5;
                }
            });
            ui.horizontal(|ui| {
                let lbl = ui.label("chamber temp");
                ui.add(
                    egui::DragValue::new(&mut s.design.chamber_temperature)
                        .speed(10.0)
                        .range(1_000.0..=4_500.0)
                        .suffix(" K"),
                )
                .labelled_by(lbl.id);
            });
            ui.horizontal(|ui| {
                let lbl = ui.label("gas γ");
                ui.add(
                    egui::DragValue::new(&mut s.design.gamma)
                        .speed(0.005)
                        .range(1.05..=1.40),
                )
                .labelled_by(lbl.id);
            });
            ui.horizontal(|ui| {
                let lbl = ui.label("molar mass");
                ui.add(
                    egui::DragValue::new(&mut s.design.molar_mass)
                        .speed(0.2)
                        .range(5.0..=40.0)
                        .suffix(" g/mol"),
                )
                .labelled_by(lbl.id);
            });
            ui.horizontal(|ui| {
                let lbl = ui.label("throat area");
                let mut cm2 = s.design.throat_area * 1.0e4;
                if ui
                    .add(
                        egui::DragValue::new(&mut cm2)
                            .speed(1.0)
                            .range(1.0..=5_000.0)
                            .suffix(" cm²"),
                    )
                    .labelled_by(lbl.id)
                    .changed()
                {
                    s.design.throat_area = cm2 * 1.0e-4;
                }
            });
            ui.horizontal(|ui| {
                let lbl = ui.label("expansion ratio ε");
                ui.add(
                    egui::DragValue::new(&mut s.design.expansion_ratio)
                        .speed(0.5)
                        .range(1.0..=200.0),
                )
                .labelled_by(lbl.id);
            });
            // Use the Slider's own `.text(..)` caption (not a separate
            // `ui.label`) so egui links BOTH the slider node and its internal
            // value DragValue (Role::SpinButton) to the caption — an external
            // `.labelled_by` would only name the slider, leaving the embedded
            // SpinButton unnamed.
            ui.add(
                egui::Slider::new(&mut s.target_margin, 1.0..=4.0).text("target cooling margin"),
            );

            // ── Combustion chemistry (predict chamber conditions) ─
            ui.add_space(4.0);
            ui.separator();
            ui.label(egui::RichText::new("Combustion — equilibrium chamber prediction").strong());
            egui::ComboBox::from_label("propellant")
                .selected_text(propellant_label(s.propellant))
                .show_ui(ui, |ui| {
                    for p in [Propellant::H2Lox, Propellant::Rp1Lox, Propellant::Ch4Lox] {
                        ui.selectable_value(&mut s.propellant, p, propellant_label(p));
                    }
                });
            ui.add(egui::Slider::new(&mut s.mixture_ratio, 1.5..=8.0).text("mixture ratio (O/F)"));
            let comb = combust(
                s.propellant,
                s.mixture_ratio,
                s.design.chamber_pressure / 1.0e5,
            );
            ui.label(
                egui::RichText::new(format!(
                    "  chamber T  : {:.0} K\n  \
                             gas γ      : {:.3}\n  \
                             molar mass : {:.1} g/mol\n  \
                             c*         : {:.0} m/s",
                    comb.chamber_temperature, comb.gamma, comb.molar_mass, comb.c_star,
                ))
                .monospace()
                .small(),
            );
            if ui
                .button("Apply combustion → engine (sets T, γ, molar mass)")
                .clicked()
            {
                s.design.chamber_temperature = comb.chamber_temperature;
                s.design.gamma = comb.gamma;
                s.design.molar_mass = comb.molar_mass;
            }
            ui.label(
                egui::RichText::new(
                    "first-order equilibrium thermochem — H₂/LOX is validated; \
                             RP-1 / CH₄ run ~10% high vs NASA CEA.",
                )
                .weak()
                .small(),
            );

            ui.add_space(4.0);
            ui.separator();

            // ── Live analysis ────────────────────────────────────
            match analyze(&s.design, &s.cooling) {
                Ok(r) => {
                    ui.label(egui::RichText::new("Performance — ideal nozzle").strong());
                    ui.label(
                        egui::RichText::new(format!(
                            "  vac thrust : {:.0} kN\n  \
                                     vac Isp    : {:.0} s\n  \
                                     SL  thrust : {:.0} kN\n  \
                                     SL  Isp    : {:.0} s\n  \
                                     c*         : {:.0} m/s\n  \
                                     exit Mach  : {:.2}",
                            r.vacuum.thrust / 1.0e3,
                            r.vacuum.isp,
                            r.sea_level.thrust / 1.0e3,
                            r.sea_level.isp,
                            r.vacuum.c_star,
                            r.vacuum.exit_mach,
                        ))
                        .monospace()
                        .small(),
                    );

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Cooling — Bartz regen balance").strong());
                    ui.label(
                        egui::RichText::new(format!(
                            "  throat flux : {:.1} MW/m²\n  \
                                     coolant ΔT  : {:.0} K\n  \
                                     margin      : {:.2}  (target {:.2})",
                            r.cooling.throat_heat_flux / 1.0e6,
                            r.cooling.coolant_temperature_rise,
                            r.cooling.cooling_margin,
                            s.target_margin,
                        ))
                        .monospace()
                        .small(),
                    );
                    let ok = r.cooling.cooling_margin >= s.target_margin;
                    let (txt, col) = if ok {
                        (
                            format!(
                                "✔ COOLED  ·  {:.0}% margin over target",
                                (r.cooling.cooling_margin / s.target_margin - 1.0) * 100.0
                            ),
                            egui::Color32::from_rgb(80, 220, 120),
                        )
                    } else {
                        (
                            "✖ OVER-FLUX  ·  raise ε / lower chamber pressure".to_string(),
                            egui::Color32::from_rgb(220, 90, 90),
                        )
                    };
                    ui.add_space(2.0);
                    ui.colored_label(col, egui::RichText::new(txt).strong());
                }
                Err(e) => {
                    ui.colored_label(
                        egui::Color32::from_rgb(220, 90, 90),
                        format!("invalid design: {e}"),
                    );
                }
            }

            // ── Staged-combustion cycle (full-flow power balance) ─
            ui.add_space(6.0);
            ui.separator();
            ui.label(
                egui::RichText::new("Staged-combustion cycle — full-flow power balance").strong(),
            );
            ui.label(
                egui::RichText::new(
                    "can twin turbopumps drive this chamber pressure? \
                             (Raptor-class methalox FFSC reference)",
                )
                .weak()
                .small(),
            );
            let mut ci = CycleInputs::raptor_methalox();
            ci.chamber_pressure = s.design.chamber_pressure;
            let cyc = solve_cycle(&ci);
            ui.label(
                egui::RichText::new(format!(
                    "  max chamber : {:.0} bar\n  \
                             ox turbine  : {:.1} MW @ {:.0} K\n  \
                             fuel turbine: {:.1} MW @ {:.0} K",
                    cyc.max_chamber_pressure / 1.0e5,
                    cyc.ox.turbine_power / 1.0e6,
                    ci.ox.turbine_inlet_temperature,
                    cyc.fuel.turbine_power / 1.0e6,
                    ci.fuel.turbine_inlet_temperature,
                ))
                .monospace()
                .small(),
            );
            let (cyc_txt, cyc_col) = if cyc.closes {
                (
                    format!(
                        "✔ CYCLE CLOSES at {:.0} bar",
                        s.design.chamber_pressure / 1.0e5
                    ),
                    egui::Color32::from_rgb(80, 220, 120),
                )
            } else {
                (
                    format!(
                        "✖ WON'T CLOSE · turbopumps top out at ~{:.0} bar",
                        cyc.max_chamber_pressure / 1.0e5
                    ),
                    egui::Color32::from_rgb(220, 90, 90),
                )
            };
            ui.colored_label(cyc_col, egui::RichText::new(cyc_txt).strong());

            // ── Optimizer ────────────────────────────────────────
            ui.add_space(6.0);
            ui.separator();
            ui.label(egui::RichText::new("AI optimizer — max sea-level Isp vs cooling").strong());
            ui.label(
                egui::RichText::new(
                    "searches chamber pressure × expansion ratio for the highest \
                             sea-level Isp whose cooling margin clears the target.",
                )
                .weak()
                .small(),
            );
            if ui
                .button(egui::RichText::new("Run engine optimizer").strong())
                .clicked()
            {
                s.opt = optimize_engine(
                    &s.design,
                    &s.cooling,
                    s.target_margin,
                    (1.0e6, 30.0e6),
                    (2.0, 120.0),
                    28,
                );
            }
            if let Some(o) = s.opt {
                ui.label(
                    egui::RichText::new(format!(
                        "best: {:.0} bar · ε {:.1} → SL Isp {:.0} s\n\
                                 margin {:.2} · ran {} sims ({} feasible)",
                        o.design.chamber_pressure / 1.0e5,
                        o.design.expansion_ratio,
                        o.sea_level.isp,
                        o.cooling.cooling_margin,
                        o.evaluations,
                        o.feasible_count,
                    ))
                    .monospace()
                    .small(),
                );
                if ui.button("Apply optimized design").clicked() {
                    s.design.chamber_pressure = o.design.chamber_pressure;
                    s.design.expansion_ratio = o.design.expansion_ratio;
                }
            }

            // ── 3-D geometry + export ────────────────────────────
            ui.add_space(6.0);
            ui.separator();
            ui.label(egui::RichText::new("Engine geometry").strong());
            if ui
                .button(egui::RichText::new("Show 3-D engine (powerhead)").strong())
                .clicked()
            {
                s.show_engine_3d_request = true;
            }
            ui.horizontal(|ui| {
                if ui.button("Show 3-D nozzle").clicked() {
                    s.show_3d_request = true;
                }
                if ui.button("Export nozzle STL").clicked() {
                    let mesh = nozzle_mesh(s.design.throat_area, s.design.expansion_ratio, 64);
                    let path = std::env::temp_dir().join("valenx_nozzle.stl");
                    s.last_export = Some(match valenx_mesh::write_stl_binary(&mesh, &path) {
                        Ok(()) => format!("exported → {}", path.display()),
                        Err(e) => format!("export error: {e}"),
                    });
                }
            });
            if let Some(msg) = &s.last_export {
                ui.label(egui::RichText::new(msg).weak().small());
            }
        });
}

/// The agent-bridge **`show_3d{kind:"engine"}`** product: the canonical
/// procedurally-detailed liquid-rocket engine (chamber + regen nozzle +
/// powerhead — the panel's primary "Show 3-D engine" geometry) paired with
/// the workbench's own ideal-nozzle + Bartz-cooling headline numbers, at a
/// fixed 3/4 camera. Registered in [`crate::products_registry`]; the per-tool
/// builder the registry dispatches to. Pure — driven off
/// [`EngineWorkbenchState::default`].
///
/// The readout rows mirror the panel's `analyze()` "Performance — ideal
/// nozzle" + "Cooling — Bartz regen balance" sections (this workbench formats
/// its readout inline rather than via a shared `compute()` string fn, so the
/// rows are assembled here from the same [`EngineReport`]).
pub(crate) fn engine_product() -> crate::WorkspaceProduct {
    let s = EngineWorkbenchState::default();
    let mesh = crate::rocket_mesh::detailed_engine_mesh();
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<engine>/detailed");
    let r =
        analyze(&s.design, &s.cooling).expect("canonical kerolox engine ⇒ analysis is physical");
    let readout = format!(
        "vac thrust : {:.0} kN\n\
         vac Isp    : {:.0} s\n\
         SL  thrust : {:.0} kN\n\
         SL  Isp    : {:.0} s\n\
         c*         : {:.0} m/s\n\
         exit Mach  : {:.2}\n\
         throat flux: {:.1} MW/m²\n\
         coolant ΔT : {:.0} K\n\
         margin     : {:.2}  (target {:.2})",
        r.vacuum.thrust / 1.0e3,
        r.vacuum.isp,
        r.sea_level.thrust / 1.0e3,
        r.sea_level.isp,
        r.vacuum.c_star,
        r.vacuum.exit_mach,
        r.cooling.throat_heat_flux / 1.0e6,
        r.cooling.coolant_temperature_rise,
        r.cooling.cooling_margin,
        s.target_margin,
    );
    let lines = crate::products_registry::lines_from_readout(&readout);
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Engine (liquid-rocket)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        animation: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_design_is_a_valid_kerolox_engine() {
        let s = EngineWorkbenchState::default();
        let r = analyze(&s.design, &s.cooling).expect("default design is physical");
        // Kerolox-class numbers: positive thrust, vacuum Isp beats sea level,
        // Isp in a sane band.
        assert!(r.vacuum.thrust > 0.0 && r.sea_level.thrust > 0.0);
        assert!(r.vacuum.isp > r.sea_level.isp);
        assert!(
            (250.0..380.0).contains(&r.vacuum.isp),
            "Isp {}",
            r.vacuum.isp
        );
        assert!(r.cooling.cooling_margin > 0.0 && r.cooling.cooling_margin.is_finite());
        assert!((s.target_margin - 1.5).abs() < 1e-12);
        assert!(s.opt.is_none());
    }

    #[test]
    fn analyze_errors_on_invalid_design() {
        let mut d = EngineWorkbenchState::default().design;
        d.gamma = 1.0; // must be > 1
        assert!(analyze(&d, &CoolingInputs::default()).is_err());
    }

    #[test]
    fn optimizer_wires_to_a_feasible_result() {
        let s = EngineWorkbenchState::default();
        let o = optimize_engine(
            &s.design,
            &s.cooling,
            s.target_margin,
            (1.0e6, 30.0e6),
            (2.0, 120.0),
            24,
        )
        .expect("a feasible engine exists");
        assert!(o.cooling.cooling_margin >= s.target_margin - 1e-9);
        assert!(o.sea_level.isp > 0.0);
    }

    #[test]
    fn combustion_prediction_is_physical_and_applies_cleanly() {
        let mut s = EngineWorkbenchState::default();
        assert_eq!(s.propellant, Propellant::Ch4Lox);
        let comb = combust(
            s.propellant,
            s.mixture_ratio,
            s.design.chamber_pressure / 1.0e5,
        );
        assert!(
            (2_500.0..4_500.0).contains(&comb.chamber_temperature),
            "Tc {}",
            comb.chamber_temperature
        );
        // Applying the prediction into the design keeps it physical.
        s.design.chamber_temperature = comb.chamber_temperature;
        s.design.gamma = comb.gamma;
        s.design.molar_mass = comb.molar_mass;
        assert!(analyze(&s.design, &s.cooling).is_ok());
    }

    #[test]
    fn cycle_readout_closes_at_a_modest_chamber_pressure() {
        // The default 97-bar design is well under the Raptor-class ceiling,
        // so the full-flow cycle must close with a healthy max-Pc headroom.
        let s = EngineWorkbenchState::default();
        let mut ci = CycleInputs::raptor_methalox();
        ci.chamber_pressure = s.design.chamber_pressure;
        let cyc = solve_cycle(&ci);
        assert!(cyc.closes, "97 bar should close on Raptor-class hardware");
        assert!(cyc.max_chamber_pressure / 1.0e5 > 250.0);
    }

    // ---- agent bridge (SetControl) ----

    #[test]
    fn agent_set_sets_representative_params_verified_via_state() {
        use crate::agent_commands::AgentValue;
        let mut s = EngineWorkbenchState::default();

        // Numeric, with the bar→Pa display-unit conversion.
        s.agent_set("chamber pressure", &AgentValue::Float(250.0))
            .expect("set chamber pressure");
        assert!(
            (s.design.chamber_pressure - 250.0 * 1.0e5).abs() < 1e-3,
            "250 bar must store as 25 MPa, got {} Pa",
            s.design.chamber_pressure
        );

        // Numeric, with the cm²→m² display-unit conversion.
        s.agent_set("throat area", &AgentValue::Float(80.0))
            .expect("set throat area");
        assert!(
            (s.design.throat_area - 80.0 * 1.0e-4).abs() < 1e-12,
            "80 cm² must store as 0.008 m², got {} m²",
            s.design.throat_area
        );

        // Plain numeric (no conversion).
        s.agent_set("expansion ratio ε", &AgentValue::Float(40.0))
            .expect("set expansion ratio");
        assert!((s.design.expansion_ratio - 40.0).abs() < 1e-12);

        s.agent_set("target cooling margin", &AgentValue::Float(2.5))
            .expect("set target margin");
        assert!((s.target_margin - 2.5).abs() < 1e-12);

        // Enum combo by name.
        s.agent_set("propellant", &AgentValue::Str("hydrolox".into()))
            .expect("set propellant");
        assert_eq!(s.propellant, Propellant::H2Lox);
        s.agent_set("propellant", &AgentValue::Str("RP-1 / LOX".into()))
            .expect("set propellant kerolox");
        assert_eq!(s.propellant, Propellant::Rp1Lox);

        s.agent_set("mixture ratio (O/F)", &AgentValue::Float(2.7))
            .expect("set mixture ratio");
        assert!((s.mixture_ratio - 2.7).abs() < 1e-12);
    }

    #[test]
    fn agent_set_unknown_control_is_err_not_panic() {
        use crate::agent_commands::AgentValue;
        let mut s = EngineWorkbenchState::default();
        assert!(
            s.agent_set("no such knob", &AgentValue::Float(1.0))
                .is_err(),
            "an unknown caption must return Err, not panic"
        );
        // Output-only / action captions are deliberately NOT settable.
        assert!(s.agent_set("vac Isp", &AgentValue::Float(330.0)).is_err());
    }

    #[test]
    fn agent_set_type_mismatch_is_err_not_panic() {
        use crate::agent_commands::AgentValue;
        let mut s = EngineWorkbenchState::default();
        // A numeric field given a string must fail-loud (no write).
        let before = s.design.chamber_pressure;
        assert!(
            s.agent_set("chamber pressure", &AgentValue::Str("lots".into()))
                .is_err(),
            "a string for a numeric field must return Err"
        );
        assert_eq!(
            s.design.chamber_pressure, before,
            "no field is written on a type-mismatch error"
        );
        // An enum field given a number must fail-loud.
        assert!(
            s.agent_set("propellant", &AgentValue::Float(3.0)).is_err(),
            "a number for the propellant enum must return Err"
        );
        // An unrecognised enum name must fail-loud.
        assert!(
            s.agent_set("propellant", &AgentValue::Str("plasma".into()))
                .is_err(),
            "an unknown propellant name must return Err"
        );
    }

    #[test]
    fn agent_control_names_are_all_settable() {
        // Every advertised caption must be accepted by agent_set (with an
        // appropriately-typed value) — no caption is listed but unhandled.
        use crate::agent_commands::AgentValue;
        for &name in EngineWorkbenchState::agent_control_names() {
            let mut s = EngineWorkbenchState::default();
            let v = if name == "propellant" {
                AgentValue::Str("methalox".into())
            } else {
                AgentValue::Float(2.0)
            };
            assert!(
                s.agent_set(name, &v).is_ok(),
                "advertised control '{name}' must be settable via agent_set"
            );
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_engine_workbench(app, ctx);
        });
    }

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_engine_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_engine_workbench);
        draw(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_engine_workbench = true;
        draw(&mut app);
    }

    #[test]
    fn workbench_draws_invalid_design_without_panic() {
        let mut app = ValenxApp::default();
        app.show_engine_workbench = true;
        app.engine.design.gamma = 1.0; // invalid → error branch must render
        draw(&mut app);
    }

    #[test]
    fn show_3d_engine_request_loads_the_powerhead_into_the_viewport() {
        // Clicking "Show 3-D engine (powerhead)" sets the request flag; the
        // deferred handler must load the detailed engine mesh into the central
        // viewport and clear the flag.
        let mut app = ValenxApp::default();
        app.show_engine_workbench = true;
        app.engine.show_engine_3d_request = true;
        draw(&mut app);
        assert!(
            app.mesh.is_some(),
            "detailed engine mesh loaded into viewport"
        );
        assert!(!app.engine.show_engine_3d_request, "request flag cleared");
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        let mut app = ValenxApp::default();
        app.show_engine_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 6,
            "expected numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        // Each of the six design-knob `DragValue`s is associated with its
        // caption via `labelled_by` (the AI-drivable name). The panel also
        // hosts two `Slider`s whose inner value field is itself a spin button
        // we do not caption, so require at least the six DragValues be
        // labelled rather than every spin button in the tree.
        let labelled = spin_buttons
            .iter()
            .filter(|n| !n.labelled_by().is_empty())
            .count();
        assert!(
            labelled >= 6,
            "every design-knob DragValue must be labelled_by its caption \
             (AI-drivable name); only {labelled} of {} spin buttons are labelled",
            spin_buttons.len()
        );
        for caption in ["chamber pressure", "chamber temp", "expansion ratio ε"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
