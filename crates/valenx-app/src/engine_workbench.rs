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
            last_export: None,
            propellant: Propellant::Ch4Lox,
            mixture_ratio: 3.6,
        }
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

/// Draw the Engine workbench right-side panel. A no-op when the
/// `show_engine_workbench` toggle is off.
pub fn draw_engine_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_engine_workbench {
        return;
    }

    egui::SidePanel::right("valenx_engine_workbench")
        .resizable(true)
        .default_width(380.0)
        .width_range(330.0..=620.0)
        .show(ctx, |ui| {
            ui.heading("Engine — design → analyze");
            ui.label(
                egui::RichText::new(
                    "ideal-nozzle performance + Bartz regen-cooling · valenx-astro",
                )
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
                        ui.label("chamber pressure");
                        let mut bar = s.design.chamber_pressure / 1.0e5;
                        if ui
                            .add(
                                egui::DragValue::new(&mut bar)
                                    .speed(1.0)
                                    .range(1.0..=400.0)
                                    .suffix(" bar"),
                            )
                            .changed()
                        {
                            s.design.chamber_pressure = bar * 1.0e5;
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.label("chamber temp");
                        ui.add(
                            egui::DragValue::new(&mut s.design.chamber_temperature)
                                .speed(10.0)
                                .range(1_000.0..=4_500.0)
                                .suffix(" K"),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("gas γ");
                        ui.add(
                            egui::DragValue::new(&mut s.design.gamma)
                                .speed(0.005)
                                .range(1.05..=1.40),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("molar mass");
                        ui.add(
                            egui::DragValue::new(&mut s.design.molar_mass)
                                .speed(0.2)
                                .range(5.0..=40.0)
                                .suffix(" g/mol"),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("throat area");
                        let mut cm2 = s.design.throat_area * 1.0e4;
                        if ui
                            .add(
                                egui::DragValue::new(&mut cm2)
                                    .speed(1.0)
                                    .range(1.0..=5_000.0)
                                    .suffix(" cm²"),
                            )
                            .changed()
                        {
                            s.design.throat_area = cm2 * 1.0e-4;
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.label("expansion ratio ε");
                        ui.add(
                            egui::DragValue::new(&mut s.design.expansion_ratio)
                                .speed(0.5)
                                .range(1.0..=200.0),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("target cooling margin");
                        ui.add(egui::Slider::new(&mut s.target_margin, 1.0..=4.0));
                    });

                    // ── Combustion chemistry (predict chamber conditions) ─
                    ui.add_space(4.0);
                    ui.separator();
                    ui.label(
                        egui::RichText::new("Combustion — equilibrium chamber prediction").strong(),
                    );
                    egui::ComboBox::from_label("propellant")
                        .selected_text(propellant_label(s.propellant))
                        .show_ui(ui, |ui| {
                            for p in [Propellant::H2Lox, Propellant::Rp1Lox, Propellant::Ch4Lox] {
                                ui.selectable_value(&mut s.propellant, p, propellant_label(p));
                            }
                        });
                    ui.horizontal(|ui| {
                        ui.label("mixture ratio (O/F)");
                        ui.add(egui::Slider::new(&mut s.mixture_ratio, 1.5..=8.0));
                    });
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
                        egui::RichText::new("Staged-combustion cycle — full-flow power balance")
                            .strong(),
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
                    ui.label(
                        egui::RichText::new("AI optimizer — max sea-level Isp vs cooling").strong(),
                    );
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

                    // ── 3-D nozzle + export ──────────────────────────────
                    ui.add_space(6.0);
                    ui.separator();
                    ui.label(egui::RichText::new("Nozzle geometry").strong());
                    ui.horizontal(|ui| {
                        if ui.button("Show 3-D nozzle").clicked() {
                            s.show_3d_request = true;
                        }
                        if ui.button("Export nozzle STL").clicked() {
                            let mesh =
                                nozzle_mesh(s.design.throat_area, s.design.expansion_ratio, 64);
                            let path = std::env::temp_dir().join("valenx_nozzle.stl");
                            s.last_export =
                                Some(match valenx_mesh::write_stl_binary(&mesh, &path) {
                                    Ok(()) => format!("exported → {}", path.display()),
                                    Err(e) => format!("export error: {e}"),
                                });
                        }
                    });
                    if let Some(msg) = &s.last_export {
                        ui.label(egui::RichText::new(msg).weak().small());
                    }
                });
        });

    // Deferred (outside the panel borrow): load the 3-D nozzle into the
    // central viewport when requested.
    if app.engine.show_3d_request {
        app.engine.show_3d_request = false;
        load_nozzle_3d(app);
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
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_engine_workbench(app, ctx);
        });
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
}
