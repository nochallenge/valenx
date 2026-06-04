//! The workflow sections of the Astro / Launch workbench panel.
//!
//! Each `draw_*_section` renders one step of the launch-simulation
//! workflow into the workbench side panel:
//!
//! - [`draw_vehicle_section`] — the stage stack + payload + reference
//!   area + launch site.
//! - [`draw_guidance_section`] — steering law, pitch program, target
//!   altitude, wind.
//! - [`draw_run_section`] — the Run button + a status line.
//! - [`draw_results_section`] — the orbit / max-Q / Δv summary cards and
//!   the flight-profile chart (altitude / velocity / Mach / dynamic
//!   pressure vs time).
//! - [`draw_planners_section`] — the closed-form mission planners
//!   (Hohmann, hoverslam, rendezvous, launch azimuth), each an
//!   input → output card with errors shown as text.
//!
//! Layout only — all the logic lives in [`super::model`] and
//! [`super::run`]. Mirrors the CFD-side [`crate::aero::panels`].

use eframe::egui;
use egui_plot::{Legend, Line, Plot, PlotPoints};
use nalgebra::Vector3;

use super::model::{self, AstroTab, GuidanceChoice, StageForm};
use crate::ValenxApp;
use valenx_astro::{landing, maneuver, rendezvous, windows};

// ---------------------------------------------------------------------------
// Shared little helpers — mirror aero/panels.rs so the two workbenches
// read the same.
// ---------------------------------------------------------------------------

/// A two-column key / value row inside a grid.
fn kv(ui: &mut egui::Ui, key: &str, value: impl Into<String>) {
    ui.label(egui::RichText::new(key).weak());
    ui.label(value.into());
    ui.end_row();
}

/// A derived-quantity line — a weak label with a computed value.
fn derived(ui: &mut egui::Ui, text: impl Into<String>) {
    ui.label(egui::RichText::new(text).weak().small());
}

/// A prominent result card — a bordered group with a small caption and
/// a large value. Used for the headline orbit / Δv cards.
fn result_card(ui: &mut egui::Ui, caption: &str, value: &str, tint: egui::Color32) {
    egui::Frame::group(ui.style())
        .fill(ui.visuals().faint_bg_color)
        .show(ui, |ui| {
            ui.vertical(|ui| {
                ui.set_min_width(82.0);
                ui.label(egui::RichText::new(caption).weak().small());
                ui.label(egui::RichText::new(value).heading().color(tint).strong());
            });
        });
}

/// A red error line, drawn only when `err` is `Some`. Sources its colour
/// from the canonical `accent::ERROR` token (matches aero).
fn error_line(ui: &mut egui::Ui, err: &Option<String>) {
    if let Some(e) = err {
        let hex = valenx_design_tokens::color::accent::ERROR;
        let c = egui::Color32::from_rgb(
            ((hex >> 16) & 0xFF) as u8,
            ((hex >> 8) & 0xFF) as u8,
            (hex & 0xFF) as u8,
        );
        ui.colored_label(c, format!("\u{2716} {e}"));
    }
}

// ---------------------------------------------------------------------------
// Tab selector
// ---------------------------------------------------------------------------

/// The Ascent / Planners tab strip at the top of the panel body.
pub fn draw_tab_selector(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        for tab in AstroTab::ALL {
            ui.selectable_value(&mut app.astro.tab, tab, tab.label());
        }
    });
    ui.separator();
}

// ---------------------------------------------------------------------------
// Ascent — 1. Vehicle / mission setup
// ---------------------------------------------------------------------------

/// Section 1 — the stage stack, payload, reference area and launch site.
pub fn draw_vehicle_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("1 \u{00b7} Vehicle")
        .id_source("astro_vehicle")
        .default_open(true)
        .show(ui, |ui| {
            let form = &mut app.astro.ascent;

            // Per-stage inputs. Index 0 fires first.
            let mut remove: Option<usize> = None;
            let n_stages = form.stages.len();
            for (i, stage) in form.stages.iter_mut().enumerate() {
                egui::CollapsingHeader::new(format!("Stage {} \u{2014} {}", i + 1, stage.name))
                    .id_source(("astro_stage", i))
                    .default_open(i == 0)
                    .show(ui, |ui| {
                        egui::Grid::new(("astro_stage_grid", i))
                            .num_columns(2)
                            .spacing([8.0, 4.0])
                            .show(ui, |ui| {
                                ui.label("Dry mass")
                                    .on_hover_text("Inert structural mass jettisoned at staging. SI unit: kg.");
                                ui.add(
                                    egui::DragValue::new(&mut stage.dry_mass)
                                        .speed(100.0)
                                        .range(1.0..=5.0e6)
                                        .suffix(" kg"),
                                );
                                ui.end_row();
                                ui.label("Propellant mass")
                                    .on_hover_text("Usable propellant in this stage. SI unit: kg.");
                                ui.add(
                                    egui::DragValue::new(&mut stage.propellant_mass)
                                        .speed(500.0)
                                        .range(0.0..=5.0e6)
                                        .suffix(" kg"),
                                );
                                ui.end_row();
                                ui.label("Thrust (vac)")
                                    .on_hover_text("Vacuum thrust. SI unit: N.");
                                ui.add(
                                    egui::DragValue::new(&mut stage.thrust_vac)
                                        .speed(10_000.0)
                                        .range(1.0..=1.0e8)
                                        .custom_formatter(|v, _| format!("{:.0} kN", v / 1_000.0)),
                                );
                                ui.end_row();
                                ui.label("Thrust (sea level)")
                                    .on_hover_text("Sea-level thrust. For an upper stage that never fires in the atmosphere set this equal to the vacuum thrust. SI unit: N.");
                                ui.add(
                                    egui::DragValue::new(&mut stage.thrust_sl)
                                        .speed(10_000.0)
                                        .range(1.0..=1.0e8)
                                        .custom_formatter(|v, _| format!("{:.0} kN", v / 1_000.0)),
                                );
                                ui.end_row();
                                ui.label("Isp (vac)")
                                    .on_hover_text("Vacuum specific impulse. SI unit: s.");
                                ui.add(
                                    egui::DragValue::new(&mut stage.isp_vac)
                                        .speed(1.0)
                                        .range(1.0..=600.0)
                                        .suffix(" s"),
                                );
                                ui.end_row();
                                ui.label("Isp (sea level)")
                                    .on_hover_text("Sea-level specific impulse. SI unit: s.");
                                ui.add(
                                    egui::DragValue::new(&mut stage.isp_sl)
                                        .speed(1.0)
                                        .range(1.0..=600.0)
                                        .suffix(" s"),
                                );
                                ui.end_row();
                            });
                        // Allow removing a stage only when more than one
                        // remains — a vehicle needs at least one stage.
                        if n_stages > 1 && ui.small_button("Remove this stage").clicked() {
                            remove = Some(i);
                        }
                    });
            }
            if let Some(i) = remove {
                form.stages.remove(i);
            }

            ui.add_space(2.0);
            if ui
                .button("+ Add stage")
                .on_hover_text("Append a stage to the top of the stack (fires last).")
                .clicked()
            {
                form.stages.push(StageForm {
                    name: format!("stage {}", form.stages.len() + 1),
                    dry_mass: 4_000.0,
                    propellant_mass: 100_000.0,
                    thrust_vac: 980_000.0,
                    thrust_sl: 980_000.0,
                    isp_vac: 348.0,
                    isp_sl: 348.0,
                });
            }

            ui.add_space(4.0);
            egui::Grid::new("astro_vehicle_grid")
                .num_columns(2)
                .spacing([8.0, 4.0])
                .show(ui, |ui| {
                    ui.label("Payload mass")
                        .on_hover_text("Payload carried all the way to orbit. SI unit: kg.");
                    ui.add(
                        egui::DragValue::new(&mut form.payload_mass)
                            .speed(100.0)
                            .range(0.0..=1.0e6)
                            .suffix(" kg"),
                    );
                    ui.end_row();
                    ui.label("Reference area")
                        .on_hover_text("Aerodynamic frontal area driving drag. SI unit: m².");
                    ui.add(
                        egui::DragValue::new(&mut form.reference_area)
                            .speed(0.1)
                            .range(0.1..=500.0)
                            .suffix(" m\u{00b2}"),
                    );
                    ui.end_row();
                    ui.label("Launch altitude")
                        .on_hover_text("Launch-site elevation above the equatorial radius. SI unit: m.");
                    ui.add(
                        egui::DragValue::new(&mut form.launch_altitude_m)
                            .speed(10.0)
                            .range(0.0..=6_000.0)
                            .suffix(" m"),
                    );
                    ui.end_row();
                });

            // Derived liftoff mass — a quick sanity read-out.
            let vehicle = form.build_vehicle();
            derived(
                ui,
                format!(
                    "Liftoff gross mass: {:.0} t \u{00b7} {} stage(s)",
                    vehicle.initial_mass() / 1_000.0,
                    vehicle.stages.len()
                ),
            );
        });
}

// ---------------------------------------------------------------------------
// Ascent — 2. Guidance
// ---------------------------------------------------------------------------

/// Section 2 — steering law, pitch program, target altitude, wind.
pub fn draw_guidance_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("2 \u{00b7} Guidance & target")
        .id_source("astro_guidance")
        .default_open(true)
        .show(ui, |ui| {
            let form = &mut app.astro.ascent;
            ui.label("Steering law");
            egui::ComboBox::from_id_source("astro_guidance_mode")
                .selected_text(form.guidance.label())
                .width(ui.available_width())
                .show_ui(ui, |ui| {
                    for g in GuidanceChoice::ALL {
                        ui.selectable_value(&mut form.guidance, g, g.label());
                    }
                });
            match form.guidance {
                GuidanceChoice::OpenLoopGravityTurn => derived(
                    ui,
                    "Burns every stage to depletion along the gravity turn — \
                     an overpowered vehicle ends on an eccentric orbit.",
                ),
                GuidanceChoice::ClosedLoopInsertion => derived(
                    ui,
                    "Ascends, coasts to apoapsis, then circularises at the \
                     target altitude — a near-circular low orbit.",
                ),
            }

            ui.add_space(4.0);
            egui::Grid::new("astro_guidance_grid")
                .num_columns(2)
                .spacing([8.0, 4.0])
                .show(ui, |ui| {
                    ui.label("Target altitude")
                        .on_hover_text("Circular-orbit altitude the closed-loop insertion targets. SI unit: km.");
                    ui.add(
                        egui::DragValue::new(&mut form.target_altitude_km)
                            .speed(5.0)
                            .range(80.0..=40_000.0)
                            .suffix(" km"),
                    );
                    ui.end_row();
                    ui.label("Vertical rise time")
                        .on_hover_text("How long to hold a vertical attitude after liftoff. SI unit: s.");
                    ui.add(
                        egui::DragValue::new(&mut form.vertical_rise_time)
                            .speed(0.5)
                            .range(0.0..=120.0)
                            .suffix(" s"),
                    );
                    ui.end_row();
                    ui.label("Pitch kick")
                        .on_hover_text("Pitch-over angle off vertical that starts the gravity turn. SI unit: deg.");
                    ui.add(
                        egui::DragValue::new(&mut form.pitch_kick_deg)
                            .speed(0.1)
                            .range(0.0..=45.0)
                            .suffix(" \u{00b0}"),
                    );
                    ui.end_row();
                    ui.label("Kick duration")
                        .on_hover_text("How long to hold the kicked attitude before the pure gravity turn. SI unit: s.");
                    ui.add(
                        egui::DragValue::new(&mut form.kick_duration)
                            .speed(0.5)
                            .range(0.0..=60.0)
                            .suffix(" s"),
                    );
                    ui.end_row();
                    ui.label("Steady wind")
                        .on_hover_text("Constant horizontal wind applied to the drag calculation (0 = calm air). SI unit: m/s.");
                    ui.add(
                        egui::DragValue::new(&mut form.wind_speed_ms)
                            .speed(1.0)
                            .range(0.0..=400.0)
                            .suffix(" m/s"),
                    );
                    ui.end_row();
                });
        });
}

// ---------------------------------------------------------------------------
// Ascent — 3. Run
// ---------------------------------------------------------------------------

/// Section 3 — the Run button + status line. The closed-form ascent is
/// fast and bounded, so it runs synchronously on click (see
/// [`super::run::run_ascent`]).
pub fn draw_run_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("3 \u{00b7} Run")
        .id_source("astro_run")
        .default_open(true)
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                let btn = egui::Button::new(
                    egui::RichText::new("\u{25b6}  Run ascent").strong(),
                );
                if ui
                    .add(btn)
                    .on_hover_text(
                        "Integrate the ascent from the pad to orbital insertion.\n\
                         The RK4 run is bounded, so it completes synchronously.\n\
                         Shortcut: Ctrl+R",
                    )
                    .clicked()
                {
                    super::run::run_ascent(app);
                }
                ui.separator();
                // Inline undo / redo — rewinds the ascent form to its
                // prior snapshot, recorded on each Run (matches aero).
                let undo = ui
                    .add_enabled(
                        app.astro.can_undo(),
                        egui::Button::new(valenx_icons::edit::UNDO),
                    )
                    .on_hover_text("Undo last form edit (Ctrl+Z)")
                    .clicked();
                let redo = ui
                    .add_enabled(
                        app.astro.can_redo(),
                        egui::Button::new(valenx_icons::edit::REDO),
                    )
                    .on_hover_text("Redo last form edit (Ctrl+Y)")
                    .clicked();
                if undo {
                    app.astro.undo_edit();
                }
                if redo {
                    app.astro.redo_edit();
                }
            });

            if !app.astro.status.is_empty() {
                ui.label(egui::RichText::new(&app.astro.status).weak());
            }
            error_line(ui, &app.astro.error);
        });
}

// ---------------------------------------------------------------------------
// Ascent — 4. Results
// ---------------------------------------------------------------------------

/// Section 4 — the orbit / max-Q / Δv summary cards + the flight-profile
/// chart.
pub fn draw_results_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("4 \u{00b7} Results")
        .id_source("astro_results")
        .default_open(true)
        .show(ui, |ui| {
            let Some(result) = app.astro.last_result.as_ref() else {
                derived(
                    ui,
                    "Run the ascent to see the orbit reached, the Δv budget, \
                     peak dynamic pressure and the flight profile.",
                );
                return;
            };

            // Headline cards.
            let orbit_tint = egui::Color32::from_rgb(120, 190, 235);
            let dv_tint = egui::Color32::from_rgb(150, 220, 150);
            let q_tint = egui::Color32::from_rgb(235, 150, 90);
            ui.horizontal_wrapped(|ui| {
                result_card(
                    ui,
                    "Apoapsis",
                    &format!("{:.0} km", result.apoapsis_km()),
                    orbit_tint,
                );
                result_card(
                    ui,
                    "Periapsis",
                    &format!("{:.0} km", result.periapsis_km()),
                    orbit_tint,
                );
                result_card(
                    ui,
                    "\u{0394}v budget",
                    &model::format_delta_v(result.ideal_delta_v),
                    dv_tint,
                );
                result_card(
                    ui,
                    "Max-Q",
                    &format!("{:.1} kPa", result.max_dynamic_pressure / 1_000.0),
                    q_tint,
                );
            });

            ui.add_space(4.0);
            ui.label(egui::RichText::new("Flight summary").strong());
            egui::Grid::new("astro_summary_grid")
                .num_columns(2)
                .spacing([8.0, 3.0])
                .show(ui, |ui| {
                    kv(ui, "outcome", format!("{:?}", result.outcome));
                    kv(
                        ui,
                        "reached space",
                        if result.reached_space { "yes (above 100 km)" } else { "no" },
                    );
                    kv(
                        ui,
                        "reached orbit",
                        if result.reached_orbit {
                            "yes (bound, periapsis > 100 km)"
                        } else {
                            "no (suborbital / re-entry arc)"
                        },
                    );
                    kv(ui, "eccentricity", format!("{:.4}", result.orbit.eccentricity));
                    kv(
                        ui,
                        "MECO time",
                        model::format_duration(result.final_time),
                    );
                    kv(
                        ui,
                        "MECO speed",
                        format!("{:.0} m/s", result.final_speed_inertial),
                    );
                    kv(
                        ui,
                        "max-Q altitude",
                        format!("{:.1} km", result.max_q_altitude_m / 1_000.0),
                    );
                    kv(ui, "peak g-load", format!("{:.1} g", result.max_acceleration_g));
                    kv(ui, "liftoff mass", format!("{:.0} t", result.liftoff_mass / 1_000.0));
                });

            // Flight milestones.
            if !result.events.is_empty() {
                egui::CollapsingHeader::new("Flight events")
                    .id_source("astro_events")
                    .default_open(false)
                    .show(ui, |ui| {
                        egui::Grid::new("astro_events_grid")
                            .num_columns(2)
                            .spacing([8.0, 3.0])
                            .show(ui, |ui| {
                                for e in &result.events {
                                    kv(
                                        ui,
                                        &model::format_duration(e.time),
                                        format!("{} ({:.0} km)", e.kind, e.altitude_m / 1_000.0),
                                    );
                                }
                            });
                    });
            }

            // Flight-profile chart — altitude / velocity / Mach /
            // dynamic-pressure vs time, all on shared, readable scales
            // (km, km/s, Mach, kPa) so one plot shows the whole story.
            // Built from plain Vecs exactly as aero/panels.rs does.
            if result.samples.len() > 1 {
                ui.add_space(4.0);
                ui.label(egui::RichText::new("Flight profile").strong());
                let altitude: PlotPoints = result
                    .samples
                    .iter()
                    .map(|s| [s.time, s.altitude_m / 1_000.0])
                    .collect();
                let velocity: PlotPoints = result
                    .samples
                    .iter()
                    .map(|s| [s.time, s.speed_inertial / 1_000.0])
                    .collect();
                let mach: PlotPoints = result
                    .samples
                    .iter()
                    .map(|s| [s.time, s.mach])
                    .collect();
                let qpress: PlotPoints = result
                    .samples
                    .iter()
                    .map(|s| [s.time, s.dynamic_pressure / 1_000.0])
                    .collect();
                Plot::new("astro_profile_plot")
                    .height(190.0)
                    .legend(Legend::default())
                    .x_axis_label("mission time (s)")
                    .show(ui, |plot_ui| {
                        plot_ui.line(Line::new(altitude).name("altitude (km)"));
                        plot_ui.line(Line::new(velocity).name("inertial speed (km/s)"));
                        plot_ui.line(Line::new(mach).name("Mach"));
                        plot_ui.line(Line::new(qpress).name("dynamic pressure (kPa)"));
                    });
                derived(
                    ui,
                    "Altitude in km, inertial speed in km/s, Mach dimensionless, \
                     dynamic pressure in kPa \u{2014} sharing one time axis.",
                );
            }
        });
}

// ---------------------------------------------------------------------------
// Planners tab
// ---------------------------------------------------------------------------

/// The Planners tab — four independent closed-form mission planners,
/// each an input → output card with errors shown as text.
pub fn draw_planners_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    draw_hohmann_planner(app, ui);
    draw_plane_change_planner(app, ui);
    draw_hoverslam_planner(app, ui);
    draw_rendezvous_planner(app, ui);
    draw_azimuth_planner(app, ui);
}

/// Hohmann Δv between two circular altitudes (`maneuver::hohmann_transfer`).
fn draw_hohmann_planner(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("Hohmann transfer \u{0394}v")
        .id_source("astro_hohmann")
        .default_open(true)
        .show(ui, |ui| {
            let form = &mut app.astro.planner;
            egui::Grid::new("astro_hohmann_grid")
                .num_columns(2)
                .spacing([8.0, 4.0])
                .show(ui, |ui| {
                    ui.label("From altitude")
                        .on_hover_text("Departure circular-orbit altitude. SI unit: km.");
                    ui.add(
                        egui::DragValue::new(&mut form.hohmann_from_km)
                            .speed(10.0)
                            .range(80.0..=400_000.0)
                            .suffix(" km"),
                    );
                    ui.end_row();
                    ui.label("To altitude")
                        .on_hover_text("Arrival circular-orbit altitude. SI unit: km.");
                    ui.add(
                        egui::DragValue::new(&mut form.hohmann_to_km)
                            .speed(10.0)
                            .range(80.0..=400_000.0)
                            .suffix(" km"),
                    );
                    ui.end_row();
                });

            let r1 = model::altitude_km_to_radius_m(form.hohmann_from_km);
            let r2 = model::altitude_km_to_radius_m(form.hohmann_to_km);
            match maneuver::hohmann_transfer(r1, r2) {
                Ok(t) => {
                    egui::Grid::new("astro_hohmann_out")
                        .num_columns(2)
                        .spacing([8.0, 3.0])
                        .show(ui, |ui| {
                            kv(ui, "burn 1 \u{0394}v", model::format_delta_v(t.delta_v1));
                            kv(ui, "burn 2 \u{0394}v", model::format_delta_v(t.delta_v2));
                            kv(ui, "total \u{0394}v", model::format_delta_v(t.total_delta_v));
                            kv(ui, "transfer time", model::format_duration(t.transfer_time));
                        });
                }
                Err(e) => derived(ui, format!("Cannot plan: {}", model::friendly_error(&e))),
            }
        });
}

/// Pure plane-change Δv on a circular orbit
/// (`maneuver::circular_plane_change_dv`).
fn draw_plane_change_planner(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("Plane change \u{0394}v")
        .id_source("astro_plane_change")
        .default_open(false)
        .show(ui, |ui| {
            let form = &mut app.astro.planner;
            egui::Grid::new("astro_plane_change_grid")
                .num_columns(2)
                .spacing([8.0, 4.0])
                .show(ui, |ui| {
                    ui.label("Orbit altitude").on_hover_text(
                        "Circular-orbit altitude where the burn happens. SI unit: km.",
                    );
                    ui.add(
                        egui::DragValue::new(&mut form.plane_change_altitude_km)
                            .speed(10.0)
                            .range(80.0..=400_000.0)
                            .suffix(" km"),
                    );
                    ui.end_row();
                    ui.label("Inclination change")
                        .on_hover_text("Pure plane rotation Δi. SI unit: deg.");
                    ui.add(
                        egui::DragValue::new(&mut form.plane_change_delta_inc_deg)
                            .speed(0.5)
                            .range(0.0..=180.0)
                            .suffix(" deg"),
                    );
                    ui.end_row();
                });

            match model::plane_change_dv(
                form.plane_change_altitude_km,
                form.plane_change_delta_inc_deg,
            ) {
                Ok(dv) => {
                    egui::Grid::new("astro_plane_change_out")
                        .num_columns(2)
                        .spacing([8.0, 3.0])
                        .show(ui, |ui| {
                            kv(ui, "plane-change \u{0394}v", model::format_delta_v(dv));
                        });
                    derived(
                        ui,
                        "Pure plane changes are costly (Δv = 2·v·sin(Δi/2)); rotate where v is low or fold into another burn."
                            .to_string(),
                    );
                }
                Err(e) => derived(ui, format!("Cannot plan: {}", model::friendly_error(&e))),
            }
        });
}

/// Hoverslam ignition altitude (`landing::ignition_altitude`).
fn draw_hoverslam_planner(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("Landing hoverslam ignition altitude")
        .id_source("astro_hoverslam")
        .default_open(false)
        .show(ui, |ui| {
            let form = &mut app.astro.planner;
            egui::Grid::new("astro_hoverslam_grid")
                .num_columns(2)
                .spacing([8.0, 4.0])
                .show(ui, |ui| {
                    ui.label("Descent speed")
                        .on_hover_text("Vehicle descent speed at the start of the burn. SI unit: m/s.");
                    ui.add(
                        egui::DragValue::new(&mut form.hoverslam_descent_speed)
                            .speed(1.0)
                            .range(0.0..=2_000.0)
                            .suffix(" m/s"),
                    );
                    ui.end_row();
                    ui.label("Net deceleration")
                        .on_hover_text("Thrust-minus-weight net deceleration T/m − g. Must be positive to land. SI unit: m/s².");
                    ui.add(
                        egui::DragValue::new(&mut form.hoverslam_net_decel)
                            .speed(0.5)
                            .range(-50.0..=100.0)
                            .suffix(" m/s\u{00b2}"),
                    );
                    ui.end_row();
                });

            match landing::ignition_altitude(form.hoverslam_descent_speed, form.hoverslam_net_decel)
            {
                Ok(h) => {
                    egui::Grid::new("astro_hoverslam_out")
                        .num_columns(2)
                        .spacing([8.0, 3.0])
                        .show(ui, |ui| {
                            kv(ui, "ignition altitude", format!("{h:.1} m"));
                            kv(ui, "h = v\u{00b2} / (2\u{00b7}a_net)", format!("{:.3} km", h / 1_000.0));
                        });
                }
                Err(e) => derived(ui, format!("Cannot land: {}", model::friendly_error(&e))),
            }
        });
}

/// Two-impulse rendezvous Δv (`rendezvous::two_impulse_rendezvous`).
fn draw_rendezvous_planner(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("Rendezvous two-impulse \u{0394}v")
        .id_source("astro_rendezvous")
        .default_open(false)
        .show(ui, |ui| {
            let form = &mut app.astro.planner;
            egui::Grid::new("astro_rdv_grid")
                .num_columns(2)
                .spacing([8.0, 4.0])
                .show(ui, |ui| {
                    ui.label("Target orbit altitude")
                        .on_hover_text("Circular-orbit altitude of the target, which sets the mean motion n. SI unit: km.");
                    ui.add(
                        egui::DragValue::new(&mut form.rdv_orbit_altitude_km)
                            .speed(5.0)
                            .range(80.0..=40_000.0)
                            .suffix(" km"),
                    );
                    ui.end_row();
                    ui.label("Radial offset")
                        .on_hover_text("Initial chaser radial offset in the target's LVLH frame. SI unit: m.");
                    ui.add(
                        egui::DragValue::new(&mut form.rdv_offset_radial)
                            .speed(50.0)
                            .range(-100_000.0..=100_000.0)
                            .suffix(" m"),
                    );
                    ui.end_row();
                    ui.label("Along-track offset")
                        .on_hover_text("Initial chaser along-track offset in LVLH. SI unit: m.");
                    ui.add(
                        egui::DragValue::new(&mut form.rdv_offset_along)
                            .speed(50.0)
                            .range(-100_000.0..=100_000.0)
                            .suffix(" m"),
                    );
                    ui.end_row();
                    ui.label("Transfer time")
                        .on_hover_text("Transfer time as a fraction of the orbital period. Avoid whole half-periods (0.5, 1.0, …) where the boundary-value problem is singular.");
                    ui.add(
                        egui::DragValue::new(&mut form.rdv_transfer_fraction)
                            .speed(0.01)
                            .range(0.05..=0.95)
                            .custom_formatter(|v, _| format!("{:.0} % of period", v * 100.0)),
                    );
                    ui.end_row();
                });

            // Build the inputs: mean motion from the target's semi-major
            // axis, transfer time as a fraction of the period.
            let a = model::altitude_km_to_radius_m(form.rdv_orbit_altitude_km);
            match rendezvous::mean_motion(a) {
                Ok(n) => {
                    let period = std::f64::consts::TAU / n;
                    let transfer_time = period * form.rdv_transfer_fraction;
                    let r0 = Vector3::new(form.rdv_offset_radial, form.rdv_offset_along, 0.0);
                    let v0 = Vector3::zeros();
                    let z = Vector3::zeros();
                    match rendezvous::two_impulse_rendezvous(n, transfer_time, r0, v0, z, z) {
                        Ok(plan) => {
                            egui::Grid::new("astro_rdv_out")
                                .num_columns(2)
                                .spacing([8.0, 3.0])
                                .show(ui, |ui| {
                                    kv(ui, "burn 1 \u{0394}v", format!("{:.2} m/s", plan.delta_v1.norm()));
                                    kv(ui, "burn 2 \u{0394}v", format!("{:.2} m/s", plan.delta_v2.norm()));
                                    kv(ui, "total \u{0394}v", model::format_delta_v(plan.total_delta_v));
                                    kv(ui, "transfer time", model::format_duration(plan.transfer_time));
                                    kv(ui, "orbital period", model::format_duration(period));
                                });
                        }
                        Err(e) => derived(
                            ui,
                            format!("Cannot solve: {}", model::friendly_error(&e)),
                        ),
                    }
                }
                Err(e) => derived(ui, format!("Cannot solve: {}", model::friendly_error(&e))),
            }
        });
}

/// Launch azimuth for a target inclination (`windows::launch_azimuth`).
fn draw_azimuth_planner(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("Launch azimuth")
        .id_source("astro_azimuth")
        .default_open(false)
        .show(ui, |ui| {
            let form = &mut app.astro.planner;
            egui::Grid::new("astro_azimuth_grid")
                .num_columns(2)
                .spacing([8.0, 4.0])
                .show(ui, |ui| {
                    ui.label("Launch latitude")
                        .on_hover_text("Geodetic latitude of the launch site. SI unit: deg.");
                    ui.add(
                        egui::DragValue::new(&mut form.azimuth_latitude_deg)
                            .speed(0.1)
                            .range(-90.0..=90.0)
                            .suffix(" \u{00b0}"),
                    );
                    ui.end_row();
                    ui.label("Target inclination")
                        .on_hover_text("Inclination of the target orbit. Must be at least |latitude| to reach directly. SI unit: deg.");
                    ui.add(
                        egui::DragValue::new(&mut form.azimuth_inclination_deg)
                            .speed(0.1)
                            .range(0.0..=180.0)
                            .suffix(" \u{00b0}"),
                    );
                    ui.end_row();
                });

            let lat = form.azimuth_latitude_deg.to_radians();
            let inc = form.azimuth_inclination_deg.to_radians();
            match windows::launch_azimuth(lat, inc) {
                Ok(az) => {
                    egui::Grid::new("astro_azimuth_out")
                        .num_columns(2)
                        .spacing([8.0, 3.0])
                        .show(ui, |ui| {
                            kv(
                                ui,
                                "ascending azimuth",
                                format!("{:.2} \u{00b0} (clockwise from north)", az.ascending.to_degrees()),
                            );
                            kv(
                                ui,
                                "descending azimuth",
                                format!("{:.2} \u{00b0}", az.descending.to_degrees()),
                            );
                        });
                }
                Err(e) => derived(
                    ui,
                    format!("Unreachable directly: {}", model::friendly_error(&e)),
                ),
            }
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn helpers_render_in_a_headless_grid() {
        // The little kv / derived / result_card helpers must render in a
        // windowless context without panicking.
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                egui::Grid::new("t").show(ui, |ui| {
                    kv(ui, "k", "v");
                });
                derived(ui, "note");
                result_card(ui, "cap", "42", ui.visuals().strong_text_color());
                error_line(ui, &Some("boom".to_string()));
                error_line(ui, &None);
            });
        });
    }
}
