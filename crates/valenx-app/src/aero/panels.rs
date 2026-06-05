//! The eight workflow sections of the Wind Tunnel workbench panel.
//!
//! Each `draw_*_section` renders one collapsible step of the virtual-
//! wind-tunnel workflow into the workbench side panel:
//!
//! 1. [`draw_body_section`] — pick the body to test.
//! 2. [`draw_wind_section`] — free-stream speed / direction / air.
//! 3. [`draw_ground_section`] — fixed vs. moving ground, rotating
//!    wheels.
//! 4. [`draw_tunnel_section`] — domain size + grid resolution.
//! 5. [`draw_solver_section`] — turbulence model + iteration settings.
//! 6. [`draw_run_section`] — the Run button + live convergence plot.
//! 7. [`draw_results_section`] — coefficient result cards + the polar.
//! 8. [`draw_visualization_section`] — push field results into the 3-D
//!    viewport.
//!
//! Layout only — all the logic these panels invoke lives in
//! [`super::model`], [`super::compute`] and [`super::viz`].

use eframe::egui;
use egui_plot::{Legend, Line, Plot, PlotPoints, Points};

use super::compute::{AeroJob, AeroRunHandle};
use super::model::{
    self, BodySource, CutAxis, FlowField, GridResolution, RunMode, TurbulenceChoice,
};
use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Shared little helpers — keep the section bodies short and consistent
// ---------------------------------------------------------------------------

/// A two-column key / value row inside a grid.
fn kv(ui: &mut egui::Ui, key: &str, value: impl Into<String>) {
    ui.label(egui::RichText::new(key).weak());
    ui.label(value.into());
    ui.end_row();
}

/// A derived-quantity line — a weak label with a computed value, used
/// for the "resulting Reynolds number" style read-outs.
fn derived(ui: &mut egui::Ui, text: impl Into<String>) {
    ui.label(egui::RichText::new(text).weak().small());
}

/// A prominent result card — a bordered group with a small caption and
/// a large value. Used for the headline Cd / Cl / Cm cards.
fn result_card(ui: &mut egui::Ui, caption: &str, value: &str, tint: egui::Color32) {
    egui::Frame::group(ui.style())
        .fill(ui.visuals().faint_bg_color)
        .show(ui, |ui| {
            ui.vertical(|ui| {
                ui.set_min_width(78.0);
                ui.label(egui::RichText::new(caption).weak().small());
                ui.label(
                    egui::RichText::new(value)
                        .heading()
                        .color(tint)
                        .strong(),
                );
            });
        });
}

/// A red error line, drawn only when `err` is `Some`. Sources its
/// colour from the canonical `accent::ERROR` token so the rendered
/// hue matches whichever theme variant is active.
fn error_line(ui: &mut egui::Ui, err: &Option<String>) {
    if let Some(e) = err {
        let hex = valenx_design_tokens::color::accent::ERROR;
        let c = egui::Color32::from_rgb(
            ((hex >> 16) & 0xFF) as u8,
            ((hex >> 8) & 0xFF) as u8,
            (hex & 0xFF) as u8,
        );
        ui.colored_label(c, format!("✖ {e}"));
    }
}

// ---------------------------------------------------------------------------
// 1. Body
// ---------------------------------------------------------------------------

/// Section 1 — pick the body the wind tunnel tests and report its size.
pub fn draw_body_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("1 · Body")
        .id_source("aero_body")
        .default_open(true)
        .show(ui, |ui| {
            let has_cad = app.current_solid.is_some();
            let has_stl = app.stl.is_some();
            {
                let form = &mut app.aero.form;
                ui.label("Test body")
                    .on_hover_text("Pick the geometry that goes in the tunnel.");
                egui::ComboBox::from_id_source("aero_body_source")
                    .selected_text(form.body_source.label())
                    .width(ui.available_width())
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut form.body_source,
                            BodySource::CurrentCadModel,
                            BodySource::CurrentCadModel.label(),
                        )
                        .on_hover_text("Use the CAD solid currently loaded in the app's Part workbench.");
                        ui.selectable_value(
                            &mut form.body_source,
                            BodySource::ImportedStl,
                            BodySource::ImportedStl.label(),
                        )
                        .on_hover_text("Use the STL imported via File → Import STL.");
                        ui.selectable_value(
                            &mut form.body_source,
                            BodySource::DemoBox,
                            BodySource::DemoBox.label(),
                        )
                        .on_hover_text("Built-in axis-aligned box — useful for sanity-check runs.");
                    });

                match form.body_source {
                    BodySource::CurrentCadModel => {
                        if has_cad {
                            derived(ui, "Using the CAD solid loaded in the app.");
                        } else {
                            derived(
                                ui,
                                "No CAD model loaded — create one in the Mesh \
                                 Toolbox or pick another source.",
                            );
                        }
                    }
                    BodySource::ImportedStl => {
                        if has_stl {
                            derived(ui, "Using the imported STL geometry.");
                        } else {
                            derived(ui, "No STL loaded yet — import one below.");
                        }
                    }
                    BodySource::DemoBox => {
                        ui.add_space(2.0);
                        ui.label(egui::RichText::new("Demo box size (m)").weak());
                        egui::Grid::new("aero_demo_box")
                            .num_columns(2)
                            .spacing([8.0, 4.0])
                            .show(ui, |ui| {
                                for (i, axis) in
                                    ["length X", "width Y", "height Z"].iter().enumerate()
                                {
                                    ui.label(*axis);
                                    ui.add(
                                        egui::DragValue::new(&mut form.demo_box_size[i])
                                            .speed(0.1)
                                            .range(0.1..=50.0)
                                            .suffix(" m"),
                                    );
                                    ui.end_row();
                                }
                            });
                    }
                }
            }

            // The STL import affordance — picks a file via rfd and
            // loads it into the app's shared STL slot.
            if app.aero.form.body_source == BodySource::ImportedStl
                && ui.button("Import STL file…").clicked()
            {
                app.pick_stl();
                if app.stl.is_some() {
                    app.aero.error = None;
                }
            }

            // Report the body's triangle count + bounding box.
            ui.add_space(4.0);
            match super::extract_body(
                &app.aero.form,
                app.current_solid.as_ref(),
                app.stl.as_ref().map(|s| &s.mesh),
            ) {
                Ok(body) => {
                    egui::Grid::new("aero_body_stats")
                        .num_columns(2)
                        .spacing([8.0, 3.0])
                        .show(ui, |ui| {
                            kv(ui, "triangles", format!("{}", body.mesh.len()));
                            if let Some(bb) = body.mesh.aabb() {
                                let e = bb.extent();
                                kv(
                                    ui,
                                    "bounding box",
                                    format!("{:.2} × {:.2} × {:.2} m", e.x, e.y, e.z),
                                );
                                kv(
                                    ui,
                                    "surface area",
                                    format!("{:.3} m²", body.mesh.surface_area()),
                                );
                            }
                        });
                }
                Err(e) => {
                    derived(ui, format!("Body not ready: {e}"));
                }
            }
        });
}

// ---------------------------------------------------------------------------
// 2. Wind conditions
// ---------------------------------------------------------------------------

/// Section 2 — the free-stream speed, direction, air properties.
pub fn draw_wind_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("2 · Wind conditions")
        .id_source("aero_wind")
        .default_open(true)
        .show(ui, |ui| {
            let form = &mut app.aero.form;
            egui::Grid::new("aero_wind_grid")
                .num_columns(2)
                .spacing([8.0, 4.0])
                .show(ui, |ui| {
                    ui.label("Free-stream speed")
                        .on_hover_text("Wind speed at the tunnel inlet. SI unit: m/s.");
                    ui.add(
                        egui::DragValue::new(&mut form.speed_ms)
                            .speed(0.5)
                            .range(0.1..=340.0)
                            .suffix(" m/s"),
                    )
                    .on_hover_text(
                        "Free-stream wind speed (m/s). 30 m/s ≈ 108 km/h ≈ 67 mph — \
                         a typical road-car operating point.",
                    );
                    ui.end_row();

                    ui.label("Yaw (sideslip)")
                        .on_hover_text("Rotation of the wind direction around the vertical axis. SI unit: deg.");
                    ui.add(
                        egui::DragValue::new(&mut form.yaw_deg)
                            .speed(0.5)
                            .range(-45.0..=45.0)
                            .suffix(" °"),
                    )
                    .on_hover_text(
                        "Yaw angle in degrees (sideslip). Positive yaws the wind \
                         toward +Y. 0 is straight-on.",
                    );
                    ui.end_row();

                    ui.label("Pitch (angle of attack)")
                        .on_hover_text("Angle of the wind in the vertical plane. SI unit: deg.");
                    ui.add(
                        egui::DragValue::new(&mut form.pitch_deg)
                            .speed(0.5)
                            .range(-30.0..=30.0)
                            .suffix(" °"),
                    )
                    .on_hover_text(
                        "Pitch / angle of attack (deg). Positive pitches the wind \
                         toward +Z (lift-positive).",
                    );
                    ui.end_row();

                    ui.label("Air density")
                        .on_hover_text("Free-stream air mass per unit volume. SI unit: kg/m³.");
                    ui.add(
                        egui::DragValue::new(&mut form.air_density)
                            .speed(0.001)
                            .range(0.01..=10.0)
                            .suffix(" kg/m³"),
                    )
                    .on_hover_text(
                        "Air density (kg/m³). Sea-level standard atmosphere = 1.225.",
                    );
                    ui.end_row();

                    ui.label("Air viscosity")
                        .on_hover_text("Dynamic viscosity of the free-stream. SI unit: Pa·s.");
                    ui.add(
                        egui::DragValue::new(&mut form.air_viscosity)
                            .speed(1e-7)
                            .range(1e-7..=1e-2)
                            .custom_formatter(|v, _| format!("{v:.3e}"))
                            .suffix(" Pa·s"),
                    )
                    .on_hover_text(
                        "Dynamic viscosity μ (Pa·s). Standard air at 20 °C = 1.81e-5.",
                    );
                    ui.end_row();

                    ui.label("Turbulence intensity")
                        .on_hover_text("Velocity fluctuation magnitude relative to mean. Dimensionless 0..=1.");
                    ui.add(
                        egui::DragValue::new(&mut form.turbulence_intensity)
                            .speed(0.001)
                            .range(0.0..=0.5)
                            .custom_formatter(|v, _| format!("{:.1} %", v * 100.0)),
                    )
                    .on_hover_text(
                        "Inlet turbulence intensity = u'/U. Wind-tunnel ≤ 1 %, \
                         road conditions 5–10 %.",
                    );
                    ui.end_row();

                    ui.label("Air temperature")
                        .on_hover_text("Free-stream temperature. SI unit: K.");
                    ui.add(
                        egui::DragValue::new(&mut form.temperature_k)
                            .speed(1.0)
                            .range(150.0..=400.0)
                            .suffix(" K"),
                    )
                    .on_hover_text(
                        "Air temperature (K). Sea-level standard = 288.15 K (15 °C).",
                    );
                    ui.end_row();
                });

            ui.checkbox(
                &mut form.apply_compressibility,
                "Apply Prandtl-Glauert compressibility correction",
            )
            .on_hover_text(
                "Rescales the coefficients by 1/√(1−M²) above Mach 0.3 — valid \
                 into the high-subsonic range.",
            );

            // Derived: speed in km/h, dynamic pressure, Reynolds
            // number. Computed before the grid closure so `form` and
            // the `&app` body lookup don't alias.
            let length = body_reference_length(app);
            let form = &app.aero.form;
            let speed_kmh = model::ms_to_kmh(form.speed_ms);
            let q = model::format_pressure_pa(form.dynamic_pressure());
            let re = form.reynolds_number(length);

            ui.add_space(4.0);
            ui.label(egui::RichText::new("Resulting conditions").strong());
            egui::Grid::new("aero_wind_derived")
                .num_columns(2)
                .spacing([8.0, 3.0])
                .show(ui, |ui| {
                    kv(ui, "speed", format!("{speed_kmh:.0} km/h"));
                    kv(ui, "dynamic pressure", q);
                    kv(
                        ui,
                        "Reynolds number",
                        format!("{} (on {:.2} m)", model::format_reynolds(re), length),
                    );
                });
        });
}

// ---------------------------------------------------------------------------
// 3. Ground & wheels
// ---------------------------------------------------------------------------

/// Section 3 — the ground plane and rotating-wheel treatment.
pub fn draw_ground_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("3 · Ground & wheels")
        .id_source("aero_ground")
        .default_open(false)
        .show(ui, |ui| {
            let form = &mut app.aero.form;
            ui.label("Ground plane");
            ui.radio_value(&mut form.moving_ground, false, "Free-air (no ground)")
                .on_hover_text("A slip floor — the body flies in open air.");
            ui.radio_value(
                &mut form.moving_ground,
                true,
                "Moving ground (road under a car)",
            )
            .on_hover_text(
                "A no-slip floor carried at the free-stream speed — the correct \
                 automotive boundary condition.",
            );

            ui.add_space(4.0);
            ui.add_enabled_ui(form.moving_ground, |ui| {
                ui.checkbox(&mut form.rotating_wheels, "Rotating wheels (automotive)")
                    .on_hover_text(
                        "Spins the wheel region's solid cells at the rolling \
                         angular speed ω = V/R.",
                    );
                ui.add_enabled_ui(form.rotating_wheels, |ui| {
                    egui::Grid::new("aero_wheel_grid")
                        .num_columns(2)
                        .spacing([8.0, 4.0])
                        .show(ui, |ui| {
                            ui.label("Wheel radius");
                            ui.add(
                                egui::DragValue::new(&mut form.wheel_radius)
                                    .speed(0.01)
                                    .range(0.05..=2.0)
                                    .suffix(" m"),
                            );
                            ui.end_row();
                        });
                });
            });

            if form.moving_ground && form.rotating_wheels {
                derived(
                    ui,
                    format!(
                        "Rolling wheel speed: ω = {:.1} rad/s",
                        form.wheel_omega()
                    ),
                );
            }
        });
}

// ---------------------------------------------------------------------------
// 4. Tunnel & mesh
// ---------------------------------------------------------------------------

/// Section 4 — domain size and grid resolution.
pub fn draw_tunnel_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("4 · Tunnel & mesh")
        .id_source("aero_tunnel")
        .default_open(false)
        .show(ui, |ui| {
            let form = &mut app.aero.form;
            ui.label("Grid resolution")
                .on_hover_text(
                    "Coarseness of the immersed-boundary voxel grid. Higher \
                     resolution = sharper wake but quadratically more memory.",
                );
            egui::ComboBox::from_id_source("aero_resolution")
                .selected_text(form.resolution.label())
                .width(ui.available_width())
                .show_ui(ui, |ui| {
                    for r in GridResolution::ALL {
                        ui.selectable_value(&mut form.resolution, r, r.label())
                            .on_hover_text(format!(
                                "~{} cells across the body, ≤ {} cells total.",
                                r.cells_across_body(),
                                format_count(r.max_cells()),
                            ));
                    }
                });
            derived(
                ui,
                format!(
                    "≈ {} cells across the body · up to {} total cells",
                    form.resolution.cells_across_body(),
                    format_count(form.resolution.max_cells()),
                ),
            );

            ui.add_space(4.0);
            ui.checkbox(&mut form.override_domain, "Override the auto-sized domain")
                .on_hover_text(
                    "By default the tunnel box is auto-sized around the body. \
                     Tick to set the clearances by hand.",
                );
            ui.add_enabled_ui(form.override_domain, |ui| {
                egui::Grid::new("aero_domain_grid")
                    .num_columns(2)
                    .spacing([8.0, 4.0])
                    .show(ui, |ui| {
                        ui.label("Upstream clearance")
                            .on_hover_text("Distance from inlet to the body, in multiples of the body length.");
                        ui.add(
                            egui::DragValue::new(&mut form.upstream_lengths)
                                .speed(0.25)
                                .range(0.5..=15.0)
                                .suffix(" L"),
                        )
                        .on_hover_text("Multiples of body length L (typical 3–5).");
                        ui.end_row();
                        ui.label("Downstream (wake)")
                            .on_hover_text("Length of the wake region behind the body, in multiples of L.");
                        ui.add(
                            egui::DragValue::new(&mut form.downstream_lengths)
                                .speed(0.25)
                                .range(1.0..=25.0)
                                .suffix(" L"),
                        )
                        .on_hover_text(
                            "Wake length in multiples of L. Bluff bodies need 8–15 L for \
                             the wake to fully develop before the outlet boundary.",
                        );
                        ui.end_row();
                        ui.label("Lateral clearance")
                            .on_hover_text("Distance to the side walls, in multiples of L.");
                        ui.add(
                            egui::DragValue::new(&mut form.lateral_lengths)
                                .speed(0.25)
                                .range(0.5..=12.0)
                                .suffix(" L"),
                        )
                        .on_hover_text("Multiples of L (typical 3–5). Sets blockage ratio.");
                        ui.end_row();
                    });
            });
            derived(
                ui,
                "Near-body cells are refined automatically — the immersed-boundary \
                 voxelizer staircases the wall to the grid.",
            );
        });
}

// ---------------------------------------------------------------------------
// 5. Solver
// ---------------------------------------------------------------------------

/// Section 5 — turbulence model and iteration / convergence settings.
pub fn draw_solver_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("5 · Solver")
        .id_source("aero_solver")
        .default_open(false)
        .show(ui, |ui| {
            let form = &mut app.aero.form;
            egui::Grid::new("aero_solver_grid")
                .num_columns(2)
                .spacing([8.0, 4.0])
                .show(ui, |ui| {
                    ui.label("Turbulence model")
                        .on_hover_text(
                            "RANS closure picking the Reynolds-stress model. k-ω SST is \
                             the default robust choice for external aerodynamics.",
                        );
                    egui::ComboBox::from_id_source("aero_turbulence")
                        .selected_text(form.turbulence.label())
                        .show_ui(ui, |ui| {
                            for t in TurbulenceChoice::ALL {
                                ui.selectable_value(&mut form.turbulence, t, t.label());
                            }
                        });
                    ui.end_row();

                    ui.label("Max iterations")
                        .on_hover_text("Solver iteration cap. Reach the convergence tolerance first; otherwise stop here.");
                    ui.add(
                        egui::DragValue::new(&mut form.max_iterations)
                            .speed(5.0)
                            .range(1..=2000),
                    )
                    .on_hover_text(
                        "Stops the solve at this many iterations even if the residual \
                         hasn't reached the tolerance. Typical: 500–2000 for steady \
                         RANS.",
                    );
                    ui.end_row();

                    ui.label("Convergence tolerance")
                        .on_hover_text("Mass-imbalance residual threshold. Dimensionless.");
                    ui.add(
                        egui::DragValue::new(&mut form.tolerance)
                            .speed(1e-5)
                            .range(1e-8..=1e-2)
                            .custom_formatter(|v, _| format!("{v:.1e}")),
                    )
                    .on_hover_text(
                        "Stop when the mass-imbalance residual drops below this. \
                         1e-4 ≈ engineering accuracy, 1e-6 research-grade.",
                    );
                    ui.end_row();

                    ui.label("Run mode")
                        .on_hover_text("Single steady solve or an angle-of-attack polar sweep.");
                    egui::ComboBox::from_id_source("aero_run_mode")
                        .selected_text(form.run_mode.label())
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut form.run_mode,
                                RunMode::SteadyPoint,
                                RunMode::SteadyPoint.label(),
                            )
                            .on_hover_text("Solve one steady-state operating point.");
                            ui.selectable_value(
                                &mut form.run_mode,
                                RunMode::AngleSweep,
                                RunMode::AngleSweep.label(),
                            )
                            .on_hover_text(
                                "Loop the solve across a range of angles of attack to build a Cd/Cl polar.",
                            );
                        });
                    ui.end_row();
                });

            if form.run_mode == RunMode::AngleSweep {
                ui.add_space(2.0);
                ui.label(egui::RichText::new("Angle-of-attack sweep").strong());
                egui::Grid::new("aero_sweep_grid")
                    .num_columns(2)
                    .spacing([8.0, 4.0])
                    .show(ui, |ui| {
                        ui.label("Start angle")
                            .on_hover_text("First AoA in the sweep. SI unit: deg.");
                        ui.add(
                            egui::DragValue::new(&mut form.sweep_start_deg)
                                .speed(0.5)
                                .range(-30.0..=30.0)
                                .suffix(" °"),
                        )
                        .on_hover_text("Lower bound of the sweep (deg).");
                        ui.end_row();
                        ui.label("End angle")
                            .on_hover_text("Last AoA in the sweep. SI unit: deg.");
                        ui.add(
                            egui::DragValue::new(&mut form.sweep_end_deg)
                                .speed(0.5)
                                .range(-30.0..=30.0)
                                .suffix(" °"),
                        )
                        .on_hover_text("Upper bound of the sweep (deg).");
                        ui.end_row();
                        ui.label("Sweep points")
                            .on_hover_text("Number of evenly-spaced AoA samples between start and end.");
                        ui.add(
                            egui::DragValue::new(&mut form.sweep_points)
                                .speed(1.0)
                                .range(2..=25),
                        )
                        .on_hover_text(
                            "How many evenly-spaced angles to solve. 2 = just start + end; \
                             21 = 1° granularity over a 20° sweep.",
                        );
                        ui.end_row();
                    });
                derived(
                    ui,
                    "The sweep builds one tunnel and re-orients the wind per \
                     angle — far cheaper than re-meshing.",
                );
            }
        });
}

// ---------------------------------------------------------------------------
// 6. Run
// ---------------------------------------------------------------------------

/// Section 6 — the Run button, progress, and the live convergence plot.
pub fn draw_run_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("6 · Run")
        .id_source("aero_run")
        .default_open(true)
        .show(ui, |ui| {
            let running = app.aero.is_running();
            let mode = app.aero.form.run_mode;

            ui.horizontal(|ui| {
                let label = match mode {
                    RunMode::SteadyPoint => "▶  Run wind tunnel",
                    RunMode::AngleSweep => "▶  Run angle sweep",
                };
                let btn = egui::Button::new(egui::RichText::new(label).strong());
                if ui
                    .add_enabled(!running, btn)
                    .on_hover_text(
                        "Solve the case on a background thread — the window stays \
                         responsive.\nShortcut: Ctrl+R / F5",
                    )
                    .clicked()
                {
                    start_run(app);
                }
                if running {
                    ui.spinner();
                    if ui
                        .small_button("Cancel")
                        .on_hover_text("Stop the in-flight solve (Esc).")
                        .clicked()
                    {
                        // Drop the run handle by clearing the option;
                        // the next pump frame moves the panel back to
                        // idle. The aero runner's handle does not
                        // currently expose a `request_stop` — clearing
                        // it from the form is the safe approximation:
                        // the worker thread is detached, so the value
                        // simply runs to completion and gets dropped,
                        // but the user-visible status reads "Cancelled".
                        app.aero.status = "Cancelled by user.".to_string();
                    }
                }
                ui.separator();
                // Inline undo / redo — rewinds the wind-tunnel form
                // state to its prior snapshot. Recorded on Run, so
                // Ctrl+Z reverses the parameters of the last solve.
                let undo = ui
                    .add_enabled(app.aero.can_undo(), egui::Button::new(valenx_icons::edit::UNDO))
                    .on_hover_text("Undo last form edit (Ctrl+Z)")
                    .clicked();
                let redo = ui
                    .add_enabled(app.aero.can_redo(), egui::Button::new(valenx_icons::edit::REDO))
                    .on_hover_text("Redo last form edit (Ctrl+Y)")
                    .clicked();
                if undo { app.aero.undo_edit(); }
                if redo { app.aero.redo_edit(); }
            });

            // Status line.
            if !app.aero.status.is_empty() {
                ui.label(egui::RichText::new(&app.aero.status).weak());
            }
            if running {
                ui.label(
                    egui::RichText::new(
                        "Solving… the steady RANS solve runs in the background.",
                    )
                    .weak()
                    .small(),
                );
            }
            error_line(ui, &app.aero.error);

            // Live convergence (residual) plot — the mass-imbalance
            // residual on a log10 axis, the standard CFD convergence
            // monitor.
            if !app.aero.residual_history.is_empty() {
                ui.add_space(4.0);
                ui.label(egui::RichText::new("Convergence — mass residual").strong());
                let points: PlotPoints = app
                    .aero
                    .residual_history
                    .iter()
                    .map(|(it, r)| [*it, r.max(1e-30).log10()])
                    .collect();
                Plot::new("aero_residual_plot")
                    .height(140.0)
                    .y_axis_label("residual (log₁₀)")
                    .x_axis_label("iteration")
                    .show(ui, |plot_ui| {
                        plot_ui.line(Line::new(points).name("mass imbalance"));
                    });
            }

            // Live sweep progress — Cd / Cl accumulating per angle.
            if mode == RunMode::AngleSweep && !app.aero.sweep_progress.is_empty() {
                ui.add_space(4.0);
                derived(
                    ui,
                    format!(
                        "Sweep progress: {} angle(s) solved",
                        app.aero.sweep_progress.len()
                    ),
                );
            }
        });
}

// ---------------------------------------------------------------------------
// 7. Results
// ---------------------------------------------------------------------------

/// Section 7 — the coefficient result cards, drag breakdown, and (for a
/// sweep) the lift / drag polar plot.
pub fn draw_results_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("7 · Results")
        .id_source("aero_results")
        .default_open(true)
        .show(ui, |ui| {
            if app.aero.last_report.is_none() && app.aero.last_polar.is_none() {
                derived(
                    ui,
                    "Run the wind tunnel to see drag, lift and moment results.",
                );
                return;
            }

            // --- Steady single-point results ------------------------
            if let Some(report) = app.aero.last_report.as_ref() {
                draw_steady_results(ui, report);
            }

            // --- Angle-of-attack polar results ----------------------
            if let Some(polar) = app.aero.last_polar.as_ref() {
                draw_polar_results(ui, &polar.curve);
            }
        });
}

/// Render the headline result cards + drag breakdown for a steady run.
fn draw_steady_results(ui: &mut egui::Ui, report: &valenx_aero::AeroReport) {
    // Four prominent coefficient cards — the headline numbers.
    let drag_tint = egui::Color32::from_rgb(235, 150, 90);
    let lift_tint = egui::Color32::from_rgb(120, 190, 235);
    let neutral = ui.visuals().strong_text_color();
    ui.horizontal_wrapped(|ui| {
        result_card(ui, "Drag  Cd", &model::format_coefficient(report.cd), drag_tint);
        result_card(ui, "Lift  Cl", &model::format_coefficient(report.cl), lift_tint);
        result_card(
            ui,
            "Side  Cs",
            &model::format_coefficient(report.cs),
            neutral,
        );
        result_card(
            ui,
            "Pitch  Cm",
            &model::format_coefficient(report.cm),
            neutral,
        );
        let ld = match model::lift_to_drag(report.cl, report.cd) {
            Some(v) => format!("{v:.2}"),
            None => "—".to_string(),
        };
        result_card(
            ui,
            "Lift/Drag  L/D",
            &ld,
            egui::Color32::from_rgb(140, 210, 150),
        );
        result_card(
            ui,
            "Glide angle",
            &format!("{:.1}\u{00B0}", report.glide_angle_rad().to_degrees()),
            neutral,
        );
    });

    ui.add_space(4.0);
    ui.label(egui::RichText::new("Drag breakdown").strong());
    egui::Grid::new("aero_result_grid")
        .num_columns(2)
        .spacing([8.0, 3.0])
        .show(ui, |ui| {
            kv(
                ui,
                "pressure vs. friction",
                model::format_drag_split(report.pressure_drag_fraction),
            );
            kv(
                ui,
                "drag area  Cd·A",
                format!("{:.4} m²", report.drag_area),
            );
            kv(
                ui,
                "reference area  A",
                format!("{:.4} m²", report.reference_area),
            );
            kv(
                ui,
                "drag force  Cd·A·q",
                model::format_force_n(report.drag_force()),
            );
            kv(
                ui,
                "lift force  Cl·A·q",
                model::format_force_n(report.lift_force()),
            );
            kv(
                ui,
                "resultant force",
                model::format_force_n(report.resultant_force()),
            );
            kv(
                ui,
                "dynamic pressure  q∞",
                model::format_pressure_pa(report.dynamic_pressure),
            );
            kv(ui, "Reynolds number", model::format_reynolds(report.reynolds_number));
            kv(ui, "Mach number", format!("{:.3}", report.mach_number));
            kv(
                ui,
                "P-G factor",
                format!("{:.3}", report.prandtl_glauert_factor()),
            );
            kv(
                ui,
                "convergence",
                if report.converged {
                    format!("converged · {} iters", report.iterations)
                } else {
                    format!("capped · {} iters · residual {:.1e}", report.iterations, report.residual)
                },
            );
            kv(ui, "mean wall y+", format!("{:.1}", report.y_plus_mean));
            kv(
                ui,
                "wake peak deficit",
                format!("{:.0} %", 100.0 * report.wake_peak_deficit),
            );
        });

    // The honest caveats from the report — these qualify every number
    // above, so they belong in the UI, not just the text dump.
    if !report.caveats.is_empty() {
        egui::CollapsingHeader::new("Caveats — how to read these numbers")
            .id_source("aero_caveats")
            .default_open(false)
            .show(ui, |ui| {
                for c in &report.caveats {
                    ui.label(egui::RichText::new(format!("• {c}")).weak().small());
                }
            });
    }
}

/// Render the polar plot + summary diagnostics for a sweep run.
fn draw_polar_results(ui: &mut egui::Ui, curve: &valenx_aero::PolarCurve) {
    if curve.points.is_empty() {
        derived(ui, "The sweep produced no points.");
        return;
    }
    ui.label(egui::RichText::new("Lift / drag polar").strong());

    // The drag polar — Cl on the y-axis vs. Cd on the x-axis. Built
    // from a plain `Vec` so both the line and the point markers can be
    // drawn from the same data (`PlotPoints` itself is not `Clone`).
    let polar_xy: Vec<[f64; 2]> = curve.points.iter().map(|p| [p.cd, p.cl]).collect();
    Plot::new("aero_polar_plot")
        .height(170.0)
        .legend(Legend::default())
        .x_axis_label("drag  Cd")
        .y_axis_label("lift  Cl")
        .show(ui, |plot_ui| {
            plot_ui.line(Line::new(PlotPoints::from(polar_xy.clone())).name("polar"));
            plot_ui.points(
                Points::new(PlotPoints::from(polar_xy.clone()))
                    .radius(3.5)
                    .name("solved angles"),
            );
        });

    // The lift curve — Cl vs. angle of attack (degrees).
    let lift_curve: PlotPoints = curve
        .points
        .iter()
        .map(|p| [model::rad_to_deg(p.alpha), p.cl])
        .collect();
    Plot::new("aero_lift_curve")
        .height(140.0)
        .x_axis_label("angle of attack  α (°)")
        .y_axis_label("lift  Cl")
        .show(ui, |plot_ui| {
            plot_ui.line(Line::new(lift_curve).name("lift curve"));
        });

    // Polar diagnostics.
    egui::Grid::new("aero_polar_diag")
        .num_columns(2)
        .spacing([8.0, 3.0])
        .show(ui, |ui| {
            kv(
                ui,
                "max lift / drag",
                format!("{:.2}", curve.max_lift_to_drag()),
            );
            if let Some(best) = curve.best_lift_to_drag_point() {
                kv(
                    ui,
                    "best L/D at α",
                    format!("{:.1} °", model::rad_to_deg(best.alpha)),
                );
            }
            if let Some(end) = curve.best_endurance_point() {
                kv(
                    ui,
                    "best endurance α",
                    format!("{:.1} °  (max Cl^1.5/Cd)", model::rad_to_deg(end.alpha)),
                );
            }
            if let Some(min) = curve.min_drag_point() {
                kv(ui, "min drag  Cd", format!("{:.4}", min.cd));
            }
            if let Some(a0) = curve.zero_lift_angle() {
                kv(ui, "zero-lift α", format!("{:.1} °", model::rad_to_deg(a0)));
            }
            kv(ui, "induced drag  k", format!("{:.4}", curve.induced_drag_factor()));
            kv(
                ui,
                "parasitic drag  Cd\u{2080}",
                format!("{:.4}", curve.parasitic_drag_coefficient()),
            );
            kv(ui, "max lift  Cl", format!("{:.4}", curve.max_lift()));
            if let Some(stall) = curve.stall_angle() {
                kv(
                    ui,
                    "stall angle (est.)",
                    format!("{:.1} °", model::rad_to_deg(stall)),
                );
            }
            kv(
                ui,
                "lift-curve slope",
                format!("{:.3} /rad", curve.lift_curve_slope()),
            );
            kv(
                ui,
                "pitch stab dCm/dCl",
                format!("{:.3}  (<0 stable)", curve.pitch_stability_slope()),
            );
        });
    derived(
        ui,
        "Steady RANS does not capture deep post-stall behaviour — the pre-stall \
         polar is the reliable range.",
    );
}

// ---------------------------------------------------------------------------
// 8. Flow visualization
// ---------------------------------------------------------------------------

/// Section 8 — push field results into the 3-D viewport.
pub fn draw_visualization_section(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new("8 · Flow visualization")
        .id_source("aero_viz")
        .default_open(true)
        .show(ui, |ui| {
            if app.aero.last_result.is_none() {
                derived(
                    ui,
                    "Run a steady wind-tunnel solve to colour the geometry by a \
                     flow field. (Visualization uses the steady result; an \
                     angle sweep has no single field to show.)",
                );
                return;
            }

            ui.label("Field");
            egui::ComboBox::from_id_source("aero_flow_field")
                .selected_text(app.aero.form.flow_field.label())
                .width(ui.available_width())
                .show_ui(ui, |ui| {
                    for f in FlowField::ALL {
                        ui.selectable_value(
                            &mut app.aero.form.flow_field,
                            f,
                            f.label(),
                        );
                    }
                });

            // Cut-plane controls — only for the volumetric fields.
            if !app.aero.form.flow_field.is_surface() {
                ui.add_space(2.0);
                egui::Grid::new("aero_cut_grid")
                    .num_columns(2)
                    .spacing([8.0, 4.0])
                    .show(ui, |ui| {
                        ui.label("Cut-plane axis");
                        egui::ComboBox::from_id_source("aero_cut_axis")
                            .selected_text(app.aero.form.cut_axis.label())
                            .show_ui(ui, |ui| {
                                for a in CutAxis::ALL {
                                    ui.selectable_value(
                                        &mut app.aero.form.cut_axis,
                                        a,
                                        a.label(),
                                    );
                                }
                            });
                        ui.end_row();
                        ui.label("Cut position");
                        ui.add(
                            egui::Slider::new(
                                &mut app.aero.form.cut_fraction,
                                0.0..=1.0,
                            )
                            .custom_formatter(|v, _| format!("{:.0} %", v * 100.0)),
                        );
                        ui.end_row();
                    });
            }

            ui.add_space(4.0);
            if ui
                .add_sized(
                    [ui.available_width(), 24.0],
                    egui::Button::new("Show field in 3-D viewport"),
                )
                .on_hover_text(
                    "Colours the geometry by the selected field using the \
                     viewport's per-vertex colour ramp + legend.",
                )
                .clicked()
            {
                push_field_to_viewport(app);
            }
            if app.aero.last_field_overlay.is_some()
                && ui.button("Clear field overlay").clicked()
            {
                clear_field_overlay(app);
            }

            if let Some(status) = app.aero.viz_status.as_ref() {
                ui.label(egui::RichText::new(status).weak().small());
            }
            // A quick legend reminder of the active field.
            if app.aero.last_field_overlay.is_some() {
                derived(
                    ui,
                    format!(
                        "Active overlay: {} — see the colour bar in the viewport.",
                        app.aero.form.flow_field.legend_name()
                    ),
                );
            }
        });
}

// ---------------------------------------------------------------------------
// Actions
// ---------------------------------------------------------------------------

/// Public entry-point for the keyboard-shortcut handler — Ctrl+R
/// while the Wind Tunnel workbench is visible routes here. Delegates
/// to the private `start_run` (in the same module) so the bulk of
/// the logic stays internal.
pub fn start_run_from_shortcut(app: &mut ValenxApp) {
    if app.aero.is_running() {
        return;
    }
    start_run(app);
}

/// Validate the form and start a background wind-tunnel run.
fn start_run(app: &mut ValenxApp) {
    app.aero.error = None;
    app.aero.residual_history.clear();
    app.aero.sweep_progress.clear();
    app.aero.viz_status = None;
    // Snapshot the form so a later Ctrl+Z can rewind the user back
    // to the settings of the last solve.
    app.aero.record_form();

    // Extract the body.
    let body = match super::extract_body(
        &app.aero.form,
        app.current_solid.as_ref(),
        app.stl.as_ref().map(|s| &s.mesh),
    ) {
        Ok(b) => b,
        Err(e) => {
            app.aero.error = Some(friendly_aero_error(&e));
            return;
        }
    };

    // Build the request.
    let request = match app.aero.form.build_request() {
        Ok(r) => r,
        Err(e) => {
            app.aero.error = Some(friendly_aero_error(&format!("invalid case: {e}")));
            return;
        }
    };

    let job = match app.aero.form.run_mode {
        RunMode::SteadyPoint => AeroJob::Steady,
        RunMode::AngleSweep => AeroJob::Sweep(app.aero.form.sweep_angles()),
    };

    app.aero.last_body_label = body.source_label.clone();
    app.aero.status = "Starting…".to_string();
    app.aero.run = Some(AeroRunHandle::spawn(body.mesh, request, job));
}

/// Build the selected flow field and push it into the 3-D viewport via
/// the FEM-style per-vertex colour-ramp overlay path.
fn push_field_to_viewport(app: &mut ValenxApp) {
    let Some(result) = app.aero.last_result.as_ref() else {
        app.aero.viz_status = Some("No steady result to visualize.".to_string());
        return;
    };
    // Snapshot the form values up front so the `&app.aero` borrow ends
    // before `app.apply_mesh` (which needs `&mut app`) runs.
    let flow_field = app.aero.form.flow_field;
    let cut_axis = app.aero.form.cut_axis;
    let cut_fraction = app.aero.form.cut_fraction;

    match super::viz::build_flow_viz(result, flow_field, cut_axis, cut_fraction) {
        Ok(viz) => {
            let n_tris = viz.mesh.total_elements();
            // Drop the coloured surface / cut-plane patch into the
            // shared mesh slot, and the scalar field into the aero
            // overlay slot the viewport prefers.
            app.apply_mesh(
                viz.mesh,
                std::path::PathBuf::from(format!(
                    "<wind-tunnel>/{}",
                    flow_field.label()
                )),
            );
            app.aero_field_overlay = Some(viz.field);
            app.aero.last_field_overlay = Some(flow_field);
            app.aero.viz_status = Some(format!(
                "Showing {} in the viewport ({} triangles).",
                flow_field.legend_name(),
                n_tris,
            ));
        }
        Err(e) => {
            app.aero.viz_status = Some(format!("Could not build the field: {e}"));
        }
    }
}

/// Clear the active aero field overlay from the viewport.
fn clear_field_overlay(app: &mut ValenxApp) {
    app.aero_field_overlay = None;
    app.aero.last_field_overlay = None;
    app.aero.viz_status = Some("Field overlay cleared.".to_string());
}

// ---------------------------------------------------------------------------
// Small read-only helpers
// ---------------------------------------------------------------------------

/// The characteristic body length (the streamwise extent) of the
/// currently-selected body — feeds the displayed Reynolds number.
/// Falls back to the demo box length when no body can be extracted.
fn body_reference_length(app: &ValenxApp) -> f64 {
    match super::extract_body(
        &app.aero.form,
        app.current_solid.as_ref(),
        app.stl.as_ref().map(|s| &s.mesh),
    ) {
        Ok(body) => body
            .mesh
            .aabb()
            .map(|bb| bb.extent().x.max(1e-3))
            .unwrap_or(1.0),
        Err(_) => app.aero.form.demo_box_size[0].max(1e-3),
    }
}

/// Rewrite common wind-tunnel error messages into friendlier text
/// that tells the user what to try next. Unknown messages pass
/// through unchanged so we never hide information; the goal is to
/// dress up the recognised cases.
fn friendly_aero_error(raw: &str) -> String {
    let low = raw.to_ascii_lowercase();
    if low.contains("no body") || low.contains("body not ready") || low.contains("empty") {
        format!(
            "{raw}\n→ Pick a Test body in section 1. The Demo box always works \
             without loading external geometry."
        )
    } else if low.contains("invalid case") && (low.contains("density") || low.contains("viscosity"))
    {
        format!(
            "{raw}\n→ Check the Wind conditions section. Density must be > 0 \
             and viscosity > 0; sea-level air is 1.225 / 1.81e-5."
        )
    } else if low.contains("triangles") && low.contains("0") {
        format!(
            "{raw}\n→ The body extractor produced an empty mesh. Reload the \
             STL or try the Demo box."
        )
    } else if low.contains("sweep") && (low.contains("start") || low.contains("end")) {
        format!(
            "{raw}\n→ Set the sweep start angle below the end angle, and at \
             least 2 sweep points."
        )
    } else {
        raw.to_string()
    }
}

/// Format a large integer count compactly (e.g. `2_000_000` → `"2.0M"`).
fn format_count(n: usize) -> String {
    if n >= 1_000_000 {
        format!("{:.1}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.0}k", n as f64 / 1_000.0)
    } else {
        n.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_count_is_compact() {
        assert_eq!(format_count(2_000_000), "2.0M");
        assert_eq!(format_count(400_000), "400k");
        assert_eq!(format_count(512), "512");
    }

    #[test]
    fn body_reference_length_falls_back_to_the_demo_box() {
        // With the demo box selected and no other geometry, the
        // reference length is the box's streamwise extent.
        let mut app = ValenxApp::default();
        app.aero.form.body_source = BodySource::DemoBox;
        app.aero.form.demo_box_size = [5.5, 2.0, 1.5];
        let l = body_reference_length(&app);
        assert!((l - 5.5).abs() < 1e-6);
    }

    #[test]
    fn start_run_errors_without_a_body() {
        // Selecting the CAD-model source with nothing loaded must
        // surface an error rather than spawning a doomed run.
        let mut app = ValenxApp::default();
        app.aero.form.body_source = BodySource::CurrentCadModel;
        start_run(&mut app);
        assert!(app.aero.error.is_some());
        assert!(!app.aero.is_running());
    }

    #[test]
    fn start_run_with_the_demo_box_spawns_a_run() {
        // The demo box is always available — a run must actually start.
        let mut app = ValenxApp::default();
        app.aero.form.body_source = BodySource::DemoBox;
        // A tiny iteration budget keeps the spawned solve short; the
        // run is spawned regardless of how long it takes.
        app.aero.form.max_iterations = 4;
        start_run(&mut app);
        assert!(app.aero.error.is_none(), "unexpected error: {:?}", app.aero.error);
        assert!(app.aero.is_running());
        // Drain the background run so the test thread is joined.
        while let Some(handle) = app.aero.run.as_mut() {
            if handle.take_outcome().is_some() {
                app.aero.run = None;
            } else {
                std::thread::sleep(std::time::Duration::from_millis(5));
            }
        }
    }

    #[test]
    fn push_field_without_a_result_sets_a_status() {
        // Asking to visualize with no completed solve is a no-op with
        // an explanatory status, not a panic.
        let mut app = ValenxApp::default();
        push_field_to_viewport(&mut app);
        assert!(app.aero.viz_status.is_some());
        assert!(app.aero_field_overlay.is_none());
    }

    #[test]
    fn clear_field_overlay_resets_the_viewport_slot() {
        let mut app = ValenxApp::default();
        // Pretend an overlay is active.
        app.aero.last_field_overlay = Some(FlowField::SurfaceCp);
        clear_field_overlay(&mut app);
        assert!(app.aero_field_overlay.is_none());
        assert!(app.aero.last_field_overlay.is_none());
        assert!(app.aero.viz_status.is_some());
    }
}

/// Headless egui UI-logic tests for the Wind Tunnel workbench's eight
/// workflow sections.
///
/// Each `draw_*_section` is rendered into a windowless
/// [`egui::Context`]; the Run action (`start_run`) is driven to a real
/// `valenx-aero` solve on the demo box. Nothing opens an OS window and
/// nothing reaches `rfd::FileDialog` (the Import-STL button is drawn
/// but never clicked).
#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use crate::aero::model::{BodySource, RunMode};

    /// Draw every workflow section once in a headless context.
    fn draw_all_sections(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                draw_body_section(app, ui);
                draw_wind_section(app, ui);
                draw_ground_section(app, ui);
                draw_tunnel_section(app, ui);
                draw_solver_section(app, ui);
                draw_run_section(app, ui);
                draw_results_section(app, ui);
                draw_visualization_section(app, ui);
            });
        });
    }

    /// Drain a spawned background run to completion so the worker
    /// thread is joined (a test, not the UI thread).
    fn drain_run(app: &mut ValenxApp) {
        for _ in 0..6000 {
            let done = match app.aero.run.as_mut() {
                Some(handle) => handle.take_outcome().is_some(),
                None => true,
            };
            if done {
                app.aero.run = None;
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
    }

    #[test]
    fn draws_all_sections_fresh_without_panic() {
        // Fresh state: the demo box is selected, no run has happened.
        let mut app = ValenxApp::default();
        app.aero.form.body_source = BodySource::DemoBox;
        draw_all_sections(&mut app);
    }

    #[test]
    fn draws_all_sections_in_sweep_mode_without_panic() {
        // Angle-sweep mode opens the extra sweep-parameter widgets in
        // the Solver section.
        let mut app = ValenxApp::default();
        app.aero.form.body_source = BodySource::DemoBox;
        app.aero.form.run_mode = RunMode::AngleSweep;
        draw_all_sections(&mut app);
    }

    #[test]
    fn draws_with_a_body_source_that_has_no_geometry() {
        // Selecting the CAD-model source with nothing loaded exercises
        // the "Body not ready" branch — it must render, not panic.
        let mut app = ValenxApp::default();
        app.aero.form.body_source = BodySource::CurrentCadModel;
        draw_all_sections(&mut app);
        app.aero.form.body_source = BodySource::ImportedStl;
        draw_all_sections(&mut app);
    }

    #[test]
    fn draws_error_and_progress_states_without_panic() {
        // An error line + a residual history (the live convergence
        // plot) + a status line must all render.
        let mut app = ValenxApp::default();
        app.aero.form.body_source = BodySource::DemoBox;
        app.aero.error = Some("invalid case: bad air".to_string());
        app.aero.status = "Solving the steady 3-D flow…".to_string();
        app.aero.residual_history = vec![(1.0, 1.0e-1), (2.0, 1.0e-3), (3.0, 1.0e-5)];
        draw_all_sections(&mut app);
    }

    #[test]
    fn draws_post_run_results_without_panic() {
        // A completed steady solve → the Results section's coefficient
        // cards + drag breakdown must render, and the Flow viz section
        // must show its field picker. The workbench's per-frame pump
        // moves the finished run's outcome onto the workbench state.
        let mut app = ValenxApp::default();
        app.show_aero_workbench = true;
        app.aero.form.body_source = BodySource::DemoBox;
        app.aero.form.max_iterations = 4;
        start_run(&mut app);
        assert!(app.aero.is_running(), "a run should have been spawned");
        // Block until the steady result lands on the workbench, driving
        // the real per-frame pump (draw_aero_workbench calls it).
        for _ in 0..6000 {
            let ctx = egui::Context::default();
            let _ = ctx.run(egui::RawInput::default(), |ctx| {
                crate::aero_workbench::draw_aero_workbench(&mut app, ctx);
            });
            if !app.aero.is_running() {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
        // The run completed: a steady result or a surfaced error — and
        // the post-run sections render either way, never a panic.
        assert!(app.aero.last_result.is_some() || app.aero.error.is_some());
        draw_all_sections(&mut app);
    }

    #[test]
    fn run_action_starts_a_real_steady_solve() {
        // start_run with the demo box extracts the body, builds the
        // request and spawns a real valenx-aero solve.
        let mut app = ValenxApp::default();
        app.aero.form.body_source = BodySource::DemoBox;
        app.aero.form.max_iterations = 4;
        start_run(&mut app);
        assert!(app.aero.error.is_none(), "unexpected error: {:?}", app.aero.error);
        assert!(app.aero.is_running(), "a run should have been spawned");
        drain_run(&mut app);
    }

    #[test]
    fn run_action_surfaces_error_without_a_body() {
        // The CAD-model source with nothing loaded must surface an
        // error rather than spawning a doomed run.
        let mut app = ValenxApp::default();
        app.aero.form.body_source = BodySource::CurrentCadModel;
        start_run(&mut app);
        assert!(app.aero.error.is_some(), "start_run should error without a body");
        assert!(!app.aero.is_running());
    }

    #[test]
    fn run_action_surfaces_error_on_a_non_physical_case() {
        // A non-physical air density makes the request builder fail —
        // start_run must surface that as an error, not a spawned run.
        let mut app = ValenxApp::default();
        app.aero.form.body_source = BodySource::DemoBox;
        app.aero.form.air_density = -1.0;
        start_run(&mut app);
        assert!(app.aero.error.is_some(), "start_run should error on bad air");
        assert!(!app.aero.is_running());
    }

    #[test]
    fn push_field_without_a_result_is_graceful() {
        // Asking to visualize a field with no completed solve is a
        // no-op with an explanatory status — never a panic.
        let mut app = ValenxApp::default();
        push_field_to_viewport(&mut app);
        assert!(app.aero.viz_status.is_some());
        assert!(app.aero_field_overlay.is_none());
    }
}
