//! The right-side **Rotor / Drone (BEMT) Workbench** panel — a native
//! front-end over the in-house `valenx-rotor` blade-element-momentum-theory
//! engine.
//!
//! Mirrors the other workbenches (`fem_workbench`, `blackhole_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_rotor_workbench`], toggled from the View menu and
//! openable by the agent bridge under the workbench id `"rotor"` (see
//! [`crate::project_tabs::TabKind`] / [`crate::agent_commands`]).
//!
//! The user sets the rotor geometry — blade count, tip & hub radius, and an
//! editable radial station table (radius / chord / twist per strip) — and the
//! operating point (rpm, axial freestream `V`, air density). **Compute**
//! builds a [`valenx_rotor::Rotor`] via [`Rotor::from_slices`] and runs
//! [`Rotor::solve`], then reads back thrust, torque, shaft power, propeller
//! efficiency, hover figure of merit, disk area and the per-element converged
//! aerodynamics (inflow angle, α, induction, Cl/Cd, dT/dr, dQ/dr).
//!
//! Honesty: every geometry / operating-point error (non-monotonic stations,
//! hub ≥ tip, zero/negative/non-finite rpm·V·ρ, a station whose inflow root
//! cannot be bracketed) surfaces as the engine's [`RotorError`] verbatim in a
//! `⚠ …` status line — the workbench never invents a number. The model is
//! research/educational-grade 1-D strip theory with the standard corrections
//! (Prandtl loss, Glauert/Buhl high-induction); absolute magnitudes are not
//! calibrated against measured propeller data. No new physics — wiring and
//! presentation over the validated engine only.

use std::fmt::Write as _;

use eframe::egui;

use valenx_rotor::{Rotor, RotorError, RotorPerformance, SEA_LEVEL_AIR_DENSITY};

use crate::ValenxApp;

/// One editable radial blade station in the workbench table: radius, chord and
/// twist (the twist is held in **degrees** for the UI and converted to radians
/// only when the [`Rotor`] is built, so the user edits a familiar unit).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StationRow {
    /// Radial position `r` of the station (m).
    pub radius_m: f64,
    /// Local chord `c(r)` (m).
    pub chord_m: f64,
    /// Local geometric twist / pitch `β(r)` (degrees, converted to radians on
    /// solve).
    pub twist_deg: f64,
}

/// Persistent state for the Rotor / Drone (BEMT) workbench.
pub struct RotorWorkbenchState {
    /// Number of blades `n` (≥ 1).
    pub blade_count: u32,
    /// Tip radius `R` (m).
    pub tip_radius_m: f64,
    /// Hub (root) radius (m), `0 ≤ hub < R`.
    pub hub_radius_m: f64,
    /// Editable radial station table (radius / chord / twist°).
    pub stations: Vec<StationRow>,
    /// Rotor speed (rev/min).
    pub rpm: f64,
    /// Axial freestream speed `V` (m/s; `0` is hover).
    pub freestream_v: f64,
    /// Air density `ρ` (kg/m³).
    pub air_density: f64,
    /// The last successful solve, if any (cleared when a solve errors).
    pub perf: Option<RotorPerformance>,
    /// Short status / error line (e.g. `⚠ <error>` when a solve fails, or a
    /// one-line success summary). Empty before the first compute.
    pub status: String,
}

impl Default for RotorWorkbenchState {
    /// A sensible 5-station, 2-blade, 0.15 m-radius tapered + twisted
    /// propeller — the same example geometry the `valenx-rotor` crate tests
    /// use (R = 0.15 m, hub = 0.02 m; twist washing out 25° → 8°). Forward
    /// flight at 5000 rpm / 5 m/s, sea-level air.
    fn default() -> Self {
        Self {
            blade_count: 2,
            tip_radius_m: 0.15,
            hub_radius_m: 0.02,
            stations: vec![
                StationRow {
                    radius_m: 0.03,
                    chord_m: 0.025,
                    twist_deg: 25.0,
                },
                StationRow {
                    radius_m: 0.06,
                    chord_m: 0.022,
                    twist_deg: 18.0,
                },
                StationRow {
                    radius_m: 0.09,
                    chord_m: 0.018,
                    twist_deg: 13.0,
                },
                StationRow {
                    radius_m: 0.12,
                    chord_m: 0.014,
                    twist_deg: 10.0,
                },
                StationRow {
                    radius_m: 0.15,
                    chord_m: 0.008,
                    twist_deg: 8.0,
                },
            ],
            rpm: 5000.0,
            freestream_v: 5.0,
            air_density: SEA_LEVEL_AIR_DENSITY,
            perf: None,
            status: String::new(),
        }
    }
}

impl RotorWorkbenchState {
    /// Build a validated [`Rotor`] from the current geometry, converting the
    /// per-station twist from degrees to radians.
    ///
    /// # Errors
    /// Propagates any [`RotorError`] from [`Rotor::from_slices`] (mismatched
    /// lengths can't occur here — the three vectors are built in lock-step —
    /// but bad geometry: too few stations, non-monotonic radii, hub ≥ tip,
    /// non-positive chord/radius, do).
    pub fn build_rotor(&self) -> Result<Rotor, RotorError> {
        let radii: Vec<f64> = self.stations.iter().map(|s| s.radius_m).collect();
        let chords: Vec<f64> = self.stations.iter().map(|s| s.chord_m).collect();
        let twists: Vec<f64> = self
            .stations
            .iter()
            .map(|s| s.twist_deg.to_radians())
            .collect();
        Rotor::from_slices(
            self.blade_count,
            self.tip_radius_m,
            self.hub_radius_m,
            &radii,
            &chords,
            &twists,
        )
    }

    /// Build + solve the rotor at the current operating point.
    ///
    /// # Errors
    /// Propagates [`RotorError`] from [`Self::build_rotor`] (geometry) or
    /// [`Rotor::solve`] (operating point / non-convergence).
    pub fn solve(&self) -> Result<RotorPerformance, RotorError> {
        let rotor = self.build_rotor()?;
        rotor.solve(self.rpm, self.freestream_v, self.air_density)
    }
}

/// A one-line success summary for the status row (the full numbers render in
/// the readout grid below). Kept terse so it reads at a glance.
fn success_summary(perf: &RotorPerformance) -> String {
    let mut s = String::new();
    let _ = write!(
        s,
        "✔ T = {:.3} N · P = {:.1} W",
        perf.thrust_n, perf.power_w
    );
    if perf.freestream_v_m_s > 0.0 {
        let _ = write!(s, " · η = {:.3}", perf.efficiency);
    } else {
        let _ = write!(s, " · FM = {:.3}", perf.figure_of_merit);
    }
    s
}

/// Draw the Rotor / Drone (BEMT) workbench (a no-op unless toggled on via
/// View → Rotor / Drone (BEMT)). Mirrors
/// [`crate::blackhole_workbench::draw_blackhole_workbench`].
pub fn draw_rotor_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_rotor_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_rotor_workbench",
        "Rotor / Drone (BEMT)",
        rotor_workbench_body,
    );
    if close {
        app.show_rotor_workbench = false;
    }
}

/// The workbench body: rotor geometry (blade count, tip/hub radius, the
/// editable station table) + operating point (rpm, V, ρ), a **Compute**
/// button, and the integrated-performance readout + per-element table. Split
/// out so [`workbench_shell`](crate::workbench_chrome::workbench_shell) (and
/// any future dockable layout) can call it with the same `(app, ui)` shape.
fn rotor_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new("propeller / rotor BEMT · valenx-rotor")
            .weak()
            .small(),
    );
    ui.separator();

    let mut do_compute = false;
    let mut add_station = false;
    let mut remove_station: Option<usize> = None;
    {
        let s = &mut app.rotor;

        // Associate each numeric `DragValue` with its caption (the grid's label
        // cell) via `labelled_by`, so the spin button carries that caption as
        // its accessibility / UI-Automation Name. egui clears a DragValue's own
        // Name, so without the association the control is anonymous to a screen
        // reader / AI driver; the hover text gives a mouse user the same caption.
        egui::Grid::new("rotor_geometry")
            .num_columns(2)
            .show(ui, |ui| {
                let n = ui.label("blade count n");
                ui.add(
                    egui::DragValue::new(&mut s.blade_count)
                        .speed(1)
                        .range(1..=12),
                )
                .labelled_by(n.id)
                .on_hover_text("Number of blades n");
                ui.end_row();

                let r = ui.label("tip radius R (m)");
                ui.add(
                    egui::DragValue::new(&mut s.tip_radius_m)
                        .speed(0.01)
                        .range(0.001..=50.0),
                )
                .labelled_by(r.id)
                .on_hover_text("Tip radius R (m)");
                ui.end_row();

                let hub = ui.label("hub radius (m)");
                ui.add(
                    egui::DragValue::new(&mut s.hub_radius_m)
                        .speed(0.005)
                        .range(0.0..=50.0),
                )
                .labelled_by(hub.id)
                .on_hover_text("Hub (root) radius (m)");
                ui.end_row();
            });

        ui.add_space(4.0);
        ui.label(egui::RichText::new("Radial stations (root → tip)").strong());
        ui.label(
            egui::RichText::new("≥ 2 rows, strictly increasing radius, all within [hub, tip].")
                .weak()
                .small(),
        );

        egui::Grid::new("rotor_stations")
            .num_columns(4)
            .striped(true)
            .show(ui, |ui| {
                // The column-header labels double as the per-cell accessibility
                // captions: every station-row `DragValue` is `labelled_by` its
                // column header (so a screen reader / AI reads "radius (m)" etc.
                // for an otherwise-anonymous spin button), and also carries an
                // explicit hover Name including the row index.
                let h_radius = ui.label(egui::RichText::new("radius (m)").strong());
                let h_chord = ui.label(egui::RichText::new("chord (m)").strong());
                let h_twist = ui.label(egui::RichText::new("twist (°)").strong());
                ui.label("");
                ui.end_row();

                for (i, st) in s.stations.iter_mut().enumerate() {
                    let n = i + 1;
                    ui.add(
                        egui::DragValue::new(&mut st.radius_m)
                            .speed(0.005)
                            .range(0.0..=50.0),
                    )
                    .labelled_by(h_radius.id)
                    .on_hover_text(format!("Station {n} radius (m)"));
                    ui.add(
                        egui::DragValue::new(&mut st.chord_m)
                            .speed(0.002)
                            .range(0.0001..=10.0),
                    )
                    .labelled_by(h_chord.id)
                    .on_hover_text(format!("Station {n} chord (m)"));
                    ui.add(
                        egui::DragValue::new(&mut st.twist_deg)
                            .speed(0.5)
                            .range(-90.0..=90.0),
                    )
                    .labelled_by(h_twist.id)
                    .on_hover_text(format!("Station {n} twist (°)"));
                    // Always offer remove; the ≥ 2-rows invariant (BEMT
                    // integration needs at least two stations) is re-checked
                    // authoritatively when the click is applied below, because
                    // the live `stations.len()` isn't cheaply visible inside
                    // this `iter_mut()` borrow.
                    if ui.small_button("✕ remove").clicked() {
                        remove_station = Some(i);
                    }
                    ui.end_row();
                }
            });

        ui.horizontal(|ui| {
            if ui.button("＋ add station").clicked() {
                add_station = true;
            }
            if ui
                .button("Reset to default prop")
                .on_hover_text("Restore the 5-station 0.15 m example propeller.")
                .clicked()
            {
                *s = RotorWorkbenchState::default();
            }
        });

        ui.add_space(6.0);
        ui.label(egui::RichText::new("Operating point").strong());
        egui::Grid::new("rotor_operating_point")
            .num_columns(2)
            .show(ui, |ui| {
                let rpm = ui.label("rotor speed (rpm)");
                ui.add(
                    egui::DragValue::new(&mut s.rpm)
                        .speed(50.0)
                        .range(0.1..=1.0e6),
                )
                .labelled_by(rpm.id)
                .on_hover_text("Rotor speed (rev/min)");
                ui.end_row();

                let v = ui.label("freestream V (m/s)");
                ui.add(
                    egui::DragValue::new(&mut s.freestream_v)
                        .speed(0.5)
                        .range(0.0..=400.0),
                )
                .labelled_by(v.id)
                .on_hover_text("Axial inflow speed. 0 = hover (reports figure of merit).");
                ui.end_row();

                let rho = ui.label("air density ρ (kg/m³)");
                ui.add(
                    egui::DragValue::new(&mut s.air_density)
                        .speed(0.01)
                        .range(0.001..=10.0),
                )
                .labelled_by(rho.id)
                .on_hover_text("Air density ρ (kg/m³)");
                ui.end_row();
            });

        ui.add_space(6.0);
        if ui
            .button(egui::RichText::new("Compute").strong())
            .on_hover_text("Build the rotor and solve the BEMT inflow at this operating point.")
            .clicked()
        {
            do_compute = true;
        }
    }

    // Apply the table edits requested above (outside the per-row borrow). A
    // new station is appended just outboard of the current last one, clamped
    // to the tip, so it lands in a plausible spot for the user to tweak.
    if let Some(i) = remove_station {
        if app.rotor.stations.len() > 2 && i < app.rotor.stations.len() {
            app.rotor.stations.remove(i);
        }
    }
    if add_station {
        let s = &mut app.rotor;
        let (last_r, last_c, last_t) = s
            .stations
            .last()
            .map(|st| (st.radius_m, st.chord_m, st.twist_deg))
            .unwrap_or((s.hub_radius_m, 0.01, 5.0));
        let new_r = (last_r + 0.01).min(s.tip_radius_m);
        s.stations.push(StationRow {
            radius_m: new_r,
            chord_m: last_c,
            twist_deg: last_t,
        });
    }

    // Run the (synchronous, cheap) solve outside the borrow above. Store the
    // performance on success, or surface the engine's error verbatim — never
    // fabricate numbers.
    if do_compute {
        match app.rotor.solve() {
            Ok(perf) => {
                app.rotor.status = success_summary(&perf);
                app.rotor.perf = Some(perf);
            }
            Err(e) => {
                app.rotor.status = format!("⚠ {e}");
                app.rotor.perf = None;
            }
        }
    }

    let s = &app.rotor;
    if !s.status.is_empty() {
        ui.add_space(6.0);
        let color = if s.status.starts_with('⚠') {
            egui::Color32::from_rgb(220, 120, 60)
        } else {
            egui::Color32::from_rgb(90, 180, 110)
        };
        ui.label(egui::RichText::new(&s.status).color(color).strong());
    }

    if let Some(perf) = &s.perf {
        ui.add_space(6.0);
        ui.label(egui::RichText::new("Integrated performance").strong());
        egui::Grid::new("rotor_results")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let row = |ui: &mut egui::Ui, k: &str, v: String| {
                    ui.label(k);
                    ui.label(v);
                    ui.end_row();
                };
                row(ui, "thrust T (N)", format!("{:.4}", perf.thrust_n));
                row(ui, "torque Q (N·m)", format!("{:.4}", perf.torque_nm));
                row(ui, "shaft power P (W)", format!("{:.3}", perf.power_w));
                row(
                    ui,
                    "efficiency η",
                    if perf.freestream_v_m_s > 0.0 {
                        format!("{:.4}", perf.efficiency)
                    } else {
                        "— (hover)".to_string()
                    },
                );
                row(
                    ui,
                    "figure of merit FM",
                    if perf.freestream_v_m_s == 0.0 {
                        format!("{:.4}", perf.figure_of_merit)
                    } else {
                        "— (forward flight)".to_string()
                    },
                );
                row(ui, "disk area A (m²)", format!("{:.5}", perf.disk_area_m2));
                row(ui, "Ω (rad/s)", format!("{:.2}", perf.omega_rad_s));
            });

        ui.add_space(6.0);
        ui.label(egui::RichText::new("Per-element results").strong());
        egui::ScrollArea::horizontal()
            .id_source("rotor_elements_scroll")
            .show(ui, |ui| {
                egui::Grid::new("rotor_elements")
                    .num_columns(7)
                    .striped(true)
                    .show(ui, |ui| {
                        for h in ["r (m)", "φ (°)", "α (°)", "a", "Cl", "Cd", "dT/dr (N/m)"] {
                            ui.label(egui::RichText::new(h).strong());
                        }
                        ui.end_row();
                        for e in &perf.elements {
                            ui.label(format!("{:.4}", e.radius_m));
                            ui.label(format!("{:.2}", e.inflow_angle_rad.to_degrees()));
                            ui.label(format!("{:.2}", e.alpha_rad.to_degrees()));
                            ui.label(format!("{:.3}", e.axial_induction));
                            ui.label(format!("{:.3}", e.cl));
                            ui.label(format!("{:.4}", e.cd));
                            ui.label(format!("{:.3}", e.dt_dr));
                            ui.end_row();
                        }
                    });
            });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_solves_to_finite_positive_loads() {
        // The default 5-station prop in forward flight must build + solve to
        // finite, positive thrust and power — the workbench's happy path.
        let s = RotorWorkbenchState::default();
        let perf = s.solve().expect("default prop should solve");
        assert!(perf.thrust_n.is_finite() && perf.thrust_n > 0.0);
        assert!(perf.power_w.is_finite() && perf.power_w > 0.0);
        assert_eq!(perf.elements.len(), s.stations.len());
        // Forward flight (V > 0): efficiency is the reported figure, physical.
        assert!(perf.efficiency > 0.0 && perf.efficiency < 1.0);
    }

    #[test]
    fn hover_reports_physical_figure_of_merit() {
        // V = 0 is hover: FM must be in (0, 1] (cannot beat the actuator-disk
        // limit) and efficiency is defined as 0.
        let mut s = RotorWorkbenchState::default();
        s.freestream_v = 0.0;
        s.rpm = 6000.0;
        let perf = s.solve().expect("hover should solve");
        assert_eq!(perf.efficiency, 0.0);
        assert!(perf.figure_of_merit > 0.0 && perf.figure_of_merit <= 1.0);
    }

    #[test]
    fn bad_geometry_surfaces_error_not_fake_numbers() {
        // Hub at/above the tip is unphysical — the engine must return an Err,
        // which the workbench renders as a `⚠ …` status (never invents loads).
        let mut s = RotorWorkbenchState::default();
        s.hub_radius_m = 0.2; // >= tip (0.15)
        assert!(s.solve().is_err());

        // Non-monotonic stations (a duplicate radius) also error.
        let mut s2 = RotorWorkbenchState::default();
        s2.stations[1].radius_m = s2.stations[0].radius_m;
        assert!(s2.solve().is_err());
    }

    #[test]
    fn out_of_domain_operating_point_errors() {
        // Zero rpm / zero density / negative V are rejected by solve().
        let mut s = RotorWorkbenchState::default();
        s.rpm = 0.0;
        assert!(s.solve().is_err());
        let mut s2 = RotorWorkbenchState::default();
        s2.air_density = 0.0;
        assert!(s2.solve().is_err());
        let mut s3 = RotorWorkbenchState::default();
        s3.freestream_v = -1.0;
        assert!(s3.solve().is_err());
    }

    #[test]
    fn success_summary_picks_efficiency_or_fm_by_regime() {
        let mut s = RotorWorkbenchState::default();
        let cruise = s.solve().unwrap();
        assert!(success_summary(&cruise).contains("η"));
        s.freestream_v = 0.0;
        let hover = s.solve().unwrap();
        assert!(success_summary(&hover).contains("FM"));
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    /// Draw the panel once in a headless egui context **with accesskit enabled**
    /// and return the emitted accessibility tree nodes — the same tree a screen
    /// reader / AI UI-Automation driver consumes. `accesskit` is re-exported by
    /// egui, so this needs no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_rotor_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn has_named_node(nodes: &[(NodeId, Node)], name: &str) -> bool {
        nodes.iter().any(|(_, n)| n.name() == Some(name))
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_rotor_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_rotor_workbench(&mut app, ctx);
        });
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_rotor_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // The geometry, station-table and operating-point DragValues are all
        // SpinButtons; each must be `labelled_by` a caption (egui clears a
        // DragValue's own Name), including the per-row station cells which are
        // associated with their column-header labels.
        let mut app = ValenxApp::default();
        app.show_rotor_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // 3 geometry + 3 operating-point + 3×(default 5 stations) = 21.
        assert!(
            spin_buttons.len() >= 12,
            "expected the rotor numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every rotor DragValue must be labelled_by a caption (AI-drivable name)"
        );

        for caption in ["blade count n", "tip radius R (m)", "rotor speed (rpm)"] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
        // The Compute button keeps its plain-text, invokable Name.
        assert!(
            nodes.iter().any(|(_, n)| n.role() == Role::Button
                && n.name().is_some_and(|s| s.contains("Compute"))),
            "the Compute button is a named, invokable node"
        );
    }
}
