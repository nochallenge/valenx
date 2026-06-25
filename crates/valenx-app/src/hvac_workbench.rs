//! HVAC workbench — duct sizing + pressure drop on `valenx-hvac`.
//!
//! A right-side calculation panel: pick a round or rectangular duct section,
//! set the run length and air velocity, and compute the **Darcy–Weisbach
//! pressure drop** plus the hydraulic diameter, cross-sectional area, and a
//! CFM-based recommended duct size. Pure calculation — fully headless-testable.

use eframe::egui;

use valenx_hvac::{duct::CrossSection, flow, pressure_drop};

use crate::ValenxApp;

#[derive(Clone, Copy, PartialEq, Eq)]
enum DuctShape {
    Round,
    Rect,
}

/// Persistent state for the HVAC workbench.
pub struct HvacWorkbenchState {
    shape: DuctShape,
    diameter_mm: f64,
    width_mm: f64,
    height_mm: f64,
    length_m: f64,
    velocity_ms: f64,
    friction: f64,
    cfm: f64,
    max_velocity_fpm: f64,
    results: Option<HvacResults>,
}

impl Default for HvacWorkbenchState {
    fn default() -> Self {
        Self {
            shape: DuctShape::Round,
            diameter_mm: 200.0,
            width_mm: 300.0,
            height_mm: 200.0,
            length_m: 10.0,
            velocity_ms: 5.0,
            friction: 0.02,
            cfm: 500.0,
            max_velocity_fpm: 700.0,
            results: None,
        }
    }
}

struct HvacResults {
    hydraulic_diameter_mm: f64,
    area_m2: f64,
    pressure_drop_pa: f64,
    duct_w_in: f64,
    duct_h_in: f64,
}

impl HvacWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`.
    /// Covers both duct shapes' dimensions (round `diameter` + rectangular
    /// `width`/`height`) so an agent can set either regardless of the selected
    /// shape; each string matches the form caption exactly.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "diameter (mm)",
            "width (mm)",
            "height (mm)",
            "length (m)",
            "velocity (m/s)",
            "friction factor",
            "flow (CFM)",
            "max vel (FPM)",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Fail-loud: an unknown caption or a value of the
    /// wrong type / out of range returns `Err(String)` — never a panic. Ranges
    /// mirror the form's `DragValue` clamps exactly.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        let ranged = |v: f64, lo: f64, hi: f64, what: &str| -> Result<f64, String> {
            if v.is_finite() && (lo..=hi).contains(&v) {
                Ok(v)
            } else {
                Err(format!("{what} must be in {lo}..={hi}, got {v}"))
            }
        };
        match name {
            "diameter (mm)" => {
                self.diameter_mm = ranged(value.as_f64()?, 10.0, 2000.0, "diameter")?
            }
            "width (mm)" => self.width_mm = ranged(value.as_f64()?, 10.0, 2000.0, "width")?,
            "height (mm)" => self.height_mm = ranged(value.as_f64()?, 10.0, 2000.0, "height")?,
            "length (m)" => self.length_m = ranged(value.as_f64()?, 0.1, 500.0, "length")?,
            "velocity (m/s)" => self.velocity_ms = ranged(value.as_f64()?, 0.1, 30.0, "velocity")?,
            "friction factor" => self.friction = ranged(value.as_f64()?, 0.001, 0.1, "friction")?,
            "flow (CFM)" => self.cfm = ranged(value.as_f64()?, 1.0, 100_000.0, "flow")?,
            "max vel (FPM)" => {
                self.max_velocity_fpm = ranged(value.as_f64()?, 50.0, 5000.0, "max vel")?
            }
            other => return Err(format!("unknown HVAC control: {other:?}")),
        }
        Ok(())
    }
}

/// Compute the duct hydraulics for the current settings.
fn run_hvac(s: &HvacWorkbenchState) -> HvacResults {
    let cs = match s.shape {
        DuctShape::Round => CrossSection::Round {
            d: s.diameter_mm.max(10.0),
        },
        DuctShape::Rect => CrossSection::Rect {
            w: s.width_mm.max(10.0),
            h: s.height_mm.max(10.0),
        },
    };
    let dh_mm = cs.hydraulic_diameter_mm();
    let area_m2 = cs.area_mm2() * 1.0e-6;
    let pressure_drop_pa = pressure_drop::darcy_weisbach(
        dh_mm / 1000.0,
        s.length_m.max(0.1),
        s.velocity_ms.max(0.1),
        s.friction.max(0.001),
    );
    let (duct_w_in, duct_h_in) =
        flow::cfm_to_duct_size(s.cfm.max(1.0), s.max_velocity_fpm.max(50.0));
    HvacResults {
        hydraulic_diameter_mm: dh_mm,
        area_m2,
        pressure_drop_pa,
        duct_w_in,
        duct_h_in,
    }
}

/// Draw the HVAC workbench (a no-op unless toggled on via View → HVAC).
pub fn draw_hvac_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_hvac_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_hvac_workbench",
        "HVAC",
        |app, ui| {
            ui.label(
                egui::RichText::new("duct sizing + pressure drop · valenx-hvac")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.hvac;
            ui.horizontal(|ui| {
                ui.selectable_value(&mut s.shape, DuctShape::Round, "Round");
                ui.selectable_value(&mut s.shape, DuctShape::Rect, "Rectangular");
            });
            // Associate each numeric `DragValue` with its caption via
            // `labelled_by`, so the spin button carries the caption as its
            // accessibility / UI-Automation Name (egui clears a DragValue's own
            // Name, leaving it anonymous to a screen reader / AI driver
            // otherwise).
            egui::Grid::new("hvac_params")
                .num_columns(2)
                .show(ui, |ui| {
                    match s.shape {
                        DuctShape::Round => {
                            let l = ui.label("diameter (mm)");
                            ui.add(
                                egui::DragValue::new(&mut s.diameter_mm)
                                    .speed(2.0)
                                    .range(10.0..=2000.0),
                            )
                            .labelled_by(l.id);
                            ui.end_row();
                        }
                        DuctShape::Rect => {
                            let l = ui.label("width (mm)");
                            ui.add(
                                egui::DragValue::new(&mut s.width_mm)
                                    .speed(2.0)
                                    .range(10.0..=2000.0),
                            )
                            .labelled_by(l.id);
                            ui.end_row();
                            let l = ui.label("height (mm)");
                            ui.add(
                                egui::DragValue::new(&mut s.height_mm)
                                    .speed(2.0)
                                    .range(10.0..=2000.0),
                            )
                            .labelled_by(l.id);
                            ui.end_row();
                        }
                    }
                    let l = ui.label("length (m)");
                    ui.add(
                        egui::DragValue::new(&mut s.length_m)
                            .speed(0.5)
                            .range(0.1..=500.0),
                    )
                    .labelled_by(l.id);
                    ui.end_row();
                    let l = ui.label("velocity (m/s)");
                    ui.add(
                        egui::DragValue::new(&mut s.velocity_ms)
                            .speed(0.1)
                            .range(0.1..=30.0),
                    )
                    .labelled_by(l.id);
                    ui.end_row();
                    let l = ui.label("friction factor");
                    ui.add(
                        egui::DragValue::new(&mut s.friction)
                            .speed(0.001)
                            .range(0.001..=0.1),
                    )
                    .labelled_by(l.id);
                    ui.end_row();
                    let l = ui.label("flow (CFM)");
                    ui.add(
                        egui::DragValue::new(&mut s.cfm)
                            .speed(10.0)
                            .range(1.0..=100000.0),
                    )
                    .labelled_by(l.id);
                    ui.end_row();
                    let l = ui.label("max vel (FPM)");
                    ui.add(
                        egui::DragValue::new(&mut s.max_velocity_fpm)
                            .speed(10.0)
                            .range(50.0..=5000.0),
                    )
                    .labelled_by(l.id);
                    ui.end_row();
                });
            ui.separator();
            if ui.button("▶ Compute").clicked() {
                let r = run_hvac(s);
                s.results = Some(r);
            }
            if let Some(r) = &s.results {
                ui.separator();
                ui.label(
                    egui::RichText::new(format!(
                        "hydraulic Ø {:.1} mm\narea {:.4} m²\npressure drop {:.2} Pa\nCFM duct size {:.1} × {:.1} in",
                        r.hydraulic_diameter_mm, r.area_m2, r.pressure_drop_pa, r.duct_w_in, r.duct_h_in,
                    ))
                    .monospace()
                    .small(),
                );
            }
        },
    );
    if close {
        app.show_hvac_workbench = false;
    }
}

/// Build the **HVAC** result card for the Workbench+Agent bridge — a DATA-ONLY
/// [`crate::WorkspaceProduct`] (`mesh: None`) whose `lines` are the genuine duct
/// hydraulics ([`run_hvac`]) for the canonical default duct (a 200 mm round duct,
/// 10 m long at 5 m/s): hydraulic diameter, area, Darcy–Weisbach pressure drop,
/// and the CFM-recommended duct size. Registered as the `"hvac"` producer in
/// [`crate::products_registry::lookup`]; the tile renders it as a text card, not
/// a 3-D view. The rows mirror the panel's own readout format.
pub(crate) fn hvac_product() -> crate::WorkspaceProduct {
    let s = HvacWorkbenchState::default();
    let r = run_hvac(&s);
    let readout = format!(
        "hydraulic Ø {:.1} mm\narea {:.4} m²\npressure drop {:.2} Pa\nCFM duct size {:.1} × {:.1} in",
        r.hydraulic_diameter_mm, r.area_m2, r.pressure_drop_pa, r.duct_w_in, r.duct_h_in,
    );
    crate::WorkspaceProduct {
        title: "HVAC".into(),
        lines: crate::products_registry::lines_from_readout(&readout),
        mesh: None,
        vertex_colors: None,
        camera: Default::default(),
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
    use crate::agent_commands::AgentValue;

    #[test]
    fn agent_set_sets_param_unknown_and_type_mismatch_err() {
        let mut s = HvacWorkbenchState::default();
        s.agent_set("velocity (m/s)", &AgentValue::Float(5.0))
            .unwrap();
        assert_eq!(s.velocity_ms, 5.0);
        // Unknown caption -> Err (no panic).
        assert!(s.agent_set("no such control", &AgentValue::Int(1)).is_err());
        // Type mismatch (bool into a numeric field) -> Err.
        assert!(s
            .agent_set("velocity (m/s)", &AgentValue::Bool(true))
            .is_err());
        // Out-of-range (velocity > 30) -> Err, field untouched.
        assert!(s
            .agent_set("velocity (m/s)", &AgentValue::Float(99.0))
            .is_err());
        assert_eq!(s.velocity_ms, 5.0, "rejected set leaves field untouched");
    }

    #[test]
    fn pressure_drop_is_positive_and_grows_with_length() {
        let r1 = run_hvac(&HvacWorkbenchState::default());
        assert!(
            r1.pressure_drop_pa > 0.0,
            "ΔP positive: {}",
            r1.pressure_drop_pa
        );
        assert!(r1.hydraulic_diameter_mm > 0.0 && r1.area_m2 > 0.0);
        let r2 = run_hvac(&HvacWorkbenchState {
            length_m: 20.0,
            ..Default::default()
        });
        assert!(
            r2.pressure_drop_pa > r1.pressure_drop_pa,
            "ΔP grows with duct length"
        );
    }

    #[test]
    fn rectangular_hydraulic_diameter_matches_4a_over_p() {
        let r = run_hvac(&HvacWorkbenchState {
            shape: DuctShape::Rect,
            width_mm: 200.0,
            height_mm: 100.0,
            ..Default::default()
        });
        // 4·A/P = 4·200·100 / (2·300) = 133.33 mm.
        assert!(
            (r.hydraulic_diameter_mm - 133.333).abs() < 0.1,
            "got {}",
            r.hydraulic_diameter_mm
        );
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    /// Render the whole workbench panel once in a headless egui context.
    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_hvac_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_hvac_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_hvac_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_hvac_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_without_panic() {
        // Populate the results card so the readout branch is exercised too.
        let mut app = ValenxApp::default();
        app.show_hvac_workbench = true;
        app.hvac.results = Some(run_hvac(&app.hvac));
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // The duct-sizing DragValues are SpinButtons; each must be `labelled_by`
        // its caption (egui clears a DragValue's own Name), so an AI / screen
        // reader can find the control by the caption text. The default Round
        // shape shows the diameter row (not width/height).
        let mut app = ValenxApp::default();
        app.show_hvac_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // diameter, length, velocity, friction, cfm, max vel.
        assert!(
            spin_buttons.len() >= 6,
            "expected the HVAC numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every HVAC DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["diameter (mm)", "length (m)", "max vel (FPM)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
        // The Compute button stays a named, invokable node.
        assert!(
            nodes.iter().any(|(_, n)| n.role() == Role::Button
                && n.name().is_some_and(|s| s.contains("Compute"))),
            "the Compute button is a named, invokable node"
        );
    }
}
