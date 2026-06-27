//! The right-side **Thermoregulation Workbench** panel — a native
//! single-node human heat balance over `valenx-thermoreg`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_thermoreg_workbench`,
//! toggled from the View menu. The form sets a body (mass, height, core and
//! skin temperature), a surrounding environment (air / radiant temperature,
//! convective film, emissivity) and a metabolic + sweat state; "Analyze"
//! closes the conceptual physiology heat balance `M - W = R + C + E + S` —
//! reporting the DuBois surface area, the radiative / convective /
//! evaporative losses, the storage residual and the resulting core-temperature
//! drift — and "Show 3-D" loads a representative human-body capsule into the
//! central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_thermoreg::{
    core_temp_change, heat_balance, Body, Environment, Metabolism, Sweat, SKIN_EMISSIVITY,
};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// A coarse activity preset that sets the metabolic heat-production rate
/// `M`. Convenience over typing a watt figure; the user can still override
/// `M` (and the external work `W`) directly afterwards.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Activity {
    /// At rest / basal-plus-resting, `M ≈ 100 W`.
    Resting,
    /// Light activity (slow walking, desk work standing), `M ≈ 160 W`.
    Light,
    /// Moderate activity (brisk walk, light cycling), `M ≈ 250 W`.
    Moderate,
    /// Vigorous activity (running, hard cycling), `M ≈ 450 W`.
    Vigorous,
}

impl Activity {
    /// The metabolic heat-production rate `M` (W) this preset selects.
    fn metabolic_w(self) -> f64 {
        match self {
            Activity::Resting => 100.0,
            Activity::Light => 160.0,
            Activity::Moderate => 250.0,
            Activity::Vigorous => 450.0,
        }
    }
}

/// Persistent form + result state for the Thermoregulation Workbench.
pub struct ThermoRegWorkbenchState {
    /// Body mass (kg).
    mass_kg: f64,
    /// Standing height (cm) — drives the DuBois surface-area estimate.
    height_cm: f64,
    /// Core (deep-body) temperature (°C).
    core_temp_c: f64,
    /// Mean skin temperature (°C).
    skin_temp_c: f64,
    /// Ambient (dry-bulb) air temperature (°C).
    air_temp_c: f64,
    /// Mean radiant temperature of the surroundings (°C).
    radiant_temp_c: f64,
    /// Convective heat-transfer coefficient `h` (W/m²·K).
    convective_coeff: f64,
    /// Long-wave skin / clothing emissivity (0–1).
    emissivity: f64,
    /// Activity preset that sets the metabolic rate `M`.
    activity: Activity,
    /// Metabolic heat-production rate `M` (W).
    metabolic_w: f64,
    /// External mechanical work rate `W` (W).
    work_w: f64,
    /// Evaporated-sweat mass rate (g/min) — converted to kg/s for the model.
    sweat_g_per_min: f64,
    /// Integration window for the core-temperature drift readout (minutes).
    window_min: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D body solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for ThermoRegWorkbenchState {
    fn default() -> Self {
        // A 70 kg, 175 cm adult (DuBois ~1.85 m²) at moderate activity
        // (M = 250 W) in a still 24 °C room, sweating lightly (2.4 g/min
        // ~ 4e-5 kg/s): losses ~272 W, so the body cools gently
        // (S ~ -22 W, ~ -0.32 °C over an hour) — near, but not at, balance.
        Self {
            mass_kg: 70.0,
            height_cm: 175.0,
            core_temp_c: 37.0,
            skin_temp_c: 34.0,
            air_temp_c: 24.0,
            radiant_temp_c: 24.0,
            convective_coeff: 3.5,
            emissivity: SKIN_EMISSIVITY,
            activity: Activity::Moderate,
            metabolic_w: 250.0,
            work_w: 0.0,
            sweat_g_per_min: 2.4,
            window_min: 60.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Thermoregulation Workbench right-side panel. A no-op when the
/// `show_thermoreg_workbench` toggle is off.
pub fn draw_thermoreg_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_thermoreg_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_thermoreg_workbench",
        "Thermoregulation",
        |app, ui| {
            ui.label(
                egui::RichText::new("single-node human heat balance · valenx-thermoreg")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.thermoreg;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Body").strong());
                    // Associate each numeric `DragValue` with its caption via
                    // `labelled_by`, so the spin button carries the caption as
                    // its accessibility / UI-Automation Name (egui clears a
                    // DragValue's own Name, leaving it anonymous to a screen
                    // reader / AI driver otherwise); the hover text mirrors the
                    // caption for a mouse user.
                    ui.horizontal(|ui| {
                        let lbl = ui.label("mass (kg)");
                        ui.add(egui::DragValue::new(&mut s.mass_kg).speed(0.5))
                            .labelled_by(lbl.id)
                            .on_hover_text("mass (kg)");
                    });
                    ui.horizontal(|ui| {
                        let lbl = ui.label("height (cm)");
                        ui.add(egui::DragValue::new(&mut s.height_cm).speed(0.5))
                            .labelled_by(lbl.id)
                            .on_hover_text("height (cm)");
                    });
                    ui.horizontal(|ui| {
                        let lbl = ui.label("core T (°C)");
                        ui.add(egui::DragValue::new(&mut s.core_temp_c).speed(0.1))
                            .labelled_by(lbl.id)
                            .on_hover_text("core T (°C)");
                    });
                    ui.horizontal(|ui| {
                        let lbl = ui.label("skin T (°C)");
                        ui.add(egui::DragValue::new(&mut s.skin_temp_c).speed(0.1))
                            .labelled_by(lbl.id)
                            .on_hover_text("skin T (°C)");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Environment").strong());
                    ui.horizontal(|ui| {
                        let lbl = ui.label("air T (°C)");
                        ui.add(egui::DragValue::new(&mut s.air_temp_c).speed(0.5))
                            .labelled_by(lbl.id)
                            .on_hover_text("air T (°C)");
                    });
                    ui.horizontal(|ui| {
                        let lbl = ui.label("radiant T (°C)");
                        ui.add(egui::DragValue::new(&mut s.radiant_temp_c).speed(0.5))
                            .labelled_by(lbl.id)
                            .on_hover_text("radiant T (°C)");
                    });
                    ui.horizontal(|ui| {
                        let lbl = ui.label("convective h (W/m²K)");
                        ui.add(egui::DragValue::new(&mut s.convective_coeff).speed(0.25))
                            .labelled_by(lbl.id)
                            .on_hover_text("convective h (W/m²K)");
                    });
                    ui.horizontal(|ui| {
                        let lbl = ui.label("emissivity");
                        ui.add(egui::DragValue::new(&mut s.emissivity).speed(0.01))
                            .labelled_by(lbl.id)
                            .on_hover_text("emissivity");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Metabolism & sweat").strong());
                    ui.horizontal(|ui| {
                        ui.label("activity");
                        if ui
                            .radio_value(&mut s.activity, Activity::Resting, "rest")
                            .clicked()
                        {
                            s.metabolic_w = s.activity.metabolic_w();
                        }
                        if ui
                            .radio_value(&mut s.activity, Activity::Light, "light")
                            .clicked()
                        {
                            s.metabolic_w = s.activity.metabolic_w();
                        }
                        if ui
                            .radio_value(&mut s.activity, Activity::Moderate, "moderate")
                            .clicked()
                        {
                            s.metabolic_w = s.activity.metabolic_w();
                        }
                        if ui
                            .radio_value(&mut s.activity, Activity::Vigorous, "vigorous")
                            .clicked()
                        {
                            s.metabolic_w = s.activity.metabolic_w();
                        }
                    });
                    ui.horizontal(|ui| {
                        let lbl = ui.label("metabolic M (W)");
                        ui.add(egui::DragValue::new(&mut s.metabolic_w).speed(5.0))
                            .labelled_by(lbl.id)
                            .on_hover_text("metabolic M (W)");
                    });
                    ui.horizontal(|ui| {
                        let lbl = ui.label("work W (W)");
                        ui.add(egui::DragValue::new(&mut s.work_w).speed(5.0))
                            .labelled_by(lbl.id)
                            .on_hover_text("work W (W)");
                    });
                    ui.horizontal(|ui| {
                        let lbl = ui.label("sweat (g/min)");
                        ui.add(egui::DragValue::new(&mut s.sweat_g_per_min).speed(0.1))
                            .labelled_by(lbl.id)
                            .on_hover_text("sweat (g/min)");
                    });
                    ui.horizontal(|ui| {
                        let lbl = ui.label("window (min)");
                        ui.add(egui::DragValue::new(&mut s.window_min).speed(1.0))
                            .labelled_by(lbl.id)
                            .on_hover_text("window (min)");
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_thermoreg(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative human-body capsule (torso + head, sized by height) as a 3-D solid and load it into the central viewport to orbit",
                        )
                        .clicked()
                    {
                        s.show_3d_request = true;
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Heat balance").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_thermoreg_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.thermoreg` borrow is
    // released here): build the body's 3-D solid and load it.
    if app.thermoreg.show_3d_request {
        app.thermoreg.show_3d_request = false;
        load_body_3d(app);
    }
}

/// Validate the form, resolve the heat balance and format the readout.
fn run_thermoreg(s: &mut ThermoRegWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated `valenx-thermoreg` objects from the form — the body,
/// environment, metabolism and sweat state shared by both the readout and
/// the 3-D gate. Extracted so it is unit-testable.
fn models(s: &ThermoRegWorkbenchState) -> Result<(Body, Environment, Metabolism, Sweat), String> {
    let body = Body::standard_adult(s.mass_kg, s.height_cm, s.core_temp_c, s.skin_temp_c)
        .map_err(|e| e.to_string())?;
    let env = Environment::new(
        s.air_temp_c,
        s.radiant_temp_c,
        s.convective_coeff,
        s.emissivity,
    )
    .map_err(|e| e.to_string())?;
    let met = Metabolism::new(s.metabolic_w, s.work_w).map_err(|e| e.to_string())?;
    // g/min -> kg/s : * 1e-3 / 60.
    let sweat = Sweat::from_rate(s.sweat_g_per_min * 1.0e-3 / 60.0).map_err(|e| e.to_string())?;
    Ok((body, env, met, sweat))
}

/// Resolve the heat balance and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &ThermoRegWorkbenchState) -> Result<String, String> {
    let (body, env, met, sweat) = models(s)?;
    let bal = heat_balance(&body, &env, &met, &sweat);
    let window_s = s.window_min * 60.0;
    let drift_c = core_temp_change(&body, bal.net_heat_w(), window_s).map_err(|e| e.to_string())?;
    let projected_core_c = body.core_temp_c + drift_c;
    let tendency = if bal.net_heat_w() > 0.0 {
        "warming"
    } else if bal.net_heat_w() < 0.0 {
        "cooling"
    } else {
        "balanced"
    };

    Ok(format!(
        "DuBois area     : {:.3} m²\n\
         metabolic M-W   : {:.1} W\n\n\
         radiation R     : {:.1} W\n\
         convection C    : {:.1} W\n\
         evaporation E   : {:.1} W\n\
         total loss      : {:.1} W\n\n\
         storage S       : {:.1} W  ({tendency})\n\
         drift over {:.0} min: {drift_c:+.3} °C\n\
         projected core  : {projected_core_c:.2} °C",
        body.surface_area_m2,
        bal.metabolic_net_w,
        bal.radiation_w,
        bal.convection_w,
        bal.evaporation_w,
        bal.total_loss_w(),
        bal.net_heat_w(),
        s.window_min,
    ))
}

/// Append a capsule (a cylinder of radius `r` capped by hemispheres) aligned
/// with the world `z` axis, spanning the cylinder centres `c0` (bottom) to
/// `c1` (top), to the triangle buffers. `seg` is the radial segment count
/// and `rings` the number of latitude rings per hemispherical cap.
#[allow(clippy::too_many_arguments)]
fn push_capsule(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c0: Vector3<f64>,
    c1: Vector3<f64>,
    r: f64,
    seg: usize,
    rings: usize,
) {
    use std::f64::consts::PI;

    // Build a stack of rings from the bottom cap, up the cylinder, to the
    // top cap. Each ring is `seg` vertices around the z axis.
    let mut ring_centres: Vec<(Vector3<f64>, f64)> = Vec::new();
    // Bottom hemisphere: latitude from -pi/2 (pole) up to 0 (equator).
    for i in 0..=rings {
        let phi = -PI / 2.0 + (PI / 2.0) * (i as f64 / rings as f64);
        ring_centres.push((c0 + Vector3::new(0.0, 0.0, r * phi.sin()), r * phi.cos()));
    }
    // Cylinder body: bottom and top rims (the equators sit at c0 / c1).
    ring_centres.push((c0, r));
    ring_centres.push((c1, r));
    // Top hemisphere: latitude from 0 (equator) up to +pi/2 (pole).
    for i in 0..=rings {
        let phi = (PI / 2.0) * (i as f64 / rings as f64);
        ring_centres.push((c1 + Vector3::new(0.0, 0.0, r * phi.sin()), r * phi.cos()));
    }

    let base = nodes.len();
    for (centre, radius) in &ring_centres {
        for k in 0..seg {
            let theta = 2.0 * PI * (k as f64 / seg as f64);
            nodes.push(centre + Vector3::new(radius * theta.cos(), radius * theta.sin(), 0.0));
        }
    }

    // Stitch adjacent rings into a quad-strip of two triangles each.
    for ring in 0..ring_centres.len() - 1 {
        let a = base + ring * seg;
        let b = base + (ring + 1) * seg;
        for k in 0..seg {
            let kn = (k + 1) % seg;
            tris.extend_from_slice(&[a + k, a + kn, b + kn, a + k, b + kn, b + k]);
        }
    }
}

/// Build a representative human-body solid as a triangle [`Mesh`] — a torso
/// capsule with a smaller head capsule above it, sized from the standing
/// height. Representative geometry (not anatomically exact; the heat-balance
/// numbers are the `valenx-thermoreg` result). `None` for an invalid
/// configuration.
fn body_solid_mesh(s: &ThermoRegWorkbenchState) -> Option<Mesh> {
    // Gate on the same validation the readout uses.
    models(s).ok()?;

    let h = (s.height_cm / 100.0).max(0.5); // metres, guarded.
    let torso_r = 0.16 * h;
    let head_r = 0.07 * h;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Torso capsule: cylinder spanning roughly hip to shoulder.
    push_capsule(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.30 * h),
        Vector3::new(0.0, 0.0, 0.78 * h),
        torso_r,
        24,
        6,
    );
    // Head capsule above the torso.
    push_capsule(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.90 * h),
        Vector3::new(0.0, 0.0, 0.96 * h),
        head_r,
        20,
        5,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-thermoreg");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D body solid and load it into the central viewport.
fn load_body_3d(app: &mut ValenxApp) {
    let Some(mesh) = body_solid_mesh(&app.thermoreg) else {
        app.thermoreg.error =
            Some("body parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<body>/valenx-thermoreg"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"thermoreg"}`** product: a DATA-ONLY text
/// card of the workbench's own heat-balance headline numbers. The heat-balance
/// result has no characteristic shape — the panel's human-body capsule is a
/// representative decorative solid, not a real object — so the bridge product
/// is right-sized to a card (`mesh: None`) carrying just the readout (the
/// confidence badge is appended centrally). The panel's "Show 3-D" button
/// still builds that representative solid into the central viewport.
/// Registered in [`crate::products_registry`]; the per-tool builder the
/// registry dispatches to. Pure — driven off [`ThermoRegWorkbenchState::default`].
///
/// The readout rows mirror the panel's `compute()` readout.
pub(crate) fn thermoreg_product() -> crate::WorkspaceProduct {
    let s = ThermoRegWorkbenchState::default();
    let readout = compute(&s).expect("default body ⇒ a valid readout");
    let lines = crate::products_registry::lines_from_readout(&readout);
    crate::WorkspaceProduct {
        title: "Thermoregulation".into(),
        lines,
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

    #[test]
    fn default_state_is_idle() {
        let s = ThermoRegWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_balance_terms() {
        let mut s = ThermoRegWorkbenchState::default();
        run_thermoreg(&mut s);
        assert!(
            s.error.is_none(),
            "default body should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("DuBois area"));
        assert!(s.result.contains("radiation R"));
        assert!(s.result.contains("storage S"));
        assert!(s.result.contains("projected core"));
        // DuBois for 70 kg / 175 cm is ~1.848 m^2 (formatted to 3 dp).
        assert!(s.result.contains("1.848"), "readout: {}", s.result);
    }

    #[test]
    fn dubois_ground_truth_70kg_175cm() {
        // Hand calc: A = 0.007184 * 70^0.425 * 175^0.725.
        // 70^0.425 = 5.9078..., 175^0.725 = 43.555..., product * const
        // = 1.8481... m^2 — the canonical ~1.85 m^2 "standard man" area.
        let a = Body::dubois_surface_area(70.0, 175.0).unwrap();
        let expected = 0.007184 * 70.0_f64.powf(0.425) * 175.0_f64.powf(0.725);
        assert!((a - expected).abs() < 1e-12);
        assert!((a - 1.848_143).abs() < 1e-5, "got {a}");
    }

    #[test]
    fn analyze_rejects_zero_mass() {
        let mut s = ThermoRegWorkbenchState {
            mass_kg: 0.0,
            ..Default::default()
        };
        run_thermoreg(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_rejects_emissivity_above_one() {
        let mut s = ThermoRegWorkbenchState {
            emissivity: 1.5,
            ..Default::default()
        };
        run_thermoreg(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn hotter_environment_raises_storage() {
        // Same body/metabolism, a hotter & more radiant room cuts the
        // sensible losses, so more heat is stored (less cooling).
        let cool = {
            let s = ThermoRegWorkbenchState::default();
            let (b, e, m, w) = models(&s).unwrap();
            heat_balance(&b, &e, &m, &w).net_heat_w()
        };
        let warm = {
            let s = ThermoRegWorkbenchState {
                air_temp_c: 33.0,
                radiant_temp_c: 33.0,
                ..Default::default()
            };
            let (b, e, m, w) = models(&s).unwrap();
            heat_balance(&b, &e, &m, &w).net_heat_w()
        };
        assert!(warm > cool, "warm {warm} should exceed cool {cool}");
    }

    #[test]
    fn body_mesh_for_default_is_nonempty_and_in_range() {
        let s = ThermoRegWorkbenchState::default();
        let mesh = body_solid_mesh(&s).expect("default body yields a solid");
        assert!(mesh.nodes.len() > 8, "expected torso + head capsules");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn body_mesh_none_for_invalid() {
        let s = ThermoRegWorkbenchState {
            height_cm: 0.0,
            ..Default::default()
        };
        assert!(body_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_thermoreg_workbench(app, ctx);
        });
    }

    /// As `draw_workbench`, but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_thermoreg_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_thermoreg_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_thermoreg_workbench = true;
        run_thermoreg(&mut app.thermoreg);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every body / environment / metabolism DragValue is a SpinButton;
        // each must be `labelled_by` its caption (egui clears a DragValue's own
        // Name), so an AI / screen reader can find the control by the caption
        // text. All twelve are visible in the default state (no conditional
        // fields in this panel).
        let mut app = ValenxApp::default();
        app.show_thermoreg_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // mass, height, core T, skin T, air T, radiant T, convective h,
        // emissivity, metabolic M, work W, sweat, window.
        assert!(
            spin_buttons.len() >= 12,
            "expected the thermoreg numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every thermoreg DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["mass (kg)", "core T (°C)", "metabolic M (W)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
        // The Analyze button stays named/invokable.
        assert!(
            nodes.iter().any(|(_, n)| n.role() == Role::Button
                && n.name().is_some_and(|s| s.contains("Analyze"))),
            "the Analyze button is a named, invokable node"
        );
    }
}
