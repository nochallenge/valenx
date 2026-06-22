//! The right-side **Power Factor Workbench** panel — native single-phase
//! AC power-triangle and shunt-capacitor correction over
//! `valenx-powerfactor`.
//!
//! Mirrors the Heat Transfer / Buckling workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_powerfactor_workbench`,
//! toggled from the View menu. The form sets a lagging load (real power,
//! present and target power factors) and a supply (voltage, frequency);
//! "Analyze" resolves the present power triangle (`S`, `P`, `Q`,
//! `PF = cos(phi)`) and sizes the shunt capacitor that lifts the power
//! factor to the target (`Qc = P(tan(phi1) - tan(phi2))`, plus the
//! physical capacitance and current), and "Show 3-D" loads a
//! power-triangle prism (legs scaled to `P` and `Q`) into the central
//! viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_powerfactor::{Correction, Phase, PowerTriangle};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which leg of the present triangle the user supplies the correction from.
///
/// Both modes correct a lagging load; they differ only in how the present
/// operating point is described to `valenx-powerfactor`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum LoadInput {
    /// Describe the load by real power `P` (watts) and a present power
    /// factor.
    RealPowerPf,
    /// Describe the load by RMS voltage, RMS current and a present power
    /// factor (the nameplate case); `P` is derived as `V*I*PF`.
    VoltAmpsPf,
}

/// Persistent form + result state for the Power Factor Workbench.
pub struct PowerFactorWorkbenchState {
    /// How the present load is specified.
    load_input: LoadInput,
    /// Real power `P` (W) — used in [`LoadInput::RealPowerPf`].
    real_w: f64,
    /// RMS supply voltage `V` (V) — also sizes the physical capacitor.
    voltage_v: f64,
    /// RMS load current `I` (A) — used in [`LoadInput::VoltAmpsPf`].
    current_a: f64,
    /// Supply frequency `f` (Hz) — sizes the physical capacitor.
    frequency_hz: f64,
    /// Present (lagging) power factor `cos(phi1)`, in `(0, 1]`.
    pf_present: f64,
    /// Target power factor `cos(phi2)`, in `(0, 1]`, strictly above the
    /// present value.
    pf_target: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D triangle solid (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for PowerFactorWorkbenchState {
    fn default() -> Self {
        // A 10 kW lagging load at 230 V / 50 Hz, raised from 0.7 to 0.95:
        // present Q1 = P*tan(acos 0.7) = 10202 var, capacitor
        // Qc = P*(tan(acos 0.7) - tan(acos 0.95)) = 6915.20 var, and the
        // apparent power drops from 14285.71 VA to 10526.32 VA. The
        // physical shunt is C = Qc/(2*pi*f*V^2) ~ 416 uF carrying ~30 A.
        Self {
            load_input: LoadInput::RealPowerPf,
            real_w: 10000.0,
            voltage_v: 230.0,
            current_a: 62.1,
            frequency_hz: 50.0,
            pf_present: 0.7,
            pf_target: 0.95,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Power Factor Workbench right-side panel. A no-op when the
/// `show_powerfactor_workbench` toggle is off.
pub fn draw_powerfactor_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_powerfactor_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_powerfactor_workbench",
        "Power Factor",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native AC power triangle + shunt-capacitor correction · valenx-powerfactor",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.powerfactor;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Load specified by").strong());
                    ui.radio_value(
                        &mut s.load_input,
                        LoadInput::RealPowerPf,
                        "real power P + power factor",
                    );
                    ui.radio_value(
                        &mut s.load_input,
                        LoadInput::VoltAmpsPf,
                        "voltage × current + power factor",
                    );

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Load").strong());
                    match s.load_input {
                        LoadInput::RealPowerPf => {
                            ui.horizontal(|ui| {
                                ui.label("real power P (W)");
                                ui.add(egui::DragValue::new(&mut s.real_w).speed(100.0));
                            });
                        }
                        LoadInput::VoltAmpsPf => {
                            ui.horizontal(|ui| {
                                ui.label("voltage V (V)");
                                ui.add(egui::DragValue::new(&mut s.voltage_v).speed(1.0));
                            });
                            ui.horizontal(|ui| {
                                ui.label("current I (A)");
                                ui.add(egui::DragValue::new(&mut s.current_a).speed(0.5));
                            });
                        }
                    }
                    ui.horizontal(|ui| {
                        ui.label("present PF cos φ₁");
                        ui.add(egui::DragValue::new(&mut s.pf_present).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Correction target").strong());
                    ui.horizontal(|ui| {
                        ui.label("target PF cos φ₂");
                        ui.add(egui::DragValue::new(&mut s.pf_target).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Supply (for capacitor sizing)").strong());
                    ui.horizontal(|ui| {
                        ui.label("voltage V (V)");
                        ui.add(egui::DragValue::new(&mut s.voltage_v).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("frequency f (Hz)");
                        ui.add(egui::DragValue::new(&mut s.frequency_hz).speed(1.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_powerfactor(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build the present power triangle as a 3-D prism (its legs scaled to the real power P and reactive power Q) and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Power factor").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_powerfactor_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.powerfactor` borrow is
    // released here): build the triangle's 3-D prism and load it.
    if app.powerfactor.show_3d_request {
        app.powerfactor.show_3d_request = false;
        load_triangle_3d(app);
    }
}

/// Validate the form, evaluate the load + correction and format the
/// readout.
fn run_powerfactor(s: &mut PowerFactorWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Resolve the present (lagging) power triangle from the active input
/// mode. The quantity both the readout and the 3-D gate need. Extracted
/// so it is unit-testable and shared.
fn present_triangle(s: &PowerFactorWorkbenchState) -> Result<PowerTriangle, String> {
    match s.load_input {
        LoadInput::RealPowerPf => {
            // Reconstruct the triangle from P and the lagging reactive
            // power Q1 = P*tan(phi1); clamp the power factor to a valid
            // cos before taking acos so a stray entry maps to a domain
            // error rather than a NaN angle.
            let pf = s.pf_present.clamp(-1.0, 1.0);
            let phi1 = pf.acos();
            let q1 = s.real_w * phi1.tan();
            PowerTriangle::from_p_q(s.real_w, q1).map_err(|e| e.to_string())
        }
        LoadInput::VoltAmpsPf => {
            PowerTriangle::from_vi_pf(s.voltage_v, s.current_a, s.pf_present, Phase::Lagging)
                .map_err(|e| e.to_string())
        }
    }
}

/// Evaluate the present triangle and the shunt-capacitor correction and
/// format the full readout, mapping any domain error to a display string.
/// Extracted so it is unit-testable.
fn compute(s: &PowerFactorWorkbenchState) -> Result<String, String> {
    let load = present_triangle(s)?;
    let corr = Correction::for_triangle(&load, s.pf_target).map_err(|e| e.to_string())?;
    let cap_f = corr
        .capacitance_farads(s.voltage_v, s.frequency_hz)
        .map_err(|e| e.to_string())?;
    let cap_a = corr
        .capacitor_current_a(s.voltage_v)
        .map_err(|e| e.to_string())?;
    let phi1_deg = load.phase_angle_rad().to_degrees();

    Ok(format!(
        "present load\n\
         apparent S      : {:.2} VA\n\
         real P          : {:.2} W\n\
         reactive Q      : {:.2} var\n\
         power factor    : {:.4}\n\
         phase angle φ₁  : {:.2} deg\n\n\
         correction → PF {:.4}\n\
         capacitor Qc    : {:.2} var\n\
         reactive after  : {:.2} var\n\
         apparent after  : {:.2} VA\n\
         capacitance C   : {:.3e} F\n\
         cap. current Ic : {:.2} A",
        load.apparent_va,
        load.real_w,
        load.reactive_var,
        load.power_factor,
        phi1_deg,
        corr.power_factor_after,
        corr.capacitor_var,
        corr.reactive_after_var,
        corr.apparent_after_va,
        cap_f,
        cap_a,
    ))
}

/// Build the present power triangle as a triangle [`Mesh`] — a right
/// triangular prism whose base leg is the real power `P` and whose
/// upright leg is the reactive power `Q`, both normalised to a unit box
/// so the shape is legible at any scale (the numbers are the
/// `valenx-powerfactor` result, not the geometry). `None` for an invalid
/// configuration.
fn triangle_prism_mesh(s: &PowerFactorWorkbenchState) -> Option<Mesh> {
    let load = present_triangle(s).ok()?;
    let p = load.real_w;
    let q = load.reactive_magnitude_var();
    // Normalise the legs to a unit-ish box; guard the degenerate all-zero
    // case (present_triangle already rejects a collapsed triangle).
    let scale = p.hypot(q);
    if scale <= 0.0 {
        return None;
    }
    let leg_p = p / scale; // along +x
    let leg_q = q / scale; // along +z
    let depth = 0.25_f64; // prism thickness along +y

    // Front face (y = 0) and back face (y = depth): right triangle with
    // the right angle at the origin.
    let nodes: Vec<Vector3<f64>> = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(leg_p, 0.0, 0.0),
        Vector3::new(0.0, 0.0, leg_q),
        Vector3::new(0.0, depth, 0.0),
        Vector3::new(leg_p, depth, 0.0),
        Vector3::new(0.0, depth, leg_q),
    ];
    // Two triangular caps (front y=0, back y=depth) plus the three
    // rectangular sides, each split into two triangles.
    let tris: Vec<usize> = vec![
        0, 2, 1, // front cap
        3, 4, 5, // back cap
        0, 1, 4, 0, 4, 3, // base (P leg)
        0, 3, 5, 0, 5, 2, // upright (Q leg)
        1, 2, 5, 1, 5, 4, // hypotenuse (S)
    ];

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-powerfactor");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D triangle prism and load it into the central viewport.
fn load_triangle_3d(app: &mut ValenxApp) {
    let Some(mesh) = triangle_prism_mesh(&app.powerfactor) else {
        app.powerfactor.error =
            Some("load parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<triangle>/valenx-powerfactor"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"powerfactor"}`** product: a DATA-ONLY
/// text card of the workbench's own `compute()` readout rows. The power
/// triangle is a 2-D vector relationship, not a physical object — the panel's
/// extruded triangular prism is a decorative chart-as-geometry — so the bridge
/// product is right-sized to a card (`mesh: None`) carrying just the readout
/// (the confidence badge is appended centrally). The panel's "Show 3-D" button
/// still builds that representative prism into the central viewport. Registered
/// in [`crate::products_registry`]; the per-tool builder the registry
/// dispatches to. Pure — driven off [`PowerFactorWorkbenchState::default`].
pub(crate) fn powerfactor_product() -> crate::WorkspaceProduct {
    let s = PowerFactorWorkbenchState::default();
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical power factor ⇒ readout computes"),
    );
    crate::WorkspaceProduct {
        title: "Power factor (correction)".into(),
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
        let s = PowerFactorWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_triangle_and_correction() {
        let mut s = PowerFactorWorkbenchState::default();
        run_powerfactor(&mut s);
        assert!(
            s.error.is_none(),
            "default load should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("apparent S"));
        assert!(s.result.contains("capacitor Qc"));
        assert!(s.result.contains("capacitance C"));
        // P=10 kW, 0.7 -> 0.95: Qc = P*(tan(acos 0.7)-tan(acos 0.95))
        // = 6915.20 var, and S1 = P/0.7 = 14285.71 VA.
        assert!(s.result.contains("6915.20"));
        assert!(s.result.contains("14285.71"));
    }

    #[test]
    fn analyze_rejects_non_improving_target() {
        // A target below the present PF (0.7) leaves nothing to correct.
        let mut s = PowerFactorWorkbenchState {
            pf_target: 0.6,
            ..Default::default()
        };
        run_powerfactor(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn power_factor_equals_cos_phi() {
        // Ground truth: PF = cos(phi). Build the present triangle and
        // confirm its power factor matches the cosine of its own phase
        // angle to full precision.
        let s = PowerFactorWorkbenchState::default();
        let load = present_triangle(&s).unwrap();
        let phi = load.phase_angle_rad();
        assert!((load.power_factor - phi.cos()).abs() < 1e-12);
        // And it reproduces the requested 0.7 within reconstruction noise.
        assert!((load.power_factor - 0.7).abs() < 1e-9);
    }

    #[test]
    fn capacitor_var_matches_hand_computed_closed_form() {
        // GROUND TRUTH: Qc = P*(tan(phi1) - tan(phi2)). Hand-compute the
        // tangents from the power factors and compare to the crate result.
        let s = PowerFactorWorkbenchState::default();
        let load = present_triangle(&s).unwrap();
        let corr = Correction::for_triangle(&load, s.pf_target).unwrap();
        let tan1: f64 = (1.0_f64 - 0.7 * 0.7).sqrt() / 0.7;
        let tan2: f64 = (1.0_f64 - 0.95 * 0.95).sqrt() / 0.95;
        let expected_qc = 10000.0 * (tan1 - tan2);
        assert!(
            (corr.capacitor_var - expected_qc).abs() < 1e-6,
            "Qc {} vs hand {expected_qc}",
            corr.capacitor_var
        );
        // Textbook magnitude: ~6915.20 var.
        assert!((corr.capacitor_var - 6915.20).abs() < 1.0);
    }

    #[test]
    fn triangle_prism_for_default_is_nonempty_and_in_range() {
        let s = PowerFactorWorkbenchState::default();
        let mesh = triangle_prism_mesh(&s).expect("default load yields a prism");
        assert!(mesh.nodes.len() >= 6, "expected a triangular prism");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn triangle_prism_none_for_invalid() {
        // A present PF above unity is rejected by from_p_q's domain (the
        // reconstructed Q becomes non-finite / the triangle invalid).
        let s = PowerFactorWorkbenchState {
            load_input: LoadInput::VoltAmpsPf,
            voltage_v: 0.0,
            ..Default::default()
        };
        assert!(triangle_prism_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_powerfactor_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_powerfactor_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_powerfactor_workbench = true;
        run_powerfactor(&mut app.powerfactor);
        draw_workbench(&mut app);
    }
}
