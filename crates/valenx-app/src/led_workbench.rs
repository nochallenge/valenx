//! The right-side **LED Workbench** panel — native series current-limiting
//! resistor sizing over `valenx-led`.
//!
//! Mirrors the MOSFET / Heat Transfer / Antenna workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_led_workbench`,
//! toggled from the View menu. The form sets a supply voltage, an LED
//! forward voltage and a target LED current, with a `LedMode` selecting a
//! single LED or `n` identical LEDs in series; "Analyze" solves the single
//! Kirchhoff voltage loop — the current-limiting resistor `R = (Vs - Vf) / I`
//! and the steady-state power split between the LED(s) and the resistor — and
//! "Show 3-D" loads a representative through-hole LED (cylindrical body plus
//! a domed lens) solid into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_led::{LedCircuit, LedString};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Whether the form drives a single LED or a series string of `n` identical
/// LEDs (their forward voltages add).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LedMode {
    /// One LED in series with the current-limiting resistor
    /// ([`LedCircuit`]).
    Single,
    /// `count` identical LEDs in series, sharing one resistor
    /// ([`LedString`]).
    Series,
}

/// Persistent form + result state for the LED Workbench.
pub struct LedWorkbenchState {
    /// Single LED or a series string.
    mode: LedMode,
    /// Number of LEDs in the series string (used only in [`LedMode::Series`]).
    count: usize,
    /// Supply (source) voltage `Vs` (V).
    supply_v: f64,
    /// LED forward voltage `Vf` of one LED (V).
    forward_v: f64,
    /// Target LED current `I` (A).
    current_a: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D LED solid (serviced after the panel
    /// draws).
    show_3d_request: bool,
}

impl Default for LedWorkbenchState {
    fn default() -> Self {
        // The canonical textbook indicator: a 5 V rail, a 2.0 V red LED and a
        // 20 mA target give R = (5 - 2) / 0.020 = 150 ohm, with P_led = 40 mW,
        // P_resistor = 60 mW and an LED efficiency of Vf/Vs = 0.4.
        Self {
            mode: LedMode::Single,
            count: 3,
            supply_v: 5.0,
            forward_v: 2.0,
            current_a: 0.020,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the LED Workbench right-side panel. A no-op when the
/// `show_led_workbench` toggle is off.
pub fn draw_led_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_led_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_led_workbench",
        "LED",
        |app, ui| {
            ui.label(
                egui::RichText::new("native series current-limiting resistor sizing · valenx-led")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.led;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Configuration").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.mode, LedMode::Single, "single LED");
                        ui.radio_value(&mut s.mode, LedMode::Series, "series string");
                    });
                    if s.mode == LedMode::Series {
                        ui.horizontal(|ui| {
                            ui.label("LED count n");
                            ui.add(egui::DragValue::new(&mut s.count).speed(1.0));
                        });
                    }

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Supply & LED").strong());
                    ui.horizontal(|ui| {
                        ui.label("supply Vs (V)");
                        ui.add(egui::DragValue::new(&mut s.supply_v).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("forward Vf (V)");
                        ui.add(egui::DragValue::new(&mut s.forward_v).speed(0.05));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Target current").strong());
                    ui.horizontal(|ui| {
                        ui.label("current I (A)");
                        ui.add(egui::DragValue::new(&mut s.current_a).speed(0.001));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_led(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative through-hole LED (cylindrical body with a domed lens) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Operating point").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_led_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.led` borrow is released
    // here): build the LED's 3-D solid and load it.
    if app.led.show_3d_request {
        app.led.show_3d_request = false;
        load_led_3d(app);
    }
}

/// Validate the form, solve the circuit and format the readout.
fn run_led(s: &mut LedWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The validated equivalent single-loop [`LedCircuit`] for the current form,
/// the quantity both the readout and the 3-D gate need. A series string is
/// collapsed to its equivalent circuit (same supply, current and *total*
/// forward voltage), so every derived scalar matches the string exactly.
/// Extracted so it is unit-testable and shared.
fn circuit(s: &LedWorkbenchState) -> Result<LedCircuit, String> {
    match s.mode {
        LedMode::Single => {
            LedCircuit::new(s.supply_v, s.forward_v, s.current_a).map_err(|e| e.to_string())
        }
        LedMode::Series => LedString::new(s.count, s.forward_v, s.supply_v, s.current_a)
            .map_err(|e| e.to_string())?
            .as_circuit()
            .map_err(|e| e.to_string()),
    }
}

/// Solve the circuit and format the full readout, mapping any domain error to
/// a display string. Extracted so it is unit-testable.
fn compute(s: &LedWorkbenchState) -> Result<String, String> {
    let c = circuit(s)?;
    let leds = match s.mode {
        LedMode::Single => 1,
        LedMode::Series => s.count,
    };

    Ok(format!(
        "configuration   : {leds} LED(s) in series\n\
         supply Vs       : {:.3} V\n\
         forward Vf      : {:.3} V each / {:.3} V total\n\
         target current  : {:.4} A\n\n\
         resistor R      : {:.2} Ω\n\
         resistor voltage: {:.3} V\n\
         LED power       : {:.4} W\n\
         resistor power  : {:.4} W\n\
         total power     : {:.4} W\n\
         LED efficiency  : {:.3}",
        c.supply_v,
        s.forward_v,
        c.total_forward_v(),
        c.current_a,
        c.resistor_ohm(),
        c.resistor_voltage(),
        c.led_power_w(),
        c.resistor_power_w(),
        c.total_power_w(),
        c.led_efficiency(),
    ))
}

/// Append an open `z`-axis cylinder (centre axis at `(cx, cy)`, base at
/// `z0`, the given `height` and `radius`, `seg` facets) to the buffers as a
/// side wall plus a flat bottom cap. The top is left open so the next stacked
/// piece (the body, then the dome) closes it.
#[allow(clippy::too_many_arguments)]
fn push_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    cx: f64,
    cy: f64,
    z0: f64,
    height: f64,
    radius: f64,
    seg: usize,
) {
    let bottom = nodes.len();
    for j in 0..seg {
        let theta = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(
            cx + radius * theta.cos(),
            cy + radius * theta.sin(),
            z0,
        ));
    }
    let top = nodes.len();
    for j in 0..seg {
        let theta = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(
            cx + radius * theta.cos(),
            cy + radius * theta.sin(),
            z0 + height,
        ));
    }
    // Side wall.
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[
            bottom + j,
            bottom + jn,
            top + jn,
            bottom + j,
            top + jn,
            top + j,
        ]);
    }
    // Flat bottom cap (outward-facing).
    let center = nodes.len();
    nodes.push(Vector3::new(cx, cy, z0));
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[center, bottom + jn, bottom + j]);
    }
}

/// Append a hemispherical dome (centre `(cx, cy)`, equator at `z0`, the given
/// `radius`, `seg` longitudinal facets, `rings` latitudinal bands) to the
/// buffers, bulging in `+z`. Used as the LED's domed lens cap.
#[allow(clippy::too_many_arguments)]
fn push_dome(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    cx: f64,
    cy: f64,
    z0: f64,
    radius: f64,
    seg: usize,
    rings: usize,
) {
    let ring_base = nodes.len();
    // `rings` latitude rings from the equator (phi = 0) up to just below the
    // pole, then a single apex vertex.
    for r in 0..rings {
        let phi = (r as f64 / rings as f64) * (TAU / 4.0);
        let z = z0 + radius * phi.sin();
        let rr = radius * phi.cos();
        for j in 0..seg {
            let theta = j as f64 / seg as f64 * TAU;
            nodes.push(Vector3::new(
                cx + rr * theta.cos(),
                cy + rr * theta.sin(),
                z,
            ));
        }
    }
    // Bands between successive rings.
    for r in 0..rings - 1 {
        let a = ring_base + r * seg;
        let b = ring_base + (r + 1) * seg;
        for j in 0..seg {
            let jn = (j + 1) % seg;
            tris.extend_from_slice(&[a + j, a + jn, b + jn, a + j, b + jn, b + j]);
        }
    }
    // Apex fan closing the top band.
    let apex = nodes.len();
    nodes.push(Vector3::new(cx, cy, z0 + radius));
    let top_ring = ring_base + (rings - 1) * seg;
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[apex, top_ring + j, top_ring + jn]);
    }
}

/// Build the LED as a triangle [`Mesh`] — a representative through-hole LED:
/// a cylindrical epoxy body capped by a hemispherical lens dome, on a short
/// base flange. Representative geometry (not to scale; the resistor and power
/// numbers are the `valenx-led` result). `None` for an invalid circuit.
fn led_solid_mesh(s: &LedWorkbenchState) -> Option<Mesh> {
    circuit(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    let seg = 24usize;

    // Base flange (a short, slightly wider disc).
    push_cylinder(&mut nodes, &mut tris, 0.0, 0.0, 0.0, 0.08, 0.62, seg);
    // Main cylindrical body.
    push_cylinder(&mut nodes, &mut tris, 0.0, 0.0, 0.08, 0.9, 0.5, seg);
    // Hemispherical lens dome on top of the body.
    push_dome(&mut nodes, &mut tris, 0.0, 0.0, 0.98, 0.5, seg, 6);

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-led");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D LED solid and load it into the central viewport.
fn load_led_3d(app: &mut ValenxApp) {
    let Some(mesh) = led_solid_mesh(&app.led) else {
        app.led.error = Some("LED parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<led>/valenx-led"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"led"}`** product: the canonical LED
/// series circuit built as a 3-D solid, paired with the workbench's own
/// `compute()` readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`LedWorkbenchState::default`].
pub(crate) fn led_product() -> crate::WorkspaceProduct {
    let s = LedWorkbenchState::default();
    let mesh = led_solid_mesh(&s).expect("canonical LED ⇒ solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<led>/valenx-led");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical LED ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "LED (series circuit)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = LedWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_resistor_and_powers() {
        let mut s = LedWorkbenchState::default();
        run_led(&mut s);
        assert!(
            s.error.is_none(),
            "default LED should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("resistor R"));
        assert!(s.result.contains("LED power"));
        assert!(s.result.contains("total power"));
        // Default red LED on 5 V at 20 mA: R = 150 ohm, efficiency Vf/Vs = 0.4.
        assert!(s.result.contains("150.00"));
        assert!(s.result.contains("0.400"));
    }

    #[test]
    fn analyze_rejects_insufficient_headroom() {
        // Vf >= Vs leaves no voltage for the resistor.
        let mut s = LedWorkbenchState {
            forward_v: 5.0,
            ..Default::default()
        };
        run_led(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn series_string_sums_forward_voltage_in_readout() {
        // Three 3.2 V LEDs on a 12 V supply at 20 mA: total Vf = 9.6 V,
        // R = (12 - 9.6) / 0.020 = 120 ohm.
        let mut s = LedWorkbenchState {
            mode: LedMode::Series,
            count: 3,
            supply_v: 12.0,
            forward_v: 3.2,
            current_a: 0.020,
            ..Default::default()
        };
        run_led(&mut s);
        assert!(
            s.error.is_none(),
            "series string should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("3 LED(s) in series"));
        assert!(s.result.contains("9.600 V total"));
        assert!(s.result.contains("120.00"));
    }

    #[test]
    fn resistor_matches_closed_form_ground_truth() {
        // Ground truth: R = (Vs - Vf) / I, hand-computed.
        // 5 V supply, 2.0 V LED, 20 mA => (5 - 2) / 0.020 = 150 ohm exactly,
        // P_led = Vf*I = 0.040 W, P_resistor = (Vs-Vf)*I = 0.060 W.
        let c = LedCircuit::new(5.0, 2.0, 0.020).expect("valid");
        let vs: f64 = 5.0;
        let vf: f64 = 2.0;
        let i: f64 = 0.020;
        let hand_r = (vs - vf) / i;
        assert!((hand_r - 150.0).abs() < 1e-9);
        assert!((c.resistor_ohm() - hand_r).abs() < 1e-9);
        assert!((c.led_power_w() - vf * i).abs() < 1e-12);
        assert!((c.resistor_power_w() - (vs - vf) * i).abs() < 1e-12);
    }

    #[test]
    fn led_mesh_for_default_is_nonempty_and_in_range() {
        let s = LedWorkbenchState::default();
        let mesh = led_solid_mesh(&s).expect("default LED yields a solid");
        assert!(mesh.nodes.len() > 8, "expected flange + body + dome");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&idx| idx < n));
        }
    }

    #[test]
    fn led_mesh_none_for_invalid() {
        let s = LedWorkbenchState {
            forward_v: 10.0,
            ..Default::default()
        };
        assert!(led_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_led_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_led_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_led_workbench = true;
        run_led(&mut app.led);
        draw_workbench(&mut app);
    }
}
