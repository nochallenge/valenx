//! The right-side **Three-Phase Workbench** panel — native balanced
//! three-phase AC power analysis over `valenx-threephase`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_threephase_workbench`,
//! toggled from the View menu. The form sets a balanced load by its wiring
//! (wye / delta), line-to-line voltage, line current and power factor;
//! "Analyze" derives the per-element phase quantities and the full power
//! triangle (real `P`, apparent `S`, reactive `|Q|`) via `valenx-threephase`,
//! and "Show 3-D" loads three transformer-style limb cylinders, spaced 120
//! degrees apart, into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_threephase::{
    apparent_power_from_line, power_from_line, reactive_power_from_line, BalancedLoad, Connection,
};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Three-Phase Workbench.
pub struct ThreePhaseWorkbenchState {
    /// How the three identical load elements are wired.
    connection: Connection,
    /// Line-to-line voltage `V_line` (V, RMS).
    v_line: f64,
    /// Line (conductor) current `I_line` (A, RMS).
    i_line: f64,
    /// Power factor `cos(phi)`, in `[-1, 1]`.
    power_factor: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D limbs (serviced after the panel
    /// draws).
    show_3d_request: bool,
}

impl Default for ThreePhaseWorkbenchState {
    fn default() -> Self {
        // A 400 V (line) wye load drawing 15 A per line at 0.8 pf:
        // P = sqrt(3) * 400 * 15 * 0.8 ~= 8.3 kW, S ~= 10.4 kVA.
        Self {
            connection: Connection::Wye,
            v_line: 400.0,
            i_line: 15.0,
            power_factor: 0.8,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Three-Phase Workbench right-side panel. A no-op when the
/// `show_threephase_workbench` toggle is off.
pub fn draw_threephase_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_threephase_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_threephase_workbench",
        "Three-Phase",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "balanced wye/delta line-phase & power triangle · valenx-threephase",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.threephase;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Connection").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.connection, Connection::Wye, "wye (Y)");
                        ui.radio_value(&mut s.connection, Connection::Delta, "delta (Δ)");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Supply").strong());
                    ui.horizontal(|ui| {
                        let cap = ui.label("line voltage (V)");
                        ui.add(egui::DragValue::new(&mut s.v_line).speed(1.0))
                            .labelled_by(cap.id);
                    });
                    ui.horizontal(|ui| {
                        let cap = ui.label("line current (A)");
                        ui.add(egui::DragValue::new(&mut s.i_line).speed(0.5))
                            .labelled_by(cap.id);
                    });
                    ui.horizontal(|ui| {
                        let cap = ui.label("power factor cos φ");
                        ui.add(egui::DragValue::new(&mut s.power_factor).speed(0.01))
                            .labelled_by(cap.id);
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_threephase(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build three transformer-style limb cylinders, spaced 120° apart, as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Power").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_threephase_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.threephase` borrow is
    // released here): build the limb solid and load it.
    if app.threephase.show_3d_request {
        app.threephase.show_3d_request = false;
        load_phases_3d(app);
    }
}

/// Validate the form, evaluate the load and format the readout.
fn run_threephase(s: &mut ThreePhaseWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`BalancedLoad`] for the current form by deriving
/// the per-element phase quantities from the line values and wiring. The
/// quantity both the readout and the 3-D gate need; extracted so it is
/// unit-testable and shared.
fn balanced_load(s: &ThreePhaseWorkbenchState) -> Result<BalancedLoad, String> {
    let v_phase = s.connection.phase_voltage(s.v_line);
    let i_phase = s.connection.phase_current(s.i_line);
    BalancedLoad::new(s.connection, v_phase, i_phase, s.power_factor).map_err(|e| e.to_string())
}

/// Evaluate the balanced load and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &ThreePhaseWorkbenchState) -> Result<String, String> {
    let load = balanced_load(s)?;
    // The power triangle straight from the validated line quantities.
    let p = power_from_line(s.v_line, s.i_line, s.power_factor).map_err(|e| e.to_string())?;
    let apparent = apparent_power_from_line(s.v_line, s.i_line).map_err(|e| e.to_string())?;
    let reactive =
        reactive_power_from_line(s.v_line, s.i_line, s.power_factor).map_err(|e| e.to_string())?;

    let wiring = match s.connection {
        Connection::Wye => "wye (Y)",
        Connection::Delta => "delta (Δ)",
    };

    Ok(format!(
        "connection      : {wiring}\n\
         line voltage    : {:.2} V\n\
         line current    : {:.2} A\n\
         phase voltage   : {:.2} V\n\
         phase current   : {:.2} A\n\
         power factor    : {:.3}\n\n\
         real power P    : {:.1} W\n\
         apparent power S: {:.1} VA\n\
         reactive |Q|    : {:.1} var\n\
         power per phase : {:.1} W",
        load.line_voltage(),
        load.line_current(),
        load.phase_voltage(),
        load.phase_current(),
        load.power_factor(),
        p,
        apparent,
        reactive,
        load.power_per_phase(),
    ))
}

/// Append a (double-sided) cylinder whose axis runs along `+z`, spanning
/// `base.z ..= base.z + length` with circle centre `(base.x, base.y)`.
fn push_cyl_z(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: Vector3<f64>,
    length: f64,
    r: f64,
    seg: usize,
) {
    let (z0, z1) = (base.z, base.z + length);
    let lo = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(base.x + r * a.cos(), base.y + r * a.sin(), z0));
    }
    let hi = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(base.x + r * a.cos(), base.y + r * a.sin(), z1));
    }
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[
            lo + j,
            hi + j,
            hi + jn,
            lo + j,
            hi + jn,
            lo + jn,
            lo + j,
            hi + jn,
            hi + j,
            lo + j,
            lo + jn,
            hi + jn,
        ]);
    }
}

/// Build the three-phase assembly as a triangle [`Mesh`] — three identical
/// vertical limb cylinders standing on a common ring, spaced 120 degrees
/// apart (the textbook symmetric layout of a balanced set / three-limb
/// core). Representative geometry (not to scale; the power numbers are the
/// `valenx-threephase` result). `None` for an invalid configuration.
fn phases_solid_mesh(s: &ThreePhaseWorkbenchState) -> Option<Mesh> {
    balanced_load(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Three limbs on a circle of this radius, each a vertical cylinder.
    let ring = 0.5;
    let limb_r = 0.16;
    let limb_h = 1.0;
    for k in 0..3 {
        let a = k as f64 / 3.0 * TAU;
        push_cyl_z(
            &mut nodes,
            &mut tris,
            Vector3::new(ring * a.cos(), ring * a.sin(), 0.1),
            limb_h,
            limb_r,
            24,
        );
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-threephase");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D limb solid and load it into the central viewport.
fn load_phases_3d(app: &mut ValenxApp) {
    let Some(mesh) = phases_solid_mesh(&app.threephase) else {
        app.threephase.error =
            Some("load parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<phases>/valenx-threephase"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"threephase"}`** product: a DATA-ONLY text
/// card of the balanced three-phase workbench's `compute()` readout rows (see
/// [`crate::products_registry`]). A line/phase voltage-current and power
/// result has no characteristic shape — the panel's three cylindrical limbs on
/// a ring are an abstract stand-in for "three phases", not a real fabricated
/// object — so the bridge product is right-sized to a card (`mesh: None`)
/// carrying just the readout (the confidence badge is appended centrally). The
/// panel's "Show 3-D" button still builds that representative schematic into
/// the central viewport. Pure — driven off
/// [`ThreePhaseWorkbenchState::default`].
pub(crate) fn threephase_product() -> crate::WorkspaceProduct {
    let s = ThreePhaseWorkbenchState::default();
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical three-phase ⇒ readout computes"),
    );
    crate::WorkspaceProduct {
        title: "Three-phase (balanced load)".into(),
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
        let s = ThreePhaseWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_power_triangle() {
        let mut s = ThreePhaseWorkbenchState::default();
        run_threephase(&mut s);
        assert!(
            s.error.is_none(),
            "default load should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("real power P"));
        assert!(s.result.contains("apparent power S"));
        assert!(s.result.contains("reactive |Q|"));
        // Ground truth: P = sqrt(3) * 400 * 15 * 0.8 = 8313.8 W (1 dp).
        assert!(s.result.contains("8313.8"), "readout was:\n{}", s.result);
    }

    #[test]
    fn analyze_rejects_zero_voltage() {
        let mut s = ThreePhaseWorkbenchState {
            v_line: 0.0,
            ..Default::default()
        };
        run_threephase(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_rejects_power_factor_out_of_range() {
        let mut s = ThreePhaseWorkbenchState {
            power_factor: 1.5,
            ..Default::default()
        };
        run_threephase(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn wye_line_voltage_is_sqrt3_times_phase_and_power_matches() {
        // Ground truth #1: for a wye load the line-to-line voltage is
        // sqrt(3) times the per-element phase voltage.
        let s = ThreePhaseWorkbenchState::default();
        let load = balanced_load(&s).expect("default load is valid");
        let sqrt3: f64 = 3.0_f64.sqrt();
        assert!(
            (load.line_voltage() - sqrt3 * load.phase_voltage()).abs() < 1e-9,
            "V_line {} != sqrt(3) * V_phase {}",
            load.line_voltage(),
            load.phase_voltage()
        );
        // Ground truth #2: P = sqrt(3) * V_line * I_line * cos(phi),
        // hand-computed for 400 V, 15 A, 0.8 pf.
        let p = power_from_line(s.v_line, s.i_line, s.power_factor).unwrap();
        let expected: f64 = sqrt3 * 400.0 * 15.0 * 0.8;
        assert!((p - expected).abs() < 1e-6, "P {p} != {expected}");
        // And it equals three times the per-phase power for a balanced set.
        assert!((p - 3.0 * load.power_per_phase()).abs() < 1e-6);
    }

    #[test]
    fn phases_mesh_for_default_is_nonempty_and_in_range() {
        let s = ThreePhaseWorkbenchState::default();
        let mesh = phases_solid_mesh(&s).expect("default load yields a solid");
        assert!(mesh.nodes.len() > 8, "expected three limb cylinders");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn phases_mesh_none_for_invalid() {
        let s = ThreePhaseWorkbenchState {
            i_line: 0.0,
            ..Default::default()
        };
        assert!(phases_solid_mesh(&s).is_none());
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
            draw_threephase_workbench(app, ctx);
        });
    }

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_threephase_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_threephase_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_threephase_workbench = true;
        run_threephase(&mut app.threephase);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        let mut app = ValenxApp::default();
        app.show_threephase_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 3,
            "expected the numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by its caption (AI-drivable name)"
        );
        for caption in ["line voltage (V)", "power factor cos φ"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
