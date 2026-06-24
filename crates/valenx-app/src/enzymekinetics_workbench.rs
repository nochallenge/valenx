//! The right-side **Enzyme Kinetics Workbench** panel — native
//! Michaelis-Menten / reversible-inhibition rate-law evaluation over
//! `valenx-enzymekinetics`.
//!
//! Mirrors the Heat Transfer / Heat Exchanger workbenches: a resizable
//! [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_enzymekinetics_workbench`, toggled from the View
//! menu. The form sets the uninhibited parameters of a
//! [`valenx_enzymekinetics::MichaelisMenten`] enzyme (`Vmax`, `Km`), a
//! substrate concentration `[S]`, an inhibition mode and (when a mode is
//! selected) an inhibitor concentration `[I]` with the relevant inhibition
//! constant; "Analyze" reports the initial velocity, fractional saturation
//! and — for an inhibited run — the apparent `Vmax` and `Km`, and "Show
//! 3-D vessel" loads a representative stirred-tank bioreactor solid into the
//! central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_enzymekinetics::inhibition::{Competitive, Noncompetitive, Uncompetitive};
use valenx_enzymekinetics::MichaelisMenten;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Reversible-inhibition mode selected in the form.
///
/// `None` is bare Michaelis-Menten; the other three are the standard
/// single-inhibitor modes, each evaluated as Michaelis-Menten with apparent
/// parameters by the matching `valenx-enzymekinetics` type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InhibitionMode {
    /// No inhibitor — bare Michaelis-Menten.
    None,
    /// Competitive: inhibitor raises apparent `Km`, leaves `Vmax`.
    Competitive,
    /// Noncompetitive: inhibitor lowers apparent `Vmax`, leaves `Km`.
    Noncompetitive,
    /// Uncompetitive: inhibitor lowers both apparent `Vmax` and `Km`.
    Uncompetitive,
}

impl InhibitionMode {
    /// Short human-readable label for the readout.
    fn label(self) -> &'static str {
        match self {
            InhibitionMode::None => "none (Michaelis-Menten)",
            InhibitionMode::Competitive => "competitive",
            InhibitionMode::Noncompetitive => "noncompetitive",
            InhibitionMode::Uncompetitive => "uncompetitive",
        }
    }
}

/// Persistent form + result state for the Enzyme Kinetics Workbench.
pub struct EnzymeKineticsWorkbenchState {
    /// Maximal velocity `Vmax` (µmol/min).
    vmax_umol_per_min: f64,
    /// Michaelis constant `Km` (mM).
    km_mm: f64,
    /// Substrate concentration `[S]` (mM).
    s_mm: f64,
    /// Inhibition mode.
    mode: InhibitionMode,
    /// Inhibitor concentration `[I]` (mM); used only when `mode != None`.
    i_mm: f64,
    /// Inhibition constant `Ki` (mM); `Ki'` for the uncompetitive mode.
    /// Used only when `mode != None`.
    ki_mm: f64,
    /// Formatted result readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D vessel solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for EnzymeKineticsWorkbenchState {
    fn default() -> Self {
        // A representative enzyme right at half-saturation: Vmax = 100
        // µmol/min, Km = 5 mM, with [S] = Km = 5 mM and no inhibitor gives
        // exactly v = Vmax/2 = 50 µmol/min and a fractional saturation of
        // 0.5 — the defining property of the Michaelis constant.
        Self {
            vmax_umol_per_min: 100.0,
            km_mm: 5.0,
            s_mm: 5.0,
            mode: InhibitionMode::None,
            i_mm: 2.0,
            ki_mm: 1.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Enzyme Kinetics Workbench right-side panel. A no-op when the
/// `show_enzymekinetics_workbench` toggle is off.
pub fn draw_enzymekinetics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_enzymekinetics_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_enzymekinetics_workbench",
        "Enzyme Kinetics",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native Michaelis-Menten / inhibition rate laws · valenx-enzymekinetics",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.enzymekinetics;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Enzyme (uninhibited)").strong());
                    ui.horizontal(|ui| {
                        let cap = ui.label("Vmax (µmol/min)");
                        ui.add(egui::DragValue::new(&mut s.vmax_umol_per_min).speed(1.0))
                            .labelled_by(cap.id);
                    });
                    ui.horizontal(|ui| {
                        let cap = ui.label("Km (mM)");
                        ui.add(egui::DragValue::new(&mut s.km_mm).speed(0.1))
                            .labelled_by(cap.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Substrate").strong());
                    ui.horizontal(|ui| {
                        let cap = ui.label("[S] (mM)");
                        ui.add(egui::DragValue::new(&mut s.s_mm).speed(0.1))
                            .labelled_by(cap.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Inhibition").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.mode, InhibitionMode::None, "None");
                        ui.radio_value(&mut s.mode, InhibitionMode::Competitive, "Comp.");
                    });
                    ui.horizontal(|ui| {
                        ui.radio_value(
                            &mut s.mode,
                            InhibitionMode::Noncompetitive,
                            "Noncomp.",
                        );
                        ui.radio_value(
                            &mut s.mode,
                            InhibitionMode::Uncompetitive,
                            "Uncomp.",
                        );
                    });
                    if s.mode != InhibitionMode::None {
                        ui.horizontal(|ui| {
                            let cap = ui.label("[I] (mM)");
                            ui.add(egui::DragValue::new(&mut s.i_mm).speed(0.1))
                                .labelled_by(cap.id);
                        });
                        let ki_label = if s.mode == InhibitionMode::Uncompetitive {
                            "Ki' (mM)"
                        } else {
                            "Ki (mM)"
                        };
                        ui.horizontal(|ui| {
                            let cap = ui.label(ki_label);
                            ui.add(egui::DragValue::new(&mut s.ki_mm).speed(0.1))
                                .labelled_by(cap.id);
                        });
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_enzyme(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D vessel").strong())
                        .on_hover_text(
                            "Build a representative stirred-tank bioreactor (vessel, stirrer shaft and a few substrate parcels) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Rate").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_enzymekinetics_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.enzymekinetics` borrow
    // is released here): build the vessel's 3-D solid and load it.
    if app.enzymekinetics.show_3d_request {
        app.enzymekinetics.show_3d_request = false;
        load_vessel_3d(app);
    }
}

/// Validate the form, evaluate the rate law and format the readout.
fn run_enzyme(s: &mut EnzymeKineticsWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The evaluated kinetics `(velocity, saturation, apparent_vmax,
/// apparent_km)` for the current form — the quantities both the readout and
/// the 3-D gate need. Extracted so it is unit-testable and shared.
///
/// The velocity and saturation come from the selected inhibition mode (bare
/// Michaelis-Menten for `None`); the apparent `Vmax` / `Km` are the
/// inhibited effective parameters (equal to the uninhibited values when
/// `mode == None`). Any domain error is mapped to a display string.
fn evaluate(s: &EnzymeKineticsWorkbenchState) -> Result<(f64, f64, f64, f64), String> {
    let base = MichaelisMenten::new(s.vmax_umol_per_min, s.km_mm).map_err(|e| e.to_string())?;
    match s.mode {
        InhibitionMode::None => {
            let v = base.velocity(s.s_mm).map_err(|e| e.to_string())?;
            let sat = base.saturation(s.s_mm).map_err(|e| e.to_string())?;
            Ok((v, sat, base.vmax(), base.km()))
        }
        InhibitionMode::Competitive => {
            let m = Competitive::new(base, s.ki_mm).map_err(|e| e.to_string())?;
            let v = m.velocity(s.s_mm, s.i_mm).map_err(|e| e.to_string())?;
            let vmax_app = m.apparent_vmax(s.i_mm).map_err(|e| e.to_string())?;
            let km_app = m.apparent_km(s.i_mm).map_err(|e| e.to_string())?;
            Ok((v, v / vmax_app, vmax_app, km_app))
        }
        InhibitionMode::Noncompetitive => {
            let m = Noncompetitive::new(base, s.ki_mm).map_err(|e| e.to_string())?;
            let v = m.velocity(s.s_mm, s.i_mm).map_err(|e| e.to_string())?;
            let vmax_app = m.apparent_vmax(s.i_mm).map_err(|e| e.to_string())?;
            let km_app = m.apparent_km(s.i_mm).map_err(|e| e.to_string())?;
            Ok((v, v / vmax_app, vmax_app, km_app))
        }
        InhibitionMode::Uncompetitive => {
            let m = Uncompetitive::new(base, s.ki_mm).map_err(|e| e.to_string())?;
            let v = m.velocity(s.s_mm, s.i_mm).map_err(|e| e.to_string())?;
            let vmax_app = m.apparent_vmax(s.i_mm).map_err(|e| e.to_string())?;
            let km_app = m.apparent_km(s.i_mm).map_err(|e| e.to_string())?;
            Ok((v, v / vmax_app, vmax_app, km_app))
        }
    }
}

/// Evaluate the rate law and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &EnzymeKineticsWorkbenchState) -> Result<String, String> {
    let (v, sat, vmax_app, km_app) = evaluate(s)?;

    let header = format!(
        "Vmax / Km       : {:.2} µmol/min / {:.3} mM\n\
         substrate [S]   : {:.3} mM\n\
         inhibition mode : {}\n",
        s.vmax_umol_per_min,
        s.km_mm,
        s.s_mm,
        s.mode.label(),
    );

    let inhibitor = if s.mode == InhibitionMode::None {
        String::new()
    } else {
        let ki_name = if s.mode == InhibitionMode::Uncompetitive {
            "Ki'"
        } else {
            "Ki "
        };
        format!(
            "inhibitor [I]   : {:.3} mM\n\
             {ki_name}            : {:.3} mM\n\
             apparent Vmax   : {:.3} µmol/min\n\
             apparent Km     : {:.3} mM\n",
            s.i_mm, s.ki_mm, vmax_app, km_app,
        )
    };

    Ok(format!(
        "{header}{inhibitor}\n\
         velocity v      : {v:.3} µmol/min\n\
         saturation v/Vm : {sat:.3}",
    ))
}

/// Append an outward-facing box (centre `c`, half-extents `h`) to the
/// buffers.
fn push_box(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    h: Vector3<f64>,
) {
    let base = nodes.len();
    let signs = [
        (-1.0, -1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, 1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0),
        (-1.0, 1.0, 1.0),
    ];
    for (sx, sy, sz) in signs {
        nodes.push(c + Vector3::new(sx * h.x, sy * h.y, sz * h.z));
    }
    let faces = [
        [1, 2, 6, 5],
        [0, 4, 7, 3],
        [3, 7, 6, 2],
        [0, 1, 5, 4],
        [4, 5, 6, 7],
        [0, 3, 2, 1],
    ];
    for f in faces {
        tris.extend_from_slice(&[
            base + f[0],
            base + f[1],
            base + f[2],
            base + f[0],
            base + f[2],
            base + f[3],
        ]);
    }
}

/// Append a (double-sided) cylinder whose axis runs along `+z`, spanning
/// `base.z ..= base.z + height` with circle centre `(base.x, base.y)`.
fn push_cyl_z(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: Vector3<f64>,
    height: f64,
    r: f64,
    seg: usize,
) {
    let (z0, z1) = (base.z, base.z + height);
    let bot = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(base.x + r * a.cos(), base.y + r * a.sin(), z0));
    }
    let top = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(base.x + r * a.cos(), base.y + r * a.sin(), z1));
    }
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[
            bot + j,
            top + j,
            top + jn,
            bot + j,
            top + jn,
            bot + jn,
            bot + j,
            top + jn,
            top + j,
            bot + j,
            bot + jn,
            top + jn,
        ]);
    }
}

/// Build the bioreactor as a triangle [`Mesh`] — an upright cylindrical
/// vessel with a stirrer shaft on the axis, two impeller blades and a few
/// suspended substrate parcels, on a base. Representative geometry (not to
/// scale; the kinetics numbers are the `valenx-enzymekinetics` result).
/// `None` for an invalid parameter set.
fn vessel_solid_mesh(s: &EnzymeKineticsWorkbenchState) -> Option<Mesh> {
    // Gate the 3-D build on a valid enzyme + a successful evaluation.
    evaluate(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Vessel wall — the large upright cylinder.
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.12),
        1.0,
        0.45,
        32,
    );
    // Vessel floor.
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.1),
        0.04,
        0.45,
        32,
    );
    // Stirrer shaft on the axis, from above the lid down into the broth.
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.2),
        1.1,
        0.035,
        12,
    );
    // Two impeller blades low on the shaft.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.32),
        Vector3::new(0.22, 0.03, 0.05),
    );
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.32),
        Vector3::new(0.03, 0.22, 0.05),
    );
    // A few suspended substrate parcels inside the broth.
    let parcels = [
        (0.2, 0.1, 0.55),
        (-0.18, 0.16, 0.7),
        (0.05, -0.22, 0.45),
        (-0.12, -0.1, 0.85),
    ];
    for (px, py, pz) in parcels {
        push_box(
            &mut nodes,
            &mut tris,
            Vector3::new(px, py, pz),
            Vector3::new(0.05, 0.05, 0.05),
        );
    }
    // Base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.05),
        Vector3::new(0.55, 0.55, 0.05),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-enzymekinetics");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D vessel solid and load it into the central viewport.
fn load_vessel_3d(app: &mut ValenxApp) {
    let Some(mesh) = vessel_solid_mesh(&app.enzymekinetics) else {
        app.enzymekinetics.error =
            Some("enzyme parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<vessel>/valenx-enzymekinetics"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"enzymekinetics"}`** product: the
/// canonical stirred-tank bioreactor solid (the panel's "Show 3-D vessel"
/// geometry) paired with the workbench's own Michaelis-Menten rate headline
/// numbers, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`EnzymeKineticsWorkbenchState::default`].
///
/// The readout rows mirror the panel's `compute()` rate readout.
pub(crate) fn enzymekinetics_product() -> crate::WorkspaceProduct {
    let s = EnzymeKineticsWorkbenchState::default();
    let mesh = vessel_solid_mesh(&s).expect("default half-saturated enzyme ⇒ a 3-D vessel");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<vessel>/valenx-enzymekinetics");
    let readout = compute(&s).expect("default half-saturated enzyme ⇒ a valid readout");
    let lines = crate::products_registry::lines_from_readout(&readout);
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Enzyme Kinetics".into(),
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
    fn default_state_is_idle() {
        let s = EnzymeKineticsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_velocity_and_saturation() {
        let mut s = EnzymeKineticsWorkbenchState::default();
        run_enzyme(&mut s);
        assert!(
            s.error.is_none(),
            "default enzyme should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("velocity v"));
        assert!(s.result.contains("saturation v/Vm"));
        // [S] = Km with no inhibitor → exactly half-saturation: v = 50,
        // v/Vmax = 0.5 (the defining property of Km).
        assert!(s.result.contains("50.000"));
        assert!(s.result.contains("0.500"));
    }

    #[test]
    fn analyze_rejects_zero_km() {
        let mut s = EnzymeKineticsWorkbenchState {
            km_mm: 0.0,
            ..Default::default()
        };
        run_enzyme(&mut s);
        assert!(s.error.is_some());
    }

    /// Ground truth: at `[S] = Km` with no inhibitor the Michaelis-Menten
    /// velocity is exactly `Vmax / 2` and the fractional saturation is
    /// exactly `0.5`, for any positive parameters.
    #[test]
    fn half_vmax_at_km_is_exact() {
        let s = EnzymeKineticsWorkbenchState {
            vmax_umol_per_min: 100.0,
            km_mm: 5.0,
            s_mm: 5.0,
            mode: InhibitionMode::None,
            ..Default::default()
        };
        let (v, sat, vmax_app, km_app) = evaluate(&s).expect("valid enzyme");
        assert!((v - 50.0).abs() < 1e-9, "v = {v}, want 50");
        assert!((sat - 0.5).abs() < 1e-9, "saturation = {sat}");
        // No inhibitor → apparent parameters equal the uninhibited ones.
        assert!((vmax_app - 100.0).abs() < 1e-9, "Vmax_app = {vmax_app}");
        assert!((km_app - 5.0).abs() < 1e-9, "Km_app = {km_app}");
    }

    /// A competitive inhibitor raises the apparent `Km` but leaves the
    /// apparent `Vmax` untouched — the canonical diagnostic signature.
    #[test]
    fn competitive_raises_apparent_km_only() {
        let s = EnzymeKineticsWorkbenchState {
            mode: InhibitionMode::Competitive,
            i_mm: 2.0,
            ki_mm: 1.0,
            ..Default::default()
        };
        let (_v, _sat, vmax_app, km_app) = evaluate(&s).expect("valid");
        // Km_app = Km*(1 + I/Ki) = 5*(1 + 2/1) = 15; Vmax unchanged at 100.
        assert!((km_app - 15.0).abs() < 1e-9, "Km_app = {km_app}");
        assert!((vmax_app - 100.0).abs() < 1e-9, "Vmax_app = {vmax_app}");
    }

    #[test]
    fn vessel_mesh_for_default_is_nonempty_and_in_range() {
        let s = EnzymeKineticsWorkbenchState::default();
        let mesh = vessel_solid_mesh(&s).expect("default enzyme yields a solid");
        assert!(
            mesh.nodes.len() > 8,
            "expected vessel + shaft + blades + parcels + base"
        );
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn vessel_mesh_none_for_invalid() {
        let s = EnzymeKineticsWorkbenchState {
            km_mm: 0.0,
            ..Default::default()
        };
        assert!(vessel_solid_mesh(&s).is_none());
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
            draw_enzymekinetics_workbench(app, ctx);
        });
    }

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_enzymekinetics_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_enzymekinetics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_enzymekinetics_workbench = true;
        run_enzyme(&mut app.enzymekinetics);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        let mut app = ValenxApp::default();
        app.show_enzymekinetics_workbench = true;
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
        for caption in ["Vmax (µmol/min)", "Km (mM)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
