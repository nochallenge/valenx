//! The right-side **Transmission Line Workbench** panel — native lossless
//! RF transmission-line reflection / standing-wave analysis over
//! `valenx-transmissionline`.
//!
//! Mirrors the Heat Transfer / Bearing workbenches: a resizable
//! [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_transmissionline_workbench`, toggled from the
//! View menu. The form sets the line's characteristic impedance `Z0` and a
//! purely resistive termination (a finite load, a short, or an idealised
//! open); "Analyze" evaluates the voltage reflection coefficient
//! `gamma = (ZL - Z0) / (ZL + Z0)` and the standard standing-wave figures
//! of merit (VSWR, return loss, mismatch loss, reflected / transmitted
//! power), and "Show 3-D" loads a representative coaxial line (inner
//! conductor inside a concentric outer shield) into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_transmissionline::{Line, Load};

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// The kind of termination presented to the line.
///
/// Maps to [`Load`]: a finite resistive load (its magnitude taken from
/// `load_ohms`), an idealised open circuit, or a short circuit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LoadKind {
    /// A finite resistive load, magnitude taken from `load_ohms`.
    Resistive,
    /// An idealised open circuit (`ZL -> infinity`), `gamma = +1`.
    Open,
    /// A short circuit (`ZL = 0`), `gamma = -1`.
    Short,
}

/// Persistent form + result state for the Transmission Line Workbench.
pub struct TransmissionLineWorkbenchState {
    /// Line characteristic impedance `Z0` (ohms). Finite and `> 0`.
    z0_ohms: f64,
    /// Termination kind selected in the form.
    load_kind: LoadKind,
    /// Resistive load magnitude `ZL` (ohms), used when `load_kind` is
    /// [`LoadKind::Resistive`].
    load_ohms: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D coax solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for TransmissionLineWorkbenchState {
    fn default() -> Self {
        // A catalogued 50 Ω line feeding a 75 Ω resistive load (the
        // classic 50→75 mismatch): gamma = (75-50)/(75+50) = 0.2,
        // VSWR = 1.5, return loss ≈ 13.98 dB.
        Self {
            z0_ohms: 50.0,
            load_kind: LoadKind::Resistive,
            load_ohms: 75.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Transmission Line Workbench right-side panel. A no-op when the
/// `show_transmissionline_workbench` toggle is off.
pub fn draw_transmissionline_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_transmissionline_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_transmissionline_workbench",
        "Transmission Line",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native lossless RF reflection / VSWR · valenx-transmissionline",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.transmissionline;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Line").strong());
                    ui.horizontal(|ui| {
                        ui.label("Z₀ (Ω)");
                        ui.add(egui::DragValue::new(&mut s.z0_ohms).speed(1.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Termination").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.load_kind, LoadKind::Resistive, "resistive");
                        ui.radio_value(&mut s.load_kind, LoadKind::Short, "short");
                        ui.radio_value(&mut s.load_kind, LoadKind::Open, "open");
                    });
                    if s.load_kind == LoadKind::Resistive {
                        ui.horizontal(|ui| {
                            ui.label("Z_L (Ω)");
                            ui.add(egui::DragValue::new(&mut s.load_ohms).speed(1.0));
                        });
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_transmissionline(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative coaxial line (an inner conductor inside a concentric outer shield) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Reflection").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_transmissionline_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.transmissionline`
    // borrow is released here): build the coax 3-D solid and load it.
    if app.transmissionline.show_3d_request {
        app.transmissionline.show_3d_request = false;
        load_line_3d(app);
    }
}

/// Validate the form, evaluate the line and format the readout.
fn run_transmissionline(s: &mut TransmissionLineWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`Load`] for the current form selection. Extracted
/// so it is shared by `compute` and the 3-D gate.
fn build_load(s: &TransmissionLineWorkbenchState) -> Result<Load, String> {
    match s.load_kind {
        LoadKind::Resistive => Load::resistive(s.load_ohms).map_err(|e| e.to_string()),
        LoadKind::Short => Ok(Load::short()),
        LoadKind::Open => Ok(Load::Open),
    }
}

/// Evaluate the line + load and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &TransmissionLineWorkbenchState) -> Result<String, String> {
    let line = Line::from_z0(s.z0_ohms).map_err(|e| e.to_string())?;
    let load = build_load(s)?;
    let r = line.reflection(load);

    // Diverging quantities (VSWR / return loss) are reported as `None` by
    // the model; render them as the limit they represent.
    let vswr = match r.vswr() {
        Some(v) => format!("{v:.4}"),
        None => "∞ (total reflection)".to_string(),
    };
    let return_loss = match r.return_loss_db() {
        Some(rl) => format!("{rl:.3} dB"),
        None => "∞ (matched)".to_string(),
    };
    let mismatch_loss = match r.mismatch_loss_db() {
        Some(ml) => format!("{ml:.3} dB"),
        None => "∞ (total reflection)".to_string(),
    };
    let load_str = match r.load_ohms() {
        Some(zl) => format!("{zl:.1} Ω"),
        None => "open (∞ Ω)".to_string(),
    };

    let gamma = r.gamma();
    let mag = r.gamma_magnitude();
    let p_refl = r.power_reflected_fraction() * 100.0;
    let p_trans = r.power_transmitted_fraction() * 100.0;

    Ok(format!(
        "Z₀              : {z0:.1} Ω\n\
         load Z_L        : {load_str}\n\n\
         gamma           : {gamma:.4}\n\
         |gamma|         : {mag:.4}\n\
         VSWR            : {vswr}\n\
         return loss     : {return_loss}\n\
         mismatch loss   : {mismatch_loss}\n\
         power reflected : {p_refl:.2} %\n\
         power to load   : {p_trans:.2} %",
        z0 = s.z0_ohms,
    ))
}

/// Append a (double-sided) cylinder whose axis runs along `+x`, spanning
/// `base.x ..= base.x + length` with circle centre `(base.y, base.z)`.
fn push_cyl_x(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: Vector3<f64>,
    length: f64,
    r: f64,
    seg: usize,
) {
    let (x0, x1) = (base.x, base.x + length);
    let lo = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x0, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    let hi = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x1, base.y + r * a.cos(), base.z + r * a.sin()));
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

/// Build the line as a triangle [`Mesh`] — a representative coaxial cable:
/// a slender inner conductor running along `x` inside a concentric, wider
/// outer shield. Representative geometry (not to scale; the reflection /
/// standing-wave numbers are the `valenx-transmissionline` result). `None`
/// for an invalid configuration (e.g. a non-positive `Z0`).
fn coax_solid_mesh(s: &TransmissionLineWorkbenchState) -> Option<Mesh> {
    // Reuse the same validation the readout uses; bail if it would fail.
    compute(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    let length = 1.6;
    let base = Vector3::new(-length * 0.5, 0.0, 0.6);
    let inner_r = 0.08;
    let outer_r = 0.26;

    // Inner conductor (slender concentric x-cylinder).
    push_cyl_x(&mut nodes, &mut tris, base, length, inner_r, 24);
    // Outer shield (wider concentric x-cylinder, sharing the axis).
    push_cyl_x(&mut nodes, &mut tris, base, length, outer_r, 48);

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-transmissionline");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D coax solid and load it into the central viewport.
fn load_line_3d(app: &mut ValenxApp) {
    let Some(mesh) = coax_solid_mesh(&app.transmissionline) else {
        app.transmissionline.error =
            Some("line parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<line>/valenx-transmissionline"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"transmissionline"}`** product: the
/// canonical coaxial line built as a 3-D solid, paired with the workbench's
/// own `compute()` readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`TransmissionLineWorkbenchState::default`].
pub(crate) fn transmissionline_product() -> crate::WorkspaceProduct {
    let s = TransmissionLineWorkbenchState::default();
    let mesh = coax_solid_mesh(&s).expect("canonical line ⇒ coax solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<line>/valenx-transmissionline");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical line ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Transmission line (coax)".into(),
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
        let s = TransmissionLineWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_gamma_vswr_and_return_loss() {
        let mut s = TransmissionLineWorkbenchState::default();
        run_transmissionline(&mut s);
        assert!(
            s.error.is_none(),
            "default line should analyze: {:?}",
            s.error
        );
        // 50→75 mismatch: gamma = 0.2, VSWR = 1.5, RL ≈ 13.979 dB.
        assert!(s.result.contains("gamma"));
        assert!(s.result.contains("0.2000"));
        assert!(s.result.contains("VSWR"));
        assert!(s.result.contains("1.5000"));
        assert!(s.result.contains("return loss"));
        assert!(s.result.contains("13.979"));
    }

    #[test]
    fn analyze_rejects_non_positive_z0() {
        let mut s = TransmissionLineWorkbenchState {
            z0_ohms: 0.0,
            ..Default::default()
        };
        run_transmissionline(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn open_and_short_report_total_reflection() {
        // Both an open and a short fully reflect: VSWR diverges (rendered
        // as the total-reflection limit, not a finite number).
        for kind in [LoadKind::Open, LoadKind::Short] {
            let mut s = TransmissionLineWorkbenchState {
                load_kind: kind,
                ..Default::default()
            };
            run_transmissionline(&mut s);
            assert!(s.error.is_none(), "{kind:?} should analyze");
            assert!(s.result.contains("total reflection"), "{kind:?}");
            // Total reflection ⇒ 0 dB return loss.
            assert!(s.result.contains("0.000 dB"), "{kind:?}");
        }
    }

    #[test]
    fn gamma_and_vswr_match_hand_computed_ground_truth() {
        // Ground truth (textbook): a 50 Ω line into a 75 Ω resistive load.
        //   gamma = (ZL - Z0)/(ZL + Z0) = (75 - 50)/(75 + 50) = 0.2
        //   VSWR  = (1 + |gamma|)/(1 - |gamma|) = 1.2/0.8 = 1.5
        let line = Line::from_z0(50.0).unwrap();
        let r = line.reflection(Load::resistive(75.0).unwrap());
        let gamma = r.gamma();
        assert!((gamma - 0.2).abs() < 1e-12, "gamma {gamma} != 0.2");
        let vswr = r.vswr().unwrap();
        assert!((vswr - 1.5).abs() < 1e-12, "VSWR {vswr} != 1.5");
        // 4 % of incident power reflected, 96 % delivered to the load.
        assert!((r.power_reflected_fraction() - 0.04).abs() < 1e-12);
        assert!((r.power_transmitted_fraction() - 0.96).abs() < 1e-12);
    }

    #[test]
    fn coax_mesh_for_default_is_nonempty_and_in_range() {
        let s = TransmissionLineWorkbenchState::default();
        let mesh = coax_solid_mesh(&s).expect("default line yields a solid");
        assert!(mesh.nodes.len() > 8, "expected inner conductor + shield");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn coax_mesh_none_for_invalid() {
        let s = TransmissionLineWorkbenchState {
            z0_ohms: -50.0,
            ..Default::default()
        };
        assert!(coax_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_transmissionline_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_transmissionline_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_transmissionline_workbench = true;
        run_transmissionline(&mut app.transmissionline);
        draw_workbench(&mut app);
    }
}
