//! The right-side **Weir Flow Workbench** panel — native sharp-crested
//! open-channel weir gauging over `valenx-weir`.
//!
//! Mirrors the Orifice Meter / Heat Transfer workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_weir_workbench`,
//! toggled from the View menu. The form picks a weir family (rectangular
//! crest or triangular V-notch), its geometry (crest length `L`, or full
//! vertex angle `theta`) and a discharge coefficient `Cd`, plus the
//! upstream head `H`; "Analyze" evaluates the textbook closed form —
//! `Q = Cd·(2/3)·√(2g)·L·H^(3/2)` for the rectangular crest and
//! `Q = Cd·(8/15)·√(2g)·tan(theta/2)·H^(5/2)` for the V-notch — reporting
//! the discharge `Q` and round-tripping it back through the inverse rating
//! curve `head_for_discharge(Q)`, and "Show 3-D" loads a representative
//! weir plate (an upright plate with a rectangular or V notch cut into its
//! crest) into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_weir::{RectangularWeir, VNotchWeir};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// The two sharp-crested weir families this workbench gauges.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WeirKind {
    /// A horizontal rectangular crest of width `L`, with `Q ∝ H^(3/2)`.
    Rectangular,
    /// A triangular V-notch of full vertex angle `theta`, with
    /// `Q ∝ H^(5/2)`.
    VNotch,
}

/// Persistent form + result state for the Weir Flow Workbench.
pub struct WeirWorkbenchState {
    /// Selected weir family.
    kind: WeirKind,
    /// Rectangular crest length (weir width) `L` (m).
    crest_length_m: f64,
    /// V-notch full vertex angle `theta` (degrees in the form; converted
    /// to radians for the solver).
    vertex_angle_deg: f64,
    /// Dimensionless discharge coefficient `Cd`.
    discharge_coefficient: f64,
    /// Upstream head `H` above the crest / vertex (m).
    head_m: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D weir plate (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for WeirWorkbenchState {
    fn default() -> Self {
        // The crate's canonical rectangular case: a 2 m-wide crest,
        // Cd = 0.62, at 0.30 m of head -> Q ~ 0.6016 m^3/s.
        Self {
            kind: WeirKind::Rectangular,
            crest_length_m: 2.0,
            vertex_angle_deg: 90.0,
            discharge_coefficient: 0.62,
            head_m: 0.30,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Weir Flow Workbench right-side panel. A no-op when the
/// `show_weir_workbench` toggle is off.
pub fn draw_weir_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_weir_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_weir_workbench",
        "Weir Flow",
        |app, ui| {
            ui.label(
                egui::RichText::new("native sharp-crested open-channel weir gauging · valenx-weir")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.weir;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Weir type").strong());
                    ui.radio_value(&mut s.kind, WeirKind::Rectangular, "rectangular crest (Q∝H³ᐟ²)");
                    ui.radio_value(&mut s.kind, WeirKind::VNotch, "triangular V-notch (Q∝H⁵ᐟ²)");

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Geometry").strong());
                    match s.kind {
                        WeirKind::Rectangular => {
                            ui.horizontal(|ui| {
                                let cap = ui.label("crest length L (m)");
                                ui.add(egui::DragValue::new(&mut s.crest_length_m).speed(0.05))
                                    .labelled_by(cap.id);
                            });
                        }
                        WeirKind::VNotch => {
                            ui.horizontal(|ui| {
                                let cap = ui.label("vertex angle θ (°)");
                                ui.add(egui::DragValue::new(&mut s.vertex_angle_deg).speed(1.0))
                                    .labelled_by(cap.id);
                            });
                        }
                    }
                    ui.horizontal(|ui| {
                        let cap = ui.label("discharge Cd");
                        ui.add(egui::DragValue::new(&mut s.discharge_coefficient).speed(0.005))
                            .labelled_by(cap.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Flow").strong());
                    ui.horizontal(|ui| {
                        let cap = ui.label("upstream head H (m)");
                        ui.add(egui::DragValue::new(&mut s.head_m).speed(0.01))
                            .labelled_by(cap.id);
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_weir(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative weir plate (an upright plate with a rectangular or V notch cut into its crest) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Discharge").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_weir_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.weir` borrow is
    // released here): build the weir's 3-D solid and load it.
    if app.weir.show_3d_request {
        app.weir.show_3d_request = false;
        load_weir_3d(app);
    }
}

/// Validate the form, evaluate the weir and format the readout.
fn run_weir(s: &mut WeirWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The weir's discharge `Q` (m³·s⁻¹) at the form's head, plus the head the
/// inverse rating curve returns for that `Q` — the quantities both the
/// readout and the 3-D gate need. Extracted so it is unit-testable and
/// shared; maps any domain error to a display string.
fn discharge_and_back(s: &WeirWorkbenchState) -> Result<(f64, f64), String> {
    match s.kind {
        WeirKind::Rectangular => {
            let weir = RectangularWeir::new(s.crest_length_m, s.discharge_coefficient)
                .map_err(|e| e.to_string())?;
            let q = weir.discharge(s.head_m).map_err(|e| e.to_string())?;
            let back = weir.head_for_discharge(q).map_err(|e| e.to_string())?;
            Ok((q, back))
        }
        WeirKind::VNotch => {
            let weir = VNotchWeir::new(s.vertex_angle_deg.to_radians(), s.discharge_coefficient)
                .map_err(|e| e.to_string())?;
            let q = weir.discharge(s.head_m).map_err(|e| e.to_string())?;
            let back = weir.head_for_discharge(q).map_err(|e| e.to_string())?;
            Ok((q, back))
        }
    }
}

/// Evaluate the weir and format the full readout. Extracted so it is
/// unit-testable.
fn compute(s: &WeirWorkbenchState) -> Result<String, String> {
    let (q, back) = discharge_and_back(s)?;
    let q_lps = q * 1000.0;

    let geometry = match s.kind {
        WeirKind::Rectangular => format!("rectangular  L = {l:.3} m", l = s.crest_length_m),
        WeirKind::VNotch => format!("V-notch      θ = {a:.1}°", a = s.vertex_angle_deg),
    };

    Ok(format!(
        "weir            : {geometry}\n\
         discharge Cd    : {cd:.3}\n\
         upstream head H : {head:.4} m\n\n\
         discharge Q     : {q:.6} m³/s  ({q_lps:.2} L/s)\n\
         inverse H(Q)    : {back:.4} m  (round-trips H)",
        cd = s.discharge_coefficient,
        head = s.head_m,
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

/// Build the weir as a triangle [`Mesh`] — an upright plate spanning the
/// channel with a notch cut into its crest: a centred rectangular slot for
/// a rectangular weir, or a centred triangular slot for a V-notch. The
/// plate is assembled from solid boxes around the opening (left jamb, right
/// jamb, and a sill below the notch), plus a thin base. Representative
/// geometry (the plate is drawn at a fixed visual size; the discharge
/// number is the `valenx-weir` result). `None` for an invalid
/// configuration.
fn weir_plate_mesh(s: &WeirWorkbenchState) -> Option<Mesh> {
    // Gate on the same validation the readout uses.
    discharge_and_back(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Fixed visual plate: half-width 1.0 in `y`, height 0..1.0 in `z`,
    // thin in the through-flow `x` direction.
    let half_w = 1.0_f64;
    let top = 1.0_f64;
    let tx = 0.06_f64;

    // The notch half-width at the crest. The V-notch opens to a point at
    // its vertex; the rectangular notch keeps full width down its slot.
    let notch_half_w = match s.kind {
        WeirKind::Rectangular => 0.35 * half_w,
        WeirKind::VNotch => {
            // Width of the V at the crest from the form's full vertex angle:
            // w = depth · tan(theta/2), capped to the plate.
            let depth = 0.6 * top;
            let half = depth * (s.vertex_angle_deg.to_radians() / 2.0).tan();
            half.clamp(0.05 * half_w, 0.9 * half_w)
        }
    };
    // Notch sill (the lowest point of the opening) measured from `z = 0`.
    let sill = match s.kind {
        WeirKind::Rectangular => 0.55 * top,
        WeirKind::VNotch => 0.4 * top,
    };

    // Left and right jambs flank the opening for the full plate height.
    let jamb_w = (half_w - notch_half_w) / 2.0;
    let jamb_cy = notch_half_w + jamb_w;
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, -jamb_cy, top / 2.0),
        Vector3::new(tx, jamb_w, top / 2.0),
    );
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, jamb_cy, top / 2.0),
        Vector3::new(tx, jamb_w, top / 2.0),
    );

    // Sill below the notch (a full-width band from the base up to `sill`).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, sill / 2.0),
        Vector3::new(tx, notch_half_w, sill / 2.0),
    );

    // For the V-notch, fill the two triangular shoulders so the opening
    // tapers to a point at the sill: stack thin full-thickness boxes whose
    // half-width shrinks linearly from the crest down to the vertex.
    if matches!(s.kind, WeirKind::VNotch) {
        let layers = 12usize;
        for i in 0..layers {
            let f0 = i as f64 / layers as f64;
            let f1 = (i + 1) as f64 / layers as f64;
            // Opening half-width at the top of this layer (linear taper:
            // full `notch_half_w` at the crest, 0 at the sill).
            let open_top = notch_half_w * (1.0 - f0);
            // Solid shoulder occupies `[open_top, notch_half_w]` on each
            // side; skip the top layer where it is degenerate.
            let shoulder_w = (notch_half_w - open_top) / 2.0;
            if shoulder_w <= 1e-6 {
                continue;
            }
            let shoulder_cy = open_top + shoulder_w;
            let z0 = sill + (top - sill) * f0;
            let z1 = sill + (top - sill) * f1;
            let hz = (z1 - z0) / 2.0;
            let cz = (z0 + z1) / 2.0;
            push_box(
                &mut nodes,
                &mut tris,
                Vector3::new(0.0, -shoulder_cy, cz),
                Vector3::new(tx, shoulder_w, hz),
            );
            push_box(
                &mut nodes,
                &mut tris,
                Vector3::new(0.0, shoulder_cy, cz),
                Vector3::new(tx, shoulder_w, hz),
            );
        }
    }

    // Thin base slab the plate stands on.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, -0.05),
        Vector3::new(0.3, half_w, 0.05),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-weir");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D weir plate and load it into the central viewport.
fn load_weir_3d(app: &mut ValenxApp) {
    let Some(mesh) = weir_plate_mesh(&app.weir) else {
        app.weir.error = Some("weir parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<weir>/valenx-weir"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"weir"}`** product: the canonical
/// sharp-crested weir plate built as a 3-D solid, paired with the workbench's
/// own `compute()` readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`WeirWorkbenchState::default`].
pub(crate) fn weir_product() -> crate::WorkspaceProduct {
    let s = WeirWorkbenchState::default();
    let mesh = weir_plate_mesh(&s).expect("canonical weir ⇒ plate solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<weir>/valenx-weir");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical weir ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Weir flow (sharp-crested)".into(),
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
    use valenx_weir::G_STANDARD;

    #[test]
    fn default_state_is_idle() {
        let s = WeirWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_discharge_and_round_trip() {
        let mut s = WeirWorkbenchState::default();
        run_weir(&mut s);
        assert!(
            s.error.is_none(),
            "default rectangular weir should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("discharge Q"));
        assert!(s.result.contains("inverse H(Q)"));
        // The crate's canonical magnitude for L=2, Cd=0.62, H=0.30:
        // Q ~ 0.601572 m^3/s; and H(Q) round-trips to 0.3000 m.
        assert!(s.result.contains("0.601572"), "got: {}", s.result);
        assert!(s.result.contains("0.3000"), "got: {}", s.result);
    }

    #[test]
    fn analyze_rejects_non_positive_head() {
        let mut s = WeirWorkbenchState {
            head_m: 0.0,
            ..Default::default()
        };
        run_weir(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn discharge_matches_weir_equation_ground_truth() {
        // Ground truth: rectangular Q = Cd·(2/3)·√(2g)·L·H^(3/2),
        // hand-computed for L=2, Cd=0.62, H=0.30, g = G_STANDARD.
        let cd = 0.62_f64;
        let l = 2.0_f64;
        let head = 0.30_f64;
        let q_expected = cd * (2.0 / 3.0) * (2.0 * G_STANDARD).sqrt() * l * head.powf(1.5);

        let s = WeirWorkbenchState::default();
        let (q, back) = discharge_and_back(&s).expect("valid weir");
        assert!(
            (q - q_expected).abs() < 1e-12,
            "Q = {q}, hand calc {q_expected}"
        );
        // The independently known magnitude, and the exact inverse.
        assert!((q - 0.601_572_0).abs() < 1e-6, "Q = {q}");
        assert!((back - head).abs() < 1e-9, "H(Q) = {back}, want {head}");
    }

    #[test]
    fn plate_mesh_for_default_is_nonempty_and_in_range() {
        let s = WeirWorkbenchState::default();
        let mesh = weir_plate_mesh(&s).expect("default rectangular weir yields a solid");
        assert!(mesh.nodes.len() > 8, "expected jambs + sill + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
        // The V-notch variant adds tapered shoulder layers -> more nodes.
        let vee = WeirWorkbenchState {
            kind: WeirKind::VNotch,
            ..Default::default()
        };
        let vee_mesh = weir_plate_mesh(&vee).expect("V-notch yields a solid");
        assert!(vee_mesh.nodes.len() > mesh.nodes.len());
    }

    #[test]
    fn plate_mesh_none_for_invalid() {
        let s = WeirWorkbenchState {
            head_m: -1.0,
            ..Default::default()
        };
        assert!(weir_plate_mesh(&s).is_none());
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
            draw_weir_workbench(app, ctx);
        });
    }

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_weir_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_weir_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_weir_workbench = true;
        run_weir(&mut app.weir);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton whose own Name egui clears;
        // each must be `labelled_by` its caption so an AI / screen reader can
        // find it by the caption text. The geometry row is weir-type-gated;
        // the always-drawn "discharge Cd" + "upstream head H (m)" are asserted.
        let mut app = ValenxApp::default();
        app.show_weir_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 2,
            "expected the numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by its caption (AI-drivable name)"
        );
        for caption in ["discharge Cd", "upstream head H (m)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
