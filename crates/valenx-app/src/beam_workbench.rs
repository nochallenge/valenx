//! The right-side **Beam Workbench** panel — native Euler-Bernoulli beam
//! bending over `valenx-beam`.
//!
//! Mirrors the DC Motor / Truss workbenches: a resizable [`egui::SidePanel`]
//! gated on `crate::ValenxApp::show_beam_workbench`, toggled from the View
//! menu. The form drives a [`valenx_beam::Beam`] with a rectangular section
//! under a chosen support condition and load; "Analyze" reports the peak
//! deflection, bending moment and bending stress from the textbook closed
//! forms, and "Show 3-D beam" loads a schematic beam-plus-supports-plus-load
//! solid into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_beam::{Beam, Load, Section, Support};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Beam Workbench.
pub struct BeamWorkbenchState {
    /// Span `L` (m).
    length: f64,
    /// Young's modulus `E` (Pa).
    youngs_modulus: f64,
    /// Rectangular-section width `b` (m).
    width: f64,
    /// Rectangular-section height `h` (m).
    height: f64,
    /// Load magnitude (N for a point load, N/m for a UDL).
    load_mag: f64,
    /// Support condition: 0 = cantilever, 1 = simply-supported, 2 = fixed-fixed.
    support_idx: usize,
    /// Apply the load as a uniformly distributed load (else a point load).
    udl: bool,
    /// Formatted readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D beam solid (serviced after the panel
    /// draws).
    show_3d_request: bool,
}

impl Default for BeamWorkbenchState {
    fn default() -> Self {
        // A 3 m simply-supported steel beam, 50 x 100 mm, central 1 kN load.
        Self {
            length: 3.0,
            youngs_modulus: 2.0e11,
            width: 0.05,
            height: 0.1,
            load_mag: 1000.0,
            support_idx: 1,
            udl: false,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// The [`Support`] the form currently selects.
fn support_of(idx: usize) -> Support {
    match idx {
        0 => Support::Cantilever,
        2 => Support::FixedFixed,
        _ => Support::SimplySupported,
    }
}

/// Draw the Beam Workbench right-side panel. A no-op when the
/// `show_beam_workbench` toggle is off.
pub fn draw_beam_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_beam_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_beam_workbench",
        "Beam",
        |app, ui| {
            ui.label(
                egui::RichText::new("native Euler-Bernoulli beam bending · valenx-beam")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.beam;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Beam + section").strong());
                    // Associate each numeric `DragValue` with its caption via
                    // `labelled_by`, so the spin button carries the caption as
                    // its accessibility / UI-Automation Name (egui clears a
                    // DragValue's own Name, leaving it anonymous to a screen
                    // reader / AI driver otherwise). The b/h pair shares the one
                    // "section b x h (m)" caption.
                    ui.horizontal(|ui| {
                        let l = ui.label("span L (m)");
                        ui.add(egui::DragValue::new(&mut s.length).speed(0.05))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("Young's E (Pa)");
                        ui.add(egui::DragValue::new(&mut s.youngs_modulus).speed(1.0e9))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("section b x h (m)");
                        ui.add(egui::DragValue::new(&mut s.width).speed(0.002))
                            .labelled_by(l.id);
                        ui.add(egui::DragValue::new(&mut s.height).speed(0.002))
                            .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Support").strong());
                    ui.radio_value(&mut s.support_idx, 0, "cantilever (fixed-free)");
                    ui.radio_value(&mut s.support_idx, 1, "simply supported");
                    ui.radio_value(&mut s.support_idx, 2, "fixed-fixed");

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Load").strong());
                    ui.checkbox(&mut s.udl, "uniformly distributed (else point)");
                    ui.horizontal(|ui| {
                        let l = ui.label(if s.udl { "w (N/m)" } else { "P (N)" });
                        ui.add(egui::DragValue::new(&mut s.load_mag).speed(10.0))
                            .labelled_by(l.id);
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_beam(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D beam").strong())
                        .on_hover_text(
                            "Build a schematic beam with its supports and load as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Bending").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_beam_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.beam` borrow is released
    // here): build the beam's 3-D solid and load it.
    if app.beam.show_3d_request {
        app.beam.show_3d_request = false;
        load_beam_3d(app);
    }
}

/// Build a validated [`Beam`] from the form, mapping the domain error to a
/// display string.
fn build_beam(s: &BeamWorkbenchState) -> Result<Beam, String> {
    let section = Section::rectangular(s.width, s.height).map_err(|e| e.to_string())?;
    Beam::new(s.length, s.youngs_modulus, section).map_err(|e| e.to_string())
}

/// Validate the form, analyse the beam and format the readout.
fn run_beam(s: &mut BeamWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build + analyse the beam and format the readout. Extracted so it is
/// unit-testable.
fn compute(s: &BeamWorkbenchState) -> Result<String, String> {
    let beam = build_beam(s)?;
    let support = support_of(s.support_idx);
    let load = if s.udl {
        Load::Udl {
            intensity: s.load_mag,
        }
    } else {
        Load::Point { force: s.load_mag }
    };
    let res = beam.analyze(support, load).map_err(|e| e.to_string())?;
    let i = beam.section.second_moment_area();
    let support_name = match s.support_idx {
        0 => "cantilever",
        2 => "fixed-fixed",
        _ => "simply supported",
    };
    let load_name = if s.udl { "UDL w =" } else { "point P =" };
    Ok(format!(
        "span L  : {:.3} m\n\
         E       : {:.2e} Pa\n\
         section : {:.3} x {:.3} m (rect)\n\
         I       : {:.3e} m⁴\n\
         support : {}\n\
         load    : {} {:.1}\n\n\
         max deflection: {:.4} mm\n\
         max moment    : {:.1} N·m\n\
         max stress    : {:.2} MPa",
        s.length,
        s.youngs_modulus,
        s.width,
        s.height,
        i,
        support_name,
        load_name,
        s.load_mag,
        res.max_deflection * 1000.0,
        res.max_moment,
        res.max_stress / 1.0e6,
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

/// Build a schematic beam as a triangle [`Mesh`] — the prismatic beam (length
/// along x), its supports for the selected condition, and a load marker.
/// `None` for an invalid beam. The deflection is reported numerically; the
/// solid is the undeformed beam-and-loading diagram.
fn beam_solid_mesh(s: &BeamWorkbenchState) -> Option<Mesh> {
    build_beam(s).ok()?;
    let (l, w, h) = (s.length, s.width.max(1e-3), s.height.max(1e-3));
    let hl = l / 2.0;
    let blk = (w + h) * 0.5; // support / marker size scale
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Beam body.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::zeros(),
        Vector3::new(hl, w / 2.0, h / 2.0),
    );

    // Supports.
    match s.support_idx {
        0 => {
            // Cantilever: a wall at the fixed (-x) end.
            push_box(
                &mut nodes,
                &mut tris,
                Vector3::new(-hl - blk * 0.5, 0.0, 0.0),
                Vector3::new(blk * 0.5, w * 0.9, h * 1.2),
            );
        }
        2 => {
            // Fixed-fixed: walls at both ends.
            for sx in [-(hl + blk * 0.5), hl + blk * 0.5] {
                push_box(
                    &mut nodes,
                    &mut tris,
                    Vector3::new(sx, 0.0, 0.0),
                    Vector3::new(blk * 0.5, w * 0.9, h * 1.2),
                );
            }
        }
        _ => {
            // Simply supported: a block under each end.
            for sx in [-hl, hl] {
                push_box(
                    &mut nodes,
                    &mut tris,
                    Vector3::new(sx, 0.0, -h / 2.0 - blk * 0.4),
                    Vector3::new(blk * 0.3, w * 0.7, blk * 0.4),
                );
            }
        }
    }

    // Load marker, above the beam.
    if s.udl {
        push_box(
            &mut nodes,
            &mut tris,
            Vector3::new(0.0, 0.0, h / 2.0 + blk * 0.25),
            Vector3::new(hl, w * 0.5, blk * 0.2),
        );
    } else {
        let load_x = if s.support_idx == 0 { hl } else { 0.0 };
        push_box(
            &mut nodes,
            &mut tris,
            Vector3::new(load_x, 0.0, h / 2.0 + blk * 0.7),
            Vector3::new(blk * 0.15, blk * 0.15, blk * 0.7),
        );
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-beam");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D beam solid and load it into the central viewport.
fn load_beam_3d(app: &mut ValenxApp) {
    let Some(mesh) = beam_solid_mesh(&app.beam) else {
        app.beam.error = Some("beam parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<beam>/valenx-beam"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"beam"}`** product: the schematic
/// beam-plus-supports-plus-load solid built from the canonical 3 m
/// simply-supported steel beam (50×100 mm section, central 1 kN point load),
/// paired with the Euler-Bernoulli bending readout rows (max deflection /
/// moment / stress), at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`BeamWorkbenchState::default`].
pub(crate) fn beam_product() -> crate::WorkspaceProduct {
    let s = BeamWorkbenchState::default();
    let mesh = beam_solid_mesh(&s).expect("canonical beam ⇒ beam-and-supports solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<beam>/valenx-beam");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical beam ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Beam (Euler-Bernoulli bending)".into(),
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
        let s = BeamWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_deflection_moment_stress() {
        let mut s = BeamWorkbenchState::default();
        run_beam(&mut s);
        assert!(s.error.is_none(), "default beam analyzes: {:?}", s.error);
        assert!(s.result.contains("max deflection"));
        assert!(s.result.contains("max moment"));
        assert!(s.result.contains("max stress"));
    }

    #[test]
    fn analyze_rejects_nonpositive_section() {
        let mut s = BeamWorkbenchState {
            height: 0.0,
            ..Default::default()
        };
        run_beam(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn beam_mesh_for_default_is_nonempty_and_in_range() {
        let s = BeamWorkbenchState::default();
        let mesh = beam_solid_mesh(&s).expect("default beam yields a solid");
        assert!(mesh.nodes.len() > 8, "expected beam + supports + load");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn beam_mesh_none_for_invalid() {
        let s = BeamWorkbenchState {
            width: 0.0,
            ..Default::default()
        };
        assert!(beam_solid_mesh(&s).is_none());
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
            draw_beam_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_beam_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_beam_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_beam_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_beam_workbench = true;
        run_beam(&mut app.beam);
        app.beam.error = Some("invalid beam parameters".to_string());
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // The beam + section + load DragValues are SpinButtons; each must be
        // `labelled_by` its caption (egui clears a DragValue's own Name), so an
        // AI / screen reader can find the control by the caption text. The width
        // and height share the one "section b x h (m)" caption.
        let mut app = ValenxApp::default();
        app.show_beam_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // span L, E, width, height, load magnitude.
        assert!(
            spin_buttons.len() >= 5,
            "expected the beam numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every beam DragValue must be labelled_by its caption (AI-drivable name)"
        );

        // Default state has udl == false, so the load caption is "P (N)".
        for caption in ["span L (m)", "Young's E (Pa)", "section b x h (m)", "P (N)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
        // The Analyze button stays a named, invokable node.
        assert!(
            nodes.iter().any(|(_, n)| n.role() == Role::Button
                && n.name().is_some_and(|s| s.contains("Analyze"))),
            "the Analyze button is a named, invokable node"
        );
    }
}
