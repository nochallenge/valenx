//! Reverse-engineering workbench — point cloud → mesh on `valenx-reverse`.
//!
//! A right-side panel that samples a demo point cloud (a sphere or torus
//! surface), runs the k-nearest-neighbour **triangulation** reconstruction,
//! and pushes the reconstructed mesh into the central 3-D viewport. Stands in
//! for a scan-to-CAD workflow. The reconstruction is headless-testable.

use eframe::egui;
use std::f64::consts::PI;
use std::path::PathBuf;

use nalgebra::Vector3;

use valenx_reverse::{triangulate, PointCloud};
use valenx_viz::stl::{StlTriangle, TriangleMesh};

use crate::types::LoadedStl;
use crate::ValenxApp;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Shape {
    Sphere,
    Torus,
}

/// Persistent state for the reverse-engineering workbench.
pub struct ReverseWorkbenchState {
    shape: Shape,
    density: usize,
    k: usize,
    status: String,
}

impl Default for ReverseWorkbenchState {
    fn default() -> Self {
        Self {
            shape: Shape::Sphere,
            density: 24,
            k: 8,
            status: String::new(),
        }
    }
}

/// Parse a demo-shape name (for the agent `SetControl` bridge) into a [`Shape`].
/// Case-insensitive; matches the two picker words. Fail-loud on an unknown name.
fn parse_shape(s: &str) -> Result<Shape, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "sphere" => Ok(Shape::Sphere),
        "torus" => Ok(Shape::Torus),
        other => Err(format!(
            "unknown shape '{other}' (expected 'sphere' or 'torus')"
        )),
    }
}

impl ReverseWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). `shape` is the
    /// Sphere/Torus selection (set by option name).
    pub fn agent_control_names() -> &'static [&'static str] {
        &["shape", "cloud density", "k neighbours"]
    }

    /// Set one labelled control by caption, for the agent `SetControl` bridge.
    /// Numeric fields read [`AgentValue::as_i64`] (the two demo counts); `shape`
    /// reads [`AgentValue::as_str`] and matches the picker names. Fail-loud: an
    /// unknown caption, wrong type, or negative count returns `Err` — never a
    /// panic, no field written on error.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        match name {
            "shape" => self.shape = parse_shape(value.as_str()?)?,
            "cloud density" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("cloud density must be >= 0, got {n}"));
                }
                self.density = n as usize;
            }
            "k neighbours" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("k neighbours must be >= 0, got {n}"));
                }
                self.k = n as usize;
            }
            other => return Err(format!("unknown Reverse Engineering control: {other:?}")),
        }
        Ok(())
    }
}

/// Sample a demo point cloud on a sphere or torus surface.
fn demo_cloud(shape: Shape, density: usize) -> Vec<Vector3<f64>> {
    let d = density.clamp(6, 64);
    let mut pts = Vec::new();
    match shape {
        Shape::Sphere => {
            for i in 0..=d {
                let theta = PI * i as f64 / d as f64;
                for j in 0..d {
                    let phi = 2.0 * PI * j as f64 / d as f64;
                    pts.push(Vector3::new(
                        theta.sin() * phi.cos(),
                        theta.sin() * phi.sin(),
                        theta.cos(),
                    ));
                }
            }
        }
        Shape::Torus => {
            let (rr, tr) = (1.0, 0.35);
            for i in 0..d {
                let u = 2.0 * PI * i as f64 / d as f64;
                for j in 0..d {
                    let v = 2.0 * PI * j as f64 / d as f64;
                    pts.push(Vector3::new(
                        (rr + tr * v.cos()) * u.cos(),
                        (rr + tr * v.cos()) * u.sin(),
                        tr * v.sin(),
                    ));
                }
            }
        }
    }
    pts
}

/// Reconstruct a mesh from the demo cloud.
fn run_reverse(s: &ReverseWorkbenchState) -> Result<valenx_mesh::Mesh, String> {
    let cloud = PointCloud::from_points(demo_cloud(s.shape, s.density));
    triangulate(&cloud, s.k.clamp(4, 20)).map_err(|e| e.to_string())
}

/// Convert a triangle (`Tri3`) mesh into the viewport's triangle soup.
fn mesh_to_triangle_soup(mesh: &valenx_mesh::Mesh) -> TriangleMesh {
    let mut out = TriangleMesh::new();
    for block in &mesh.element_blocks {
        for tri in block.connectivity.chunks_exact(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            let mut t = StlTriangle {
                normal: [0.0, 0.0, 0.0],
                vertices: [
                    [a.x as f32, a.y as f32, a.z as f32],
                    [b.x as f32, b.y as f32, b.z as f32],
                    [c.x as f32, c.y as f32, c.z as f32],
                ],
            };
            t.normal = t.computed_normal();
            out.triangles.push(t);
        }
    }
    out
}

/// Draw the reverse-engineering workbench (a no-op unless toggled on via
/// View → Reverse Engineering).
pub fn draw_reverse_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_reverse_workbench {
        return;
    }
    let mut reconstruct = false;
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_reverse_workbench",
        "Reverse Engineering",
        |app, ui| {
            ui.label(
                egui::RichText::new("point cloud → mesh · valenx-reverse")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.reverse;
            ui.horizontal(|ui| {
                ui.selectable_value(&mut s.shape, Shape::Sphere, "Sphere");
                ui.selectable_value(&mut s.shape, Shape::Torus, "Torus");
            });
            egui::Grid::new("reverse_params")
                .num_columns(2)
                .show(ui, |ui| {
                    // Associate each DragValue with its caption via
                    // `labelled_by`, so the spin button is named for a screen
                    // reader / AI driver (egui clears a DragValue's own Name).
                    let density = ui.label("cloud density");
                    ui.add(
                        egui::DragValue::new(&mut s.density)
                            .speed(0.5)
                            .range(6..=64),
                    )
                    .labelled_by(density.id);
                    ui.end_row();
                    let k = ui.label("k neighbours");
                    ui.add(egui::DragValue::new(&mut s.k).speed(0.2).range(4..=20))
                        .labelled_by(k.id);
                    ui.end_row();
                });
            ui.separator();
            if ui.button("▶ Reconstruct → 3D viewport").clicked() {
                reconstruct = true;
            }
            if !s.status.is_empty() {
                ui.label(egui::RichText::new(&s.status).small().weak());
            }
            ui.label(
                egui::RichText::new("Reconstructed mesh renders in the central 3D viewport.")
                    .small()
                    .weak(),
            );
        },
    );
    if close {
        app.show_reverse_workbench = false;
    }

    if reconstruct {
        let n_pts = demo_cloud(app.reverse.shape, app.reverse.density).len();
        match run_reverse(&app.reverse) {
            Ok(mesh) => {
                let soup = mesh_to_triangle_soup(&mesh);
                let n = soup.triangles.len();
                if n == 0 {
                    app.reverse.status = "no triangles reconstructed".to_string();
                } else {
                    app.mesh = None;
                    app.stl = Some(LoadedStl::new(
                        PathBuf::from("<reverse>/reconstruction"),
                        soup,
                    ));
                    app.frame_current_stl();
                    app.reverse.status = format!("{n_pts} points → {n} triangles");
                }
            }
            Err(e) => app.reverse.status = format!("error: {e}"),
        }
    }
}

/// The agent-bridge product for the reverse-engineering workbench
/// (`show_3d{kind="reverse"}`).
///
/// Runs the **default** demo reconstruction — the canonical sphere point cloud
/// (density 24 ⇒ ~600 points) reconstructed by the k-NN
/// [`valenx_reverse::triangulate`] (k = 8) into a `Tri3`
/// [`valenx_mesh::Mesh`]. Pure and app-state-free: it builds a fresh
/// [`ReverseWorkbenchState::default`] so it yields exactly the surface the
/// workbench shows on first reconstruct. The readout reports the point and
/// triangle counts.
pub(crate) fn reverse_product() -> crate::WorkspaceProduct {
    let state = ReverseWorkbenchState::default();
    let n_pts = demo_cloud(state.shape, state.density).len();
    let (mesh, lines) = match run_reverse(&state) {
        Ok(mesh) => {
            let tris = mesh.total_elements();
            let lines = vec![
                "scan-to-CAD: sphere point cloud → k-NN surface".to_string(),
                format!("{n_pts} sampled points · k = {} neighbours", state.k),
                format!("reconstructed surface: {tris} triangles"),
            ];
            (mesh, lines)
        }
        Err(e) => {
            // Theoretically unreachable for the canonical cloud; degrade to a
            // tiny placeholder triangle + a note rather than panicking.
            let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
            block.connectivity = vec![0, 1, 2];
            let mut placeholder = valenx_mesh::Mesh::new("valenx-reverse-surface");
            placeholder.nodes = vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.0, 1.0, 0.0),
            ];
            placeholder.element_blocks.push(block);
            placeholder.recompute_stats();
            (
                placeholder,
                vec![
                    "point-cloud surface reconstruction".to_string(),
                    format!("reconstruction unavailable — showing placeholder ({e})"),
                ],
            )
        }
    };
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<reverse>/reconstruction");
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Reverse-engineered surface (k-NN)".into(),
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
    fn agent_set_sets_param_unknown_and_type_mismatch_err() {
        use crate::agent_commands::AgentValue;
        let mut s = ReverseWorkbenchState::default();
        // A representative integer set lands in state.
        s.agent_set("cloud density", &AgentValue::Int(32)).unwrap();
        assert_eq!(s.density, 32);
        // The shape enum is set by option name.
        s.agent_set("shape", &AgentValue::Str("torus".into()))
            .unwrap();
        assert_eq!(s.shape, Shape::Torus);
        // Unknown caption -> Err.
        assert!(s.agent_set("no such control", &AgentValue::Int(1)).is_err());
        // Type mismatch (string into the integer count) -> Err, field untouched.
        assert!(s
            .agent_set("cloud density", &AgentValue::Str("many".into()))
            .is_err());
        assert_eq!(s.density, 32, "rejected set leaves field untouched");
    }

    #[test]
    fn sphere_cloud_reconstructs_to_a_mesh() {
        let mesh = run_reverse(&ReverseWorkbenchState::default()).expect("triangulate");
        assert!(!mesh.nodes.is_empty(), "mesh has nodes");
        let soup = mesh_to_triangle_soup(&mesh);
        assert!(
            !soup.triangles.is_empty(),
            "reconstruction yields triangles"
        );
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_reverse_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_reverse_workbench(&mut app, ctx);
        });
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        use egui::accesskit::{Node, NodeId, Role};

        // Render with accesskit enabled and read the emitted a11y tree — the
        // same tree a screen reader / AI UI-Automation driver consumes. Every
        // DragValue (Role::SpinButton) must carry a caption via `labelled_by`,
        // since egui clears a DragValue's own Name.
        let mut app = ValenxApp::default();
        app.show_reverse_workbench = true;
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_reverse_workbench(&mut app, ctx);
        });
        let nodes: Vec<(NodeId, Node)> = out
            .platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes;

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // The reconstruction params draw cloud-density + k-neighbours.
        assert!(
            spin_buttons.len() >= 2,
            "expected the reverse numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every reverse DragValue must be labelled_by its caption"
        );
        for caption in ["cloud density", "k neighbours"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
