//! Concrete-reinforcement workbench — rebar cages on `valenx-reinforcement`.
//!
//! A right-side panel that parametrically generates a beam or column **rebar
//! cage**, tessellates it (`cage::to_mesh`), and pushes the result into the
//! central 3-D viewport (the same `ValenxApp::stl` slot the STL importer and
//! the molecule viewer fill). Orbit the viewport to inspect the cage. The cage
//! geometry is headless-testable; the viewport is the interactive view.

use eframe::egui;
use std::path::PathBuf;

use valenx_reinforcement::cage;
use valenx_viz::stl::{StlTriangle, TriangleMesh};

use crate::types::LoadedStl;
use crate::ValenxApp;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Section {
    Beam,
    Column,
}

/// Persistent state for the reinforcement workbench.
pub struct ReinforcementWorkbenchState {
    section: Section,
    width: f64,
    depth: f64,
    length: f64,
    n_bars: usize,
    hoop_spacing: f64,
    status: String,
}

impl Default for ReinforcementWorkbenchState {
    fn default() -> Self {
        Self {
            section: Section::Beam,
            width: 0.3,
            depth: 0.5,
            length: 4.0,
            n_bars: 4,
            hoop_spacing: 0.2,
            status: String::new(),
        }
    }
}

/// Generate the rebar-cage mesh for the current settings.
fn run_reinforcement(s: &ReinforcementWorkbenchState) -> Result<valenx_mesh::Mesh, String> {
    let hoop = s.hoop_spacing.max(0.02);
    let n = s.n_bars.max(2);
    let (w, d, l) = (s.width.max(0.05), s.depth.max(0.05), s.length.max(0.1));
    let cage = match s.section {
        Section::Beam => cage::generate_beam(w, d, l, n, hoop),
        Section::Column => cage::generate_column(w, d, l, n, hoop),
    }
    .map_err(|e| e.to_string())?;
    Ok(cage::to_mesh(&cage))
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

/// Draw the reinforcement workbench (a no-op unless toggled on via
/// View → Reinforcement).
pub fn draw_reinforcement_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_reinforcement_workbench {
        return;
    }
    let mut generate = false;
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_reinforcement_workbench",
        "Concrete Reinforcement",
        |app, ui| {
            ui.label(
                egui::RichText::new("rebar cages · valenx-reinforcement")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.reinforcement;
            ui.horizontal(|ui| {
                ui.selectable_value(&mut s.section, Section::Beam, "Beam");
                ui.selectable_value(&mut s.section, Section::Column, "Column");
            });
            egui::Grid::new("reinf_params")
                .num_columns(2)
                .show(ui, |ui| {
                    ui.label("width (m)");
                    ui.add(
                        egui::DragValue::new(&mut s.width)
                            .speed(0.01)
                            .range(0.05..=3.0),
                    );
                    ui.end_row();
                    ui.label("depth (m)");
                    ui.add(
                        egui::DragValue::new(&mut s.depth)
                            .speed(0.01)
                            .range(0.05..=3.0),
                    );
                    ui.end_row();
                    ui.label(if s.section == Section::Beam {
                        "length (m)"
                    } else {
                        "height (m)"
                    });
                    ui.add(
                        egui::DragValue::new(&mut s.length)
                            .speed(0.05)
                            .range(0.1..=20.0),
                    );
                    ui.end_row();
                    ui.label("bars");
                    ui.add(egui::DragValue::new(&mut s.n_bars).speed(0.2));
                    ui.end_row();
                    ui.label("hoop spacing (m)");
                    ui.add(
                        egui::DragValue::new(&mut s.hoop_spacing)
                            .speed(0.01)
                            .range(0.02..=1.0),
                    );
                    ui.end_row();
                });
            ui.separator();
            if ui.button("▶ Generate cage → 3D viewport").clicked() {
                generate = true;
            }
            if !s.status.is_empty() {
                ui.label(egui::RichText::new(&s.status).small().weak());
            }
            ui.label(
                egui::RichText::new("Cage renders in the central 3D viewport — orbit to inspect.")
                    .small()
                    .weak(),
            );
        },
    );
    if close {
        app.show_reinforcement_workbench = false;
    }

    if generate {
        match run_reinforcement(&app.reinforcement) {
            Ok(mesh) => {
                let soup = mesh_to_triangle_soup(&mesh);
                let n = soup.triangles.len();
                if n == 0 {
                    app.reinforcement.status = "cage produced no geometry".to_string();
                } else {
                    app.mesh = None;
                    app.stl = Some(LoadedStl {
                        path: PathBuf::from("<reinforcement>/cage"),
                        mesh: soup,
                    });
                    app.frame_current_stl();
                    app.reinforcement.status = format!("{n} triangles in the viewport");
                }
            }
            Err(e) => app.reinforcement.status = format!("error: {e}"),
        }
    }
}

/// The agent-bridge product for the concrete-reinforcement workbench
/// (`show_3d{kind="reinforcement"}`).
///
/// Generates the **default** canonical rebar cage — a 0.3 × 0.5 × 4.0 m beam
/// with 4 longitudinal bars per face (Ø16) and Ø10 U-stirrups at 0.2 m centres —
/// via [`valenx_reinforcement::cage::generate_beam`] and tessellates it to a
/// `Tri3` [`valenx_mesh::Mesh`] with `cage::to_mesh`. Pure and app-state-free:
/// it builds a fresh [`ReinforcementWorkbenchState::default`] so it yields
/// exactly the cage the workbench shows on first generate. The readout reports
/// the beam dimensions and triangle count.
pub(crate) fn reinforcement_product() -> crate::WorkspaceProduct {
    let state = ReinforcementWorkbenchState::default();
    let (mesh, lines) = match run_reinforcement(&state) {
        Ok(mesh) => {
            let tris = mesh.total_elements();
            let lines = vec![
                format!(
                    "RC beam {:.1} × {:.1} × {:.1} m · {} longitudinal bars/face",
                    state.width, state.depth, state.length, state.n_bars
                ),
                format!(
                    "Ø10 U-stirrups at {:.2} m centres · 30 mm cover",
                    state.hoop_spacing
                ),
                format!("rebar cage: {tris} triangles"),
            ];
            (mesh, lines)
        }
        Err(e) => {
            // Theoretically unreachable for the canonical beam; degrade to a
            // tiny placeholder triangle + a note rather than panicking.
            let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
            block.connectivity = vec![0, 1, 2];
            let mut placeholder = valenx_mesh::Mesh::new("valenx-rebar-cage");
            placeholder.nodes = vec![
                nalgebra::Vector3::new(0.0, 0.0, 0.0),
                nalgebra::Vector3::new(1.0, 0.0, 0.0),
                nalgebra::Vector3::new(0.0, 1.0, 0.0),
            ];
            placeholder.element_blocks.push(block);
            placeholder.recompute_stats();
            (
                placeholder,
                vec![
                    "parametric rebar cage".to_string(),
                    format!("cage generation unavailable — showing placeholder ({e})"),
                ],
            )
        }
    };
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<reinforcement>/cage");
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Rebar cage (RC beam)".into(),
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
    fn beam_cage_tessellates_to_triangles() {
        let s = ReinforcementWorkbenchState::default();
        let mesh = run_reinforcement(&s).expect("beam cage");
        assert!(!mesh.nodes.is_empty(), "cage has nodes");
        let soup = mesh_to_triangle_soup(&mesh);
        assert!(!soup.triangles.is_empty(), "cage tessellates to triangles");
    }

    #[test]
    fn column_cage_generates() {
        let s = ReinforcementWorkbenchState {
            section: Section::Column,
            ..Default::default()
        };
        let mesh = run_reinforcement(&s).expect("column cage");
        assert!(!mesh.nodes.is_empty());
    }
}
