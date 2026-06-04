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

#[derive(Clone, Copy, PartialEq, Eq)]
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
        Self { shape: Shape::Sphere, density: 24, k: 8, status: String::new() }
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
    egui::SidePanel::right("valenx_reverse_workbench")
        .resizable(true)
        .default_width(330.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            ui.heading("Reverse Engineering");
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
            egui::Grid::new("reverse_params").num_columns(2).show(ui, |ui| {
                ui.label("cloud density");
                ui.add(egui::DragValue::new(&mut s.density).speed(0.5).range(6..=64));
                ui.end_row();
                ui.label("k neighbours");
                ui.add(egui::DragValue::new(&mut s.k).speed(0.2).range(4..=20));
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
        });

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
                    app.stl = Some(LoadedStl {
                        path: PathBuf::from("<reverse>/reconstruction"),
                        mesh: soup,
                    });
                    app.frame_current_stl();
                    app.reverse.status = format!("{n_pts} points → {n} triangles");
                }
            }
            Err(e) => app.reverse.status = format!("error: {e}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sphere_cloud_reconstructs_to_a_mesh() {
        let mesh = run_reverse(&ReverseWorkbenchState::default()).expect("triangulate");
        assert!(!mesh.nodes.is_empty(), "mesh has nodes");
        let soup = mesh_to_triangle_soup(&mesh);
        assert!(!soup.triangles.is_empty(), "reconstruction yields triangles");
    }
}
