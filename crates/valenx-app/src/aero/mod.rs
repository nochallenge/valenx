//! The Aerodynamics / Wind Tunnel workbench — implementation modules.
//!
//! The workbench wraps the native [`valenx_aero`] 3-D external-
//! aerodynamics CFD engine in a polished egui side panel. It is
//! structured the same way as the Genetics workbench
//! ([`crate::genetics`]): a thin `*_workbench.rs` host module mounts the
//! panel; the real work lives here, split by concern:
//!
//! - [`model`] — the form state, the `valenx-aero` request builder, the
//!   unit conversions and result formatters (pure, fully `#[test]`-ed).
//! - [`compute`] — background-thread orchestration of a steady solve or
//!   an angle-of-attack sweep, so the egui thread never blocks.
//! - [`viz`] — turning a converged result into coloured geometry for
//!   the app's 3-D viewport field-overlay path.
//! - [`panels`] — the eight workflow sections the side panel renders.
//!
//! The body-extraction helper that turns whatever geometry the app has
//! loaded into a [`valenx_aero::TriMesh`] lives in this module.

pub mod compute;
pub mod model;
pub mod panels;
pub mod viz;

use valenx_aero::{geometry::box_body, TriMesh};

use model::{BodySource, WindTunnelForm};

/// A body extracted and ready to drop into the virtual wind tunnel.
#[derive(Debug)]
pub struct ExtractedBody {
    /// The triangle-soup geometry.
    pub mesh: TriMesh,
    /// A short human-readable description of where it came from.
    pub source_label: String,
}

/// Extract the test body the form selects from the app's current state.
///
/// - [`BodySource::CurrentCadModel`] — tessellates the app's loaded CAD
///   solid (operand A).
/// - [`BodySource::ImportedStl`] — converts an already-loaded STL
///   triangle soup.
/// - [`BodySource::DemoBox`] — builds the parametric demo box; always
///   available.
///
/// Returns an error string when the chosen source has nothing loaded or
/// the geometry cannot be voxelized.
pub fn extract_body(
    form: &WindTunnelForm,
    current_solid: Option<&valenx_cad::Solid>,
    stl: Option<&valenx_viz::TriangleMesh>,
) -> Result<ExtractedBody, String> {
    match form.body_source {
        BodySource::CurrentCadModel => {
            let solid = current_solid.ok_or_else(|| {
                "no CAD model is loaded — create one in the Mesh Toolbox's Part \
                 section, import an STL, or pick the demo box"
                    .to_string()
            })?;
            // Tessellate the BRep at a chord tolerance scaled to the
            // solid — fine enough for a voxelized aero body.
            let mesh = TriMesh::from_solid(solid, 0.2)
                .map_err(|e| format!("tessellating the CAD model failed: {e}"))?;
            check_nonempty(&mesh)?;
            Ok(ExtractedBody {
                mesh,
                source_label: "current CAD model".to_string(),
            })
        }
        BodySource::ImportedStl => {
            let stl = stl.ok_or_else(|| {
                "no STL is loaded — use the Import STL button, or pick the demo box".to_string()
            })?;
            let mesh = trimesh_from_stl(stl)?;
            check_nonempty(&mesh)?;
            Ok(ExtractedBody {
                mesh,
                source_label: "imported STL".to_string(),
            })
        }
        BodySource::DemoBox => {
            let [sx, sy, sz] = form.demo_box_size;
            let mesh = box_body(
                nalgebra::Vector3::new(0.0, 0.0, 0.0),
                nalgebra::Vector3::new(sx.max(0.05), sy.max(0.05), sz.max(0.05)),
            );
            Ok(ExtractedBody {
                mesh,
                source_label: format!("demo box {sx:.1}×{sy:.1}×{sz:.1} m"),
            })
        }
    }
}

/// Convert an STL triangle soup into a [`valenx_aero::TriMesh`].
///
/// `valenx_viz::TriangleMesh` carries `f32` triangle vertices; the aero
/// solver works in `f64`, so the vertices are widened here.
fn trimesh_from_stl(stl: &valenx_viz::TriangleMesh) -> Result<TriMesh, String> {
    use valenx_aero::geometry::Triangle;
    let v = |p: [f32; 3]| nalgebra::Vector3::new(p[0] as f64, p[1] as f64, p[2] as f64);
    let tris: Vec<Triangle> = stl
        .triangles
        .iter()
        .map(|t| Triangle::new(v(t.vertices[0]), v(t.vertices[1]), v(t.vertices[2])))
        .filter(|t| !t.is_degenerate())
        .collect();
    if tris.is_empty() {
        return Err("the STL has no non-degenerate triangles".to_string());
    }
    Ok(TriMesh::from_triangles(tris))
}

/// Confirm an extracted body has a non-trivial, finite bounding box.
fn check_nonempty(mesh: &TriMesh) -> Result<(), String> {
    if mesh.is_empty() {
        return Err("the extracted body has no triangles".to_string());
    }
    match mesh.aabb() {
        Some(bb) => {
            let e = bb.extent();
            if !(e.x.is_finite() && e.y.is_finite() && e.z.is_finite())
                || e.x <= 0.0
                || e.y <= 0.0
                || e.z <= 0.0
            {
                Err("the body's bounding box is degenerate (zero or non-finite extent)".to_string())
            } else {
                Ok(())
            }
        }
        None => Err("the body has an empty bounding box".to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn demo_box_extracts_without_any_loaded_geometry() {
        // The demo box must always be available — no CAD model, no STL.
        let form = WindTunnelForm {
            body_source: BodySource::DemoBox,
            demo_box_size: [4.0, 2.0, 1.5],
            ..WindTunnelForm::default()
        };
        let body = extract_body(&form, None, None).expect("demo box");
        assert!(!body.mesh.is_empty());
        // A box is 12 triangles.
        assert_eq!(body.mesh.len(), 12);
        assert!(body.source_label.contains("demo box"));
        // The bounding box matches the requested size.
        let bb = body.mesh.aabb().expect("aabb");
        let e = bb.extent();
        assert!((e.x - 4.0).abs() < 1e-9);
        assert!((e.y - 2.0).abs() < 1e-9);
        assert!((e.z - 1.5).abs() < 1e-9);
    }

    #[test]
    fn demo_box_clamps_a_zero_dimension() {
        // A zero-size dimension would give a degenerate body; the
        // extractor clamps it to a small positive value so the case
        // is still well-posed.
        let form = WindTunnelForm {
            body_source: BodySource::DemoBox,
            demo_box_size: [0.0, 2.0, 1.5],
            ..WindTunnelForm::default()
        };
        let body = extract_body(&form, None, None).expect("clamped box");
        let bb = body.mesh.aabb().expect("aabb");
        assert!(bb.extent().x > 0.0);
    }

    #[test]
    fn current_cad_model_errors_when_none_is_loaded() {
        let form = WindTunnelForm {
            body_source: BodySource::CurrentCadModel,
            ..WindTunnelForm::default()
        };
        let err = extract_body(&form, None, None).unwrap_err();
        assert!(err.contains("no CAD model"));
    }

    #[test]
    fn imported_stl_errors_when_none_is_loaded() {
        let form = WindTunnelForm {
            body_source: BodySource::ImportedStl,
            ..WindTunnelForm::default()
        };
        let err = extract_body(&form, None, None).unwrap_err();
        assert!(err.contains("no STL"));
    }

    #[test]
    fn trimesh_from_stl_converts_a_triangle_soup() {
        use valenx_viz::stl::{StlTriangle, TriangleMesh};
        // A two-triangle STL square.
        let mesh = TriangleMesh {
            format: None,
            name: Some("square".to_string()),
            triangles: vec![
                StlTriangle {
                    normal: [0.0, 0.0, 1.0],
                    vertices: [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
                },
                StlTriangle {
                    normal: [0.0, 0.0, 1.0],
                    vertices: [[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
                },
            ],
        };
        let tm = trimesh_from_stl(&mesh).expect("converted");
        assert_eq!(tm.len(), 2);
        // Total area of the unit square is 1.
        assert!((tm.surface_area() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn trimesh_from_stl_rejects_an_empty_soup() {
        use valenx_viz::stl::TriangleMesh;
        let empty = TriangleMesh {
            format: None,
            name: None,
            triangles: Vec::new(),
        };
        assert!(trimesh_from_stl(&empty).is_err());
    }

    #[test]
    fn imported_stl_extracts_a_loaded_soup() {
        use valenx_viz::stl::{StlTriangle, TriangleMesh};
        // A closed-ish tetrahedron STL — four triangles.
        let stl = TriangleMesh {
            format: None,
            name: Some("tet".to_string()),
            triangles: vec![
                StlTriangle {
                    normal: [0.0, 0.0, 0.0],
                    vertices: [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                },
                StlTriangle {
                    normal: [0.0, 0.0, 0.0],
                    vertices: [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
                },
                StlTriangle {
                    normal: [0.0, 0.0, 0.0],
                    vertices: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                },
                StlTriangle {
                    normal: [0.0, 0.0, 0.0],
                    vertices: [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                },
            ],
        };
        let form = WindTunnelForm {
            body_source: BodySource::ImportedStl,
            ..WindTunnelForm::default()
        };
        let body = extract_body(&form, None, Some(&stl)).expect("stl body");
        assert_eq!(body.mesh.len(), 4);
        assert!(body.source_label.contains("STL"));
        // The body has a non-degenerate bounding box.
        assert!(check_nonempty(&body.mesh).is_ok());
    }

    #[test]
    fn check_nonempty_rejects_a_degenerate_body() {
        // A flat (zero-thickness) body has a degenerate bounding box.
        use valenx_aero::geometry::Triangle;
        let flat = TriMesh::from_triangles(vec![Triangle::new(
            nalgebra::Vector3::new(0.0, 0.0, 0.0),
            nalgebra::Vector3::new(1.0, 0.0, 0.0),
            nalgebra::Vector3::new(0.0, 1.0, 0.0),
        )]);
        // A single triangle has zero z-extent.
        assert!(check_nonempty(&flat).is_err());
    }
}
