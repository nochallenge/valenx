//! Expanded OpenSCAD primitive set.
//!
//! Phase 45's evaluator covered `cube` / `sphere` / `cylinder` /
//! `union` / `difference` / `intersection`.  Phase 52 adds:
//!
//! - **3-D**: `polyhedron` (explicit vertex/face mesh).
//! - **2-D**: `square`, `circle`, `polygon`, `text` — every 2-D
//!   primitive is emitted as a degenerate zero-height prism so it
//!   round-trips through the BRep boolean pipeline alongside the 3-D
//!   solids. v1 emits the polygon as a [`valenx_mesh::Mesh`] of fan
//!   triangles and wraps it as a mesh-backed [`valenx_cad::Solid`].
//! - **CSG**: `hull` (convex hull) — v1 is a pass-through union;
//!   true 3-D hulls land in Phase 58 (CGAL).
//! - **CSG**: `minkowski` — typed
//!   [`crate::OpenScadCsgError::NotImplemented`] stub. Minkowski sums
//!   require offset surfaces that the truck kernel doesn't expose;
//!   Phase 58 / 59 will provide the algorithm.
//!
//! All primitives accept their parameters in SI millimetres to match
//! the rest of the Valenx kernel.

use nalgebra::Vector3;

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::OpenScadCsgError;

/// Explicit-mesh `polyhedron(points, faces)`.
///
/// `points` is a list of 3-D vertex positions and `faces` is a list
/// of vertex-index loops (any face with >3 indices is fan-triangulated
/// from the first vertex).
pub fn polyhedron(
    points: &[[f64; 3]],
    faces: &[Vec<usize>],
) -> Result<Solid, OpenScadCsgError> {
    if points.len() < 4 {
        return Err(OpenScadCsgError::BadParameter {
            name: "points",
            reason: format!("need ≥4 vertices, got {}", points.len()),
        });
    }
    if faces.is_empty() {
        return Err(OpenScadCsgError::BadParameter {
            name: "faces",
            reason: "at least one face required".into(),
        });
    }
    let mut mesh = Mesh::new("polyhedron");
    mesh.nodes = points
        .iter()
        .map(|p| Vector3::new(p[0], p[1], p[2]))
        .collect();
    let mut block = ElementBlock {
        element_type: ElementType::Tri3,
        connectivity: Vec::new(),
    };
    for face in faces {
        if face.len() < 3 {
            return Err(OpenScadCsgError::BadParameter {
                name: "faces",
                reason: "face has <3 vertices".into(),
            });
        }
        // Fan-triangulate. v0 → vi → vi+1.
        for i in 1..face.len() - 1 {
            block.connectivity.extend_from_slice(&[
                face[0] as u32,
                face[i] as u32,
                face[i + 1] as u32,
            ]);
        }
    }
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

/// 2-D `square([sx, sy])`. Emitted as a degenerate zero-height prism
/// in the XY plane.
pub fn square(sx: f64, sy: f64) -> Result<Solid, OpenScadCsgError> {
    if sx <= 0.0 || sy <= 0.0 {
        return Err(OpenScadCsgError::BadParameter {
            name: "size",
            reason: format!("non-positive ({sx}, {sy})"),
        });
    }
    let mut mesh = Mesh::new("square");
    mesh.nodes = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(sx, 0.0, 0.0),
        Vector3::new(sx, sy, 0.0),
        Vector3::new(0.0, sy, 0.0),
    ];
    mesh.element_blocks.push(ElementBlock {
        element_type: ElementType::Tri3,
        connectivity: vec![0u32, 1, 2, 0, 2, 3],
    });
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

/// 2-D `circle(r)` — 32-sided regular polygon.
pub fn circle(r: f64) -> Result<Solid, OpenScadCsgError> {
    if r <= 0.0 {
        return Err(OpenScadCsgError::BadParameter {
            name: "r",
            reason: format!("non-positive radius {r}"),
        });
    }
    let n = 32u32;
    let mut mesh = Mesh::new("circle");
    // Centre vertex.
    mesh.nodes.push(Vector3::new(0.0, 0.0, 0.0));
    for k in 0..n {
        let theta = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
        mesh.nodes
            .push(Vector3::new(r * theta.cos(), r * theta.sin(), 0.0));
    }
    let mut conn = Vec::<u32>::with_capacity((n as usize) * 3);
    for k in 0..n {
        let a = k + 1;
        let b = if k + 1 < n { k + 2 } else { 1 };
        conn.extend_from_slice(&[0u32, a, b]);
    }
    mesh.element_blocks.push(ElementBlock {
        element_type: ElementType::Tri3,
        connectivity: conn,
    });
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

/// 2-D `polygon(points)` — explicit list of `(x, y)`. Fan-triangulated
/// from `points[0]`.
pub fn polygon(points: &[[f64; 2]]) -> Result<Solid, OpenScadCsgError> {
    if points.len() < 3 {
        return Err(OpenScadCsgError::BadParameter {
            name: "points",
            reason: format!("need ≥3 vertices, got {}", points.len()),
        });
    }
    let mut mesh = Mesh::new("polygon");
    mesh.nodes = points
        .iter()
        .map(|p| Vector3::new(p[0], p[1], 0.0))
        .collect();
    let mut conn = Vec::<u32>::with_capacity((points.len() - 2) * 3);
    for i in 1..points.len() - 1 {
        conn.extend_from_slice(&[0u32, i as u32, (i + 1) as u32]);
    }
    mesh.element_blocks.push(ElementBlock {
        element_type: ElementType::Tri3,
        connectivity: conn,
    });
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

/// 2-D `text(str, size)` — v1 emits a placeholder rectangle the size
/// of the bounding box that a 1-em glyph would occupy.  Real font
/// rendering ships in Phase 52.5 once `valenx-fonts` glyph outlines
/// are wired through.
pub fn text(s: &str, size: f64) -> Result<Solid, OpenScadCsgError> {
    if size <= 0.0 {
        return Err(OpenScadCsgError::BadParameter {
            name: "size",
            reason: format!("non-positive size {size}"),
        });
    }
    // Width is a rough average — 0.55 × size per glyph (works for
    // most sans-serif faces).  Empty strings yield a 0-width box so
    // the caller still gets a valid (degenerate) solid.
    let w = (s.chars().count() as f64) * 0.55 * size;
    let h = size;
    square(w.max(1e-3), h)
}

/// `hull(children)` — v1 is a pass-through union.  True convex hull
/// lands in Phase 58 via `valenx-cgal-port::convex_hull::hull_3d`.
pub fn hull(children: &[Solid]) -> Result<Solid, OpenScadCsgError> {
    if children.is_empty() {
        return Err(OpenScadCsgError::BadParameter {
            name: "children",
            reason: "hull(): no children".into(),
        });
    }
    let mut acc = children[0].clone();
    for c in &children[1..] {
        acc = valenx_cad::union(&acc, c).map_err(|e| OpenScadCsgError::Cad {
            op: "hull",
            reason: e.to_string(),
        })?;
    }
    Ok(acc)
}

/// `minkowski(children)` — typed `NotImplemented`.
pub fn minkowski(_children: &[Solid]) -> Result<Solid, OpenScadCsgError> {
    Err(OpenScadCsgError::NotImplemented { op: "minkowski" })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn square_two_tris() {
        let s = square(2.0, 3.0).expect("ok");
        let m = match s {
            Solid::Mesh(m) => m,
            _ => panic!("expected mesh"),
        };
        assert_eq!(m.nodes.len(), 4);
        assert_eq!(m.element_blocks[0].connectivity.len(), 6);
    }

    #[test]
    fn circle_32_segments() {
        let s = circle(1.0).expect("ok");
        let m = match s {
            Solid::Mesh(m) => m,
            _ => panic!("expected mesh"),
        };
        assert_eq!(m.nodes.len(), 33); // centre + 32 rim
        assert_eq!(m.element_blocks[0].connectivity.len(), 32 * 3);
    }

    #[test]
    fn polygon_triangle() {
        let s = polygon(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]).expect("ok");
        let m = match s {
            Solid::Mesh(m) => m,
            _ => panic!("expected mesh"),
        };
        assert_eq!(m.element_blocks[0].connectivity, vec![0, 1, 2]);
    }

    #[test]
    fn polyhedron_tetrahedron() {
        let pts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let faces = vec![
            vec![0, 1, 2],
            vec![0, 1, 3],
            vec![0, 2, 3],
            vec![1, 2, 3],
        ];
        let s = polyhedron(&pts, &faces).expect("ok");
        if let Solid::Mesh(m) = s {
            assert_eq!(m.nodes.len(), 4);
            assert_eq!(m.element_blocks[0].connectivity.len(), 12);
        } else {
            panic!("expected mesh");
        }
    }

    #[test]
    fn text_width_proportional() {
        let s = text("hi", 10.0).expect("ok");
        // Just need to confirm it produces a valid Solid; layout
        // tweaks ship later.
        assert!(matches!(s, Solid::Mesh(_)));
    }

    #[test]
    fn minkowski_not_implemented() {
        let r = minkowski(&[]);
        assert!(matches!(
            r,
            Err(OpenScadCsgError::NotImplemented { op: "minkowski" })
        ));
    }

    #[test]
    fn polyhedron_rejects_short_face() {
        let pts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let faces = vec![vec![0, 1]]; // only 2 verts
        let r = polyhedron(&pts, &faces);
        assert!(r.is_err());
    }
}
