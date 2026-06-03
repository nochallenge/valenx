//! Section cuts.
//!
//! [`cut`] slices a [`valenx_cad::Solid`] by an arbitrary plane and
//! returns the cross-section as line segments plus the
//! visible/hidden edges of what's still left of the solid.
//! [`hatch`] generates parallel-line hatch patterns inside the
//! section polygons, following the standard engineering convention
//! (45° at 2 mm spacing).
//!
//! Heavy lifting delegates to
//! [`valenx_mesh::cut::intersect_plane_triangles`] — we tessellate
//! the BRep first when needed, then hand it a triangle soup.

use nalgebra::Vector3;

use crate::error::TechDrawError;
use crate::hlr;

/// Result of a [`cut`] call.
#[derive(Clone, Debug)]
pub struct SectionResult {
    /// Cross-section outline line segments in **world space**. Each
    /// segment is two `(x, y, z)` endpoints lying on the cutting
    /// plane. Use [`SectionResult::project_to_2d`] to flatten via the
    /// view's camera matrix.
    pub cross_section_segments_world: Vec<[Vector3<f64>; 2]>,

    /// Visible edges of the remaining solid (in drawing-plane 2D),
    /// for the same camera/view.
    pub visible_edges: Vec<[(f64, f64); 2]>,
    /// Hidden edges of the remaining solid (in drawing-plane 2D).
    pub hidden_edges: Vec<[(f64, f64); 2]>,
}

impl SectionResult {
    /// Project the cross-section segments through a view camera
    /// matrix (same matrix used by the visible / hidden edge fields).
    pub fn project_to_2d(&self, camera: &nalgebra::Matrix4<f64>) -> Vec<[(f64, f64); 2]> {
        self.cross_section_segments_world
            .iter()
            .map(|seg| {
                [
                    {
                        let p = crate::projection::project_point(seg[0], camera);
                        (p[0], p[1])
                    },
                    {
                        let p = crate::projection::project_point(seg[1], camera);
                        (p[0], p[1])
                    },
                ]
            })
            .collect()
    }
}

/// Intersect `solid` with the plane defined by `plane_origin` and
/// `plane_normal`, classify what's left through `camera`.
///
/// The cross-section segments live in **world space** (the
/// drawing/section code projects them at render time).
pub fn cut(
    solid: &valenx_cad::Solid,
    plane_origin: Vector3<f64>,
    plane_normal: Vector3<f64>,
    camera: &nalgebra::Matrix4<f64>,
) -> Result<SectionResult, TechDrawError> {
    // Pull triangles from the solid (tessellate BRep if needed).
    let mesh = match solid {
        valenx_cad::Solid::Brep(_) => valenx_cad::tessellate::solid_to_mesh(solid, 0.1)
            .map_err(|e| TechDrawError::ExportFailed(format!("section tessellation: {e}")))?,
        valenx_cad::Solid::Mesh(m) => m.clone(),
    };
    if mesh.nodes.is_empty() {
        return Err(TechDrawError::EmptySolid);
    }

    // Build a flat triangle soup for valenx_mesh::cut::intersect_plane_triangles.
    let mut soup: Vec<[[f64; 3]; 3]> = Vec::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::element::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            // Bounds-check the connectivity: a loaded/corrupt mesh can
            // carry a node index past `nodes.len()`. Index by `.get()`
            // and skip the triangle on an out-of-range value rather than
            // panicking (mirrors projection.rs / hlr.rs).
            let (Some(&a), Some(&b), Some(&c)) = (
                mesh.nodes.get(tri[0] as usize),
                mesh.nodes.get(tri[1] as usize),
                mesh.nodes.get(tri[2] as usize),
            ) else {
                continue;
            };
            let v = [a, b, c];
            soup.push([
                [v[0].x, v[0].y, v[0].z],
                [v[1].x, v[1].y, v[1].z],
                [v[2].x, v[2].y, v[2].z],
            ]);
        }
    }

    let segs = valenx_mesh::cut::intersect_plane_triangles(&soup, plane_origin, plane_normal);
    let cross_section_segments_world: Vec<[Vector3<f64>; 2]> =
        segs.into_iter().map(|s| [s.a, s.b]).collect();

    // Visible / hidden edges of the original solid, classified by the
    // existing HLR pass. Trimming the solid to one side of the plane
    // before HLR would be the "proper" thing to do (currently
    // valenx_mesh::cut::slice() does centroid-keep slicing) — Phase
    // 5.5 will revisit. For v1 we report the full-solid classification
    // so the section overlay is decoupled from the trimmed-solid
    // outline (which matches typical CAD workflows).
    let (visible_edges, hidden_edges) = hlr::classify_edges(solid, camera)?;
    Ok(SectionResult {
        cross_section_segments_world,
        visible_edges,
        hidden_edges,
    })
}

/// Generate a parallel-line hatch pattern that fills the rectangular
/// 2D bbox of `polygon_segments_2d`. `spacing` is the distance
/// between adjacent hatch lines in mm, `angle_rad` is measured from
/// the +X axis.
///
/// v1 simplification: we hatch the **bounding box** of the section
/// polygons, not the polygons themselves. True polygon-inset hatching
/// needs the segments to form a closed loop with consistent
/// orientation — which the mesh-cut output doesn't guarantee (it
/// emits unordered line segments). Phase 5.5 will add a proper
/// segment-to-loop stitcher.
///
/// Returns each hatch line as a `[(x0, y0), (x1, y1)]` segment in
/// drawing-plane millimeters.
pub fn hatch(
    polygon_segments_2d: &[[(f64, f64); 2]],
    spacing: f64,
    angle_rad: f64,
) -> Vec<[(f64, f64); 2]> {
    if spacing <= 0.0 || polygon_segments_2d.is_empty() {
        return Vec::new();
    }
    // Bounding box of all segment endpoints.
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for seg in polygon_segments_2d {
        for p in seg {
            if p.0 < min_x {
                min_x = p.0;
            }
            if p.1 < min_y {
                min_y = p.1;
            }
            if p.0 > max_x {
                max_x = p.0;
            }
            if p.1 > max_y {
                max_y = p.1;
            }
        }
    }
    if !min_x.is_finite() {
        return Vec::new();
    }
    // Build a rotated coordinate frame: u-axis = direction of the
    // hatch lines, v-axis = perpendicular. We sweep `t` along v from
    // bbox_v_min to bbox_v_max in `spacing` increments, emitting one
    // hatch segment per step.
    let dir = (angle_rad.cos(), angle_rad.sin());
    let nrm = (-dir.1, dir.0);
    // Project the four bbox corners onto u and v to find the (umin,
    // umax, vmin, vmax) extent in the rotated frame.
    let corners = [
        (min_x, min_y),
        (max_x, min_y),
        (max_x, max_y),
        (min_x, max_y),
    ];
    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;
    for c in &corners {
        let u = c.0 * dir.0 + c.1 * dir.1;
        let v = c.0 * nrm.0 + c.1 * nrm.1;
        if u < u_min {
            u_min = u;
        }
        if u > u_max {
            u_max = u;
        }
        if v < v_min {
            v_min = v;
        }
        if v > v_max {
            v_max = v;
        }
    }
    let mut out = Vec::new();
    let mut v = v_min;
    while v <= v_max {
        // Endpoints in rotated frame, then back-rotate to world.
        let a_u = u_min;
        let a_v = v;
        let b_u = u_max;
        let b_v = v;
        let ax = a_u * dir.0 + a_v * nrm.0;
        let ay = a_u * dir.1 + a_v * nrm.1;
        let bx = b_u * dir.0 + b_v * nrm.0;
        let by = b_u * dir.1 + b_v * nrm.1;
        out.push([(ax, ay), (bx, by)]);
        v += spacing;
    }
    out
}

/// Hatch a section polygon using a named pattern from
/// [`crate::hatch_lib`] (Phase 18G).
///
/// For each angle in the pattern, runs the existing bbox-fill
/// [`hatch`] engine and concatenates the lines. Dot patterns
/// additionally emit short 0.1 mm "tick" segments at every grid
/// intersection so the SVG / PDF / DXF exporters can render them as
/// dots without a new primitive type.
///
/// Returns an empty vec when the pattern name is unknown.
pub fn hatch_with_pattern(
    polygon_segments_2d: &[[(f64, f64); 2]],
    pattern_name: &str,
) -> Vec<[(f64, f64); 2]> {
    let Some(pat) = crate::hatch_lib::by_name(pattern_name) else {
        return Vec::new();
    };
    let mut out: Vec<[(f64, f64); 2]> = Vec::new();
    for angle in &pat.angles {
        out.extend(hatch(polygon_segments_2d, pat.spacing.max(0.25), *angle));
    }
    if let Some(d) = pat.dot_spacing {
        // Bounding box of the polygon segments.
        let mut min_x = f64::INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        for seg in polygon_segments_2d {
            for p in seg {
                if p.0 < min_x {
                    min_x = p.0;
                }
                if p.1 < min_y {
                    min_y = p.1;
                }
                if p.0 > max_x {
                    max_x = p.0;
                }
                if p.1 > max_y {
                    max_y = p.1;
                }
            }
        }
        if min_x.is_finite() && d > 0.0 {
            let mut y = min_y;
            while y <= max_y {
                let mut x = min_x;
                while x <= max_x {
                    // 0.2 mm horizontal tick as a "dot" approximation.
                    out.push([(x, y), (x + 0.2, y)]);
                    x += d;
                }
                y += d;
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::view::ViewKind;
    use valenx_cad::primitives::box_solid;

    /// Task 15 — a cube cut by a plane through its center produces a
    /// square section polygon + hatch lines.
    ///
    /// We cut a 2×2×2 cube (centered around (1,1,1)) by the z=1
    /// plane. The cross-section is a 2×2 square in world XY at z=1.
    /// The cut yields multiple line segments (mesh triangles aren't
    /// stitched), and the hatch fills the bbox.
    #[test]
    fn cube_section_through_center_yields_segments_and_hatch() {
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        let cam = ViewKind::Top.camera_matrix();
        let res = cut(
            &cube,
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 0.0, 1.0),
            &cam,
        )
        .unwrap();
        // The cut must produce at least one segment. A z-cut through
        // a tessellated cube hits at least two horizontal-face edges
        // per side of the cube's "belt" — at least 4 segments total.
        assert!(
            !res.cross_section_segments_world.is_empty(),
            "expected non-empty cross-section, got 0"
        );
        // Project to 2D and verify every point is in the [0, 2]
        // square in drawing coords (top-view of the cube).
        let segs_2d = res.project_to_2d(&cam);
        for seg in &segs_2d {
            for (x, y) in seg {
                assert!(*x >= -1e-6 && *x <= 2.0 + 1e-6);
                assert!(*y >= -1e-6 && *y <= 2.0 + 1e-6);
            }
        }
        // Hatch the 2D segments. spacing=2mm at 45° → at least one
        // line spanning the 2-mm-wide bbox.
        let hatch_lines = hatch(&segs_2d, 2.0, std::f64::consts::FRAC_PI_4);
        assert!(
            !hatch_lines.is_empty(),
            "hatch should emit at least one line"
        );
    }

    #[test]
    fn cut_out_of_range_connectivity_does_not_panic() {
        // R32 M1: sibling of the R31 projection.rs fix. cut() built its
        // triangle soup via `mesh.nodes[tri[N] as usize]` directly, so a
        // loaded/corrupt Tri3 mesh whose connectivity references a node
        // past `nodes.len()` panicked ("index out of bounds"). The bad
        // triangle must be skipped gracefully.
        use valenx_mesh::element::{ElementBlock, ElementType};
        let mut mesh = valenx_mesh::Mesh::new("corrupt");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        // First triangle valid; second references node 99 (only 3 exist).
        block.connectivity = vec![0, 1, 2, 0, 1, 99];
        mesh.element_blocks.push(block);
        let s = valenx_cad::Solid::from_mesh(mesh);
        let cam = ViewKind::Top.camera_matrix();
        // Must not panic. A cut through z=0.0 touches the valid triangle.
        let _ = cut(
            &s,
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            &cam,
        );
    }

    #[test]
    fn hatch_empty_input_returns_empty() {
        assert!(hatch(&[], 2.0, 0.0).is_empty());
    }

    #[test]
    fn hatch_zero_spacing_returns_empty() {
        let segs = vec![[(0.0, 0.0), (1.0, 1.0)]];
        assert!(hatch(&segs, 0.0, 0.0).is_empty());
    }

    #[test]
    fn hatch_at_45_degrees_fills_bbox() {
        let segs = vec![[(0.0, 0.0), (10.0, 0.0)], [(10.0, 0.0), (10.0, 10.0)]];
        let h = hatch(&segs, 2.0, std::f64::consts::FRAC_PI_4);
        assert!(!h.is_empty());
    }

    /// Phase 18G — pattern hatch emits more segments than a single
    /// hatch for a crossed pattern.
    #[test]
    fn hatch_with_pattern_steel_doubles_single_pattern() {
        let segs = vec![[(0.0, 0.0), (10.0, 0.0)], [(10.0, 0.0), (10.0, 10.0)]];
        let h_single = hatch(&segs, 2.5, std::f64::consts::FRAC_PI_4);
        let h_pattern = hatch_with_pattern(&segs, "ANSI32");
        assert!(
            h_pattern.len() > h_single.len(),
            "crossed pattern should add the perpendicular set"
        );
    }

    #[test]
    fn hatch_with_pattern_concrete_includes_dot_ticks() {
        let segs = vec![[(0.0, 0.0), (10.0, 0.0)], [(10.0, 0.0), (10.0, 10.0)]];
        let h = hatch_with_pattern(&segs, "AR-CONC");
        // Dots show up as tiny 0.2-mm-long horizontal ticks.
        let n_dots = h
            .iter()
            .filter(|seg| {
                (seg[0].1 - seg[1].1).abs() < 1e-9 && (seg[1].0 - seg[0].0 - 0.2).abs() < 1e-9
            })
            .count();
        assert!(n_dots > 0, "concrete pattern should emit dot ticks");
    }

    #[test]
    fn hatch_with_pattern_unknown_returns_empty() {
        let segs = vec![[(0.0, 0.0), (10.0, 0.0)]];
        assert!(hatch_with_pattern(&segs, "NOT_A_PATTERN").is_empty());
    }
}
