//! Phase 169 — `Graphic3d_ClipPlane` — interactive clipping planes
//! (1-6 simultaneous planes per view).
//!
//! ## What OCCT does
//!
//! `V3d_View::AddClipPlane(Graphic3d_ClipPlane)` configures an
//! arbitrary half-space `n · p + d ≥ 0` that the OpenGl pipeline
//! evaluates per fragment via `gl_ClipDistance[i]`. Hardware supports
//! 6 simultaneous planes (`GL_MAX_CLIP_DISTANCES`). Each plane carries
//! its own capping style (visible cap drawn at the clip boundary or
//! transparent), highlight colour, and on/off toggle. Use case: cut-
//! away inspection of internal geometry without re-meshing.
//!
//! ## v1 status — real implementation (CPU clip-plane set + clipper)
//!
//! The GPU `gl_ClipDistance` fast path needs a live `wgpu::Device`
//! with the `CLIP_DISTANCES` feature, which the app layer owns. What
//! this module ships is the renderer-independent half that OCCT's
//! `Graphic3d_ClipPlane` *is*: a managed set of up to
//! [`MAX_CLIP_PLANES`] half-space planes, plus a **real geometric
//! clipper** ([`clip_mesh`]) that cuts a triangle mesh against the
//! active planes — the same Sutherland-Hodgman half-space clip the
//! GPU does per fragment, done on the CPU so the result is usable for
//! cutaway export, capping-cross-section generation, and software
//! rendering.
//!
//! [`v3d_viewer_clipping_plane`] installs a plane into a
//! [`ClipPlaneSet`]; the app passes that set to the renderer (which
//! uploads it as a `vec4 planes[6]` uniform) and/or to [`clip_mesh`].

use nalgebra::Vector3;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::OcctVizError;

/// Maximum number of simultaneously-active clipping planes — matches
/// `GL_MAX_CLIP_DISTANCES` on baseline hardware.
pub const MAX_CLIP_PLANES: usize = 6;

/// Plane defined as `n · p + d ≥ 0` keeps the points on the `+n` side;
/// the other half-space is clipped.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ClipPlane {
    /// Outward-pointing normal (need not be unit — [`clip_mesh`]
    /// normalises internally).
    pub normal: [f32; 3],
    /// Signed offset along `normal`.
    pub offset: f32,
}

impl ClipPlane {
    /// Signed distance from `p` to this plane. Positive ⇒ kept,
    /// negative ⇒ clipped. The normal is normalised so the value is a
    /// true Euclidean distance.
    pub fn signed_distance(&self, p: Vector3<f64>) -> f64 {
        let n = Vector3::new(
            self.normal[0] as f64,
            self.normal[1] as f64,
            self.normal[2] as f64,
        );
        let len = n.norm();
        if len < 1e-12 {
            return 0.0;
        }
        (n.dot(&p) + self.offset as f64) / len
    }
}

/// An ordered set of up to [`MAX_CLIP_PLANES`] clipping planes —
/// Valenx's analogue of OCCT's per-view `Graphic3d_SequenceOfHClipPlane`.
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct ClipPlaneSet {
    /// Slot `i` holds `Some(plane)` when plane `i` is active.
    planes: [Option<ClipPlane>; MAX_CLIP_PLANES],
}

impl ClipPlaneSet {
    /// Empty set — no planes installed.
    pub fn new() -> Self {
        Self::default()
    }

    /// Install (or replace) the plane at `index`.
    pub fn set(&mut self, index: usize, plane: ClipPlane) {
        if index < MAX_CLIP_PLANES {
            self.planes[index] = Some(plane);
        }
    }

    /// Remove the plane at `index`, if any.
    pub fn clear(&mut self, index: usize) {
        if index < MAX_CLIP_PLANES {
            self.planes[index] = None;
        }
    }

    /// The active planes, in slot order.
    pub fn active(&self) -> impl Iterator<Item = &ClipPlane> {
        self.planes.iter().flatten()
    }

    /// Number of active planes.
    pub fn len(&self) -> usize {
        self.planes.iter().filter(|p| p.is_some()).count()
    }

    /// True when no planes are installed.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Pack the active planes into a `vec4[6]` uniform layout
    /// (`xyz` = normalised normal, `w` = offset). Inactive slots are
    /// zero so a fragment shader's `dot(plane, vec4(world, 1))` is a
    /// no-op (>= 0 always passes). This is exactly what the renderer
    /// uploads when the GPU `CLIP_DISTANCES` path is wired.
    pub fn to_uniform(&self) -> [[f32; 4]; MAX_CLIP_PLANES] {
        let mut out = [[0.0f32; 4]; MAX_CLIP_PLANES];
        for (slot, plane) in self.planes.iter().enumerate() {
            if let Some(p) = plane {
                let len = (p.normal[0] * p.normal[0]
                    + p.normal[1] * p.normal[1]
                    + p.normal[2] * p.normal[2])
                    .sqrt()
                    .max(1e-12);
                out[slot] = [
                    p.normal[0] / len,
                    p.normal[1] / len,
                    p.normal[2] / len,
                    p.offset / len,
                ];
            }
        }
        out
    }
}

/// Install `plane` at the given `index` (0..=5) into `set`.
///
/// This is the [`ClipPlaneSet`]-mutating analogue of OCCT's
/// `V3d_View::AddClipPlane`.
///
/// # Errors
///
/// [`OcctVizError::BadInput`] if `index >= MAX_CLIP_PLANES`, the
/// normal is zero or non-finite, or the offset is non-finite.
pub fn v3d_viewer_clipping_plane(
    set: &mut ClipPlaneSet,
    index: usize,
    plane: ClipPlane,
) -> Result<(), OcctVizError> {
    if index >= MAX_CLIP_PLANES {
        return Err(OcctVizError::bad_input(
            "index",
            format!("must be < {MAX_CLIP_PLANES} (got {index})"),
        ));
    }
    for (i, v) in plane.normal.iter().enumerate() {
        if !v.is_finite() {
            return Err(OcctVizError::bad_input(
                "normal",
                format!("component {i} is non-finite"),
            ));
        }
    }
    let n2 = plane.normal[0] * plane.normal[0]
        + plane.normal[1] * plane.normal[1]
        + plane.normal[2] * plane.normal[2];
    if n2 < 1e-12 {
        return Err(OcctVizError::bad_input("normal", "must be non-zero"));
    }
    if !plane.offset.is_finite() {
        return Err(OcctVizError::bad_input("offset", "must be finite"));
    }
    set.set(index, plane);
    Ok(())
}

/// Clip `mesh`'s `Tri3` geometry against every active plane in `set`,
/// returning the kept (`n · p + d ≥ 0` for all planes) portion.
///
/// Each triangle is intersected with each plane's half-space via
/// Sutherland-Hodgman polygon clipping; the resulting convex polygon
/// is fan-triangulated. A triangle wholly inside survives unchanged;
/// one wholly outside is dropped; one straddling a plane is cut along
/// the plane. With no active planes the input mesh is returned as-is.
///
/// This is the CPU equivalent of the GPU's per-fragment clip — useful
/// for generating a cutaway mesh for export or a capping cross-section.
pub fn clip_mesh(mesh: &Mesh, set: &ClipPlaneSet) -> Mesh {
    if set.is_empty() {
        return mesh.clone();
    }
    let planes: Vec<&ClipPlane> = set.active().collect();

    let mut out = Mesh::new(format!("{}-clipped", mesh.id));
    let mut conn: Vec<u32> = Vec::new();

    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            // Start with the triangle as a 3-vertex polygon.
            let mut poly: Vec<Vector3<f64>> = vec![
                mesh.nodes[tri[0] as usize],
                mesh.nodes[tri[1] as usize],
                mesh.nodes[tri[2] as usize],
            ];
            // Clip against each plane in turn.
            for plane in &planes {
                poly = clip_polygon_to_plane(&poly, plane);
                if poly.len() < 3 {
                    break;
                }
            }
            if poly.len() < 3 {
                continue;
            }
            // Fan-triangulate the clipped convex polygon, appending
            // new vertices to the output mesh.
            let base = out.nodes.len() as u32;
            for p in &poly {
                out.nodes.push(*p);
            }
            for k in 1..poly.len() - 1 {
                conn.extend_from_slice(&[base, base + k as u32, base + k as u32 + 1]);
            }
        }
    }

    if !conn.is_empty() {
        out.element_blocks.push(ElementBlock {
            element_type: ElementType::Tri3,
            connectivity: conn,
        });
    }
    out.recompute_stats();
    out
}

/// Sutherland-Hodgman clip of a convex polygon against one plane's
/// keep half-space (`signed_distance ≥ 0`).
fn clip_polygon_to_plane(poly: &[Vector3<f64>], plane: &ClipPlane) -> Vec<Vector3<f64>> {
    if poly.is_empty() {
        return Vec::new();
    }
    let mut out: Vec<Vector3<f64>> = Vec::with_capacity(poly.len() + 1);
    let n = poly.len();
    for i in 0..n {
        let cur = poly[i];
        let prev = poly[(i + n - 1) % n];
        let dc = plane.signed_distance(cur);
        let dp = plane.signed_distance(prev);
        let cur_in = dc >= 0.0;
        let prev_in = dp >= 0.0;
        if cur_in {
            if !prev_in {
                // Entering — add the crossing point.
                let t = dp / (dp - dc);
                out.push(prev + (cur - prev) * t);
            }
            out.push(cur);
        } else if prev_in {
            // Leaving — add the crossing point only.
            let t = dp / (dp - dc);
            out.push(prev + (cur - prev) * t);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn valid_plane() -> ClipPlane {
        ClipPlane {
            normal: [0.0, 0.0, 1.0],
            offset: 0.0,
        }
    }

    /// A single triangle straddling the z=0 plane: two verts below,
    /// one above.
    fn straddling_mesh() -> Mesh {
        let mut m = Mesh::new("strad");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, -1.0),
            Vector3::new(2.0, 0.0, -1.0),
            Vector3::new(1.0, 0.0, 1.0),
        ];
        m.element_blocks.push(ElementBlock {
            element_type: ElementType::Tri3,
            connectivity: vec![0, 1, 2],
        });
        m
    }

    #[test]
    fn rejects_out_of_range_index() {
        let mut set = ClipPlaneSet::new();
        let err = v3d_viewer_clipping_plane(&mut set, 99, valid_plane()).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_zero_normal() {
        let mut set = ClipPlaneSet::new();
        let p = ClipPlane {
            normal: [0.0; 3],
            offset: 0.0,
        };
        let err = v3d_viewer_clipping_plane(&mut set, 0, p).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_nan_offset() {
        let mut set = ClipPlaneSet::new();
        let p = ClipPlane {
            normal: [0.0, 0.0, 1.0],
            offset: f32::NAN,
        };
        let err = v3d_viewer_clipping_plane(&mut set, 0, p).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn valid_input_installs_plane() {
        let mut set = ClipPlaneSet::new();
        v3d_viewer_clipping_plane(&mut set, 0, valid_plane()).unwrap();
        assert_eq!(set.len(), 1);
        v3d_viewer_clipping_plane(&mut set, 3, valid_plane()).unwrap();
        assert_eq!(set.len(), 2);
        set.clear(0);
        assert_eq!(set.len(), 1);
    }

    #[test]
    fn signed_distance_normalises() {
        // A non-unit normal must still report a true Euclidean distance.
        let p = ClipPlane {
            normal: [0.0, 0.0, 10.0],
            offset: 0.0,
        };
        let d = p.signed_distance(Vector3::new(0.0, 0.0, 3.0));
        assert!((d - 3.0).abs() < 1e-9, "d={d}");
    }

    #[test]
    fn empty_set_clip_is_identity() {
        let m = straddling_mesh();
        let set = ClipPlaneSet::new();
        let out = clip_mesh(&m, &set);
        assert_eq!(out.nodes.len(), 3);
    }

    #[test]
    fn clip_keeps_only_the_positive_side() {
        // Clip the straddling triangle by z >= 0. The kept region is
        // the part above z=0: a smaller triangle near the apex. Every
        // output vertex must have z >= -EPS.
        let m = straddling_mesh();
        let mut set = ClipPlaneSet::new();
        set.set(0, valid_plane());
        let out = clip_mesh(&m, &set);
        assert!(!out.nodes.is_empty(), "clip should keep the apex region");
        for n in &out.nodes {
            assert!(n.z >= -1e-9, "vertex below the clip plane: z={}", n.z);
        }
    }

    #[test]
    fn clip_drops_a_fully_outside_triangle() {
        // A triangle entirely below z=0, clipped by z >= 0, vanishes.
        let mut m = Mesh::new("below");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, -2.0),
            Vector3::new(1.0, 0.0, -2.0),
            Vector3::new(0.0, 1.0, -3.0),
        ];
        m.element_blocks.push(ElementBlock {
            element_type: ElementType::Tri3,
            connectivity: vec![0, 1, 2],
        });
        let mut set = ClipPlaneSet::new();
        set.set(0, valid_plane());
        let out = clip_mesh(&m, &set);
        assert_eq!(out.element_blocks.len(), 0, "fully-clipped mesh is empty");
    }

    #[test]
    fn uniform_layout_normalises_and_zeros_inactive() {
        let mut set = ClipPlaneSet::new();
        set.set(
            2,
            ClipPlane {
                normal: [0.0, 0.0, 4.0],
                offset: 8.0,
            },
        );
        let u = set.to_uniform();
        // Slot 0 inactive → all zeros.
        assert_eq!(u[0], [0.0, 0.0, 0.0, 0.0]);
        // Slot 2 active → unit normal, offset scaled by 1/len.
        assert!((u[2][2] - 1.0).abs() < 1e-6, "nz={}", u[2][2]);
        assert!((u[2][3] - 2.0).abs() < 1e-6, "offset={}", u[2][3]);
    }
}
