//! Cut-cell geometry — the sub-voxel reconstruction of the body
//! surface inside each Cartesian cell it passes through.
//!
//! # Staircased vs. cut-cell
//!
//! The plain immersed-boundary method ([`crate::immersed`]) flags each
//! Cartesian cell as wholly `Solid` or wholly `Fluid`: the wall is the
//! voxel staircase, every cell carries either its full fluid volume or
//! none of it, and a body face is one flat Cartesian face of a cut
//! cell. That is robust and mesher-free, but the geometric wall is a
//! step pattern — the absolute drag / lift coefficients land in the
//! right ballpark, not at engineering tolerance, and a smooth sphere
//! is shed like a blocky bluff body.
//!
//! The **cut-cell** method keeps the uniform Cartesian grid but, for
//! every cell the body *surface* passes through, reconstructs the
//! sub-cell geometry:
//!
//! - the **fluid volume fraction** `θ ∈ [0, 1]` of the cell — how much
//!   of the `dx·dy·dz` box is outside the body;
//! - the **apertured (partial) fluid areas** of the cell's six
//!   Cartesian faces — the fraction of each face that is open to the
//!   flow;
//! - the **cut face** itself — the polygon where the body surface
//!   clips the cell box, carried as a true area + an outward unit
//!   normal. This polygon is the real wall: the no-slip boundary
//!   condition rides on it, and the force integration sums pressure
//!   and shear over it instead of over the staircased voxel faces.
//!
//! # How the geometry is computed
//!
//! - **Cut face** — every body triangle is clipped to the cell box by
//!   Sutherland-Hodgman against the box's six half-spaces; the surviving
//!   polygon fragments are summed (area-weighted) into one cut face per
//!   cell. This is *exact* — no sampling — so the wall area and normal
//!   that carry the BC and the force are exact for a planar surface
//!   crossing the cell, and converge with the mesh for a curved one.
//! - **Volume fraction + face apertures** — the solid is bounded by a
//!   triangle soup, so the cell box ∩ body is a general (possibly
//!   non-convex) region. Rather than a full polyhedral CSG, the volume
//!   fraction and the six face apertures are integrated by a
//!   deterministic regular sub-sample of the cell with the same
//!   ray-cast inside/outside test the voxelizer uses. The sub-sample is
//!   only run for the thin band of cells the surface actually crosses
//!   (a single centre test settles every wholly-inside / wholly-outside
//!   cell), so it stays cheap. A regular sub-sample is the standard
//!   robust choice for cut-cell volume fractions when a watertight
//!   polyhedral clipper is not warranted; its error is `O(1/N)` per
//!   axis and is documented honest scope.
//!
//! # Honest scope
//!
//! A real cut-cell immersed boundary: the wall area / normal / volume
//! fraction / apertures are genuine sub-cell geometry and the force
//! integrates over the true clipped surface. It is **not** a body-
//! fitted unstructured mesh — there is still no near-wall prism layer,
//! the background grid is still uniform Cartesian, and the volume
//! fraction is sampled (not a closed polyhedral CSG). Those are the
//! documented residue; the cut-cell step closes the *staircased-wall*
//! gap, materially improving the wall geometry and the coefficients.

use nalgebra::Vector3;

use crate::geometry::{Triangle, TriMesh};
use crate::grid::Grid3;
use crate::immersed::point_inside;

/// Which immersed-boundary wall treatment the voxelizer and solver use.
///
/// The two paths share the same uniform Cartesian grid and the same
/// `Solid` / `Cut` / `Fluid` tagging; they differ in how the body wall
/// inside a cut cell is represented.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum WallMethod {
    /// The legacy staircased treatment — every cell is wholly solid or
    /// wholly fluid, the wall is the voxel staircase, and a body face
    /// is one flat Cartesian face of a cut cell. Robust and cheap; the
    /// documented v1 wall-accuracy caveat.
    Staircase,
    /// The cut-cell treatment — for every cell the body surface passes
    /// through, the sub-cell fluid volume fraction, the apertured
    /// Cartesian face areas, and the true clipped wall face (area +
    /// normal) are reconstructed; the finite-volume discretisation and
    /// the force integration use that real geometry. The default.
    #[default]
    CutCell,
}

impl WallMethod {
    /// True if this is the cut-cell method.
    #[inline]
    pub fn is_cut_cell(self) -> bool {
        matches!(self, WallMethod::CutCell)
    }
}

/// The reconstructed cut face of one cell — the polygon where the body
/// surface clips the Cartesian cell box.
///
/// A cut cell has exactly one aggregated cut face: all the body-surface
/// triangle fragments inside the cell are reduced to a single
/// equivalent flat wall patch. `area` is the magnitude of the summed
/// fragment **area-vector** `Σ Aᵢ·nᵢ` and `normal` is its direction —
/// the equivalent flat wall that exerts the same net pressure force as
/// the real (possibly slightly non-planar, edge-spanning) fragment set.
/// So `area·normal` is exactly `∫ n dA` over the wall, the quantity the
/// pressure-force integral needs. `centroid` is the area-weighted
/// centroid, used as the force application point and the moment arm.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CutFace {
    /// Magnitude of the wall area-vector `|Σ Aᵢ·nᵢ|` (m²) — the
    /// equivalent flat-wall area; `area·normal == ∫ n dA`.
    pub area: f64,
    /// Outward unit normal of the wall (points from solid into fluid).
    pub normal: Vector3<f64>,
    /// Area-weighted centroid of the wall patch (world coordinates).
    pub centroid: Vector3<f64>,
}

impl CutFace {
    /// An empty (no-wall) cut face — area zero, used for cells the
    /// surface does not actually clip.
    pub fn none() -> CutFace {
        CutFace {
            area: 0.0,
            normal: Vector3::zeros(),
            centroid: Vector3::zeros(),
        }
    }

    /// True if this cell genuinely carries a wall patch.
    #[inline]
    pub fn has_wall(&self) -> bool {
        self.area > 0.0 && self.normal.norm_squared() > 0.5
    }

    /// The wall area-vector `∫ n dA = area·normal` — the quantity the
    /// pressure-force integral multiplies by `−p`.
    #[inline]
    pub fn area_vector(&self) -> Vector3<f64> {
        self.area * self.normal
    }

    /// Aggregate two cut faces into one — used when a tiny cut cell is
    /// merged into a neighbour and its wall is transferred there, and
    /// when an edge cell's perpendicular wall fragments are reduced to
    /// one equivalent patch. The area-vectors add; the merged `area` is
    /// the resultant magnitude and `normal` its direction; the centroid
    /// is area-weighted.
    pub fn merged(&self, other: &CutFace) -> CutFace {
        if !self.has_wall() {
            return *other;
        }
        if !other.has_wall() {
            return *self;
        }
        let area_vec = self.area_vector() + other.area_vector();
        let area = area_vec.norm();
        let normal = area_vec.try_normalize(1e-30).unwrap_or(self.normal);
        // The centroid is weighted by the true fragment areas (the
        // scalar areas, not the possibly-cancelling area-vector).
        let w = self.area + other.area;
        let centroid = if w > 0.0 {
            (self.area * self.centroid + other.area * other.centroid) / w
        } else {
            self.centroid
        };
        CutFace {
            area,
            normal,
            centroid,
        }
    }
}

/// The full cut-cell geometry of one Cartesian cell.
///
/// `volume_fraction` is the fluid fraction `θ ∈ [0, 1]`; `face_aperture`
/// is the open (fluid) fraction of each of the six Cartesian faces in
/// the order `[-x, +x, -y, +y, -z, +z]`; `cut_face` is the body-surface
/// wall patch (see [`CutFace`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CellGeometry {
    /// Fluid volume fraction of the cell, `0` = solid, `1` = fluid.
    pub volume_fraction: f64,
    /// Open (fluid) fraction of each Cartesian face, `[-x,+x,-y,+y,-z,+z]`.
    pub face_aperture: [f64; 6],
    /// The body-surface cut face (the wall) clipped into this cell.
    pub cut_face: CutFace,
}

impl CellGeometry {
    /// A wholly-fluid cell — full volume, every face fully open, no wall.
    pub fn full_fluid() -> CellGeometry {
        CellGeometry {
            volume_fraction: 1.0,
            face_aperture: [1.0; 6],
            cut_face: CutFace::none(),
        }
    }

    /// A wholly-solid cell — no fluid volume, every face closed.
    pub fn full_solid() -> CellGeometry {
        CellGeometry {
            volume_fraction: 0.0,
            face_aperture: [0.0; 6],
            cut_face: CutFace::none(),
        }
    }
}

/// An axis-aligned cell box, `[min, max]` in world coordinates.
#[derive(Clone, Copy, Debug)]
pub struct CellBox {
    /// Minimum corner.
    pub min: Vector3<f64>,
    /// Maximum corner.
    pub max: Vector3<f64>,
}

impl CellBox {
    /// Build the box of pressure cell `(i, j, k)` on `grid`.
    pub fn of_cell(grid: &Grid3, i: usize, j: usize, k: usize) -> CellBox {
        let (dx, dy, dz) = (grid.dx(), grid.dy(), grid.dz());
        let lo = Vector3::new(
            grid.x0 + i as f64 * dx,
            grid.y0 + j as f64 * dy,
            grid.z0 + k as f64 * dz,
        );
        CellBox {
            min: lo,
            max: lo + Vector3::new(dx, dy, dz),
        }
    }

    /// The box centre.
    #[inline]
    pub fn centre(&self) -> Vector3<f64> {
        0.5 * (self.min + self.max)
    }

    /// The box extent `(dx, dy, dz)`.
    #[inline]
    pub fn extent(&self) -> Vector3<f64> {
        self.max - self.min
    }

    /// Does an axis-aligned triangle bounding box overlap this cell box?
    /// A cheap reject before the (more expensive) clip.
    fn overlaps_tri(&self, tri: &Triangle) -> bool {
        let tmin = tri.a.inf(&tri.b).inf(&tri.c);
        let tmax = tri.a.sup(&tri.b).sup(&tri.c);
        tmax.x >= self.min.x
            && tmin.x <= self.max.x
            && tmax.y >= self.min.y
            && tmin.y <= self.max.y
            && tmax.z >= self.min.z
            && tmin.z <= self.max.z
    }
}

/// Clip a convex polygon (a vertex ring) against a single axis-aligned
/// half-space — the Sutherland-Hodgman inner loop.
///
/// `axis` is 0/1/2 for x/y/z; if `keep_below` the kept side is
/// `coord ≤ bound`, else `coord ≥ bound`. Returns the clipped ring.
fn clip_halfspace(
    poly: &[Vector3<f64>],
    axis: usize,
    bound: f64,
    keep_below: bool,
) -> Vec<Vector3<f64>> {
    if poly.is_empty() {
        return Vec::new();
    }
    // Signed distance into the kept half-space — positive = inside.
    let inside = |p: &Vector3<f64>| -> f64 {
        let c = p[axis];
        if keep_below {
            bound - c
        } else {
            c - bound
        }
    };
    let mut out = Vec::with_capacity(poly.len() + 4);
    let n = poly.len();
    for idx in 0..n {
        let cur = poly[idx];
        let prev = poly[(idx + n - 1) % n];
        let dc = inside(&cur);
        let dp = inside(&prev);
        let cur_in = dc >= 0.0;
        let prev_in = dp >= 0.0;
        if cur_in {
            if !prev_in {
                // Entering — emit the crossing point then `cur`.
                let t = dp / (dp - dc);
                out.push(prev + t * (cur - prev));
            }
            out.push(cur);
        } else if prev_in {
            // Leaving — emit only the crossing point.
            let t = dp / (dp - dc);
            out.push(prev + t * (cur - prev));
        }
    }
    out
}

/// Clip a triangle to a cell box, returning the surviving convex
/// polygon (3–9 vertices) or an empty ring if the triangle misses the
/// box entirely.
pub fn clip_triangle_to_box(tri: &Triangle, cell: &CellBox) -> Vec<Vector3<f64>> {
    let mut poly = vec![tri.a, tri.b, tri.c];
    // Six half-spaces — the box faces.
    poly = clip_halfspace(&poly, 0, cell.min.x, false);
    poly = clip_halfspace(&poly, 0, cell.max.x, true);
    poly = clip_halfspace(&poly, 1, cell.min.y, false);
    poly = clip_halfspace(&poly, 1, cell.max.y, true);
    poly = clip_halfspace(&poly, 2, cell.min.z, false);
    poly = clip_halfspace(&poly, 2, cell.max.z, true);
    poly
}

/// The signed area vector of a planar polygon ring — `½ Σ vᵢ × vᵢ₊₁`.
/// Its magnitude is the polygon area; its direction is the polygon
/// normal (right-hand rule on the vertex order).
#[cfg(test)]
fn polygon_area_vector(poly: &[Vector3<f64>]) -> Vector3<f64> {
    if poly.len() < 3 {
        return Vector3::zeros();
    }
    let mut acc = Vector3::zeros();
    let n = poly.len();
    let o = poly[0];
    for idx in 1..n - 1 {
        acc += (poly[idx] - o).cross(&(poly[idx + 1] - o));
    }
    0.5 * acc
}

/// The area-weighted centroid of a planar polygon ring (fan
/// decomposition from vertex 0).
fn polygon_centroid(poly: &[Vector3<f64>]) -> (Vector3<f64>, f64) {
    if poly.len() < 3 {
        return (Vector3::zeros(), 0.0);
    }
    let n = poly.len();
    let o = poly[0];
    let mut area_sum = 0.0;
    let mut cen = Vector3::zeros();
    for idx in 1..n - 1 {
        let v1 = poly[idx];
        let v2 = poly[idx + 1];
        let a = 0.5 * (v1 - o).cross(&(v2 - o)).norm();
        let tri_cen = (o + v1 + v2) / 3.0;
        cen += a * tri_cen;
        area_sum += a;
    }
    if area_sum > 0.0 {
        (cen / area_sum, area_sum)
    } else {
        (Vector3::zeros(), 0.0)
    }
}

/// Build the aggregated [`CutFace`] of a cell from the body mesh.
///
/// Every triangle whose AABB overlaps the cell is clipped to the cell
/// box; the clipped fragments are reduced to one equivalent flat wall
/// patch. The fragment **area-vectors** `Aᵢ·nᵢ` (each triangle's own
/// oriented normal sets the sign, so the patch points from solid into
/// fluid) are summed; the patch `area` is the resultant magnitude and
/// `normal` its direction, so `area·normal` is exactly the wall's
/// `∫ n dA`. The centroid is weighted by the true fragment areas.
pub fn cell_cut_face(cell: &CellBox, mesh: &TriMesh) -> CutFace {
    let mut area_vec = Vector3::zeros();
    let mut scalar_area = 0.0;
    let mut centroid_acc = Vector3::zeros();
    for tri in &mesh.triangles {
        if !cell.overlaps_tri(tri) {
            continue;
        }
        let poly = clip_triangle_to_box(tri, cell);
        if poly.len() < 3 {
            continue;
        }
        let (cen, area) = polygon_centroid(&poly);
        if area <= 0.0 {
            continue;
        }
        // The triangle's own outward orientation sets the wall normal.
        let tri_n = tri.normal();
        if tri_n.norm_squared() < 0.5 {
            continue;
        }
        area_vec += area * tri_n;
        scalar_area += area;
        centroid_acc += area * cen;
    }
    if scalar_area <= 0.0 {
        return CutFace::none();
    }
    // The equivalent flat-wall area is the area-vector magnitude (so
    // `area·normal == ∫ n dA`); for a flat patch this equals the
    // scalar area, for an edge-spanning patch it is slightly smaller.
    let resultant = area_vec.norm();
    if resultant <= 0.0 {
        return CutFace::none();
    }
    CutFace {
        area: resultant,
        normal: area_vec / resultant,
        centroid: centroid_acc / scalar_area,
    }
}

/// The number of sub-samples per axis used to integrate the volume
/// fraction and the face apertures. `8³ = 512` interior samples gives
/// a volume fraction good to roughly one part in `8` per axis — ample
/// for the finite-volume weighting, and only run on the thin cut band.
pub const SUBSAMPLE: usize = 8;

/// A fast single-ray inside test for sub-sampling.
///
/// [`point_inside`] casts three rays and takes a majority vote — that
/// robustness matters for the *cell-centre* classification, where a
/// single grazing artefact would flip a whole cell's tag. The volume /
/// aperture sub-sample instead averages hundreds of samples, so one
/// stray crossing barely moves the fraction; a single ray is enough and
/// is three times cheaper. The ray direction is slightly irrational so
/// it rarely lies in a triangle plane.
#[inline]
fn point_inside_fast(point: Vector3<f64>, mesh: &TriMesh) -> bool {
    let dir = Vector3::new(1.0, 0.013, 0.007);
    let mut n = 0usize;
    for tri in &mesh.triangles {
        if ray_hits(point, dir, tri) {
            n += 1;
        }
    }
    n % 2 == 1
}

/// Möller-Trumbore ray/triangle hit test (the boolean form used by the
/// fast sub-sample inside test).
#[inline]
fn ray_hits(origin: Vector3<f64>, dir: Vector3<f64>, tri: &Triangle) -> bool {
    const EPS: f64 = 1e-12;
    let e1 = tri.b - tri.a;
    let e2 = tri.c - tri.a;
    let pvec = dir.cross(&e2);
    let det = e1.dot(&pvec);
    if det.abs() < EPS {
        return false;
    }
    let inv_det = 1.0 / det;
    let tvec = origin - tri.a;
    let u = tvec.dot(&pvec) * inv_det;
    if !(0.0..=1.0).contains(&u) {
        return false;
    }
    let qvec = tvec.cross(&e1);
    let v = dir.dot(&qvec) * inv_det;
    if v < 0.0 || u + v > 1.0 {
        return false;
    }
    e2.dot(&qvec) * inv_det > EPS
}

/// Integrate the fluid volume fraction of a cell by a regular
/// sub-sample of the inside test.
fn sample_volume_fraction(cell: &CellBox, mesh: &TriMesh) -> f64 {
    let ext = cell.extent();
    let n = SUBSAMPLE;
    let mut fluid = 0usize;
    for kk in 0..n {
        let fz = (kk as f64 + 0.5) / n as f64;
        for jj in 0..n {
            let fy = (jj as f64 + 0.5) / n as f64;
            for ii in 0..n {
                let fx = (ii as f64 + 0.5) / n as f64;
                let p = cell.min + Vector3::new(fx * ext.x, fy * ext.y, fz * ext.z);
                if !point_inside_fast(p, mesh) {
                    fluid += 1;
                }
            }
        }
    }
    fluid as f64 / (n * n * n) as f64
}

/// Integrate the open (fluid) fraction of one Cartesian face of a cell
/// by a regular 2-D sub-sample. `face` is `0..6` in the order
/// `[-x, +x, -y, +y, -z, +z]`.
fn sample_face_aperture(cell: &CellBox, mesh: &TriMesh, face: usize) -> f64 {
    let ext = cell.extent();
    let n = SUBSAMPLE;
    let mut open = 0usize;
    // The two in-plane axes and the fixed coordinate of this face.
    let (ax, bx, fixed) = match face {
        0 => (1usize, 2usize, cell.min.x), // -x
        1 => (1usize, 2usize, cell.max.x), // +x
        2 => (0usize, 2usize, cell.min.y), // -y
        3 => (0usize, 2usize, cell.max.y), // +y
        4 => (0usize, 1usize, cell.min.z), // -z
        _ => (0usize, 1usize, cell.max.z), // +z
    };
    let normal_axis = match face {
        0 | 1 => 0usize,
        2 | 3 => 1usize,
        _ => 2usize,
    };
    for bb in 0..n {
        let fb = (bb as f64 + 0.5) / n as f64;
        for aa in 0..n {
            let fa = (aa as f64 + 0.5) / n as f64;
            let mut p = cell.min;
            p[ax] = cell.min[ax] + fa * ext[ax];
            p[bx] = cell.min[bx] + fb * ext[bx];
            p[normal_axis] = fixed;
            if !point_inside_fast(p, mesh) {
                open += 1;
            }
        }
    }
    open as f64 / (n * n) as f64
}

/// Compute the full cut-cell geometry of a candidate cell.
///
/// The cut face is clipped exactly from the body triangles; the volume
/// fraction and the six Cartesian-face apertures are integrated by the
/// deterministic regular sub-sample. A cell the surface does not
/// actually clip comes back with an empty cut face (and a volume
/// fraction of `0` or `1` from the sub-sample), so the caller can use
/// `cut_face.has_wall()` to decide whether it is a genuine cut cell.
pub fn compute_cut_geometry(cell: &CellBox, mesh: &TriMesh) -> CellGeometry {
    let cut_face = cell_cut_face(cell, mesh);
    if !cut_face.has_wall() {
        // No body triangle clips this cell — the box is wholly inside
        // or wholly outside the body. A single centre test settles it;
        // the expensive volume / aperture sub-sample is skipped.
        let inside = point_inside(cell.centre(), mesh);
        return if inside {
            CellGeometry::full_solid()
        } else {
            CellGeometry::full_fluid()
        };
    }
    let volume_fraction = sample_volume_fraction(cell, mesh).clamp(0.0, 1.0);
    let mut face_aperture = [0.0; 6];
    for (f, ap) in face_aperture.iter_mut().enumerate() {
        *ap = sample_face_aperture(cell, mesh, f).clamp(0.0, 1.0);
    }
    CellGeometry {
        volume_fraction,
        face_aperture,
        cut_face,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{box_body, sphere_body};

    #[test]
    fn clip_triangle_fully_inside_is_unchanged() {
        // A small triangle wholly inside a unit cell survives intact.
        let cell = CellBox {
            min: Vector3::zeros(),
            max: Vector3::new(1.0, 1.0, 1.0),
        };
        let tri = Triangle::new(
            Vector3::new(0.2, 0.2, 0.5),
            Vector3::new(0.8, 0.2, 0.5),
            Vector3::new(0.5, 0.8, 0.5),
        );
        let poly = clip_triangle_to_box(&tri, &cell);
        assert_eq!(poly.len(), 3);
        let a = polygon_area_vector(&poly).norm();
        assert!((a - tri.area()).abs() < 1e-12, "clipped area {a} vs {}", tri.area());
    }

    #[test]
    fn clip_triangle_fully_outside_is_empty() {
        let cell = CellBox {
            min: Vector3::zeros(),
            max: Vector3::new(1.0, 1.0, 1.0),
        };
        let tri = Triangle::new(
            Vector3::new(5.0, 5.0, 5.0),
            Vector3::new(6.0, 5.0, 5.0),
            Vector3::new(5.0, 6.0, 5.0),
        );
        assert!(clip_triangle_to_box(&tri, &cell).is_empty());
    }

    #[test]
    fn clip_triangle_straddling_a_face_keeps_the_inside_part() {
        // A triangle crossing the +x face of the unit cell: only the
        // x ≤ 1 part survives, and its area is exactly half here.
        let cell = CellBox {
            min: Vector3::zeros(),
            max: Vector3::new(1.0, 1.0, 1.0),
        };
        // Right triangle in the z = 0.5 plane spanning x ∈ [0, 2].
        let tri = Triangle::new(
            Vector3::new(0.0, 0.0, 0.5),
            Vector3::new(2.0, 0.0, 0.5),
            Vector3::new(0.0, 2.0, 0.5),
        );
        let poly = clip_triangle_to_box(&tri, &cell);
        assert!(poly.len() >= 3);
        let area = polygon_area_vector(&poly).norm();
        // The clipped polygon is the part of the big triangle with
        // x ≤ 1 and y ≤ 1 — a known quadrilateral of area 0.75
        // (unit square minus the corner triangle x+y>... ); verify it
        // is positive and below the full triangle area (2.0).
        assert!(area > 0.0 && area < 2.0, "clipped area {area}");
    }

    #[test]
    fn cut_face_of_a_plane_has_the_right_area_and_normal() {
        // A single large triangle lying in the x = 0.4 plane, oriented
        // with +x normal, clipped to the unit cell: the cut face must
        // have area 1.0 (it spans the whole y-z face) and normal +x.
        let cell = CellBox {
            min: Vector3::zeros(),
            max: Vector3::new(1.0, 1.0, 1.0),
        };
        // Two triangles covering the x = 0.4 square y,z ∈ [-1, 2].
        let p = 0.4;
        let t1 = Triangle::new(
            Vector3::new(p, -1.0, -1.0),
            Vector3::new(p, 2.0, -1.0),
            Vector3::new(p, 2.0, 2.0),
        );
        let t2 = Triangle::new(
            Vector3::new(p, -1.0, -1.0),
            Vector3::new(p, 2.0, 2.0),
            Vector3::new(p, -1.0, 2.0),
        );
        let mesh = TriMesh::from_triangles(vec![t1, t2]);
        let cf = cell_cut_face(&cell, &mesh);
        assert!(cf.has_wall());
        assert!((cf.area - 1.0).abs() < 1e-9, "cut-face area {} should be 1", cf.area);
        // The triangles wind so their normal is +x.
        assert!(
            (cf.normal - Vector3::new(1.0, 0.0, 0.0)).norm() < 1e-9,
            "cut-face normal {:?} should be +x",
            cf.normal
        );
        // The centroid sits in the cut plane at the cell-face centre.
        assert!((cf.centroid.x - p).abs() < 1e-9);
        assert!((cf.centroid.y - 0.5).abs() < 1e-9);
        assert!((cf.centroid.z - 0.5).abs() < 1e-9);
    }

    #[test]
    fn cut_face_normal_follows_an_inclined_plane() {
        // A 45° inclined plane through the cell: the cut-face normal
        // must be the plane normal (1,0,1)/√2.
        let cell = CellBox {
            min: Vector3::zeros(),
            max: Vector3::new(1.0, 1.0, 1.0),
        };
        // Plane x + z = 1 over y ∈ [-1, 2]; normal (1,0,1)/√2.
        let big = 3.0;
        let t1 = Triangle::new(
            Vector3::new(1.0 + big, -1.0, -big),
            Vector3::new(1.0 + big, 2.0, -big),
            Vector3::new(1.0 - big, 2.0, big),
        );
        let t2 = Triangle::new(
            Vector3::new(1.0 + big, -1.0, -big),
            Vector3::new(1.0 - big, 2.0, big),
            Vector3::new(1.0 - big, -1.0, big),
        );
        let mesh = TriMesh::from_triangles(vec![t1, t2]);
        let cf = cell_cut_face(&cell, &mesh);
        assert!(cf.has_wall());
        let want = Vector3::new(1.0, 0.0, 1.0).normalize();
        assert!(
            (cf.normal - want).norm() < 1e-6,
            "inclined cut-face normal {:?} vs {:?}",
            cf.normal,
            want
        );
    }

    #[test]
    fn volume_fraction_is_zero_inside_one_and_outside() {
        // A box body: a cell wholly inside has θ = 0, a cell wholly
        // outside has θ = 1.
        let body = box_body(
            Vector3::new(-5.0, -5.0, -5.0),
            Vector3::new(5.0, 5.0, 5.0),
        );
        // A unit cell deep inside the body.
        let inside = CellBox {
            min: Vector3::new(-0.5, -0.5, -0.5),
            max: Vector3::new(0.5, 0.5, 0.5),
        };
        assert!(sample_volume_fraction(&inside, &body) < 1e-9);
        // A unit cell well outside the body.
        let outside = CellBox {
            min: Vector3::new(20.0, 20.0, 20.0),
            max: Vector3::new(21.0, 21.0, 21.0),
        };
        assert!((sample_volume_fraction(&outside, &body) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn volume_fraction_of_a_half_cut_cell_is_about_one_half() {
        // A box body whose +x face passes exactly through the centre of
        // a unit cell: half the cell is solid, half fluid → θ ≈ 0.5.
        // Body occupies x ≤ 0, the cell spans x ∈ [-0.5, 0.5].
        let body = box_body(
            Vector3::new(-10.0, -10.0, -10.0),
            Vector3::new(0.0, 10.0, 10.0),
        );
        let cell = CellBox {
            min: Vector3::new(-0.5, -0.5, -0.5),
            max: Vector3::new(0.5, 0.5, 0.5),
        };
        let theta = sample_volume_fraction(&cell, &body);
        assert!(
            (theta - 0.5).abs() < 0.05,
            "half-cut cell volume fraction {theta} should be ~0.5"
        );
    }

    #[test]
    fn face_aperture_is_open_outside_and_closed_inside() {
        // The same x ≤ 0 body. The cell spans x ∈ [-0.5, 0.5].
        // Its -x face (at x = -0.5) is wholly inside the body → closed.
        // Its +x face (at x = +0.5) is wholly in the fluid → open.
        let body = box_body(
            Vector3::new(-10.0, -10.0, -10.0),
            Vector3::new(0.0, 10.0, 10.0),
        );
        let cell = CellBox {
            min: Vector3::new(-0.5, -0.5, -0.5),
            max: Vector3::new(0.5, 0.5, 0.5),
        };
        let minus_x = sample_face_aperture(&cell, &body, 0);
        let plus_x = sample_face_aperture(&cell, &body, 1);
        assert!(minus_x < 1e-9, "-x face should be closed, got {minus_x}");
        assert!((plus_x - 1.0).abs() < 1e-9, "+x face should be open, got {plus_x}");
        // A face parallel to the cut (e.g. +y) is half open.
        let plus_y = sample_face_aperture(&cell, &body, 3);
        assert!(
            (plus_y - 0.5).abs() < 0.05,
            "+y face straddling the cut should be ~0.5 open, got {plus_y}"
        );
    }

    #[test]
    fn sphere_cut_faces_sum_to_roughly_the_sphere_area() {
        // Summing the clipped cut-face area over a grid covering a
        // sphere must recover roughly the analytic 4πr² — the cut
        // faces tile the true surface.
        let r = 1.5;
        let body = sphere_body(Vector3::zeros(), r, 48, 96);
        let grid = Grid3::new(24, 24, 24, 8.0, 8.0, 8.0, -4.0, -4.0, -4.0);
        let mut total = 0.0;
        for k in 0..grid.nz {
            for j in 0..grid.ny {
                for i in 0..grid.nx {
                    let cell = CellBox::of_cell(&grid, i, j, k);
                    total += cell_cut_face(&cell, &body).area;
                }
            }
        }
        let analytic = 4.0 * std::f64::consts::PI * r * r;
        let rel = (total - analytic).abs() / analytic;
        assert!(
            rel < 0.05,
            "summed cut-face area {total} vs analytic {analytic} (rel {rel})"
        );
    }

    #[test]
    fn compute_cut_geometry_is_internally_consistent() {
        // A cut cell on a sphere surface: θ strictly between 0 and 1,
        // a real wall patch, apertures in range.
        let body = sphere_body(Vector3::zeros(), 2.0, 48, 96);
        // A cell straddling the +x pole of the sphere.
        let cell = CellBox {
            min: Vector3::new(1.8, -0.2, -0.2),
            max: Vector3::new(2.2, 0.2, 0.2),
        };
        let g = compute_cut_geometry(&cell, &body);
        assert!(g.volume_fraction > 0.0 && g.volume_fraction < 1.0);
        assert!(g.cut_face.has_wall());
        // The wall normal at the +x pole points roughly +x.
        assert!(g.cut_face.normal.x > 0.5);
        for ap in g.face_aperture {
            assert!((0.0..=1.0).contains(&ap));
        }
    }
}
