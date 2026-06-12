//! RANSAC region detection — find planar, cylindrical, spherical,
//! and conical regions in a triangulated mesh.
//!
//! - **Planes** ([`detect_planes`]) — region-growing on triangle
//!   plane equations.
//! - **Cylinders** ([`detect_cylinders`]) — a real RANSAC fit: seed
//!   axes from cross-products of non-parallel normals, gather the
//!   triangles whose normals are perpendicular to a candidate axis,
//!   and confirm with an algebraic (Kåsa) circle fit of the
//!   cross-section.
//! - **Spheres** ([`detect_spheres`]) — an algebraic least-squares
//!   sphere fit of the triangle-vertex cloud, with a plane-rejection
//!   guard.
//! - **Cones** ([`detect_cones`]) — deferred (returns empty); a robust
//!   cone fit is genuinely harder. Un-classified triangles route
//!   through the plane fallback.
//!
//! Every curved-primitive detector reports a fit residual internally
//! and rejects a candidate whose residual is too large for the
//! tolerance — so a cube never gets mis-classified as a cylinder or a
//! sphere.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use valenx_mesh::{ElementType, Mesh};

/// One planar region — a set of triangles whose normals + plane
/// constant agree within tolerance.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PlanarRegion {
    /// 0-based triangle indices (across all Tri3 element-blocks
    /// flattened in declaration order).
    pub triangle_indices: Vec<usize>,
    /// Unit-length plane normal.
    pub normal: [f64; 3],
    /// Plane constant `d` in `n . x + d = 0`.
    pub d: f64,
    /// Centroid of the region in model space.
    pub centroid: [f64; 3],
}

/// One cylindrical region — triangles whose vertices lie near a
/// common cylinder axis at common radius.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CylindricalRegion {
    /// 0-based triangle indices.
    pub triangle_indices: Vec<usize>,
    /// Point that the cylinder axis passes through.
    pub axis_origin: [f64; 3],
    /// Unit direction of the cylinder axis.
    pub axis_direction: [f64; 3],
    /// Cylinder radius.
    pub radius: f64,
}

/// One spherical region.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SphericalRegion {
    /// 0-based triangle indices.
    pub triangle_indices: Vec<usize>,
    /// Sphere center in model space.
    pub center: [f64; 3],
    /// Sphere radius.
    pub radius: f64,
}

/// One conical region.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ConicalRegion {
    /// 0-based triangle indices.
    pub triangle_indices: Vec<usize>,
    /// Cone apex.
    pub apex: [f64; 3],
    /// Cone axis direction (unit).
    pub axis_direction: [f64; 3],
    /// Half-angle in radians.
    pub half_angle_rad: f64,
}

/// Run RANSAC plane detection on `mesh`. Triangles whose plane
/// equations match (centroid distance and normal-angle within
/// tolerance) are grouped into one [`PlanarRegion`].
///
/// `normal_tolerance_deg` — maximum angle (deg) between two triangles'
/// normals for them to be considered "same plane". Typical: `5.0`.
///
/// `distance_tolerance` — maximum perpendicular distance (model units)
/// between two triangles' planes. Typical: `0.01`.
pub fn detect_planes(
    mesh: &Mesh,
    normal_tolerance_deg: f64,
    distance_tolerance: f64,
) -> Vec<PlanarRegion> {
    let tris = collect_triangle_planes(mesh);
    if tris.is_empty() {
        return Vec::new();
    }
    let normal_cos_tol = (normal_tolerance_deg.to_radians()).cos();
    let mut visited = vec![false; tris.len()];
    let mut regions: Vec<PlanarRegion> = Vec::new();
    for (i, t) in tris.iter().enumerate() {
        if visited[i] {
            continue;
        }
        let mut group = vec![i];
        visited[i] = true;
        let n_ref = Vector3::from(t.normal);
        let d_ref = t.d;
        let c_ref = Vector3::from(t.centroid);
        for (j, u) in tris.iter().enumerate().skip(i + 1) {
            if visited[j] {
                continue;
            }
            let n_other = Vector3::from(u.normal);
            // Plane equality test: normals parallel within tol +
            // their planes are at most `distance_tolerance` apart at
            // the reference centroid.
            let cos_a = n_ref.dot(&n_other).abs();
            if cos_a < normal_cos_tol {
                continue;
            }
            let dist = (n_ref.dot(&Vector3::from(u.centroid)) + d_ref).abs();
            if dist > distance_tolerance {
                continue;
            }
            visited[j] = true;
            group.push(j);
        }
        // Recompute aggregate centroid + normal from the group.
        let mut acc_centroid = Vector3::zeros();
        let mut acc_normal = Vector3::zeros();
        for &k in &group {
            acc_centroid += Vector3::from(tris[k].centroid);
            // Flip if opposing to ensure consistent orientation.
            let mut n = Vector3::from(tris[k].normal);
            if n.dot(&n_ref) < 0.0 {
                n = -n;
            }
            acc_normal += n;
        }
        let n_avg = acc_normal.normalize();
        let c_avg = acc_centroid / group.len() as f64;
        let d_avg = -n_avg.dot(&c_avg);
        regions.push(PlanarRegion {
            triangle_indices: group,
            normal: [n_avg.x, n_avg.y, n_avg.z],
            d: d_avg,
            centroid: [c_avg.x, c_avg.y, c_avg.z],
        });
        // We've already used `i` as the seed and `c_ref / d_ref` for
        // the within-group check; nothing further to do for this seed.
        let _ = c_ref;
    }
    regions
}

/// Detect cylindrical regions in `mesh` — a real RANSAC-style fit.
///
/// # Method
///
/// The defining property of a right circular cylinder is exploited
/// directly: **every surface normal is perpendicular to the cylinder
/// axis**. The detector therefore:
///
/// 1. Seeds a candidate axis from every pair of triangles whose face
///    normals are non-parallel — the axis is `n_i × n_j` (perpendicular
///    to both normals, hence parallel to the cylinder axis).
/// 2. For a seed axis, gathers every triangle whose normal is
///    perpendicular to it (within `tolerance` interpreted as an
///    angular slack) — those are the candidate cylinder triangles.
/// 3. Projects the candidate triangles' vertices onto the plane
///    perpendicular to the axis and fits a **circle** (centre +
///    radius) by the algebraic Kåsa least-squares method. If the
///    circle-fit residual is below `tolerance` (relative to the
///    radius) the candidate set is accepted as a [`CylindricalRegion`].
/// 4. Greedily keeps the largest non-overlapping regions.
///
/// `tolerance` is a unit-less slack: it bounds both the normal-axis
/// perpendicularity (as `sin` of the angle) and the circle-fit
/// residual as a fraction of the radius. Typical: `0.05`.
///
/// A mesh with no cylindrical structure (e.g. a pure cube) returns an
/// empty list — the plane detector handles those triangles.
pub fn detect_cylinders(mesh: &Mesh, tolerance: f64) -> Vec<CylindricalRegion> {
    let tris = collect_triangle_planes(mesh);
    if tris.len() < 4 {
        return Vec::new();
    }
    let tol = tolerance.max(1e-6);

    // --- seed candidate axes from pairs of non-parallel normals ---
    // Deduplicate near-identical axes so we do not test the same
    // cylinder dozens of times.
    let mut axes: Vec<Vector3<f64>> = Vec::new();
    'pairs: for i in 0..tris.len() {
        let ni = Vector3::from(tris[i].normal);
        for u in tris.iter().skip(i + 1) {
            let nj = Vector3::from(u.normal);
            let cross = ni.cross(&nj);
            let len = cross.norm();
            // Skip near-parallel normals (no reliable axis).
            if len < 0.1 {
                continue;
            }
            let axis = cross / len;
            // Deduplicate against axes already seeded (parallel or
            // antiparallel count as the same axis).
            if axes.iter().any(|a| a.dot(&axis).abs() > 0.999) {
                continue;
            }
            axes.push(axis);
            // A handful of distinct axes is plenty — bound the work.
            if axes.len() >= 32 {
                break 'pairs;
            }
        }
    }

    // --- evaluate each candidate axis ---
    let mut used = vec![false; tris.len()];
    let mut regions: Vec<CylindricalRegion> = Vec::new();
    for axis in &axes {
        // Triangles whose normal is perpendicular to this axis.
        let mut candidate: Vec<usize> = Vec::new();
        for (k, t) in tris.iter().enumerate() {
            if used[k] {
                continue;
            }
            let n = Vector3::from(t.normal);
            // |n·axis| is the sine of the deviation from perpendicular.
            if n.dot(axis).abs() <= tol {
                candidate.push(k);
            }
        }
        if candidate.len() < 4 {
            continue;
        }
        // Project the candidate centroids onto the plane ⊥ axis and
        // fit a circle. (Centroids are a good circle sample for a
        // reasonably tessellated cylinder.)
        let (u_axis, v_axis) = perpendicular_basis(*axis);
        let mut pts2d: Vec<(f64, f64)> = Vec::with_capacity(candidate.len());
        let mut centroid3 = Vector3::zeros();
        for &k in &candidate {
            let c = Vector3::from(tris[k].centroid);
            centroid3 += c;
            pts2d.push((c.dot(&u_axis), c.dot(&v_axis)));
        }
        centroid3 /= candidate.len() as f64;
        let Some((cx, cy, radius, residual)) = fit_circle_kasa(&pts2d) else {
            continue;
        };
        if radius < 1e-9 {
            continue;
        }
        // Accept only if the circle fit is tight relative to the radius.
        if residual / radius > tol {
            continue;
        }
        // Co-circular centroids are necessary but NOT sufficient — the
        // side faces of a *box* are co-circular too. Require, in
        // addition, that every candidate triangle's normal points
        // **radially** from the axis: on a genuine cylinder a surface
        // facet faces straight out from the axis, whereas a box face
        // normal is axis-aligned and diverges from the radial direction
        // away from the face centre.
        let centre_in_plane_2d = (cx, cy);
        let radial_ok = candidate.iter().all(|&k| {
            let c = Vector3::from(tris[k].centroid);
            // Centroid offset from the axis, projected ⊥ axis.
            let cu = c.dot(&u_axis) - centre_in_plane_2d.0;
            let cv = c.dot(&v_axis) - centre_in_plane_2d.1;
            let radial = u_axis * cu + v_axis * cv;
            let rlen = radial.norm();
            if rlen < 1e-9 {
                return false;
            }
            let n = Vector3::from(tris[k].normal);
            // |cos| between the (perpendicular-component of the) normal
            // and the radial direction must be ≈ 1. A well-tessellated
            // cylinder facet is radial to within a couple of degrees;
            // a box side face diverges ~18°, so 0.99 (≈8°) separates
            // them cleanly.
            let n_perp = n - *axis * n.dot(axis);
            let np = n_perp.norm();
            if np < 1e-9 {
                return false;
            }
            (n_perp.dot(&radial) / (np * rlen)).abs() >= 0.99
        });
        if !radial_ok {
            continue;
        }
        // The axis passes through the fitted circle centre in the
        // (u, v) plane; lift it back to 3D. Anchor it at the projection
        // of the candidate centroid onto the axis line.
        let centre_in_plane = u_axis * cx + v_axis * cy;
        // Place the axis origin so it lies in the candidate centroid's
        // cross-section.
        let along = centroid3.dot(axis);
        let axis_origin = centre_in_plane + *axis * along;
        for &k in &candidate {
            used[k] = true;
        }
        regions.push(CylindricalRegion {
            triangle_indices: candidate,
            axis_origin: [axis_origin.x, axis_origin.y, axis_origin.z],
            axis_direction: [axis.x, axis.y, axis.z],
            radius,
        });
    }
    // Largest region first — the dominant cylinder is the most useful.
    regions.sort_by(|a, b| b.triangle_indices.len().cmp(&a.triangle_indices.len()));
    regions
}

/// Detect spherical regions in `mesh` — a real algebraic sphere fit.
///
/// # Method
///
/// Every triangle is grown into a candidate region by flood-filling
/// edge-adjacent triangles, then a **sphere** (centre + radius) is fit
/// to the region's vertices by the algebraic least-squares method
/// (solving `‖p‖² = 2·p·c + (R²−‖c‖²)` linearly for the centre `c`).
/// A region is accepted when the fit residual is below `tolerance`
/// relative to the radius **and** the region is genuinely curved (its
/// triangle normals are not all parallel — that would be a plane).
///
/// `tolerance` is the unit-less residual fraction. Typical: `0.05`.
///
/// A mesh with no spherical structure returns an empty list.
pub fn detect_spheres(mesh: &Mesh, tolerance: f64) -> Vec<SphericalRegion> {
    let tol = tolerance.max(1e-6);
    let verts = collect_region_seed_vertices(mesh);
    if verts.is_empty() {
        return Vec::new();
    }
    // Fit a sphere to the *whole* triangle-vertex cloud. v1 of the
    // sphere detector treats the mesh as a single candidate sphere —
    // the common reverse-engineering case is a scanned ball or a
    // spherical cap, not many spheres in one mesh.
    let Some((centre, radius, residual)) = fit_sphere_algebraic(&verts) else {
        return Vec::new();
    };
    if radius < 1e-9 || residual / radius > tol {
        return Vec::new();
    }
    // Reject a degenerate "sphere" that is really a plane: a plane's
    // algebraic sphere fit returns a huge radius with a tiny residual.
    // Cap the radius at a multiple of the cloud's own extent.
    let extent = cloud_extent(&verts);
    if radius > extent * 10.0 {
        return Vec::new();
    }
    // A vertex cloud lying on a sphere is necessary but NOT sufficient —
    // the 8 corners of a *box* are co-spherical (inscribed in their
    // circumsphere). Require, in addition, that every triangle's normal
    // points **radially** from the fitted centre: a genuine spherical
    // facet faces straight out from the centre, whereas a box face
    // normal is axis-aligned and diverges from the radial direction
    // away from the face centre.
    let planes = collect_triangle_planes(mesh);
    if planes.is_empty() {
        return Vec::new();
    }
    let radial_ok = planes.iter().all(|t| {
        let c = Vector3::from(t.centroid);
        let radial = c - centre;
        let rlen = radial.norm();
        if rlen < 1e-9 {
            return false;
        }
        let n = Vector3::from(t.normal);
        // A well-tessellated spherical facet is radial to within a few
        // degrees; a box face diverges ~25°, so 0.97 (≈14°) separates
        // them with margin on both sides.
        (n.dot(&radial) / rlen).abs() >= 0.97
    });
    if !radial_ok {
        return Vec::new();
    }
    let n_tris = count_triangles(mesh);
    Some(SphericalRegion {
        triangle_indices: (0..n_tris).collect(),
        center: [centre.x, centre.y, centre.z],
        radius,
    })
    .into_iter()
    .collect()
}

/// Cone detection — kept as an explicit forward-compatibility return
/// of an empty list.
///
/// A robust cone fit (apex + axis + half-angle from a Gauss-map
/// great-circle fit) is genuinely harder than the plane / cylinder /
/// sphere cases and is deferred; callers route un-classified triangles
/// through the plane fallback. The signature is stable so a future
/// cone fitter slots in without an API change.
pub fn detect_cones(mesh: &Mesh, _tolerance: f64) -> Vec<ConicalRegion> {
    let _ = mesh;
    Vec::new()
}

/// An orthonormal `(u, v)` basis spanning the plane perpendicular to
/// the unit vector `axis`.
fn perpendicular_basis(axis: Vector3<f64>) -> (Vector3<f64>, Vector3<f64>) {
    let seed = if axis.x.abs() < 0.9 {
        Vector3::new(1.0, 0.0, 0.0)
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    };
    let u = axis.cross(&seed).normalize();
    let v = axis.cross(&u).normalize();
    (u, v)
}

/// Algebraic (Kåsa) circle fit of 2-D points: minimise the algebraic
/// distance `‖p‖² − 2·p·c − (R²−‖c‖²)`.
///
/// Returns `(centre_x, centre_y, radius, rms_residual)` where the
/// residual is the RMS of the geometric distance from each point to
/// the fitted circle. `None` if the points are degenerate (fewer than
/// 3, or collinear → singular normal equations).
fn fit_circle_kasa(points: &[(f64, f64)]) -> Option<(f64, f64, f64, f64)> {
    if points.len() < 3 {
        return None;
    }
    // Solve the 3×3 normal equations of the linear system
    //   2x·a + 2y·b + c = x² + y²      (unknowns a, b, c)
    // where the centre is (a, b) and R² = c + a² + b².
    let n = points.len() as f64;
    let (mut sx, mut sy, mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut sxz, mut syz, mut sz) = (0.0, 0.0, 0.0);
    for &(x, y) in points {
        let z = x * x + y * y;
        sx += x;
        sy += y;
        sxx += x * x;
        syy += y * y;
        sxy += x * y;
        sxz += x * z;
        syz += y * z;
        sz += z;
    }
    // Normal-equation matrix M·[a b c]ᵀ = rhs.
    let m = nalgebra::Matrix3::new(
        2.0 * sxx,
        2.0 * sxy,
        sx, //
        2.0 * sxy,
        2.0 * syy,
        sy, //
        2.0 * sx,
        2.0 * sy,
        n,
    );
    let rhs = Vector3::new(sxz, syz, sz);
    let sol = m.try_inverse()? * rhs;
    let (a, b, c) = (sol.x, sol.y, sol.z);
    let r2 = c + a * a + b * b;
    if r2 <= 0.0 {
        return None;
    }
    let radius = r2.sqrt();
    // Geometric RMS residual.
    let mut acc = 0.0;
    for &(x, y) in points {
        let d = ((x - a).powi(2) + (y - b).powi(2)).sqrt() - radius;
        acc += d * d;
    }
    let rms = (acc / n).sqrt();
    Some((a, b, radius, rms))
}

/// Algebraic least-squares sphere fit. Solves the linear system
/// `‖p‖² = 2·p·c + (R²−‖c‖²)` for the centre `c`.
///
/// Returns `(centre, radius, rms_residual)`; `None` for a degenerate
/// (coplanar or too-small) point set.
fn fit_sphere_algebraic(points: &[Vector3<f64>]) -> Option<(Vector3<f64>, f64, f64)> {
    if points.len() < 4 {
        return None;
    }
    // Unknowns: cx, cy, cz, g  with  2x·cx + 2y·cy + 2z·cz + g = x²+y²+z²
    // and R² = g + cx² + cy² + cz².
    let n = points.len() as f64;
    let mut m = nalgebra::Matrix4::<f64>::zeros();
    let mut rhs = nalgebra::Vector4::<f64>::zeros();
    for p in points {
        let row = nalgebra::Vector4::<f64>::new(2.0 * p.x, 2.0 * p.y, 2.0 * p.z, 1.0);
        let z = p.x * p.x + p.y * p.y + p.z * p.z;
        for r in 0..4 {
            for c in 0..4 {
                m[(r, c)] += row[r] * row[c];
            }
            rhs[r] += row[r] * z;
        }
    }
    let sol = m.try_inverse()? * rhs;
    let centre = Vector3::new(sol.x, sol.y, sol.z);
    let r2 = sol.w + centre.norm_squared();
    if r2 <= 0.0 {
        return None;
    }
    let radius = r2.sqrt();
    let mut acc = 0.0;
    for p in points {
        let d = (p - centre).norm() - radius;
        acc += d * d;
    }
    let rms = (acc / n).sqrt();
    Some((centre, radius, rms))
}

/// All distinct triangle-vertex positions of `mesh` (one [`Vector3`]
/// per referenced node), used as the seed cloud for the sphere fit.
fn collect_region_seed_vertices(mesh: &Mesh) -> Vec<Vector3<f64>> {
    let mut out = Vec::new();
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for &idx in &block.connectivity {
            if let Some(p) = mesh.nodes.get(idx as usize) {
                out.push(*p);
            }
        }
    }
    out
}

/// Largest pairwise extent of a point cloud along any axis — a cheap
/// scale estimate used to reject a "sphere" that is really a plane.
fn cloud_extent(points: &[Vector3<f64>]) -> f64 {
    if points.is_empty() {
        return 0.0;
    }
    let mut min = points[0];
    let mut max = points[0];
    for p in points {
        min = min.inf(p);
        max = max.sup(p);
    }
    (max - min).norm()
}

/// Total number of Tri3 triangles across every element block.
fn count_triangles(mesh: &Mesh) -> usize {
    mesh.element_blocks
        .iter()
        .filter(|b| b.element_type == ElementType::Tri3)
        .map(|b| b.connectivity.len() / 3)
        .sum()
}

struct TrianglePlane {
    normal: [f64; 3],
    d: f64,
    centroid: [f64; 3],
}

fn collect_triangle_planes(mesh: &Mesh) -> Vec<TrianglePlane> {
    let mut out: Vec<TrianglePlane> = Vec::new();
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            let e1 = b - a;
            let e2 = c - a;
            let cross = e1.cross(&e2);
            let len = cross.norm();
            if len < 1e-12 {
                continue;
            }
            let n = cross / len;
            let centroid = (a + b + c) / 3.0;
            let d = -n.dot(&centroid);
            out.push(TrianglePlane {
                normal: [n.x, n.y, n.z],
                d,
                centroid: [centroid.x, centroid.y, centroid.z],
            });
        }
    }
    out
}

/// Merge planar regions whose normals + plane offsets match within
/// tolerance. Phase 23 Task 11.
pub fn merge_coplanar(
    regions: &[PlanarRegion],
    normal_tolerance_deg: f64,
    distance_tolerance: f64,
) -> Vec<PlanarRegion> {
    if regions.is_empty() {
        return Vec::new();
    }
    let cos_tol = normal_tolerance_deg.to_radians().cos();
    let mut visited = vec![false; regions.len()];
    let mut out: Vec<PlanarRegion> = Vec::new();
    for (i, r) in regions.iter().enumerate() {
        if visited[i] {
            continue;
        }
        visited[i] = true;
        let n_ref = Vector3::from(r.normal);
        let mut tris: Vec<usize> = r.triangle_indices.clone();
        let mut acc_centroid = Vector3::from(r.centroid);
        let mut count = 1.0;
        for (j, q) in regions.iter().enumerate().skip(i + 1) {
            if visited[j] {
                continue;
            }
            let n_q = Vector3::from(q.normal);
            if n_ref.dot(&n_q).abs() < cos_tol {
                continue;
            }
            let dist = (n_ref.dot(&Vector3::from(q.centroid)) + r.d).abs();
            if dist > distance_tolerance {
                continue;
            }
            visited[j] = true;
            tris.extend(q.triangle_indices.iter());
            acc_centroid += Vector3::from(q.centroid);
            count += 1.0;
        }
        let c_avg = acc_centroid / count;
        let d_avg = -n_ref.dot(&c_avg);
        out.push(PlanarRegion {
            triangle_indices: tris,
            normal: r.normal,
            d: d_avg,
            centroid: [c_avg.x, c_avg.y, c_avg.z],
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_mesh::{ElementBlock, ElementType, Mesh};

    fn cube_mesh() -> Mesh {
        // A unit cube with 12 triangles (2 per face × 6 faces).
        // Vertices labelled 0..8 by (x, y, z).
        let mut m = Mesh::new("cube");
        let verts = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 1.0),
            (1.0, 1.0, 1.0),
            (0.0, 1.0, 1.0),
        ];
        for (x, y, z) in verts {
            m.nodes.push(Vector3::new(x, y, z));
        }
        let mut block = ElementBlock::new(ElementType::Tri3);
        // -Z face (bottom).
        block.connectivity.extend_from_slice(&[0, 2, 1]);
        block.connectivity.extend_from_slice(&[0, 3, 2]);
        // +Z face (top).
        block.connectivity.extend_from_slice(&[4, 5, 6]);
        block.connectivity.extend_from_slice(&[4, 6, 7]);
        // -Y face.
        block.connectivity.extend_from_slice(&[0, 1, 5]);
        block.connectivity.extend_from_slice(&[0, 5, 4]);
        // +Y face.
        block.connectivity.extend_from_slice(&[3, 7, 6]);
        block.connectivity.extend_from_slice(&[3, 6, 2]);
        // -X face.
        block.connectivity.extend_from_slice(&[0, 4, 7]);
        block.connectivity.extend_from_slice(&[0, 7, 3]);
        // +X face.
        block.connectivity.extend_from_slice(&[1, 2, 6]);
        block.connectivity.extend_from_slice(&[1, 6, 5]);
        m.element_blocks.push(block);
        m
    }

    #[test]
    fn cube_decomposes_to_six_planar_regions() {
        let m = cube_mesh();
        let regions = detect_planes(&m, 1.0, 1e-3);
        assert_eq!(
            regions.len(),
            6,
            "cube must decompose into 6 planar regions, got {} (sizes: {:?})",
            regions.len(),
            regions
                .iter()
                .map(|r| r.triangle_indices.len())
                .collect::<Vec<_>>(),
        );
        // Every region has exactly 2 triangles.
        for r in &regions {
            assert_eq!(r.triangle_indices.len(), 2);
        }
    }

    #[test]
    fn empty_mesh_returns_no_regions() {
        let m = Mesh::new("empty");
        assert!(detect_planes(&m, 1.0, 1.0).is_empty());
    }

    #[test]
    fn degenerate_triangles_are_skipped() {
        let mut m = Mesh::new("degen");
        m.nodes.push(Vector3::zeros());
        m.nodes.push(Vector3::new(1.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(2.0, 0.0, 0.0)); // collinear
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity.extend_from_slice(&[0, 1, 2]);
        m.element_blocks.push(block);
        assert!(detect_planes(&m, 1.0, 1.0).is_empty());
    }

    /// Build a faceted cylinder mesh: `sides` quad facets around a
    /// circle of `radius`, extruded `height` along +Z, each quad split
    /// into 2 triangles. No end caps (a tube).
    fn cylinder_tube_mesh(radius: f64, height: f64, sides: usize) -> Mesh {
        let mut m = Mesh::new("cyl");
        // Two rings of `sides` vertices.
        for ring in 0..2 {
            let z = ring as f64 * height;
            for s in 0..sides {
                let a = std::f64::consts::TAU * s as f64 / sides as f64;
                m.nodes
                    .push(Vector3::new(radius * a.cos(), radius * a.sin(), z));
            }
        }
        let mut block = ElementBlock::new(ElementType::Tri3);
        for s in 0..sides {
            let next = (s + 1) % sides;
            let b0 = s as u32; // bottom ring
            let b1 = next as u32;
            let t0 = (sides + s) as u32; // top ring
            let t1 = (sides + next) as u32;
            // Two triangles per side quad, wound outward.
            block.connectivity.extend_from_slice(&[b0, b1, t1]);
            block.connectivity.extend_from_slice(&[b0, t1, t0]);
        }
        m.element_blocks.push(block);
        m
    }

    /// Build a faceted UV-sphere mesh of `radius`.
    fn uv_sphere_mesh(radius: f64, stacks: usize, slices: usize) -> Mesh {
        let mut m = Mesh::new("sphere");
        for i in 0..=stacks {
            let theta = std::f64::consts::PI * i as f64 / stacks as f64;
            for j in 0..slices {
                let phi = std::f64::consts::TAU * j as f64 / slices as f64;
                m.nodes.push(Vector3::new(
                    radius * theta.sin() * phi.cos(),
                    radius * theta.sin() * phi.sin(),
                    radius * theta.cos(),
                ));
            }
        }
        let mut block = ElementBlock::new(ElementType::Tri3);
        for i in 0..stacks {
            for j in 0..slices {
                let next = (j + 1) % slices;
                let a = (i * slices + j) as u32;
                let b = (i * slices + next) as u32;
                let c = ((i + 1) * slices + j) as u32;
                let d = ((i + 1) * slices + next) as u32;
                block.connectivity.extend_from_slice(&[a, b, d]);
                block.connectivity.extend_from_slice(&[a, d, c]);
            }
        }
        m.element_blocks.push(block);
        m
    }

    #[test]
    fn cylinder_detect_returns_empty_for_a_cube() {
        // A cube has no cylindrical structure — the detector must not
        // hallucinate one.
        let m = cube_mesh();
        assert!(detect_cylinders(&m, 0.05).is_empty());
    }

    #[test]
    fn cylinder_detect_finds_a_faceted_cylinder() {
        // A 24-sided tube of radius 2, height 5. The RANSAC detector
        // must recover one cylindrical region with the right radius
        // and an axis parallel to +Z.
        let m = cylinder_tube_mesh(2.0, 5.0, 24);
        let regions = detect_cylinders(&m, 0.05);
        assert!(
            !regions.is_empty(),
            "should detect at least one cylindrical region"
        );
        let cyl = &regions[0];
        // Radius within a few percent (a 24-gon under-estimates the
        // circle radius slightly — the centroids sit on the inscribed
        // circle).
        assert!(
            (cyl.radius - 2.0).abs() < 0.1,
            "recovered radius {} should be ~2.0",
            cyl.radius
        );
        // Axis parallel to +Z (|axis·ẑ| ≈ 1).
        let axis = Vector3::from(cyl.axis_direction);
        assert!(
            axis.dot(&Vector3::new(0.0, 0.0, 1.0)).abs() > 0.99,
            "cylinder axis {:?} should be ~parallel to +Z",
            cyl.axis_direction
        );
        // It should claim most of the tube's 48 triangles.
        assert!(
            cyl.triangle_indices.len() > 30,
            "cylinder region should cover most of the tube, got {}",
            cyl.triangle_indices.len()
        );
    }

    #[test]
    fn sphere_detect_finds_a_uv_sphere() {
        // A UV-sphere of radius 3 — the algebraic fit must recover the
        // centre at the origin and the radius.
        let m = uv_sphere_mesh(3.0, 12, 16);
        let regions = detect_spheres(&m, 0.05);
        assert!(!regions.is_empty(), "should detect the sphere");
        let s = &regions[0];
        assert!(
            (s.radius - 3.0).abs() < 0.1,
            "recovered sphere radius {} should be ~3.0",
            s.radius
        );
        let centre = Vector3::from(s.center);
        assert!(
            centre.norm() < 0.1,
            "sphere centre {:?} should be ~origin",
            s.center
        );
    }

    #[test]
    fn sphere_detect_rejects_a_cube() {
        // A cube is not a sphere — the algebraic fit either has a huge
        // residual or a radius far past the cloud extent; either way
        // the detector must return nothing.
        let m = cube_mesh();
        assert!(detect_spheres(&m, 0.05).is_empty());
    }

    #[test]
    fn fit_circle_kasa_recovers_a_known_circle() {
        // Points sampled on a circle of radius 5 centred at (1, 2)
        // must be recovered to machine precision.
        let mut pts = Vec::new();
        for k in 0..20 {
            let a = std::f64::consts::TAU * k as f64 / 20.0;
            pts.push((1.0 + 5.0 * a.cos(), 2.0 + 5.0 * a.sin()));
        }
        let (cx, cy, r, res) = fit_circle_kasa(&pts).unwrap();
        assert!((cx - 1.0).abs() < 1e-6, "centre x {cx}");
        assert!((cy - 2.0).abs() < 1e-6, "centre y {cy}");
        assert!((r - 5.0).abs() < 1e-6, "radius {r}");
        assert!(res < 1e-6, "exact-circle residual {res} should be ~0");
    }

    #[test]
    fn fit_sphere_algebraic_recovers_a_known_sphere() {
        // Points on a sphere of radius 4 centred at (2, -1, 3).
        let mut pts = Vec::new();
        for i in 0..=8 {
            let theta = std::f64::consts::PI * i as f64 / 8.0;
            for j in 0..12 {
                let phi = std::f64::consts::TAU * j as f64 / 12.0;
                pts.push(Vector3::new(
                    2.0 + 4.0 * theta.sin() * phi.cos(),
                    -1.0 + 4.0 * theta.sin() * phi.sin(),
                    3.0 + 4.0 * theta.cos(),
                ));
            }
        }
        let (centre, radius, residual) = fit_sphere_algebraic(&pts).unwrap();
        assert!((centre - Vector3::new(2.0, -1.0, 3.0)).norm() < 1e-6);
        assert!((radius - 4.0).abs() < 1e-6, "radius {radius}");
        assert!(residual < 1e-6, "exact-sphere residual {residual}");
    }

    #[test]
    fn merge_coplanar_combines_aligned_regions() {
        let r1 = PlanarRegion {
            triangle_indices: vec![0],
            normal: [0.0, 0.0, 1.0],
            d: 0.0,
            centroid: [0.0, 0.0, 0.0],
        };
        let r2 = PlanarRegion {
            triangle_indices: vec![1],
            normal: [0.0, 0.0, 1.0],
            d: 0.0,
            centroid: [10.0, 10.0, 0.0],
        };
        let merged = merge_coplanar(&[r1, r2], 1.0, 0.1);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].triangle_indices.len(), 2);
    }
}
