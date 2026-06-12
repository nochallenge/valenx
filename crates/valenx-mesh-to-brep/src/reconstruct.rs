//! Region → NURBS surface fitting + assembly into a solid.

use thiserror::Error;

use nalgebra::Vector3;
use valenx_cad::Solid;
use valenx_mesh::{ElementBlock, ElementType, Mesh};
use valenx_surface::fit::surface_through_scattered;
use valenx_surface::NurbsSurface;

use crate::feature_detect::{detect_planes, CylindricalRegion, PlanarRegion};

/// Errors during BRep reconstruction.
#[derive(Debug, Error)]
pub enum ReconstructError {
    /// The mesh had no Tri3 triangles to work with.
    #[error("mesh has no Tri3 triangles")]
    EmptyMesh,
    /// The mesh decomposed into more regions than v1 can handle.
    #[error("too many regions: {0} (v1 limit: {1})")]
    TooManyRegions(usize, usize),
    /// Surface fit failed.
    #[error("surface fit: {0}")]
    SurfaceFit(String),
}

/// Maximum number of regions v1's [`brep_from_mesh`] will attempt to
/// reconstruct from. Above this, the function falls back to a
/// pass-through (return the original mesh wrapped as `Solid::Mesh`).
pub const V1_MAX_REGIONS: usize = 10;

/// Fit a NURBS surface to a [`PlanarRegion`].
///
/// v1 returns a flat tessellated quad patch that lies in the region's
/// plane — true rational NURBS fitting per region is queued behind
/// `valenx-surface`'s public API integration.
pub fn nurbs_from_region(mesh: &Mesh, region: &PlanarRegion) -> Result<Mesh, ReconstructError> {
    // v1: collect the bounding box of the region's triangles within
    // the plane, then emit a single quad's worth of triangles aligned
    // to the plane.
    let normal = Vector3::from(region.normal).normalize();
    // Build two in-plane basis vectors `u`, `v` perpendicular to
    // `normal`.
    let arbitrary = if normal.x.abs() < 0.9 {
        Vector3::new(1.0, 0.0, 0.0)
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    };
    let u = normal.cross(&arbitrary).normalize();
    let v = normal.cross(&u).normalize();
    let centroid = Vector3::from(region.centroid);
    // Project the region's triangles' vertices into (u, v).
    let mut min_u = f64::INFINITY;
    let mut max_u = f64::NEG_INFINITY;
    let mut min_v = f64::INFINITY;
    let mut max_v = f64::NEG_INFINITY;
    for &tri_idx in &region.triangle_indices {
        let block = mesh
            .element_blocks
            .iter()
            .find(|b| b.element_type == ElementType::Tri3);
        if let Some(block) = block {
            if tri_idx * 3 + 2 < block.connectivity.len() {
                for k in 0..3 {
                    let node_idx = block.connectivity[tri_idx * 3 + k] as usize;
                    if node_idx >= mesh.nodes.len() {
                        continue;
                    }
                    let p = mesh.nodes[node_idx] - centroid;
                    let uu = p.dot(&u);
                    let vv = p.dot(&v);
                    min_u = min_u.min(uu);
                    max_u = max_u.max(uu);
                    min_v = min_v.min(vv);
                    max_v = max_v.max(vv);
                }
            }
        }
    }
    if !min_u.is_finite() {
        return Err(ReconstructError::SurfaceFit(
            "region had no recoverable triangles".into(),
        ));
    }
    let mut out = Mesh::new("nurbs_planar_fit");
    let corner_a = centroid + u * min_u + v * min_v;
    let corner_b = centroid + u * max_u + v * min_v;
    let corner_c = centroid + u * max_u + v * max_v;
    let corner_d = centroid + u * min_u + v * max_v;
    out.nodes.push(corner_a);
    out.nodes.push(corner_b);
    out.nodes.push(corner_c);
    out.nodes.push(corner_d);
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity.extend_from_slice(&[0, 1, 2]);
    block.connectivity.extend_from_slice(&[0, 2, 3]);
    out.element_blocks.push(block);
    Ok(out)
}

/// The fitted NURBS surface for one region plus a tolerance report.
///
/// Production reverse-engineering needs to *quantify* how well a fitted
/// surface matches the scan it came from — a fit is only trustworthy
/// if its deviation is within the scanner's accuracy. This carries the
/// fitted surface alongside that deviation.
#[derive(Clone, Debug)]
pub struct RegionFit {
    /// The fitted NURBS surface for the region.
    pub surface: NurbsSurface,
    /// Root-mean-squared distance (model units) from the region's
    /// sample points to the fitted surface. The headline tolerance
    /// number — compare it against the scanner accuracy.
    pub rms_error: f64,
    /// Worst single-point deviation (model units) — the largest
    /// distance from any region sample to the fitted surface. A fit
    /// can have a small RMS but a large local spike; this catches it.
    pub max_error: f64,
    /// Number of mesh sample points the fit was computed from.
    pub sample_count: usize,
}

impl RegionFit {
    /// Whether the fit meets a caller-supplied tolerance: both the RMS
    /// and the worst-case deviation are within `tolerance` model units.
    pub fn within_tolerance(&self, tolerance: f64) -> bool {
        self.rms_error <= tolerance && self.max_error <= tolerance
    }
}

/// Fit a real NURBS surface to a planar region **with a tolerance
/// report** — the production-quality path (Phase 23 v2).
///
/// Where [`nurbs_from_region`] emits a flat tessellated quad, this
/// gathers every distinct vertex of the region's triangles and fits a
/// genuine tensor-product NURBS surface through them with
/// [`valenx_surface::fit::surface_through_scattered`], then measures
/// the deviation of every sample from the fitted surface and returns
/// it in the [`RegionFit`] report.
///
/// `cps_u` / `cps_v` are the target control-point counts of the
/// fitted surface (a `3×3` bicubic-ish patch is a reasonable default
/// for a single planar region). `degree` is the surface degree in both
/// parameters.
///
/// # Errors
///
/// [`ReconstructError::SurfaceFit`] if the region has too few
/// recoverable vertices for the requested control grid, or the
/// underlying surface fit fails.
pub fn fit_planar_region_nurbs(
    mesh: &Mesh,
    region: &PlanarRegion,
    cps_u: usize,
    cps_v: usize,
    degree: usize,
) -> Result<RegionFit, ReconstructError> {
    let points = region_vertices(mesh, &region.triangle_indices);
    fit_region_points(&points, cps_u, cps_v, degree)
}

/// Fit a real NURBS surface to a **cylindrical** region with a
/// tolerance report.
///
/// The cylinder's lateral surface is developable, so a scattered NURBS
/// fit through the region's vertices reproduces it well. This shares
/// the scattered-fit machinery with the planar path and returns the
/// same [`RegionFit`] tolerance report — so a caller can hold the
/// cylinder fit to the same scanner-accuracy bar as a plane.
///
/// # Errors
///
/// [`ReconstructError::SurfaceFit`] as for [`fit_planar_region_nurbs`].
pub fn fit_cylindrical_region_nurbs(
    mesh: &Mesh,
    region: &CylindricalRegion,
    cps_u: usize,
    cps_v: usize,
    degree: usize,
) -> Result<RegionFit, ReconstructError> {
    let points = region_vertices(mesh, &region.triangle_indices);
    fit_region_points(&points, cps_u, cps_v, degree)
}

/// Shared core: fit a NURBS surface through a scattered point set and
/// build the [`RegionFit`] tolerance report.
fn fit_region_points(
    points: &[Vector3<f64>],
    cps_u: usize,
    cps_v: usize,
    degree: usize,
) -> Result<RegionFit, ReconstructError> {
    let needed = cps_u.max(degree + 1) * cps_v.max(degree + 1);
    if points.len() < needed {
        return Err(ReconstructError::SurfaceFit(format!(
            "region has {} sample points; a {cps_u}×{cps_v} NURBS fit needs at least {needed}",
            points.len()
        )));
    }
    let fit = surface_through_scattered(points, degree, degree, cps_u, cps_v)
        .map_err(|e| ReconstructError::SurfaceFit(format!("{e}")))?;

    // `valenx-surface`'s reported RMS is measured in its *internal
    // parameter space* (it evaluates the surface at assumed-uniform
    // (u,v) and differences against the binned grid) — a v1
    // simplification that over-reports whenever the scattered binning
    // warps the parameterisation. `RegionFit` promises a *geometric*
    // RMS (true model-unit distance of each sample to the surface), so
    // we measure both the RMS and the worst-case deviation ourselves by
    // projecting every sample onto the fitted surface.
    let (rms_error, max_error) = surface_deviation(&fit.surface, points);
    Ok(RegionFit {
        surface: fit.surface,
        rms_error,
        max_error,
        sample_count: points.len(),
    })
}

/// Geometric deviation of `points` from `surface`: returns
/// `(rms, max)` — the root-mean-squared and worst-case model-unit
/// distance from each sample to the fitted surface.
///
/// Each sample is projected onto the surface by [`closest_point`] — a
/// grid-seeded Gauss-Newton search for the true perpendicular foot, so
/// the distance is the genuine point-to-surface distance, not a
/// pre-sampled-grid approximation. That matters: a perfect fit then
/// reports a deviation of ~0 (a sample *on* the surface projects onto
/// itself), as a tolerance report must.
fn surface_deviation(surface: &NurbsSurface, points: &[Vector3<f64>]) -> (f64, f64) {
    let mut worst = 0.0_f64;
    let mut sumsq = 0.0_f64;
    for p in points {
        let foot = closest_point(surface, *p);
        let d = (p - foot).norm();
        sumsq += d * d;
        if d > worst {
            worst = d;
        }
    }
    let rms = (sumsq / points.len().max(1) as f64).sqrt();
    (rms, worst)
}

/// Project `p` onto `surface`, returning the nearest surface point.
///
/// A coarse grid seeds the nearest `(u, v)`, then a Gauss-Newton
/// iteration with central-difference tangents drives the residual
/// perpendicular to the surface — quadratic convergence near the
/// solution, so a sample that genuinely lies on the surface projects
/// back onto itself to machine precision.
fn closest_point(surface: &NurbsSurface, p: Vector3<f64>) -> Vector3<f64> {
    let (u_min, u_max) = surface.u_range();
    let (v_min, v_max) = surface.v_range();
    // Seed: nearest node of a coarse grid.
    let mut best_u = 0.5 * (u_min + u_max);
    let mut best_v = 0.5 * (v_min + v_max);
    let mut best_d = f64::INFINITY;
    const SEED: usize = 12;
    for i in 0..=SEED {
        for j in 0..=SEED {
            let u = u_min + i as f64 * (u_max - u_min) / SEED as f64;
            let v = v_min + j as f64 * (v_max - v_min) / SEED as f64;
            let d = (surface.evaluate(u, v) - p).norm();
            if d < best_d {
                best_d = d;
                best_u = u;
                best_v = v;
            }
        }
    }
    // Gauss-Newton refinement.
    let h = ((u_max - u_min) + (v_max - v_min)) * 1.0e-6;
    let (mut u, mut v) = (best_u, best_v);
    for _ in 0..32 {
        let q = surface.evaluate(u, v);
        let r = q - p;
        if r.norm() < 1.0e-13 {
            break;
        }
        let u_lo = (u - h).max(u_min);
        let u_hi = (u + h).min(u_max);
        let v_lo = (v - h).max(v_min);
        let v_hi = (v + h).min(v_max);
        let du_span = (u_hi - u_lo).max(1.0e-12);
        let dv_span = (v_hi - v_lo).max(1.0e-12);
        let tu = (surface.evaluate(u_hi, v) - surface.evaluate(u_lo, v)) / du_span;
        let tv = (surface.evaluate(u, v_hi) - surface.evaluate(u, v_lo)) / dv_span;
        let a11 = tu.dot(&tu);
        let a12 = tu.dot(&tv);
        let a22 = tv.dot(&tv);
        let b1 = -tu.dot(&r);
        let b2 = -tv.dot(&r);
        let det = a11 * a22 - a12 * a12;
        if det.abs() < 1.0e-15 {
            break;
        }
        u = (u + (a22 * b1 - a12 * b2) / det).clamp(u_min, u_max);
        v = (v + (-a12 * b1 + a11 * b2) / det).clamp(v_min, v_max);
    }
    surface.evaluate(u, v)
}

/// Collect the distinct vertex positions of a region's triangles.
///
/// `triangle_indices` are 0-based indices into the flattened Tri3
/// element blocks (the convention every `*Region` type uses).
fn region_vertices(mesh: &Mesh, triangle_indices: &[usize]) -> Vec<Vector3<f64>> {
    // Flatten the Tri3 connectivity once.
    let mut tri_conn: Vec<[u32; 3]> = Vec::new();
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            tri_conn.push([tri[0], tri[1], tri[2]]);
        }
    }
    // Gather vertices, de-duplicating by node index so a shared vertex
    // is not weighted twice in the fit.
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    for &t in triangle_indices {
        if let Some(tri) = tri_conn.get(t) {
            for &node in tri {
                if seen.insert(node) {
                    if let Some(p) = mesh.nodes.get(node as usize) {
                        out.push(*p);
                    }
                }
            }
        }
    }
    out
}

/// Top-level reverse-engineer entry point. Detects planar regions,
/// fits each to a NURBS patch, glues the patches into a single mesh-
/// backed [`Solid`].
///
/// v1 short-circuits to "wrap the input as a `Solid::Mesh`" if the
/// region count exceeds [`V1_MAX_REGIONS`].
pub fn brep_from_mesh(
    mesh: &Mesh,
    normal_tolerance_deg: f64,
    distance_tolerance: f64,
) -> Result<Solid, ReconstructError> {
    if mesh
        .element_blocks
        .iter()
        .all(|b| b.element_type != ElementType::Tri3)
    {
        return Err(ReconstructError::EmptyMesh);
    }
    let regions = detect_planes(mesh, normal_tolerance_deg, distance_tolerance);
    if regions.is_empty() {
        return Ok(Solid::from_mesh(mesh.clone()));
    }
    if regions.len() > V1_MAX_REGIONS {
        tracing::warn!(
            target: "valenx-mesh-to-brep",
            "{} regions exceeds v1 limit ({}); returning pass-through mesh",
            regions.len(),
            V1_MAX_REGIONS,
        );
        return Ok(Solid::from_mesh(mesh.clone()));
    }
    let mut combined = Mesh::new("reverse_engineered");
    let mut block = ElementBlock::new(ElementType::Tri3);
    for region in &regions {
        let fit = nurbs_from_region(mesh, region)?;
        let base = combined.nodes.len() as u32;
        combined.nodes.extend(fit.nodes.iter().copied());
        for fb in &fit.element_blocks {
            if fb.element_type == ElementType::Tri3 {
                for c in &fb.connectivity {
                    block.connectivity.push(*c + base);
                }
            }
        }
    }
    if !block.connectivity.is_empty() {
        combined.element_blocks.push(block);
    }
    Ok(Solid::from_mesh(combined))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature_detect::detect_planes;
    use valenx_mesh::{ElementBlock, ElementType, Mesh};

    fn cube_mesh() -> Mesh {
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
        block.connectivity.extend_from_slice(&[0, 2, 1, 0, 3, 2]);
        block.connectivity.extend_from_slice(&[4, 5, 6, 4, 6, 7]);
        block.connectivity.extend_from_slice(&[0, 1, 5, 0, 5, 4]);
        block.connectivity.extend_from_slice(&[3, 7, 6, 3, 6, 2]);
        block.connectivity.extend_from_slice(&[0, 4, 7, 0, 7, 3]);
        block.connectivity.extend_from_slice(&[1, 2, 6, 1, 6, 5]);
        m.element_blocks.push(block);
        m
    }

    #[test]
    fn brep_from_cube_returns_solid_mesh() {
        let m = cube_mesh();
        let s = brep_from_mesh(&m, 1.0, 1e-3).unwrap();
        match s {
            Solid::Mesh(out) => {
                assert!(!out.nodes.is_empty(), "should produce reconstructed mesh");
            }
            _ => panic!("expected mesh solid"),
        }
    }

    #[test]
    fn empty_mesh_returns_empty_error() {
        let m = Mesh::new("empty");
        let err = brep_from_mesh(&m, 1.0, 1e-3).unwrap_err();
        assert!(matches!(err, ReconstructError::EmptyMesh));
    }

    #[test]
    fn nurbs_from_region_handles_cube_face() {
        let m = cube_mesh();
        let regions = detect_planes(&m, 1.0, 1e-3);
        let fit = nurbs_from_region(&m, &regions[0]).unwrap();
        // Four corners + two triangles.
        assert_eq!(fit.nodes.len(), 4);
        assert_eq!(fit.element_blocks[0].connectivity.len(), 6);
    }

    /// A flat n×n grid of triangles in the z=0 plane spanning the unit
    /// square — a dense planar region for the NURBS fit tests.
    fn flat_grid_mesh(n: usize) -> Mesh {
        let mut m = Mesh::new("grid");
        for j in 0..=n {
            for i in 0..=n {
                m.nodes
                    .push(Vector3::new(i as f64 / n as f64, j as f64 / n as f64, 0.0));
            }
        }
        let stride = (n + 1) as u32;
        let mut block = ElementBlock::new(ElementType::Tri3);
        for j in 0..n as u32 {
            for i in 0..n as u32 {
                let a = j * stride + i;
                let b = a + 1;
                let c = a + stride;
                let d = c + 1;
                block.connectivity.extend_from_slice(&[a, b, d]);
                block.connectivity.extend_from_slice(&[a, d, c]);
            }
        }
        m.element_blocks.push(block);
        m
    }

    #[test]
    fn fit_planar_region_nurbs_fits_a_flat_grid_tightly() {
        // A genuinely flat region must fit a NURBS surface with a tiny
        // deviation — the headline tolerance report should be ~0.
        let m = flat_grid_mesh(6);
        let regions = detect_planes(&m, 1.0, 1e-3);
        assert_eq!(regions.len(), 1, "flat grid is one planar region");
        let fit = fit_planar_region_nurbs(&m, &regions[0], 3, 3, 2).unwrap();
        // A flat region fits a flat NURBS surface essentially exactly.
        assert!(
            fit.rms_error < 1e-6,
            "flat-region RMS {} should be ~0",
            fit.rms_error
        );
        assert!(
            fit.max_error < 1e-3,
            "flat-region worst deviation {} should be tiny",
            fit.max_error
        );
        assert!(fit.sample_count > 0, "fit should report a sample count");
        // The tolerance gate should pass at a loose tolerance.
        assert!(fit.within_tolerance(0.01));
    }

    #[test]
    fn fit_planar_region_nurbs_rejects_a_too_sparse_region() {
        // A single cube face has only 4 vertices — far too few for a
        // 3×3 NURBS control grid. The fit must error, not panic.
        let m = cube_mesh();
        let regions = detect_planes(&m, 1.0, 1e-3);
        let err = fit_planar_region_nurbs(&m, &regions[0], 3, 3, 2).unwrap_err();
        assert!(matches!(err, ReconstructError::SurfaceFit(_)));
    }

    #[test]
    fn region_fit_within_tolerance_checks_both_errors() {
        // within_tolerance must require BOTH the RMS and the worst-case
        // deviation to be inside the bound.
        let surface = {
            // Reuse a real flat fit to get a valid NurbsSurface.
            let m = flat_grid_mesh(6);
            let regions = detect_planes(&m, 1.0, 1e-3);
            fit_planar_region_nurbs(&m, &regions[0], 3, 3, 2)
                .unwrap()
                .surface
        };
        let good = RegionFit {
            surface: surface.clone(),
            rms_error: 0.001,
            max_error: 0.005,
            sample_count: 49,
        };
        assert!(good.within_tolerance(0.01));
        // A small RMS but a large spike must fail.
        let spiky = RegionFit {
            surface,
            rms_error: 0.001,
            max_error: 0.5,
            sample_count: 49,
        };
        assert!(
            !spiky.within_tolerance(0.01),
            "a local spike must fail the gate"
        );
    }
}
