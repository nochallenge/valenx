//! Generator dispatch — produce a [`Vec<Placement>`] from a lattice
//! recipe.
//!
//! ## Per-placement orientation (Phase 28.5)
//!
//! Every generator derives a real instance orientation from local
//! geometry rather than leaving it as the identity:
//!
//! - **Grid** — identity (an axis-aligned grid has no intrinsic
//!   per-cell direction).
//! - **Polar** — rotation about the array axis by the per-step angle.
//! - **Bezier / OnCurve** — the source solid's local **+Z** is aligned
//!   with the curve **tangent** at the sample parameter (the tangent
//!   is the analytic NURBS derivative for `OnCurve`, the de Casteljau
//!   hodograph for `Bezier`).
//! - **OnSurface** — local +Z is aligned with the surface **normal**
//!   (cross product of the finite-difference u- and v-tangents).
//! - **OnMesh** — vertices use the area-weighted vertex normal; face
//!   centroids use the triangle's face normal.
//!
//! Orientations are computed with the internal `orient_z_to` helper,
//! which yields the minimal rotation carrying +Z onto the target
//! direction.

use nalgebra::{Unit, UnitQuaternion, Vector3};

use crate::error::LatticeError;
use crate::lattice::{Lattice, MeshSamplingMode};
use crate::placement::Placement;

/// The source-frame reference axis that gets aligned with a derived
/// tangent / normal. Solids are conventionally modelled "up" the +Z
/// axis, so +Z is the natural axis to steer.
const REFERENCE_AXIS: Vector3<f64> = Vector3::new(0.0, 0.0, 1.0);

/// Round-10 fix: cap on the per-recipe placement count. A hostile
/// `.ron` (or a typo) like
/// `Lattice::Grid { rows: usize::MAX, cols: usize::MAX, levels: 1, ... }`
/// pre-fix flowed into `Vec::with_capacity(rows * cols * levels)`,
/// which wraps `usize` and then back-fills with billions of placements
/// — the process OOMs before the caller sees an error.
///
/// 10 million placements is well past any real lattice users hand-author
/// (a real CAD assembly tops out around 10⁴–10⁵ instances) but still
/// finishes in seconds on a modern box. Past the cap the call returns
/// [`LatticeError::TooManyPlacements`] before any allocation.
pub const MAX_LATTICE_PLACEMENTS: usize = 10_000_000;

/// Multiply `factors` with `checked_mul` and reject anything past
/// [`MAX_LATTICE_PLACEMENTS`]. Returns the validated count.
///
/// Used by every generator that derives its placement count from a
/// product of caller-supplied counts (Grid: `rows·cols·levels`;
/// OnSurface: `n_u·n_v`). The Vec::with_capacity that follows now
/// has a value the allocator can honour without wrapping.
fn checked_placement_count(factors: &[usize]) -> Result<usize, LatticeError> {
    let mut count: usize = 1;
    for &f in factors {
        count = count
            .checked_mul(f)
            .ok_or(LatticeError::TooManyPlacements {
                // Promote to u128 so the error message survives the
                // overflow — `count * f` already overflowed `usize`,
                // but `u128` can still represent it for the report.
                count: (count as u128).saturating_mul(f as u128),
                max: MAX_LATTICE_PLACEMENTS,
            })?;
    }
    if count > MAX_LATTICE_PLACEMENTS {
        return Err(LatticeError::TooManyPlacements {
            count: count as u128,
            max: MAX_LATTICE_PLACEMENTS,
        });
    }
    Ok(count)
}

/// Minimal rotation that carries the source [`REFERENCE_AXIS`] (+Z)
/// onto `target`. Returns the identity when `target` is degenerate
/// (near-zero length) — there is no defined direction to steer to.
fn orient_z_to(target: Vector3<f64>) -> UnitQuaternion<f64> {
    let len = target.norm();
    if len < 1e-12 {
        return UnitQuaternion::identity();
    }
    let dir = target / len;
    UnitQuaternion::rotation_between(&REFERENCE_AXIS, &dir)
        .unwrap_or_else(UnitQuaternion::identity)
}

/// Dispatch the recipe to the matching per-variant generator.
pub fn generate(lattice: &Lattice) -> Result<Vec<Placement>, LatticeError> {
    match lattice {
        Lattice::Grid {
            rows,
            cols,
            levels,
            spacing,
        } => grid(*rows, *cols, *levels, *spacing),
        Lattice::Polar {
            center,
            axis,
            count,
            total_angle,
        } => polar(*center, *axis, *count, *total_angle),
        Lattice::Bezier {
            control_points,
            n_samples,
        } => bezier(control_points, *n_samples),
        Lattice::OnCurve { curve, n_samples } => on_curve(curve, *n_samples),
        Lattice::OnSurface { surface, n_u, n_v } => on_surface(surface, *n_u, *n_v),
        Lattice::OnMesh { mesh, mode } => on_mesh(mesh, *mode),
    }
}

fn grid(
    rows: usize,
    cols: usize,
    levels: usize,
    spacing: Vector3<f64>,
) -> Result<Vec<Placement>, LatticeError> {
    if rows == 0 || cols == 0 || levels == 0 {
        return Err(LatticeError::BadParameter {
            name: "rows/cols/levels",
            reason: "must all be > 0".into(),
        });
    }
    // Round-10 fix: rows · cols · levels is three Deserialize-fed
    // usize values multiplied with no check. Past the cap (or on
    // overflow), bail before `Vec::with_capacity` wraps and asks
    // the allocator for an exabyte.
    let capacity = checked_placement_count(&[rows, cols, levels])?;
    let mut out = Vec::with_capacity(capacity);
    for x in 0..rows {
        for y in 0..cols {
            for z in 0..levels {
                let p = Vector3::new(
                    x as f64 * spacing.x,
                    y as f64 * spacing.y,
                    z as f64 * spacing.z,
                );
                out.push(Placement::at(p));
            }
        }
    }
    Ok(out)
}

fn polar(
    center: Vector3<f64>,
    axis: Vector3<f64>,
    count: usize,
    total_angle: f64,
) -> Result<Vec<Placement>, LatticeError> {
    if count == 0 {
        return Err(LatticeError::BadParameter {
            name: "count",
            reason: "must be > 0".into(),
        });
    }
    let axis_norm = axis.norm();
    if axis_norm < f64::EPSILON {
        return Err(LatticeError::Degenerate("axis has zero length".into()));
    }
    let unit_axis = Unit::new_normalize(axis);
    // Round-10 fix: clamp the placement count even on the single-
    // factor Polar / Bezier / OnCurve variants. `count = usize::MAX`
    // pre-fix flowed straight into `Vec::with_capacity` (which the
    // allocator then refuses, panicking the process).
    let capacity = checked_placement_count(&[count])?;
    let mut out = Vec::with_capacity(capacity);
    for i in 0..count {
        let t = i as f64 / count.max(1) as f64;
        let angle = t * total_angle;
        let rot = UnitQuaternion::from_axis_angle(&unit_axis, angle);
        out.push(Placement {
            position: center,
            orientation: rot,
        });
    }
    Ok(out)
}

/// de Casteljau on a degree-`n` Bezier at parameter `t`.
fn bezier_at(cps: &[Vector3<f64>], t: f64) -> Vector3<f64> {
    let mut work: Vec<Vector3<f64>> = cps.to_vec();
    let n = work.len();
    for k in 1..n {
        for i in 0..n - k {
            work[i] = work[i] * (1.0 - t) + work[i + 1] * t;
        }
    }
    work[0]
}

/// Derivative (hodograph) of a degree-`n` Bezier at parameter `t`.
/// The Bezier derivative is itself a degree-(n-1) Bezier whose
/// control points are `n · (P_{i+1} - P_i)`; we de Casteljau that.
fn bezier_tangent(cps: &[Vector3<f64>], t: f64) -> Vector3<f64> {
    let n = cps.len();
    if n < 2 {
        return Vector3::zeros();
    }
    let deg = (n - 1) as f64;
    let hodograph: Vec<Vector3<f64>> = (0..n - 1)
        .map(|i| (cps[i + 1] - cps[i]) * deg)
        .collect();
    bezier_at(&hodograph, t)
}

fn bezier(cps: &[Vector3<f64>], n_samples: usize) -> Result<Vec<Placement>, LatticeError> {
    if cps.len() < 2 {
        return Err(LatticeError::Degenerate(
            "bezier needs at least 2 control points".into(),
        ));
    }
    if n_samples < 2 {
        return Err(LatticeError::BadParameter {
            name: "n_samples",
            reason: "need >= 2".into(),
        });
    }
    let capacity = checked_placement_count(&[n_samples])?;
    let mut out = Vec::with_capacity(capacity);
    for i in 0..n_samples {
        let t = i as f64 / (n_samples - 1) as f64;
        // Orient each instance so its +Z follows the curve tangent.
        let tangent = bezier_tangent(cps, t);
        out.push(Placement {
            position: bezier_at(cps, t),
            orientation: orient_z_to(tangent),
        });
    }
    Ok(out)
}

fn on_curve(
    curve: &valenx_surface::NurbsCurve,
    n_samples: usize,
) -> Result<Vec<Placement>, LatticeError> {
    if n_samples < 2 {
        return Err(LatticeError::BadParameter {
            name: "n_samples",
            reason: "need >= 2".into(),
        });
    }
    let u_min = curve.knots[curve.degree];
    let u_max = curve.knots[curve.knots.len() - 1 - curve.degree];
    let capacity = checked_placement_count(&[n_samples])?;
    let mut out = Vec::with_capacity(capacity);
    for i in 0..n_samples {
        let t = i as f64 / (n_samples - 1) as f64;
        let u = u_min + t * (u_max - u_min);
        // Orient each instance so its +Z follows the analytic curve
        // tangent (the first NURBS derivative).
        let tangent = curve.derivative(u, 1);
        out.push(Placement {
            position: curve.evaluate(u),
            orientation: orient_z_to(tangent),
        });
    }
    Ok(out)
}

fn on_surface(
    surface: &valenx_surface::NurbsSurface,
    n_u: usize,
    n_v: usize,
) -> Result<Vec<Placement>, LatticeError> {
    if n_u < 2 || n_v < 2 {
        return Err(LatticeError::BadParameter {
            name: "n_u/n_v",
            reason: "need >= 2".into(),
        });
    }
    let u_min = surface.u_knots[surface.u_degree];
    let u_max = surface.u_knots[surface.u_knots.len() - 1 - surface.u_degree];
    let v_min = surface.v_knots[surface.v_degree];
    let v_max = surface.v_knots[surface.v_knots.len() - 1 - surface.v_degree];
    // Round-10 fix: `n_u * n_v` is two Deserialize-fed usize values;
    // check the multiply before `Vec::with_capacity` honours the
    // wrapped-around value.
    let capacity = checked_placement_count(&[n_u, n_v])?;
    let mut out = Vec::with_capacity(capacity);
    for i in 0..n_u {
        let tu = i as f64 / (n_u - 1) as f64;
        let u = u_min + tu * (u_max - u_min);
        for j in 0..n_v {
            let tv = j as f64 / (n_v - 1) as f64;
            let v = v_min + tv * (v_max - v_min);
            // Orient each instance so its +Z follows the surface
            // normal at the sample (u, v).
            let normal = surface_normal(surface, u, v, (u_min, u_max), (v_min, v_max));
            out.push(Placement {
                position: surface.evaluate(u, v),
                orientation: orient_z_to(normal),
            });
        }
    }
    Ok(out)
}

/// Finite-difference surface normal at `(u, v)`: the cross product of
/// the central-difference u- and v-tangents. The step is a small
/// fraction of each parameter span, clamped inside the domain.
fn surface_normal(
    surface: &valenx_surface::NurbsSurface,
    u: f64,
    v: f64,
    (u_min, u_max): (f64, f64),
    (v_min, v_max): (f64, f64),
) -> Vector3<f64> {
    let hu = ((u_max - u_min) * 1e-4).max(1e-9);
    let hv = ((v_max - v_min) * 1e-4).max(1e-9);
    let up = (u + hu).min(u_max);
    let um = (u - hu).max(u_min);
    let vp = (v + hv).min(v_max);
    let vm = (v - hv).max(v_min);
    let du = surface.evaluate(up, v) - surface.evaluate(um, v);
    let dv = surface.evaluate(u, vp) - surface.evaluate(u, vm);
    du.cross(&dv)
}

fn on_mesh(
    mesh: &valenx_mesh::Mesh,
    mode: MeshSamplingMode,
) -> Result<Vec<Placement>, LatticeError> {
    match mode {
        MeshSamplingMode::Vertices => {
            if mesh.nodes.is_empty() {
                return Err(LatticeError::Degenerate("mesh has no vertices".into()));
            }
            // Area-weighted vertex normals: a raw face cross product
            // is already area-weighted, so summing them per vertex
            // gives the standard area-weighted normal.
            let mut normals = vec![Vector3::<f64>::zeros(); mesh.nodes.len()];
            for block in &mesh.element_blocks {
                if block.element_type != valenx_mesh::element::ElementType::Tri3 {
                    continue;
                }
                for chunk in block.connectivity.chunks_exact(3) {
                    let (i0, i1, i2) =
                        (chunk[0] as usize, chunk[1] as usize, chunk[2] as usize);
                    if i0 >= mesh.nodes.len()
                        || i1 >= mesh.nodes.len()
                        || i2 >= mesh.nodes.len()
                    {
                        continue;
                    }
                    let fn_ = (mesh.nodes[i1] - mesh.nodes[i0])
                        .cross(&(mesh.nodes[i2] - mesh.nodes[i0]));
                    normals[i0] += fn_;
                    normals[i1] += fn_;
                    normals[i2] += fn_;
                }
            }
            Ok(mesh
                .nodes
                .iter()
                .zip(normals.iter())
                .map(|(p, n)| Placement {
                    position: *p,
                    // Isolated vertices have a zero normal → identity.
                    orientation: orient_z_to(*n),
                })
                .collect())
        }
        MeshSamplingMode::FaceCentroids => {
            let mut out = Vec::new();
            for block in &mesh.element_blocks {
                if block.element_type != valenx_mesh::element::ElementType::Tri3 {
                    continue;
                }
                for chunk in block.connectivity.chunks_exact(3) {
                    let a = mesh
                        .nodes
                        .get(chunk[0] as usize)
                        .ok_or_else(|| LatticeError::Degenerate("idx oob".into()))?;
                    let b = mesh
                        .nodes
                        .get(chunk[1] as usize)
                        .ok_or_else(|| LatticeError::Degenerate("idx oob".into()))?;
                    let c = mesh
                        .nodes
                        .get(chunk[2] as usize)
                        .ok_or_else(|| LatticeError::Degenerate("idx oob".into()))?;
                    let centroid = (a + b + c) / 3.0;
                    // Orient each instance to the triangle's face normal.
                    let face_normal = (b - a).cross(&(c - a));
                    out.push(Placement {
                        position: centroid,
                        orientation: orient_z_to(face_normal),
                    });
                }
            }
            if out.is_empty() {
                return Err(LatticeError::Degenerate(
                    "mesh has no Tri3 face blocks".into(),
                ));
            }
            Ok(out)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grid_2x2x1_returns_4() {
        let p = generate(&Lattice::Grid {
            rows: 2,
            cols: 2,
            levels: 1,
            spacing: Vector3::new(1.0, 1.0, 1.0),
        })
        .unwrap();
        assert_eq!(p.len(), 4);
    }

    #[test]
    fn polar_4_quarter_circle() {
        let p = generate(&Lattice::Polar {
            center: Vector3::zeros(),
            axis: Vector3::z(),
            count: 4,
            total_angle: std::f64::consts::PI * 2.0,
        })
        .unwrap();
        assert_eq!(p.len(), 4);
    }

    #[test]
    fn polar_zero_axis_errors() {
        let r = generate(&Lattice::Polar {
            center: Vector3::zeros(),
            axis: Vector3::zeros(),
            count: 4,
            total_angle: 1.0,
        });
        assert!(matches!(r.err(), Some(LatticeError::Degenerate(_))));
    }

    #[test]
    fn bezier_n3_passes() {
        let p = generate(&Lattice::Bezier {
            control_points: vec![
                Vector3::zeros(),
                Vector3::new(1.0, 1.0, 0.0),
                Vector3::new(2.0, 0.0, 0.0),
            ],
            n_samples: 5,
        })
        .unwrap();
        assert_eq!(p.len(), 5);
    }

    #[test]
    fn bezier_under_2_cps_errors() {
        let r = generate(&Lattice::Bezier {
            control_points: vec![Vector3::zeros()],
            n_samples: 3,
        });
        assert!(matches!(r.err(), Some(LatticeError::Degenerate(_))));
    }

    #[test]
    fn on_mesh_vertices_counts_match() {
        let mut m = valenx_mesh::Mesh::new("t");
        m.nodes = vec![Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0)];
        let p = generate(&Lattice::OnMesh {
            mesh: m,
            mode: MeshSamplingMode::Vertices,
        })
        .unwrap();
        assert_eq!(p.len(), 2);
    }

    // --- Phase 28.5 orientation tests ---

    #[test]
    fn orient_z_to_aligns_with_target() {
        // The quaternion must carry +Z onto the requested direction.
        let q = orient_z_to(Vector3::new(1.0, 0.0, 0.0));
        let rotated = q * Vector3::z();
        assert!((rotated - Vector3::new(1.0, 0.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn orient_z_to_degenerate_is_identity() {
        let q = orient_z_to(Vector3::zeros());
        assert!((q.angle()).abs() < 1e-12);
    }

    #[test]
    fn bezier_orientation_follows_tangent() {
        // A straight Bezier along +X: every instance's +Z must rotate
        // onto +X (the constant tangent).
        let p = bezier(
            &[Vector3::zeros(), Vector3::new(10.0, 0.0, 0.0)],
            5,
        )
        .unwrap();
        for pl in &p {
            let local_z = pl.orientation * Vector3::z();
            assert!(
                (local_z - Vector3::new(1.0, 0.0, 0.0)).norm() < 1e-9,
                "bezier instance +Z should follow the +X tangent"
            );
        }
    }

    // --- Round-10 H1 RED→GREEN tests: lattice DoS via
    //     Vec::with_capacity(rows * cols * levels) ---

    /// Round-10 H1: pre-fix `Lattice::Grid { rows: usize::MAX, cols:
    /// usize::MAX, levels: 1, ... }` flowed straight into
    /// `Vec::with_capacity(rows * cols * levels)`, which wrapped
    /// `usize` and either panicked the allocator or back-filled with
    /// billions of placements. The generator now rejects past
    /// `MAX_LATTICE_PLACEMENTS` before any allocation.
    #[test]
    fn grid_rejects_overflowing_placement_count() {
        let r = grid(usize::MAX, usize::MAX, 1, Vector3::new(1.0, 1.0, 1.0));
        assert!(
            matches!(r.as_ref().err(), Some(LatticeError::TooManyPlacements { .. })),
            "expected TooManyPlacements, got: {:?}",
            r.as_ref().err()
        );
    }

    #[test]
    fn grid_rejects_count_past_cap_without_overflow() {
        // Past the cap but well inside usize — checks the > MAX branch
        // separately from the checked_mul overflow branch.
        let r = grid(2_001, 2_001, 3, Vector3::new(1.0, 1.0, 1.0));
        assert!(matches!(
            r.err(),
            Some(LatticeError::TooManyPlacements { .. })
        ));
    }

    #[test]
    fn on_surface_rejects_overflowing_placement_count() {
        // OnSurface multiplies n_u * n_v — give it the same overflow
        // pattern as grid. Build the simplest degree-1×1 flat patch.
        let surface = valenx_surface::NurbsSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![
                vec![Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0)],
                vec![Vector3::new(0.0, 1.0, 0.0), Vector3::new(1.0, 1.0, 0.0)],
            ],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        )
        .expect("flat unit patch");
        let r = on_surface(&surface, usize::MAX, usize::MAX);
        assert!(matches!(
            r.err(),
            Some(LatticeError::TooManyPlacements { .. })
        ));
    }

    #[test]
    fn bezier_rejects_oversized_sample_count() {
        // Single-factor variants get the cap too — `n_samples =
        // usize::MAX` pre-fix flowed straight into Vec::with_capacity.
        let r = bezier(
            &[Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0)],
            usize::MAX,
        );
        assert!(matches!(
            r.err(),
            Some(LatticeError::TooManyPlacements { .. })
        ));
    }

    #[test]
    fn on_mesh_face_centroid_orientation_follows_face_normal() {
        // One triangle in the z=0 plane wound CCW → face normal +Z →
        // orientation is the identity (source +Z already matches).
        let mut m = valenx_mesh::Mesh::new("tri");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        m.element_blocks.push(valenx_mesh::element::ElementBlock {
            element_type: valenx_mesh::element::ElementType::Tri3,
            connectivity: vec![0, 1, 2],
        });
        let p = generate(&Lattice::OnMesh {
            mesh: m,
            mode: MeshSamplingMode::FaceCentroids,
        })
        .unwrap();
        assert_eq!(p.len(), 1);
        let local_z = p[0].orientation * Vector3::z();
        assert!(
            (local_z - Vector3::z()).norm() < 1e-9,
            "+Z-facing triangle should leave the instance +Z unrotated"
        );
    }
}
