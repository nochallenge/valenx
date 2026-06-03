//! Structural member — a profile swept along a linear-segment path.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::FramesError;
use crate::profile::{cross_section_polygon, Profile};

/// One structural member.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Member {
    /// Path centreline, world coordinates (mm), as a polyline. A
    /// curved member is represented by a finely-sampled polyline —
    /// see [`Member::from_nurbs_curve`].
    pub path: Vec<Vector3<f64>>,
    /// Cross-section profile.
    pub profile: Profile,
    /// Rotation of the profile about the member's tangent axis,
    /// in degrees. 0 means the profile's local +v aligns with the
    /// world +Z (or with the closest projection if path tangent is
    /// near-vertical).
    pub orientation_angle_deg: f64,
}

impl Member {
    /// Convenience: straight member from `a` to `b`.
    pub fn straight(a: Vector3<f64>, b: Vector3<f64>, profile: Profile) -> Self {
        Self {
            path: vec![a, b],
            profile,
            orientation_angle_deg: 0.0,
        }
    }

    /// Build a **curved** member by sweeping `profile` along a
    /// `valenx_surface::NurbsCurve` (Phase 38.5).
    ///
    /// The NURBS curve is sampled at `n_samples` equally-spaced
    /// parameter values across its valid knot range; the resulting
    /// points become the member's polyline `path`. [`to_solid`] then
    /// sweeps the section through that path with the existing
    /// per-station-tangent framing, so a NURBS-defined arc, helix or
    /// spline produces a smoothly curved structural member.
    ///
    /// # Errors
    ///
    /// - [`FramesError::BadParameter`] when `n_samples < 2`.
    pub fn from_nurbs_curve(
        curve: &valenx_surface::NurbsCurve,
        n_samples: usize,
        profile: Profile,
    ) -> Result<Self, FramesError> {
        if n_samples < 2 {
            return Err(FramesError::BadParameter {
                name: "n_samples",
                reason: format!("need >= 2 samples to form a path, got {n_samples}"),
            });
        }
        // Valid parameter range of a clamped NURBS curve.
        let u_min = curve.knots[curve.degree];
        let u_max = curve.knots[curve.knots.len() - 1 - curve.degree];
        let mut path = Vec::with_capacity(n_samples);
        for i in 0..n_samples {
            let t = i as f64 / (n_samples - 1) as f64;
            let u = u_min + t * (u_max - u_min);
            path.push(curve.evaluate(u));
        }
        Ok(Self {
            path,
            profile,
            orientation_angle_deg: 0.0,
        })
    }

    /// Total cumulative length (mm).
    pub fn length_mm(&self) -> f64 {
        let mut total = 0.0;
        for w in self.path.windows(2) {
            total += (w[1] - w[0]).norm();
        }
        total
    }
}

/// Pick a local frame `(u, v, t)` for the tangent direction `t`. The
/// algorithm grabs an arbitrary perpendicular for `u` (preferring +X
/// when the tangent is near-vertical) then takes `v = t × u`.
fn local_frame(tangent: Vector3<f64>) -> (Vector3<f64>, Vector3<f64>, Vector3<f64>) {
    let t = tangent.normalize();
    // Reference up axis = world +Z unless tangent is near-parallel.
    let world_up = Vector3::new(0.0, 0.0, 1.0);
    let dot = t.dot(&world_up).abs();
    let ref_up = if dot > 0.99 {
        Vector3::new(0.0, 1.0, 0.0)
    } else {
        world_up
    };
    let u = ref_up.cross(&t).normalize();
    let v = t.cross(&u).normalize();
    (u, v, t)
}

/// Apply the optional roll rotation about the tangent axis.
fn rotate_about_tangent(u: Vector3<f64>, v: Vector3<f64>, angle_deg: f64) -> (Vector3<f64>, Vector3<f64>) {
    if angle_deg == 0.0 {
        return (u, v);
    }
    let a = angle_deg.to_radians();
    let cs = a.cos();
    let sn = a.sin();
    let u2 = u * cs + v * sn;
    let v2 = -u * sn + v * cs;
    (u2, v2)
}

/// Sweep the cross-section polygon along the member's polyline path.
///
/// Each path vertex gets a cross-section ring framed by the local
/// tangent (averaged at interior vertices); consecutive rings are
/// stitched into tube walls. A curved member built with
/// [`Member::from_nurbs_curve`] is a finely-sampled polyline, so the
/// same sweep produces a smooth curved member — there is no
/// linear-only restriction.
///
/// Returns a [`Solid::Mesh`] with one triangulated tube per path
/// segment. End caps fan to the polygon centroid.
pub fn to_solid(m: &Member) -> Result<Solid, FramesError> {
    if m.path.len() < 2 {
        return Err(FramesError::DegeneratePath(format!(
            "path needs >= 2 vertices, got {}",
            m.path.len()
        )));
    }
    let poly2d = cross_section_polygon(m.profile);
    if poly2d.len() < 3 {
        return Err(FramesError::BadParameter {
            name: "profile",
            reason: format!("polygon has {} vertices", poly2d.len()),
        });
    }

    let mut mesh = Mesh::new(format!("frames_member_{}", m.profile.label()));
    let mut block = ElementBlock::new(ElementType::Tri3);

    let n = poly2d.len();
    let stations = m.path.len();
    // Per-station vertices: (station × poly_index).
    let mut ring_base: Vec<u32> = Vec::with_capacity(stations);
    for (i, station) in m.path.iter().enumerate() {
        // Pick the local frame from the tangent at this station.
        let tangent = if i == 0 {
            m.path[1] - m.path[0]
        } else if i == stations - 1 {
            m.path[i] - m.path[i - 1]
        } else {
            (m.path[i + 1] - m.path[i - 1]).normalize()
        };
        let (mut u, mut v, _t) = local_frame(tangent);
        let (uu, vv) = rotate_about_tangent(u, v, m.orientation_angle_deg);
        u = uu;
        v = vv;
        let base = mesh.nodes.len() as u32;
        ring_base.push(base);
        for p in &poly2d {
            let world = station + u * p[0] + v * p[1];
            mesh.nodes.push(world);
        }
    }

    // Side walls: stitch each adjacent pair of rings.
    for s in 0..(stations - 1) {
        let a = ring_base[s];
        let b = ring_base[s + 1];
        for i in 0..n {
            let j = (i + 1) % n;
            let p00 = a + i as u32;
            let p01 = a + j as u32;
            let p10 = b + i as u32;
            let p11 = b + j as u32;
            block.connectivity.extend_from_slice(&[p00, p01, p11]);
            block.connectivity.extend_from_slice(&[p00, p11, p10]);
        }
    }

    // End caps: fan to centroid at first + last station.
    let mut cap_centroid =
        |station: Vector3<f64>, base: u32, _flip: bool| -> u32 {
            // Compute centroid from the ring nodes already pushed for
            // this station — they're consecutive in `nodes` starting
            // at `base`.
            let mut c = Vector3::zeros();
            for i in 0..n {
                c += mesh.nodes[base as usize + i];
            }
            c /= n as f64;
            // Replace one component with the station to keep cap
            // co-planar with the section (best-effort: caps live in
            // the section plane already).
            let _ = station;
            mesh.nodes.push(c);
            mesh.nodes.len() as u32 - 1
        };
    let c_first = cap_centroid(m.path[0], ring_base[0], true);
    let c_last = cap_centroid(*m.path.last().unwrap(), *ring_base.last().unwrap(), false);
    for i in 0..n {
        let j = (i + 1) % n;
        // First cap fan — orient opposite to side walls.
        block.connectivity.extend_from_slice(&[
            c_first,
            ring_base[0] + j as u32,
            ring_base[0] + i as u32,
        ]);
        // Last cap fan.
        let base_last = *ring_base.last().unwrap();
        block.connectivity.extend_from_slice(&[
            c_last,
            base_last + i as u32,
            base_last + j as u32,
        ]);
    }

    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn straight_ibeam_member_meshes() {
        let m = Member::straight(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(2000.0, 0.0, 0.0),
            Profile::default_ipe200(),
        );
        let s = to_solid(&m).expect("mesh");
        match s {
            Solid::Mesh(mesh) => {
                assert!(!mesh.nodes.is_empty());
                assert!(mesh.total_elements() > 0);
            }
            _ => panic!("expected mesh-backed solid"),
        }
        assert!((m.length_mm() - 2000.0).abs() < 1e-6);
    }

    #[test]
    fn degenerate_single_point_errors() {
        let m = Member {
            path: vec![Vector3::zeros()],
            profile: Profile::default_ipe200(),
            orientation_angle_deg: 0.0,
        };
        assert!(matches!(to_solid(&m), Err(FramesError::DegeneratePath(_))));
    }

    #[test]
    fn two_segment_path_produces_two_tubes() {
        let m = Member {
            path: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(1000.0, 0.0, 0.0),
                Vector3::new(1000.0, 1000.0, 0.0),
            ],
            profile: Profile::RhsRect {
                h: 80.0,
                b: 40.0,
                t: 3.0,
            },
            orientation_angle_deg: 0.0,
        };
        let s = to_solid(&m).unwrap();
        match s {
            Solid::Mesh(mesh) => {
                // 3 rings * 4 verts + 2 cap centroids = 14 nodes.
                assert_eq!(mesh.nodes.len(), 14);
            }
            _ => panic!(),
        }
    }

    // --- Phase 38.5 NurbsCurve-path tests ---

    /// A degree-2 NURBS arc-ish curve (3 control points).
    fn sample_nurbs() -> valenx_surface::NurbsCurve {
        valenx_surface::NurbsCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(500.0, 800.0, 0.0),
                Vector3::new(1000.0, 0.0, 0.0),
            ],
            vec![1.0, 1.0, 1.0],
        )
        .expect("valid NURBS curve")
    }

    #[test]
    fn from_nurbs_curve_samples_the_path() {
        let m = Member::from_nurbs_curve(&sample_nurbs(), 20, Profile::default_ipe200())
            .expect("curved member builds");
        assert_eq!(m.path.len(), 20, "path should have n_samples vertices");
        // Curve endpoints are the first and last control points.
        assert!((m.path[0] - Vector3::new(0.0, 0.0, 0.0)).norm() < 1e-6);
        assert!(
            (m.path[19] - Vector3::new(1000.0, 0.0, 0.0)).norm() < 1e-6
        );
        // The curve bulges toward +Y — at least one interior sample
        // must have y > 0 (a straight line could not).
        assert!(m.path.iter().any(|p| p.y > 1.0), "path should be curved");
    }

    #[test]
    fn from_nurbs_curve_member_sweeps_to_a_solid() {
        let m = Member::from_nurbs_curve(&sample_nurbs(), 16, Profile::default_ipe200())
            .unwrap();
        let s = to_solid(&m).expect("curved member meshes");
        match s {
            Solid::Mesh(mesh) => {
                assert!(!mesh.nodes.is_empty());
                assert!(mesh.total_elements() > 0);
            }
            _ => panic!("expected mesh-backed solid"),
        }
    }

    #[test]
    fn from_nurbs_curve_rejects_too_few_samples() {
        let err = Member::from_nurbs_curve(&sample_nurbs(), 1, Profile::default_ipe200())
            .unwrap_err();
        assert!(matches!(err, FramesError::BadParameter { .. }));
    }
}
