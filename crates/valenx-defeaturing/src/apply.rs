//! Tessellation-based defeature pipeline.

use nalgebra::Vector3;

use valenx_cad::{solid_to_mesh, Solid};
use valenx_mesh::element::ElementType;
use valenx_mesh::Mesh;

use crate::defeature::Defeature;
use crate::error::DefeatureError;

/// Default tessellation tolerance (mm).
pub const DEFAULT_TOLERANCE_MM: f64 = 0.5;

/// Apply a sequence of defeatures to a solid. Returns a
/// [`Solid::Mesh`] with the rejected triangles removed.
pub fn apply(solid: &Solid, defeatures: &[Defeature]) -> Result<Solid, DefeatureError> {
    let mut mesh = match solid {
        Solid::Mesh(m) => m.clone(),
        Solid::Brep(_) => solid_to_mesh(solid, DEFAULT_TOLERANCE_MM)
            .map_err(|e| DefeatureError::Tessellation(e.to_string()))?,
    };

    for d in defeatures {
        apply_one(&mut mesh, *d)?;
    }
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

fn apply_one(mesh: &mut Mesh, d: Defeature) -> Result<(), DefeatureError> {
    match d {
        Defeature::SliverRemove { max_aspect } => {
            if !(0.0..1.0).contains(&max_aspect) {
                return Err(DefeatureError::BadParameter {
                    name: "max_aspect",
                    reason: format!(
                        "must be in [0, 1), got {max_aspect}"
                    ),
                });
            }
            filter_tri_blocks(mesh, |a, b, c| {
                let ab = (b - a).norm();
                let bc = (c - b).norm();
                let ca = (a - c).norm();
                let max = ab.max(bc).max(ca);
                if max < 1e-12 {
                    return false;
                }
                // Sliver metric = minimum altitude / longest edge.
                // The min-edge / max-edge ratio does NOT detect a
                // "cap" sliver — one vertex sitting on the opposite
                // edge gives a near-zero-area triangle whose three
                // edge lengths are still well-proportioned. The
                // smallest altitude (= 2·area / longest edge) collapses
                // to ~0 for exactly that degenerate case, so the area-
                // normalised ratio is the correct thinness measure.
                let area = 0.5 * (b - a).cross(&(c - a)).norm();
                let min_altitude = 2.0 * area / max;
                let aspect = min_altitude / max;
                aspect >= max_aspect
            });
            Ok(())
        }
        Defeature::FilletRemove { max_radius_mm } => {
            if max_radius_mm <= 0.0 {
                return Err(DefeatureError::BadParameter {
                    name: "max_radius_mm",
                    reason: format!(
                        "must be > 0, got {max_radius_mm}"
                    ),
                });
            }
            // v1 proxy: drop triangles whose minimum circumscribed
            // circle radius is below `max_radius_mm` and which have
            // a sharp dihedral with their neighbours (approximated by
            // the per-triangle "blendiness" — small + angular).
            filter_tri_blocks(mesh, |a, b, c| {
                let ab = (b - a).norm();
                let bc = (c - b).norm();
                let ca = (a - c).norm();
                let perim = ab + bc + ca;
                let area = 0.5 * (b - a).cross(&(c - a)).norm();
                if area < 1e-12 {
                    return false;
                }
                let circumradius = ab * bc * ca / (4.0 * area);
                // Keep triangles bigger than the threshold.
                circumradius > max_radius_mm || perim > max_radius_mm * 6.0
            });
            Ok(())
        }
        Defeature::HoleRemove { max_diameter_mm } => {
            if max_diameter_mm <= 0.0 {
                return Err(DefeatureError::BadParameter {
                    name: "max_diameter_mm",
                    reason: format!(
                        "must be > 0, got {max_diameter_mm}"
                    ),
                });
            }
            // v1 proxy: triangles whose centroid sits in a region
            // with a strong inward-pointing radial concavity are
            // candidates. Without a proper feature recognizer we
            // approximate by dropping triangles whose normal is
            // close to -z and whose bounding-box diagonal is small.
            let r = max_diameter_mm * 0.5;
            filter_tri_blocks(mesh, |a, b, c| {
                let n = (b - a).cross(&(c - a)).try_normalize(1e-12);
                let centroid = (a + b + c) / 3.0;
                let max_edge = ((b - a).norm())
                    .max((c - b).norm())
                    .max((a - c).norm());
                if let Some(nn) = n {
                    if nn.z < -0.95 && max_edge < r && centroid.z.abs() < r {
                        return false;
                    }
                }
                true
            });
            Ok(())
        }
        Defeature::TextRemove { max_depth_mm } => {
            if max_depth_mm <= 0.0 {
                return Err(DefeatureError::BadParameter {
                    name: "max_depth_mm",
                    reason: format!("must be > 0, got {max_depth_mm}"),
                });
            }
            // v1 proxy: drop triangles whose surface normal does NOT
            // match the global outward direction AND which sit within
            // `max_depth_mm` of the bounding-box surface. Without a
            // patch detector this collapses to "drop tiny inward-
            // facing triangles" which approximates engraved text
            // removal on simple cases.
            // Compute mesh bbox.
            let mut bb_min = Vector3::repeat(f64::INFINITY);
            let mut bb_max = Vector3::repeat(f64::NEG_INFINITY);
            for p in &mesh.nodes {
                bb_min = bb_min.zip_map(p, f64::min);
                bb_max = bb_max.zip_map(p, f64::max);
            }
            filter_tri_blocks(mesh, |a, b, c| {
                let centroid = (a + b + c) / 3.0;
                let max_edge = ((b - a).norm())
                    .max((c - b).norm())
                    .max((a - c).norm());
                if max_edge < max_depth_mm
                    && (centroid.z - bb_min.z).abs() < max_depth_mm
                {
                    return false;
                }
                true
            });
            Ok(())
        }
    }
}

fn filter_tri_blocks(
    mesh: &mut Mesh,
    mut keep: impl FnMut(Vector3<f64>, Vector3<f64>, Vector3<f64>) -> bool,
) {
    for block in mesh.element_blocks.iter_mut() {
        if !matches!(block.element_type, ElementType::Tri3) {
            continue;
        }
        let mut new = Vec::with_capacity(block.connectivity.len());
        for tri in block.connectivity.chunks(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            if keep(a, b, c) {
                new.extend_from_slice(tri);
            }
        }
        block.connectivity = new;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_mesh::element::ElementBlock;

    fn two_tri_mesh_with_sliver() -> Mesh {
        let mut m = Mesh::new("test");
        m.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(10.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(0.0, 10.0, 0.0));
        // Sliver: nearly co-linear.
        m.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(10.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(5.0, 0.0001, 0.0));
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity.extend_from_slice(&[0, 1, 2, 3, 4, 5]);
        m.element_blocks.push(block);
        m.recompute_stats();
        m
    }

    #[test]
    fn sliver_remove_drops_thin_triangle() {
        let s = Solid::from_mesh(two_tri_mesh_with_sliver());
        let s2 = apply(
            &s,
            &[Defeature::SliverRemove { max_aspect: 0.1 }],
        )
        .unwrap();
        match s2 {
            Solid::Mesh(m) => {
                // One triangle remaining.
                assert_eq!(m.total_elements(), 1);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn sliver_remove_rejects_aspect_outside_range() {
        let s = Solid::from_mesh(two_tri_mesh_with_sliver());
        let err = apply(
            &s,
            &[Defeature::SliverRemove { max_aspect: 1.5 }],
        )
        .unwrap_err();
        assert!(matches!(err, DefeatureError::BadParameter { .. }));
    }

    #[test]
    fn fillet_remove_keeps_large_triangles() {
        let s = Solid::from_mesh(two_tri_mesh_with_sliver());
        let s2 = apply(
            &s,
            &[Defeature::FilletRemove {
                max_radius_mm: 0.5,
            }],
        )
        .unwrap();
        match s2 {
            Solid::Mesh(m) => assert!(m.total_elements() >= 1),
            _ => panic!(),
        }
    }

    #[test]
    fn empty_defeature_list_is_identity() {
        let s = Solid::from_mesh(two_tri_mesh_with_sliver());
        let s2 = apply(&s, &[]).unwrap();
        match s2 {
            Solid::Mesh(m) => assert_eq!(m.total_elements(), 2),
            _ => panic!(),
        }
    }
}
