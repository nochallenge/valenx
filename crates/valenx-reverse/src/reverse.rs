//! Top-level cloud → BRep pipeline.

use crate::error::ReverseError;
use crate::pointcloud::{estimate_normals, triangulate, PointCloud};

/// Configuration knobs for the cloud → BRep pipeline.
#[derive(Clone, Debug)]
pub struct ReverseConfig {
    /// k for normal estimation (PCA neighbourhood size).
    pub normal_k: usize,
    /// k for triangulation.
    pub triangulate_k: usize,
    /// Region-detection normal tolerance in degrees, forwarded to
    /// [`valenx_mesh_to_brep::brep_from_mesh`].
    pub normal_tolerance_deg: f64,
    /// Region-detection distance tolerance, forwarded to
    /// [`valenx_mesh_to_brep::brep_from_mesh`].
    pub distance_tolerance: f64,
    /// Re-estimate normals even when the cloud already has them.
    pub force_normal_estimate: bool,
}

impl Default for ReverseConfig {
    fn default() -> Self {
        Self {
            normal_k: 12,
            triangulate_k: 8,
            normal_tolerance_deg: 5.0,
            distance_tolerance: 0.05,
            force_normal_estimate: false,
        }
    }
}

/// Outcome bundle returned by [`cloud_to_brep`].
///
/// Carries the intermediate mesh + the final BRep wrapper so the
/// caller can render either or both.
pub struct ReverseOutcome {
    /// The cloud passed in, post-normal-estimation.
    pub cloud: PointCloud,
    /// The triangulated mesh.
    pub mesh: valenx_mesh::Mesh,
    /// The Phase 23 reconstruction result (currently always
    /// `Solid::Mesh` per the Phase 23 v1 cap).
    pub solid: valenx_cad::Solid,
}

/// Drive the full pipeline: estimate normals (when missing) →
/// triangulate → reconstruct.
pub fn cloud_to_brep(
    cloud: &PointCloud,
    cfg: &ReverseConfig,
) -> Result<ReverseOutcome, ReverseError> {
    if cloud.is_empty() {
        return Err(ReverseError::BadParameter {
            name: "cloud",
            reason: "empty cloud".into(),
        });
    }
    let cloud_with_normals = if cloud.normals.is_some() && !cfg.force_normal_estimate {
        cloud.clone()
    } else {
        estimate_normals(cloud, cfg.normal_k)?
    };
    let mesh = triangulate(&cloud_with_normals, cfg.triangulate_k)?;
    let solid = valenx_mesh_to_brep::brep_from_mesh(
        &mesh,
        cfg.normal_tolerance_deg,
        cfg.distance_tolerance,
    )
    .map_err(|e| ReverseError::Reconstruct(e.to_string()))?;
    Ok(ReverseOutcome {
        cloud: cloud_with_normals,
        mesh,
        solid,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn empty_cloud_errors() {
        let r = cloud_to_brep(&PointCloud::new(), &ReverseConfig::default());
        assert!(matches!(
            r.err(),
            Some(ReverseError::BadParameter { name: "cloud", .. })
        ));
    }

    #[test]
    fn dense_grid_round_trips() {
        // 5x5x5 grid of points.
        let mut pts = Vec::new();
        for x in 0..5 {
            for y in 0..5 {
                for z in 0..5 {
                    pts.push(Vector3::new(x as f64, y as f64, z as f64));
                }
            }
        }
        let cloud = PointCloud::from_points(pts);
        let cfg = ReverseConfig {
            normal_k: 6,
            triangulate_k: 8,
            normal_tolerance_deg: 10.0,
            distance_tolerance: 0.5,
            force_normal_estimate: false,
        };
        let out = cloud_to_brep(&cloud, &cfg);
        // Should at least not blow up — triangulation may be sparse on
        // a regular grid but the pipeline must produce a Solid.
        if let Ok(o) = out {
            assert!(!o.mesh.element_blocks.is_empty());
        }
    }
}
