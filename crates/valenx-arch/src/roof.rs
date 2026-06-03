//! Roof entity.
//!
//! v1 ships four roof types — Flat (flat slab), Gable (two pitched
//! planes meeting at a ridge), Hip (four pitched planes meeting at a
//! ridge or point), and Shed (single pitched plane). Tessellation is
//! deliberately simple: each roof type emits the topologically
//! simplest mesh covering the boundary, treating the boundary as an
//! axis-aligned rectangle when it's not — Phase 15.5 will swap in a
//! proper roof solver.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::error::ArchError;

/// Roof shape.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum RoofType {
    /// Flat roof — the boundary is extruded slightly along +Z to give
    /// the roof some thickness.
    Flat,
    /// Gable roof — two sloped planes meeting at a ridge that runs
    /// parallel to the boundary's long axis.
    Gable,
    /// Hip roof — four sloped planes meeting at a horizontal ridge
    /// (or a peak for square footprints).
    Hip,
    /// Shed roof — a single sloped plane.
    Shed,
}

impl RoofType {
    /// Human label.
    pub fn label(&self) -> &'static str {
        match self {
            RoofType::Flat => "Flat",
            RoofType::Gable => "Gable",
            RoofType::Hip => "Hip",
            RoofType::Shed => "Shed",
        }
    }
}

/// Parameters describing a roof.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RoofParams {
    /// Closed boundary (footprint) in world space. Z taken from the
    /// first point (treated as eave height).
    pub boundary: Vec<Vector3<f64>>,
    /// Distance from the boundary plane to the highest point (ridge /
    /// peak / top of the shed). Ignored for [`RoofType::Flat`] (use a
    /// small ε internally to give the roof some thickness).
    pub peak_height: f64,
    /// Roof shape.
    pub roof_type: RoofType,
}

impl RoofParams {
    /// Validate dimensions.
    pub fn validate(&self) -> Result<(), ArchError> {
        if self.boundary.len() < 3 {
            return Err(ArchError::BadDimension {
                name: "boundary",
                reason: format!("need at least 3 points (got {})", self.boundary.len()),
            });
        }
        if !self.peak_height.is_finite() {
            return Err(ArchError::BadDimension {
                name: "peak_height",
                reason: format!("must be finite (got {})", self.peak_height),
            });
        }
        if !matches!(self.roof_type, RoofType::Flat) && self.peak_height <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "peak_height",
                reason: format!(
                    "must be > 0 for {:?} (got {})",
                    self.roof_type, self.peak_height
                ),
            });
        }
        Ok(())
    }

    /// Footprint axis-aligned bounding box (min, max, eave-z).
    fn footprint_aabb(&self) -> (f64, f64, f64, f64, f64) {
        let mut xmin = f64::INFINITY;
        let mut ymin = f64::INFINITY;
        let mut xmax = f64::NEG_INFINITY;
        let mut ymax = f64::NEG_INFINITY;
        let z = self.boundary[0].z;
        for p in &self.boundary {
            xmin = xmin.min(p.x);
            ymin = ymin.min(p.y);
            xmax = xmax.max(p.x);
            ymax = ymax.max(p.y);
        }
        (xmin, ymin, xmax, ymax, z)
    }

    /// Tessellate to a [`valenx_mesh::Mesh`].
    pub fn tessellate_mesh(&self) -> Result<Mesh, ArchError> {
        self.validate()?;
        let mut mesh = Mesh::new("roof");
        let mut block = ElementBlock::new(ElementType::Tri3);
        let (xmin, ymin, xmax, ymax, z0) = self.footprint_aabb();
        let cx = (xmin + xmax) * 0.5;
        let cy = (ymin + ymax) * 0.5;
        let zp = z0 + self.peak_height;

        match self.roof_type {
            RoofType::Flat => {
                // Slab-like thin extrusion above the boundary.
                let thickness = self.peak_height.max(0.05);
                // Reuse the slab tessellator by building a temp slab.
                let slab = crate::slab::SlabParams {
                    boundary: self.boundary.clone(),
                    thickness,
                    material: "Roof".into(),
                    structural: None,
                };
                return slab.tessellate_mesh();
            }
            RoofType::Gable => {
                // Ridge runs along the longer axis. For an
                // axis-aligned rectangle: ridge from
                // (xmin, cy, zp) to (xmax, cy, zp) when xmax-xmin >
                // ymax-ymin, else from (cx, ymin, zp) to (cx, ymax, zp).
                let long_x = (xmax - xmin) >= (ymax - ymin);
                let ridge = if long_x {
                    (Vector3::new(xmin, cy, zp), Vector3::new(xmax, cy, zp))
                } else {
                    (Vector3::new(cx, ymin, zp), Vector3::new(cx, ymax, zp))
                };
                // 6 nodes: 4 eaves + 2 ridge ends.
                mesh.nodes.push(Vector3::new(xmin, ymin, z0)); // 0
                mesh.nodes.push(Vector3::new(xmax, ymin, z0)); // 1
                mesh.nodes.push(Vector3::new(xmax, ymax, z0)); // 2
                mesh.nodes.push(Vector3::new(xmin, ymax, z0)); // 3
                mesh.nodes.push(ridge.0); // 4
                mesh.nodes.push(ridge.1); // 5
                if long_x {
                    // South slope (0-1-5-4) and north slope (3-2-5-4 reversed)
                    block.connectivity.extend_from_slice(&[0, 1, 5]);
                    block.connectivity.extend_from_slice(&[0, 5, 4]);
                    block.connectivity.extend_from_slice(&[3, 4, 5]);
                    block.connectivity.extend_from_slice(&[3, 5, 2]);
                    // Gable triangles on the short walls.
                    block.connectivity.extend_from_slice(&[0, 4, 3]);
                    block.connectivity.extend_from_slice(&[1, 2, 5]);
                } else {
                    block.connectivity.extend_from_slice(&[0, 1, 4]);
                    block.connectivity.extend_from_slice(&[1, 5, 4]);
                    block.connectivity.extend_from_slice(&[3, 4, 5]);
                    block.connectivity.extend_from_slice(&[3, 5, 2]);
                    block.connectivity.extend_from_slice(&[0, 4, 3]);
                    block.connectivity.extend_from_slice(&[1, 2, 5]);
                }
            }
            RoofType::Hip => {
                // Four sloped planes meeting at a single peak.
                mesh.nodes.push(Vector3::new(xmin, ymin, z0)); // 0
                mesh.nodes.push(Vector3::new(xmax, ymin, z0)); // 1
                mesh.nodes.push(Vector3::new(xmax, ymax, z0)); // 2
                mesh.nodes.push(Vector3::new(xmin, ymax, z0)); // 3
                mesh.nodes.push(Vector3::new(cx, cy, zp)); // 4
                block.connectivity.extend_from_slice(&[0, 1, 4]);
                block.connectivity.extend_from_slice(&[1, 2, 4]);
                block.connectivity.extend_from_slice(&[2, 3, 4]);
                block.connectivity.extend_from_slice(&[3, 0, 4]);
            }
            RoofType::Shed => {
                // Single plane sloping from low end (xmin) to high end (xmax).
                mesh.nodes.push(Vector3::new(xmin, ymin, z0)); // 0
                mesh.nodes.push(Vector3::new(xmax, ymin, zp)); // 1
                mesh.nodes.push(Vector3::new(xmax, ymax, zp)); // 2
                mesh.nodes.push(Vector3::new(xmin, ymax, z0)); // 3
                block.connectivity.extend_from_slice(&[0, 1, 2]);
                block.connectivity.extend_from_slice(&[0, 2, 3]);
            }
        }
        mesh.element_blocks.push(block);
        mesh.recompute_stats();
        Ok(mesh)
    }

    /// Wrap [`Self::tessellate_mesh`] in a mesh-backed solid.
    pub fn tessellate(&self) -> Result<valenx_cad::Solid, ArchError> {
        let m = self.tessellate_mesh()?;
        Ok(valenx_cad::Solid::from_mesh(m))
    }

    /// Hint points for [`crate::ArchDocument::bbox`].
    pub fn bbox_hint_points(&self) -> impl Iterator<Item = Vector3<f64>> + '_ {
        let mut v = Vec::with_capacity(self.boundary.len() + 1);
        for p in &self.boundary {
            v.push(*p);
        }
        if let Some(first) = self.boundary.first() {
            v.push(Vector3::new(first.x, first.y, first.z + self.peak_height));
        }
        v.into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rect_boundary() -> Vec<Vector3<f64>> {
        vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(8.0, 0.0, 0.0),
            Vector3::new(8.0, 5.0, 0.0),
            Vector3::new(0.0, 5.0, 0.0),
        ]
    }

    #[test]
    fn all_four_types_tessellate() {
        for rt in [
            RoofType::Flat,
            RoofType::Gable,
            RoofType::Hip,
            RoofType::Shed,
        ] {
            let r = RoofParams {
                boundary: rect_boundary(),
                peak_height: 2.0,
                roof_type: rt.clone(),
            };
            let m = r.tessellate_mesh().unwrap();
            assert!(m.total_elements() > 0, "{rt:?} produced empty mesh");
        }
    }

    #[test]
    fn rejects_short_boundary() {
        let r = RoofParams {
            boundary: vec![Vector3::zeros()],
            peak_height: 2.0,
            roof_type: RoofType::Gable,
        };
        assert!(matches!(
            r.validate(),
            Err(ArchError::BadDimension {
                name: "boundary",
                ..
            })
        ));
    }

    #[test]
    fn rejects_zero_peak_for_pitched() {
        let r = RoofParams {
            boundary: rect_boundary(),
            peak_height: 0.0,
            roof_type: RoofType::Gable,
        };
        assert!(matches!(
            r.validate(),
            Err(ArchError::BadDimension {
                name: "peak_height",
                ..
            })
        ));
    }

    #[test]
    fn flat_roof_accepts_zero_peak() {
        let r = RoofParams {
            boundary: rect_boundary(),
            peak_height: 0.0,
            roof_type: RoofType::Flat,
        };
        r.validate().unwrap();
    }
}
