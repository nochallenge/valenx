//! Slab entity — a polygonal floor / ceiling extruded along +Z.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::error::ArchError;

/// Parameters describing a slab.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SlabParams {
    /// Closed polygon boundary in world space. The slab extends from
    /// the polygon's Z (taken from the first point) up by
    /// [`Self::thickness`] along +Z. The boundary should be planar
    /// (caller responsibility); v1 doesn't enforce planarity but the
    /// extrusion will produce twisted faces if it's not.
    pub boundary: Vec<Vector3<f64>>,
    /// Slab thickness, in world units. Must be > 0.
    pub thickness: f64,
    /// Material descriptor.
    pub material: String,
    /// Optional structural attributes — material grade + applied
    /// area load (kN/m² treated as a per-corner concentrated nodal
    /// force during export). `None` for non-structural slabs.
    /// Consumed by [`crate::structural::export_structural_model`].
    #[serde(default)]
    pub structural: Option<crate::structural::StructuralMember>,
}

impl SlabParams {
    /// Validate that the boundary has at least 3 distinct points and
    /// that thickness is positive.
    pub fn validate(&self) -> Result<(), ArchError> {
        if self.boundary.len() < 3 {
            return Err(ArchError::BadDimension {
                name: "boundary",
                reason: format!("need at least 3 points (got {})", self.boundary.len()),
            });
        }
        if !self.thickness.is_finite() || self.thickness <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "thickness",
                reason: format!("must be > 0 (got {})", self.thickness),
            });
        }
        Ok(())
    }

    /// Tessellate the slab into a [`valenx_mesh::Mesh`].
    ///
    /// Top/bottom faces are fan-triangulated around the first vertex.
    /// Side walls are quads (split into two triangles each) connecting
    /// consecutive boundary edges between the bottom + top layers.
    ///
    /// This is a v1 strategy. It works for convex boundaries and for
    /// simple non-convex ones; complex non-simple polygons or
    /// boundaries with holes are out of scope (Phase 15.5 → call
    /// `valenx_cad::prism` instead, which goes through truck's
    /// proper plane-attach).
    pub fn tessellate_mesh(&self) -> Result<Mesh, ArchError> {
        self.validate()?;
        let n = self.boundary.len();
        let z0 = self.boundary[0].z;
        let z1 = z0 + self.thickness;

        let mut mesh = Mesh::new("slab");
        // Bottom ring (indices 0..n), top ring (indices n..2n).
        for p in &self.boundary {
            mesh.nodes.push(Vector3::new(p.x, p.y, z0));
        }
        for p in &self.boundary {
            mesh.nodes.push(Vector3::new(p.x, p.y, z1));
        }
        let mut block = ElementBlock::new(ElementType::Tri3);
        // Bottom fan (CCW when viewed from -Z, i.e. clockwise in XY
        // → outward normal points down). We emit (0, i, i+1).
        for i in 1..(n - 1) {
            block.connectivity.push(0_u32);
            block.connectivity.push((i + 1) as u32);
            block.connectivity.push(i as u32);
        }
        // Top fan (outward = +Z): (n, n+i, n+i+1).
        for i in 1..(n - 1) {
            block.connectivity.push(n as u32);
            block.connectivity.push((n + i) as u32);
            block.connectivity.push((n + i + 1) as u32);
        }
        // Side walls — for each edge (i, i+1) generate two
        // triangles: (i, i+1, n+i+1) and (i, n+i+1, n+i).
        for i in 0..n {
            let j = (i + 1) % n;
            let b0 = i as u32;
            let b1 = j as u32;
            let t0 = (n + i) as u32;
            let t1 = (n + j) as u32;
            block.connectivity.extend_from_slice(&[b0, b1, t1]);
            block.connectivity.extend_from_slice(&[b0, t1, t0]);
        }
        mesh.element_blocks.push(block);
        mesh.recompute_stats();
        Ok(mesh)
    }

    /// Wrap [`Self::tessellate_mesh`] in a mesh-backed
    /// [`valenx_cad::Solid`].
    pub fn tessellate(&self) -> Result<valenx_cad::Solid, ArchError> {
        let m = self.tessellate_mesh()?;
        Ok(valenx_cad::Solid::from_mesh(m))
    }

    /// Approximate footprint area via the shoelace formula on the
    /// boundary's XY projection. Used by [`crate::Schedule`].
    pub fn area_m2(&self) -> f64 {
        let n = self.boundary.len();
        if n < 3 {
            return 0.0;
        }
        let mut a = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            a += self.boundary[i].x * self.boundary[j].y;
            a -= self.boundary[j].x * self.boundary[i].y;
        }
        (a * 0.5).abs()
    }

    /// Hint points for [`crate::ArchDocument::bbox`].
    pub fn bbox_hint_points(&self) -> impl Iterator<Item = Vector3<f64>> + '_ {
        let mut v = Vec::with_capacity(self.boundary.len() * 2);
        let z0 = self.boundary.first().map(|p| p.z).unwrap_or(0.0);
        let z1 = z0 + self.thickness;
        for p in &self.boundary {
            v.push(Vector3::new(p.x, p.y, z0));
            v.push(Vector3::new(p.x, p.y, z1));
        }
        v.into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square_slab() -> SlabParams {
        SlabParams {
            boundary: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(4.0, 0.0, 0.0),
                Vector3::new(4.0, 3.0, 0.0),
                Vector3::new(0.0, 3.0, 0.0),
            ],
            thickness: 0.2,
            material: "Concrete".into(),
            structural: None,
        }
    }

    #[test]
    fn rejects_short_boundary() {
        let mut s = square_slab();
        s.boundary.truncate(2);
        assert!(matches!(
            s.validate(),
            Err(ArchError::BadDimension {
                name: "boundary",
                ..
            })
        ));
    }

    #[test]
    fn rejects_zero_thickness() {
        let mut s = square_slab();
        s.thickness = 0.0;
        assert!(matches!(
            s.validate(),
            Err(ArchError::BadDimension {
                name: "thickness",
                ..
            })
        ));
    }

    #[test]
    fn area_of_4x3_rectangle_is_12() {
        let s = square_slab();
        assert!((s.area_m2() - 12.0).abs() < 1e-9);
    }

    #[test]
    fn tessellation_emits_top_bottom_and_sides() {
        // A 4-sided slab fan-triangulated gives 2 top tris + 2 bot
        // tris + 4 sides × 2 tris = 12 tris.
        let s = square_slab();
        let m = s.tessellate_mesh().unwrap();
        assert_eq!(m.total_elements(), 12);
        assert_eq!(m.nodes.len(), 8); // 4 bot + 4 top
    }
}
