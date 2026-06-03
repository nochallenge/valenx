//! Space entity — a named room enclosing a floor footprint up to a
//! ceiling height. Spaces are zero-thickness logical volumes — the
//! IFC writer maps them to `IfcSpace`, the schedule reports area +
//! volume, and the tessellator emits a wireframe-ish thin box (top +
//! bottom slabs only) so the user can see the space outline in the
//! viewport without it occluding adjacent geometry.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::error::ArchError;

/// Parameters describing a named space.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SpaceParams {
    /// Closed boundary in world space. Z taken from the first point.
    pub boundary: Vec<Vector3<f64>>,
    /// Ceiling height (along +Z from the boundary z).
    pub ceiling_height: f64,
    /// Room / space name (`"Living"`, `"Kitchen"`, …).
    pub space_name: String,
}

impl SpaceParams {
    /// Validate dimensions.
    pub fn validate(&self) -> Result<(), ArchError> {
        if self.boundary.len() < 3 {
            return Err(ArchError::BadDimension {
                name: "boundary",
                reason: format!("need at least 3 points (got {})", self.boundary.len()),
            });
        }
        if !self.ceiling_height.is_finite() || self.ceiling_height <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "ceiling_height",
                reason: format!("must be > 0 (got {})", self.ceiling_height),
            });
        }
        Ok(())
    }

    /// Floor area via the shoelace formula.
    pub fn floor_area(&self) -> f64 {
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

    /// Volume (floor area × ceiling height).
    pub fn volume(&self) -> f64 {
        self.floor_area() * self.ceiling_height
    }

    /// Tessellate to a [`valenx_mesh::Mesh`] — top + bottom fans only.
    pub fn tessellate_mesh(&self) -> Result<Mesh, ArchError> {
        self.validate()?;
        let n = self.boundary.len();
        let z0 = self.boundary[0].z;
        let z1 = z0 + self.ceiling_height;
        let mut mesh = Mesh::new("space");
        for p in &self.boundary {
            mesh.nodes.push(Vector3::new(p.x, p.y, z0));
        }
        for p in &self.boundary {
            mesh.nodes.push(Vector3::new(p.x, p.y, z1));
        }
        let mut block = ElementBlock::new(ElementType::Tri3);
        for i in 1..(n - 1) {
            block
                .connectivity
                .extend_from_slice(&[0_u32, (i + 1) as u32, i as u32]);
        }
        for i in 1..(n - 1) {
            block
                .connectivity
                .extend_from_slice(&[n as u32, (n + i) as u32, (n + i + 1) as u32]);
        }
        mesh.element_blocks.push(block);
        mesh.recompute_stats();
        Ok(mesh)
    }

    /// Wrap in a mesh-backed solid.
    pub fn tessellate(&self) -> Result<valenx_cad::Solid, ArchError> {
        let m = self.tessellate_mesh()?;
        Ok(valenx_cad::Solid::from_mesh(m))
    }

    /// Hint points for [`crate::ArchDocument::bbox`].
    pub fn bbox_hint_points(&self) -> impl Iterator<Item = Vector3<f64>> + '_ {
        let mut v = Vec::with_capacity(self.boundary.len() * 2);
        let z0 = self.boundary.first().map(|p| p.z).unwrap_or(0.0);
        let z1 = z0 + self.ceiling_height;
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

    fn living_room() -> SpaceParams {
        SpaceParams {
            boundary: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(5.0, 0.0, 0.0),
                Vector3::new(5.0, 4.0, 0.0),
                Vector3::new(0.0, 4.0, 0.0),
            ],
            ceiling_height: 2.7,
            space_name: "Living".into(),
        }
    }

    #[test]
    fn area_and_volume() {
        let s = living_room();
        assert!((s.floor_area() - 20.0).abs() < 1e-9);
        assert!((s.volume() - 54.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_bad_dims() {
        let mut s = living_room();
        s.ceiling_height = 0.0;
        assert!(matches!(
            s.validate(),
            Err(ArchError::BadDimension {
                name: "ceiling_height",
                ..
            })
        ));
    }

    #[test]
    fn tessellation_top_and_bottom_fans() {
        let s = living_room();
        let m = s.tessellate_mesh().unwrap();
        // 4-gon × 2 (top + bot) = 4 triangles.
        assert_eq!(m.total_elements(), 4);
        assert_eq!(m.nodes.len(), 8);
    }
}
