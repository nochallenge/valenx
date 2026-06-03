//! Wall entity — a vertical, prismatic segment between two world-space
//! points with explicit height and thickness.
//!
//! ## Geometry
//!
//! The wall is treated as a rectangular box of `length × thickness ×
//! height` where:
//! - `length = ‖end - start‖`
//! - the long axis lies along the `start → end` direction in the XY
//!   plane (z component of `end - start` is ignored for direction —
//!   wall is assumed vertical),
//! - the thickness extends symmetrically perpendicular to the long
//!   axis in the XY plane,
//! - height extends along +Z from the start point.
//!
//! Tessellation builds the box as 6 quad faces → 12 triangles (the
//! same topology a [`valenx_cad::box_solid`] would produce, but we
//! build the mesh directly to keep the wall agnostic to the truck
//! kernel's coordinate frame).

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::error::ArchError;

/// Parameters describing a wall segment.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct WallParams {
    /// Wall start point (one end of the long axis), in world units.
    pub start: Vector3<f64>,
    /// Wall end point (other end of the long axis), in world units.
    pub end: Vector3<f64>,
    /// Wall height (along +Z), in world units.
    pub height: f64,
    /// Wall thickness (perpendicular to long axis, in XY), in world
    /// units. Distributed symmetrically (half on each side of the
    /// centreline).
    pub thickness: f64,
    /// Material descriptor (free-form string — e.g. `"Brick"`,
    /// `"Concrete"`, `"CLT"`). Used by [`crate::Schedule`] grouping
    /// and the IFC writer's material name.
    pub material: String,
}

impl WallParams {
    /// Validate that height + thickness are strictly positive and
    /// that the long-axis length is non-zero.
    pub fn validate(&self) -> Result<(), ArchError> {
        if !self.height.is_finite() || self.height <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "height",
                reason: format!("must be > 0 (got {})", self.height),
            });
        }
        if !self.thickness.is_finite() || self.thickness <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "thickness",
                reason: format!("must be > 0 (got {})", self.thickness),
            });
        }
        if self.length() < 1e-9 {
            return Err(ArchError::BadDimension {
                name: "length",
                reason: format!("start and end coincide (length {})", self.length()),
            });
        }
        Ok(())
    }

    /// Wall length along the start→end axis (XY projection if `end - start`
    /// has a Z component).
    pub fn length(&self) -> f64 {
        let d = self.end - self.start;
        (d.x * d.x + d.y * d.y).sqrt()
    }

    /// Unit vector along the wall's long axis (XY plane).
    /// Returns `(1, 0, 0)` for degenerate inputs.
    pub fn axis_xy(&self) -> Vector3<f64> {
        let d = self.end - self.start;
        let l = (d.x * d.x + d.y * d.y).sqrt();
        if l < 1e-12 {
            Vector3::new(1.0, 0.0, 0.0)
        } else {
            Vector3::new(d.x / l, d.y / l, 0.0)
        }
    }

    /// Perpendicular-in-XY unit vector (90° CCW from the long axis).
    pub fn perp_xy(&self) -> Vector3<f64> {
        let a = self.axis_xy();
        Vector3::new(-a.y, a.x, 0.0)
    }

    /// Compute the 8 corner positions of the wall box in world space.
    ///
    /// Ordering (used by [`Self::tessellate_mesh`]):
    /// `[bl_bot, br_bot, tr_bot, tl_bot, bl_top, br_top, tr_top, tl_top]`
    /// where `bl` is the start side / perpendicular-minus side, `br`
    /// is start side / perpendicular-plus, etc.
    pub fn corners(&self) -> [Vector3<f64>; 8] {
        let a = self.axis_xy();
        let p = self.perp_xy();
        let half_t = self.thickness * 0.5;
        let bl_bot = self.start - p * half_t;
        let br_bot = self.start + p * half_t;
        let tl_bot = self.start + a * self.length() - p * half_t;
        let tr_bot = self.start + a * self.length() + p * half_t;
        let up = Vector3::new(0.0, 0.0, self.height);
        [
            bl_bot,
            br_bot,
            tr_bot,
            tl_bot,
            bl_bot + up,
            br_bot + up,
            tr_bot + up,
            tl_bot + up,
        ]
    }

    /// Tessellate to a [`valenx_mesh::Mesh`] (12 triangles, 8 nodes).
    ///
    /// Faces (CCW outward normals):
    /// - bottom (z-min): 0-1-2-3
    /// - top (z-max): 4-7-6-5
    /// - +perp side (y-plus-ish): 1-5-6-2
    /// - -perp side: 0-3-7-4
    /// - start cap: 0-4-5-1
    /// - end cap: 3-2-6-7
    pub fn tessellate_mesh(&self) -> Result<Mesh, ArchError> {
        self.validate()?;
        let c = self.corners();
        let mut mesh = Mesh::new("wall");
        for v in &c {
            mesh.nodes.push(*v);
        }
        let mut block = ElementBlock::new(ElementType::Tri3);
        let quads: [[u32; 4]; 6] = [
            [0, 1, 2, 3], // bottom
            [4, 7, 6, 5], // top
            [1, 5, 6, 2], // +perp
            [0, 3, 7, 4], // -perp
            [0, 4, 5, 1], // start cap
            [3, 2, 6, 7], // end cap
        ];
        for q in quads {
            block.connectivity.extend_from_slice(&[q[0], q[1], q[2]]);
            block.connectivity.extend_from_slice(&[q[0], q[2], q[3]]);
        }
        mesh.element_blocks.push(block);
        mesh.recompute_stats();
        Ok(mesh)
    }

    /// Return a CAD [`valenx_cad::Solid`] for the wall — a
    /// mesh-backed solid wrapping [`Self::tessellate_mesh`].
    ///
    /// Mesh-backed because a true BRep `box_solid` would need to be
    /// translated + rotated to place the wall, and `valenx_cad`'s
    /// `Solid::translated`/`rotated` paths work but the resulting
    /// truck solid would still need extra surgery for opening cuts
    /// (windows / doors) — we prefer to do all of that in the mesh
    /// domain so opening boolean cuts stay simple.
    pub fn tessellate(&self) -> Result<valenx_cad::Solid, ArchError> {
        let mesh = self.tessellate_mesh()?;
        Ok(valenx_cad::Solid::from_mesh(mesh))
    }

    /// Iterator over the wall's 8 corners — used by
    /// [`crate::ArchDocument::bbox`].
    pub fn bbox_hint_points(&self) -> impl Iterator<Item = Vector3<f64>> + '_ {
        // Materialise into a small Vec so we don't have to plumb the
        // closure ownership all the way back through the iter chain.
        let c = self.corners();
        c.into_iter().collect::<Vec<_>>().into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn wall_3x025() -> WallParams {
        WallParams {
            start: Vector3::zeros(),
            end: Vector3::new(3.0, 0.0, 0.0),
            height: 2.7,
            thickness: 0.25,
            material: "Brick".into(),
        }
    }

    #[test]
    fn length_and_axes() {
        let w = wall_3x025();
        assert!((w.length() - 3.0).abs() < 1e-12);
        let a = w.axis_xy();
        assert!((a - Vector3::new(1.0, 0.0, 0.0)).norm() < 1e-12);
        let p = w.perp_xy();
        assert!((p - Vector3::new(0.0, 1.0, 0.0)).norm() < 1e-12);
    }

    #[test]
    fn rejects_zero_or_negative_dims() {
        let mut w = wall_3x025();
        w.height = -1.0;
        assert!(matches!(
            w.validate(),
            Err(ArchError::BadDimension { name: "height", .. })
        ));

        let mut w = wall_3x025();
        w.thickness = 0.0;
        assert!(matches!(
            w.validate(),
            Err(ArchError::BadDimension {
                name: "thickness",
                ..
            })
        ));

        let mut w = wall_3x025();
        w.end = w.start;
        assert!(matches!(
            w.validate(),
            Err(ArchError::BadDimension { name: "length", .. })
        ));
    }

    #[test]
    fn tessellation_has_12_triangles_and_8_nodes() {
        let w = wall_3x025();
        let m = w.tessellate_mesh().unwrap();
        assert_eq!(m.nodes.len(), 8);
        assert_eq!(m.total_elements(), 12);
    }

    #[test]
    fn tessellate_returns_solid() {
        let w = wall_3x025();
        let s = w.tessellate().unwrap();
        // mesh-backed, so the cached mesh round-trips through
        // `valenx_cad::solid_to_mesh`.
        let m = valenx_cad::solid_to_mesh(&s, 0.5).unwrap();
        assert_eq!(m.total_elements(), 12);
    }
}
