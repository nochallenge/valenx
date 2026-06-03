//! Door entity — hosted on a wall, opening cut into the wall during
//! tessellation. Mirrors [`crate::WindowParams`] but starts at floor
//! level by convention.

use serde::{Deserialize, Serialize};
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::error::ArchError;
use crate::wall::WallParams;

/// Door style — descriptive.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum DoorStyle {
    /// Single-leaf swing door.
    Single,
    /// Double-leaf swing door.
    Double,
    /// Sliding door.
    Sliding,
    /// Bi-fold door.
    Bifold,
}

impl DoorStyle {
    /// Human label.
    pub fn label(self) -> &'static str {
        match self {
            DoorStyle::Single => "Single",
            DoorStyle::Double => "Double",
            DoorStyle::Sliding => "Sliding",
            DoorStyle::Bifold => "Bifold",
        }
    }
}

/// Hinge side of a swing door.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Side {
    /// Hinge on the left jamb (when looking at the door from
    /// outside).
    Left,
    /// Hinge on the right jamb.
    Right,
}

impl Side {
    /// Human label.
    pub fn label(self) -> &'static str {
        match self {
            Side::Left => "Left",
            Side::Right => "Right",
        }
    }
}

/// Parameters describing a door.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DoorParams {
    /// Id of the host wall.
    pub host: usize,
    /// Position along the wall's long axis (centre of the door).
    pub position_along_wall: f64,
    /// Door width (along wall axis).
    pub width: f64,
    /// Door height (along +Z from the wall's bottom).
    pub height: f64,
    /// Style.
    pub style: DoorStyle,
    /// Hinge side (only meaningful for swing styles).
    pub hinge_side: Side,
}

impl DoorParams {
    /// Validate dimensions.
    pub fn validate(&self) -> Result<(), ArchError> {
        if !self.width.is_finite() || self.width <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "width",
                reason: format!("must be > 0 (got {})", self.width),
            });
        }
        if !self.height.is_finite() || self.height <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "height",
                reason: format!("must be > 0 (got {})", self.height),
            });
        }
        if !self.position_along_wall.is_finite() {
            return Err(ArchError::BadDimension {
                name: "position_along_wall",
                reason: format!("must be finite (got {})", self.position_along_wall),
            });
        }
        Ok(())
    }

    /// Door void — extends from the floor (host wall's bottom Z) up
    /// by `height`. Overpenetrates the wall thickness so the
    /// subtraction is clean.
    pub fn void_mesh(&self, host: &WallParams) -> Result<Mesh, ArchError> {
        self.validate()?;
        let axis = host.axis_xy();
        let perp = host.perp_xy();
        let pen_thickness = host.thickness * 1.5;
        let cx = host.start + axis * self.position_along_wall;
        let half_w = self.width * 0.5;
        let z0 = host.start.z;
        let z1 = z0 + self.height;
        let half_t = pen_thickness * 0.5;
        let corners = [
            cx - axis * half_w - perp * half_t,
            cx + axis * half_w - perp * half_t,
            cx + axis * half_w + perp * half_t,
            cx - axis * half_w + perp * half_t,
        ];
        let mut mesh = Mesh::new("door_void");
        for c in &corners {
            mesh.nodes.push(nalgebra::Vector3::new(c.x, c.y, z0));
        }
        for c in &corners {
            mesh.nodes.push(nalgebra::Vector3::new(c.x, c.y, z1));
        }
        let mut block = ElementBlock::new(ElementType::Tri3);
        let quads: [[u32; 4]; 6] = [
            [0, 1, 2, 3],
            [4, 7, 6, 5],
            [1, 5, 6, 2],
            [0, 3, 7, 4],
            [0, 4, 5, 1],
            [3, 2, 6, 7],
        ];
        for q in quads {
            block.connectivity.extend_from_slice(&[q[0], q[1], q[2]]);
            block.connectivity.extend_from_slice(&[q[0], q[2], q[3]]);
        }
        mesh.element_blocks.push(block);
        mesh.recompute_stats();
        Ok(mesh)
    }

    /// Visible door leaf — thin panel positioned in the wall's plane.
    pub fn tessellate_in_wall(&self, host: &WallParams) -> Result<Mesh, ArchError> {
        self.validate()?;
        let axis = host.axis_xy();
        let perp = host.perp_xy();
        let cx = host.start + axis * self.position_along_wall;
        let half_w = self.width * 0.5;
        let z0 = host.start.z;
        let z1 = z0 + self.height;
        let leaf_t = 0.04;
        let half_t = leaf_t * 0.5;
        let corners = [
            cx - axis * half_w - perp * half_t,
            cx + axis * half_w - perp * half_t,
            cx + axis * half_w + perp * half_t,
            cx - axis * half_w + perp * half_t,
        ];
        let mut mesh = Mesh::new("door_leaf");
        for c in &corners {
            mesh.nodes.push(nalgebra::Vector3::new(c.x, c.y, z0));
        }
        for c in &corners {
            mesh.nodes.push(nalgebra::Vector3::new(c.x, c.y, z1));
        }
        let mut block = ElementBlock::new(ElementType::Tri3);
        let quads: [[u32; 4]; 6] = [
            [0, 1, 2, 3],
            [4, 7, 6, 5],
            [1, 5, 6, 2],
            [0, 3, 7, 4],
            [0, 4, 5, 1],
            [3, 2, 6, 7],
        ];
        for q in quads {
            block.connectivity.extend_from_slice(&[q[0], q[1], q[2]]);
            block.connectivity.extend_from_slice(&[q[0], q[2], q[3]]);
        }
        mesh.element_blocks.push(block);
        mesh.recompute_stats();
        Ok(mesh)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn sample_wall() -> WallParams {
        WallParams {
            start: Vector3::zeros(),
            end: Vector3::new(4.0, 0.0, 0.0),
            height: 2.5,
            thickness: 0.15,
            material: "Brick".into(),
        }
    }

    fn sample_door() -> DoorParams {
        DoorParams {
            host: 1,
            position_along_wall: 2.0,
            width: 0.9,
            height: 2.1,
            style: DoorStyle::Single,
            hinge_side: Side::Left,
        }
    }

    #[test]
    fn rejects_bad_dims() {
        let mut d = sample_door();
        d.width = 0.0;
        assert!(matches!(
            d.validate(),
            Err(ArchError::BadDimension { name: "width", .. })
        ));

        let mut d = sample_door();
        d.height = -1.0;
        assert!(matches!(
            d.validate(),
            Err(ArchError::BadDimension { name: "height", .. })
        ));
    }

    #[test]
    fn leaf_emits_12_tris() {
        let d = sample_door();
        let m = d.tessellate_in_wall(&sample_wall()).unwrap();
        assert_eq!(m.total_elements(), 12);
    }

    #[test]
    fn void_emits_12_tris() {
        let d = sample_door();
        let v = d.void_mesh(&sample_wall()).unwrap();
        assert_eq!(v.total_elements(), 12);
    }

    #[test]
    fn style_and_side_labels() {
        assert_eq!(DoorStyle::Single.label(), "Single");
        assert_eq!(DoorStyle::Double.label(), "Double");
        assert_eq!(Side::Left.label(), "Left");
        assert_eq!(Side::Right.label(), "Right");
    }
}
