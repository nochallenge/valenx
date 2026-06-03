//! Window entity — hosted on a wall and cut out of that wall during
//! tessellation.

use serde::{Deserialize, Serialize};
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::error::ArchError;
use crate::wall::WallParams;

/// Window style — purely descriptive in v1 (no per-style geometry).
/// Used by [`crate::Schedule`] grouping and the IFC writer.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum WindowStyle {
    /// Side-hinged casement.
    Casement,
    /// Horizontal sliding window.
    Sliding,
    /// Top-hinged awning.
    Awning,
    /// Non-operable fixed.
    Fixed,
}

impl WindowStyle {
    /// Human label.
    pub fn label(self) -> &'static str {
        match self {
            WindowStyle::Casement => "Casement",
            WindowStyle::Sliding => "Sliding",
            WindowStyle::Awning => "Awning",
            WindowStyle::Fixed => "Fixed",
        }
    }
}

/// Parameters describing a window.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct WindowParams {
    /// Id of the host wall (a [`crate::ArchEntity::Wall`] in the same
    /// [`crate::ArchDocument`]).
    pub host: usize,
    /// Position along the wall's long axis, measured from
    /// `wall.start` toward `wall.end` (centre of the window).
    pub position_along_wall: f64,
    /// Height of the window's sill above the wall's bottom edge.
    pub position_height: f64,
    /// Window width (along the wall's long axis).
    pub width: f64,
    /// Window height (along +Z).
    pub height: f64,
    /// Frame thickness. Recorded for IFC + schedule; v1 tessellation
    /// treats the window as a simple void in the wall (no frame
    /// geometry).
    pub frame_thickness: f64,
    /// Style — descriptive (no geometric effect in v1).
    pub style: WindowStyle,
}

impl WindowParams {
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
        if !self.position_height.is_finite() || self.position_height < 0.0 {
            return Err(ArchError::BadDimension {
                name: "position_height",
                reason: format!("must be ≥ 0 (got {})", self.position_height),
            });
        }
        Ok(())
    }

    /// Render the window's *void* — a tiny solid that, when subtracted
    /// from the host wall via [`crate::opening`], leaves the window
    /// opening. The void is sized slightly larger than the wall
    /// thickness so it fully penetrates.
    pub fn void_mesh(&self, host: &WallParams) -> Result<Mesh, ArchError> {
        self.validate()?;
        let axis = host.axis_xy();
        let perp = host.perp_xy();
        let pen_thickness = host.thickness * 1.5; // overpenetrate.
        let mut mesh = Mesh::new("window_void");
        let mut block = ElementBlock::new(ElementType::Tri3);
        let cx = host.start + axis * self.position_along_wall;
        let half_w = self.width * 0.5;
        let z0 = host.start.z + self.position_height;
        let z1 = z0 + self.height;
        let half_t = pen_thickness * 0.5;
        // 8 corners.
        let corners = [
            cx - axis * half_w - perp * half_t,
            cx + axis * half_w - perp * half_t,
            cx + axis * half_w + perp * half_t,
            cx - axis * half_w + perp * half_t,
        ];
        for c in &corners {
            mesh.nodes.push(nalgebra::Vector3::new(c.x, c.y, z0));
        }
        for c in &corners {
            mesh.nodes.push(nalgebra::Vector3::new(c.x, c.y, z1));
        }
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

    /// Render the window's visible geometry — a thin pane occupying
    /// the window's footprint. Caller controls placement; if no host
    /// wall is supplied the pane is centred on the world origin.
    pub fn tessellate_in_wall(&self, host: &WallParams) -> Result<Mesh, ArchError> {
        self.validate()?;
        let axis = host.axis_xy();
        let perp = host.perp_xy();
        let cx = host.start + axis * self.position_along_wall;
        let half_w = self.width * 0.5;
        let z0 = host.start.z + self.position_height;
        let z1 = z0 + self.height;
        // Pane thickness — frame_thickness if positive, else 0.02.
        let pane_t = if self.frame_thickness > 0.0 {
            self.frame_thickness
        } else {
            0.02
        };
        let half_t = pane_t * 0.5;
        let mut mesh = Mesh::new("window_pane");
        let corners = [
            cx - axis * half_w - perp * half_t,
            cx + axis * half_w - perp * half_t,
            cx + axis * half_w + perp * half_t,
            cx - axis * half_w + perp * half_t,
        ];
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
            end: Vector3::new(5.0, 0.0, 0.0),
            height: 2.7,
            thickness: 0.2,
            material: "Brick".into(),
        }
    }

    fn sample_window() -> WindowParams {
        WindowParams {
            host: 1,
            position_along_wall: 2.5,
            position_height: 1.0,
            width: 1.2,
            height: 1.0,
            frame_thickness: 0.05,
            style: WindowStyle::Casement,
        }
    }

    #[test]
    fn rejects_zero_size() {
        let mut w = sample_window();
        w.width = 0.0;
        assert!(matches!(
            w.validate(),
            Err(ArchError::BadDimension { name: "width", .. })
        ));

        let mut w = sample_window();
        w.height = -1.0;
        assert!(matches!(
            w.validate(),
            Err(ArchError::BadDimension { name: "height", .. })
        ));
    }

    #[test]
    fn pane_emits_12_tris() {
        let w = sample_window();
        let mesh = w.tessellate_in_wall(&sample_wall()).unwrap();
        assert_eq!(mesh.total_elements(), 12);
        assert_eq!(mesh.nodes.len(), 8);
    }

    #[test]
    fn void_overpenetrates_thickness() {
        let w = sample_window();
        let v = w.void_mesh(&sample_wall()).unwrap();
        assert_eq!(v.total_elements(), 12);
        assert_eq!(v.nodes.len(), 8);
    }

    #[test]
    fn label_strings_stable() {
        assert_eq!(WindowStyle::Casement.label(), "Casement");
        assert_eq!(WindowStyle::Sliding.label(), "Sliding");
        assert_eq!(WindowStyle::Awning.label(), "Awning");
        assert_eq!(WindowStyle::Fixed.label(), "Fixed");
    }
}
