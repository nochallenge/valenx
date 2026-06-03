//! Staircase entity — a straight stair with uniform-rise / uniform-run
//! steps.
//!
//! v1 ships straight stairs only. Helical / L-shaped / dog-leg stairs
//! are Phase 15.5+ — the parametric API picks up additional variants
//! (`StairKind::Straight` / `Helical { ... }`).

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::error::ArchError;

/// Parameters describing a stair.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StairParams {
    /// Bottom-of-stair base point (front-left corner of the first
    /// tread), in world space.
    pub base: Vector3<f64>,
    /// Direction the stair travels (horizontal). Must be non-zero;
    /// the implementation normalises and ignores the Z component.
    pub direction: Vector3<f64>,
    /// Total height climbed (top - bottom).
    pub total_rise: f64,
    /// Total horizontal distance covered.
    pub total_run: f64,
    /// Number of equal-size steps.
    pub num_steps: u32,
    /// Stair width perpendicular to [`Self::direction`].
    pub width: f64,
}

impl StairParams {
    /// Validate dimensions.
    pub fn validate(&self) -> Result<(), ArchError> {
        if self.num_steps < 1 {
            return Err(ArchError::BadDimension {
                name: "num_steps",
                reason: format!("must be ≥ 1 (got {})", self.num_steps),
            });
        }
        if !self.total_rise.is_finite() || self.total_rise <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "total_rise",
                reason: format!("must be > 0 (got {})", self.total_rise),
            });
        }
        if !self.total_run.is_finite() || self.total_run <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "total_run",
                reason: format!("must be > 0 (got {})", self.total_run),
            });
        }
        if !self.width.is_finite() || self.width <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "width",
                reason: format!("must be > 0 (got {})", self.width),
            });
        }
        if self.dir_xy().norm() < 1e-9 {
            return Err(ArchError::BadDimension {
                name: "direction",
                reason: "horizontal component is zero".into(),
            });
        }
        Ok(())
    }

    /// Direction vector projected onto XY and normalised.
    pub fn dir_xy(&self) -> Vector3<f64> {
        let d = Vector3::new(self.direction.x, self.direction.y, 0.0);
        let n = d.norm();
        if n < 1e-12 {
            Vector3::new(1.0, 0.0, 0.0)
        } else {
            d / n
        }
    }

    /// Perpendicular (CCW in XY) to the stair direction.
    pub fn perp_xy(&self) -> Vector3<f64> {
        let d = self.dir_xy();
        Vector3::new(-d.y, d.x, 0.0)
    }

    /// Per-step rise (height of each step).
    pub fn rise_per_step(&self) -> f64 {
        self.total_rise / self.num_steps.max(1) as f64
    }

    /// Per-step run (depth of each step).
    pub fn run_per_step(&self) -> f64 {
        self.total_run / self.num_steps.max(1) as f64
    }

    /// Tessellate to a [`valenx_mesh::Mesh`] — each step is a box
    /// with the appropriate rise/run/width.
    pub fn tessellate_mesh(&self) -> Result<Mesh, ArchError> {
        self.validate()?;
        let dir = self.dir_xy();
        let perp = self.perp_xy();
        let rise = self.rise_per_step();
        let run = self.run_per_step();
        let mut mesh = Mesh::new("stair");
        let mut block = ElementBlock::new(ElementType::Tri3);
        for i in 0..self.num_steps as usize {
            // Each step's bounding box: extends from cumulative
            // (run * i, 0) to (total_run, rise * (i+1)).
            let step_base = self.base + dir * (run * i as f64);
            let step_top_z = (i as f64 + 1.0) * rise;
            let step_bot_z = 0.0; // anchored to base.z for now (treat
                                  // each step like a stack of boxes
                                  // grounded on the base).
                                  // Use the start of each step as the bottom of a single
                                  // box that runs from base.z to base.z + step_top_z and
                                  // extends `run × width` horizontally from the step's
                                  // start.
            let z0 = self.base.z + step_bot_z;
            let z1 = self.base.z + step_top_z;
            let half_w = self.width * 0.5;
            let bl_bot = step_base - perp * half_w;
            let br_bot = step_base + perp * half_w;
            let tl_bot = step_base + dir * run - perp * half_w;
            let tr_bot = step_base + dir * run + perp * half_w;
            let n0 = mesh.nodes.len() as u32;
            mesh.nodes.push(Vector3::new(bl_bot.x, bl_bot.y, z0));
            mesh.nodes.push(Vector3::new(br_bot.x, br_bot.y, z0));
            mesh.nodes.push(Vector3::new(tr_bot.x, tr_bot.y, z0));
            mesh.nodes.push(Vector3::new(tl_bot.x, tl_bot.y, z0));
            mesh.nodes.push(Vector3::new(bl_bot.x, bl_bot.y, z1));
            mesh.nodes.push(Vector3::new(br_bot.x, br_bot.y, z1));
            mesh.nodes.push(Vector3::new(tr_bot.x, tr_bot.y, z1));
            mesh.nodes.push(Vector3::new(tl_bot.x, tl_bot.y, z1));
            let quads: [[u32; 4]; 6] = [
                [n0, n0 + 1, n0 + 2, n0 + 3],
                [n0 + 4, n0 + 7, n0 + 6, n0 + 5],
                [n0 + 1, n0 + 5, n0 + 6, n0 + 2],
                [n0, n0 + 3, n0 + 7, n0 + 4],
                [n0, n0 + 4, n0 + 5, n0 + 1],
                [n0 + 3, n0 + 2, n0 + 6, n0 + 7],
            ];
            for q in quads {
                block.connectivity.extend_from_slice(&[q[0], q[1], q[2]]);
                block.connectivity.extend_from_slice(&[q[0], q[2], q[3]]);
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

    /// Hint points for [`crate::ArchDocument::bbox`] — extreme
    /// corners of the stair's bounding box.
    pub fn bbox_hint_points(&self) -> impl Iterator<Item = Vector3<f64>> + '_ {
        let dir = self.dir_xy();
        let perp = self.perp_xy();
        let half_w = self.width * 0.5;
        let bl = self.base - perp * half_w;
        let br = self.base + perp * half_w;
        let top_bl = bl + dir * self.total_run;
        let top_br = br + dir * self.total_run;
        let z0 = self.base.z;
        let z1 = self.base.z + self.total_rise;
        vec![
            Vector3::new(bl.x, bl.y, z0),
            Vector3::new(br.x, br.y, z0),
            Vector3::new(top_bl.x, top_bl.y, z0),
            Vector3::new(top_br.x, top_br.y, z0),
            Vector3::new(bl.x, bl.y, z1),
            Vector3::new(br.x, br.y, z1),
            Vector3::new(top_bl.x, top_bl.y, z1),
            Vector3::new(top_br.x, top_br.y, z1),
        ]
        .into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_stair() -> StairParams {
        StairParams {
            base: Vector3::zeros(),
            direction: Vector3::new(1.0, 0.0, 0.0),
            total_rise: 3.0,
            total_run: 4.0,
            num_steps: 12,
            width: 1.2,
        }
    }

    #[test]
    fn rise_run_per_step() {
        let s = sample_stair();
        assert!((s.rise_per_step() - 0.25).abs() < 1e-9);
        assert!((s.run_per_step() - 4.0 / 12.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_zero_steps_and_bad_dims() {
        let mut s = sample_stair();
        s.num_steps = 0;
        assert!(matches!(
            s.validate(),
            Err(ArchError::BadDimension {
                name: "num_steps",
                ..
            })
        ));

        let mut s = sample_stair();
        s.total_rise = 0.0;
        assert!(matches!(
            s.validate(),
            Err(ArchError::BadDimension {
                name: "total_rise",
                ..
            })
        ));

        let mut s = sample_stair();
        s.width = -1.0;
        assert!(matches!(
            s.validate(),
            Err(ArchError::BadDimension { name: "width", .. })
        ));
    }

    #[test]
    fn tessellation_has_steps_times_12_triangles() {
        let s = sample_stair();
        let m = s.tessellate_mesh().unwrap();
        // 12 steps × 12 tris each = 144 tris.
        assert_eq!(m.total_elements(), 144);
    }
}
