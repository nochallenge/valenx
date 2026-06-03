//! Local grid-refinement guidance — near the body and in the wake.
//!
//! The core solver runs on a *uniform* Cartesian grid. The accuracy
//! that grid delivers is set by its cell size: the body's shear layer
//! and the wake need fine cells, but the far-field does not. A true
//! adaptive-mesh-refinement (AMR) octree solver re-grids on the fly;
//! that is a major subsystem and is **not** what this v1 ships.
//!
//! What this module does ship is the *guidance* an AMR pass would
//! use, and a practical substitute: it identifies the cells that need
//! refinement (the cut-cell band around the body and the wake region
//! behind it), reports a [`RefinementPlan`], and computes the
//! **recommended uniform cell size** to resolve the body at a target
//! cell count — which a caller feeds straight back into
//! [`crate::TunnelSizing`]. In other words: rather than a non-uniform
//! grid, the v1 sizes the *uniform* grid correctly for the body, and
//! tells you where a future AMR pass should add cells.
//!
//! # Honest scope
//!
//! Refinement *guidance* and uniform-grid sizing, not an adaptive
//! solver. A genuine octree AMR with hanging-node flux matching is a
//! documented major extension. The plan this module returns is real
//! and useful — it is exactly the cell set an AMR pass would target —
//! but acting on it (a non-uniform grid) is future work.

use crate::domain::WindTunnel;
use crate::immersed::CellTag;

/// A region flagged for refinement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RefinementZone {
    /// The band of cut cells wrapped around the body — the shear layer
    /// lives here and needs the finest cells.
    NearBody,
    /// The wake region downstream of the body — vortex shedding and
    /// the velocity deficit need resolution here.
    Wake,
}

/// A refinement plan — which cells a future adaptive pass should
/// refine, and the recommended uniform cell size for the current run.
#[derive(Clone, Debug)]
pub struct RefinementPlan {
    /// Linear cell indices flagged near the body.
    pub near_body_cells: Vec<usize>,
    /// Linear cell indices flagged in the wake.
    pub wake_cells: Vec<usize>,
    /// The current uniform cell size (m).
    pub current_cell_size: f64,
    /// The recommended uniform cell size (m) to resolve the body at
    /// the target cell count — feed this back into the tunnel sizing.
    pub recommended_cell_size: f64,
}

impl RefinementPlan {
    /// Total cell count flagged for refinement.
    pub fn flagged_count(&self) -> usize {
        self.near_body_cells.len() + self.wake_cells.len()
    }

    /// True if the recommended cell size is meaningfully finer than
    /// the current one — i.e. the run would benefit from a finer grid.
    pub fn needs_finer_grid(&self) -> bool {
        self.recommended_cell_size < 0.85 * self.current_cell_size
    }
}

/// Build a refinement plan for a wind-tunnel case.
///
/// `wake_length` is how far downstream (in body lengths) the wake
/// region extends; `target_cells_across_body` is the cell count the
/// body's smallest dimension should be resolved with. The returned
/// plan flags the near-body cut-cell band and the wake cells, and
/// reports the cell size that would hit the target resolution.
pub fn plan_refinement(
    tunnel: &WindTunnel,
    wake_length: f64,
    target_cells_across_body: usize,
) -> RefinementPlan {
    let grid = tunnel.grid;
    let h = grid.dx().max(grid.dy()).max(grid.dz());

    // Body bounding box from the solid cells.
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let mut min_z = f64::INFINITY;
    let mut max_z = f64::NEG_INFINITY;
    let mut any_solid = false;
    for k in 0..grid.nz {
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                if tunnel.body.is_solid(i, j, k) {
                    let (cx, cy, cz) = grid.cell_centre(i, j, k);
                    min_x = min_x.min(cx);
                    max_x = max_x.max(cx);
                    min_y = min_y.min(cy);
                    max_y = max_y.max(cy);
                    min_z = min_z.min(cz);
                    max_z = max_z.max(cz);
                    any_solid = true;
                }
            }
        }
    }

    let mut near_body_cells = Vec::new();
    let mut wake_cells = Vec::new();

    if any_solid {
        let body_len = (max_x - min_x).max(max_y - min_y).max(max_z - min_z);
        let wake_end = max_x + wake_length * body_len.max(h);
        for k in 0..grid.nz {
            for j in 0..grid.ny {
                for i in 0..grid.nx {
                    let idx = i + grid.nx * (j + grid.ny * k);
                    let tag = tunnel.body.tag(i, j, k);
                    if tag == CellTag::Cut {
                        near_body_cells.push(idx);
                        continue;
                    }
                    if tag == CellTag::Fluid {
                        let (cx, cy, cz) = grid.cell_centre(i, j, k);
                        // Wake: downstream of the body, within its
                        // lateral footprint (slightly widened).
                        let pad = 0.5 * body_len.max(h);
                        if cx > max_x
                            && cx < wake_end
                            && cy > min_y - pad
                            && cy < max_y + pad
                            && cz > min_z - pad
                            && cz < max_z + pad
                        {
                            wake_cells.push(idx);
                        }
                    }
                }
            }
        }
    }

    // Recommended cell size: resolve the body's smallest dimension
    // with the target cell count.
    let recommended_cell_size = if any_solid {
        let smallest = (max_x - min_x).min(max_y - min_y).min(max_z - min_z);
        (smallest / target_cells_across_body.max(2) as f64).max(1e-9)
    } else {
        h
    };

    RefinementPlan {
        near_body_cells,
        wake_cells,
        current_cell_size: h,
        recommended_cell_size,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::box_body;
    use crate::wind::Wind;
    use nalgebra::Vector3;

    #[test]
    fn plan_flags_the_cut_band_and_the_wake() {
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let tunnel = WindTunnel::build(&body, Wind::straight(20.0).unwrap()).unwrap();
        let plan = plan_refinement(&tunnel, 5.0, 24);
        // The body must produce a near-body cut band.
        assert!(
            !plan.near_body_cells.is_empty(),
            "the body should flag near-body cells"
        );
        // And a wake region behind it.
        assert!(!plan.wake_cells.is_empty(), "the body should flag a wake");
        assert_eq!(
            plan.flagged_count(),
            plan.near_body_cells.len() + plan.wake_cells.len()
        );
    }

    #[test]
    fn near_body_cells_are_exactly_the_cut_cells() {
        let body = box_body(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let tunnel = WindTunnel::build(&body, Wind::straight(20.0).unwrap()).unwrap();
        let plan = plan_refinement(&tunnel, 4.0, 16);
        let (_, cut, _) = tunnel.body.tag_counts();
        assert_eq!(
            plan.near_body_cells.len(),
            cut,
            "near-body set should be exactly the cut cells"
        );
    }

    #[test]
    fn recommended_cell_size_targets_the_resolution() {
        // A box whose smallest dimension is 1 m, target 20 cells →
        // recommended cell size ≈ 0.05 m.
        let body = box_body(Vector3::zeros(), Vector3::new(4.0, 2.0, 1.0));
        let tunnel = WindTunnel::build(&body, Wind::straight(20.0).unwrap()).unwrap();
        let plan = plan_refinement(&tunnel, 5.0, 20);
        assert!(
            plan.recommended_cell_size > 0.02 && plan.recommended_cell_size < 0.1,
            "recommended cell size {} should be ~0.05",
            plan.recommended_cell_size
        );
        assert!(plan.recommended_cell_size.is_finite());
    }

    #[test]
    fn empty_tunnel_flags_nothing() {
        // A tunnel whose body has been stripped flags no cells.
        let body = box_body(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let mut tunnel = WindTunnel::build(&body, Wind::straight(20.0).unwrap()).unwrap();
        for t in tunnel.body.tags.iter_mut() {
            *t = CellTag::Fluid;
        }
        let plan = plan_refinement(&tunnel, 4.0, 16);
        assert_eq!(plan.flagged_count(), 0);
    }
}
