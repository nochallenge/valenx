//! [`Animation`] — playback wrapper around the Phase 17 voxel sim.
//!
//! Each frame replays the first `t * total_moves` moves of the
//! toolpath into a fresh voxel grid and meshes the result. Replaying
//! from scratch each frame keeps the API simple (no implicit history),
//! and `to_mesh()` on the voxel grid already returns a complete
//! triangle mesh suitable for the viewport.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::CamoticsError;
use crate::report::{FrameMetadata, MaterialRemovalReport};
use valenx_cam::stock::Stock;
use valenx_cam::toolpath::Toolpath;
use valenx_cam::voxel::Voxel;
use valenx_mesh::Mesh;

/// Animated material-removal simulation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Animation {
    /// Block of starting material.
    pub stock: Stock,
    /// Sequence of moves applied to the stock.
    pub toolpath: Toolpath,
    /// Tool radius in mm (cylindrical bit assumed). Must be > 0.
    pub tool_radius_mm: f64,
    /// Number of frames in the animation. `frame(0)` = stock, `frame(n-1)` = post-cut.
    pub n_frames: usize,
    /// Voxel resolution `(nx, ny, nz)` used for every frame.
    pub voxel_resolution: (u32, u32, u32),
}

impl Animation {
    /// Build a validated [`Animation`].
    pub fn new(
        stock: Stock,
        toolpath: Toolpath,
        tool_radius_mm: f64,
        n_frames: usize,
        voxel_resolution: (u32, u32, u32),
    ) -> Result<Self, CamoticsError> {
        if !tool_radius_mm.is_finite() || tool_radius_mm <= 0.0 {
            return Err(CamoticsError::BadParameter {
                name: "tool_radius_mm",
                reason: format!("must be > 0, got {tool_radius_mm}"),
            });
        }
        if n_frames < 2 {
            return Err(CamoticsError::BadParameter {
                name: "n_frames",
                reason: format!("must be >= 2, got {n_frames}"),
            });
        }
        if voxel_resolution.0 == 0 || voxel_resolution.1 == 0 || voxel_resolution.2 == 0 {
            return Err(CamoticsError::BadParameter {
                name: "voxel_resolution",
                reason: format!("each axis must be > 0, got {voxel_resolution:?}"),
            });
        }
        Ok(Self {
            stock,
            toolpath,
            tool_radius_mm,
            n_frames,
            voxel_resolution,
        })
    }

    /// Mesh of the stock at `t ∈ [0, 1]`. `t=0` = original stock,
    /// `t=1` = fully cut. Re-runs the simulator from scratch — for a
    /// 100-frame playback this is `O(100 * total_moves)` cuts which is
    /// fine for the v1 scope.
    ///
    /// # Errors
    ///
    /// Returns [`CamoticsError::Cam`] if the voxel grid construction
    /// hits the round-12 `MAX_VOXEL_CELLS` cap.
    pub fn frame(&self, t: f64) -> Result<Mesh, CamoticsError> {
        let t = t.clamp(0.0, 1.0);
        let mut grid = self.fresh_grid()?;
        let moves = &self.toolpath.moves;
        if moves.is_empty() {
            return Ok(grid.to_mesh());
        }
        // Number of moves to replay = round(t * (moves.len() - 1))
        let n_moves = ((t * (moves.len() - 1) as f64).round() as usize).min(moves.len() - 1);
        for window in moves[..=n_moves].windows(2) {
            grid.cut_segment(window[0].position, window[1].position, self.tool_radius_mm);
        }
        Ok(grid.to_mesh())
    }

    /// Mesh of the stock at `t ∈ [0, 1]`, extracted with **Surface
    /// Nets** — the smooth-surface counterpart of [`Animation::frame`]
    /// (Phase 56.5).
    ///
    /// [`Animation::frame`] returns the blocky axis-aligned voxel
    /// boundary; this returns the smoothed dual-contour surface
    /// ([`Voxel::to_mesh_surface_nets`]), which rounds off the voxel
    /// stair-stepping and is much closer to the real machined finish
    /// CAMotics renders. Use it for the final-quality viewport /
    /// export; the faceted [`Animation::frame`] stays available for
    /// the cheapest interactive scrub.
    ///
    /// # Errors
    ///
    /// Returns [`CamoticsError::Cam`] if the voxel grid construction
    /// or surface-net extraction hits the round-12
    /// `MAX_VOXEL_CELLS` cap.
    pub fn frame_smooth(&self, t: f64) -> Result<Mesh, CamoticsError> {
        let t = t.clamp(0.0, 1.0);
        let mut grid = self.fresh_grid()?;
        let moves = &self.toolpath.moves;
        if moves.is_empty() {
            return grid
                .to_mesh_surface_nets()
                .map_err(|e| CamoticsError::Cam(e.to_string()));
        }
        let n_moves = ((t * (moves.len() - 1) as f64).round() as usize).min(moves.len() - 1);
        for window in moves[..=n_moves].windows(2) {
            grid.cut_segment(window[0].position, window[1].position, self.tool_radius_mm);
        }
        grid.to_mesh_surface_nets()
            .map_err(|e| CamoticsError::Cam(e.to_string()))
    }

    /// All frames at uniformly-spaced `t` values.
    ///
    /// # Errors
    ///
    /// Returns [`CamoticsError::Cam`] on the first frame that fails
    /// to construct.
    pub fn frames(&self) -> Result<Vec<Mesh>, CamoticsError> {
        (0..self.n_frames)
            .map(|i| {
                let t = i as f64 / (self.n_frames - 1) as f64;
                self.frame(t)
            })
            .collect()
    }

    /// One frame at index `i`. Returns `FrameOutOfRange` if `i >= n_frames`.
    pub fn frame_at(&self, i: usize) -> Result<Mesh, CamoticsError> {
        if i >= self.n_frames {
            return Err(CamoticsError::FrameOutOfRange(i, self.n_frames));
        }
        let t = i as f64 / (self.n_frames - 1) as f64;
        self.frame(t)
    }

    /// Metadata for frame `i` — time, position, MRR, swarf volume.
    pub fn frame_metadata(&self, i: usize) -> Result<FrameMetadata, CamoticsError> {
        if i >= self.n_frames {
            return Err(CamoticsError::FrameOutOfRange(i, self.n_frames));
        }
        let t = i as f64 / (self.n_frames - 1) as f64;
        let moves = &self.toolpath.moves;
        let pos = if moves.is_empty() {
            self.stock.origin
        } else {
            let idx = ((t * (moves.len() - 1) as f64).round() as usize).min(moves.len() - 1);
            moves[idx].position
        };
        // Solid voxel count → cut volume = (initial - now) * cell_size.
        let grid_now = {
            let mut g = self.fresh_grid()?;
            if !moves.is_empty() {
                let n_moves =
                    ((t * (moves.len() - 1) as f64).round() as usize).min(moves.len() - 1);
                for window in moves[..=n_moves].windows(2) {
                    g.cut_segment(window[0].position, window[1].position, self.tool_radius_mm);
                }
            }
            g
        };
        let initial = self.fresh_grid()?.solid_count() as i64;
        let now = grid_now.solid_count() as i64;
        let cs = grid_now.cell_size();
        let cell_vol = cs.x * cs.y * cs.z;
        let removed_vol = (initial - now).max(0) as f64 * cell_vol;
        // MRR (per-frame) = (removed since previous frame) / Δt.
        // Δt is symbolic (1.0 / n_frames) — we don't track real machining time.
        let mrr = if i > 0 {
            let dt = 1.0 / (self.n_frames - 1) as f64;
            let prev_meta = self.frame_metadata(i - 1)?;
            (removed_vol - prev_meta.cut_volume_mm3) / dt
        } else {
            0.0
        };
        Ok(FrameMetadata {
            index: i,
            t,
            position: pos,
            cut_volume_mm3: removed_vol,
            swarf_volume_mm3: removed_vol,
            mrr_mm3_per_unit_t: mrr,
        })
    }

    /// Build a full material-removal report by walking every frame.
    pub fn material_removal_report(&self) -> Result<MaterialRemovalReport, CamoticsError> {
        let mut frames = Vec::with_capacity(self.n_frames);
        for i in 0..self.n_frames {
            frames.push(self.frame_metadata(i)?);
        }
        let total_cut = frames.last().map(|f| f.cut_volume_mm3).unwrap_or(0.0);
        Ok(MaterialRemovalReport {
            frames,
            total_cut_volume_mm3: total_cut,
            tool_deflection_mm: 0.0,
        })
    }

    fn fresh_grid(&self) -> Result<Voxel, CamoticsError> {
        Voxel::from_aabb(
            self.stock.origin,
            self.stock.origin + self.stock.size,
            self.voxel_resolution,
        )
        .map_err(|e| CamoticsError::Cam(e.to_string()))
    }

    /// Side-by-side viewport payload — `(target_mesh, current_mesh)`.
    /// The target is the stock at `t=0`; downstream UI overlays it
    /// alongside the simulated material. v1 returns the stock mesh
    /// for the target since the "ideal target solid" is upstream of
    /// this crate; callers who own a target Solid mesh should swap
    /// the first element themselves.
    ///
    /// # Errors
    ///
    /// Returns [`CamoticsError::Cam`] if either frame's voxel grid
    /// construction fails.
    pub fn side_by_side(&self, t: f64) -> Result<(Mesh, Mesh), CamoticsError> {
        let target = self.frame(0.0)?;
        let current = self.frame(t)?;
        Ok((target, current))
    }

    /// Stock vector helpers — wraps [`Vector3`] arithmetic so callers
    /// don't need to import nalgebra just to inspect a stock corner.
    pub fn stock_origin(&self) -> Vector3<f64> {
        self.stock.origin
    }

    /// Stock extent.
    pub fn stock_size(&self) -> Vector3<f64> {
        self.stock.size
    }
}
