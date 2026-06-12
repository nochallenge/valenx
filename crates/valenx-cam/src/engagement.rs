//! Engagement-angle analysis — the modern HSM cutter-load proxy.
//!
//! In high-feedrate machining, the **tool engagement angle** — the
//! arc of the cutter circumference actually in contact with the
//! remaining stock — is the dominant tool-load proxy. Mastercam
//! Dynamic Motion, HSMWorks / Fusion 360 Adaptive Clearing, and
//! SolidCAM iMachining all design their toolpaths around bounding
//! it: a small, *constant* engagement lets the spindle hold a much
//! higher chip-load and feedrate without the localised force spikes
//! that snap cutters at corners.
//!
//! This module computes the engagement angle along an existing
//! toolpath by **rasterising the remaining stock into a 2D
//! occupancy grid** and querying how much of each cut move's cutter
//! disc intersects that grid. The grid representation handles
//! arbitrary pocket shapes (square pockets, pockets with islands,
//! pockets with corners that would over-engage on a naive concentric
//! ring) without needing exact polygon CSG.
//!
//! ## What "engagement angle" means here
//!
//! At a cutter centre `c` moving in direction `d`, the engagement
//! angle is the angular extent of the arc on the cutter's outline
//! whose outward-facing material is *solid*. A full-slot cut
//! engages 180° (the leading half-disc is in material); a finishing
//! pass along an open side engages near 0°; a corner over-engages
//! to near 360° before the corner is partially cleared.
//!
//! We sample the angle by walking `N` rays from `c` outward at
//! evenly-spaced angles `θ ∈ [0, 2π)`. For each ray of length
//! `tool.radius`, the grid is sampled at the endpoint; rays that
//! land in **solid stock** count as engaged. The engagement angle is
//! `2π · (engaged / N)` radians.
//!
//! ## v1 scope (honest)
//!
//! - **2D in the XY plane** — the analysis is per Z step-down, just
//!   like the toolpath generator. The cutter is modelled as a circle
//!   of `tool.radius`.
//! - **Binary occupancy** — the cell either is or isn't solid; no
//!   sub-cell coverage. Cell size is `cell_size_mm`, defaults to a
//!   fraction of the tool radius (`tool.radius / 8` is the
//!   commercial-CAM rule of thumb).
//! - **Discrete angular sampling** — the engagement angle is
//!   bucketed to `2π/n_samples`. With the default 64 samples that
//!   resolves to `~5.6°`, which is finer than typical tool-load
//!   plots in Mastercam Verify / HSMWorks ToolPath Verifier.
//! - **No tool-removal-rate optimisation** — engagement is computed
//!   post-hoc, then bounded by the path generator. The path
//!   generator (`adaptive_constant_engagement`) regenerates the path
//!   so the engagement stays under the bound everywhere.
//!
//! See `crate::op::adaptive_constant_engagement` for the toolpath
//! generator that *uses* this analysis to keep the engagement
//! bounded.

use nalgebra::Vector3;

/// Binary occupancy grid in the XY plane — `true` = solid stock.
///
/// The grid origin (cell `(0, 0)`) sits at `min`; cells are
/// `cell_size_mm` square and there are `n_x × n_y` of them. The
/// grid is a Z slice — Z is implicit in the caller's choice of slice.
#[derive(Clone, Debug)]
pub struct StockGrid {
    /// XY origin of cell `(0, 0)`'s minimum corner (mm).
    pub min: nalgebra::Vector2<f64>,
    /// Side length of each square cell (mm).
    pub cell_size_mm: f64,
    /// Number of cells along X.
    pub n_x: usize,
    /// Number of cells along Y.
    pub n_y: usize,
    /// Row-major `n_x · n_y` solidity flags.
    pub solid: Vec<bool>,
}

impl StockGrid {
    /// Allocate a grid covering `[min, max]` at the requested cell
    /// size. Every cell starts solid (the full stock outline is
    /// material).
    pub fn new(
        min: nalgebra::Vector2<f64>,
        max: nalgebra::Vector2<f64>,
        cell_size_mm: f64,
    ) -> Self {
        let span = max - min;
        let n_x = (span.x / cell_size_mm).ceil().max(1.0) as usize;
        let n_y = (span.y / cell_size_mm).ceil().max(1.0) as usize;
        Self {
            min,
            cell_size_mm,
            n_x,
            n_y,
            solid: vec![true; n_x * n_y],
        }
    }

    /// Map an XY point to its cell index. Out-of-grid points
    /// return `None`.
    pub fn cell_of(&self, x: f64, y: f64) -> Option<(usize, usize)> {
        let ix = ((x - self.min.x) / self.cell_size_mm).floor() as isize;
        let iy = ((y - self.min.y) / self.cell_size_mm).floor() as isize;
        if ix < 0 || iy < 0 {
            return None;
        }
        let (ix, iy) = (ix as usize, iy as usize);
        if ix >= self.n_x || iy >= self.n_y {
            return None;
        }
        Some((ix, iy))
    }

    /// True if `(x, y)` lies in a solid cell. Off-grid points are
    /// treated as *non-solid* (we don't have a notion of stock past
    /// the grid extent).
    pub fn is_solid(&self, x: f64, y: f64) -> bool {
        match self.cell_of(x, y) {
            Some((ix, iy)) => self.solid[iy * self.n_x + ix],
            None => false,
        }
    }

    /// Total solid-cell count. Useful for area-removed estimates +
    /// progress reports.
    pub fn solid_count(&self) -> usize {
        self.solid.iter().filter(|b| **b).count()
    }

    /// Carve the disc of radius `r` centred on `c` — every solid
    /// cell whose centre falls inside the disc becomes non-solid.
    ///
    /// Used to "consume" the cutter's footprint along a toolpath
    /// when simulating engagement-angle history.
    pub fn carve_disc(&mut self, c: nalgebra::Vector2<f64>, r: f64) {
        let r2 = r * r;
        let cs = self.cell_size_mm;
        // Compute the AABB of the disc in cell index space.
        let ix_min = (((c.x - r) - self.min.x) / cs).floor() as isize;
        let ix_max = (((c.x + r) - self.min.x) / cs).ceil() as isize;
        let iy_min = (((c.y - r) - self.min.y) / cs).floor() as isize;
        let iy_max = (((c.y + r) - self.min.y) / cs).ceil() as isize;
        let ix_min = ix_min.max(0) as usize;
        let iy_min = iy_min.max(0) as usize;
        let ix_max = (ix_max.max(0) as usize).min(self.n_x);
        let iy_max = (iy_max.max(0) as usize).min(self.n_y);
        for iy in iy_min..iy_max {
            for ix in ix_min..ix_max {
                let cx = self.min.x + (ix as f64 + 0.5) * cs;
                let cy = self.min.y + (iy as f64 + 0.5) * cs;
                let dx = cx - c.x;
                let dy = cy - c.y;
                if dx * dx + dy * dy <= r2 {
                    self.solid[iy * self.n_x + ix] = false;
                }
            }
        }
    }

    /// Carve the swept disc from `a` to `b` (inclusive) at radius
    /// `r`. Walks the segment in steps of `cell_size_mm / 2` and
    /// carves a disc at each step — fine-enough to leave no holes.
    pub fn carve_segment(&mut self, a: nalgebra::Vector2<f64>, b: nalgebra::Vector2<f64>, r: f64) {
        let d = b - a;
        let len = d.norm();
        if len < 1e-12 {
            self.carve_disc(a, r);
            return;
        }
        let step = (self.cell_size_mm * 0.5).min(r * 0.5).max(1e-6);
        let n_steps = (len / step).ceil() as usize;
        for i in 0..=n_steps {
            let t = (i as f64) / (n_steps as f64);
            self.carve_disc(a + d * t, r);
        }
    }
}

/// One engagement-angle sample at a single point on the toolpath.
#[derive(Clone, Copy, Debug)]
pub struct EngagementSample {
    /// Cutter centre when this sample was taken (XY only).
    pub centre: nalgebra::Vector2<f64>,
    /// Engagement angle in radians (`0` = no contact, `π` = half
    /// engagement / full slot, `2π` = totally surrounded — only seen
    /// at sharp interior corners before they're cleared).
    pub engagement_rad: f64,
}

/// Compute the engagement angle of a cutter at `centre` against
/// `grid` — sample `n_samples` outward rays of length `tool_radius`
/// and count how many hit solid stock.
///
/// Returns radians.
pub fn engagement_at(
    centre: nalgebra::Vector2<f64>,
    tool_radius: f64,
    grid: &StockGrid,
    n_samples: usize,
) -> f64 {
    let n = n_samples.max(8);
    let mut hit = 0_usize;
    // Sample slightly inside the cutter outline (0.98 R) so a
    // freshly-carved cell on the cutter outline isn't falsely
    // counted as not-in-contact.
    let r_sample = tool_radius * 0.98;
    for i in 0..n {
        let theta = (i as f64) * std::f64::consts::TAU / (n as f64);
        let x = centre.x + r_sample * theta.cos();
        let y = centre.y + r_sample * theta.sin();
        if grid.is_solid(x, y) {
            hit += 1;
        }
    }
    (hit as f64) * std::f64::consts::TAU / (n as f64)
}

/// Walk an XY toolpath sample-by-sample and report the engagement
/// angle at every cut move's endpoint, **carving** the grid as the
/// cutter moves so each query sees the *remaining* stock.
///
/// `cut_points` is the sequence of cutter-centre XY positions —
/// typically the `Cut` moves projected to XY (the caller has already
/// filtered out rapids + plunges).
pub fn engagement_along(
    cut_points: &[nalgebra::Vector2<f64>],
    tool_radius: f64,
    grid: &mut StockGrid,
    n_samples: usize,
) -> Vec<EngagementSample> {
    let mut out = Vec::with_capacity(cut_points.len());
    if cut_points.is_empty() {
        return out;
    }
    // Initial cutter footprint at the entry point.
    grid.carve_disc(cut_points[0], tool_radius);
    out.push(EngagementSample {
        centre: cut_points[0],
        engagement_rad: 0.0, // entry point — already carved
    });
    for w in cut_points.windows(2) {
        // Sample engagement at the *destination* using the grid
        // state *before* this move's swept volume is removed.
        let centre = w[1];
        let eng = engagement_at(centre, tool_radius, grid, n_samples);
        out.push(EngagementSample {
            centre,
            engagement_rad: eng,
        });
        // Now carve the swept volume so the next sample sees the
        // post-move stock state.
        grid.carve_segment(w[0], w[1], tool_radius);
    }
    out
}

/// Convenience: project a 3D toolpath move sequence to a flat XY
/// trace, keeping only `Cut` moves. Used by the engagement-along
/// pipeline to drop Z and to skip rapids / plunges that don't
/// participate in side engagement.
pub fn extract_cut_xy(moves: &[crate::toolpath::Move]) -> Vec<nalgebra::Vector2<f64>> {
    moves
        .iter()
        .filter(|m| m.kind == crate::toolpath::MoveKind::Cut)
        .map(|m| nalgebra::Vector2::new(m.position.x, m.position.y))
        .collect()
}

/// Convenience: project a single 3D position to XY.
pub fn xy(p: Vector3<f64>) -> nalgebra::Vector2<f64> {
    nalgebra::Vector2::new(p.x, p.y)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector2;

    #[test]
    fn grid_starts_fully_solid() {
        let g = StockGrid::new(Vector2::new(0.0, 0.0), Vector2::new(10.0, 10.0), 1.0);
        assert_eq!(g.solid_count(), g.n_x * g.n_y);
        assert!(g.is_solid(5.0, 5.0));
        assert!(!g.is_solid(-1.0, 5.0));
    }

    #[test]
    fn carve_disc_clears_circle() {
        let mut g = StockGrid::new(Vector2::new(-5.0, -5.0), Vector2::new(5.0, 5.0), 0.5);
        g.carve_disc(Vector2::new(0.0, 0.0), 2.0);
        // Centre is hollow.
        assert!(!g.is_solid(0.0, 0.0));
        assert!(!g.is_solid(1.0, 0.0));
        // Far corner still solid.
        assert!(g.is_solid(4.5, 4.5));
    }

    #[test]
    fn engagement_unobstructed_disc_zero() {
        let g = StockGrid::new(Vector2::new(0.0, 0.0), Vector2::new(10.0, 10.0), 0.5);
        // A cutter at the grid centre with the whole disc surface
        // *over solid* — every outward sample hits material → full
        // engagement (2π).
        let eng = engagement_at(Vector2::new(5.0, 5.0), 1.0, &g, 64);
        assert!((eng - std::f64::consts::TAU).abs() < 1e-6, "got {eng}");
    }

    #[test]
    fn engagement_void_zero() {
        let mut g = StockGrid::new(Vector2::new(0.0, 0.0), Vector2::new(10.0, 10.0), 0.5);
        // Carve a big disc around the test centre so the sample ring
        // sees no solid.
        g.carve_disc(Vector2::new(5.0, 5.0), 3.0);
        let eng = engagement_at(Vector2::new(5.0, 5.0), 1.0, &g, 64);
        assert!(eng < 0.05, "got {eng}");
    }

    #[test]
    fn engagement_half_slot_about_pi() {
        // Carve the LEFT half-plane and check that a cutter sitting
        // on the freshly-cleared boundary engages ~π (the right
        // half is still solid).
        let mut g = StockGrid::new(Vector2::new(-5.0, -5.0), Vector2::new(5.0, 5.0), 0.25);
        for iy in 0..g.n_y {
            for ix in 0..g.n_x {
                let cx = g.min.x + (ix as f64 + 0.5) * g.cell_size_mm;
                if cx < 0.0 {
                    g.solid[iy * g.n_x + ix] = false;
                }
            }
        }
        // Place a 1-mm-radius cutter centred at x=0 (the boundary).
        let eng = engagement_at(Vector2::new(0.0, 0.0), 1.0, &g, 128);
        // Should be near π. Bucket resolution is 2π/128 ≈ 0.05 rad.
        assert!(
            (eng - std::f64::consts::PI).abs() < 0.1,
            "expected ~π, got {eng}"
        );
    }

    #[test]
    fn engagement_along_carves_progressively() {
        // 20x20 grid; cutter walks across it left→right. After the
        // walk completes, every sample after the first should report
        // engagement, and the grid should be partially hollow along
        // the strip.
        let mut g = StockGrid::new(Vector2::new(0.0, 0.0), Vector2::new(20.0, 20.0), 0.5);
        let cut_xy: Vec<_> = (0..20).map(|i| Vector2::new(i as f64, 10.0)).collect();
        let report = engagement_along(&cut_xy, 1.5, &mut g, 64);
        assert_eq!(report.len(), 20);
        // First sample is the entry-carve marker, engagement=0.
        assert!(report[0].engagement_rad < 1e-6);
        // Subsequent samples — cutter is mostly slotting (~π for a
        // straight cut along virgin material on each side). Some
        // engagement, well above 0.
        let mid = report[10].engagement_rad;
        assert!(
            mid > 1.0,
            "mid-walk engagement should be substantial, got {mid}"
        );
    }
}
