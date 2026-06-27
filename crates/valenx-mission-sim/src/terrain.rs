//! # `terrain` — an in-house **procedural elevation heightfield** for terrain-aware
//! mission planning.
//!
//! A small, dependency-free **2.5-D terrain model**: a rectangular grid of
//! per-cell ground **elevation** in metres ([`HeightGrid`]) plus the analysis it
//! enables. It is the capstone that ties valenx's tactical [`crate::routing`]
//! (M3) and [`crate::los`] (M4) together: routing becomes **slope-aware** (gentle
//! ground is cheap, a steep ridge flank is expensive, and an impassably steep
//! slope is a wall) and line-of-sight becomes **terrain-masked** in 2.5-D (a ridge
//! between the observer and target casts "dead ground" behind it).
//!
//! This is terrain analysis for **movement and sensor planning** only — a
//! defensive / analysis capability ("where can a vehicle go, and what can a sensor
//! position see across this landscape?"). There is no targeting, no adversary
//! model, and no weapon semantics anywhere.
//!
//! ## Model
//!
//! - The grid is `w × h` cells in **row-major** order (`idx = y · w + x`), each an
//!   `f32` **elevation in metres**. `(0, 0)` is the NW corner; the same row-major
//!   convention as [`crate::routing::CostGrid`], so the two grids overlay
//!   cell-for-cell.
//! - [`demo_terrain`] synthesises a deterministic landscape from a handful of
//!   summed 2-D **Gaussian** hills plus a long diagonal **ridge**, on a gently
//!   tilted plane. No data files, no network, fully reproducible — the same
//!   `(w, h)` always yields the same field.
//! - [`HeightGrid::slope_at`] is the gradient magnitude estimated from the four
//!   orthogonal neighbours by central differences (a dimensionless rise-over-run,
//!   using the cell spacing). It is ~`0` on flat ground and large on a ridge flank.
//! - [`CostGrid::from_terrain`](crate::routing::CostGrid::from_terrain) derives a
//!   routing cost field from slope, and [`crate::los::line_of_sight_terrain`]
//!   ray-marches the elevation for 2.5-D dead-ground masking.
//!
//! Everything here is **pure** (no I/O, no clock, no randomness): deterministic in
//! the grid dimensions and query coordinates.
//!
//! ## DEM hook (future)
//!
//! The procedural [`demo_terrain`] is the in-house default (reliable, offline). A
//! real digital-elevation-model (DEM) raster could be loaded later by constructing
//! a [`HeightGrid`] directly from the sampled elevations (`HeightGrid { w, h, elev
//! }`); every analysis here operates on that struct, so no consumer would change.
//! We deliberately do **not** depend on `gdal` (a heavy C binding) for the
//! built-in model.
//!
//! ## Example
//!
//! ```
//! use valenx_mission_sim::terrain::demo_terrain;
//!
//! let t = demo_terrain(64, 48);
//! let (lo, hi) = t.minmax();
//! assert!(hi > lo, "the demo landscape has relief");
//! // Slope is a non-negative rise-over-run.
//! assert!(t.slope_at(10, 10) >= 0.0);
//! ```

/// Horizontal spacing between adjacent grid cells, in metres. The demo landscape
/// is sized in real metres so slopes come out as realistic rise-over-run values;
/// a `30 m` post spacing matches common DEM resolutions (e.g. SRTM-1).
pub const CELL_SPACING_M: f32 = 30.0;

/// A rectangular grid of per-cell ground **elevation** in metres.
///
/// Cells are stored **row-major**: cell `(x, y)` lives at `elev[y * w + x]`, with
/// `x` in `0..w` (column) and `y` in `0..h` (row), and `(0, 0)` the NW corner.
/// Mirrors the layout of [`crate::routing::CostGrid`] so a terrain grid and a cost
/// grid of the same dimensions index identically.
#[derive(Debug, Clone, PartialEq)]
pub struct HeightGrid {
    /// Width in cells (number of columns; `x` ranges `0..w`).
    pub w: usize,
    /// Height in cells (number of rows; `y` ranges `0..h`).
    pub h: usize,
    /// Row-major ground elevation per cell in **metres** (`len == w * h`).
    pub elev: Vec<f32>,
}

impl HeightGrid {
    /// A `w × h` grid with every cell at the same `elevation_m` (a flat plateau).
    /// Useful as a baseline and in tests (slope is exactly `0` everywhere).
    pub fn flat(w: usize, h: usize, elevation_m: f32) -> Self {
        Self {
            w,
            h,
            elev: vec![elevation_m; w.saturating_mul(h)],
        }
    }

    /// Whether `(x, y)` is inside the grid bounds.
    #[inline]
    pub fn in_bounds(&self, x: usize, y: usize) -> bool {
        x < self.w && y < self.h
    }

    /// The ground **elevation** in metres at cell `(x, y)`.
    ///
    /// An out-of-bounds cell **clamps** to the nearest edge cell (so neighbour
    /// queries near a border are well-defined), and returns `0.0` only for a
    /// genuinely empty grid (`w == 0 || h == 0`).
    #[inline]
    pub fn elevation_at(&self, x: usize, y: usize) -> f32 {
        if self.w == 0 || self.h == 0 {
            return 0.0;
        }
        let cx = x.min(self.w - 1);
        let cy = y.min(self.h - 1);
        self.elev[cy * self.w + cx]
    }

    /// The terrain **slope** at cell `(x, y)` as a dimensionless gradient magnitude
    /// (rise over run): `sqrt((dz/dx)² + (dz/dy)²)`.
    ///
    /// Estimated by **central differences** over the four orthogonal neighbours
    /// (each one [`CELL_SPACING_M`] away), with edge cells clamped (see
    /// [`elevation_at`](Self::elevation_at)). On flat ground this is ~`0`; on a
    /// steep ridge flank it is large (e.g. a 45° slope is `1.0`). Always
    /// non-negative and finite for a finite field.
    pub fn slope_at(&self, x: usize, y: usize) -> f32 {
        if self.w == 0 || self.h == 0 {
            return 0.0;
        }
        // Central difference: (z[x+1] - z[x-1]) / (2·spacing). Saturating_sub keeps
        // the left/down neighbour in-bounds; elevation_at clamps the right/up one,
        // so a border cell uses a one-sided-ish estimate without panicking.
        let zx1 = self.elevation_at(x + 1, y);
        let zx0 = self.elevation_at(x.saturating_sub(1), y);
        let zy1 = self.elevation_at(x, y + 1);
        let zy0 = self.elevation_at(x, y.saturating_sub(1));
        let dzdx = (zx1 - zx0) / (2.0 * CELL_SPACING_M);
        let dzdy = (zy1 - zy0) / (2.0 * CELL_SPACING_M);
        (dzdx * dzdx + dzdy * dzdy).sqrt()
    }

    /// The minimum and maximum elevation `(min_m, max_m)` over the whole grid (both
    /// `0.0` for an empty grid). Used to scale the elevation colour ramp in the
    /// map viz and to report the terrain relief in the agent readout.
    pub fn minmax(&self) -> (f32, f32) {
        let mut lo = f32::INFINITY;
        let mut hi = f32::NEG_INFINITY;
        for &z in &self.elev {
            if z < lo {
                lo = z;
            }
            if z > hi {
                hi = z;
            }
        }
        if lo.is_finite() && hi.is_finite() {
            (lo, hi)
        } else {
            (0.0, 0.0)
        }
    }
}

/// A single 2-D **Gaussian hill**: a peak of height `amp` metres centred at grid
/// cell `(cx, cy)` with a `1/e` radius of `sigma` cells. Summing several of these
/// (plus a ridge and a tilt) builds the demo landscape.
fn gaussian_hill(x: f32, y: f32, cx: f32, cy: f32, sigma: f32, amp: f32) -> f32 {
    let dx = x - cx;
    let dy = y - cy;
    let r2 = dx * dx + dy * dy;
    amp * (-r2 / (2.0 * sigma * sigma)).exp()
}

/// Synthesise a **deterministic procedural landscape** as a `w × h` [`HeightGrid`]
/// of elevations in metres.
///
/// The field is the sum of:
/// - a gently **tilted plane** (a low regional slope rising to the NE) so the
///   ground is never perfectly flat;
/// - a few rounded **Gaussian hills** of differing height and width; and
/// - a long, narrow diagonal **ridge** (a Gaussian ridge-line) that a slope-aware
///   route must climb over or skirt, and that masks line-of-sight behind it,
///   **broken by a single low saddle (a mountain pass)** so a cross-ridge route
///   can thread the pass rather than being walled off — the classic terrain
///   routing problem (find the pass), while dead ground persists elsewhere.
///
/// Heights are clamped to be non-negative (sea level `0 m` floor). The same
/// `(w, h)` always produces the same terrain — no randomness, no data files, no
/// network. Dimensions are taken as-is; a `0`-sized grid yields an empty field.
pub fn demo_terrain(w: usize, h: usize) -> HeightGrid {
    let mut elev = vec![0.0f32; w.saturating_mul(h)];
    if w == 0 || h == 0 {
        return HeightGrid { w, h, elev };
    }
    let wf = w as f32;
    let hf = h as f32;
    // Scale hill placement / width to the grid so the landscape looks the same at
    // any resolution.
    let s = (wf.min(hf)).max(1.0);

    for y in 0..h {
        for x in 0..w {
            let xf = x as f32;
            let yf = y as f32;

            // Regional tilt: a gentle plane rising toward the NE corner (a few
            // hundred metres of relief across the whole field).
            let tilt = 120.0 * (xf / wf) + 80.0 * (1.0 - yf / hf);

            // A handful of rounded hills of varying height / spread.
            let mut z = tilt;
            z += gaussian_hill(xf, yf, 0.28 * wf, 0.32 * hf, 0.16 * s, 420.0);
            z += gaussian_hill(xf, yf, 0.68 * wf, 0.28 * hf, 0.11 * s, 300.0);
            z += gaussian_hill(xf, yf, 0.50 * wf, 0.70 * hf, 0.20 * s, 260.0);
            z += gaussian_hill(xf, yf, 0.82 * wf, 0.78 * hf, 0.09 * s, 200.0);

            // A long diagonal RIDGE: distance from the line y = x (NW→SE), made
            // narrow so its flanks are steep — the classic terrain feature a route
            // must cross and that creates dead ground for line-of-sight.
            // Perpendicular distance from the cell to the diagonal, in cells.
            let ridge_dist = ((xf - yf) * std::f32::consts::FRAC_1_SQRT_2).abs();
            let ridge_sigma = 0.045 * s;
            let ridge_peak = 520.0;
            // A single low SADDLE (mountain pass): a Gaussian dip in the ridge
            // height centred partway along the crest, so the ridge is fully
            // traversable at one place. Position it ~40% along the diagonal and
            // make it a few cells wide; it removes most of the ridge height there,
            // opening a gap a slope-aware route can thread.
            let along = 0.5 * (xf + yf); // position along the NW→SE diagonal (cells)
            let pass_center = 0.40 * s;
            let pass_sigma = 0.10 * s;
            let pass_dip = (-(along - pass_center) * (along - pass_center)
                / (2.0 * pass_sigma * pass_sigma))
                .exp(); // 1 at the pass, ~0 elsewhere
            let ridge_here = ridge_peak * (1.0 - 0.92 * pass_dip);
            z += ridge_here * (-ridge_dist * ridge_dist / (2.0 * ridge_sigma * ridge_sigma)).exp();

            elev[y * w + x] = z.max(0.0);
        }
    }
    HeightGrid { w, h, elev }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flat_grid_has_zero_slope_everywhere() {
        let t = HeightGrid::flat(8, 6, 250.0);
        let (lo, hi) = t.minmax();
        assert_eq!((lo, hi), (250.0, 250.0));
        for y in 0..t.h {
            for x in 0..t.w {
                assert!(
                    t.slope_at(x, y) < 1e-6,
                    "flat ground must have ~zero slope at ({x},{y})"
                );
            }
        }
    }

    #[test]
    fn demo_terrain_is_deterministic() {
        let a = demo_terrain(40, 30);
        let b = demo_terrain(40, 30);
        assert_eq!(a, b, "the procedural terrain must be reproducible");
    }

    #[test]
    fn demo_terrain_has_relief_and_is_nonnegative() {
        let t = demo_terrain(64, 48);
        let (lo, hi) = t.minmax();
        assert!(lo >= 0.0, "elevations are clamped to a sea-level floor");
        assert!(
            hi - lo > 200.0,
            "the demo landscape should have real relief (hills + ridge), got {lo}..{hi}"
        );
        assert_eq!(t.elev.len(), 64 * 48);
    }

    #[test]
    fn ridge_flank_is_steeper_than_open_ground() {
        // The diagonal ridge runs along y = x. A cell just OFF the ridge line sits
        // on its steep flank; a cell far from the ridge (and from any hill centre)
        // is comparatively gentle. The flank slope must exceed the open-ground one.
        let t = demo_terrain(80, 80);
        // On the ridge flank: near the diagonal but a few cells to one side.
        let flank = t.slope_at(40, 36);
        // Open ground: a corner region away from the ridge and the hills.
        let open = t.slope_at(4, 74);
        assert!(
            flank > open,
            "ridge-flank slope ({flank}) must exceed open-ground slope ({open})"
        );
        assert!(flank > 0.05, "the ridge flank should be genuinely steep");
    }

    #[test]
    fn elevation_at_clamps_out_of_bounds() {
        let t = demo_terrain(10, 10);
        // Out-of-bounds queries clamp to the nearest edge cell (no panic).
        assert_eq!(t.elevation_at(100, 100), t.elevation_at(9, 9));
        assert_eq!(t.elevation_at(0, 100), t.elevation_at(0, 9));
    }

    #[test]
    fn slope_is_finite_and_nonnegative_across_the_field() {
        let t = demo_terrain(50, 40);
        for y in 0..t.h {
            for x in 0..t.w {
                let s = t.slope_at(x, y);
                assert!(s.is_finite() && s >= 0.0, "slope at ({x},{y}) = {s}");
            }
        }
    }

    #[test]
    fn empty_grid_is_well_defined() {
        let t = demo_terrain(0, 0);
        assert_eq!(t.minmax(), (0.0, 0.0));
        assert_eq!(t.elevation_at(0, 0), 0.0);
        assert_eq!(t.slope_at(0, 0), 0.0);
        assert!(t.elev.is_empty());
    }

    #[test]
    fn field_is_smooth_between_neighbours() {
        // A summed-Gaussian field is continuous: adjacent cells differ by a bounded
        // amount (no discontinuities). Check the max neighbour jump is modest
        // relative to the total relief.
        let t = demo_terrain(60, 60);
        let (lo, hi) = t.minmax();
        let relief = hi - lo;
        let mut max_jump = 0.0f32;
        for y in 0..t.h {
            for x in 0..t.w {
                let z = t.elevation_at(x, y);
                if x + 1 < t.w {
                    max_jump = max_jump.max((z - t.elevation_at(x + 1, y)).abs());
                }
                if y + 1 < t.h {
                    max_jump = max_jump.max((z - t.elevation_at(x, y + 1)).abs());
                }
            }
        }
        assert!(
            max_jump < 0.25 * relief,
            "neighbour jump {max_jump} should be small vs relief {relief} (smooth field)"
        );
    }
}
