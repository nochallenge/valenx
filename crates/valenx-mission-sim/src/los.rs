//! # `los` — in-house **2-D line-of-sight** over the routing cost grid.
//!
//! A small, dependency-free **supercover ray-march** that answers a single
//! defensive / analysis question: from one grid cell, can an observer *see*
//! another cell, or does the terrain / an obstacle **mask** the view? It reuses
//! the very same [`crate::routing::CostGrid`] the [`crate::routing`] A\*
//! planner runs over — open ground is transparent, an `f32::INFINITY` cell is an
//! opaque obstacle that blocks the sight line.
//!
//! This is sensor / visibility *analysis* only — a defensive capability ("what
//! can this observer position see, and what is dead ground behind a ridge?") with
//! no targeting, no adversary model, and no weapon semantics. It complements
//! [`crate::routing`]'s path planning: where routing finds a way **around**
//! obstacles, this reports what is **visible past** them.
//!
//! ## Model
//!
//! - The grid is the routing [`CostGrid`]: `w × h` row-major cells, each an `f32`
//!   traversal cost. For line-of-sight a cell is **opaque** exactly when it is
//!   **not** [`CostGrid::passable`] — i.e. an `f32::INFINITY` obstacle (or out of
//!   bounds). Every finite, non-negative cell is **transparent** (the terrain is
//!   treated as flat occupancy; there is no per-cell height field yet).
//! - A sight line is the segment between two **cell centres**. We walk every cell
//!   the segment passes through — a **supercover** traversal (unlike Bresenham,
//!   which can skip a corner) so a wall that only clips the corner of the ray
//!   still blocks it. The standard [Amanatides–Woo] voxel DDA drives the walk.
//! - The two **endpoints are excluded** from the occlusion test: an observer
//!   standing on (or looking at) a marginal cell can still "see" out of / into it.
//!   Only an opaque cell strictly **between** the endpoints blocks the view.
//!
//! Everything here is **pure** (no I/O, no clock, no randomness): the same grid
//! and endpoints always yield the same answer.
//!
//! [Amanatides–Woo]: https://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm)
//!
//! ## Example
//!
//! ```
//! use valenx_mission_sim::routing::CostGrid;
//! use valenx_mission_sim::los::line_of_sight;
//!
//! // 5×1 open corridor: both ends see each other.
//! let mut grid = CostGrid::uniform(5, 1, 1.0);
//! assert!(line_of_sight(&grid, (0, 0), (4, 0)));
//!
//! // Drop a wall in the middle and the view is masked.
//! grid.cost[2] = f32::INFINITY;
//! assert!(!line_of_sight(&grid, (0, 0), (4, 0)));
//! ```

use crate::routing::CostGrid;

/// Whether the cell `(x, y)` **blocks** a sight line — i.e. it is opaque. A cell
/// is opaque exactly when it is not passable on the grid: an `f32::INFINITY`
/// obstacle, or anything out of bounds. Mirrors the routing notion of an
/// impassable cell so LoS and routing agree on what an obstacle is.
#[inline]
fn opaque(grid: &CostGrid, x: usize, y: usize) -> bool {
    !grid.passable(x, y)
}

/// Whether the observer at cell `from` has a clear **line of sight** to the cell
/// `to`, across the occupancy of `grid`.
///
/// Returns `true` when no **opaque** cell (an `f32::INFINITY` obstacle) lies
/// strictly between the two cell centres, and `false` when the sight line crosses
/// any obstacle. The traversal is a **supercover** DDA, so a wall that merely
/// clips the corner of the ray still blocks it (no diagonal "leak" between two
/// touching obstacles).
///
/// Endpoint handling:
/// - `from == to` is always `true` (a cell sees itself).
/// - Adjacent cells (8-neighbours) are always mutually visible — there is no cell
///   strictly between them to occlude.
/// - The **endpoints themselves are not tested** for opacity: an observer on, or
///   looking at, an obstacle cell can still establish the line (only cells
///   *between* them mask it). Callers that require clear endpoints should check
///   [`CostGrid::passable`] separately.
/// - An **out-of-bounds** endpoint yields `false` (nothing to see / off-map).
///
/// Pure — no I/O, deterministic in the grid and endpoints.
pub fn line_of_sight(grid: &CostGrid, from: (usize, usize), to: (usize, usize)) -> bool {
    let (x0, y0) = from;
    let (x1, y1) = to;

    // Off-map endpoints can never establish a line of sight.
    if !grid.in_bounds(x0, y0) || !grid.in_bounds(x1, y1) {
        return false;
    }
    // A cell always sees itself; nothing lies between.
    if from == to {
        return true;
    }

    // Amanatides–Woo grid traversal between the two cell *centres*. Work in cell
    // units: the ray starts at the centre of `from` (+0.5, +0.5) and ends at the
    // centre of `to`. `step_*` is the direction of travel along each axis;
    // `t_max_*` is the ray parameter `t` at which the next cell boundary is
    // crossed; `t_delta_*` is the `t` increment per whole cell on that axis.
    let mut x = x0 as isize;
    let mut y = y0 as isize;
    let tx = x1 as isize;
    let ty = y1 as isize;

    let dx = (x1 as f64) - (x0 as f64);
    let dy = (y1 as f64) - (y0 as f64);

    let step_x = dx.signum() as isize; // -1, 0, or +1
    let step_y = dy.signum() as isize;

    // Distance in `t` (0..=1 over the whole segment) to the first boundary, and
    // per-cell thereafter. A zero component means the ray is axis-aligned on that
    // axis: it never crosses a boundary there, so push its `t_max` to +inf.
    let (mut t_max_x, t_delta_x) = if step_x != 0 {
        // From a centre, the first boundary on this axis is half a cell away.
        (0.5 / dx.abs(), 1.0 / dx.abs())
    } else {
        (f64::INFINITY, f64::INFINITY)
    };
    let (mut t_max_y, t_delta_y) = if step_y != 0 {
        (0.5 / dy.abs(), 1.0 / dy.abs())
    } else {
        (f64::INFINITY, f64::INFINITY)
    };

    // Step until we arrive at the target cell. Each iteration advances to the
    // next cell the segment enters; we test that cell (unless it is the target)
    // for opacity. The loop is bounded by the Manhattan distance + a small margin
    // so it can never spin (defensive guard against FP edge cases).
    let max_steps = (tx - x).unsigned_abs() + (ty - y).unsigned_abs() + 2;
    for _ in 0..max_steps {
        // Advance one cell along whichever axis hits its next boundary first.
        // On an exact tie (the ray crosses a lattice corner) step BOTH axes in
        // the same iteration: that visits the shared corner without slipping
        // diagonally between two obstacles, which is the supercover property.
        if (t_max_x - t_max_y).abs() < 1e-9 {
            x += step_x;
            y += step_y;
            t_max_x += t_delta_x;
            t_max_y += t_delta_y;
        } else if t_max_x < t_max_y {
            x += step_x;
            t_max_x += t_delta_x;
        } else {
            y += step_y;
            t_max_y += t_delta_y;
        }

        // Reached the target cell: no obstacle was found strictly between → clear.
        if x == tx && y == ty {
            return true;
        }

        // An opaque cell strictly between the endpoints masks the sight line.
        // (x, y are in-bounds here: they march from one in-bounds cell toward
        // another, and `opaque` treats any out-of-bounds cell as blocking too.)
        if x < 0 || y < 0 || opaque(grid, x as usize, y as usize) {
            return false;
        }
    }

    // Should be unreachable (the target is always reached within max_steps), but
    // fail safe: if the walk did not arrive, report no clear line.
    false
}

/// Whether the observer at cell `from` has a clear **2.5-D terrain-masked** line
/// of sight to the cell `to`, across the elevation of a
/// [`HeightGrid`](crate::terrain::HeightGrid).
///
/// Unlike the flat-occupancy [`line_of_sight`] (which treats `f32::INFINITY` cost
/// cells as opaque walls), this samples the **ground elevation** along the ray and
/// checks it against the straight 3-D **sight line** from the observer's eye to the
/// target. The sight line runs from `(elev[from] + obs_h)` at the observer to
/// `(elev[to] + tgt_h)` at the target; at each intermediate cell the line's height
/// is interpolated, and the view is **BLOCKED** if the terrain there rises above
/// it. This models true **dead ground** behind a ridge: an intervening crest taller
/// than the line of sight hides the target, even though nothing is an "obstacle".
///
/// `obs_h` and `tgt_h` are the observer and target heights **above the ground** in
/// the same metres as the elevation (e.g. a `2 m` standing observer, a `0 m`
/// ground target). A taller observer (`obs_h` large) can see *over* a low hill.
///
/// Endpoint handling mirrors [`line_of_sight`]: `from == to` is `true`, adjacent
/// cells are mutually visible (nothing strictly between), the **endpoint cells
/// themselves are not occlusion-tested**, and an **out-of-bounds** endpoint yields
/// `false`. The traversal is the same supercover DDA, so a ridge the ideal ray only
/// clips still masks the view. Pure — deterministic in the terrain and endpoints.
pub fn line_of_sight_terrain(
    height: &crate::terrain::HeightGrid,
    from: (usize, usize),
    to: (usize, usize),
    obs_h: f32,
    tgt_h: f32,
) -> bool {
    let (x0, y0) = from;
    let (x1, y1) = to;

    if !height.in_bounds(x0, y0) || !height.in_bounds(x1, y1) {
        return false;
    }
    if from == to {
        return true;
    }

    // Eye height at the observer and aim height at the target (ground + offset).
    let z_from = height.elevation_at(x0, y0) as f64 + obs_h as f64;
    let z_to = height.elevation_at(x1, y1) as f64 + tgt_h as f64;

    // Total horizontal span, used to interpolate the sight-line height by fraction
    // of distance covered (param t in 0..=1 along the segment).
    let dx = (x1 as f64) - (x0 as f64);
    let dy = (y1 as f64) - (y0 as f64);
    let span2 = dx * dx + dy * dy;
    if span2 == 0.0 {
        return true; // same cell (already handled, but guards the divide below).
    }

    let mut x = x0 as isize;
    let mut y = y0 as isize;
    let tx = x1 as isize;
    let ty = y1 as isize;

    let step_x = dx.signum() as isize;
    let step_y = dy.signum() as isize;
    let (mut t_max_x, t_delta_x) = if step_x != 0 {
        (0.5 / dx.abs(), 1.0 / dx.abs())
    } else {
        (f64::INFINITY, f64::INFINITY)
    };
    let (mut t_max_y, t_delta_y) = if step_y != 0 {
        (0.5 / dy.abs(), 1.0 / dy.abs())
    } else {
        (f64::INFINITY, f64::INFINITY)
    };

    let max_steps = (tx - x).unsigned_abs() + (ty - y).unsigned_abs() + 2;
    for _ in 0..max_steps {
        if (t_max_x - t_max_y).abs() < 1e-9 {
            x += step_x;
            y += step_y;
            t_max_x += t_delta_x;
            t_max_y += t_delta_y;
        } else if t_max_x < t_max_y {
            x += step_x;
            t_max_x += t_delta_x;
        } else {
            y += step_y;
            t_max_y += t_delta_y;
        }

        if x == tx && y == ty {
            return true; // reached target; nothing rose above the sight line.
        }
        if x < 0 || y < 0 {
            return false;
        }

        // Fraction of the horizontal distance covered to this cell centre, and the
        // sight-line height there (linear interpolation between eye and aim).
        let cdx = (x as f64) - (x0 as f64);
        let cdy = (y as f64) - (y0 as f64);
        let t = ((cdx * cdx + cdy * cdy) / span2).sqrt().clamp(0.0, 1.0);
        let line_z = z_from + (z_to - z_from) * t;

        // Dead ground: terrain here taller than the sight line masks the target.
        let ground_z = height.elevation_at(x as usize, y as usize) as f64;
        if ground_z > line_z {
            return false;
        }
    }

    false // fail safe — target should have been reached within max_steps.
}

/// For an `observer` cell, whether it has line of sight to **each** target in
/// `targets`, returned in the same order (`out[i] == line_of_sight(grid,
/// observer, targets[i])`).
///
/// A convenience batch wrapper over [`line_of_sight`] for the common "what can
/// this sensor position see among these contacts?" query. Pure; the result has
/// the same length as `targets`.
pub fn visible_from(
    grid: &CostGrid,
    observer: (usize, usize),
    targets: &[(usize, usize)],
) -> Vec<bool> {
    targets
        .iter()
        .map(|&t| line_of_sight(grid, observer, t))
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    // These tests index grids as the readable `row * width + col`; with a literal
    // row of 0 or 1 that form trips identity_op (`1 * w`, `+ 0`) and erasing_op
    // (`0 * w`). The explicit arithmetic documents the (row, col) cell being poked,
    // so keep it and silence the style lints for the test module only.
    #![allow(clippy::identity_op, clippy::erasing_op)]
    use super::*;

    #[test]
    fn clear_on_open_grid() {
        // Nothing blocks across a fully open field, in any direction.
        let grid = CostGrid::uniform(10, 8, 1.0);
        assert!(line_of_sight(&grid, (0, 0), (9, 7)));
        assert!(line_of_sight(&grid, (9, 7), (0, 0)));
        assert!(line_of_sight(&grid, (0, 7), (9, 0)));
        assert!(line_of_sight(&grid, (5, 0), (5, 7))); // vertical
        assert!(line_of_sight(&grid, (0, 3), (9, 3))); // horizontal
    }

    #[test]
    fn from_equals_to_is_visible() {
        let grid = CostGrid::uniform(5, 5, 1.0);
        assert!(line_of_sight(&grid, (2, 2), (2, 2)));
        // Even an obstacle cell sees "itself" (endpoints are not occlusion-tested).
        let mut g = CostGrid::uniform(5, 5, 1.0);
        g.cost[2 * 5 + 2] = f32::INFINITY;
        assert!(line_of_sight(&g, (2, 2), (2, 2)));
    }

    #[test]
    fn adjacent_cells_are_always_visible() {
        // Even ringed by obstacles, two touching cells see each other (nothing
        // lies strictly between 8-neighbours).
        let mut grid = CostGrid::uniform(4, 4, 1.0);
        grid.block_rect(0, 0, 4, 4); // everything opaque...
        grid.cost[1 * 4 + 1] = 1.0; // ...but re-open two adjacent cells.
        grid.cost[1 * 4 + 2] = 1.0;
        assert!(line_of_sight(&grid, (1, 1), (2, 1))); // orthogonal neighbour
                                                       // A diagonal neighbour too (reopen one more cell).
        let mut g2 = CostGrid::uniform(4, 4, 1.0);
        assert!(line_of_sight(&g2, (1, 1), (2, 2)));
        g2.cost[1 * 4 + 1] = f32::INFINITY; // even if the source cell is opaque
        assert!(line_of_sight(&g2, (1, 1), (2, 2)));
    }

    #[test]
    fn blocked_by_a_wall_between() {
        // A vertical wall between the observer and target masks the horizontal
        // sight line.
        let mut grid = CostGrid::uniform(7, 5, 1.0);
        grid.block_rect(3, 0, 4, 5); // full-height wall at column x=3
        assert!(!line_of_sight(&grid, (0, 2), (6, 2)));
        // Removing the wall restores the view.
        let open = CostGrid::uniform(7, 5, 1.0);
        assert!(line_of_sight(&open, (0, 2), (6, 2)));
    }

    #[test]
    fn partial_wall_masks_only_behind_it() {
        // A wall spanning the top rows blocks a line that passes behind it but
        // not one that goes under the gap.
        let mut grid = CostGrid::uniform(7, 5, 1.0);
        grid.block_rect(3, 0, 4, 3); // column x=3, rows 0..3 (gap at y=3,4)
                                     // A line through row 1 crosses the wall → blocked.
        assert!(!line_of_sight(&grid, (0, 1), (6, 1)));
        // A line through the open bottom row is clear.
        assert!(line_of_sight(&grid, (0, 4), (6, 4)));
    }

    #[test]
    fn endpoint_on_obstacle_is_handled() {
        // The far endpoint sits on an obstacle, but the cells *between* are open:
        // the observer can still see the (occupied) target cell.
        let mut grid = CostGrid::uniform(6, 1, 1.0);
        grid.cost[5] = f32::INFINITY; // target (5,0) is an obstacle
        assert!(line_of_sight(&grid, (0, 0), (5, 0)));
        // But an obstacle *one short* of the endpoint masks it.
        let mut g2 = CostGrid::uniform(6, 1, 1.0);
        g2.cost[4] = f32::INFINITY; // (4,0) between observer and target
        assert!(!line_of_sight(&g2, (0, 0), (5, 0)));
    }

    #[test]
    fn supercover_does_not_leak_past_a_clipped_corner() {
        // The supercover property: a wall the *ideal* ray only clips the CORNER
        // of must still block it (Bresenham could skip that cell). Observer
        // (0,1) -> target (6,4) over a wall at column x=3, rows 0..3: the line
        // enters cell (3,2), which is opaque, so the view is masked even though a
        // naive line might thread just past the wall's lower corner.
        let mut grid = CostGrid::uniform(7, 5, 1.0);
        grid.block_rect(3, 0, 4, 3);
        assert!(
            !line_of_sight(&grid, (0, 1), (6, 4)),
            "a wall the sight line clips the corner of must block it"
        );
        // With that wall cell re-opened the same line is clear (proving it was
        // exactly the clipped cell that mattered).
        grid.cost[2 * grid.w + 3] = 1.0; // re-open (3,2)
        assert!(line_of_sight(&grid, (0, 1), (6, 4)));
    }

    #[test]
    fn exact_diagonal_through_open_pinhole_is_visible() {
        // A pure 45° ray (0,0)->(2,2) steps the shared corner (1,1) directly. If
        // (1,1) is open the line is clear, even with obstacles flanking the seam
        // at (1,0) and (0,1) (they only touch the ray at a lattice point).
        let mut grid = CostGrid::uniform(3, 3, 1.0);
        grid.cost[0 * 3 + 1] = f32::INFINITY; // (1,0)
        grid.cost[1 * 3 + 0] = f32::INFINITY; // (0,1)
        assert!(line_of_sight(&grid, (0, 0), (2, 2)));
        // Close the on-ray corner cell and the diagonal is masked.
        grid.cost[1 * 3 + 1] = f32::INFINITY; // (1,1)
        assert!(!line_of_sight(&grid, (0, 0), (2, 2)));
    }

    #[test]
    fn out_of_bounds_endpoint_is_not_visible() {
        let grid = CostGrid::uniform(4, 4, 1.0);
        assert!(!line_of_sight(&grid, (9, 0), (0, 0)));
        assert!(!line_of_sight(&grid, (0, 0), (0, 9)));
    }

    #[test]
    fn symmetry_holds_for_random_ish_obstacles() {
        // LoS is symmetric: A sees B iff B sees A. Check it over a scattered field.
        let mut grid = CostGrid::uniform(9, 9, 1.0);
        grid.block_rect(2, 2, 3, 7);
        grid.block_rect(5, 1, 6, 5);
        grid.block_rect(6, 6, 8, 7);
        for ay in 0..9 {
            for ax in 0..9 {
                for by in 0..9 {
                    for bx in 0..9 {
                        assert_eq!(
                            line_of_sight(&grid, (ax, ay), (bx, by)),
                            line_of_sight(&grid, (bx, by), (ax, ay)),
                            "LoS must be symmetric for ({ax},{ay}) <-> ({bx},{by})"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn terrain_los_clear_over_flat_ground() {
        use crate::terrain::HeightGrid;
        // A flat plateau: every line of sight is clear (no terrain rises to mask).
        let t = HeightGrid::flat(12, 8, 300.0);
        assert!(line_of_sight_terrain(&t, (0, 0), (11, 7), 2.0, 0.0));
        assert!(line_of_sight_terrain(&t, (11, 7), (0, 0), 2.0, 0.0));
        assert!(line_of_sight_terrain(&t, (0, 4), (11, 4), 2.0, 2.0)); // horizontal
    }

    #[test]
    fn terrain_los_blocked_by_a_ridge_between() {
        use crate::terrain::HeightGrid;
        // Build a wall-like ridge: a tall north-south crest in the middle column.
        // Two ground-level observers either side of it cannot see each other —
        // the crest rises above the (flat) sight line between them = dead ground.
        let w = 11;
        let h = 5;
        let mut elev = vec![100.0f32; w * h];
        for y in 0..h {
            elev[y * w + 5] = 900.0; // a sheer ridge at column x=5
        }
        let t = HeightGrid { w, h, elev };
        assert!(
            !line_of_sight_terrain(&t, (0, 2), (10, 2), 2.0, 2.0),
            "a tall ridge between the two cells must mask the view (dead ground)"
        );
        // Without the ridge (flat) the same line is clear.
        let flat = HeightGrid::flat(w, h, 100.0);
        assert!(line_of_sight_terrain(&flat, (0, 2), (10, 2), 2.0, 2.0));
    }

    #[test]
    fn tall_observer_sees_over_a_low_hill() {
        use crate::terrain::HeightGrid;
        // A modest hill (a single raised cell) sits between observer and target.
        // A ground-level observer is masked by it, but raising the observer high
        // enough lifts the sight line above the hill crest → target becomes visible.
        let w = 9;
        let h = 3;
        let mut elev = vec![100.0f32; w * h];
        elev[1 * w + 4] = 130.0; // a low hill (+30 m) at the midpoint (4,1)
        let t = HeightGrid { w, h, elev };
        // Low observer & target (1 m): the hill rises above the flat sight line.
        assert!(
            !line_of_sight_terrain(&t, (0, 1), (8, 1), 1.0, 1.0),
            "a low observer should be masked by the intervening hill"
        );
        // A tall observer (100 m mast) lifts the sight line above the +30 m hill.
        assert!(
            line_of_sight_terrain(&t, (0, 1), (8, 1), 100.0, 1.0),
            "a tall-enough observer should see over the low hill"
        );
    }

    #[test]
    fn terrain_los_endpoints_and_bounds() {
        use crate::terrain::HeightGrid;
        let t = HeightGrid::flat(5, 5, 200.0);
        // A cell sees itself; adjacent cells are mutually visible.
        assert!(line_of_sight_terrain(&t, (2, 2), (2, 2), 2.0, 0.0));
        assert!(line_of_sight_terrain(&t, (2, 2), (3, 2), 2.0, 0.0));
        // Out-of-bounds endpoints are never visible.
        assert!(!line_of_sight_terrain(&t, (9, 0), (0, 0), 2.0, 0.0));
        assert!(!line_of_sight_terrain(&t, (0, 0), (0, 9), 2.0, 0.0));
    }

    #[test]
    fn terrain_and_flat_los_are_independent_models() {
        use crate::terrain::demo_terrain;
        // The demo landscape has a steep diagonal ridge; a ground observer in the
        // NW corner looking to the SE corner crosses the ridge → masked in 2.5-D.
        let t = demo_terrain(48, 48);
        assert!(
            !line_of_sight_terrain(&t, (2, 2), (45, 45), 2.0, 2.0),
            "the SE target is in dead ground behind the demo ridge"
        );
    }

    #[test]
    fn visible_from_batches_in_order() {
        // A wall at column x=3 (rows 0..3) blocks the targets behind it on rows
        // 0..2 but not one along the fully-open bottom row.
        let mut grid = CostGrid::uniform(7, 5, 1.0);
        grid.block_rect(3, 0, 4, 3);
        let observer = (0, 1);
        let targets = [(6, 1), (6, 4), (1, 1)];
        let got = visible_from(&grid, observer, &targets);
        assert_eq!(got.len(), targets.len());
        assert!(!got[0], "(6,1) is behind the wall on the same row → masked");
        // (0,1)->(6,4) clips the wall's lower cell (3,2) → also masked.
        assert!(!got[1], "(6,4)'s line clips the wall corner → masked");
        assert!(got[2], "(1,1) is adjacent/open → visible");
        // A target reachable along the fully open bottom row IS visible.
        assert!(
            line_of_sight(&grid, (0, 4), (6, 4)),
            "open bottom row is clear"
        );
        // Equivalent to calling line_of_sight per target.
        for (i, &t) in targets.iter().enumerate() {
            assert_eq!(got[i], line_of_sight(&grid, observer, t));
        }
    }

    #[test]
    fn empty_targets_yield_empty() {
        let grid = CostGrid::uniform(4, 4, 1.0);
        assert!(visible_from(&grid, (0, 0), &[]).is_empty());
    }
}
