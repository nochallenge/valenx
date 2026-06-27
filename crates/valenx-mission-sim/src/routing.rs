//! # `routing` — in-house **A\* tactical routing** over a 2-D grid cost field.
//!
//! A small, dependency-free **A\*** path planner over a rectangular grid of
//! per-cell traversal costs ([`CostGrid`]). It finds a least-cost 8-connected
//! path from a start cell to a goal cell across a **cost field** — open ground is
//! cheap, difficult ground is expensive, and an obstacle is a cell of
//! `f32::INFINITY` cost that the search never enters.
//!
//! This is route *planning* only — a defensive / analysis capability (find a way
//! **around** obstacles), with no targeting, no adversary model, and no weapon
//! semantics. It complements [`crate::planner`]'s waypoint routes: where the
//! planner moves entities along pre-set legs, this computes a leg-by-leg path
//! through a cost field for the user to inspect.
//!
//! ## Model
//!
//! - The grid is `w × h` cells in **row-major** order (`idx = y · w + x`), each
//!   holding an `f32` *traversal cost* (the cost to step **into** that cell).
//!   `f32::INFINITY` marks an impassable **obstacle**; any finite `>= 0` value is
//!   passable terrain (`1.0` = nominal).
//! - The search is **8-connected** (4 orthogonal + 4 diagonal neighbours). An
//!   orthogonal step costs the destination cell's cost; a **diagonal** step costs
//!   `√2 ×` the destination cost (so diagonals are not "free"), and a diagonal is
//!   only taken when it does **not** cut the corner of an obstacle (both shared
//!   orthogonal neighbours must be passable).
//! - The heuristic is the admissible **octile distance** scaled by the minimum
//!   passable cell cost, so A\* stays optimal on a non-uniform cost field.
//!
//! Everything here is **pure** (no I/O, no clock, no randomness): the same grid
//! and endpoints always yield the same path.
//!
//! ## Example
//!
//! ```
//! use valenx_mission_sim::routing::{astar, CostGrid};
//!
//! // 5×1 open corridor: the straight path is start..=goal inclusive.
//! let grid = CostGrid::uniform(5, 1, 1.0);
//! let path = astar(&grid, (0, 0), (4, 0)).expect("reachable");
//! assert_eq!(path.first(), Some(&(0, 0)));
//! assert_eq!(path.last(), Some(&(4, 0)));
//! ```

use std::collections::BinaryHeap;

/// A rectangular grid of per-cell **traversal costs** for [`astar`].
///
/// Cells are stored **row-major**: cell `(x, y)` lives at `cost[y * w + x]`,
/// with `x` in `0..w` (column) and `y` in `0..h` (row). Each entry is the cost
/// to *enter* that cell; `f32::INFINITY` is an impassable **obstacle** the search
/// never steps into. Finite values should be `>= 0` (a nominal traversable cell
/// is `1.0`); larger values model difficult ground the planner will route around
/// when cheaper.
#[derive(Debug, Clone, PartialEq)]
pub struct CostGrid {
    /// Width in cells (number of columns; `x` ranges `0..w`).
    pub w: usize,
    /// Height in cells (number of rows; `y` ranges `0..h`).
    pub h: usize,
    /// Row-major traversal cost per cell (`len == w * h`). `f32::INFINITY` is an
    /// obstacle.
    pub cost: Vec<f32>,
}

impl CostGrid {
    /// A `w × h` grid with every cell set to the same finite `cost` (no
    /// obstacles). Useful for tests and as a baseline field.
    pub fn uniform(w: usize, h: usize, cost: f32) -> Self {
        Self {
            w,
            h,
            cost: vec![cost; w.saturating_mul(h)],
        }
    }

    /// Whether `(x, y)` is inside the grid bounds.
    #[inline]
    pub fn in_bounds(&self, x: usize, y: usize) -> bool {
        x < self.w && y < self.h
    }

    /// The traversal cost of cell `(x, y)` (`f32::INFINITY` for an out-of-bounds
    /// cell or an obstacle).
    #[inline]
    pub fn cost_at(&self, x: usize, y: usize) -> f32 {
        if self.in_bounds(x, y) {
            self.cost[y * self.w + x]
        } else {
            f32::INFINITY
        }
    }

    /// Whether cell `(x, y)` is passable (in bounds and finite, non-negative
    /// cost). An `f32::INFINITY` cell — or anything out of bounds — is blocked.
    #[inline]
    pub fn passable(&self, x: usize, y: usize) -> bool {
        let c = self.cost_at(x, y);
        c.is_finite() && c >= 0.0
    }

    /// Mark the rectangle `[x0, x1) × [y0, y1)` (clamped to the grid) as an
    /// `f32::INFINITY` obstacle wall. A convenience used by [`demo_field`].
    pub fn block_rect(&mut self, x0: usize, y0: usize, x1: usize, y1: usize) {
        let x1 = x1.min(self.w);
        let y1 = y1.min(self.h);
        for y in y0..y1 {
            for x in x0..x1 {
                self.cost[y * self.w + x] = f32::INFINITY;
            }
        }
    }

    /// Derive a routing cost field from a terrain [`HeightGrid`](crate::terrain::HeightGrid)
    /// by its **slope**: gentle ground is cheap, steep ground is expensive, and an
    /// impassably steep slope is an obstacle. This is what makes A\* routing
    /// **terrain-aware** — the planner prefers gentle terrain and routes around
    /// steep ridges.
    ///
    /// For each cell the cost is `base + slope_k · slope(cell)` where
    /// `slope(cell)` is the dimensionless rise-over-run from
    /// [`HeightGrid::slope_at`](crate::terrain::HeightGrid::slope_at). When that
    /// slope **exceeds** `impassable_slope` (and `impassable_slope > 0`) the cell
    /// is set to `f32::INFINITY` — a wall the route never enters (e.g. a cliff). A
    /// non-positive `impassable_slope` disables the cutoff (every cell is
    /// passable). `base` is floored at a small positive value so every passable
    /// cell has a positive cost (keeps the A\* heuristic well-scaled).
    ///
    /// The returned grid has the **same dimensions** as `height`, indexing
    /// cell-for-cell, so a derived cost grid and the terrain overlay align exactly.
    /// Pure — deterministic in the terrain and parameters.
    pub fn from_terrain(
        height: &crate::terrain::HeightGrid,
        base: f32,
        slope_k: f32,
        impassable_slope: f32,
    ) -> Self {
        let (w, h) = (height.w, height.h);
        let base = base.max(1e-3);
        let mut cost = vec![base; w.saturating_mul(h)];
        for y in 0..h {
            for x in 0..w {
                let s = height.slope_at(x, y);
                cost[y * w + x] = if impassable_slope > 0.0 && s > impassable_slope {
                    f32::INFINITY
                } else {
                    base + slope_k.max(0.0) * s
                };
            }
        }
        Self { w, h, cost }
    }

    /// The smallest **finite, positive** cell cost in the grid (used to scale the
    /// admissible heuristic). Falls back to `1.0` when the grid has no positive
    /// finite cell (all-zero or all-obstacle), keeping the heuristic valid.
    fn min_step_cost(&self) -> f32 {
        let mut m = f32::INFINITY;
        for &c in &self.cost {
            if c.is_finite() && c > 0.0 && c < m {
                m = c;
            }
        }
        if m.is_finite() {
            m
        } else {
            1.0
        }
    }
}

/// A demonstration cost field: a `w × h` grid of nominal `1.0` ground with a
/// couple of obstacle **walls** that a straight start→goal line would hit, so a
/// route must detour around them. The walls each leave a gap, so a path always
/// exists between opposite corners.
///
/// Wall layout (proportional to the grid, so it scales): a vertical wall at
/// `x ≈ w/3` spanning the **top** with a gap near the bottom, and a second
/// vertical wall at `x ≈ 2w/3` spanning the **bottom** with a gap near the top —
/// the classic offset-slit detour.
pub fn demo_field(w: usize, h: usize) -> CostGrid {
    let mut grid = CostGrid::uniform(w, h, 1.0);
    if w < 4 || h < 4 {
        return grid; // too small to place meaningful walls; leave it open.
    }
    let x1 = w / 3;
    let x2 = (2 * w) / 3;
    let gap = (h / 4).max(1); // gap height left open at one end of each wall.

    // First wall: column x1, blocked from the top down, gap at the bottom.
    grid.block_rect(x1, 0, x1 + 1, h - gap);
    // Second wall: column x2, blocked from the bottom up, gap at the top.
    grid.block_rect(x2, gap, x2 + 1, h);
    grid
}

/// Open-set node for the A\* min-heap. `BinaryHeap` is a *max*-heap, so [`Ord`]
/// is reversed on `f` (then `g`) to pop the **lowest** `f`-score first.
#[derive(Copy, Clone)]
struct Node {
    /// Estimated total cost `f = g + h` through this cell.
    f: f32,
    /// Cost-so-far from the start to this cell.
    g: f32,
    /// Flat row-major index of the cell.
    idx: usize,
}

impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.f == other.f && self.g == other.g
    }
}
impl Eq for Node {}
impl Ord for Node {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Reverse so the BinaryHeap yields the smallest f (ties: smallest g)
        // first. NaN cannot arise (all f/g are finite when pushed).
        other
            .f
            .partial_cmp(&self.f)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                other
                    .g
                    .partial_cmp(&self.g)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    }
}
impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// The 8 neighbour offsets `(dx, dy)`; the first four are orthogonal, the last
/// four diagonal. Diagonal steps cost [`DIAG`]× and must not cut an obstacle
/// corner.
const NEIGHBORS: [(isize, isize); 8] = [
    (1, 0),
    (-1, 0),
    (0, 1),
    (0, -1),
    (1, 1),
    (1, -1),
    (-1, 1),
    (-1, -1),
];

/// Diagonal step-cost multiplier `√2` (so a diagonal move is geometrically
/// longer than an orthogonal one).
const DIAG: f32 = std::f32::consts::SQRT_2;

/// Admissible **octile** heuristic between two cells, scaled by `min_cost` (the
/// cheapest passable cell). Octile distance is the exact 8-connected grid
/// distance ignoring obstacles: `(dmax − dmin) + √2 · dmin`, i.e. travel
/// diagonally while both axes have distance left, then straight. Scaling by the
/// minimum step cost keeps it a lower bound on the true cost over a non-uniform
/// field, so A\* remains optimal.
#[inline]
fn octile(ax: usize, ay: usize, bx: usize, by: usize, min_cost: f32) -> f32 {
    let dx = (ax as isize - bx as isize).unsigned_abs() as f32;
    let dy = (ay as isize - by as isize).unsigned_abs() as f32;
    let (dmin, dmax) = if dx < dy { (dx, dy) } else { (dy, dx) };
    (dmax - dmin + DIAG * dmin) * min_cost
}

/// Compute the least-cost **8-connected A\*** path across `grid` from `start` to
/// `goal`, both `(x, y)` cells.
///
/// Returns the contiguous cell path **including both endpoints**
/// (`start..=goal`), where each consecutive pair is 8-adjacent, or `None` when
/// the goal is unreachable (walled off), either endpoint is out of bounds, or
/// either endpoint is an obstacle (`f32::INFINITY`). When `start == goal` (and it
/// is passable) the path is the single cell `[start]`.
///
/// Cost model: stepping into a cell costs that cell's [`CostGrid`] value; a
/// diagonal step costs `√2 ×` that value and is forbidden when it would cut the
/// corner of an obstacle (either shared orthogonal neighbour is blocked).
/// `f32::INFINITY` cells are never entered. Pure — no I/O.
pub fn astar(
    grid: &CostGrid,
    start: (usize, usize),
    goal: (usize, usize),
) -> Option<Vec<(usize, usize)>> {
    let (sx, sy) = start;
    let (gx, gy) = goal;

    // Reject degenerate grids and impassable / out-of-bounds endpoints up front.
    if grid.w == 0 || grid.h == 0 {
        return None;
    }
    if !grid.passable(sx, sy) || !grid.passable(gx, gy) {
        return None;
    }

    let w = grid.w;
    let n = w * grid.h;
    let start_idx = sy * w + sx;
    let goal_idx = gy * w + gx;

    // start == goal: trivial single-cell path (endpoint already verified passable).
    if start_idx == goal_idx {
        return Some(vec![start]);
    }

    let min_cost = grid.min_step_cost();

    // g-score per cell (INFINITY = unvisited) and predecessor for path rebuild.
    let mut g_score = vec![f32::INFINITY; n];
    let mut came_from = vec![usize::MAX; n];
    let mut closed = vec![false; n];

    g_score[start_idx] = 0.0;
    let mut open = BinaryHeap::new();
    open.push(Node {
        f: octile(sx, sy, gx, gy, min_cost),
        g: 0.0,
        idx: start_idx,
    });

    while let Some(Node { g, idx, .. }) = open.pop() {
        if idx == goal_idx {
            return Some(reconstruct(&came_from, goal_idx, w));
        }
        // Stale heap entry (a cheaper path to this cell was already expanded).
        if closed[idx] {
            continue;
        }
        closed[idx] = true;

        let cx = idx % w;
        let cy = idx / w;

        for (k, &(dx, dy)) in NEIGHBORS.iter().enumerate() {
            let nx = cx as isize + dx;
            let ny = cy as isize + dy;
            if nx < 0 || ny < 0 {
                continue;
            }
            let (nx, ny) = (nx as usize, ny as usize);
            if !grid.passable(nx, ny) {
                continue;
            }
            let diagonal = k >= 4;
            if diagonal {
                // Forbid cutting an obstacle corner: both orthogonal cells that
                // the diagonal "slips between" must be passable.
                let ortho_a = grid.passable((cx as isize + dx) as usize, cy);
                let ortho_b = grid.passable(cx, (cy as isize + dy) as usize);
                if !ortho_a || !ortho_b {
                    continue;
                }
            }
            let nidx = ny * w + nx;
            if closed[nidx] {
                continue;
            }
            let step = grid.cost_at(nx, ny) * if diagonal { DIAG } else { 1.0 };
            let tentative = g + step;
            if tentative < g_score[nidx] {
                g_score[nidx] = tentative;
                came_from[nidx] = idx;
                open.push(Node {
                    f: tentative + octile(nx, ny, gx, gy, min_cost),
                    g: tentative,
                    idx: nidx,
                });
            }
        }
    }

    None // open set exhausted: goal unreachable.
}

/// Walk the `came_from` predecessor chain back from `goal_idx` to the start and
/// return the path in forward order (`start..=goal`) as `(x, y)` cells.
fn reconstruct(came_from: &[usize], goal_idx: usize, w: usize) -> Vec<(usize, usize)> {
    let mut path = Vec::new();
    let mut cur = goal_idx;
    loop {
        path.push((cur % w, cur / w));
        let prev = came_from[cur];
        if prev == usize::MAX {
            break; // reached the start (no predecessor).
        }
        cur = prev;
    }
    path.reverse();
    path
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Every consecutive pair in a path must be 8-adjacent (the path is
    /// contiguous, no teleporting between cells).
    fn is_contiguous(path: &[(usize, usize)]) -> bool {
        path.windows(2).all(|w| {
            let dx = (w[0].0 as isize - w[1].0 as isize).abs();
            let dy = (w[0].1 as isize - w[1].1 as isize).abs();
            dx <= 1 && dy <= 1 && (dx + dy) != 0
        })
    }

    #[test]
    fn open_grid_is_near_straight() {
        // On an open uniform field, the diagonal-then-straight path from corner
        // to corner has length max(dx, dy) + 1 cells (octile-optimal).
        let grid = CostGrid::uniform(10, 6, 1.0);
        let path = astar(&grid, (0, 0), (9, 5)).expect("open grid is reachable");
        assert_eq!(path.first(), Some(&(0, 0)));
        assert_eq!(path.last(), Some(&(9, 5)));
        // dx = 9, dy = 5 → 9 steps → 10 cells inclusive.
        assert_eq!(path.len(), 10, "octile-optimal corner path length");
        assert!(is_contiguous(&path));
    }

    #[test]
    fn routes_around_an_obstacle_wall() {
        // A vertical wall spanning all but the bottom row forces a detour; the
        // straight horizontal line through the middle is blocked.
        let mut grid = CostGrid::uniform(7, 5, 1.0);
        grid.block_rect(3, 0, 4, 4); // column x=3, rows 0..4 blocked (gap at y=4).
        let path = astar(&grid, (0, 2), (6, 2)).expect("a detour around the wall exists");
        assert_eq!(path.first(), Some(&(0, 2)));
        assert_eq!(path.last(), Some(&(6, 2)));
        assert!(is_contiguous(&path));
        // The path must never step onto the wall column above the gap.
        assert!(
            path.iter().all(|&(x, y)| !(x == 3 && y < 4)),
            "route must go around, not through, the wall: {path:?}"
        );
    }

    #[test]
    fn fully_walled_goal_is_unreachable() {
        // Ring the goal cell with obstacles on all sides → no path.
        let mut grid = CostGrid::uniform(7, 7, 1.0);
        let (gx, gy) = (3usize, 3usize);
        grid.block_rect(gx - 1, gy - 1, gx + 2, gy + 2); // 3×3 wall...
        grid.cost[gy * grid.w + gx] = 1.0; // ...but re-open the goal itself.
        assert!(
            astar(&grid, (0, 0), (gx, gy)).is_none(),
            "a goal enclosed by obstacles must be unreachable"
        );
    }

    #[test]
    fn start_equals_goal_is_single_cell() {
        let grid = CostGrid::uniform(5, 5, 1.0);
        let path = astar(&grid, (2, 2), (2, 2)).expect("trivial path");
        assert_eq!(path, vec![(2, 2)]);
    }

    #[test]
    fn infinite_cost_cells_are_never_entered() {
        // Scatter obstacles; any returned path must avoid every INFINITY cell.
        let mut grid = CostGrid::uniform(8, 8, 1.0);
        grid.block_rect(2, 2, 3, 6);
        grid.block_rect(5, 0, 6, 4);
        let path = astar(&grid, (0, 0), (7, 7)).expect("a route threads the gaps");
        assert!(is_contiguous(&path));
        for &(x, y) in &path {
            assert!(
                grid.cost_at(x, y).is_finite(),
                "path entered an obstacle cell at ({x},{y})"
            );
        }
    }

    #[test]
    fn obstacle_endpoint_is_rejected() {
        let mut grid = CostGrid::uniform(5, 5, 1.0);
        grid.cost[0] = f32::INFINITY; // start (0,0) is an obstacle.
        assert!(astar(&grid, (0, 0), (4, 4)).is_none());
        let mut grid2 = CostGrid::uniform(5, 5, 1.0);
        grid2.cost[4 * 5 + 4] = f32::INFINITY; // goal (4,4) is an obstacle.
        assert!(astar(&grid2, (0, 0), (4, 4)).is_none());
    }

    #[test]
    fn out_of_bounds_endpoints_are_rejected() {
        let grid = CostGrid::uniform(4, 4, 1.0);
        assert!(astar(&grid, (9, 0), (0, 0)).is_none());
        assert!(astar(&grid, (0, 0), (0, 9)).is_none());
    }

    #[test]
    fn prefers_cheaper_terrain_over_a_cost_ridge() {
        // A vertical ridge of very expensive (but passable) cells across the
        // middle row band; the optimal path should detour through the cheap gap
        // rather than plough straight through the ridge.
        let mut grid = CostGrid::uniform(7, 5, 1.0);
        for y in 0..4 {
            grid.cost[y * 7 + 3] = 50.0; // costly column, gap at y=4.
        }
        let path = astar(&grid, (0, 2), (6, 2)).expect("reachable");
        let ridge_cells = path.iter().filter(|&&(x, y)| x == 3 && y < 4).count();
        assert_eq!(
            ridge_cells, 0,
            "A* should avoid the costly ridge via the cheap gap: {path:?}"
        );
        assert!(is_contiguous(&path));
    }

    #[test]
    fn demo_field_has_walls_but_is_solvable() {
        let grid = demo_field(30, 20);
        // It must contain at least some obstacles...
        assert!(
            grid.cost.iter().any(|c| c.is_infinite()),
            "demo field should place obstacle walls"
        );
        // ...yet a corner-to-corner route must still exist around them.
        let path = astar(&grid, (0, 0), (29, 19)).expect("demo field is solvable");
        assert!(is_contiguous(&path));
        for &(x, y) in &path {
            assert!(grid.cost_at(x, y).is_finite());
        }
    }

    #[test]
    fn from_terrain_makes_routing_slope_aware() {
        use crate::terrain::demo_terrain;
        // Derive a cost field from the demo landscape. The long diagonal ridge has
        // very steep FLANKS (rise-over-run ~5-6) while open valleys are gentle
        // (~0.9), so a slope cutoff of 4.0 turns the steep ridge flanks into
        // impassable walls and leaves the rest of the field traversable. A
        // least-cost route from the NW to the SE corner must then go AROUND those
        // steep flanks rather than climbing them.
        let terrain = demo_terrain(40, 40);
        let cutoff = 4.0;
        let grid = CostGrid::from_terrain(&terrain, 1.0, 5.0, cutoff);
        // The derived grid matches the terrain dimensions cell-for-cell.
        assert_eq!((grid.w, grid.h), (terrain.w, terrain.h));
        // Some cells (the steep ridge flanks) must have been marked impassable...
        let blocked = grid.cost.iter().filter(|c| c.is_infinite()).count();
        assert!(blocked > 0, "steep ridge flanks should be impassable");
        // ...but not so many that the whole map is walled off.
        assert!(
            blocked < grid.cost.len() / 2,
            "most of the field should remain traversable, blocked {blocked}/{}",
            grid.cost.len()
        );
        // Routing across the field still succeeds (it skirts the steep flanks).
        let path = astar(&grid, (0, 0), (39, 39)).expect("a route around the steep ridge exists");
        assert!(is_contiguous(&path));
        // The route never steps onto a too-steep (impassable) cell, and every cell
        // it uses is below the impassability slope — i.e. it prefers gentle ground.
        for &(x, y) in &path {
            assert!(
                grid.cost_at(x, y).is_finite(),
                "slope-aware route entered an impassable ridge cell at ({x},{y})"
            );
            assert!(
                terrain.slope_at(x, y) <= cutoff,
                "route used a slope steeper than the cutoff at ({x},{y})"
            );
        }
    }

    #[test]
    fn from_terrain_cutoff_disabled_keeps_all_cells_passable() {
        use crate::terrain::demo_terrain;
        // A non-positive impassable_slope disables the cutoff: no cell is INFINITY,
        // but steep cells are still more expensive than gentle ones.
        let terrain = demo_terrain(30, 30);
        let grid = CostGrid::from_terrain(&terrain, 1.0, 500.0, 0.0);
        assert!(
            grid.cost.iter().all(|c| c.is_finite()),
            "no cutoff => every cell passable"
        );
        // A cell on the steep ridge costs more than a gentle corner cell.
        let ridge = grid.cost_at(15, 15);
        let corner = grid.cost_at(2, 27);
        assert!(
            ridge > corner,
            "steep ridge cost {ridge} should exceed gentle corner cost {corner}"
        );
    }

    #[test]
    fn diagonal_does_not_cut_obstacle_corners() {
        // Two obstacles meeting at a corner; a diagonal that slips between them
        // must be forbidden, so the path is strictly longer than a naive diagonal.
        let mut grid = CostGrid::uniform(3, 3, 1.0);
        grid.cost[1] = f32::INFINITY; // (1,0)
        grid.cost[3] = f32::INFINITY; // (0,1)
                                      // Going (0,0)->(1,1) diagonally would cut between the two walls; it must
                                      // be blocked, leaving (1,1) reachable only the long way (which here is
                                      // also blocked off), so corner (0,0) cannot reach (1,1) diagonally.
        let path = astar(&grid, (0, 0), (1, 1));
        if let Some(p) = path {
            // If a path exists it must NOT be the illegal single diagonal step.
            assert!(
                !(p.len() == 2 && p[0] == (0, 0) && p[1] == (1, 1)),
                "must not cut the obstacle corner diagonally"
            );
        }
    }
}
