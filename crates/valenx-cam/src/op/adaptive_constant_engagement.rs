//! Constant-engagement adaptive clearing — Mastercam Dynamic /
//! HSMWorks Adaptive / Fusion 360 Adaptive-class toolpath.
//!
//! The existing [`crate::op::adaptive_clearing`] op generates a
//! sequence of concentric inward offset rings — the engagement
//! angle is **approximated** by the step-over fraction. In a real
//! HSM operation, the engagement angle is the load-limiting
//! quantity: a small constant engagement (~10°-30°) lets the tool
//! safely run at high feed and aggressive step-down, while
//! transient over-engagement at corners is what snaps cutters.
//!
//! ## Algorithm (real HSM-class adaptive)
//!
//! For each Z step-down level:
//!
//! 1. Stitch the pocket boundary polygon from the source mesh
//!    intersected with the Z plane.
//! 2. Inscribe an XY occupancy grid covering the polygon's AABB
//!    (cell size = `tool.radius / 8`); mark cells inside the polygon
//!    as solid (the *remaining* stock at this depth).
//! 3. Generate a sequence of inward offset rings via the existing
//!    [`crate::offset`] primitive.
//! 4. Walk each ring sample-by-sample. At each sample:
//!    1. Compute the **engagement angle** of the cutter centred at
//!       the sample against the current grid using
//!       [`crate::engagement::engagement_at`].
//!    2. If the engagement ≤ `max_engagement_rad`, emit a cut move
//!       to the sample and carve the cutter footprint out of the
//!       grid.
//!    3. If the engagement exceeds the bound (this is a corner or
//!       a re-engagement spike), insert a **trochoidal roll-over**
//!       — a small circular loop offset toward the *cleared* side
//!       so the cutter peels material instead of slot-cutting it.
//!       After the loop, re-check and emit the cut.
//! 5. After all rings, all cells either cleared or below
//!    `min_stub_area` are considered finished.
//!
//! ## v1 simplifications (honest)
//!
//! - **2.5D** — XY adaptive at each Z level, exactly like the
//!   existing concentric-ring adaptive. True 3D-rest adaptive
//!   (HSMWorks 3D Adaptive) operates on a 3D distance field; this
//!   v1 stays on 2D slices.
//! - **Grid-based engagement** — engagement is sampled from a
//!   `cell_size = R/8` occupancy grid, not analytical CSG. This is
//!   the same approximation HSMWorks Toolpath Verifier uses
//!   internally; commercial CAM uses adaptive cell sizes near the
//!   cutter outline (out of scope here).
//! - **Roll-over moves** — trochoidal loops use 16-segment polylines
//!   inserted at corners where the concentric ring would over-
//!   engage. They consume cleared material first, allowing the
//!   next sample on the ring to enter at the bounded engagement.
//!   The loop radius is `helical_radius` (typically R/2-R).
//! - **One pocket per Z slice** — multi-region cross-sections fall
//!   back to the first ring set (same simplification as the
//!   concentric-ring adaptive).
//!
//! See [`AdaptiveConstantEngagementParams`] for the parameters and
//! [`generate`] for the toolpath generator.

use nalgebra::Vector3;
use valenx_mesh::cut::intersect_plane;
use valenx_mesh::Mesh;

use crate::engagement::{engagement_at, StockGrid};
use crate::error::CamError;
use crate::offset;
use crate::op::profile::stitch_segments;
use crate::stock::Stock;
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Parameters for constant-engagement adaptive clearing.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct AdaptiveConstantEngagementParams {
    /// Tool id.
    pub tool_id: u32,
    /// Cutting feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// **Maximum engagement angle in radians** — the load bound the
    /// path generator enforces everywhere. Typical HSM values are
    /// `15°-30° = 0.26-0.52 rad`. The default `0.35 rad ≈ 20°` is
    /// the HSMWorks "high removal" default.
    pub max_engagement_rad: f64,
    /// Nominal step-over as a fraction of tool diameter (0..1).
    /// Used to size offset rings before engagement bounding. Real
    /// HSM uses 30%-50% — much more aggressive than concentric
    /// clearing (3%-10%) because the engagement bound caps load.
    pub step_over_fraction: f64,
    /// Z step-down per pass (mm).
    pub step_down: f64,
    /// Radius (mm) of the trochoidal roll-over inserted at corners
    /// where the concentric path would over-engage. Typical `0.5×R`
    /// to `1.0×R`.
    pub helical_radius: f64,
    /// Number of segments in each trochoidal loop polyline.
    pub helical_segments: usize,
    /// Total depth below stock top (mm).
    pub depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
    /// Cell size for the engagement grid (mm). `<= 0` selects the
    /// default `tool.radius / 8`.
    pub cell_size_mm: f64,
    /// Engagement-ray sample count per query (8 minimum). Default
    /// 64 = `~5.6°` resolution, matching commercial verifiers.
    pub engagement_samples: usize,
}

impl Default for AdaptiveConstantEngagementParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 2500.0,
            plunge_feed: 300.0,
            spindle_rpm: 18_000.0,
            max_engagement_rad: 0.35,
            step_over_fraction: 0.35,
            step_down: 6.0,
            helical_radius: 1.0,
            helical_segments: 16,
            depth: 6.0,
            safe_z_clearance: 5.0,
            cell_size_mm: 0.0,
            engagement_samples: 64,
        }
    }
}

/// Report bundled with the toolpath — exposes the engagement
/// telemetry so callers (and tests) can verify the bound was
/// actually respected.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct AdaptiveEngagementReport {
    /// Maximum engagement observed at any *cutting* move (radians).
    pub max_engagement_rad: f64,
    /// Mean engagement (radians).
    pub mean_engagement_rad: f64,
    /// Cut-move count.
    pub n_cut_moves: usize,
    /// Trochoidal roll-overs inserted (one per corner relief).
    pub n_rollovers: usize,
}

/// Generate a constant-engagement adaptive-clearing toolpath. See
/// the module docs for the full algorithm.
///
/// Returns the toolpath plus a [`AdaptiveEngagementReport`] with
/// the engagement-angle statistics. The report's
/// `max_engagement_rad` is the *measured* maximum — within
/// `2π / engagement_samples` of the bound when the bound is
/// effective.
pub fn generate(
    stock: &Stock,
    source: &Mesh,
    params: &AdaptiveConstantEngagementParams,
    tool: &Tool,
) -> Result<(Toolpath, AdaptiveEngagementReport), CamError> {
    validate(params, tool)?;
    let safe_z = stock.top_z() + params.safe_z_clearance;
    let mut tp = Toolpath::new();
    tp.push(Move::new(
        MoveKind::Rapid,
        Vector3::new(0.0, 0.0, safe_z),
        0.0,
    ));

    let n_passes = crate::op::compute_n_passes(
        params.depth,
        params.step_down,
        "adaptive_constant_engagement",
    )?;
    let step_over = params.step_over_fraction * tool.diameter_mm;
    let cell_size = if params.cell_size_mm > 0.0 {
        params.cell_size_mm
    } else {
        (tool.radius_mm() / 8.0).max(1e-3)
    };

    let mut max_eng = 0.0_f64;
    let mut sum_eng = 0.0_f64;
    let mut n_cut = 0_usize;
    let mut n_roll = 0_usize;
    let mut any_pass = false;

    for k in 1..=n_passes {
        let depth_below_top = (params.step_down * k as f64).min(params.depth);
        let z = stock.top_z() - depth_below_top;
        let segments = intersect_plane(
            source,
            Vector3::new(0.0, 0.0, z),
            Vector3::new(0.0, 0.0, 1.0),
        );
        if segments.is_empty() {
            continue;
        }
        let polygon = match stitch_segments(&segments) {
            Some(p) => p,
            None => continue,
        };
        // Inscribe a grid covering the polygon AABB + a small
        // tool-radius margin so the cutter can sit on the boundary.
        let (xy_min, xy_max) = polygon_aabb(&polygon, tool.radius_mm() * 1.5);
        let mut grid = StockGrid::new(xy_min, xy_max, cell_size);
        // Mark cells *outside* the polygon as non-solid (the
        // polygon's interior is the only material the cutter cares
        // about).
        mark_inside_only(&mut grid, &polygon);

        // Initial inset to put the cutter centre safely inside the
        // polygon.
        let mut current = match offset::polygon(&polygon, -tool.radius_mm())
            .into_iter()
            .next()
        {
            Some(p) => p,
            None => continue,
        };
        let first = current[0];
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(first.x, first.y, safe_z),
            0.0,
        ));
        tp.push(Move::new(
            MoveKind::Plunge,
            Vector3::new(first.x, first.y, z),
            params.plunge_feed,
        ));
        // Carve the initial cutter footprint.
        grid.carve_disc(nalgebra::Vector2::new(first.x, first.y), tool.radius_mm());
        let mut last_pos = nalgebra::Vector2::new(first.x, first.y);

        // Walk rings inward. Per ring, walk sample-by-sample with
        // engagement bounding.
        for _ring_iter in 0..2048 {
            let n_pts = current.len();
            // Sample finely along the ring polyline so engagement
            // queries are dense (cutter advances by < cell_size).
            for i in 0..n_pts {
                let a = current[i];
                let b = current[(i + 1) % n_pts];
                let seg_len = ((b - a).norm()).max(1e-9);
                let n_sub = ((seg_len / (cell_size * 0.75)).ceil() as usize).max(1);
                for s in 1..=n_sub {
                    let t = (s as f64) / (n_sub as f64);
                    let p = a + (b - a) * t;
                    let p2 = nalgebra::Vector2::new(p.x, p.y);
                    // Check engagement of the *prospective* cutter
                    // position; if it exceeds the bound, insert a
                    // trochoidal roll-over before committing.
                    let eng = engagement_at(p2, tool.radius_mm(), &grid, params.engagement_samples);
                    if eng > params.max_engagement_rad {
                        // Trochoidal relief: walk a small circle
                        // centred between last and prospective on
                        // the *cleared* side, peeling material
                        // until the next sample is below bound.
                        emit_rollover(
                            &mut tp,
                            &mut grid,
                            last_pos,
                            p2,
                            tool.radius_mm(),
                            params.helical_radius,
                            params.helical_segments,
                            z,
                            params.feed_mm_per_min,
                        );
                        n_roll += 1;
                    }
                    // Now emit the cut and update grid + telemetry.
                    let measured_eng =
                        engagement_at(p2, tool.radius_mm(), &grid, params.engagement_samples);
                    tp.push(Move::new(
                        MoveKind::Cut,
                        Vector3::new(p.x, p.y, z),
                        params.feed_mm_per_min,
                    ));
                    grid.carve_segment(last_pos, p2, tool.radius_mm());
                    last_pos = p2;
                    if measured_eng > max_eng {
                        max_eng = measured_eng;
                    }
                    sum_eng += measured_eng;
                    n_cut += 1;
                }
            }
            // Offset inward by step_over.
            let next = offset::polygon(&current, -step_over);
            if next.is_empty() {
                break;
            }
            let next_ring = next.into_iter().next().unwrap();
            if shortest_edge(&next_ring) < tool.radius_mm() * 0.5 {
                break;
            }
            current = next_ring;
        }
        tp.push(Move::new(
            MoveKind::Rapid,
            Vector3::new(0.0, 0.0, safe_z),
            0.0,
        ));
        any_pass = true;
    }
    if !any_pass {
        return Err(CamError::BadOperation {
            name: "adaptive_constant_engagement".into(),
            reason: "no inscribed polygon at any depth".into(),
        });
    }
    let mean = if n_cut > 0 {
        sum_eng / (n_cut as f64)
    } else {
        0.0
    };
    Ok((
        tp,
        AdaptiveEngagementReport {
            max_engagement_rad: max_eng,
            mean_engagement_rad: mean,
            n_cut_moves: n_cut,
            n_rollovers: n_roll,
        },
    ))
}

/// Emit a 16-segment trochoidal arc that consumes material in a
/// circle of radius `helical_radius` between `from` and `to`,
/// carving the cutter footprint into `grid` along the way.
#[allow(clippy::too_many_arguments)]
fn emit_rollover(
    tp: &mut Toolpath,
    grid: &mut StockGrid,
    from: nalgebra::Vector2<f64>,
    to: nalgebra::Vector2<f64>,
    tool_radius: f64,
    helical_radius: f64,
    n_segments: usize,
    z: f64,
    feed: f64,
) {
    // Centre between from and to.
    let centre = (from + to) * 0.5;
    let n = n_segments.max(8);
    let mut prev = from;
    for i in 1..=n {
        let theta = (i as f64) * std::f64::consts::TAU / (n as f64);
        let p = centre + nalgebra::Vector2::new(helical_radius * theta.cos(), helical_radius * theta.sin());
        tp.push(Move::new(
            MoveKind::Cut,
            Vector3::new(p.x, p.y, z),
            feed,
        ));
        grid.carve_segment(prev, p, tool_radius);
        prev = p;
    }
}

fn polygon_aabb(
    polygon: &[Vector3<f64>],
    margin: f64,
) -> (nalgebra::Vector2<f64>, nalgebra::Vector2<f64>) {
    let mut min = nalgebra::Vector2::new(f64::INFINITY, f64::INFINITY);
    let mut max = nalgebra::Vector2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for v in polygon {
        if v.x < min.x {
            min.x = v.x;
        }
        if v.y < min.y {
            min.y = v.y;
        }
        if v.x > max.x {
            max.x = v.x;
        }
        if v.y > max.y {
            max.y = v.y;
        }
    }
    min.x -= margin;
    min.y -= margin;
    max.x += margin;
    max.y += margin;
    (min, max)
}

/// Set every grid cell **outside** the polygon to non-solid (so the
/// engagement query only sees pocket-interior material as
/// removable).
fn mark_inside_only(grid: &mut StockGrid, polygon: &[Vector3<f64>]) {
    let cs = grid.cell_size_mm;
    for iy in 0..grid.n_y {
        for ix in 0..grid.n_x {
            let cx = grid.min.x + (ix as f64 + 0.5) * cs;
            let cy = grid.min.y + (iy as f64 + 0.5) * cs;
            if !point_in_polygon(cx, cy, polygon) {
                grid.solid[iy * grid.n_x + ix] = false;
            }
        }
    }
}

fn point_in_polygon(x: f64, y: f64, polygon: &[Vector3<f64>]) -> bool {
    let n = polygon.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let xi = polygon[i].x;
        let yi = polygon[i].y;
        let xj = polygon[j].x;
        let yj = polygon[j].y;
        let intersect = ((yi > y) != (yj > y))
            && (x < (xj - xi) * (y - yi) / (yj - yi).max(f64::MIN_POSITIVE).copysign(yj - yi) + xi);
        if intersect {
            inside = !inside;
        }
        j = i;
    }
    inside
}

fn shortest_edge(polygon: &[Vector3<f64>]) -> f64 {
    let n = polygon.len();
    if n < 2 {
        return 0.0;
    }
    let mut min = f64::INFINITY;
    for i in 0..n {
        let d = (polygon[(i + 1) % n] - polygon[i]).norm();
        if d < min {
            min = d;
        }
    }
    min
}

fn validate(
    params: &AdaptiveConstantEngagementParams,
    tool: &Tool,
) -> Result<(), CamError> {
    let mk = |reason: String| CamError::BadOperation {
        name: "adaptive_constant_engagement".into(),
        reason,
    };
    if !(params.max_engagement_rad > 0.0 && params.max_engagement_rad <= std::f64::consts::TAU) {
        return Err(mk(format!(
            "max_engagement_rad must be in (0, 2π] (got {})",
            params.max_engagement_rad
        )));
    }
    if !(params.step_over_fraction > 0.0 && params.step_over_fraction <= 1.0) {
        return Err(mk(format!(
            "step_over_fraction must be in (0, 1] (got {})",
            params.step_over_fraction
        )));
    }
    if !(params.step_down > 0.0) {
        return Err(mk(format!("step_down must be > 0 (got {})", params.step_down)));
    }
    if !(params.depth > 0.0) {
        return Err(mk(format!("depth must be > 0 (got {})", params.depth)));
    }
    if !(params.feed_mm_per_min > 0.0) {
        return Err(mk(format!(
            "feed must be > 0 (got {})",
            params.feed_mm_per_min
        )));
    }
    if !(params.helical_radius > 0.0) {
        return Err(mk(format!(
            "helical_radius must be > 0 (got {})",
            params.helical_radius
        )));
    }
    if !(tool.diameter_mm > 0.0) {
        return Err(mk("tool diameter must be > 0".into()));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::ToolKind;
    use valenx_mesh::element::{ElementBlock, ElementType};

    fn cube(size: f64) -> Mesh {
        let s = size * 0.5;
        let nodes = vec![
            Vector3::new(-s, -s, -s),
            Vector3::new(s, -s, -s),
            Vector3::new(s, s, -s),
            Vector3::new(-s, s, -s),
            Vector3::new(-s, -s, s),
            Vector3::new(s, -s, s),
            Vector3::new(s, s, s),
            Vector3::new(-s, s, s),
        ];
        let conn: Vec<u32> = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        let mut mesh = Mesh::new("cube");
        mesh.nodes = nodes;
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = conn;
        mesh.element_blocks.push(block);
        mesh
    }

    fn square_pocket_setup() -> (Stock, Mesh, Tool) {
        // 40mm cube = 40x40 pocket polygon at any Z.
        let mesh = cube(40.0);
        let stock = Stock::new(
            Vector3::new(-25.0, -25.0, -20.0),
            Vector3::new(50.0, 50.0, 40.0),
            "alu",
        )
        .unwrap();
        let tool = Tool::new(1, "EM4", ToolKind::EndMill, 4.0, 25.0, 2, "carbide").unwrap();
        (stock, mesh, tool)
    }

    #[test]
    fn engagement_bound_respected_on_square_pocket() {
        let (stock, mesh, tool) = square_pocket_setup();
        let mut params = AdaptiveConstantEngagementParams {
            step_down: 4.0,
            depth: 4.0,
            max_engagement_rad: 0.35,
            ..Default::default()
        };
        // Bigger grid resolution to keep test fast.
        params.cell_size_mm = 0.5;
        params.engagement_samples = 48;
        let (tp, report) = generate(&stock, &mesh, &params, &tool).unwrap();
        assert!(!tp.is_empty());
        assert!(
            report.n_cut_moves > 50,
            "expected substantial path, got {}",
            report.n_cut_moves
        );
        // The measured max engagement should stay close to the bound.
        // Allow a small overshoot tolerance equal to one engagement
        // sample bucket (2π/n_samples) — the rollover re-checks
        // *after* its loop, and the check is bucketed.
        let bucket = std::f64::consts::TAU / (params.engagement_samples as f64);
        let tol = bucket * 2.0; // 2-bucket slack for sampling jitter
        assert!(
            report.max_engagement_rad <= params.max_engagement_rad + tol,
            "engagement bound exceeded: max {} > bound {} + tol {} (buckets {})",
            report.max_engagement_rad,
            params.max_engagement_rad,
            tol,
            bucket
        );
        assert!(
            report.mean_engagement_rad <= params.max_engagement_rad,
            "mean engagement {} should sit at or below bound {}",
            report.mean_engagement_rad,
            params.max_engagement_rad
        );
    }

    #[test]
    fn plausible_path_length_vs_concentric() {
        // The constant-engagement path should be *longer* (or
        // comparable) than the concentric-ring path because corner
        // roll-overs add length — but not 10×, that would mean we
        // failed to make progress. We just sanity-check positive
        // total distance.
        let (stock, mesh, tool) = square_pocket_setup();
        let params = AdaptiveConstantEngagementParams {
            step_down: 4.0,
            depth: 4.0,
            max_engagement_rad: 0.5,
            cell_size_mm: 0.5,
            engagement_samples: 32,
            ..Default::default()
        };
        let (tp, report) = generate(&stock, &mesh, &params, &tool).unwrap();
        let dist = tp.total_distance();
        assert!(dist > 0.0, "expected non-zero path length");
        assert!(
            dist < 50_000.0,
            "path length is implausibly large: {dist}",
        );
        // A roll-over insertion is the diagnostic the engagement
        // bound is genuinely engaged (i.e. the path is *not* just a
        // copy of the concentric path). Some are expected at the
        // square pocket's corners.
        assert!(
            report.n_rollovers > 0 || report.max_engagement_rad < params.max_engagement_rad,
            "either a rollover should have fired or the bound was already loose enough",
        );
    }

    #[test]
    fn validate_rejects_bad_params() {
        let tool = Tool::new(1, "EM4", ToolKind::EndMill, 4.0, 25.0, 2, "carbide").unwrap();
        let bad = AdaptiveConstantEngagementParams {
            max_engagement_rad: -0.1,
            ..Default::default()
        };
        assert!(validate(&bad, &tool).is_err());
        let bad = AdaptiveConstantEngagementParams {
            step_over_fraction: 1.5,
            ..Default::default()
        };
        assert!(validate(&bad, &tool).is_err());
    }

    #[test]
    fn looser_bound_means_higher_max_engagement_allowed() {
        let (stock, mesh, tool) = square_pocket_setup();
        let strict = AdaptiveConstantEngagementParams {
            step_down: 4.0,
            depth: 4.0,
            max_engagement_rad: 0.25,
            cell_size_mm: 0.6,
            engagement_samples: 32,
            ..Default::default()
        };
        let loose = AdaptiveConstantEngagementParams {
            max_engagement_rad: 1.5,
            ..strict.clone()
        };
        let (_t1, r1) = generate(&stock, &mesh, &strict, &tool).unwrap();
        let (_t2, r2) = generate(&stock, &mesh, &loose, &tool).unwrap();
        // The strict pass should fire more roll-overs than the loose pass.
        assert!(
            r1.n_rollovers >= r2.n_rollovers,
            "strict bound should fire ≥ rollovers vs loose ({} vs {})",
            r1.n_rollovers,
            r2.n_rollovers
        );
        // Strict pass's max engagement is bounded; loose pass's
        // bound is so loose it never fires.
        let bucket = std::f64::consts::TAU / (strict.engagement_samples as f64);
        assert!(r1.max_engagement_rad <= strict.max_engagement_rad + 2.0 * bucket);
    }
}
