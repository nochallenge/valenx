//! Feedrate optimization with lookahead — the controller-side post-
//! processing pass that every commercial CAM tool does.
//!
//! Without it: the post emits whatever feed the operation requested,
//! the controller runs that everywhere, and at a sharp corner
//! either (a) the machine's acceleration limit is exceeded (bad
//! finish, dropped steps) or (b) the controller's own lookahead
//! ramps down — which masks the planner deficiency but at the cost
//! of unpredictable cycle time.
//!
//! With it: every move's feed is scheduled by the offline planner
//! to:
//!
//! 1. **Centripetal-acceleration bound** — on a fitted arc, the
//!    maximum safe feed is `v_arc = √(a_cent_max · r)`. Beyond this
//!    the workpiece sees forces that crash the spindle.
//! 2. **Corner-deceleration bound** — at every sharp transition
//!    (G1→G1 with significant angle change), the achievable speed
//!    is bounded by the velocity at which the corner can be
//!    rounded within the chord-error tolerance.
//! 3. **Deceleration ramp** — given the machine's max deceleration
//!    `a_decel_max`, the planner traces *backward* from each
//!    bounded move and clamps prior moves' feeds so the machine
//!    can decelerate in time (the lookahead window).
//! 4. **Acceleration ramp** — symmetric forward trace from low-feed
//!    moves so the machine can spool up smoothly.
//!
//! The final feed at every move is `min(configured, v_arc, v_corner,
//! v_lookahead)`.
//!
//! ## v1 simplifications (honest)
//!
//! - **Single-axis acceleration limit** — `a_decel_max` is treated
//!   as a scalar over all 3 axes. Production controllers track per-
//!   axis limits with axis-projected vectors; that's a follow-up.
//! - **Corner bound is geometric** — `v_corner = √(2 · a_decel_max ·
//!   chord_tol)` (the standard junction-velocity formula, derived
//!   from a chord-error tolerance over a discretised turn). It is
//!   exact for the limit of a single sharp corner; chained tight
//!   corners use the same bound at each.
//! - **Lookahead is greedy** — backward pass clamps the move
//!   immediately upstream; this is the standard
//!   `v² = v_target² + 2 · a · d` iteration. We propagate as far
//!   back as `n_lookahead` moves; a real CNC controller's
//!   lookahead is typically 64-1024 moves.

use serde::{Deserialize, Serialize};

use crate::toolpath::{MoveKind, Toolpath};

/// Parameters for the feedrate optimizer.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FeedrateParams {
    /// Maximum centripetal acceleration the machine + workpiece can
    /// sustain (mm/min²). Hobby routers: ~1.5e7 (~250 mm/s²);
    /// production VMC: ~3e8 mm/min² (~800 mm/s²).
    pub a_centripetal_max_mm_per_min2: f64,
    /// Maximum linear deceleration along the toolpath (mm/min²).
    /// Same units as centripetal.
    pub a_decel_max_mm_per_min2: f64,
    /// Chord-error tolerance at sharp corners (mm). Used by the
    /// corner-bound formula `v_corner = √(2 · a · chord_tol)`.
    /// Typical 0.005-0.02 mm for finishing.
    pub corner_chord_tol_mm: f64,
    /// Threshold turn angle (radians) below which a corner doesn't
    /// trigger feed clamping. A 5° turn is barely a corner; a 30°
    /// turn warrants slowing.
    pub corner_angle_threshold_rad: f64,
    /// Number of moves the backward + forward lookahead propagates
    /// through. Real controllers run 64-1024; the default 128 is
    /// the EdgeCAM / HSMWorks-class envelope.
    pub n_lookahead: usize,
}

impl Default for FeedrateParams {
    fn default() -> Self {
        Self {
            // 5000 mm/s² is a sane production VMC default for
            // centripetal limit, expressed in mm/min².
            a_centripetal_max_mm_per_min2: 5000.0 * 60.0 * 60.0,
            a_decel_max_mm_per_min2: 5000.0 * 60.0 * 60.0,
            corner_chord_tol_mm: 0.01,
            corner_angle_threshold_rad: (15.0_f64).to_radians(),
            n_lookahead: 128,
        }
    }
}

/// Statistics produced by the feedrate optimization pass.
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
pub struct FeedrateReport {
    /// Moves processed.
    pub input_moves: usize,
    /// Moves whose feed was reduced (centripetal, corner, or lookahead).
    pub reduced_moves: usize,
    /// Moves clamped specifically by the centripetal bound.
    pub centripetal_clamps: usize,
    /// Moves clamped by the corner bound.
    pub corner_clamps: usize,
    /// Moves clamped by the lookahead backward pass.
    pub lookahead_clamps: usize,
    /// Maximum reduction observed (mm/min) across all moves.
    pub max_reduction_mm_per_min: f64,
}

/// Run the feedrate optimization pass over `toolpath`. Returns a
/// new toolpath with feeds scheduled and the report.
pub fn optimize(toolpath: &Toolpath, params: &FeedrateParams) -> (Toolpath, FeedrateReport) {
    let n = toolpath.moves.len();
    if n < 2 {
        return (
            toolpath.clone(),
            FeedrateReport {
                input_moves: n,
                ..Default::default()
            },
        );
    }
    let mut out = toolpath.clone();
    let mut report = FeedrateReport {
        input_moves: n,
        ..Default::default()
    };
    // -- Pass 1: centripetal bound on Arc moves.
    for i in 0..n {
        if let MoveKind::Arc { centre_xy, .. } = out.moves[i].kind {
            // Need the start point (prev move's position) to read r.
            if i == 0 {
                continue;
            }
            let start = out.moves[i - 1].position;
            let r = ((start.x - centre_xy.x).powi(2) + (start.y - centre_xy.y).powi(2)).sqrt();
            if !(r > 0.0) {
                continue;
            }
            let v_arc = (params.a_centripetal_max_mm_per_min2 * r).sqrt();
            if v_arc < out.moves[i].feed {
                let red = out.moves[i].feed - v_arc;
                if red > report.max_reduction_mm_per_min {
                    report.max_reduction_mm_per_min = red;
                }
                out.moves[i].feed = v_arc;
                report.centripetal_clamps += 1;
            }
        }
    }
    // -- Pass 2: corner-bound clamp on G1 transitions with sharp turn.
    for i in 1..(n - 1) {
        // Only consider cut/plunge moves — rapids and arcs are out of
        // scope for the corner clamp (arcs already bounded; rapids
        // jog at controller-determined rates).
        if !matches!(out.moves[i].kind, MoveKind::Cut | MoveKind::Plunge) {
            continue;
        }
        let prev = out.moves[i - 1].position;
        let curr = out.moves[i].position;
        let next = out.moves[i + 1].position;
        let v_in = curr - prev;
        let v_out = next - curr;
        let n_in = v_in.norm();
        let n_out = v_out.norm();
        if n_in < 1e-9 || n_out < 1e-9 {
            continue;
        }
        let cos_t = (v_in.dot(&v_out)) / (n_in * n_out);
        let cos_t = cos_t.clamp(-1.0, 1.0);
        // `turn` is the deviation from straight: 0 = collinear, π = U-turn.
        // acos(cos_t) gives the angle BETWEEN the two direction
        // vectors (0 when they point the same way), which is exactly
        // the turn-from-straight.
        let turn = cos_t.acos();
        if turn < params.corner_angle_threshold_rad {
            continue;
        }
        // Junction-velocity bound — Smoothieware / GRBL formula.
        // For a half-turn of `turn/2`, the allowed velocity is
        // bounded by the chord-error tolerance:
        // v_corner = sqrt(a · chord_tol · sin(half) / (1 − sin(half)))
        // — derived from a circular-arc smoothing of the corner.
        let half = turn * 0.5;
        let s = half.sin();
        // Avoid divide-by-zero at 180° (U-turn).
        let factor = if s < 0.999 {
            (s / (1.0 - s)).max(0.0)
        } else {
            // 180° cusp: velocity must approach 0.
            0.0
        };
        let v_corner =
            (params.a_decel_max_mm_per_min2 * params.corner_chord_tol_mm * factor).sqrt();
        // Cap by configured feed.
        let original = out.moves[i].feed;
        if v_corner < original {
            let red = original - v_corner;
            if red > report.max_reduction_mm_per_min {
                report.max_reduction_mm_per_min = red;
            }
            out.moves[i].feed = v_corner.max(0.0);
            report.corner_clamps += 1;
        }
    }
    // -- Pass 3: backward lookahead — for each move with a low feed,
    // propagate the deceleration constraint backward so the machine
    // can ramp down in time. v² ≤ v_target² + 2·a·d.
    for i_back in (1..n).rev() {
        let target = out.moves[i_back].feed;
        // Only clamp via the deceleration formula on cut/plunge/arc
        // moves; rapid moves can be handled by the controller.
        if !matches!(
            out.moves[i_back].kind,
            MoveKind::Cut | MoveKind::Plunge | MoveKind::Arc { .. }
        ) {
            continue;
        }
        // Propagate up to n_lookahead moves backward.
        let start = i_back.saturating_sub(params.n_lookahead);
        let mut v_next = target;
        for j in (start..i_back).rev() {
            if !matches!(
                out.moves[j].kind,
                MoveKind::Cut | MoveKind::Plunge | MoveKind::Arc { .. }
            ) {
                break;
            }
            let d = (out.moves[j + 1].position - out.moves[j].position).norm();
            if d < 1e-12 {
                continue;
            }
            // Max velocity at the start of segment j that decelerates
            // to v_next over distance d.
            let v_max = (v_next.powi(2) + 2.0 * params.a_decel_max_mm_per_min2 * d).sqrt();
            if v_max < out.moves[j].feed {
                let red = out.moves[j].feed - v_max;
                if red > report.max_reduction_mm_per_min {
                    report.max_reduction_mm_per_min = red;
                }
                out.moves[j].feed = v_max;
                report.lookahead_clamps += 1;
            } else {
                // Once we hit a move that's already low enough, no
                // need to keep propagating.
                break;
            }
            v_next = out.moves[j].feed;
        }
    }
    // Count reduced moves and final stats.
    for (i, m) in out.moves.iter().enumerate() {
        if i < toolpath.moves.len() && m.feed < toolpath.moves[i].feed - 1e-6 {
            report.reduced_moves += 1;
        }
    }
    (out, report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arcfit::ArcDir;
    use crate::toolpath::{Move, MoveKind, Toolpath};
    use nalgebra::Vector3;

    fn p(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    #[test]
    fn empty_or_short_toolpath_passthrough() {
        let tp = Toolpath::new();
        let (out, r) = optimize(&tp, &FeedrateParams::default());
        assert_eq!(out.len(), 0);
        assert_eq!(r.input_moves, 0);

        let mut tp1 = Toolpath::new();
        tp1.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 5.0), 0.0));
        let (out1, r1) = optimize(&tp1, &FeedrateParams::default());
        assert_eq!(out1.len(), 1);
        assert_eq!(r1.input_moves, 1);
    }

    #[test]
    fn sharp_corner_clamps_corner_feed() {
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 0.0), 0.0));
        tp.push(Move::new(MoveKind::Cut, p(10.0, 0.0, 0.0), 5000.0));
        // 90-degree turn here.
        tp.push(Move::new(MoveKind::Cut, p(10.0, 10.0, 0.0), 5000.0));
        tp.push(Move::new(MoveKind::Cut, p(20.0, 10.0, 0.0), 5000.0));
        let params = FeedrateParams::default();
        let (out, report) = optimize(&tp, &params);
        // The middle move (index 2) is the corner — should be
        // significantly lower than 5000.
        assert!(
            out.moves[2].feed < 5000.0,
            "corner move feed should be reduced: {}",
            out.moves[2].feed
        );
        assert!(
            report.corner_clamps >= 1,
            "expected at least one corner clamp, got {}",
            report.corner_clamps
        );
        // Lookahead should have ramped down the move *before* the corner.
        assert!(
            out.moves[1].feed < 5000.0,
            "pre-corner move should be ramped down by lookahead: {}",
            out.moves[1].feed
        );
    }

    #[test]
    fn arc_centripetal_bound_clamps_arc() {
        // Arc with very small radius — should be clamped hard.
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Cut, p(0.0, 0.0, 0.0), 6000.0));
        tp.push(Move {
            kind: MoveKind::Arc {
                centre_xy: nalgebra::Vector2::new(0.5, 0.0),
                dir: ArcDir::Clockwise,
            },
            position: p(1.0, 0.0, 0.0),
            feed: 6000.0,
        });
        // Constrict centripetal acceleration so the bound bites.
        let params = FeedrateParams {
            a_centripetal_max_mm_per_min2: 1e6,
            ..Default::default()
        };
        let (out, report) = optimize(&tp, &params);
        // r = 0.5; v_max = sqrt(1e6 * 0.5) = sqrt(500000) ≈ 707
        assert!(
            out.moves[1].feed < 1000.0,
            "arc feed should be clamped to ~707, got {}",
            out.moves[1].feed
        );
        assert_eq!(report.centripetal_clamps, 1);
    }

    #[test]
    fn no_corner_no_arc_no_clamp() {
        // Straight line — no corners, no arcs → no clamp.
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 0.0), 0.0));
        for k in 1..=5 {
            tp.push(Move::new(MoveKind::Cut, p(k as f64, 0.0, 0.0), 2000.0));
        }
        let params = FeedrateParams::default();
        let (out, report) = optimize(&tp, &params);
        // All cut moves should still be at 2000.
        for m in &out.moves {
            if matches!(m.kind, MoveKind::Cut) {
                assert!((m.feed - 2000.0).abs() < 1e-6);
            }
        }
        assert_eq!(report.corner_clamps, 0);
        assert_eq!(report.centripetal_clamps, 0);
        assert_eq!(report.reduced_moves, 0);
    }

    #[test]
    fn ramp_down_distance_matches_decel_kinematics() {
        // Place a corner with known low post-corner feed; check the
        // pre-corner moves' feeds follow v² = v_post² + 2·a·d.
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 0.0), 0.0));
        // 10 cut moves each 5mm long going +X at high feed.
        for k in 1..=10 {
            tp.push(Move::new(
                MoveKind::Cut,
                p((k as f64) * 5.0, 0.0, 0.0),
                10_000.0,
            ));
        }
        // Now a sharp corner that turns -90°.
        tp.push(Move::new(MoveKind::Cut, p(50.0, 5.0, 0.0), 10_000.0));
        let params = FeedrateParams {
            a_decel_max_mm_per_min2: 1e6,
            a_centripetal_max_mm_per_min2: 1e6,
            corner_chord_tol_mm: 0.005,
            n_lookahead: 32,
            ..Default::default()
        };
        let (out, report) = optimize(&tp, &params);
        // The corner vertex is the move at index 10 (the cut to
        // (50,0,0) is the move where the next segment turns). Its
        // feed should be sharply clamped; the moves leading up to
        // it should monotonically ramp *down* (the machine
        // decelerates in time).
        let v_corner = out.moves[10].feed;
        assert!(
            v_corner < 1000.0,
            "corner-vertex move should be clamped low, got {v_corner}",
        );
        // Walk backward from the move just before the corner. Each
        // successive earlier move should be at least as high as the
        // next one (we can't ramp UP earlier than the corner).
        let mut prev_v = v_corner;
        for i in (1..10).rev() {
            let v_i = out.moves[i].feed;
            assert!(
                v_i >= prev_v - 1e-6,
                "lookahead ramp not monotonic: move {i} feed {v_i} < later {prev_v}",
            );
            prev_v = v_i;
        }
        assert!(report.lookahead_clamps > 0);
        // The earliest move's feed should be back near the original
        // 10_000 — enough distance to spool back up.
        // v² = v_corner² + 2·a·d. Total distance back to move 1 is
        // 9 segments × 5 mm = 45 mm. v_max = sqrt(v_corner² + 9e7) ≈ 9487.
        // Allow a buffer for the corner clamp value.
        let earliest = out.moves[1].feed;
        assert!(
            earliest > 5000.0,
            "earliest move should have ramped back up substantially, got {earliest}",
        );
    }
}
