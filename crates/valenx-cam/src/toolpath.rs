//! Toolpath: the canonical output of every CAM operation.
//!
//! A toolpath is a sequence of [`Move`]s. Each move records a single
//! transition from the previous position to a new `position`, with a
//! [`MoveKind`] tag (rapid / cut / plunge) and a feed-rate.
//!
//! - The first move in a toolpath is interpreted as the absolute
//!   starting point — its `kind` is honoured (typically `Rapid`).
//! - All subsequent moves are deltas from the previous position.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::arcfit::ArcDir;

/// Why this move exists. Drives postprocessor formatting (`G0`/`G1`)
/// and the simulation overlay's colour.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum MoveKind {
    /// Non-cutting rapid traverse (G0). Tool is above the stock and
    /// not in contact with material.
    Rapid,
    /// Cutting move along the part surface (G1). Honours `feed`.
    Cut,
    /// Vertical plunge down into material (G1). Honours `feed`
    /// (typically slower than `Cut`).
    Plunge,
    /// Circular-arc cut in the XY plane (G2 / G3). The arc's centre
    /// is fixed and the move's `position` is the **end-point**; the
    /// start-point is the previous move's position. The
    /// post-processor emits `G2`/`G3 X Y I J F` lines (or equivalent
    /// dialect). Inserted by [`crate::arcfit::fit_arcs`].
    Arc {
        /// Arc centre in the XY plane (mm).
        centre_xy: nalgebra::Vector2<f64>,
        /// Arc direction — clockwise (G2) or counter-clockwise (G3).
        dir: ArcDir,
    },
}

impl Eq for MoveKind {}

impl MoveKind {
    /// Short label for panels / debug.
    pub fn label(self) -> &'static str {
        match self {
            MoveKind::Rapid => "Rapid",
            MoveKind::Cut => "Cut",
            MoveKind::Plunge => "Plunge",
            MoveKind::Arc { .. } => "Arc",
        }
    }
}

/// One discrete tool motion. `feed` is mm/min; ignored by `Rapid`
/// moves at postprocessor time.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Move {
    /// What kind of motion this is.
    pub kind: MoveKind,
    /// Absolute target position (mm) in stock-local coordinates.
    pub position: Vector3<f64>,
    /// Feed rate (mm/min). Ignored for `Rapid`.
    pub feed: f64,
}

impl Move {
    /// Convenience constructor — `feed` defaults to 0.
    pub fn new(kind: MoveKind, position: Vector3<f64>, feed: f64) -> Self {
        Self {
            kind,
            position,
            feed,
        }
    }
}

/// Sequence of moves describing a complete CAM operation (or a chain
/// of concatenated operations).
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Toolpath {
    /// All moves, in execution order.
    pub moves: Vec<Move>,
}

impl Toolpath {
    /// Empty toolpath — append with [`Toolpath::push`].
    pub fn new() -> Self {
        Self { moves: Vec::new() }
    }

    /// Push a single move onto the end.
    pub fn push(&mut self, m: Move) {
        self.moves.push(m);
    }

    /// `true` if the toolpath has no moves.
    pub fn is_empty(&self) -> bool {
        self.moves.is_empty()
    }

    /// Number of moves.
    pub fn len(&self) -> usize {
        self.moves.len()
    }

    /// Sum of segment lengths (mm). The first move's start point is
    /// itself — so distance is 0 for a single-move toolpath.
    pub fn total_distance(&self) -> f64 {
        let mut total = 0.0;
        for w in self.moves.windows(2) {
            total += (w[1].position - w[0].position).norm();
        }
        total
    }

    /// Axis-aligned bounding box over every move position. Returns
    /// `None` if the toolpath is empty.
    pub fn bounding_box(&self) -> Option<(Vector3<f64>, Vector3<f64>)> {
        let first = self.moves.first()?.position;
        let mut min = first;
        let mut max = first;
        for m in self.moves.iter().skip(1) {
            min.x = min.x.min(m.position.x);
            min.y = min.y.min(m.position.y);
            min.z = min.z.min(m.position.z);
            max.x = max.x.max(m.position.x);
            max.y = max.y.max(m.position.y);
            max.z = max.z.max(m.position.z);
        }
        Some((min, max))
    }

    /// Append every move from `other` to `self` in order.
    ///
    /// Used by multi-op sequencing — e.g. Face → Pocket → Drill chains
    /// in [`crate::simulate`]'s estimator.
    pub fn concatenate(&mut self, other: &Toolpath) {
        self.moves.extend_from_slice(&other.moves);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    #[test]
    fn empty_toolpath() {
        let t = Toolpath::new();
        assert!(t.is_empty());
        assert_eq!(t.len(), 0);
        assert!(t.bounding_box().is_none());
        assert!((t.total_distance() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn total_distance_sums_segments() {
        let mut t = Toolpath::new();
        t.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 10.0), 0.0));
        t.push(Move::new(MoveKind::Plunge, p(0.0, 0.0, 0.0), 200.0));
        t.push(Move::new(MoveKind::Cut, p(10.0, 0.0, 0.0), 500.0));
        // 10 down + 10 across = 20
        assert!((t.total_distance() - 20.0).abs() < 1e-9);
    }

    #[test]
    fn bounding_box_envelopes_positions() {
        let mut t = Toolpath::new();
        t.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 5.0), 0.0));
        t.push(Move::new(MoveKind::Cut, p(10.0, 20.0, -3.0), 500.0));
        let (min, max) = t.bounding_box().unwrap();
        assert_eq!(min, p(0.0, 0.0, -3.0));
        assert_eq!(max, p(10.0, 20.0, 5.0));
    }

    #[test]
    fn concatenate_appends() {
        let mut a = Toolpath::new();
        a.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 10.0), 0.0));
        let mut b = Toolpath::new();
        b.push(Move::new(MoveKind::Cut, p(5.0, 5.0, 0.0), 200.0));
        b.push(Move::new(MoveKind::Cut, p(10.0, 5.0, 0.0), 200.0));
        a.concatenate(&b);
        assert_eq!(a.len(), 3);
        assert_eq!(a.moves[2].position, p(10.0, 5.0, 0.0));
    }
}
