//! Fixed-axis indexing — set rotary axes to a fixed angle, then run a
//! 3-axis op on the resulting tilted plane.
//!
//! v1 produces a 5-axis toolpath by:
//!
//! 1. Emitting one positioning move at `(a_deg, b_deg)` (rapid).
//! 2. Lifting every 3-axis move from the source toolpath into the
//!    5-axis toolpath at the same `(a_deg, b_deg)` orientation.
//!
//! No interpolation of the rotary axes — a single index pose covers
//! the entire 3-axis op. Use [`crate::op::tcp_5ax_contour`] for true
//! continuous 5-axis machining.

use serde::{Deserialize, Serialize};

use crate::axis::{Move5Ax, Toolpath5Ax};
use crate::toolpath::Toolpath;

/// Parameters for a fixed-axis indexed op.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FixedAxisIndexingParams {
    /// A-axis index (degrees).
    pub a_deg: f64,
    /// B-axis index (degrees).
    pub b_deg: f64,
}

/// Lift a 3-axis [`Toolpath`] into a [`Toolpath5Ax`] at the supplied
/// rotary pose.
pub fn lift_to_indexed_pose(source: &Toolpath, params: &FixedAxisIndexingParams) -> Toolpath5Ax {
    let mut tp = Toolpath5Ax::new();
    for m in &source.moves {
        tp.push(Move5Ax::new(
            m.kind,
            m.position,
            params.a_deg,
            params.b_deg,
            m.feed,
        ));
    }
    tp
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::toolpath::{Move, MoveKind};
    use nalgebra::Vector3;

    #[test]
    fn lift_preserves_move_count_and_pose() {
        let mut t = Toolpath::new();
        t.push(Move::new(MoveKind::Rapid, Vector3::new(0.0, 0.0, 5.0), 0.0));
        t.push(Move::new(
            MoveKind::Cut,
            Vector3::new(10.0, 0.0, 0.0),
            500.0,
        ));
        let params = FixedAxisIndexingParams {
            a_deg: 30.0,
            b_deg: 45.0,
        };
        let tp5 = lift_to_indexed_pose(&t, &params);
        assert_eq!(tp5.len(), 2);
        assert_eq!(tp5.moves[0].a_deg, 30.0);
        assert_eq!(tp5.moves[0].b_deg, 45.0);
    }
}
