//! Geometric mate constraints between two parts.
//!
//! Each [`Mate`] contributes one or more scalar *residuals* to the
//! solver's system. The solver drives all residuals to zero
//! simultaneously to find a transform configuration that satisfies
//! every mate. See [`crate::solver`] for the Newton-Raphson driver.
//!
//! ## Residual conventions
//!
//! - **Coincident** — 3 residuals: `(R_b·p_b + t_b) - (R_a·p_a + t_a)`,
//!   one per axis. Solved → the two anchor points share world position.
//! - **Distance** — 1 residual: `||p_b_world - p_a_world|| - target`.
//! - **Angle** — 1 residual: `angle(v_a_world, v_b_world) - target`
//!   (radians).
//! - **Parallel** — 2 residuals: the two components of
//!   `v_a_world × v_b_world` in the plane orthogonal to `v_a_world`.
//!   Solved → the two direction vectors are parallel.
//! - **Perpendicular** — 1 residual: `v_a_world · v_b_world`.
//! - **Tangent** — 1 residual: gap between the two axes minus the sum
//!   of radii.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// One mate variant — the structural constraint between two parts.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum MateKind {
    /// Two world-space points (one on each part, expressed in each
    /// part's local frame) must coincide.
    Coincident {
        /// Source part id.
        part_a: usize,
        /// Anchor point in `part_a`'s local frame.
        point_a: Vector3<f64>,
        /// Target part id.
        part_b: usize,
        /// Anchor point in `part_b`'s local frame.
        point_b: Vector3<f64>,
    },
    /// Two points must be at a fixed Euclidean distance.
    Distance {
        /// Source part id.
        part_a: usize,
        /// Anchor point in `part_a`'s local frame.
        point_a: Vector3<f64>,
        /// Target part id.
        part_b: usize,
        /// Anchor point in `part_b`'s local frame.
        point_b: Vector3<f64>,
        /// Target distance (world units).
        target: f64,
    },
    /// Two direction vectors must meet at a fixed angle (radians).
    Angle {
        /// Source part id.
        part_a: usize,
        /// Direction vector in `part_a`'s local frame.
        vec_a: Vector3<f64>,
        /// Target part id.
        part_b: usize,
        /// Direction vector in `part_b`'s local frame.
        vec_b: Vector3<f64>,
        /// Target angle (radians).
        target: f64,
    },
    /// Two direction vectors must be parallel.
    Parallel {
        /// Source part id.
        part_a: usize,
        /// Direction vector in `part_a`'s local frame.
        vec_a: Vector3<f64>,
        /// Target part id.
        part_b: usize,
        /// Direction vector in `part_b`'s local frame.
        vec_b: Vector3<f64>,
    },
    /// Two direction vectors must be perpendicular.
    Perpendicular {
        /// Source part id.
        part_a: usize,
        /// Direction vector in `part_a`'s local frame.
        vec_a: Vector3<f64>,
        /// Target part id.
        part_b: usize,
        /// Direction vector in `part_b`'s local frame.
        vec_b: Vector3<f64>,
    },
    /// Two cylindrical surfaces (axis + radius) must be tangent
    /// (axis-to-axis distance equals the sum of radii).
    Tangent {
        /// Source part id.
        part_a: usize,
        /// Axis origin on `part_a` in local frame.
        axis_a_origin: Vector3<f64>,
        /// Axis direction on `part_a` in local frame.
        axis_a_dir: Vector3<f64>,
        /// Radius of `part_a`'s cylinder.
        radius_a: f64,
        /// Target part id.
        part_b: usize,
        /// Axis origin on `part_b` in local frame.
        axis_b_origin: Vector3<f64>,
        /// Axis direction on `part_b` in local frame.
        axis_b_dir: Vector3<f64>,
        /// Radius of `part_b`'s cylinder.
        radius_b: f64,
    },
}

impl MateKind {
    /// Number of scalar residual equations this mate contributes to
    /// the solver's system.
    pub fn n_residuals(&self) -> usize {
        match self {
            MateKind::Coincident { .. } => 3,
            MateKind::Distance { .. } => 1,
            MateKind::Angle { .. } => 1,
            MateKind::Parallel { .. } => 2,
            MateKind::Perpendicular { .. } => 1,
            MateKind::Tangent { .. } => 1,
        }
    }

    /// Return the two part ids this mate constrains, in `(part_a,
    /// part_b)` order. Used by the solver to figure out which pose
    /// columns are touched without exhaustively pattern-matching.
    pub fn parts(&self) -> (usize, usize) {
        match self {
            MateKind::Coincident { part_a, part_b, .. }
            | MateKind::Distance { part_a, part_b, .. }
            | MateKind::Angle { part_a, part_b, .. }
            | MateKind::Parallel { part_a, part_b, .. }
            | MateKind::Perpendicular { part_a, part_b, .. }
            | MateKind::Tangent { part_a, part_b, .. } => (*part_a, *part_b),
        }
    }
}

/// One mate in the assembly — a [`MateKind`] payload plus stable id
/// and a `suppressed` toggle that lets users disable the mate without
/// deleting it (the solver skips suppressed mates).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Mate {
    /// Stable id assigned by [`crate::Assembly::add_mate`].
    pub id: usize,
    /// The constraint payload.
    pub kind: MateKind,
    /// When `true` the mate is skipped by the solver. Useful for
    /// debugging an over-constrained system or for staged solving.
    pub suppressed: bool,
}

impl Mate {
    /// Build a fresh mate, not suppressed.
    pub fn new(id: usize, kind: MateKind) -> Self {
        Self {
            id,
            kind,
            suppressed: false,
        }
    }

    /// Convenience pass-through to [`MateKind::n_residuals`].
    pub fn n_residuals(&self) -> usize {
        self.kind.n_residuals()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residual_counts_per_mate_kind() {
        let coin = MateKind::Coincident {
            part_a: 0,
            point_a: Vector3::zeros(),
            part_b: 1,
            point_b: Vector3::zeros(),
        };
        assert_eq!(coin.n_residuals(), 3);

        let dist = MateKind::Distance {
            part_a: 0,
            point_a: Vector3::zeros(),
            part_b: 1,
            point_b: Vector3::zeros(),
            target: 5.0,
        };
        assert_eq!(dist.n_residuals(), 1);

        let par = MateKind::Parallel {
            part_a: 0,
            vec_a: Vector3::x(),
            part_b: 1,
            vec_b: Vector3::x(),
        };
        assert_eq!(par.n_residuals(), 2);
    }

    #[test]
    fn parts_returns_both_part_ids() {
        let m = MateKind::Distance {
            part_a: 3,
            point_a: Vector3::zeros(),
            part_b: 7,
            point_b: Vector3::zeros(),
            target: 1.0,
        };
        assert_eq!(m.parts(), (3, 7));
    }

    #[test]
    fn new_mate_is_not_suppressed() {
        let m = Mate::new(
            0,
            MateKind::Perpendicular {
                part_a: 0,
                vec_a: Vector3::x(),
                part_b: 1,
                vec_b: Vector3::y(),
            },
        );
        assert!(!m.suppressed);
        assert_eq!(m.id, 0);
        assert_eq!(m.n_residuals(), 1);
    }
}
