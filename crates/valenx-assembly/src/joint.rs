//! Mechanical joints — constrained motion between two parts.
//!
//! A [`Joint`] differs from a [`crate::Mate`] in that it carries a
//! *parameter* (current state — rotation angle for Revolute,
//! translation distance for Prismatic, etc.) plus an applier that
//! turns that parameter into a relative transform between the two
//! parts.
//!
//! The kinematics live in [`crate::kinematics`]. Joints are *not*
//! consumed by the constraint solver — instead, the UI's "joint
//! slider" drives the parameter and re-applies the kinematics each
//! frame, producing a preview of the constrained motion.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// One joint variant.
///
/// Every joint has a `part_a` (the "ground" side) and a `part_b` (the
/// "driven" side). When the joint's parameter changes, `part_b`'s
/// transform is recomputed relative to `part_a`'s. Setting `part_a`'s
/// `fixed` flag on the [`crate::Part`] turns this into a one-way
/// driven motion; leaving both un-fixed lets the user drag either
/// side and watch the other follow.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum JointKind {
    /// Rigid bond — part_b's transform is locked to part_a's. The
    /// parameter is unused.
    Fixed {
        /// Source part id.
        part_a: usize,
        /// Constrained part id.
        part_b: usize,
    },
    /// Hinge — part_b rotates about a fixed axis on part_a. The
    /// parameter is the rotation angle in radians.
    Revolute {
        /// Source part id.
        part_a: usize,
        /// Constrained part id.
        part_b: usize,
        /// Axis origin in part_a's local frame.
        axis_origin: Vector3<f64>,
        /// Axis direction in part_a's local frame (normalized inside
        /// [`crate::kinematics`]).
        axis_dir: Vector3<f64>,
    },
    /// Slider — part_b translates along a fixed axis on part_a. The
    /// parameter is the translation distance in world units.
    Prismatic {
        /// Source part id.
        part_a: usize,
        /// Constrained part id.
        part_b: usize,
        /// Slide direction in part_a's local frame.
        axis_dir: Vector3<f64>,
    },
    /// Slider + hinge sharing one axis — part_b can both rotate around
    /// and slide along the axis. The parameter is the rotation angle
    /// (the slide distance is *not* parameterized in v1 — Phase 6.5
    /// will add the dual-parameter form).
    Cylindrical {
        /// Source part id.
        part_a: usize,
        /// Constrained part id.
        part_b: usize,
        /// Axis origin in part_a's local frame.
        axis_origin: Vector3<f64>,
        /// Axis direction in part_a's local frame.
        axis_dir: Vector3<f64>,
    },
    /// Ball joint — part_b pivots freely about a single point. The
    /// parameter is unused for v1 (the rotation is read from
    /// part_b's transform).
    Spherical {
        /// Source part id.
        part_a: usize,
        /// Constrained part id.
        part_b: usize,
        /// Pivot point in part_a's local frame.
        point: Vector3<f64>,
    },
    /// Planar / "face" joint — part_b can translate within a plane on
    /// part_a but cannot leave the plane. The parameter is unused for
    /// v1.
    Planar {
        /// Source part id.
        part_a: usize,
        /// Constrained part id.
        part_b: usize,
        /// Plane origin in part_a's local frame.
        plane_origin: Vector3<f64>,
        /// Plane normal in part_a's local frame.
        plane_normal: Vector3<f64>,
    },
}

impl JointKind {
    /// Return the two part ids this joint connects, `(part_a, part_b)`.
    pub fn parts(&self) -> (usize, usize) {
        match self {
            JointKind::Fixed { part_a, part_b, .. }
            | JointKind::Revolute { part_a, part_b, .. }
            | JointKind::Prismatic { part_a, part_b, .. }
            | JointKind::Cylindrical { part_a, part_b, .. }
            | JointKind::Spherical { part_a, part_b, .. }
            | JointKind::Planar { part_a, part_b, .. } => (*part_a, *part_b),
        }
    }

    /// Short identifier used by the UI dropdowns and the toolbox
    /// status line.
    pub fn label(&self) -> &'static str {
        match self {
            JointKind::Fixed { .. } => "Fixed",
            JointKind::Revolute { .. } => "Revolute",
            JointKind::Prismatic { .. } => "Prismatic",
            JointKind::Cylindrical { .. } => "Cylindrical",
            JointKind::Spherical { .. } => "Spherical",
            JointKind::Planar { .. } => "Planar",
        }
    }
}

/// One joint in the assembly — kind, parameter (current state), and
/// suppress flag.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Joint {
    /// Stable id assigned by [`crate::Assembly::add_joint`].
    pub id: usize,
    /// The constrained-motion payload.
    pub kind: JointKind,
    /// Current parameter — meaning depends on `kind` (angle for
    /// Revolute / Cylindrical, distance for Prismatic, unused for
    /// Fixed / Spherical / Planar).
    pub parameter: f64,
    /// When `true` the kinematics applier skips this joint (the
    /// preview pretends the joint isn't there).
    pub suppressed: bool,
}

impl Joint {
    /// Build a fresh joint with parameter = 0 and not suppressed.
    pub fn new(id: usize, kind: JointKind) -> Self {
        Self {
            id,
            kind,
            parameter: 0.0,
            suppressed: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parts_returns_pair() {
        let j = JointKind::Revolute {
            part_a: 4,
            part_b: 9,
            axis_origin: Vector3::zeros(),
            axis_dir: Vector3::z(),
        };
        assert_eq!(j.parts(), (4, 9));
    }

    #[test]
    fn labels_are_distinct() {
        let kinds = [
            JointKind::Fixed {
                part_a: 0,
                part_b: 1,
            },
            JointKind::Revolute {
                part_a: 0,
                part_b: 1,
                axis_origin: Vector3::zeros(),
                axis_dir: Vector3::z(),
            },
            JointKind::Prismatic {
                part_a: 0,
                part_b: 1,
                axis_dir: Vector3::x(),
            },
            JointKind::Cylindrical {
                part_a: 0,
                part_b: 1,
                axis_origin: Vector3::zeros(),
                axis_dir: Vector3::z(),
            },
            JointKind::Spherical {
                part_a: 0,
                part_b: 1,
                point: Vector3::zeros(),
            },
            JointKind::Planar {
                part_a: 0,
                part_b: 1,
                plane_origin: Vector3::zeros(),
                plane_normal: Vector3::z(),
            },
        ];
        let labels: Vec<&str> = kinds.iter().map(|k| k.label()).collect();
        let unique: std::collections::HashSet<_> = labels.iter().collect();
        assert_eq!(unique.len(), kinds.len(), "labels collided: {labels:?}");
    }

    #[test]
    fn new_joint_defaults_zero_param_unsuppressed() {
        let j = Joint::new(
            0,
            JointKind::Fixed {
                part_a: 0,
                part_b: 1,
            },
        );
        assert_eq!(j.parameter, 0.0);
        assert!(!j.suppressed);
    }
}
