//! FEM constraint kinds.

use serde::{Deserialize, Serialize};

/// One boundary condition applied to a face / plane / face-pair.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum FemConstraint {
    /// Fully fixed (encastré) constraint — all 6 DOFs zero.
    Fixed {
        /// 0-based face index in the source mesh tessellation.
        face_id: usize,
    },
    /// Prescribed displacement [ux, uy, uz] in metres.
    Displacement {
        /// 0-based face index.
        face_id: usize,
        /// Imposed displacement vector.
        value: [f64; 3],
    },
    /// Symmetry plane (constrains normal displacement + rotations
    /// about the in-plane axes).
    Symmetry {
        /// Plane normal (need not be unit length — will be normalised).
        plane_normal: [f64; 3],
        /// Point on the symmetry plane.
        plane_origin: [f64; 3],
    },
    /// Contact between two faces (small-deformation, no friction in v1).
    Contact {
        /// 0-based index of face 1.
        face1: usize,
        /// 0-based index of face 2.
        face2: usize,
    },
}

impl FemConstraint {
    /// Short label for the FEM panel's constraint list.
    pub fn kind_label(&self) -> &'static str {
        match self {
            FemConstraint::Fixed { .. } => "Fixed",
            FemConstraint::Displacement { .. } => "Displacement",
            FemConstraint::Symmetry { .. } => "Symmetry",
            FemConstraint::Contact { .. } => "Contact",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_variant_has_kind_label() {
        let cases = [
            (FemConstraint::Fixed { face_id: 0 }, "Fixed"),
            (
                FemConstraint::Displacement {
                    face_id: 0,
                    value: [0.0; 3],
                },
                "Displacement",
            ),
            (
                FemConstraint::Symmetry {
                    plane_normal: [0.0, 0.0, 1.0],
                    plane_origin: [0.0; 3],
                },
                "Symmetry",
            ),
            (FemConstraint::Contact { face1: 0, face2: 1 }, "Contact"),
        ];
        for (c, expected) in cases {
            assert_eq!(c.kind_label(), expected);
        }
    }

    #[test]
    fn ron_round_trip() {
        let c = FemConstraint::Displacement {
            face_id: 7,
            value: [0.001, 0.0, 0.0],
        };
        let j = serde_json::to_string(&c).unwrap();
        let c2: FemConstraint = serde_json::from_str(&j).unwrap();
        match c2 {
            FemConstraint::Displacement { face_id, value } => {
                assert_eq!(face_id, 7);
                assert!((value[0] - 0.001).abs() < 1e-12);
            }
            _ => panic!("wrong variant"),
        }
    }
}
