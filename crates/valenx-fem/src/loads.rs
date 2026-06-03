//! FEM load kinds.

use serde::{Deserialize, Serialize};

/// One applied load.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum FemLoad {
    /// Concentrated force at a single vertex.
    PointForce {
        /// 0-based vertex index in the source mesh tessellation.
        vertex: usize,
        /// Force vector [Fx, Fy, Fz] in newtons.
        vector: [f64; 3],
    },
    /// Pressure normal to a face.
    PressureForce {
        /// 0-based face index.
        face: usize,
        /// Pressure magnitude in Pa (positive = pushes into the face).
        magnitude: f64,
    },
    /// Force-per-unit-length along an edge.
    LineForce {
        /// 0-based edge index in the source mesh.
        edge: usize,
        /// Force vector per metre.
        vector: [f64; 3],
    },
    /// Body force per unit volume (e.g. centrifugal).
    BodyForce {
        /// Vector in N/m^3.
        vector: [f64; 3],
    },
    /// Earth gravity (-9.81 m/s^2 along the load's direction).
    Gravity {
        /// Direction of gravity (need not be unit-length).
        direction: [f64; 3],
    },
    /// Imposed temperature on a face.
    Temperature {
        /// 0-based face index.
        face: usize,
        /// Temperature in K.
        value: f64,
    },
}

impl FemLoad {
    /// Short label for the FEM panel's load list.
    pub fn kind_label(&self) -> &'static str {
        match self {
            FemLoad::PointForce { .. } => "Point Force",
            FemLoad::PressureForce { .. } => "Pressure",
            FemLoad::LineForce { .. } => "Line Force",
            FemLoad::BodyForce { .. } => "Body Force",
            FemLoad::Gravity { .. } => "Gravity",
            FemLoad::Temperature { .. } => "Temperature",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_variant_has_kind_label() {
        let cases = [
            (
                FemLoad::PointForce {
                    vertex: 0,
                    vector: [0.0; 3],
                },
                "Point Force",
            ),
            (
                FemLoad::PressureForce {
                    face: 0,
                    magnitude: 1.0,
                },
                "Pressure",
            ),
            (
                FemLoad::LineForce {
                    edge: 0,
                    vector: [0.0; 3],
                },
                "Line Force",
            ),
            (FemLoad::BodyForce { vector: [0.0; 3] }, "Body Force"),
            (
                FemLoad::Gravity {
                    direction: [0.0, 0.0, -1.0],
                },
                "Gravity",
            ),
            (
                FemLoad::Temperature {
                    face: 0,
                    value: 293.15,
                },
                "Temperature",
            ),
        ];
        for (l, expected) in cases {
            assert_eq!(l.kind_label(), expected);
        }
    }
}
