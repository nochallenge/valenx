//! Geometric measurements.
//!
//! Each [`Measurement`] variant is a self-contained recipe for
//! computing a single scalar quantity. The [`compute`] dispatcher
//! routes by variant and returns the value (or an [`InspectError`]
//! when inputs are insufficient).

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::InspectError;

/// Things this workbench can measure. All length units are mm; angles
/// are radians; area is mm²; volume is mm³.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Measurement {
    /// Euclidean distance between two 3D points.
    Distance {
        /// Origin point.
        from: Vector3<f64>,
        /// Destination point.
        to: Vector3<f64>,
    },
    /// Angle at vertex `v` formed by rays to `a` and `b`. Returns the
    /// (0..π) unsigned angle (radians).
    Angle {
        /// Vertex.
        v: Vector3<f64>,
        /// First ray endpoint.
        a: Vector3<f64>,
        /// Second ray endpoint.
        b: Vector3<f64>,
    },
    /// Radius of a fitted circle defined by `center` and a sample point.
    Radius {
        /// Circle centre.
        center: Vector3<f64>,
        /// Any point on the circle.
        sample: Vector3<f64>,
    },
    /// Total length of a polyline (sum of segment lengths). Requires
    /// at least 2 points.
    LinearLength {
        /// Polyline vertices in order.
        polyline: Vec<Vector3<f64>>,
    },
    /// Signed planar polygon area via the 2D shoelace formula (z is
    /// dropped). Requires at least 3 points.
    Area {
        /// Polygon vertices in CCW order in the xy plane.
        polygon: Vec<[f64; 2]>,
    },
    /// Closed-mesh volume via the signed-tetrahedron sum. The mesh
    /// **must** be a closed surface for the answer to be meaningful;
    /// a non-closed mesh still returns a number (just not the volume).
    Volume {
        /// Vertex positions.
        vertices: Vec<Vector3<f64>>,
        /// Triangle indices (3 per triangle).
        triangles: Vec<[usize; 3]>,
    },
    /// Length of the axis-aligned bounding-box diagonal of a point set.
    /// Useful as a single "size" scalar for a part. Requires at least
    /// 1 point.
    BoundingBoxDiagonal {
        /// Points to bound.
        points: Vec<Vector3<f64>>,
    },
}

impl Measurement {
    /// Short label used in error messages and CSV headers.
    pub fn label(&self) -> &'static str {
        match self {
            Self::Distance { .. } => "Distance",
            Self::Angle { .. } => "Angle",
            Self::Radius { .. } => "Radius",
            Self::LinearLength { .. } => "LinearLength",
            Self::Area { .. } => "Area",
            Self::Volume { .. } => "Volume",
            Self::BoundingBoxDiagonal { .. } => "BoundingBoxDiagonal",
        }
    }
}

/// Compute the scalar value for one measurement.
pub fn compute(m: &Measurement) -> Result<f64, InspectError> {
    match m {
        Measurement::Distance { from, to } => Ok((to - from).norm()),
        Measurement::Angle { v, a, b } => {
            let va = a - v;
            let vb = b - v;
            let na = va.norm();
            let nb = vb.norm();
            if na < f64::EPSILON || nb < f64::EPSILON {
                return Err(InspectError::NotEnoughGeometry {
                    kind: "Angle",
                    reason: "ray length is zero".into(),
                });
            }
            let cos = (va.dot(&vb) / (na * nb)).clamp(-1.0, 1.0);
            Ok(cos.acos())
        }
        Measurement::Radius { center, sample } => Ok((sample - center).norm()),
        Measurement::LinearLength { polyline } => {
            if polyline.len() < 2 {
                return Err(InspectError::NotEnoughGeometry {
                    kind: "LinearLength",
                    reason: "needs at least 2 points".into(),
                });
            }
            let mut s = 0.0;
            for w in polyline.windows(2) {
                s += (w[1] - w[0]).norm();
            }
            Ok(s)
        }
        Measurement::Area { polygon } => {
            if polygon.len() < 3 {
                return Err(InspectError::NotEnoughGeometry {
                    kind: "Area",
                    reason: "needs at least 3 points".into(),
                });
            }
            let mut s = 0.0;
            let n = polygon.len();
            for i in 0..n {
                let j = (i + 1) % n;
                s += polygon[i][0] * polygon[j][1] - polygon[j][0] * polygon[i][1];
            }
            Ok(s.abs() * 0.5)
        }
        Measurement::Volume {
            vertices,
            triangles,
        } => {
            if triangles.is_empty() {
                return Err(InspectError::NotEnoughGeometry {
                    kind: "Volume",
                    reason: "needs at least one triangle".into(),
                });
            }
            let mut v6 = 0.0;
            for tri in triangles {
                let a = vertices.get(tri[0]).ok_or(InspectError::BadParameter {
                    name: "triangle.idx0",
                    reason: format!("out of bounds: {}", tri[0]),
                })?;
                let b = vertices.get(tri[1]).ok_or(InspectError::BadParameter {
                    name: "triangle.idx1",
                    reason: format!("out of bounds: {}", tri[1]),
                })?;
                let c = vertices.get(tri[2]).ok_or(InspectError::BadParameter {
                    name: "triangle.idx2",
                    reason: format!("out of bounds: {}", tri[2]),
                })?;
                // signed tetra volume / 6
                v6 += a.dot(&b.cross(c));
            }
            Ok((v6 / 6.0).abs())
        }
        Measurement::BoundingBoxDiagonal { points } => {
            if points.is_empty() {
                return Err(InspectError::NotEnoughGeometry {
                    kind: "BoundingBoxDiagonal",
                    reason: "needs at least 1 point".into(),
                });
            }
            let mut min = points[0];
            let mut max = points[0];
            for p in points.iter().skip(1) {
                for i in 0..3 {
                    if p[i] < min[i] {
                        min[i] = p[i];
                    }
                    if p[i] > max[i] {
                        max[i] = p[i];
                    }
                }
            }
            Ok((max - min).norm())
        }
    }
}

/// Triangle-soup form of a mesh suitable for
/// [`Measurement::Volume`] inputs.
pub type TriangleSoup = (Vec<Vector3<f64>>, Vec<[usize; 3]>);

/// Tessellate a [`valenx_mesh::Mesh`]'s triangle-only blocks into the
/// (vertices, triangles) shape [`Measurement::Volume`] expects.
///
/// Returns `Err(InspectError::NotEnoughGeometry)` when the mesh has no
/// triangle blocks. Non-triangle blocks (tets, hexes, …) are skipped
/// silently — those don't carry a meaningful surface for the
/// signed-tetra method without an extra surface-extraction step.
pub fn triangles_from_mesh(mesh: &valenx_mesh::Mesh) -> Result<TriangleSoup, InspectError> {
    use valenx_mesh::element::ElementType;
    let mut tris: Vec<[usize; 3]> = Vec::new();
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        // 3 connectivity entries per triangle, stored flat.
        for chunk in block.connectivity.chunks_exact(3) {
            tris.push([chunk[0] as usize, chunk[1] as usize, chunk[2] as usize]);
        }
    }
    if tris.is_empty() {
        return Err(InspectError::NotEnoughGeometry {
            kind: "Volume",
            reason: "mesh has no Tri3 blocks".into(),
        });
    }
    Ok((mesh.nodes.clone(), tris))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distance_unit_x() {
        let m = Measurement::Distance {
            from: Vector3::zeros(),
            to: Vector3::new(3.0, 4.0, 0.0),
        };
        assert!((compute(&m).unwrap() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn right_angle_is_pi_over_2() {
        let m = Measurement::Angle {
            v: Vector3::zeros(),
            a: Vector3::new(1.0, 0.0, 0.0),
            b: Vector3::new(0.0, 1.0, 0.0),
        };
        let got = compute(&m).unwrap();
        assert!((got - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    }

    #[test]
    fn polyline_length_three_segments() {
        let m = Measurement::LinearLength {
            polyline: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(1.0, 1.0, 0.0),
                Vector3::new(0.0, 1.0, 0.0),
            ],
        };
        assert!((compute(&m).unwrap() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn unit_square_area() {
        let m = Measurement::Area {
            polygon: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        };
        assert!((compute(&m).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn unit_tetra_volume() {
        // Tetrahedron with 4 triangles
        let vs = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let tris = vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        let m = Measurement::Volume {
            vertices: vs,
            triangles: tris,
        };
        let v = compute(&m).unwrap();
        assert!((v - 1.0 / 6.0).abs() < 1e-12);
    }

    #[test]
    fn bbox_diag_unit_cube() {
        let m = Measurement::BoundingBoxDiagonal {
            points: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(1.0, 1.0, 1.0),
            ],
        };
        assert!((compute(&m).unwrap() - 3f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn empty_polyline_errors() {
        let m = Measurement::LinearLength {
            polyline: vec![Vector3::zeros()],
        };
        assert!(compute(&m).is_err());
    }
}
