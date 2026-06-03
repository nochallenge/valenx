//! PCA-based shape descriptor.

use nalgebra::{DMatrix, Vector3};

use crate::error::LibiglError;
use crate::triangle::TriMesh;

/// Returns the 3 eigenvalues of the vertex covariance matrix in
/// descending order — a coarse shape descriptor (sphericity /
/// flatness / linearity ratios fall out of these).
pub fn shape_descriptor(mesh: &TriMesh) -> Result<Vec<f64>, LibiglError> {
    if mesh.vertices.len() < 3 {
        return Err(LibiglError::NotEnough {
            what: "vertices",
            needed: 3,
            given: mesh.vertices.len(),
        });
    }
    let centroid = mesh.centroid();
    // Build 3x3 covariance.
    let mut c = [[0.0_f64; 3]; 3];
    for v in &mesh.vertices {
        let d = v - centroid;
        for i in 0..3 {
            for j in 0..3 {
                c[i][j] += d[i] * d[j];
            }
        }
    }
    let n = mesh.vertices.len() as f64;
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] /= n;
        }
    }
    let m = DMatrix::<f64>::from_row_slice(
        3,
        3,
        &[
            c[0][0], c[0][1], c[0][2], c[1][0], c[1][1], c[1][2], c[2][0], c[2][1], c[2][2],
        ],
    );
    let eig = m.symmetric_eigen();
    let mut vals: Vec<f64> = eig.eigenvalues.iter().copied().collect();
    vals.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let _ = Vector3::<f64>::zeros();
    Ok(vals)
}
