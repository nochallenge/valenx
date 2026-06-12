//! Optimal structural superposition and RMSD.
//!
//! Given two equal-length, one-to-one corresponded point sets this
//! module finds the rigid transform (rotation + translation) that
//! minimises the root-mean-square deviation between them.
//!
//! Two solvers are provided, both standard:
//!
//! - [`kabsch`] — the Kabsch (1976) algorithm: covariance matrix →
//!   SVD → rotation, with the reflection correction.
//! - [`quaternion_superpose`] — the Kearsley / quaternion variant:
//!   build the 4×4 key matrix, take its largest-eigenvalue
//!   eigenvector, convert to a rotation.
//!
//! Both yield the same optimal rotation; the quaternion route is
//! offered because it never needs the reflection fix-up and is
//! numerically convenient.

use crate::error::{BiostructError, Result};
use nalgebra::{Matrix3, Matrix4, Point3, Vector3};

/// A rigid-body transform: rotate about the origin, then translate.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidTransform {
    /// 3×3 rotation matrix (proper, `det = +1`).
    pub rotation: Matrix3<f64>,
    /// Translation vector applied after rotation.
    pub translation: Vector3<f64>,
}

impl RigidTransform {
    /// The identity transform.
    pub fn identity() -> Self {
        RigidTransform {
            rotation: Matrix3::identity(),
            translation: Vector3::zeros(),
        }
    }

    /// Apply the transform to a point: `x' = R·x + t`.
    pub fn apply(&self, p: &Point3<f64>) -> Point3<f64> {
        Point3::from(self.rotation * p.coords + self.translation)
    }

    /// Apply the transform to every point of a slice.
    pub fn apply_all(&self, points: &[Point3<f64>]) -> Vec<Point3<f64>> {
        points.iter().map(|p| self.apply(p)).collect()
    }
}

/// The outcome of a superposition: the transform that maps the
/// *mobile* set onto the *reference* set, and the resulting RMSD.
#[derive(Clone, Debug, PartialEq)]
pub struct Superposition {
    /// The transform mapping mobile → reference.
    pub transform: RigidTransform,
    /// Root-mean-square deviation after superposition, ångström.
    pub rmsd: f64,
    /// Number of corresponded point pairs.
    pub n: usize,
}

/// Plain RMSD of two equal-length point sets *without* any
/// superposition — measures them where they sit.
pub fn rmsd(a: &[Point3<f64>], b: &[Point3<f64>]) -> Result<f64> {
    check_pair(a, b)?;
    let mut sum = 0.0;
    for (p, q) in a.iter().zip(b) {
        sum += (p - q).norm_squared();
    }
    Ok((sum / a.len() as f64).sqrt())
}

/// Per-point deviation of two equal-length sets (no superposition).
pub fn per_point_deviation(a: &[Point3<f64>], b: &[Point3<f64>]) -> Result<Vec<f64>> {
    check_pair(a, b)?;
    Ok(a.iter().zip(b).map(|(p, q)| (p - q).norm()).collect())
}

/// Kabsch optimal superposition of `mobile` onto `reference`.
///
/// Steps: centre both sets on their centroids, form the
/// cross-covariance `H = Σ mobileᵢ · referenceᵢᵀ`, take its SVD
/// `H = U S Vᵀ`, and set `R = V · D · Uᵀ` where `D` corrects a
/// reflection when `det(V·Uᵀ) < 0`.
pub fn kabsch(mobile: &[Point3<f64>], reference: &[Point3<f64>]) -> Result<Superposition> {
    check_pair(mobile, reference)?;
    let n = mobile.len();

    let cm = centroid(mobile);
    let cr = centroid(reference);

    // Covariance H = Σ (mobile - cm)(reference - cr)^T.
    let mut h = Matrix3::zeros();
    for (p, q) in mobile.iter().zip(reference) {
        let pm = p.coords - cm;
        let qr = q.coords - cr;
        h += pm * qr.transpose();
    }

    let svd = h.svd(true, true);
    let u = svd
        .u
        .ok_or_else(|| BiostructError::invalid("svd", "covariance SVD produced no U"))?;
    let v_t = svd
        .v_t
        .ok_or_else(|| BiostructError::invalid("svd", "covariance SVD produced no Vᵀ"))?;
    let v = v_t.transpose();

    // Reflection correction: ensure a proper rotation.
    let d = (v * u.transpose()).determinant();
    let mut diag = Matrix3::identity();
    if d < 0.0 {
        diag[(2, 2)] = -1.0;
    }
    let rotation = v * diag * u.transpose();

    // Translation: t = cr - R·cm.
    let translation = cr - rotation * cm;
    let transform = RigidTransform {
        rotation,
        translation,
    };

    let moved = transform.apply_all(mobile);
    let result_rmsd = rmsd(&moved, reference)?;
    Ok(Superposition {
        transform,
        rmsd: result_rmsd,
        n,
    })
}

/// Quaternion (Kearsley) optimal superposition.
///
/// Builds the symmetric 4×4 Kearsley matrix from the centred
/// coordinate differences and sums; its eigenvector for the *largest*
/// eigenvalue is the optimal rotation quaternion. The resulting
/// rotation is identical to [`kabsch`]'s.
pub fn quaternion_superpose(
    mobile: &[Point3<f64>],
    reference: &[Point3<f64>],
) -> Result<Superposition> {
    check_pair(mobile, reference)?;
    let n = mobile.len();
    let cm = centroid(mobile);
    let cr = centroid(reference);

    // Cross-covariance R.
    let mut r = Matrix3::zeros();
    for (p, q) in mobile.iter().zip(reference) {
        let pm = p.coords - cm;
        let qr = q.coords - cr;
        r += pm * qr.transpose();
    }

    // Kearsley 4x4 key matrix from the covariance components.
    let (r11, r12, r13) = (r[(0, 0)], r[(0, 1)], r[(0, 2)]);
    let (r21, r22, r23) = (r[(1, 0)], r[(1, 1)], r[(1, 2)]);
    let (r31, r32, r33) = (r[(2, 0)], r[(2, 1)], r[(2, 2)]);

    let mut k = Matrix4::zeros();
    k[(0, 0)] = r11 + r22 + r33;
    k[(0, 1)] = r23 - r32;
    k[(0, 2)] = r31 - r13;
    k[(0, 3)] = r12 - r21;
    k[(1, 0)] = k[(0, 1)];
    k[(1, 1)] = r11 - r22 - r33;
    k[(1, 2)] = r12 + r21;
    k[(1, 3)] = r13 + r31;
    k[(2, 0)] = k[(0, 2)];
    k[(2, 1)] = k[(1, 2)];
    k[(2, 2)] = -r11 + r22 - r33;
    k[(2, 3)] = r23 + r32;
    k[(3, 0)] = k[(0, 3)];
    k[(3, 1)] = k[(1, 3)];
    k[(3, 2)] = k[(2, 3)];
    k[(3, 3)] = -r11 - r22 + r33;

    // The optimal quaternion is the eigenvector of the LARGEST
    // eigenvalue of this matrix.
    let eig = nalgebra::SymmetricEigen::new(k);
    let mut best = 0usize;
    for i in 1..4 {
        if eig.eigenvalues[i] > eig.eigenvalues[best] {
            best = i;
        }
    }
    let q = eig.eigenvectors.column(best);
    let (qw, qx, qy, qz) = (q[0], q[1], q[2], q[3]);

    // Quaternion -> rotation matrix.
    let rotation = quat_to_matrix(qw, qx, qy, qz);
    let translation = cr - rotation * cm;
    let transform = RigidTransform {
        rotation,
        translation,
    };
    let moved = transform.apply_all(mobile);
    let result_rmsd = rmsd(&moved, reference)?;
    Ok(Superposition {
        transform,
        rmsd: result_rmsd,
        n,
    })
}

/// Convert a unit quaternion `(w, x, y, z)` to a 3×3 rotation matrix.
fn quat_to_matrix(w: f64, x: f64, y: f64, z: f64) -> Matrix3<f64> {
    let norm = (w * w + x * x + y * y + z * z).sqrt().max(1e-12);
    let (w, x, y, z) = (w / norm, x / norm, y / norm, z / norm);
    Matrix3::new(
        1.0 - 2.0 * (y * y + z * z),
        2.0 * (x * y - w * z),
        2.0 * (x * z + w * y),
        2.0 * (x * y + w * z),
        1.0 - 2.0 * (x * x + z * z),
        2.0 * (y * z - w * x),
        2.0 * (x * z - w * y),
        2.0 * (y * z + w * x),
        1.0 - 2.0 * (x * x + y * y),
    )
}

/// Per-residue RMSD after a superposition: applies `transform` to the
/// mobile points and returns the deviation of each pair.
pub fn per_residue_rmsd(
    transform: &RigidTransform,
    mobile: &[Point3<f64>],
    reference: &[Point3<f64>],
) -> Result<Vec<f64>> {
    check_pair(mobile, reference)?;
    let moved = transform.apply_all(mobile);
    per_point_deviation(&moved, reference)
}

/// Centroid of a point set (caller guarantees non-empty).
fn centroid(points: &[Point3<f64>]) -> Vector3<f64> {
    let mut acc = Vector3::zeros();
    for p in points {
        acc += p.coords;
    }
    acc / points.len() as f64
}

/// Validate that two point sets are non-empty and equal length.
fn check_pair(a: &[Point3<f64>], b: &[Point3<f64>]) -> Result<()> {
    if a.is_empty() || b.is_empty() {
        return Err(BiostructError::invalid("points", "empty point set"));
    }
    if a.len() != b.len() {
        return Err(BiostructError::invalid(
            "points",
            format!(
                "point sets must be equal length: {} vs {}",
                a.len(),
                b.len()
            ),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cloud() -> Vec<Point3<f64>> {
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
            Point3::new(0.0, 0.0, 3.0),
            Point3::new(1.5, 1.5, 1.5),
        ]
    }

    #[test]
    fn rmsd_of_identical_sets_is_zero() {
        let c = cloud();
        assert!(rmsd(&c, &c).unwrap() < 1e-12);
    }

    #[test]
    fn kabsch_recovers_a_known_rotation() {
        // Rotate the cloud 90 deg about z, then translate.
        let c = cloud();
        let rot = Matrix3::new(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let t = Vector3::new(5.0, -3.0, 2.0);
        let moved: Vec<_> = c.iter().map(|p| Point3::from(rot * p.coords + t)).collect();
        // kabsch(mobile=c, reference=moved) should find (rot, t).
        let sup = kabsch(&c, &moved).unwrap();
        assert!(sup.rmsd < 1e-9, "rmsd after fit was {}", sup.rmsd);
        assert!((sup.transform.rotation - rot).norm() < 1e-9);
        assert!((sup.transform.translation - t).norm() < 1e-9);
    }

    #[test]
    fn quaternion_matches_kabsch() {
        let c = cloud();
        // some arbitrary rotation about (1,1,1)
        let axis = nalgebra::Unit::new_normalize(Vector3::new(1.0, 1.0, 1.0));
        let angle = 0.7_f64;
        let rot: Matrix3<f64> = nalgebra::Rotation3::from_axis_angle(&axis, angle).into_inner();
        let moved: Vec<_> = c.iter().map(|p| Point3::from(rot * p.coords)).collect();
        let k = kabsch(&c, &moved).unwrap();
        let q = quaternion_superpose(&c, &moved).unwrap();
        assert!((k.transform.rotation - q.transform.rotation).norm() < 1e-7);
        assert!(q.rmsd < 1e-7);
    }

    #[test]
    fn superposition_reduces_rmsd() {
        // Two clouds that differ by a rotation: the post-fit RMSD must
        // be far below the as-placed RMSD.
        let c = cloud();
        let rot = Matrix3::new(0.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0);
        let moved: Vec<_> = c.iter().map(|p| Point3::from(rot * p.coords)).collect();
        let raw = rmsd(&c, &moved).unwrap();
        let fitted = kabsch(&c, &moved).unwrap().rmsd;
        assert!(fitted < raw);
        assert!(fitted < 1e-9);
    }

    #[test]
    fn reflection_is_corrected() {
        // A mirrored cloud must still produce a proper rotation
        // (det = +1), not an improper one.
        let c = cloud();
        let mirror: Vec<_> = c.iter().map(|p| Point3::new(-p.x, p.y, p.z)).collect();
        let sup = kabsch(&c, &mirror).unwrap();
        let det = sup.transform.rotation.determinant();
        assert!((det - 1.0).abs() < 1e-6, "rotation det was {det}");
    }

    #[test]
    fn per_residue_rmsd_lengths() {
        let c = cloud();
        let sup = kabsch(&c, &c).unwrap();
        let dev = per_residue_rmsd(&sup.transform, &c, &c).unwrap();
        assert_eq!(dev.len(), c.len());
        assert!(dev.iter().all(|d| *d < 1e-9));
    }

    #[test]
    fn mismatched_lengths_error() {
        let a = cloud();
        let b = vec![Point3::origin()];
        assert!(rmsd(&a, &b).is_err());
        assert!(kabsch(&a, &b).is_err());
    }
}
