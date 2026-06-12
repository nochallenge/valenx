//! RMSD & RMSF with optimal superposition — **roadmap feature 29**.
//!
//! Comparing two structures (or a trajectory against a reference)
//! requires first removing the rigid-body translation and rotation —
//! otherwise a molecule that merely drifted would look "different".
//!
//! - **Kabsch superposition** ([`kabsch`]) finds the rotation that
//!   minimises the RMSD between two equal-length point sets. Both sets
//!   are first centred on their centroids; the optimal rotation comes
//!   from the singular-value decomposition of the cross-covariance
//!   matrix, with a reflection guard so the result is a proper
//!   rotation (det = +1).
//!
//! - **RMSD** ([`rmsd`]) is the root-mean-square atom deviation
//!   *after* the optimal superposition — the standard measure of how
//!   far one structure is from another.
//!
//! - **RMSF** ([`rmsf`]) is the per-atom root-mean-square *fluctuation*
//!   over a trajectory: each frame is Kabsch-fitted to a reference,
//!   then the spread of each atom about its mean position is reported.
//!   It pinpoints the flexible regions of a molecule (loops vs the
//!   rigid core) and maps directly onto crystallographic B-factors.
//!
//! The SVD is `nalgebra`'s; everything else is hand-rolled.

use nalgebra::{Matrix3, Vector3};

use crate::error::{MdError, Result};

/// Centroid (mean position) of a point set.
fn centroid(points: &[Vector3<f64>]) -> Vector3<f64> {
    if points.is_empty() {
        return Vector3::zeros();
    }
    points.iter().sum::<Vector3<f64>>() / points.len() as f64
}

/// The result of an optimal Kabsch superposition of `mobile` onto
/// `reference`.
#[derive(Clone, Debug, PartialEq)]
pub struct Superposition {
    /// The proper rotation matrix (det = +1).
    pub rotation: Matrix3<f64>,
    /// Translation: `reference_centroid − rotation·mobile_centroid`.
    pub translation: Vector3<f64>,
    /// RMSD after the optimal superposition (nm).
    pub rmsd: f64,
}

impl Superposition {
    /// Applies this superposition to a point: `rotation·p + translation`.
    pub fn apply(&self, p: Vector3<f64>) -> Vector3<f64> {
        self.rotation * p + self.translation
    }

    /// Applies this superposition to a whole point set.
    pub fn apply_all(&self, points: &[Vector3<f64>]) -> Vec<Vector3<f64>> {
        points.iter().map(|p| self.apply(*p)).collect()
    }
}

/// Computes the optimal rigid-body superposition of `mobile` onto
/// `reference` (the **Kabsch** algorithm).
///
/// # Errors
/// [`MdError::DimensionMismatch`] if the two sets differ in length;
/// [`MdError::Invalid`] if either is empty.
pub fn kabsch(mobile: &[Vector3<f64>], reference: &[Vector3<f64>]) -> Result<Superposition> {
    if mobile.len() != reference.len() {
        return Err(MdError::dimension(format!(
            "{} mobile atoms vs {} reference atoms",
            mobile.len(),
            reference.len()
        )));
    }
    if mobile.is_empty() {
        return Err(MdError::invalid("kabsch", "needs at least one atom"));
    }
    let n = mobile.len();
    let c_mob = centroid(mobile);
    let c_ref = centroid(reference);

    // Cross-covariance matrix H = Σ (mobile−c)·(ref−c)ᵀ.
    let mut h = Matrix3::zeros();
    for i in 0..n {
        let p = mobile[i] - c_mob;
        let q = reference[i] - c_ref;
        h += p * q.transpose();
    }

    // Optimal rotation from the SVD of H.
    let svd = h.svd(true, true);
    let u = svd
        .u
        .ok_or_else(|| MdError::invalid("kabsch", "SVD failed"))?;
    let v_t = svd
        .v_t
        .ok_or_else(|| MdError::invalid("kabsch", "SVD failed"))?;
    // R = V·Uᵀ; fix a reflection by flipping the sign of the last
    // column of V when det(V·Uᵀ) is negative.
    let mut rotation = v_t.transpose() * u.transpose();
    if rotation.determinant() < 0.0 {
        let mut v = v_t.transpose();
        v.column_mut(2).neg_mut();
        rotation = v * u.transpose();
    }
    let translation = c_ref - rotation * c_mob;

    // RMSD after superposition.
    let mut sumsq = 0.0;
    for i in 0..n {
        let fitted = rotation * (mobile[i] - c_mob) + c_ref;
        sumsq += (fitted - reference[i]).norm_squared();
    }
    let rmsd = (sumsq / n as f64).sqrt();

    Ok(Superposition {
        rotation,
        translation,
        rmsd,
    })
}

/// The minimal RMSD between two structures after optimal Kabsch
/// superposition (nm).
///
/// # Errors
/// As [`kabsch`].
pub fn rmsd(mobile: &[Vector3<f64>], reference: &[Vector3<f64>]) -> Result<f64> {
    Ok(kabsch(mobile, reference)?.rmsd)
}

/// The RMSD between two structures *without* any superposition — a
/// plain root-mean-square coordinate difference.
///
/// # Errors
/// [`MdError::DimensionMismatch`] for unequal lengths;
/// [`MdError::Invalid`] for empty input.
pub fn rmsd_no_fit(a: &[Vector3<f64>], b: &[Vector3<f64>]) -> Result<f64> {
    if a.len() != b.len() {
        return Err(MdError::dimension("structures differ in atom count"));
    }
    if a.is_empty() {
        return Err(MdError::invalid("rmsd", "needs at least one atom"));
    }
    let sumsq: f64 = a.iter().zip(b).map(|(x, y)| (x - y).norm_squared()).sum();
    Ok((sumsq / a.len() as f64).sqrt())
}

/// Per-atom root-mean-square fluctuation over a trajectory.
///
/// Every frame is Kabsch-fitted onto `reference`; then for each atom
/// the RMS spread about its mean fitted position is returned. The
/// output vector is one value per atom (nm).
///
/// # Errors
/// [`MdError::Invalid`] for an empty trajectory; [`MdError::DimensionMismatch`]
/// if a frame or the reference disagrees on atom count.
pub fn rmsf(frames: &[Vec<Vector3<f64>>], reference: &[Vector3<f64>]) -> Result<Vec<f64>> {
    if frames.is_empty() {
        return Err(MdError::invalid("rmsf", "needs at least one frame"));
    }
    let natoms = reference.len();
    if natoms == 0 {
        return Err(MdError::invalid("rmsf", "reference has no atoms"));
    }
    if frames.iter().any(|f| f.len() != natoms) {
        return Err(MdError::dimension(
            "a frame disagrees with the reference on atom count",
        ));
    }
    // Fit every frame onto the reference.
    let fitted: Vec<Vec<Vector3<f64>>> = frames
        .iter()
        .map(|f| {
            let sup = kabsch(f, reference)?;
            Ok(sup.apply_all(f))
        })
        .collect::<Result<_>>()?;

    // Per-atom mean position over the fitted frames.
    let n_frames = fitted.len() as f64;
    let mut mean = vec![Vector3::zeros(); natoms];
    for frame in &fitted {
        for a in 0..natoms {
            mean[a] += frame[a];
        }
    }
    for m in &mut mean {
        *m /= n_frames;
    }
    // Per-atom RMS fluctuation about the mean.
    let mut fluct = vec![0.0; natoms];
    for frame in &fitted {
        for a in 0..natoms {
            fluct[a] += (frame[a] - mean[a]).norm_squared();
        }
    }
    for f in &mut fluct {
        *f = (*f / n_frames).sqrt();
    }
    Ok(fluct)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A small L-shaped point set.
    fn shape() -> Vec<Vector3<f64>> {
        vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 2.0, 0.0),
        ]
    }

    #[test]
    fn identical_structures_have_zero_rmsd() {
        let s = shape();
        assert!(rmsd(&s, &s).unwrap() < 1e-9);
    }

    #[test]
    fn pure_translation_is_removed() {
        let s = shape();
        let shifted: Vec<Vector3<f64>> =
            s.iter().map(|p| p + Vector3::new(5.0, -3.0, 2.0)).collect();
        // Kabsch removes the translation -> RMSD ~ 0.
        assert!(rmsd(&shifted, &s).unwrap() < 1e-9);
        // ...but the no-fit RMSD is large.
        assert!(rmsd_no_fit(&shifted, &s).unwrap() > 1.0);
    }

    #[test]
    fn pure_rotation_is_removed() {
        let s = shape();
        // Rotate 90 degrees about z.
        let rot = Matrix3::new(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let rotated: Vec<Vector3<f64>> = s.iter().map(|p| rot * p).collect();
        let result = kabsch(&rotated, &s).unwrap();
        // RMSD after optimal superposition should vanish.
        assert!(result.rmsd < 1e-9, "rmsd = {}", result.rmsd);
        // The recovered rotation is a proper rotation.
        assert!((result.rotation.determinant() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn superposition_apply_round_trips() {
        let s = shape();
        let rot = Matrix3::new(0.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0);
        let moved: Vec<Vector3<f64>> = s
            .iter()
            .map(|p| rot * p + Vector3::new(1.0, 2.0, 3.0))
            .collect();
        let sup = kabsch(&moved, &s).unwrap();
        let back = sup.apply_all(&moved);
        for (a, b) in back.iter().zip(&s) {
            assert!((a - b).norm() < 1e-9);
        }
    }

    #[test]
    fn rmsd_detects_a_real_difference() {
        let mut s = shape();
        let r = s.clone();
        // Genuinely move one atom.
        s[2] += Vector3::new(0.0, 0.5, 0.0);
        let d = rmsd(&s, &r).unwrap();
        assert!(d > 0.05 && d < 0.5, "rmsd = {d}");
    }

    #[test]
    fn rejects_size_mismatch() {
        let s = shape();
        let short = &s[..3];
        assert!(rmsd(short, &s).is_err());
        assert!(kabsch(&[], &[]).is_err());
    }

    #[test]
    fn rmsf_flags_the_mobile_atom() {
        // A trajectory where atom 2 wiggles and the rest are fixed.
        let base = shape();
        let mut frames = Vec::new();
        for f in 0..20 {
            let mut frame = base.clone();
            let wiggle = 0.1 * ((f as f64) * 0.7).sin();
            frame[2] += Vector3::new(0.0, wiggle, 0.0);
            frames.push(frame);
        }
        let fluct = rmsf(&frames, &base).unwrap();
        assert_eq!(fluct.len(), base.len());
        // Atom 2 fluctuates much more than the others.
        let others_max = fluct
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != 2)
            .map(|(_, v)| *v)
            .fold(0.0, f64::max);
        assert!(
            fluct[2] > others_max,
            "atom 2 RMSF {} not the largest",
            fluct[2]
        );
    }

    #[test]
    fn rmsf_is_zero_for_a_static_trajectory() {
        let base = shape();
        let frames = vec![base.clone(); 10];
        let fluct = rmsf(&frames, &base).unwrap();
        for f in fluct {
            assert!(f < 1e-9);
        }
    }
}
