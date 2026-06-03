//! Phase 191 — exploded-view animation — parts move outward then back
//! to reveal the assembly's internal structure.
//!
//! ## What OCCT does
//!
//! Real OCCT exposes `AIS_Animation_Object` chained per part with a
//! pre-computed "explode vector" pointing from the assembly centroid
//! to each part's centroid. The animation lerp is on the magnitude:
//! `t=0` keeps parts at rest, `t=1` moves them along their explode
//! vectors by `factor × |explode_vec|`. Standard explode factor is 2x
//! (parts end up at twice the centroid distance), with linear or
//! ease-in-out timing.
//!
//! ## v1 status
//!
//! **Honest v1.** Given a list of per-part centroids + an
//! `explode_factor`, returns each part's lerp'd translation at
//! parameter `t ∈ [0, 1]`. The caller composes the result with the
//! part's resting pose to produce the per-frame model matrix. Linear
//! timing only (ease-in-out is a 1-line `t = t*t*(3-2*t)` change the
//! caller can apply pre-call if they want it; Phase 191.5 will expose
//! the timing function as a parameter).

use nalgebra::Vector3;

use crate::error::OcctVizError;

/// Compute the explode-offset translation for each part at parameter
/// `t`. Returns one [`Vector3`] per input centroid — add it to the
/// part's resting world position to get the exploded position.
///
/// `centroids[i]` — world-space centre of part `i`.
/// `assembly_centroid` — world-space centre of the whole assembly
///   (typically the volume-weighted mean of per-part centroids; the
///   caller pre-computes this).
/// `factor` — explode magnitude; 0 = no explode, 1 = parts at twice
///   the centroid distance.
/// `t` — animation parameter in `[0, 1]` (clamped).
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if `t` or `factor` is non-finite, or
///   any centroid component is non-finite.
pub fn view_animation_explode(
    centroids: &[Vector3<f32>],
    assembly_centroid: Vector3<f32>,
    factor: f32,
    t: f32,
) -> Result<Vec<Vector3<f32>>, OcctVizError> {
    if !t.is_finite() {
        return Err(OcctVizError::bad_input("t", "must be finite"));
    }
    if !factor.is_finite() {
        return Err(OcctVizError::bad_input("factor", "must be finite"));
    }
    for (i, c) in centroids.iter().enumerate() {
        if !c.x.is_finite() || !c.y.is_finite() || !c.z.is_finite() {
            return Err(OcctVizError::bad_input(
                "centroids",
                format!("centroids[{i}] contains non-finite component"),
            ));
        }
    }
    if !assembly_centroid.x.is_finite()
        || !assembly_centroid.y.is_finite()
        || !assembly_centroid.z.is_finite()
    {
        return Err(OcctVizError::bad_input(
            "assembly_centroid",
            "must be finite",
        ));
    }
    let t = t.clamp(0.0, 1.0);
    Ok(centroids
        .iter()
        .map(|c| (c - assembly_centroid) * factor * t)
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_nan_t() {
        let err = view_animation_explode(&[], Vector3::zeros(), 1.0, f32::NAN).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn t_zero_keeps_parts_at_rest() {
        let cs = vec![Vector3::new(10.0, 0.0, 0.0), Vector3::new(-10.0, 0.0, 0.0)];
        let out = view_animation_explode(&cs, Vector3::zeros(), 1.0, 0.0).unwrap();
        assert!(out[0].norm() < 1e-4);
        assert!(out[1].norm() < 1e-4);
    }

    #[test]
    fn t_one_full_explode() {
        let cs = vec![Vector3::new(10.0, 0.0, 0.0), Vector3::new(-10.0, 0.0, 0.0)];
        let out = view_animation_explode(&cs, Vector3::zeros(), 1.0, 1.0).unwrap();
        // Part at (10,0,0) moves +10 along x (away from centroid).
        assert!((out[0].x - 10.0).abs() < 1e-4);
        // Part at (-10,0,0) moves -10 along x (away from centroid).
        assert!((out[1].x - (-10.0)).abs() < 1e-4);
    }

    #[test]
    fn factor_scales_magnitude() {
        let cs = vec![Vector3::new(10.0, 0.0, 0.0)];
        let half = view_animation_explode(&cs, Vector3::zeros(), 0.5, 1.0).unwrap();
        let double = view_animation_explode(&cs, Vector3::zeros(), 2.0, 1.0).unwrap();
        assert!((half[0].x - 5.0).abs() < 1e-4);
        assert!((double[0].x - 20.0).abs() < 1e-4);
    }

    #[test]
    fn empty_input_is_valid() {
        let out = view_animation_explode(&[], Vector3::zeros(), 1.0, 0.5).unwrap();
        assert!(out.is_empty());
    }
}
