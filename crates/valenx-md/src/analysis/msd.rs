//! Mean-squared displacement & diffusion — **roadmap feature 28**.
//!
//! The **mean-squared displacement** measures how far atoms wander
//! from where they started:
//!
//! ```text
//! MSD(τ) = ⟨ |rᵢ(t+τ) − rᵢ(t)|² ⟩
//! ```
//!
//! averaged over atoms `i` and over all time origins `t`. For a
//! diffusing system the MSD grows *linearly* with the lag `τ`, and the
//! slope gives the **self-diffusion coefficient** via the Einstein
//! relation:
//!
//! ```text
//! D = MSD(τ) / (2·d·τ)        d = dimensionality
//! ```
//!
//! — `d = 3` for normal three-dimensional diffusion, so `D = slope/6`.
//!
//! ## Unwrapped coordinates
//!
//! The MSD must be computed from **unwrapped** trajectories — if an
//! atom is wrapped back across a periodic boundary its displacement
//! would jump by a box length. [`unwrap_trajectory`] removes those
//! jumps by tracking each atom frame-to-frame; feed raw wrapped
//! frames through it first.

use nalgebra::Vector3;

use crate::error::{MdError, Result};
use crate::pbc::SimBox;

/// A mean-squared-displacement curve.
#[derive(Clone, Debug, PartialEq)]
pub struct MsdCurve {
    /// Lag times τ (ps).
    pub lag_times: Vec<f64>,
    /// MSD at each lag (nm²).
    pub msd: Vec<f64>,
}

impl MsdCurve {
    /// The Einstein self-diffusion coefficient (nm²/ps) from a linear
    /// fit of `MSD` vs `τ` over the central portion of the curve.
    ///
    /// The first and last fifths of the curve are dropped: the start
    /// is dominated by ballistic (non-diffusive) motion and the tail
    /// by poor statistics. `D = slope / (2·d)` with `d = 3`.
    pub fn diffusion_coefficient(&self) -> f64 {
        self.diffusion_coefficient_dim(3)
    }

    /// The diffusion coefficient for an arbitrary dimensionality `d`
    /// (`D = slope/(2·d)`).
    pub fn diffusion_coefficient_dim(&self, dimensionality: usize) -> f64 {
        let n = self.lag_times.len();
        if n < 4 || dimensionality == 0 {
            return 0.0;
        }
        let lo = n / 5;
        let hi = n - n / 5;
        let slope = linear_slope(&self.lag_times[lo..hi], &self.msd[lo..hi]);
        slope / (2.0 * dimensionality as f64)
    }
}

/// Least-squares slope of `y` against `x` (through a free intercept).
fn linear_slope(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    if n < 2.0 {
        return 0.0;
    }
    let sx: f64 = x.iter().sum();
    let sy: f64 = y.iter().sum();
    let sxx: f64 = x.iter().map(|v| v * v).sum();
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| a * b).sum();
    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-18 {
        0.0
    } else {
        (n * sxy - sx * sy) / denom
    }
}

/// Removes periodic-wrap jumps from a trajectory.
///
/// Given frames of *wrapped* positions, returns frames of continuous
/// (unwrapped) positions by adding, frame to frame, the minimum-image
/// displacement of every atom. The first frame is returned unchanged.
///
/// # Errors
/// [`MdError::Invalid`] for an empty trajectory or inconsistent frame
/// sizes.
pub fn unwrap_trajectory(
    frames: &[Vec<Vector3<f64>>],
    cell: &SimBox,
) -> Result<Vec<Vec<Vector3<f64>>>> {
    if frames.is_empty() {
        return Err(MdError::invalid("frames", "needs at least one frame"));
    }
    let natoms = frames[0].len();
    if frames.iter().any(|f| f.len() != natoms) {
        return Err(MdError::dimension("trajectory frames differ in atom count"));
    }
    let mut out = Vec::with_capacity(frames.len());
    out.push(frames[0].clone());
    for f in 1..frames.len() {
        let prev_unwrapped = &out[f - 1];
        let prev_wrapped = &frames[f - 1];
        let mut current = Vec::with_capacity(natoms);
        for a in 0..natoms {
            // Minimum-image step from the previous wrapped frame.
            let step = cell.min_image(frames[f][a] - prev_wrapped[a]);
            current.push(prev_unwrapped[a] + step);
        }
        out.push(current);
    }
    Ok(out)
}

/// Computes the MSD curve from an *unwrapped* trajectory.
///
/// All time origins are averaged (the standard "sliding-window"
/// estimator): the MSD at lag `k` averages over every pair of frames
/// `k` apart. `dt` is the time between stored frames (ps).
///
/// # Errors
/// [`MdError::Invalid`] for fewer than two frames, inconsistent frame
/// sizes, or a non-positive `dt`.
pub fn mean_squared_displacement(frames: &[Vec<Vector3<f64>>], dt: f64) -> Result<MsdCurve> {
    if frames.len() < 2 {
        return Err(MdError::invalid("frames", "needs at least two frames"));
    }
    if !(dt.is_finite() && dt > 0.0) {
        return Err(MdError::invalid("dt", "must be finite and positive"));
    }
    let natoms = frames[0].len();
    if natoms == 0 {
        return Err(MdError::invalid("frames", "frames have no atoms"));
    }
    if frames.iter().any(|f| f.len() != natoms) {
        return Err(MdError::dimension("trajectory frames differ in atom count"));
    }
    let n_frames = frames.len();
    let max_lag = n_frames - 1;
    let mut lag_times = Vec::with_capacity(max_lag);
    let mut msd = Vec::with_capacity(max_lag);
    for lag in 1..=max_lag {
        let mut sum = 0.0;
        let mut count = 0u64;
        for origin in 0..(n_frames - lag) {
            for a in 0..natoms {
                let d = frames[origin + lag][a] - frames[origin][a];
                sum += d.norm_squared();
                count += 1;
            }
        }
        lag_times.push(lag as f64 * dt);
        msd.push(if count > 0 { sum / count as f64 } else { 0.0 });
    }
    Ok(MsdCurve { lag_times, msd })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A trajectory of `n_frames` frames in which every atom moves at a
    /// constant velocity `v` per frame (ballistic — MSD ∝ τ²) or, when
    /// `diffusive_step` is used, by a fixed displacement per frame.
    fn ballistic(n_atoms: usize, n_frames: usize, v: Vector3<f64>) -> Vec<Vec<Vector3<f64>>> {
        (0..n_frames)
            .map(|f| {
                (0..n_atoms)
                    .map(|a| Vector3::new(a as f64 * 0.3, 0.0, 0.0) + v * f as f64)
                    .collect()
            })
            .collect()
    }

    #[test]
    fn rejects_bad_input() {
        let one = vec![vec![Vector3::zeros()]];
        assert!(mean_squared_displacement(&one, 0.01).is_err());
        let two = ballistic(2, 2, Vector3::new(0.1, 0.0, 0.0));
        assert!(mean_squared_displacement(&two, 0.0).is_err());
    }

    #[test]
    fn stationary_trajectory_has_zero_msd() {
        let frames = ballistic(10, 20, Vector3::zeros());
        let curve = mean_squared_displacement(&frames, 0.01).unwrap();
        for m in &curve.msd {
            assert!(m.abs() < 1e-12);
        }
    }

    #[test]
    fn ballistic_msd_grows_quadratically() {
        // Constant velocity -> displacement at lag k is k·v ->
        // MSD(k) = (k·|v|)².
        let v = Vector3::new(0.05, 0.0, 0.0);
        let frames = ballistic(5, 30, v);
        let curve = mean_squared_displacement(&frames, 1.0).unwrap();
        // MSD at lag 1 should be |v|² = 0.0025.
        assert!(
            (curve.msd[0] - 0.0025).abs() < 1e-9,
            "MSD(1) = {}",
            curve.msd[0]
        );
        // MSD at lag 4 should be (4|v|)² = 16·0.0025.
        assert!((curve.msd[3] - 16.0 * 0.0025).abs() < 1e-9);
    }

    #[test]
    fn unwrap_removes_boundary_jumps() {
        // An atom crossing a periodic boundary: wrapped coords jump,
        // unwrapped coords stay continuous.
        let cell = SimBox::cubic(1.0).unwrap();
        // Atom marches in +x by 0.4 nm per frame, wrapped into [0,1).
        let mut wrapped = Vec::new();
        for f in 0..5 {
            let x = (0.4 * f as f64) % 1.0;
            wrapped.push(vec![Vector3::new(x, 0.0, 0.0)]);
        }
        let unwrapped = unwrap_trajectory(&wrapped, &cell).unwrap();
        // Frame 4 unwrapped x should be ~1.6 nm, not wrapped 0.6.
        assert!(
            (unwrapped[4][0].x - 1.6).abs() < 1e-9,
            "x = {}",
            unwrapped[4][0].x
        );
    }

    #[test]
    fn diffusion_coefficient_recovers_known_slope() {
        // Build an MSD curve that is exactly linear: MSD = 6·D·τ with
        // D = 0.5 nm²/ps. The fit should recover D.
        let d_true = 0.5;
        let lag_times: Vec<f64> = (1..=50).map(|k| k as f64 * 0.1).collect();
        let msd: Vec<f64> = lag_times.iter().map(|t| 6.0 * d_true * t).collect();
        let curve = MsdCurve { lag_times, msd };
        let d = curve.diffusion_coefficient();
        assert!((d - d_true).abs() < 1e-6, "recovered D = {d}");
    }

    #[test]
    fn diffusion_coefficient_handles_short_curve() {
        let curve = MsdCurve {
            lag_times: vec![0.1, 0.2],
            msd: vec![0.1, 0.2],
        };
        assert_eq!(curve.diffusion_coefficient(), 0.0);
    }
}
