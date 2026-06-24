//! Reconstruction-quality tooling: metrics, outlier filtering, and export
//! (the **quality** stage that rounds out the in-house SfM port).
//!
//! Stages 1–5 plus the [`crate::mapper`] produce a [`Reconstruction`] — a set
//! of registered [`CameraPose`](crate::pnp::CameraPose)s and a sparse 3-D point
//! cloud tied to its 2-D observations by [`Track`](crate::mapper::Track)s. This
//! module is what you reach for *after* the
//! solver has run: it answers "**how good is this reconstruction, and can I
//! clean and export it?**". It mirrors the post-reconstruction utilities of a
//! production SfM tool (COLMAP's model statistics, observation filtering, and
//! its PLY / `cameras.txt` exporters), but is a small, self-contained,
//! dependency-light in-house implementation.
//!
//! Three capabilities:
//!
//! 1. **Quality metrics** ([`reconstruction_quality`] → [`QualityMetrics`]).
//!    The mean and median **reprojection error** (the pixel distance between
//!    each observed feature and the projection of its 3-D point through the
//!    observing camera), the **track-length** distribution (mean / max number
//!    of observations per triangulated point), the count of **registered
//!    cameras** and of **3-D points**, and the **camera connectivity** (how
//!    many image pairs share at least `N` triangulated points). These are the
//!    numbers you scan to decide whether a reconstruction is healthy: a low
//!    mean reprojection error and long, well-connected tracks indicate a strong
//!    model.
//!
//! 2. **Reprojection-outlier filtering** ([`filter_outliers`] →
//!    `(Reconstruction, FilterReport)`). Gross mismatches survive RANSAC and
//!    pollute a model with high-residual observations; this removes the worst
//!    of them. Two rules are offered ([`OutlierRule`]): a fixed pixel
//!    **threshold**, or a **robust** `k·MAD` rule (median absolute deviation —
//!    a heavy-tail-resistant analogue of "k standard deviations"). An
//!    observation whose reprojection error exceeds the cutoff is dropped; a
//!    track left with fewer than two observations (no longer triangulable) is
//!    dropped wholesale along with its 3-D point. The filtered
//!    [`Reconstruction`] is returned beside a [`FilterReport`] recording
//!    exactly what was removed and the before/after mean error.
//!
//! 3. **Export** ([`write_ply`] / `point_cloud_ply`, [`write_camera_poses`] /
//!    `camera_poses_text`). A standard **PLY** point cloud (ASCII, `x y z` with
//!    optional per-point `red green blue`) that any mesh viewer (MeshLab,
//!    CloudCompare, Blender) opens, and a plain-text **camera-pose** dump
//!    (per registered image: the rotation `R`, the translation `t`, and the
//!    derived camera centre `C = −Rᵀ t`).
//!
//! ## Reprojection-error convention
//!
//! Every metric and filter here is built on a single per-observation quantity:
//! for a triangulated track point `X` observed as pixel `x` by registered
//! camera `(R, t)`, the **reprojection error** is the Euclidean pixel distance
//! `‖ x − π(K[R|t]X) ‖` where `π` is the perspective projection
//! ([`crate::project_point`]). Observations of un-triangulated tracks, of
//! unregistered cameras, or with an out-of-range keypoint index contribute no
//! residual (there is nothing to compare against) and are skipped — the same
//! convention the mapper's own RMSE helper uses.
//!
//! ## Honesty notes
//!
//! - These are **diagnostic** metrics, not a calibrated accuracy certificate.
//!   A low reprojection error means the model is *self-consistent* (the points
//!   reproject onto their features); it does **not** by itself certify metric
//!   accuracy against ground-truth survey coordinates — the reconstruction is
//!   still only defined up to a global similarity (see [`crate::mapper`]).
//! - The `k·MAD` rule assumes the bulk of the residuals are unimodal; with a
//!   very small number of observations the MAD estimate is noisy (as for any
//!   robust scale estimate), so prefer a fixed threshold when you have only a
//!   handful of points.
//! - Colours are passed in by the caller (this crate is grayscale-only and
//!   does not sample image colour); the PLY writer simply serializes whatever
//!   RGB you supply, or writes an uncoloured cloud when you supply none.

use crate::error::{PhotogrammetryError, Result};
use crate::mapper::Reconstruction;
use crate::pnp::project_point;
use crate::twoview::CameraIntrinsics;
use crate::ImageFeatures;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::path::Path;

/// Pixel `(u, v)` of `images[image_idx].keypoints[kp]`, if both indices are in
/// range. Mirrors the mapper's private helper (kept local to this module).
fn pixel_of(images: &[ImageFeatures], image_idx: usize, kp: usize) -> Option<(f64, f64)> {
    let p = images.get(image_idx)?.keypoints.get(kp)?;
    Some((f64::from(p.x), f64::from(p.y)))
}

/// Per-observation reprojection errors (px) over every observation of every
/// **triangulated** track that is seen by a **registered** camera.
///
/// Each returned tuple is `(image_idx, keypoint_idx, error_px)`. Observations
/// that cannot yield a residual — un-triangulated track, unregistered camera,
/// out-of-range keypoint, or a point on the camera's principal plane
/// (`z ≈ 0`, where the projection is undefined) — are skipped. The output is
/// the raw material for [`QualityMetrics`] and [`filter_outliers`].
///
/// # Errors
///
/// Returns [`PhotogrammetryError::EmptyReconstruction`] when the reconstruction
/// has no registered cameras *and* no points (there is nothing to measure), and
/// [`PhotogrammetryError::NonFiniteValue`] if any stored 3-D point coordinate
/// is non-finite (`NaN`/`±∞`) — a corrupt model is reported loudly rather than
/// silently producing `NaN` statistics.
fn observation_errors(
    recon: &Reconstruction,
    images: &[ImageFeatures],
    k: &CameraIntrinsics,
) -> Result<Vec<(usize, usize, f64)>> {
    // Fail loud on a wholly empty reconstruction.
    if recon.num_registered() == 0 && recon.points.is_empty() {
        return Err(PhotogrammetryError::EmptyReconstruction);
    }
    // Fail loud on a corrupt (non-finite) point cloud before any averaging.
    for p in &recon.points {
        if !p.iter().all(|c| c.is_finite()) {
            return Err(PhotogrammetryError::NonFiniteValue);
        }
    }

    let mut out = Vec::new();
    for tr in &recon.tracks {
        let Some(pi) = tr.point_idx else { continue };
        let Some(&x) = recon.points.get(pi) else {
            continue;
        };
        for &(im, kp) in &tr.observations {
            let Some(Some(pose)) = recon.cameras.get(im) else {
                continue; // unregistered or out-of-range camera
            };
            let Some((ou, ov)) = pixel_of(images, im, kp) else {
                continue; // out-of-range keypoint
            };
            let Some((pu, pv)) = project_point(k, &pose.rotation, &pose.translation, &x) else {
                continue; // z ≈ 0, no finite image
            };
            let err = ((pu - ou).powi(2) + (pv - ov).powi(2)).sqrt();
            out.push((im, kp, err));
        }
    }
    Ok(out)
}

/// Quality metrics summarizing a [`Reconstruction`].
///
/// Produced by [`reconstruction_quality`]. All reprojection errors are in
/// pixels; track lengths are observation counts. The error/length fields are
/// computed only over observations that yield a residual (triangulated track,
/// registered camera) — see [`observation_errors`](self).
#[derive(Debug, Clone, PartialEq)]
pub struct QualityMetrics {
    /// Number of registered (posed) cameras.
    pub num_registered_cameras: usize,
    /// Number of triangulated 3-D points in the cloud.
    pub num_points: usize,
    /// Total number of valid 2-D observations contributing a reprojection
    /// residual (the sum of triangulated track lengths over registered views).
    pub num_observations: usize,
    /// Mean reprojection error over all valid observations, in pixels. `0.0`
    /// when there are no valid observations (guarded divide).
    pub mean_reprojection_error: f64,
    /// Median reprojection error over all valid observations, in pixels. `0.0`
    /// when there are no valid observations.
    pub median_reprojection_error: f64,
    /// Maximum reprojection error over all valid observations, in pixels.
    /// `0.0` when there are no valid observations.
    pub max_reprojection_error: f64,
    /// Mean track length: the average number of (registered-view) observations
    /// per triangulated point. `0.0` when there are no triangulated points.
    pub mean_track_length: f64,
    /// Maximum track length (most observations on any one triangulated point).
    /// `0` when there are no triangulated points.
    pub max_track_length: usize,
    /// Camera connectivity: the number of unordered image pairs that share at
    /// least [`QualityParams::min_shared_points`] triangulated points. A
    /// well-connected model has many such pairs.
    pub connected_pairs: usize,
}

/// Parameters for [`reconstruction_quality`].
#[derive(Debug, Clone, Copy)]
pub struct QualityParams {
    /// Two images count as a "connected pair" (toward
    /// [`QualityMetrics::connected_pairs`]) when they jointly observe at least
    /// this many **triangulated** points. A typical value is the minimum match
    /// support you would trust for a verified edge (e.g. 15–30).
    pub min_shared_points: usize,
}

impl Default for QualityParams {
    /// Default connectivity threshold of 15 shared triangulated points per
    /// image pair.
    fn default() -> Self {
        Self {
            min_shared_points: 15,
        }
    }
}

/// The median of a slice of errors. Returns `0.0` for an empty slice (the
/// callers treat "no observations" as a zero metric, guarded explicitly).
///
/// Sorts a copy by total order (NaNs cannot occur here — the inputs come from
/// [`observation_errors`], which has already rejected non-finite geometry —
/// but `total_cmp` keeps the sort total regardless).
fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut v: Vec<f64> = values.to_vec();
    v.sort_by(|a, b| a.total_cmp(b));
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        // Average the two central order statistics.
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// Compute [`QualityMetrics`] for a reconstruction.
///
/// `images` are the same per-image [`ImageFeatures`] the reconstruction was
/// built from (so observed pixels can be looked up), and `k` the shared
/// [`CameraIntrinsics`]. `params` tunes only the connectivity threshold.
///
/// Every average is **guarded**: the mean/median/max reprojection error are
/// `0.0` when there are no valid observations, and the mean track length is
/// `0.0` when there are no triangulated points — there is never a divide by a
/// zero count.
///
/// # Errors
///
/// - [`PhotogrammetryError::EmptyReconstruction`] when the reconstruction has
///   no registered cameras and no points.
/// - [`PhotogrammetryError::NonFiniteValue`] when a stored 3-D point is
///   non-finite.
pub fn reconstruction_quality(
    recon: &Reconstruction,
    images: &[ImageFeatures],
    k: &CameraIntrinsics,
    params: &QualityParams,
) -> Result<QualityMetrics> {
    let errors = observation_errors(recon, images, k)?;

    // Reprojection-error statistics (guarded divides).
    let num_observations = errors.len();
    let err_only: Vec<f64> = errors.iter().map(|&(_, _, e)| e).collect();
    let (mean_reprojection_error, median_reprojection_error, max_reprojection_error) =
        if num_observations == 0 {
            (0.0, 0.0, 0.0)
        } else {
            let sum: f64 = err_only.iter().sum();
            let mean = sum / num_observations as f64;
            let med = median(&err_only);
            let max = err_only.iter().copied().fold(0.0_f64, f64::max);
            (mean, med, max)
        };

    // Track-length distribution over TRIANGULATED tracks, counting only
    // observations by registered cameras (the same observations that can carry
    // a residual — a track's "useful length" in the posed model).
    let mut track_lengths: Vec<usize> = Vec::new();
    // Per-track list of the images that registered-observe it, for connectivity.
    let mut track_images: Vec<Vec<usize>> = Vec::new();
    for tr in &recon.tracks {
        if tr.point_idx.is_none() {
            continue;
        }
        let mut imgs: Vec<usize> = Vec::new();
        for &(im, kp) in &tr.observations {
            let registered = matches!(recon.cameras.get(im), Some(Some(_)));
            let kp_in_range = images.get(im).is_some_and(|f| kp < f.keypoints.len());
            if registered && kp_in_range {
                imgs.push(im);
            }
        }
        if imgs.is_empty() {
            continue;
        }
        track_lengths.push(imgs.len());
        track_images.push(imgs);
    }

    let num_triangulated = track_lengths.len();
    let (mean_track_length, max_track_length) = if num_triangulated == 0 {
        (0.0, 0)
    } else {
        let total: usize = track_lengths.iter().sum();
        let mean = total as f64 / num_triangulated as f64;
        let max = track_lengths.iter().copied().max().unwrap_or(0);
        (mean, max)
    };

    // Camera connectivity: count unordered image pairs sharing >= N points.
    // Accumulate per-pair shared-point counts over the triangulated tracks.
    let mut pair_counts: HashMap<(usize, usize), usize> = HashMap::new();
    for imgs in &track_images {
        // Each unordered pair of distinct images observing this track shares
        // this one point.
        for a in 0..imgs.len() {
            for b in (a + 1)..imgs.len() {
                let (lo, hi) = if imgs[a] < imgs[b] {
                    (imgs[a], imgs[b])
                } else {
                    (imgs[b], imgs[a])
                };
                *pair_counts.entry((lo, hi)).or_insert(0) += 1;
            }
        }
    }
    let connected_pairs = pair_counts
        .values()
        .filter(|&&c| c >= params.min_shared_points)
        .count();

    Ok(QualityMetrics {
        num_registered_cameras: recon.num_registered(),
        num_points: recon.points.len(),
        num_observations,
        mean_reprojection_error,
        median_reprojection_error,
        max_reprojection_error,
        mean_track_length,
        max_track_length,
        connected_pairs,
    })
}

/// The rule [`filter_outliers`] uses to decide which observations are outliers.
#[derive(Debug, Clone, Copy)]
pub enum OutlierRule {
    /// Drop every observation whose reprojection error exceeds this fixed pixel
    /// threshold. The threshold must be finite and non-negative (a negative or
    /// `NaN` threshold is rejected loudly by [`filter_outliers`]).
    Threshold(f64),
    /// Robust cutoff at `median + k · MAD`, where `MAD` is the median absolute
    /// deviation of the reprojection errors (scaled to be a consistent estimate
    /// of the standard deviation under normality by the usual `1.4826` factor).
    /// `k` is the multiplier (e.g. `3.0` ≈ "3 sigma"). `k` must be finite and
    /// non-negative.
    RobustMad(f64),
}

/// A record of what [`filter_outliers`] removed.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterReport {
    /// The effective pixel cutoff that was applied (for [`OutlierRule::Threshold`]
    /// this is the threshold itself; for [`OutlierRule::RobustMad`] it is the
    /// computed `median + k·MAD`).
    pub cutoff_px: f64,
    /// Number of individual observations removed (their reprojection error
    /// exceeded [`Self::cutoff_px`]).
    pub observations_removed: usize,
    /// Number of whole tracks (and their 3-D points) removed because filtering
    /// left them with fewer than two observations — no longer triangulable.
    pub tracks_removed: usize,
    /// Number of 3-D points in the cloud before filtering.
    pub points_before: usize,
    /// Number of 3-D points in the cloud after filtering.
    pub points_after: usize,
    /// Mean reprojection error (px) before filtering.
    pub mean_error_before: f64,
    /// Mean reprojection error (px) after filtering.
    pub mean_error_after: f64,
}

/// Remove high-reprojection-error observations (and any track they leave with
/// fewer than two surviving observations) from a reconstruction.
///
/// Returns the filtered [`Reconstruction`] **rebuilt with a compacted point
/// cloud** (surviving points are re-indexed contiguously and the tracks'
/// `point_idx` / `point3d` updated to match) together with a [`FilterReport`].
/// Cameras are left untouched (filtering removes structure, not poses).
///
/// The rule is chosen by [`OutlierRule`]: a fixed pixel `Threshold`, or the
/// robust `RobustMad(k)` cutoff `median + k·MAD`. An observation is removed
/// when its reprojection error is **strictly greater** than the cutoff. A
/// surviving track with at least two observations keeps its (unchanged) 3-D
/// point; a track reduced below two observations is dropped entirely (its point
/// leaves the cloud).
///
/// # Errors
///
/// - [`PhotogrammetryError::EmptyReconstruction`] when the reconstruction has
///   no registered cameras and no points.
/// - [`PhotogrammetryError::NonFiniteValue`] when a stored 3-D point is
///   non-finite, or when the rule parameter (`Threshold`'s threshold or
///   `RobustMad`'s `k`) is non-finite.
/// - [`PhotogrammetryError::NegativeThreshold`] when the rule parameter is
///   negative.
pub fn filter_outliers(
    recon: &Reconstruction,
    images: &[ImageFeatures],
    k: &CameraIntrinsics,
    rule: OutlierRule,
) -> Result<(Reconstruction, FilterReport)> {
    // Validate the rule parameter loudly (fail-loud on NaN / negative).
    let param = match rule {
        OutlierRule::Threshold(t) => t,
        OutlierRule::RobustMad(kk) => kk,
    };
    if !param.is_finite() {
        return Err(PhotogrammetryError::NonFiniteValue);
    }
    if param < 0.0 {
        return Err(PhotogrammetryError::NegativeThreshold);
    }

    // Per-observation errors (also validates emptiness / finite geometry).
    let errors = observation_errors(recon, images, k)?;
    let mean_error_before = if errors.is_empty() {
        0.0
    } else {
        errors.iter().map(|&(_, _, e)| e).sum::<f64>() / errors.len() as f64
    };

    // Resolve the rule to a concrete pixel cutoff.
    let cutoff_px = match rule {
        OutlierRule::Threshold(t) => t,
        OutlierRule::RobustMad(kk) => {
            let errs: Vec<f64> = errors.iter().map(|&(_, _, e)| e).collect();
            let med = median(&errs);
            // MAD = median(|e_i − median|), scaled by 1.4826 for σ-consistency.
            let abs_dev: Vec<f64> = errs.iter().map(|e| (e - med).abs()).collect();
            let mad = median(&abs_dev);
            med + kk * 1.4826 * mad
        }
    };

    // Build a fast lookup of (image, keypoint) -> error so we can test each
    // observation. A track observes any image at most once, so (im, kp) keys
    // are unique per residual.
    let mut err_of: HashMap<(usize, usize), f64> = HashMap::with_capacity(errors.len());
    for &(im, kp, e) in &errors {
        err_of.insert((im, kp), e);
    }

    // Rebuild the reconstruction: keep cameras as-is; for each track drop the
    // outlier observations, then keep the track (with its point) only if >= 2
    // observations survive.
    let mut new_points: Vec<nalgebra::Vector3<f64>> = Vec::new();
    let mut new_tracks: Vec<crate::mapper::Track> = Vec::with_capacity(recon.tracks.len());
    let mut observations_removed = 0usize;
    let mut tracks_removed = 0usize;

    for tr in &recon.tracks {
        // Untriangulated tracks are passed through unchanged (no residual to
        // test against; they carry no point).
        let Some(_pi) = tr.point_idx else {
            new_tracks.push(crate::mapper::Track {
                observations: tr.observations.clone(),
                point3d: None,
                point_idx: None,
            });
            continue;
        };

        // Keep observations whose error is within the cutoff. An observation
        // with no residual (e.g. its camera is unregistered) is retained: it
        // was never an outlier candidate, and dropping it could needlessly
        // dissolve a track.
        let mut kept_obs: Vec<(usize, usize)> = Vec::with_capacity(tr.observations.len());
        for &(im, kp) in &tr.observations {
            match err_of.get(&(im, kp)) {
                Some(&e) if e > cutoff_px => {
                    observations_removed += 1; // an outlier: drop it
                }
                _ => kept_obs.push((im, kp)),
            }
        }

        if kept_obs.len() < 2 {
            // No longer triangulable: drop the whole track and its point.
            tracks_removed += 1;
            // Any observations it still held that were NOT counted as outliers
            // are effectively removed too, but the report counts only the
            // reprojection-outlier removals above (a dropped sub-2 track is
            // reported via tracks_removed); this keeps observations_removed
            // meaning "removed for exceeding the cutoff".
            continue;
        }

        // Survivor: re-home its 3-D point into the compacted cloud.
        let point = tr.point3d.unwrap_or_else(|| {
            // point_idx was Some, so points[pi] is the canonical value; fall
            // back to it if the cached point3d was somehow absent.
            recon.points[_pi]
        });
        let new_idx = new_points.len();
        new_points.push(point);
        new_tracks.push(crate::mapper::Track {
            observations: kept_obs,
            point3d: Some(point),
            point_idx: Some(new_idx),
        });
    }

    let filtered = Reconstruction {
        cameras: recon.cameras.clone(),
        points: new_points,
        tracks: new_tracks,
    };

    // Mean error after, recomputed on the filtered model.
    let after_errors = observation_errors(&filtered, images, k)?;
    let mean_error_after = if after_errors.is_empty() {
        0.0
    } else {
        after_errors.iter().map(|&(_, _, e)| e).sum::<f64>() / after_errors.len() as f64
    };

    let report = FilterReport {
        cutoff_px,
        observations_removed,
        tracks_removed,
        points_before: recon.points.len(),
        points_after: filtered.points.len(),
        mean_error_before,
        mean_error_after,
    };
    Ok((filtered, report))
}

/// Serialize a reconstruction's point cloud to an ASCII **PLY** string.
///
/// Emits one vertex per triangulated 3-D point in [`Reconstruction::points`],
/// in their stored order. When `colors` is `Some`, it must have exactly one
/// `[r, g, b]` (`u8`) entry per point and the PLY gains `red`/`green`/`blue`
/// properties; when `None`, an uncoloured `x y z` cloud is written.
///
/// The output is the standard ASCII PLY any mesh viewer reads, e.g.
///
/// ```text
/// ply
/// format ascii 1.0
/// element vertex 2
/// property float x
/// property float y
/// property float z
/// end_header
/// 0.10000000 0.20000000 0.30000000
/// 0.40000000 0.50000000 0.60000000
/// ```
///
/// Coordinates are written with full `f64` round-trip precision (`{:.17}`), so
/// parsing them back recovers the originals to within floating-point epsilon.
///
/// # Errors
///
/// - [`PhotogrammetryError::EmptyReconstruction`] when the point cloud is
///   empty (nothing to export).
/// - [`PhotogrammetryError::NonFiniteValue`] when any point coordinate is
///   non-finite.
/// - [`PhotogrammetryError::ColorCountMismatch`] when `colors` is `Some` but
///   its length does not equal the number of points.
pub fn point_cloud_ply(recon: &Reconstruction, colors: Option<&[[u8; 3]]>) -> Result<String> {
    let n = recon.points.len();
    if n == 0 {
        return Err(PhotogrammetryError::EmptyReconstruction);
    }
    for p in &recon.points {
        if !p.iter().all(|c| c.is_finite()) {
            return Err(PhotogrammetryError::NonFiniteValue);
        }
    }
    if let Some(c) = colors {
        if c.len() != n {
            return Err(PhotogrammetryError::ColorCountMismatch {
                expected: n,
                actual: c.len(),
            });
        }
    }

    let mut s = String::new();
    s.push_str("ply\n");
    s.push_str("format ascii 1.0\n");
    s.push_str("comment generated by valenx-photogrammetry\n");
    let _ = writeln!(s, "element vertex {n}");
    s.push_str("property float x\n");
    s.push_str("property float y\n");
    s.push_str("property float z\n");
    if colors.is_some() {
        s.push_str("property uchar red\n");
        s.push_str("property uchar green\n");
        s.push_str("property uchar blue\n");
    }
    s.push_str("end_header\n");

    for (i, p) in recon.points.iter().enumerate() {
        // 17 significant digits round-trips an f64 exactly.
        let _ = write!(s, "{:.17} {:.17} {:.17}", p.x, p.y, p.z);
        if let Some(c) = colors {
            let rgb = c[i];
            let _ = write!(s, " {} {} {}", rgb[0], rgb[1], rgb[2]);
        }
        s.push('\n');
    }
    Ok(s)
}

/// Write a reconstruction's point cloud to a PLY file at `path`.
///
/// A thin file-writing wrapper over [`point_cloud_ply`]; see it for the format,
/// the `colors` contract, and the validation errors.
///
/// # Errors
///
/// The same model/colour validation errors as [`point_cloud_ply`], plus
/// [`PhotogrammetryError::Io`] if the file cannot be written.
pub fn write_ply(
    recon: &Reconstruction,
    colors: Option<&[[u8; 3]]>,
    path: impl AsRef<Path>,
) -> Result<()> {
    let body = point_cloud_ply(recon, colors)?;
    std::fs::write(path, body).map_err(|e| PhotogrammetryError::Io(e.to_string()))
}

/// Serialize the registered camera poses to a plain-text string.
///
/// One block per registered camera (unregistered images are skipped), giving
/// the image index, the rotation matrix `R` (row-major, `X_cam = R·X_world + t`),
/// the translation `t`, and the world-frame camera centre `C = −Rᵀ t`. The
/// format is human-readable and stable, intended for inspection / simple
/// downstream parsing rather than as a standardized interchange format.
///
/// # Errors
///
/// - [`PhotogrammetryError::EmptyReconstruction`] when no camera is registered.
pub fn camera_poses_text(recon: &Reconstruction) -> Result<String> {
    if recon.num_registered() == 0 {
        return Err(PhotogrammetryError::EmptyReconstruction);
    }
    let mut s = String::new();
    s.push_str("# valenx-photogrammetry camera poses\n");
    s.push_str("# convention: X_cam = R * X_world + t ; centre C = -R^T t\n");
    let _ = writeln!(s, "# registered_cameras {}", recon.num_registered());
    for (im, cam) in recon.cameras.iter().enumerate() {
        let Some(pose) = cam else { continue };
        let r = &pose.rotation;
        let t = &pose.translation;
        let c = -(r.transpose() * t);
        let _ = writeln!(s, "camera {im}");
        let _ = writeln!(s, "R {:.17} {:.17} {:.17}", r[(0, 0)], r[(0, 1)], r[(0, 2)]);
        let _ = writeln!(s, "R {:.17} {:.17} {:.17}", r[(1, 0)], r[(1, 1)], r[(1, 2)]);
        let _ = writeln!(s, "R {:.17} {:.17} {:.17}", r[(2, 0)], r[(2, 1)], r[(2, 2)]);
        let _ = writeln!(s, "t {:.17} {:.17} {:.17}", t.x, t.y, t.z);
        let _ = writeln!(s, "C {:.17} {:.17} {:.17}", c.x, c.y, c.z);
    }
    Ok(s)
}

/// Write the registered camera poses to a text file at `path`.
///
/// A thin file-writing wrapper over [`camera_poses_text`]; see it for the
/// format.
///
/// # Errors
///
/// [`PhotogrammetryError::EmptyReconstruction`] when no camera is registered,
/// or [`PhotogrammetryError::Io`] if the file cannot be written.
pub fn write_camera_poses(recon: &Reconstruction, path: impl AsRef<Path>) -> Result<()> {
    let body = camera_poses_text(recon)?;
    std::fs::write(path, body).map_err(|e| PhotogrammetryError::Io(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mapper::{Reconstruction, Track};
    use crate::pnp::CameraPose;
    use crate::Keypoint;
    use nalgebra::{Matrix3, Vector3};

    fn test_intrinsics() -> CameraIntrinsics {
        CameraIntrinsics::new(800.0, 800.0, 320.0, 240.0)
    }

    /// Project a world point through `K, R, t` to a pixel (test-local truth).
    fn project(k: &CameraIntrinsics, pose: &CameraPose, x: &Vector3<f64>) -> (f64, f64) {
        let cam = pose.rotation * x + pose.translation;
        (k.fx * cam.x / cam.z + k.cx, k.fy * cam.y / cam.z + k.cy)
    }

    /// A small known camera rig (identity + two offset poses), all looking down
    /// +z so the scene points (positive z) are in front of every camera.
    fn rig() -> Vec<CameraPose> {
        vec![
            CameraPose {
                rotation: Matrix3::identity(),
                translation: Vector3::zeros(),
            },
            CameraPose {
                rotation: Matrix3::identity(),
                translation: Vector3::new(-1.0, 0.0, 0.0),
            },
            CameraPose {
                rotation: Matrix3::identity(),
                translation: Vector3::new(0.6, -0.2, 0.1),
            },
        ]
    }

    /// Eight known 3-D points spanning depth, all in front of the rig.
    fn points() -> Vec<Vector3<f64>> {
        [
            (-0.4, -0.3, 5.0),
            (0.3, -0.2, 6.0),
            (0.1, 0.4, 4.5),
            (-0.2, 0.25, 7.0),
            (0.45, 0.35, 5.5),
            (-0.35, -0.15, 6.5),
            (0.2, -0.35, 4.8),
            (-0.1, 0.1, 8.0),
        ]
        .iter()
        .map(|&(x, y, z)| Vector3::new(x, y, z))
        .collect()
    }

    /// Build a SYNTHETIC, noise-free reconstruction: every camera observes
    /// every point, keypoint index == point index, points triangulated exactly
    /// at truth, every track linked. Returns the reconstruction and the images.
    fn synth_reconstruction(
        k: &CameraIntrinsics,
        cams: &[CameraPose],
        pts: &[Vector3<f64>],
    ) -> (Reconstruction, Vec<ImageFeatures>) {
        // Per-image keypoints (one per point, in point order).
        let images: Vec<ImageFeatures> = cams
            .iter()
            .map(|cam| {
                let kps = pts
                    .iter()
                    .map(|x| {
                        let (u, v) = project(k, cam, x);
                        Keypoint::new(u as f32, v as f32, 1.0)
                    })
                    .collect();
                ImageFeatures::new(kps)
            })
            .collect();

        // One track per point, observed by every camera, triangulated at truth.
        let tracks: Vec<Track> = pts
            .iter()
            .enumerate()
            .map(|(p, &x)| Track {
                observations: (0..cams.len()).map(|im| (im, p)).collect(),
                point3d: Some(x),
                point_idx: Some(p),
            })
            .collect();

        let recon = Reconstruction {
            cameras: cams.iter().copied().map(Some).collect(),
            points: pts.to_vec(),
            tracks,
        };
        (recon, images)
    }

    // ===== BENCHMARK-PIN 1: synthetic noise-free metrics are exact. =========
    #[test]
    fn synthetic_metrics_are_exact() {
        let k = test_intrinsics();
        let cams = rig();
        let pts = points();
        let (recon, images) = synth_reconstruction(&k, &cams, &pts);

        let metrics =
            reconstruction_quality(&recon, &images, &k, &QualityParams::default()).unwrap();

        // Mean reprojection error is ~0 on a noise-free synthetic scene. The
        // residual floor here is set by [`Keypoint`] storing pixels as `f32`
        // (the algorithm itself is exact): an `f32` ulp at a pixel magnitude of
        // a few hundred is ~3e-5 px, so a noise-free scene reprojects to within
        // that quantization, not literally 0. We pin tightly below 1e-4 — far
        // below any real measurement noise — and use the strict <1e-6 round-trip
        // pin where the value is genuinely `f64` (the PLY export test).
        assert!(
            metrics.mean_reprojection_error < 1e-4,
            "noise-free mean reprojection error must be ~0 (f32-keypoint floor), got {}",
            metrics.mean_reprojection_error
        );
        assert!(metrics.median_reprojection_error < 1e-4);
        assert!(metrics.max_reprojection_error < 1e-4);

        // Counts are EXACT.
        assert_eq!(metrics.num_registered_cameras, cams.len());
        assert_eq!(metrics.num_points, pts.len());
        // Every point seen by every camera => observations = points * cameras.
        assert_eq!(metrics.num_observations, pts.len() * cams.len());

        // Track lengths: every track has exactly `cams.len()` observations.
        assert!((metrics.mean_track_length - cams.len() as f64).abs() < 1e-12);
        assert_eq!(metrics.max_track_length, cams.len());

        // Connectivity: all 3 image pairs share all 8 points (>= 15? no — 8).
        // With the default threshold (15) NO pair qualifies; with threshold 8
        // all C(3,2)=3 pairs qualify. Check both to pin the semantics.
        assert_eq!(metrics.connected_pairs, 0, "8 < 15 default threshold");
        let m8 = reconstruction_quality(
            &recon,
            &images,
            &k,
            &QualityParams {
                min_shared_points: 8,
            },
        )
        .unwrap();
        assert_eq!(m8.connected_pairs, 3, "all 3 pairs share 8 >= 8 points");
        // And with threshold 9, none (only 8 shared).
        let m9 = reconstruction_quality(
            &recon,
            &images,
            &k,
            &QualityParams {
                min_shared_points: 9,
            },
        )
        .unwrap();
        assert_eq!(m9.connected_pairs, 0);
    }

    // ===== BENCHMARK-PIN 2: injecting outliers then filtering removes exactly
    //        those and lowers the mean reprojection error. ====================
    #[test]
    fn outlier_injection_then_filter() {
        let k = test_intrinsics();
        let cams = rig();
        let pts = points();
        let (recon, mut images) = synth_reconstruction(&k, &cams, &pts);

        // Baseline: clean model, ~0 mean error (f32-keypoint floor; see
        // `synthetic_metrics_are_exact` for why this is ~1e-5, not literally 0).
        let clean = reconstruction_quality(&recon, &images, &k, &QualityParams::default()).unwrap();
        assert!(clean.mean_reprojection_error < 1e-4);

        // Inject 3 GROSS outliers: perturb specific observed pixels far from
        // where their 3-D point reprojects. We corrupt keypoint p in image im.
        // Track p observes (im, p); move that keypoint by a big offset.
        let gross: [(usize, usize); 3] = [(1, 2), (2, 5), (0, 7)];
        for &(im, p) in &gross {
            let kp = &mut images[im].keypoints[p];
            kp.x += 50.0; // 50 px gross error, way above any noise
            kp.y -= 40.0;
        }

        // Now the mean error is clearly raised.
        let dirty = reconstruction_quality(&recon, &images, &k, &QualityParams::default()).unwrap();
        assert!(
            dirty.mean_reprojection_error > 1.0,
            "gross outliers must raise the mean error, got {}",
            dirty.mean_reprojection_error
        );

        // Filter with a fixed threshold (10 px): the 3 gross observations
        // (~64 px error) exceed it; all clean observations (~0) survive.
        let (filtered, report) =
            filter_outliers(&recon, &images, &k, OutlierRule::Threshold(10.0)).unwrap();
        assert_eq!(
            report.observations_removed, 3,
            "exactly the 3 injected outliers must be removed, got {}",
            report.observations_removed
        );
        // Each corrupted track still has 2 clean observations (3 cams − 1),
        // so NO track is dropped and all 8 points remain.
        assert_eq!(report.tracks_removed, 0);
        assert_eq!(report.points_after, pts.len());

        // The mean error after filtering drops back to ~0.
        assert!(
            report.mean_error_after < report.mean_error_before,
            "filtering must lower the mean error ({} -> {})",
            report.mean_error_before,
            report.mean_error_after
        );
        assert!(
            report.mean_error_after < 1e-4,
            "post-filter mean error must be ~0 (f32-keypoint floor), got {}",
            report.mean_error_after
        );

        // The robust k*MAD rule also catches the same 3 (median/MAD of the
        // residuals are ~0, so any gross error blows past median + k*MAD).
        let (_f2, r2) = filter_outliers(&recon, &images, &k, OutlierRule::RobustMad(5.0)).unwrap();
        assert_eq!(
            r2.observations_removed, 3,
            "MAD rule removes the 3 outliers"
        );

        // Sanity: re-running the metrics through filtered recon agrees with the
        // report's after-error.
        let post =
            reconstruction_quality(&filtered, &images, &k, &QualityParams::default()).unwrap();
        assert!((post.mean_reprojection_error - report.mean_error_after).abs() < 1e-12);
    }

    // ===== A filter that strips a track below 2 observations drops it. =======
    #[test]
    fn filter_drops_undersupported_track() {
        let k = test_intrinsics();
        // Two cameras only, so a track has exactly 2 observations; corrupting
        // ONE of them and filtering leaves 1 -> the whole track is dropped.
        let cams = vec![rig()[0], rig()[1]];
        let pts = points();
        let (recon, mut images) = synth_reconstruction(&k, &cams, &pts);

        // Corrupt one observation of track 3 (image 0, keypoint 3).
        images[0].keypoints[3].x += 80.0;

        let (filtered, report) =
            filter_outliers(&recon, &images, &k, OutlierRule::Threshold(10.0)).unwrap();
        assert_eq!(report.observations_removed, 1);
        assert_eq!(report.tracks_removed, 1, "track 3 drops below 2 obs");
        assert_eq!(report.points_before, pts.len());
        assert_eq!(report.points_after, pts.len() - 1);
        // The filtered point cloud is compacted (indices contiguous, all valid).
        for tr in &filtered.tracks {
            if let Some(pi) = tr.point_idx {
                assert!(pi < filtered.points.len(), "compacted index in range");
            }
        }
    }

    // ===== BENCHMARK-PIN 3: PLY export round-trips the point coordinates. ====
    #[test]
    fn ply_roundtrips_coordinates() {
        let k = test_intrinsics();
        let cams = rig();
        let pts = points();
        let (recon, _images) = synth_reconstruction(&k, &cams, &pts);

        let ply = point_cloud_ply(&recon, None).unwrap();

        // Parse the vertices back out of the ASCII PLY and compare < 1e-6.
        let parsed = parse_ply_xyz(&ply);
        assert_eq!(parsed.len(), pts.len(), "vertex count must match");
        for (orig, got) in pts.iter().zip(parsed.iter()) {
            assert!(
                (orig.x - got[0]).abs() < 1e-6
                    && (orig.y - got[1]).abs() < 1e-6
                    && (orig.z - got[2]).abs() < 1e-6,
                "PLY xyz round-trip mismatch: {orig:?} vs {got:?}"
            );
        }

        // The coloured variant emits the rgb columns and still round-trips xyz.
        let colors: Vec<[u8; 3]> = (0..pts.len())
            .map(|i| [(i * 10) as u8, (i * 20) as u8, (i * 30) as u8])
            .collect();
        let ply_rgb = point_cloud_ply(&recon, Some(&colors)).unwrap();
        assert!(ply_rgb.contains("property uchar red"));
        let (parsed_xyz, parsed_rgb) = parse_ply_xyz_rgb(&ply_rgb);
        assert_eq!(parsed_xyz.len(), pts.len());
        for (orig, got) in pts.iter().zip(parsed_xyz.iter()) {
            assert!((orig.x - got[0]).abs() < 1e-6);
            assert!((orig.y - got[1]).abs() < 1e-6);
            assert!((orig.z - got[2]).abs() < 1e-6);
        }
        assert_eq!(parsed_rgb, colors, "rgb columns must round-trip exactly");
    }

    /// Minimal ASCII-PLY vertex parser (xyz only) for the round-trip test.
    fn parse_ply_xyz(ply: &str) -> Vec<[f64; 3]> {
        let mut in_body = false;
        let mut out = Vec::new();
        for line in ply.lines() {
            if in_body {
                let nums: Vec<f64> = line
                    .split_whitespace()
                    .filter_map(|t| t.parse::<f64>().ok())
                    .collect();
                if nums.len() >= 3 {
                    out.push([nums[0], nums[1], nums[2]]);
                }
            } else if line.trim() == "end_header" {
                in_body = true;
            }
        }
        out
    }

    /// ASCII-PLY parser returning xyz and rgb columns.
    fn parse_ply_xyz_rgb(ply: &str) -> (Vec<[f64; 3]>, Vec<[u8; 3]>) {
        let mut in_body = false;
        let mut xyz = Vec::new();
        let mut rgb = Vec::new();
        for line in ply.lines() {
            if in_body {
                let toks: Vec<&str> = line.split_whitespace().collect();
                if toks.len() >= 6 {
                    xyz.push([
                        toks[0].parse().unwrap(),
                        toks[1].parse().unwrap(),
                        toks[2].parse().unwrap(),
                    ]);
                    rgb.push([
                        toks[3].parse().unwrap(),
                        toks[4].parse().unwrap(),
                        toks[5].parse().unwrap(),
                    ]);
                }
            } else if line.trim() == "end_header" {
                in_body = true;
            }
        }
        (xyz, rgb)
    }

    // ===== Camera-pose text export contains every registered camera. ========
    #[test]
    fn camera_pose_text_lists_all() {
        let k = test_intrinsics();
        let cams = rig();
        let pts = points();
        let (recon, _images) = synth_reconstruction(&k, &cams, &pts);

        let txt = camera_poses_text(&recon).unwrap();
        for im in 0..cams.len() {
            assert!(txt.contains(&format!("camera {im}")), "missing camera {im}");
        }
        // The export reports the right number of registered cameras and a line
        // per camera's rotation (3), translation (1) and centre (1).
        assert!(txt.contains(&format!("# registered_cameras {}", cams.len())));
        assert_eq!(
            txt.lines().filter(|l| l.starts_with("camera ")).count(),
            cams.len()
        );
        // Camera 0's centre is the origin (identity pose, zero translation):
        // C = -R^T t = 0 in every coordinate (modulo the sign of zero).
        let c_lines: Vec<&str> = txt.lines().filter(|l| l.starts_with("C ")).collect();
        let first_c: Vec<f64> = c_lines[0]
            .split_whitespace()
            .skip(1)
            .map(|t| t.parse::<f64>().unwrap())
            .collect();
        assert_eq!(first_c.len(), 3);
        assert!(
            first_c.iter().all(|&v| v.abs() < 1e-12),
            "cam0 centre is origin"
        );
    }

    // ===== FAIL-LOUD: empty reconstruction, NaN, negative threshold. =========
    #[test]
    fn fail_loud_paths() {
        let k = test_intrinsics();

        // Empty reconstruction (no cameras, no points) -> EmptyReconstruction.
        let empty = Reconstruction {
            cameras: vec![None, None],
            points: vec![],
            tracks: vec![],
        };
        let images = vec![ImageFeatures::new(vec![]), ImageFeatures::new(vec![])];
        assert!(matches!(
            reconstruction_quality(&empty, &images, &k, &QualityParams::default()),
            Err(PhotogrammetryError::EmptyReconstruction)
        ));
        assert!(matches!(
            point_cloud_ply(&empty, None),
            Err(PhotogrammetryError::EmptyReconstruction)
        ));
        assert!(matches!(
            camera_poses_text(&empty),
            Err(PhotogrammetryError::EmptyReconstruction)
        ));

        // NaN point coordinate -> NonFiniteValue.
        let cams = rig();
        let pts = points();
        let (mut recon, imgs) = synth_reconstruction(&k, &cams, &pts);
        recon.points[0].x = f64::NAN;
        assert!(matches!(
            reconstruction_quality(&recon, &imgs, &k, &QualityParams::default()),
            Err(PhotogrammetryError::NonFiniteValue)
        ));
        assert!(matches!(
            point_cloud_ply(&recon, None),
            Err(PhotogrammetryError::NonFiniteValue)
        ));

        // Negative threshold -> NegativeThreshold.
        let (good, gimgs) = synth_reconstruction(&k, &cams, &pts);
        assert!(matches!(
            filter_outliers(&good, &gimgs, &k, OutlierRule::Threshold(-1.0)),
            Err(PhotogrammetryError::NegativeThreshold)
        ));
        // NaN threshold -> NonFiniteValue.
        assert!(matches!(
            filter_outliers(&good, &gimgs, &k, OutlierRule::RobustMad(f64::NAN)),
            Err(PhotogrammetryError::NonFiniteValue)
        ));

        // Colour-count mismatch on PLY export.
        assert!(matches!(
            point_cloud_ply(&good, Some(&[[0, 0, 0]])),
            Err(PhotogrammetryError::ColorCountMismatch { .. })
        ));
    }
}
