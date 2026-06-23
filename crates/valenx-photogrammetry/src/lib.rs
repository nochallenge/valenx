//! # valenx-photogrammetry — structure-from-motion building blocks (Stages 1–4)
//!
//! An **in-house, clean-room Rust** reimplementation of the front end of a
//! classical Structure-from-Motion (SfM) pipeline, of the kind made famous
//! by COLMAP. The full pipeline (built out over later stages) is:
//!
//! ```text
//! images → feature detection & description   ← STAGE 1
//!        → descriptor matching                ← STAGE 2
//!        → two-view geometry verification     ← STAGE 2
//!        → relative pose (R, t) + triangulation   ← STAGE 3
//!        → incremental mapper: register new view (PnP)  ← STAGE 4 (this stage)
//!        → bundle adjustment
//!        → sparse point cloud + camera poses
//! ```
//!
//! **Stage 1** covers the first arrow: the grayscale [`GrayImage`] input
//! type, FAST-9 corner detection (the [`fast`] module), and an ORB-style
//! oriented binary descriptor (the [`descriptor`] module).
//!
//! **Stage 2** covers descriptor *matching* and two-view geometric
//! *verification*:
//!
//! - [`matching`] — brute-force Hamming nearest-neighbour matching with the
//!   Lowe ratio test and a mutual cross-check ([`match_descriptors`],
//!   [`Match`]).
//! - [`geometry`] — RANSAC estimation of the **fundamental matrix** via the
//!   Hartley-normalized 8-point algorithm, scoring inliers by Sampson
//!   distance ([`verify_two_view`], [`RansacParams`], [`TwoViewResult`]).
//!   This stage recovers `F` and the inlier set **only**.
//!
//! **Stage 3 (added here)** is the [`twoview`] module — the **calibrated**
//! two-view geometry step. Given the camera [`CameraIntrinsics`] it upgrades
//! the fundamental matrix to the **essential matrix**
//! ([`essential_from_fundamental`]), decomposes it into the four candidate
//! relative poses ([`decompose_essential`]), triangulates correspondences by
//! the linear DLT ([`triangulate_point`]), and selects the physically correct
//! pose `(R, t)` by the cheirality test, returning a
//! [`TwoViewReconstruction`] ([`recover_pose`]). Two honest caveats stated in
//! the module docs: it is the **calibrated path** (intrinsics must be
//! supplied — no auto-calibration), and the recovered translation (and hence
//! the triangulated points) is fixed only **up to a global scale**.
//!
//! **Stage 4 (added here)** is the [`pnp`] module — the **incremental
//! mapper's** core operation: registering a *new* view into an existing
//! reconstruction by **camera resectioning** from 2D–3D correspondences
//! (already-triangulated scene points matched to the new image's features).
//! It solves the **Perspective-n-Point (PnP)** problem with the linear
//! **DLT** ([`solve_pnp`]): normalize pixels to calibrated rays, build the
//! `2n×12` homogeneous system from `x̂ᵢ × ([R|t]Xᵢ) = 0`, SVD it, reshape to a
//! `3×4` `[R|t]`, then orthonormalize the rotation block and rescale the
//! translation to stay metrically consistent. [`solve_pnp_ransac`] wraps that
//! in RANSAC over **reprojection error** for robustness to mismatches,
//! returning a [`PnpResult`] (pose + inlier set). Because the 3-D points come
//! in a fixed world frame, the recovered [`CameraPose`] is at **metric
//! scale** — there is no free scale here (contrast Stage 3). Honest caveats in
//! the module docs: the linear DLT is a strong *initializer* (below EPnP /
//! iterative accuracy on noisy data; bundle adjustment polishes it later), it
//! needs `n ≥ 6` points, and a **coplanar** point configuration is degenerate.
//!
//! ## Provenance / licensing
//!
//! This is a *clean-room* implementation written from the published
//! computer-vision literature — FAST (Rosten & Drummond), ORB (Rublee et
//! al.), the Lowe ratio test (Lowe 2004), the Hartley normalized 8-point
//! algorithm (Hartley 1997), RANSAC (Fischler & Bolles 1981), and DLT camera
//! resectioning / PnP (Hartley & Zisserman ch. 7) are all textbook
//! algorithms. **No COLMAP (or OpenCV) source code is used or
//! copied.** COLMAP is BSD-3-Clause and is credited purely as the *method
//! reference* for the overall SfM design in the crate's
//! `THIRD-PARTY-NOTICES` file. valenx itself is `MIT OR Apache-2.0`.
//!
//! ## What's implemented
//!
//! - [`GrayImage`] — a validated, row-major 8-bit grayscale buffer (no
//!   dependency on the `image` crate).
//! - [`Keypoint`] — a detected feature point `(x, y, score, angle)`.
//! - [`detect_fast`] — FAST-9 corner detection (16-px Bresenham circle,
//!   contiguous-arc ≥ 9 test, score-based non-maximum suppression).
//! - [`describe_keypoint`] / [`detect_and_describe`] — intensity-centroid
//!   orientation plus a 256-bit **steered-BRIEF** descriptor (`[u8; 32]`).
//!
//! ## Honesty note on the descriptor
//!
//! The descriptor is **oriented-FAST + steered-BRIEF**, not the *full*
//! learned `rBRIEF` of ORB: it steers a deterministic Gaussian test
//! pattern by the patch orientation but does **not** perform ORB's learned
//! test-selection / decorrelation pass. It is therefore a faithful steered
//! BRIEF and a good ORB approximation, but is neither bit-compatible with
//! OpenCV's `ORB` nor quite as discriminative as true rBRIEF. See the
//! [`descriptor`] module docs for the full explanation.
//!
//! ## Scope / limitations
//!
//! - Research/educational grade. Not a calibrated, survey-grade
//!   metrology product.
//! - Single-scale: there is no image pyramid yet, so detection is not
//!   scale-invariant. (A pyramid is a natural Stage-1.5 addition.)
//! - Operates on grayscale only; colour-to-gray conversion and image
//!   decoding are out of scope for this dependency-light stage.

#![forbid(unsafe_code)]

pub mod descriptor;
pub mod fast;
pub mod geometry;
pub mod image;
pub mod matching;
pub mod pnp;
pub mod twoview;

mod error;
mod keypoint;

pub use descriptor::{
    describe_keypoint, hamming_distance, DESCRIPTOR_BITS, DESCRIPTOR_BYTES, PATCH_RADIUS,
};
pub use error::{PhotogrammetryError, Result};
pub use fast::detect_fast;
pub use geometry::{verify_two_view, RansacParams, TwoViewResult};
pub use image::GrayImage;
pub use keypoint::Keypoint;
pub use matching::{match_descriptors, Match};
pub use pnp::{project_point, solve_pnp, solve_pnp_ransac, CameraPose, PnpRansacParams, PnpResult};
pub use twoview::{
    decompose_essential, essential_from_fundamental, recover_pose, triangulate_point,
    CameraIntrinsics, TwoViewReconstruction,
};

/// Detect FAST-9 corners in `img` and compute an ORB-style oriented
/// 256-bit (`[u8; 32]`) descriptor for each, returning the oriented
/// keypoints paired with their descriptors.
///
/// `threshold` is the FAST brightness margin `t` (see [`detect_fast`]);
/// larger values yield fewer, stronger corners. Each returned [`Keypoint`]
/// carries its FAST score and its intensity-centroid orientation; the
/// companion `[u8; 32]` is the steered-BRIEF signature (see the
/// [`descriptor`] module — this is steered BRIEF, **not** the full learned
/// rBRIEF; an honest distinction documented there).
///
/// This is the one-call Stage-1 entry point: `image in, (keypoint,
/// descriptor) pairs out`, ready to feed Stage-2 descriptor matching.
#[must_use]
pub fn detect_and_describe(img: &GrayImage, threshold: u8) -> Vec<(Keypoint, [u8; 32])> {
    detect_fast(img, threshold)
        .iter()
        .map(|kp| describe_keypoint(img, kp))
        .collect()
}
