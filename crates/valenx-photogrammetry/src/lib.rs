//! # valenx-photogrammetry — structure-from-motion building blocks (Stage 1)
//!
//! An **in-house, clean-room Rust** reimplementation of the front end of a
//! classical Structure-from-Motion (SfM) pipeline, of the kind made famous
//! by COLMAP. The full pipeline (built out over later stages) is:
//!
//! ```text
//! images → feature detection & description   ← THIS STAGE
//!        → descriptor matching
//!        → two-view geometry (essential/fundamental matrix, RANSAC)
//!        → incremental mapper (PnP, triangulation)
//!        → bundle adjustment
//!        → sparse point cloud + camera poses
//! ```
//!
//! **Stage 1 (this crate as shipped) covers the first arrow only:** the
//! grayscale [`GrayImage`] input type, FAST-9 corner detection (the
//! [`fast`] module), and an ORB-style oriented binary descriptor (the
//! [`descriptor`] module).
//!
//! ## Provenance / licensing
//!
//! This is a *clean-room* implementation written from the published
//! computer-vision literature — FAST (Rosten & Drummond) and ORB (Rublee
//! et al.) are textbook algorithms. **No COLMAP (or OpenCV) source code is
//! used or copied.** COLMAP is BSD-3-Clause and is credited purely as the
//! *method reference* for the overall SfM design in the crate's
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
pub mod image;

mod error;
mod keypoint;

pub use descriptor::{
    describe_keypoint, hamming_distance, DESCRIPTOR_BITS, DESCRIPTOR_BYTES, PATCH_RADIUS,
};
pub use error::{PhotogrammetryError, Result};
pub use fast::detect_fast;
pub use image::GrayImage;
pub use keypoint::Keypoint;

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
