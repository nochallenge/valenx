//! FAST-9 corner detection.
//!
//! Clean-room implementation of the FAST ("Features from Accelerated
//! Segment Test") corner detector of Rosten & Drummond, in its FAST-9
//! variant: a pixel is a corner when at least **9 contiguous** pixels on
//! the radius-3, 16-pixel Bresenham circle around it are all *brighter*
//! than the centre by a threshold `t`, or all *darker* by `t`.
//!
//! The detector here follows the textbook algorithm directly — no source
//! is taken from any existing library. Steps:
//!
//! 1. **High-speed rejection.** Examine the four compass pixels (circle
//!    indices 0, 4, 8, 12). A 9-arc can only exist if at least three of
//!    these four are uniformly brighter (or uniformly darker) than the
//!    centre ± `t`. This rejects the overwhelming majority of pixels with
//!    four comparisons.
//! 2. **Full segment test.** For survivors, classify all 16 circle pixels
//!    as brighter / darker / similar and look for a contiguous run of ≥ 9
//!    of one polarity, treating the circle as cyclic.
//! 3. **Score.** A corner's strength is the largest threshold at which it
//!    still passes the segment test, computed as the best (over the two
//!    polarities and all valid arcs) minimum absolute centre-to-arc
//!    intensity difference. This is the standard FAST score used to drive
//!    non-maximum suppression.
//! 4. **Non-maximum suppression.** A corner is kept only if its score is
//!    `>=` that of every corner in its 8-pixel neighbourhood, so dense
//!    clusters collapse to their local peak.
//!
//! Pixels within the circle radius (3 px) of the image border are skipped
//! so the 16 circle taps are always in bounds — no out-of-range access is
//! possible.

use crate::image::GrayImage;
use crate::keypoint::Keypoint;

/// Radius of the Bresenham sampling circle, in pixels. Pixels closer than
/// this to any image border are never detected as corners (their circle
/// would leave the image), so this is also the detector's border margin.
pub const FAST_RADIUS: usize = 3;

/// Number of points on the sampling circle.
pub(crate) const CIRCLE_LEN: usize = 16;

/// Minimum contiguous-arc length for the FAST-9 corner test.
const FAST_N: usize = 9;

/// The 16 `(dx, dy)` offsets of the radius-3 Bresenham circle, ordered
/// clockwise starting at the topmost point (12 o'clock). Index 0 is
/// straight up `(0, -3)`; indices advance clockwise. The four compass
/// points used by the high-speed test are at indices 0, 4, 8, 12.
pub(crate) const CIRCLE: [(i32, i32); CIRCLE_LEN] = [
    (0, -3),
    (1, -3),
    (2, -2),
    (3, -1),
    (3, 0),
    (3, 1),
    (2, 2),
    (1, 3),
    (0, 3),
    (-1, 3),
    (-2, 2),
    (-3, 1),
    (-3, 0),
    (-3, -1),
    (-2, -2),
    (-1, -3),
];

/// Classification of one circle pixel relative to the centre.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Cls {
    /// Circle pixel is brighter than `centre + t`.
    Brighter,
    /// Circle pixel is darker than `centre - t`.
    Darker,
    /// Within `[centre - t, centre + t]` — neither.
    Similar,
}

/// Detect FAST-9 corners in `img` with brightness threshold `threshold`,
/// returning the surviving keypoints after non-maximum suppression.
///
/// `threshold` is the per-pixel intensity margin `t`: a circle pixel
/// counts as brighter when its value exceeds `centre + t` and darker when
/// it is below `centre - t`. Larger values yield fewer, stronger corners.
///
/// Border pixels closer than [`FAST_RADIUS`] to any edge are never
/// reported (their circle would leave the image), so very small images
/// simply yield no corners rather than panicking. The returned keypoints
/// carry their FAST score in [`Keypoint::score`]; orientation is left at
/// `0.0` for the descriptor stage to fill in.
#[must_use]
pub fn detect_fast(img: &GrayImage, threshold: u8) -> Vec<Keypoint> {
    let w = img.width;
    let h = img.height;

    // No interior exists if the image is too small for a full circle.
    if w <= 2 * FAST_RADIUS || h <= 2 * FAST_RADIUS {
        return Vec::new();
    }

    // First pass: score every interior pixel (0 = not a corner). A dense
    // score grid over the interior makes the non-max-suppression neighbour
    // lookups O(1).
    let mut scores = vec![0u16; w * h];
    let t = i32::from(threshold);

    for y in FAST_RADIUS..(h - FAST_RADIUS) {
        for x in FAST_RADIUS..(w - FAST_RADIUS) {
            if let Some(score) = corner_score(img, x, y, t) {
                scores[y * w + x] = score;
            }
        }
    }

    // Second pass: non-maximum suppression over the 8-neighbourhood. A
    // candidate survives unless a neighbour scores strictly higher; ties on
    // an equal-score plateau are broken deterministically by keeping only
    // the pixel with the lowest (y, x) raster index, so a flat ridge of
    // equal scores collapses to a single representative rather than every
    // pixel on it.
    let mut keypoints = Vec::new();
    for y in FAST_RADIUS..(h - FAST_RADIUS) {
        for x in FAST_RADIUS..(w - FAST_RADIUS) {
            let s = scores[y * w + x];
            if s == 0 {
                continue;
            }
            if is_local_max(&scores, w, x, y, s) {
                keypoints.push(Keypoint::new(x as f32, y as f32, f32::from(s)));
            }
        }
    }

    keypoints
}

/// Return the FAST score at `(x, y)` if it is a corner, else `None`.
///
/// `(x, y)` must already be at least [`FAST_RADIUS`] from every border so
/// the 16 circle taps are in bounds (the caller guarantees this).
fn corner_score(img: &GrayImage, x: usize, y: usize, t: i32) -> Option<u16> {
    let centre = i32::from(img.at(x, y));

    // --- High-speed rejection on the four compass points. ---
    // The four points at indices 0/4/8/12 are spaced every 90° around the
    // circle, so any window of 9 contiguous indices (the FAST-9 arc length)
    // contains AT LEAST 2 of them. A corner therefore always has >= 2
    // compass points of a single polarity; if neither polarity reaches 2,
    // no 9-arc can exist and we reject with four comparisons.
    //
    // (Note: the familiar "3 of 4" high-speed test is only valid for the
    // FAST-12 variant, where a 12-arc must cover 3 compass points. Using 3
    // here would wrongly discard genuine FAST-9 corners — e.g. the apex of
    // a right-angle corner, which has exactly 2 bright + 2 dark compass
    // points.)
    let mut compass_bright = 0u8;
    let mut compass_dark = 0u8;
    for &idx in &[0usize, 4, 8, 12] {
        let (dx, dy) = CIRCLE[idx];
        let v = i32::from(sample(img, x, y, dx, dy));
        if v > centre + t {
            compass_bright += 1;
        } else if v < centre - t {
            compass_dark += 1;
        }
    }
    if compass_bright < 2 && compass_dark < 2 {
        return None;
    }

    // --- Full segment test. Classify all 16 circle pixels. ---
    let mut cls = [Cls::Similar; CIRCLE_LEN];
    let mut diff = [0i32; CIRCLE_LEN];
    for (i, &(dx, dy)) in CIRCLE.iter().enumerate() {
        let v = i32::from(sample(img, x, y, dx, dy));
        diff[i] = v - centre;
        cls[i] = if v > centre + t {
            Cls::Brighter
        } else if v < centre - t {
            Cls::Darker
        } else {
            Cls::Similar
        };
    }

    // Best contiguous arc of each polarity (cyclic). The score for an arc
    // is the minimum, over the arc, of the absolute centre difference — the
    // largest threshold that still keeps the whole arc on one side.
    let bright = best_arc_score(&cls, &diff, Cls::Brighter);
    let dark = best_arc_score(&cls, &diff, Cls::Darker);

    let score = bright.max(dark);
    if score > 0 {
        // Clamp to u16: diffs are at most 255, so this always fits.
        Some(score as u16)
    } else {
        None
    }
}

/// Length and strength of the best contiguous arc of class `target` on the
/// cyclic 16-pixel circle. Returns the arc *score* (min |diff| along the
/// longest qualifying arc) if some arc reaches length [`FAST_N`], else `0`.
///
/// Because the circle is cyclic, an arc may wrap from index 15 back to 0;
/// we walk `2 * CIRCLE_LEN` positions to cover every wrapped run exactly
/// once.
fn best_arc_score(cls: &[Cls; CIRCLE_LEN], diff: &[i32; CIRCLE_LEN], target: Cls) -> i32 {
    let mut best = 0i32;
    let mut run_len = 0usize;
    // Track the minimum absolute diff within the current run.
    let mut run_min = i32::MAX;

    for step in 0..(2 * CIRCLE_LEN) {
        let i = step % CIRCLE_LEN;
        if cls[i] == target {
            run_len += 1;
            let mag = diff[i].abs();
            if mag < run_min {
                run_min = mag;
            }
            // Any contiguous run of >= FAST_N qualifies. `run_min` only ever
            // decreases as a run grows, so a longer run can never raise the
            // score; that also makes the "whole circle is one class" case
            // (run_len up to 2*CIRCLE_LEN over the doubled walk) safe to
            // treat identically — no special-casing or double counting.
            if run_len >= FAST_N && run_min > best {
                best = run_min;
            }
        } else {
            run_len = 0;
            run_min = i32::MAX;
        }
    }
    best
}

/// Sample the pixel at `(x + dx, y + dy)`. The caller guarantees the
/// coordinate is in bounds (interior pixels only), so the signed offset
/// can be applied and cast back to `usize` safely.
#[inline]
fn sample(img: &GrayImage, x: usize, y: usize, dx: i32, dy: i32) -> u8 {
    let sx = (x as i32 + dx) as usize;
    let sy = (y as i32 + dy) as usize;
    img.at(sx, sy)
}

/// Is the score at `(x, y)` a local maximum over its 8-neighbourhood?
///
/// A candidate survives non-maximum suppression unless a neighbour has a
/// strictly greater score; ties are resolved deterministically by keeping
/// the candidate with the lowest `(y, x)` raster order (a later neighbour
/// of equal score does not suppress an earlier one).
fn is_local_max(scores: &[u16], w: usize, x: usize, y: usize, s: u16) -> bool {
    for ny in (y - 1)..=(y + 1) {
        for nx in (x - 1)..=(x + 1) {
            if nx == x && ny == y {
                continue;
            }
            let n = scores[ny * w + nx];
            if n > s {
                return false;
            }
            // Equal score: suppress the later (greater raster index) one.
            if n == s {
                let here = y * w + x;
                let there = ny * w + nx;
                if there < here {
                    return false;
                }
            }
        }
    }
    true
}
