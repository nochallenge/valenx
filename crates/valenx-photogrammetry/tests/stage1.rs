//! Stage-1 correctness gate for `valenx-photogrammetry`.
//!
//! Exercises the public API end-to-end against synthetic images with
//! known ground truth:
//!
//! 1. FAST finds corners near the four corners of a white square on black,
//!    and none deep in the flat interior.
//! 2. A uniform (flat) image yields zero corners.
//! 3. A pixel-count / dimension mismatch is rejected by the constructor.
//! 4. The descriptor is deterministic, with ~0 Hamming distance for the
//!    same patch and a large distance for a clearly different patch.
//! 5. No panic on tiny images or keypoints near the border.

use valenx_photogrammetry::{
    describe_keypoint, detect_and_describe, detect_fast, hamming_distance, GrayImage, Keypoint,
    PhotogrammetryError, DESCRIPTOR_BITS,
};

/// Build a `w x h` image filled with `bg`, then paint an axis-aligned solid
/// square of intensity `fg` whose top-left corner is at `(sx, sy)` and
/// whose side length is `side`.
fn white_square_on_black(
    w: usize,
    h: usize,
    sx: usize,
    sy: usize,
    side: usize,
    bg: u8,
    fg: u8,
) -> GrayImage {
    let mut px = vec![bg; w * h];
    for y in sy..(sy + side) {
        for x in sx..(sx + side) {
            px[y * w + x] = fg;
        }
    }
    GrayImage::new(w, h, px).expect("valid synthetic image")
}

/// Smallest distance from `(px, py)` to any keypoint, or `f32::INFINITY`
/// if there are none.
fn nearest_kp_dist(kps: &[Keypoint], px: f32, py: f32) -> f32 {
    kps.iter()
        .map(|k| ((k.x - px).powi(2) + (k.y - py).powi(2)).sqrt())
        .fold(f32::INFINITY, f32::min)
}

// ---------------------------------------------------------------------------
// Test 1 — corners of a white square are found; flat interior is not.
// ---------------------------------------------------------------------------

#[test]
fn finds_corners_of_white_square_not_interior() {
    // A generous 40x40 square so each corner is a clean convex corner and
    // the interior has a large flat region.
    let w = 80;
    let h = 80;
    let sx = 20;
    let sy = 20;
    let side = 40;
    let img = white_square_on_black(w, h, sx, sy, side, 0, 255);

    let kps = detect_fast(&img, 30);
    assert!(!kps.is_empty(), "expected some corners on the square");

    // The four geometric corners (top-left, top-right, bottom-left,
    // bottom-right) of the filled square.
    let corners = [
        (sx as f32, sy as f32),
        ((sx + side - 1) as f32, sy as f32),
        (sx as f32, (sy + side - 1) as f32),
        ((sx + side - 1) as f32, (sy + side - 1) as f32),
    ];
    // FAST may fire a pixel or two off the exact apex; allow a small radius.
    let tol = 3.0;
    for (cx, cy) in corners {
        let d = nearest_kp_dist(&kps, cx, cy);
        assert!(
            d <= tol,
            "no keypoint within {tol} px of square corner ({cx}, {cy}); nearest was {d}"
        );
    }

    // No keypoint should land deep inside the flat (uniform) interior.
    // Use a margin well inside the FAST border so edge effects can't count.
    let in_margin = 6.0;
    let lo_x = sx as f32 + in_margin;
    let hi_x = (sx + side) as f32 - in_margin;
    let lo_y = sy as f32 + in_margin;
    let hi_y = (sy + side) as f32 - in_margin;
    for k in &kps {
        let inside = k.x > lo_x && k.x < hi_x && k.y > lo_y && k.y < hi_y;
        assert!(
            !inside,
            "unexpected corner in the flat interior at ({}, {})",
            k.x, k.y
        );
    }
}

// ---------------------------------------------------------------------------
// Test 2 — a uniform image yields zero corners.
// ---------------------------------------------------------------------------

#[test]
fn uniform_image_has_no_corners() {
    // Mid-grey, a black, and a white flat field — none should produce a
    // corner regardless of the constant level.
    for level in [0u8, 128, 255] {
        let img = GrayImage::new(64, 48, vec![level; 64 * 48]).unwrap();
        let kps = detect_fast(&img, 20);
        assert!(
            kps.is_empty(),
            "flat field at level {level} produced {} corners",
            kps.len()
        );
    }
}

// ---------------------------------------------------------------------------
// Test 3 — pixel-count / dimension validation.
// ---------------------------------------------------------------------------

#[test]
fn constructor_rejects_bad_inputs() {
    // Wrong buffer length.
    let err = GrayImage::new(4, 4, vec![0u8; 15]).unwrap_err();
    assert!(
        matches!(
            err,
            PhotogrammetryError::PixelCountMismatch {
                expected: 16,
                actual: 15,
                ..
            }
        ),
        "expected PixelCountMismatch, got {err:?}"
    );

    // Too many pixels also fails.
    assert!(matches!(
        GrayImage::new(2, 2, vec![0u8; 5]).unwrap_err(),
        PhotogrammetryError::PixelCountMismatch {
            expected: 4,
            actual: 5,
            ..
        }
    ));

    // Zero dimensions fail.
    assert!(matches!(
        GrayImage::new(0, 10, vec![]).unwrap_err(),
        PhotogrammetryError::ZeroDimension { .. }
    ));
    assert!(matches!(
        GrayImage::new(10, 0, vec![]).unwrap_err(),
        PhotogrammetryError::ZeroDimension { .. }
    ));

    // A correctly sized buffer succeeds.
    let ok = GrayImage::new(4, 4, vec![7u8; 16]).unwrap();
    assert_eq!(ok.len(), 16);
    assert_eq!(ok.at(1, 1), 7);
    assert_eq!(ok.get(4, 0), None);
}

// ---------------------------------------------------------------------------
// Test 4 — descriptor determinism + discrimination.
// ---------------------------------------------------------------------------

/// A textured patch so the centroid orientation and BRIEF tests are well
/// defined (a flat patch would give a degenerate, all-zero-ish descriptor).
fn textured_image(w: usize, h: usize, seed: u32) -> GrayImage {
    let mut px = vec![0u8; w * h];
    for y in 0..h {
        for x in 0..w {
            // Deterministic pseudo-texture: a couple of mixed frequencies.
            let v = ((x as u32)
                .wrapping_mul(37)
                .wrapping_add((y as u32).wrapping_mul(101))
                ^ seed.wrapping_mul(2654435761))
            .wrapping_add((x as u32 * y as u32).wrapping_mul(13));
            px[y * w + x] = (v % 256) as u8;
        }
    }
    GrayImage::new(w, h, px).unwrap()
}

#[test]
fn descriptor_is_deterministic_and_discriminative() {
    let img_a = textured_image(64, 64, 1);
    let kp = Keypoint::new(32.0, 32.0, 1.0);

    // Determinism: same image + same keypoint → identical descriptor and
    // identical orientation, every time.
    let (k1, d1) = describe_keypoint(&img_a, &kp);
    let (k2, d2) = describe_keypoint(&img_a, &kp);
    assert_eq!(d1, d2, "descriptor must be deterministic");
    assert_eq!(
        k1.angle.to_bits(),
        k2.angle.to_bits(),
        "angle must be deterministic"
    );
    assert_eq!(
        hamming_distance(&d1, &d2),
        0,
        "identical patch must have zero Hamming distance"
    );

    // Descriptor must carry real information (not all-zero / all-one).
    let ones: u32 = d1.iter().map(|b| b.count_ones()).sum();
    assert!(
        ones > 16 && ones < (DESCRIPTOR_BITS as u32 - 16),
        "descriptor looks degenerate: {ones} set bits of {DESCRIPTOR_BITS}"
    );

    // Discrimination: a clearly different patch (different texture seed)
    // must be far in Hamming space.
    let img_b = textured_image(64, 64, 999);
    let (_, d_diff) = describe_keypoint(&img_b, &kp);
    let far = hamming_distance(&d1, &d_diff);
    assert!(
        far > 40,
        "different patches should be far apart in Hamming space, got {far} of {DESCRIPTOR_BITS}"
    );

    // And the same-vs-same distance is strictly smaller than same-vs-diff.
    assert!(hamming_distance(&d1, &d2) < far);
}

// ---------------------------------------------------------------------------
// Test 5 — no panic on tiny images / near borders.
// ---------------------------------------------------------------------------

#[test]
fn tiny_images_and_borders_do_not_panic() {
    // Images smaller than / around the FAST circle diameter (7 px): the
    // detector must return cleanly (empty) rather than index out of bounds.
    for (w, h) in [(1usize, 1usize), (3, 3), (6, 6), (7, 7), (7, 9), (9, 7)] {
        let img = GrayImage::new(w, h, vec![123u8; w * h]).unwrap();
        let kps = detect_fast(&img, 10);
        // No assertion on count — just that it ran without panicking. Tiny
        // flat images legitimately have no corners.
        let _ = kps.len();

        // The full pipeline must also be panic-free on tiny images.
        let pairs = detect_and_describe(&img, 10);
        let _ = pairs.len();
    }

    // Describing a keypoint sitting right at a corner of a small image must
    // not panic, even though its 31x31 BRIEF patch reaches far outside the
    // image (sampling is edge-clamped).
    let small = textured_image(10, 10, 5);
    for &(x, y) in &[
        (0.0f32, 0.0f32),
        (9.0, 0.0),
        (0.0, 9.0),
        (9.0, 9.0),
        (1.0, 1.0),
    ] {
        let kp = Keypoint::new(x, y, 1.0);
        let (_oriented, _desc) = describe_keypoint(&small, &kp);
    }
}

// ---------------------------------------------------------------------------
// Bonus — the one-call entry point produces 32-byte descriptors.
// ---------------------------------------------------------------------------

#[test]
fn detect_and_describe_shape() {
    let img = white_square_on_black(80, 80, 20, 20, 40, 0, 255);
    let results = detect_and_describe(&img, 30);
    assert!(!results.is_empty());
    for (kp, desc) in &results {
        // Descriptor is exactly 32 bytes = 256 bits.
        assert_eq!(desc.len(), 32);
        // Orientation is a finite angle in (-pi, pi].
        assert!(kp.angle.is_finite());
        assert!(kp.angle > -std::f32::consts::PI - 1e-3 && kp.angle <= std::f32::consts::PI + 1e-3);
        // Score is positive (it passed the FAST test).
        assert!(kp.score > 0.0);
    }
}
