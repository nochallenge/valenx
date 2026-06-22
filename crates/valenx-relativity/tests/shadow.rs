//! Ground-truth validation of the shadow ray-tracer: the Schwarzschild shadow
//! is a centred disk of radius 3√3 M; the Kerr shadow is asymmetric.

use valenx_relativity::{equatorial_shadow_edges, kerr, render_shadow, schwarzschild};

/// The Schwarzschild shadow is symmetric with edges at ∓3√3 M.
#[test]
fn schwarzschild_shadow_is_symmetric_sqrt27() {
    let b_crit = 27.0_f64.sqrt(); // 3√3 ≈ 5.196
    let (left, right) = equatorial_shadow_edges(&schwarzschild(1.0), 10.0).unwrap();
    assert!(
        (right - b_crit).abs() / b_crit < 0.01,
        "right edge {right} vs {b_crit}"
    );
    assert!(
        (-left - b_crit).abs() / b_crit < 0.01,
        "left edge {left} vs {}",
        -b_crit
    );
    assert!(
        (right + left).abs() < 0.05,
        "should be symmetric, offset = {}",
        right + left
    );
}

/// The Kerr shadow is shifted/asymmetric: the two edges are not mirror images,
/// because frame dragging moves the apparent rim toward the co-rotating side.
#[test]
fn kerr_shadow_is_asymmetric() {
    let (left, right) = equatorial_shadow_edges(&kerr(1.0, 0.9), 10.0).unwrap();
    assert!(
        left < 0.0 && right > 0.0,
        "edges bracket the origin: {left}, {right}"
    );
    let offset = right + left; // zero if symmetric
    assert!(
        offset.abs() > 0.3,
        "Kerr shadow should be visibly shifted, offset = {offset}"
    );
    // Still a sensible overall size.
    let width = right - left;
    assert!(
        (6.0..=12.0).contains(&width),
        "shadow width {width} out of range"
    );
}

/// A rendered image has a dark disk in the middle and sky at the corners.
#[test]
fn rendered_image_has_shadow_disk() {
    let img = render_shadow(
        &schwarzschild(1.0),
        200.0,
        std::f64::consts::FRAC_PI_2,
        8.0,
        21,
        21,
    )
    .unwrap();
    assert_eq!(img.width, 21);
    assert_eq!(img.height, 21);
    // Centre pixel is in shadow; corners see the sky.
    assert!(img.is_shadow(10, 10), "centre should be shadow");
    assert!(!img.is_shadow(0, 0), "corner should be sky");
    assert!(!img.is_shadow(20, 20), "corner should be sky");
    // The shadow disk (radius ≈5.2 in a 16-wide plane) covers a sensible area.
    let frac = img.shadow_fraction();
    assert!(
        (0.15..0.5).contains(&frac),
        "shadow fraction {frac} out of range"
    );
}
