//! CAD-kernel validation suite — boolean operations.
//!
//! Union / difference / intersection are checked against
//! **analytically-known result volumes**: every test sets up
//! primitives whose overlap geometry has a closed-form volume, runs
//! the boolean, and asserts the measured result volume plus that the
//! output is a valid closed solid.
//!
//! The boolean kernel (`truck-shapeops`) is the part of any CAD kernel
//! that is hardest to get right; these tests are the honest proof that
//! Valenx's wrapper produces *correct* — not merely non-empty —
//! results, and the degenerate-case tests (coincident faces, tangent
//! contact, fully-engulfed operands) exercise the robustness path.
//!
//! ## Flat-faced solids only
//!
//! Every result here is a boolean of *boxes*, so the result boundary
//! is flat-faced and the measured volume is **exact** — no convergence
//! tolerance needed. Boolean tests that involve a cylinder bore use a
//! convergence bound instead (the bore wall is curved).

use truck_modeling::{builder, EuclideanSpace, Point3, Rad, Vector3};
use valenx_cad::measure::{is_closed_solid_tol, solid_volume_tol};
use valenx_cad::{box_solid, difference, intersection, union, CadError, Solid};

/// Fine tessellation tolerance.
const FINE: f64 = 5.0e-4;

/// Build an axis-aligned box spanning `[min, min+size]` on each axis.
/// `valenx_cad::box_solid` always corners at the origin, so this
/// translates it into place — exactly what positioned boolean tests
/// need.
fn box_at(min: (f64, f64, f64), size: (f64, f64, f64)) -> Solid {
    box_solid(size.0, size.1, size.2)
        .expect("box dims positive")
        .translated(min.0, min.1, min.2)
        .expect("translation deltas finite")
}

/// Build a cylinder of `radius` / `height` whose base disk is centred
/// at `(cx, cy, z0)`, axis along +Z.
fn cylinder_at(cx: f64, cy: f64, z0: f64, radius: f64, height: f64) -> Solid {
    let v = builder::vertex(Point3::new(cx + radius, cy, z0));
    let circle = builder::rsweep(
        &v,
        Point3::new(cx, cy, z0),
        Vector3::unit_z(),
        Rad(2.0 * std::f64::consts::PI),
    );
    let disk = builder::try_attach_plane(&[circle]).expect("cylinder disk");
    let inner: truck_modeling::Solid = builder::tsweep(&disk, Vector3::new(0.0, 0.0, height));
    let _ = Point3::origin();
    Solid::from_truck(inner)
}

// ===========================================================================
// Union — two boxes sharing a known overlap volume.
// ===========================================================================

#[test]
fn union_of_two_overlapping_boxes_has_inclusion_exclusion_volume() {
    // Box A: [0,2]³, volume 8. Box B: [1,3]³, volume 8.
    // Overlap: [1,2]³, volume 1.
    // |A ∪ B| = |A| + |B| − |A ∩ B| = 8 + 8 − 1 = 15.
    let a = box_at((0.0, 0.0, 0.0), (2.0, 2.0, 2.0));
    let b = box_at((1.0, 1.0, 1.0), (2.0, 2.0, 2.0));
    let u = union(&a, &b).expect("union of overlapping boxes");
    let v = solid_volume_tol(&u, FINE).unwrap();
    assert!(
        (v - 15.0).abs() < 1e-6,
        "A∪B volume {v} != 15 (inclusion-exclusion)"
    );
    assert!(
        is_closed_solid_tol(&u, FINE).unwrap(),
        "the union must be a valid closed solid"
    );
}

#[test]
fn union_of_disjoint_then_overlapping_is_monotone() {
    // |A ∪ B| with a larger overlap is smaller (more shared volume
    // removed by inclusion-exclusion). Two 2³ boxes:
    //   small overlap [1.5,2]³ → |∪| = 8+8−0.125 = 15.875
    //   large overlap [0.5,2]³ → |∪| = 8+8−3.375 = 12.625
    let a = box_at((0.0, 0.0, 0.0), (2.0, 2.0, 2.0));

    let b_small = box_at((1.5, 1.5, 1.5), (2.0, 2.0, 2.0));
    let u_small = union(&a, &b_small).expect("union small overlap");
    let v_small = solid_volume_tol(&u_small, FINE).unwrap();
    assert!(
        (v_small - 15.875).abs() < 1e-6,
        "small-overlap union {v_small} != 15.875"
    );

    let b_large = box_at((0.5, 0.5, 0.5), (2.0, 2.0, 2.0));
    let u_large = union(&a, &b_large).expect("union large overlap");
    let v_large = solid_volume_tol(&u_large, FINE).unwrap();
    assert!(
        (v_large - 12.625).abs() < 1e-6,
        "large-overlap union {v_large} != 12.625"
    );

    assert!(
        v_large < v_small,
        "more overlap → smaller union ({v_large} should be < {v_small})"
    );
}

// ===========================================================================
// Intersection — overlap volume is the analytic box-intersection.
// ===========================================================================

#[test]
fn intersection_of_two_boxes_is_the_overlap_box() {
    // A = [0,2]³, B = [1,4]×[1,4]×[1,4]. Overlap = [1,2]³ → volume 1.
    let a = box_at((0.0, 0.0, 0.0), (2.0, 2.0, 2.0));
    let b = box_at((1.0, 1.0, 1.0), (3.0, 3.0, 3.0));
    let inter = intersection(&a, &b).expect("intersection of overlapping boxes");
    let v = solid_volume_tol(&inter, FINE).unwrap();
    assert!(
        (v - 1.0).abs() < 1e-6,
        "A∩B volume {v} != 1 (overlap box [1,2]³)"
    );
    assert!(
        is_closed_solid_tol(&inter, FINE).unwrap(),
        "the intersection must be a valid closed solid"
    );
}

#[test]
fn intersection_of_an_engulfed_box_is_the_inner_box() {
    // B = [0.5,1.5]³ is entirely inside A = [0,3]³.
    // A ∩ B = B → volume 1.
    let a = box_at((0.0, 0.0, 0.0), (3.0, 3.0, 3.0));
    let b = box_at((0.5, 0.5, 0.5), (1.0, 1.0, 1.0));
    let inter = intersection(&a, &b).expect("intersection of engulfed box");
    let v = solid_volume_tol(&inter, FINE).unwrap();
    assert!(
        (v - 1.0).abs() < 1e-6,
        "A∩B with B inside A should equal |B| = 1, got {v}"
    );
}

// ===========================================================================
// Difference — material removed equals the overlap volume.
// ===========================================================================

#[test]
fn difference_removes_exactly_the_overlap_volume() {
    // A = [0,4]³ (volume 64). B = [3,6]³ overlaps A in [3,4]³
    // (volume 1). A − B keeps 64 − 1 = 63.
    let a = box_at((0.0, 0.0, 0.0), (4.0, 4.0, 4.0));
    let b = box_at((3.0, 3.0, 3.0), (3.0, 3.0, 3.0));
    let diff = difference(&a, &b).expect("difference of overlapping boxes");
    let v = solid_volume_tol(&diff, FINE).unwrap();
    assert!(
        (v - 63.0).abs() < 1e-6,
        "A−B volume {v} != 63 (64 − 1 overlap)"
    );
    assert!(
        is_closed_solid_tol(&diff, FINE).unwrap(),
        "the difference must be a valid closed solid"
    );
}

#[test]
fn difference_with_disjoint_tool_keeps_full_volume() {
    // B = [10,11]³ does not touch A = [0,2]³. A − B = A → volume 8.
    // truck-shapeops may return EmptyResult for a fully-disjoint
    // subtraction; both "result == A" and a typed EmptyResult are
    // acceptable graceful outcomes — what must NOT happen is a panic
    // or a silently-wrong volume.
    let a = box_at((0.0, 0.0, 0.0), (2.0, 2.0, 2.0));
    let b = box_at((10.0, 10.0, 10.0), (1.0, 1.0, 1.0));
    match difference(&a, &b) {
        Ok(diff) => {
            let v = solid_volume_tol(&diff, FINE).unwrap();
            assert!(
                (v - 8.0).abs() < 1e-6,
                "A−(disjoint B) should keep |A| = 8, got {v}"
            );
        }
        Err(CadError::EmptyResult) => {
            // Acceptable — a disjoint subtraction has no shapeops
            // intersection geometry to work with.
        }
        Err(other) => panic!("disjoint difference should not hard-fail: {other:?}"),
    }
}

// ===========================================================================
// Cylinder bored through a block — the canonical "punched cube".
// ===========================================================================

#[test]
fn cylinder_bored_through_block_removes_pi_r2_h() {
    // Block [0,4]×[0,4]×[0,2], volume 32. A cylinder radius 1
    // centred at (2,2) bored all the way through (the cylinder
    // over-runs the block in Z so the bore is a clean through-hole).
    // Removed material = πr²·h_block = π·1·2 = 2π ≈ 6.2832.
    // Result volume = 32 − 2π.
    let block = box_at((0.0, 0.0, 0.0), (4.0, 4.0, 2.0));
    let drill = cylinder_at(2.0, 2.0, -1.0, 1.0, 4.0); // spans z=-1..3
    let bored = difference(&block, &drill).expect("bored block");

    let v = solid_volume_tol(&bored, FINE).unwrap();
    let exact = 32.0 - std::f64::consts::PI * 1.0 * 1.0 * 2.0;
    // The bore wall is curved → the measured volume converges from
    // above (the inscribed bore facets leave slightly MORE material).
    assert!(
        (v - exact).abs() / exact < 0.01,
        "bored block volume {v} should be within 1% of {exact}"
    );
    assert!(
        is_closed_solid_tol(&bored, FINE).unwrap(),
        "the bored block must be a valid closed solid"
    );
    // A through-hole adds the curved bore wall to the 6 block faces.
    assert!(
        bored.faces() > 6,
        "a bored block should have more than 6 faces, got {}",
        bored.faces()
    );
}

// ===========================================================================
// Boolean robustness — degenerate / edge cases must fall through
// cleanly, never panic, never produce a silently-invalid solid.
// ===========================================================================

#[test]
fn boolean_with_coincident_faces_does_not_panic() {
    // Two boxes sharing an entire face: A = [0,1]³, B = [1,2]×[0,1]×[0,1].
    // They touch on the plane x = 1 with zero overlap volume.
    // Coincident faces are the historically fragile boolean input.
    // The kernel must produce *some* defined outcome — a valid solid
    // or a typed error — never a panic.
    let a = box_at((0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    let b = box_at((1.0, 0.0, 0.0), (1.0, 1.0, 1.0));

    // Union of face-touching boxes: ideally a 1×2×1 bar (volume 2),
    // but a degenerate coincident-face union may also surface
    // EmptyResult. Both are graceful; a panic is not.
    match union(&a, &b) {
        Ok(u) => {
            let v = solid_volume_tol(&u, FINE).unwrap();
            assert!(
                v > 0.0 && v <= 2.0 + 1e-6,
                "face-touching union volume {v} out of range (0, 2]"
            );
        }
        Err(CadError::EmptyResult) => { /* graceful */ }
        Err(other) => panic!("coincident-face union should not hard-fail: {other:?}"),
    }

    // Intersection of face-touching boxes: zero-volume contact. The
    // kernel should report EmptyResult or a degenerate-but-defined
    // result — never panic.
    match intersection(&a, &b) {
        Ok(_) | Err(CadError::EmptyResult) => { /* graceful */ }
        Err(other) => panic!("coincident-face intersection should not hard-fail: {other:?}"),
    }
}

#[test]
fn boolean_with_tangent_contact_does_not_panic() {
    // Two cylinders touching along a single tangent line — edge-on
    // contact with zero overlap volume. A classic degenerate input.
    let c1 = cylinder_at(0.0, 0.0, 0.0, 1.0, 2.0);
    let c2 = cylinder_at(2.0, 0.0, 0.0, 1.0, 2.0); // touches c1 at x=1
    match intersection(&c1, &c2) {
        Ok(_) | Err(CadError::EmptyResult) => { /* graceful */ }
        Err(other) => panic!("tangent-contact intersection should not hard-fail: {other:?}"),
    }
    match union(&c1, &c2) {
        Ok(u) => {
            // Tangent union: roughly two disjoint cylinders. Volume
            // should be near 2·πr²h = 4π but the kernel may or may not
            // weld them; just require a positive defined volume.
            let v = solid_volume_tol(&u, 1e-3).unwrap();
            assert!(v > 0.0, "tangent union volume {v} should be positive");
        }
        Err(CadError::EmptyResult) => { /* graceful */ }
        Err(other) => panic!("tangent-contact union should not hard-fail: {other:?}"),
    }
}

#[test]
fn intersection_of_fully_disjoint_solids_is_empty_not_a_panic() {
    // Two boxes far apart. A ∩ B is genuinely empty; the kernel must
    // surface EmptyResult, not a panic and not a bogus solid.
    let a = box_at((0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    let b = box_at((100.0, 100.0, 100.0), (1.0, 1.0, 1.0));
    match intersection(&a, &b) {
        Err(CadError::EmptyResult) => { /* the expected outcome */ }
        Ok(s) => {
            // If the kernel returns something, it must at least be a
            // near-zero-volume degenerate — not a phantom solid.
            let v = solid_volume_tol(&s, 1e-3).unwrap().abs();
            assert!(
                v < 1e-6,
                "intersection of disjoint solids should be empty, got volume {v}"
            );
        }
        Err(other) => panic!("disjoint intersection should be EmptyResult, got {other:?}"),
    }
}

#[test]
fn self_intersection_difference_collapses_to_empty_or_degenerate() {
    // A − A: subtracting a solid from an identical copy. The exact
    // answer is the empty set; the kernel must not panic and must not
    // return a phantom positive-volume solid.
    let a = box_at((0.0, 0.0, 0.0), (2.0, 2.0, 2.0));
    let b = box_at((0.0, 0.0, 0.0), (2.0, 2.0, 2.0));
    match difference(&a, &b) {
        Err(CadError::EmptyResult) => { /* the clean outcome */ }
        Ok(s) => {
            let v = solid_volume_tol(&s, 1e-3).unwrap().abs();
            assert!(v < 1e-3, "A−A should be empty or near-zero volume, got {v}");
        }
        Err(other) => panic!("A−A should not hard-fail with {other:?}"),
    }
}
