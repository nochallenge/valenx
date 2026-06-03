//! CAD-kernel validation suite — primitives.
//!
//! Each primitive (box, cylinder, sphere, cone, torus, prism) is built
//! and checked against its **exact analytic ground truth**:
//!
//! - **Volume** — measured by the divergence-theorem integral and
//!   compared to the closed-form formula.
//! - **Surface area** — measured by summing boundary triangles and
//!   compared to the closed-form formula.
//! - **Closed-solid validity** — the boundary must be a closed
//!   2-manifold.
//! - **Topology** — the Euler characteristic `V−E+F` must match the
//!   genus the primitive's shape demands (2 for a genus-0 solid, 0 for
//!   the genus-1 torus).
//!
//! This is the core of the commercial-depth correctness pass: a CAD
//! kernel that cannot prove its primitives is not trustworthy. Curved
//! primitives are checked with a *convergence* assertion — the
//! tessellated measurement undershoots a convex curved boundary, so
//! the test asserts the measured value approaches the exact value from
//! below within a stated tolerance, never that it equals it.

use std::f64::consts::PI;

use valenx_cad::measure::{
    euler_characteristic, is_closed_solid_tol, solid_area_tol, solid_volume_tol,
};
use valenx_cad::{box_solid, cone, cylinder, prism, sphere, torus};

/// Fine tessellation tolerance for curved-solid convergence checks.
const FINE: f64 = 5.0e-4;

// ===========================================================================
// Box — flat-faced, every measurement is exact.
// ===========================================================================

#[test]
fn box_volume_area_topology_are_exact() {
    let (dx, dy, dz) = (2.0, 3.0, 5.0);
    let b = box_solid(dx, dy, dz).unwrap();

    // Volume = dx·dy·dz, exact for a flat-faced solid.
    let v = solid_volume_tol(&b, FINE).unwrap();
    assert!(
        (v - dx * dy * dz).abs() < 1e-9,
        "box volume {v} != {}",
        dx * dy * dz
    );

    // Surface area = 2(dx·dy + dy·dz + dz·dx), exact.
    let a = solid_area_tol(&b, FINE).unwrap();
    let exact_a = 2.0 * (dx * dy + dy * dz + dz * dx);
    assert!((a - exact_a).abs() < 1e-9, "box area {a} != {exact_a}");

    // Topology: a box is genus-0 → V−E+F = 8−12+6 = 2.
    assert_eq!(b.vertices(), 8, "box has 8 vertices");
    assert_eq!(b.edges(), 12, "box has 12 edges");
    assert_eq!(b.faces(), 6, "box has 6 faces");
    assert_eq!(euler_characteristic(&b), Some(2), "box Euler χ must be 2");

    assert!(
        is_closed_solid_tol(&b, FINE).unwrap(),
        "a box must be a closed 2-manifold"
    );
}

// ===========================================================================
// Cylinder — πr²h volume, 2πr(r+h) area. Curved → converges from below.
// ===========================================================================

#[test]
fn cylinder_volume_converges_to_analytic() {
    let (r, h) = (1.5, 4.0);
    let cyl = cylinder(r, h).unwrap();
    let exact = PI * r * r * h;
    let v = solid_volume_tol(&cyl, FINE).unwrap();
    assert!(v > 0.0, "cylinder volume must be positive, got {v}");
    // An inscribed-facet tessellation undershoots a convex curved
    // boundary — the measured volume must not exceed the exact value
    // (a small tessellation-curvature slop is allowed).
    assert!(
        v <= exact * (1.0 + 1e-4),
        "cylinder volume {v} should not exceed exact {exact}"
    );
    assert!(
        (exact - v) / exact < 0.01,
        "cylinder volume {v} should converge within 1% of {exact} (rel err {})",
        (exact - v) / exact
    );
}

#[test]
fn cylinder_surface_area_converges_to_analytic() {
    let (r, h) = (1.5, 4.0);
    let cyl = cylinder(r, h).unwrap();
    // Total area = 2 caps (2·πr²) + side (2πr·h).
    let exact = 2.0 * PI * r * r + 2.0 * PI * r * h;
    let a = solid_area_tol(&cyl, FINE).unwrap();
    assert!(
        a <= exact * (1.0 + 1e-4),
        "cylinder area {a} should not exceed exact {exact}"
    );
    assert!(
        (exact - a) / exact < 0.01,
        "cylinder area {a} should converge within 1% of {exact}"
    );
}

#[test]
fn cylinder_is_a_closed_solid() {
    let cyl = cylinder(1.0, 2.0).unwrap();
    assert!(
        is_closed_solid_tol(&cyl, FINE).unwrap(),
        "a cylinder must be a closed 2-manifold"
    );
}

#[test]
fn cylinder_euler_characteristic_is_genus_zero() {
    // A cylinder is topologically a sphere (genus 0) → χ = 2,
    // regardless of how truck splits its closed side sweep.
    let cyl = cylinder(1.0, 2.0).unwrap();
    assert_eq!(
        euler_characteristic(&cyl),
        Some(2),
        "a cylinder is genus-0; V−E+F must be 2"
    );
}

// ===========================================================================
// Sphere — 4/3·πr³ volume, 4πr² area. Curved → converges from below.
// ===========================================================================

#[test]
fn sphere_volume_converges_to_analytic() {
    let r = 2.0;
    let s = sphere(r).unwrap();
    let exact = 4.0 / 3.0 * PI * r * r * r;
    let v = solid_volume_tol(&s, FINE).unwrap();
    assert!(v > 0.0, "sphere volume must be positive, got {v}");
    assert!(
        v <= exact * (1.0 + 1e-4),
        "sphere volume {v} should not exceed exact {exact}"
    );
    assert!(
        (exact - v) / exact < 0.01,
        "sphere volume {v} should converge within 1% of {exact} (rel err {})",
        (exact - v) / exact
    );
}

#[test]
fn sphere_surface_area_converges_to_analytic() {
    let r = 2.0;
    let s = sphere(r).unwrap();
    let exact = 4.0 * PI * r * r;
    let a = solid_area_tol(&s, FINE).unwrap();
    assert!(
        a <= exact * (1.0 + 1e-4),
        "sphere area {a} should not exceed exact {exact}"
    );
    assert!(
        (exact - a) / exact < 0.01,
        "sphere area {a} should converge within 1% of {exact}"
    );
}

#[test]
fn sphere_is_a_closed_solid() {
    let s = sphere(1.0).unwrap();
    assert!(
        is_closed_solid_tol(&s, FINE).unwrap(),
        "a sphere must be a closed 2-manifold"
    );
}

// ===========================================================================
// Cone — ⅓πr²h volume (pointed); frustum ⅓πh(R²+Rr+r²).
// ===========================================================================

#[test]
fn pointed_cone_volume_converges_to_analytic() {
    let (r, h) = (2.0, 3.0);
    let c = cone(r, 0.0, h).unwrap();
    let exact = PI * r * r * h / 3.0;
    let v = solid_volume_tol(&c, FINE).unwrap();
    assert!(v > 0.0, "cone volume must be positive, got {v}");
    assert!(
        v <= exact * (1.0 + 1e-4),
        "cone volume {v} should not exceed exact {exact}"
    );
    assert!(
        (exact - v) / exact < 0.01,
        "pointed cone volume {v} should converge within 1% of {exact} (rel err {})",
        (exact - v) / exact
    );
}

#[test]
fn frustum_volume_converges_to_analytic() {
    // Truncated cone: V = ⅓πh(R² + R·r + r²).
    let (base_r, top_r, h) = (3.0, 1.5, 4.0);
    let frustum = cone(base_r, top_r, h).unwrap();
    let exact = PI * h / 3.0 * (base_r * base_r + base_r * top_r + top_r * top_r);
    let v = solid_volume_tol(&frustum, FINE).unwrap();
    assert!(
        v <= exact * (1.0 + 1e-4),
        "frustum volume {v} should not exceed exact {exact}"
    );
    assert!(
        (exact - v) / exact < 0.01,
        "frustum volume {v} should converge within 1% of {exact} (rel err {})",
        (exact - v) / exact
    );
}

#[test]
fn cone_is_a_closed_solid() {
    let pointed = cone(1.0, 0.0, 2.0).unwrap();
    assert!(
        is_closed_solid_tol(&pointed, FINE).unwrap(),
        "a pointed cone must be a closed 2-manifold"
    );
    let frustum = cone(2.0, 1.0, 3.0).unwrap();
    assert!(
        is_closed_solid_tol(&frustum, FINE).unwrap(),
        "a frustum must be a closed 2-manifold"
    );
}

// ===========================================================================
// Torus — 2π²·R·r² volume, 4π²·R·r area. Genus 1.
// ===========================================================================

#[test]
fn torus_volume_converges_to_analytic() {
    // V = 2π²·R·r² for major radius R, minor radius r.
    let (major, minor) = (3.0, 1.0);
    let t = torus(major, minor).unwrap();
    let exact = 2.0 * PI * PI * major * minor * minor;
    let v = solid_volume_tol(&t, FINE).unwrap();
    assert!(v > 0.0, "torus volume must be positive, got {v}");
    assert!(
        v <= exact * (1.0 + 1e-4),
        "torus volume {v} should not exceed exact {exact}"
    );
    assert!(
        (exact - v) / exact < 0.02,
        "torus volume {v} should converge within 2% of {exact} (rel err {})",
        (exact - v) / exact
    );
}

#[test]
fn torus_surface_area_converges_to_analytic() {
    // A = 4π²·R·r.
    let (major, minor) = (3.0, 1.0);
    let t = torus(major, minor).unwrap();
    let exact = 4.0 * PI * PI * major * minor;
    let a = solid_area_tol(&t, FINE).unwrap();
    assert!(
        a <= exact * (1.0 + 1e-4),
        "torus area {a} should not exceed exact {exact}"
    );
    assert!(
        (exact - a) / exact < 0.02,
        "torus area {a} should converge within 2% of {exact}"
    );
}

#[test]
fn torus_is_a_closed_solid() {
    let t = torus(3.0, 1.0).unwrap();
    assert!(
        is_closed_solid_tol(&t, FINE).unwrap(),
        "a torus must be a closed 2-manifold"
    );
}

#[test]
fn torus_euler_characteristic_is_genus_one() {
    // A torus is genus-1: the Euler–Poincaré formula gives
    // χ = 2 − 2g = 0. This is the topology check that a sphere-genus
    // assumption would get wrong.
    let t = torus(3.0, 1.0).unwrap();
    assert_eq!(
        euler_characteristic(&t),
        Some(0),
        "a torus is genus-1; V−E+F must be 0"
    );
}

// ===========================================================================
// Prism — flat-faced, base-area·height volume is exact.
// ===========================================================================

#[test]
fn triangular_prism_volume_and_topology_are_exact() {
    // Right-triangle base (legs 3 and 4 → area 6), height 5 → V = 30.
    let tri = prism(&[(0.0, 0.0), (3.0, 0.0), (0.0, 4.0)], 5.0).unwrap();
    let v = solid_volume_tol(&tri, FINE).unwrap();
    assert!((v - 30.0).abs() < 1e-9, "triangular prism volume {v} != 30");

    // 2 triangular ends + 3 rectangular sides = 5 faces; 6 vertices;
    // 9 edges → χ = 6 − 9 + 5 = 2.
    assert_eq!(tri.faces(), 5, "triangular prism has 5 faces");
    assert_eq!(tri.vertices(), 6, "triangular prism has 6 vertices");
    assert_eq!(tri.edges(), 9, "triangular prism has 9 edges");
    assert_eq!(euler_characteristic(&tri), Some(2));

    assert!(
        is_closed_solid_tol(&tri, FINE).unwrap(),
        "a prism must be a closed 2-manifold"
    );
}

#[test]
fn pentagonal_prism_volume_matches_shoelace_area() {
    // A regular pentagon of circumradius 2, extruded by height 3. The
    // base area is computed by the shoelace formula on the same
    // vertices, so the prism volume must equal base_area·3.
    use std::f64::consts::TAU;
    let r = 2.0;
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for i in 0..5 {
        let a = i as f64 / 5.0 * TAU;
        pts.push((r * a.cos(), r * a.sin()));
    }
    // Shoelace area of the closed polygon.
    let mut shoelace = 0.0;
    for i in 0..pts.len() {
        let (x0, y0) = pts[i];
        let (x1, y1) = pts[(i + 1) % pts.len()];
        shoelace += x0 * y1 - x1 * y0;
    }
    let base_area = shoelace.abs() * 0.5;

    let height = 3.0;
    let p = prism(&pts, height).unwrap();
    let v = solid_volume_tol(&p, FINE).unwrap();
    assert!(
        (v - base_area * height).abs() < 1e-9,
        "pentagonal prism volume {v} != base_area·height {}",
        base_area * height
    );
    assert!(
        is_closed_solid_tol(&p, FINE).unwrap(),
        "a pentagonal prism must be a closed solid"
    );
}

// ===========================================================================
// Cross-primitive — outward orientation: every primitive has a
// positive measured volume (a flipped boundary would report negative).
// ===========================================================================

#[test]
fn every_primitive_has_outward_facing_boundary() {
    // A correctly built solid has outward-facing normals → the
    // divergence-theorem volume integral is positive. A negative value
    // would mean the boundary is inside-out — a real construction bug.
    let solids: Vec<(&str, valenx_cad::Solid)> = vec![
        ("box", box_solid(1.0, 1.0, 1.0).unwrap()),
        ("cylinder", cylinder(1.0, 2.0).unwrap()),
        ("sphere", sphere(1.0).unwrap()),
        ("pointed cone", cone(1.0, 0.0, 2.0).unwrap()),
        ("frustum", cone(2.0, 1.0, 3.0).unwrap()),
        ("torus", torus(3.0, 1.0).unwrap()),
        (
            "prism",
            prism(&[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)], 1.0).unwrap(),
        ),
    ];
    for (name, s) in solids {
        let v = solid_volume_tol(&s, 1e-3).unwrap();
        assert!(
            v > 0.0,
            "{name} has a non-positive volume {v} — boundary is inside-out"
        );
    }
}
