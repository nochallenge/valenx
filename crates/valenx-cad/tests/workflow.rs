//! End-to-end integration tests for `valenx-cad`.
//!
//! Exercises the full pipeline a desktop user takes — build a
//! primitive, position it, boolean-combine it, then tessellate to a
//! triangle mesh suitable for the viewport — through the public API
//! surface only. If any of these break, the Part workbench is
//! broken.

use valenx_cad::{
    box_solid, cone, cylinder, difference, fillet_edges, intersection, prism, solid_to_mesh,
    sphere, torus, union, CadError, DEFAULT_BOOL_TOLERANCE, DEFAULT_TESS_TOLERANCE,
};

#[test]
fn defaults_match_module_constants() {
    // Cheap sanity check that nothing has redefined the constants
    // under us — the toolbox UI hardcodes these so users will get a
    // surprise if the kernel and the UI drift apart.
    assert_eq!(DEFAULT_TESS_TOLERANCE, 0.5);
    assert_eq!(DEFAULT_BOOL_TOLERANCE, 0.05);
}

#[test]
fn full_workflow_box_minus_cylinder_then_mesh() {
    // Build a unit cube at the origin, then a cylinder centred at
    // (0.5, 0.5) and offset along -Z so it punches through the cube.
    // This is the canonical "punched cube" workflow the Part toolbox
    // advertises — keep it green.
    let cube = box_solid(1.0, 1.0, 1.0).expect("cube");
    let drill = cylinder(0.25, 2.0)
        .expect("cylinder")
        .translated(0.5, 0.5, -0.5)
        .unwrap();
    let punched = difference(&cube, &drill).expect("difference");
    assert!(
        punched.faces() > cube.faces(),
        "punched cube should have more faces than the original"
    );

    let mesh = solid_to_mesh(&punched, DEFAULT_TESS_TOLERANCE).expect("mesh");
    assert!(
        mesh.total_elements() > 0,
        "the punched-cube tessellation should produce triangles"
    );
}

#[test]
fn primitives_each_round_trip_through_tessellation() {
    let primitives: Vec<(&str, valenx_cad::Solid)> = vec![
        ("box", box_solid(2.0, 1.0, 1.0).unwrap()),
        ("cylinder", cylinder(1.0, 2.0).unwrap()),
        ("sphere", sphere(1.0).unwrap()),
        ("pointed cone", cone(1.0, 0.0, 2.0).unwrap()),
        ("frustum", cone(2.0, 1.0, 2.0).unwrap()),
        ("torus", torus(2.0, 0.5).unwrap()),
        (
            "triangular prism",
            prism(&[(0.0, 0.0), (2.0, 0.0), (1.0, 1.0)], 1.0).unwrap(),
        ),
    ];
    for (name, solid) in primitives {
        let mesh = solid_to_mesh(&solid, 0.1)
            .unwrap_or_else(|e| panic!("{name}: tessellation failed: {e}"));
        assert!(
            !mesh.nodes.is_empty(),
            "{name}: tessellated mesh should have nodes"
        );
        assert!(
            mesh.total_elements() > 0,
            "{name}: tessellated mesh should have triangles"
        );
    }
}

#[test]
fn fillet_returns_typed_not_implemented_not_a_silent_pass() {
    // The whole point of returning NotImplemented from fillet_edges
    // is so callers can detect the unsupported case at runtime. Make
    // sure that contract holds.
    let cube = box_solid(1.0, 1.0, 1.0).unwrap();
    match fillet_edges(&cube, 0.05) {
        Err(CadError::NotImplemented { op, reason }) => {
            assert_eq!(op, "fillet_edges");
            assert!(!reason.is_empty());
        }
        other => panic!("expected NotImplemented, got {other:?}"),
    }
}

#[test]
fn box_dimensions_round_trip_through_tessellation() {
    // Tessellate a 3×2×1 box and check the resulting AABB matches
    // — sanity that we didn't scramble axes anywhere.
    let b = box_solid(3.0, 2.0, 1.0).unwrap();
    let mesh = solid_to_mesh(&b, 0.1).unwrap();
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for n in &mesh.nodes {
        for i in 0..3 {
            if n[i] < min[i] {
                min[i] = n[i];
            }
            if n[i] > max[i] {
                max[i] = n[i];
            }
        }
    }
    assert!((max[0] - min[0] - 3.0).abs() < 1e-6, "X span 3.0");
    assert!((max[1] - min[1] - 2.0).abs() < 1e-6, "Y span 2.0");
    assert!((max[2] - min[2] - 1.0).abs() < 1e-6, "Z span 1.0");
}

#[test]
fn translation_moves_aabb() {
    // Cube at origin, translated by (10, 20, 30). The resulting mesh
    // bounding box should be exactly offset.
    let a = box_solid(1.0, 1.0, 1.0).unwrap();
    let b = a.translated(10.0, 20.0, 30.0).unwrap();
    let mesh = solid_to_mesh(&b, 0.1).unwrap();
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for n in &mesh.nodes {
        for i in 0..3 {
            if n[i] < min[i] {
                min[i] = n[i];
            }
            if n[i] > max[i] {
                max[i] = n[i];
            }
        }
    }
    assert!((min[0] - 10.0).abs() < 1e-6);
    assert!((min[1] - 20.0).abs() < 1e-6);
    assert!((min[2] - 30.0).abs() < 1e-6);
}

#[test]
fn union_yields_more_or_equal_faces_than_either_input() {
    let a = box_solid(1.0, 1.0, 1.0).unwrap();
    // Offset the drill so it doesn't share faces with cube boundary.
    let b = cylinder(0.25, 2.0)
        .unwrap()
        .translated(0.5, 0.5, -0.5)
        .unwrap();
    let u = union(&a, &b).expect("union");
    // After welding two overlapping solids, the face count should be
    // at least as large as either operand.
    assert!(u.faces() >= a.faces().min(b.faces()));
}

#[test]
fn intersection_face_count_is_positive() {
    let a = box_solid(1.0, 1.0, 1.0).unwrap();
    let b = cylinder(0.25, 2.0)
        .unwrap()
        .translated(0.5, 0.5, -0.5)
        .unwrap();
    let inter = intersection(&a, &b).expect("intersection");
    assert!(inter.faces() > 0);
}

#[test]
fn rotation_preserves_topology() {
    // Rotating a box around the X axis by 90 degrees produces a
    // solid with the same face / edge / vertex counts.
    let a = box_solid(1.0, 2.0, 3.0).unwrap();
    let b = a
        .rotated(
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            std::f64::consts::FRAC_PI_2,
        )
        .unwrap();
    assert_eq!(a.faces(), b.faces());
    assert_eq!(a.edges(), b.edges());
    assert_eq!(a.vertices(), b.vertices());
}
