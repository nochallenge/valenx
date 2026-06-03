//! Phase 12E Task 46 — op regression: each sketch-level op preserves
//! solvability when applied to constrained entities.

use crate::constraint::Constraint;
use crate::ops::{copy, linear_array, mirror, polar_array, r#move, rotate};
use crate::sketch::Sketch;
use crate::solver::{solve, SolverStatus};

fn build_constrained_pair() -> Sketch {
    // Two points + a line between them with a Distance(2) constraint.
    let mut s = Sketch::new();
    let a = s.add_point(0.0, 0.0);
    let b = s.add_point(2.0, 0.0);
    let _ab = s.add_line(a, b).unwrap();
    s.add_constraint(Constraint::Distance { a, b, target: 2.0 });
    s
}

#[test]
fn copy_preserves_solvability() {
    let mut s = build_constrained_pair();
    // Copy the line at offset (5, 0).
    let new = copy::copy(&mut s, &[crate::geom::EntityId(3)], (5.0, 0.0));
    assert_eq!(new.len(), 1);
    let r = solve(&mut s, Default::default()).unwrap();
    assert!(matches!(r.status, SolverStatus::Converged));
}

#[test]
fn translate_preserves_solvability() {
    let mut s = build_constrained_pair();
    r#move::translate(&mut s, &[crate::geom::EntityId(3)], (1.0, 0.0));
    let r = solve(&mut s, Default::default()).unwrap();
    assert!(matches!(r.status, SolverStatus::Converged));
}

#[test]
fn rotate_preserves_solvability() {
    let mut s = build_constrained_pair();
    rotate::rotate(&mut s, &[crate::geom::EntityId(3)], (0.0, 0.0), 0.5);
    let r = solve(&mut s, Default::default()).unwrap();
    assert!(matches!(r.status, SolverStatus::Converged));
}

#[test]
fn mirror_preserves_solvability() {
    let mut s = build_constrained_pair();
    let line = mirror::MirrorLine {
        point: (0.0, 0.0),
        direction: (0.0, 1.0),
    };
    let _ = mirror::mirror(&mut s, &[crate::geom::EntityId(3)], &line);
    let r = solve(&mut s, Default::default()).unwrap();
    assert!(matches!(r.status, SolverStatus::Converged));
}

#[test]
fn linear_array_preserves_solvability() {
    let mut s = build_constrained_pair();
    let _ = linear_array::linear_array(&mut s, &[crate::geom::EntityId(3)], (1.0, 0.0), 3, 5.0);
    let r = solve(&mut s, Default::default()).unwrap();
    assert!(matches!(r.status, SolverStatus::Converged));
}

#[test]
fn polar_array_preserves_solvability() {
    let mut s = build_constrained_pair();
    let _ = polar_array::polar_array(
        &mut s,
        &[crate::geom::EntityId(3)],
        (10.0, 0.0),
        4,
        std::f64::consts::TAU,
    );
    let r = solve(&mut s, Default::default()).unwrap();
    assert!(matches!(r.status, SolverStatus::Converged));
}
