//! Boolean set-ops on solids — union, difference, intersection.
//!
//! Thin wrappers over [`truck_shapeops::and`] / [`truck_shapeops::or`].
//! All three flavours route through the same crate so callers don't
//! have to remember which truck symbol matches which set operation:
//!
//! - [`union`] (`A ∪ B`)         → `truck_shapeops::or`
//! - [`intersection`] (`A ∩ B`)  → `truck_shapeops::and`
//! - [`difference`] (`A − B`)    → `A ∩ ¬B`, implemented by inverting
//!   face orientations on `B` (`Solid::not`) and calling the AND.
//!
//! Tolerance
//! ---------
//!
//! truck-shapeops needs an explicit linear tolerance to decide when
//! two coincident-looking edges should be merged. We default to
//! [`crate::DEFAULT_BOOL_TOLERANCE`] = 0.05 model units, which works
//! well for primitives sized in the 1.0–10.0 range. Callers can pass
//! a tighter tolerance with [`union_tol`] / [`difference_tol`] /
//! [`intersection_tol`] if they're modelling at a smaller scale.
//!
//! Robustness
//! ----------
//!
//! `truck-shapeops` is not a hardened boolean kernel. On degenerate
//! input — coincident faces, tangent contact, a self-intersection, a
//! disjoint difference — it has two failure modes Valenx must contain:
//!
//! 1. **It panics.** A non-simple intermediate wire trips a
//!    `panic!` deep inside `truck-topology`. Left unguarded that
//!    unwinds the caller's thread. Every boolean here runs inside
//!    [`std::panic::catch_unwind`] so the panic is converted into a
//!    clean [`CadError::EmptyResult`].
//! 2. **It returns a phantom empty solid.** A disjoint difference
//!    (`A − B` with `B` not touching `A`) comes back as `Some(solid)`
//!    where the "solid" has zero boundary shells — a result that
//!    silently measures to volume 0. A shell-less / face-less result
//!    is caught and converted to [`CadError::EmptyResult`] too, so a
//!    boolean never returns an `Ok` solid that is not a real solid.
//!
//! The contract is therefore: a boolean either returns a genuine
//! non-empty [`Solid::Brep`], or a typed [`CadError`] — it never
//! panics and never returns a silently-invalid solid.

use std::panic::{catch_unwind, AssertUnwindSafe};

use crate::solid::{CadError, Solid};
use crate::DEFAULT_BOOL_TOLERANCE;

/// Whether a truck boolean result is a degenerate non-solid.
///
/// `truck-shapeops` can return `Some(solid)` for a result that is not
/// actually a solid — most often a disjoint difference, which comes
/// back with **no boundary shells at all**. A solid with zero shells
/// (or shells that contain zero faces) has no volume and no surface;
/// treating it as a valid `Ok` result would let a silently-wrong
/// volume-0 answer escape. This predicate flags that case so the
/// boolean wrappers can convert it to [`CadError::EmptyResult`].
fn is_degenerate_result(solid: &truck_modeling::Solid) -> bool {
    let shells = solid.boundaries();
    if shells.is_empty() {
        return true;
    }
    // Every shell empty ⇒ no boundary surface ⇒ not a real solid.
    shells.iter().all(|shell| shell.is_empty())
}

/// Run a `truck-shapeops` boolean closure with full robustness
/// containment.
///
/// Wraps the call in [`catch_unwind`] (a `truck-topology` `panic!` on
/// a non-simple wire is converted to [`CadError::EmptyResult`]) and
/// post-filters the result through `is_degenerate_result` (a phantom
/// shell-less solid becomes [`CadError::EmptyResult`] too).
///
/// `AssertUnwindSafe` is sound here: the closure only *reads* its
/// borrowed `truck` operands and `truck-shapeops` builds a fresh
/// `Option<Solid>` — a panic cannot leave a caller-visible value in a
/// torn state, because the operands are not mutated through the shared
/// borrow.
fn run_boolean<F>(op: F) -> Result<Solid, CadError>
where
    F: FnOnce() -> Option<truck_modeling::Solid>,
{
    let outcome = catch_unwind(AssertUnwindSafe(op));
    match outcome {
        // truck panicked inside (non-simple wire, degenerate topology).
        // A panic is a hard failure for that geometry — surface it as
        // an empty result, never let it unwind the caller.
        Err(_) => Err(CadError::EmptyResult),
        // truck returned no solid.
        Ok(None) => Err(CadError::EmptyResult),
        // truck returned a "solid" — accept it only if it is a real
        // one. A shell-less phantom is treated as an empty result.
        Ok(Some(solid)) => {
            if is_degenerate_result(&solid) {
                Err(CadError::EmptyResult)
            } else {
                Ok(Solid::from_inner(solid))
            }
        }
    }
}

/// Union (A ∪ B). Equivalent to "weld both solids into one".
pub fn union(a: &Solid, b: &Solid) -> Result<Solid, CadError> {
    union_tol(a, b, DEFAULT_BOOL_TOLERANCE)
}

/// Union with an explicit tolerance. Pass a smaller value for
/// tighter coincidence detection at the cost of more solver work.
///
/// Both operands must be [`Solid::Brep`] — mesh-backed solids surface
/// [`CadError::MeshBackedSolid`].
pub fn union_tol(a: &Solid, b: &Solid, tol: f64) -> Result<Solid, CadError> {
    check_tol(tol)?;
    let a_brep = a.try_inner().map_err(|_| CadError::MeshBackedSolid {
        op: "union",
        reason: "left operand is mesh-backed; rebuild it without the fillet/chamfer first"
            .to_string(),
    })?;
    let b_brep = b.try_inner().map_err(|_| CadError::MeshBackedSolid {
        op: "union",
        reason: "right operand is mesh-backed; rebuild it without the fillet/chamfer first"
            .to_string(),
    })?;
    run_boolean(|| truck_shapeops::or(a_brep, b_brep, tol))
}

/// Intersection (A ∩ B). Empty if the solids don't overlap.
pub fn intersection(a: &Solid, b: &Solid) -> Result<Solid, CadError> {
    intersection_tol(a, b, DEFAULT_BOOL_TOLERANCE)
}

/// Intersection with an explicit tolerance.
///
/// Both operands must be [`Solid::Brep`].
pub fn intersection_tol(a: &Solid, b: &Solid, tol: f64) -> Result<Solid, CadError> {
    check_tol(tol)?;
    let a_brep = a.try_inner().map_err(|_| CadError::MeshBackedSolid {
        op: "intersection",
        reason: "left operand is mesh-backed".to_string(),
    })?;
    let b_brep = b.try_inner().map_err(|_| CadError::MeshBackedSolid {
        op: "intersection",
        reason: "right operand is mesh-backed".to_string(),
    })?;
    run_boolean(|| truck_shapeops::and(a_brep, b_brep, tol))
}

/// Difference (A − B). Implemented as `A ∩ ¬B` after flipping the
/// face orientations on a clone of `B` — that's the documented
/// truck-shapeops pattern, e.g. `tests/punched_cube.rs`.
pub fn difference(a: &Solid, b: &Solid) -> Result<Solid, CadError> {
    difference_tol(a, b, DEFAULT_BOOL_TOLERANCE)
}

/// Difference with an explicit tolerance.
///
/// Both operands must be [`Solid::Brep`].
pub fn difference_tol(a: &Solid, b: &Solid, tol: f64) -> Result<Solid, CadError> {
    check_tol(tol)?;
    // Clone B and invert it. We could mutate the original in-place
    // but that'd surprise callers and rule out reusing `b` for a
    // subsequent op.
    let a_brep = a.try_inner().map_err(|_| CadError::MeshBackedSolid {
        op: "difference",
        reason: "left operand is mesh-backed".to_string(),
    })?;
    let mut b_inverted = b.clone();
    b_inverted
        .try_inner_mut()
        .map_err(|_| CadError::MeshBackedSolid {
            op: "difference",
            reason: "right operand is mesh-backed".to_string(),
        })?
        .not();
    let b_brep = b_inverted.try_inner().expect("just verified above");
    run_boolean(|| truck_shapeops::and(a_brep, b_brep, tol))
}

fn check_tol(tol: f64) -> Result<(), CadError> {
    if !tol.is_finite() {
        return Err(CadError::InvalidParam(format!(
            "boolean tolerance must be finite, got {tol}"
        )));
    }
    if tol <= 0.0 {
        return Err(CadError::InvalidParam(format!(
            "boolean tolerance must be > 0, got {tol}"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use truck_modeling::{builder, EuclideanSpace, Point3, Rad, Vector3};

    use super::*;
    use crate::primitives::box_solid;

    /// Helper: build a cylinder at an explicit `(cx, cy)` centre with
    /// the given radius and height, ranging Z=`z0..z0+height`. truck's
    /// `builder::rsweep` doesn't need its swept centre to be the origin
    /// so this is straightforward to position. Used by the boolean
    /// tests to set up overlapping primitives without re-implementing
    /// the move/translate plumbing.
    fn cylinder_at(cx: f64, cy: f64, z0: f64, radius: f64, height: f64) -> Solid {
        let v = builder::vertex(Point3::new(cx + radius, cy, z0));
        let circle = builder::rsweep(
            &v,
            Point3::new(cx, cy, z0),
            Vector3::unit_z(),
            Rad(2.0 * PI),
        );
        let disk = builder::try_attach_plane(&[circle]).unwrap();
        let inner: truck_modeling::Solid = builder::tsweep(&disk, Vector3::new(0.0, 0.0, height));
        let _ = Point3::origin(); // keep the EuclideanSpace import live
        Solid::from_inner(inner)
    }

    #[test]
    fn union_overlapping_cube_and_cylinder_yields_solid() {
        // Reproduces truck-shapeops' own `punched_cube` setup but
        // with OR instead of AND-NOT — the unit cube plus a cylinder
        // straddling its mid-face is a known-good shapeops input.
        let a = box_solid(1.0, 1.0, 1.0).unwrap();
        let b = cylinder_at(0.5, 0.5, -0.5, 0.25, 2.0);
        let u = union(&a, &b).expect("union of overlapping solids succeeds");
        assert!(u.faces() > 0, "union should have faces");
    }

    #[test]
    fn intersection_overlapping_cube_and_cylinder_yields_solid() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let cyl = cylinder_at(0.5, 0.5, -0.5, 0.25, 2.0);
        let inter = intersection(&cube, &cyl).expect("intersection succeeds");
        assert!(inter.faces() > 0);
    }

    #[test]
    fn difference_punches_cylinder_out_of_cube() {
        // The canonical "punched cube" — cube minus a cylinder. Lifted
        // straight from truck-shapeops' own integration test so we
        // know shapeops handles it. Used here to verify our
        // difference() wrapper matches truck's expectations.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let cyl = cylinder_at(0.5, 0.5, -0.5, 0.25, 2.0);
        let punched = difference(&cube, &cyl).expect("difference succeeds");
        // A punched cube has the cube's 6 outer faces plus the inner
        // hole walls — strictly more than 6.
        assert!(
            punched.faces() > 6,
            "punched cube should have more faces than a plain cube ({} vs 6)",
            punched.faces()
        );
    }

    #[test]
    fn boolean_rejects_bad_tolerance() {
        let a = box_solid(1.0, 1.0, 1.0).unwrap();
        let b = box_solid(1.0, 1.0, 1.0).unwrap();
        assert!(matches!(
            union_tol(&a, &b, -0.01),
            Err(CadError::InvalidParam(_))
        ));
        assert!(matches!(
            union_tol(&a, &b, f64::NAN),
            Err(CadError::InvalidParam(_))
        ));
    }

    #[test]
    fn union_with_mesh_backed_operand_rejects() {
        // Phase 3: mesh-backed Solids can't participate in booleans.
        // Each one must surface as MeshBackedSolid so the UI can tell
        // the user to either reorder or rebuild without the fillet.
        let brep = box_solid(1.0, 1.0, 1.0).unwrap();
        let mesh = Solid::from_mesh(valenx_mesh::Mesh::new("dummy"));
        assert!(matches!(
            union(&brep, &mesh),
            Err(CadError::MeshBackedSolid { op: "union", .. })
        ));
        assert!(matches!(
            union(&mesh, &brep),
            Err(CadError::MeshBackedSolid { op: "union", .. })
        ));
        assert!(matches!(
            difference(&mesh, &brep),
            Err(CadError::MeshBackedSolid {
                op: "difference",
                ..
            })
        ));
        assert!(matches!(
            intersection(&brep, &mesh),
            Err(CadError::MeshBackedSolid {
                op: "intersection",
                ..
            })
        ));
    }

    #[test]
    fn self_difference_does_not_panic_and_surfaces_empty() {
        // A − A trips a "non-simple wire" panic deep inside
        // truck-topology. `run_boolean`'s catch_unwind must convert
        // that into a clean EmptyResult — never let it unwind the
        // caller's thread.
        let a = box_solid(2.0, 2.0, 2.0).unwrap();
        let b = box_solid(2.0, 2.0, 2.0).unwrap();
        match difference(&a, &b) {
            Err(CadError::EmptyResult) => {}
            other => panic!("A−A should surface EmptyResult, got {other:?}"),
        }
    }

    #[test]
    fn disjoint_difference_does_not_return_a_phantom_solid() {
        // truck-shapeops returns Some(shell-less solid) for a disjoint
        // A − B. is_degenerate_result must catch the empty-shell
        // phantom and convert it to EmptyResult so a silently-wrong
        // volume-0 "solid" never escapes as Ok.
        let a = box_solid(1.0, 1.0, 1.0).unwrap();
        let b = box_solid(1.0, 1.0, 1.0)
            .unwrap()
            .translated(50.0, 50.0, 50.0)
            .unwrap();
        match difference(&a, &b) {
            Err(CadError::EmptyResult) => {}
            Ok(_) => panic!("disjoint A−B must not return an Ok phantom solid"),
            Err(other) => panic!("disjoint A−B should be EmptyResult, got {other:?}"),
        }
    }

    #[test]
    fn is_degenerate_result_flags_a_shell_less_solid() {
        // A truck Solid with no boundary shells is not a real solid.
        let empty = truck_modeling::Solid::new(Vec::new());
        assert!(
            is_degenerate_result(&empty),
            "a zero-shell solid must be flagged degenerate"
        );
        // A genuine primitive is not degenerate.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        assert!(
            !is_degenerate_result(cube.try_inner().unwrap()),
            "a real box must not be flagged degenerate"
        );
    }
}
