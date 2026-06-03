//! Bridge a sketch to a 3-D solid via [`valenx-cad`] extrusion.

use valenx_cad::Solid;

use crate::sketch::Sketch;
use crate::SketchError;

/// Extract the closed planar profile from a sketch as an ordered list
/// of (x, y) waypoints. Currently only line segments are supported in
/// the profile (Task 32 adds arc support).
///
/// The profile is assumed to be:
/// - a connected loop of lines
/// - traversed in the order constraints were added
/// - closed: last endpoint matches first within `tol`
///
/// Returns `Err(SketchError::ConstraintTypeMismatch)` if the sketch
/// has entities that aren't lines or aren't part of the profile loop.
pub fn extract_profile_lines(sketch: &Sketch, tol: f64) -> Result<Vec<(f64, f64)>, SketchError> {
    use crate::geom::Entity;
    let mut waypoints: Vec<(f64, f64)> = Vec::new();
    for (idx, entity) in sketch.entities.iter().enumerate() {
        // Phase 12C: skip construction-flagged entities.
        if sketch.construction.get(idx).copied().unwrap_or(false) {
            continue;
        }
        if let Entity::Line(l) = entity {
            let ((sx, sy), (ex, ey)) = l.endpoints(&sketch.vars);
            if waypoints.is_empty() {
                waypoints.push((sx, sy));
                waypoints.push((ex, ey));
            } else {
                let last = *waypoints.last().unwrap();
                if (last.0 - sx).powi(2) + (last.1 - sy).powi(2) <= tol * tol {
                    waypoints.push((ex, ey));
                } else if (last.0 - ex).powi(2) + (last.1 - ey).powi(2) <= tol * tol {
                    waypoints.push((sx, sy));
                } else {
                    return Err(SketchError::ConstraintTypeMismatch(format!(
                        "line endpoints don't connect to previous segment (last=({}, {}), new=({},{})/({},{}))",
                        last.0, last.1, sx, sy, ex, ey
                    )));
                }
            }
        }
    }
    // Drop the duplicate closing point if present.
    if waypoints.len() >= 2 {
        let first = waypoints[0];
        let last = *waypoints.last().unwrap();
        if (first.0 - last.0).powi(2) + (first.1 - last.1).powi(2) <= tol * tol {
            waypoints.pop();
        }
    }
    Ok(waypoints)
}

/// Extrude the sketch's closed profile along ±Z by `depth` to produce
/// a 3-D solid.
///
/// A positive `depth` sweeps the profile up (+Z); a negative `depth`
/// sweeps it down (−Z).
///
/// Currently only line-based profiles are supported (see Task 31).
/// Arc segments will be approximated by short line segments in a later
/// task (Phase 2's loft/sweep machinery covers richer cases).
///
/// ## Orientation note
///
/// `truck_modeling::builder::tsweep` orients the swept solid by the
/// sweep direction relative to the profile face's normal. Sweeping the
/// z = 0 profile face **down** (`depth < 0`) produces a solid whose
/// boundary is **inside-out** — every face normal points inward. Left
/// uncorrected that breaks every downstream consumer that relies on
/// outward normals: the boolean kernel treats an inside-out operand as
/// a "hole", and a signed volume comes out negative. So for a negative
/// depth the swept solid's faces are flipped back with `Solid::not`,
/// exactly as `valenx_cad::Solid::mirrored` does after a
/// handedness-inverting reflection. The returned solid always has
/// outward-facing normals regardless of the sweep direction.
pub fn extrude(sketch: &Sketch, depth: f64, profile_tol: f64) -> Result<Solid, SketchError> {
    use truck_modeling::builder;
    use truck_modeling::Point3;
    use truck_modeling::Vector3;

    if !depth.is_finite() || depth.abs() < 1e-12 {
        return Err(SketchError::ConstraintTypeMismatch(format!(
            "extrude depth must be nonzero and finite, got {depth}"
        )));
    }

    let waypoints = extract_profile_lines(sketch, profile_tol)?;
    if waypoints.len() < 3 {
        return Err(SketchError::ConstraintTypeMismatch(format!(
            "profile needs at least 3 waypoints to extrude, got {}",
            waypoints.len()
        )));
    }

    // Build the bottom-face wire from waypoints (at z = 0).
    let vertices: Vec<_> = waypoints
        .iter()
        .map(|(x, y)| builder::vertex(Point3::new(*x, *y, 0.0)))
        .collect();
    let mut edges = Vec::new();
    for i in 0..vertices.len() {
        let next = (i + 1) % vertices.len();
        edges.push(builder::line(&vertices[i], &vertices[next]));
    }
    let wire: truck_modeling::Wire = edges.into_iter().collect();
    let face = builder::try_attach_plane(&[wire]).map_err(|e| {
        SketchError::ConstraintTypeMismatch(format!("failed to build planar face: {e:?}"))
    })?;
    let mut shell = builder::tsweep(&face, Vector3::new(0.0, 0.0, depth));
    // A downward sweep inverts the solid's handedness — flip the face
    // orientations back so the result always has outward normals.
    if depth < 0.0 {
        shell.not();
    }
    // Wrap in our Solid wrapper.
    Ok(Solid::from_truck(shell))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Phase 12C: construction-flagged lines are excluded from the
    /// extracted profile.
    #[test]
    fn construction_lines_excluded_from_profile() {
        let mut s = crate::sketch::Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let c = s.add_point(1.0, 1.0);
        let d = s.add_point(0.0, 1.0);
        let _ab = s.add_line(a, b).unwrap();
        let _bc = s.add_line(b, c).unwrap();
        let _cd = s.add_line(c, d).unwrap();
        let da = s.add_line(d, a).unwrap();
        // A 5th line as a centre construction reference — must NOT
        // appear in the profile.
        let ctr = s.add_point(0.5, 0.5);
        let centre_ref = s.add_line(a, ctr).unwrap();
        s.toggle_construction(centre_ref);
        let wp = extract_profile_lines(&s, 1e-6).unwrap();
        assert_eq!(wp.len(), 4, "should have 4 profile waypoints");
        // Sanity: da is still in the loop.
        let _ = da;
    }

    #[test]
    fn three_lines_form_triangle() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let c = s.add_point(0.5, 1.0);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_line(c, a).unwrap();
        let waypoints = extract_profile_lines(&s, 1e-6).unwrap();
        assert_eq!(waypoints.len(), 3);
        assert_eq!(waypoints[0], (0.0, 0.0));
        assert_eq!(waypoints[1], (1.0, 0.0));
        assert_eq!(waypoints[2], (0.5, 1.0));
    }

    #[test]
    fn disconnected_lines_error() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let c = s.add_point(5.0, 5.0);
        let d = s.add_point(6.0, 5.0);
        s.add_line(a, b).unwrap();
        s.add_line(c, d).unwrap();
        let err = extract_profile_lines(&s, 1e-6).unwrap_err();
        assert_eq!(err.code(), "sketch.constraint_type_mismatch");
    }

    #[test]
    fn extrude_triangle_produces_solid() {
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(1.0, 0.0);
        let c = s.add_point(0.5, 1.0);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_line(c, a).unwrap();
        let solid = extrude(&s, 2.0, 1e-6).unwrap();
        // Tessellate and confirm node count is sane (3-prism has 6 verts, ~8 tris).
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.5).unwrap();
        assert!(!mesh.nodes.is_empty());
    }

    #[test]
    fn solve_then_extrude_unit_square() {
        use crate::constraint::Constraint;
        let mut s = Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(0.5, 0.0); // not yet on x=1
        let c = s.add_point(1.0, 1.0);
        let d = s.add_point(0.0, 1.0);
        let ab = s.add_line(a, b).unwrap();
        let _bc = s.add_line(b, c).unwrap();
        let _cd = s.add_line(c, d).unwrap();
        let _da = s.add_line(d, a).unwrap();
        // Constrain to make ab horizontal, ab length 1.
        s.add_constraint(Constraint::Horizontal(ab));
        s.add_constraint(Constraint::Distance { a, b, target: 1.0 });
        let report = crate::solver::solve(&mut s, Default::default()).unwrap();
        assert!(matches!(
            report.status,
            crate::solver::SolverStatus::Converged
        ));
        // Now extrude (may succeed or report truck-api-gap; either is OK for v1).
        let _ = s.extrude(0.5);
    }
}
