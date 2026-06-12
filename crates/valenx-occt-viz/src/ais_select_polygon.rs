//! Phase 174 — `AIS_InteractiveContext::SelectPolygon()` — lasso /
//! freehand polygon selection.
//!
//! ## What OCCT does
//!
//! The user clicks a series of polygon vertices (or drags a freehand
//! curve that's sampled to vertices). OCCT then tests each selectable
//! object's screen-projected centroid against polygon-inclusion via
//! a winding-number test. Objects whose centroid lies inside the
//! polygon are promoted to `Selected`.
//!
//! ## v1 status — real lasso selection
//!
//! With the picking substrate from Phase 171 in place, this is a real
//! selector. Each object's geometry centroid is projected to screen
//! pixels with the installed camera; a ray-crossing point-in-polygon
//! test decides inclusion; matching objects are promoted to
//! `Selected`.
//!
//! Objects with no registered geometry, or whose centroid is behind
//! the camera, are skipped. With no camera installed nothing is
//! selected.

use crate::ais_interactive_context::{InteractiveContext, ObjectState};
use crate::error::OcctVizError;

/// Select every object in `ctx` whose screen-projected centroid falls
/// inside `polygon`. Returns the list of newly-selected IDs.
///
/// `polygon` is an ordered list of 2D screen-pixel vertices; the
/// polygon closes implicitly between `polygon.last()` and
/// `polygon[0]`.
///
/// # Errors
///
/// [`OcctVizError::BadInput`] if `polygon.len() < 3` or any vertex is
/// non-finite.
pub fn ais_select_polygon(
    ctx: &mut InteractiveContext,
    polygon: &[[f32; 2]],
) -> Result<Vec<usize>, OcctVizError> {
    if polygon.len() < 3 {
        return Err(OcctVizError::bad_input(
            "polygon",
            format!("need ≥ 3 vertices; got {}", polygon.len()),
        ));
    }
    for (i, v) in polygon.iter().enumerate() {
        if !v[0].is_finite() || !v[1].is_finite() {
            return Err(OcctVizError::bad_input(
                "polygon",
                format!("vertex {i} contains non-finite coordinate"),
            ));
        }
    }

    let mut selected: Vec<usize> = Vec::new();
    for id in ctx.geometry_ids() {
        if ctx.state(id) == Some(ObjectState::Hidden) {
            continue;
        }
        let Some(geom) = ctx.geometry(id) else {
            continue;
        };
        let centroid = geom.centroid();
        let Some((sx, sy)) = ctx.project(centroid) else {
            continue; // behind the camera
        };
        if point_in_polygon(sx, sy, polygon) {
            ctx.set_state(id, ObjectState::Selected);
            selected.push(id);
        }
    }
    selected.sort();
    Ok(selected)
}

/// Ray-crossing (even-odd) point-in-polygon test.
fn point_in_polygon(px: f32, py: f32, poly: &[[f32; 2]]) -> bool {
    let mut inside = false;
    let n = poly.len();
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = (poly[i][0], poly[i][1]);
        let (xj, yj) = (poly[j][0], poly[j][1]);
        // Does the horizontal ray from (px, py) cross edge i-j?
        let crosses = (yi > py) != (yj > py) && px < (xj - xi) * (py - yi) / (yj - yi) + xi;
        if crosses {
            inside = !inside;
        }
        j = i;
    }
    inside
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ais_interactive_context::{ais_interactive_context, Pickable};
    use nalgebra::{Matrix4, Point3, Vector3};

    fn test_view() -> Matrix4<f64> {
        let view = Matrix4::look_at_rh(
            &Point3::new(0.0, 0.0, 10.0),
            &Point3::origin(),
            &Vector3::y(),
        );
        let proj = Matrix4::new_perspective(1.0, 60f64.to_radians(), 0.1, 100.0);
        proj * view
    }

    /// A point primitive at the world origin → projects to the screen
    /// centre.
    fn origin_point() -> Pickable {
        Pickable::Point {
            position: [0.0, 0.0, 0.0],
        }
    }

    #[test]
    fn rejects_too_few_vertices() {
        let mut ctx = ais_interactive_context().unwrap();
        let err = ais_select_polygon(&mut ctx, &[[0.0, 0.0], [10.0, 0.0]]).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_non_finite_vertex() {
        let mut ctx = ais_interactive_context().unwrap();
        let err =
            ais_select_polygon(&mut ctx, &[[0.0, 0.0], [10.0, 0.0], [f32::NAN, 10.0]]).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn point_in_polygon_basic() {
        let square = [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]];
        assert!(point_in_polygon(5.0, 5.0, &square));
        assert!(!point_in_polygon(15.0, 5.0, &square));
        assert!(!point_in_polygon(-1.0, 5.0, &square));
    }

    #[test]
    fn lasso_around_centre_selects_origin_object() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(origin_point());
        // A polygon covering the screen centre (origin projects to
        // ~(50, 50)).
        let poly = [[20.0, 20.0], [80.0, 20.0], [80.0, 80.0], [20.0, 80.0]];
        let sel = ais_select_polygon(&mut ctx, &poly).unwrap();
        assert_eq!(sel, vec![id]);
        assert_eq!(ctx.state(id), Some(ObjectState::Selected));
    }

    #[test]
    fn lasso_off_to_the_side_selects_nothing() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        ctx.display_geometry(origin_point());
        let poly = [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0]];
        let sel = ais_select_polygon(&mut ctx, &poly).unwrap();
        assert!(sel.is_empty());
    }

    #[test]
    fn no_camera_selects_nothing() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.display_geometry(origin_point());
        let poly = [[0.0, 0.0], [100.0, 0.0], [100.0, 100.0], [0.0, 100.0]];
        let sel = ais_select_polygon(&mut ctx, &poly).unwrap();
        assert!(sel.is_empty());
    }
}
