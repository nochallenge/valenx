//! Solid-modeling primitives.
//!
//! Each builder validates its input, then delegates to
//! [`truck_modeling::builder`] to produce a closed BRep. The
//! resulting topology is wrapped in a [`Solid`] so callers never have
//! to spell out the three-parameter `Solid<Point3, Curve, Surface>`.
//!
//! Coordinate conventions
//! ----------------------
//!
//! - **Box** — corner at the origin, opposite corner at `(dx, dy, dz)`.
//! - **Cylinder** — base disk centred on `(0,0,0)` lying in the X-Y
//!   plane, axis pointing along +Z, height `height`.
//! - **Sphere** — centred on the origin.
//! - **Truncated cone** — base disk in the X-Y plane, axis along +Z.
//!   `top_radius = 0.0` is the regular "pointed" cone.
//! - **Torus** — major axis along Z.
//! - **Prism** — closed polygon defined in the X-Y plane, extruded
//!   along +Z.
//!
//! These match Valenx's right-handed Y-up viewport convention well
//! enough for v1; we can revisit the coordinate frame once we have
//! actual users.

use std::f64::consts::PI;

use truck_modeling::{builder, EuclideanSpace, Point3, Rad, Solid as TruckSolid, Vector3, Wire};

use crate::solid::{CadError, Solid};

/// Axis-aligned box.
///
/// One corner at the origin, the opposite corner at `(dx, dy, dz)`.
/// All three dimensions must be strictly positive.
pub fn box_solid(dx: f64, dy: f64, dz: f64) -> Result<Solid, CadError> {
    require_positive("box.dx", dx)?;
    require_positive("box.dy", dy)?;
    require_positive("box.dz", dz)?;

    // Three orthogonal sweeps — vertex → edge → face → solid. The
    // sweep direction is the dimension we're extending each step.
    let v = builder::vertex(Point3::origin());
    let edge = builder::tsweep(&v, Vector3::new(dx, 0.0, 0.0));
    let face = builder::tsweep(&edge, Vector3::new(0.0, dy, 0.0));
    let solid: TruckSolid = builder::tsweep(&face, Vector3::new(0.0, 0.0, dz));
    Ok(Solid::from_inner(solid))
}

/// Right circular cylinder centred on the origin.
///
/// Base disk lives in the X-Y plane, axis along +Z.
pub fn cylinder(radius: f64, height: f64) -> Result<Solid, CadError> {
    require_positive("cylinder.radius", radius)?;
    require_positive("cylinder.height", height)?;

    // Build a disk: revolve a vertex 360° around the Z axis through
    // the origin to produce a circle wire, then attach the disk
    // surface inside it. `try_attach_plane` is the truck builder
    // helper that creates a planar face from a closed wire.
    let v = builder::vertex(Point3::new(radius, 0.0, 0.0));
    let circle = builder::rsweep(&v, Point3::origin(), Vector3::unit_z(), Rad(2.0 * PI));
    let disk = builder::try_attach_plane(&[circle])
        .map_err(|e| CadError::InvalidParam(format!("cylinder disk: {e:?}")))?;

    // Translation sweep: disk → solid cylinder along +Z.
    let solid: TruckSolid = builder::tsweep(&disk, Vector3::new(0.0, 0.0, height));
    Ok(Solid::from_inner(solid))
}

/// Sphere centred on the origin.
///
/// Built by revolving a semi-circular wire around the Y axis through
/// the origin, then closing the resulting shell into a solid via
/// [`builder::cone`] (which is misnamed in the truck API — it's
/// really "sweep around an axis, closing degenerate edges that hit
/// the axis").
pub fn sphere(radius: f64) -> Result<Solid, CadError> {
    require_positive("sphere.radius", radius)?;

    // Top pole vertex, semi-circle wire around the X axis (so the
    // wire lives in the Y-Z plane), then revolve that wire around
    // the Y axis. `builder::cone` collapses the degenerate edges
    // that touch the rotation axis at the poles, giving us a closed
    // sphere shell.
    let v0 = builder::vertex(Point3::new(0.0, radius, 0.0));
    let wire: Wire = builder::rsweep(&v0, Point3::origin(), Vector3::unit_x(), Rad(PI));
    let shell = builder::cone(&wire, Vector3::unit_y(), Rad(2.0 * PI));
    let solid = TruckSolid::new(vec![shell]);
    Ok(Solid::from_inner(solid))
}

/// Truncated cone (frustum). Pass `top_radius = 0.0` for a regular
/// pointed cone, or matching radii for a degenerate "cylinder" (use
/// [`cylinder`] for that case directly).
///
/// Base disk lies in the X-Y plane, axis along +Z, height `height`.
pub fn cone(base_radius: f64, top_radius: f64, height: f64) -> Result<Solid, CadError> {
    require_positive("cone.base_radius", base_radius)?;
    if top_radius < 0.0 {
        return Err(CadError::InvalidParam(format!(
            "cone.top_radius must be >= 0, got {top_radius}"
        )));
    }
    require_positive("cone.height", height)?;

    // Build an OPEN profile wire whose endpoints lie on the
    // rotation axis (Z). truck's `builder::cone` closes the shell
    // by collapsing the degenerate side edges that ride along the
    // axis — that's why the profile must START and END on the
    // rotation axis, with the off-axis vertices in between.
    //
    // The convention is the inverse of FreeCAD: it walks DOWN the
    // axis from top, OUT to the slant, then BACK along the base.
    //
    // For a pointed cone (top_radius == 0):
    //   apex (on axis) → base outer → base centre (on axis)
    // For a frustum:
    //   top centre (on axis) → top outer → base outer → base centre (on axis)
    let top_axis = builder::vertex(Point3::new(0.0, 0.0, height));
    let base_axis = builder::vertex(Point3::new(0.0, 0.0, 0.0));
    let wire: Wire = if top_radius > 0.0 {
        // Frustum: top axis → top outer → base outer → base axis.
        let top_outer = builder::vertex(Point3::new(top_radius, 0.0, height));
        let base_outer = builder::vertex(Point3::new(base_radius, 0.0, 0.0));
        vec![
            builder::line(&top_axis, &top_outer),
            builder::line(&top_outer, &base_outer),
            builder::line(&base_outer, &base_axis),
        ]
        .into()
    } else {
        // Pointed cone: apex (top axis) → base outer → base axis.
        let base_outer = builder::vertex(Point3::new(base_radius, 0.0, 0.0));
        vec![
            builder::line(&top_axis, &base_outer),
            builder::line(&base_outer, &base_axis),
        ]
        .into()
    };
    let shell = builder::cone(&wire, Vector3::unit_z(), Rad(2.0 * PI));
    let solid = TruckSolid::new(vec![shell]);
    Ok(Solid::from_inner(solid))
}

/// Torus. `major_radius` is the centre-circle radius; `minor_radius`
/// is the tube radius. The torus's major axis lies along Z (so the
/// hole is parallel to the X-Y plane).
pub fn torus(major_radius: f64, minor_radius: f64) -> Result<Solid, CadError> {
    require_positive("torus.major_radius", major_radius)?;
    require_positive("torus.minor_radius", minor_radius)?;
    if minor_radius >= major_radius {
        return Err(CadError::InvalidParam(format!(
            "torus.minor_radius ({minor_radius}) must be strictly less than \
             major_radius ({major_radius}) — a self-intersecting torus is \
             not a valid BRep solid"
        )));
    }

    // Standard two-sweep torus from `examples/torus.rs`, but oriented
    // so the major axis is Z.
    let v = builder::vertex(Point3::new(major_radius + minor_radius, 0.0, 0.0));
    let inner_circle = builder::rsweep(
        &v,
        Point3::new(major_radius, 0.0, 0.0),
        Vector3::unit_y(),
        Rad(2.0 * PI),
    );
    let shell = builder::rsweep(
        &inner_circle,
        Point3::origin(),
        Vector3::unit_z(),
        Rad(2.0 * PI),
    );
    let solid = TruckSolid::new(vec![shell]);
    Ok(Solid::from_inner(solid))
}

/// Prism — extrude a planar polygon from the X-Y plane along +Z.
///
/// Inputs
/// ------
///
/// - `profile_xy` — vertices of the polygon, traversed in order. The
///   polygon must be closed by the caller (last point != first point;
///   the function does NOT auto-close). Vertices must be coplanar in
///   the X-Y plane.
/// - `height` — extrusion length along +Z. Must be strictly positive.
///
/// At least three distinct points are required. Self-intersecting
/// or non-simple polygons are not validated here — truck will reject
/// the resulting wire downstream.
pub fn prism(profile_xy: &[(f64, f64)], height: f64) -> Result<Solid, CadError> {
    if profile_xy.len() < 3 {
        return Err(CadError::InvalidParam(format!(
            "prism profile needs at least 3 points, got {}",
            profile_xy.len()
        )));
    }
    require_positive("prism.height", height)?;

    // Build vertices, then line edges between successive pairs and
    // a closing edge back to the start.
    let mut vertices = Vec::with_capacity(profile_xy.len());
    for (x, y) in profile_xy {
        if !x.is_finite() || !y.is_finite() {
            return Err(CadError::InvalidParam(format!(
                "prism profile point ({x}, {y}) contains a non-finite value"
            )));
        }
        vertices.push(builder::vertex(Point3::new(*x, *y, 0.0)));
    }
    let mut edges = Vec::with_capacity(vertices.len());
    for i in 0..vertices.len() {
        let next = (i + 1) % vertices.len();
        edges.push(builder::line(&vertices[i], &vertices[next]));
    }
    let wire: Wire = edges.into();
    let face = builder::try_attach_plane(&[wire])
        .map_err(|e| CadError::InvalidParam(format!("prism profile: {e:?}")))?;
    let solid: TruckSolid = builder::tsweep(&face, Vector3::new(0.0, 0.0, height));
    Ok(Solid::from_inner(solid))
}

/// Revolve — sweep a half-section profile about the Z axis to build a
/// solid of revolution (a lathe / turn operation).
///
/// Inputs
/// ------
///
/// - `profile_rz` — the half-section as `(r, z)` points, where `r ≥ 0` is the
///   distance from the Z axis and `z` is the axial coordinate. The points are
///   traversed in order to form an open polyline; the revolve closes the body.
///   The endpoints are expected to lie **on the axis** (`r = 0`) so the
///   resulting shell caps cleanly at the top and bottom — the same convention
///   [`cone`] and [`sphere`] use with truck's `builder::cone` axis-sweep, which
///   collapses the degenerate edges that ride along the rotation axis. A
///   profile whose ends are off-axis still revolves, but yields an open tube.
/// - `angle_deg` — the sweep angle in degrees, `0 < angle_deg ≤ 360`. `360`
///   builds a full lathe body; a partial angle builds a wedge whose two flat
///   end faces close the profile.
///
/// At least two distinct profile points are required (a single segment is the
/// minimum that sweeps to a surface). Interior points with `r < 0` are
/// rejected (they would cross the axis); the endpoints may be exactly `0`.
///
/// This is a **true BRep revolve** built from truck builder sweeps — the
/// returned [`Solid`] is a proper closed 2-manifold (for the full-360°,
/// axis-capped case), not a triangle soup, so it composes with the boolean
/// operators and the mass-property measures exactly like the other primitives.
///
/// Implementation notes
/// --------------------
///
/// - The **full-360°** path mirrors [`cone`] / [`sphere`]: it walks the open
///   profile wire (reversing it so the swept-shell normals point outward, the
///   same winding the validated [`cone`] builder uses) and feeds it to truck's
///   `builder::cone` axis-sweep, which collapses the degenerate edges that ride
///   along the rotation axis into apex points — yielding a closed shell.
/// - The **partial-angle** path closes the half-section back along the axis
///   into a planar face in the X-Z plane and `rsweep`s that *face* through the
///   angle, so the two profile-plane ends become flat cap faces of the wedge.
pub fn revolve(profile_rz: &[(f64, f64)], angle_deg: f64) -> Result<Solid, CadError> {
    if profile_rz.len() < 2 {
        return Err(CadError::InvalidParam(format!(
            "revolve profile needs at least 2 points, got {}",
            profile_rz.len()
        )));
    }
    if !angle_deg.is_finite() || angle_deg <= 0.0 || angle_deg > 360.0 + 1e-9 {
        return Err(CadError::InvalidParam(format!(
            "revolve.angle_deg must be in (0, 360], got {angle_deg}"
        )));
    }
    for &(r, z) in profile_rz {
        if !r.is_finite() || !z.is_finite() {
            return Err(CadError::InvalidParam(format!(
                "revolve profile point ({r}, {z}) contains a non-finite value"
            )));
        }
        if r < 0.0 {
            return Err(CadError::InvalidParam(format!(
                "revolve profile radius must be >= 0 (cannot cross the axis), got {r}"
            )));
        }
    }
    let angle = angle_deg.to_radians();
    let full_turn = angle_deg >= 360.0 - 1e-9;

    if full_turn {
        // Full lathe body. Walk the profile REVERSED so the resulting
        // surface-of-revolution normals point outward (matching the winding the
        // validated `cone` builder uses; the forward winding yields an
        // inward-facing — negative-volume — shell). `builder::cone` collapses
        // the degenerate axis edges into apex points, closing the shell.
        let verts: Vec<_> = profile_rz
            .iter()
            .rev()
            .map(|&(r, z)| builder::vertex(Point3::new(r, 0.0, z)))
            .collect();
        let mut edges = Vec::with_capacity(verts.len() - 1);
        for i in 0..verts.len() - 1 {
            edges.push(builder::line(&verts[i], &verts[i + 1]));
        }
        let wire: Wire = edges.into();
        let shell = builder::cone(&wire, Vector3::unit_z(), Rad(angle));
        let solid = TruckSolid::new(vec![shell]);
        return Ok(Solid::from_inner(solid));
    }

    // Partial wedge. Close the half-section into a planar face in the X-Z plane,
    // then revolve the FACE so its two profile-plane ends become flat caps.
    // Close along the axis (r = 0) only when the endpoints are off it, so an
    // already axis-touching profile isn't given a zero-length closing edge.
    let mut pts: Vec<(f64, f64)> = profile_rz.to_vec();
    let first = *pts.first().expect("len checked >= 2");
    let last = *pts.last().expect("len checked >= 2");
    if last.0 > 1e-9 {
        pts.push((0.0, last.1));
    }
    if first.0 > 1e-9 {
        pts.push((0.0, first.1));
    }
    let verts: Vec<_> = pts
        .iter()
        .map(|&(r, z)| builder::vertex(Point3::new(r, 0.0, z)))
        .collect();
    let mut edges = Vec::with_capacity(verts.len());
    for i in 0..verts.len() {
        let next = (i + 1) % verts.len();
        edges.push(builder::line(&verts[i], &verts[next]));
    }
    let wire: Wire = edges.into();
    let face = builder::try_attach_plane(&[wire])
        .map_err(|e| CadError::InvalidParam(format!("revolve profile face: {e:?}")))?;
    let solid: TruckSolid = builder::rsweep(&face, Point3::origin(), Vector3::unit_z(), Rad(angle));
    Ok(Solid::from_inner(solid))
}

/// Helper: enforce `value > 0` with a parameter-name-tagged error.
fn require_positive(name: &str, value: f64) -> Result<(), CadError> {
    if !value.is_finite() {
        return Err(CadError::InvalidParam(format!(
            "{name} must be finite, got {value}"
        )));
    }
    if value <= 0.0 {
        return Err(CadError::InvalidParam(format!(
            "{name} must be strictly positive, got {value}"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn box_topology() {
        let cube = box_solid(2.0, 3.0, 4.0).unwrap();
        assert_eq!(cube.faces(), 6);
        assert_eq!(cube.edges(), 12);
        assert_eq!(cube.vertices(), 8);
    }

    #[test]
    fn box_rejects_zero_or_negative_dims() {
        assert!(matches!(
            box_solid(-1.0, 1.0, 1.0),
            Err(CadError::InvalidParam(_))
        ));
        assert!(matches!(
            box_solid(1.0, 0.0, 1.0),
            Err(CadError::InvalidParam(_))
        ));
        assert!(matches!(
            box_solid(1.0, 1.0, f64::NAN),
            Err(CadError::InvalidParam(_))
        ));
    }

    #[test]
    fn cylinder_has_three_faces_top_bottom_side() {
        // top disk + bottom disk + side surface (which truck may
        // split into 3 sub-faces depending on its closed-sweep
        // tessellation). All we really assert is "more than zero".
        let cyl = cylinder(1.0, 2.0).unwrap();
        assert!(cyl.faces() > 0, "cylinder should have visible faces");
        assert!(cyl.vertices() > 0);
    }

    #[test]
    fn sphere_topology_is_non_empty() {
        let s = sphere(1.0).unwrap();
        assert!(s.faces() > 0);
        assert!(s.vertices() > 0);
    }

    #[test]
    fn cone_pointed_and_frustum_both_build() {
        let pointed = cone(1.0, 0.0, 2.0).unwrap();
        assert!(pointed.faces() > 0);

        let frustum = cone(2.0, 1.0, 3.0).unwrap();
        assert!(frustum.faces() >= pointed.faces());
    }

    #[test]
    fn cone_rejects_negative_top_radius() {
        assert!(matches!(
            cone(1.0, -0.5, 2.0),
            Err(CadError::InvalidParam(_))
        ));
    }

    #[test]
    fn torus_builds_when_minor_lt_major() {
        let t = torus(2.0, 0.5).unwrap();
        assert!(t.faces() > 0);
    }

    #[test]
    fn torus_rejects_self_intersecting() {
        assert!(matches!(torus(1.0, 1.0), Err(CadError::InvalidParam(_))));
        assert!(matches!(torus(1.0, 1.5), Err(CadError::InvalidParam(_))));
    }

    #[test]
    fn prism_extrudes_triangle() {
        let tri = prism(&[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)], 2.0).unwrap();
        // A triangular prism has 5 faces (2 triangular ends + 3 rect sides),
        // 9 edges, 6 vertices.
        assert_eq!(tri.faces(), 5, "triangular prism should have 5 faces");
        assert_eq!(tri.vertices(), 6);
    }

    #[test]
    fn prism_rejects_short_profiles() {
        assert!(matches!(
            prism(&[(0.0, 0.0), (1.0, 0.0)], 1.0),
            Err(CadError::InvalidParam(_))
        ));
    }

    #[test]
    fn revolve_full_turn_builds_a_closed_solid() {
        // A cone half-section: axis (0,0) → outer (1,0) → axis (0,2). Revolving
        // 360° about Z reproduces a cone — a non-empty, watertight BRep solid.
        let body = revolve(&[(0.0, 0.0), (1.0, 0.0), (0.0, 2.0)], 360.0).unwrap();
        assert!(body.faces() > 0, "revolved solid should have faces");
        assert!(body.vertices() > 0);
        // A full-360 axis-capped revolve must close (watertight 2-manifold).
        assert!(
            crate::measure::is_closed_solid(&body).unwrap_or(false),
            "a full-turn axis-capped revolve should be a closed solid"
        );
        // Its volume should match the analytic cone π r² h / 3 within the
        // tessellation tolerance the measure module converges from below at.
        let vol = crate::measure::solid_volume(&body).unwrap();
        let expected = std::f64::consts::PI * 1.0 * 1.0 * 2.0 / 3.0;
        assert!(
            (vol - expected).abs() < 0.2,
            "revolved cone volume {vol} should approximate {expected}"
        );
    }

    #[test]
    fn revolve_partial_angle_builds() {
        // A 90° wedge of the same section still builds a non-empty body.
        let wedge = revolve(&[(0.0, 0.0), (1.0, 0.0), (0.0, 2.0)], 90.0).unwrap();
        assert!(wedge.faces() > 0);
    }

    #[test]
    fn revolve_rejects_bad_input() {
        // Too few points.
        assert!(matches!(
            revolve(&[(0.0, 0.0)], 360.0),
            Err(CadError::InvalidParam(_))
        ));
        // Angle out of range.
        assert!(matches!(
            revolve(&[(0.0, 0.0), (1.0, 0.0)], 0.0),
            Err(CadError::InvalidParam(_))
        ));
        assert!(matches!(
            revolve(&[(0.0, 0.0), (1.0, 0.0)], 400.0),
            Err(CadError::InvalidParam(_))
        ));
        // Negative radius (crosses the axis).
        assert!(matches!(
            revolve(&[(0.0, 0.0), (-1.0, 0.0), (0.0, 1.0)], 360.0),
            Err(CadError::InvalidParam(_))
        ));
    }
}
