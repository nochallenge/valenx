//! Parametric **L-bracket** producer for the Workbench+Agent 3-D workspace
//! tile.
//!
//! Builds a recognizable structural L-bracket entirely in-house over the
//! native `valenx-cad` BRep kernel — an extruded L-profile with one rounded
//! re-entrant (concave) corner and two M5 clearance holes drilled through —
//! then tessellates it to a [`crate::types::LoadedMesh`] and exports a binary
//! STL alongside. The single source of truth for the agent-bridge bracket
//! product (see [`crate::agent_commands::AgentCommand::Show3d`]
//! `kind:"bracket"`).
//!
//! ## Geometry
//!
//! - Footprint **80 × 60 mm**, **5 mm** thick (extruded along +z).
//! - Two **20 mm**-wide legs forming the L (outer bottom leg 80×20, outer
//!   vertical leg 20×60).
//! - The inner re-entrant corner is **rounded** by baking a short tangent arc
//!   (`FILLET_R` mm, [`FILLET_PTS`] points) straight into the 2-D profile —
//!   an honest *profile round*, not a 3-D edge blend. (`valenx_cad::fillet_edges`
//!   is a documented NotImplemented stub and is deliberately **not** called.)
//! - Two **M5** holes (Ø5 → r 2.5 mm), one centred in each leg, cut with
//!   [`valenx_cad::difference`] of a through-cylinder.
//!
//! ## STL
//!
//! [`valenx_mesh::write_stl_binary`] writes the tessellated surface to a
//! fixed temp path (`<temp>/valenx_bracket.stl`). A write error is reported in
//! the readout rows rather than panicking.

use std::f64::consts::FRAC_PI_2;
use std::path::PathBuf;

use valenx_viz::OrbitCamera;

use crate::types::LoadedMesh;

/// Bracket footprint along +x (mm).
const LEN_X: f64 = 80.0;
/// Bracket footprint along +y (mm).
const LEN_Y: f64 = 60.0;
/// Leg width (mm) — both the horizontal and vertical leg are this wide.
const LEG_W: f64 = 20.0;
/// Plate thickness, extruded along +z (mm).
const THICK: f64 = 5.0;
/// Radius of the rounded re-entrant (inner) corner (mm).
const FILLET_R: f64 = 6.0;
/// Number of arc vertices baked in for the rounded inner corner.
const FILLET_PTS: usize = 6;
/// M5 clearance-hole radius (mm).
const HOLE_R: f64 = 2.5;

/// Build the closed L-profile polygon in the X-Y plane (CCW), with the inner
/// re-entrant corner replaced by a short tangent arc of [`FILLET_PTS`] points.
///
/// Outer outline (CCW):
/// `(0,0) → (LEN_X,0) → (LEN_X,LEG_W) → ⟪arc⟫ → (LEG_W,LEN_Y) → (0,LEN_Y) →`
/// close. The re-entrant corner the arc rounds is the concave vertex at
/// `(LEG_W, LEG_W)`; the arc is centred at `(LEG_W+r, LEG_W+r)` and runs from
/// `(LEG_W+r, LEG_W)` (tangent to the bottom-leg top edge) to
/// `(LEG_W, LEG_W+r)` (tangent to the vertical-leg inner edge).
fn l_profile() -> Vec<(f64, f64)> {
    let mut p: Vec<(f64, f64)> = Vec::with_capacity(5 + FILLET_PTS);
    p.push((0.0, 0.0));
    p.push((LEN_X, 0.0));
    p.push((LEN_X, LEG_W));
    // Rounded re-entrant corner: sweep the arc CCW (as seen from +z) from the
    // bottom-edge tangent point to the left-edge tangent point. Centre at
    // (cx, cy); start angle -90° (pointing down to (cx, LEG_W)), end angle
    // 180° (pointing left to (LEG_W, cy)). Going CCW from -90° to 180° is a
    // +270° sweep that bulges away from the centre — wrong. We instead sweep
    // the *concave* quarter: from angle 0° (point (cx+? )) ... simplest is to
    // parametrize the inner quarter-circle directly.
    let cx = LEG_W + FILLET_R;
    let cy = LEG_W + FILLET_R;
    // Angles measured from +x at the arc centre. The bottom tangent point
    // (cx, LEG_W) is at angle -90° (i.e. 270°); the left tangent point
    // (LEG_W, cy) is at angle 180°. Walk from 270° down to 180° (a -90° / 90°
    // CCW-in-profile sweep) so the inserted vertices replace the sharp corner
    // with a concave round that stays inside the material.
    for k in 0..=FILLET_PTS {
        let t = k as f64 / FILLET_PTS as f64; // 0..=1
        let ang = 1.5 * std::f64::consts::PI - t * FRAC_PI_2; // 270° → 180°
        p.push((cx + FILLET_R * ang.cos(), cy + FILLET_R * ang.sin()));
    }
    p.push((LEG_W, LEN_Y));
    p.push((0.0, LEN_Y));
    p
}

/// Build the L-bracket solid, tessellate it to a [`LoadedMesh`], write a
/// binary STL beside it, and return the mesh paired with the readout rows
/// (dims, hole spec, STL path / write status).
///
/// Infallible at the API level: the profile + primitives are valid by
/// construction, so the prism / difference / tessellation succeed. A failed
/// STL write degrades to a note in the returned rows (no panic).
pub(crate) fn bracket_loaded_mesh() -> (LoadedMesh, Vec<String>) {
    // Extrude the rounded L-profile to a 5 mm plate.
    let profile = l_profile();
    let body = valenx_cad::prism(&profile, THICK).expect("valid L-profile prism");

    // Two M5 through-holes, one centred in each leg. Cylinders are built on
    // the X-Y plane about the origin (axis +z); translate each to its hole
    // centre, then difference it from the body. Extend the cutter slightly
    // beyond both faces so the boolean cleanly punches through.
    let bottom_leg_hole = (LEN_X - LEG_W * 0.5, LEG_W * 0.5); // (70, 10)
    let vert_leg_hole = (LEG_W * 0.5, LEN_Y - LEG_W * 0.5); // (10, 50)
    let mut solid = body;
    for (hx, hy) in [bottom_leg_hole, vert_leg_hole] {
        let cutter = valenx_cad::cylinder(HOLE_R, THICK * 2.0)
            .expect("valid hole cylinder")
            .translated(hx, hy, -THICK * 0.5)
            .expect("finite hole translation");
        solid = valenx_cad::difference(&solid, &cutter).expect("hole boolean succeeds");
    }

    let mesh = valenx_cad::solid_to_mesh(&solid, valenx_cad::DEFAULT_TESS_TOLERANCE)
        .expect("bracket solid tessellates");

    // Write a binary STL next to the model (fixed name under the temp dir).
    let stl_path = std::env::temp_dir().join("valenx_bracket.stl");
    let stl_line = match valenx_mesh::write_stl_binary(&mesh, &stl_path) {
        Ok(()) => format!("STL: {}", stl_path.display()),
        Err(e) => format!("STL write failed: {e}"),
    };

    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    let loaded = LoadedMesh {
        path: PathBuf::from("<bracket>/valenx-l-bracket"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    };

    let lines = vec![
        format!("L-bracket: {LEN_X:.0} × {LEN_Y:.0} mm footprint, {THICK:.0} mm thick"),
        format!("legs: 2 × {LEG_W:.0} mm wide"),
        format!("rounded inner corner: r {FILLET_R:.0} mm (profile round)"),
        format!("holes: 2 × M5 (Ø5.0 mm) — one per leg"),
        stl_line,
    ];
    (loaded, lines)
}

/// A fixed 3/4-view [`OrbitCamera`] framing the bracket `mesh` (same
/// `frame_bounds` fit + hero angle as [`crate::rocket_workbench::lv1_camera`]),
/// for the Workbench+Agent bracket product's per-tile 3-D view.
pub(crate) fn bracket_camera(mesh: &valenx_mesh::Mesh) -> OrbitCamera {
    let mut camera = OrbitCamera::default();
    if let Some((min, max)) = crate::mesh_loader::mesh_bounding_box(mesh) {
        camera.frame_bounds(min, max);
    }
    camera.azimuth_deg = 35.0;
    camera.elevation_deg = 22.0;
    camera
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profile_is_closed_and_has_the_arc() {
        let p = l_profile();
        // 5 corner vertices of the L (minus the rounded one) + (FILLET_PTS+1)
        // arc vertices.
        assert_eq!(p.len(), 5 + (FILLET_PTS + 1));
        // First/last are distinct (prism auto-closes; caller must not repeat).
        assert_ne!(p.first(), p.last());
        // Arc endpoints are tangent to the two leg edges.
        let arc_start = (LEG_W + FILLET_R, LEG_W);
        let arc_end = (LEG_W, LEG_W + FILLET_R);
        assert!(p
            .iter()
            .any(|&(x, y)| (x - arc_start.0).abs() < 1e-9 && (y - arc_start.1).abs() < 1e-9));
        assert!(p
            .iter()
            .any(|&(x, y)| (x - arc_end.0).abs() < 1e-9 && (y - arc_end.1).abs() < 1e-9));
    }

    #[test]
    fn bracket_builds_a_nonempty_surface_mesh() {
        let (loaded, lines) = bracket_loaded_mesh();
        assert!(!loaded.mesh.nodes.is_empty(), "bracket has vertices");
        assert!(loaded.mesh.total_elements() > 0, "bracket has triangles");
        assert!(
            lines.iter().any(|l| l.contains("M5")),
            "hole spec line present"
        );
    }
}
