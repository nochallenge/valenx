//! Flow-field visualization — turning a converged [`AeroResult`] into
//! coloured geometry the app's 3-D viewport can render.
//!
//! The viewport already ships a per-vertex scalar colour-ramp path
//! (used by the FEM stress visualisation): it draws a
//! [`valenx_mesh::Mesh`] tinted by a paired [`valenx_fields::Field`]
//! through the cool-to-warm colormap, with a colour-bar legend. This
//! module produces exactly that pair from aero data:
//!
//! - **Surface fields** (Cp, skin-friction) — one small quad per
//!   immersed-body surface face, the scalar value carried on its four
//!   vertices, so the body shell is painted with the field.
//! - **Cut-plane fields** (velocity magnitude, static pressure,
//!   Q-criterion) — a flat rectangular grid patch sampled from a
//!   [`valenx_aero::FieldSlice`], one cell per slice sample.
//!
//! The result is a `(Mesh, Field)` pair; the workbench drops the mesh
//! into `ValenxApp::mesh` and the field into the dedicated
//! `aero_field_overlay` slot the viewport picks up.

use nalgebra::Vector3;

use valenx_aero::{slice_field, surface_field, AeroResult, FieldSlice, SliceAxis};
use valenx_fields::{Field, FieldKind, Location, RegionRef, TimeKey, Units};
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use super::model::{CutAxis, FlowField};

/// A built flow-visualization payload — a triangle mesh plus the scalar
/// field that colours it.
pub struct FlowVizMesh {
    /// The geometry — body surface shell or a cut-plane patch.
    pub mesh: Mesh,
    /// The per-node scalar field, ready for the viewport colour ramp.
    pub field: Field,
}

/// Build the flow-visualization payload for a chosen field.
///
/// `result` is the converged solve; `field` selects which scalar to
/// paint; `cut_axis` / `cut_fraction` place the cut plane for the
/// volumetric fields (ignored for the surface fields). Returns an
/// error string if the field cannot be built (e.g. an empty body).
pub fn build_flow_viz(
    result: &AeroResult,
    field: FlowField,
    cut_axis: CutAxis,
    cut_fraction: f64,
) -> Result<FlowVizMesh, String> {
    match field {
        FlowField::SurfaceCp => surface_viz(result, true),
        FlowField::SkinFriction => surface_viz(result, false),
        FlowField::VelocityMagnitude => cutplane_viz(result, field, cut_axis, cut_fraction),
        FlowField::PressureSlice => cutplane_viz(result, field, cut_axis, cut_fraction),
        FlowField::VortexQ => cutplane_viz(result, field, cut_axis, cut_fraction),
    }
}

// ---------------------------------------------------------------------------
// Surface fields — Cp / Cf painted on the body shell
// ---------------------------------------------------------------------------

/// Build a body-surface mesh coloured by `Cp` (when `use_cp`) or by the
/// skin-friction coefficient `Cf`.
///
/// `valenx-aero`'s [`surface_field`] returns one [`SurfacePoint`] per
/// exposed immersed-body face — a centre, an outward normal, the scalar
/// values. Each becomes a small square quad (two triangles) lying in
/// the body surface, with the scalar carried on all four corners so the
/// viewport's per-vertex ramp paints it flat per face.
fn surface_viz(result: &AeroResult, use_cp: bool) -> Result<FlowVizMesh, String> {
    let points = surface_field(&result.tunnel, &result.flow);
    if points.is_empty() {
        return Err("the body voxelized to no surface faces — try a finer grid".to_string());
    }
    // Each face's quad spans one grid cell. Use the grid cell size as
    // the patch half-extent so the quads tile the staircased surface.
    let g = result.tunnel.grid;
    let half = 0.5 * g.dx().min(g.dy()).min(g.dz());

    let mut mesh = Mesh::new("aero-surface-field");
    let mut block = ElementBlock::new(ElementType::Tri3);
    let mut values: Vec<f64> = Vec::with_capacity(points.len() * 4);

    for p in &points {
        let (u, v) = perpendicular_basis(p.normal);
        // Four corners of the face quad, in the body surface plane.
        let base = mesh.nodes.len() as u32;
        for &(su, sv) in &[(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)] {
            let corner = p.position + (su * half) * u + (sv * half) * v;
            mesh.nodes.push(corner);
        }
        // Two CCW triangles for the quad.
        block
            .connectivity
            .extend_from_slice(&[base, base + 1, base + 2, base, base + 2, base + 3]);
        let scalar = if use_cp { p.cp } else { p.cf };
        // The same value on all four corners → a flat-shaded face.
        for _ in 0..4 {
            values.push(scalar);
        }
    }
    mesh.element_blocks.push(block);
    mesh.recompute_stats();

    let (name, units) = if use_cp {
        ("Cp", valenx_fields::units::DIMENSIONLESS)
    } else {
        ("Cf", valenx_fields::units::DIMENSIONLESS)
    };
    let field = make_node_field(name, units, values);
    Ok(FlowVizMesh { mesh, field })
}

// ---------------------------------------------------------------------------
// Cut-plane fields — velocity / pressure / Q on an axis-aligned slice
// ---------------------------------------------------------------------------

/// Build a cut-plane patch coloured by a volumetric field.
fn cutplane_viz(
    result: &AeroResult,
    field: FlowField,
    cut_axis: CutAxis,
    cut_fraction: f64,
) -> Result<FlowVizMesh, String> {
    let g = result.tunnel.grid;
    let axis = cut_axis.slice_axis();
    let frac = cut_fraction.clamp(0.0, 1.0);
    // The world coordinate of the cut plane along its normal axis.
    let coordinate = match axis {
        SliceAxis::X => g.x0 + frac * g.lx,
        SliceAxis::Y => g.y0 + frac * g.ly,
        SliceAxis::Z => g.z0 + frac * g.lz,
    };

    // The cell-centred scalar field to slice.
    let (cell_field, name, units) = match field {
        FlowField::VelocityMagnitude => {
            let mut f = g.scalar_field();
            for k in 0..g.nz {
                for j in 0..g.ny {
                    for i in 0..g.nx {
                        f.set(i, j, k, result.flow.speed_at_cell(i, j, k));
                    }
                }
            }
            (f, "|U|", valenx_fields::units::METER_PER_SECOND)
        }
        FlowField::PressureSlice => (
            result.flow.pressure.clone(),
            "p",
            valenx_fields::units::PASCAL,
        ),
        FlowField::VortexQ => {
            // The Q-criterion field marks vortex cores; slice it.
            let q = vorticity_q(result);
            (q, "Q", valenx_fields::units::DIMENSIONLESS)
        }
        // The surface fields never reach here.
        FlowField::SurfaceCp | FlowField::SkinFriction => {
            return Err("internal: a surface field routed to the cut-plane builder".into())
        }
    };

    let slice = slice_field(&g, &cell_field, axis, coordinate);
    let (mesh, values) = slice_to_mesh(&slice, &g, coordinate);
    let field = make_node_field(name, units, values);
    Ok(FlowVizMesh { mesh, field })
}

/// The Q-criterion field for a result — a thin wrapper so `cutplane_viz`
/// stays readable.
fn vorticity_q(result: &AeroResult) -> valenx_aero::Field3 {
    valenx_aero::q_criterion(&result.flow)
}

/// Turn a [`FieldSlice`] into a flat rectangular grid mesh.
///
/// The slice is a `width × height` array of cell values; this builds a
/// `width × height` quad grid in world space at the cut plane, one cell
/// quad per slice sample. Returns the mesh plus a node-aligned scalar
/// `Vec`: each cell's value is carried onto all **four** of its corner
/// nodes, so the colour reads flat per cell (matching the slice's
/// cell-centred sampling) and `field.data.len() == mesh.nodes.len()`.
fn slice_to_mesh(slice: &FieldSlice, g: &valenx_aero::Grid3, coordinate: f64) -> (Mesh, Vec<f64>) {
    let mut mesh = Mesh::new("aero-cut-plane");
    let mut block = ElementBlock::new(ElementType::Tri3);
    let mut values = Vec::new();
    let (w, h) = (slice.width.max(1), slice.height.max(1));

    // The in-plane cell sizes + the world position of cell (a, b).
    let cell = |a: usize, b: usize| -> Vector3<f64> {
        match slice.axis {
            SliceAxis::X => Vector3::new(
                coordinate,
                g.y0 + (a as f64 + 0.5) * g.dy(),
                g.z0 + (b as f64 + 0.5) * g.dz(),
            ),
            SliceAxis::Y => Vector3::new(
                g.x0 + (a as f64 + 0.5) * g.dx(),
                coordinate,
                g.z0 + (b as f64 + 0.5) * g.dz(),
            ),
            SliceAxis::Z => Vector3::new(
                g.x0 + (a as f64 + 0.5) * g.dx(),
                g.y0 + (b as f64 + 0.5) * g.dy(),
                coordinate,
            ),
        }
    };
    // Half-extents of one cell quad along the two in-plane axes.
    let (du, dv) = match slice.axis {
        SliceAxis::X => (g.dy(), g.dz()),
        SliceAxis::Y => (g.dx(), g.dz()),
        SliceAxis::Z => (g.dx(), g.dy()),
    };
    let (au, av): (Vector3<f64>, Vector3<f64>) = match slice.axis {
        SliceAxis::X => (Vector3::new(0.0, 1.0, 0.0), Vector3::new(0.0, 0.0, 1.0)),
        SliceAxis::Y => (Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0)),
        SliceAxis::Z => (Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)),
    };

    for b in 0..h {
        for a in 0..w {
            let c = cell(a, b);
            let base = mesh.nodes.len() as u32;
            for &(su, sv) in &[(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)] {
                mesh.nodes.push(c + (su * du) * au + (sv * dv) * av);
            }
            // The cell-centred slice value, carried onto all four corner
            // nodes so the field stays node-aligned with the mesh.
            let v = slice
                .values
                .get(a + b * slice.width)
                .copied()
                .unwrap_or(0.0);
            values.extend_from_slice(&[v, v, v, v]);
            block.connectivity.extend_from_slice(&[
                base,
                base + 1,
                base + 2,
                base,
                base + 2,
                base + 3,
            ]);
        }
    }
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    (mesh, values)
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build an OnNode scalar [`Field`] with a recomputed value range.
///
/// `values` already carries one entry per mesh node (four per face quad
/// or per cell quad), so the field is a node-located scalar the
/// viewport's per-vertex ramp consumes directly.
fn make_node_field(name: &str, units: Units, values: Vec<f64>) -> Field {
    let mut field = Field {
        name: name.to_string(),
        kind: FieldKind::Scalar,
        location: Location::OnNode,
        region: RegionRef("aero".to_string()),
        units,
        time: TimeKey::Steady,
        data: values,
        range: None,
    };
    field.recompute_range();
    field
}

/// An orthonormal basis pair perpendicular to `normal`.
///
/// Used to lay a face quad flat in the body surface plane. The seed
/// axis is picked to be the least aligned with `normal` so the cross
/// product is well-conditioned; a degenerate (zero) normal falls back
/// to the world x/y plane.
fn perpendicular_basis(normal: Vector3<f64>) -> (Vector3<f64>, Vector3<f64>) {
    let n = normal.try_normalize(1e-12).unwrap_or_else(Vector3::z);
    let seed = if n.x.abs() < 0.9 {
        Vector3::new(1.0, 0.0, 0.0)
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    };
    let u = n
        .cross(&seed)
        .try_normalize(1e-12)
        .unwrap_or_else(Vector3::x);
    let v = n.cross(&u).try_normalize(1e-12).unwrap_or_else(Vector3::y);
    (u, v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_aero::{geometry::box_body, run_windtunnel, AeroRequest, TurbulenceModel};

    /// A converged solve over a small box — the fixture every viz test
    /// builds against.
    fn solved_box() -> AeroResult {
        let body = box_body(Vector3::new(-0.5, -0.5, -0.5), Vector3::new(0.5, 0.5, 0.5));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_max_iterations(12);
        run_windtunnel(&body, &req).expect("box solve")
    }

    #[test]
    fn perpendicular_basis_is_orthonormal() {
        for n in [
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(0.3, 0.7, -0.2),
        ] {
            let (u, v) = perpendicular_basis(n);
            let nn = n.normalize();
            assert!((u.norm() - 1.0).abs() < 1e-9);
            assert!((v.norm() - 1.0).abs() < 1e-9);
            // u and v are perpendicular to the normal and to each other.
            assert!(u.dot(&nn).abs() < 1e-9);
            assert!(v.dot(&nn).abs() < 1e-9);
            assert!(u.dot(&v).abs() < 1e-9);
        }
    }

    #[test]
    fn perpendicular_basis_survives_a_zero_normal() {
        // A degenerate normal must not produce NaNs.
        let (u, v) = perpendicular_basis(Vector3::zeros());
        assert!(u.iter().all(|c| c.is_finite()));
        assert!(v.iter().all(|c| c.is_finite()));
    }

    #[test]
    fn surface_cp_viz_paints_the_body_shell() {
        let result = solved_box();
        let viz =
            build_flow_viz(&result, FlowField::SurfaceCp, CutAxis::Y, 0.5).expect("surface Cp viz");
        // The mesh has triangles and the field has one value per node.
        assert!(!viz.mesh.nodes.is_empty());
        assert_eq!(viz.field.data.len(), viz.mesh.nodes.len());
        assert_eq!(viz.field.location, Location::OnNode);
        assert_eq!(viz.field.kind, FieldKind::Scalar);
        // The node count is a multiple of four (one quad per face).
        assert_eq!(viz.mesh.nodes.len() % 4, 0);
        // Every Cp value is finite.
        assert!(viz.field.data.iter().all(|v| v.is_finite()));
        // The field has a populated range.
        assert!(viz.field.range.is_some());
    }

    #[test]
    fn skin_friction_viz_builds_a_distinct_field() {
        let result = solved_box();
        let cp = build_flow_viz(&result, FlowField::SurfaceCp, CutAxis::Y, 0.5).unwrap();
        let cf = build_flow_viz(&result, FlowField::SkinFriction, CutAxis::Y, 0.5).unwrap();
        // Same geometry topology, different field name.
        assert_eq!(cp.mesh.nodes.len(), cf.mesh.nodes.len());
        assert_eq!(cf.field.name, "Cf");
        assert!(cf.field.data.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn velocity_cutplane_viz_builds_a_grid_patch() {
        let result = solved_box();
        let viz = build_flow_viz(&result, FlowField::VelocityMagnitude, CutAxis::Y, 0.5)
            .expect("velocity cut-plane viz");
        assert!(!viz.mesh.nodes.is_empty());
        assert_eq!(viz.field.data.len(), viz.mesh.nodes.len());
        // A cut plane is a grid of cell quads — node count is a
        // multiple of four.
        assert_eq!(viz.mesh.nodes.len() % 4, 0);
        assert_eq!(viz.field.name, "|U|");
        assert!(viz.field.data.iter().all(|v| v.is_finite() && *v >= 0.0));
    }

    #[test]
    fn pressure_and_q_cutplanes_build_for_every_axis() {
        let result = solved_box();
        for axis in CutAxis::ALL {
            for field in [FlowField::PressureSlice, FlowField::VortexQ] {
                let viz = build_flow_viz(&result, field, axis, 0.5)
                    .unwrap_or_else(|e| panic!("{field:?} on {axis:?}: {e}"));
                assert!(!viz.mesh.nodes.is_empty());
                assert_eq!(viz.field.data.len(), viz.mesh.nodes.len());
                assert!(viz.field.data.iter().all(|v| v.is_finite()));
            }
        }
    }

    #[test]
    fn cut_fraction_is_clamped() {
        // An out-of-range cut fraction must not panic — it is clamped
        // into [0, 1] before placing the plane.
        let result = solved_box();
        let lo = build_flow_viz(&result, FlowField::PressureSlice, CutAxis::X, -3.0);
        let hi = build_flow_viz(&result, FlowField::PressureSlice, CutAxis::X, 9.0);
        assert!(lo.is_ok());
        assert!(hi.is_ok());
    }

    #[test]
    fn surface_field_data_is_flat_per_face() {
        // The four corners of every face quad must carry the same value
        // (a flat-shaded face) — so consecutive groups of four are
        // equal.
        let result = solved_box();
        let viz = build_flow_viz(&result, FlowField::SurfaceCp, CutAxis::Y, 0.5).unwrap();
        for chunk in viz.field.data.chunks_exact(4) {
            assert!((chunk[0] - chunk[1]).abs() < 1e-12);
            assert!((chunk[0] - chunk[2]).abs() < 1e-12);
            assert!((chunk[0] - chunk[3]).abs() < 1e-12);
        }
    }
}
