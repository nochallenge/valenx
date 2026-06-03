//! Shared machinery for the feature-based `feat_make_*` family
//! (Phases 97-100).
//!
//! OCCT's `BRepFeat_*` builders all share the same two-step shape:
//! (1) build the feature body, oriented onto a sketch plane / spine in
//! world space; (2) **fuse** it onto (boss) or **subtract** it from
//! (pocket) an existing base solid. This module provides:
//!
//! - [`orient_z_to`] — the rigid transform carrying a `+Z`-canonical
//!   body onto an arbitrary `(origin, normal)` frame.
//! - [`feature_combine`] — the topology-aware add/subtract: a real
//!   BRep boolean when both operands are BReps, with a co-refinement
//!   mesh-CSG fallback (via `valenx-cgal-port`) when an operand is
//!   mesh-backed or the BRep boolean fails.
//!
//! The fallback means the feature ops always return *real* carved
//! geometry — never a stub — degrading from exact-BRep to
//! mesh-domain only when the BRep kernel cannot proceed.

use valenx_cad::Solid;
use valenx_cgal_port::mesh_boolean::{self, Mesh3};
use valenx_cgal_port::Triangle3;

use crate::error::OcctSurfaceError;

/// Tessellation chord tolerance for the mesh-CSG fallback.
const FEAT_TESS_TOLERANCE: f64 = 0.05;

/// Rigidly transform a body built canonically along `+Z` (base at the
/// origin) so its axis lies along `normal` and its base sits at
/// `origin`.
///
/// `normal` need not be unit. Returns the body unchanged when
/// `normal` is already `+Z`; handles the antiparallel degenerate
/// case. Topology is preserved (rigid transform).
pub fn orient_z_to(body: Solid, origin: [f64; 3], normal: [f64; 3]) -> Solid {
    // Callers (feat_make_prism, prim_api_cylinder, …) pre-validate
    // that `origin` and `normal` are finite. The translated/rotated
    // round-6 fallible variants are therefore safe to `expect()` here:
    // an infinite/NaN input means the caller forgot to validate, which
    // is a programming error worth surfacing loudly.
    let len = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    if len < 1e-12 {
        return body
            .translated(origin[0], origin[1], origin[2])
            .expect("orient_z_to: origin must be finite (caller invariant)");
    }
    let t = [normal[0] / len, normal[1] / len, normal[2] / len];
    let dot = t[2].clamp(-1.0, 1.0); // +Z · t
    let oriented = if dot > 1.0 - 1e-12 {
        body
    } else if dot < -1.0 + 1e-12 {
        body.rotated((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), std::f64::consts::PI)
            .expect("orient_z_to: 180° about +X is always valid")
    } else {
        // Rotation axis = (+Z) × t, angle = acos(dot).
        let axis = (-t[1], t[0], 0.0);
        let al = (axis.0 * axis.0 + axis.1 * axis.1 + axis.2 * axis.2).sqrt();
        let unit = (axis.0 / al, axis.1 / al, axis.2 / al);
        body.rotated((0.0, 0.0, 0.0), unit, dot.acos())
            .expect("orient_z_to: rotation parameters derived from validated normal")
    };
    oriented
        .translated(origin[0], origin[1], origin[2])
        .expect("orient_z_to: origin must be finite (caller invariant)")
}

/// Combine a feature body with a base solid: `fuse = true` adds the
/// feature (boss), `fuse = false` removes it (pocket).
///
/// Prefers the exact BRep boolean (`valenx_cad::union` /
/// `difference`); if either operand is mesh-backed, or the BRep
/// boolean returns an empty / failed result, falls back to the
/// `valenx-cgal-port` co-refinement mesh CSG and returns a
/// mesh-backed [`Solid`].
///
/// # Errors
///
/// [`OcctSurfaceError::TruckLimit`] only when even the mesh fallback
/// cannot tessellate an operand.
pub fn feature_combine(
    base: &Solid,
    feature: &Solid,
    fuse: bool,
) -> Result<Solid, OcctSurfaceError> {
    // Fast path: both BRep → try the exact BRep boolean.
    if matches!(base, Solid::Brep(_)) && matches!(feature, Solid::Brep(_)) {
        let brep_result = if fuse {
            valenx_cad::union(base, feature)
        } else {
            valenx_cad::difference(base, feature)
        };
        if let Ok(s) = brep_result {
            return Ok(s);
        }
        // BRep boolean failed (EmptyResult / degenerate) — fall
        // through to the mesh CSG rather than erroring out.
    }
    // Mesh-domain fallback: tessellate both, run the co-refinement CSG.
    let base_soup = solid_to_soup(base)?;
    let feat_soup = solid_to_soup(feature)?;
    let result = if fuse {
        mesh_boolean::union(&base_soup, &feat_soup)
    } else {
        mesh_boolean::difference(&base_soup, &feat_soup)
    };
    Ok(Solid::from_mesh(soup_to_mesh(&result)))
}

/// Tessellate a solid into a `valenx-cgal-port` triangle soup.
fn solid_to_soup(solid: &Solid) -> Result<Mesh3, OcctSurfaceError> {
    let mesh = valenx_cad::solid_to_mesh(solid, FEAT_TESS_TOLERANCE)
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("feature: tessellate: {e:?}")))?;
    let mut soup = Mesh3::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let mut v = [nalgebra::Vector3::zeros(); 3];
            let mut ok = true;
            for (k, &idx) in tri.iter().enumerate() {
                match mesh.nodes.get(idx as usize) {
                    Some(p) => v[k] = *p,
                    None => {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                soup.triangles.push(Triangle3 { v });
            }
        }
    }
    Ok(soup)
}

/// Convert a `valenx-cgal-port` triangle soup back into a welded
/// [`valenx_mesh::Mesh`].
fn soup_to_mesh(soup: &Mesh3) -> valenx_mesh::Mesh {
    let (verts, faces) = mesh_boolean::mesh3_to_indexed(soup, 1e-6);
    let mut mesh = valenx_mesh::Mesh::new("feature");
    mesh.nodes = verts;
    let mut conn: Vec<u32> = Vec::new();
    for f in &faces {
        conn.extend_from_slice(&[f[0] as u32, f[1] as u32, f[2] as u32]);
    }
    mesh.element_blocks.push(valenx_mesh::element::ElementBlock {
        element_type: valenx_mesh::ElementType::Tri3,
        connectivity: conn,
    });
    mesh.recompute_stats();
    mesh
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn orient_z_to_plus_z_is_just_translation() {
        let body = box_solid(1.0, 1.0, 1.0).unwrap();
        let moved = orient_z_to(body, [5.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        // Rigid transform preserves topology.
        assert_eq!(moved.faces(), 6);
    }

    #[test]
    fn feature_combine_fuse_of_two_brep_boxes() {
        // Two overlapping BRep boxes — the fuse should yield a single
        // valid solid.
        let a = box_solid(2.0, 2.0, 2.0).unwrap();
        let b = box_solid(2.0, 2.0, 2.0)
            .unwrap()
            .translated(1.0, 0.0, 0.0)
            .unwrap();
        let fused = feature_combine(&a, &b, true).unwrap();
        // The result tessellates to non-empty geometry.
        let mesh = valenx_cad::solid_to_mesh(&fused, 0.2).unwrap();
        assert!(!mesh.nodes.is_empty());
    }

    #[test]
    fn feature_combine_subtract_carves_the_base() {
        let a = box_solid(4.0, 4.0, 4.0).unwrap();
        // A bar poking all the way through.
        let b = box_solid(1.0, 1.0, 8.0)
            .unwrap()
            .translated(1.0, 1.0, -2.0)
            .unwrap();
        let carved = feature_combine(&a, &b, false).unwrap();
        let mesh = valenx_cad::solid_to_mesh(&carved, 0.2).unwrap();
        assert!(!mesh.nodes.is_empty());
    }
}
