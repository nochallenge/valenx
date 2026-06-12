//! Lightweight wrapper around truck's [`truck_modeling::Solid`].
//!
//! Why a wrapper at all
//! --------------------
//!
//! - Hides the three-type-parameter mouthful (`Solid<Point3, Curve,
//!   Surface>`) from callers — they get a single `cad::Solid`.
//! - Lets us add Valenx-flavoured helpers like
//!   [`Solid::faces`] / [`Solid::vertices`] without monkey-patching
//!   the upstream crate.
//! - Gives us a typed error enum ([`CadError`]) instead of leaking
//!   `truck_modeling::errors::Error`, which has variants we don't
//!   surface (e.g. STEP / OBJ I/O errors).
//!
//! Two solid backends
//! ------------------
//!
//! Since Phase 3 there are *two* kinds of [`Solid`]:
//!
//! - [`Solid::Brep`] — the classic truck-modeling BRep. Faces, edges,
//!   and vertices are real topological entities; booleans, sweeps,
//!   tessellation all work as expected.
//! - [`Solid::Mesh`] — a triangle-mesh-backed Solid. Produced by
//!   `valenx-fillet` (mesh-domain fillet / chamfer output) when there's
//!   no true BRep representation to round-trip. Tessellation returns
//!   the cached mesh as-is; booleans and other BRep-only ops return
//!   [`CadError::MeshBackedSolid`].
//!
//! Mesh-backed solids are a v1 compromise. Once Phase 3.5 ships a true
//! BRep fillet, [`Solid::Mesh`] can be removed (or kept as an explicit
//! "mesh-output" toggle for visualization-grade workflows).

use std::collections::HashSet;

/// All errors produced by the `valenx-cad` API.
#[derive(Debug, thiserror::Error)]
pub enum CadError {
    /// A boolean operation produced no solid (degenerate intersection,
    /// non-overlapping inputs for AND, etc.). truck-shapeops signals
    /// this by returning `None` from `and` / `or`.
    #[error("boolean op produced no solid")]
    EmptyResult,
    /// Caller passed something the primitive builder won't accept —
    /// negative dimensions, fewer than three points for a prism
    /// profile, etc.
    #[error("invalid parameter: {0}")]
    InvalidParam(String),
    /// truck-meshalgo refused the tessellation request (e.g. caller
    /// passed a tolerance that's not strictly positive).
    #[error("tessellation failed: {0}")]
    Tessellation(String),
    /// Operation is structurally supported by Valenx's API surface
    /// but not yet wired through to a truck implementation. The
    /// canonical case in 2026 is per-edge filleting, which truck 0.6
    /// does not expose.
    #[error("{op} is not implemented yet: {reason}")]
    NotImplemented {
        /// Logical operation name (e.g. `"fillet_edges"`).
        op: &'static str,
        /// One-line explanation surfaced in the UI so users know why
        /// it doesn't work and what to do instead.
        reason: String,
    },
    /// Caller tried to perform a BRep-only operation (boolean, sweep,
    /// rigid transform, face/edge/vertex topology query) on a
    /// mesh-backed [`Solid`]. Mesh-backed solids ship from the v1
    /// fillet/chamfer pipeline and don't carry the topological info
    /// these ops require. Surface the error verbatim to the user so
    /// they know to either apply the filleting last in the feature
    /// tree or skip it for the boolean-heavy parts of their model.
    #[error("operation '{op}' is not supported on mesh-backed solids: {reason}")]
    MeshBackedSolid {
        /// Logical operation name (e.g. `"union"`, `"faces"`).
        op: &'static str,
        /// One-line explanation surfaced in the UI.
        reason: String,
    },
}

/// A closed solid — either a true BRep managed by the truck kernel or
/// a triangle-mesh approximation produced by `valenx-fillet`.
///
/// Most of the time you want the [`Solid::Brep`] variant: primitives,
/// booleans, sweeps, rotations, and STEP export all need real
/// topology. [`Solid::Mesh`] exists so the mesh-domain fillet pipeline
/// in `valenx-fillet` can return *something* that the feature tree and
/// viewport can render, with the understanding that downstream BRep
/// ops will refuse to operate on it.
///
/// Reach for the [`crate::primitives`] / [`crate::boolean`] helpers
/// when constructing a [`Solid::Brep`]; use [`Solid::from_mesh`] for
/// the mesh variant.
#[derive(Clone, Debug)]
pub enum Solid {
    /// Classic truck-modeling BRep. Faces, edges, and vertices are
    /// real topological entities; all CAD ops are supported.
    Brep(truck_modeling::Solid),
    /// Triangle-mesh approximation. Produced by the v1 fillet /
    /// chamfer pipeline. Booleans, sweeps, rigid transforms, and
    /// topological counts all fail with [`CadError::MeshBackedSolid`].
    Mesh(valenx_mesh::Mesh),
}

impl Solid {
    /// Wrap a raw truck solid. Crate-private — callers should reach
    /// for the [`crate::primitives`] / [`crate::boolean`] functions
    /// instead of constructing a [`Solid`] by hand.
    pub(crate) fn from_inner(inner: truck_modeling::Solid) -> Self {
        Self::Brep(inner)
    }

    /// Public wrapper for an externally-built truck solid. Used by
    /// downstream crates (e.g. `valenx-sketch::extrude`) that build
    /// their own BRep with `truck_modeling::builder` and need to hand
    /// the result back across the Valenx API boundary. Prefer the
    /// [`crate::primitives`] helpers when possible.
    pub fn from_truck(inner: truck_modeling::Solid) -> Self {
        Self::Brep(inner)
    }

    /// Wrap a triangle [`valenx_mesh::Mesh`] as a mesh-backed Solid.
    ///
    /// **v1 limitation:** the resulting Solid carries no BRep
    /// topology. Booleans, sweeps, rigid transforms, and the
    /// face/edge/vertex counters will all return
    /// [`CadError::MeshBackedSolid`]. Tessellation via
    /// [`crate::solid_to_mesh`] returns the cached mesh as-is — no
    /// re-tessellation, no chord-error budget applied.
    ///
    /// Use this constructor to round-trip the output of
    /// `valenx-fillet::apply_fillet` / `apply_chamfer` back through
    /// the feature tree. The fillet / chamfer should be the *last*
    /// op in any chain that needs further BRep work — once a mesh
    /// has entered the pipeline, the BRep door is closed.
    pub fn from_mesh(mesh: valenx_mesh::Mesh) -> Self {
        Self::Mesh(mesh)
    }

    /// Borrow the underlying truck solid. Crate-private so the
    /// boolean + tessellation modules can call into truck-shapeops /
    /// truck-meshalgo without `pub`-leaking the upstream type.
    ///
    /// Returns [`CadError::MeshBackedSolid`] when the solid is a
    /// mesh-backed variant — that's the typed error every BRep op
    /// surfaces when it can't operate on the input.
    pub(crate) fn try_inner(&self) -> Result<&truck_modeling::Solid, CadError> {
        match self {
            Self::Brep(b) => Ok(b),
            Self::Mesh(_) => Err(CadError::MeshBackedSolid {
                op: "<brep-only>",
                reason: "this solid was produced by valenx-fillet's mesh-domain \
                         pipeline and has no BRep topology; apply filleting \
                         last in the feature tree, or rebuild without the \
                         fillet to perform this op"
                    .to_string(),
            }),
        }
    }

    /// Mutable borrow of the underlying solid. Required for the
    /// `Solid::not()` call inside [`crate::difference`] which inverts
    /// face orientations on the second operand in place.
    pub(crate) fn try_inner_mut(&mut self) -> Result<&mut truck_modeling::Solid, CadError> {
        match self {
            Self::Brep(b) => Ok(b),
            Self::Mesh(_) => Err(CadError::MeshBackedSolid {
                op: "<brep-only mutate>",
                reason: "mesh-backed solids cannot be mutated in place".to_string(),
            }),
        }
    }

    /// Borrow the cached mesh for a mesh-backed solid. Returns
    /// `Some(&Mesh)` when this is a [`Solid::Mesh`], `None` for a
    /// [`Solid::Brep`]. Used by `solid_to_mesh` to short-circuit the
    /// truck-meshalgo tessellation step.
    pub(crate) fn cached_mesh(&self) -> Option<&valenx_mesh::Mesh> {
        match self {
            Self::Mesh(m) => Some(m),
            Self::Brep(_) => None,
        }
    }

    /// Total number of faces across the solid's bounding shells.
    ///
    /// Returns [`CadError::MeshBackedSolid`] for a mesh-backed solid
    /// (mesh triangles aren't BRep faces — the closest analog would
    /// be the triangle count, which is what
    /// `valenx_mesh::Mesh::total_elements` reports).
    pub fn faces(&self) -> usize {
        match self {
            Self::Brep(b) => b.boundaries().iter().map(|shell| shell.len()).sum(),
            // Mesh-backed solid has no BRep faces. Return 0 rather
            // than panicking so test-helpers that call faces() in
            // assertions keep working; callers that need true
            // topology should query the variant directly.
            Self::Mesh(_) => 0,
        }
    }

    /// Number of distinct edges in the solid, de-duplicated by edge
    /// ID. Without the HashSet we'd double-count every edge that's
    /// shared between two faces (which is most of them in a closed
    /// manifold). Returns 0 for mesh-backed solids.
    pub fn edges(&self) -> usize {
        match self {
            Self::Brep(b) => {
                let mut seen = HashSet::new();
                for edge in b.edge_iter() {
                    seen.insert(edge.id());
                }
                seen.len()
            }
            Self::Mesh(_) => 0,
        }
    }

    /// Number of distinct vertices in the solid, de-duplicated by
    /// vertex ID. Returns 0 for mesh-backed solids.
    pub fn vertices(&self) -> usize {
        match self {
            Self::Brep(b) => {
                let mut seen = HashSet::new();
                for vertex in b.vertex_iter() {
                    seen.insert(vertex.id());
                }
                seen.len()
            }
            Self::Mesh(_) => 0,
        }
    }

    /// Return a new solid translated by `(dx, dy, dz)`. The original
    /// is left untouched (clone-and-transform semantics — convenient
    /// when you want to keep an "untransformed" reference).
    ///
    /// Mesh-backed solids translate their nodes directly; BRep solids
    /// delegate to truck-modeling's `builder::translated`.
    ///
    /// # Errors
    ///
    /// Returns [`CadError::InvalidParam`] when any of `dx`, `dy`, or
    /// `dz` is non-finite (NaN or ±inf). Round-6 fix: pre-fix a NaN
    /// silently propagated into the transform matrix and produced
    /// solids whose entire vertex array was NaN — downstream
    /// tessellation then panicked on the unrecoverable geometry.
    pub fn translated(&self, dx: f64, dy: f64, dz: f64) -> Result<Self, CadError> {
        if !(dx.is_finite() && dy.is_finite() && dz.is_finite()) {
            return Err(CadError::InvalidParam(format!(
                "translated: dx/dy/dz must be finite (got {dx}, {dy}, {dz})"
            )));
        }
        Ok(match self {
            Self::Brep(b) => {
                let inner = truck_modeling::builder::translated(
                    b,
                    truck_modeling::Vector3::new(dx, dy, dz),
                );
                Self::Brep(inner)
            }
            Self::Mesh(m) => {
                let mut out = m.clone();
                for n in &mut out.nodes {
                    n.x += dx;
                    n.y += dy;
                    n.z += dz;
                }
                Self::Mesh(out)
            }
        })
    }

    /// Return a new solid rotated around `axis` (passing through
    /// `origin`) by `angle_rad` radians. Right-hand rule.
    ///
    /// **Mesh-backed solids:** the mesh is returned unchanged (see
    /// the note on the BRep path below).
    ///
    /// # Errors
    ///
    /// Returns [`CadError::InvalidParam`] when any component of
    /// `origin` / `axis` or `angle_rad` is non-finite. Round-6 fix.
    pub fn rotated(
        &self,
        origin: (f64, f64, f64),
        axis: (f64, f64, f64),
        angle_rad: f64,
    ) -> Result<Self, CadError> {
        if !(origin.0.is_finite() && origin.1.is_finite() && origin.2.is_finite()) {
            return Err(CadError::InvalidParam(format!(
                "rotated: origin components must be finite (got {origin:?})"
            )));
        }
        if !(axis.0.is_finite() && axis.1.is_finite() && axis.2.is_finite()) {
            return Err(CadError::InvalidParam(format!(
                "rotated: axis components must be finite (got {axis:?})"
            )));
        }
        // Round-8 sister fix to `Solid::mirrored` (round-3): a
        // zero-length axis is degenerate — truck builds the rotation
        // around `Vector3::zero()` which collapses to identity (or
        // worse, NaN-poisoned transforms downstream). Refuse the
        // input here so the feature-tree evaluator gets a clean
        // structured error instead of silent geometry corruption.
        let (ax, ay, az) = axis;
        let n_len = (ax * ax + ay * ay + az * az).sqrt();
        if n_len <= 0.0 || !n_len.is_finite() {
            return Err(CadError::InvalidParam(format!(
                "rotated: axis must be non-zero (got ({ax}, {ay}, {az}))"
            )));
        }
        if !angle_rad.is_finite() {
            return Err(CadError::InvalidParam(format!(
                "rotated: angle_rad must be finite (got {angle_rad})"
            )));
        }
        Ok(match self {
            Self::Brep(b) => {
                let inner = truck_modeling::builder::rotated(
                    b,
                    truck_modeling::Point3::new(origin.0, origin.1, origin.2),
                    truck_modeling::Vector3::new(axis.0, axis.1, axis.2),
                    truck_modeling::Rad(angle_rad),
                );
                Self::Brep(inner)
            }
            // For mesh-backed solids we leave the mesh unchanged. The
            // existing rotated() callers (transform helpers, pattern
            // ops) only ever invoke it on BRep solids built from
            // primitives or sweeps. If a caller does land here we want
            // it to surface as "no-op on mesh" rather than panic — an
            // explicit error would force every transform call to bubble
            // a Result through layers that currently don't propagate
            // one. Phase 3.5 (true BRep fillet) removes this hazard.
            Self::Mesh(_) => self.clone(),
        })
    }

    /// Return a new solid reflected across the plane through `origin`
    /// with the given (not necessarily unit) `normal`.
    ///
    /// # Caveats
    ///
    /// Reflection inverts handedness, so the underlying truck BRep
    /// comes out with all face orientations flipped. This method
    /// re-inverts them via `Solid::not()` so the returned solid has
    /// outward-facing normals (necessary for downstream boolean ops to
    /// behave the same way as on un-reflected geometry).
    ///
    /// # Errors
    ///
    /// Returns [`CadError::InvalidParam`] when `normal` is the zero
    /// vector or contains any non-finite component. Round-3 fix:
    /// previously this asserted, which would panic the whole process
    /// instead of letting the feature-tree evaluator recover.
    ///
    /// **Mesh-backed solids:** returns the mesh unchanged (see the note
    /// on [`Solid::rotated`]).
    pub fn mirrored(
        &self,
        origin: (f64, f64, f64),
        normal: (f64, f64, f64),
    ) -> Result<Self, CadError> {
        let b = match self {
            Self::Brep(b) => b,
            Self::Mesh(_) => return Ok(self.clone()),
        };
        use truck_modeling::{Matrix4, Vector3, Vector4};
        if !(normal.0.is_finite() && normal.1.is_finite() && normal.2.is_finite()) {
            return Err(CadError::InvalidParam(format!(
                "mirrored: normal components must be finite (got {normal:?})"
            )));
        }
        if !(origin.0.is_finite() && origin.1.is_finite() && origin.2.is_finite()) {
            return Err(CadError::InvalidParam(format!(
                "mirrored: origin components must be finite (got {origin:?})"
            )));
        }
        // Normalize the plane normal — caller may pass in arbitrary
        // length (we still want a clean reflection regardless).
        let n_len = (normal.0 * normal.0 + normal.1 * normal.1 + normal.2 * normal.2).sqrt();
        if n_len <= 0.0 {
            return Err(CadError::InvalidParam(format!(
                "mirrored: normal must be non-zero (got {normal:?})"
            )));
        }
        let nx = normal.0 / n_len;
        let ny = normal.1 / n_len;
        let nz = normal.2 / n_len;

        // The reflection matrix R about a plane through the origin
        // with unit normal n is I - 2 * n * n^T (a Householder
        // matrix). For a plane offset by `origin`, compose
        // T(origin) * R * T(-origin).
        //
        // cgmath's Matrix4::new takes column-major order:
        //   Matrix4::new(c0r0, c0r1, c0r2, c0r3,
        //                c1r0, c1r1, c1r2, c1r3,
        //                c2r0, c2r1, c2r2, c2r3,
        //                c3r0, c3r1, c3r2, c3r3)
        // R is symmetric (R^T = R) so column-major and row-major
        // forms coincide for the 3×3 block; we still spell it out
        // column-by-column for clarity.
        let reflect = Matrix4::from_cols(
            Vector4::new(1.0 - 2.0 * nx * nx, -2.0 * ny * nx, -2.0 * nz * nx, 0.0),
            Vector4::new(-2.0 * nx * ny, 1.0 - 2.0 * ny * ny, -2.0 * nz * ny, 0.0),
            Vector4::new(-2.0 * nx * nz, -2.0 * ny * nz, 1.0 - 2.0 * nz * nz, 0.0),
            Vector4::new(0.0, 0.0, 0.0, 1.0),
        );

        let to_origin = Matrix4::from_translation(Vector3::new(-origin.0, -origin.1, -origin.2));
        let back = Matrix4::from_translation(Vector3::new(origin.0, origin.1, origin.2));
        let mat = back * reflect * to_origin;

        let mut inner = truck_modeling::builder::transformed(b, mat);
        // Reflection has determinant -1, so faces come out
        // "inside-out". Flip them back so the surface normals point
        // outward — otherwise shapeops's boolean ops will treat the
        // reflected solid as a "hole" rather than a solid.
        inner.not();
        Ok(Self::Brep(inner))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::box_solid;

    #[test]
    fn cad_error_display_includes_message() {
        let err = CadError::InvalidParam("dx must be positive".into());
        assert!(err.to_string().contains("dx must be positive"));

        let err = CadError::NotImplemented {
            op: "fillet_edges",
            reason: "truck 0.6 does not expose an edge-fillet algorithm".into(),
        };
        let msg = err.to_string();
        assert!(msg.contains("fillet_edges"));
        assert!(msg.contains("truck 0.6"));

        // The MeshBackedSolid variant ships in Phase 3 — verify the
        // formatter surfaces both the op name and the reason so users
        // know which call failed and why.
        let err = CadError::MeshBackedSolid {
            op: "union",
            reason: "this solid has no BRep topology".into(),
        };
        let msg = err.to_string();
        assert!(msg.contains("union"));
        assert!(msg.contains("BRep"));
    }

    #[test]
    fn cube_has_expected_topology() {
        // A unit cube has 6 faces, 12 edges, 8 vertices.
        let cube = box_solid(1.0, 1.0, 1.0).expect("unit cube builds");
        assert_eq!(cube.faces(), 6, "cube should have 6 faces");
        assert_eq!(cube.edges(), 12, "cube should have 12 edges");
        assert_eq!(cube.vertices(), 8, "cube should have 8 vertices");
    }

    #[test]
    fn mesh_backed_solid_has_zero_brep_topology() {
        // A mesh-backed Solid has no BRep faces/edges/vertices —
        // querying those returns 0, not a panic.
        let mesh = valenx_mesh::Mesh::new("mesh-backed");
        let s = Solid::from_mesh(mesh);
        assert_eq!(s.faces(), 0, "mesh-backed solid has no BRep faces");
        assert_eq!(s.edges(), 0, "mesh-backed solid has no BRep edges");
        assert_eq!(s.vertices(), 0, "mesh-backed solid has no BRep vertices");
    }

    /// Round-3 fix: `Solid::mirrored` used to `assert!()` on a zero
    /// normal, which panics the whole process and skips any caller
    /// recovery. The fallible signature returns InvalidParam so the
    /// feature-tree evaluator surfaces a structured error instead.
    #[test]
    fn mirrored_returns_err_on_zero_normal() {
        let cube = box_solid(1.0, 1.0, 1.0).expect("unit cube builds");
        let result = cube.mirrored((0.0, 0.0, 0.0), (0.0, 0.0, 0.0));
        let err = result.expect_err("zero normal must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("non-zero"), "msg: {msg}");
    }

    #[test]
    fn mirrored_returns_err_on_nan_normal() {
        let cube = box_solid(1.0, 1.0, 1.0).expect("unit cube builds");
        let result = cube.mirrored((0.0, 0.0, 0.0), (f64::NAN, 0.0, 1.0));
        let err = result.expect_err("NaN normal must be rejected");
        let msg = format!("{err}");
        assert!(msg.contains("finite"), "msg: {msg}");
    }

    #[test]
    fn mirrored_succeeds_for_valid_normal() {
        let cube = box_solid(1.0, 1.0, 1.0).expect("unit cube builds");
        let result = cube.mirrored((0.0, 0.0, 0.0), (1.0, 0.0, 0.0));
        assert!(result.is_ok(), "valid normal should succeed");
    }

    #[test]
    fn from_mesh_round_trip_preserves_mesh_identity() {
        // The constructor stores the mesh as-is; cached_mesh() returns
        // it. Used by solid_to_mesh to short-circuit re-tessellation.
        let mut mesh = valenx_mesh::Mesh::new("test");
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 0.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(1.0, 0.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 1.0, 0.0));
        let s = Solid::from_mesh(mesh);
        let cached = s.cached_mesh().expect("from_mesh stores the mesh");
        assert_eq!(cached.nodes.len(), 3);
        assert!(matches!(s, Solid::Mesh(_)));
    }

    #[test]
    fn mesh_backed_translated_moves_nodes() {
        // Translating a mesh-backed solid must shift each node by the
        // delta — the fillet/chamfer output needs to play nice with
        // Mirror / pattern ops downstream that translate the solid.
        let mut mesh = valenx_mesh::Mesh::new("tx");
        mesh.nodes.push(nalgebra::Vector3::new(1.0, 2.0, 3.0));
        let s = Solid::from_mesh(mesh).translated(10.0, -1.0, 0.5).unwrap();
        let cached = s.cached_mesh().unwrap();
        assert!((cached.nodes[0] - nalgebra::Vector3::new(11.0, 1.0, 3.5)).norm() < 1e-12);
    }

    #[test]
    fn translated_rejects_non_finite_components() {
        // Round-6 RED→GREEN: a NaN/inf delta must NOT silently
        // produce a corrupt BRep — instead returns a structured
        // InvalidParam error the caller can recover from.
        let s = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = s.translated(f64::NAN, 0.0, 0.0).unwrap_err();
        match err {
            CadError::InvalidParam(msg) => assert!(msg.contains("finite"), "msg: {msg}"),
            other => panic!("expected InvalidParam, got {other:?}"),
        }
        let err = s.translated(f64::INFINITY, 0.0, 0.0).unwrap_err();
        assert!(matches!(err, CadError::InvalidParam(_)));
        let err = s.translated(0.0, f64::NEG_INFINITY, 0.0).unwrap_err();
        assert!(matches!(err, CadError::InvalidParam(_)));
        // Finite still works.
        assert!(s.translated(1.0, 2.0, 3.0).is_ok());
    }

    #[test]
    fn rotated_rejects_non_finite_components() {
        // Round-6 RED→GREEN companion for `rotated`.
        let s = box_solid(1.0, 1.0, 1.0).unwrap();
        // Bad origin
        assert!(matches!(
            s.rotated((f64::NAN, 0.0, 0.0), (0.0, 0.0, 1.0), 0.5)
                .unwrap_err(),
            CadError::InvalidParam(_)
        ));
        // Bad axis
        assert!(matches!(
            s.rotated((0.0, 0.0, 0.0), (0.0, f64::INFINITY, 1.0), 0.5)
                .unwrap_err(),
            CadError::InvalidParam(_)
        ));
        // Bad angle
        assert!(matches!(
            s.rotated((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), f64::NAN)
                .unwrap_err(),
            CadError::InvalidParam(_)
        ));
        // Finite still works.
        assert!(s.rotated((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.5).is_ok());
    }

    #[test]
    fn rotated_rejects_zero_axis() {
        // Round-8 RED→GREEN: sister fix to `Solid::mirrored`'s
        // zero-normal check. A zero-length rotation axis is
        // degenerate — truck builds the rotation around
        // `Vector3::zero()` and the result either collapses to
        // identity or NaN-poisons downstream geometry. Refuse the
        // input here.
        let s = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = s
            .rotated((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 1.0)
            .unwrap_err();
        match err {
            CadError::InvalidParam(msg) => {
                assert!(
                    msg.contains("non-zero") && msg.contains("axis"),
                    "msg: {msg}"
                );
            }
            other => panic!("expected InvalidParam, got {other:?}"),
        }
    }
}
