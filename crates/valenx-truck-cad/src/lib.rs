//! `valenx-truck-cad` — a focused **B-Rep solid-modeling V1** built
//! directly on the [`truck`](https://crates.io/crates/truck-modeling)
//! Rust CAD kernel (RICOS Co. Ltd, Apache-2.0).
//!
//! This crate IS Valenx's `truck` integration. It exposes the three
//! things a real boundary-representation modeller needs, behind one
//! small AI-drivable facade ([`BrepKernel`]):
//!
//! 1. **Primitives** — [`BrepKernel::primitive`] builds a `box`,
//!    `cylinder`, or `sphere` as a genuine closed [`Solid`] BRep (faces
//!    / edges / vertices forming a closed 2-manifold), via truck's
//!    `builder` sweeps — *not* a triangle soup. The curved surfaces are
//!    truck's NURBS surfaces of revolution, so the solid carries exact
//!    geometry that boolean ops and tessellation both consume.
//! 2. **Booleans** — [`BrepKernel::boolean`] performs `union`,
//!    `difference`, or `intersection` of two solids (`truck-shapeops`
//!    `or` / `and`, with difference = `A ∩ ¬B`).
//! 3. **Tessellation** — [`BrepKernel::tessellate`] meshes a solid into
//!    a flat [`TriMesh`] (positions + triangle indices + an axis-aligned
//!    bounding box) suitable for the egui viewport / STL export.
//!
//! Design note
//! ===========
//!
//! The primitive builders and the boolean operators are delegated to
//! the in-house [`valenx_cad`] kernel, which itself sits on `truck`.
//! That crate already contains the hardening this V1 would otherwise
//! have to duplicate: `truck-shapeops` can `panic!` on degenerate input
//! (a self-intersecting intermediate wire) and can return a *phantom*
//! shell-less "solid" for a disjoint difference. `valenx-cad` contains
//! both failure modes (a `catch_unwind` guard + a degenerate-result
//! filter), so a boolean here either yields a genuine non-empty solid or
//! a typed [`BrepError`] — it never unwinds the caller and never returns
//! a silently-invalid solid. Building the V1 facade on top of that,
//! rather than re-implementing the guards, is the honest CRATE-FIRST
//! choice.
//!
//! The crate re-exports the concrete truck [`Solid`] so callers never
//! have to spell out the three-parameter
//! `truck_topology::Solid<Point3, Curve, Surface>`.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use thiserror::Error;

/// The concrete truck B-Rep solid type — `Solid<Point3, Curve, Surface>`
/// — re-exported so downstream code (and this crate's own API) can name
/// it without the three generic parameters. A [`Solid`] is a closed
/// boundary representation: a set of oriented shells of NURBS faces.
pub use truck_modeling::Solid;

/// Default chord-error tolerance for [`BrepKernel::tessellate`] (model
/// units). Mirrors [`valenx_cad::DEFAULT_TESS_TOLERANCE`].
pub const DEFAULT_TESS_TOLERANCE: f64 = valenx_cad::DEFAULT_TESS_TOLERANCE;

/// Default linear tolerance for the boolean set-ops (model units).
/// Mirrors [`valenx_cad::DEFAULT_BOOL_TOLERANCE`].
pub const DEFAULT_BOOL_TOLERANCE: f64 = valenx_cad::DEFAULT_BOOL_TOLERANCE;

/// Errors from the B-Rep kernel facade.
///
/// Wraps the underlying [`valenx_cad::CadError`] with a thin,
/// crate-local variant set so callers depend on this crate's error type
/// rather than the kernel's. The [`From`] impl preserves the original
/// message.
#[derive(Debug, Error)]
pub enum BrepError {
    /// An input parameter was invalid (non-positive size, non-finite
    /// value, a self-intersecting torus, …).
    #[error("invalid parameter: {0}")]
    InvalidParam(String),
    /// A boolean produced no real solid — an empty intersection, a
    /// disjoint difference, or a degenerate input that tripped (and was
    /// contained from) a `truck-shapeops` panic.
    #[error("boolean produced an empty / degenerate result")]
    EmptyResult,
    /// Tessellation failed (bad tolerance, or a degenerate solid that
    /// meshed to zero triangles).
    #[error("tessellation failed: {0}")]
    Tessellation(String),
    /// Any other kernel error, carried through verbatim.
    #[error("{0}")]
    Kernel(String),
}

impl From<valenx_cad::CadError> for BrepError {
    fn from(e: valenx_cad::CadError) -> Self {
        use valenx_cad::CadError as C;
        match e {
            C::InvalidParam(m) => BrepError::InvalidParam(m),
            C::EmptyResult => BrepError::EmptyResult,
            C::Tessellation(m) => BrepError::Tessellation(m),
            other => BrepError::Kernel(other.to_string()),
        }
    }
}

/// A primitive solid the kernel can build. Sizes are in model units and
/// must be strictly positive; the builders reject anything else.
///
/// Each variant maps onto a validated [`valenx_cad`] primitive built
/// from truck `builder` sweeps, so the result is always a proper closed
/// BRep. Coordinate conventions match `valenx-cad`:
/// box corner at the origin, cylinder/sphere centred on the origin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Primitive {
    /// Axis-aligned box: one corner at the origin, opposite corner at
    /// `(dx, dy, dz)`.
    Box {
        /// Extent along X.
        dx: f64,
        /// Extent along Y.
        dy: f64,
        /// Extent along Z.
        dz: f64,
    },
    /// Right circular cylinder centred on the origin, axis along +Z.
    Cylinder {
        /// Base-disk radius.
        radius: f64,
        /// Height along +Z.
        height: f64,
    },
    /// Sphere centred on the origin.
    Sphere {
        /// Sphere radius.
        radius: f64,
    },
}

impl Primitive {
    /// Stable lowercase id for this primitive's *shape*
    /// (`"box"`/`"cylinder"`/`"sphere"`), independent of its sizes —
    /// handy for AI-drivable selectors and readouts.
    pub fn kind_id(&self) -> &'static str {
        match self {
            Primitive::Box { .. } => "box",
            Primitive::Cylinder { .. } => "cylinder",
            Primitive::Sphere { .. } => "sphere",
        }
    }
}

/// A boolean set operation between two solids.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoolOp {
    /// `A ∪ B` — weld both solids into one.
    Union,
    /// `A − B` — carve `B` out of `A` (`A ∩ ¬B`).
    Difference,
    /// `A ∩ B` — keep only the overlap.
    Intersection,
}

impl BoolOp {
    /// Stable lowercase id (`"union"`/`"difference"`/`"intersection"`).
    pub fn id(&self) -> &'static str {
        match self {
            BoolOp::Union => "union",
            BoolOp::Difference => "difference",
            BoolOp::Intersection => "intersection",
        }
    }

    /// Parse a (case-insensitive, whitespace-tolerant) id / common alias
    /// into a [`BoolOp`]. Accepts `or`/`and`/`sub`/`subtract`/`cut`/
    /// `minus`/`diff`/`intersect` as friendly synonyms. `None` for an
    /// unknown string.
    pub fn from_id(s: &str) -> Option<BoolOp> {
        match s.trim().to_ascii_lowercase().as_str() {
            "union" | "or" | "add" | "fuse" | "weld" => Some(BoolOp::Union),
            "difference" | "diff" | "sub" | "subtract" | "cut" | "minus" => {
                Some(BoolOp::Difference)
            }
            "intersection" | "intersect" | "and" | "common" => Some(BoolOp::Intersection),
            _ => None,
        }
    }
}

/// A flat triangle mesh — the tessellation output the viewport / STL
/// export consume.
///
/// `positions` is a flat list of `[x, y, z]` vertices; `triangles` is a
/// flat list of vertex indices, three per triangle (so
/// `triangles.len()` is always a multiple of 3). `bounds` is the
/// axis-aligned bounding box `(min, max)` of the positions.
///
/// NB: truck tessellates each face with its own vertex array, so
/// vertices on shared edges appear more than once. That is fine for
/// rendering and STL; weld downstream if you need a manifold index.
#[derive(Clone, Debug, Default)]
pub struct TriMesh {
    /// Flat `[x, y, z]` vertex positions.
    pub positions: Vec<[f64; 3]>,
    /// Flat triangle vertex indices (3 per triangle).
    pub triangles: Vec<u32>,
    /// Axis-aligned bounding box `(min, max)`, each `[x, y, z]`. The
    /// zero box for an empty mesh.
    pub bounds: ([f64; 3], [f64; 3]),
}

impl TriMesh {
    /// Number of triangles (`triangles.len() / 3`).
    pub fn triangle_count(&self) -> usize {
        self.triangles.len() / 3
    }

    /// Number of (un-welded) vertices.
    pub fn vertex_count(&self) -> usize {
        self.positions.len()
    }

    /// Per-axis size of the bounding box `(max − min)`.
    pub fn bound_extents(&self) -> [f64; 3] {
        let (min, max) = self.bounds;
        [max[0] - min[0], max[1] - min[1], max[2] - min[2]]
    }
}

/// The B-Rep modeling facade: build primitives, combine them with
/// boolean ops, and tessellate the result.
///
/// Stateless — every method is a pure function of its inputs. Held as a
/// unit struct so the API reads as a kernel object (`BrepKernel::new()
/// .primitive(...)`) and so future state (a tolerance default, a unit
/// system) has a home without an API break.
#[derive(Clone, Copy, Debug, Default)]
pub struct BrepKernel;

impl BrepKernel {
    /// A fresh kernel handle.
    pub fn new() -> Self {
        BrepKernel
    }

    /// Build a [`Primitive`] into a closed truck [`Solid`] BRep.
    ///
    /// # Errors
    /// [`BrepError::InvalidParam`] if any size is non-positive or
    /// non-finite, or the primitive is otherwise geometrically invalid.
    pub fn primitive(&self, prim: Primitive) -> Result<Solid, BrepError> {
        let solid = match prim {
            Primitive::Box { dx, dy, dz } => valenx_cad::primitives::box_solid(dx, dy, dz)?,
            Primitive::Cylinder { radius, height } => {
                valenx_cad::primitives::cylinder(radius, height)?
            }
            Primitive::Sphere { radius } => valenx_cad::primitives::sphere(radius)?,
        };
        // valenx-cad's primitive builders always return a BRep-backed
        // Solid; unwrap the inner truck solid for this crate's API.
        into_truck(solid)
    }

    /// Combine two solids with a [`BoolOp`] at [`DEFAULT_BOOL_TOLERANCE`].
    ///
    /// # Errors
    /// [`BrepError::EmptyResult`] for an empty / disjoint result (and
    /// for any degenerate input that tripped — and was contained from —
    /// a `truck-shapeops` panic). [`BrepError::InvalidParam`] for a bad
    /// tolerance (never, here — the default is valid).
    pub fn boolean(&self, op: BoolOp, a: &Solid, b: &Solid) -> Result<Solid, BrepError> {
        self.boolean_tol(op, a, b, DEFAULT_BOOL_TOLERANCE)
    }

    /// Combine two solids with an explicit linear tolerance.
    ///
    /// # Errors
    /// As [`BrepKernel::boolean`], plus [`BrepError::InvalidParam`] if
    /// `tol` is not finite and strictly positive.
    pub fn boolean_tol(
        &self,
        op: BoolOp,
        a: &Solid,
        b: &Solid,
        tol: f64,
    ) -> Result<Solid, BrepError> {
        // Wrap the borrowed truck solids in valenx-cad's Solid enum
        // (a cheap clone — truck topology is Rc-shared internally) so we
        // can call the hardened, panic-guarded boolean wrappers.
        let av = valenx_cad::Solid::from_truck(a.clone());
        let bv = valenx_cad::Solid::from_truck(b.clone());
        let out = match op {
            BoolOp::Union => valenx_cad::boolean::union_tol(&av, &bv, tol),
            BoolOp::Difference => valenx_cad::boolean::difference_tol(&av, &bv, tol),
            BoolOp::Intersection => valenx_cad::boolean::intersection_tol(&av, &bv, tol),
        }?;
        into_truck(out)
    }

    /// Tessellate a solid into a flat [`TriMesh`] at
    /// [`DEFAULT_TESS_TOLERANCE`].
    ///
    /// # Errors
    /// [`BrepError::Tessellation`] if the solid meshes to zero triangles.
    pub fn tessellate(&self, solid: &Solid) -> Result<TriMesh, BrepError> {
        self.tessellate_tol(solid, DEFAULT_TESS_TOLERANCE)
    }

    /// Tessellate a solid into a flat [`TriMesh`] at an explicit
    /// chord-error tolerance (smaller = denser).
    ///
    /// # Errors
    /// [`BrepError::Tessellation`] if `tol` is not finite and strictly
    /// positive, or the solid meshes to zero triangles.
    pub fn tessellate_tol(&self, solid: &Solid, tol: f64) -> Result<TriMesh, BrepError> {
        // Reuse valenx-cad's BRep→Mesh path (truck-meshalgo
        // constrained-Delaunay + quad/n-gon splitting), then flatten its
        // canonical `valenx_mesh::Mesh` into the simple TriMesh this
        // crate exposes.
        let cad_solid = valenx_cad::Solid::from_truck(solid.clone());
        let mesh = valenx_cad::solid_to_mesh(&cad_solid, tol)?;

        let mut positions = Vec::with_capacity(mesh.nodes.len());
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for n in &mesh.nodes {
            let p = [n.x, n.y, n.z];
            for k in 0..3 {
                min[k] = min[k].min(p[k]);
                max[k] = max[k].max(p[k]);
            }
            positions.push(p);
        }

        // Flatten every element block's connectivity. valenx-cad only
        // ever emits Tri3 from a solid tessellation, but stride-guard so
        // a stray non-triangle block can't desync the index stream.
        let mut triangles = Vec::new();
        for block in &mesh.element_blocks {
            let stride = block.element_type.nodes_per_element();
            if stride != 3 {
                continue;
            }
            triangles.extend(block.connectivity.iter().copied());
        }

        if triangles.is_empty() {
            return Err(BrepError::Tessellation(
                "tessellation produced zero triangles — the solid may be degenerate".into(),
            ));
        }
        if positions.is_empty() {
            // Defensive: triangles without positions is malformed.
            min = [0.0; 3];
            max = [0.0; 3];
        }

        Ok(TriMesh {
            positions,
            triangles,
            bounds: (min, max),
        })
    }
}

// ---------------------------------------------------------------------------
// Re-exported measurement helpers (validation ground truth)
// ---------------------------------------------------------------------------

/// Signed volume of a solid (cubic model units), at the kernel's default
/// measurement tolerance. Positive for a correctly-oriented closed
/// solid. Thin re-export of [`valenx_cad::solid_volume`] — the
/// validation suite checks constructed solids against analytic ground
/// truth with it.
///
/// # Errors
/// [`BrepError::Tessellation`] if the solid cannot be measured.
pub fn solid_volume(solid: &Solid) -> Result<f64, BrepError> {
    let cad = valenx_cad::Solid::from_truck(solid.clone());
    Ok(valenx_cad::measure::solid_volume(&cad)?)
}

/// Whether a solid's boundary is a closed orientable 2-manifold
/// (watertight). Thin re-export of [`valenx_cad::measure::is_closed_solid`].
///
/// # Errors
/// [`BrepError::Tessellation`] if the solid cannot be tessellated for
/// the check.
pub fn is_closed_solid(solid: &Solid) -> Result<bool, BrepError> {
    let cad = valenx_cad::Solid::from_truck(solid.clone());
    Ok(valenx_cad::measure::is_closed_solid(&cad)?)
}

/// Topology counts `(faces, edges, vertices)` of a solid's boundary —
/// the BRep cell counts (not the tessellation's). Useful for an
/// Euler-characteristic sanity check (`V − E + F = 2` for a solid
/// topologically equivalent to a sphere).
pub fn topology_counts(solid: &Solid) -> (usize, usize, usize) {
    let cad = valenx_cad::Solid::from_truck(solid.clone());
    (cad.faces(), cad.edges(), cad.vertices())
}

/// Unwrap a [`valenx_cad::Solid`] into the inner truck [`Solid`]. The
/// kernel's primitive/boolean builders always return a BRep-backed
/// solid, so the mesh-backed arm is unreachable in this crate; it is
/// mapped to a typed error rather than panicking, for total safety.
///
/// Matches on the public [`valenx_cad::Solid::Brep`] variant (the inner
/// truck solid is a public field) — no crate-private accessor needed.
fn into_truck(s: valenx_cad::Solid) -> Result<Solid, BrepError> {
    match s {
        valenx_cad::Solid::Brep(inner) => Ok(inner),
        valenx_cad::Solid::Mesh(_) => Err(BrepError::Kernel(
            "expected a BRep-backed solid from the kernel, got a mesh-backed one".into(),
        )),
    }
}

// ===========================================================================
// Tests — V1 validation
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const CUBE: Primitive = Primitive::Box {
        dx: 2.0,
        dy: 2.0,
        dz: 2.0,
    };

    fn k() -> BrepKernel {
        BrepKernel::new()
    }

    #[test]
    fn box_primitive_is_a_valid_closed_solid_with_cube_topology() {
        let solid = k().primitive(CUBE).expect("box builds");
        // A box BRep: 6 faces, 12 edges, 8 vertices.
        let (f, e, v) = topology_counts(&solid);
        assert_eq!((f, e, v), (6, 12, 8), "box topology");
        // Euler characteristic V − E + F = 2 (genus-0 solid).
        assert_eq!(v as i64 - e as i64 + f as i64, 2, "Euler V-E+F=2");
        // Watertight.
        assert!(
            is_closed_solid(&solid).expect("closed check runs"),
            "a box must be a closed solid"
        );
        // Volume == 2·2·2 = 8 (flat-faced ⇒ exact).
        let vol = solid_volume(&solid).expect("volume");
        assert!((vol - 8.0).abs() < 1e-6, "box volume {vol} != 8");
    }

    #[test]
    fn cylinder_and_sphere_build_non_empty_solids() {
        let cyl = k()
            .primitive(Primitive::Cylinder {
                radius: 1.0,
                height: 2.0,
            })
            .expect("cylinder builds");
        let (f, _, v) = topology_counts(&cyl);
        assert!(f > 0 && v > 0, "cylinder has faces + vertices");

        let sph = k()
            .primitive(Primitive::Sphere { radius: 1.0 })
            .expect("sphere builds");
        let (sf, _, sv) = topology_counts(&sph);
        assert!(sf > 0 && sv > 0, "sphere has faces + vertices");
    }

    #[test]
    fn primitive_rejects_bad_sizes() {
        assert!(matches!(
            k().primitive(Primitive::Box {
                dx: -1.0,
                dy: 1.0,
                dz: 1.0
            }),
            Err(BrepError::InvalidParam(_))
        ));
        assert!(matches!(
            k().primitive(Primitive::Cylinder {
                radius: 0.0,
                height: 1.0
            }),
            Err(BrepError::InvalidParam(_))
        ));
        assert!(matches!(
            k().primitive(Primitive::Sphere { radius: f64::NAN }),
            Err(BrepError::InvalidParam(_))
        ));
    }

    /// A unit cube ∪ a cylinder straddling its mid-face — the
    /// documented-good `truck-shapeops` input (the "punched cube"
    /// geometry, but OR'd). Avoids the coincident-face case two
    /// axis-aligned boxes overlapping along a shared plane would create,
    /// which `truck-shapeops` is not robust to.
    fn cube_and_overlapping_cylinder() -> (Solid, Solid) {
        let cube = into_truck(valenx_cad::primitives::box_solid(1.0, 1.0, 1.0).expect("cube"))
            .expect("cube inner");
        let cyl = into_truck(
            valenx_cad::primitives::cylinder(0.25, 2.0)
                .expect("cyl")
                .translated(0.5, 0.5, -0.5)
                .expect("translate cyl"),
        )
        .expect("cyl inner");
        (cube, cyl)
    }

    #[test]
    fn union_of_overlapping_cube_and_cylinder_is_a_valid_solid() {
        let (cube, cyl) = cube_and_overlapping_cylinder();
        let u = k()
            .boolean(BoolOp::Union, &cube, &cyl)
            .expect("union builds");
        // The union must be a real (non-empty) solid with faces.
        let (f, _, _) = topology_counts(&u);
        assert!(f > 0, "union has faces");
        // Volume: the cube (1.0) plus the part of the cylinder sticking
        // out above z=1 and below z=0 (each a 0.5-tall plug of the
        // r=0.25 cylinder), so strictly MORE than the cube alone.
        let vol = solid_volume(&u).expect("union volume");
        assert!(
            vol > 1.0,
            "union of cube + protruding cylinder should exceed the cube's volume 1.0, got {vol}"
        );
    }

    #[test]
    fn difference_box_minus_cylinder_punches_a_hole() {
        // The canonical "punched cube": a unit cube minus a cylinder
        // straddling its mid-face. Lifted from truck-shapeops' own
        // integration setup, so shapeops is known to handle it.
        let cube = valenx_cad::primitives::box_solid(1.0, 1.0, 1.0).expect("cube");
        let cube = into_truck(cube).expect("cube inner");
        // Cylinder centred at (0.5, 0.5), radius 0.25, tall enough to
        // pierce the cube (z from −0.5 to 1.5).
        let cyl = valenx_cad::primitives::cylinder(0.25, 2.0)
            .expect("cyl")
            .translated(0.5, 0.5, -0.5)
            .expect("translate cyl");
        let cyl = into_truck(cyl).expect("cyl inner");

        let punched = k()
            .boolean(BoolOp::Difference, &cube, &cyl)
            .expect("difference builds");
        // A punched cube has the cube's 6 outer faces PLUS the inner
        // hole wall — strictly more than 6.
        let (f, _, _) = topology_counts(&punched);
        assert!(
            f > 6,
            "punched cube should have more faces than a plain cube, got {f}"
        );
        // Volume: 1·1·1 − π·0.25²·1 ≈ 1 − 0.19635 = 0.8037 (the cylinder
        // pierces a unit-height cube, so the removed plug is one unit
        // tall). Allow a loose band for boolean/tess discretisation.
        let vol = solid_volume(&punched).expect("punched volume");
        let expected = 1.0 - std::f64::consts::PI * 0.25 * 0.25;
        assert!(
            (vol - expected).abs() < 0.05,
            "punched volume {vol} should be ~{expected}"
        );
    }

    #[test]
    fn intersection_of_overlapping_cube_and_cylinder_is_the_overlap() {
        let (cube, cyl) = cube_and_overlapping_cylinder();
        let inter = k()
            .boolean(BoolOp::Intersection, &cube, &cyl)
            .expect("intersection builds");
        let (f, _, _) = topology_counts(&inter);
        assert!(f > 0, "intersection has faces");
        // The overlap is the cylinder's r=0.25 column clipped to the
        // unit cube's z∈[0,1]: a 1-tall plug, volume ≈ π·0.25²·1 ≈ 0.196.
        let vol = solid_volume(&inter).expect("intersection volume");
        let expected = std::f64::consts::PI * 0.25 * 0.25;
        assert!(
            (vol - expected).abs() < 0.05,
            "intersection volume {vol} should be ~{expected} (the clipped cylinder plug)"
        );
    }

    #[test]
    fn disjoint_difference_surfaces_empty_not_a_phantom() {
        // A − B with B far away must NOT return an Ok phantom solid —
        // valenx-cad's degenerate-result filter converts the shell-less
        // truck result into an empty error, which we map to EmptyResult.
        let a = k().primitive(CUBE).expect("a");
        let b_solid = valenx_cad::primitives::box_solid(2.0, 2.0, 2.0)
            .expect("b")
            .translated(50.0, 50.0, 50.0)
            .expect("translate b");
        let b = into_truck(b_solid).expect("b inner");
        match k().boolean(BoolOp::Difference, &a, &b) {
            Err(BrepError::EmptyResult) => {}
            other => panic!("disjoint A−B should be EmptyResult, got {other:?}"),
        }
    }

    #[test]
    fn tessellate_box_yields_nonempty_mesh_with_sane_bounds() {
        let solid = k().primitive(CUBE).expect("box");
        let mesh = k().tessellate(&solid).expect("tessellation");
        // A box meshes to >= 12 triangles (2 per face).
        assert!(
            mesh.triangle_count() >= 12,
            "box tessellation should have >=12 triangles, got {}",
            mesh.triangle_count()
        );
        assert!(mesh.vertex_count() > 0, "mesh has vertices");
        // triangles flat-list length is a multiple of 3.
        assert_eq!(mesh.triangles.len() % 3, 0, "triangle indices come in 3s");
        // Every index is in range.
        assert!(
            mesh.triangles
                .iter()
                .all(|&i| (i as usize) < mesh.positions.len()),
            "all triangle indices reference a real vertex"
        );
        // Bounds: the 2×2×2 box spans [0,0,0]..[2,2,2] (corner at origin).
        let (min, max) = mesh.bounds;
        for ax in 0..3 {
            assert!(min[ax] <= 1e-9, "min[{ax}]={} should be ~0", min[ax]);
            assert!(
                (max[ax] - 2.0).abs() < 1e-6,
                "max[{ax}]={} should be ~2",
                max[ax]
            );
        }
        let ext = mesh.bound_extents();
        assert!(
            ext.iter().all(|&e| (e - 2.0).abs() < 1e-6),
            "box extents should all be ~2, got {ext:?}"
        );
    }

    #[test]
    fn tessellate_sphere_is_denser_with_finer_tolerance() {
        let solid = k()
            .primitive(Primitive::Sphere { radius: 1.0 })
            .expect("sphere");
        let coarse = k().tessellate_tol(&solid, 0.3).expect("coarse mesh");
        let fine = k().tessellate_tol(&solid, 0.02).expect("fine mesh");
        assert!(coarse.triangle_count() > 0);
        assert!(
            fine.triangle_count() > coarse.triangle_count(),
            "finer tolerance ({} tris) should beat coarser ({} tris)",
            fine.triangle_count(),
            coarse.triangle_count()
        );
        // Sphere of radius 1 spans [-1,-1,-1]..[1,1,1]; the fine mesh's
        // bounds should approach that (inscribed facets undershoot a
        // little, so allow a small inward slack).
        let (min, max) = fine.bounds;
        for ax in 0..3 {
            assert!(min[ax] < -0.9 && min[ax] >= -1.001, "min[{ax}]={}", min[ax]);
            assert!(max[ax] > 0.9 && max[ax] <= 1.001, "max[{ax}]={}", max[ax]);
        }
    }

    #[test]
    fn boolean_id_round_trips_and_aliases_parse() {
        for op in [BoolOp::Union, BoolOp::Difference, BoolOp::Intersection] {
            assert_eq!(BoolOp::from_id(op.id()), Some(op));
            assert_eq!(BoolOp::from_id(&op.id().to_uppercase()), Some(op));
        }
        assert_eq!(BoolOp::from_id("or"), Some(BoolOp::Union));
        assert_eq!(BoolOp::from_id("cut"), Some(BoolOp::Difference));
        assert_eq!(BoolOp::from_id("and"), Some(BoolOp::Intersection));
        assert_eq!(BoolOp::from_id("nonsense"), None);
    }

    #[test]
    fn primitive_kind_ids() {
        assert_eq!(CUBE.kind_id(), "box");
        assert_eq!(
            Primitive::Cylinder {
                radius: 1.0,
                height: 1.0
            }
            .kind_id(),
            "cylinder"
        );
        assert_eq!(Primitive::Sphere { radius: 1.0 }.kind_id(), "sphere");
    }
}
