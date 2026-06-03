//! Feature variants — the six operations a Phase-2 part can perform.
//!
//! Each variant carries its parameter struct (`PadParams`, `PocketParams`,
//! etc.) plus a reference to whatever sketches or earlier features it
//! consumes. The [`Feature`] enum is the in-memory representation; the
//! RON-based [`crate::persist`] layer serializes the same shape.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use valenx_spreadsheet::{CellRef, Spreadsheet, SpreadsheetError};

use crate::threads::ThreadSpec;

/// Either a literal `f64` or a spreadsheet formula source string that
/// resolves to one (Phase 16).
///
/// Feature parameters that hold a [`Value`] can be either a hard-coded
/// number (the common case) or a parametric expression that references
/// a [`Spreadsheet`] cell. The [`crate::tree::FeatureTree`] resolves
/// every [`Value::Expression`] via
/// [`crate::tree::FeatureTree::replay_with_spreadsheet`] before
/// dispatching to the per-op evaluator.
///
/// `Default` is `Value::Literal(0.0)` so a freshly-constructed param
/// struct still typechecks; callers should always assign a meaningful
/// value before use.
///
/// ## Constructing from literals
///
/// `From<f64>` is implemented so existing call sites continue to work
/// with a small `.into()` sprinkle:
///
/// ```
/// use valenx_feature_tree::feature::Value;
///
/// let v: Value = 10.0.into();
/// assert_eq!(v.literal(), Some(10.0));
/// ```
///
/// ## Resolving against a spreadsheet
///
/// ```
/// use valenx_feature_tree::feature::Value;
/// use valenx_spreadsheet::{Cell, CellRef, Spreadsheet};
///
/// let mut ss = Spreadsheet::new();
/// ss.add_sheet("S");
/// ss.set_cell(&CellRef::parse("S.A1").unwrap(), Cell::Number(7.0)).unwrap();
///
/// let v = Value::Expression("S.A1 * 2".into());
/// assert_eq!(v.resolve(&ss).unwrap(), 14.0);
/// ```
///
/// `PartialEq` is derived structurally — `Literal(f64)` uses IEEE 754
/// semantics so a NaN literal never compares equal to itself, and an
/// undo snapshot containing a NaN-valued literal fails to dedupe.
/// The Part Design panel's drag-value widgets never let a NaN reach
/// here in practice.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Value {
    /// Concrete numeric value — most callers use this.
    Literal(f64),
    /// Spreadsheet formula source. Resolved against a
    /// [`Spreadsheet`] by [`Value::resolve`].
    Expression(String),
}

impl Default for Value {
    fn default() -> Self {
        Value::Literal(0.0)
    }
}

impl From<f64> for Value {
    fn from(n: f64) -> Self {
        Value::Literal(n)
    }
}

impl Value {
    /// Return the literal numeric value if this is a
    /// [`Value::Literal`]; `None` for an expression.
    pub fn literal(&self) -> Option<f64> {
        match self {
            Value::Literal(n) => Some(*n),
            Value::Expression(_) => None,
        }
    }

    /// Resolve this value against a [`Spreadsheet`].
    ///
    /// Literals pass through; expressions are parsed + evaluated.
    /// An [`Value::Expression`] referencing a cell that doesn't
    /// exist in `ss` returns `0.0` (consistent with
    /// [`Spreadsheet::cell`] returning [`valenx_spreadsheet::Cell::Empty`]
    /// for missing entries).
    ///
    /// # Errors
    ///
    /// Returns the underlying [`SpreadsheetError`] from parsing or
    /// evaluation. Callers usually map this into a
    /// [`crate::FeatureError::BadParameter`] with the offending
    /// parameter name.
    pub fn resolve(&self, ss: &Spreadsheet) -> Result<f64, SpreadsheetError> {
        match self {
            Value::Literal(n) => Ok(*n),
            Value::Expression(src) => {
                let expr = valenx_spreadsheet::parser::parse(src)?;
                valenx_spreadsheet::evaluator::evaluate(&expr, ss)
            }
        }
    }

    /// Construct a [`Value::Expression`] from a [`CellRef`]
    /// (convenience for "just point this param at one cell").
    pub fn cell(r: &CellRef) -> Self {
        Value::Expression(r.to_string())
    }
}

/// Identifier for a sketch stored alongside the feature tree.
///
/// `SketchRef(i)` indexes into [`crate::tree::FeatureTree::sketches`].
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct SketchRef(
    /// Zero-based index into the tree's sketches vector.
    pub usize,
);

/// Identifier for a feature within a tree.
///
/// `FeatureId(i)` indexes into [`crate::tree::FeatureTree::features`].
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub struct FeatureId(
    /// Zero-based index into the tree's features vector.
    pub usize,
);

/// Parameters for a Pad (extrude a sketch profile into a solid).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PadParams {
    /// Sketch whose closed profile is extruded.
    pub sketch: SketchRef,
    /// Extrusion length in model units (must be nonzero, finite).
    ///
    /// Phase 16: this is a [`Value`] so callers can plug in a
    /// spreadsheet expression (`Value::Expression("S.A1 * 2".into())`)
    /// instead of a literal. The standard [`crate::replay::replay`]
    /// resolves literals only; expressions require
    /// [`crate::FeatureTree::replay_with_spreadsheet`].
    pub depth: Value,
    /// `true` = extrude +Z, `false` = extrude -Z. Sketch working plane is XY.
    pub direction_positive: bool,
}

/// Parameters for a Pocket (boolean-subtract an extruded profile from the
/// preceding solid).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PocketParams {
    /// Sketch whose closed profile is extruded into the cutting tool.
    pub sketch: SketchRef,
    /// Pocket depth in model units (must be nonzero, finite).
    ///
    /// Phase 16: see [`PadParams::depth`] for the [`Value`] semantics.
    pub depth: Value,
    /// `true` = pocket goes +Z, `false` = -Z.
    pub direction_positive: bool,
}

/// Parameters for a Revolve (rotational sweep of a sketch profile about
/// an axis).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RevolveParams {
    /// Sketch whose profile is swept.
    pub sketch: SketchRef,
    /// Point that the rotation axis passes through (world coordinates).
    pub axis_origin: Vector3<f64>,
    /// Axis direction (need not be unit-length; will be normalized).
    pub axis_direction: Vector3<f64>,
    /// Revolution angle in radians (full sweep = 2*pi).
    ///
    /// Phase 16: see [`PadParams::depth`] for the [`Value`] semantics.
    pub angle: Value,
}

/// Parameters for a Mirror (planar reflection of an earlier feature's
/// solid across a plane).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MirrorParams {
    /// Earlier feature in the tree to mirror (must be solid-producing).
    pub target: FeatureId,
    /// A point on the mirror plane (world coordinates).
    pub plane_origin: Vector3<f64>,
    /// Plane normal (need not be unit-length; will be normalized).
    pub plane_normal: Vector3<f64>,
    /// `true` = combine the original and the mirrored copy into one
    /// solid (union); `false` = output only the mirrored copy.
    pub keep_original: bool,
}

/// Parameters for a LinearPattern (translate-and-union N instances of an
/// earlier feature along a direction).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LinearPatternParams {
    /// Earlier feature in the tree to pattern (must be solid-producing).
    pub target: FeatureId,
    /// Direction the pattern advances in (need not be unit-length).
    pub direction: Vector3<f64>,
    /// Number of total instances (including the original; minimum 1).
    pub count: u32,
    /// Distance between consecutive instances along `direction`.
    pub spacing: f64,
}

/// Parameters for a CircularPattern (rotate-and-union N instances of an
/// earlier feature about an axis).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CircularPatternParams {
    /// Earlier feature in the tree to pattern (must be solid-producing).
    pub target: FeatureId,
    /// Point that the rotation axis passes through.
    pub axis_origin: Vector3<f64>,
    /// Axis direction (need not be unit-length).
    pub axis_direction: Vector3<f64>,
    /// Number of total instances (including the original; minimum 1).
    pub count: u32,
    /// Total swept angle in radians; typical full circle = 2π.
    pub total_angle: f64,
}

/// Parameters for a Fillet (round every sharp convex edge of a target
/// feature's solid with a cylindrical strip).
///
/// **Phase 14 dispatch:** when the target is a `Solid::Brep`, the
/// evaluator tries `valenx-fillet-brep` first; on a soft error
/// (non-planar adjacency, non-convex edge, mesh-backed, or the
/// pending Phase 14.5 truck-substitution gap) it falls through to
/// the Phase 3 mesh-domain pipeline. The result type is
/// `Solid::Brep` only when the BRep path succeeds end-to-end;
/// otherwise it's `Solid::Mesh` as in Phase 3.
///
/// **Output caveat:** if the fillet falls through to mesh-domain,
/// downstream BRep ops (booleans, sweeps) against the result will
/// fail with `CadError::MeshBackedSolid`. Apply Fillets last in the
/// tree to avoid this hazard.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FilletParams {
    /// Earlier feature in the tree to fillet (must be solid-producing).
    pub target: FeatureId,
    /// Fillet radius in model units; must be > 0.
    pub radius: f64,
    /// Dihedral-angle threshold in degrees — edges sharper than this
    /// qualify for filleting. Typical value: 45°. 0° fillets every
    /// interior edge, including imperceptibly flat ones.
    pub threshold_deg: f64,
    /// Optional explicit list of 0-based edge indices to fillet
    /// (Phase 14). `None` (the default for backward compatibility)
    /// means "auto-select by `threshold_deg`". When `Some`, the
    /// indices are interpreted against `valenx_fillet_brep::bridge::
    /// unique_edges(brep)` for BRep targets and ignored for
    /// mesh-backed targets (mesh-domain has no concept of indexed
    /// BRep edges).
    #[serde(default)]
    pub edge_indices: Option<Vec<usize>>,
}

/// Parameters for a Chamfer (replace every sharp convex edge of a
/// target feature's solid with a flat bevel of constant width).
///
/// **Phase 14 dispatch:** same as [`FilletParams`] — BRep path first,
/// fall through to mesh-domain on soft errors.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ChamferParams {
    /// Earlier feature in the tree to chamfer (must be solid-producing).
    pub target: FeatureId,
    /// Chamfer width in model units — how far inward from the original
    /// edge each face is offset before the bevel connects them. Must
    /// be > 0.
    pub distance: f64,
    /// Dihedral-angle threshold in degrees (same semantics as Fillet).
    pub threshold_deg: f64,
    /// Optional explicit list of 0-based edge indices to chamfer
    /// (Phase 14). See [`FilletParams::edge_indices`] for semantics.
    #[serde(default)]
    pub edge_indices: Option<Vec<usize>>,
}

/// Parameters for an ImportedSolid (Phase 8 — STEP / IGES bring-in).
///
/// Holds the path the solid was loaded from so the project can re-read
/// it on next open. The solid itself is **not** serialised — too large
/// and lossy — so callers replay the tree, which re-reads from disk via
/// `valenx_step_iges::import`.
///
/// If the original file moved or was deleted, replay returns
/// [`crate::FeatureError::Io`].
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ImportedSolidParams {
    /// Absolute or relative path to the original STEP / IGES file.
    /// Persisted in the `.valenx` project envelope; re-read on load.
    pub source_path: String,
}

/// Parameters for an ImportedAdvanced (Phase 20 — STEP AP242 / IGES
/// trimmed surface bring-in). Same shape as [`ImportedSolidParams`],
/// plus optional AP242 metadata captured at import time.
///
/// Replay re-reads the source file via [`valenx_step_iges::import`]
/// (same dispatch as [`ImportedSolidParams`]) and additionally surfaces
/// the AP242 metadata (product hierarchy, materials, colors, GD&T)
/// captured at original import time.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ImportedAdvancedParams {
    /// Absolute or relative path to the original AP242 STEP or trimmed-
    /// surface IGES file. Persisted; re-read on load.
    pub source_path: String,
    /// Captured AP242 product-path components, if any.
    #[serde(default)]
    pub product_path: Vec<String>,
    /// Captured feature-history hints (BASED_FEATURE names).
    #[serde(default)]
    pub feature_hints: Vec<String>,
    /// Captured parametric (key, value) pairs.
    #[serde(default)]
    pub parametric_values: Vec<(String, String)>,
    /// Captured GD&T tolerance strings.
    #[serde(default)]
    pub tolerances: Vec<String>,
    /// Material names attached to the import.
    #[serde(default)]
    pub material_names: Vec<String>,
    /// True if the source file uses the AP242 schema (vs plain AP203/
    /// AP214 or vanilla IGES wireframe). Lets the feature-tree view
    /// surface a small "AP242" badge.
    #[serde(default)]
    pub is_ap242: bool,
}

// =============================================================================
// Phase 13 — new feature variants (10 additions)
// =============================================================================

/// Counterbore add-on for a [`HoleParams`] — a flat-bottomed cylindrical
/// recess at the top of the hole sized to recess a socket head.
///
/// `diameter` is the recess outer diameter (always > the drill
/// diameter); `depth` is how far it cuts into the part from the top
/// face.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CounterboreParams {
    /// Counterbore outer diameter in mm. Must be > drill diameter.
    pub diameter: f64,
    /// Counterbore depth (how deep the recess goes from the top face).
    pub depth: f64,
}

/// Countersink add-on for a [`HoleParams`] — a conical recess at the top
/// of the hole sized for a flat-head screw.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CountersinkParams {
    /// Outer (top) diameter of the conical recess.
    pub diameter: f64,
    /// Included angle of the cone in degrees (typical 82°, 90°, 100°,
    /// or 120° per fastener standard).
    pub angle_deg: f64,
}

/// Depth-control variants for a [`HoleParams`].
///
/// v1 supports the three common modes; `UpToFace` falls back to a long
/// blind cut when the referenced face cannot be resolved (typical for
/// mesh-backed solids).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum HoleDepthMode {
    /// Drill a fixed depth into the part.
    Blind {
        /// Depth from the entry face along the drill direction (mm).
        depth: f64,
    },
    /// Punch all the way through the base solid.
    Through,
    /// Drill down until reaching the referenced face (v1: emulated as
    /// a long blind cut tagged with the face hint).
    UpToFace {
        /// Free-form reference label for the limiting face (e.g.
        /// `"face #3"`). v1 does not resolve this against the BRep —
        /// stored for downstream consumers.
        face_ref: String,
    },
}

/// Parameters for a Hole — drills 1+ cylindrical pockets at the points
/// of a sketch, optionally with counterbore / countersink modifiers and
/// thread metadata.
///
/// v1: thread *geometry* is not modelled (no helical thread); the
/// [`ThreadSpec`] is metadata attached to the resulting solid for
/// downstream callouts. See module-level docs in [`crate::threads`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct HoleParams {
    /// Sketch whose Point entities define hole positions on the XY
    /// working plane. v1 reads all `Entity::Point` entries; lines /
    /// circles in the sketch are ignored.
    pub sketch: SketchRef,
    /// Depth-control variant.
    pub depth_mode: HoleDepthMode,
    /// Drill (through-bore) diameter in mm. Must be > 0.
    pub drill_diameter: f64,
    /// `true` = drill downward (-Z) into the base solid (typical);
    /// `false` = drill upward (+Z). Matches Pocket's flip semantics.
    pub direction_negative: bool,
    /// Optional counterbore (flat-bottom recess for a socket head).
    pub counterbore: Option<CounterboreParams>,
    /// Optional countersink (conical recess for a flat head).
    pub countersink: Option<CountersinkParams>,
    /// Optional thread specification — metadata only in v1.
    pub thread: Option<ThreadSpec>,
}

/// Parameters for a Loft — interpolates a surface between 2+ profile
/// sketches, optionally guided by additional rail/guide curves.
///
/// v1 limitation: uses tessellation-based intermediate cross-section
/// generation; the output is a mesh-backed [`valenx_cad::Solid`].
/// `guide_curves` are stored but not yet applied to the interpolation
/// (Phase 14+).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct LoftParams {
    /// Ordered list of profile sketches to interpolate between (2+).
    pub profile_sketches: Vec<SketchRef>,
    /// Optional guide curves that steer the loft. v1: stored for
    /// future use, not applied.
    pub guide_curves: Vec<SketchRef>,
    /// Close the loft into a periodic shape (first ↔ last profile).
    pub closed: bool,
    /// `true` = straight ruled connections between corresponding
    /// vertices on adjacent profiles; `false` = smoothed.
    pub ruled: bool,
}

/// Parameters for a Sweep — sweeps a profile sketch along a path
/// sketch, optionally with a twist.
///
/// v1: discrete sampling along the path; mesh-backed output.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SweepParams {
    /// 2D profile to be swept.
    pub profile_sketch: SketchRef,
    /// Path / centerline the profile follows.
    pub path_sketch: SketchRef,
    /// Total twist along the path in radians (0 = no twist).
    pub twist_angle: f64,
    /// `true` = the profile rotates with the path tangent (frenet);
    /// `false` = the profile keeps its world-space orientation.
    pub keep_profile_orientation: bool,
}

/// Parameters for a Pipe — a specialization of Sweep with a fixed
/// circular / polygonal cross-section and optional bend filleting.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PipeParams {
    /// Cross-section sketch (typically a circle / polygon).
    pub cross_section_sketch: SketchRef,
    /// Path the cross-section follows.
    pub centerline_sketch: SketchRef,
    /// Radius applied to corners of the centerline when filleting
    /// path bends. 0.0 = sharp corners (no fillet).
    pub bend_radius: f64,
}

/// Parameters for a Helix — sweeps a profile along a parametric helix
/// path (axis + pitch + turns + taper).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct HelixParams {
    /// Profile to sweep (the cross-section that gets coiled).
    pub profile_sketch: SketchRef,
    /// Distance between adjacent turns (mm).
    pub pitch: f64,
    /// Number of turns. May be fractional.
    pub turns: f64,
    /// World-space point the helix axis passes through.
    pub axis_origin: Vector3<f64>,
    /// Direction of the helix axis (normalised internally).
    pub axis_direction: Vector3<f64>,
    /// Taper half-angle in degrees: 0 = constant-radius cylinder,
    /// >0 grows the radius along +axis (conical helix).
    pub taper_angle: f64,
    /// `true` = left-handed (counter-clockwise viewed along +axis).
    pub left_handed: bool,
}

/// One operation in a [`MultiTransformParams`] — the recipe for one
/// instance of the patterned solid.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum TransformOp {
    /// Translate by `delta`.
    Translate {
        /// Translation vector applied to the cloned target.
        delta: Vector3<f64>,
    },
    /// Rotate by `angle_rad` around `axis` through the world origin.
    Rotate {
        /// Rotation axis (need not be unit-length).
        axis: Vector3<f64>,
        /// Rotation angle in radians (right-hand rule).
        angle_rad: f64,
    },
    /// Uniform scale by `factor` (applied per-vertex relative to
    /// world origin; v1 is mesh-domain only).
    Scale {
        /// Uniform scale factor (1.0 = no change).
        factor: f64,
    },
    /// Reflect across a plane through the world origin with the given
    /// normal (need not be unit-length).
    Mirror {
        /// Plane normal (need not be unit-length).
        plane_normal: Vector3<f64>,
    },
}

/// Parameters for a MultiTransform — apply N arbitrary transforms to a
/// target feature and union all instances.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MultiTransformParams {
    /// Earlier feature in the tree to clone.
    pub target: FeatureId,
    /// Sequence of transforms — one cloned + transformed copy per
    /// entry, plus the original. Empty list = the original alone.
    pub transforms: Vec<TransformOp>,
}

/// Parameters for a DraftAngle — tilt selected faces of a target solid
/// by `draft_angle_deg` about an axis in the neutral plane.
///
/// v1: tessellation-based approximation. Face indices are interpreted
/// against the *triangle index* of a mesh-backed proxy of the target
/// — for BRep targets we tessellate first and re-wrap.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DraftAngleParams {
    /// Earlier feature to draft.
    pub target: FeatureId,
    /// 0-based indices of triangles (post-tessellation) to draft. For
    /// real FreeCAD-like face selection, see Phase 14+.
    pub face_indices: Vec<usize>,
    /// Normal of the neutral plane (the axis about which to tilt).
    pub neutral_plane_normal: Vector3<f64>,
    /// Draft angle in degrees (>0 tilts faces outward).
    pub draft_angle_deg: f64,
}

/// Which side of the input solid a Shell op hollows toward.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ShellSide {
    /// Keep the original outer boundary, push the inner wall toward
    /// the centre by `thickness`. (Typical "hollow box" use.)
    Inward,
    /// Keep the original boundary as the inner wall, grow outward by
    /// `thickness`. (Less common — "shell around" semantics.)
    Outward,
}

/// Parameters for a Shell — hollow out the target solid leaving a
/// thin-walled shell with `face_indices_to_remove` removed entirely
/// (open faces).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ShellParams {
    /// Earlier feature to shell.
    pub target: FeatureId,
    /// 0-based triangle indices to *remove* (leave open). Empty list
    /// = closed shell on all faces.
    pub face_indices_to_remove: Vec<usize>,
    /// Wall thickness (mm). Must be > 0.
    pub thickness: f64,
    /// Which side the offset goes.
    pub inward_or_outward: ShellSide,
}

/// Parameters for a Thickness — add wall thickness to a single face of
/// a target solid (turn a "sheet" of geometry into a thin slab).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ThicknessParams {
    /// Earlier feature to thicken.
    pub target: FeatureId,
    /// 0-based triangle index of the face to thicken.
    pub face_index: usize,
    /// Slab thickness (mm). Positive = grow along face normal.
    pub thickness: f64,
}

/// Boolean operation kinds supported by a [`BooleanHistoryParams`].
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum BoolKind {
    /// Union (A ∪ B ∪ …).
    Union,
    /// Difference (A − B − …).
    Difference,
    /// Intersection (A ∩ B ∩ …).
    Intersection,
    /// Cross-section curve (v1: returns the intersection mesh as a
    /// thin slab — true curve extraction is Phase 14+).
    Section,
}

/// Parameters for a BooleanHistory — apply an N-way boolean operation
/// to the listed targets in order.
///
/// Differs from Pocket (subtract one extrusion from "the last solid")
/// in that BooleanHistory can union / intersect / difference *any* set
/// of evaluated features, not just the immediately preceding one.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct BooleanHistoryParams {
    /// Which operation to apply across the targets.
    pub operation: BoolKind,
    /// Ordered list of target features. For `Difference`, the first
    /// entry is the base and subsequent entries are subtracted; for
    /// `Union`/`Intersection`/`Section`, order doesn't matter
    /// semantically but is preserved for replay determinism.
    pub targets: Vec<FeatureId>,
}

/// One operation in a feature tree. Each variant wraps its
/// strongly-typed parameter struct so [`crate::replay::replay`] can
/// dispatch without juggling untyped maps.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Feature {
    /// Extrude a sketch profile into a fresh solid.
    Pad(PadParams),
    /// Subtract an extruded profile from the preceding solid.
    Pocket(PocketParams),
    /// Sweep a sketch profile about an axis.
    Revolve(RevolveParams),
    /// Reflect an earlier feature across a plane (optionally combine).
    Mirror(MirrorParams),
    /// Translate-and-union N instances of an earlier feature.
    LinearPattern(LinearPatternParams),
    /// Rotate-and-union N instances of an earlier feature.
    CircularPattern(CircularPatternParams),
    /// Round every sharp convex edge of an earlier feature with a
    /// cylindrical strip (v1: mesh-domain output).
    Fillet(FilletParams),
    /// Replace every sharp convex edge of an earlier feature with a
    /// flat bevel (v1: mesh-domain output).
    Chamfer(ChamferParams),
    /// A Solid imported from a STEP or IGES file (Phase 8). The solid
    /// is re-read from `source_path` during replay; the path is what's
    /// persisted in the `.valenx` project envelope.
    ImportedSolid(ImportedSolidParams),
    /// Drill 1+ holes at the points of a sketch (Phase 13A).
    Hole(HoleParams),
    /// Interpolate a surface between 2+ profile sketches (Phase 13B).
    Loft(LoftParams),
    /// Sweep a profile along a path sketch (Phase 13B).
    Sweep(SweepParams),
    /// Pipe — sweep specialised for tube / piping with bend radius
    /// (Phase 13B).
    Pipe(PipeParams),
    /// Coil a profile along a parametric helix (Phase 13C).
    Helix(HelixParams),
    /// Apply N arbitrary transforms to a target and union the results
    /// (Phase 13C).
    MultiTransform(MultiTransformParams),
    /// Draft (tilt) selected faces of a target solid (Phase 13D).
    DraftAngle(DraftAngleParams),
    /// Hollow out a solid into a thin-walled shell (Phase 13D).
    Shell(ShellParams),
    /// Add thickness to a single face of a target (Phase 13D).
    Thickness(ThicknessParams),
    /// General-purpose N-way boolean op (Phase 13E).
    BooleanHistory(BooleanHistoryParams),
    /// A Solid imported from a STEP AP242 or IGES trimmed-surface file
    /// (Phase 20). Like [`Feature::ImportedSolid`] but additionally
    /// carries the AP242 metadata (product hierarchy, materials,
    /// colors, parametric history hints).
    ImportedAdvanced(ImportedAdvancedParams),
}

impl Feature {
    /// Short label for UI display (tree-view kind column, status bar, etc).
    pub fn kind_label(&self) -> &'static str {
        match self {
            Feature::Pad(_) => "Pad",
            Feature::Pocket(_) => "Pocket",
            Feature::Revolve(_) => "Revolve",
            Feature::Mirror(_) => "Mirror",
            Feature::LinearPattern(_) => "Linear Pattern",
            Feature::CircularPattern(_) => "Circular Pattern",
            Feature::Fillet(_) => "Fillet",
            Feature::Chamfer(_) => "Chamfer",
            Feature::ImportedSolid(_) => "Imported",
            Feature::Hole(_) => "Hole",
            Feature::Loft(_) => "Loft",
            Feature::Sweep(_) => "Sweep",
            Feature::Pipe(_) => "Pipe",
            Feature::Helix(_) => "Helix",
            Feature::MultiTransform(_) => "Multi-Transform",
            Feature::DraftAngle(_) => "Draft Angle",
            Feature::Shell(_) => "Shell",
            Feature::Thickness(_) => "Thickness",
            Feature::BooleanHistory(_) => "Boolean History",
            Feature::ImportedAdvanced(_) => "Imported (AP242)",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_variant_constructs_and_has_kind_label() {
        let pad = Feature::Pad(PadParams {
            sketch: SketchRef(0),
            depth: 10.0.into(),
            direction_positive: true,
        });
        assert_eq!(pad.kind_label(), "Pad");

        let pocket = Feature::Pocket(PocketParams {
            sketch: SketchRef(1),
            depth: 5.0.into(),
            direction_positive: false,
        });
        assert_eq!(pocket.kind_label(), "Pocket");

        let revolve = Feature::Revolve(RevolveParams {
            sketch: SketchRef(2),
            axis_origin: Vector3::new(0.0, 0.0, 0.0),
            axis_direction: Vector3::new(0.0, 1.0, 0.0),
            angle: std::f64::consts::TAU.into(),
        });
        assert_eq!(revolve.kind_label(), "Revolve");

        let mirror = Feature::Mirror(MirrorParams {
            target: FeatureId(0),
            plane_origin: Vector3::zeros(),
            plane_normal: Vector3::new(1.0, 0.0, 0.0),
            keep_original: true,
        });
        assert_eq!(mirror.kind_label(), "Mirror");

        let lp = Feature::LinearPattern(LinearPatternParams {
            target: FeatureId(0),
            direction: Vector3::new(1.0, 0.0, 0.0),
            count: 4,
            spacing: 2.0,
        });
        assert_eq!(lp.kind_label(), "Linear Pattern");

        let cp = Feature::CircularPattern(CircularPatternParams {
            target: FeatureId(0),
            axis_origin: Vector3::zeros(),
            axis_direction: Vector3::new(0.0, 0.0, 1.0),
            count: 6,
            total_angle: std::f64::consts::TAU,
        });
        assert_eq!(cp.kind_label(), "Circular Pattern");

        let fillet = Feature::Fillet(FilletParams {
            target: FeatureId(0),
            radius: 0.5,
            threshold_deg: 45.0,
            edge_indices: None,
        });
        assert_eq!(fillet.kind_label(), "Fillet");

        let chamfer = Feature::Chamfer(ChamferParams {
            target: FeatureId(0),
            distance: 0.5,
            threshold_deg: 45.0,
            edge_indices: None,
        });
        assert_eq!(chamfer.kind_label(), "Chamfer");
    }

    #[test]
    fn phase_13_variants_have_kind_labels() {
        // Build one instance of each Phase 13 variant and confirm
        // kind_label() returns the expected display string.
        use crate::threads::iso_metric_table;
        let m6 = iso_metric_table()
            .into_iter()
            .find(|s| s.designation == "M6")
            .unwrap();
        let cases: Vec<(Feature, &'static str)> = vec![
            (
                Feature::Hole(HoleParams {
                    sketch: SketchRef(0),
                    depth_mode: HoleDepthMode::Through,
                    drill_diameter: 1.0,
                    direction_negative: true,
                    counterbore: None,
                    countersink: None,
                    thread: Some(m6),
                }),
                "Hole",
            ),
            (
                Feature::Loft(LoftParams {
                    profile_sketches: vec![SketchRef(0), SketchRef(1)],
                    guide_curves: vec![],
                    closed: false,
                    ruled: true,
                }),
                "Loft",
            ),
            (
                Feature::Sweep(SweepParams {
                    profile_sketch: SketchRef(0),
                    path_sketch: SketchRef(1),
                    twist_angle: 0.0,
                    keep_profile_orientation: true,
                }),
                "Sweep",
            ),
            (
                Feature::Pipe(PipeParams {
                    cross_section_sketch: SketchRef(0),
                    centerline_sketch: SketchRef(1),
                    bend_radius: 0.5,
                }),
                "Pipe",
            ),
            (
                Feature::Helix(HelixParams {
                    profile_sketch: SketchRef(0),
                    pitch: 1.0,
                    turns: 2.0,
                    axis_origin: Vector3::zeros(),
                    axis_direction: Vector3::new(0.0, 0.0, 1.0),
                    taper_angle: 0.0,
                    left_handed: false,
                }),
                "Helix",
            ),
            (
                Feature::MultiTransform(MultiTransformParams {
                    target: FeatureId(0),
                    transforms: vec![],
                }),
                "Multi-Transform",
            ),
            (
                Feature::DraftAngle(DraftAngleParams {
                    target: FeatureId(0),
                    face_indices: vec![],
                    neutral_plane_normal: Vector3::new(0.0, 0.0, 1.0),
                    draft_angle_deg: 5.0,
                }),
                "Draft Angle",
            ),
            (
                Feature::Shell(ShellParams {
                    target: FeatureId(0),
                    face_indices_to_remove: vec![],
                    thickness: 0.1,
                    inward_or_outward: ShellSide::Inward,
                }),
                "Shell",
            ),
            (
                Feature::Thickness(ThicknessParams {
                    target: FeatureId(0),
                    face_index: 0,
                    thickness: 0.5,
                }),
                "Thickness",
            ),
            (
                Feature::BooleanHistory(BooleanHistoryParams {
                    operation: BoolKind::Union,
                    targets: vec![FeatureId(0)],
                }),
                "Boolean History",
            ),
        ];
        for (f, expected) in cases {
            assert_eq!(f.kind_label(), expected, "kind_label mismatch for {f:?}");
        }
    }

    #[test]
    fn id_newtypes_are_distinct() {
        // SketchRef and FeatureId are both `pub usize` newtypes but are
        // type-distinct so callers can't accidentally pass one for the
        // other.
        let s = SketchRef(3);
        let f = FeatureId(3);
        assert_eq!(s.0, f.0);
        // Inequality is a type check, not a value check — these don't
        // compare directly. Just confirm copy/clone work.
        let s2 = s;
        let f2 = f;
        assert_eq!(s, s2);
        assert_eq!(f, f2);
    }

    #[test]
    fn value_literal_round_trip() {
        // Phase 16: Value::from(f64) produces a Literal; literal()
        // unwraps. Default is Literal(0.0).
        let v: Value = 2.5_f64.into();
        assert_eq!(v.literal(), Some(2.5));
        let d = Value::default();
        assert_eq!(d.literal(), Some(0.0));
    }

    #[test]
    fn value_expression_resolves_against_spreadsheet() {
        // Phase 16 Task 26: a Value::Expression that references a cell
        // resolves to the cell's value via the spreadsheet evaluator.
        use valenx_spreadsheet::{Cell, CellRef};
        let mut ss = Spreadsheet::new();
        ss.add_sheet("Sheet1");
        ss.set_cell(&CellRef::parse("Sheet1.A1").unwrap(), Cell::Number(5.0))
            .unwrap();
        let v = Value::Expression("Sheet1.A1 * 2".into());
        assert_eq!(v.resolve(&ss).unwrap(), 10.0);
    }

    #[test]
    fn value_expression_resolves_with_literal_passthrough() {
        // Literal Value passes through resolve() unchanged even when the
        // spreadsheet is empty.
        let v: Value = 7.5.into();
        assert_eq!(v.resolve(&Spreadsheet::new()).unwrap(), 7.5);
    }

    #[test]
    fn value_cell_constructor() {
        // Value::cell(CellRef) → Value::Expression("Sheet.A1").
        let r = CellRef::parse("S.B7").unwrap();
        let v = Value::cell(&r);
        match v {
            Value::Expression(s) => assert_eq!(s, "S.B7"),
            other => panic!("expected Expression, got {other:?}"),
        }
    }

    #[test]
    fn pad_with_expression_depth_uses_spreadsheet() {
        // Phase 16 Task 26 integration: Pad with depth =
        // Value::Expression("Sheet1.A1 * 2") + spreadsheet with
        // Sheet1.A1 = 5.0 -> resolved depth = 10.0. We don't run the
        // full pad evaluator here (it needs a real sketch + truck);
        // resolving the Value is the unit under test.
        use valenx_spreadsheet::{Cell, CellRef};
        let mut ss = Spreadsheet::new();
        ss.add_sheet("Sheet1");
        ss.set_cell(&CellRef::parse("Sheet1.A1").unwrap(), Cell::Number(5.0))
            .unwrap();
        let params = PadParams {
            sketch: SketchRef(0),
            depth: Value::Expression("Sheet1.A1 * 2".into()),
            direction_positive: true,
        };
        assert_eq!(params.depth.resolve(&ss).unwrap(), 10.0);
    }
}
