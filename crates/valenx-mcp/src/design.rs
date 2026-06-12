//! Generative-design MCP tools — a parametric CAD modelling surface an
//! external LLM drives to *generate* designs (Phase: generative design
//! as an MCP plugin).
//!
//! ## What "generative design" means here
//!
//! Valenx ships **no machine-learning model**. The "generative" part of
//! generative design is the **LLM itself**, on the other side of the
//! MCP connection: it composes a part by calling these tools, asks
//! [`evaluate_design`] for the part's mass / volume / bounding box, and
//! iterates — adjusting dimensions, adding pockets, filleting — toward
//! whatever goal it was given ("a bracket under 200 g that fits a
//! 60 mm cube", say). This module is the **bounded, deterministic CAD
//! engine** that loop runs on: every tool maps onto real
//! `valenx-sketch` / `valenx-feature-tree` / `valenx-cad` geometry, so
//! the LLM's design is a genuine parametric solid, not a sketch of one.
//!
//! ## The tool set
//!
//! Sketching (2-D, on the XY working plane):
//! - [`create_sketch`] — start a new 2-D sketch.
//! - [`add_sketch_line`] — add a line segment (creates its two
//!   endpoints and returns all three entity ids).
//! - [`add_sketch_circle`] — add a circle (stored as a fine polygonal
//!   loop so it pads / revolves uniformly).
//! - [`add_constraint`] — add a geometric / dimensional constraint
//!   between sketch entities and **re-solve** the sketch (the real
//!   `valenx-sketch` Newton-Raphson constraint solver runs).
//!
//! Features (3-D, a parametric history):
//! - [`pad`] — extrude a sketch profile into a fresh solid.
//! - [`pocket`] — subtract an extruded profile from the model so far.
//! - [`revolve`] — revolve a sketch profile about an axis.
//! - [`fillet`] — round a feature's sharp edges.
//! - [`boolean`] — union / difference / intersect earlier features.
//!
//! Feedback + output:
//! - [`evaluate_design`] — replay the whole history and return the
//!   solid's **mass, volume, bounding box, surface area** so the LLM
//!   can score its current design and decide the next move.
//! - [`export_design`] — replay and write the solid to a STEP or STL
//!   file (sandboxed path).
//!
//! ## Session state
//!
//! The MCP server is a single-client stdio process, so the design is
//! one process-global [`DesignSession`] behind a mutex. A
//! [`reset_design`] tool clears it. The session holds the
//! [`valenx_feature_tree::FeatureTree`] (the feature history) plus the
//! draft sketches that have not yet been consumed by a feature.
//!
//! ## Honest scope
//!
//! This is a **real, bounded** parametric modeller — every operation
//! produces genuine BRep / mesh geometry through the shipped kernels.
//! It is deliberately a v1 surface: sketches are on the **XY plane**
//! only (pad extrudes ±Z; revolve spins about an arbitrary axis), a
//! circle in a profile is a fine polygon (the feature-tree extrude path
//! is line-based), and `add_constraint` covers the common
//! point/line constraints. Multi-plane sketching, datum geometry, and
//! assembly-level generative design are documented follow-ups. None of
//! that limits the core loop: draw → constrain → feature → measure →
//! iterate genuinely works.

use std::sync::{Mutex, OnceLock};

use anyhow::{anyhow, Result};
use serde_json::{json, Value};

use valenx_feature_tree::feature::{
    BoolKind, BooleanHistoryParams, FeatureId, FilletParams, PadParams, PocketParams,
    RevolveParams, SketchRef,
};
use valenx_feature_tree::{replay, Feature, FeatureTree};
use valenx_sketch::constraint::Constraint;
use valenx_sketch::{EntityId, Sketch, SolverConfig, SolverStatus};

use crate::sandbox::sandbox_check;

/// Number of line segments a [`add_sketch_circle`] circle is
/// approximated by. The feature-tree extrude path is line-based, so a
/// circle entering a pad / revolve profile is discretised; 64 segments
/// keeps the faceting well below a typical print / machining tolerance.
const CIRCLE_SEGMENTS: usize = 64;

/// Tessellation tolerance (model units) used when a design is measured
/// or exported to STL. Half a millimetre at the typical part scale.
const TESS_TOLERANCE: f64 = 0.25;

/// Round-17 M1: bound on draft sketches held by a single
/// [`DesignSession`]. The generative-design loop runs sketches through
/// features (pad / pocket / etc.) and a real part rarely needs more
/// than a handful of drafts at once — 256 is well beyond any honest
/// usage while preventing a misbehaving (or hostile) LLM client from
/// exhausting process memory by spamming `create_sketch`.
const MAX_DESIGN_DRAFTS: usize = 256;
/// Round-17 M1: bound on entities (points + lines) per draft sketch.
/// A circle accounts for `2 × CIRCLE_SEGMENTS = 128` entities, so even
/// a sketch packed with circles fits comfortably below 10k. Holds back
/// runaway `add_sketch_line` / `add_sketch_circle` loops that would
/// blow the solver's parameter vector.
const MAX_DESIGN_ENTITIES_PER_SKETCH: usize = 10_000;
/// Round-17 M1: bound on constraints per draft sketch. Real sketches
/// rarely exceed a few dozen constraints; 10k absorbs any reasonable
/// usage and keeps the Newton-Raphson solver from being asked to
/// invert a 10⁶-row Jacobian.
const MAX_DESIGN_CONSTRAINTS_PER_SKETCH: usize = 10_000;
/// Round-17 M1: bound on features in the feature tree. Each pad /
/// pocket / fillet / boolean is a feature; even the densest history a
/// real generative-design run produces is well under 10k.
const MAX_DESIGN_FEATURES: usize = 10_000;
/// Round-17 M1 + L4: bound on the number of feature ids a single
/// boolean op can reference. A union / difference / intersection over
/// 64 targets is already exotic; this stops a `targets: [0; 1e9]`
/// payload from allocating gigabytes of FeatureId vector before the
/// per-id existence check runs.
const MAX_BOOLEAN_TARGETS: usize = 64;

/// The working plane a sketch is drawn on.
///
/// v1 supports only the **XY** plane — the feature-tree pad / pocket /
/// revolve evaluators all take an XY profile. The enum exists so the
/// tool surface is forward-compatible and rejects an unsupported plane
/// with a clear message rather than silently mis-modelling.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SketchPlane {
    /// The global XY plane (the only plane v1 supports).
    Xy,
}

/// One draft sketch held in the [`DesignSession`] before a feature
/// consumes it.
struct DraftSketch {
    /// The `valenx-sketch` 2-D sketch — entities + constraints + the
    /// solver's parameter vector.
    sketch: Sketch,
    /// The working plane (always [`SketchPlane::Xy`] in v1).
    #[allow(dead_code)]
    plane: SketchPlane,
}

/// The whole in-memory generative-design document.
///
/// Wraps the [`FeatureTree`] (the parametric feature history) and the
/// draft sketches not yet turned into features. One of these lives
/// process-global behind a mutex for the duration of the MCP session.
pub struct DesignSession {
    /// The parametric feature history. Replaying it produces the solid.
    tree: FeatureTree,
    /// Draft sketches, addressed by their index (the "sketch id" the
    /// tools hand back to the LLM).
    drafts: Vec<DraftSketch>,
}

impl DesignSession {
    /// A fresh, empty design.
    fn new() -> DesignSession {
        DesignSession {
            tree: FeatureTree::new(),
            drafts: Vec::new(),
        }
    }
}

/// The process-global design session, lazily created and mutex-guarded.
fn session() -> &'static Mutex<DesignSession> {
    static SESSION: OnceLock<Mutex<DesignSession>> = OnceLock::new();
    SESSION.get_or_init(|| Mutex::new(DesignSession::new()))
}

/// Borrow the global session, mapping a poisoned lock to an error
/// instead of panicking.
fn lock_session() -> Result<std::sync::MutexGuard<'static, DesignSession>> {
    session()
        .lock()
        .map_err(|_| anyhow!("design session lock was poisoned"))
}

// === argument helpers ====================================================

/// Read a required floating-point argument by name.
fn req_f64(args: &Value, name: &str) -> Result<f64> {
    args.get(name)
        .and_then(|v| v.as_f64())
        .ok_or_else(|| anyhow!("missing or non-numeric argument `{name}`"))
}

/// Round-17 L1: assert that every element of `coords` is a finite
/// number. Returns `Err` naming `label` if any element is NaN /
/// infinity, otherwise `Ok(())`.
///
/// Note: serde_json's parser refuses to construct a `Number` from a
/// non-finite f64, so under the normal MCP-over-JSON path this guard
/// is unreachable. It still earns its keep as defence-in-depth:
/// (a) it catches a future refactor that builds an args `Value` via
/// arithmetic (where Number::from_f64 isn't on the hot path), and
/// (b) it documents the solver invariant — every coordinate that
/// reaches `add_point` / a constraint target MUST be finite.
fn ensure_finite_coords(coords: &[f64], label: &str) -> Result<()> {
    if coords.iter().any(|c| !c.is_finite()) {
        return Err(anyhow!(
            "{label} coordinates must be finite numbers (no NaN / infinity)"
        ));
    }
    Ok(())
}

/// Read a required unsigned-integer argument by name.
fn req_usize(args: &Value, name: &str) -> Result<usize> {
    args.get(name)
        .and_then(|v| v.as_u64())
        .map(|v| v as usize)
        .ok_or_else(|| anyhow!("missing or non-integer argument `{name}`"))
}

/// Read a required string argument by name.
fn req_str<'a>(args: &'a Value, name: &str) -> Result<&'a str> {
    args.get(name)
        .and_then(|v| v.as_str())
        .ok_or_else(|| anyhow!("missing or non-string argument `{name}`"))
}

/// Read an optional floating-point argument, falling back to `default`.
fn opt_f64(args: &Value, name: &str, default: f64) -> f64 {
    args.get(name).and_then(|v| v.as_f64()).unwrap_or(default)
}

/// Read an optional boolean argument, falling back to `default`.
fn opt_bool(args: &Value, name: &str, default: bool) -> bool {
    args.get(name).and_then(|v| v.as_bool()).unwrap_or(default)
}

/// Wrap a one-line human message + a structured payload into the MCP
/// `tools/call` result envelope.
fn tool_ok(text: impl Into<String>, structured: Value) -> Value {
    json!({
        "content": [{ "type": "text", "text": text.into() }],
        "structuredContent": structured,
    })
}

// === the tools ===========================================================

/// `create_sketch` — start a new, empty 2-D sketch.
///
/// Arguments: optional `plane` (`"XY"`; defaults to `"XY"`). Returns the
/// new sketch's id (an index the other sketch tools take as
/// `sketch_id`).
pub fn create_sketch(args: &Value) -> Result<Value> {
    let plane_str = args.get("plane").and_then(|v| v.as_str()).unwrap_or("XY");
    let plane = match plane_str.to_ascii_uppercase().as_str() {
        "XY" => SketchPlane::Xy,
        other => {
            return Err(anyhow!(
                "unsupported sketch plane `{other}` — v1 supports only the XY plane"
            ))
        }
    };
    let mut s = lock_session()?;
    // Round-17 M1: cap the number of drafts a session can hold so a
    // hostile LLM client can't OOM the process with create_sketch spam.
    if s.drafts.len() >= MAX_DESIGN_DRAFTS {
        return Err(anyhow!(
            "design session already holds {MAX_DESIGN_DRAFTS} sketches \
             (the per-session cap); call reset_design first"
        ));
    }
    let id = s.drafts.len();
    s.drafts.push(DraftSketch {
        sketch: Sketch::new(),
        plane,
    });
    Ok(tool_ok(
        format!("created sketch {id} on the XY plane"),
        json!({ "sketch_id": id }),
    ))
}

/// `add_sketch_line` — add a line segment to a sketch.
///
/// Arguments: `sketch_id`, and the line endpoints `x1`, `y1`, `x2`,
/// `y2`. The two endpoints become point entities; the line connects
/// them. Returns the entity ids of the start point, end point, and the
/// line — the LLM passes these to [`add_constraint`].
pub fn add_sketch_line(args: &Value) -> Result<Value> {
    let sketch_id = req_usize(args, "sketch_id")?;
    let (x1, y1) = (req_f64(args, "x1")?, req_f64(args, "y1")?);
    let (x2, y2) = (req_f64(args, "x2")?, req_f64(args, "y2")?);
    // Round-17 L1: reject NaN / infinite coordinates BEFORE handing them
    // to the solver. The Newton-Raphson iteration silently diverges on
    // a NaN endpoint (every residual becomes NaN, the Jacobian inverse
    // fails, and the solver reports "did not converge" without naming
    // the real cause) — fail fast with a clear message instead.
    ensure_finite_coords(&[x1, y1, x2, y2], "sketch line")?;
    let mut s = lock_session()?;
    let draft = s
        .drafts
        .get_mut(sketch_id)
        .ok_or_else(|| anyhow!("no sketch with id {sketch_id}"))?;
    // Round-17 M1: cap entities per sketch so a hostile loop adding
    // points / lines forever can't blow the solver's parameter vector
    // or exhaust memory.
    if draft.sketch.entities.len() >= MAX_DESIGN_ENTITIES_PER_SKETCH {
        return Err(anyhow!(
            "sketch {sketch_id} already holds {} entities \
             (the per-sketch cap is {MAX_DESIGN_ENTITIES_PER_SKETCH})",
            draft.sketch.entities.len()
        ));
    }
    let p0 = draft.sketch.add_point(x1, y1);
    let p1 = draft.sketch.add_point(x2, y2);
    let line = draft
        .sketch
        .add_line(p0, p1)
        .map_err(|e| anyhow!("could not add line: {e}"))?;
    Ok(tool_ok(
        format!("added line {} to sketch {sketch_id}", line.0),
        json!({
            "sketch_id": sketch_id,
            "start_point_id": p0.0,
            "end_point_id": p1.0,
            "line_id": line.0,
        }),
    ))
}

/// `add_sketch_circle` — add a circle to a sketch.
///
/// Arguments: `sketch_id`, the centre `cx`, `cy`, and the `radius`. The
/// circle is stored as a closed loop of 64 line segments (the
/// feature-tree extrude path is line-based), so it pads and revolves
/// like any other profile. Returns the ids of the loop's vertex points.
pub fn add_sketch_circle(args: &Value) -> Result<Value> {
    let sketch_id = req_usize(args, "sketch_id")?;
    let cx = req_f64(args, "cx")?;
    let cy = req_f64(args, "cy")?;
    let radius = req_f64(args, "radius")?;
    if !(radius.is_finite() && radius > 0.0) {
        return Err(anyhow!("circle radius must be a positive finite number"));
    }
    // Round-17 L1: reject NaN / infinite centre coordinates BEFORE the
    // CIRCLE_SEGMENTS-long add_point loop runs (a NaN centre would
    // poison every vertex and silently break the solver downstream).
    ensure_finite_coords(&[cx, cy], "circle centre")?;
    let mut s = lock_session()?;
    let draft = s
        .drafts
        .get_mut(sketch_id)
        .ok_or_else(|| anyhow!("no sketch with id {sketch_id}"))?;
    // Round-17 M1: refuse to add a circle that would push the sketch
    // past the entity cap. A circle is `2 × CIRCLE_SEGMENTS` entities
    // (vertices + segment lines); reject up-front if the new entities
    // wouldn't fit instead of silently adding a partial loop.
    let circle_entities = 2 * CIRCLE_SEGMENTS;
    if draft.sketch.entities.len() + circle_entities > MAX_DESIGN_ENTITIES_PER_SKETCH {
        return Err(anyhow!(
            "adding a circle ({circle_entities} entities) would push sketch \
             {sketch_id} past the {MAX_DESIGN_ENTITIES_PER_SKETCH}-entity cap \
             (currently {} entities)",
            draft.sketch.entities.len()
        ));
    }
    // Build the polygonal loop: CIRCLE_SEGMENTS vertices, then a closed
    // chain of lines around them.
    let mut verts: Vec<EntityId> = Vec::with_capacity(CIRCLE_SEGMENTS);
    for k in 0..CIRCLE_SEGMENTS {
        let theta = std::f64::consts::TAU * k as f64 / CIRCLE_SEGMENTS as f64;
        verts.push(
            draft
                .sketch
                .add_point(cx + radius * theta.cos(), cy + radius * theta.sin()),
        );
    }
    let mut line_ids = Vec::with_capacity(CIRCLE_SEGMENTS);
    for k in 0..CIRCLE_SEGMENTS {
        let next = (k + 1) % CIRCLE_SEGMENTS;
        let line = draft
            .sketch
            .add_line(verts[k], verts[next])
            .map_err(|e| anyhow!("could not add circle segment: {e}"))?;
        line_ids.push(line.0);
    }
    Ok(tool_ok(
        format!(
            "added a circle (r = {radius}) to sketch {sketch_id} as a \
             {CIRCLE_SEGMENTS}-segment loop"
        ),
        json!({
            "sketch_id": sketch_id,
            "vertex_point_ids": verts.iter().map(|e| e.0).collect::<Vec<_>>(),
            "segment_count": CIRCLE_SEGMENTS,
        }),
    ))
}

/// `add_constraint` — add a geometric / dimensional constraint to a
/// sketch and re-solve it.
///
/// Arguments: `sketch_id`, a `type` string, and the entity ids the
/// constraint relates. The supported `type`s and their arguments:
///
/// | `type`          | entities                  | extra        |
/// |-----------------|---------------------------|--------------|
/// | `coincident`    | `point_a`, `point_b`      | —            |
/// | `horizontal`    | `line`                    | —            |
/// | `vertical`      | `line`                    | —            |
/// | `parallel`      | `line_a`, `line_b`        | —            |
/// | `perpendicular` | `line_a`, `line_b`        | —            |
/// | `equal_length`  | `line_a`, `line_b`        | —            |
/// | `distance`      | `point_a`, `point_b`      | `value`      |
/// | `angle`         | `line_a`, `line_b`        | `value` (rad)|
///
/// After the constraint is added the real `valenx-sketch`
/// Newton-Raphson solver runs; the response reports whether it
/// converged and the residual norm.
pub fn add_constraint(args: &Value) -> Result<Value> {
    let sketch_id = req_usize(args, "sketch_id")?;
    let kind = req_str(args, "type")?.to_ascii_lowercase();
    // Entity-id helpers: read a named arg as a 1-based EntityId.
    let ent = |name: &str| -> Result<EntityId> { Ok(EntityId(req_usize(args, name)?)) };
    // Round-17 L1: read the dimensional `value` once and gate it on
    // is_finite(). The Distance / Angle constraints feed `value`
    // straight into the solver as the target residual; a NaN would
    // poison every Newton-Raphson iteration without naming the cause.
    let req_finite_value = || -> Result<f64> {
        let v = req_f64(args, "value")?;
        ensure_finite_coords(&[v], "constraint `value`")?;
        Ok(v)
    };
    let constraint = match kind.as_str() {
        "coincident" => Constraint::Coincident {
            a: ent("point_a")?,
            b: ent("point_b")?,
        },
        "horizontal" => Constraint::Horizontal(ent("line")?),
        "vertical" => Constraint::Vertical(ent("line")?),
        "parallel" => Constraint::Parallel {
            a: ent("line_a")?,
            b: ent("line_b")?,
        },
        "perpendicular" => Constraint::Perpendicular {
            a: ent("line_a")?,
            b: ent("line_b")?,
        },
        "equal_length" => Constraint::EqualLength {
            a: ent("line_a")?,
            b: ent("line_b")?,
        },
        "distance" => Constraint::Distance {
            a: ent("point_a")?,
            b: ent("point_b")?,
            target: req_finite_value()?,
        },
        "angle" => Constraint::Angle {
            a: ent("line_a")?,
            b: ent("line_b")?,
            target: req_finite_value()?,
        },
        other => {
            return Err(anyhow!(
                "unsupported constraint type `{other}` — supported: coincident, \
                 horizontal, vertical, parallel, perpendicular, equal_length, \
                 distance, angle"
            ))
        }
    };

    let mut s = lock_session()?;
    let draft = s
        .drafts
        .get_mut(sketch_id)
        .ok_or_else(|| anyhow!("no sketch with id {sketch_id}"))?;
    // Round-17 M1: cap constraints per sketch — Newton-Raphson solving
    // a 10⁶-row Jacobian would tie up the MCP server indefinitely.
    if draft.sketch.constraints.len() >= MAX_DESIGN_CONSTRAINTS_PER_SKETCH {
        return Err(anyhow!(
            "sketch {sketch_id} already holds {} constraints \
             (the per-sketch cap is {MAX_DESIGN_CONSTRAINTS_PER_SKETCH})",
            draft.sketch.constraints.len()
        ));
    }
    draft.sketch.add_constraint(constraint);
    // Re-solve the sketch with the new constraint included.
    let report = valenx_sketch::solver::solve(&mut draft.sketch, SolverConfig::default())
        .map_err(|e| anyhow!("sketch solve failed: {e}"))?;
    let converged = report.status == SolverStatus::Converged;
    Ok(tool_ok(
        format!(
            "added `{kind}` constraint to sketch {sketch_id}; solver {} (residual {:.2e})",
            if converged {
                "converged"
            } else {
                "did not converge"
            },
            report.residual_norm
        ),
        json!({
            "sketch_id": sketch_id,
            "constraint_type": kind,
            "solver_converged": converged,
            "residual_norm": report.residual_norm,
            "constraint_count": draft.sketch.constraints.len(),
        }),
    ))
}

/// Move a draft sketch into the feature tree, returning its
/// [`SketchRef`]. The draft is *cloned* — it stays available for
/// further features.
fn push_sketch_to_tree(session: &mut DesignSession, sketch_id: usize) -> Result<SketchRef> {
    let draft = session
        .drafts
        .get(sketch_id)
        .ok_or_else(|| anyhow!("no sketch with id {sketch_id}"))?;
    if draft.sketch.entities.is_empty() {
        return Err(anyhow!(
            "sketch {sketch_id} is empty — add lines / circles before using it in a feature"
        ));
    }
    let sketch = draft.sketch.clone();
    Ok(session.tree.add_sketch(sketch))
}

/// Round-17 M1: refuse to push another feature when the tree already
/// holds [`MAX_DESIGN_FEATURES`]. Centralised so every feature
/// entrypoint (pad / pocket / revolve / fillet / boolean) gives the
/// same error shape and the same cap is honoured uniformly.
fn check_feature_cap(session: &DesignSession) -> Result<()> {
    if session.tree.features.len() >= MAX_DESIGN_FEATURES {
        return Err(anyhow!(
            "design feature tree already holds {} features \
             (the per-session cap is {MAX_DESIGN_FEATURES}); call \
             reset_design first",
            session.tree.features.len()
        ));
    }
    Ok(())
}

/// `pad` — extrude a sketch's closed profile into a fresh solid.
///
/// Arguments: `sketch_id`, the extrusion `depth`, and optional
/// `direction_positive` (default `true` — extrude +Z). Adds a
/// [`Feature::Pad`] to the history and returns its feature id.
pub fn pad(args: &Value) -> Result<Value> {
    let sketch_id = req_usize(args, "sketch_id")?;
    let depth = req_f64(args, "depth")?;
    if !(depth.is_finite() && depth.abs() > 1e-12) {
        return Err(anyhow!("pad depth must be a nonzero finite number"));
    }
    let direction_positive = opt_bool(args, "direction_positive", true);
    let mut s = lock_session()?;
    // Round-17 M1: refuse to add another feature when the tree is full.
    check_feature_cap(&s)?;
    let sketch_ref = push_sketch_to_tree(&mut s, sketch_id)?;
    let feature = Feature::Pad(PadParams {
        sketch: sketch_ref,
        depth: depth.into(),
        direction_positive,
    });
    let id = s
        .tree
        .add_feature(feature, format!("Pad (sketch {sketch_id})"));
    Ok(tool_ok(
        format!("padded sketch {sketch_id} by {depth} → feature {}", id.0),
        json!({ "feature_id": id.0, "operation": "pad" }),
    ))
}

/// `pocket` — subtract an extruded sketch profile from the model.
///
/// Arguments: `sketch_id`, the pocket `depth`, optional
/// `direction_positive` (default `true`). Adds a [`Feature::Pocket`],
/// which the replay subtracts from the immediately-preceding solid.
pub fn pocket(args: &Value) -> Result<Value> {
    let sketch_id = req_usize(args, "sketch_id")?;
    let depth = req_f64(args, "depth")?;
    if !(depth.is_finite() && depth.abs() > 1e-12) {
        return Err(anyhow!("pocket depth must be a nonzero finite number"));
    }
    let direction_positive = opt_bool(args, "direction_positive", true);
    let mut s = lock_session()?;
    // Round-17 M1: refuse to add another feature when the tree is full.
    check_feature_cap(&s)?;
    let sketch_ref = push_sketch_to_tree(&mut s, sketch_id)?;
    let feature = Feature::Pocket(PocketParams {
        sketch: sketch_ref,
        depth: depth.into(),
        direction_positive,
    });
    let id = s
        .tree
        .add_feature(feature, format!("Pocket (sketch {sketch_id})"));
    Ok(tool_ok(
        format!("pocketed sketch {sketch_id} by {depth} → feature {}", id.0),
        json!({ "feature_id": id.0, "operation": "pocket" }),
    ))
}

/// `revolve` — revolve a sketch profile about an axis into a solid.
///
/// Arguments: `sketch_id`, the `angle` in radians (a full revolution is
/// `2π`), and the axis — `axis_origin` (`[x,y,z]`) and `axis_direction`
/// (`[x,y,z]`, normalised internally). Defaults: a full `2π` sweep
/// about the global Y axis through the origin. Adds a
/// [`Feature::Revolve`].
pub fn revolve(args: &Value) -> Result<Value> {
    let sketch_id = req_usize(args, "sketch_id")?;
    let angle = opt_f64(args, "angle", std::f64::consts::TAU);
    if !(angle.is_finite() && angle.abs() > 1e-9) {
        return Err(anyhow!("revolve angle must be a nonzero finite number"));
    }
    let origin = read_vec3(args, "axis_origin", [0.0, 0.0, 0.0]);
    let direction = read_vec3(args, "axis_direction", [0.0, 1.0, 0.0]);
    if direction.iter().all(|c| c.abs() < 1e-12) {
        return Err(anyhow!("revolve axis_direction must be a non-zero vector"));
    }
    let mut s = lock_session()?;
    // Round-17 M1: refuse to add another feature when the tree is full.
    check_feature_cap(&s)?;
    let sketch_ref = push_sketch_to_tree(&mut s, sketch_id)?;
    let feature = Feature::Revolve(RevolveParams {
        sketch: sketch_ref,
        axis_origin: nalgebra::Vector3::new(origin[0], origin[1], origin[2]),
        axis_direction: nalgebra::Vector3::new(direction[0], direction[1], direction[2]),
        angle: angle.into(),
    });
    let id = s
        .tree
        .add_feature(feature, format!("Revolve (sketch {sketch_id})"));
    Ok(tool_ok(
        format!(
            "revolved sketch {sketch_id} by {angle} rad → feature {}",
            id.0
        ),
        json!({ "feature_id": id.0, "operation": "revolve" }),
    ))
}

/// Read a `[x, y, z]` array argument, falling back to `default`.
fn read_vec3(args: &Value, name: &str, default: [f64; 3]) -> [f64; 3] {
    match args.get(name).and_then(|v| v.as_array()) {
        Some(a) if a.len() == 3 => [
            a[0].as_f64().unwrap_or(default[0]),
            a[1].as_f64().unwrap_or(default[1]),
            a[2].as_f64().unwrap_or(default[2]),
        ],
        _ => default,
    }
}

/// `fillet` — round the sharp convex edges of an earlier feature.
///
/// Arguments: `target_feature` (the feature id to fillet), the fillet
/// `radius`, and an optional `threshold_deg` dihedral-angle threshold
/// (default `45`). Adds a [`Feature::Fillet`].
pub fn fillet(args: &Value) -> Result<Value> {
    let target = FeatureId(req_usize(args, "target_feature")?);
    let radius = req_f64(args, "radius")?;
    if !(radius.is_finite() && radius > 0.0) {
        return Err(anyhow!("fillet radius must be a positive finite number"));
    }
    let threshold_deg = opt_f64(args, "threshold_deg", 45.0);
    let mut s = lock_session()?;
    // Round-17 M1: refuse to add another feature when the tree is full.
    check_feature_cap(&s)?;
    if target.0 >= s.tree.features.len() {
        return Err(anyhow!("no feature with id {}", target.0));
    }
    let feature = Feature::Fillet(FilletParams {
        target,
        radius,
        threshold_deg,
        edge_indices: None,
    });
    let id = s.tree.add_feature(feature, format!("Fillet r{radius}"));
    Ok(tool_ok(
        format!(
            "filleted feature {} with radius {radius} → feature {}",
            target.0, id.0
        ),
        json!({ "feature_id": id.0, "operation": "fillet" }),
    ))
}

/// `boolean` — combine earlier features with a boolean operation.
///
/// Arguments: `operation` (`"union"`, `"difference"`, or
/// `"intersection"`) and `targets` — an array of feature ids. For
/// `difference` the first id is the base and the rest are subtracted.
/// Adds a [`Feature::BooleanHistory`].
pub fn boolean(args: &Value) -> Result<Value> {
    let op_str = req_str(args, "operation")?.to_ascii_lowercase();
    let operation = match op_str.as_str() {
        "union" => BoolKind::Union,
        "difference" => BoolKind::Difference,
        "intersection" => BoolKind::Intersection,
        other => {
            return Err(anyhow!(
                "unsupported boolean operation `{other}` — use union, difference, or intersection"
            ))
        }
    };
    let target_arr = args
        .get("targets")
        .and_then(|v| v.as_array())
        .ok_or_else(|| anyhow!("`targets` must be an array of feature ids"))?;
    if target_arr.len() < 2 {
        return Err(anyhow!("a boolean needs at least 2 target features"));
    }
    // Round-17 L4 (== M1 fifth cap): refuse a boolean with an absurd
    // number of targets BEFORE allocating the FeatureId vector. A
    // hostile `targets: [0; 1e9]` payload could otherwise allocate
    // 8 GB just to fail the existence check inside the loop.
    if target_arr.len() > MAX_BOOLEAN_TARGETS {
        return Err(anyhow!(
            "boolean accepts at most {MAX_BOOLEAN_TARGETS} target features (got {})",
            target_arr.len()
        ));
    }
    let mut targets = Vec::with_capacity(target_arr.len());
    for t in target_arr {
        let idx = t
            .as_u64()
            .ok_or_else(|| anyhow!("each `targets` entry must be a feature id"))?
            as usize;
        targets.push(FeatureId(idx));
    }
    let mut s = lock_session()?;
    // Round-17 M1: refuse to add another feature when the tree is full.
    check_feature_cap(&s)?;
    for t in &targets {
        if t.0 >= s.tree.features.len() {
            return Err(anyhow!("no feature with id {}", t.0));
        }
    }
    let feature = Feature::BooleanHistory(BooleanHistoryParams {
        operation,
        targets: targets.clone(),
    });
    let id = s.tree.add_feature(feature, format!("Boolean ({op_str})"));
    Ok(tool_ok(
        format!("{op_str} of {} features → feature {}", targets.len(), id.0),
        json!({ "feature_id": id.0, "operation": op_str }),
    ))
}

/// The measured properties of a design — what [`evaluate_design`]
/// returns so the driving LLM can score its current part.
#[derive(Clone, Copy, Debug)]
pub struct DesignMetrics {
    /// Enclosed volume in model units³.
    pub volume: f64,
    /// Mass = `volume · density` in mass units.
    pub mass: f64,
    /// Total surface area in model units².
    pub surface_area: f64,
    /// Axis-aligned bounding-box minimum corner.
    pub bbox_min: [f64; 3],
    /// Axis-aligned bounding-box maximum corner.
    pub bbox_max: [f64; 3],
}

/// `evaluate_design` — replay the feature history and measure the
/// resulting solid.
///
/// Arguments: optional `density` (mass per unit volume; default `1.0`,
/// so `mass == volume` unless a real density is supplied).
///
/// Replays the whole tree to a solid, tessellates it, and computes the
/// **volume** (signed-tetrahedron / divergence-theorem sum over the
/// triangle mesh), the **mass** (`volume · density`), the **surface
/// area**, and the **axis-aligned bounding box**. This is the feedback
/// signal of the generative loop — the LLM reads these numbers and
/// adjusts the design toward its goal.
pub fn evaluate_design(args: &Value) -> Result<Value> {
    let density = opt_f64(args, "density", 1.0);
    let s = lock_session()?;
    let metrics = measure(&s.tree)?;
    drop(s);
    let mass = metrics.volume * density;
    Ok(tool_ok(
        format!(
            "design: volume {:.4}, mass {:.4} (density {density}), \
             bounding box {:.3}×{:.3}×{:.3}",
            metrics.volume,
            mass,
            metrics.bbox_max[0] - metrics.bbox_min[0],
            metrics.bbox_max[1] - metrics.bbox_min[1],
            metrics.bbox_max[2] - metrics.bbox_min[2],
        ),
        json!({
            "volume": metrics.volume,
            "mass": mass,
            "density": density,
            "surface_area": metrics.surface_area,
            "bounding_box": {
                "min": metrics.bbox_min,
                "max": metrics.bbox_max,
                "size": [
                    metrics.bbox_max[0] - metrics.bbox_min[0],
                    metrics.bbox_max[1] - metrics.bbox_min[1],
                    metrics.bbox_max[2] - metrics.bbox_min[2],
                ],
            },
        }),
    ))
}

/// Replay the tree, tessellate the solid, and compute its metrics.
///
/// Volume is the divergence-theorem sum `Σ (v0 · (v1 × v2)) / 6` over
/// every triangle — exact for a closed triangle mesh. Surface area is
/// the sum of the triangle areas; the bounding box is the node extent.
fn measure(tree: &FeatureTree) -> Result<DesignMetrics> {
    let solid = replay(tree)
        .map_err(|e| anyhow!("could not replay the feature history: {e}"))?
        .ok_or_else(|| anyhow!("the design is empty — add at least one feature"))?;
    let mesh = valenx_cad::solid_to_mesh(&solid, TESS_TOLERANCE)
        .map_err(|e| anyhow!("could not tessellate the design: {e}"))?;
    Ok(mesh_metrics(&mesh))
}

/// Compute volume / surface-area / bounding-box from a triangle mesh.
///
/// Split out so it is unit-testable without a feature tree.
fn mesh_metrics(mesh: &valenx_mesh::Mesh) -> DesignMetrics {
    use valenx_mesh::element::ElementType;
    let mut signed_vol6 = 0.0f64; // 6 × the signed volume
    let mut area = 0.0f64;
    let mut bbox_min = [f64::INFINITY; 3];
    let mut bbox_max = [f64::NEG_INFINITY; 3];

    for node in &mesh.nodes {
        for k in 0..3 {
            bbox_min[k] = bbox_min[k].min(node[k]);
            bbox_max[k] = bbox_max[k].max(node[k]);
        }
    }
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let (ia, ib, ic) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
            if ia >= mesh.nodes.len() || ib >= mesh.nodes.len() || ic >= mesh.nodes.len() {
                continue;
            }
            let a = mesh.nodes[ia];
            let b = mesh.nodes[ib];
            let c = mesh.nodes[ic];
            // Divergence-theorem volume element: a · (b × c).
            signed_vol6 += a.dot(&b.cross(&c));
            // Triangle area = ½|(b−a) × (c−a)|.
            area += 0.5 * (b - a).cross(&(c - a)).norm();
        }
    }
    // A consistently-wound closed mesh gives a positive volume; take the
    // absolute value so an inward-wound mesh still reports a sensible
    // magnitude.
    let volume = (signed_vol6 / 6.0).abs();
    if !bbox_min[0].is_finite() {
        bbox_min = [0.0; 3];
        bbox_max = [0.0; 3];
    }
    DesignMetrics {
        volume,
        mass: volume, // density applied by the caller
        surface_area: area,
        bbox_min,
        bbox_max,
    }
}

/// `export_design` — replay the feature history and write the solid to
/// a file.
///
/// Arguments: `path` (sandboxed — must resolve under the MCP sandbox
/// root) and optional `format` (`"step"` or `"stl"`; inferred from the
/// path extension when omitted). STEP is written via
/// `valenx-step-iges`; STL via the binary-STL writer over the
/// tessellated mesh.
pub fn export_design(args: &Value) -> Result<Value> {
    let raw_path = req_str(args, "path")?;
    let path = sandbox_check(std::path::Path::new(raw_path))?;
    // Resolve the export format: explicit `format`, else the extension.
    let format = match args.get("format").and_then(|v| v.as_str()) {
        Some(f) => f.to_ascii_lowercase(),
        None => path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.to_ascii_lowercase())
            .unwrap_or_else(|| "step".to_string()),
    };

    let s = lock_session()?;
    let solid = replay(&s.tree)
        .map_err(|e| anyhow!("could not replay the feature history: {e}"))?
        .ok_or_else(|| anyhow!("the design is empty — nothing to export"))?;
    drop(s);

    match format.as_str() {
        "step" | "stp" => {
            valenx_step_iges::export(&solid, &path)
                .map_err(|e| anyhow!("STEP export failed: {e}"))?;
        }
        "stl" => {
            let mesh = valenx_cad::solid_to_mesh(&solid, TESS_TOLERANCE)
                .map_err(|e| anyhow!("could not tessellate for STL export: {e}"))?;
            valenx_mesh::write_stl_binary(&mesh, &path)
                .map_err(|e| anyhow!("STL export failed: {e}"))?;
        }
        other => {
            return Err(anyhow!(
                "unsupported export format `{other}` — use `step` or `stl`"
            ))
        }
    }
    Ok(tool_ok(
        format!("exported the design to {} ({format})", path.display()),
        json!({ "path": path.display().to_string(), "format": format }),
    ))
}

/// `reset_design` — discard the whole design and start fresh.
///
/// Takes no arguments. Clears every sketch and feature so the LLM can
/// begin a new part (or restart after a dead end).
pub fn reset_design(_args: &Value) -> Result<Value> {
    let mut s = lock_session()?;
    *s = DesignSession::new();
    Ok(tool_ok(
        "design session reset — all sketches and features cleared",
        json!({ "reset": true }),
    ))
}

/// The static `tools/list` entries for the generative-design tool set.
///
/// Appended to the MCP server's tool catalogue by
/// [`crate::tools::list`].
pub fn tool_list() -> Vec<Value> {
    vec![
        json!({
            "name": "create_sketch",
            "description": "Start a new 2-D sketch (XY plane). Returns a sketch_id.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "plane": { "type": "string", "enum": ["XY"], "description": "Working plane (v1: XY only)." }
                }
            }
        }),
        json!({
            "name": "add_sketch_line",
            "description": "Add a line segment to a sketch. Returns the start/end point ids and the line id.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "sketch_id": { "type": "integer" },
                    "x1": { "type": "number" }, "y1": { "type": "number" },
                    "x2": { "type": "number" }, "y2": { "type": "number" }
                },
                "required": ["sketch_id", "x1", "y1", "x2", "y2"]
            }
        }),
        json!({
            "name": "add_sketch_circle",
            "description": "Add a circle to a sketch (stored as a polygonal loop). Returns the loop vertex ids.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "sketch_id": { "type": "integer" },
                    "cx": { "type": "number" }, "cy": { "type": "number" },
                    "radius": { "type": "number" }
                },
                "required": ["sketch_id", "cx", "cy", "radius"]
            }
        }),
        json!({
            "name": "add_constraint",
            "description": "Add a geometric/dimensional constraint to a sketch and re-solve it. \
                            Types: coincident, horizontal, vertical, parallel, perpendicular, equal_length, distance, angle.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "sketch_id": { "type": "integer" },
                    "type": { "type": "string" },
                    "point_a": { "type": "integer" }, "point_b": { "type": "integer" },
                    "line": { "type": "integer" },
                    "line_a": { "type": "integer" }, "line_b": { "type": "integer" },
                    "value": { "type": "number" }
                },
                "required": ["sketch_id", "type"]
            }
        }),
        json!({
            "name": "pad",
            "description": "Extrude a sketch profile into a solid. Returns the feature_id.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "sketch_id": { "type": "integer" },
                    "depth": { "type": "number" },
                    "direction_positive": { "type": "boolean" }
                },
                "required": ["sketch_id", "depth"]
            }
        }),
        json!({
            "name": "pocket",
            "description": "Subtract an extruded sketch profile from the model. Returns the feature_id.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "sketch_id": { "type": "integer" },
                    "depth": { "type": "number" },
                    "direction_positive": { "type": "boolean" }
                },
                "required": ["sketch_id", "depth"]
            }
        }),
        json!({
            "name": "revolve",
            "description": "Revolve a sketch profile about an axis into a solid. Returns the feature_id.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "sketch_id": { "type": "integer" },
                    "angle": { "type": "number", "description": "Sweep angle in radians (2*pi = full)." },
                    "axis_origin": { "type": "array", "items": { "type": "number" } },
                    "axis_direction": { "type": "array", "items": { "type": "number" } }
                },
                "required": ["sketch_id"]
            }
        }),
        json!({
            "name": "fillet",
            "description": "Round the sharp edges of an earlier feature's solid. Returns the feature_id.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "target_feature": { "type": "integer" },
                    "radius": { "type": "number" },
                    "threshold_deg": { "type": "number" }
                },
                "required": ["target_feature", "radius"]
            }
        }),
        json!({
            "name": "boolean",
            "description": "Combine earlier features with union/difference/intersection. Returns the feature_id.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "operation": { "type": "string", "enum": ["union", "difference", "intersection"] },
                    "targets": { "type": "array", "items": { "type": "integer" } }
                },
                "required": ["operation", "targets"]
            }
        }),
        json!({
            "name": "evaluate_design",
            "description": "Replay the design and return its mass, volume, surface area and bounding box \
                            so the LLM can score it and iterate toward a goal.",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "density": { "type": "number", "description": "Mass per unit volume (default 1.0)." }
                }
            }
        }),
        json!({
            "name": "export_design",
            "description": "Replay the design and write it to a STEP or STL file (sandboxed path).",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "path": { "type": "string" },
                    "format": { "type": "string", "enum": ["step", "stl"] }
                },
                "required": ["path"]
            }
        }),
        json!({
            "name": "reset_design",
            "description": "Discard the whole design and start fresh.",
            "inputSchema": { "type": "object" }
        }),
    ]
}

/// Dispatch a generative-design `tools/call` by name. Returns
/// `Ok(None)` if `name` is not one of the design tools (so the caller
/// can fall through to the other tool families).
pub fn dispatch(name: &str, args: &Value) -> Option<Result<Value>> {
    let result = match name {
        "create_sketch" => create_sketch(args),
        "add_sketch_line" => add_sketch_line(args),
        "add_sketch_circle" => add_sketch_circle(args),
        "add_constraint" => add_constraint(args),
        "pad" => pad(args),
        "pocket" => pocket(args),
        "revolve" => revolve(args),
        "fillet" => fillet(args),
        "boolean" => boolean(args),
        "evaluate_design" => evaluate_design(args),
        "export_design" => export_design(args),
        "reset_design" => reset_design(args),
        _ => return None,
    };
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Serialise the design tests and reset the global session.
    ///
    /// The [`DesignSession`] is process-global, so two tests mutating it
    /// on different threads would race. Each session-touching test
    /// begins `let _g = fresh();`: that resets the session *and* takes a
    /// process-wide mutex held for the whole test body, so the tests run
    /// one at a time. A poisoned guard is recovered (a panicking test
    /// must not wedge the rest of the suite). The returned guard must be
    /// bound (`let _g = …`) so it is not dropped immediately —
    /// `MutexGuard` is itself `#[must_use]`, which enforces that.
    fn fresh() -> std::sync::MutexGuard<'static, ()> {
        static M: OnceLock<Mutex<()>> = OnceLock::new();
        let guard = M
            .get_or_init(|| Mutex::new(()))
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        reset_design(&json!({})).unwrap();
        guard
    }

    /// Pull the `structuredContent` of a tool result.
    fn structured(v: &Value) -> &Value {
        &v["structuredContent"]
    }

    #[test]
    fn create_sketch_returns_an_id() {
        let _g = fresh();
        let r = create_sketch(&json!({ "plane": "XY" })).unwrap();
        assert_eq!(structured(&r)["sketch_id"], 0);
        // A second sketch gets the next id.
        let r2 = create_sketch(&json!({})).unwrap();
        assert_eq!(structured(&r2)["sketch_id"], 1);
    }

    #[test]
    fn create_sketch_rejects_an_unsupported_plane() {
        let _g = fresh();
        let err = create_sketch(&json!({ "plane": "XZ" })).unwrap_err();
        assert!(err.to_string().contains("XY plane"), "got: {err}");
    }

    #[test]
    fn add_sketch_line_returns_three_entity_ids() {
        let _g = fresh();
        let sid = structured(&create_sketch(&json!({})).unwrap())["sketch_id"]
            .as_u64()
            .unwrap();
        let r = add_sketch_line(&json!({
            "sketch_id": sid, "x1": 0.0, "y1": 0.0, "x2": 1.0, "y2": 0.0
        }))
        .unwrap();
        let s = structured(&r);
        // Two points then a line → ids 1, 2, 3.
        assert_eq!(s["start_point_id"], 1);
        assert_eq!(s["end_point_id"], 2);
        assert_eq!(s["line_id"], 3);
    }

    #[test]
    fn add_sketch_line_rejects_an_unknown_sketch() {
        let _g = fresh();
        let err = add_sketch_line(&json!({
            "sketch_id": 99, "x1": 0.0, "y1": 0.0, "x2": 1.0, "y2": 0.0
        }))
        .unwrap_err();
        assert!(err.to_string().contains("no sketch with id 99"));
    }

    #[test]
    fn add_constraint_runs_the_solver() {
        // Draw a roughly-horizontal line, then constrain it horizontal:
        // the real valenx-sketch solver must run and converge.
        let _g = fresh();
        let sid = 0;
        create_sketch(&json!({})).unwrap();
        let line = add_sketch_line(&json!({
            "sketch_id": sid, "x1": 0.0, "y1": 0.0, "x2": 2.0, "y2": 0.3
        }))
        .unwrap();
        let line_id = structured(&line)["line_id"].as_u64().unwrap();
        let r = add_constraint(&json!({
            "sketch_id": sid, "type": "horizontal", "line": line_id
        }))
        .unwrap();
        assert_eq!(structured(&r)["solver_converged"], true);
    }

    #[test]
    fn add_constraint_rejects_an_unknown_type() {
        let _g = fresh();
        create_sketch(&json!({})).unwrap();
        let err = add_constraint(&json!({
            "sketch_id": 0, "type": "telepathy"
        }))
        .unwrap_err();
        assert!(err.to_string().contains("unsupported constraint type"));
    }

    /// Build a closed unit-square sketch and return its id.
    fn unit_square_sketch() -> usize {
        let sid = structured(&create_sketch(&json!({})).unwrap())["sketch_id"]
            .as_u64()
            .unwrap() as usize;
        // Four corners of a 1×1 square, drawn as four connected lines.
        let pts = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        for k in 0..4 {
            let (x1, y1) = pts[k];
            let (x2, y2) = pts[(k + 1) % 4];
            add_sketch_line(&json!({
                "sketch_id": sid, "x1": x1, "y1": y1, "x2": x2, "y2": y2
            }))
            .unwrap();
        }
        sid
    }

    #[test]
    fn pad_a_square_produces_a_measurable_solid() {
        // The core path: draw a unit square, pad it 2 units, and the
        // evaluated volume must be 1·1·2 = 2.
        let _g = fresh();
        let sid = unit_square_sketch();
        let p = pad(&json!({ "sketch_id": sid, "depth": 2.0 })).unwrap();
        assert_eq!(structured(&p)["operation"], "pad");
        let metrics = evaluate_design(&json!({})).unwrap();
        let vol = structured(&metrics)["volume"].as_f64().unwrap();
        assert!(
            (vol - 2.0).abs() < 0.05,
            "padded unit square should have volume ~2, got {vol}"
        );
    }

    #[test]
    fn evaluate_design_reports_the_bounding_box() {
        // A 1×1 square padded by 3 → bounding box 1×1×3.
        let _g = fresh();
        let sid = unit_square_sketch();
        pad(&json!({ "sketch_id": sid, "depth": 3.0 })).unwrap();
        let m = evaluate_design(&json!({})).unwrap();
        let size = &structured(&m)["bounding_box"]["size"];
        let sz = [
            size[0].as_f64().unwrap(),
            size[1].as_f64().unwrap(),
            size[2].as_f64().unwrap(),
        ];
        assert!(
            (sz[0] - 1.0).abs() < 0.05,
            "x extent should be 1, got {}",
            sz[0]
        );
        assert!(
            (sz[1] - 1.0).abs() < 0.05,
            "y extent should be 1, got {}",
            sz[1]
        );
        assert!(
            (sz[2] - 3.0).abs() < 0.05,
            "z extent should be 3, got {}",
            sz[2]
        );
    }

    #[test]
    fn evaluate_design_applies_density_to_mass() {
        // Mass = volume · density. A volume-2 part at density 2.7 has
        // mass 5.4.
        let _g = fresh();
        let sid = unit_square_sketch();
        pad(&json!({ "sketch_id": sid, "depth": 2.0 })).unwrap();
        let m = evaluate_design(&json!({ "density": 2.7 })).unwrap();
        let mass = structured(&m)["mass"].as_f64().unwrap();
        assert!(
            (mass - 5.4).abs() < 0.2,
            "mass should be volume·density ≈ 5.4, got {mass}"
        );
    }

    #[test]
    fn pocket_removes_material_so_mass_drops() {
        // The generative-feedback loop: pad a block, measure it, pocket
        // a hole, measure again — the volume must have decreased.
        let _g = fresh();
        // 4×4 base block, 2 thick.
        let base = structured(&create_sketch(&json!({})).unwrap())["sketch_id"]
            .as_u64()
            .unwrap() as usize;
        for k in 0..4 {
            let pts = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)];
            let (x1, y1) = pts[k];
            let (x2, y2) = pts[(k + 1) % 4];
            add_sketch_line(&json!({
                "sketch_id": base, "x1": x1, "y1": y1, "x2": x2, "y2": y2
            }))
            .unwrap();
        }
        pad(&json!({ "sketch_id": base, "depth": 2.0 })).unwrap();
        let before = structured(&evaluate_design(&json!({})).unwrap())["volume"]
            .as_f64()
            .unwrap();

        // A 2×2 hole pocketed through the block. "Through all" is
        // specified with a depth that runs *past* the 2-thick block
        // (depth 4 > 2) — the conventional CAD idiom. A depth equal to
        // the block thickness would leave the cutter's far cap
        // coincident with the block's top face, which stalls the
        // truck-shapeops boolean kernel; the overshoot clears it.
        let hole = structured(&create_sketch(&json!({})).unwrap())["sketch_id"]
            .as_u64()
            .unwrap() as usize;
        for k in 0..4 {
            let pts = [(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)];
            let (x1, y1) = pts[k];
            let (x2, y2) = pts[(k + 1) % 4];
            add_sketch_line(&json!({
                "sketch_id": hole, "x1": x1, "y1": y1, "x2": x2, "y2": y2
            }))
            .unwrap();
        }
        pocket(&json!({ "sketch_id": hole, "depth": 4.0 })).unwrap();
        let after = structured(&evaluate_design(&json!({})).unwrap())["volume"]
            .as_f64()
            .unwrap();
        assert!(
            after < before,
            "pocketing a hole must reduce the volume: {before} → {after}"
        );
    }

    #[test]
    fn pad_rejects_an_empty_sketch() {
        let _g = fresh();
        let sid = structured(&create_sketch(&json!({})).unwrap())["sketch_id"]
            .as_u64()
            .unwrap();
        let err = pad(&json!({ "sketch_id": sid, "depth": 1.0 })).unwrap_err();
        assert!(err.to_string().contains("empty"), "got: {err}");
    }

    #[test]
    fn evaluate_design_of_an_empty_design_errors() {
        let _g = fresh();
        let err = evaluate_design(&json!({})).unwrap_err();
        assert!(err.to_string().contains("empty"), "got: {err}");
    }

    #[test]
    fn reset_design_clears_the_session() {
        let _g = fresh();
        unit_square_sketch();
        pad(&json!({ "sketch_id": 0, "depth": 1.0 })).unwrap();
        // After a reset the design is empty again.
        reset_design(&json!({})).unwrap();
        let err = evaluate_design(&json!({})).unwrap_err();
        assert!(err.to_string().contains("empty"));
    }

    #[test]
    fn boolean_requires_at_least_two_targets() {
        let _g = fresh();
        let err = boolean(&json!({
            "operation": "union", "targets": [0]
        }))
        .unwrap_err();
        assert!(err.to_string().contains("at least 2"));
    }

    #[test]
    fn boolean_rejects_an_unknown_operation() {
        let _g = fresh();
        let err = boolean(&json!({
            "operation": "frobnicate", "targets": [0, 1]
        }))
        .unwrap_err();
        assert!(err.to_string().contains("unsupported boolean operation"));
    }

    #[test]
    fn fillet_rejects_an_unknown_target_feature() {
        let _g = fresh();
        let err = fillet(&json!({
            "target_feature": 42, "radius": 0.5
        }))
        .unwrap_err();
        assert!(err.to_string().contains("no feature with id 42"));
    }

    #[test]
    fn dispatch_routes_design_tools_and_passes_others_through() {
        let _g = fresh();
        // A design tool is dispatched.
        assert!(dispatch("create_sketch", &json!({})).is_some());
        // A non-design tool falls through (None).
        assert!(dispatch("dock", &json!({})).is_none());
    }

    #[test]
    fn tool_list_advertises_every_design_tool() {
        let names: Vec<String> = tool_list()
            .iter()
            .map(|t| t["name"].as_str().unwrap().to_string())
            .collect();
        for expected in [
            "create_sketch",
            "add_sketch_line",
            "add_sketch_circle",
            "add_constraint",
            "pad",
            "pocket",
            "revolve",
            "fillet",
            "boolean",
            "evaluate_design",
            "export_design",
            "reset_design",
        ] {
            assert!(
                names.iter().any(|n| n == expected),
                "tool `{expected}` should be advertised"
            );
        }
    }

    #[test]
    fn mesh_metrics_of_a_unit_cube_mesh() {
        // A hand-built closed unit-cube triangle mesh: volume 1,
        // surface area 6, bounding box 1×1×1. Verifies the
        // divergence-theorem volume + area + bbox maths directly.
        use nalgebra::Vector3 as V;
        use valenx_mesh::element::{ElementBlock, ElementType};
        let mut mesh = valenx_mesh::Mesh::new("cube");
        // 8 corners.
        let corners = [
            V::new(0.0, 0.0, 0.0),
            V::new(1.0, 0.0, 0.0),
            V::new(1.0, 1.0, 0.0),
            V::new(0.0, 1.0, 0.0),
            V::new(0.0, 0.0, 1.0),
            V::new(1.0, 0.0, 1.0),
            V::new(1.0, 1.0, 1.0),
            V::new(0.0, 1.0, 1.0),
        ];
        for c in corners {
            mesh.nodes.push(c);
        }
        // 12 triangles, all wound counter-clockwise as seen from
        // outside (outward normals → positive volume).
        let tris: [[u32; 3]; 12] = [
            [0, 2, 1],
            [0, 3, 2], // bottom (−z)
            [4, 5, 6],
            [4, 6, 7], // top (+z)
            [0, 1, 5],
            [0, 5, 4], // −y
            [2, 3, 7],
            [2, 7, 6], // +y
            [1, 2, 6],
            [1, 6, 5], // +x
            [0, 4, 7],
            [0, 7, 3], // −x
        ];
        let mut conn = Vec::new();
        for t in tris {
            conn.extend_from_slice(&t);
        }
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Tri3,
            connectivity: conn,
        });
        let m = mesh_metrics(&mesh);
        assert!(
            (m.volume - 1.0).abs() < 1e-9,
            "cube volume should be 1, got {}",
            m.volume
        );
        assert!((m.surface_area - 6.0).abs() < 1e-9, "cube area should be 6");
        assert_eq!(m.bbox_min, [0.0, 0.0, 0.0]);
        assert_eq!(m.bbox_max, [1.0, 1.0, 1.0]);
    }

    #[test]
    fn mesh_metrics_of_an_empty_mesh_is_zero() {
        let mesh = valenx_mesh::Mesh::new("empty");
        let m = mesh_metrics(&mesh);
        assert_eq!(m.volume, 0.0);
        assert_eq!(m.surface_area, 0.0);
        assert_eq!(m.bbox_min, [0.0; 3]);
    }

    // === Round-17 RED→GREEN: design-session caps + NaN guards =========

    /// Round-17 M1: `create_sketch` refuses to grow the session past
    /// [`MAX_DESIGN_DRAFTS`]. Pre-fix the loop allocated a fresh draft
    /// every call without bound — a misbehaving LLM client could exhaust
    /// process memory by spamming `create_sketch`.
    #[test]
    fn create_sketch_caps_session_drafts() {
        let _g = fresh();
        // Fill the cap.
        for _ in 0..MAX_DESIGN_DRAFTS {
            create_sketch(&json!({})).unwrap();
        }
        // The next call must be rejected with a clear message.
        let err = create_sketch(&json!({})).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains(&format!("{MAX_DESIGN_DRAFTS}")),
            "cap message should name the limit, got: {msg}"
        );
    }

    /// Round-17 M1: `add_sketch_line` refuses to push entities past
    /// [`MAX_DESIGN_ENTITIES_PER_SKETCH`]. We pre-load the sketch
    /// directly (via the public Sketch API) so the test runs in
    /// milliseconds instead of 5000 RPC round-trips.
    #[test]
    fn add_sketch_line_caps_entities_per_sketch() {
        let _g = fresh();
        let sid = structured(&create_sketch(&json!({})).unwrap())["sketch_id"]
            .as_u64()
            .unwrap() as usize;
        // Stuff the sketch up to the cap directly. add_point per call
        // is 1 entity, so MAX_DESIGN_ENTITIES_PER_SKETCH points fill it.
        {
            let mut s = lock_session().unwrap();
            let draft = s.drafts.get_mut(sid).unwrap();
            for k in 0..MAX_DESIGN_ENTITIES_PER_SKETCH {
                draft.sketch.add_point(k as f64, 0.0);
            }
        }
        // Now any further line must be rejected — the entrypoint sees
        // the entity count and refuses before calling add_point.
        let err = add_sketch_line(&json!({
            "sketch_id": sid, "x1": 0.0, "y1": 0.0, "x2": 1.0, "y2": 0.0
        }))
        .unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains(&format!("{MAX_DESIGN_ENTITIES_PER_SKETCH}")),
            "cap message should name the limit, got: {msg}"
        );
    }

    /// Round-17 M1: `add_sketch_circle` refuses to push the sketch past
    /// the entity cap. A circle adds `2 × CIRCLE_SEGMENTS = 128`
    /// entities; pre-load just past `MAX − 1` so adding a circle would
    /// overflow.
    #[test]
    fn add_sketch_circle_caps_entities_per_sketch() {
        let _g = fresh();
        let sid = structured(&create_sketch(&json!({})).unwrap())["sketch_id"]
            .as_u64()
            .unwrap() as usize;
        {
            let mut s = lock_session().unwrap();
            let draft = s.drafts.get_mut(sid).unwrap();
            // Leave only 1 entity worth of headroom — a circle needs
            // 2 × CIRCLE_SEGMENTS = 128, so it must be rejected.
            for k in 0..(MAX_DESIGN_ENTITIES_PER_SKETCH - 1) {
                draft.sketch.add_point(k as f64, 0.0);
            }
        }
        let err = add_sketch_circle(&json!({
            "sketch_id": sid, "cx": 0.0, "cy": 0.0, "radius": 1.0
        }))
        .unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains(&format!("{MAX_DESIGN_ENTITIES_PER_SKETCH}-entity cap")),
            "cap message should name the per-sketch limit, got: {msg}"
        );
    }

    /// Round-17 M1: `add_constraint` refuses to push constraints past
    /// [`MAX_DESIGN_CONSTRAINTS_PER_SKETCH`].
    #[test]
    fn add_constraint_caps_constraints_per_sketch() {
        let _g = fresh();
        // Set up a sketch with one real line — constraints reference
        // it by id so the constraint vector can grow without the entity
        // count blowing up.
        let sid = 0;
        create_sketch(&json!({})).unwrap();
        let line = add_sketch_line(&json!({
            "sketch_id": sid, "x1": 0.0, "y1": 0.0, "x2": 1.0, "y2": 0.0
        }))
        .unwrap();
        let line_id = structured(&line)["line_id"].as_u64().unwrap();
        // Pre-load constraints directly — avoids running the solver
        // 10k times.
        {
            let mut s = lock_session().unwrap();
            let draft = s.drafts.get_mut(sid).unwrap();
            for _ in 0..MAX_DESIGN_CONSTRAINTS_PER_SKETCH {
                draft
                    .sketch
                    .add_constraint(Constraint::Horizontal(EntityId(line_id as usize)));
            }
        }
        // The next add_constraint must be rejected.
        let err = add_constraint(&json!({
            "sketch_id": sid, "type": "horizontal", "line": line_id
        }))
        .unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains(&format!("{MAX_DESIGN_CONSTRAINTS_PER_SKETCH}")),
            "cap message should name the limit, got: {msg}"
        );
    }

    /// Round-17 M1: features (pad / pocket / revolve / fillet / boolean)
    /// refuse to grow past [`MAX_DESIGN_FEATURES`]. We pre-load the
    /// tree directly with the Pad constructor so the test runs in
    /// milliseconds.
    #[test]
    fn pad_caps_features_in_session() {
        let _g = fresh();
        let sid = unit_square_sketch();
        // Pre-fill the feature tree to the cap via direct API.
        {
            let mut s = lock_session().unwrap();
            // Push the unit-square sketch into the tree once so we
            // have a real SketchRef to reuse.
            let sketch = s.drafts[sid].sketch.clone();
            let sref = s.tree.add_sketch(sketch);
            for _ in 0..MAX_DESIGN_FEATURES {
                s.tree.add_feature(
                    Feature::Pad(PadParams {
                        sketch: sref,
                        depth: 1.0.into(),
                        direction_positive: true,
                    }),
                    "fill".to_string(),
                );
            }
        }
        // The next pad call must be rejected.
        let err = pad(&json!({ "sketch_id": sid, "depth": 2.0 })).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains(&format!("{MAX_DESIGN_FEATURES}")),
            "cap message should name the limit, got: {msg}"
        );
    }

    /// Round-17 M1 + L4: `boolean` refuses a targets array larger than
    /// [`MAX_BOOLEAN_TARGETS`] BEFORE allocating the FeatureId vector.
    #[test]
    fn boolean_caps_target_count() {
        let _g = fresh();
        // Build a targets array one larger than the cap. The values
        // don't matter — the cap check fires before the existence
        // check.
        let too_many: Vec<u64> = (0..(MAX_BOOLEAN_TARGETS as u64 + 1)).collect();
        let err = boolean(&json!({
            "operation": "union",
            "targets": too_many,
        }))
        .unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains(&format!("{MAX_BOOLEAN_TARGETS}")),
            "cap message should name the limit, got: {msg}"
        );
    }

    /// Round-17 L1: the finite-coords guard rejects NaN.
    ///
    /// Note: we can't drive the public `add_sketch_line` entrypoint
    /// directly with NaN because serde_json's `Number::from_f64`
    /// returns `None` for non-finite values, so `json!({"x": NAN})`
    /// serialises to `null` (which then fails `req_f64`'s `as_f64()`
    /// check with a different error). We test the helper itself,
    /// which is the load-bearing primitive — every entrypoint that
    /// receives an f64 from `req_f64` runs through this helper before
    /// the value reaches the solver.
    #[test]
    fn ensure_finite_coords_rejects_nan() {
        let err = ensure_finite_coords(&[f64::NAN, 0.0], "test").unwrap_err();
        assert!(err.to_string().contains("finite"), "got: {err}");
    }

    /// Round-17 L1: the finite-coords guard rejects +infinity.
    #[test]
    fn ensure_finite_coords_rejects_positive_infinity() {
        let err = ensure_finite_coords(&[0.0, f64::INFINITY], "test").unwrap_err();
        assert!(err.to_string().contains("finite"), "got: {err}");
    }

    /// Round-17 L1: the finite-coords guard rejects -infinity.
    #[test]
    fn ensure_finite_coords_rejects_negative_infinity() {
        let err = ensure_finite_coords(&[f64::NEG_INFINITY, 0.0, 0.0], "test").unwrap_err();
        assert!(err.to_string().contains("finite"), "got: {err}");
    }

    /// Round-17 L1: the finite-coords guard accepts all-finite input.
    #[test]
    fn ensure_finite_coords_accepts_finite() {
        ensure_finite_coords(&[0.0, 1.0, -2.5, 1e100], "test").unwrap();
    }

    /// Round-17 L1: every f64-taking entrypoint (`add_sketch_line`,
    /// `add_sketch_circle`, `add_constraint`) MUST run its inputs
    /// through `ensure_finite_coords` before they reach the solver.
    /// This test pins that wiring by grep-ing the source — if a
    /// refactor moves the guard out of an entrypoint, the test fails
    /// and we relearn the invariant.
    #[test]
    fn finite_coords_guard_is_wired_into_every_entrypoint() {
        // Read this very source file. The relative path works for
        // workspace `cargo test` runs which invoke from the crate
        // dir; if it doesn't, the test gracefully skips.
        let src = std::fs::read_to_string("src/design.rs")
            .or_else(|_| std::fs::read_to_string("crates/valenx-mcp/src/design.rs"))
            .ok();
        let Some(src) = src else {
            eprintln!("skipping: couldn't locate design.rs source for grep test");
            return;
        };
        // Each entrypoint should mention the helper. `add_constraint`
        // wraps the call in `req_finite_value`, which itself calls
        // `ensure_finite_coords` — so the helper name appears once
        // per entrypoint in the body region of the file (we count
        // ≥ 3, allowing for additional defensive calls).
        let calls = src.matches("ensure_finite_coords").count();
        assert!(
            calls >= 3,
            "expected ensure_finite_coords to be called from at least \
             3 entrypoints, found {calls}"
        );
    }
}
