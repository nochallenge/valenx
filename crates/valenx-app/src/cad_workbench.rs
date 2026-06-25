//! Parametric-CAD workbench — named parameters drive sketch geometry *and* a
//! CSG feature tree.
//!
//! A right-side panel over `valenx-solvespace-3d`. One shared table of **named
//! parameters** (Fusion's "Change Parameters") feeds two consumers:
//!
//! 1. **Sketch** — pick a parameter to drive a circle's radius and Solve; the
//!    constraint solver lands the circle on the parameter-driven radius.
//! 2. **Feature tree (CSG)** — an ordered list of steps, each placing a
//!    primitive (box / cylinder) and combining it with the running body via
//!    New / Join / Cut / Intersect. Rebuild folds the tree into one solid,
//!    tessellates it, and pushes it into the central 3-D viewport. Edit a
//!    parameter, rebuild, and the whole model re-drives — a hole moves, a boss
//!    grows.
//!
//! Compute is synchronous: parameter resolution is sub-millisecond and a
//! handful of boolean ops on primitives is well under a frame.

use eframe::egui;

use valenx_solvespace_3d::{
    timeline_from_ron, timeline_to_ron, Constraint3D, Feature, FeatureTimeline, Op, ParameterTable,
    Sketch3D, Step,
};

use crate::agent_commands::AgentValue;
use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which drawing tool the sketch canvas is in. `Line` drops straight segments
/// (one click each); `Arc` builds a 3-point circular arc (start, end, a point on
/// the arc); `Spline` builds a smooth Catmull-Rom curve through clicked points.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
enum SketchTool {
    #[default]
    Line,
    Arc,
    Spline,
}

/// One ordered edge of the drawn sketch profile. Each segment knows its own end
/// anchor (`to`) so the profile is an ordered chain `start → seg₀.to → seg₁.to …`.
/// Tessellated by [`CadWorkbenchState::sketch_profile_polyline`] into the closed
/// polyline that [`CadWorkbenchState::add_extrude_from_sketch`] sweeps into a solid.
#[derive(Clone, Debug, PartialEq)]
enum SketchSeg {
    /// Straight segment from the previous anchor to `to`.
    Line { to: [f64; 2] },
    /// Circular arc from the previous anchor to `to`, passing through `via` (a
    /// third point ON the arc). The circle is the unique one through the three
    /// points; if they are collinear it degenerates to a straight line.
    Arc { to: [f64; 2], via: [f64; 2] },
    /// Smooth Catmull-Rom spline from the previous anchor through every point in
    /// `thru` (the last of which is this segment's end anchor).
    Spline { thru: Vec<[f64; 2]> },
}

impl SketchSeg {
    /// The end anchor of this segment (where the next segment starts).
    fn end(&self) -> [f64; 2] {
        match self {
            SketchSeg::Line { to } | SketchSeg::Arc { to, .. } => *to,
            // A Spline always carries ≥1 `thru` point; its end is the last one.
            SketchSeg::Spline { thru } => *thru.last().expect("spline has ≥1 point"),
        }
    }
}

/// Which primitive a feature-tree step builds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FeatureKind {
    Box,
    Cylinder,
    Sphere,
    Cone,
    Torus,
    Extrude,
    Revolve,
}

/// One UI-editable feature-tree step. Carries both the box and cylinder
/// dimension fields so toggling `kind` preserves whatever the user typed.
#[derive(Clone)]
struct UiStep {
    op: Op,
    kind: FeatureKind,
    /// Extrude profile — `(x, y)` points (literal coords, not expressions). For
    /// a Revolve step the same field holds the `(r, z)` half-section.
    profile: Vec<(f64, f64)>,
    dx: String,
    dy: String,
    dz: String,
    radius: String,
    height: String,
    top_radius: String,
    major: String,
    minor: String,
    /// Revolve sweep angle in degrees (expression). `360` = full lathe body.
    angle: String,
    x: String,
    y: String,
    z: String,
    rx: String,
    ry: String,
    rz: String,
}

impl UiStep {
    /// A step with default dimension fields for every primitive; callers
    /// override the ones the chosen `kind` actually uses.
    fn base(op: Op, kind: FeatureKind) -> Self {
        Self {
            op,
            kind,
            profile: vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
            dx: "1".into(),
            dy: "1".into(),
            dz: "1".into(),
            radius: "0.5".into(),
            height: "1".into(),
            top_radius: "0".into(),
            major: "1".into(),
            minor: "0.25".into(),
            angle: "360".into(),
            x: "0".into(),
            y: "0".into(),
            z: "0".into(),
            rx: "0".into(),
            ry: "0".into(),
            rz: "0".into(),
        }
    }

    fn new_box() -> Self {
        Self::base(Op::Join, FeatureKind::Box)
    }

    fn new_cylinder() -> Self {
        let mut s = Self::base(Op::Cut, FeatureKind::Cylinder);
        s.radius = "0.25".into();
        s.height = "2".into();
        s
    }

    fn new_sphere() -> Self {
        Self::base(Op::Join, FeatureKind::Sphere)
    }

    fn new_cone() -> Self {
        Self::base(Op::Join, FeatureKind::Cone)
    }

    fn new_torus() -> Self {
        Self::base(Op::Join, FeatureKind::Torus)
    }

    fn new_extrude() -> Self {
        Self::base(Op::Join, FeatureKind::Extrude)
    }

    fn new_revolve() -> Self {
        let mut s = Self::base(Op::Join, FeatureKind::Revolve);
        // Default half-section is a cone profile that starts and ends on the
        // axis (r = 0) so a full-turn revolve caps cleanly: axis → outer → axis.
        s.profile = vec![(0.0, 0.0), (1.0, 0.0), (0.0, 2.0)];
        s.angle = "360".into();
        s
    }

    /// Translate into a solver-crate [`Step`].
    fn to_step(&self) -> Step {
        let feature = match self.kind {
            FeatureKind::Box => Feature::Box {
                dx: self.dx.clone(),
                dy: self.dy.clone(),
                dz: self.dz.clone(),
            },
            FeatureKind::Cylinder => Feature::Cylinder {
                radius: self.radius.clone(),
                height: self.height.clone(),
            },
            FeatureKind::Sphere => Feature::Sphere {
                radius: self.radius.clone(),
            },
            FeatureKind::Cone => Feature::Cone {
                base_radius: self.radius.clone(),
                top_radius: self.top_radius.clone(),
                height: self.height.clone(),
            },
            FeatureKind::Torus => Feature::Torus {
                major_radius: self.major.clone(),
                minor_radius: self.minor.clone(),
            },
            FeatureKind::Extrude => Feature::Extrude {
                profile: self.profile.clone(),
                height: self.height.clone(),
            },
            FeatureKind::Revolve => Feature::Revolve {
                profile: self.profile.clone(),
                angle_deg: self.angle.clone(),
            },
        };
        let mut step = Step::placed(
            self.op,
            feature,
            self.x.clone(),
            self.y.clone(),
            self.z.clone(),
        );
        step.rotate_deg = [self.rx.clone(), self.ry.clone(), self.rz.clone()];
        step
    }
}

/// The default feature tree: a unit box with a cylinder punched through it —
/// `valenx-cad`'s proven "punched cube" geometry, so the seed always rebuilds.
fn default_steps() -> Vec<UiStep> {
    let mut base = UiStep::base(Op::New, FeatureKind::Box);
    base.dx = "size".into();
    base.dy = "size".into();
    base.dz = "size".into();
    let mut hole = UiStep::base(Op::Cut, FeatureKind::Cylinder);
    hole.radius = "hole_r".into();
    hole.height = "hole_h".into();
    hole.x = "size / 2".into();
    hole.y = "size / 2".into();
    hole.z = "-0.5".into();
    vec![base, hole]
}

/// Persistent state for the parametric-CAD workbench.
pub struct CadWorkbenchState {
    /// Editable named parameters as (name, expression) rows — shared by the
    /// sketch demo and the feature tree.
    params: Vec<(String, String)>,
    /// Name of the parameter that drives the circle's radius.
    radius_param: String,
    results: Option<CadResults>,
    /// Feature-tree steps, in build order.
    steps: Vec<UiStep>,
    /// Last rebuild outcome: `Ok(status)` or `Err(message)`.
    tree_status: Option<Result<String, String>>,
    /// Tessellated body waiting to be pushed into the viewport (deferred out
    /// of the panel borrow).
    rebuilt_mesh: Option<valenx_mesh::Mesh>,
    /// Set when a fresh rebuild needs pushing into the viewport.
    push_rebuild: bool,
    /// Snapshots from the last rebuild — `history[k]` is the set of bodies
    /// present after step k. Drives the scrubber.
    history: Option<Vec<Vec<valenx_cad::Solid>>>,
    /// 1-based step the history scrubber is showing (`1..=history.len()`).
    scrub: usize,
    /// Per-body visibility for the final body set (multi-body). Empty = all
    /// visible; `body_visible[i] == false` hides body i in the viewport.
    body_visible: Vec<bool>,
    /// Material density (mass per unit volume) — drives the mass readout.
    density: f64,
    /// Set by [`CadWorkbenchState::request_rebuild`] (e.g. the top-bar Part
    /// Design menu) to ask `draw_cad_workbench` to rebuild the tree into the
    /// viewport on this frame — the programmatic equivalent of clicking the
    /// panel's "Rebuild → viewport" button.
    rebuild_request: bool,
    /// First anchor of the mouse-drawn sketch profile (model units), or `None`
    /// when the sketch is empty. The profile is the ordered chain
    /// `sketch_start → sketch_segs[0].end → sketch_segs[1].end → …`.
    sketch_start: Option<[f64; 2]>,
    /// Ordered edges of the mouse-drawn profile (Line / Arc / Spline), in draw
    /// order. Feeds [`CadWorkbenchState::sketch_profile_polyline`].
    sketch_segs: Vec<SketchSeg>,
    /// Which drawing tool clicks add (Line / Arc / Spline).
    sketch_tool: SketchTool,
    /// Clicks staged for a multi-click segment: 3 for an Arc (start, end, via),
    /// ≥1 for a Spline run (committed on **Finish curve**). Empty between segments.
    sketch_pending: Vec<[f64; 2]>,
    /// Whether the sketched profile's loop is closed (the user clicked back near
    /// the first anchor). A closed loop renders a translucent fill.
    sketch_closed: bool,
    /// Snap canvas clicks to the [`SKETCH_SNAP`]-unit grid.
    sketch_grid_snap: bool,
    /// Extrude depth (model units) the **Extrude sketch** button applies to the
    /// drawn profile.
    sketch_extrude_height: f64,
    /// Set by the top-bar Part Design **Sketch** menu item to ask the panel to
    /// scroll the sketch canvas into view and flash its header on this frame.
    sketch_focus_request: bool,
}

/// Grid spacing of the sketch canvas, in model units — also the snap step when
/// **snap-to-grid** is enabled.
const SKETCH_SNAP: f64 = 0.25;

/// Half-extent of the sketch canvas in model units: the canvas shows
/// `[-SKETCH_VIEW, +SKETCH_VIEW]` on each axis (a `2·SKETCH_VIEW`-unit window).
const SKETCH_VIEW: f64 = 4.0;

/// Click/close pick radius on the sketch canvas, in screen pixels — a click
/// within this distance of the first vertex closes the polygon loop.
const SKETCH_PICK_PX: f32 = 8.0;

/// Number of straight chords each circular Arc segment is sampled into (both for
/// the live canvas render and the extruded profile). 24 chords keeps the
/// faceting under ~0.5° per chord on a quarter-turn arc — visually smooth and a
/// faithful extrusion cross-section. The arc's endpoints are always included.
const SKETCH_ARC_SAMPLES: usize = 24;

/// Number of straight chords sampled **per Catmull-Rom span** (one span between
/// each pair of consecutive control points). The spline's control points are
/// always included, so a curve through N points yields N + (N−1)·(SAMPLES−1)
/// polyline vertices.
const SKETCH_SPLINE_SAMPLES: usize = 16;

impl Default for CadWorkbenchState {
    fn default() -> Self {
        Self {
            params: vec![
                ("base".to_string(), "4".to_string()),
                ("radius".to_string(), "base + 1".to_string()),
                ("size".to_string(), "1".to_string()),
                ("hole_r".to_string(), "0.25".to_string()),
                ("hole_h".to_string(), "2".to_string()),
            ],
            radius_param: "radius".to_string(),
            results: None,
            steps: default_steps(),
            tree_status: None,
            rebuilt_mesh: None,
            push_rebuild: false,
            history: None,
            scrub: 1,
            body_visible: Vec::new(),
            density: 1.0,
            rebuild_request: false,
            sketch_start: None,
            sketch_segs: Vec::new(),
            sketch_tool: SketchTool::default(),
            sketch_pending: Vec::new(),
            sketch_closed: false,
            sketch_grid_snap: true,
            sketch_extrude_height: 1.0,
            sketch_focus_request: false,
        }
    }
}

/// Solve the centre and radius of the circle through three points. Returns
/// `None` when the points are (near-)collinear — the circumradius diverges, so
/// the caller falls back to a straight chord. Closed form from the perpendicular
/// bisectors: with `d = 2·(aₓ(bᵧ−cᵧ)+bₓ(cᵧ−aᵧ)+cₓ(aᵧ−bᵧ))`, the centre is the
/// standard circumcentre and `r = |centre − a|`.
fn circle_through_3(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Option<([f64; 2], f64)> {
    let d = 2.0 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
    if d.abs() < 1e-9 {
        return None; // collinear → no finite circle
    }
    let a2 = a[0] * a[0] + a[1] * a[1];
    let b2 = b[0] * b[0] + b[1] * b[1];
    let c2 = c[0] * c[0] + c[1] * c[1];
    let ux = (a2 * (b[1] - c[1]) + b2 * (c[1] - a[1]) + c2 * (a[1] - b[1])) / d;
    let uy = (a2 * (c[0] - b[0]) + b2 * (a[0] - c[0]) + c2 * (b[0] - a[0])) / d;
    let centre = [ux, uy];
    let r = ((a[0] - ux).powi(2) + (a[1] - uy).powi(2)).sqrt();
    Some((centre, r))
}

/// Sample the circular arc that starts at `start`, ends at `end`, and passes
/// through `via`, into a polyline. The first point is `start`; the returned
/// vector ends just *before* `end` (the caller appends the shared end anchor, so
/// chained segments don't duplicate vertices). The swept direction and span are
/// chosen so the arc actually passes through `via` (i.e. the minor-or-major arc
/// that contains the via point). Collinear points fall back to the single
/// straight point `[start]` (a plain chord to `end`).
fn sample_arc(start: [f64; 2], via: [f64; 2], end: [f64; 2]) -> Vec<[f64; 2]> {
    let Some((centre, r)) = circle_through_3(start, via, end) else {
        return vec![start]; // collinear → straight chord start→end
    };
    let ang = |p: [f64; 2]| (p[1] - centre[1]).atan2(p[0] - centre[0]);
    let a0 = ang(start);
    let a1 = ang(end);
    let av = ang(via);
    let tau = std::f64::consts::TAU;
    // Normalise a positive (CCW) sweep start→end, and the via offset within it.
    let mut sweep = a1 - a0;
    while sweep <= 0.0 {
        sweep += tau;
    }
    let mut via_off = av - a0;
    while via_off <= 0.0 {
        via_off += tau;
    }
    // If `via` is not inside the CCW start→end arc, the arc through it is the
    // complementary one — sweep clockwise (negative) instead.
    if via_off > sweep {
        sweep -= tau;
    }
    let mut out = Vec::with_capacity(SKETCH_ARC_SAMPLES);
    // Sample [0, SAMPLES): includes start (t=0), excludes end (t=1).
    for i in 0..SKETCH_ARC_SAMPLES {
        let t = i as f64 / SKETCH_ARC_SAMPLES as f64;
        let a = a0 + sweep * t;
        out.push([centre[0] + r * a.cos(), centre[1] + r * a.sin()]);
    }
    out
}

/// One uniform (centripetal-free, standard) Catmull-Rom interpolation at
/// parameter `t ∈ [0,1]` on the span `p1 → p2`, using neighbours `p0`/`p3` for
/// the tangents. The curve passes through `p1` at `t=0` and `p2` at `t=1`.
fn catmull_rom(p0: [f64; 2], p1: [f64; 2], p2: [f64; 2], p3: [f64; 2], t: f64) -> [f64; 2] {
    let t2 = t * t;
    let t3 = t2 * t;
    let comp = |a: f64, b: f64, c: f64, d: f64| {
        0.5 * ((2.0 * b)
            + (-a + c) * t
            + (2.0 * a - 5.0 * b + 4.0 * c - d) * t2
            + (-a + 3.0 * b - 3.0 * c + d) * t3)
    };
    [
        comp(p0[0], p1[0], p2[0], p3[0]),
        comp(p0[1], p1[1], p2[1], p3[1]),
    ]
}

/// Sample a Catmull-Rom spline through `ctrl` (≥2 control points) into a
/// polyline. The first point is `ctrl[0]`; the returned vector ends just
/// *before* the final control point (the caller appends the shared end anchor).
/// End tangents use a reflected phantom point so the curve still passes through
/// the endpoints. A single control point returns `[ctrl[0]]`.
fn sample_spline(ctrl: &[[f64; 2]]) -> Vec<[f64; 2]> {
    if ctrl.len() < 2 {
        return ctrl.first().map(|&p| vec![p]).unwrap_or_default();
    }
    let n = ctrl.len();
    // Phantom endpoints (reflect the neighbour about the endpoint) so the spline
    // interpolates ctrl[0] and ctrl[n-1] cleanly.
    let pt = |i: isize| -> [f64; 2] {
        if i < 0 {
            [2.0 * ctrl[0][0] - ctrl[1][0], 2.0 * ctrl[0][1] - ctrl[1][1]]
        } else if i as usize >= n {
            [
                2.0 * ctrl[n - 1][0] - ctrl[n - 2][0],
                2.0 * ctrl[n - 1][1] - ctrl[n - 2][1],
            ]
        } else {
            ctrl[i as usize]
        }
    };
    let mut out = Vec::with_capacity((n - 1) * SKETCH_SPLINE_SAMPLES);
    // One span per consecutive control pair; sample [0,SAMPLES) so spans join
    // without duplicating the shared control point. The final endpoint is left
    // for the caller to append.
    for k in 0..n - 1 {
        let p0 = pt(k as isize - 1);
        let p1 = pt(k as isize);
        let p2 = pt(k as isize + 1);
        let p3 = pt(k as isize + 2);
        for i in 0..SKETCH_SPLINE_SAMPLES {
            let t = i as f64 / SKETCH_SPLINE_SAMPLES as f64;
            out.push(catmull_rom(p0, p1, p2, p3, t));
        }
    }
    out
}

/// Public feature-tree mutators — the **single code path** shared by the
/// workbench's own `+ Box / + Cylinder / …` buttons and the top-bar
/// **Part Design** menu, so a feature added from the menu behaves exactly like
/// one added from the panel. Each `add_*` appends a step; [`Self::set_last_op`]
/// sets the boolean op of the last step; [`Self::request_rebuild`] flags a
/// rebuild-to-viewport on the next frame (mirroring the panel's "Rebuild →
/// viewport" button); [`Self::reset_to_base_solid`] starts a fresh single-solid
/// part.
impl CadWorkbenchState {
    /// Append a Box feature (the panel's `+ Box`).
    pub fn add_box(&mut self) {
        self.steps.push(UiStep::new_box());
    }

    /// Append a Cylinder feature (the panel's `+ Cylinder`).
    pub fn add_cylinder(&mut self) {
        self.steps.push(UiStep::new_cylinder());
    }

    /// Append a Sphere feature (the panel's `+ Sphere`).
    pub fn add_sphere(&mut self) {
        self.steps.push(UiStep::new_sphere());
    }

    /// Append a Cone feature (the panel's `+ Cone`).
    pub fn add_cone(&mut self) {
        self.steps.push(UiStep::new_cone());
    }

    /// Append a Torus feature (the panel's `+ Torus`).
    pub fn add_torus(&mut self) {
        self.steps.push(UiStep::new_torus());
    }

    /// Append an Extrude feature — a 2-D sketch profile swept along +Z (the
    /// panel's `+ Extrude`). This is the discoverable *sketch → extrude → solid*
    /// path: the step's editable `(x, y)` profile **is** the sketch.
    pub fn add_extrude(&mut self) {
        self.steps.push(UiStep::new_extrude());
    }

    /// Append a Revolve feature — a half-section profile revolved about the Z
    /// axis (the panel's `+ Revolve`), the lathe/turn operation.
    pub fn add_revolve(&mut self) {
        self.steps.push(UiStep::new_revolve());
    }

    /// Append an Extrude feature whose `(x, y)` profile is the supplied polygon
    /// `points` (in model units), swept along +Z by `height`. This is the
    /// **draw-then-extrude** path: the mouse-drawn sketch canvas hands its
    /// vertices straight to the feature tree. The op is `New` when the tree is
    /// empty (the first step must start a body) and `Join` otherwise (weld the
    /// prism onto the running body, matching the panel's `+ Extrude`). A profile
    /// of fewer than 3 points or a non-finite `height` is ignored. Caller should
    /// follow with [`Self::request_rebuild`] (or `perform_rebuild`) to fold the
    /// new prism into the viewport.
    pub fn add_extrude_from_sketch(&mut self, points: &[[f64; 2]], height: f64) {
        if points.len() < 3 || !height.is_finite() {
            return;
        }
        let op = if self.steps.is_empty() {
            Op::New
        } else {
            Op::Join
        };
        let mut st = UiStep::base(op, FeatureKind::Extrude);
        st.profile = points.iter().map(|p| (p[0], p[1])).collect();
        // Literal extrude depth — the timeline resolves the height as an
        // expression, so a plain number string drives it parametrically-compatibly.
        st.height = fmt_num(height);
        self.steps.push(st);
    }

    /// Append a straight-line vertex (Line tool, and the canvas's plain click).
    /// The first vertex becomes `sketch_start`; each later vertex appends a
    /// `Line` segment from the previous anchor. A no-op once the loop is closed.
    pub fn sketch_add_point(&mut self, x: f64, y: f64) {
        if self.sketch_closed {
            return;
        }
        let p = [x, y];
        match self.sketch_start {
            None => self.sketch_start = Some(p),
            Some(_) => self.sketch_segs.push(SketchSeg::Line { to: p }),
        }
    }

    /// Commit a 3-point circular arc (Arc tool): `start → end` passing through
    /// `via`. The previous anchor is the implicit start of the first arc; if the
    /// sketch is empty the supplied `start` seeds `sketch_start`. A no-op once the
    /// loop is closed.
    pub fn sketch_add_arc(&mut self, start: [f64; 2], via: [f64; 2], end: [f64; 2]) {
        if self.sketch_closed {
            return;
        }
        if self.sketch_start.is_none() {
            self.sketch_start = Some(start);
        }
        self.sketch_segs.push(SketchSeg::Arc { to: end, via });
    }

    /// Commit a Catmull-Rom spline (Spline tool) through `pts`. The previous
    /// anchor is the implicit start of the curve; if the sketch is empty `pts[0]`
    /// seeds `sketch_start` and the spline runs through the remaining points. A
    /// no-op once the loop is closed, or when fewer than 2 usable points result.
    pub fn sketch_add_spline(&mut self, pts: &[[f64; 2]]) {
        if self.sketch_closed || pts.is_empty() {
            return;
        }
        let thru: Vec<[f64; 2]> = if self.sketch_start.is_none() {
            self.sketch_start = Some(pts[0]);
            pts[1..].to_vec()
        } else {
            pts.to_vec()
        };
        if thru.is_empty() {
            return; // a lone seed point is just the start anchor, no curve yet
        }
        self.sketch_segs.push(SketchSeg::Spline { thru });
    }

    /// Remove the last sketched segment (canvas **Undo**). Dropping the only
    /// anchor empties the sketch. Always re-opens a closed loop.
    pub fn sketch_undo(&mut self) {
        if self.sketch_segs.pop().is_none() {
            self.sketch_start = None;
        }
        self.sketch_pending.clear();
        self.sketch_closed = false;
    }

    /// Discard the whole mouse-drawn sketch (canvas **Clear sketch**).
    pub fn sketch_clear(&mut self) {
        self.sketch_start = None;
        self.sketch_segs.clear();
        self.sketch_pending.clear();
        self.sketch_closed = false;
    }

    /// The ordered **anchor** vertices of the drawn profile (model units): the
    /// start anchor followed by each segment's end anchor. For a Line-only sketch
    /// these are exactly the clicked points; for Arc/Spline they are the segment
    /// endpoints (the curve interior lives between them — see
    /// [`Self::sketch_profile_polyline`]).
    pub fn sketch_points(&self) -> Vec<[f64; 2]> {
        let mut v = Vec::with_capacity(1 + self.sketch_segs.len());
        if let Some(p) = self.sketch_start {
            v.push(p);
            v.extend(self.sketch_segs.iter().map(|s| s.end()));
        }
        v
    }

    /// Number of anchor vertices currently in the sketch (start + one per
    /// segment). `0` when empty.
    pub fn sketch_anchor_count(&self) -> usize {
        if self.sketch_start.is_some() {
            1 + self.sketch_segs.len()
        } else {
            0
        }
    }

    /// Whether the mouse-drawn sketch loop is closed.
    pub fn sketch_is_closed(&self) -> bool {
        self.sketch_closed
    }

    /// Walk the segment chain and produce the closed **polyline** (model units)
    /// that [`Self::add_extrude_from_sketch`] sweeps into a solid: line endpoints
    /// are kept as-is, each Arc is sampled into `SKETCH_ARC_SAMPLES` chords
    /// along the circle through its 3 points, and each Spline into
    /// `SKETCH_SPLINE_SAMPLES` chords per control span. The result starts at the
    /// start anchor and ends at the last segment's end anchor (the extruder closes
    /// the loop back to the start). Returns an empty vector for an empty sketch.
    pub fn sketch_profile_polyline(&self) -> Vec<[f64; 2]> {
        let Some(start) = self.sketch_start else {
            return Vec::new();
        };
        let mut poly = vec![start];
        let mut cur = start;
        for seg in &self.sketch_segs {
            match seg {
                SketchSeg::Line { to } => {
                    poly.push(*to);
                    cur = *to;
                }
                SketchSeg::Arc { to, via } => {
                    // sample_arc yields [start .. just-before-end]; skip its
                    // first point (== cur, already in poly), then add the end.
                    let arc = sample_arc(cur, *via, *to);
                    poly.extend(arc.into_iter().skip(1));
                    poly.push(*to);
                    cur = *to;
                }
                SketchSeg::Spline { thru } => {
                    // Control points = previous anchor + the run-through points.
                    let mut ctrl = Vec::with_capacity(1 + thru.len());
                    ctrl.push(cur);
                    ctrl.extend_from_slice(thru);
                    let end = *thru.last().expect("spline thru is non-empty");
                    let curve = sample_spline(&ctrl);
                    poly.extend(curve.into_iter().skip(1)); // skip cur (already in poly)
                    poly.push(end);
                    cur = end;
                }
            }
        }
        poly
    }

    /// Ask the panel to scroll the sketch canvas into view (the top-bar Part
    /// Design → **Sketch** menu item). `draw_cad_workbench` consumes the flag.
    pub fn focus_sketch(&mut self) {
        self.sketch_focus_request = true;
    }

    /// Set the boolean [`Op`] of the **last** step (the next feature combines
    /// with the running body this way). No-op on an empty tree. Used by the
    /// Part Design → Boolean submenu (Union / Cut / Intersect).
    pub fn set_last_op(&mut self, op: Op) {
        if let Some(last) = self.steps.last_mut() {
            last.op = op;
        }
    }

    /// Part Design → Boolean → Union: set the last step to weld (boolean union).
    pub fn set_last_op_union(&mut self) {
        self.set_last_op(Op::Join);
    }

    /// Part Design → Boolean → Cut: set the last step to subtract (difference).
    pub fn set_last_op_cut(&mut self) {
        self.set_last_op(Op::Cut);
    }

    /// Part Design → Boolean → Intersect: set the last step to keep the overlap.
    pub fn set_last_op_intersect(&mut self) {
        self.set_last_op(Op::Intersect);
    }

    /// Whether the feature tree has at least one step (so the Boolean submenu
    /// can disable itself when there's nothing to combine).
    pub fn has_steps(&self) -> bool {
        !self.steps.is_empty()
    }

    /// Reset the feature tree to a single base solid (a unit box) — Part Design
    /// → "New part". Clears any cached rebuild/history so the viewport reflects
    /// the fresh part on the next rebuild.
    pub fn reset_to_base_solid(&mut self) {
        let mut base = UiStep::base(Op::New, FeatureKind::Box);
        base.dx = "size".into();
        base.dy = "size".into();
        base.dz = "size".into();
        self.steps = vec![base];
        self.history = None;
        self.rebuilt_mesh = None;
        self.body_visible = Vec::new();
        self.tree_status = None;
    }

    /// Flag that the tree should be rebuilt and pushed into the central 3-D
    /// viewport on the next frame — the same effect as the panel's
    /// "Rebuild → viewport" button. `draw_cad_workbench` performs the rebuild.
    pub fn request_rebuild(&mut self) {
        self.rebuild_request = true;
    }

    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`
    /// so an agent can discover the name space. Matches the workspace-wide
    /// `fn() -> &'static [&'static str]` contract (a fixed, discoverable list).
    ///
    /// Two groups:
    /// * **Singletons** — the always-present model / sketch controls
    ///   (`Material density`, `Extrude height`, `circle radius =`,
    ///   `Snap to grid`, `Sketch tool`).
    /// * **Per feature-tree step** — the feature tree is a `Vec<UiStep>`, so its
    ///   step fields are addressed by a **1-based index**: `step {N} <field>`.
    ///   Because the list must be `&'static`, the canonical entries below cover
    ///   the representative **first** step (`step 1 …`) across *all* per-kind
    ///   dimension captions (`dx`/`dy`/`dz`, `r`/`h`, `base r`/`top r`,
    ///   `major`/`minor`, `height`, `angle°`) plus `op` / `shape` and the
    ///   always-present placement (`x`/`y`/`z`) and rotation (`rx°`/`ry°`/`rz°`).
    ///   [`agent_set`](Self::agent_set) generalises the **same** field names to
    ///   any existing step — `step 2 dx`, `step 3 angle°`, … — and validates the
    ///   index against the live tree.
    ///
    /// NOTE — live-canvas sketch point placement (`sketch_start` /
    /// `sketch_segs` / `sketch_pending`) is intentionally **not** here: it is
    /// already bridge-driven via the dedicated `add_sketch_point` path, so a
    /// settable caption would duplicate that channel. Pure actions (Rebuild,
    /// Undo, + Box, …) and read-only outputs (resolved values, mass, the
    /// per-Extrude/Revolve `P{j}` profile points) are out of scope for
    /// `SetControl`.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            // -- Singletons --
            "Material density",
            "Extrude height",
            "circle radius =",
            "Snap to grid",
            "Sketch tool",
            // -- Representative first step (the same field names work for any
            //    `step {N}`; see agent_set) --
            "step 1 op",
            "step 1 shape",
            "step 1 dx",
            "step 1 dy",
            "step 1 dz",
            "step 1 r",
            "step 1 h",
            "step 1 base r",
            "step 1 top r",
            "step 1 major",
            "step 1 minor",
            "step 1 height",
            "step 1 angle\u{00B0}",
            "step 1 x",
            "step 1 y",
            "step 1 z",
            "step 1 rx\u{00B0}",
            "step 1 ry\u{00B0}",
            "step 1 rz\u{00B0}",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Caption strings match exactly what the workbench
    /// draws (and what each control is `labelled_by`), so an agent can set a
    /// parameter by the same name a user reads.
    ///
    /// Fail-loud: an unknown caption or a value of the wrong type returns
    /// `Err(String)` (the bridge turns it into a `warn` feed note) — never a
    /// panic, and no field is written on error. Numeric singletons read
    /// [`AgentValue::as_f64`]; the `Snap to grid` checkbox reads
    /// [`AgentValue::as_bool`]; the enum / expression-string controls (the
    /// `Sketch tool`, the parameter name in `circle radius =`, and every per-step
    /// `op` / `shape` / dimension field) read [`AgentValue::as_str`]. Per-step
    /// dimension fields are CAD **expressions** (e.g. `"size / 2"`, `"hole_r"`),
    /// so they are stored verbatim as strings — a bare number is accepted as the
    /// string form of that number.
    pub fn agent_set(&mut self, name: &str, value: &AgentValue) -> Result<(), String> {
        // A per-step dimension / placement / rotation field is a CAD
        // **expression** stored verbatim as a string (e.g. `"size / 2"`,
        // `"hole_r"`). Per the documented contract, a bare number is accepted as
        // the string form of that number (`4` → `"4"`, `2.5` → `"2.5"`), so an
        // agent may set a constant dimension without quoting it; a string is
        // taken verbatim. A `bool` is the wrong shape and fails loudly.
        let expr = |v: &AgentValue| -> Result<String, String> {
            match v {
                AgentValue::Str(s) => Ok(s.clone()),
                AgentValue::Int(i) => Ok(i.to_string()),
                AgentValue::Float(f) => Ok(f.to_string()),
                other => Err(format!(
                    "expected an expression string or number, got {other:?}"
                )),
            }
        };

        // -- Singletons --------------------------------------------------------
        match name {
            "Material density" => {
                self.density = value.as_f64()?;
                return Ok(());
            }
            "Extrude height" => {
                self.sketch_extrude_height = value.as_f64()?;
                return Ok(());
            }
            "circle radius =" => {
                // The parameter *name* that drives the sketch circle's radius.
                self.radius_param = value.as_str()?.to_string();
                return Ok(());
            }
            "Snap to grid" => {
                self.sketch_grid_snap = value.as_bool()?;
                return Ok(());
            }
            "Sketch tool" => {
                self.sketch_tool = parse_sketch_tool(value.as_str()?)?;
                return Ok(());
            }
            _ => {}
        }

        // -- Per feature-tree step (`step {N} <field>`, 1-based) ---------------
        if let Some(rest) = name.strip_prefix("step ") {
            // Split "<N> <field>" on the first space.
            let (num, field) = rest
                .split_once(' ')
                .ok_or_else(|| format!("malformed step control caption: {name:?}"))?;
            let idx1: usize = num
                .parse()
                .map_err(|_| format!("step control caption has a non-numeric index: {name:?}"))?;
            if idx1 == 0 {
                return Err(format!("step index is 1-based, got 0 in {name:?}"));
            }
            // Compute the length up front so the out-of-range message can borrow
            // it (the `get_mut` below holds a mutable borrow of `self.steps`).
            let n_steps = self.steps.len();
            let st = self.steps.get_mut(idx1 - 1).ok_or_else(|| {
                format!("step {idx1} is out of range (the feature tree has {n_steps} step(s))")
            })?;
            match field {
                // Enums (combos) — read by option name (string only).
                "op" => st.op = parse_op(value.as_str()?)?,
                "shape" => st.kind = parse_feature_kind(value.as_str()?)?,
                // Dimension / placement / rotation expressions (stored verbatim;
                // a bare number is taken as its string form, see `expr`).
                "dx" => st.dx = expr(value)?,
                "dy" => st.dy = expr(value)?,
                "dz" => st.dz = expr(value)?,
                "r" => st.radius = expr(value)?,
                "h" => st.height = expr(value)?,
                "base r" => st.radius = expr(value)?,
                "top r" => st.top_radius = expr(value)?,
                "major" => st.major = expr(value)?,
                "minor" => st.minor = expr(value)?,
                "height" => st.height = expr(value)?,
                "angle\u{00B0}" => st.angle = expr(value)?,
                "x" => st.x = expr(value)?,
                "y" => st.y = expr(value)?,
                "z" => st.z = expr(value)?,
                "rx\u{00B0}" => st.rx = expr(value)?,
                "ry\u{00B0}" => st.ry = expr(value)?,
                "rz\u{00B0}" => st.rz = expr(value)?,
                other => {
                    return Err(format!(
                        "unknown step field {other:?} in CAD control {name:?}"
                    ))
                }
            }
            return Ok(());
        }

        Err(format!("unknown CAD control: {name:?}"))
    }
}

struct CadResults {
    /// Each parameter's resolved value or error message.
    resolved: Vec<(String, Result<f64, String>)>,
    /// The solved circle radius, if the sketch solved.
    solved_radius: Option<f64>,
    /// Solver / status message.
    status: String,
}

/// Build a [`ParameterTable`] from the editable rows, skipping blank names.
fn build_table(params: &[(String, String)]) -> ParameterTable {
    let mut table = ParameterTable::new();
    for (n, e) in params {
        let n = n.trim();
        if !n.is_empty() {
            table.set(n, e);
        }
    }
    table
}

/// Resolve the parameters and solve a circle whose radius is driven by the
/// chosen parameter.
fn run_cad(s: &CadWorkbenchState) -> CadResults {
    let table = build_table(&s.params);
    let resolved: Vec<(String, Result<f64, String>)> = s
        .params
        .iter()
        .filter(|(n, _)| !n.trim().is_empty())
        .map(|(n, _)| (n.clone(), table.value(n.trim()).map_err(|e| e.to_string())))
        .collect();

    let (solved_radius, status) = match table.value(s.radius_param.trim()) {
        Ok(r) => {
            let mut sk = Sketch3D::new();
            let c = sk.add_point(0.0, 0.0, 0.0);
            let circle = sk
                .add_circle(c, 1.0, 0.0, 0.0, 1.0)
                .expect("centre is a point");
            sk.add_constraint(Constraint3D::CircleRadius { circle, target: r });
            match sk.solve() {
                Ok(rep) => {
                    let solved = sk.circle_radius(circle);
                    (
                        Some(solved),
                        format!("{:?} — circle radius = {solved:.4}", rep.status),
                    )
                }
                Err(e) => (None, format!("solve error: {e}")),
            }
        }
        Err(e) => (
            None,
            format!("radius parameter '{}': {e}", s.radius_param.trim()),
        ),
    };

    CadResults {
        resolved,
        solved_radius,
        status,
    }
}

/// Build a solver-crate [`FeatureTimeline`] from the UI steps.
fn steps_to_timeline(steps: &[UiStep]) -> FeatureTimeline {
    let mut tl = FeatureTimeline::new();
    for st in steps {
        tl.push(st.to_step());
    }
    tl
}

/// Reconstruct a UI step from a persisted solver-crate [`Step`]. Only box and
/// cylinder primitives are editable here; an extrude feature in a hand-authored
/// file is reported as an error rather than silently mangled.
fn ui_step_from(step: &Step) -> Result<UiStep, String> {
    let mut us = UiStep::new_box();
    us.op = step.op;
    us.x = step.at[0].clone();
    us.y = step.at[1].clone();
    us.z = step.at[2].clone();
    us.rx = step.rotate_deg[0].clone();
    us.ry = step.rotate_deg[1].clone();
    us.rz = step.rotate_deg[2].clone();
    match &step.feature {
        Feature::Box { dx, dy, dz } => {
            us.kind = FeatureKind::Box;
            us.dx = dx.clone();
            us.dy = dy.clone();
            us.dz = dz.clone();
        }
        Feature::Cylinder { radius, height } => {
            us.kind = FeatureKind::Cylinder;
            us.radius = radius.clone();
            us.height = height.clone();
        }
        Feature::Sphere { radius } => {
            us.kind = FeatureKind::Sphere;
            us.radius = radius.clone();
        }
        Feature::Cone {
            base_radius,
            top_radius,
            height,
        } => {
            us.kind = FeatureKind::Cone;
            us.radius = base_radius.clone();
            us.top_radius = top_radius.clone();
            us.height = height.clone();
        }
        Feature::Torus {
            major_radius,
            minor_radius,
        } => {
            us.kind = FeatureKind::Torus;
            us.major = major_radius.clone();
            us.minor = minor_radius.clone();
        }
        Feature::Extrude { profile, height } => {
            us.kind = FeatureKind::Extrude;
            us.profile = profile.clone();
            us.height = height.clone();
        }
        Feature::Revolve { profile, angle_deg } => {
            us.kind = FeatureKind::Revolve;
            us.profile = profile.clone();
            us.angle = angle_deg.clone();
        }
    }
    Ok(us)
}

/// Serialise the workbench's parameters + feature tree to a RON string.
fn save_string(s: &CadWorkbenchState) -> Result<String, String> {
    timeline_to_ron(&s.params, &steps_to_timeline(&s.steps)).map_err(|e| e.to_string())
}

/// Editable parameters (`(name, expression)` rows) plus UI steps, parsed from
/// a saved document.
type LoadedTree = (Vec<(String, String)>, Vec<UiStep>);

/// Parse a RON document into editable parameters + UI steps.
fn load_from_string(txt: &str) -> Result<LoadedTree, String> {
    let doc = timeline_from_ron(txt).map_err(|e| e.to_string())?;
    let steps = doc
        .timeline
        .steps
        .iter()
        .map(ui_step_from)
        .collect::<Result<Vec<_>, _>>()?;
    Ok((doc.parameters, steps))
}

/// Total solid volume of a set of bodies, in cubic model units. Sums each
/// body's best-effort volume (`valenx_cad::solid_volume`); a body whose volume
/// can't be evaluated is skipped rather than poisoning the total. Exact for
/// flat-faced solids; curved solids converge from below at the measure
/// tolerance, so this slightly under-reports them.
fn total_volume(bodies: &[valenx_cad::Solid]) -> f64 {
    bodies
        .iter()
        .filter_map(|b| valenx_cad::solid_volume(b).ok())
        .sum()
}

/// Total boundary surface area of a set of bodies, in square model units. Sums
/// each body's best-effort area (`valenx_cad::solid_area`), skipping any that
/// can't be evaluated. Exact for flat-faced solids; curved solids converge from
/// below at the measure tolerance, so this slightly under-reports them.
fn total_area(bodies: &[valenx_cad::Solid]) -> f64 {
    bodies
        .iter()
        .filter_map(|b| valenx_cad::solid_area(b).ok())
        .sum()
}

/// The **B-rep Euler–Poincaré characteristic** `χ = V − E + F` summed over a set of
/// solid bodies, threading [`valenx_cad::euler_characteristic`] (each body's exact
/// boundary-representation `V − E + F`). Unlike the *tessellation* Euler
/// characteristic [`mesh_euler_characteristic`] — which counts the triangle mesh —
/// this reads the solid's true CAD topology directly: a single closed orientable body
/// has `χ = 2 − 2g` (`2` for a box or sphere, genus 0; `0` for a torus), and the
/// value is additive over disjoint bodies. Because the Euler characteristic is a
/// topological invariant the two must agree (a B-rep-vs-mesh consistency check). A
/// mesh-backed solid contributes `0` (no B-rep topology — its `euler_characteristic`
/// is `None`); `0` for an empty body list.
fn brep_euler_characteristic(bodies: &[valenx_cad::Solid]) -> i64 {
    bodies
        .iter()
        .filter_map(valenx_cad::euler_characteristic)
        .sum()
}

/// Overall bounding-box dimensions `[dx, dy, dz]` (model units) of a mesh, from
/// its axis-aligned min/max corner. `None` for an empty mesh.
fn mesh_dimensions(mesh: &valenx_mesh::Mesh) -> Option<[f32; 3]> {
    let (min, max) = crate::mesh_loader::mesh_bounding_box(mesh)?;
    Some([max[0] - min[0], max[1] - min[1], max[2] - min[2]])
}

/// Fill fraction (rectangularity): solid `volume` divided by the bounding-box
/// volume implied by `dims`. `1.0` for a box that fills its bounding box, `~π/4`
/// for a cylinder, `~π/6` for a sphere. `None` when the bounding box has
/// effectively zero volume.
fn fill_fraction(volume: f64, dims: [f32; 3]) -> Option<f64> {
    let bbox_volume = dims[0] as f64 * dims[1] as f64 * dims[2] as f64;
    if bbox_volume > 1e-12 {
        Some(volume / bbox_volume)
    } else {
        None
    }
}

/// Wadell sphericity `ψ = π^(1/3)·(6V)^(2/3) / A` — the dimensionless
/// compactness of a solid of `volume` and surface `area`, normalised so a
/// perfect sphere is `1.0` (a cube is ≈ 0.806; more elongated or branched
/// shapes score lower). Independent of orientation and scale. `None` when the
/// area or volume is non-positive. Because curved-solid `volume`/`area`
/// under-report slightly at the measure tolerance, a tessellated sphere reads a
/// hair under 1.
fn sphericity(volume: f64, area: f64) -> Option<f64> {
    if area > 1e-12 && volume > 0.0 {
        Some(std::f64::consts::PI.cbrt() * (6.0 * volume).powf(2.0 / 3.0) / area)
    } else {
        None
    }
}

/// The space diagonal `√(dx² + dy² + dz²)` of a bounding box of dimensions
/// `dims` (model units) — the part's overall reach, equal to the diameter of
/// the smallest sphere that encloses the box. Exact from the axis-aligned
/// extents, so it is independent of tessellation density.
fn bbox_diagonal(dims: [f32; 3]) -> f64 {
    let (dx, dy, dz) = (dims[0] as f64, dims[1] as f64, dims[2] as f64);
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// The exact **mesh diameter** — the largest straight-line distance between any
/// two vertices (the longest chord, the caliper extent of the part). Unlike the
/// axis-aligned `bbox_diagonal`, which measures the box *around* the part and so
/// overestimates the true span of any shape not aligned to the axes, this is the
/// genuine maximum reach of the geometry; and unlike the *approximate* Ritter
/// enclosing-sphere radius from `mesh_bounding_sphere_radius`, it is the exact
/// extremal pair — the diameter satisfies `d ≤ 2·R_enc`, with equality only when
/// the farthest pair is antipodal through the sphere centre. Computed by an
/// exact all-pairs scan, `O(n²)` in the vertex count, which is comfortably fast
/// for the workshop-scale tessellations the workbench produces. `None` for a
/// mesh with no vertices; a lone vertex gives `0`. Model units.
fn mesh_diameter(mesh: &valenx_mesh::Mesh) -> Option<f64> {
    let pts = &mesh.nodes;
    if pts.is_empty() {
        return None;
    }
    let mut max_sq = 0.0_f64;
    for (i, a) in pts.iter().enumerate() {
        for b in &pts[i + 1..] {
            max_sq = max_sq.max((a - b).norm_squared());
        }
    }
    Some(max_sq.sqrt())
}

/// The number of **sharp (feature) edges** of a `Tri3` mesh — interior edges
/// (shared by exactly two triangles) where the surface *folds* by more than
/// `angle_threshold_rad`, the angle between the two triangles' outward face
/// normals. These are the crease / feature edges a CAD kernel walks for
/// filleting, chamfering and crisp shaded display: unlike the volume, inertia,
/// topology and quality measures — which describe how *much* material there is
/// or how the mesh connects — this finds *where the surface bends*. A flat or
/// smoothly curved patch contributes none; a hard edge (a box rim, a Boolean
/// seam) contributes one each. The in-facet diagonals a tessellator adds to
/// split a flat polygon are coplanar (zero fold) and so are never counted.
/// Degenerate (zero-area) triangles and boundary edges (only one incident
/// triangle) are skipped.
fn mesh_sharp_edge_count(mesh: &valenx_mesh::Mesh, angle_threshold_rad: f64) -> usize {
    use std::collections::HashMap;
    let mut edge_normals: HashMap<(u32, u32), Vec<nalgebra::Vector3<f64>>> = HashMap::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            let normal = (b - a).cross(&(c - a));
            let len = normal.norm();
            if len <= 0.0 {
                continue; // degenerate triangle has no defined normal
            }
            let unit = normal / len;
            for &(i, j) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if i <= j { (i, j) } else { (j, i) };
                edge_normals.entry(key).or_default().push(unit);
            }
        }
    }
    edge_normals
        .values()
        .filter(|normals| normals.len() == 2)
        .filter(|normals| normals[0].dot(&normals[1]).clamp(-1.0, 1.0).acos() > angle_threshold_rad)
        .count()
}

/// The largest **dihedral (crease) angle** anywhere on the surface (rad) — the
/// sharpest fold, the maximum over every interior manifold edge (shared by
/// exactly two triangles) of the angle between the two triangles' outward face
/// normals (`0` where the surface is flat, `π` where it folds back on itself).
/// Where [`mesh_sharp_edge_count`] *counts* how many edges fold past a fixed
/// threshold, this reports *how sharp the sharpest fold is* — the single number
/// that sizes the smallest fillet/chamfer the part needs and flags its crispest
/// feature. A box rim folds at exactly `π/2`; a smoothly tessellated curved
/// patch folds only gently. The in-facet diagonals a tessellator adds to split a
/// flat polygon are coplanar (zero fold), and degenerate (zero-area) triangles
/// and boundary edges (one incident triangle) are ignored. `None` for a mesh
/// with no interior manifold edge.
fn mesh_max_dihedral_angle(mesh: &valenx_mesh::Mesh) -> Option<f64> {
    use std::collections::HashMap;
    let mut edge_normals: HashMap<(u32, u32), Vec<nalgebra::Vector3<f64>>> = HashMap::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            let normal = (b - a).cross(&(c - a));
            let len = normal.norm();
            if len <= 0.0 {
                continue; // degenerate triangle has no defined normal
            }
            let unit = normal / len;
            for &(i, j) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if i <= j { (i, j) } else { (j, i) };
                edge_normals.entry(key).or_default().push(unit);
            }
        }
    }
    edge_normals
        .values()
        .filter(|normals| normals.len() == 2)
        .map(|normals| normals[0].dot(&normals[1]).clamp(-1.0, 1.0).acos())
        .reduce(f64::max)
}

/// The volume enclosed by a closed `Tri3` surface mesh `V = ⅙·|Σ aᵢ·(bᵢ×cᵢ)|`
/// (cubic model units), via the divergence theorem: each triangle spans a signed
/// tetrahedron with the origin contributing `a·(b×c)/6`, the parts outside the
/// surface cancelling. This is the **zeroth mass-moment** that the volume centroid
/// [`mesh_centroid`], radius of gyration [`mesh_radius_of_gyration`] and principal
/// moments [`mesh_principal_moments`] are all weighted by — surfaced here in its own
/// right to lead that suite. It is integrated over the **tessellated boundary**, so
/// comparing it to the CSG [`total_volume`] of the source solid gauges tessellation
/// fidelity: the two agree for an exact watertight mesh and diverge as faceting
/// coarsens. The magnitude is returned, so the result is independent of the global
/// winding direction. `None` for a mesh with no triangles or effectively zero
/// volume (open / degenerate).
fn mesh_volume(mesh: &valenx_mesh::Mesh) -> Option<f64> {
    let mut total_vol = 0.0;
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            total_vol += a.dot(&b.cross(&c)) / 6.0;
        }
    }
    if total_vol.abs() > 1e-12 {
        Some(total_vol.abs())
    } else {
        None
    }
}

/// The volume centroid (centre of mass at uniform density) `[x, y, z]` of a
/// closed `Tri3` surface mesh, via the divergence theorem: each triangle forms
/// a signed tetrahedron with the origin, contributing its signed volume
/// `a·(b×c)/6` weighted by the tet centroid `(a+b+c)/4`. This is the true
/// volume centroid — distinct from the vertex average for any non-symmetric
/// tessellation. `None` for a mesh with no triangles or effectively zero signed
/// volume (open / degenerate). Model units.
fn mesh_centroid(mesh: &valenx_mesh::Mesh) -> Option<[f64; 3]> {
    let mut total_vol = 0.0;
    let mut acc = [0.0_f64; 3];
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            let v = a.dot(&b.cross(&c)) / 6.0;
            total_vol += v;
            acc[0] += v * (a.x + b.x + c.x) / 4.0;
            acc[1] += v * (a.y + b.y + c.y) / 4.0;
            acc[2] += v * (a.z + b.z + c.z) / 4.0;
        }
    }
    if total_vol.abs() > 1e-12 {
        Some([acc[0] / total_vol, acc[1] / total_vol, acc[2] / total_vol])
    } else {
        None
    }
}

/// The radius of gyration about the centroid `k = √⟨r²⟩` (model units) of a
/// closed `Tri3` mesh — the root-mean-square distance of the solid's volume from
/// its centre of mass, the length that sets its rotational inertia (`I ≈ m·k²`).
/// Extends the divergence-theorem volume integral with the diagonal second
/// moments `∫x²,∫y²,∫z² dV` (tetrahedron covariance `(V/20)(Σpᵢ² + (Σpᵢ)²)` per
/// axis), then removes the centroid offset via the parallel-axis theorem. For an
/// `Lx×Ly×Lz` box it is `√((Lx²+Ly²+Lz²)/12)`. `None` for a mesh with no
/// triangles or effectively zero signed volume.
fn mesh_radius_of_gyration(mesh: &valenx_mesh::Mesh) -> Option<f64> {
    let mut total_vol = 0.0;
    let mut first = [0.0_f64; 3]; // ∫x, ∫y, ∫z dV (→ centroid)
    let mut second = [0.0_f64; 3]; // ∫x², ∫y², ∫z² dV about the origin
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            let v = a.dot(&b.cross(&c)) / 6.0;
            total_vol += v;
            first[0] += v * (a.x + b.x + c.x) / 4.0;
            first[1] += v * (a.y + b.y + c.y) / 4.0;
            first[2] += v * (a.z + b.z + c.z) / 4.0;
            second[0] += v / 20.0 * (a.x * a.x + b.x * b.x + c.x * c.x + (a.x + b.x + c.x).powi(2));
            second[1] += v / 20.0 * (a.y * a.y + b.y * b.y + c.y * c.y + (a.y + b.y + c.y).powi(2));
            second[2] += v / 20.0 * (a.z * a.z + b.z * b.z + c.z * c.z + (a.z + b.z + c.z).powi(2));
        }
    }
    if total_vol.abs() < 1e-12 {
        return None;
    }
    // ⟨r²⟩ about the origin, minus the centroid offset (parallel-axis theorem).
    let mean_sq_origin = (second[0] + second[1] + second[2]) / total_vol;
    let centroid_sq =
        (first[0] * first[0] + first[1] * first[1] + first[2] * first[2]) / (total_vol * total_vol);
    Some((mean_sq_origin - centroid_sq).max(0.0).sqrt())
}

/// The three **principal moments of inertia** about the centroid (unit density)
/// of a closed `Tri3` mesh — the eigenvalues of the centroidal inertia tensor,
/// returned sorted descending (`I₁ ≥ I₂ ≥ I₃`, model units `u⁵`). This completes
/// the mass-property suite after the volume, centroid and radius of gyration:
/// it extends the divergence-theorem second moments with the **products of
/// inertia** `∫xy, ∫xz, ∫yz dV` (same tetrahedron covariance
/// `(V/20)(Σpᵢpⱼ + (Σpᵢ)(Σpⱼ))` as the diagonal terms), shifts every second
/// moment to the centroid by the parallel-axis theorem, assembles the symmetric
/// inertia tensor `I = [[∫y²+z², −∫xy, −∫xz], [−∫xy, ∫x²+z², −∫yz], [−∫xz, −∫yz,
/// ∫x²+y²]]`, and takes its symmetric eigendecomposition (the eigenvectors are
/// the principal axes). For an `Lx×Ly×Lz` box they are `(m/12)(Lⱼ²+Lₖ²)` with
/// `m = V`. `None` for a mesh with no triangles or effectively zero volume.
fn mesh_principal_moments(mesh: &valenx_mesh::Mesh) -> Option<[f64; 3]> {
    let mut vol = 0.0;
    let mut first = [0.0_f64; 3]; // ∫x, ∫y, ∫z dV (→ centroid)
    let mut m = [0.0_f64; 6]; // second moments about the origin: xx, yy, zz, xy, xz, yz
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            let v = a.dot(&b.cross(&c)) / 6.0;
            vol += v;
            let s = a + b + c; // Σ pᵢ
            first[0] += v * s.x / 4.0;
            first[1] += v * s.y / 4.0;
            first[2] += v * s.z / 4.0;
            // ∫xᵢxⱼ contribution = (v/20)·(aᵢaⱼ + bᵢbⱼ + cᵢcⱼ + sᵢsⱼ).
            m[0] += v / 20.0 * (a.x * a.x + b.x * b.x + c.x * c.x + s.x * s.x);
            m[1] += v / 20.0 * (a.y * a.y + b.y * b.y + c.y * c.y + s.y * s.y);
            m[2] += v / 20.0 * (a.z * a.z + b.z * b.z + c.z * c.z + s.z * s.z);
            m[3] += v / 20.0 * (a.x * a.y + b.x * b.y + c.x * c.y + s.x * s.y);
            m[4] += v / 20.0 * (a.x * a.z + b.x * b.z + c.x * c.z + s.x * s.z);
            m[5] += v / 20.0 * (a.y * a.z + b.y * b.z + c.y * c.z + s.y * s.z);
        }
    }
    if vol.abs() < 1e-12 {
        return None;
    }
    let cx = first[0] / vol;
    let cy = first[1] / vol;
    let cz = first[2] / vol;
    // Parallel-axis shift of every second moment to the centroid: M_ij,c = M_ij,o − V·cᵢcⱼ.
    let mxx = m[0] - vol * cx * cx;
    let myy = m[1] - vol * cy * cy;
    let mzz = m[2] - vol * cz * cz;
    let mxy = m[3] - vol * cx * cy;
    let mxz = m[4] - vol * cx * cz;
    let myz = m[5] - vol * cy * cz;
    let tensor = nalgebra::Matrix3::new(
        myy + mzz,
        -mxy,
        -mxz,
        -mxy,
        mxx + mzz,
        -myz,
        -mxz,
        -myz,
        mxx + myy,
    );
    let eig = tensor.symmetric_eigen().eigenvalues;
    let mut vals = [eig[0], eig[1], eig[2]];
    // Sort descending: I₁ ≥ I₂ ≥ I₃.
    vals.sort_by(|p, q| q.partial_cmp(p).unwrap_or(std::cmp::Ordering::Equal));
    Some(vals)
}

/// Whether a `Tri3` surface mesh is **watertight** (a closed 2-manifold) — every
/// triangle edge is shared by *exactly two* triangles, so the surface bounds a
/// solid with no holes or cracks. This is the printability / validity gate a
/// slicer applies before it can decide what is "inside": an open edge (count 1)
/// is a boundary hole and a non-manifold edge (count ≥ 3) is a self-intersection
/// or T-junction. Edges are counted undirected via a sorted `(min, max)`
/// node-index key, so triangle winding is ignored. Returns `false` for a mesh
/// with no triangles.
fn mesh_is_watertight(mesh: &valenx_mesh::Mesh) -> bool {
    let mut edges: std::collections::HashMap<(u32, u32), u32> = std::collections::HashMap::new();
    let mut any = false;
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            any = true;
            for &(i, j) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if i <= j { (i, j) } else { (j, i) };
                *edges.entry(key).or_insert(0) += 1;
            }
        }
    }
    any && edges.values().all(|&c| c == 2)
}

/// The number of **open (boundary) edges** of a `Tri3` surface mesh — edges
/// belonging to exactly *one* triangle, i.e. the combined perimeter of every
/// hole and crack in the surface. Where [`mesh_is_watertight`] answers the
/// yes/no "is the surface closed?", this answers *how* open it is: a watertight
/// solid reports `0`, and a larger count means bigger or more numerous holes —
/// the number a mesh-repair or slicer step reports to size the defect before
/// deciding whether the mesh can be printed or Boolean-ed. Edges are counted
/// undirected via the same sorted `(min, max)` node-index key as the watertight
/// check, so triangle winding is ignored. Returns `0` for a mesh with no
/// triangles.
fn mesh_open_edge_count(mesh: &valenx_mesh::Mesh) -> usize {
    let mut edges: std::collections::HashMap<(u32, u32), u32> = std::collections::HashMap::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            for &(i, j) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if i <= j { (i, j) } else { (j, i) };
                *edges.entry(key).or_insert(0) += 1;
            }
        }
    }
    edges.values().filter(|&&c| c == 1).count()
}

/// The number of **non-manifold edges** of a `Tri3` surface mesh — edges shared by
/// *three or more* triangles, the self-intersection / T-junction defect a
/// well-formed 2-manifold never has (every interior edge of a clean surface joins
/// exactly two faces). It is the quantitative companion to
/// [`mesh_open_edge_count`] on the *other* side of the manifold ideal: where the
/// open-edge count sizes the holes (edges with too *few* — one — incident
/// triangles), this sizes the self-intersections (edges with too *many*), the two
/// failure modes the boolean [`mesh_is_watertight`] collapses into a single
/// yes/no. Non-manifold edges are the harder mesh-repair case — a slicer cannot
/// decide inside-from-outside across them. Edges are counted undirected via the
/// same sorted `(min, max)` node-index key. Returns `0` for a clean or empty mesh.
fn mesh_nonmanifold_edge_count(mesh: &valenx_mesh::Mesh) -> usize {
    let mut edges: std::collections::HashMap<(u32, u32), u32> = std::collections::HashMap::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            for &(i, j) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if i <= j { (i, j) } else { (j, i) };
                *edges.entry(key).or_insert(0) += 1;
            }
        }
    }
    edges.values().filter(|&&c| c >= 3).count()
}

/// The **mean edge length** — the characteristic element size `h̄` (model units)
/// of a `Tri3` mesh, the average length of its *distinct* undirected edges. This
/// is the single "mesh resolution" scalar: the typical element size that sets the
/// interpolation and quadrature accuracy of any analysis run on the surface
/// (discretisation error usually scaling as a power of `h`). Where the shape
/// metrics ([`mesh_min_triangle_quality`], [`mesh_max_aspect_ratio`]) grade *how
/// well-formed* each element is, and the [`mesh_diameter`] / bounding box give
/// the part's *overall* extent, this reports the *local* discretisation scale.
/// Each edge is counted once however many triangles share it (the in-facet
/// diagonals a tessellator adds are included). `None` for a mesh with no edges.
fn mesh_mean_edge_length(mesh: &valenx_mesh::Mesh) -> Option<f64> {
    let mut edges: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            for &(i, j) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                edges.insert(if i <= j { (i, j) } else { (j, i) });
            }
        }
    }
    if edges.is_empty() {
        return None;
    }
    let total: f64 = edges
        .iter()
        .map(|&(i, j)| (mesh.nodes[i as usize] - mesh.nodes[j as usize]).norm())
        .sum();
    Some(total / edges.len() as f64)
}

/// The **total edge length** — the cumulative length of all *distinct* undirected
/// edges of a `Tri3` mesh (model units), i.e. the length of its 1-skeleton
/// (wireframe). It is the 1-D member of the tessellation's dimensional-measure
/// family — alongside the 2-D [`mesh_surface_area`] and the 3-D [`mesh_volume`] —
/// and the aggregate companion to the *average* [`mesh_mean_edge_length`]
/// (`total = mean · edge_count`). Each edge is counted once however many triangles
/// share it (the in-facet diagonals a tessellator adds are included). `None` for a
/// mesh with no edges.
fn mesh_total_edge_length(mesh: &valenx_mesh::Mesh) -> Option<f64> {
    let mut edges: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            for &(i, j) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                edges.insert(if i <= j { (i, j) } else { (j, i) });
            }
        }
    }
    if edges.is_empty() {
        return None;
    }
    let total: f64 = edges
        .iter()
        .map(|&(i, j)| (mesh.nodes[i as usize] - mesh.nodes[j as usize]).norm())
        .sum();
    Some(total)
}

/// The **maximum vertex valence** (graph-degree) of a `Tri3` surface mesh — the
/// largest number of DISTINCT neighbouring vertices joined to any single vertex by
/// a triangle edge. It is a topological mesh-quality indicator: a regular
/// triangulation has interior valence 6, and an unusually high maximum flags an
/// irregular "pole"/fan vertex (many triangles crowding one point) — the kind of
/// singular vertex that degrades subdivision and FEM conditioning. It is the
/// connectivity companion to the geometric worst-case indicators
/// [`mesh_max_aspect_ratio`] and [`mesh_max_dihedral_angle`], and is distinct from
/// [`mesh_sharp_edge_count`] / [`mesh_open_edge_count`] (which count EDGES, not
/// vertex degree). Returns `None` for an empty mesh with no `Tri3` faces.
fn mesh_max_vertex_valence(mesh: &valenx_mesh::Mesh) -> Option<usize> {
    use std::collections::{HashMap, HashSet};
    let mut neighbours: HashMap<u32, HashSet<u32>> = HashMap::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            for &(i, j) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                neighbours.entry(i).or_default().insert(j);
                neighbours.entry(j).or_default().insert(i);
            }
        }
    }
    neighbours.values().map(HashSet::len).max()
}

/// The **mean vertex valence** (average graph-degree) of a `Tri3` surface mesh —
/// the number of incident edges per vertex averaged over all vertices, `2·E/V`. It
/// is the regularity hallmark of a triangle mesh: a well-formed closed surface
/// tends toward a mean valence of **6** (Euler's relation forces `2E/V = 6 − 12/V`
/// for a closed triangulation), and a mean far from 6 signals an irregular or open
/// mesh. It is the *average*-case companion to the worst-case
/// [`mesh_max_vertex_valence`]: the max flags the single most irregular "pole"
/// vertex, the mean reports the typical connectivity. Returns `None` for an empty
/// mesh with no `Tri3` faces.
fn mesh_mean_vertex_valence(mesh: &valenx_mesh::Mesh) -> Option<f64> {
    use std::collections::{HashMap, HashSet};
    let mut neighbours: HashMap<u32, HashSet<u32>> = HashMap::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            for &(i, j) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                neighbours.entry(i).or_default().insert(j);
                neighbours.entry(j).or_default().insert(i);
            }
        }
    }
    if neighbours.is_empty() {
        return None;
    }
    let total: usize = neighbours.values().map(HashSet::len).sum();
    Some(total as f64 / neighbours.len() as f64)
}

/// The number of **boundary loops** (holes) in a `Tri3` surface mesh — the
/// connected components of the boundary-edge graph (the multiplicity-1 edges that
/// bound an open patch). It answers the headline mesh-repair question — *how many
/// distinct holes are there to fill* — and is distinct from
/// [`mesh_open_edge_count`], which counts the boundary EDGES (a triangular hole is
/// 3 open edges but **one** loop), from [`mesh_shell_count`] (disjoint solids),
/// and from the global [`mesh_euler_characteristic`]. A closed (watertight) mesh
/// has `0`; each missing patch contributes one loop. Returns `0` for a closed or
/// empty mesh.
fn mesh_boundary_loop_count(mesh: &valenx_mesh::Mesh) -> usize {
    use std::collections::{HashMap, HashSet};
    // Edge multiplicities — a boundary edge belongs to exactly one triangle.
    let mut mult: HashMap<(u32, u32), u32> = HashMap::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            for &(i, j) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if i <= j { (i, j) } else { (j, i) };
                *mult.entry(key).or_insert(0) += 1;
            }
        }
    }
    // Adjacency among the endpoints of the boundary (multiplicity-1) edges.
    let mut adj: HashMap<u32, Vec<u32>> = HashMap::new();
    for (&(i, j), &c) in &mult {
        if c == 1 {
            adj.entry(i).or_default().push(j);
            adj.entry(j).or_default().push(i);
        }
    }
    // Count the connected components — each closed boundary loop is one hole.
    let mut visited: HashSet<u32> = HashSet::new();
    let mut loops = 0;
    for &start in adj.keys() {
        if !visited.insert(start) {
            continue; // already part of a counted loop
        }
        loops += 1;
        let mut stack = vec![start];
        while let Some(v) = stack.pop() {
            if let Some(neighbours) = adj.get(&v) {
                for &n in neighbours {
                    if visited.insert(n) {
                        stack.push(n);
                    }
                }
            }
        }
    }
    loops
}

/// The **Euler characteristic** `χ = V − E + F` of a `Tri3` surface mesh — `V`
/// the distinct *referenced* vertices, `E` the distinct undirected edges, `F`
/// the triangles. This is the fundamental topological invariant of the surface:
/// a closed genus-`g` orientable surface has `χ = 2 − 2g`, so a sphere or box
/// gives `χ = 2` (genus 0), a torus `χ = 0` (genus 1), and every extra handle
/// drops `χ` by 2. Unlike the geometric and mass measures it is invariant to any
/// bending or scaling — it sees only the connectivity — and so reports *which*
/// topology the watertight check has merely confirmed is closed. `None` for a
/// mesh with no triangles.
fn mesh_euler_characteristic(mesh: &valenx_mesh::Mesh) -> Option<i64> {
    let mut verts: std::collections::HashSet<u32> = std::collections::HashSet::new();
    let mut edges: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
    let mut faces = 0_i64;
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            faces += 1;
            for &v in tri {
                verts.insert(v);
            }
            for &(i, j) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                edges.insert(if i <= j { (i, j) } else { (j, i) });
            }
        }
    }
    if faces == 0 {
        return None;
    }
    Some(verts.len() as i64 - edges.len() as i64 + faces)
}

/// The **total angle defect** `Σ_v (2π − Σθ)` over a `Tri3` surface mesh (rad) —
/// the discrete **total Gaussian curvature**. At each vertex the *angle defect* is
/// `2π` minus the sum of the incident triangles' interior angles there (positive at
/// a convex corner, negative at a saddle, zero where the surface is locally flat),
/// and this is their sum over every vertex. By the **Gauss–Bonnet theorem** it
/// equals `2π·χ` for a closed surface — `4π` for a sphere or box (genus 0), `0` for
/// a torus — so it is the *geometric* (angle-based) reading of the very topology
/// that [`mesh_euler_characteristic`] reads combinatorially from `V − E + F`. Being
/// built from angles it is invariant to translation and uniform scaling. `None` for
/// a mesh with no triangles.
fn mesh_total_angle_defect(mesh: &valenx_mesh::Mesh) -> Option<f64> {
    let mut angle_sum: std::collections::HashMap<u32, f64> = std::collections::HashMap::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let p = [
                mesh.nodes[tri[0] as usize],
                mesh.nodes[tri[1] as usize],
                mesh.nodes[tri[2] as usize],
            ];
            for k in 0..3 {
                let u = p[(k + 1) % 3] - p[k];
                let v = p[(k + 2) % 3] - p[k];
                *angle_sum.entry(tri[k]).or_insert(0.0) += u.cross(&v).norm().atan2(u.dot(&v));
            }
        }
    }
    if angle_sum.is_empty() {
        return None;
    }
    let two_pi = 2.0 * std::f64::consts::PI;
    let total: f64 = angle_sum.values().map(|&s| two_pi - s).sum();
    Some(total)
}

/// The number of **connected shells** (disjoint solids) in a `Tri3` mesh — how
/// many separate connected surfaces it holds, by union-find over the triangle
/// vertices (every triangle ties its three nodes into one component). This is
/// the count the watertight check (is the surface closed) and the Euler
/// characteristic (which topology) leave implicit: a Cut that severed the body,
/// or a Boolean that left two pieces, shows up here as more than one shell
/// (`χ` is the sum over them). Returns `0` for a mesh with no triangles.
fn mesh_shell_count(mesh: &valenx_mesh::Mesh) -> usize {
    use std::collections::{HashMap, HashSet};
    // Union-find (path-halving find) over the referenced vertices.
    fn find(parent: &mut HashMap<u32, u32>, mut x: u32) -> u32 {
        loop {
            let p = parent.get(&x).copied().unwrap_or(x);
            if p == x {
                return x;
            }
            let gp = parent.get(&p).copied().unwrap_or(p);
            parent.insert(x, gp); // path halving
            x = gp;
        }
    }
    let mut parent: HashMap<u32, u32> = HashMap::new();
    let mut any = false;
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            any = true;
            for &v in tri {
                parent.entry(v).or_insert(v);
            }
            let r0 = find(&mut parent, tri[0]);
            let r1 = find(&mut parent, tri[1]);
            let r2 = find(&mut parent, tri[2]);
            if r1 != r0 {
                parent.insert(r1, r0);
            }
            if r2 != r0 {
                parent.insert(r2, r0);
            }
        }
    }
    if !any {
        return 0;
    }
    let verts: Vec<u32> = parent.keys().copied().collect();
    let mut roots: HashSet<u32> = HashSet::new();
    for v in verts {
        let r = find(&mut parent, v);
        roots.insert(r);
    }
    roots.len()
}

/// The radius of the **minimal enclosing sphere** (model units) of the mesh
/// vertices, by **Ritter's two-pass algorithm**. Unlike the axis-aligned
/// bounding box, the enclosing sphere is *orientation-independent* — the single
/// tightest sphere that contains the whole part — which is why it is the
/// standard bound for broad-phase collision tests and view-frustum culling.
/// Ritter seeds the sphere on the most-separated vertex pair, then grows it
/// minimally to swallow any outlier; the result is within a few percent of the
/// true minimum (and never smaller). For the 8 corners of an `Lx×Ly×Lz` box it
/// is exactly the half space-diagonal `½√(Lx²+Ly²+Lz²)`. `None` for a mesh with
/// no vertices.
fn mesh_bounding_sphere_radius(mesh: &valenx_mesh::Mesh) -> Option<f64> {
    let pts = &mesh.nodes;
    if pts.is_empty() {
        return None;
    }
    // The vertex farthest from `q` (squared distance, to avoid the sqrt).
    let farthest = |q: nalgebra::Vector3<f64>| -> nalgebra::Vector3<f64> {
        pts.iter()
            .copied()
            .max_by(|a, b| {
                (a - q)
                    .norm_squared()
                    .partial_cmp(&(b - q).norm_squared())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(q)
    };
    // Pass 1: seed on the most-separated pair found by the two farthest-point hops.
    let y = farthest(pts[0]);
    let z = farthest(y);
    let mut centre = (y + z) / 2.0;
    let mut radius = (y - z).norm() / 2.0;
    // Pass 2: grow the sphere just enough to cover any vertex left outside.
    for &p in pts {
        let d = (p - centre).norm();
        if d > radius {
            // Shift the centre toward `p` by (d−r)/2 (old radius), then update r.
            centre += (p - centre) * ((d - radius) / (2.0 * d));
            radius = (radius + d) / 2.0;
        }
    }
    Some(radius)
}

/// The **minimum normalized triangle shape quality** over every `Tri3` face —
/// `Q = 4√3·A / (ℓ₁² + ℓ₂² + ℓ₃²)`, the triangle area over its sum of squared
/// edge lengths, scaled so an *equilateral* triangle scores exactly `Q = 1` and
/// a degenerate sliver scores `Q → 0`. Where every other measure here describes
/// the *shape the mesh represents* (its volume, inertia, topology), this grades
/// the *mesh itself*: the worst-shaped element caps interpolation accuracy and
/// solver conditioning, so the minimum `Q` over all faces is the one number a
/// mesher checks before handing a surface to analysis (a long, thin "sliver"
/// triangle is the classic culprit). Returns `None` for a mesh with no
/// non-degenerate triangles.
fn mesh_min_triangle_quality(mesh: &valenx_mesh::Mesh) -> Option<f64> {
    let scale = 4.0 * 3.0_f64.sqrt();
    let mut worst: Option<f64> = None;
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            let sum_sq_edges =
                (b - a).norm_squared() + (c - b).norm_squared() + (a - c).norm_squared();
            if sum_sq_edges <= 0.0 {
                continue; // a zero-perimeter (coincident-vertex) triangle is not an element
            }
            let area = 0.5 * (b - a).cross(&(c - a)).norm();
            let q = scale * area / sum_sq_edges;
            worst = Some(worst.map_or(q, |w| w.min(q)));
        }
    }
    worst
}

/// The largest **triangle aspect ratio** over every `Tri3` face — the worst
/// longest-edge-to-shortest-edge ratio `ℓ_max/ℓ_min`, which is exactly `1` for an
/// equilateral triangle and grows without bound as the triangle is stretched into
/// a sliver. It is the *elongation* companion to [`mesh_min_triangle_quality`]:
/// that radius-ratio score grades overall triangle shape, while this isolates pure
/// stretch, and the two together are the standard mesh-quality pair a mesher
/// reports — a high aspect ratio warns of the anisotropy that skews interpolation
/// gradients and solver conditioning even when the area-based quality still looks
/// acceptable. Degenerate (zero-shortest-edge) triangles are skipped; returns
/// `None` for a mesh with no usable triangles. Dimensionless.
fn mesh_max_aspect_ratio(mesh: &valenx_mesh::Mesh) -> Option<f64> {
    let mut worst: Option<f64> = None;
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            let edges = [(b - a).norm(), (c - b).norm(), (a - c).norm()];
            let shortest = edges[0].min(edges[1]).min(edges[2]);
            if shortest <= 0.0 {
                continue; // a coincident-vertex edge is not a real triangle
            }
            let longest = edges[0].max(edges[1]).max(edges[2]);
            worst = Some(worst.map_or(longest / shortest, |w| w.max(longest / shortest)));
        }
    }
    worst
}

/// The **smallest triangle (facet) area** in a `Tri3` surface mesh (model units²)
/// — a degeneracy / resolution indicator: a near-zero-area sliver triangle wrecks
/// finite-element conditioning and shading, and the minimum facet area sets the
/// finest spatial scale the tessellation resolves. It is the *area* companion to
/// the *shape* score [`mesh_min_triangle_quality`] (a unitless 0–1 measure) and
/// the *length* scale [`mesh_mean_edge_length`], measuring magnitude rather than
/// shape or edge length. Each facet's area is `½‖(b−a)×(c−a)‖`. Returns `None` for
/// an empty mesh with no `Tri3` faces.
fn mesh_min_triangle_area(mesh: &valenx_mesh::Mesh) -> Option<f64> {
    let mut smallest: Option<f64> = None;
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            let area = 0.5 * (b - a).cross(&(c - a)).norm();
            smallest = Some(smallest.map_or(area, |s| s.min(area)));
        }
    }
    smallest
}

/// The **largest triangle (facet) area** in a `Tri3` surface mesh (model units²) —
/// the coarsest facet, the area-extremes companion to [`mesh_min_triangle_area`].
/// Together the two bound the mesh's *area grading*: a large max/min spread marks a
/// non-uniformly refined tessellation (fine where it matters, coarse elsewhere),
/// while the max alone caps how finely the surface can resolve curvature. Each
/// facet's area is `½‖(b−a)×(c−a)‖`. Returns `None` for an empty mesh with no
/// `Tri3` faces.
fn mesh_max_triangle_area(mesh: &valenx_mesh::Mesh) -> Option<f64> {
    let mut largest: Option<f64> = None;
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            let area = 0.5 * (b - a).cross(&(c - a)).norm();
            largest = Some(largest.map_or(area, |s| s.max(area)));
        }
    }
    largest
}

/// The **total surface area** of a `Tri3` surface mesh (model units²) — the sum of
/// every facet area `Σ ½‖(b−a)×(c−a)‖`. Integrated over the *tessellation*, it is
/// the 2-D companion to the divergence-theorem [`mesh_volume`]: comparing it to the
/// CSG [`total_area`] of the source solid gauges tessellation fidelity (the two
/// agree for an exact mesh and diverge as faceting coarsens), and together with
/// `mesh_volume` it yields a mesh-only Wadell sphericity. It necessarily lies
/// between `n·`[`mesh_min_triangle_area`] and `n·`[`mesh_max_triangle_area`] for `n`
/// facets. Returns `None` for an empty mesh with no `Tri3` faces.
fn mesh_surface_area(mesh: &valenx_mesh::Mesh) -> Option<f64> {
    let mut total: Option<f64> = None;
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let a = mesh.nodes[tri[0] as usize];
            let b = mesh.nodes[tri[1] as usize];
            let c = mesh.nodes[tri[2] as usize];
            let area = 0.5 * (b - a).cross(&(c - a)).norm();
            total = Some(total.map_or(area, |s| s + area));
        }
    }
    total
}

/// Surface-area-to-volume ratio `S/V` (per model unit) — the compactness /
/// heat-exchange scaling figure: a sphere of radius `r` gives `3/r`, a cube of
/// side `a` gives `6/a`, and it falls as a body grows (the square–cube law).
/// `None` when the volume is non-positive.
fn surface_area_to_volume_ratio(area: f64, volume: f64) -> Option<f64> {
    if volume > 1e-12 {
        Some(area / volume)
    } else {
        None
    }
}

/// Bounding-box aspect ratio (slenderness) — the longest box edge over the
/// shortest, `max(dims) / min(dims)`. `1.0` for a cube; large for a slender rod
/// or a thin plate, where it flags poor machinability / 3-D-print stability and
/// a buckling tendency. `None` when the shortest extent is effectively zero (a
/// degenerate flat box).
fn bbox_aspect_ratio(dims: [f32; 3]) -> Option<f64> {
    let (mut lo, mut hi) = (f64::INFINITY, 0.0_f64);
    for d in dims {
        let d = d as f64;
        lo = lo.min(d);
        hi = hi.max(d);
    }
    if lo > 1e-12 {
        Some(hi / lo)
    } else {
        None
    }
}

/// Rebuild the feature tree against the parameters. Returns the per-step
/// snapshots (`snapshots[k]` = the set of bodies after step k) plus a one-line
/// status, or an error message.
fn rebuild_tree(s: &CadWorkbenchState) -> Result<(Vec<Vec<valenx_cad::Solid>>, String), String> {
    let table = build_table(&s.params);
    let tl = steps_to_timeline(&s.steps);
    let model = tl.rebuild(&table).map_err(|e| e.to_string())?;
    let nbodies = model.bodies.len();
    let faces: usize = model.bodies.iter().map(|b| b.faces()).sum();
    // B-rep Euler–Poincaré characteristic χ = V−E+F — the solid's true CAD topology
    // (distinct from the tessellation euler); shown only when there are B-rep faces.
    let brep_str = if faces > 0 {
        format!(" · brep χ {}", brep_euler_characteristic(&model.bodies))
    } else {
        String::new()
    };
    let volume = total_volume(&model.bodies);
    let area = total_area(&model.bodies);
    // Mass = density × solid volume (model mass units).
    let mass = s.density * volume;
    // Overall bounding-box size + fill fraction of the final model (best-effort:
    // these terms are dropped from the status if tessellation fails).
    let mesh = tessellate_step(&model.snapshots, model.snapshots.len(), &[]).ok();
    let dims = mesh.as_ref().and_then(mesh_dimensions);
    let bbox_str = match dims {
        Some(d @ [dx, dy, dz]) => {
            let ar = bbox_aspect_ratio(d)
                .map(|r| format!(" · AR {r:.1}:1"))
                .unwrap_or_default();
            format!(
                " · bbox {dx:.3}×{dy:.3}×{dz:.3} u (diag {:.3}{ar})",
                bbox_diagonal(d)
            )
        }
        None => String::new(),
    };
    let fill_str = dims
        .and_then(|d| fill_fraction(volume, d))
        .map(|f| format!(" · fill {:.0}%", 100.0 * f))
        .unwrap_or_default();
    // Wadell sphericity (1.0 = sphere) — a shape-compactness readout alongside fill.
    let sphericity_str = sphericity(volume, area)
        .map(|psi| format!(" · ψ {psi:.2}"))
        .unwrap_or_default();
    // Surface-area-to-volume ratio (square–cube law) — paired onto the area term.
    let sv_str = surface_area_to_volume_ratio(area, volume)
        .map(|sv| format!(" (S/V {sv:.3}/u)"))
        .unwrap_or_default();
    // Enclosed volume of the tessellated boundary (divergence theorem) — the zeroth
    // mass-moment leading the centroid/gyration/inertia suite, and a tessellation-
    // fidelity check against the solid `volume` above (equal for an exact mesh).
    let meshvol_str = mesh
        .as_ref()
        .and_then(mesh_volume)
        .map(|v| format!(" · mesh vol {v:.4} u³"))
        .unwrap_or_default();
    // Tessellated total surface area (Σ facet areas) — the 2-D companion to mesh vol
    // and a tessellation-fidelity check against the solid `area` above.
    let mesharea_str = mesh
        .as_ref()
        .and_then(mesh_surface_area)
        .map(|a| format!(" · mesh area {a:.3} u²"))
        .unwrap_or_default();
    // Volume centroid (centre of mass) of the tessellated model.
    let centroid_str = mesh
        .as_ref()
        .and_then(mesh_centroid)
        .map(|c| format!(" · centroid ({:.2}, {:.2}, {:.2}) u", c[0], c[1], c[2]))
        .unwrap_or_default();
    // Radius of gyration about the centroid — the inertia length scale (I ≈ m·k²).
    let gyration_str = mesh
        .as_ref()
        .and_then(mesh_radius_of_gyration)
        .map(|k| format!(" · k_gyr {k:.3} u"))
        .unwrap_or_default();
    // Principal moments of inertia about the centroid (I₁≥I₂≥I₃, unit density) —
    // the rotational-inertia completion of the mass-property suite.
    let moments_str = mesh
        .as_ref()
        .and_then(mesh_principal_moments)
        .map(|[i1, i2, i3]| format!(" · I ({i1:.3}, {i2:.3}, {i3:.3}) u⁵"))
        .unwrap_or_default();
    // Watertightness — a closed 2-manifold (every edge shared by exactly two
    // triangles) is the "is this a printable solid" gate. Dropped if tessellation failed.
    let watertight_str = mesh
        .as_ref()
        .map(|m| {
            if mesh_is_watertight(m) {
                " · watertight \u{2713}"
            } else {
                " · OPEN \u{2717}"
            }
        })
        .unwrap_or_default();
    // Boundary-edge count — quantifies the holes the OPEN flag warns of (the
    // perimeter to repair); shown only when the surface is actually open.
    let open_str = mesh
        .as_ref()
        .map(mesh_open_edge_count)
        .filter(|&n| n > 0)
        .map(|n| format!(" · {n} open edges"))
        .unwrap_or_default();
    // Distinct holes (boundary loops) — how many gaps to fill, vs the open-edge
    // count above (one triangular hole is 3 open edges); shown only when > 0.
    let hole_str = mesh
        .as_ref()
        .map(mesh_boundary_loop_count)
        .filter(|&n| n > 0)
        .map(|n| format!(" · {n} hole{}", if n == 1 { "" } else { "s" }))
        .unwrap_or_default();
    // Non-manifold edges (shared by ≥3 triangles) — the self-intersection / T-junction
    // defect, the other manifold failure beside the holes above; shown only when > 0.
    let nonmanifold_str = mesh
        .as_ref()
        .map(mesh_nonmanifold_edge_count)
        .filter(|&n| n > 0)
        .map(|n| format!(" · {n} non-manifold edge{}", if n == 1 { "" } else { "s" }))
        .unwrap_or_default();
    // Euler characteristic χ = V−E+F (and genus for a closed orientable surface)
    // — the topology the watertight flag has confirmed is closed.
    let euler_str = mesh
        .as_ref()
        .and_then(mesh_euler_characteristic)
        .map(|chi| {
            // genus g = (2−χ)/2 is meaningful only for a single closed orientable
            // surface (χ ≤ 2, even); a multi-shell model (χ > 2) just shows χ.
            if chi <= 2 && chi % 2 == 0 {
                format!(" · χ {chi} (genus {})", (2 - chi) / 2)
            } else {
                format!(" · χ {chi}")
            }
        })
        .unwrap_or_default();
    // Total angle defect ∑(2π−Σθ) = discrete total Gaussian curvature — by Gauss–
    // Bonnet it equals 2π·χ for a closed surface (4π for a genus-0 box), the
    // angle-based reading of the topology χ reports combinatorially.
    let defect_str = mesh
        .as_ref()
        .and_then(mesh_total_angle_defect)
        .map(|d| format!(" · ∑κ {d:.3}"))
        .unwrap_or_default();
    // Connected-shell count — how many disjoint solids the tessellation holds
    // (1 for a single connected body). Dropped if tessellation failed or empty.
    let shells_str = mesh
        .as_ref()
        .map(mesh_shell_count)
        .filter(|&n| n >= 1)
        .map(|n| format!(" · {n} shell{}", if n == 1 { "" } else { "s" }))
        .unwrap_or_default();
    // Minimal enclosing-sphere radius (Ritter) — the orientation-independent
    // tightest sphere around the part, the broad-phase culling/collision bound.
    let encl_str = mesh
        .as_ref()
        .and_then(mesh_bounding_sphere_radius)
        .map(|r| format!(" · r_enc {r:.3} u"))
        .unwrap_or_default();
    // Exact mesh diameter (longest chord between vertices) — the true caliper
    // extent, distinct from the axis-aligned bbox diagonal and the Ritter sphere.
    let diam_str = mesh
        .as_ref()
        .and_then(mesh_diameter)
        .map(|d| format!(" · diam {d:.3} u"))
        .unwrap_or_default();
    // Mean edge length — the characteristic element size (mesh resolution), the
    // typical discretisation scale distinct from the caliper extent above.
    let mean_edge_str = mesh
        .as_ref()
        .and_then(mesh_mean_edge_length)
        .map(|h| format!(" · mean edge {h:.3} u"))
        .unwrap_or_default();
    // Total edge (wireframe / 1-skeleton) length — the 1-D tessellation measure
    // completing the trio with mesh area (2-D) and mesh vol (3-D).
    let wire_str = mesh
        .as_ref()
        .and_then(mesh_total_edge_length)
        .map(|w| format!(" · wire {w:.3} u"))
        .unwrap_or_default();
    // Worst triangle shape quality (4√3·A/Σℓ²; 1 = equilateral, →0 = sliver) —
    // the mesh-quality gate, distinct from the shape/topology measures above.
    let quality_str = mesh
        .as_ref()
        .and_then(mesh_min_triangle_quality)
        .map(|q| format!(" · Qmin {q:.2}"))
        .unwrap_or_default();
    // Smallest facet area — the area-scale resolution/degeneracy flag (a near-zero
    // sliver wrecks conditioning & shading), vs the unitless Qmin shape score.
    let min_tri_area_str = mesh
        .as_ref()
        .and_then(mesh_min_triangle_area)
        .map(|a| format!(" · min tri area {a:.3} u²"))
        .unwrap_or_default();
    // Largest facet area — the coarsest element; with the min above it brackets the
    // mesh's area grading (a wide spread = non-uniform refinement).
    let max_tri_area_str = mesh
        .as_ref()
        .and_then(mesh_max_triangle_area)
        .map(|a| format!(" · max tri area {a:.3} u²"))
        .unwrap_or_default();
    // Worst triangle aspect ratio (longest/shortest edge; 1 = equilateral, ↑ for a
    // stretched sliver) — the elongation companion to the Qmin shape score.
    let aspect_str = mesh
        .as_ref()
        .and_then(mesh_max_aspect_ratio)
        .map(|r| format!(" · ARtri {r:.2}"))
        .unwrap_or_default();
    // Sharp (feature) edges — interior edges that fold past a 30° crease angle,
    // the chamfer/fillet candidates, distinct from the topology/quality measures.
    let sharp_str = mesh
        .as_ref()
        .map(|m| mesh_sharp_edge_count(m, std::f64::consts::FRAC_PI_6))
        .filter(|&n| n > 0)
        .map(|n| format!(" · {n} sharp edges"))
        .unwrap_or_default();
    // Sharpest fold — the maximum dihedral (crease) angle on the surface; the
    // single number that sizes the smallest fillet, vs the sharp-edge count.
    let crease_str = mesh
        .as_ref()
        .and_then(mesh_max_dihedral_angle)
        .map(|theta| format!(" · max crease {:.0}°", theta.to_degrees()))
        .unwrap_or_default();
    // Maximum vertex valence (graph-degree) — the connectivity-irregularity
    // indicator; a high value flags a "pole"/fan vertex, vs the geometric crease.
    let valence_str = mesh
        .as_ref()
        .and_then(mesh_max_vertex_valence)
        .map(|n| format!(" · max valence {n}"))
        .unwrap_or_default();
    // Mean vertex valence (average graph-degree) — the typical connectivity, → 6 for a
    // well-formed closed triangle mesh; the average-case companion to the max above.
    let mean_valence_str = mesh
        .as_ref()
        .and_then(mesh_mean_vertex_valence)
        .map(|v| format!(" · mean valence {v:.2}"))
        .unwrap_or_default();
    // B-rep shape descriptors of the final body (best-effort; blank if unavailable):
    // volume-equivalent sphere radius, the two Zingg bounding-box shape ratios, and the
    // rotation-invariant inertia anisotropy I_max/I_min.
    let shape_str = model
        .bodies
        .last()
        .map(|body| {
            let r_eq = valenx_cad::solid_equivalent_sphere_radius(body).unwrap_or(0.0);
            let elong = valenx_cad::solid_bounding_box_elongation(body).unwrap_or(0.0);
            let flat = valenx_cad::solid_bounding_box_flatness(body).unwrap_or(0.0);
            let i_aniso = valenx_cad::solid_inertia_anisotropy(body).unwrap_or(0.0);
            format!(
                " · r_eq {r_eq:.4} u · elongation {elong:.3} · flatness {flat:.3} · I_aniso {i_aniso:.3}"
            )
        })
        .unwrap_or_default();
    let mut status = format!(
        "{nbodies} bodies · {faces} faces{brep_str} · {volume:.4} u³ · {mass:.4} mass · {area:.4} u²{sv_str}{bbox_str}{fill_str}{sphericity_str}{meshvol_str}{mesharea_str}{centroid_str}{gyration_str}{moments_str}{watertight_str}{open_str}{hole_str}{nonmanifold_str}{euler_str}{defect_str}{shells_str}{encl_str}{diam_str}{mean_edge_str}{wire_str}{quality_str}{min_tri_area_str}{max_tri_area_str}{aspect_str}{sharp_str}{crease_str}{valence_str}{mean_valence_str} · {} steps",
        s.steps.len()
    );
    status.push_str(&shape_str);
    Ok((model.snapshots, status))
}

/// Tessellate the set of bodies at a given 1-based step into one display mesh
/// (the bodies are concatenated — they're separate solids drawn together).
fn tessellate_step(
    history: &[Vec<valenx_cad::Solid>],
    step_1based: usize,
    visible: &[bool],
) -> Result<valenx_mesh::Mesh, String> {
    let idx = step_1based
        .saturating_sub(1)
        .min(history.len().saturating_sub(1));
    let bodies = history.get(idx).ok_or_else(|| "no such step".to_string())?;
    let mut merged: Option<valenx_mesh::Mesh> = None;
    for (i, solid) in bodies.iter().enumerate() {
        // Empty `visible` ⇒ all bodies shown; otherwise honour the flag.
        if !visible.get(i).copied().unwrap_or(true) {
            continue;
        }
        let m = valenx_cad::solid_to_mesh(solid, valenx_cad::DEFAULT_TESS_TOLERANCE)
            .map_err(|e| e.to_string())?;
        merged = Some(match merged {
            None => m,
            Some(acc) => valenx_mesh::boolean::concatenate(&acc, &m),
        });
    }
    merged.ok_or_else(|| "no visible bodies to display".to_string())
}

/// Export the body set at a given 1-based step to a binary STL file. The actual
/// write is `valenx_mesh::stl_write::write_stl_binary` (an allowlisted streaming
/// export — see valenx-core's no-raw-fs-write guard).
fn export_stl(
    history: &[Vec<valenx_cad::Solid>],
    step_1based: usize,
    path: &std::path::Path,
) -> Result<(), String> {
    let mesh = tessellate_step(history, step_1based, &[])?;
    valenx_mesh::stl_write::write_stl_binary(&mesh, path).map_err(|e| e.to_string())
}

fn op_label(op: Op) -> &'static str {
    match op {
        Op::New => "New",
        Op::Join => "Join",
        Op::Cut => "Cut",
        Op::Intersect => "Intersect",
    }
}

fn kind_label(kind: FeatureKind) -> &'static str {
    match kind {
        FeatureKind::Box => "Box",
        FeatureKind::Cylinder => "Cylinder",
        FeatureKind::Sphere => "Sphere",
        FeatureKind::Cone => "Cone",
        FeatureKind::Torus => "Torus",
        FeatureKind::Extrude => "Extrude",
        FeatureKind::Revolve => "Revolve",
    }
}

/// Parse a boolean-op name (for the agent `SetControl` bridge) into an [`Op`].
/// Case-insensitive; accepts the combo's menu words. Fail-loud on an
/// unrecognised name so a typo is a `warn` note, not a silent no-op.
fn parse_op(s: &str) -> Result<Op, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "new" => Ok(Op::New),
        "join" | "union" | "add" => Ok(Op::Join),
        "cut" | "subtract" | "difference" => Ok(Op::Cut),
        "intersect" | "intersection" => Ok(Op::Intersect),
        other => Err(format!(
            "unknown op '{other}' (expected 'New', 'Join', 'Cut', or 'Intersect')"
        )),
    }
}

/// Parse a primitive-shape name (for the agent `SetControl` bridge) into a
/// [`FeatureKind`]. Case-insensitive; accepts the combo's menu words. Fail-loud.
fn parse_feature_kind(s: &str) -> Result<FeatureKind, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "box" | "cube" => Ok(FeatureKind::Box),
        "cylinder" | "cyl" => Ok(FeatureKind::Cylinder),
        "sphere" => Ok(FeatureKind::Sphere),
        "cone" => Ok(FeatureKind::Cone),
        "torus" => Ok(FeatureKind::Torus),
        "extrude" => Ok(FeatureKind::Extrude),
        "revolve" => Ok(FeatureKind::Revolve),
        other => Err(format!(
            "unknown shape '{other}' (expected 'Box', 'Cylinder', 'Sphere', \
             'Cone', 'Torus', 'Extrude', or 'Revolve')"
        )),
    }
}

/// Parse a sketch-tool name (for the agent `SetControl` bridge) into a
/// [`SketchTool`]. Case-insensitive; accepts the tool-row words. Fail-loud.
fn parse_sketch_tool(s: &str) -> Result<SketchTool, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "line" => Ok(SketchTool::Line),
        "arc" => Ok(SketchTool::Arc),
        "spline" => Ok(SketchTool::Spline),
        other => Err(format!(
            "unknown sketch tool '{other}' (expected 'Line', 'Arc', or 'Spline')"
        )),
    }
}

/// A narrow single-line editor for a dimension / placement expression, with a
/// compact leading `caption` rendered as its label and associated via
/// `labelled_by` so each field carries a unique accessible name (a bare
/// `TextEdit` has no own-name, leaving it unaddressable by a screen reader / AI
/// driver — closing AI-drivability gap 2).
fn dim_edit(ui: &mut egui::Ui, caption: &str, v: &mut String) {
    let cap = ui.label(egui::RichText::new(caption).weak().small());
    ui.add(egui::TextEdit::singleline(v).desired_width(52.0))
        .labelled_by(cap.id);
}

/// Format a float as a compact expression string for a feature field — up to 4
/// decimal places with trailing zeros trimmed (so `1.0 → "1"`, `0.25 → "0.25"`).
fn fmt_num(v: f64) -> String {
    let s = format!("{v:.4}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s.is_empty() || s == "-" {
        "0".to_string()
    } else {
        s.to_string()
    }
}

/// Rebuild the feature tree and stage the final solid's mesh for the central
/// 3-D viewport. This is the shared implementation of the panel's
/// "Rebuild → viewport" button and the programmatic [`CadWorkbenchState::
/// request_rebuild`] path (Part Design menu), so both behave identically:
/// on success it caches the per-step history, sets the scrubber to the last
/// step, and flags `push_rebuild`; on failure it records the error in
/// `tree_status` and clears any stale rebuild.
fn perform_rebuild(s: &mut CadWorkbenchState) {
    match rebuild_tree(s) {
        Ok((history, status)) => {
            let k = history.len();
            let nbodies = history.last().map_or(0, |b| b.len());
            let visible = vec![true; nbodies];
            match tessellate_step(&history, k, &visible) {
                Ok(mesh) => {
                    s.rebuilt_mesh = Some(mesh);
                    s.push_rebuild = true;
                    s.tree_status = Some(Ok(status));
                    s.scrub = k;
                    s.history = Some(history);
                    s.body_visible = visible;
                }
                Err(e) => {
                    s.history = None;
                    s.rebuilt_mesh = None;
                    s.push_rebuild = false;
                    s.tree_status = Some(Err(e));
                }
            }
        }
        Err(e) => {
            s.history = None;
            s.rebuilt_mesh = None;
            s.push_rebuild = false;
            s.tree_status = Some(Err(e));
        }
    }
}

/// Draw the interactive 2-D **sketch canvas** — a Fusion-style draw-then-extrude
/// surface supporting **Line / Arc / Spline** tools. Clicks (optionally snapped
/// to the grid) drop straight segments, 3-point circular arcs, or Catmull-Rom
/// splines; clicking back near the first anchor closes the loop. The in-progress
/// profile renders live (straight segments, smooth tessellated arcs/splines,
/// vertex + via/control-point dots, a translucent fill when closed, and a
/// rubber-band preview to the cursor). **Undo** / **Clear sketch** / a
/// **snap-to-grid** checkbox edit the sketch; **Extrude sketch** (enabled once
/// the profile encloses an area) tessellates the curved profile and sweeps it +
/// the height into the feature tree, rebuilding the viewport. Every control is a
/// named, AI-drivable egui widget. All state lives on [`CadWorkbenchState`];
/// nothing touches the app god-struct.
fn draw_sketch_canvas(s: &mut CadWorkbenchState, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(match s.sketch_tool {
            SketchTool::Line => "Line: click to drop points · click the first point to close",
            SketchTool::Arc => "Arc: click start, click end, click a 3rd point on the arc",
            SketchTool::Spline => "Spline: click points · Finish curve to commit the spline",
        })
        .weak()
        .small(),
    );

    // ── Tool selector (Line | Arc | Spline) ─────────────────────────────────
    ui.horizontal(|ui| {
        ui.label("Tool:");
        for (tool, name) in [
            (SketchTool::Line, "Line"),
            (SketchTool::Arc, "Arc"),
            (SketchTool::Spline, "Spline"),
        ] {
            if ui
                .selectable_label(s.sketch_tool == tool, name)
                .on_hover_text(match tool {
                    SketchTool::Line => "Straight segments — one click per vertex.",
                    SketchTool::Arc => "3-point circular arc — start, end, a point on the arc.",
                    SketchTool::Spline => "Smooth Catmull-Rom curve through clicked points.",
                })
                .clicked()
                && s.sketch_tool != tool
            {
                // Switching tool abandons any half-entered multi-click segment.
                s.sketch_tool = tool;
                s.sketch_pending.clear();
            }
        }
    });

    // ── Controls row ────────────────────────────────────────────────────────
    ui.horizontal_wrapped(|ui| {
        ui.checkbox(&mut s.sketch_grid_snap, "Snap to grid")
            .on_hover_text(format!("Snap clicks to the {SKETCH_SNAP}-unit grid."));
        if ui
            .button("Undo")
            .on_hover_text("Remove the most recently placed segment (or pending click).")
            .clicked()
            && s.sketch_pending.pop().is_none()
        {
            // No staged click to drop → remove the last committed segment.
            s.sketch_undo();
        }
        if ui
            .button("Clear sketch")
            .on_hover_text("Discard the whole drawn profile.")
            .clicked()
        {
            s.sketch_clear();
        }
        // Spline runs are open-ended: an explicit button commits the pending
        // points as one smooth curve (≥1 staged point on top of an anchor, or
        // ≥2 staged points when the sketch is still empty).
        let can_finish_spline = s.sketch_tool == SketchTool::Spline
            && !s.sketch_closed
            && (s.sketch_pending.len() >= 2
                || (s.sketch_start.is_some() && !s.sketch_pending.is_empty()));
        if ui
            .add_enabled(can_finish_spline, egui::Button::new("Finish curve"))
            .on_hover_text("Commit the staged points as one Catmull-Rom spline segment.")
            .clicked()
        {
            let pts = std::mem::take(&mut s.sketch_pending);
            s.sketch_add_spline(&pts);
        }
    });

    // ── Canvas ──────────────────────────────────────────────────────────────
    let side = ui.available_width().clamp(160.0, 360.0);
    let (resp, painter) = ui.allocate_painter(egui::vec2(side, side), egui::Sense::click());
    let rect = resp.rect;
    painter.rect_filled(rect, 3.0, egui::Color32::from_gray(18));

    // Model (y-up, centred, ±SKETCH_VIEW window) ↔ screen (y-down) mapping.
    let half = rect.width().min(rect.height()) * 0.5;
    let c = rect.center();
    let scale = (half / SKETCH_VIEW as f32).max(1.0);
    let to_screen = |p: [f64; 2]| -> egui::Pos2 {
        egui::pos2(c.x + p[0] as f32 * scale, c.y - p[1] as f32 * scale)
    };
    let to_model = |pos: egui::Pos2| -> [f64; 2] {
        [
            ((pos.x - c.x) / scale) as f64,
            ((c.y - pos.y) / scale) as f64,
        ]
    };
    let snap = |v: [f64; 2]| -> [f64; 2] {
        if s.sketch_grid_snap {
            [
                (v[0] / SKETCH_SNAP).round() * SKETCH_SNAP,
                (v[1] / SKETCH_SNAP).round() * SKETCH_SNAP,
            ]
        } else {
            v
        }
    };

    // Light grid on the SKETCH_SNAP step, brighter X/Y axes.
    let grid_stroke = egui::Stroke::new(0.5, egui::Color32::from_gray(40));
    let axis_stroke = egui::Stroke::new(1.0, egui::Color32::from_gray(90));
    let mut g = -SKETCH_VIEW;
    while g <= SKETCH_VIEW + 1e-9 {
        let sx = to_screen([g, 0.0]).x;
        let sy = to_screen([0.0, g]).y;
        painter.line_segment(
            [egui::pos2(sx, rect.top()), egui::pos2(sx, rect.bottom())],
            grid_stroke,
        );
        painter.line_segment(
            [egui::pos2(rect.left(), sy), egui::pos2(rect.right(), sy)],
            grid_stroke,
        );
        g += SKETCH_SNAP;
    }
    let origin = to_screen([0.0, 0.0]);
    painter.line_segment(
        [
            egui::pos2(rect.left(), origin.y),
            egui::pos2(rect.right(), origin.y),
        ],
        axis_stroke,
    );
    painter.line_segment(
        [
            egui::pos2(origin.x, rect.top()),
            egui::pos2(origin.x, rect.bottom()),
        ],
        axis_stroke,
    );

    // Handle a click. Closing the loop (near the first anchor) takes precedence;
    // otherwise the active tool consumes the (snapped) point.
    if resp.clicked() {
        if let Some(pos) = resp.interact_pointer_pos() {
            let first_screen = s.sketch_start.map(to_screen);
            let near_first = first_screen
                .map(|f| f.distance(pos) <= SKETCH_PICK_PX)
                .unwrap_or(false);
            let m = snap(to_model(pos));
            if near_first && s.sketch_anchor_count() >= 3 && s.sketch_pending.is_empty() {
                s.sketch_closed = true;
            } else if !s.sketch_closed {
                match s.sketch_tool {
                    SketchTool::Line => s.sketch_add_point(m[0], m[1]),
                    SketchTool::Arc => {
                        s.sketch_pending.push(m);
                        // For the first arc we need 3 clicks (start, end, via);
                        // once a chain exists the previous anchor is the start,
                        // so only end + via are needed.
                        let need = if s.sketch_start.is_none() { 3 } else { 2 };
                        if s.sketch_pending.len() >= need {
                            let p = std::mem::take(&mut s.sketch_pending);
                            // p = [start, end, via] (fresh) or [end, via] (chained).
                            let (start, end, via) = if need == 3 {
                                (p[0], p[1], p[2])
                            } else {
                                let start = s.sketch_points().last().copied().unwrap_or(p[0]);
                                (start, p[0], p[1])
                            };
                            s.sketch_add_arc(start, via, end);
                        }
                    }
                    SketchTool::Spline => {
                        // Accumulate; the user commits the run with "Finish curve".
                        s.sketch_pending.push(m);
                    }
                }
            }
        }
    }

    // ── Render ────────────────────────────────────────────────────────────────
    let line_stroke = egui::Stroke::new(1.5, egui::Color32::from_rgb(120, 200, 255));
    // The fully tessellated committed profile (curves smoothed). Drives both the
    // outline and the closed-loop fill, so arcs/splines show as smooth curves.
    let profile: Vec<egui::Pos2> = s
        .sketch_profile_polyline()
        .iter()
        .map(|&p| to_screen(p))
        .collect();
    if s.sketch_closed && profile.len() >= 3 {
        // Translucent fill for the closed loop (the tessellated outline may be
        // concave, so fill via a triangle fan around the centroid).
        let cx = profile.iter().map(|p| p.x).sum::<f32>() / profile.len() as f32;
        let cy = profile.iter().map(|p| p.y).sum::<f32>() / profile.len() as f32;
        let centroid = egui::pos2(cx, cy);
        let fill = egui::Color32::from_rgba_unmultiplied(120, 200, 255, 40);
        for w in profile.windows(2) {
            painter.add(egui::Shape::convex_polygon(
                vec![centroid, w[0], w[1]],
                fill,
                egui::Stroke::NONE,
            ));
        }
        painter.add(egui::Shape::convex_polygon(
            vec![centroid, *profile.last().unwrap(), profile[0]],
            fill,
            egui::Stroke::NONE,
        ));
    }
    // Outline (smooth curves) + the closing edge or a rubber-band preview.
    for w in profile.windows(2) {
        painter.line_segment([w[0], w[1]], line_stroke);
    }
    if s.sketch_closed && profile.len() >= 3 {
        painter.line_segment([*profile.last().unwrap(), profile[0]], line_stroke);
    } else {
        // Rubber-band preview from the last committed anchor to the cursor (only
        // while no multi-click segment is mid-entry).
        let band = (!s.sketch_closed && s.sketch_pending.is_empty())
            .then(|| profile.last().copied().zip(resp.hover_pos()))
            .flatten();
        if let Some((last, cursor)) = band {
            painter.line_segment(
                [last, cursor],
                egui::Stroke::new(1.0, egui::Color32::from_gray(120)),
            );
        }
    }

    // Pending (un-committed) clicks: a dashed-ish preview chain in amber so the
    // user sees the arc-via / spline points building up.
    let pending_stroke = egui::Stroke::new(1.0, egui::Color32::from_rgb(255, 180, 80));
    if !s.sketch_pending.is_empty() {
        // Chain from the current end anchor (if any) through the pending points.
        let mut chain: Vec<egui::Pos2> = Vec::new();
        if let Some(&anchor) = s.sketch_points().last() {
            chain.push(to_screen(anchor));
        }
        chain.extend(s.sketch_pending.iter().map(|&p| to_screen(p)));
        for w in chain.windows(2) {
            painter.line_segment([w[0], w[1]], pending_stroke);
        }
        for &p in &s.sketch_pending {
            painter.circle_filled(to_screen(p), 3.0, egui::Color32::from_rgb(255, 180, 80));
        }
        // Live spline preview through anchor + pending points.
        if s.sketch_tool == SketchTool::Spline && chain.len() >= 2 {
            let ctrl: Vec<[f64; 2]> = {
                let mut v = Vec::new();
                if let Some(&anchor) = s.sketch_points().last() {
                    v.push(anchor);
                }
                v.extend(s.sketch_pending.iter().copied());
                v
            };
            let mut curve = sample_spline(&ctrl);
            curve.push(*ctrl.last().unwrap());
            let cs: Vec<egui::Pos2> = curve.iter().map(|&p| to_screen(p)).collect();
            for w in cs.windows(2) {
                painter.line_segment([w[0], w[1]], pending_stroke);
            }
        }
    }

    // Anchor dots: the first is highlighted amber as the loop-closing target.
    for (i, &p) in s.sketch_points().iter().enumerate() {
        let col = if i == 0 {
            egui::Color32::from_rgb(255, 210, 90)
        } else {
            egui::Color32::from_rgb(120, 200, 255)
        };
        painter.circle_filled(to_screen(p), 3.0, col);
    }

    // ── Status + extrude controls ───────────────────────────────────────────
    let n = s.sketch_anchor_count();
    let pend = s.sketch_pending.len();
    ui.label(
        egui::RichText::new(format!(
            "{n} anchor{}{} · {}",
            if n == 1 { "" } else { "s" },
            if pend > 0 {
                format!(" (+{pend} pending)")
            } else {
                String::new()
            },
            if s.sketch_closed { "closed" } else { "open" }
        ))
        .small()
        .monospace(),
    );
    ui.horizontal(|ui| {
        // Associate the DragValue with its caption via `labelled_by`, so the
        // spin button is named for a screen reader / AI driver (egui clears a
        // DragValue's own Name).
        let cap = ui.label("Extrude height");
        ui.add(
            egui::DragValue::new(&mut s.sketch_extrude_height)
                .speed(0.1)
                .range(0.01..=1.0e6)
                .suffix(" u"),
        )
        .labelled_by(cap.id);
    });
    // A profile encloses area once its tessellated polyline has ≥3 points; arcs
    // and splines reach that with as few as 2 anchors.
    let can_extrude = s.sketch_profile_polyline().len() >= 3;
    if ui
        .add_enabled(can_extrude, egui::Button::new("Extrude sketch"))
        .on_hover_text(
            "Tessellate the drawn profile (arcs/splines sampled into a smooth \
             polyline), sweep it along +Z by the height above into a solid, add \
             it to the feature tree, and rebuild the viewport.",
        )
        .clicked()
    {
        let pts = s.sketch_profile_polyline();
        let h = s.sketch_extrude_height;
        s.add_extrude_from_sketch(&pts, h);
        perform_rebuild(s);
    }
}

/// Draw the parametric-CAD workbench (a no-op unless toggled on via
/// View → Parametric CAD).
pub fn draw_cad_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_cad_workbench {
        return;
    }
    // A programmatic rebuild requested this frame (e.g. from the top-bar Part
    // Design menu, which sets the flag before this panel draws). Run it before
    // the panel so the staged mesh is pushed into the viewport this frame.
    if std::mem::take(&mut app.cad.rebuild_request) {
        perform_rebuild(&mut app.cad);
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_cad_workbench",
        "Parametric CAD",
        |app, ui| {
            ui.label(
                egui::RichText::new("named parameters · valenx-solvespace-3d")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.cad;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    // ── Robotics: native in-house URDF import ────────────
                    ui.label(egui::RichText::new("Robotics — import URDF (native)").strong());
                    if ui
                        .button(egui::RichText::new("Load 16-DOF robot hand").strong())
                        .clicked()
                    {
                        match valenx_assembly::import_urdf(&valenx_assembly::demo_hand_urdf()) {
                            Ok(robot) => {
                                let mesh =
                                    valenx_assembly::assembly_to_mesh(&robot.assembly, 0.005);
                                s.tree_status = Some(Ok(format!(
                                    "loaded {} — {} parts, {} joints (no adapter)",
                                    robot.name,
                                    robot.assembly.parts.len(),
                                    robot.assembly.joints.len(),
                                )));
                                s.rebuilt_mesh = Some(mesh);
                                s.push_rebuild = true;
                            }
                            Err(e) => {
                                s.tree_status = Some(Err(format!("URDF import failed: {e}")));
                            }
                        }
                    }
                    ui.label(
                        egui::RichText::new(
                            "parses a URDF robot natively (no adapter) into the 3-D viewport",
                        )
                        .weak()
                        .small(),
                    );
                    ui.add_space(6.0);
                    ui.separator();

                    ui.label(egui::RichText::new("Parameters (name = expression)").strong());
                    let mut remove: Option<usize> = None;
                    for (i, (name, expr)) in s.params.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            // Per-row captions so the name / expr editors carry
                            // unique accessible names (a bare TextEdit has no
                            // own-name) — AI-drivable / screen-reader addressable.
                            let name_cap = ui.label(
                                egui::RichText::new(format!("param {} name", i + 1))
                                    .weak()
                                    .small(),
                            );
                            ui.add(
                                egui::TextEdit::singleline(name)
                                    .desired_width(80.0)
                                    .hint_text("name"),
                            )
                            .labelled_by(name_cap.id);
                            let expr_cap = ui.label("=");
                            ui.add(
                                egui::TextEdit::singleline(expr)
                                    .desired_width(130.0)
                                    .hint_text("expr"),
                            )
                            .labelled_by(expr_cap.id);
                            if ui.small_button("✕").clicked() {
                                remove = Some(i);
                            }
                        });
                    }
                    if let Some(i) = remove {
                        s.params.remove(i);
                    }
                    if ui.button("+ parameter").clicked() {
                        s.params.push((String::new(), String::new()));
                    }

                    // ---- Material: density → mass readout ----
                    ui.add_space(6.0);
                    ui.horizontal(|ui| {
                        // Associate the DragValue with its caption via
                        // `labelled_by` so it is named for a screen reader / AI
                        // driver (egui clears a DragValue's own Name).
                        let cap = ui.label("Material density").on_hover_text(
                            "Mass per unit volume; the rebuild status shows \
                             mass = density × solid volume.",
                        );
                        ui.add(
                            egui::DragValue::new(&mut s.density)
                                .speed(0.1)
                                .range(0.0..=1.0e9)
                                .suffix(" /u³"),
                        )
                        .labelled_by(cap.id);
                    });

                    // ---- Sketch: parameter-driven circle radius ----
                    ui.add_space(6.0);
                    ui.label(egui::RichText::new("Sketch").strong());
                    ui.horizontal(|ui| {
                        let radius_cap = ui.label("circle radius =");
                        ui.add(
                            egui::TextEdit::singleline(&mut s.radius_param)
                                .desired_width(100.0)
                                .hint_text("parameter"),
                        )
                        .labelled_by(radius_cap.id);
                    });
                    if ui.button("▶ Solve").clicked() {
                        let res = run_cad(s);
                        s.results = Some(res);
                    }
                    if let Some(r) = &s.results {
                        ui.label(egui::RichText::new("Resolved").strong());
                        for (name, val) in &r.resolved {
                            match val {
                                Ok(v) => ui.label(
                                    egui::RichText::new(format!("{name} = {v:.4}"))
                                        .monospace()
                                        .small(),
                                ),
                                Err(e) => ui.colored_label(
                                    egui::Color32::from_rgb(220, 120, 80),
                                    egui::RichText::new(format!("{name}: {e}")).small(),
                                ),
                            };
                        }
                        let color = if r.solved_radius.is_some() {
                            egui::Color32::from_rgb(80, 220, 120)
                        } else {
                            egui::Color32::from_rgb(220, 120, 80)
                        };
                        ui.colored_label(color, &r.status);
                    }

                    // ---- Sketch canvas: draw a polygon → extrude → solid ----
                    ui.separator();
                    let sketch_header =
                        ui.label(egui::RichText::new("Sketch canvas (draw → extrude)").strong());
                    // The top-bar Part Design → Sketch menu item flags a focus
                    // request; scroll the canvas header into view this frame.
                    if std::mem::take(&mut s.sketch_focus_request) {
                        sketch_header.scroll_to_me(Some(egui::Align::TOP));
                    }
                    draw_sketch_canvas(s, ui);

                    // ---- Feature tree (CSG) ----
                    ui.separator();
                    ui.label(egui::RichText::new("Feature tree (CSG)").strong());
                    ui.label(
                        egui::RichText::new(
                            "each step places a primitive and combines it with the running body",
                        )
                        .weak()
                        .small(),
                    );

                    let mut remove_step: Option<usize> = None;
                    for (i, st) in s.steps.iter_mut().enumerate() {
                        ui.group(|ui| {
                            ui.horizontal(|ui| {
                                // Step-number label doubles as the boolean-op
                                // combo's accessible caption; the kind combo gets
                                // its own so both are addressable by name.
                                let step_cap = ui.label(format!("step {} op", i + 1));
                                egui::ComboBox::from_id_source(("cad_op", i))
                                    .selected_text(op_label(st.op))
                                    .width(92.0)
                                    .show_ui(ui, |ui| {
                                        for op in [Op::New, Op::Join, Op::Cut, Op::Intersect] {
                                            ui.selectable_value(&mut st.op, op, op_label(op));
                                        }
                                    })
                                    .response
                                    .labelled_by(step_cap.id);
                                let kind_cap =
                                    ui.label(egui::RichText::new("shape").weak().small());
                                egui::ComboBox::from_id_source(("cad_kind", i))
                                    .selected_text(kind_label(st.kind))
                                    .width(92.0)
                                    .show_ui(ui, |ui| {
                                        ui.selectable_value(&mut st.kind, FeatureKind::Box, "Box");
                                        ui.selectable_value(
                                            &mut st.kind,
                                            FeatureKind::Cylinder,
                                            "Cylinder",
                                        );
                                        ui.selectable_value(
                                            &mut st.kind,
                                            FeatureKind::Sphere,
                                            "Sphere",
                                        );
                                        ui.selectable_value(
                                            &mut st.kind,
                                            FeatureKind::Cone,
                                            "Cone",
                                        );
                                        ui.selectable_value(
                                            &mut st.kind,
                                            FeatureKind::Torus,
                                            "Torus",
                                        );
                                        ui.selectable_value(
                                            &mut st.kind,
                                            FeatureKind::Extrude,
                                            "Extrude",
                                        );
                                        ui.selectable_value(
                                            &mut st.kind,
                                            FeatureKind::Revolve,
                                            "Revolve",
                                        );
                                    })
                                    .response
                                    .labelled_by(kind_cap.id);
                                if ui.small_button("✕").clicked() {
                                    remove_step = Some(i);
                                }
                            });
                            match st.kind {
                                FeatureKind::Box => {
                                    ui.horizontal(|ui| {
                                        dim_edit(ui, "dx", &mut st.dx);
                                        dim_edit(ui, "dy", &mut st.dy);
                                        dim_edit(ui, "dz", &mut st.dz);
                                    });
                                }
                                FeatureKind::Cylinder => {
                                    ui.horizontal(|ui| {
                                        dim_edit(ui, "r", &mut st.radius);
                                        dim_edit(ui, "h", &mut st.height);
                                    });
                                }
                                FeatureKind::Sphere => {
                                    ui.horizontal(|ui| {
                                        dim_edit(ui, "r", &mut st.radius);
                                    });
                                }
                                FeatureKind::Cone => {
                                    ui.horizontal(|ui| {
                                        dim_edit(ui, "base r", &mut st.radius);
                                        dim_edit(ui, "top r", &mut st.top_radius);
                                        dim_edit(ui, "h", &mut st.height);
                                    });
                                }
                                FeatureKind::Torus => {
                                    ui.horizontal(|ui| {
                                        dim_edit(ui, "major", &mut st.major);
                                        dim_edit(ui, "minor", &mut st.minor);
                                    });
                                }
                                FeatureKind::Extrude => {
                                    ui.horizontal(|ui| {
                                        dim_edit(ui, "height", &mut st.height);
                                    });
                                    ui.label(
                                        egui::RichText::new("profile (x, y) — ≥3 points")
                                            .weak()
                                            .small(),
                                    );
                                    let mut rm_pt: Option<usize> = None;
                                    for (j, pt) in st.profile.iter_mut().enumerate() {
                                        ui.horizontal(|ui| {
                                            // Per-point captions associated via
                                            // `labelled_by` so each profile
                                            // spin button is uniquely named for
                                            // a screen reader / AI driver (egui
                                            // clears a DragValue's own Name).
                                            let xcap = ui.label(format!("P{j} x"));
                                            ui.add(
                                                egui::DragValue::new(&mut pt.0)
                                                    .speed(0.1)
                                                    .prefix("x "),
                                            )
                                            .labelled_by(xcap.id);
                                            let ycap = ui.label(format!("P{j} y"));
                                            ui.add(
                                                egui::DragValue::new(&mut pt.1)
                                                    .speed(0.1)
                                                    .prefix("y "),
                                            )
                                            .labelled_by(ycap.id);
                                            if ui.small_button("✕").clicked() {
                                                rm_pt = Some(j);
                                            }
                                        });
                                    }
                                    if let Some(j) = rm_pt {
                                        if st.profile.len() > 3 {
                                            st.profile.remove(j);
                                        }
                                    }
                                    if ui.small_button("+ point").clicked() {
                                        st.profile.push((0.0, 0.0));
                                    }
                                }
                                FeatureKind::Revolve => {
                                    ui.horizontal(|ui| {
                                        dim_edit(ui, "angle°", &mut st.angle);
                                    });
                                    ui.label(
                                        egui::RichText::new(
                                            "half-section (r, z) — ends on axis (r=0) cap the body",
                                        )
                                        .weak()
                                        .small(),
                                    );
                                    let mut rm_pt: Option<usize> = None;
                                    for (j, pt) in st.profile.iter_mut().enumerate() {
                                        ui.horizontal(|ui| {
                                            // Per-point captions associated via
                                            // `labelled_by` so each half-section
                                            // spin button is uniquely named for
                                            // a screen reader / AI driver (egui
                                            // clears a DragValue's own Name).
                                            let rcap = ui.label(format!("P{j} r"));
                                            ui.add(
                                                egui::DragValue::new(&mut pt.0)
                                                    .speed(0.1)
                                                    .prefix("r ")
                                                    .range(0.0..=f64::INFINITY),
                                            )
                                            .labelled_by(rcap.id);
                                            let zcap = ui.label(format!("P{j} z"));
                                            ui.add(
                                                egui::DragValue::new(&mut pt.1)
                                                    .speed(0.1)
                                                    .prefix("z "),
                                            )
                                            .labelled_by(zcap.id);
                                            if ui.small_button("✕").clicked() {
                                                rm_pt = Some(j);
                                            }
                                        });
                                    }
                                    if let Some(j) = rm_pt {
                                        if st.profile.len() > 2 {
                                            st.profile.remove(j);
                                        }
                                    }
                                    if ui.small_button("+ point").clicked() {
                                        st.profile.push((0.0, 0.0));
                                    }
                                }
                            }
                            ui.horizontal(|ui| {
                                ui.label("at");
                                dim_edit(ui, "x", &mut st.x);
                                dim_edit(ui, "y", &mut st.y);
                                dim_edit(ui, "z", &mut st.z);
                            });
                            ui.horizontal(|ui| {
                                ui.label("rot");
                                dim_edit(ui, "rx°", &mut st.rx);
                                dim_edit(ui, "ry°", &mut st.ry);
                                dim_edit(ui, "rz°", &mut st.rz);
                            });
                        });
                    }
                    if let Some(i) = remove_step {
                        s.steps.remove(i);
                    }
                    ui.horizontal_wrapped(|ui| {
                        // Each button calls the same public mutator the top-bar
                        // Part Design menu uses, so panel and menu share one path.
                        if ui.button("+ Box").clicked() {
                            s.add_box();
                        }
                        if ui.button("+ Cylinder").clicked() {
                            s.add_cylinder();
                        }
                        if ui.button("+ Sphere").clicked() {
                            s.add_sphere();
                        }
                        if ui.button("+ Cone").clicked() {
                            s.add_cone();
                        }
                        if ui.button("+ Torus").clicked() {
                            s.add_torus();
                        }
                        if ui.button("+ Extrude").clicked() {
                            s.add_extrude();
                        }
                        if ui.button("+ Revolve").clicked() {
                            s.add_revolve();
                        }
                    });
                    ui.horizontal(|ui| {
                        if ui.button("💾 Save…").clicked() {
                            match save_string(s) {
                                Ok(txt) => {
                                    if let Some(path) = rfd::FileDialog::new()
                                        .add_filter("valenx CAD tree", &["ron"])
                                        .set_file_name("feature-tree.ron")
                                        .save_file()
                                    {
                                        // Durable small-state save → crash-safe
                                        // atomic write (sidecar + fsync + rename),
                                        // per the no-raw-fs-write guard.
                                        let res =
                                            valenx_core::io_caps::atomic_write_str(&path, &txt);
                                        s.tree_status = Some(match res {
                                            Ok(()) => Ok(format!("saved {}", path.display())),
                                            Err(e) => Err(format!("save failed: {e}")),
                                        });
                                    }
                                }
                                Err(e) => {
                                    s.tree_status = Some(Err(format!("serialize failed: {e}")));
                                }
                            }
                        }
                        if ui.button("📂 Load…").clicked() {
                            if let Some(path) = rfd::FileDialog::new()
                                .add_filter("valenx CAD tree", &["ron"])
                                .pick_file()
                            {
                                let loaded = std::fs::read_to_string(&path)
                                    .map_err(|e| e.to_string())
                                    .and_then(|txt| load_from_string(&txt));
                                match loaded {
                                    Ok((params, steps)) => {
                                        s.params = params;
                                        s.steps = steps;
                                        s.history = None;
                                        s.body_visible = Vec::new();
                                        s.tree_status = Some(Ok(format!(
                                            "loaded {} ({} steps)",
                                            path.display(),
                                            s.steps.len()
                                        )));
                                    }
                                    Err(e) => {
                                        s.tree_status = Some(Err(format!("load failed: {e}")));
                                    }
                                }
                            }
                        }
                        if ui.button("⭳ Export STL…").clicked() {
                            match rebuild_tree(s) {
                                Ok((history, _status)) => {
                                    if let Some(path) = rfd::FileDialog::new()
                                        .add_filter("STL mesh", &["stl"])
                                        .set_file_name("valenx-cad.stl")
                                        .save_file()
                                    {
                                        let res = export_stl(&history, history.len(), &path);
                                        s.tree_status = Some(match res {
                                            Ok(()) => Ok(format!("exported {}", path.display())),
                                            Err(e) => Err(format!("export failed: {e}")),
                                        });
                                    }
                                }
                                Err(e) => {
                                    s.tree_status = Some(Err(format!("rebuild failed: {e}")));
                                }
                            }
                        }
                    });
                    if ui.button("▶ Rebuild → viewport").clicked() {
                        perform_rebuild(s);
                    }
                    if let Some(res) = &s.tree_status {
                        match res {
                            Ok(status) => {
                                ui.colored_label(egui::Color32::from_rgb(80, 220, 120), status)
                            }
                            Err(e) => ui.colored_label(
                                egui::Color32::from_rgb(220, 120, 80),
                                format!("rebuild failed: {e}"),
                            ),
                        };
                    }

                    // History scrubber — roll the model back/forward through
                    // the per-step snapshots from the last rebuild, pushing the
                    // selected step into the viewport.
                    let n = s.history.as_ref().map_or(0, |h| h.len());
                    if n > 1 {
                        ui.add_space(2.0);
                        ui.label(egui::RichText::new("History").strong());
                        let resp = ui.add(
                            egui::Slider::new(&mut s.scrub, 1..=n)
                                .integer()
                                .text("step"),
                        );
                        if resp.changed() {
                            let scrub = s.scrub;
                            let mesh = s
                                .history
                                .as_ref()
                                .and_then(|h| tessellate_step(h, scrub, &s.body_visible).ok());
                            if let Some(mesh) = mesh {
                                s.rebuilt_mesh = Some(mesh);
                                s.push_rebuild = true;
                            }
                        }
                        let label = s
                            .history
                            .as_ref()
                            .and_then(|h| h.get(s.scrub.saturating_sub(1)))
                            .map(|bodies| {
                                let faces: usize = bodies.iter().map(|b| b.faces()).sum();
                                let vol = total_volume(bodies);
                                format!(
                                    "step {} / {n} — {} bodies, {faces} faces, {vol:.4} u³",
                                    s.scrub,
                                    bodies.len()
                                )
                            });
                        if let Some(label) = label {
                            ui.label(egui::RichText::new(label).small().monospace());
                        }
                    }

                    // Per-body visibility (multi-body): toggle which bodies of
                    // the final set are shown in the viewport.
                    if s.body_visible.len() > 1 {
                        ui.add_space(2.0);
                        ui.label(egui::RichText::new("Bodies").strong());
                        let mut changed = false;
                        ui.horizontal_wrapped(|ui| {
                            for (i, vis) in s.body_visible.iter_mut().enumerate() {
                                if ui.checkbox(vis, format!("{}", i + 1)).changed() {
                                    changed = true;
                                }
                            }
                        });
                        if changed {
                            let scrub = s.scrub;
                            let mesh = s
                                .history
                                .as_ref()
                                .and_then(|h| tessellate_step(h, scrub, &s.body_visible).ok());
                            if let Some(mesh) = mesh {
                                s.rebuilt_mesh = Some(mesh);
                                s.push_rebuild = true;
                            }
                        }
                    }
                });
        },
    );
    if close {
        app.show_cad_workbench = false;
    }

    // Deferred (outside the panel borrow): push the rebuilt solid's mesh into
    // the central 3-D viewport.
    if app.cad.push_rebuild {
        app.cad.push_rebuild = false;
        if let Some(mesh) = app.cad.rebuilt_mesh.take() {
            let quality = valenx_mesh::quality_report(&mesh);
            let aspect_hist =
                valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
            let skew_hist =
                valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
            app.stl = None;
            app.aero_field_overlay = None;
            app.mesh = Some(LoadedMesh {
                path: std::path::PathBuf::from("<cad>/feature-tree"),
                mesh,
                quality,
                aspect_hist,
                skew_hist,
            });
            app.frame_current_mesh();
        }
    }
}

/// The agent-bridge product for the parametric CSG feature-tree workbench
/// (`show_3d{kind="cad"}`).
///
/// Rebuilds the **default feature tree** — valenx-cad's canonical punched cube
/// (a unit box with a Ø0.5 cylinder cut through it) — into BRep solid(s) via
/// [`rebuild_tree`], then tessellates the final step to a `Tri3`
/// [`valenx_mesh::Mesh`] with [`tessellate_step`] (which calls
/// `valenx_cad::solid_to_mesh`). Pure and app-state-free: it constructs a fresh
/// [`CadWorkbenchState::default`], so it produces exactly the geometry the
/// workbench shows on first open. The readout rows report the solid's volume,
/// surface area and triangle count.
pub(crate) fn cad_product() -> crate::WorkspaceProduct {
    let state = CadWorkbenchState::default();
    let built = (|| -> Result<(valenx_mesh::Mesh, f64, f64), String> {
        let (snapshots, _status) = rebuild_tree(&state)?;
        let bodies = snapshots.last().cloned().unwrap_or_default();
        let volume = total_volume(&bodies);
        let area = total_area(&bodies);
        let mesh = tessellate_step(&snapshots, snapshots.len(), &[])?;
        Ok((mesh, volume, area))
    })();
    let (mesh, lines) = match built {
        Ok((mesh, volume, area)) => {
            let tris = mesh.total_elements();
            let lines = vec![
                "parametric CSG: unit box − Ø0.5 through-cylinder".to_string(),
                format!("volume {volume:.4} u³ · surface area {area:.4} u²"),
                format!("BRep tessellation: {tris} triangles"),
            ];
            (mesh, lines)
        }
        Err(e) => {
            // Theoretically unreachable for the canonical default tree; degrade
            // to a tiny placeholder triangle + a note rather than panicking.
            let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
            block.connectivity = vec![0, 1, 2];
            let mut placeholder = valenx_mesh::Mesh::new("valenx-cad-part");
            placeholder.nodes = vec![
                nalgebra::Vector3::new(0.0, 0.0, 0.0),
                nalgebra::Vector3::new(1.0, 0.0, 0.0),
                nalgebra::Vector3::new(0.0, 1.0, 0.0),
            ];
            placeholder.element_blocks.push(block);
            placeholder.recompute_stats();
            (
                placeholder,
                vec![
                    "parametric CSG feature tree".to_string(),
                    format!("rebuild unavailable — showing placeholder ({e})"),
                ],
            )
        }
    };
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<cad>/feature-tree");
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "CAD part (parametric CSG)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        animation: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parameter_drives_circle_radius() {
        // Defaults: base = 4, radius = base + 1, circle radius driven by `radius`.
        let s = CadWorkbenchState::default();
        let r = run_cad(&s);
        let solved = r.solved_radius.expect("sketch solved");
        assert!((solved - 5.0).abs() < 1e-4, "radius = {solved}");
        let radius = r.resolved.iter().find(|(n, _)| n == "radius").unwrap();
        assert_eq!(
            radius.1.as_ref().ok().map(|v| (v * 1e4).round() / 1e4),
            Some(5.0)
        );
    }

    #[test]
    fn cyclic_parameters_report_an_error_not_a_panic() {
        let s = CadWorkbenchState {
            params: vec![
                ("a".to_string(), "b + 1".to_string()),
                ("b".to_string(), "a + 1".to_string()),
            ],
            radius_param: "a".to_string(),
            ..CadWorkbenchState::default()
        };
        let r = run_cad(&s);
        assert!(r.solved_radius.is_none());
        assert!(
            r.status.to_lowercase().contains("cyclic"),
            "status: {}",
            r.status
        );
    }

    #[test]
    fn feature_tree_rebuilds_punched_cube_with_history() {
        // The default tree is New box + Cut cylinder — the punched cube.
        let s = CadWorkbenchState::default();
        let (history, status) = rebuild_tree(&s).expect("default tree rebuilds");
        assert_eq!(history.len(), 2, "two steps → two snapshots");
        assert!(status.contains("faces"), "status: {status}");
        // The scrubber's intermediate history: step 1 is the bare box (6
        // faces), step 2 is the punched cube (more than 6).
        assert_eq!(history[0][0].faces(), 6, "first snapshot is the bare box");
        assert!(history[1][0].faces() > 6, "second snapshot is punched");
        // Every step tessellates to a non-empty viewport mesh.
        for k in 1..=history.len() {
            let mesh = tessellate_step(&history, k, &[]).expect("tessellate step");
            assert!(
                crate::mesh_loader::mesh_bounding_box(&mesh).is_some(),
                "step {k} should tessellate to a non-empty mesh"
            );
        }
    }

    #[test]
    fn cut_step_reduces_the_per_step_volume() {
        // The default tree is a box (step 1) with a cylinder cut (step 2). The
        // per-step volume the scrubber now shows must drop when the hole is
        // punched.
        let s = CadWorkbenchState::default();
        let (history, _) = rebuild_tree(&s).expect("default tree rebuilds");
        assert_eq!(history.len(), 2, "two steps → two snapshots");
        let v_box = total_volume(&history[0]);
        let v_punched = total_volume(&history[1]);
        assert!(v_box > 0.0, "the bare box has positive volume");
        assert!(
            v_punched < v_box,
            "cutting the hole reduces the volume: {v_box} → {v_punched}"
        );
    }

    #[test]
    fn revolve_step_rebuilds_to_a_nonempty_solid() {
        // A feature tree whose sole step is a New Revolve of a cone half-section
        // must rebuild to a non-empty body that tessellates to a real mesh —
        // proving the Revolve FeatureKind threads through the timeline + kernel.
        let mut s = CadWorkbenchState::default();
        let mut rev = UiStep::new_revolve();
        rev.op = Op::New; // start a fresh body from the revolve
        s.steps = vec![rev];
        let (history, status) = rebuild_tree(&s).expect("revolve tree rebuilds");
        assert_eq!(history.len(), 1, "one step → one snapshot");
        assert!(
            history[0][0].faces() > 0,
            "the revolved body has faces; status: {status}"
        );
        let mesh = tessellate_step(&history, 1, &[]).expect("tessellate revolve");
        assert!(
            crate::mesh_loader::mesh_bounding_box(&mesh).is_some(),
            "the revolved solid tessellates to a non-empty viewport mesh"
        );
        // The default half-section is a cone (r 1, h 2) → analytic volume
        // π r² h / 3 ≈ 2.094; the BRep volume converges from below.
        let vol = total_volume(&history[0]);
        assert!(
            vol > 1.5,
            "revolved cone volume {vol} is in the right range"
        );
    }

    #[test]
    fn add_extrude_from_sketch_rebuilds_a_triangular_prism() {
        // The draw-then-extrude path: a 3-point triangle drawn on the sketch
        // canvas + a height must yield a feature tree that rebuilds to a
        // non-empty solid carrying that exact profile.
        let mut s = CadWorkbenchState::default();
        s.steps.clear(); // start from an empty tree → the extrude is the base body
        let triangle = [[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]];
        let height = 3.0;
        s.add_extrude_from_sketch(&triangle, height);

        // One Extrude step, op New (first step starts the body), carrying the
        // drawn profile and the requested height.
        assert_eq!(s.steps.len(), 1, "one extrude step appended");
        let st = &s.steps[0];
        assert_eq!(st.kind, FeatureKind::Extrude);
        assert_eq!(st.op, Op::New, "an empty tree's first step must be New");
        assert_eq!(
            st.profile,
            vec![(0.0, 0.0), (2.0, 0.0), (0.0, 2.0)],
            "the step carries the drawn triangle profile"
        );
        assert_eq!(st.height, "3", "the step carries the requested height");

        // It rebuilds to a non-empty solid.
        let (history, status) = rebuild_tree(&s).expect("extruded sketch rebuilds");
        assert_eq!(history.len(), 1, "one step → one snapshot");
        assert!(
            history[0][0].faces() > 0,
            "the extruded prism has faces; status: {status}"
        );
        let mesh = tessellate_step(&history, 1, &[]).expect("tessellate extrude");
        assert!(
            crate::mesh_loader::mesh_bounding_box(&mesh).is_some(),
            "the extruded prism tessellates to a non-empty viewport mesh"
        );
        // Triangular prism volume = base area · height = (½·2·2) · 3 = 6 u³.
        let vol = total_volume(&history[0]);
        assert!(
            (vol - 6.0).abs() < 1e-6,
            "triangular-prism volume should be 6, got {vol}"
        );

        // On a non-empty tree the helper welds with Join instead of New.
        let mut t = CadWorkbenchState::default();
        let before = t.steps.len();
        t.add_extrude_from_sketch(&triangle, height);
        assert_eq!(t.steps.len(), before + 1, "appended onto the existing tree");
        assert_eq!(
            t.steps.last().unwrap().op,
            Op::Join,
            "a non-empty tree welds the extrude with Join"
        );

        // Degenerate inputs (fewer than 3 points, or a non-finite height) are
        // ignored rather than producing a broken step.
        let mut u = CadWorkbenchState::default();
        u.steps.clear();
        u.add_extrude_from_sketch(&[[0.0, 0.0], [1.0, 0.0]], 1.0);
        assert!(u.steps.is_empty(), "a 2-point profile is rejected");
        u.add_extrude_from_sketch(&triangle, f64::NAN);
        assert!(u.steps.is_empty(), "a non-finite height is rejected");
    }

    #[test]
    fn sketch_canvas_state_edits_build_and_close_a_polygon() {
        // Exercise the canvas's state mutators (the click handler / Undo / Clear
        // / close-loop logic call these): points accumulate, Undo pops + re-opens
        // the loop, Clear empties it, and a closed loop refuses further points.
        let mut s = CadWorkbenchState::default();
        assert!(s.sketch_points().is_empty());
        assert!(!s.sketch_is_closed());

        s.sketch_add_point(0.0, 0.0);
        s.sketch_add_point(1.0, 0.0);
        s.sketch_add_point(1.0, 1.0);
        assert_eq!(s.sketch_points().len(), 3);

        // Closing the loop (the canvas does this on a click near the first
        // vertex) then refuses further points.
        s.sketch_closed = true;
        s.sketch_add_point(2.0, 2.0);
        assert_eq!(
            s.sketch_points().len(),
            3,
            "a closed loop rejects new points"
        );

        // Undo pops the last vertex and re-opens the loop.
        s.sketch_undo();
        assert_eq!(s.sketch_points().len(), 2);
        assert!(!s.sketch_is_closed(), "undo re-opens a closed loop");

        // Clear empties the whole sketch.
        s.sketch_clear();
        assert!(s.sketch_points().is_empty());
        assert!(!s.sketch_is_closed());

        // focus_sketch raises the flag the panel consumes to scroll the canvas.
        assert!(!s.sketch_focus_request);
        s.focus_sketch();
        assert!(
            s.sketch_focus_request,
            "focus_sketch flags a scroll-to-canvas"
        );
    }

    #[test]
    fn circle_through_3_recovers_a_known_circle() {
        // Three points on the unit circle centred at (2, 1): the solver must
        // recover that centre + radius. Collinear points return None.
        let (centre, r) = circle_through_3([3.0, 1.0], [2.0, 2.0], [1.0, 1.0])
            .expect("non-collinear points define a circle");
        assert!((centre[0] - 2.0).abs() < 1e-9 && (centre[1] - 1.0).abs() < 1e-9);
        assert!((r - 1.0).abs() < 1e-9, "radius {r}");
        assert!(
            circle_through_3([0.0, 0.0], [1.0, 1.0], [2.0, 2.0]).is_none(),
            "collinear points have no finite circle"
        );
    }

    #[test]
    fn arc_segment_tessellates_onto_its_circle() {
        // A quarter-ish arc start (1,0) → end (0,1) through via (√½, √½): every
        // sampled polyline point must lie on the circle through the 3 points
        // (here the unit circle, centre origin, r = 1).
        let start = [1.0, 0.0];
        let via = [
            std::f64::consts::FRAC_1_SQRT_2,
            std::f64::consts::FRAC_1_SQRT_2,
        ];
        let end = [0.0, 1.0];
        let (centre, r) = circle_through_3(start, via, end).expect("circle");
        assert!(centre[0].abs() < 1e-9 && centre[1].abs() < 1e-9 && (r - 1.0).abs() < 1e-9);

        let poly = sample_arc(start, via, end);
        assert!(poly.len() >= 2, "arc sampled into a polyline");
        // First sample is exactly the start anchor.
        assert!((poly[0][0] - start[0]).abs() < 1e-9 && (poly[0][1] - start[1]).abs() < 1e-9);
        for p in &poly {
            let d = ((p[0] - centre[0]).powi(2) + (p[1] - centre[1]).powi(2)).sqrt();
            assert!(
                (d - r).abs() < 1e-6,
                "sample {p:?} is on the circle (d={d})"
            );
        }

        // The full profile through one Arc segment keeps the end anchor and every
        // interior point still rides the circle.
        let mut s = CadWorkbenchState::default();
        s.sketch_clear();
        s.sketch_add_arc(start, via, end);
        let profile = s.sketch_profile_polyline();
        assert!(
            (profile[0][0] - start[0]).abs() < 1e-9 && (profile[0][1] - start[1]).abs() < 1e-9,
            "profile starts at the arc start"
        );
        let last = *profile.last().unwrap();
        assert!(
            (last[0] - end[0]).abs() < 1e-9 && (last[1] - end[1]).abs() < 1e-9,
            "profile ends at the arc end anchor"
        );
        for p in &profile {
            let d = ((p[0] - centre[0]).powi(2) + (p[1] - centre[1]).powi(2)).sqrt();
            assert!(
                (d - r).abs() < 1e-6,
                "profile point {p:?} on circle (d={d})"
            );
        }
    }

    #[test]
    fn collinear_arc_falls_back_to_a_straight_chord() {
        // Three collinear "arc" points must not panic or NaN — the segment
        // degrades to the straight chord start→end (just the two endpoints).
        let mut s = CadWorkbenchState::default();
        s.sketch_clear();
        s.sketch_add_arc([0.0, 0.0], [1.0, 0.0], [2.0, 0.0]);
        let profile = s.sketch_profile_polyline();
        assert_eq!(
            profile,
            vec![[0.0, 0.0], [2.0, 0.0]],
            "a collinear arc is a plain chord"
        );
        assert!(
            profile.iter().all(|p| p[0].is_finite() && p[1].is_finite()),
            "no NaN/inf from the degenerate arc"
        );
    }

    #[test]
    fn catmull_rom_spline_passes_through_its_control_points() {
        // The tessellated curve must interpolate every control point (Catmull-Rom
        // is an interpolating spline). Sample an S-shaped 4-point curve and assert
        // each control point appears (within tol) somewhere on the polyline.
        let ctrl = [[0.0, 0.0], [1.0, 1.5], [2.0, -0.5], [3.0, 1.0]];
        let mut curve = sample_spline(&ctrl);
        // sample_spline omits the final endpoint (caller appends it); add it so we
        // can check the last control point too.
        curve.push(*ctrl.last().unwrap());
        for cp in &ctrl {
            let hit = curve
                .iter()
                .any(|p| (p[0] - cp[0]).abs() < 1e-6 && (p[1] - cp[1]).abs() < 1e-6);
            assert!(hit, "spline passes through control point {cp:?}");
        }
        // The first sample is exactly the first control point.
        assert!((curve[0][0]).abs() < 1e-9 && (curve[0][1]).abs() < 1e-9);

        // The same through the public segment path: a Spline segment's profile
        // hits the seed anchor and every run-through point.
        let mut s = CadWorkbenchState::default();
        s.sketch_clear();
        s.sketch_add_spline(&ctrl);
        let profile = s.sketch_profile_polyline();
        for cp in &ctrl {
            let hit = profile
                .iter()
                .any(|p| (p[0] - cp[0]).abs() < 1e-6 && (p[1] - cp[1]).abs() < 1e-6);
            assert!(hit, "profile spline passes through {cp:?}");
        }
    }

    #[test]
    fn arc_profile_extrudes_to_a_nonempty_solid() {
        // A closed profile containing an Arc must tessellate and extrude into a
        // real solid via add_extrude_from_sketch — the headline curved-section
        // capability. Build: line (0,0)→(2,0), arc up-and-over to (0,2) bulging
        // through (2.2, 1.0), then back to start.
        let mut s = CadWorkbenchState::default();
        s.sketch_clear();
        s.sketch_add_point(0.0, 0.0); // start anchor
        s.sketch_add_point(2.0, 0.0); // straight base edge
        s.sketch_add_arc([2.0, 0.0], [2.2, 1.0], [0.0, 2.0]); // bulging arc side

        let poly = s.sketch_profile_polyline();
        assert!(
            poly.len() > 3,
            "the arc inflates the profile past a bare triangle ({} pts)",
            poly.len()
        );
        assert!(
            poly.iter().all(|p| p[0].is_finite() && p[1].is_finite()),
            "tessellated profile is finite"
        );

        // Extrude the curved profile into the (empty) feature tree.
        let mut cad = CadWorkbenchState::default();
        cad.steps.clear();
        cad.add_extrude_from_sketch(&poly, 1.5);
        assert_eq!(
            cad.steps.len(),
            1,
            "one extrude step from the curved profile"
        );
        assert_eq!(cad.steps[0].kind, FeatureKind::Extrude);
        assert_eq!(
            cad.steps[0].profile.len(),
            poly.len(),
            "the step carries every tessellated profile point"
        );

        let (history, status) = rebuild_tree(&cad).expect("curved profile rebuilds");
        assert_eq!(history.len(), 1, "one step → one snapshot");
        assert!(
            history[0][0].faces() > 0,
            "the extruded curved solid has faces; status: {status}"
        );
        let mesh = tessellate_step(&history, 1, &[]).expect("tessellate curved extrude");
        assert!(
            crate::mesh_loader::mesh_bounding_box(&mesh).is_some(),
            "the curved-section solid tessellates to a non-empty viewport mesh"
        );
        let vol = total_volume(&history[0]);
        assert!(
            vol > 0.0,
            "curved-section solid has positive volume ({vol})"
        );
    }

    #[test]
    fn sketch_tool_switch_clears_pending_and_undo_drops_segments() {
        // The segment model: Line/Arc/Spline build an ordered segment chain;
        // switching tool abandons a half-entered segment; Undo drops whole
        // committed segments (and the lone start anchor last).
        let mut s = CadWorkbenchState::default();
        s.sketch_clear();
        assert_eq!(s.sketch_anchor_count(), 0);

        // Two line clicks → start + 1 Line segment → 2 anchors.
        s.sketch_add_point(0.0, 0.0);
        s.sketch_add_point(1.0, 0.0);
        assert_eq!(s.sketch_anchor_count(), 2);

        // An arc adds one segment (3 anchors total).
        s.sketch_add_arc([1.0, 0.0], [1.5, 0.5], [1.0, 1.0]);
        assert_eq!(s.sketch_anchor_count(), 3);
        assert_eq!(s.sketch_segs.len(), 2, "one Line + one Arc segment");

        // Undo drops the Arc, then the Line, then the start anchor.
        s.sketch_undo();
        assert_eq!(s.sketch_segs.len(), 1);
        s.sketch_undo();
        assert_eq!(s.sketch_segs.len(), 0);
        assert_eq!(s.sketch_anchor_count(), 1, "start anchor remains");
        s.sketch_undo();
        assert_eq!(
            s.sketch_anchor_count(),
            0,
            "undo drops the lone start anchor"
        );
        assert!(s.sketch_start.is_none());

        // A degenerate 1-point spline seed (no run-through points) adds no segment.
        s.sketch_add_spline(&[[0.5, 0.5]]);
        assert_eq!(
            s.sketch_anchor_count(),
            1,
            "seed-only spline is just the start"
        );
        assert!(
            s.sketch_segs.is_empty(),
            "no spline segment from a lone seed"
        );
    }

    #[test]
    fn fmt_num_trims_trailing_zeros() {
        assert_eq!(fmt_num(1.0), "1");
        assert_eq!(fmt_num(0.25), "0.25");
        assert_eq!(fmt_num(1.5), "1.5");
        assert_eq!(fmt_num(-2.0), "-2");
        assert_eq!(fmt_num(0.0), "0");
    }

    #[test]
    fn part_design_menu_mutators_share_the_panel_path() {
        // The top-bar Part Design menu drives these same `CadWorkbenchState`
        // methods the panel's + buttons call. Exercise each and assert the
        // feature tree + flags land as expected.
        let mut s = CadWorkbenchState::default();

        // Populate a real rebuild history first (via the shared rebuild path),
        // so we can prove the reset clears it rather than injecting a fake one.
        perform_rebuild(&mut s);
        assert!(s.history.is_some(), "a rebuild populates the history");

        // New part → a single base solid (a New box), cached rebuild cleared.
        s.reset_to_base_solid();
        assert_eq!(s.steps.len(), 1, "reset leaves exactly one base step");
        assert_eq!(s.steps[0].op, Op::New, "the base step starts a new body");
        assert_eq!(s.steps[0].kind, FeatureKind::Box);
        assert!(
            s.history.is_none(),
            "reset clears the stale rebuild history"
        );

        // Each add_* appends one step of the expected kind.
        s.add_extrude();
        assert_eq!(s.steps.last().unwrap().kind, FeatureKind::Extrude);
        s.add_revolve();
        assert_eq!(s.steps.last().unwrap().kind, FeatureKind::Revolve);
        s.add_box();
        assert_eq!(s.steps.last().unwrap().kind, FeatureKind::Box);
        s.add_cylinder();
        assert_eq!(s.steps.last().unwrap().kind, FeatureKind::Cylinder);
        s.add_sphere();
        assert_eq!(s.steps.last().unwrap().kind, FeatureKind::Sphere);
        s.add_cone();
        assert_eq!(s.steps.last().unwrap().kind, FeatureKind::Cone);
        s.add_torus();
        assert_eq!(s.steps.last().unwrap().kind, FeatureKind::Torus);
        assert_eq!(s.steps.len(), 8, "7 adds on top of the 1 base step");

        // Boolean submenu sets the last step's op.
        s.set_last_op_cut();
        assert_eq!(s.steps.last().unwrap().op, Op::Cut);
        s.set_last_op_intersect();
        assert_eq!(s.steps.last().unwrap().op, Op::Intersect);
        s.set_last_op_union();
        assert_eq!(s.steps.last().unwrap().op, Op::Join);

        // request_rebuild raises the flag draw_cad_workbench consumes.
        assert!(!s.rebuild_request);
        s.request_rebuild();
        assert!(
            s.rebuild_request,
            "request_rebuild flags a viewport rebuild"
        );

        // "New part (reset)" leaves a lone base box that rebuilds to a real
        // solid on its own (the menu's New-part → rebuild path).
        let mut t = CadWorkbenchState::default();
        t.reset_to_base_solid();
        let (history, _status) = rebuild_tree(&t).expect("reset base solid rebuilds");
        assert_eq!(history.len(), 1, "the reset part is a single base step");
        assert_eq!(history[0][0].faces(), 6, "the base solid is a box");
    }

    #[test]
    fn rebuild_reports_total_solid_volume() {
        // A single 2×3×4 box — flat-faced, so its measured volume is exact.
        let mut bx = UiStep::base(Op::New, FeatureKind::Box);
        bx.dx = "2".into();
        bx.dy = "3".into();
        bx.dz = "4".into();
        let s = CadWorkbenchState {
            steps: vec![bx],
            ..CadWorkbenchState::default()
        };
        let (history, status) = rebuild_tree(&s).expect("box rebuilds");
        let v = total_volume(history.last().expect("one snapshot"));
        assert!(
            (v - 24.0).abs() < 1e-6,
            "2×3×4 box volume should be 24, got {v}"
        );
        // The status line surfaces the volume in cubic model units.
        assert!(
            status.contains("u³"),
            "status should report volume: {status}"
        );

        // The default tree punches a cylinder through a unit (size = 1) box:
        // the cut removes material, so the result has positive volume strictly
        // below the 1×1×1 box it started from.
        let punched = CadWorkbenchState::default();
        let (ph, _) = rebuild_tree(&punched).expect("punched cube rebuilds");
        let punched_vol = total_volume(ph.last().expect("snapshot"));
        assert!(
            punched_vol > 0.0 && punched_vol < 1.0,
            "punched unit cube volume should be in (0, 1): {punched_vol}"
        );
    }

    #[test]
    fn rebuild_reports_total_surface_area() {
        // A single 2×3×4 box — flat-faced, so its measured area is exact:
        // 2·(2·3 + 3·4 + 2·4) = 2·26 = 52 u².
        let mut bx = UiStep::base(Op::New, FeatureKind::Box);
        bx.dx = "2".into();
        bx.dy = "3".into();
        bx.dz = "4".into();
        let s = CadWorkbenchState {
            steps: vec![bx],
            ..CadWorkbenchState::default()
        };
        let (history, status) = rebuild_tree(&s).expect("box rebuilds");
        let a = total_area(history.last().expect("one snapshot"));
        assert!(
            (a - 52.0).abs() < 1e-6,
            "2×3×4 box area should be 52, got {a}"
        );
        // The status line surfaces the area in square model units.
        assert!(status.contains("u²"), "status should report area: {status}");
    }

    #[test]
    fn rebuild_reports_bounding_box_dimensions() {
        // A 2×3×4 box — flat faces ⇒ exact tessellated extents.
        let mut bx = UiStep::base(Op::New, FeatureKind::Box);
        bx.dx = "2".into();
        bx.dy = "3".into();
        bx.dz = "4".into();
        let s = CadWorkbenchState {
            steps: vec![bx],
            ..CadWorkbenchState::default()
        };
        let (history, status) = rebuild_tree(&s).expect("box rebuilds");
        let mesh = tessellate_step(&history, history.len(), &[]).expect("tessellate");
        let dims = mesh_dimensions(&mesh).expect("non-empty mesh");
        assert!((dims[0] - 2.0).abs() < 1e-4, "dx {}", dims[0]);
        assert!((dims[1] - 3.0).abs() < 1e-4, "dy {}", dims[1]);
        assert!((dims[2] - 4.0).abs() < 1e-4, "dz {}", dims[2]);
        // The status line surfaces the bounding box.
        assert!(
            status.contains("bbox"),
            "status should report bbox: {status}"
        );
        // A box fills its bounding box exactly — fill fraction ≈ 1.0.
        let vol = total_volume(history.last().expect("snapshot"));
        let fill = fill_fraction(vol, dims).expect("fill");
        assert!(
            (fill - 1.0).abs() < 1e-3,
            "box should fill its bbox: {fill}"
        );
        assert!(
            status.contains("fill"),
            "status should report fill: {status}"
        );
    }

    #[test]
    fn fill_fraction_relates_volume_to_bounding_box() {
        // A box exactly fills its bounding box.
        assert!((fill_fraction(24.0, [2.0, 3.0, 4.0]).unwrap() - 1.0).abs() < 1e-12);
        // A unit-radius cylinder of height 2 (V = 2π) inside its 2×2×2 bounding
        // box fills exactly π/4.
        let cyl = fill_fraction(2.0 * std::f64::consts::PI, [2.0, 2.0, 2.0]).unwrap();
        assert!((cyl - std::f64::consts::FRAC_PI_4).abs() < 1e-12);
        // A degenerate (zero-volume) bounding box yields None, not a divide by zero.
        assert!(fill_fraction(1.0, [0.0, 1.0, 1.0]).is_none());
    }

    #[test]
    fn sphericity_is_one_for_a_sphere_and_known_for_a_cube() {
        let pi = std::f64::consts::PI;
        // The sphere is the maximally compact solid → ψ = 1 exactly.
        let r = 2.5;
        let v_sphere = 4.0 / 3.0 * pi * r * r * r;
        let a_sphere = 4.0 * pi * r * r;
        assert!((sphericity(v_sphere, a_sphere).unwrap() - 1.0).abs() < 1e-12);
        // A cube has the classic Wadell sphericity ≈ 0.8060 (scale-independent).
        let a = 2.0;
        assert!((sphericity(a * a * a, 6.0 * a * a).unwrap() - 0.8060).abs() < 1e-3);
        assert!((sphericity(1.0, 6.0).unwrap() - 0.8060).abs() < 1e-3);
        // Non-positive area or volume → None, never a divide blow-up.
        assert!(sphericity(1.0, 0.0).is_none());
        assert!(sphericity(0.0, 1.0).is_none());
    }

    #[test]
    fn bbox_diagonal_is_the_euclidean_norm_of_the_extents() {
        // A 3-4-12 box has space diagonal 13 (3²+4²+12² = 169).
        assert!((bbox_diagonal([3.0, 4.0, 12.0]) - 13.0).abs() < 1e-6);
        // A 1-2-2 box → 3.
        assert!((bbox_diagonal([1.0, 2.0, 2.0]) - 3.0).abs() < 1e-6);
        // A flat box collapses to its in-plane diagonal (3-4 → 5).
        assert!((bbox_diagonal([3.0, 4.0, 0.0]) - 5.0).abs() < 1e-6);
    }

    #[test]
    fn surface_area_to_volume_ratio_follows_the_square_cube_law() {
        let pi = std::f64::consts::PI;
        // A sphere of radius r: S/V = 4πr² / (4/3·πr³) = 3/r.
        let r = 2.0;
        let v = 4.0 / 3.0 * pi * r * r * r;
        let a = 4.0 * pi * r * r;
        assert!((surface_area_to_volume_ratio(a, v).unwrap() - 3.0 / r).abs() < 1e-12);
        // A cube of side s: S/V = 6s² / s³ = 6/s.
        let s = 4.0;
        let cube = surface_area_to_volume_ratio(6.0 * s * s, s * s * s).unwrap();
        assert!((cube - 6.0 / s).abs() < 1e-12, "cube S/V {cube}");
        // Doubling every length halves S/V (the square–cube law).
        let big =
            surface_area_to_volume_ratio(6.0 * (2.0 * s) * (2.0 * s), (2.0 * s).powi(3)).unwrap();
        assert!((big - cube / 2.0).abs() < 1e-12, "scaling {big}");
        // Non-positive volume → None, never a divide blow-up.
        assert!(surface_area_to_volume_ratio(1.0, 0.0).is_none());
    }

    #[test]
    fn bbox_aspect_ratio_is_longest_over_shortest_extent() {
        // A 3×4×12 box → 12/3 = 4.
        assert!((bbox_aspect_ratio([3.0, 4.0, 12.0]).unwrap() - 4.0).abs() < 1e-9);
        // A cube → 1.0.
        assert!((bbox_aspect_ratio([5.0, 5.0, 5.0]).unwrap() - 1.0).abs() < 1e-9);
        // Order-independent.
        assert!((bbox_aspect_ratio([12.0, 3.0, 4.0]).unwrap() - 4.0).abs() < 1e-9);
        // A degenerate flat box (zero shortest extent) → None.
        assert!(bbox_aspect_ratio([0.0, 4.0, 12.0]).is_none());
    }

    #[test]
    fn mesh_centroid_is_the_volume_centroid_not_the_vertex_average() {
        use nalgebra::Vector3;
        // A square pyramid: 2×2 base at z=0, apex at z=4. Its solid centroid sits
        // at h/4 = 1.0 up the axis — distinct from the 5-vertex average (z = 0.8),
        // so this discriminates the true volume centroid from a vertex mean.
        let mut mesh = valenx_mesh::Mesh::new("pyramid");
        mesh.nodes = vec![
            Vector3::new(-1.0, -1.0, 0.0),
            Vector3::new(1.0, -1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(-1.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 4.0),
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = vec![
            0, 2, 1, 0, 3, 2, // base (outward normal −z)
            0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4, // four sides
        ];
        mesh.element_blocks.push(block);
        let c = mesh_centroid(&mesh).expect("closed mesh has a centroid");
        assert!(c[0].abs() < 1e-9, "cx {}", c[0]);
        assert!(c[1].abs() < 1e-9, "cy {}", c[1]);
        assert!(
            (c[2] - 1.0).abs() < 1e-9,
            "cz {} (should be h/4 = 1.0)",
            c[2]
        );
        // An empty mesh has no centroid.
        assert!(mesh_centroid(&valenx_mesh::Mesh::new("empty")).is_none());
    }

    #[test]
    fn mesh_radius_of_gyration_matches_the_box_formula() {
        use nalgebra::Vector3;
        // A 1×2×3 box (corner at the origin). The radius of gyration about the
        // centroid is √((Lx²+Ly²+Lz²)/12), independent of where the box sits.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = vec![
            0, 2, 1, 0, 3, 2, // −z bottom
            4, 5, 6, 4, 6, 7, // +z top
            0, 1, 5, 0, 5, 4, // −y front
            3, 7, 6, 3, 6, 2, // +y back
            0, 4, 7, 0, 7, 3, // −x left
            1, 2, 6, 1, 6, 5, // +x right
        ];
        mesh.element_blocks.push(block);
        let k = mesh_radius_of_gyration(&mesh).expect("closed box");
        let expected = ((lx * lx + ly * ly + lz * lz) / 12.0_f64).sqrt();
        assert!((k - expected).abs() < 1e-9, "k {k} vs {expected}");
        // An empty mesh has no radius of gyration.
        assert!(mesh_radius_of_gyration(&valenx_mesh::Mesh::new("empty")).is_none());
    }

    #[test]
    fn mesh_volume_matches_the_box_volume() {
        use nalgebra::Vector3;
        // The canonical 1×2×3 box (corner at the origin). The divergence-theorem
        // enclosed volume must equal the elementary Lx·Ly·Lz = 6, reached by a
        // completely different route (signed-tetrahedron sum vs length×width×height).
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = vec![
            0, 2, 1, 0, 3, 2, // −z bottom
            4, 5, 6, 4, 6, 7, // +z top
            0, 1, 5, 0, 5, 4, // −y front
            3, 7, 6, 3, 6, 2, // +y back
            0, 4, 7, 0, 7, 3, // −x left
            1, 2, 6, 1, 6, 5, // +x right
        ];
        mesh.element_blocks.push(block);
        let v = mesh_volume(&mesh).expect("closed box has a volume");
        assert!((v - 6.0).abs() < 1e-9, "box volume Lx·Ly·Lz = 6, got {v}");
        // Winding-direction independence: reversing every triangle flips the signed
        // sum, but the returned magnitude is unchanged (validates the abs()).
        for tri in mesh.element_blocks[0].connectivity.chunks_exact_mut(3) {
            tri.swap(1, 2);
        }
        let v_rev = mesh_volume(&mesh).expect("reversed box");
        assert!(
            (v_rev - 6.0).abs() < 1e-9,
            "winding-independent magnitude: {v_rev}"
        );
        // Translation invariance: a closed surface's signed volume ignores position.
        for n in mesh.nodes.iter_mut() {
            *n += Vector3::new(10.0, -5.0, 7.0);
        }
        let v_shift = mesh_volume(&mesh).expect("translated box");
        assert!(
            (v_shift - 6.0).abs() < 1e-9,
            "translation-invariant: {v_shift}"
        );
        // Cubic scaling: doubling every coordinate scales the volume by 2³ = 8.
        for n in mesh.nodes.iter_mut() {
            *n *= 2.0;
        }
        let v_scaled = mesh_volume(&mesh).expect("scaled box");
        assert!(
            (v_scaled - 48.0).abs() < 1e-9,
            "V ∝ length³: 6·8 = 48, got {v_scaled}"
        );
        // An empty / triangle-free mesh has no enclosed volume.
        assert!(mesh_volume(&valenx_mesh::Mesh::new("empty")).is_none());
    }

    #[test]
    fn mesh_principal_moments_match_the_box_formula() {
        use nalgebra::Vector3;
        // The same 1×2×3 box. For a box the principal moments about its centre
        // are I = (m/12)(Lⱼ²+Lₖ²) with m = V = 6 (unit density):
        //   I_x = (6/12)(2²+3²) = 6.5, I_y = (6/12)(1²+3²) = 5.0, I_z = (6/12)(1²+2²) = 2.5.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let nodes = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let connectivity = vec![
            0, 2, 1, 0, 3, 2, // −z
            4, 5, 6, 4, 6, 7, // +z
            0, 1, 5, 0, 5, 4, // −y
            3, 7, 6, 3, 6, 2, // +y
            0, 4, 7, 0, 7, 3, // −x
            1, 2, 6, 1, 6, 5, // +x
        ];
        let make = |offset: Vector3<f64>| {
            let mut mesh = valenx_mesh::Mesh::new("box");
            mesh.nodes = nodes.iter().map(|n| *n + offset).collect();
            let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
            block.connectivity = connectivity.clone();
            mesh.element_blocks.push(block);
            mesh
        };

        let [i1, i2, i3] = mesh_principal_moments(&make(Vector3::zeros())).expect("closed box");
        let m = lx * ly * lz; // unit-density "mass" = volume = 6
        assert!(
            (i1 - m / 12.0 * (ly * ly + lz * lz)).abs() < 1e-9,
            "I1 {i1}"
        ); // 6.5
        assert!(
            (i2 - m / 12.0 * (lx * lx + lz * lz)).abs() < 1e-9,
            "I2 {i2}"
        ); // 5.0
        assert!(
            (i3 - m / 12.0 * (lx * lx + ly * ly)).abs() < 1e-9,
            "I3 {i3}"
        ); // 2.5

        // Principal moments are centroid-relative, so translating the box leaves
        // them unchanged — this validates the parallel-axis shift.
        let shifted =
            mesh_principal_moments(&make(Vector3::new(10.0, -5.0, 7.0))).expect("shifted box");
        assert!((shifted[0] - i1).abs() < 1e-6, "translation invariance I1");
        assert!((shifted[1] - i2).abs() < 1e-6, "translation invariance I2");
        assert!((shifted[2] - i3).abs() < 1e-6, "translation invariance I3");

        // An empty mesh has no inertia tensor.
        assert!(mesh_principal_moments(&valenx_mesh::Mesh::new("empty")).is_none());
    }

    #[test]
    fn mesh_is_watertight_distinguishes_a_closed_box_from_an_open_one() {
        use nalgebra::Vector3;
        // The same closed 1×2×3 box: 12 triangles, every edge shared by exactly two.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let closed = vec![
            0, 2, 1, 0, 3, 2, // −z bottom
            4, 5, 6, 4, 6, 7, // +z top
            0, 1, 5, 0, 5, 4, // −y front
            3, 7, 6, 3, 6, 2, // +y back
            0, 4, 7, 0, 7, 3, // −x left
            1, 2, 6, 1, 6, 5, // +x right
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = closed.clone();
        mesh.element_blocks.push(block);
        assert!(mesh_is_watertight(&mesh), "a closed box is watertight");

        // Drop one triangle (the last) → its three edges now belong to a single
        // triangle each, so the surface has a hole and is no longer watertight.
        let mut open = mesh.clone();
        open.element_blocks[0]
            .connectivity
            .truncate(closed.len() - 3);
        assert!(
            !mesh_is_watertight(&open),
            "a box missing a face triangle is open"
        );

        // An empty mesh (no triangles) is not watertight.
        assert!(!mesh_is_watertight(&valenx_mesh::Mesh::new("empty")));
    }

    #[test]
    fn mesh_max_aspect_ratio_finds_the_most_stretched_triangle() {
        use nalgebra::Vector3;
        // The canonical 1×2×3 box: each rectangular face splits into two right
        // triangles with legs (a,b) and hypotenuse √(a²+b²), so that triangle's
        // aspect ratio is √(a²+b²)/min(a,b). The most stretched faces are the 1×3
        // pair → √10/1 = √10 ≈ 3.162 (vs √5 on the 1×2 faces, √13/2 on the 2×3).
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let connectivity = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = connectivity;
        mesh.element_blocks.push(block);
        let ar = mesh_max_aspect_ratio(&mesh).expect("box has triangles");
        assert!(
            (ar - 10.0_f64.sqrt()).abs() < 1e-9,
            "box max aspect ratio = √10, got {ar}"
        );

        // A single equilateral triangle is the ideal: AR = 1 exactly.
        let mut eq = valenx_mesh::Mesh::new("equilateral");
        eq.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, 3.0_f64.sqrt() / 2.0, 0.0),
        ];
        let mut eb = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        eb.connectivity = vec![0, 1, 2];
        eq.element_blocks.push(eb);
        assert!(
            (mesh_max_aspect_ratio(&eq).unwrap() - 1.0).abs() < 1e-9,
            "equilateral → 1"
        );

        // A 3-4-5 right triangle → longest/shortest = 5/3.
        let mut r345 = valenx_mesh::Mesh::new("3-4-5");
        r345.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(0.0, 4.0, 0.0),
        ];
        let mut rb = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        rb.connectivity = vec![0, 1, 2];
        r345.element_blocks.push(rb);
        assert!(
            (mesh_max_aspect_ratio(&r345).unwrap() - 5.0 / 3.0).abs() < 1e-9,
            "3-4-5 → 5/3"
        );

        // No triangles → None.
        assert_eq!(
            mesh_max_aspect_ratio(&valenx_mesh::Mesh::new("empty")),
            None
        );
    }

    #[test]
    fn mesh_max_dihedral_angle_finds_the_sharpest_fold() {
        use nalgebra::Vector3;
        use std::f64::consts::FRAC_PI_2;
        // The canonical 1×2×3 box: every one of the 12 cube-rim edges folds at a
        // right angle (adjacent faces' outward normals are perpendicular), while
        // the 6 in-face diagonals are flat (coplanar). So the sharpest fold is π/2.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let connectivity = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = connectivity;
        mesh.element_blocks.push(block);
        let theta = mesh_max_dihedral_angle(&mesh).expect("box has interior edges");
        assert!(
            (theta - FRAC_PI_2).abs() < 1e-9,
            "box sharpest fold = π/2, got {theta}"
        );

        // Two coplanar triangles sharing an edge: the surface is flat → 0 fold.
        let mut flat = valenx_mesh::Mesh::new("flat");
        flat.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut fb = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        fb.connectivity = vec![0, 1, 2, 0, 2, 3];
        flat.element_blocks.push(fb);
        assert!(
            mesh_max_dihedral_angle(&flat).unwrap().abs() < 1e-9,
            "a coplanar patch has zero fold"
        );

        // Two flaps meeting at a right angle along their shared edge → π/2.
        let mut fold = valenx_mesh::Mesh::new("fold");
        fold.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0), // flap in the z=0 plane
            Vector3::new(0.0, 0.0, 1.0), // flap in the y=0 plane
        ];
        let mut gb = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        gb.connectivity = vec![0, 1, 2, 0, 1, 3];
        fold.element_blocks.push(gb);
        assert!(
            (mesh_max_dihedral_angle(&fold).unwrap() - FRAC_PI_2).abs() < 1e-9,
            "a right-angle fold = π/2"
        );

        // No interior manifold edge → None.
        assert_eq!(
            mesh_max_dihedral_angle(&valenx_mesh::Mesh::new("empty")),
            None
        );
    }

    #[test]
    fn mesh_mean_edge_length_is_the_characteristic_element_size() {
        use nalgebra::Vector3;
        // A single 3-4-5 right triangle: distinct edges {3, 4, 5}, mean = 12/3 = 4.
        let mut tri = valenx_mesh::Mesh::new("3-4-5");
        tri.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(0.0, 4.0, 0.0),
        ];
        let mut tb = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        tb.connectivity = vec![0, 1, 2];
        tri.element_blocks.push(tb);
        assert!(
            (mesh_mean_edge_length(&tri).unwrap() - 4.0).abs() < 1e-9,
            "3-4-5 mean edge = 4"
        );

        // The 1×2×3 box: 18 distinct edges — 12 cube (4×1 + 4×2 + 4×3 = 24) and 6
        // face diagonals (2√5 + 2√10 + 2√13) — so the mean is
        // (24 + 2√5 + 2√10 + 2√13) / 18 ≈ 2.334.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let connectivity = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = connectivity;
        mesh.element_blocks.push(block);
        let expected =
            (24.0 + 2.0 * 5.0_f64.sqrt() + 2.0 * 10.0_f64.sqrt() + 2.0 * 13.0_f64.sqrt()) / 18.0;
        let mean = mesh_mean_edge_length(&mesh).expect("box has edges");
        assert!(
            (mean - expected).abs() < 1e-9,
            "box mean edge ≈ 2.334, got {mean}"
        );

        // No edges → None.
        assert_eq!(
            mesh_mean_edge_length(&valenx_mesh::Mesh::new("empty")),
            None
        );
    }

    #[test]
    fn mesh_total_edge_length_sums_the_box_wireframe() {
        use nalgebra::Vector3;
        // A single 3-4-5 right triangle: distinct edges {3, 4, 5} → wireframe 12.
        let mut tri = valenx_mesh::Mesh::new("3-4-5");
        tri.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(0.0, 4.0, 0.0),
        ];
        let mut tb = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        tb.connectivity = vec![0, 1, 2];
        tri.element_blocks.push(tb);
        assert!(
            (mesh_total_edge_length(&tri).unwrap() - 12.0).abs() < 1e-9,
            "3-4-5 wireframe = 12"
        );

        // The canonical 1×2×3 box: 18 distinct edges — 12 cuboid (4×1 + 4×2 + 4×3 =
        // 24) and 6 face diagonals (2√5 + 2√10 + 2√13) → total ≈ 42.008.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        mesh.element_blocks.push(block);
        let expected = 24.0 + 2.0 * 5.0_f64.sqrt() + 2.0 * 10.0_f64.sqrt() + 2.0 * 13.0_f64.sqrt();
        let total = mesh_total_edge_length(&mesh).expect("box has edges");
        assert!(
            (total - expected).abs() < 1e-9,
            "box wireframe ≈ 42.008, got {total}"
        );
        // STRONG cross-check: total = mean × edge_count, threading mesh_mean_edge_length
        // (the box has 18 distinct edges, counted independently here).
        let mean = mesh_mean_edge_length(&mesh).expect("box has edges");
        assert!(
            (total - mean * 18.0).abs() < 1e-9,
            "total = mean × 18 distinct edges"
        );
        // No edges → None.
        assert_eq!(
            mesh_total_edge_length(&valenx_mesh::Mesh::new("empty")),
            None
        );
    }

    #[test]
    fn mesh_max_vertex_valence_counts_incident_edges() {
        use nalgebra::Vector3;
        use std::collections::{HashMap, HashSet};
        // The canonical 1×2×3 box: 8 corners, 12 Tri3 faces, 18 distinct edges.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let connectivity = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = connectivity.clone();
        mesh.element_blocks.push(block);

        // In this triangulation the two diagonal "pole" corners (vertices 0 and 6)
        // each touch 6 of the other 7 corners; the remaining six corners have
        // valence 4. So the maximum vertex valence is 6 — a hand-derived fact.
        assert_eq!(
            mesh_max_vertex_valence(&mesh),
            Some(6),
            "box max valence = 6"
        );

        // STRONG non-tautological cross-check via the HANDSHAKING LEMMA: build the
        // distinct-edge set independently (the global edge-dedup path, as in
        // mesh_mean_edge_length) and a per-vertex degree map, then assert
        // Σ(degree) == 2·|edges|. Two different aggregations of the same topology.
        let mut edges: HashSet<(u32, u32)> = HashSet::new();
        let mut nbr: HashMap<u32, HashSet<u32>> = HashMap::new();
        for tri in connectivity.chunks_exact(3) {
            for &(i, j) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                edges.insert(if i <= j { (i, j) } else { (j, i) });
                nbr.entry(i).or_default().insert(j);
                nbr.entry(j).or_default().insert(i);
            }
        }
        let sum_deg: usize = nbr.values().map(HashSet::len).sum();
        assert_eq!(edges.len(), 18, "box has 18 distinct edges");
        assert_eq!(
            sum_deg,
            2 * edges.len(),
            "handshaking lemma: Σdeg = 2E = 36"
        );
        // Mean valence = 2E/V = 36/8 = 4.5 (the average matches; the max exceeds it).
        assert!(
            (sum_deg as f64 / mesh.nodes.len() as f64 - 4.5).abs() < 1e-12,
            "mean valence 4.5"
        );

        // No Tri3 faces → None.
        assert_eq!(
            mesh_max_vertex_valence(&valenx_mesh::Mesh::new("empty")),
            None
        );
    }

    #[test]
    fn mesh_mean_vertex_valence_averages_the_degree() {
        use nalgebra::Vector3;
        use std::collections::HashSet;
        // The canonical 1×2×3 box: 8 corners, 12 Tri3 faces, 18 distinct edges.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let connectivity = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = connectivity.clone();
        mesh.element_blocks.push(block);

        // Mean valence = 2E/V = 2·18/8 = 4.5 exactly.
        let mean = mesh_mean_vertex_valence(&mesh).expect("box has vertices");
        assert!(
            (mean - 4.5).abs() < 1e-12,
            "box mean valence = 4.5, got {mean}"
        );
        // NON-TAUTOLOGICAL handshaking-lemma cross-check: independently count the
        // distinct edges and vertices; the mean degree is 2·|edges|/|V|. (The helper
        // sums per-vertex degrees; this uses the global edge set — different paths.)
        let mut edges: HashSet<(u32, u32)> = HashSet::new();
        let mut verts: HashSet<u32> = HashSet::new();
        for tri in connectivity.chunks_exact(3) {
            for &(i, j) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                edges.insert(if i <= j { (i, j) } else { (j, i) });
                verts.insert(i);
                verts.insert(j);
            }
        }
        assert_eq!(
            (edges.len(), verts.len()),
            (18, 8),
            "box has 18 edges, 8 vertices"
        );
        assert!(
            (mean - 2.0 * edges.len() as f64 / verts.len() as f64).abs() < 1e-12,
            "mean = 2E/V"
        );
        // The mean never exceeds the worst-case max (#201): 4.5 ≤ 6.
        let max = mesh_max_vertex_valence(&mesh).expect("box has vertices") as f64;
        assert!(mean <= max, "mean {mean} ≤ max {max}");
        // No Tri3 faces → None.
        assert_eq!(
            mesh_mean_vertex_valence(&valenx_mesh::Mesh::new("empty")),
            None
        );
    }

    #[test]
    fn mesh_min_triangle_area_finds_the_smallest_facet() {
        use nalgebra::Vector3;
        // The canonical 1×2×3 box: 12 Tri3 facets — four of area 1 (the 1×2 faces),
        // four of area 1.5 (the 1×3 faces), four of area 3 (the 2×3 faces).
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let connectivity = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = connectivity.clone();
        mesh.element_blocks.push(block);

        // The smallest facet is one of the 1×2-face half-triangles, area 1.0.
        let min = mesh_min_triangle_area(&mesh).expect("box has facets");
        assert!(
            (min - 1.0).abs() < 1e-9,
            "box min facet area = 1.0, got {min}"
        );

        // NON-TAUTOLOGICAL cross-checks: independently recompute all 12 facet areas;
        // the returned min must be ≤ every one, AND the areas must SUM to the box's
        // closed-form surface area 2·(1·2 + 1·3 + 2·3) = 22 (an independent total
        // tying the facet areas to the known geometry).
        let areas: Vec<f64> = connectivity
            .chunks_exact(3)
            .map(|t| {
                let a = mesh.nodes[t[0] as usize];
                let b = mesh.nodes[t[1] as usize];
                let c = mesh.nodes[t[2] as usize];
                0.5 * (b - a).cross(&(c - a)).norm()
            })
            .collect();
        assert!(
            areas.iter().all(|&a| min <= a + 1e-12),
            "min ≤ every facet area"
        );
        let total: f64 = areas.iter().sum();
        assert!(
            (total - 22.0).abs() < 1e-9,
            "Σ facet areas = box surface area 22, got {total}"
        );

        // No Tri3 faces → None.
        assert_eq!(
            mesh_min_triangle_area(&valenx_mesh::Mesh::new("empty")),
            None
        );
    }

    #[test]
    fn mesh_max_triangle_area_finds_the_largest_facet() {
        use nalgebra::Vector3;
        // The canonical 1×2×3 box: 12 Tri3 facets of areas {1×4, 1.5×4, 3×4}.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let connectivity = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = connectivity.clone();
        mesh.element_blocks.push(block);

        // The largest facet is one of the 2×3-face half-triangles, area 3.0.
        let max = mesh_max_triangle_area(&mesh).expect("box has facets");
        assert!(
            (max - 3.0).abs() < 1e-9,
            "box max facet area = 3.0, got {max}"
        );

        // NON-TAUTOLOGICAL cross-checks: independently recompute all 12 facet areas;
        // the returned max must be ≥ every one, it brackets the min #207 (3 ≥ 1), and
        // the areas sum to the box surface area 22.
        let areas: Vec<f64> = connectivity
            .chunks_exact(3)
            .map(|t| {
                let a = mesh.nodes[t[0] as usize];
                let b = mesh.nodes[t[1] as usize];
                let c = mesh.nodes[t[2] as usize];
                0.5 * (b - a).cross(&(c - a)).norm()
            })
            .collect();
        assert!(
            areas.iter().all(|&a| max >= a - 1e-12),
            "max ≥ every facet area"
        );
        let min = mesh_min_triangle_area(&mesh).expect("box has facets");
        assert!(
            max >= min && (min - 1.0).abs() < 1e-9,
            "extremes bracket: min 1 ≤ max 3"
        );
        let total: f64 = areas.iter().sum();
        assert!(
            (total - 22.0).abs() < 1e-9,
            "Σ facet areas = 22, got {total}"
        );

        // No Tri3 faces → None.
        assert_eq!(
            mesh_max_triangle_area(&valenx_mesh::Mesh::new("empty")),
            None
        );
    }

    #[test]
    fn mesh_surface_area_sums_the_box_facets() {
        use nalgebra::Vector3;
        // The canonical 1×2×3 box: 12 Tri3 facets of areas {1×4, 1.5×4, 3×4} → Σ = 22.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let connectivity = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = connectivity.clone();
        mesh.element_blocks.push(block);

        // The 1×2×3 box has surface area 2(1·2 + 1·3 + 2·3) = 22 exactly.
        let area = mesh_surface_area(&mesh).expect("box has facets");
        assert!(
            (area - 22.0).abs() < 1e-9,
            "box surface area = 22, got {area}"
        );

        // STRONG non-tautological cross-checks: independently re-sum all 12 facet
        // areas == the returned total; it is bracketed by n·min #207 ≤ area ≤ n·max
        // #225 (n = 12 facets) and bounded below by the single largest facet.
        let areas: Vec<f64> = connectivity
            .chunks_exact(3)
            .map(|t| {
                let a = mesh.nodes[t[0] as usize];
                let b = mesh.nodes[t[1] as usize];
                let c = mesh.nodes[t[2] as usize];
                0.5 * (b - a).cross(&(c - a)).norm()
            })
            .collect();
        let resum: f64 = areas.iter().sum();
        assert!(
            (area - resum).abs() < 1e-12,
            "area == Σ facet areas: {area} vs {resum}"
        );
        let n = areas.len() as f64;
        let min = mesh_min_triangle_area(&mesh).expect("box has facets");
        let max = mesh_max_triangle_area(&mesh).expect("box has facets");
        assert!(
            n * min <= area + 1e-12 && area <= n * max + 1e-12,
            "n·min ≤ area ≤ n·max"
        );
        assert!(area >= max, "a single facet cannot exceed the total");

        // No Tri3 faces → None.
        assert_eq!(mesh_surface_area(&valenx_mesh::Mesh::new("empty")), None);
    }

    #[test]
    fn mesh_open_edge_count_measures_the_hole_perimeter() {
        use nalgebra::Vector3;
        // The canonical closed 1×2×3 box: every edge is shared by two triangles.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let closed = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = closed.clone();
        mesh.element_blocks.push(block);
        assert_eq!(
            mesh_open_edge_count(&mesh),
            0,
            "a closed box has no open edges"
        );

        // Drop one triangle → its 3 edges fall from multiplicity 2 to 1.
        let mut open = mesh.clone();
        open.element_blocks[0]
            .connectivity
            .truncate(closed.len() - 3);
        assert_eq!(
            mesh_open_edge_count(&open),
            3,
            "a box missing one triangle has 3 open edges"
        );

        // A lone triangle: all three edges are boundary.
        let mut tri = valenx_mesh::Mesh::new("tri");
        tri.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
        ];
        let mut tblock = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        tblock.connectivity = vec![0, 1, 2];
        tri.element_blocks.push(tblock);
        assert_eq!(
            mesh_open_edge_count(&tri),
            3,
            "a lone triangle has a 3-edge boundary"
        );

        // Two coplanar triangles sharing edge 1–2: 5 edges, the shared one closed → 4 open.
        let mut sheet = tri.clone();
        sheet.element_blocks[0].connectivity = vec![0, 1, 2, 1, 3, 2];
        assert_eq!(
            mesh_open_edge_count(&sheet),
            4,
            "a 2-triangle sheet has a 4-edge boundary"
        );

        // An empty mesh has no edges.
        assert_eq!(mesh_open_edge_count(&valenx_mesh::Mesh::new("empty")), 0);
    }

    #[test]
    fn mesh_nonmanifold_edge_count_finds_self_intersections() {
        use nalgebra::Vector3;
        // The canonical closed 1×2×3 box: a clean 2-manifold — every edge shared by
        // exactly two triangles, so zero non-manifold edges.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        mesh.element_blocks.push(block);
        assert_eq!(
            mesh_nonmanifold_edge_count(&mesh),
            0,
            "a clean box is 2-manifold"
        );

        // A "fan" of three triangles all sharing the edge 0–1 (a T-junction): edge
        // 0–1 has multiplicity 3 → exactly one non-manifold edge; the other six edges
        // each belong to a single triangle (boundary, multiplicity 1).
        let mut fan = valenx_mesh::Mesh::new("fan");
        fan.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, -1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let mut fblock = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        fblock.connectivity = vec![0, 1, 2, 0, 1, 3, 0, 1, 4];
        fan.element_blocks.push(fblock);
        assert_eq!(
            mesh_nonmanifold_edge_count(&fan),
            1,
            "edge 0–1 is shared by 3 triangles"
        );
        // A genuinely different defect than a hole: the over-shared edge is invisible
        // to the open-edge (under-shared) count, yet both break watertightness.
        assert!(mesh_open_edge_count(&fan) > 0 && !mesh_is_watertight(&fan));

        // An empty mesh reports zero.
        assert_eq!(
            mesh_nonmanifold_edge_count(&valenx_mesh::Mesh::new("empty")),
            0
        );
    }

    #[test]
    fn mesh_boundary_loop_count_counts_the_holes() {
        use nalgebra::Vector3;
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let closed: Vec<u32> = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        let make = |conn: Vec<u32>| {
            let mut m = valenx_mesh::Mesh::new("box");
            m.nodes = nodes.clone();
            let mut b = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
            b.connectivity = conn;
            m.element_blocks.push(b);
            m
        };
        // The closed (watertight) box has no holes.
        assert_eq!(
            mesh_boundary_loop_count(&make(closed.clone())),
            0,
            "closed box → 0 holes"
        );

        // Drop the last triangle → one triangular gap (3 open edges, but ONE loop).
        let mut one = closed.clone();
        one.truncate(closed.len() - 3);
        assert_eq!(
            mesh_boundary_loop_count(&make(one)),
            1,
            "one missing triangle → 1 hole"
        );

        // Omit two vertex-disjoint triangles — T0=[0,2,1] (verts {0,1,2}) and
        // T3=[4,6,7] (verts {4,6,7}) — leaving two separate triangular holes.
        // (mesh_open_edge_count would report 6 open edges here; the loop count is 2.)
        let two: Vec<u32> = vec![
            0, 3, 2, 4, 5, 6, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6,
            5,
        ];
        assert_eq!(
            mesh_boundary_loop_count(&make(two)),
            2,
            "two disjoint gaps → 2 holes"
        );

        // An empty mesh has no boundary loops.
        assert_eq!(
            mesh_boundary_loop_count(&valenx_mesh::Mesh::new("empty")),
            0
        );
    }

    #[test]
    fn mesh_euler_characteristic_distinguishes_closed_from_open_topology() {
        use nalgebra::Vector3;
        // The closed 1×2×3 box: 8 vertices, 18 edges, 12 triangles →
        // χ = 8 − 18 + 12 = 2 (genus 0, the topology of a sphere).
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let closed: Vec<u32> = vec![
            0, 2, 1, 0, 3, 2, // −z
            4, 5, 6, 4, 6, 7, // +z
            0, 1, 5, 0, 5, 4, // −y
            3, 7, 6, 3, 6, 2, // +y
            0, 4, 7, 0, 7, 3, // −x
            1, 2, 6, 1, 6, 5, // +x
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = closed.clone();
        mesh.element_blocks.push(block);
        assert_eq!(
            mesh_euler_characteristic(&mesh),
            Some(2),
            "closed box: χ = 2 (genus 0)"
        );

        // Drop one triangle → an open surface: its vertices and edges survive on
        // the neighbouring faces, so only F drops by 1 → χ = 8 − 18 + 11 = 1.
        let mut open = mesh.clone();
        open.element_blocks[0]
            .connectivity
            .truncate(closed.len() - 3);
        assert_eq!(
            mesh_euler_characteristic(&open),
            Some(1),
            "box minus a face: χ = 1"
        );

        // Two disjoint boxes → χ is additive over components: 2 + 2 = 4.
        let mut two = mesh.clone();
        let n = two.nodes.len() as u32;
        let shifted: Vec<Vector3<f64>> = mesh
            .nodes
            .iter()
            .map(|p| *p + Vector3::new(10.0, 0.0, 0.0))
            .collect();
        two.nodes.extend(shifted);
        let mut block2 = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block2.connectivity = closed.iter().map(|&i| i + n).collect();
        two.element_blocks.push(block2);
        assert_eq!(
            mesh_euler_characteristic(&two),
            Some(4),
            "two disjoint boxes: χ = 4"
        );

        // An empty mesh has no Euler characteristic.
        assert_eq!(
            mesh_euler_characteristic(&valenx_mesh::Mesh::new("empty")),
            None
        );
    }

    #[test]
    fn brep_euler_characteristic_sums_the_solid_topology() {
        use nalgebra::Vector3;
        // A box B-rep has 6 faces, 12 edges, 8 vertices → χ = 8 − 12 + 6 = 2, the
        // genus-0 (sphere-topology) invariant read straight from the CAD topology.
        let cube = valenx_cad::box_solid(1.0, 2.0, 3.0).expect("box builds");
        assert_eq!(cube.faces(), 6);
        assert_eq!(cube.edges(), 12);
        assert_eq!(cube.vertices(), 8);
        let brep_chi = brep_euler_characteristic(std::slice::from_ref(&cube));
        assert_eq!(brep_chi, 2, "box B-rep χ = 8 − 12 + 6 = 2");

        // STRONG cross-check against an INDEPENDENT topology computation: the
        // tessellation Euler characteristic of a welded triangle mesh of the same
        // box. B-rep V−E+F (CAD topology) and triangulation V−E+F must agree — χ is a
        // topological invariant, blind to how the boundary is discretised.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut bm = valenx_mesh::Mesh::new("box");
        bm.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        bm.element_blocks.push(block);
        let mesh_chi = mesh_euler_characteristic(&bm).expect("box mesh has triangles");
        assert_eq!(mesh_chi, 2);
        assert_eq!(
            brep_chi, mesh_chi,
            "B-rep χ must equal the tessellation χ (topological invariance)"
        );

        // Additivity over disjoint bodies: two boxes → 2 + 2 = 4.
        let two = [
            valenx_cad::box_solid(1.0, 2.0, 3.0).unwrap(),
            valenx_cad::box_solid(2.0, 2.0, 2.0).unwrap(),
        ];
        assert_eq!(
            brep_euler_characteristic(&two),
            4,
            "χ is additive over bodies"
        );

        // No bodies → 0.
        assert_eq!(brep_euler_characteristic(&[]), 0);
    }

    #[test]
    fn mesh_total_angle_defect_satisfies_gauss_bonnet() {
        use nalgebra::Vector3;
        use std::f64::consts::PI;
        // The canonical 1×2×3 box: 8 corners, each a 3-right-angle corner, so each
        // angle defect is 2π − 3·(π/2) = π/2, summing to 4π = 2π·χ (χ = 2, genus 0).
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        mesh.element_blocks.push(block);

        let defect = mesh_total_angle_defect(&mesh).expect("box has triangles");
        assert!(
            (defect - 4.0 * PI).abs() < 1e-9,
            "box ∑defect = 4π, got {defect}"
        );
        // STRONG cross-check (Gauss–Bonnet): ∑defect == 2π·χ, threading the
        // combinatorial mesh_euler_characteristic (geometric vs topological paths).
        let chi = mesh_euler_characteristic(&mesh).expect("box has triangles");
        assert!(
            (defect - 2.0 * PI * chi as f64).abs() < 1e-9,
            "Gauss–Bonnet ∑defect = 2πχ"
        );
        // STRONG count identity: ∑defect == 2π·V − π·F (V=8 distinct verts, F=12).
        assert!(
            (defect - (2.0 * PI * 8.0 - PI * 12.0)).abs() < 1e-9,
            "∑defect = 2πV − πF"
        );
        // Scale-invariant — built from angles, so a ×10 zoom leaves it unchanged.
        for n in mesh.nodes.iter_mut() {
            *n *= 10.0;
        }
        let scaled = mesh_total_angle_defect(&mesh).expect("scaled box");
        assert!(
            (scaled - 4.0 * PI).abs() < 1e-9,
            "angle defect is scale-invariant: {scaled}"
        );
        // No triangles → None.
        assert_eq!(
            mesh_total_angle_defect(&valenx_mesh::Mesh::new("empty")),
            None
        );
    }

    #[test]
    fn mesh_shell_count_counts_disjoint_solids() {
        use nalgebra::Vector3;
        // One closed box → a single connected shell.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let closed: Vec<u32> = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = closed.clone();
        mesh.element_blocks.push(block);
        assert_eq!(mesh_shell_count(&mesh), 1, "one box is a single shell");

        // Two disjoint boxes (a second block on fresh node indices) → two shells.
        let mut two = mesh.clone();
        let n = two.nodes.len() as u32;
        let shifted: Vec<Vector3<f64>> = mesh
            .nodes
            .iter()
            .map(|p| *p + Vector3::new(10.0, 0.0, 0.0))
            .collect();
        two.nodes.extend(shifted);
        let mut block2 = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block2.connectivity = closed.iter().map(|&i| i + n).collect();
        two.element_blocks.push(block2);
        assert_eq!(
            mesh_shell_count(&two),
            2,
            "two disjoint boxes are two shells"
        );

        // An empty mesh has no shells.
        assert_eq!(mesh_shell_count(&valenx_mesh::Mesh::new("empty")), 0);
    }

    #[test]
    fn mesh_bounding_sphere_radius_is_the_half_space_diagonal_of_a_box() {
        use nalgebra::Vector3;
        // The minimal enclosing sphere of a box's 8 corners is centred at the box
        // centre with radius = half the space diagonal ½√(Lx²+Ly²+Lz²).
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let corners = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let mut mesh = valenx_mesh::Mesh::new("box");
        mesh.nodes = corners.to_vec();
        let r = mesh_bounding_sphere_radius(&mesh).expect("box has vertices");
        let expected = 0.5 * (lx * lx + ly * ly + lz * lz).sqrt();
        assert!((r - expected).abs() < 1e-9, "radius {r} vs {expected}");

        // Translating the part leaves the enclosing radius unchanged.
        let mut shifted = mesh.clone();
        for n in &mut shifted.nodes {
            *n += Vector3::new(10.0, -5.0, 7.0);
        }
        assert!(
            (mesh_bounding_sphere_radius(&shifted).unwrap() - expected).abs() < 1e-9,
            "translation-invariant"
        );

        // A degenerate or empty point set: one vertex → radius 0; no vertices → None.
        let mut one = valenx_mesh::Mesh::new("one");
        one.nodes = vec![Vector3::new(3.0, 4.0, 5.0)];
        assert_eq!(mesh_bounding_sphere_radius(&one), Some(0.0));
        assert!(mesh_bounding_sphere_radius(&valenx_mesh::Mesh::new("empty")).is_none());
    }

    #[test]
    fn mesh_min_triangle_quality_grades_element_shape() {
        use nalgebra::Vector3;
        // A one-triangle mesh from three corners.
        let tri_mesh = |a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>| {
            let mut mesh = valenx_mesh::Mesh::new("tri");
            mesh.nodes = vec![a, b, c];
            let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
            block.connectivity = vec![0, 1, 2];
            mesh.element_blocks.push(block);
            mesh
        };

        // An equilateral triangle is the perfect element: Q = 1.
        let equi = tri_mesh(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, 3.0_f64.sqrt() / 2.0, 0.0),
        );
        assert!(
            (mesh_min_triangle_quality(&equi).unwrap() - 1.0).abs() < 1e-12,
            "equilateral Q=1"
        );

        // A right isosceles triangle (legs 1,1) scores Q = √3/2 ≈ 0.866.
        let right_iso = tri_mesh(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        );
        assert!(
            (mesh_min_triangle_quality(&right_iso).unwrap() - 3.0_f64.sqrt() / 2.0).abs() < 1e-12,
            "right isosceles Q=√3/2"
        );

        // The canonical 1×2×3 box: each face splits into right triangles, and the
        // worst sits on the most-elongated (1×3) face → Q = 3√3/10 ≈ 0.5196.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut box_mesh = valenx_mesh::Mesh::new("box");
        box_mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        box_mesh.element_blocks.push(block);
        let qmin = mesh_min_triangle_quality(&box_mesh).expect("box has triangles");
        assert!(
            (qmin - 3.0 * 3.0_f64.sqrt() / 10.0).abs() < 1e-12,
            "box Qmin = 3√3/10, got {qmin}"
        );
        // No element can beat the equilateral ideal, and a real face is positive.
        assert!(qmin > 0.0 && qmin <= 1.0 + 1e-12, "0 < Qmin ≤ 1");

        // A near-collinear sliver collapses toward zero quality.
        let sliver = tri_mesh(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, 1e-6, 0.0),
        );
        assert!(
            mesh_min_triangle_quality(&sliver).unwrap() < 1e-3,
            "sliver Q→0"
        );

        // An empty mesh has no triangles to grade.
        assert!(mesh_min_triangle_quality(&valenx_mesh::Mesh::new("empty")).is_none());
    }

    #[test]
    fn mesh_diameter_is_the_longest_chord_not_the_bbox_diagonal() {
        use nalgebra::Vector3;
        // The 1×2×3 box: the longest chord is the space diagonal between opposite
        // corners (0,0,0) and (1,2,3) → √(1+4+9) = √14.
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut box_mesh = valenx_mesh::Mesh::new("box");
        box_mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let d = mesh_diameter(&box_mesh).expect("box has vertices");
        assert!(
            (d - 14.0_f64.sqrt()).abs() < 1e-12,
            "box diameter √14, got {d}"
        );

        // A non-axis-aligned flat triangle: the diameter is the longest side (the
        // base, 6), strictly less than the bbox diagonal √(6²+4²)=√52≈7.21 —
        // proving it is the true caliper extent, NOT the axis-aligned box diagonal.
        let mut tri = valenx_mesh::Mesh::new("tri");
        tri.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(3.0, 4.0, 0.0),
            Vector3::new(6.0, 0.0, 0.0),
        ];
        let dt = mesh_diameter(&tri).expect("triangle has vertices");
        assert!(
            (dt - 6.0).abs() < 1e-12,
            "triangle diameter = base 6, got {dt}"
        );
        assert!(
            dt < bbox_diagonal([6.0, 4.0, 0.0]) - 1e-9,
            "diameter is strictly below the bbox diagonal"
        );

        // A single vertex spans nothing; an empty mesh has no diameter.
        let mut one = valenx_mesh::Mesh::new("one");
        one.nodes = vec![Vector3::new(3.0, 4.0, 5.0)];
        assert_eq!(mesh_diameter(&one), Some(0.0));
        assert!(mesh_diameter(&valenx_mesh::Mesh::new("empty")).is_none());
    }

    #[test]
    fn mesh_sharp_edge_count_finds_the_box_feature_edges() {
        use nalgebra::Vector3;
        // The canonical 1×2×3 box: 12 cube edges (adjacent faces meet at 90°) plus
        // 6 in-face triangulation diagonals (coplanar, 0° fold).
        let (lx, ly, lz) = (1.0, 2.0, 3.0);
        let mut box_mesh = valenx_mesh::Mesh::new("box");
        box_mesh.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(lx, 0.0, 0.0),
            Vector3::new(lx, ly, 0.0),
            Vector3::new(0.0, ly, 0.0),
            Vector3::new(0.0, 0.0, lz),
            Vector3::new(lx, 0.0, lz),
            Vector3::new(lx, ly, lz),
            Vector3::new(0.0, ly, lz),
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7,
            3, 1, 2, 6, 1, 6, 5,
        ];
        box_mesh.element_blocks.push(block);
        // At a 30° crease angle the 12 right-angle cube edges are sharp; the 6
        // coplanar in-face diagonals (0° fold) are not.
        assert_eq!(mesh_sharp_edge_count(&box_mesh, 30.0_f64.to_radians()), 12);
        // Raise the threshold past the 90° cube fold and even those drop out → 0.
        assert_eq!(mesh_sharp_edge_count(&box_mesh, 100.0_f64.to_radians()), 0);

        // A coplanar two-triangle sheet folds nowhere; its lone boundary edges have
        // a single incident triangle each → no sharp edges.
        let mut sheet = valenx_mesh::Mesh::new("sheet");
        sheet.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
        ];
        let mut sblock = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        sblock.connectivity = vec![0, 1, 2, 1, 3, 2]; // two coplanar triangles sharing edge 1–2
        sheet.element_blocks.push(sblock);
        assert_eq!(
            mesh_sharp_edge_count(&sheet, 30.0_f64.to_radians()),
            0,
            "a flat sheet has no folds"
        );

        // An empty mesh has no edges to fold.
        assert_eq!(
            mesh_sharp_edge_count(&valenx_mesh::Mesh::new("empty"), 30.0_f64.to_radians()),
            0
        );
    }

    #[test]
    fn rebuild_reports_mass_from_density_and_volume() {
        // A 2×3×4 box (volume 24) at density 2.5 → mass 60.
        let mut bx = UiStep::base(Op::New, FeatureKind::Box);
        bx.dx = "2".into();
        bx.dy = "3".into();
        bx.dz = "4".into();
        let s = CadWorkbenchState {
            steps: vec![bx],
            density: 2.5,
            ..CadWorkbenchState::default()
        };
        let (_history, status) = rebuild_tree(&s).expect("box rebuilds");
        // mass = density × volume = 2.5 × 24 = 60.
        assert!(
            status.contains("mass"),
            "status should report mass: {status}"
        );
        assert!(
            status.contains("60.0000 mass"),
            "mass = 2.5×24 = 60: {status}"
        );
    }

    #[test]
    fn feature_tree_reports_a_no_base_body_error() {
        // A lone Cut step has no body to cut from — surfaces as an error,
        // not a panic, and pushes nothing.
        let s = CadWorkbenchState {
            steps: vec![UiStep::new_cylinder()], // a single Cut step
            ..CadWorkbenchState::default()
        };
        let err = rebuild_tree(&s).expect_err("a lone Cut must fail");
        assert!(!err.is_empty());
    }

    #[test]
    fn feature_tree_save_load_round_trips() {
        let s = CadWorkbenchState::default();
        let txt = save_string(&s).expect("serialize");
        let (params, steps) = load_from_string(&txt).expect("deserialize");
        assert_eq!(params.len(), s.params.len());
        assert_eq!(steps.len(), s.steps.len());
        // The reloaded tree rebuilds to the same punched cube.
        let loaded = CadWorkbenchState {
            params,
            steps,
            ..CadWorkbenchState::default()
        };
        let (history, _) = rebuild_tree(&loaded).expect("loaded tree rebuilds");
        assert!(
            history.last().unwrap()[0].faces() > 6,
            "reloaded tree still punches a hole"
        );
    }

    #[test]
    fn feature_tree_builds_new_primitives_standalone() {
        // Each new UI primitive, as a lone New step, rebuilds to a real body —
        // exercises to_step → rebuild for sphere / cone / torus.
        for mut step in [
            UiStep::new_sphere(),
            UiStep::new_cone(),
            UiStep::new_torus(),
        ] {
            step.op = Op::New;
            let s = CadWorkbenchState {
                steps: vec![step],
                ..CadWorkbenchState::default()
            };
            let (history, _) = rebuild_tree(&s).expect("primitive rebuilds");
            assert!(
                history.last().unwrap()[0].faces() > 0,
                "primitive has faces"
            );
        }
    }

    #[test]
    fn save_load_preserves_a_torus_step() {
        // Round-trip a new primitive through to_step → RON → ui_step_from.
        let torus = {
            let mut t = UiStep::new_torus();
            t.op = Op::New;
            t
        };
        let s = CadWorkbenchState {
            steps: vec![torus],
            ..CadWorkbenchState::default()
        };
        let txt = save_string(&s).expect("serialize");
        let (_params, steps) = load_from_string(&txt).expect("deserialize");
        assert_eq!(steps.len(), 1);
        assert_eq!(steps[0].kind, FeatureKind::Torus);
    }

    #[test]
    fn rotation_round_trips_and_rebuilds() {
        let mut step = UiStep::new_box();
        step.op = Op::New;
        step.rx = "30".into();
        step.rz = "45".into();
        let s = CadWorkbenchState {
            steps: vec![step],
            ..CadWorkbenchState::default()
        };
        // A rotated box rebuilds and tessellates to a non-empty mesh.
        let (history, _) = rebuild_tree(&s).expect("rotated box rebuilds");
        let mesh = tessellate_step(&history, 1, &[]).expect("tessellate");
        assert!(crate::mesh_loader::mesh_bounding_box(&mesh).is_some());
        // Rotation survives save → load.
        let txt = save_string(&s).expect("serialize");
        let (_params, steps) = load_from_string(&txt).expect("deserialize");
        assert_eq!(steps[0].rx, "30");
        assert_eq!(steps[0].rz, "45");
    }

    #[test]
    fn extrude_profile_round_trips_and_rebuilds() {
        let mut step = UiStep::new_extrude();
        step.op = Op::New;
        step.profile = vec![(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)];
        step.height = "1".into();
        let s = CadWorkbenchState {
            steps: vec![step],
            ..CadWorkbenchState::default()
        };
        // An extruded profile rebuilds and tessellates to a non-empty mesh.
        let (history, _) = rebuild_tree(&s).expect("extrude rebuilds");
        let mesh = tessellate_step(&history, 1, &[]).expect("tessellate");
        assert!(crate::mesh_loader::mesh_bounding_box(&mesh).is_some());
        // The profile survives save → load.
        let txt = save_string(&s).expect("serialize");
        let (_params, steps) = load_from_string(&txt).expect("deserialize");
        assert_eq!(steps[0].kind, FeatureKind::Extrude);
        assert_eq!(steps[0].profile.len(), 4);
    }

    #[test]
    fn two_new_steps_produce_two_bodies() {
        let mut a = UiStep::new_box();
        a.op = Op::New;
        let mut b = UiStep::new_box();
        b.op = Op::New; // a second New keeps a separate body
        b.x = "3".into();
        let s = CadWorkbenchState {
            steps: vec![a, b],
            ..CadWorkbenchState::default()
        };
        let (history, status) = rebuild_tree(&s).expect("two bodies rebuild");
        assert_eq!(
            history.last().unwrap().len(),
            2,
            "final step has two bodies"
        );
        assert!(status.contains("2 bodies"), "status: {status}");
        // Both bodies concatenate into one non-empty display mesh.
        let mesh = tessellate_step(&history, history.len(), &[]).expect("merge");
        assert!(crate::mesh_loader::mesh_bounding_box(&mesh).is_some());
    }

    #[test]
    fn body_visibility_filters_the_display() {
        let mut a = UiStep::new_box();
        a.op = Op::New;
        let mut b = UiStep::new_box();
        b.op = Op::New;
        b.x = "3".into();
        let s = CadWorkbenchState {
            steps: vec![a, b],
            ..CadWorkbenchState::default()
        };
        let (history, _) = rebuild_tree(&s).expect("two bodies");
        let last = history.len();
        // Both / one visible → a real mesh; none visible → an error.
        assert!(tessellate_step(&history, last, &[true, true]).is_ok());
        assert!(tessellate_step(&history, last, &[true, false]).is_ok());
        assert!(
            tessellate_step(&history, last, &[false, false]).is_err(),
            "hiding all bodies is an error, not an empty mesh"
        );
    }

    #[test]
    fn export_writes_a_binary_stl() {
        let s = CadWorkbenchState::default();
        let (history, _) = rebuild_tree(&s).expect("rebuild");
        let path =
            std::env::temp_dir().join(format!("valenx_cad_export_{}.stl", std::process::id()));
        export_stl(&history, history.len(), &path).expect("export");
        let len = std::fs::metadata(&path).expect("stl file exists").len();
        assert!(
            len > 84,
            "binary STL = 80B header + 4B count + triangles ({len} B)"
        );
        let _ = std::fs::remove_file(&path);
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    /// Draw the whole CAD panel once with accesskit enabled and return the
    /// emitted accessibility tree — the same tree a screen reader / AI
    /// UI-Automation driver consumes. `accesskit` is re-exported by egui, so no
    /// extra dependency is needed.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_cad_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn has_named_node(nodes: &[(NodeId, Node)], name: &str) -> bool {
        nodes.iter().any(|(_, n)| n.name() == Some(name))
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_cad_workbench);
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // (1) Default panel: the always-visible standalone DragValues (sketch
        // extrude height + material density) must be associated with their
        // captions — egui clears a DragValue's own Name, so without
        // `labelled_by` they are anonymous to a screen reader / AI driver.
        let mut app = ValenxApp::default();
        app.show_cad_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 2,
            "expected the standalone CAD numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every CAD DragValue must be labelled_by its caption (AI-drivable name)"
        );
        for caption in ["Extrude height", "Material density"] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }

        // (2) A feature tree containing an Extrude step renders the per-point
        // profile-edit loop, whose x/y spin buttons live inside a `for` and are
        // captioned per point (`P{j} x` / `P{j} y`). Verify those, too, are all
        // named — the loop-gated case the prompt called out.
        let mut app = ValenxApp::default();
        app.show_cad_workbench = true;
        app.cad.steps = vec![UiStep::new_extrude()];
        let nodes = draw_and_collect_nodes(&mut app);
        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // standalone (2) + a 4-point Extrude profile (4 × x,y = 8) = 10.
        assert!(
            spin_buttons.len() >= 10,
            "expected the Extrude profile-loop spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every CAD profile-loop DragValue must be labelled_by its per-point caption"
        );
        for caption in ["P0 x", "P0 y"] {
            assert!(
                has_named_node(&nodes, caption),
                "per-point caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }

    /// The agent `SetControl` bridge: each representative control sets the
    /// matching state field by its user-visible caption; an unknown caption, a
    /// wrong-type value, and an out-of-range step index all fail loud (Err, no
    /// panic) and leave state untouched.
    #[test]
    fn agent_set_sets_controls_and_rejects_bad_input() {
        let mut s = CadWorkbenchState::default();

        // -- Singletons: numeric, enum (by option name), bool ----------------
        s.agent_set("Material density", &AgentValue::Float(7.85))
            .expect("set material density");
        assert_eq!(s.density, 7.85);

        s.agent_set("Extrude height", &AgentValue::Float(3.5))
            .expect("set extrude height");
        assert_eq!(s.sketch_extrude_height, 3.5);

        s.agent_set("Sketch tool", &AgentValue::Str("Arc".into()))
            .expect("set sketch tool by name");
        assert_eq!(s.sketch_tool, SketchTool::Arc);

        s.agent_set("Snap to grid", &AgentValue::Bool(false))
            .expect("set snap-to-grid");
        assert!(!s.sketch_grid_snap);

        // `circle radius =` is the parameter *name* (a string), not a number.
        s.agent_set("circle radius =", &AgentValue::Str("base".into()))
            .expect("set radius-driving parameter name");
        assert_eq!(s.radius_param, "base");

        // -- Per-step fields (1-based index) ---------------------------------
        // The default tree is [Box(New), Cylinder(Cut)]; step 1 is the Box.
        s.agent_set("step 1 dx", &AgentValue::Str("size / 2".into()))
            .expect("set step 1 dx expression");
        assert_eq!(s.steps[0].dx, "size / 2");

        s.agent_set("step 1 op", &AgentValue::Str("Join".into()))
            .expect("set step 1 op by name");
        assert_eq!(s.steps[0].op, Op::Join);

        s.agent_set("step 2 shape", &AgentValue::Str("Sphere".into()))
            .expect("set step 2 shape by name");
        assert_eq!(s.steps[1].kind, FeatureKind::Sphere);

        // A bare number is accepted as the string form for an expression field.
        s.agent_set("step 1 dz", &AgentValue::Int(4))
            .expect("integer coerces to expression string");
        assert_eq!(s.steps[0].dz, "4");

        // -- Fail-loud cases (each returns Err; none panics) -----------------
        // Unknown caption.
        assert!(s
            .agent_set("does not exist", &AgentValue::Float(1.0))
            .is_err());
        // Unknown per-step field.
        assert!(s
            .agent_set("step 1 wat", &AgentValue::Str("x".into()))
            .is_err());
        // Type mismatch: a string into a numeric singleton.
        let before = s.density;
        assert!(s
            .agent_set("Material density", &AgentValue::Str("heavy".into()))
            .is_err());
        assert_eq!(s.density, before, "rejected set left density unchanged");
        // Type mismatch: a non-string into an enum control.
        assert!(s.agent_set("Sketch tool", &AgentValue::Float(2.0)).is_err());
        // Unknown enum option name.
        assert!(s
            .agent_set("step 1 shape", &AgentValue::Str("Dodecahedron".into()))
            .is_err());
        // Out-of-range step index (tree has 2 steps).
        assert!(s
            .agent_set("step 9 dx", &AgentValue::Str("1".into()))
            .is_err());
        // 0 is rejected (the index is 1-based).
        assert!(s
            .agent_set("step 0 dx", &AgentValue::Str("1".into()))
            .is_err());

        // -- Discovery list advertises the representative name space ----------
        let names = CadWorkbenchState::agent_control_names();
        for expect in ["Material density", "Sketch tool", "step 1 op", "step 1 dx"] {
            assert!(
                names.contains(&expect),
                "agent_control_names should advertise {expect:?}"
            );
        }
    }
}
