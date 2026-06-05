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

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which primitive a feature-tree step builds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FeatureKind {
    Box,
    Cylinder,
    Sphere,
    Cone,
    Torus,
    Extrude,
}

/// One UI-editable feature-tree step. Carries both the box and cylinder
/// dimension fields so toggling `kind` preserves whatever the user typed.
#[derive(Clone)]
struct UiStep {
    op: Op,
    kind: FeatureKind,
    /// Extrude profile — `(x, y)` points (literal coords, not expressions).
    profile: Vec<(f64, f64)>,
    dx: String,
    dy: String,
    dz: String,
    radius: String,
    height: String,
    top_radius: String,
    major: String,
    minor: String,
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
        };
        let mut step =
            Step::placed(self.op, feature, self.x.clone(), self.y.clone(), self.z.clone());
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
}

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
        }
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
            let circle = sk.add_circle(c, 1.0, 0.0, 0.0, 1.0).expect("centre is a point");
            sk.add_constraint(Constraint3D::CircleRadius { circle, target: r });
            match sk.solve() {
                Ok(rep) => {
                    let solved = sk.circle_radius(circle);
                    (Some(solved), format!("{:?} — circle radius = {solved:.4}", rep.status))
                }
                Err(e) => (None, format!("solve error: {e}")),
            }
        }
        Err(e) => (None, format!("radius parameter '{}': {e}", s.radius_param.trim())),
    };

    CadResults { resolved, solved_radius, status }
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
        Feature::Cone { base_radius, top_radius, height } => {
            us.kind = FeatureKind::Cone;
            us.radius = base_radius.clone();
            us.top_radius = top_radius.clone();
            us.height = height.clone();
        }
        Feature::Torus { major_radius, minor_radius } => {
            us.kind = FeatureKind::Torus;
            us.major = major_radius.clone();
            us.minor = minor_radius.clone();
        }
        Feature::Extrude { profile, height } => {
            us.kind = FeatureKind::Extrude;
            us.profile = profile.clone();
            us.height = height.clone();
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
            second[0] +=
                v / 20.0 * (a.x * a.x + b.x * b.x + c.x * c.x + (a.x + b.x + c.x).powi(2));
            second[1] +=
                v / 20.0 * (a.y * a.y + b.y * b.y + c.y * c.y + (a.y + b.y + c.y).powi(2));
            second[2] +=
                v / 20.0 * (a.z * a.z + b.z * b.z + c.z * c.z + (a.z + b.z + c.z).powi(2));
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
    let status = format!(
        "{nbodies} bodies · {faces} faces · {volume:.4} u³ · {mass:.4} mass · {area:.4} u²{sv_str}{bbox_str}{fill_str}{sphericity_str}{centroid_str}{gyration_str}{watertight_str} · {} steps",
        s.steps.len()
    );
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
    }
}

/// A narrow single-line editor for a dimension / placement expression.
fn dim_edit(ui: &mut egui::Ui, v: &mut String) {
    ui.add(egui::TextEdit::singleline(v).desired_width(52.0));
}

/// Draw the parametric-CAD workbench (a no-op unless toggled on via
/// View → Parametric CAD).
pub fn draw_cad_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_cad_workbench {
        return;
    }
    egui::SidePanel::right("valenx_cad_workbench")
        .resizable(true)
        .default_width(340.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            ui.heading("Parametric CAD");
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
                    ui.label(egui::RichText::new("Parameters (name = expression)").strong());
                    let mut remove: Option<usize> = None;
                    for (i, (name, expr)) in s.params.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            ui.add(
                                egui::TextEdit::singleline(name)
                                    .desired_width(80.0)
                                    .hint_text("name"),
                            );
                            ui.label("=");
                            ui.add(
                                egui::TextEdit::singleline(expr)
                                    .desired_width(130.0)
                                    .hint_text("expr"),
                            );
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
                        ui.label("Material density").on_hover_text(
                            "Mass per unit volume; the rebuild status shows \
                             mass = density × solid volume.",
                        );
                        ui.add(
                            egui::DragValue::new(&mut s.density)
                                .speed(0.1)
                                .range(0.0..=1.0e9)
                                .suffix(" /u³"),
                        );
                    });

                    // ---- Sketch: parameter-driven circle radius ----
                    ui.add_space(6.0);
                    ui.label(egui::RichText::new("Sketch").strong());
                    ui.horizontal(|ui| {
                        ui.label("circle radius =");
                        ui.add(
                            egui::TextEdit::singleline(&mut s.radius_param)
                                .desired_width(100.0)
                                .hint_text("parameter"),
                        );
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
                                ui.label(format!("{}.", i + 1));
                                egui::ComboBox::from_id_source(("cad_op", i))
                                    .selected_text(op_label(st.op))
                                    .width(92.0)
                                    .show_ui(ui, |ui| {
                                        for op in [Op::New, Op::Join, Op::Cut, Op::Intersect] {
                                            ui.selectable_value(&mut st.op, op, op_label(op));
                                        }
                                    });
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
                                    });
                                if ui.small_button("✕").clicked() {
                                    remove_step = Some(i);
                                }
                            });
                            match st.kind {
                                FeatureKind::Box => {
                                    ui.horizontal(|ui| {
                                        ui.label("dx,dy,dz");
                                        dim_edit(ui, &mut st.dx);
                                        dim_edit(ui, &mut st.dy);
                                        dim_edit(ui, &mut st.dz);
                                    });
                                }
                                FeatureKind::Cylinder => {
                                    ui.horizontal(|ui| {
                                        ui.label("r, h");
                                        dim_edit(ui, &mut st.radius);
                                        dim_edit(ui, &mut st.height);
                                    });
                                }
                                FeatureKind::Sphere => {
                                    ui.horizontal(|ui| {
                                        ui.label("r");
                                        dim_edit(ui, &mut st.radius);
                                    });
                                }
                                FeatureKind::Cone => {
                                    ui.horizontal(|ui| {
                                        ui.label("base r, top r, h");
                                        dim_edit(ui, &mut st.radius);
                                        dim_edit(ui, &mut st.top_radius);
                                        dim_edit(ui, &mut st.height);
                                    });
                                }
                                FeatureKind::Torus => {
                                    ui.horizontal(|ui| {
                                        ui.label("major, minor");
                                        dim_edit(ui, &mut st.major);
                                        dim_edit(ui, &mut st.minor);
                                    });
                                }
                                FeatureKind::Extrude => {
                                    ui.horizontal(|ui| {
                                        ui.label("height");
                                        dim_edit(ui, &mut st.height);
                                    });
                                    ui.label(
                                        egui::RichText::new("profile (x, y) — ≥3 points")
                                            .weak()
                                            .small(),
                                    );
                                    let mut rm_pt: Option<usize> = None;
                                    for (j, pt) in st.profile.iter_mut().enumerate() {
                                        ui.horizontal(|ui| {
                                            ui.add(
                                                egui::DragValue::new(&mut pt.0)
                                                    .speed(0.1)
                                                    .prefix("x "),
                                            );
                                            ui.add(
                                                egui::DragValue::new(&mut pt.1)
                                                    .speed(0.1)
                                                    .prefix("y "),
                                            );
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
                            }
                            ui.horizontal(|ui| {
                                ui.label("at x,y,z");
                                dim_edit(ui, &mut st.x);
                                dim_edit(ui, &mut st.y);
                                dim_edit(ui, &mut st.z);
                            });
                            ui.horizontal(|ui| {
                                ui.label("rot x,y,z°");
                                dim_edit(ui, &mut st.rx);
                                dim_edit(ui, &mut st.ry);
                                dim_edit(ui, &mut st.rz);
                            });
                        });
                    }
                    if let Some(i) = remove_step {
                        s.steps.remove(i);
                    }
                    ui.horizontal_wrapped(|ui| {
                        if ui.button("+ Box").clicked() {
                            s.steps.push(UiStep::new_box());
                        }
                        if ui.button("+ Cylinder").clicked() {
                            s.steps.push(UiStep::new_cylinder());
                        }
                        if ui.button("+ Sphere").clicked() {
                            s.steps.push(UiStep::new_sphere());
                        }
                        if ui.button("+ Cone").clicked() {
                            s.steps.push(UiStep::new_cone());
                        }
                        if ui.button("+ Torus").clicked() {
                            s.steps.push(UiStep::new_torus());
                        }
                        if ui.button("+ Extrude").clicked() {
                            s.steps.push(UiStep::new_extrude());
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
                                        let res = valenx_core::io_caps::atomic_write_str(
                                            &path, &txt,
                                        );
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
                    if let Some(res) = &s.tree_status {
                        match res {
                            Ok(status) => ui.colored_label(
                                egui::Color32::from_rgb(80, 220, 120),
                                status,
                            ),
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
                        let resp = ui
                            .add(egui::Slider::new(&mut s.scrub, 1..=n).integer().text("step"));
                        if resp.changed() {
                            let scrub = s.scrub;
                            let mesh =
                                s.history.as_ref().and_then(|h| {
                                    tessellate_step(h, scrub, &s.body_visible).ok()
                                });
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
        });

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
        assert_eq!(radius.1.as_ref().ok().map(|v| (v * 1e4).round() / 1e4), Some(5.0));
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
        assert!(r.status.to_lowercase().contains("cyclic"), "status: {}", r.status);
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
        assert!((v - 24.0).abs() < 1e-6, "2×3×4 box volume should be 24, got {v}");
        // The status line surfaces the volume in cubic model units.
        assert!(status.contains("u³"), "status should report volume: {status}");

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
        assert!((a - 52.0).abs() < 1e-6, "2×3×4 box area should be 52, got {a}");
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
        assert!(status.contains("bbox"), "status should report bbox: {status}");
        // A box fills its bounding box exactly — fill fraction ≈ 1.0.
        let vol = total_volume(history.last().expect("snapshot"));
        let fill = fill_fraction(vol, dims).expect("fill");
        assert!((fill - 1.0).abs() < 1e-3, "box should fill its bbox: {fill}");
        assert!(status.contains("fill"), "status should report fill: {status}");
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
        assert!((c[2] - 1.0).abs() < 1e-9, "cz {} (should be h/4 = 1.0)", c[2]);
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
        assert!(status.contains("mass"), "status should report mass: {status}");
        assert!(status.contains("60.0000 mass"), "mass = 2.5×24 = 60: {status}");
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
        for mut step in [UiStep::new_sphere(), UiStep::new_cone(), UiStep::new_torus()] {
            step.op = Op::New;
            let s = CadWorkbenchState {
                steps: vec![step],
                ..CadWorkbenchState::default()
            };
            let (history, _) = rebuild_tree(&s).expect("primitive rebuilds");
            assert!(history.last().unwrap()[0].faces() > 0, "primitive has faces");
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
        assert_eq!(history.last().unwrap().len(), 2, "final step has two bodies");
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
        assert!(len > 84, "binary STL = 80B header + 4B count + triangles ({len} B)");
        let _ = std::fs::remove_file(&path);
    }
}
