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
    // Worst triangle shape quality (4√3·A/Σℓ²; 1 = equilateral, →0 = sliver) —
    // the mesh-quality gate, distinct from the shape/topology measures above.
    let quality_str = mesh
        .as_ref()
        .and_then(mesh_min_triangle_quality)
        .map(|q| format!(" · Qmin {q:.2}"))
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
    let status = format!(
        "{nbodies} bodies · {faces} faces · {volume:.4} u³ · {mass:.4} mass · {area:.4} u²{sv_str}{bbox_str}{fill_str}{sphericity_str}{centroid_str}{gyration_str}{moments_str}{watertight_str}{open_str}{hole_str}{euler_str}{shells_str}{encl_str}{diam_str}{mean_edge_str}{quality_str}{aspect_str}{sharp_str}{crease_str} · {} steps",
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
        assert!((i1 - m / 12.0 * (ly * ly + lz * lz)).abs() < 1e-9, "I1 {i1}"); // 6.5
        assert!((i2 - m / 12.0 * (lx * lx + lz * lz)).abs() < 1e-9, "I2 {i2}"); // 5.0
        assert!((i3 - m / 12.0 * (lx * lx + ly * ly)).abs() < 1e-9, "I3 {i3}"); // 2.5

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
        assert!((ar - 10.0_f64.sqrt()).abs() < 1e-9, "box max aspect ratio = √10, got {ar}");

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
        assert!((mesh_max_aspect_ratio(&eq).unwrap() - 1.0).abs() < 1e-9, "equilateral → 1");

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
        assert!((mesh_max_aspect_ratio(&r345).unwrap() - 5.0 / 3.0).abs() < 1e-9, "3-4-5 → 5/3");

        // No triangles → None.
        assert_eq!(mesh_max_aspect_ratio(&valenx_mesh::Mesh::new("empty")), None);
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
        assert!((theta - FRAC_PI_2).abs() < 1e-9, "box sharpest fold = π/2, got {theta}");

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
        assert_eq!(mesh_max_dihedral_angle(&valenx_mesh::Mesh::new("empty")), None);
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
        assert!((mesh_mean_edge_length(&tri).unwrap() - 4.0).abs() < 1e-9, "3-4-5 mean edge = 4");

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
        assert!((mean - expected).abs() < 1e-9, "box mean edge ≈ 2.334, got {mean}");

        // No edges → None.
        assert_eq!(mesh_mean_edge_length(&valenx_mesh::Mesh::new("empty")), None);
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
        assert_eq!(mesh_open_edge_count(&mesh), 0, "a closed box has no open edges");

        // Drop one triangle → its 3 edges fall from multiplicity 2 to 1.
        let mut open = mesh.clone();
        open.element_blocks[0].connectivity.truncate(closed.len() - 3);
        assert_eq!(mesh_open_edge_count(&open), 3, "a box missing one triangle has 3 open edges");

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
        assert_eq!(mesh_open_edge_count(&tri), 3, "a lone triangle has a 3-edge boundary");

        // Two coplanar triangles sharing edge 1–2: 5 edges, the shared one closed → 4 open.
        let mut sheet = tri.clone();
        sheet.element_blocks[0].connectivity = vec![0, 1, 2, 1, 3, 2];
        assert_eq!(mesh_open_edge_count(&sheet), 4, "a 2-triangle sheet has a 4-edge boundary");

        // An empty mesh has no edges.
        assert_eq!(mesh_open_edge_count(&valenx_mesh::Mesh::new("empty")), 0);
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
        assert_eq!(mesh_boundary_loop_count(&make(closed.clone())), 0, "closed box → 0 holes");

        // Drop the last triangle → one triangular gap (3 open edges, but ONE loop).
        let mut one = closed.clone();
        one.truncate(closed.len() - 3);
        assert_eq!(mesh_boundary_loop_count(&make(one)), 1, "one missing triangle → 1 hole");

        // Omit two vertex-disjoint triangles — T0=[0,2,1] (verts {0,1,2}) and
        // T3=[4,6,7] (verts {4,6,7}) — leaving two separate triangular holes.
        // (mesh_open_edge_count would report 6 open edges here; the loop count is 2.)
        let two: Vec<u32> = vec![
            0, 3, 2, 4, 5, 6, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6,
            5,
        ];
        assert_eq!(mesh_boundary_loop_count(&make(two)), 2, "two disjoint gaps → 2 holes");

        // An empty mesh has no boundary loops.
        assert_eq!(mesh_boundary_loop_count(&valenx_mesh::Mesh::new("empty")), 0);
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
        assert_eq!(mesh_euler_characteristic(&open), Some(1), "box minus a face: χ = 1");

        // Two disjoint boxes → χ is additive over components: 2 + 2 = 4.
        let mut two = mesh.clone();
        let n = two.nodes.len() as u32;
        let shifted: Vec<Vector3<f64>> =
            mesh.nodes.iter().map(|p| *p + Vector3::new(10.0, 0.0, 0.0)).collect();
        two.nodes.extend(shifted);
        let mut block2 = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block2.connectivity = closed.iter().map(|&i| i + n).collect();
        two.element_blocks.push(block2);
        assert_eq!(mesh_euler_characteristic(&two), Some(4), "two disjoint boxes: χ = 4");

        // An empty mesh has no Euler characteristic.
        assert_eq!(mesh_euler_characteristic(&valenx_mesh::Mesh::new("empty")), None);
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
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0,
            7, 3, 1, 2, 6, 1, 6, 5,
        ];
        let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block.connectivity = closed.clone();
        mesh.element_blocks.push(block);
        assert_eq!(mesh_shell_count(&mesh), 1, "one box is a single shell");

        // Two disjoint boxes (a second block on fresh node indices) → two shells.
        let mut two = mesh.clone();
        let n = two.nodes.len() as u32;
        let shifted: Vec<Vector3<f64>> =
            mesh.nodes.iter().map(|p| *p + Vector3::new(10.0, 0.0, 0.0)).collect();
        two.nodes.extend(shifted);
        let mut block2 = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
        block2.connectivity = closed.iter().map(|&i| i + n).collect();
        two.element_blocks.push(block2);
        assert_eq!(mesh_shell_count(&two), 2, "two disjoint boxes are two shells");

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
        assert!((d - 14.0_f64.sqrt()).abs() < 1e-12, "box diameter √14, got {d}");

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
        assert!((dt - 6.0).abs() < 1e-12, "triangle diameter = base 6, got {dt}");
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
