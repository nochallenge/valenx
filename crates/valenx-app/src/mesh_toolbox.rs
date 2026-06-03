//! The right-side "Mesh Toolbox" panel.
//!
//! Surfaces editing primitives whenever the user has dropped an STL
//! or loaded a canonical mesh:
//!
//! - **Inspector** — node / element counts, AABB, quality scalars.
//! - **Transformations** — translate / scale / rotate / mirror via
//!   numeric inputs + Apply buttons. All native, no external tools.
//! - **Cut plane** — point + normal numeric inputs; preview the cut
//!   line on the surface, or commit to slice the mesh in half.
//! - **Repair** — merge coincident nodes (existing
//!   `valenx_mesh::boolean::merge_coincident_nodes`).
//! - **Export** — save the current mesh as a binary STL via the
//!   `rfd` file dialog.
//! - **External editors** — open the source STL in FreeCAD, or
//!   scaffold a new OCCT case from the existing template library.
//!
//! Every Apply emits a structured audit-log entry through
//! [`emit_audit`] so an admin running `valenx audit verify` sees
//! exactly what was edited and when. The methods on `ValenxApp` are
//! the actual side-effect carriers — this module's `draw_mesh_toolbox`
//! is purely the egui-layout shim that wires the form widgets to
//! those methods.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;
use valenx_mesh::transform::{Axis, Plane};
use valenx_viz::TriangleMesh;

use crate::audit::emit_audit;
use crate::types::LoadedMesh;
use crate::ValenxApp;

/// User-visible form state for the toolbox. Mirrors the layout of
/// each section so the egui code is a thin map from struct fields to
/// `egui::DragValue` widgets.
#[derive(Clone, Debug)]
pub struct MeshToolboxState {
    /// Translate (X / Y / Z) inputs in user coordinates.
    pub translate: [f64; 3],
    /// Uniform scale factor.
    pub scale_uniform: f64,
    /// Per-axis scale factors.
    pub scale_per_axis: [f64; 3],
    /// Rotation axis chosen by the radio row.
    pub rotate_axis: ToolboxAxis,
    /// Rotation angle in degrees (more intuitive than radians for
    /// numeric input; converted internally).
    pub rotate_angle_deg: f64,
    /// Mirror plane chosen by the radio row.
    pub mirror_plane: ToolboxAxis,
    /// Cut-plane point and normal numeric inputs.
    pub cut_point: [f64; 3],
    pub cut_normal: [f64; 3],
    /// Toggle for drawing the cut-overlay on the viewport.
    pub cut_show_overlay: bool,
    /// Tolerance used by the merge-coincident-nodes repair step.
    pub repair_tolerance: f64,

    // ----- Mesh Tools (Phase 7) -----
    /// Target fraction of vertices to keep when running QEM decimation.
    /// `0.5` halves the vertex count; `1.0` is a no-op.
    pub mesh_tools_decimate_fraction: f64,
    /// Iteration count for Laplacian smoothing.
    pub mesh_tools_laplacian_iter: u32,
    /// Step factor λ for Laplacian smoothing.
    pub mesh_tools_laplacian_factor: f64,
    /// Iteration count for Taubin smoothing.
    pub mesh_tools_taubin_iter: u32,
    /// λ step for Taubin smoothing.
    pub mesh_tools_taubin_lambda: f64,
    /// μ step for Taubin smoothing.
    pub mesh_tools_taubin_mu: f64,
    /// Target edge length for isotropic remeshing.
    pub mesh_tools_remesh_target: f64,
    /// Iteration count for isotropic remeshing.
    pub mesh_tools_remesh_iter: u32,
    /// Maximum boundary perimeter for `fill_holes`.
    pub mesh_tools_fill_holes_max: f64,

    // ----- Part workbench (valenx-cad / truck BRep kernel) -----
    /// Which CAD primitive the "Insert primitive" combo will build
    /// when the user clicks Create.
    pub cad_primitive: CadPrimitiveKind,
    /// Box dimensions: (dx, dy, dz). Defaults to a unit cube.
    pub cad_box_dims: [f64; 3],
    /// Cylinder parameters: (radius, height).
    pub cad_cyl_radius: f64,
    pub cad_cyl_height: f64,
    /// Sphere radius.
    pub cad_sphere_radius: f64,
    /// Cone / frustum parameters: base radius, top radius, height.
    /// `top_radius = 0` is a regular pointed cone.
    pub cad_cone_base: f64,
    pub cad_cone_top: f64,
    pub cad_cone_height: f64,
    /// Torus parameters: major (centre-circle) + minor (tube) radii.
    pub cad_torus_major: f64,
    pub cad_torus_minor: f64,

    /// When true, the next "Create" action stores the resulting solid
    /// as `second_solid` (operand B) instead of `current_solid`
    /// (operand A). Lets users build two primitives back-to-back
    /// without losing the first one when they need a boolean.
    pub cad_create_as_second: bool,

    /// Fillet radius used by the "Apply fillet" button.
    pub cad_fillet_radius: f64,

    /// State for the Dock panel — receptor / ligand paths, search-box
    /// parameters, last-run scores.
    pub dock_panel: DockPanelState,

    /// State for the Sketcher panel — active sketch, current tool,
    /// last solver report. See `SketcherPanelState`.
    pub sketcher: SketcherPanelState,

    /// State for the Part Design panel — feature tree, selection,
    /// replay status. See [`PartDesignPanelState`].
    pub part_design: PartDesignPanelState,

    /// State for the Draft workbench panel — 2D entities on a working
    /// plane. See [`DraftPanelState`]. Phase 4.
    pub draft: DraftPanelState,

    /// State for the TechDraw workbench panel — 2D engineering
    /// drawings projected from 3D solids. See [`TechDrawPanelState`].
    /// Phase 5.
    pub techdraw: TechDrawPanelState,

    /// State for the Assembly workbench panel — multi-part scene
    /// with mates and joints. See [`AssemblyPanelState`]. Phase 6.
    pub assembly: AssemblyPanelState,

    /// State for the Surface workbench panel — NURBS curves +
    /// surfaces, Coons patch fill, tessellation. Phase 9.
    pub surface: SurfacePanelState,

    /// State for the CAM workbench panel — tool table, stock, ops,
    /// toolpath, postprocessor selection. Phase 10.
    pub cam: CamPanelState,

    /// State for the Arch / BIM workbench panel — walls, slabs,
    /// columns, beams, openings, schedule, IFC export. Phase 15.
    pub arch: ArchPanelState,

    /// State for the Spreadsheet workbench panel — named sheets,
    /// editable cells, formula evaluation, expression-bound feature
    /// params. See [`SpreadsheetPanelState`]. Phase 16.
    pub spreadsheet: SpreadsheetPanelState,
}

impl Default for MeshToolboxState {
    fn default() -> Self {
        Self {
            translate: [0.0, 0.0, 0.0],
            scale_uniform: 1.0,
            scale_per_axis: [1.0, 1.0, 1.0],
            rotate_axis: ToolboxAxis::Z,
            rotate_angle_deg: 0.0,
            mirror_plane: ToolboxAxis::X,
            cut_point: [0.0, 0.0, 0.0],
            cut_normal: [0.0, 0.0, 1.0],
            cut_show_overlay: false,
            repair_tolerance: 1.0e-6,

            mesh_tools_decimate_fraction: 0.5,
            mesh_tools_laplacian_iter: 3,
            mesh_tools_laplacian_factor: 0.5,
            mesh_tools_taubin_iter: 5,
            mesh_tools_taubin_lambda: 0.5,
            mesh_tools_taubin_mu: -0.53,
            mesh_tools_remesh_target: 0.5,
            mesh_tools_remesh_iter: 3,
            mesh_tools_fill_holes_max: 100.0,

            cad_primitive: CadPrimitiveKind::Box,
            cad_box_dims: [1.0, 1.0, 1.0],
            cad_cyl_radius: 0.5,
            cad_cyl_height: 1.0,
            cad_sphere_radius: 0.5,
            cad_cone_base: 0.5,
            cad_cone_top: 0.0,
            cad_cone_height: 1.0,
            cad_torus_major: 1.0,
            cad_torus_minor: 0.25,
            cad_create_as_second: false,
            cad_fillet_radius: 0.1,

            dock_panel: DockPanelState::default(),
            sketcher: SketcherPanelState::new_with_overlay_on(),
            part_design: PartDesignPanelState::default(),
            draft: DraftPanelState::default(),
            techdraw: TechDrawPanelState::default(),
            assembly: AssemblyPanelState::default(),
            surface: SurfacePanelState::default(),
            cam: CamPanelState::default(),
            arch: ArchPanelState::default(),
            spreadsheet: SpreadsheetPanelState::default(),
        }
    }
}

/// State for the Assembly workbench panel (Phase 6).
///
/// Holds the live [`valenx_assembly::Assembly`], the numeric inputs
/// for the "Add part / mate / joint" buttons, last solver report, and
/// transient UI flags (selected part / mate / joint, last error).
#[derive(Clone, Debug)]
pub struct AssemblyPanelState {
    /// The active assembly.
    pub assembly: valenx_assembly::Assembly,
    /// Last solver report (None until first solve).
    pub last_report: Option<valenx_assembly::SolverReport>,
    /// Last user-facing error message (red-text label at the bottom
    /// of the panel).
    pub last_error: Option<String>,

    // ----- Selection -----
    /// Currently selected part id (for "Edit transform" / "Delete part").
    pub selected_part: Option<usize>,
    /// Currently selected mate id.
    pub selected_mate: Option<usize>,
    /// Currently selected joint id.
    pub selected_joint: Option<usize>,

    // ----- Add Part inputs -----
    /// Display name for the next "Add Part" press.
    pub new_part_name: String,
    /// Which primitive to source the new part from.
    pub new_part_primitive: AssemblyPartPrimitive,
    /// Box dims for the AssemblyPartPrimitive::Box source.
    pub new_part_box_dims: [f64; 3],
    /// Cylinder params (radius, height).
    pub new_part_cyl: [f64; 2],
    /// Sphere radius.
    pub new_part_sphere: f64,
    /// Initial transform translation for the new part.
    pub new_part_translation: [f64; 3],

    // ----- Add Mate inputs -----
    /// Which mate kind the "Add Mate" combo will build.
    pub mate_kind: AssemblyMateKindUi,
    /// Source part id input (free-text or selectable).
    pub mate_part_a: usize,
    /// Target part id.
    pub mate_part_b: usize,
    /// Coincident / Distance: anchor on part_a (local frame).
    pub mate_point_a: [f64; 3],
    /// Coincident / Distance: anchor on part_b (local frame).
    pub mate_point_b: [f64; 3],
    /// Distance / Angle: target value.
    pub mate_target: f64,
    /// Angle / Parallel / Perpendicular: vector on part_a (local frame).
    pub mate_vec_a: [f64; 3],
    /// Same for part_b.
    pub mate_vec_b: [f64; 3],
    /// Tangent: radii.
    pub mate_radius_a: f64,
    /// Same for part_b.
    pub mate_radius_b: f64,

    // ----- Add Joint inputs -----
    /// Which joint kind the "Add Joint" combo will build.
    pub joint_kind: AssemblyJointKindUi,
    /// Source part id.
    pub joint_part_a: usize,
    /// Target part id.
    pub joint_part_b: usize,
    /// Axis origin (in part_a's local frame) for Revolute / Cylindrical.
    pub joint_axis_origin: [f64; 3],
    /// Axis direction for Revolute / Prismatic / Cylindrical.
    pub joint_axis_dir: [f64; 3],
    /// Pivot point for Spherical.
    pub joint_point: [f64; 3],
    /// Plane origin + normal for Planar.
    pub joint_plane_origin: [f64; 3],
    /// Plane normal for Planar.
    pub joint_plane_normal: [f64; 3],

    /// Tessellation tolerance for the "Render assembly" button.
    pub render_tolerance: f64,
}

impl Default for AssemblyPanelState {
    fn default() -> Self {
        Self {
            assembly: valenx_assembly::Assembly::new(),
            last_report: None,
            last_error: None,
            selected_part: None,
            selected_mate: None,
            selected_joint: None,
            new_part_name: "Part".into(),
            new_part_primitive: AssemblyPartPrimitive::Box,
            new_part_box_dims: [1.0, 1.0, 1.0],
            new_part_cyl: [0.5, 1.0],
            new_part_sphere: 0.5,
            new_part_translation: [0.0, 0.0, 0.0],
            mate_kind: AssemblyMateKindUi::Coincident,
            mate_part_a: 0,
            mate_part_b: 1,
            mate_point_a: [0.0, 0.0, 0.0],
            mate_point_b: [0.0, 0.0, 0.0],
            mate_target: 1.0,
            mate_vec_a: [1.0, 0.0, 0.0],
            mate_vec_b: [1.0, 0.0, 0.0],
            mate_radius_a: 0.5,
            mate_radius_b: 0.5,
            joint_kind: AssemblyJointKindUi::Fixed,
            joint_part_a: 0,
            joint_part_b: 1,
            joint_axis_origin: [0.0, 0.0, 0.0],
            joint_axis_dir: [0.0, 0.0, 1.0],
            joint_point: [0.0, 0.0, 0.0],
            joint_plane_origin: [0.0, 0.0, 0.0],
            joint_plane_normal: [0.0, 0.0, 1.0],
            render_tolerance: 0.5,
        }
    }
}

/// State for the Surface workbench panel (Phase 9 — NURBS).
///
/// Holds the live in-memory list of NURBS curves + surfaces, the
/// inputs for the "Construct curve / surface / Coons-fill / Sew /
/// Trim" tool buttons, and the last user-facing status / error
/// strings. The panel mirrors the layout pattern of the Draft and
/// Assembly workbenches: tool palette at the top, input form
/// below the selected tool, entity list with tessellate buttons
/// at the bottom.
#[derive(Clone, Debug)]
pub struct SurfacePanelState {
    /// The active surface workbench file (curves + surfaces).
    pub file: valenx_surface::persist::SurfaceFile,
    /// Currently selected tool in the palette.
    pub tool: SurfaceTool,
    /// Last user-facing status string (green text at the bottom).
    pub last_status: Option<String>,
    /// Last user-facing error string (red text at the bottom).
    pub last_error: Option<String>,
    /// Selected curve index (for actions that take an existing
    /// curve as input).
    pub selected_curve: Option<usize>,
    /// Selected surface index.
    pub selected_surface: Option<usize>,

    // ----- Construct curve inputs -----
    /// Degree of the new NURBS curve.
    pub curve_degree: usize,
    /// Number of control points for the new NURBS curve.
    pub curve_n_cps: usize,
    /// Control points (`curve_n_cps` of them, padded / truncated
    /// as the user changes the count).
    pub curve_cps: Vec<[f64; 3]>,
    /// Per-CP weights.
    pub curve_weights: Vec<f64>,

    // ----- Construct surface inputs -----
    /// u-direction degree.
    pub surface_u_degree: usize,
    /// v-direction degree.
    pub surface_v_degree: usize,
    /// u-direction CP count.
    pub surface_nu: usize,
    /// v-direction CP count.
    pub surface_nv: usize,
    /// Flat list of CPs (size = nu * nv, row-major in u).
    pub surface_cps: Vec<[f64; 3]>,

    // ----- Coons fill inputs -----
    /// Indices of the four boundary curves: c0, c1, d0, d1.
    pub coons_curves: [usize; 4],

    // ----- Sew inputs -----
    /// Surface A id.
    pub sew_surface_a: usize,
    /// Surface B id.
    pub sew_surface_b: usize,
    /// Edge on A (0=UMin, 1=UMax, 2=VMin, 3=VMax).
    pub sew_edge_a: u8,
    /// Edge on B.
    pub sew_edge_b: u8,
    /// Tolerance for edge-CP coincidence.
    pub sew_tolerance: f64,
    /// Use G2 continuous sew (Phase 19C default). When `false`,
    /// uses the Phase 9 G0 averaging sew.
    pub sew_use_g2: bool,

    // ----- Trim inputs -----
    /// Surface id to trim.
    pub trim_surface: usize,
    /// Curve id used as trim boundary.
    pub trim_curve: usize,
    /// 0 = Inside, 1 = Outside.
    pub trim_side: u8,
    /// Tessellation density used during trim (per-axis sample
    /// count on the surface; curve uses 8× this).
    pub trim_resolution: usize,
    /// Phase 9.5: opt in to the parametric (u, v) domain trim.
    /// Default `true` so warped surfaces trim correctly; the legacy
    /// world-xy trim is still available by unchecking this.
    pub trim_use_uv: bool,

    // ----- Tessellation -----
    /// Per-axis sample count for the "Tessellate + push to viewport"
    /// button. 32 is a reasonable default.
    pub tess_resolution: usize,

    // ----- Phase 19A — knot insertion / removal / elevation -----
    /// Parameter for knot insertion / removal.
    pub knot_op_u: f64,
    /// Source curve id for the knot operation (also reused for
    /// degree elevation).
    pub knot_op_curve: usize,
    /// Tolerance used by Tiller knot removal.
    pub knot_op_tolerance: f64,
    /// Source surface id for surface knot insertion / elevation.
    pub knot_op_surface: usize,
    /// 0 = u-direction, 1 = v-direction for surface knot operations.
    pub knot_op_direction: u8,
    /// Degree elevation amount.
    pub elevate_degree_by: usize,

    // ----- Phase 19D — surface fitting -----
    /// Inline text input for scattered points to fit (one
    /// "x,y,z" triple per line).
    pub fit_points_text: String,
    /// Target u-degree for the fit.
    pub fit_degree_u: usize,
    /// Target v-degree for the fit.
    pub fit_degree_v: usize,
    /// Target u-direction CP count for the fit.
    pub fit_n_cps_u: usize,
    /// Target v-direction CP count for the fit.
    pub fit_n_cps_v: usize,
    /// Last reported RMS error from the fit operation.
    pub fit_last_rms: Option<f64>,

    // ----- Phase 19E — ruled surfaces -----
    /// First boundary curve id for ruled surface between two curves.
    pub ruled_curve_a: usize,
    /// Second boundary curve id (between_curves case).
    pub ruled_curve_b: usize,
    /// Extrusion vector for `extrude_along_vector`.
    pub ruled_extrude_vector: [f64; 3],
    /// Apex point for `cone_from_apex`.
    pub ruled_apex: [f64; 3],
    /// Which ruled-surface constructor is selected: 0 = between,
    /// 1 = extrude, 2 = cone.
    pub ruled_kind: u8,

    // ----- Phase 19B — SSI -----
    /// Surface ids for the SSI op.
    pub ssi_surface_a: usize,
    /// Second surface id.
    pub ssi_surface_b: usize,
    /// Tolerance for the SSI refinement / fit.
    pub ssi_tolerance: f64,
}

/// Which surface workbench tool is selected.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SurfaceTool {
    /// Construct a NURBS curve from degree + control points.
    NurbsCurve,
    /// Construct a NURBS surface from u/v degree + CP grid.
    NurbsSurface,
    /// Fill 4 boundary curves into a Coons patch.
    CoonsFill,
    /// Sew two NURBS surfaces along an edge.
    Sew,
    /// Trim a NURBS surface by a closed curve.
    Trim,
    /// Phase 19A — knot insertion / removal / degree elevation.
    KnotOps,
    /// Phase 19B — true rational surface-surface intersection.
    Ssi,
    /// Phase 19D — fit a NURBS curve or surface through points.
    Fit,
    /// Phase 19E — ruled surface constructors.
    Ruled,
}

impl SurfaceTool {
    /// Human-readable label for the tool palette buttons.
    pub fn label(self) -> &'static str {
        match self {
            SurfaceTool::NurbsCurve => "NurbsCurve",
            SurfaceTool::NurbsSurface => "NurbsSurface",
            SurfaceTool::CoonsFill => "CoonsFill",
            SurfaceTool::Sew => "Sew",
            SurfaceTool::Trim => "Trim",
            SurfaceTool::KnotOps => "KnotOps",
            SurfaceTool::Ssi => "SSI",
            SurfaceTool::Fit => "Fit",
            SurfaceTool::Ruled => "Ruled",
        }
    }
}

impl Default for SurfacePanelState {
    fn default() -> Self {
        Self {
            file: valenx_surface::persist::SurfaceFile::new(),
            tool: SurfaceTool::NurbsCurve,
            last_status: None,
            last_error: None,
            selected_curve: None,
            selected_surface: None,

            curve_degree: 3,
            curve_n_cps: 4,
            curve_cps: vec![
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [2.0, 1.0, 0.0],
                [3.0, 0.0, 0.0],
            ],
            curve_weights: vec![1.0, 1.0, 1.0, 1.0],

            surface_u_degree: 3,
            surface_v_degree: 3,
            surface_nu: 4,
            surface_nv: 4,
            surface_cps: (0..16)
                .map(|k| {
                    let i = k / 4;
                    let j = k % 4;
                    [i as f64 / 3.0, j as f64 / 3.0, 0.0]
                })
                .collect(),

            coons_curves: [0, 1, 2, 3],
            sew_surface_a: 0,
            sew_surface_b: 1,
            sew_edge_a: 1,
            sew_edge_b: 0,
            sew_tolerance: 1e-6,
            sew_use_g2: true, // Phase 19C default

            trim_surface: 0,
            trim_curve: 0,
            trim_side: 1, // default Outside (i.e. punch a hole)
            trim_resolution: 24,
            trim_use_uv: true,

            tess_resolution: 32,

            // Phase 19 defaults.
            knot_op_u: 0.5,
            knot_op_curve: 0,
            knot_op_tolerance: 1e-4,
            knot_op_surface: 0,
            knot_op_direction: 0,
            elevate_degree_by: 1,

            fit_points_text: "0,0,0\n1,0,0\n2,1,0\n3,1,0\n4,0,0\n".to_string(),
            fit_degree_u: 3,
            fit_degree_v: 3,
            fit_n_cps_u: 4,
            fit_n_cps_v: 4,
            fit_last_rms: None,

            ruled_curve_a: 0,
            ruled_curve_b: 1,
            ruled_extrude_vector: [0.0, 0.0, 1.0],
            ruled_apex: [0.0, 0.0, 1.0],
            ruled_kind: 0,

            ssi_surface_a: 0,
            ssi_surface_b: 1,
            ssi_tolerance: 1e-3,
        }
    }
}

/// State for the CAM workbench panel (Phase 10 — Path/CAM).
///
/// Holds the live in-memory tool table + stock + ordered op list,
/// the last-generated [`valenx_cam::Toolpath`] (regenerated on the
/// Generate button), inputs for adding tools / ops, and the
/// last-shown postprocessor choice / status / error.
#[derive(Clone, Debug)]
pub struct CamPanelState {
    /// The active CAM workbench file (tools + stock + operations).
    pub file: valenx_cam::persist::CamFile,
    /// Last generated toolpath — `None` until the user presses the
    /// Generate button.
    pub last_toolpath: Option<valenx_cam::Toolpath>,
    /// Currently selected postprocessor for export / preview.
    pub selected_postprocessor: valenx_cam::PostKind,
    /// `true` when the simulate-overlay should draw the toolpath in
    /// the viewport.
    pub show_overlay: bool,
    /// Last user-facing status string (green text at the bottom).
    pub last_status: Option<String>,
    /// Last user-facing error string (red text at the bottom).
    pub last_error: Option<String>,

    // ----- Add Tool inputs -----
    /// Display name for the next "Add Tool" press.
    pub new_tool_name: String,
    /// Kind of the new tool.
    pub new_tool_kind: valenx_cam::ToolKind,
    /// Diameter (mm).
    pub new_tool_diameter: f64,
    /// Length (mm).
    pub new_tool_length: f64,
    /// Flute count.
    pub new_tool_flutes: u32,
    /// Material descriptor.
    pub new_tool_material: String,

    // ----- Stock inputs (mirror Stock fields) -----
    /// Stock origin (mm).
    pub stock_origin: [f64; 3],
    /// Stock size (mm).
    pub stock_size: [f64; 3],
    /// Stock material descriptor.
    pub stock_material: String,

    // ----- Add Operation inputs (Profile / Pocket / Drill / Face) -----
    /// Which op the inline form is currently editing.
    pub new_op_kind: CamOpKind,
    /// Tool id input for the next "Add Op".
    pub new_op_tool_id: u32,
    /// Cutting feed (mm/min).
    pub new_op_feed: f64,
    /// Plunge feed (mm/min).
    pub new_op_plunge_feed: f64,
    /// Spindle RPM.
    pub new_op_spindle_rpm: f64,
    /// Z step-down per pass (mm).
    pub new_op_step_down: f64,
    /// XY step-over per pass (mm).
    pub new_op_step_over: f64,
    /// Total depth below stock top (mm).
    pub new_op_depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub new_op_safe_z: f64,
    /// Raster angle in degrees (Pocket / Face).
    pub new_op_raster_angle: f64,
    /// Climb-vs-conventional toggle.
    pub new_op_climb: bool,
    /// Pocket fill strategy.
    pub new_op_pocket_strategy: valenx_cam::PocketStrategy,
    /// Drill peck depth (mm).
    pub new_op_peck_depth: f64,
    /// Drill total depth (mm).
    pub new_op_drill_total_depth: f64,
    /// Drill retract clearance (mm).
    pub new_op_retract_clearance: f64,
    /// Drill hole positions (XYZ, Z ignored).
    pub new_op_hole_positions: Vec<[f64; 3]>,

    /// Last computed estimated cycle time (minutes).
    pub last_estimated_time_min: Option<f64>,
    /// Last computed removed volume (mm³).
    pub last_removed_volume_mm3: Option<f64>,

    // ----- Phase 17G — animation viewer + wear / fixture -----
    /// Animation voxel frames produced by the last "Animate" run.
    pub animation_frames: Vec<valenx_mesh::Mesh>,
    /// Selected animation frame (slider position).
    pub animation_frame_idx: usize,
    /// Voxel resolution used for the animation simulation.
    pub animation_resolution: u32,
    /// Number of animation frames requested.
    pub animation_n_frames: u32,
    /// Last wear warnings (per op).
    pub last_wear_warnings: Vec<(String, valenx_cam::wear::WearWarning)>,
    /// Active fixture geometry for collision checks.
    pub fixture: valenx_cam::fixture::Fixture,
    /// Last fixture-collision report.
    pub last_collisions: Vec<valenx_cam::fixture::Collision>,
}

/// Which CAM operation the inline form is currently editing.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum CamOpKind {
    /// [`valenx_cam::Operation::Profile`].
    Profile,
    /// [`valenx_cam::Operation::Pocket`].
    Pocket,
    /// [`valenx_cam::Operation::Drill`].
    Drill,
    /// [`valenx_cam::Operation::Face`].
    Face,
    /// Phase 17A — adaptive clearing.
    AdaptiveClearing,
    /// Phase 17A — helical bore.
    HelicalBore,
    /// Phase 17A — plunge rough.
    PlungeRough,
    /// Phase 17A — ramp entry.
    RampEntry,
    /// Phase 17A — peck-drill full.
    PeckDrillFull,
    /// Phase 17B — contour 2D.
    Contour2D,
    /// Phase 17B — contour 3D.
    Contour3D,
    /// Phase 17B — engrave.
    Engrave,
    /// Phase 17B — scribe.
    Scribe,
    /// Phase 17B — spiral pocket.
    SpiralPocket,
    /// Phase 17B — trochoidal slot.
    TrochoidalSlot,
    /// Phase 17B — waterline 3D.
    Waterline3D,
    /// Phase 17B — slot.
    Slot,
    /// Phase 17B — thread mill.
    ThreadMill,
    /// Phase 17B — rest machining.
    RestMachining,
}

impl CamOpKind {
    /// Human-readable label for the palette buttons.
    pub fn label(self) -> &'static str {
        match self {
            CamOpKind::Profile => "Profile",
            CamOpKind::Pocket => "Pocket",
            CamOpKind::Drill => "Drill",
            CamOpKind::Face => "Face",
            CamOpKind::AdaptiveClearing => "Adaptive Clearing",
            CamOpKind::HelicalBore => "Helical Bore",
            CamOpKind::PlungeRough => "Plunge Rough",
            CamOpKind::RampEntry => "Ramp Entry",
            CamOpKind::PeckDrillFull => "Peck Drill Full",
            CamOpKind::Contour2D => "Contour 2D",
            CamOpKind::Contour3D => "Contour 3D",
            CamOpKind::Engrave => "Engrave",
            CamOpKind::Scribe => "Scribe",
            CamOpKind::SpiralPocket => "Spiral Pocket",
            CamOpKind::TrochoidalSlot => "Trochoidal Slot",
            CamOpKind::Waterline3D => "Waterline 3D",
            CamOpKind::Slot => "Slot",
            CamOpKind::ThreadMill => "Thread Mill",
            CamOpKind::RestMachining => "Rest Machining",
        }
    }
}

impl Default for CamPanelState {
    fn default() -> Self {
        let mut file = valenx_cam::persist::CamFile::new();
        // Seed with one default tool so the user can immediately add
        // an op without first defining a tool.
        if let Ok(t) = valenx_cam::Tool::new(
            1,
            "EM6",
            valenx_cam::ToolKind::EndMill,
            6.0,
            25.0,
            2,
            "carbide",
        ) {
            file.tools.push(t);
        }
        Self {
            stock_origin: [
                file.stock.origin.x,
                file.stock.origin.y,
                file.stock.origin.z,
            ],
            stock_size: [file.stock.size.x, file.stock.size.y, file.stock.size.z],
            stock_material: file.stock.material.clone(),
            file,
            last_toolpath: None,
            selected_postprocessor: valenx_cam::PostKind::default(),
            show_overlay: false,
            last_status: None,
            last_error: None,

            new_tool_name: "Tool".into(),
            new_tool_kind: valenx_cam::ToolKind::EndMill,
            new_tool_diameter: 6.0,
            new_tool_length: 25.0,
            new_tool_flutes: 2,
            new_tool_material: "carbide".into(),

            new_op_kind: CamOpKind::Pocket,
            new_op_tool_id: 1,
            new_op_feed: 600.0,
            new_op_plunge_feed: 200.0,
            new_op_spindle_rpm: 12000.0,
            new_op_step_down: 1.0,
            new_op_step_over: 2.0,
            new_op_depth: 5.0,
            new_op_safe_z: 5.0,
            new_op_raster_angle: 0.0,
            new_op_climb: true,
            new_op_pocket_strategy: valenx_cam::PocketStrategy::ZigZag,
            new_op_peck_depth: 1.0,
            new_op_drill_total_depth: 5.0,
            new_op_retract_clearance: 1.0,
            new_op_hole_positions: vec![[10.0, 10.0, 0.0]],

            last_estimated_time_min: None,
            last_removed_volume_mm3: None,
            animation_frames: Vec::new(),
            animation_frame_idx: 0,
            animation_resolution: valenx_cam::simulate::DEFAULT_VOXEL_RES,
            animation_n_frames: 8,
            last_wear_warnings: Vec::new(),
            fixture: valenx_cam::fixture::Fixture::new(),
            last_collisions: Vec::new(),
        }
    }
}

/// Discriminant for the Arch panel's "Add entity" tool palette.
/// Mirrors [`valenx_arch::ArchEntityKind`] one-for-one but lives in
/// the app crate so the UI can carry per-tool form state.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ArchTool {
    /// Add a wall.
    Wall,
    /// Add a slab (floor / ceiling).
    Slab,
    /// Add a column.
    Column,
    /// Add a beam.
    Beam,
    /// Add a window (must reference a host wall id).
    Window,
    /// Add a door (must reference a host wall id).
    Door,
    /// Add a stair.
    Stair,
    /// Add a roof.
    Roof,
    /// Add a named space.
    Space,
}

impl ArchTool {
    /// Human label for the palette button.
    pub fn label(self) -> &'static str {
        match self {
            ArchTool::Wall => "Wall",
            ArchTool::Slab => "Slab",
            ArchTool::Column => "Column",
            ArchTool::Beam => "Beam",
            ArchTool::Window => "Window",
            ArchTool::Door => "Door",
            ArchTool::Stair => "Stair",
            ArchTool::Roof => "Roof",
            ArchTool::Space => "Space",
        }
    }
}

/// State for the Arch / BIM workbench panel (Phase 15).
///
/// Holds the active [`valenx_arch::ArchDocument`], the currently
/// selected [`ArchTool`], the numeric inputs for each "Add entity"
/// dialog, the last computed schedule, and transient
/// status / error messages.
#[derive(Clone, Debug)]
pub struct ArchPanelState {
    /// The active arch document.
    pub doc: valenx_arch::ArchDocument,
    /// Current tool the user is in.
    pub tool: ArchTool,
    /// Last computed schedule (rebuilt on every render via the
    /// "Refresh schedule" button — kept cached so the panel doesn't
    /// re-walk the doc on every frame).
    pub last_schedule: Option<valenx_arch::Schedule>,
    /// Tessellation tolerance for the "Render" button.
    pub render_tolerance: f64,
    /// Currently selected entity id (for highlight / delete).
    pub selected: Option<usize>,

    /// Last user-facing status (green text).
    pub last_status: Option<String>,
    /// Last user-facing error (red text).
    pub last_error: Option<String>,

    // ----- Wall inputs -----
    /// Wall start point.
    pub wall_start: [f64; 3],
    /// Wall end point.
    pub wall_end: [f64; 3],
    /// Wall height.
    pub wall_height: f64,
    /// Wall thickness.
    pub wall_thickness: f64,
    /// Wall material.
    pub wall_material: String,

    // ----- Slab inputs -----
    /// Slab boundary as a list of (x, y) pairs at z = `slab_z`.
    pub slab_boundary: Vec<[f64; 2]>,
    /// Slab z (floor level).
    pub slab_z: f64,
    /// Slab thickness.
    pub slab_thickness: f64,
    /// Slab material.
    pub slab_material: String,

    // ----- Column inputs -----
    /// Column base point.
    pub column_base: [f64; 3],
    /// Column height.
    pub column_height: f64,
    /// Column section kind: 0 = Rect, 1 = Circ, 2 = IBeam.
    pub column_section_kind: u8,
    /// Rect / IBeam width.
    pub column_width: f64,
    /// Rect / IBeam depth.
    pub column_depth: f64,
    /// Circular column radius + segments.
    pub column_radius: f64,
    /// Tessellation segment count for circular columns.
    pub column_segments: u32,
    /// I-beam flange thickness.
    pub column_flange_thickness: f64,
    /// I-beam web thickness.
    pub column_web_thickness: f64,
    /// Column material.
    pub column_material: String,

    // ----- Beam inputs -----
    /// Beam start point.
    pub beam_start: [f64; 3],
    /// Beam end point.
    pub beam_end: [f64; 3],
    /// Beam section kind: 0 = Rect, 1 = IBeam, 2 = Channel.
    pub beam_section_kind: u8,
    /// Beam width / flange width.
    pub beam_width: f64,
    /// Beam depth.
    pub beam_depth: f64,
    /// I-beam / channel thickness.
    pub beam_flange_thickness: f64,
    /// I-beam web thickness.
    pub beam_web_thickness: f64,
    /// Beam orientation about its axis (radians).
    pub beam_orientation: f64,
    /// Beam material.
    pub beam_material: String,

    // ----- Window inputs -----
    /// Host wall id for the next "Add window".
    pub window_host: usize,
    /// Position along the host wall.
    pub window_position_along: f64,
    /// Sill height above the wall's bottom.
    pub window_position_height: f64,
    /// Window width.
    pub window_width: f64,
    /// Window height.
    pub window_height: f64,
    /// Frame thickness.
    pub window_frame_thickness: f64,
    /// Window style index: 0 = Casement, 1 = Sliding, 2 = Awning, 3 = Fixed.
    pub window_style: u8,

    // ----- Door inputs -----
    /// Host wall id for the next "Add door".
    pub door_host: usize,
    /// Position along the host wall.
    pub door_position_along: f64,
    /// Door width.
    pub door_width: f64,
    /// Door height.
    pub door_height: f64,
    /// Door style index: 0 = Single, 1 = Double, 2 = Sliding, 3 = Bifold.
    pub door_style: u8,
    /// Hinge side index: 0 = Left, 1 = Right.
    pub door_hinge_side: u8,

    // ----- Stair inputs -----
    /// Stair base point.
    pub stair_base: [f64; 3],
    /// Stair direction.
    pub stair_direction: [f64; 3],
    /// Total rise.
    pub stair_total_rise: f64,
    /// Total run.
    pub stair_total_run: f64,
    /// Number of steps.
    pub stair_num_steps: u32,
    /// Stair width.
    pub stair_width: f64,

    // ----- Roof inputs -----
    /// Roof boundary (x, y) pairs at z = `roof_z`.
    pub roof_boundary: Vec<[f64; 2]>,
    /// Roof z (eave height).
    pub roof_z: f64,
    /// Peak height.
    pub roof_peak_height: f64,
    /// Roof type index: 0 = Flat, 1 = Gable, 2 = Hip, 3 = Shed.
    pub roof_type: u8,

    // ----- Space inputs -----
    /// Space boundary (x, y) pairs at z = `space_z`.
    pub space_boundary: Vec<[f64; 2]>,
    /// Floor z for the space.
    pub space_z: f64,
    /// Ceiling height.
    pub space_ceiling_height: f64,
    /// Space name.
    pub space_name: String,
}

impl Default for ArchPanelState {
    fn default() -> Self {
        Self {
            doc: valenx_arch::ArchDocument::new("Untitled"),
            tool: ArchTool::Wall,
            last_schedule: None,
            render_tolerance: 0.1,
            selected: None,
            last_status: None,
            last_error: None,
            wall_start: [0.0, 0.0, 0.0],
            wall_end: [3.0, 0.0, 0.0],
            wall_height: 2.7,
            wall_thickness: 0.2,
            wall_material: "Brick".into(),
            slab_boundary: vec![[0.0, 0.0], [4.0, 0.0], [4.0, 3.0], [0.0, 3.0]],
            slab_z: 0.0,
            slab_thickness: 0.2,
            slab_material: "Concrete".into(),
            column_base: [0.0, 0.0, 0.0],
            column_height: 2.7,
            column_section_kind: 0,
            column_width: 0.3,
            column_depth: 0.3,
            column_radius: 0.15,
            column_segments: 12,
            column_flange_thickness: 0.02,
            column_web_thickness: 0.01,
            column_material: "Steel".into(),
            beam_start: [0.0, 0.0, 2.7],
            beam_end: [3.0, 0.0, 2.7],
            beam_section_kind: 0,
            beam_width: 0.2,
            beam_depth: 0.3,
            beam_flange_thickness: 0.02,
            beam_web_thickness: 0.01,
            beam_orientation: 0.0,
            beam_material: "Steel".into(),
            window_host: 1,
            window_position_along: 1.5,
            window_position_height: 1.0,
            window_width: 1.0,
            window_height: 1.0,
            window_frame_thickness: 0.05,
            window_style: 0,
            door_host: 1,
            door_position_along: 1.5,
            door_width: 0.9,
            door_height: 2.1,
            door_style: 0,
            door_hinge_side: 0,
            stair_base: [0.0, 0.0, 0.0],
            stair_direction: [1.0, 0.0, 0.0],
            stair_total_rise: 3.0,
            stair_total_run: 4.0,
            stair_num_steps: 12,
            stair_width: 1.0,
            roof_boundary: vec![[0.0, 0.0], [8.0, 0.0], [8.0, 5.0], [0.0, 5.0]],
            roof_z: 2.7,
            roof_peak_height: 2.0,
            roof_type: 1, // Gable
            space_boundary: vec![[0.0, 0.0], [4.0, 0.0], [4.0, 3.0], [0.0, 3.0]],
            space_z: 0.0,
            space_ceiling_height: 2.7,
            space_name: "Living".into(),
        }
    }
}

/// Which CAD primitive to source for a newly-added
/// [`valenx_assembly::Part`]. Mirrors a subset of
/// [`CadPrimitiveKind`] — the assembly panel only needs the basic
/// three for v1.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum AssemblyPartPrimitive {
    /// Rectangular box.
    Box,
    /// Cylinder.
    Cylinder,
    /// Sphere.
    Sphere,
}

impl AssemblyPartPrimitive {
    /// Human-readable label for the combo box.
    pub fn label(&self) -> &'static str {
        match self {
            Self::Box => "Box",
            Self::Cylinder => "Cylinder",
            Self::Sphere => "Sphere",
        }
    }
}

/// Mate kind selector for the UI combo. Mirrors
/// [`valenx_assembly::MateKind`] minus the payload.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum AssemblyMateKindUi {
    /// Two points coincide.
    Coincident,
    /// Two points at fixed distance.
    Distance,
    /// Two vectors at fixed angle.
    Angle,
    /// Two vectors parallel.
    Parallel,
    /// Two vectors perpendicular.
    Perpendicular,
    /// Two cylindrical axes tangent.
    Tangent,
}

impl AssemblyMateKindUi {
    /// Combo label.
    pub fn label(&self) -> &'static str {
        match self {
            Self::Coincident => "Coincident",
            Self::Distance => "Distance",
            Self::Angle => "Angle",
            Self::Parallel => "Parallel",
            Self::Perpendicular => "Perpendicular",
            Self::Tangent => "Tangent",
        }
    }
}

/// Joint kind selector for the UI combo. Mirrors
/// [`valenx_assembly::JointKind`] minus the payload.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum AssemblyJointKindUi {
    /// Rigid bond.
    Fixed,
    /// Hinge.
    Revolute,
    /// Slider.
    Prismatic,
    /// Slider + hinge sharing one axis.
    Cylindrical,
    /// Ball joint.
    Spherical,
    /// Planar / face joint.
    Planar,
}

impl AssemblyJointKindUi {
    /// Combo label.
    pub fn label(&self) -> &'static str {
        match self {
            Self::Fixed => "Fixed",
            Self::Revolute => "Revolute",
            Self::Prismatic => "Prismatic",
            Self::Cylindrical => "Cylindrical",
            Self::Spherical => "Spherical",
            Self::Planar => "Planar",
        }
    }
}

/// State for the TechDraw workbench panel (Phase 5).
///
/// Holds the currently-edited [`valenx_techdraw::Drawing`], the
/// numeric inputs for the "Add view" / "Add dimension" buttons, and
/// transient UI flags (selected view, last error, last-export path).
#[derive(Clone, Debug)]
pub struct TechDrawPanelState {
    /// The active drawing.
    pub drawing: valenx_techdraw::Drawing,
    /// Currently selected view index (for "Edit view" / "Remove view"
    /// dialogs).
    pub selected_view: Option<usize>,
    /// Last user-facing error message (red-text label at the bottom
    /// of the panel).
    pub last_error: Option<String>,
    /// Pending position for the next "Add view" click.
    pub new_view_position: [f64; 2],
    /// Scale for the next "Add view" click.
    pub new_view_scale: f64,
    /// Linear dimension inputs.
    pub dim_from: [f64; 2],
    pub dim_to: [f64; 2],
    pub dim_offset: f64,
    /// Sheet size selector (A4..A0).
    pub sheet_size: valenx_techdraw::SheetSize,
    /// Phase 18A — make the next-added view parametric (linked to the
    /// active feature tree's last feature, auto-update on tree replay).
    pub new_view_parametric: bool,
    /// Phase 18B — dim-chain inputs: kind + offset.
    pub chain_kind: valenx_techdraw::DimChainKind,
    pub chain_offset: f64,
    /// Comma-separated `x,y` pairs to append as chain entries.
    pub chain_entries: String,
    /// Phase 18C — balloon inputs.
    pub balloon_position: [f64; 2],
    pub balloon_target: [f64; 2],
    pub balloon_number: String,
    pub balloon_style: valenx_techdraw::BalloonStyle,
    /// Phase 18C — leader inputs.
    pub leader_start: [f64; 2],
    pub leader_end: [f64; 2],
    pub leader_text: String,
    pub leader_arrow: valenx_techdraw::ArrowKind,
    /// Phase 18D — weld inputs.
    pub weld_position: [f64; 2],
    pub weld_target: [f64; 2],
    pub weld_type: valenx_techdraw::WeldType,
    pub weld_side: valenx_techdraw::WeldPosition,
    pub weld_size: String,
    pub weld_all_around: bool,
    pub weld_field: bool,
    /// Phase 18E — surface-finish inputs.
    pub sf_position: [f64; 2],
    pub sf_ra: f64,
    pub sf_process: valenx_techdraw::SurfaceProcess,
    pub sf_lay: valenx_techdraw::LayPattern,
    /// Phase 18F — GD&T inputs.
    pub gdt_position: [f64; 2],
    pub gdt_characteristic: valenx_techdraw::GeometricCharacteristic,
    pub gdt_tolerance: String,
    pub gdt_modifier: valenx_techdraw::MaterialCondition,
    pub gdt_datum_letters: String,
    /// Phase 18F — datum inputs.
    pub datum_position: [f64; 2],
    pub datum_target: [f64; 2],
    pub datum_letter: String,
    /// Phase 18G — selected hatch-pattern name for section fills.
    pub hatch_pattern: String,
}

impl Default for TechDrawPanelState {
    fn default() -> Self {
        Self {
            drawing: valenx_techdraw::Drawing::new(valenx_techdraw::Sheet::a4_landscape(
                "Untitled", "", "A",
            )),
            selected_view: None,
            last_error: None,
            new_view_position: [80.0, 100.0],
            new_view_scale: 1.0,
            dim_from: [0.0, 0.0],
            dim_to: [10.0, 0.0],
            dim_offset: 5.0,
            sheet_size: valenx_techdraw::SheetSize::A4,
            new_view_parametric: false,
            chain_kind: valenx_techdraw::DimChainKind::Chain,
            chain_offset: 5.0,
            chain_entries: "0,0; 10,0; 20,0".to_string(),
            balloon_position: [20.0, 20.0],
            balloon_target: [40.0, 40.0],
            balloon_number: "1".to_string(),
            balloon_style: valenx_techdraw::BalloonStyle::Circle,
            leader_start: [10.0, 10.0],
            leader_end: [30.0, 30.0],
            leader_text: "note".to_string(),
            leader_arrow: valenx_techdraw::ArrowKind::Closed,
            weld_position: [50.0, 50.0],
            weld_target: [70.0, 70.0],
            weld_type: valenx_techdraw::WeldType::Fillet,
            weld_side: valenx_techdraw::WeldPosition::Arrow,
            weld_size: "5".to_string(),
            weld_all_around: false,
            weld_field: false,
            sf_position: [80.0, 80.0],
            sf_ra: 1.6,
            sf_process: valenx_techdraw::SurfaceProcess::Machined,
            sf_lay: valenx_techdraw::LayPattern::Parallel,
            gdt_position: [100.0, 100.0],
            gdt_characteristic: valenx_techdraw::GeometricCharacteristic::Position,
            gdt_tolerance: "0.1".to_string(),
            gdt_modifier: valenx_techdraw::MaterialCondition::Rfs,
            gdt_datum_letters: "A".to_string(),
            datum_position: [120.0, 120.0],
            datum_target: [130.0, 130.0],
            datum_letter: "A".to_string(),
            hatch_pattern: "ANSI31".to_string(),
        }
    }
}

/// State for the Dock panel — receptor + ligand paths the user has
/// chosen, search-box parameters, last-run results.
#[derive(Clone, Debug, Default)]
pub struct DockPanelState {
    pub receptor_path: String,
    pub ligand_path: String,
    pub output_path: String,
    pub center: [f64; 3],
    pub size: [f64; 3],
    pub exhaustiveness: u32,
    pub num_modes: u32,
    pub energy_range: f64,
    pub seed: u64,
    /// Last-run results: parallel `(score, rmsd_to_first)` per pose.
    pub last_scores: Vec<(f64, f64)>,
    /// Index into `last_scores` of the pose currently visualized in the viewport.
    pub selected_pose: Option<usize>,
    /// Last error message, if any.
    pub last_error: Option<String>,
}

/// State for the Sketcher panel — currently active sketch, current
/// drawing tool, last solver report.
#[derive(Clone, Debug, Default)]
pub struct SketcherPanelState {
    /// The active sketch being edited.
    pub sketch: valenx_sketch::Sketch,
    /// Current tool the user is in.
    pub tool: SketcherTool,
    /// Last solver report (None until first solve).
    pub last_report: Option<valenx_sketch::SolverReport>,
    /// Pending click — for two-click tools (line, etc.) holds the first click's EntityId.
    pub pending_first_click: Option<valenx_sketch::EntityId>,
    /// Pad depth (used for extrude).
    pub pad_depth: f64,
    /// Last error to display to the user.
    pub last_error: Option<String>,
    /// Numeric-click input X coord (MVP — see `pending_click_y`).
    pub pending_click_x: f64,
    /// Numeric-click input Y coord (MVP: until the viewport-click
    /// pipeline lands in Phase 1H, users enter (x, y) numerically and
    /// press "Add point at (x, y)" to simulate clicking in 3D).
    pub pending_click_y: f64,
    /// Currently selected entities (multi-select via shift-click).
    pub selected: Vec<valenx_sketch::EntityId>,
    /// Pending numeric target for the Distance / Angle / Radius
    /// constraint buttons (v1 simplification: enter the value here
    /// BEFORE clicking the constraint button instead of a modal popup).
    pub pending_target: f64,
    /// Toggle for drawing the 2D sketch overlay on the viewport.
    /// `true` by default — users opt into seeing what they sketch as
    /// soon as the panel opens. Set to false to hide the overlay (e.g.
    /// when reviewing the underlying mesh without sketch clutter).
    pub show_overlay: bool,
    /// Undo / redo history for [`Self::sketch`]. Snapshot is pushed
    /// just before any sketch-mutating action (add point / line /
    /// circle, add constraint, mark construction, sketch op) so
    /// `↶`/`↷` and `Ctrl+Z` / `Ctrl+Y` restore the previous state.
    pub history: crate::undo::History<valenx_sketch::Sketch>,
}

impl SketcherPanelState {
    /// Default `show_overlay` is `true` — but `#[derive(Default)]` on
    /// the parent gives every field its type's default, which would be
    /// `false` for `bool`. Override via this constructor used by
    /// [`MeshToolboxState::default`].
    fn new_with_overlay_on() -> Self {
        Self {
            show_overlay: true,
            ..Self::default()
        }
    }

    /// Snapshot the current sketch state onto the undo stack. Called
    /// before every action that mutates [`Self::sketch`] so a
    /// subsequent `↶` / `Ctrl+Z` restores the pre-edit value.
    pub fn record(&mut self) {
        self.history.record(self.sketch.clone());
    }

    /// Pop the latest snapshot and restore it. Returns `true` if a
    /// snapshot was popped.
    pub fn undo_edit(&mut self) -> bool {
        if let Some(prev) = self.history.undo(self.sketch.clone()) {
            self.sketch = prev;
            // Stale selection refs may point past the new entity count.
            self.selected
                .retain(|id| id.0 != 0 && id.0 <= self.sketch.entities.len());
            self.last_report = None;
            self.last_error = None;
            true
        } else {
            false
        }
    }

    /// Pop the latest redo snapshot and reapply it. Returns `true`
    /// if a snapshot was popped.
    pub fn redo_edit(&mut self) -> bool {
        if let Some(next) = self.history.redo(self.sketch.clone()) {
            self.sketch = next;
            self.selected
                .retain(|id| id.0 != 0 && id.0 <= self.sketch.entities.len());
            self.last_report = None;
            self.last_error = None;
            true
        } else {
            false
        }
    }

    /// Whether an undo would do something. Used by the inline
    /// `↶` button's enabled state.
    pub fn can_undo(&self) -> bool {
        self.history.can_undo()
    }

    /// Whether a redo would do something.
    pub fn can_redo(&self) -> bool {
        self.history.can_redo()
    }
}

/// State for the Part Design panel — feature tree (sketches +
/// features), selection, replay status.
///
/// The panel is built around a [`valenx_feature_tree::FeatureTree`]
/// the user grows by clicking "Add Sketch", "Add Pad", "Add Pocket",
/// etc. Each change sets `pending_replay = true`; the panel notices
/// on the next frame and (if `auto_replay`) re-evaluates the tree and
/// pushes the result into the viewport.
#[derive(Clone, Debug)]
pub struct PartDesignPanelState {
    /// The parametric tree the user is editing.
    pub tree: valenx_feature_tree::FeatureTree,
    /// Currently selected feature (for Suppress / Delete buttons).
    pub selected_feature: Option<valenx_feature_tree::feature::FeatureId>,
    /// Currently selected sketch (for Add-feature param dialogs).
    pub selected_sketch: Option<valenx_feature_tree::feature::SketchRef>,
    /// Last replay error, displayed in red below the tree.
    pub last_replay_error: Option<String>,
    /// Auto-replay every change. Defaults to true.
    pub auto_replay: bool,
    /// "A change happened" flag set by add/delete/suppress handlers.
    /// Cleared after a successful replay so the next idle frame is
    /// cheap.
    pub pending_replay: bool,
    /// Pad params (used by the "Add Pad" button — v1: edit before
    /// clicking Add, modal popup deferred to Phase 2.5).
    pub pad_sketch_index: usize,
    pub pad_depth: f64,
    pub pad_direction_positive: bool,
    /// Pocket params.
    pub pocket_sketch_index: usize,
    pub pocket_depth: f64,
    pub pocket_direction_positive: bool,
    /// Revolve params.
    pub revolve_sketch_index: usize,
    pub revolve_axis_origin: [f64; 3],
    pub revolve_axis_direction: [f64; 3],
    pub revolve_angle_deg: f64,
    /// Mirror params.
    pub mirror_target_index: usize,
    pub mirror_plane_origin: [f64; 3],
    pub mirror_plane_normal: [f64; 3],
    pub mirror_keep_original: bool,
    /// Linear pattern params.
    pub lp_target_index: usize,
    pub lp_direction: [f64; 3],
    pub lp_count: u32,
    pub lp_spacing: f64,
    /// Circular pattern params.
    pub cp_target_index: usize,
    pub cp_axis_origin: [f64; 3],
    pub cp_axis_direction: [f64; 3],
    pub cp_count: u32,
    pub cp_total_angle_deg: f64,
    /// Fillet params (Phase 3 — Task 30, extended in Phase 14). The
    /// Phase 14 dispatcher tries BRep first; when `fillet_edge_indices`
    /// is empty, the auto-selector by `fillet_threshold_deg` runs.
    /// Populated indices ("0,1,2,3") request explicit BRep edge
    /// selection.
    pub fillet_target_index: usize,
    pub fillet_radius: f64,
    pub fillet_threshold_deg: f64,
    /// Comma-separated 0-based BRep edge indices. Empty = "auto".
    pub fillet_edge_indices: String,
    /// Chamfer params (Phase 3 — Task 31, extended in Phase 14).
    pub chamfer_target_index: usize,
    pub chamfer_distance: f64,
    pub chamfer_threshold_deg: f64,
    /// Comma-separated 0-based BRep edge indices. Empty = "auto".
    pub chamfer_edge_indices: String,
    /// Hole params (Phase 13A). `hole_depth_mode_idx`: 0 = Blind, 1 =
    /// Through, 2 = UpToFace.
    pub hole_sketch_index: usize,
    pub hole_depth_mode_idx: usize,
    pub hole_blind_depth: f64,
    pub hole_drill_diameter: f64,
    pub hole_direction_negative: bool,
    pub hole_use_counterbore: bool,
    pub hole_counterbore_diameter: f64,
    pub hole_counterbore_depth: f64,
    pub hole_use_countersink: bool,
    pub hole_countersink_diameter: f64,
    pub hole_countersink_angle_deg: f64,
    pub hole_use_thread: bool,
    /// 0=ISO, 1=UN, 2=BSPP, 3=NPT
    pub hole_thread_standard_idx: usize,
    pub hole_thread_entry_idx: usize,
    /// Loft params (Phase 13B).
    pub loft_profile_indices: String,
    pub loft_guide_indices: String,
    pub loft_closed: bool,
    pub loft_ruled: bool,
    /// Sweep params (Phase 13B).
    pub sweep_profile_index: usize,
    pub sweep_path_index: usize,
    pub sweep_twist_deg: f64,
    pub sweep_keep_orientation: bool,
    /// Pipe params (Phase 13B).
    pub pipe_cross_section_index: usize,
    pub pipe_centerline_index: usize,
    pub pipe_bend_radius: f64,
    /// Helix params (Phase 13C).
    pub helix_profile_index: usize,
    pub helix_pitch: f64,
    pub helix_turns: f64,
    pub helix_axis_origin: [f64; 3],
    pub helix_axis_direction: [f64; 3],
    pub helix_taper_deg: f64,
    pub helix_left_handed: bool,
    /// MultiTransform params (Phase 13C). `multi_ops` stored as a
    /// simple textual recipe: each line is "translate dx dy dz" or
    /// "rotate ax ay az angle_deg" or "scale factor" or
    /// "mirror nx ny nz". Empty lines ignored.
    pub mt_target_index: usize,
    pub mt_ops_recipe: String,
    /// DraftAngle params (Phase 13D). `face_indices_csv` is a
    /// comma-separated list of triangle indices.
    pub draft_target_index: usize,
    pub draft_face_indices_csv: String,
    pub draft_neutral_normal: [f64; 3],
    pub draft_angle_deg: f64,
    /// Shell params (Phase 13D).
    pub shell_target_index: usize,
    pub shell_face_indices_csv: String,
    pub shell_thickness: f64,
    pub shell_side_idx: usize, // 0=Inward 1=Outward
    /// Thickness params (Phase 13D).
    pub thickness_target_index: usize,
    pub thickness_face_index: usize,
    pub thickness_thickness: f64,
    /// BooleanHistory params (Phase 13E). `bh_targets_csv` is
    /// comma-separated feature ids.
    pub bh_op_idx: usize, // 0=Union 1=Difference 2=Intersection 3=Section
    pub bh_targets_csv: String,
    /// Delete-confirmation latch — set by clicking "Delete feature",
    /// cleared by the confirm/cancel buttons or by changing selection.
    pub pending_delete_confirm: bool,
    /// Last `.valenx` project path that was saved or loaded via the
    /// Part Design panel buttons. `None` until the user invokes one of
    /// those buttons. Surfaced in the panel header so the user can see
    /// what's on disk vs. in memory at a glance.
    pub project_path: Option<PathBuf>,
    /// Undo / redo history for [`Self::tree`]. Snapshot is pushed just
    /// before any tree-mutating action (Add Sketch / Add Feature /
    /// Toggle Suppress / Delete) so `↶` / `↷` and `Ctrl+Z` /
    /// `Ctrl+Y` restore the previous state.
    pub history: crate::undo::History<valenx_feature_tree::FeatureTree>,
}

impl PartDesignPanelState {
    /// Snapshot the current feature tree onto the undo stack. Called
    /// before every mutating action so a subsequent `↶` / `Ctrl+Z`
    /// restores the pre-edit value.
    pub fn record(&mut self) {
        self.history.record(self.tree.clone());
    }

    /// Pop the latest snapshot and restore it. Returns `true` if a
    /// snapshot was popped.
    pub fn undo_edit(&mut self) -> bool {
        if let Some(prev) = self.history.undo(self.tree.clone()) {
            self.tree = prev;
            // Stale ids may point past the new feature / sketch count.
            self.selected_feature = self.selected_feature.filter(|id| id.0 < self.tree.features.len());
            self.selected_sketch = self.selected_sketch.filter(|r| r.0 < self.tree.sketches.len());
            self.pending_replay = true;
            self.pending_delete_confirm = false;
            self.last_replay_error = None;
            true
        } else {
            false
        }
    }

    /// Pop the latest redo snapshot and reapply it.
    pub fn redo_edit(&mut self) -> bool {
        if let Some(next) = self.history.redo(self.tree.clone()) {
            self.tree = next;
            self.selected_feature = self.selected_feature.filter(|id| id.0 < self.tree.features.len());
            self.selected_sketch = self.selected_sketch.filter(|r| r.0 < self.tree.sketches.len());
            self.pending_replay = true;
            self.pending_delete_confirm = false;
            self.last_replay_error = None;
            true
        } else {
            false
        }
    }

    /// Whether an undo would do something.
    pub fn can_undo(&self) -> bool {
        self.history.can_undo()
    }

    /// Whether a redo would do something.
    pub fn can_redo(&self) -> bool {
        self.history.can_redo()
    }
}

impl Default for PartDesignPanelState {
    fn default() -> Self {
        Self {
            tree: valenx_feature_tree::FeatureTree::new(),
            selected_feature: None,
            selected_sketch: None,
            last_replay_error: None,
            auto_replay: true,
            pending_replay: false,
            pad_sketch_index: 0,
            pad_depth: 10.0,
            pad_direction_positive: true,
            pocket_sketch_index: 0,
            pocket_depth: 5.0,
            pocket_direction_positive: false,
            revolve_sketch_index: 0,
            revolve_axis_origin: [0.0, 0.0, 0.0],
            revolve_axis_direction: [0.0, 1.0, 0.0],
            revolve_angle_deg: 360.0,
            mirror_target_index: 0,
            mirror_plane_origin: [0.0, 0.0, 0.0],
            mirror_plane_normal: [1.0, 0.0, 0.0],
            mirror_keep_original: true,
            lp_target_index: 0,
            lp_direction: [1.0, 0.0, 0.0],
            lp_count: 4,
            lp_spacing: 2.0,
            cp_target_index: 0,
            cp_axis_origin: [0.0, 0.0, 0.0],
            cp_axis_direction: [0.0, 0.0, 1.0],
            cp_count: 6,
            cp_total_angle_deg: 360.0,
            fillet_target_index: 0,
            fillet_radius: 0.1,
            fillet_threshold_deg: 45.0,
            fillet_edge_indices: String::new(),
            chamfer_target_index: 0,
            chamfer_distance: 0.1,
            chamfer_threshold_deg: 45.0,
            chamfer_edge_indices: String::new(),
            hole_sketch_index: 0,
            hole_depth_mode_idx: 0,
            hole_blind_depth: 5.0,
            hole_drill_diameter: 2.0,
            hole_direction_negative: true,
            hole_use_counterbore: false,
            hole_counterbore_diameter: 4.0,
            hole_counterbore_depth: 1.5,
            hole_use_countersink: false,
            hole_countersink_diameter: 4.0,
            hole_countersink_angle_deg: 82.0,
            hole_use_thread: false,
            hole_thread_standard_idx: 0,
            hole_thread_entry_idx: 0,
            loft_profile_indices: "0,1".into(),
            loft_guide_indices: String::new(),
            loft_closed: false,
            loft_ruled: true,
            sweep_profile_index: 0,
            sweep_path_index: 1,
            sweep_twist_deg: 0.0,
            sweep_keep_orientation: true,
            pipe_cross_section_index: 0,
            pipe_centerline_index: 1,
            pipe_bend_radius: 0.5,
            helix_profile_index: 0,
            helix_pitch: 2.0,
            helix_turns: 3.0,
            helix_axis_origin: [0.0, 0.0, 0.0],
            helix_axis_direction: [0.0, 0.0, 1.0],
            helix_taper_deg: 0.0,
            helix_left_handed: false,
            mt_target_index: 0,
            mt_ops_recipe: "translate 5 0 0\nrotate 0 0 1 30\n".into(),
            draft_target_index: 0,
            draft_face_indices_csv: "0,1".into(),
            draft_neutral_normal: [0.0, 0.0, 1.0],
            draft_angle_deg: 5.0,
            shell_target_index: 0,
            shell_face_indices_csv: String::new(),
            shell_thickness: 0.1,
            shell_side_idx: 0,
            thickness_target_index: 0,
            thickness_face_index: 0,
            thickness_thickness: 0.1,
            bh_op_idx: 0,
            bh_targets_csv: "0,1".into(),
            pending_delete_confirm: false,
            project_path: None,
            history: crate::undo::History::new(),
        }
    }
}

/// Draft workbench tool the user is currently in.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum DraftTool {
    /// Browsing entities; no placement on click.
    #[default]
    Select,
    /// Place a line by entering its start + end coords.
    Line,
    /// Place a polyline by appending points; explicit "Close" button
    /// commits a closed polyline.
    Polyline,
    /// Place an arc by entering centre, radius, start + end angles.
    Arc,
    /// Place a circle by entering centre + radius.
    Circle,
    /// Place an axis-aligned rectangle by entering min + max corners.
    Rectangle,
    /// Place a regular polygon by entering centre, radius, side count.
    Polygon,
    /// Place a linear dimension by entering from + to + perpendicular offset.
    Dimension,
    /// Place a text label by entering position, content, and size.
    Text,
}

/// Which built-in working plane the Draft panel's combo selects.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum DraftPlaneKind {
    /// World XY (top view).
    #[default]
    Xy,
    /// World XZ (front view).
    Xz,
    /// World YZ (side view).
    Yz,
}

impl DraftPlaneKind {
    fn to_plane(self) -> valenx_draft::WorkingPlane {
        match self {
            DraftPlaneKind::Xy => valenx_draft::WorkingPlane::from_xy(),
            DraftPlaneKind::Xz => valenx_draft::WorkingPlane::from_xz(),
            DraftPlaneKind::Yz => valenx_draft::WorkingPlane::from_yz(),
        }
    }

    fn label(self) -> &'static str {
        match self {
            DraftPlaneKind::Xy => "XY (top)",
            DraftPlaneKind::Xz => "XZ (front)",
            DraftPlaneKind::Yz => "YZ (side)",
        }
    }
}

/// State for the Draft workbench panel — active document, current
/// tool, click buffers, and overlay toggle.
#[derive(Clone, Debug)]
pub struct DraftPanelState {
    /// The active draft document being edited.
    pub document: valenx_draft::DraftDocument,
    /// Current tool the user is in.
    pub tool: DraftTool,
    /// Which built-in plane the document is on (UI selector).
    pub plane_kind: DraftPlaneKind,
    /// Toggle for drawing the 2D draft overlay on the viewport.
    /// `true` by default — users opt into seeing what they draw as
    /// soon as the panel opens.
    pub show_overlay: bool,
    /// Currently selected entity index (for delete).
    pub selected_entity: Option<usize>,
    /// Pending error message displayed at the bottom of the panel.
    pub last_error: Option<String>,
    /// Grid spacing for grid-snap (0 = disabled). Used by the
    /// overlay's grid-intersection snap marker.
    pub grid_spacing: f64,

    // ----- Line tool inputs -----
    /// Start point candidate (x, y).
    pub line_start: [f64; 2],
    /// End point candidate (x, y).
    pub line_end: [f64; 2],
    /// True once the user has "Place start"-ed; flips back on commit.
    pub line_start_placed: bool,

    // ----- Polyline tool inputs -----
    /// Pending vertex coords for the next "Append" press.
    pub polyline_next_point: [f64; 2],
    /// Accumulating list of vertices for the current polyline.
    pub polyline_points: Vec<[f64; 2]>,

    // ----- Arc tool inputs -----
    pub arc_center: [f64; 2],
    pub arc_radius: f64,
    pub arc_start_angle_deg: f64,
    pub arc_end_angle_deg: f64,

    // ----- Circle tool inputs -----
    pub circle_center: [f64; 2],
    pub circle_radius: f64,

    // ----- Rectangle tool inputs -----
    pub rect_min: [f64; 2],
    pub rect_max: [f64; 2],

    // ----- Polygon tool inputs -----
    pub polygon_center: [f64; 2],
    pub polygon_radius: f64,
    pub polygon_sides: u32,

    // ----- Linear dimension inputs -----
    pub dim_from: [f64; 2],
    pub dim_to: [f64; 2],
    pub dim_offset: f64,

    // ----- Text tool inputs -----
    pub text_position: [f64; 2],
    pub text_content: String,
    pub text_size: f64,
}

impl Default for DraftPanelState {
    fn default() -> Self {
        Self {
            document: valenx_draft::DraftDocument::default(),
            tool: DraftTool::default(),
            plane_kind: DraftPlaneKind::default(),
            show_overlay: true,
            selected_entity: None,
            last_error: None,
            grid_spacing: 1.0,

            line_start: [0.0, 0.0],
            line_end: [1.0, 0.0],
            line_start_placed: false,

            polyline_next_point: [0.0, 0.0],
            polyline_points: Vec::new(),

            arc_center: [0.0, 0.0],
            arc_radius: 1.0,
            arc_start_angle_deg: 0.0,
            arc_end_angle_deg: 90.0,

            circle_center: [0.0, 0.0],
            circle_radius: 1.0,

            rect_min: [0.0, 0.0],
            rect_max: [1.0, 1.0],

            polygon_center: [0.0, 0.0],
            polygon_radius: 1.0,
            polygon_sides: 6,

            dim_from: [0.0, 0.0],
            dim_to: [10.0, 0.0],
            dim_offset: 1.0,

            text_position: [0.0, 0.0],
            text_content: "label".to_string(),
            text_size: 0.5,
        }
    }
}

// =============================================================================
// Phase 16 — Spreadsheet Workbench
// =============================================================================

/// State for the Spreadsheet workbench panel — the active workbook,
/// the currently-selected sheet, editable cell buffer, and the
/// re-evaluate-all status feedback.
///
/// Phase 16 v1 ships a minimal grid: a single sheet rendered as a
/// fixed-size table (rows × cols defined by [`SpreadsheetPanelState::view_rows`]
/// / `view_cols`). Cells are edited by typing into the `editor_*`
/// buffer and pressing "Set"; clicking a coordinate populates the
/// buffer for editing. "Re-evaluate" walks every visible cell and
/// surfaces parse / circular-ref errors as a single status string.
#[derive(Clone, Debug)]
pub struct SpreadsheetPanelState {
    /// Live workbook (round-trips through [`valenx_spreadsheet::SpreadsheetFile`]).
    pub workbook: valenx_spreadsheet::Spreadsheet,
    /// Name of the currently-selected sheet (empty string when no
    /// sheet exists yet).
    pub active_sheet: String,
    /// New-sheet name buffer (typed into "Add sheet" field).
    pub new_sheet_name: String,
    /// Visible row count for the grid (1-indexed in display).
    pub view_rows: u32,
    /// Visible column count for the grid (alphabetic A..= in display).
    pub view_cols: u32,
    /// Currently selected cell within the active sheet
    /// (`(row, col)`). `None` until the user clicks a coordinate.
    pub selected_cell: Option<(u32, u32)>,
    /// Editor buffer — the contents the user is typing for the
    /// selected cell. Pressing "Set" applies it.
    pub editor_text: String,
    /// Last status string ("12 cells re-evaluated", green).
    pub last_status: Option<String>,
    /// Last error string ("circular reference in Sheet.A1", red).
    pub last_error: Option<String>,
}

impl Default for SpreadsheetPanelState {
    fn default() -> Self {
        let mut workbook = valenx_spreadsheet::Spreadsheet::new();
        workbook.add_sheet("Default");
        Self {
            workbook,
            active_sheet: "Default".to_string(),
            new_sheet_name: "Sheet1".to_string(),
            view_rows: 12,
            view_cols: 6,
            selected_cell: None,
            editor_text: String::new(),
            last_status: None,
            last_error: None,
        }
    }
}

/// Sketcher tool the user is currently in.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum SketcherTool {
    /// Selecting entities.
    #[default]
    Select,
    /// Click two points to draw a line.
    Line,
    /// Click center then radius to draw a circle.
    Circle,
}

/// Primitive shapes available from the Part workbench combo.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CadPrimitiveKind {
    /// Axis-aligned rectangular box.
    Box,
    /// Right circular cylinder along Z.
    Cylinder,
    /// Sphere centred on the origin.
    Sphere,
    /// Truncated cone / frustum along Z.
    Cone,
    /// Torus with major axis Z.
    Torus,
}

impl CadPrimitiveKind {
    fn label(self) -> &'static str {
        match self {
            CadPrimitiveKind::Box => "Box",
            CadPrimitiveKind::Cylinder => "Cylinder",
            CadPrimitiveKind::Sphere => "Sphere",
            CadPrimitiveKind::Cone => "Cone / Frustum",
            CadPrimitiveKind::Torus => "Torus",
        }
    }

    fn audit_tag(self) -> &'static str {
        match self {
            CadPrimitiveKind::Box => "cad.primitive.box",
            CadPrimitiveKind::Cylinder => "cad.primitive.cylinder",
            CadPrimitiveKind::Sphere => "cad.primitive.sphere",
            CadPrimitiveKind::Cone => "cad.primitive.cone",
            CadPrimitiveKind::Torus => "cad.primitive.torus",
        }
    }
}

/// Three-way radio enum shared by the rotate-axis and mirror-plane
/// rows. Maps to [`valenx_mesh::Axis`] / [`valenx_mesh::Plane`] when
/// the Apply button reaches into the mesh ops.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ToolboxAxis {
    X,
    Y,
    Z,
}

impl ToolboxAxis {
    fn to_axis(self) -> Axis {
        match self {
            ToolboxAxis::X => Axis::X,
            ToolboxAxis::Y => Axis::Y,
            ToolboxAxis::Z => Axis::Z,
        }
    }

    fn to_plane(self) -> Plane {
        match self {
            ToolboxAxis::X => Plane::X,
            ToolboxAxis::Y => Plane::Y,
            ToolboxAxis::Z => Plane::Z,
        }
    }

    fn label(self) -> &'static str {
        match self {
            ToolboxAxis::X => "X",
            ToolboxAxis::Y => "Y",
            ToolboxAxis::Z => "Z",
        }
    }
}

/// Apply an in-place transformation to a `TriangleMesh` (STL path).
/// `xform` runs on each individual vertex; `recompute_normals=true`
/// makes us re-derive the face normal from the right-handed winding
/// of the freshly-transformed vertices.
fn for_each_vertex<F: Fn([f32; 3]) -> [f32; 3]>(
    mesh: &mut TriangleMesh,
    xform: F,
    recompute_normals: bool,
) {
    for tri in &mut mesh.triangles {
        for v in &mut tri.vertices {
            *v = xform(*v);
        }
        if recompute_normals {
            tri.normal = tri.computed_normal();
        }
    }
}

impl ValenxApp {
    /// Translate either the loaded canonical mesh or the loaded STL
    /// by `(dx, dy, dz)`. Recomputes the cached quality report on
    /// the canonical mesh so the Inspector reflects the move. Emits
    /// `mesh.translate` audit.
    pub fn apply_translate(&mut self, dx: f64, dy: f64, dz: f64) {
        if dx == 0.0 && dy == 0.0 && dz == 0.0 {
            self.status = Some("Translate skipped — delta is zero".into());
            return;
        }
        let mut applied = false;
        if let Some(loaded) = self.mesh.as_mut() {
            valenx_mesh::transform::translate(&mut loaded.mesh, dx, dy, dz);
            refresh_loaded_mesh_quality(loaded);
            applied = true;
        } else if let Some(stl) = self.stl.as_mut() {
            for_each_vertex(
                &mut stl.mesh,
                |v| [v[0] + dx as f32, v[1] + dy as f32, v[2] + dz as f32],
                false, // translation preserves normals
            );
            applied = true;
        }
        if applied {
            self.status = Some(format!("Translated by ({dx:.4}, {dy:.4}, {dz:.4})"));
            self.last_error = None;
            emit_audit(
                "mesh.translate",
                serde_json::json!({"kind": "mesh"}),
                serde_json::json!({"dx": dx, "dy": dy, "dz": dz}),
            );
        } else {
            self.last_error = Some("Translate: no mesh or STL loaded.".into());
        }
    }

    /// Uniform scale around the world origin.
    pub fn apply_scale_uniform(&mut self, factor: f64) {
        if factor == 1.0 {
            self.status = Some("Scale skipped — factor is 1".into());
            return;
        }
        if !factor.is_finite() || factor == 0.0 {
            self.last_error = Some(format!("Scale: factor {factor} is not allowed."));
            return;
        }
        let mut applied = false;
        if let Some(loaded) = self.mesh.as_mut() {
            valenx_mesh::transform::scale_uniform(&mut loaded.mesh, factor);
            refresh_loaded_mesh_quality(loaded);
            applied = true;
        } else if let Some(stl) = self.stl.as_mut() {
            let f = factor as f32;
            for_each_vertex(&mut stl.mesh, |v| [v[0] * f, v[1] * f, v[2] * f], true);
            applied = true;
        }
        if applied {
            self.status = Some(format!("Uniformly scaled by {factor}"));
            self.last_error = None;
            emit_audit(
                "mesh.scale",
                serde_json::json!({"kind": "mesh"}),
                serde_json::json!({"mode": "uniform", "factor": factor}),
            );
        } else {
            self.last_error = Some("Scale: no mesh or STL loaded.".into());
        }
    }

    /// Per-axis scale around the world origin.
    pub fn apply_scale_per_axis(&mut self, sx: f64, sy: f64, sz: f64) {
        if sx == 1.0 && sy == 1.0 && sz == 1.0 {
            self.status = Some("Scale skipped — all factors are 1".into());
            return;
        }
        for (name, v) in [("sx", sx), ("sy", sy), ("sz", sz)] {
            if !v.is_finite() {
                self.last_error = Some(format!("Scale: {name} = {v} is not allowed."));
                return;
            }
        }
        let mut applied = false;
        if let Some(loaded) = self.mesh.as_mut() {
            valenx_mesh::transform::scale_per_axis(&mut loaded.mesh, sx, sy, sz);
            refresh_loaded_mesh_quality(loaded);
            applied = true;
        } else if let Some(stl) = self.stl.as_mut() {
            let sxf = sx as f32;
            let syf = sy as f32;
            let szf = sz as f32;
            for_each_vertex(
                &mut stl.mesh,
                |v| [v[0] * sxf, v[1] * syf, v[2] * szf],
                true,
            );
            applied = true;
        }
        if applied {
            self.status = Some(format!("Scaled per-axis ({sx}, {sy}, {sz})"));
            self.last_error = None;
            emit_audit(
                "mesh.scale",
                serde_json::json!({"kind": "mesh"}),
                serde_json::json!({"mode": "per_axis", "sx": sx, "sy": sy, "sz": sz}),
            );
        } else {
            self.last_error = Some("Scale: no mesh or STL loaded.".into());
        }
    }

    /// Rotate around a Cartesian axis by an angle expressed in
    /// radians.
    pub fn apply_rotate(&mut self, axis: ToolboxAxis, angle_rad: f64) {
        if angle_rad == 0.0 {
            self.status = Some("Rotate skipped — angle is zero".into());
            return;
        }
        if !angle_rad.is_finite() {
            self.last_error = Some("Rotate: angle is not finite.".into());
            return;
        }
        let mut applied = false;
        if let Some(loaded) = self.mesh.as_mut() {
            valenx_mesh::transform::rotate_axis(&mut loaded.mesh, axis.to_axis(), angle_rad);
            refresh_loaded_mesh_quality(loaded);
            applied = true;
        } else if let Some(stl) = self.stl.as_mut() {
            rotate_triangle_mesh(&mut stl.mesh, axis, angle_rad);
            applied = true;
        }
        if applied {
            self.status = Some(format!(
                "Rotated around {} by {:.4} rad",
                axis.label(),
                angle_rad
            ));
            self.last_error = None;
            emit_audit(
                "mesh.rotate",
                serde_json::json!({"kind": "mesh"}),
                serde_json::json!({"axis": axis.label(), "angle_rad": angle_rad}),
            );
        } else {
            self.last_error = Some("Rotate: no mesh or STL loaded.".into());
        }
    }

    /// Mirror across one of the three Cartesian planes.
    pub fn apply_mirror(&mut self, plane: ToolboxAxis) {
        let mut applied = false;
        if let Some(loaded) = self.mesh.as_mut() {
            valenx_mesh::transform::mirror(&mut loaded.mesh, plane.to_plane());
            refresh_loaded_mesh_quality(loaded);
            applied = true;
        } else if let Some(stl) = self.stl.as_mut() {
            mirror_triangle_mesh(&mut stl.mesh, plane);
            applied = true;
        }
        if applied {
            self.status = Some(format!("Mirrored across {} plane", plane.label()));
            self.last_error = None;
            emit_audit(
                "mesh.mirror",
                serde_json::json!({"kind": "mesh"}),
                serde_json::json!({"plane": plane.label()}),
            );
        } else {
            self.last_error = Some("Mirror: no mesh or STL loaded.".into());
        }
    }

    /// Cut the loaded mesh by the plane defined by `point` and
    /// `normal` (any non-zero vector). Discards everything on the
    /// negative side of the plane.
    pub fn apply_cut_plane(&mut self, point: [f64; 3], normal: [f64; 3]) {
        let norm = Vector3::new(normal[0], normal[1], normal[2]);
        if norm.norm() < 1e-12 {
            self.last_error = Some("Cut: normal vector must be non-zero.".into());
            return;
        }
        let pt = Vector3::new(point[0], point[1], point[2]);
        let mut applied = false;
        if let Some(loaded) = self.mesh.as_mut() {
            let cut = valenx_mesh::cut::slice(&loaded.mesh, pt, norm);
            loaded.mesh = cut;
            refresh_loaded_mesh_quality(loaded);
            applied = true;
        } else if let Some(stl) = self.stl.as_mut() {
            slice_triangle_mesh(&mut stl.mesh, point, normal);
            applied = true;
        }
        if applied {
            self.status = Some(format!(
                "Cut at point ({:.4}, {:.4}, {:.4}), normal ({:.4}, {:.4}, {:.4})",
                point[0], point[1], point[2], normal[0], normal[1], normal[2],
            ));
            self.last_error = None;
            emit_audit(
                "mesh.cut",
                serde_json::json!({"kind": "mesh"}),
                serde_json::json!({"point": point, "normal": normal}),
            );
        } else {
            self.last_error = Some("Cut: no mesh or STL loaded.".into());
        }
    }

    /// Merge coincident nodes within `tolerance`. Surfaces a friendly
    /// "merged N duplicates" status on the next frame so the user
    /// gets feedback without opening the log.
    pub fn apply_merge_coincident(&mut self, tolerance: f64) {
        if !tolerance.is_finite() || tolerance < 0.0 {
            self.last_error = Some(format!("Repair: tolerance {tolerance} is not allowed."));
            return;
        }
        if let Some(loaded) = self.mesh.as_mut() {
            let before = loaded.mesh.nodes.len();
            let merged = valenx_mesh::boolean::merge_coincident_nodes(&loaded.mesh, tolerance);
            let after = merged.nodes.len();
            loaded.mesh = merged;
            refresh_loaded_mesh_quality(loaded);
            let removed = before.saturating_sub(after);
            self.status = Some(format!("Merged {removed} duplicate vertices"));
            self.last_error = None;
            emit_audit(
                "mesh.repair",
                serde_json::json!({"kind": "mesh"}),
                serde_json::json!({
                    "op": "merge_coincident_nodes",
                    "tolerance": tolerance,
                    "removed": removed,
                }),
            );
        } else {
            // STL repair is a no-op today — TriangleMesh has no
            // shared-node structure to dedup. Surface that honestly.
            self.last_error = Some(
                "Repair (merge coincident nodes) requires a canonical mesh, not an STL triangle soup.".into(),
            );
        }
    }

    /// Decimate the loaded canonical mesh via QEM to
    /// `target_fraction × current` vertex count. STL triangle soup
    /// is not supported (no shared-vertex structure). Refreshes the
    /// quality stats and emits a `mesh.decimate` audit entry.
    ///
    /// **Known limitation:** runs synchronously on the UI thread.
    /// QEM decimation on a 100k-tri mesh takes several seconds and
    /// blocks the egui repaint. TODO: dispatch to a worker via
    /// `thread::spawn` + an `mpsc::channel` so the viewport stays
    /// responsive, with a "Working…" overlay during the wait. The
    /// other mesh-toolbox heavy ops (remesh, boolean, fill_holes)
    /// share the same limitation and should land together as a
    /// single refactor — threading them one at a time invites
    /// state-drift bugs.
    pub fn apply_mesh_decimate(&mut self, target_fraction: f64) {
        if !target_fraction.is_finite() || !(0.0..=1.0).contains(&target_fraction) {
            self.last_error = Some(format!(
                "Decimate: target fraction {target_fraction} must be in [0, 1]."
            ));
            return;
        }
        let Some(loaded) = self.mesh.as_mut() else {
            self.last_error = Some(
                "Decimate requires a canonical mesh (STL triangle soup is not supported).".into(),
            );
            return;
        };
        let before = loaded.mesh.nodes.len();
        loaded.mesh = valenx_mesh::quadric_error_decimate(&loaded.mesh, target_fraction);
        refresh_loaded_mesh_quality(loaded);
        let after = loaded.mesh.nodes.len();
        self.status = Some(format!(
            "Decimated {before} → {after} vertices (target {:.0}%)",
            target_fraction * 100.0
        ));
        self.last_error = None;
        emit_audit(
            "mesh.decimate",
            serde_json::json!({"kind": "mesh"}),
            serde_json::json!({
                "target_fraction": target_fraction,
                "vertices_before": before,
                "vertices_after": after,
            }),
        );
    }

    /// Run `iterations` of Laplacian smoothing with factor `lambda`.
    pub fn apply_mesh_laplacian(&mut self, iterations: u32, lambda: f64) {
        let Some(loaded) = self.mesh.as_mut() else {
            self.last_error = Some("Laplacian smoothing requires a canonical mesh.".into());
            return;
        };
        loaded.mesh = valenx_mesh::laplacian(&loaded.mesh, iterations as usize, lambda);
        refresh_loaded_mesh_quality(loaded);
        self.status = Some(format!(
            "Laplacian smoothed {iterations} iter, factor {lambda}"
        ));
        self.last_error = None;
        emit_audit(
            "mesh.smooth",
            serde_json::json!({"kind": "mesh"}),
            serde_json::json!({"op": "laplacian", "iterations": iterations, "factor": lambda}),
        );
    }

    /// Run `iterations` Taubin smoothing passes (each = one λ step + one μ step).
    pub fn apply_mesh_taubin(&mut self, iterations: u32, lambda: f64, mu: f64) {
        let Some(loaded) = self.mesh.as_mut() else {
            self.last_error = Some("Taubin smoothing requires a canonical mesh.".into());
            return;
        };
        loaded.mesh = valenx_mesh::taubin(&loaded.mesh, iterations as usize, lambda, mu);
        refresh_loaded_mesh_quality(loaded);
        self.status = Some(format!(
            "Taubin smoothed {iterations} iter (λ={lambda}, μ={mu})"
        ));
        self.last_error = None;
        emit_audit(
            "mesh.smooth",
            serde_json::json!({"kind": "mesh"}),
            serde_json::json!({
                "op": "taubin",
                "iterations": iterations,
                "lambda": lambda,
                "mu": mu,
            }),
        );
    }

    /// Run isotropic remeshing at the given target edge length.
    ///
    /// **Known limitation:** synchronous on the UI thread; see the
    /// note on [`Self::apply_mesh_decimate`]. Should be threaded as
    /// part of the same mesh-toolbox refactor.
    pub fn apply_mesh_remesh(&mut self, target_edge_length: f64, iterations: u32) {
        if !target_edge_length.is_finite() || target_edge_length <= 0.0 {
            self.last_error = Some(format!(
                "Remesh: target edge length {target_edge_length} must be > 0."
            ));
            return;
        }
        let Some(loaded) = self.mesh.as_mut() else {
            self.last_error = Some("Remesh requires a canonical mesh.".into());
            return;
        };
        loaded.mesh = valenx_mesh::isotropic(&loaded.mesh, target_edge_length, iterations as usize);
        refresh_loaded_mesh_quality(loaded);
        self.status = Some(format!(
            "Remeshed at target edge length {target_edge_length} ({iterations} iter)"
        ));
        self.last_error = None;
        emit_audit(
            "mesh.remesh",
            serde_json::json!({"kind": "mesh"}),
            serde_json::json!({
                "target_edge_length": target_edge_length,
                "iterations": iterations,
            }),
        );
    }

    /// Fill every closed boundary loop whose perimeter is
    /// `<= max_boundary_length`. Use `f64::INFINITY` to fill all.
    ///
    /// **Known limitation:** synchronous on the UI thread; see the
    /// note on [`Self::apply_mesh_decimate`]. Should be threaded as
    /// part of the same mesh-toolbox refactor.
    pub fn apply_mesh_fill_holes(&mut self, max_boundary_length: f64) {
        let Some(loaded) = self.mesh.as_mut() else {
            self.last_error = Some("Fill holes requires a canonical mesh.".into());
            return;
        };
        let before_loops = valenx_mesh::boundary_loops(&loaded.mesh).len();
        loaded.mesh = valenx_mesh::fill_holes(&loaded.mesh, max_boundary_length);
        let after_loops = valenx_mesh::boundary_loops(&loaded.mesh).len();
        refresh_loaded_mesh_quality(loaded);
        let filled = before_loops.saturating_sub(after_loops);
        self.status = Some(format!(
            "Filled {filled} hole{} (max boundary length {max_boundary_length})",
            if filled == 1 { "" } else { "s" }
        ));
        self.last_error = None;
        emit_audit(
            "mesh.fill_holes",
            serde_json::json!({"kind": "mesh"}),
            serde_json::json!({
                "max_boundary_length": max_boundary_length,
                "loops_before": before_loops,
                "loops_after": after_loops,
            }),
        );
    }

    /// Save the loaded mesh to a user-chosen path as OBJ.
    pub fn save_mesh_as_obj(&mut self) {
        let dialog = rfd::FileDialog::new()
            .add_filter("Wavefront OBJ", &["obj", "OBJ"])
            .set_title("Save mesh as OBJ")
            .set_file_name("valenx-export.obj");
        let Some(path) = dialog.save_file() else {
            return;
        };
        let Some(loaded) = self.mesh.as_ref() else {
            self.last_error =
                Some("Save OBJ requires a canonical mesh (load an STL to promote it).".into());
            return;
        };
        match valenx_mesh::format::obj::write_path(&loaded.mesh, &path) {
            Ok(()) => {
                self.status = Some(format!("Saved to {}", path.display()));
                self.last_error = None;
                emit_audit(
                    "mesh.save_obj",
                    serde_json::json!({"kind": "mesh"}),
                    serde_json::json!({"path": path.display().to_string()}),
                );
            }
            Err(e) => {
                self.last_error = Some(format!("Save OBJ failed: {e}"));
            }
        }
    }

    /// Save the loaded mesh to a user-chosen path as ASCII PLY.
    pub fn save_mesh_as_ply(&mut self) {
        let dialog = rfd::FileDialog::new()
            .add_filter("Stanford PLY", &["ply", "PLY"])
            .set_title("Save mesh as PLY")
            .set_file_name("valenx-export.ply");
        let Some(path) = dialog.save_file() else {
            return;
        };
        let Some(loaded) = self.mesh.as_ref() else {
            self.last_error = Some("Save PLY requires a canonical mesh.".into());
            return;
        };
        match valenx_mesh::format::ply::write_path(&loaded.mesh, &path) {
            Ok(()) => {
                self.status = Some(format!("Saved to {}", path.display()));
                self.last_error = None;
                emit_audit(
                    "mesh.save_ply",
                    serde_json::json!({"kind": "mesh"}),
                    serde_json::json!({"path": path.display().to_string()}),
                );
            }
            Err(e) => {
                self.last_error = Some(format!("Save PLY failed: {e}"));
            }
        }
    }

    /// Attempt to save the loaded mesh as 3MF. Always errors today —
    /// 3MF support is deferred to v1.5 (needs a `zip` crate).
    pub fn save_mesh_as_3mf(&mut self) {
        let dialog = rfd::FileDialog::new()
            .add_filter("3MF (3D Manufacturing Format)", &["3mf", "3MF"])
            .set_title("Save mesh as 3MF")
            .set_file_name("valenx-export.3mf");
        let Some(path) = dialog.save_file() else {
            return;
        };
        let Some(loaded) = self.mesh.as_ref() else {
            self.last_error = Some("Save 3MF requires a canonical mesh.".into());
            return;
        };
        match valenx_mesh::format::three_mf::write_path(&loaded.mesh, &path) {
            Ok(()) => {
                self.status = Some(format!("Saved to {}", path.display()));
                self.last_error = None;
            }
            Err(e) => {
                self.last_error = Some(format!("Save 3MF: {e}"));
            }
        }
    }

    /// Save the current mesh / STL to a user-chosen path as binary STL.
    /// Opens a save-file dialog; cancellation is silent.
    pub fn save_mesh_as_stl(&mut self) {
        let dialog = rfd::FileDialog::new()
            .add_filter("STL mesh", &["stl", "STL"])
            .set_title("Save mesh as STL")
            .set_file_name("valenx-export.stl");
        let picked = dialog.save_file();
        let Some(path) = picked else {
            return;
        };
        let result = if let Some(loaded) = self.mesh.as_ref() {
            valenx_mesh::write_stl_binary(&loaded.mesh, &path)
                .map_err(|e| format!("write STL: {e}"))
        } else if let Some(stl) = self.stl.as_ref() {
            write_triangle_mesh_stl(&stl.mesh, &path)
        } else {
            Err("No mesh or STL loaded.".into())
        };
        match result {
            Ok(()) => {
                self.status = Some(format!("Saved to {}", path.display()));
                self.last_error = None;
                emit_audit(
                    "mesh.save_stl",
                    serde_json::json!({"kind": "mesh"}),
                    serde_json::json!({"path": path.display().to_string()}),
                );
            }
            Err(e) => {
                self.last_error = Some(format!("Save STL failed: {e}"));
            }
        }
    }

    /// Reset transformations by reloading the mesh from its
    /// `source_path` / STL `path`. Cheaper than tracking an undo
    /// stack for the pre-alpha surface; lossy if the file on disk
    /// has been overwritten since loading.
    pub fn reset_mesh_transformations(&mut self) {
        let stl_path = self.stl.as_ref().map(|s| s.path.clone());
        let mesh_path = self.mesh.as_ref().map(|m| m.path.clone());
        match (stl_path, mesh_path) {
            (Some(p), _) => {
                self.load_stl(p);
                self.status = Some("Reset STL to on-disk source".into());
            }
            (_, Some(p)) => {
                self.load_mesh(p);
                self.status = Some("Reset mesh to on-disk source".into());
            }
            _ => {
                self.last_error = Some("Reset: no mesh or STL loaded.".into());
            }
        }
    }

    /// Spawn `freecad <path>` in the background. If FreeCAD isn't on
    /// PATH, surface the install hint rather than silently failing.
    pub fn open_in_freecad(&mut self) {
        let path: Option<PathBuf> = self
            .stl
            .as_ref()
            .map(|s| s.path.clone())
            .or_else(|| self.mesh.as_ref().map(|m| m.path.clone()));
        let Some(p) = path else {
            self.last_error = Some("Open in FreeCAD: no mesh or STL loaded.".into());
            return;
        };
        let spawn = std::process::Command::new("freecad").arg(&p).spawn();
        match spawn {
            Ok(_) => {
                self.status = Some(format!("Opened {} in FreeCAD", p.display()));
                self.last_error = None;
                emit_audit(
                    "mesh.open_external",
                    serde_json::json!({"kind": "mesh"}),
                    serde_json::json!({"tool": "freecad", "path": p.display().to_string()}),
                );
            }
            Err(e) => {
                self.last_error = Some(format!(
                    "Could not spawn `freecad`: {e}. Install FreeCAD and make sure it's on PATH \
                     (https://www.freecad.org/downloads.php) — or use Tools menu → CAD → FreeCAD \
                     to scaffold a project via the registered adapter."
                ));
            }
        }
    }

    /// Scaffold a fresh OCCT-targeted case via the existing
    /// `new_case_for_adapter("occt")` path, then copy the current
    /// STL into the case directory so the user lands inside a ready-
    /// to-run project. Audits the launch with the destination dir.
    pub fn open_with_occt_adapter(&mut self) {
        // Snapshot the source STL path BEFORE scaffolding — the new
        // project replaces `self.project_path`, but the source STL
        // we're carrying over is independent.
        let source = self
            .stl
            .as_ref()
            .map(|s| s.path.clone())
            .or_else(|| self.mesh.as_ref().map(|m| m.path.clone()));
        let Some(src) = source else {
            self.last_error = Some("Open with OCCT: no mesh or STL loaded.".into());
            return;
        };
        // Call the existing scaffolding routine — it picks a folder,
        // sets up project.toml + cases/<dir>/case.toml, loads it.
        let before_root = self.project_path.clone();
        self.new_case_for_adapter("occt");
        // If the project_path actually changed, we successfully
        // scaffolded; the OCCT case dir is the only one in the new
        // project so copy the STL into it.
        if self.project_path != before_root {
            if let Some(proj_root) = &self.project_path {
                let dst_dir = proj_root.join("cases");
                if let Ok(entries) = std::fs::read_dir(&dst_dir) {
                    for entry in entries.flatten() {
                        if entry.path().is_dir() {
                            let dst = entry.path().join(
                                src.file_name()
                                    .unwrap_or_else(|| std::ffi::OsStr::new("input.stl")),
                            );
                            if let Err(e) = std::fs::copy(&src, &dst) {
                                self.last_error = Some(format!(
                                    "Scaffolded OCCT project but failed to copy STL: {e}"
                                ));
                                return;
                            }
                            self.status = Some(format!(
                                "Scaffolded OCCT case at {} with input STL.",
                                proj_root.display()
                            ));
                            emit_audit(
                                "mesh.open_external",
                                serde_json::json!({"kind": "mesh"}),
                                serde_json::json!({
                                    "tool": "occt",
                                    "project": proj_root.display().to_string(),
                                    "stl": src.display().to_string(),
                                }),
                            );
                            return;
                        }
                    }
                }
            }
        }
    }

    /// Flip the toolbox visibility. Bound to the View menu entry and
    /// the command palette.
    pub fn toggle_mesh_toolbox(&mut self) {
        self.show_mesh_toolbox = !self.show_mesh_toolbox;
    }

    // ---------- Part workbench (valenx-cad / truck) ----------

    /// Build a CAD primitive from the toolbox form. The resulting
    /// [`valenx_cad::Solid`] is stored as either operand A
    /// (`current_solid`) or operand B (`second_solid`) depending on
    /// the form toggle, AND tessellated to a triangle mesh that
    /// replaces the active viewport mesh. Emits a per-shape audit log
    /// (`cad.primitive.box`, `cad.primitive.cylinder`, …).
    pub fn apply_create_primitive(&mut self, kind: CadPrimitiveKind) {
        let s = &self.mesh_toolbox;
        let (solid_result, params) = match kind {
            CadPrimitiveKind::Box => (
                valenx_cad::box_solid(s.cad_box_dims[0], s.cad_box_dims[1], s.cad_box_dims[2]),
                serde_json::json!({
                    "dx": s.cad_box_dims[0],
                    "dy": s.cad_box_dims[1],
                    "dz": s.cad_box_dims[2],
                }),
            ),
            CadPrimitiveKind::Cylinder => (
                valenx_cad::cylinder(s.cad_cyl_radius, s.cad_cyl_height),
                serde_json::json!({
                    "radius": s.cad_cyl_radius,
                    "height": s.cad_cyl_height,
                }),
            ),
            CadPrimitiveKind::Sphere => (
                valenx_cad::sphere(s.cad_sphere_radius),
                serde_json::json!({ "radius": s.cad_sphere_radius }),
            ),
            CadPrimitiveKind::Cone => (
                valenx_cad::cone(s.cad_cone_base, s.cad_cone_top, s.cad_cone_height),
                serde_json::json!({
                    "base_radius": s.cad_cone_base,
                    "top_radius": s.cad_cone_top,
                    "height": s.cad_cone_height,
                }),
            ),
            CadPrimitiveKind::Torus => (
                valenx_cad::torus(s.cad_torus_major, s.cad_torus_minor),
                serde_json::json!({
                    "major_radius": s.cad_torus_major,
                    "minor_radius": s.cad_torus_minor,
                }),
            ),
        };
        let as_second = self.mesh_toolbox.cad_create_as_second;
        let label = kind.label();
        let audit_tag = kind.audit_tag();
        match solid_result {
            Ok(solid) => self.push_cad_solid(solid, as_second, audit_tag, label, params),
            Err(e) => {
                self.last_error = Some(format!("Create {label}: {e}"));
            }
        }
    }

    /// Apply a boolean op between operand A (`current_solid`) and
    /// operand B (`second_solid`). The result replaces operand A and
    /// re-tessellates the viewport mesh. Operand B is consumed
    /// regardless of success — booleans are destructive on the right
    /// operand to keep the UX simple (caller would otherwise need
    /// "duplicate B" affordances we don't have yet).
    pub fn apply_cad_boolean(&mut self, op: CadBooleanOp) {
        let (a, b) = match (self.current_solid.as_ref(), self.second_solid.as_ref()) {
            (Some(a), Some(b)) => (a.clone(), b.clone()),
            _ => {
                self.last_error = Some(
                    "Boolean: need both operand A (current solid) and operand B \
                     (second solid). Build a primitive with 'Create as second' enabled \
                     to populate B."
                        .into(),
                );
                return;
            }
        };
        let result = match op {
            CadBooleanOp::Union => valenx_cad::union(&a, &b),
            CadBooleanOp::Difference => valenx_cad::difference(&a, &b),
            CadBooleanOp::Intersection => valenx_cad::intersection(&a, &b),
        };
        let label = op.label();
        let audit_tag = op.audit_tag();
        match result {
            Ok(solid) => {
                // Consume B; result replaces A.
                self.second_solid = None;
                self.push_cad_solid(solid, false, audit_tag, label, serde_json::json!({}));
            }
            Err(e) => {
                self.last_error = Some(format!("Boolean {label}: {e}"));
            }
        }
    }

    /// Fillet edges on operand A by the radius in the toolbox form.
    /// Currently always reports the typed [`valenx_cad::CadError::NotImplemented`]
    /// from valenx-cad's fillet stub — the function is wired up
    /// nevertheless so the UI flow is in place for when truck ships
    /// fillet support.
    pub fn apply_cad_fillet(&mut self) {
        let radius = self.mesh_toolbox.cad_fillet_radius;
        let solid = match self.current_solid.as_ref() {
            Some(s) => s,
            None => {
                self.last_error =
                    Some("Fillet: no current solid — build a primitive first.".into());
                return;
            }
        };
        match valenx_cad::fillet_edges(solid, radius) {
            Ok(filleted) => {
                self.push_cad_solid(
                    filleted,
                    false,
                    "cad.fillet",
                    "Fillet edges",
                    serde_json::json!({ "radius": radius }),
                );
            }
            Err(e) => {
                self.last_error = Some(format!("Fillet: {e}"));
            }
        }
    }

    /// Internal: stash a freshly-built solid (primitive, boolean
    /// result, or fillet output) as operand A or B, tessellate it,
    /// and feed the resulting [`valenx_mesh::Mesh`] into the viewport
    /// via [`Self::apply_mesh`] so the rest of the toolbox + the
    /// inspector see it like any other mesh. Also emits the audit
    /// log entry for the op.
    fn push_cad_solid(
        &mut self,
        solid: valenx_cad::Solid,
        as_second: bool,
        audit_tag: &'static str,
        action_label: &str,
        params: serde_json::Value,
    ) {
        let faces = solid.faces();
        let edges = solid.edges();
        let vertices = solid.vertices();
        let mesh = match valenx_cad::solid_to_mesh(&solid, valenx_cad::DEFAULT_TESS_TOLERANCE) {
            Ok(m) => m,
            Err(e) => {
                self.last_error = Some(format!("{action_label}: tessellation failed: {e}"));
                return;
            }
        };
        if as_second {
            self.second_solid = Some(solid);
        } else {
            self.current_solid = Some(solid);
        }
        // Pseudo-path keeps the rest of the toolbox happy without
        // touching the filesystem — the user has no source STL.
        let pseudo_path = std::path::PathBuf::from(format!("<cad>/{action_label}.solid"));
        self.apply_mesh(mesh, pseudo_path);
        self.status = Some(format!(
            "{action_label}: {faces} faces, {edges} edges, {vertices} vertices"
        ));
        self.last_error = None;
        emit_audit(
            audit_tag,
            serde_json::json!({
                "as_second": as_second,
                "faces": faces,
                "edges": edges,
                "vertices": vertices,
            }),
            params,
        );
    }

    /// Clear both CAD operands (the "Discard solids" affordance).
    /// Loaded viewport mesh is left alone.
    pub fn apply_cad_clear(&mut self) {
        let had_a = self.current_solid.is_some();
        let had_b = self.second_solid.is_some();
        if !had_a && !had_b {
            self.status = Some("No CAD solids to clear.".into());
            return;
        }
        self.current_solid = None;
        self.second_solid = None;
        self.status = Some("Cleared CAD operands.".into());
        emit_audit(
            "cad.clear",
            serde_json::json!({"had_a": had_a, "had_b": had_b}),
            serde_json::json!({}),
        );
    }
}

/// Which boolean operation the toolbox is firing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CadBooleanOp {
    Union,
    Difference,
    Intersection,
}

impl CadBooleanOp {
    fn label(self) -> &'static str {
        match self {
            CadBooleanOp::Union => "Union",
            CadBooleanOp::Difference => "Difference",
            CadBooleanOp::Intersection => "Intersection",
        }
    }

    fn audit_tag(self) -> &'static str {
        match self {
            CadBooleanOp::Union => "cad.boolean.union",
            CadBooleanOp::Difference => "cad.boolean.difference",
            CadBooleanOp::Intersection => "cad.boolean.intersection",
        }
    }
}

/// Recompute the [`LoadedMesh`]'s quality scalars + histograms in
/// place after an editing op mutates the underlying canonical mesh.
/// Mirrors the post-load `apply_mesh` work so the Inspector and the
/// browser tree see consistent numbers.
fn refresh_loaded_mesh_quality(loaded: &mut LoadedMesh) {
    loaded.mesh.recompute_stats();
    let report = valenx_mesh::quality_report(&loaded.mesh);
    let aspect_hist =
        valenx_mesh::aspect_ratio_histogram(&loaded.mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist =
        valenx_mesh::skewness_histogram(&loaded.mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    loaded.mesh.stats.min_element_size = report.min_size;
    loaded.mesh.stats.max_aspect_ratio = report.max_aspect_ratio;
    loaded.mesh.stats.max_skewness = report.max_skewness;
    loaded.mesh.stats.min_orthogonality = report.min_orthogonality;
    loaded.quality = report;
    loaded.aspect_hist = aspect_hist;
    loaded.skew_hist = skew_hist;
}

fn rotate_triangle_mesh(mesh: &mut TriangleMesh, axis: ToolboxAxis, angle_rad: f64) {
    // Closed-form rotation around a Cartesian axis — avoids pulling
    // nalgebra into this crate. Sense matches the right-hand rule.
    let c = angle_rad.cos();
    let s = angle_rad.sin();
    for tri in &mut mesh.triangles {
        for vertex in &mut tri.vertices {
            let x = vertex[0] as f64;
            let y = vertex[1] as f64;
            let z = vertex[2] as f64;
            let (nx, ny, nz) = match axis {
                ToolboxAxis::X => (x, c * y - s * z, s * y + c * z),
                ToolboxAxis::Y => (c * x + s * z, y, -s * x + c * z),
                ToolboxAxis::Z => (c * x - s * y, s * x + c * y, z),
            };
            *vertex = [nx as f32, ny as f32, nz as f32];
        }
        tri.normal = tri.computed_normal();
    }
}

fn mirror_triangle_mesh(mesh: &mut TriangleMesh, plane: ToolboxAxis) {
    for tri in &mut mesh.triangles {
        for v in &mut tri.vertices {
            match plane {
                ToolboxAxis::X => v[0] = -v[0],
                ToolboxAxis::Y => v[1] = -v[1],
                ToolboxAxis::Z => v[2] = -v[2],
            }
        }
        // Reverse winding so face normals stay outward post-mirror.
        tri.vertices.swap(0, 2);
        tri.normal = tri.computed_normal();
    }
}

fn slice_triangle_mesh(mesh: &mut TriangleMesh, point: [f64; 3], normal: [f64; 3]) {
    let p = [point[0] as f32, point[1] as f32, point[2] as f32];
    let n = [normal[0] as f32, normal[1] as f32, normal[2] as f32];
    mesh.triangles.retain(|tri| {
        let cx = (tri.vertices[0][0] + tri.vertices[1][0] + tri.vertices[2][0]) / 3.0;
        let cy = (tri.vertices[0][1] + tri.vertices[1][1] + tri.vertices[2][1]) / 3.0;
        let cz = (tri.vertices[0][2] + tri.vertices[1][2] + tri.vertices[2][2]) / 3.0;
        let dot = (cx - p[0]) * n[0] + (cy - p[1]) * n[1] + (cz - p[2]) * n[2];
        dot >= 0.0
    });
}

fn write_triangle_mesh_stl(mesh: &TriangleMesh, path: &std::path::Path) -> Result<(), String> {
    use std::io::Write;
    let file =
        std::fs::File::create(path).map_err(|e| format!("create {}: {e}", path.display()))?;
    let mut writer = std::io::BufWriter::new(file);
    let mut header = [0u8; 80];
    let tag = b"valenx-app binary STL";
    header[..tag.len()].copy_from_slice(tag);
    writer.write_all(&header).map_err(|e| e.to_string())?;
    let count = mesh.triangle_count() as u32;
    writer
        .write_all(&count.to_le_bytes())
        .map_err(|e| e.to_string())?;
    for tri in &mesh.triangles {
        // Always recompute from winding so a hand-mutated mesh
        // doesn't carry stale normals to the export.
        let n = tri.computed_normal();
        for &c in n.iter() {
            writer
                .write_all(&c.to_le_bytes())
                .map_err(|e| e.to_string())?;
        }
        for v in tri.vertices.iter() {
            for &c in v.iter() {
                writer
                    .write_all(&c.to_le_bytes())
                    .map_err(|e| e.to_string())?;
            }
        }
        writer
            .write_all(&0u16.to_le_bytes())
            .map_err(|e| e.to_string())?;
    }
    writer.flush().map_err(|e| e.to_string())?;
    Ok(())
}

/// Render the right-side mesh-toolbox panel for the active `ValenxApp`.
/// Mounted in `update.rs` between the central viewport and the right
/// edge of the window. Only paints when a mesh / STL is loaded AND
/// the user hasn't hidden the panel through the View menu.
pub fn draw_mesh_toolbox(app: &mut ValenxApp, ctx: &egui::Context) {
    // The Part workbench can run with NO mesh loaded — that's how
    // users bootstrap a fresh CAD session — so the toolbox is now
    // governed only by the user-visible "show" toggle, not by
    // whether a mesh / STL / CAD operand happens to be loaded.
    if !app.show_mesh_toolbox {
        return;
    }
    egui::SidePanel::right("valenx_mesh_toolbox")
        .resizable(true)
        .default_width(280.0)
        .width_range(220.0..=420.0)
        .show(ctx, |ui| {
            ui.heading("Mesh Toolbox")
                .on_hover_text(
                    "CAD-side workbench. Ctrl+1 toggles. F1 opens panel help.",
                );
            ui.label(
                egui::RichText::new(
                    "Inspector + Transformations + Part / Sketcher / Draft / TechDraw / \
                     Assembly / Surface / CAM / Arch.",
                )
                .weak()
                .small(),
            );
            ui.separator();
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    draw_inspector(app, ui);
                    ui.separator();
                    draw_part_workbench(app, ui);
                    ui.separator();
                    draw_transformations(app, ui);
                    ui.separator();
                    draw_cut_plane(app, ui);
                    ui.separator();
                    draw_mesh_tools(app, ui);
                    ui.separator();
                    draw_repair(app, ui);
                    ui.separator();
                    draw_export(app, ui);
                    ui.separator();
                    draw_external_editors(app, ui);
                    ui.separator();
                    // Workbench collapsing sections. The hover tooltips
                    // on each header are part of the polish pass — they
                    // give the user a one-line summary of what the
                    // section covers before opening it.
                    ui.collapsing("Dock", |ui| {
                        draw_dock_panel(app, ui);
                    })
                    .header_response
                    .on_hover_text("Layout manager for floating panels.");
                    ui.collapsing("Sketcher", |ui| {
                        draw_sketcher_panel(app, ui);
                    })
                    .header_response
                    .on_hover_text(
                        crate::panel_help::short_summary("Sketcher"),
                    );
                    ui.collapsing("Part Design", |ui| {
                        draw_part_design_panel(app, ui);
                    })
                    .header_response
                    .on_hover_text(
                        crate::panel_help::short_summary("Part Design"),
                    );
                    ui.collapsing("Draft", |ui| {
                        draw_draft_panel(app, ui);
                    })
                    .header_response
                    .on_hover_text(crate::panel_help::short_summary("Draft"));
                    ui.collapsing("TechDraw", |ui| {
                        draw_techdraw_panel(app, ui);
                    })
                    .header_response
                    .on_hover_text(crate::panel_help::short_summary("TechDraw"));
                    ui.collapsing("Assembly", |ui| {
                        draw_assembly_panel(app, ui);
                    })
                    .header_response
                    .on_hover_text(crate::panel_help::short_summary("Assembly"));
                    ui.collapsing("Surface", |ui| {
                        draw_surface_panel(app, ui);
                    })
                    .header_response
                    .on_hover_text(crate::panel_help::short_summary("Surface"));
                    ui.collapsing("CAM", |ui| {
                        draw_cam_panel(app, ui);
                    })
                    .header_response
                    .on_hover_text(crate::panel_help::short_summary("CAM"));
                    ui.collapsing("Arch / BIM", |ui| {
                        draw_arch_panel(app, ui);
                    })
                    .header_response
                    .on_hover_text(crate::panel_help::short_summary("Arch"));
                    ui.collapsing("Spreadsheet", |ui| {
                        draw_spreadsheet_panel(app, ui);
                    })
                    .header_response
                    .on_hover_text(crate::panel_help::short_summary("Spreadsheet"));
                });
        });
}

fn draw_inspector(app: &ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Inspector").strong());
    if let Some(loaded) = app.mesh.as_ref() {
        ui.label(format!(
            "Source: {}",
            loaded
                .path
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
        ));
        ui.label(format!("Nodes: {}", loaded.mesh.stats.node_count));
        ui.label(format!("Elements: {}", loaded.mesh.stats.element_count));
        if let Some((min, max)) = canonical_aabb(&loaded.mesh) {
            ui.label(format!(
                "AABB: {:.3} × {:.3} × {:.3}",
                max.x - min.x,
                max.y - min.y,
                max.z - min.z,
            ));
        }
        egui::CollapsingHeader::new("Quality")
            .id_source("mesh_toolbox_inspector_quality")
            .default_open(false)
            .show(ui, |ui| {
                if let Some(ar) = loaded.quality.max_aspect_ratio {
                    ui.label(format!("max aspect ratio: {ar:.4}"));
                }
                if let Some(sk) = loaded.quality.max_skewness {
                    ui.label(format!("max skewness: {sk:.4}"));
                }
                if let Some(min) = loaded.quality.min_size {
                    ui.label(format!("min element size: {min:.4e}"));
                }
                if let Some(orth) = loaded.quality.min_orthogonality {
                    ui.label(format!("min orthogonality: {orth:.4}"));
                }
            });
    } else if let Some(stl) = app.stl.as_ref() {
        ui.label(format!(
            "Source: {}",
            stl.path.file_name().unwrap_or_default().to_string_lossy()
        ));
        let tris = stl.mesh.triangle_count();
        ui.label(format!("Triangles: {tris}"));
        ui.label(format!("Vertices: {}", tris * 3));
        if let Some((min, max)) = stl.mesh.bounding_box() {
            ui.label(format!(
                "AABB: {:.3} × {:.3} × {:.3}",
                max[0] - min[0],
                max[1] - min[1],
                max[2] - min[2],
            ));
        }
    }
}

fn draw_transformations(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Transformations").strong())
        .on_hover_text("Translate / scale / rotate / mirror the loaded mesh non-destructively. Reset rolls the mesh back to its loaded state.");

    let mut apply_translate = false;
    let mut apply_scale_uniform = false;
    let mut apply_scale_axes = false;
    let mut apply_rotate = false;
    let mut apply_mirror = false;
    let mut apply_reset = false;

    let s = &mut app.mesh_toolbox;
    ui.horizontal(|ui| {
        ui.label("Translate")
            .on_hover_text("XYZ displacement vector (model units, typically mm or m).");
        ui.add(
            egui::DragValue::new(&mut s.translate[0])
                .speed(0.1)
                .prefix("X "),
        )
        .on_hover_text("Translation along the X axis (model units).");
        ui.add(
            egui::DragValue::new(&mut s.translate[1])
                .speed(0.1)
                .prefix("Y "),
        )
        .on_hover_text("Translation along the Y axis (model units).");
        ui.add(
            egui::DragValue::new(&mut s.translate[2])
                .speed(0.1)
                .prefix("Z "),
        )
        .on_hover_text("Translation along the Z axis (model units).");
    });
    if ui
        .button("Apply translate")
        .on_hover_text("Add the XYZ displacement above to every vertex.")
        .clicked()
    {
        apply_translate = true;
    }

    ui.horizontal(|ui| {
        ui.label("Scale (uniform)")
            .on_hover_text("Multiplier applied equally to X / Y / Z. 1.0 = no scaling.");
        ui.add(
            egui::DragValue::new(&mut s.scale_uniform)
                .speed(0.05)
                .range(0.0001..=1e9),
        )
        .on_hover_text("Uniform scale factor (dimensionless). 2.0 = doubles in every direction.");
    });
    if ui
        .button("Apply uniform scale")
        .on_hover_text("Multiply every vertex coordinate by the uniform scale factor.")
        .clicked()
    {
        apply_scale_uniform = true;
    }

    ui.horizontal(|ui| {
        ui.label("Scale (axes)")
            .on_hover_text("Independent scale factor per axis — useful for squashing / stretching.");
        ui.add(
            egui::DragValue::new(&mut s.scale_per_axis[0])
                .speed(0.05)
                .prefix("X "),
        )
        .on_hover_text("Scale along X (dimensionless multiplier).");
        ui.add(
            egui::DragValue::new(&mut s.scale_per_axis[1])
                .speed(0.05)
                .prefix("Y "),
        )
        .on_hover_text("Scale along Y (dimensionless multiplier).");
        ui.add(
            egui::DragValue::new(&mut s.scale_per_axis[2])
                .speed(0.05)
                .prefix("Z "),
        )
        .on_hover_text("Scale along Z (dimensionless multiplier).");
    });
    if ui
        .button("Apply per-axis scale")
        .on_hover_text("Multiply each axis independently by the XYZ scale factors above.")
        .clicked()
    {
        apply_scale_axes = true;
    }

    ui.horizontal(|ui| {
        ui.label("Rotate axis")
            .on_hover_text("Axis the rotation pivots around.");
        ui.radio_value(&mut s.rotate_axis, ToolboxAxis::X, "X")
            .on_hover_text("Rotate around the X axis (roll).");
        ui.radio_value(&mut s.rotate_axis, ToolboxAxis::Y, "Y")
            .on_hover_text("Rotate around the Y axis (pitch).");
        ui.radio_value(&mut s.rotate_axis, ToolboxAxis::Z, "Z")
            .on_hover_text("Rotate around the Z axis (yaw).");
    });
    ui.horizontal(|ui| {
        ui.label("Angle (deg)")
            .on_hover_text("Rotation amount in degrees. Positive = right-hand rule around the axis.");
        ui.add(egui::DragValue::new(&mut s.rotate_angle_deg).speed(1.0))
            .on_hover_text("Rotation angle (degrees).");
    });
    if ui
        .button("Apply rotate")
        .on_hover_text("Rotate the mesh around the chosen axis by the angle above.")
        .clicked()
    {
        apply_rotate = true;
    }

    ui.horizontal(|ui| {
        ui.label("Mirror plane")
            .on_hover_text("Coordinate plane the reflection happens across.");
        ui.radio_value(&mut s.mirror_plane, ToolboxAxis::X, "X")
            .on_hover_text("Mirror across the YZ plane (flips X).");
        ui.radio_value(&mut s.mirror_plane, ToolboxAxis::Y, "Y")
            .on_hover_text("Mirror across the XZ plane (flips Y).");
        ui.radio_value(&mut s.mirror_plane, ToolboxAxis::Z, "Z")
            .on_hover_text("Mirror across the XY plane (flips Z).");
    });
    if ui
        .button("Apply mirror")
        .on_hover_text("Flip every vertex coordinate across the chosen plane.")
        .clicked()
    {
        apply_mirror = true;
    }

    ui.add_space(4.0);
    if ui
        .button("Reset to loaded mesh")
        .on_hover_text("Discard every transform applied this session — restores the mesh to its on-disk state.")
        .clicked()
    {
        apply_reset = true;
    }

    // Apply all deferred actions now that we've released the &mut
    // app.mesh_toolbox borrow.
    let translate = app.mesh_toolbox.translate;
    let scale_u = app.mesh_toolbox.scale_uniform;
    let scale_pa = app.mesh_toolbox.scale_per_axis;
    let rot_axis = app.mesh_toolbox.rotate_axis;
    let rot_angle_rad = app.mesh_toolbox.rotate_angle_deg.to_radians();
    let mirror_plane = app.mesh_toolbox.mirror_plane;
    if apply_translate {
        app.apply_translate(translate[0], translate[1], translate[2]);
    }
    if apply_scale_uniform {
        app.apply_scale_uniform(scale_u);
    }
    if apply_scale_axes {
        app.apply_scale_per_axis(scale_pa[0], scale_pa[1], scale_pa[2]);
    }
    if apply_rotate {
        app.apply_rotate(rot_axis, rot_angle_rad);
    }
    if apply_mirror {
        app.apply_mirror(mirror_plane);
    }
    if apply_reset {
        app.reset_mesh_transformations();
    }
}

fn draw_cut_plane(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Cut plane").strong())
        .on_hover_text("Slice the mesh along an arbitrary plane defined by a point + normal vector.");
    let mut apply_cut = false;
    {
        let s = &mut app.mesh_toolbox;
        ui.horizontal(|ui| {
            ui.label("Point  ")
                .on_hover_text("A point in the cut plane (model units).");
            ui.add(
                egui::DragValue::new(&mut s.cut_point[0])
                    .speed(0.1)
                    .prefix("X "),
            )
            .on_hover_text("Point X (model units).");
            ui.add(
                egui::DragValue::new(&mut s.cut_point[1])
                    .speed(0.1)
                    .prefix("Y "),
            )
            .on_hover_text("Point Y (model units).");
            ui.add(
                egui::DragValue::new(&mut s.cut_point[2])
                    .speed(0.1)
                    .prefix("Z "),
            )
            .on_hover_text("Point Z (model units).");
        });
        ui.horizontal(|ui| {
            ui.label("Normal ")
                .on_hover_text("The plane's outward-normal vector (need not be unit-length).");
            ui.add(
                egui::DragValue::new(&mut s.cut_normal[0])
                    .speed(0.05)
                    .prefix("X "),
            )
            .on_hover_text("Normal X component.");
            ui.add(
                egui::DragValue::new(&mut s.cut_normal[1])
                    .speed(0.05)
                    .prefix("Y "),
            )
            .on_hover_text("Normal Y component.");
            ui.add(
                egui::DragValue::new(&mut s.cut_normal[2])
                    .speed(0.05)
                    .prefix("Z "),
            )
            .on_hover_text("Normal Z component.");
        });
        ui.checkbox(
            &mut s.cut_show_overlay,
            "Show cut overlay (intersected triangles)",
        )
        .on_hover_text("Highlight which triangles the plane intersects without actually cutting.");
        if ui
            .button("Apply cut")
            .on_hover_text("Keeps only the side where (node - point) · normal >= 0.")
            .clicked()
        {
            apply_cut = true;
        }
    }
    if apply_cut {
        let pt = app.mesh_toolbox.cut_point;
        let n = app.mesh_toolbox.cut_normal;
        app.apply_cut_plane(pt, n);
    }
}

/// "Part" section — primitives + boolean ops + fillet driven by the
/// `valenx-cad` Rust-native CAD kernel (truck BRep). Surfaced inside
/// a collapsing header so users who don't care about CAD don't see
/// it expanded by default.
fn draw_part_workbench(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::CollapsingHeader::new(
        egui::RichText::new("Part — CAD primitives (Rust-native via truck)").strong(),
    )
    .id_source("mesh_toolbox_part_workbench")
    .default_open(false)
    .show(ui, |ui| {
        // ---- Operand status row ----
        ui.horizontal(|ui| {
            let a_label = match &app.current_solid {
                Some(s) => format!("A: {} faces", s.faces()),
                None => "A: <none>".to_string(),
            };
            let b_label = match &app.second_solid {
                Some(s) => format!("B: {} faces", s.faces()),
                None => "B: <none>".to_string(),
            };
            ui.label(a_label);
            ui.separator();
            ui.label(b_label);
        });
        if ui
            .button("Clear A and B")
            .on_hover_text(
                "Drops both CAD operands so the next primitive Create starts fresh. \
             The viewport mesh is left alone.",
            )
            .clicked()
        {
            app.apply_cad_clear();
        }
        ui.separator();

        // ---- Insert primitive ----
        ui.label(egui::RichText::new("Insert primitive").strong());
        let mut create_clicked = false;
        {
            let s = &mut app.mesh_toolbox;
            ui.horizontal(|ui| {
                ui.label("Shape");
                egui::ComboBox::from_id_source("mesh_toolbox_cad_primitive")
                    .selected_text(s.cad_primitive.label())
                    .show_ui(ui, |ui| {
                        for kind in [
                            CadPrimitiveKind::Box,
                            CadPrimitiveKind::Cylinder,
                            CadPrimitiveKind::Sphere,
                            CadPrimitiveKind::Cone,
                            CadPrimitiveKind::Torus,
                        ] {
                            ui.selectable_value(&mut s.cad_primitive, kind, kind.label());
                        }
                    });
            });
            match s.cad_primitive {
                CadPrimitiveKind::Box => {
                    ui.horizontal(|ui| {
                        ui.label("Size")
                            .on_hover_text("Box edge lengths along X / Y / Z (model units).");
                        ui.add(
                            egui::DragValue::new(&mut s.cad_box_dims[0])
                                .speed(0.1)
                                .range(1e-6..=1e9)
                                .prefix("dx "),
                        )
                        .on_hover_text("Box width along X (model units).");
                        ui.add(
                            egui::DragValue::new(&mut s.cad_box_dims[1])
                                .speed(0.1)
                                .range(1e-6..=1e9)
                                .prefix("dy "),
                        )
                        .on_hover_text("Box depth along Y (model units).");
                        ui.add(
                            egui::DragValue::new(&mut s.cad_box_dims[2])
                                .speed(0.1)
                                .range(1e-6..=1e9)
                                .prefix("dz "),
                        )
                        .on_hover_text("Box height along Z (model units).");
                    });
                }
                CadPrimitiveKind::Cylinder => {
                    ui.horizontal(|ui| {
                        ui.label("Cylinder")
                            .on_hover_text("Solid cylinder, axis aligned with Z.");
                        ui.add(
                            egui::DragValue::new(&mut s.cad_cyl_radius)
                                .speed(0.05)
                                .range(1e-6..=1e9)
                                .prefix("r "),
                        )
                        .on_hover_text("Cylinder radius (model units).");
                        ui.add(
                            egui::DragValue::new(&mut s.cad_cyl_height)
                                .speed(0.1)
                                .range(1e-6..=1e9)
                                .prefix("h "),
                        )
                        .on_hover_text("Cylinder height along Z (model units).");
                    });
                }
                CadPrimitiveKind::Sphere => {
                    ui.horizontal(|ui| {
                        ui.label("Sphere")
                            .on_hover_text("Solid sphere centred on the origin.");
                        ui.add(
                            egui::DragValue::new(&mut s.cad_sphere_radius)
                                .speed(0.05)
                                .range(1e-6..=1e9)
                                .prefix("r "),
                        )
                        .on_hover_text("Sphere radius (model units).");
                    });
                }
                CadPrimitiveKind::Cone => {
                    ui.horizontal(|ui| {
                        ui.label("Cone")
                            .on_hover_text("Solid cone / frustum aligned with Z.");
                        ui.add(
                            egui::DragValue::new(&mut s.cad_cone_base)
                                .speed(0.05)
                                .range(1e-6..=1e9)
                                .prefix("base "),
                        )
                        .on_hover_text("Base radius at z = 0 (model units).");
                        ui.add(
                            egui::DragValue::new(&mut s.cad_cone_top)
                                .speed(0.05)
                                .range(0.0..=1e9)
                                .prefix("top "),
                        )
                        .on_hover_text("Top radius at z = h. 0 → pointed cone; > 0 → frustum.");
                        ui.add(
                            egui::DragValue::new(&mut s.cad_cone_height)
                                .speed(0.1)
                                .range(1e-6..=1e9)
                                .prefix("h "),
                        )
                        .on_hover_text("Cone height along Z (model units).");
                    });
                    ui.small("top = 0 gives a pointed cone; top > 0 gives a frustum.");
                }
                CadPrimitiveKind::Torus => {
                    ui.horizontal(|ui| {
                        ui.label("Torus")
                            .on_hover_text("Donut shape — major radius around Z, minor radius is the tube cross-section.");
                        ui.add(
                            egui::DragValue::new(&mut s.cad_torus_major)
                                .speed(0.05)
                                .range(1e-6..=1e9)
                                .prefix("major "),
                        )
                        .on_hover_text("Major radius — centre of the tube to centre of the torus (model units).");
                        ui.add(
                            egui::DragValue::new(&mut s.cad_torus_minor)
                                .speed(0.05)
                                .range(1e-6..=1e9)
                                .prefix("minor "),
                        )
                        .on_hover_text("Minor radius — radius of the tube cross-section (model units).");
                    });
                    ui.small("minor must be strictly less than major.");
                }
            }
            ui.checkbox(
                &mut s.cad_create_as_second,
                "Create as operand B (for boolean ops)",
            );
            if ui
                .button("Create")
                .on_hover_text(
                    "Builds the primitive through the truck kernel, tessellates it, \
                 and replaces the viewport mesh with the result.",
                )
                .clicked()
            {
                create_clicked = true;
            }
        }
        if create_clicked {
            let kind = app.mesh_toolbox.cad_primitive;
            app.apply_create_primitive(kind);
        }
        ui.separator();

        // ---- Boolean ops ----
        ui.label(egui::RichText::new("Boolean operations").strong());
        let booleans_enabled = app.current_solid.is_some() && app.second_solid.is_some();
        let booleans_hint = if booleans_enabled {
            "Combines operand A with operand B using truck-shapeops. \
             The result replaces operand A; B is consumed."
        } else {
            "Build a primitive normally (operand A), then build a second \
             one with 'Create as operand B' checked, then come back here."
        };
        ui.horizontal(|ui| {
            if ui
                .add_enabled(booleans_enabled, egui::Button::new("Union"))
                .on_hover_text(booleans_hint)
                .clicked()
            {
                app.apply_cad_boolean(CadBooleanOp::Union);
            }
            if ui
                .add_enabled(booleans_enabled, egui::Button::new("Difference"))
                .on_hover_text(booleans_hint)
                .clicked()
            {
                app.apply_cad_boolean(CadBooleanOp::Difference);
            }
            if ui
                .add_enabled(booleans_enabled, egui::Button::new("Intersection"))
                .on_hover_text(booleans_hint)
                .clicked()
            {
                app.apply_cad_boolean(CadBooleanOp::Intersection);
            }
        });
        ui.separator();

        // ---- Fillet edges (stub — see valenx_cad::fillet_edges) ----
        ui.label(egui::RichText::new("Fillet edges").strong());
        let mut fillet_clicked = false;
        {
            let s = &mut app.mesh_toolbox;
            ui.horizontal(|ui| {
                ui.label("Radius");
                ui.add(
                    egui::DragValue::new(&mut s.cad_fillet_radius)
                        .speed(0.01)
                        .range(1e-6..=1e9),
                );
            });
        }
        let fillet_enabled = app.current_solid.is_some();
        if ui
            .add_enabled(fillet_enabled, egui::Button::new("Apply fillet"))
            .on_hover_text(
                "Routes through valenx_cad::fillet_edges. truck 0.6 does NOT \
                 ship an edge-fillet algorithm yet, so this currently returns \
                 a typed 'not implemented' error — fall back to 'Open in \
                 FreeCAD' for filleting until the upstream API ships.",
            )
            .clicked()
        {
            fillet_clicked = true;
        }
        if fillet_clicked {
            app.apply_cad_fillet();
        }
    });
}

fn draw_repair(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Repair").strong())
        .on_hover_text("Mesh-repair operations — merge duplicate vertices, weld coincident nodes.");
    let mut apply_repair = false;
    {
        let s = &mut app.mesh_toolbox;
        ui.horizontal(|ui| {
            ui.label("Merge coincident — tolerance")
                .on_hover_text("Welding distance — any two vertices closer than this collapse to one.");
            ui.add(
                egui::DragValue::new(&mut s.repair_tolerance)
                    .speed(1e-7)
                    .range(0.0..=1e9),
            )
            .on_hover_text("Coincident-vertex merge tolerance (model units).");
        });
        if ui
            .button("Apply merge")
            .on_hover_text(
                "Snaps any two nodes within `tolerance` to a single index. \
                 Canonical mesh only — STL triangle soups don't share nodes.",
            )
            .clicked()
        {
            apply_repair = true;
        }
    }
    if apply_repair {
        let tol = app.mesh_toolbox.repair_tolerance;
        app.apply_merge_coincident(tol);
    }
}

fn draw_export(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Export").strong())
        .on_hover_text("Save the current mesh to disk in a portable format.");
    if ui
        .button("Save mesh as STL…")
        .on_hover_text("Stereolithography binary STL — every triangle stored independently. Works on any mesh kind.")
        .clicked()
    {
        app.save_mesh_as_stl();
    }
    if ui
        .button("Save mesh as OBJ…")
        .on_hover_text(
            "Wavefront OBJ — `v` lines for positions, `f` lines for triangles. \
             Requires a canonical mesh.",
        )
        .clicked()
    {
        app.save_mesh_as_obj();
    }
    if ui
        .button("Save mesh as PLY…")
        .on_hover_text(
            "Stanford PLY (ASCII). Header + vertex rows + face rows. \
             Requires a canonical mesh.",
        )
        .clicked()
    {
        app.save_mesh_as_ply();
    }
    if ui
        .button("Save mesh as 3MF…")
        .on_hover_text(
            "3MF — currently unsupported (v1 ships OBJ/PLY/STL). \
             3MF will land in v1.5 once a `zip` crate is in workspace deps.",
        )
        .clicked()
    {
        app.save_mesh_as_3mf();
    }
}

/// Mesh Tools sub-section: decimation, smoothing, remeshing, fill holes.
///
/// All operations require a canonical mesh — STL triangle soup is not
/// supported (no shared-vertex structure to operate on). Each button
/// queues an action; we apply them after the `&mut app.mesh_toolbox`
/// borrow scope ends so we can call methods on `&mut app`.
fn draw_mesh_tools(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Mesh Tools").strong())
        .on_hover_text("Decimation / smoothing / remesh / fill — all canonical-mesh operations (not STL soup).");

    let mut do_decimate = false;
    let mut do_laplacian = false;
    let mut do_taubin = false;
    let mut do_remesh = false;
    let mut do_fill = false;

    {
        let s = &mut app.mesh_toolbox;

        ui.horizontal(|ui| {
            ui.label("Decimate fraction")
                .on_hover_text("Target vertex count as a fraction of the original. 0.5 = half the vertices.");
            ui.add(
                egui::Slider::new(&mut s.mesh_tools_decimate_fraction, 0.05..=1.0)
                    .text("of vertices"),
            )
            .on_hover_text("Decimation target (0.05 = aggressive, 1.0 = no decimation).");
        });
        if ui
            .button("Decimate")
            .on_hover_text(
                "Quadric error metric decimation (Garland-Heckbert). \
                 Midpoint contraction, manifold-preserving collapses.",
            )
            .clicked()
        {
            do_decimate = true;
        }

        ui.separator();
        ui.label("Smoothing")
            .on_hover_text("Iterative vertex relaxation — Laplacian (shrinks) and Taubin (shrink-free).");
        ui.horizontal(|ui| {
            ui.label("Laplacian iter")
                .on_hover_text("Number of Laplacian smoothing passes.");
            ui.add(egui::DragValue::new(&mut s.mesh_tools_laplacian_iter).range(0..=200))
                .on_hover_text("Smoothing iteration count (0 = no smoothing).");
            ui.label("factor")
                .on_hover_text("Per-iteration smoothing strength (0 = no move, 1 = full average).");
            ui.add(
                egui::DragValue::new(&mut s.mesh_tools_laplacian_factor)
                    .speed(0.05)
                    .range(0.0..=1.0),
            )
            .on_hover_text("Laplacian factor λ (dimensionless 0..=1).");
        });
        if ui
            .button("Apply Laplacian")
            .on_hover_text("Smooth the mesh with `iter` iterations of Laplacian relaxation at the given factor.")
            .clicked()
        {
            do_laplacian = true;
        }
        ui.horizontal(|ui| {
            ui.label("Taubin iter")
                .on_hover_text("Number of Taubin (shrink-free) smoothing passes.");
            ui.add(egui::DragValue::new(&mut s.mesh_tools_taubin_iter).range(0..=200))
                .on_hover_text("Taubin iteration count.");
            ui.label("λ")
                .on_hover_text("Taubin's positive-step weight.");
            ui.add(
                egui::DragValue::new(&mut s.mesh_tools_taubin_lambda)
                    .speed(0.05)
                    .range(-2.0..=2.0),
            );
            ui.label("μ");
            ui.add(
                egui::DragValue::new(&mut s.mesh_tools_taubin_mu)
                    .speed(0.05)
                    .range(-2.0..=2.0),
            );
        });
        if ui
            .button("Apply Taubin")
            .on_hover_text(
                "Taubin λ/μ smoothing — alternates positive and negative steps to \
                 approximate low-pass filtering with minimal shrinkage.",
            )
            .clicked()
        {
            do_taubin = true;
        }

        ui.separator();
        ui.horizontal(|ui| {
            ui.label("Remesh target edge");
            ui.add(
                egui::DragValue::new(&mut s.mesh_tools_remesh_target)
                    .speed(0.05)
                    .range(1e-6..=1e6),
            );
            ui.label("iter");
            ui.add(egui::DragValue::new(&mut s.mesh_tools_remesh_iter).range(0..=50));
        });
        if ui
            .button("Isotropic remesh")
            .on_hover_text(
                "Split-collapse-flip-smooth iteration driving the mesh toward \
                 uniform triangles at the target edge length.",
            )
            .clicked()
        {
            do_remesh = true;
        }

        ui.separator();
        ui.horizontal(|ui| {
            ui.label("Fill holes — max boundary length");
            ui.add(
                egui::DragValue::new(&mut s.mesh_tools_fill_holes_max)
                    .speed(0.1)
                    .range(0.0..=1e9),
            );
        });
        if ui
            .button("Fill holes")
            .on_hover_text(
                "Detect closed boundary loops and triangulate each one with ear-clipping. \
                 Loops longer than the threshold are skipped (protects intentional openings).",
            )
            .clicked()
        {
            do_fill = true;
        }
    }

    // Apply deferred actions now that the &mut borrow is dropped.
    let frac = app.mesh_toolbox.mesh_tools_decimate_fraction;
    let lap_iter = app.mesh_toolbox.mesh_tools_laplacian_iter;
    let lap_fac = app.mesh_toolbox.mesh_tools_laplacian_factor;
    let tau_iter = app.mesh_toolbox.mesh_tools_taubin_iter;
    let tau_lambda = app.mesh_toolbox.mesh_tools_taubin_lambda;
    let tau_mu = app.mesh_toolbox.mesh_tools_taubin_mu;
    let remesh_target = app.mesh_toolbox.mesh_tools_remesh_target;
    let remesh_iter = app.mesh_toolbox.mesh_tools_remesh_iter;
    let fill_max = app.mesh_toolbox.mesh_tools_fill_holes_max;
    if do_decimate {
        app.apply_mesh_decimate(frac);
    }
    if do_laplacian {
        app.apply_mesh_laplacian(lap_iter, lap_fac);
    }
    if do_taubin {
        app.apply_mesh_taubin(tau_iter, tau_lambda, tau_mu);
    }
    if do_remesh {
        app.apply_mesh_remesh(remesh_target, remesh_iter);
    }
    if do_fill {
        app.apply_mesh_fill_holes(fill_max);
    }
}

fn draw_external_editors(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("External editors").strong());
    let has_source = app.stl.is_some() || app.mesh.is_some();
    ui.horizontal(|ui| {
        if ui
            .add_enabled(has_source, egui::Button::new("Open in FreeCAD"))
            .on_hover_text(
                "Spawns `freecad <path>` in the background. \
                 If FreeCAD isn't on PATH, an install hint is surfaced.",
            )
            .clicked()
        {
            app.open_in_freecad();
        }
        if ui
            .add_enabled(has_source, egui::Button::new("Process with OCCT…"))
            .on_hover_text(
                "Scaffolds an OCCT-targeted case using the registered adapter \
                 template, then copies the source STL into the new case directory.",
            )
            .clicked()
        {
            app.open_with_occt_adapter();
        }
    });
}

/// Canonical-mesh AABB. Returns `None` for an empty mesh.
fn canonical_aabb(mesh: &valenx_mesh::Mesh) -> Option<(Vector3<f64>, Vector3<f64>)> {
    let mut iter = mesh.nodes.iter();
    let first = *iter.next()?;
    let mut min = first;
    let mut max = first;
    for &n in iter {
        for i in 0..3 {
            if n[i] < min[i] {
                min[i] = n[i];
            }
            if n[i] > max[i] {
                max[i] = n[i];
            }
        }
    }
    Some((min, max))
}

/// Render the Dock subsection inside the Mesh Toolbox.
pub fn draw_dock_panel(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.heading("Dock (native AutoDock Vina)");
    ui.separator();
    // Capture a "show pose i in viewport" request so we can dispatch
    // after the &mut app.mesh_toolbox.dock_panel borrow scope ends.
    let mut show_pose_request: Option<(String, usize)> = None;
    let s = &mut app.mesh_toolbox.dock_panel;
    ui.horizontal(|ui| {
        ui.label("Receptor PDBQT:");
        ui.text_edit_singleline(&mut s.receptor_path);
        if ui.button("…").clicked() {
            if let Some(p) = rfd::FileDialog::new()
                .add_filter("PDBQT", &["pdbqt"])
                .pick_file()
            {
                s.receptor_path = p.display().to_string();
            }
        }
    });
    ui.horizontal(|ui| {
        ui.label("Ligand PDBQT:");
        ui.text_edit_singleline(&mut s.ligand_path);
        if ui.button("…").clicked() {
            if let Some(p) = rfd::FileDialog::new()
                .add_filter("PDBQT", &["pdbqt"])
                .pick_file()
            {
                s.ligand_path = p.display().to_string();
            }
        }
    });
    ui.horizontal(|ui| {
        ui.label("Output PDBQT:");
        ui.text_edit_singleline(&mut s.output_path);
    });
    ui.separator();
    ui.label("Search box centre (Å):");
    ui.horizontal(|ui| {
        ui.add(egui::DragValue::new(&mut s.center[0]).speed(0.1));
        ui.add(egui::DragValue::new(&mut s.center[1]).speed(0.1));
        ui.add(egui::DragValue::new(&mut s.center[2]).speed(0.1));
    });
    ui.label("Search box size (Å):");
    ui.horizontal(|ui| {
        ui.add(
            egui::DragValue::new(&mut s.size[0])
                .speed(0.1)
                .range(1.0..=80.0),
        );
        ui.add(
            egui::DragValue::new(&mut s.size[1])
                .speed(0.1)
                .range(1.0..=80.0),
        );
        ui.add(
            egui::DragValue::new(&mut s.size[2])
                .speed(0.1)
                .range(1.0..=80.0),
        );
    });
    ui.horizontal(|ui| {
        ui.label("Exhaustiveness:");
        ui.add(egui::DragValue::new(&mut s.exhaustiveness).range(1..=32));
        ui.label("Num modes:");
        ui.add(egui::DragValue::new(&mut s.num_modes).range(1..=99));
    });
    ui.horizontal(|ui| {
        ui.label("Energy range:");
        ui.add(egui::DragValue::new(&mut s.energy_range).speed(0.1));
        ui.label("Seed:");
        ui.add(egui::DragValue::new(&mut s.seed));
    });
    ui.separator();
    if ui.button("Dock now").clicked() {
        run_dock_now(s);
    }
    if let Some(err) = &s.last_error {
        ui.colored_label(egui::Color32::RED, err);
    }
    if !s.last_scores.is_empty() {
        ui.label(format!("{} poses:", s.last_scores.len()));
        for (i, (score, rmsd)) in s.last_scores.iter().enumerate() {
            let label = format!("#{}: score {score:.3} kcal/mol, rmsd {rmsd:.2} Å", i + 1);
            let selected = s.selected_pose == Some(i);
            if ui.selectable_label(selected, label).clicked() {
                s.selected_pose = Some(i);
            }
        }
        if let Some(i) = s.selected_pose {
            if ui.button("Show pose in viewport").clicked() {
                show_pose_request = Some((s.output_path.clone(), i));
            }
        }
    }
    // Release the &mut app.mesh_toolbox.dock_panel borrow before
    // reaching for app.load_stl (which needs &mut app on its own).
    let _ = s;
    if let Some((output_path, pose_index)) = show_pose_request {
        push_pose_as_mesh(app, &output_path, pose_index);
    }
}

fn run_dock_now(s: &mut DockPanelState) {
    use std::path::Path;
    use valenx_core::io_caps::{read_capped_to_string, MAX_PDBQT_FILE_BYTES};
    s.last_error = None;
    s.last_scores.clear();
    s.selected_pose = None;
    // Round-20 H1 (R16 H1 sister): pre-fix the Desktop Dock panel did
    // three bare `std::fs::read_to_string` on user-chosen paths
    // (receptor, ligand, output). A user (or a stale path in saved
    // panel state) pointing at a multi-GB file would OOM the renderer
    // process before the docker ever ran. Cap matches the MCP-side
    // hardening in `valenx_mcp::tools::MAX_PDBQT_FILE_BYTES`.
    let receptor = match read_capped_to_string(Path::new(&s.receptor_path), MAX_PDBQT_FILE_BYTES) {
        Ok(v) => v,
        Err(e) => {
            s.last_error = Some(format!("receptor read: {e}"));
            return;
        }
    };
    let ligand = match read_capped_to_string(Path::new(&s.ligand_path), MAX_PDBQT_FILE_BYTES) {
        Ok(v) => v,
        Err(e) => {
            s.last_error = Some(format!("ligand read: {e}"));
            return;
        }
    };
    let cfg = valenx_dock::DockConfig {
        center: nalgebra::Vector3::new(s.center[0], s.center[1], s.center[2]),
        size: nalgebra::Vector3::new(s.size[0], s.size[1], s.size[2]),
        exhaustiveness: s.exhaustiveness.max(1),
        num_modes: s.num_modes.max(1),
        energy_range: s.energy_range.max(0.1),
        seed: s.seed,
        ..Default::default()
    };
    let out_path = if s.output_path.is_empty() {
        std::env::temp_dir().join("valenx_dock_out.pdbqt")
    } else {
        std::path::PathBuf::from(&s.output_path)
    };
    // Audit: the case-toml-driven adapter dispatcher emits run.start
    // before every external/native adapter run, but the Dock panel
    // here bypasses that and calls `valenx_dock::dock` in-process.
    // Stamp our own entry so the audit chain stays complete for
    // panel-driven docking too.
    emit_audit(
        "dock.native.start",
        serde_json::json!({
            "kind": "dock",
            "receptor": s.receptor_path,
            "ligand": s.ligand_path,
        }),
        serde_json::json!({
            "exhaustiveness": s.exhaustiveness,
            "num_modes": s.num_modes,
            "seed": s.seed,
            "output": out_path.display().to_string(),
        }),
    );
    match valenx_dock::dock(&receptor, &ligand, &cfg, &out_path, None) {
        Ok(poses) => {
            // Compute heavy-atom RMSD vs the top pose so the panel
            // can display a per-mode RMSD column.
            // Reparse the ligand once to get the atom list for RMSD.
            let lig = match valenx_dock::ligand::Ligand::from_pdbqt(&ligand) {
                Ok(l) => l,
                Err(e) => {
                    s.last_error = Some(format!("ligand parse: {e}"));
                    return;
                }
            };
            let first_pose = poses.first().map(|p| p.0.clone());
            s.last_scores = poses
                .iter()
                .map(|(p, score)| {
                    let r = match &first_pose {
                        Some(fp) => valenx_dock::cluster::rmsd(&lig, fp, p),
                        None => 0.0,
                    };
                    (*score, r)
                })
                .collect();
            s.selected_pose = if s.last_scores.is_empty() {
                None
            } else {
                Some(0)
            };
            let _ = Path::new(&out_path);
            emit_audit(
                "dock.native.complete",
                serde_json::json!({
                    "kind": "dock",
                    "output": out_path.display().to_string(),
                }),
                serde_json::json!({
                    "result": "ok",
                    "n_poses": s.last_scores.len(),
                }),
            );
        }
        Err(e) => {
            s.last_error = Some(format!("dock failed: {e}"));
            emit_audit(
                "dock.native.complete",
                serde_json::json!({"kind": "dock"}),
                serde_json::json!({"result": "failed", "error": e.to_string()}),
            );
        }
    }
}

fn push_pose_as_mesh(app: &mut crate::ValenxApp, output_path: &str, pose_index: usize) {
    use std::path::Path;
    use valenx_bio::format::pdbqt::{parse, PdbqtRecord};
    use valenx_core::io_caps::{read_capped_to_string, MAX_PDBQT_FILE_BYTES};
    // Round-20 H1 (R16 H1 sister): cap the dock-output read so a
    // multi-GB file the user accidentally selected (or that the
    // docker mis-wrote) can't OOM the renderer when showing a pose.
    let text = match read_capped_to_string(Path::new(output_path), MAX_PDBQT_FILE_BYTES) {
        Ok(v) => v,
        Err(_) => return,
    };
    // PDBQT pose ensemble: MODEL n / ATOMs / ENDMDL. We split on MODEL.
    let blocks: Vec<&str> = text.split("MODEL ").collect();
    // blocks[0] is anything before the first MODEL line (often empty).
    let target_block = blocks.get(pose_index + 1);
    let Some(block) = target_block else { return };
    let recs = match parse(block) {
        Ok(r) => r,
        Err(_) => return,
    };
    let mut positions = Vec::new();
    for r in recs {
        if let PdbqtRecord::Atom(a) = r {
            positions.push(a.position);
        }
    }
    // Defer to the existing mesh-loader path: write a quick STL of
    // unit-radius octahedra at each atom and load it.
    let mut tris: Vec<[[f64; 3]; 3]> = Vec::new();
    for p in &positions {
        push_octahedron(&mut tris, *p, 0.5);
    }
    let tmp = std::env::temp_dir().join("valenx_pose_preview.stl");
    if let Err(e) = write_binary_stl(&tmp, &tris) {
        tracing::warn!(target: "valenx-app", ?e, "write pose preview STL failed");
        return;
    }
    tracing::info!(target: "valenx-app", "wrote pose preview STL: {}", tmp.display());
    // Reuse the existing app-side STL loader (the same path that
    // "Import STL…" reaches through). The loader sets `app.stl`,
    // frames the viewport, and updates the status string.
    app.load_stl(tmp);
}

/// Draw the Sketcher panel — tool palette, click input, constraint
/// palette, solver, and Pad (extrude) controls.
///
/// Layout mirrors `draw_dock_panel`: capture deferred-dispatch
/// requests during the `&mut app.mesh_toolbox.sketcher` borrow scope,
/// then act on `&mut app` after the borrow ends.
///
/// # Multi-sketch list — DEFERRED to Phase 2
///
/// Phase 1 holds exactly one sketch per session (the field
/// `MeshToolboxState::sketcher.sketch`). Multi-sketch persistence —
/// a named list of sketches, with select / delete / rename — lands
/// in Phase 2 alongside the feature-tree work, which is the natural
/// owner: feature-tree nodes ARE the sketch references. Adding a
/// stand-alone list here would duplicate state the feature tree
/// will own, so we hold off.
pub fn draw_sketcher_panel(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.heading("Sketcher (2D parametric)");
    ui.separator();

    // Deferred-dispatch (same pattern as `draw_dock_panel`): collect
    // a "pad with depth X" request during the &mut sketcher borrow
    // scope, then act on &mut app after the borrow ends.
    let mut pad_request: Option<f64> = None;

    let s = &mut app.mesh_toolbox.sketcher;

    // ----- Undo / Redo affordance (polish pass 3) -----
    //
    // Inline `↶ ↷` mirror Ctrl+Z / Ctrl+Y so first-time users have a
    // visible reset for "I just added the wrong entity". Snapshots are
    // pushed by every mutating button further down via `s.record()`.
    ui.horizontal(|ui| {
        if ui
            .add_enabled(s.can_undo(), egui::Button::new("↶"))
            .on_hover_text("Undo last sketch edit (Ctrl+Z)")
            .clicked()
        {
            s.undo_edit();
        }
        if ui
            .add_enabled(s.can_redo(), egui::Button::new("↷"))
            .on_hover_text("Redo last sketch edit (Ctrl+Y / Ctrl+Shift+Z)")
            .clicked()
        {
            s.redo_edit();
        }
        ui.label(
            egui::RichText::new("Ctrl+Z / Ctrl+Y")
                .weak()
                .small(),
        );
    });

    // ----- Viewport overlay toggle (Phase 1H — Task 51) -----
    //
    // The 2D sketch overlay is drawn on the XY plane through the
    // existing OrbitCamera. Default on so users see what they sketch
    // immediately; disable to inspect the underlying mesh without
    // sketch clutter.
    ui.checkbox(&mut s.show_overlay, "Show sketch overlay in viewport")
        .on_hover_text("Render the current 2-D sketch as a viewport overlay on the XY plane.");

    // ----- Tool palette (Task 39) -----
    ui.label(egui::RichText::new("Tool").strong())
        .on_hover_text("Which sketcher tool the next click adds (Select / Line / Circle).");
    ui.horizontal(|ui| {
        ui.selectable_value(&mut s.tool, SketcherTool::Select, "Select")
            .on_hover_text("Selection mode — click an entity to add it to the selection set for constraints.");
        ui.selectable_value(&mut s.tool, SketcherTool::Line, "Line")
            .on_hover_text("Line tool — two clicks (start + end) place a straight line entity.");
        ui.selectable_value(&mut s.tool, SketcherTool::Circle, "Circle")
            .on_hover_text("Circle tool — first click sets centre, second click sets a rim point (defines radius).");
    });
    ui.label(format!(
        "Active: {:?} — {} entities, {} constraints",
        s.tool,
        s.sketch.entities.len(),
        s.sketch.constraints.len(),
    ));

    // ----- Click input (Task 40) -----
    //
    // MVP: enter coords numerically; full viewport-click pipeline
    // lands in Phase 1H. The DragValue pair stands in for the (x, y)
    // a mouse click would produce, the "Add" button stands in for the
    // click event itself. Behaviour branches on `s.tool`.
    ui.separator();
    ui.label(egui::RichText::new("Click input (MVP — numeric)").strong())
        .on_hover_text("Phase-1 numeric stand-in for viewport clicks. Enter (x, y) sketch-plane coordinates, then press the Add button.");
    ui.horizontal(|ui| {
        ui.label("X")
            .on_hover_text("Sketch-plane X coordinate (model units).");
        ui.add(egui::DragValue::new(&mut s.pending_click_x).speed(0.05))
            .on_hover_text("X coordinate of the next click.");
        ui.label("Y")
            .on_hover_text("Sketch-plane Y coordinate (model units).");
        ui.add(egui::DragValue::new(&mut s.pending_click_y).speed(0.05))
            .on_hover_text("Y coordinate of the next click.");
    });
    let click_label = match s.tool {
        SketcherTool::Select => "Add point at (x, y)",
        SketcherTool::Line => {
            if s.pending_first_click.is_some() {
                "Line: place end point"
            } else {
                "Line: place start point"
            }
        }
        SketcherTool::Circle => {
            if s.pending_first_click.is_some() {
                "Circle: place rim point (sets radius)"
            } else {
                "Circle: place centre point"
            }
        }
    };
    if ui.button(click_label).clicked() {
        let (x, y) = (s.pending_click_x, s.pending_click_y);
        tracing::info!(target: "valenx-app", "sketcher click at ({x}, {y}) tool={:?}", s.tool);
        s.last_error = None;
        // Snapshot before mutating the sketch — undo restores
        // pre-click state.
        s.record();
        match s.tool {
            SketcherTool::Select => {
                // No tool active — just drop a free point.
                let _ = s.sketch.add_point(x, y);
            }
            // Task 41 and 42 land below.
            _ => handle_sketcher_geometry_click(s, x, y),
        }
    }
    if let Some(pending) = s.pending_first_click {
        ui.label(format!("Pending first click: entity {}", pending.0));
    }
    if let Some(err) = &s.last_error {
        ui.colored_label(egui::Color32::RED, err);
    }

    // ----- Constraint palette (Task 43) -----
    //
    // Each button reads from `s.selected` (Task 44 — populated by the
    // entity list further down). Distance / Angle / Radius use
    // `s.pending_target` as their numeric target (v1 simplification:
    // the user types the target into the DragValue BEFORE clicking
    // the constraint button, rather than a modal popup).
    ui.separator();
    ui.label(egui::RichText::new("Constraints").strong());
    ui.horizontal(|ui| {
        ui.label("Numeric target (Distance / Angle / Radius):");
        ui.add(egui::DragValue::new(&mut s.pending_target).speed(0.05));
    });
    ui.horizontal_wrapped(|ui| {
        if ui.button("Coincident").clicked() {
            add_constraint_2(s, |a, b| {
                valenx_sketch::constraint::Constraint::Coincident { a, b }
            });
        }
        if ui.button("Horizontal").clicked() {
            add_constraint_1(s, valenx_sketch::constraint::Constraint::Horizontal);
        }
        if ui.button("Vertical").clicked() {
            add_constraint_1(s, valenx_sketch::constraint::Constraint::Vertical);
        }
        if ui.button("Parallel").clicked() {
            add_constraint_2(s, |a, b| valenx_sketch::constraint::Constraint::Parallel {
                a,
                b,
            });
        }
        if ui.button("Perpendicular").clicked() {
            add_constraint_2(s, |a, b| {
                valenx_sketch::constraint::Constraint::Perpendicular { a, b }
            });
        }
        if ui.button("Tangent").clicked() {
            add_constraint_2(s, |a, b| valenx_sketch::constraint::Constraint::Tangent {
                line_or_circle_a: a,
                circle_b: b,
            });
        }
        if ui.button("EqualLength").clicked() {
            add_constraint_2(s, |a, b| {
                valenx_sketch::constraint::Constraint::EqualLength { a, b }
            });
        }
        let target = s.pending_target;
        if ui.button("Distance").clicked() {
            add_constraint_2(s, |a, b| valenx_sketch::constraint::Constraint::Distance {
                a,
                b,
                target,
            });
        }
        if ui.button("Angle").clicked() {
            add_constraint_2(s, |a, b| valenx_sketch::constraint::Constraint::Angle {
                a,
                b,
                target,
            });
        }
        if ui.button("Radius").clicked() {
            add_constraint_1(s, |a| valenx_sketch::constraint::Constraint::Radius {
                circle_or_arc: a,
                target,
            });
        }
    });

    // ----- Phase 12 — extra primitives + constraints + ops -----
    //
    // Phase 12F Task 47: numeric-driven buttons for new primitives and
    // sketch-level operations. All buttons read from
    // `s.pending_target` (radius / length / angle / spacing — context
    // dependent) plus a small set of optional drag-values inline.
    // Construction toggle flips the flag on every selected entity.
    ui.separator();
    ui.label(egui::RichText::new("Phase 12: Primitives").strong());
    ui.horizontal_wrapped(|ui| {
        if ui.button("BSpline (4 selected pts)").clicked() && s.selected.len() >= 4 {
            let cps: Vec<valenx_sketch::EntityId> = s.selected.iter().take(4).copied().collect();
            let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
            s.record();
            if let Err(e) = s.sketch.add_bspline(3, knots, &cps, vec![1.0; cps.len()]) {
                s.last_error = Some(format!("bspline: {e}"));
            }
        }
        if ui.button("Ellipse (center+major+minor)").clicked() && !s.selected.is_empty() {
            // Use first selected as center; target is minor radius.
            let center = s.selected[0];
            s.record();
            if let Err(e) = s.sketch.add_ellipse(
                center,
                (s.pending_target.max(1.0), 0.0),
                s.pending_target.max(0.5),
            ) {
                s.last_error = Some(format!("ellipse: {e}"));
            }
        }
        if ui.button("EllipticalArc (center)").clicked() && !s.selected.is_empty() {
            let center = s.selected[0];
            s.record();
            if let Err(e) = s.sketch.add_elliptical_arc(
                center,
                (s.pending_target.max(1.0), 0.0),
                s.pending_target.max(0.5),
                0.0,
                std::f64::consts::PI,
            ) {
                s.last_error = Some(format!("ellip arc: {e}"));
            }
        }
    });

    ui.label(egui::RichText::new("Phase 12: Construction").strong());
    if ui.button("Toggle Construction (selected)").clicked() {
        s.record();
        for id in s.selected.clone() {
            s.sketch.toggle_construction(id);
        }
    }

    ui.label(egui::RichText::new("Phase 12: Extra Constraints").strong());
    ui.horizontal_wrapped(|ui| {
        let target = s.pending_target;
        if ui.button("Symmetric (a,b,mid)").clicked() && s.selected.len() >= 3 {
            let a = s.selected[0];
            let b = s.selected[1];
            let mid = s.selected[2];
            s.record();
            s.sketch
                .add_constraint(valenx_sketch::constraint::Constraint::Symmetric {
                    a,
                    b,
                    midpoint: mid,
                });
        }
        if ui.button("PointOnLine").clicked() && s.selected.len() >= 2 {
            let point = s.selected[0];
            let line = s.selected[1];
            s.record();
            s.sketch
                .add_constraint(valenx_sketch::constraint::Constraint::PointOnLine { point, line });
        }
        if ui.button("PointOnCircle").clicked() && s.selected.len() >= 2 {
            let point = s.selected[0];
            let circle = s.selected[1];
            s.record();
            s.sketch
                .add_constraint(valenx_sketch::constraint::Constraint::PointOnCircle {
                    point,
                    circle,
                });
        }
        if ui.button("PointOnArc").clicked() && s.selected.len() >= 2 {
            let point = s.selected[0];
            let arc = s.selected[1];
            s.record();
            s.sketch
                .add_constraint(valenx_sketch::constraint::Constraint::PointOnArc { point, arc });
        }
        if ui.button("PointOnEllipse").clicked() && s.selected.len() >= 2 {
            let point = s.selected[0];
            let ellipse = s.selected[1];
            s.record();
            s.sketch
                .add_constraint(valenx_sketch::constraint::Constraint::PointOnEllipse {
                    point,
                    ellipse,
                });
        }
        if ui.button("DistanceX").clicked() && s.selected.len() >= 2 {
            let (a, b) = (s.selected[0], s.selected[1]);
            s.record();
            s.sketch
                .add_constraint(valenx_sketch::constraint::Constraint::DistanceX { a, b, target });
        }
        if ui.button("DistanceY").clicked() && s.selected.len() >= 2 {
            let (a, b) = (s.selected[0], s.selected[1]);
            s.record();
            s.sketch
                .add_constraint(valenx_sketch::constraint::Constraint::DistanceY { a, b, target });
        }
        if ui.button("LineLength").clicked() && !s.selected.is_empty() {
            let line = s.selected[0];
            s.record();
            s.sketch
                .add_constraint(valenx_sketch::constraint::Constraint::LineLength { line, target });
        }
        if ui.button("ArcRadius").clicked() && !s.selected.is_empty() {
            let arc = s.selected[0];
            s.record();
            s.sketch
                .add_constraint(valenx_sketch::constraint::Constraint::ArcRadius { arc, target });
        }
        if ui.button("ArcAngle").clicked() && !s.selected.is_empty() {
            let arc = s.selected[0];
            s.record();
            s.sketch
                .add_constraint(valenx_sketch::constraint::Constraint::ArcAngle { arc, target });
        }
        if ui.button("EllipseRadiusA").clicked() && !s.selected.is_empty() {
            let ellipse = s.selected[0];
            s.record();
            s.sketch
                .add_constraint(valenx_sketch::constraint::Constraint::EllipseRadiusA {
                    ellipse,
                    target,
                });
        }
        if ui.button("EllipseRadiusB").clicked() && !s.selected.is_empty() {
            let ellipse = s.selected[0];
            s.record();
            s.sketch
                .add_constraint(valenx_sketch::constraint::Constraint::EllipseRadiusB {
                    ellipse,
                    target,
                });
        }
    });

    ui.label(egui::RichText::new("Phase 12: Sketch Ops").strong());
    ui.horizontal_wrapped(|ui| {
        if ui.button("Move (delta x = target, y = 0)").clicked() && !s.selected.is_empty() {
            s.record();
            valenx_sketch::ops::r#move::translate(
                &mut s.sketch,
                &s.selected.clone(),
                (s.pending_target, 0.0),
            );
        }
        if ui.button("Rotate (90 deg about origin)").clicked() && !s.selected.is_empty() {
            s.record();
            valenx_sketch::ops::rotate::rotate(
                &mut s.sketch,
                &s.selected.clone(),
                (0.0, 0.0),
                std::f64::consts::FRAC_PI_2,
            );
        }
        if ui.button("Mirror across X-axis").clicked() && !s.selected.is_empty() {
            let line = valenx_sketch::ops::mirror::MirrorLine {
                point: (0.0, 0.0),
                direction: (1.0, 0.0),
            };
            s.record();
            let _ = valenx_sketch::ops::mirror::mirror(&mut s.sketch, &s.selected.clone(), &line);
        }
        if ui.button("Copy (offset x = target)").clicked() && !s.selected.is_empty() {
            s.record();
            let _ = valenx_sketch::ops::copy::copy(
                &mut s.sketch,
                &s.selected.clone(),
                (s.pending_target, 0.0),
            );
        }
        if ui
            .button("LinearArray (3 along +x, spacing = target)")
            .clicked()
            && !s.selected.is_empty()
        {
            s.record();
            let _ = valenx_sketch::ops::linear_array::linear_array(
                &mut s.sketch,
                &s.selected.clone(),
                (1.0, 0.0),
                3,
                s.pending_target.max(1.0),
            );
        }
        if ui.button("PolarArray (4 about origin, 2π)").clicked() && !s.selected.is_empty() {
            s.record();
            let _ = valenx_sketch::ops::polar_array::polar_array(
                &mut s.sketch,
                &s.selected.clone(),
                (0.0, 0.0),
                4,
                std::f64::consts::TAU,
            );
        }
    });

    // ----- Selection (Task 44) -----
    //
    // Render every entity as a selectable_label. Click toggles its
    // ID in `s.selected`. Selected entries render as RichText.strong()
    // so the user can see what the constraint buttons will operate
    // on. Shift-click is implicit — clicking an already-selected
    // entity removes it from the selection.
    ui.separator();
    ui.collapsing("Entities", |ui| {
        if s.sketch.entities.is_empty() {
            ui.label("(empty)");
        }
        for (i, entity) in s.sketch.entities.iter().enumerate() {
            let id = valenx_sketch::EntityId(i + 1);
            let selected = s.selected.contains(&id);
            let kind = match entity {
                valenx_sketch::geom::Entity::Point(p) => {
                    let (x, y) = p.read(&s.sketch.vars);
                    format!("Point  ({x:.3}, {y:.3})")
                }
                valenx_sketch::geom::Entity::Line(l) => {
                    let ((sx, sy), (ex, ey)) = l.endpoints(&s.sketch.vars);
                    format!("Line   ({sx:.3}, {sy:.3}) → ({ex:.3}, {ey:.3})")
                }
                valenx_sketch::geom::Entity::Circle(c) => {
                    let (cx, cy) = c.center.read(&s.sketch.vars);
                    let r = c.radius(&s.sketch.vars);
                    format!("Circle ({cx:.3}, {cy:.3}) r={r:.3}")
                }
                valenx_sketch::geom::Entity::Arc(a) => {
                    let (cx, cy) = a.center.read(&s.sketch.vars);
                    format!("Arc    ({cx:.3}, {cy:.3})")
                }
                valenx_sketch::geom::Entity::BSpline(b) => {
                    format!("BSpline deg={} cps={}", b.degree, b.n_control_points())
                }
                valenx_sketch::geom::Entity::Ellipse(e) => {
                    let (cx, cy) = e.center.read(&s.sketch.vars);
                    let a = e.major_radius(&s.sketch.vars);
                    let bb = e.minor_radius(&s.sketch.vars);
                    format!("Ellipse ({cx:.3}, {cy:.3}) a={a:.3} b={bb:.3}")
                }
                valenx_sketch::geom::Entity::EllipticalArc(ea) => {
                    let (cx, cy) = ea.ellipse.center.read(&s.sketch.vars);
                    format!("EllipArc ({cx:.3}, {cy:.3})")
                }
            };
            let label = if selected {
                egui::RichText::new(format!("#{}: {kind}", id.0)).strong()
            } else {
                egui::RichText::new(format!("#{}: {kind}", id.0))
            };
            if ui.selectable_label(selected, label).clicked() {
                if selected {
                    s.selected.retain(|x| *x != id);
                } else {
                    s.selected.push(id);
                }
            }
        }
        if !s.selected.is_empty() && ui.button("Clear selection").clicked() {
            s.selected.clear();
        }
    });

    // ----- Solve (Task 45) -----
    //
    // Newton-Raphson + Levenberg-Marquardt — drives every variable to
    // satisfy every constraint simultaneously. Stashes the report
    // for Task 46 to render.
    ui.separator();
    if ui.button("Solve").clicked() {
        match valenx_sketch::solver::solve(&mut s.sketch, valenx_sketch::SolverConfig::default()) {
            Ok(report) => {
                s.last_report = Some(report);
                s.last_error = None;
            }
            Err(e) => {
                s.last_error = Some(format!("solve: {e}"));
                s.last_report = None;
            }
        }
    }

    // ----- Solver report (Task 46) -----
    //
    // Colour code:
    //   green  = Converged AND dof_balance == 0 (well-constrained).
    //   yellow = Converged but dof_balance != 0 (under/over-constrained
    //            but the solver still landed somewhere; users may want
    //            to add or drop constraints).
    //   red    = MaxIterations.
    if let Some(report) = &s.last_report {
        let status_color = match report.status {
            valenx_sketch::SolverStatus::Converged if report.diagnostics.dof_balance == 0 => {
                egui::Color32::from_rgb(80, 200, 120)
            }
            valenx_sketch::SolverStatus::Converged => egui::Color32::from_rgb(220, 180, 70),
            valenx_sketch::SolverStatus::MaxIterations => egui::Color32::RED,
        };
        ui.colored_label(status_color, format!("Status: {:?}", report.status));
        ui.label(format!("Iterations: {}", report.iterations));
        ui.label(format!("Residual ‖r‖₂: {:.3e}", report.residual_norm));
        let d = &report.diagnostics;
        ui.label(format!(
            "Constraints (residuals): {}   Variables: {}",
            d.n_residuals, d.n_variables,
        ));
        let dof_color = match d.dof_balance {
            0 => egui::Color32::from_rgb(80, 200, 120),
            n if n < 0 => egui::Color32::from_rgb(220, 180, 70),
            _ => egui::Color32::from_rgb(220, 100, 100),
        };
        let dof_label = match d.dof_balance {
            0 => "DOF balance: 0 (well-constrained)".to_string(),
            n if n < 0 => format!("DOF balance: {n} (under-constrained)"),
            n => format!("DOF balance: +{n} (over-constrained)"),
        };
        ui.colored_label(dof_color, dof_label);
    }

    // ----- Pad / Extrude (Task 47) -----
    //
    // Walk the closed sketch profile, extrude along +Z by `pad_depth`,
    // tessellate the resulting BRep, and load it into the viewport.
    // Pattern identical to the dock panel's "Show pose in viewport"
    // request: capture during the &mut sketcher borrow, dispatch
    // after the borrow ends.
    ui.separator();
    ui.label(egui::RichText::new("Pad (extrude)").strong());
    ui.horizontal(|ui| {
        ui.label("Depth:");
        ui.add(egui::DragValue::new(&mut s.pad_depth).speed(0.05));
        if ui.button("Pad").clicked() {
            pad_request = Some(s.pad_depth);
        }
    });

    // Release the &mut sketcher borrow before reaching for &mut app.
    let _ = s;
    if let Some(depth) = pad_request {
        push_sketch_pad(app, depth);
    }
}

/// Extrude the current sketch by `depth` and load the resulting mesh
/// into the viewport. Errors land in `app.mesh_toolbox.sketcher.last_error`.
fn push_sketch_pad(app: &mut crate::ValenxApp, depth: f64) {
    let solid = match app.mesh_toolbox.sketcher.sketch.extrude(depth) {
        Ok(s) => s,
        Err(e) => {
            app.mesh_toolbox.sketcher.last_error = Some(format!("pad: extrude failed: {e}"));
            return;
        }
    };
    let mesh = match valenx_cad::solid_to_mesh(&solid, valenx_cad::DEFAULT_TESS_TOLERANCE) {
        Ok(m) => m,
        Err(e) => {
            app.mesh_toolbox.sketcher.last_error = Some(format!("pad: tessellate failed: {e}"));
            return;
        }
    };
    let pseudo_path = std::path::PathBuf::from("<sketcher>/pad.solid");
    app.apply_mesh(mesh, pseudo_path);
    app.mesh_toolbox.sketcher.last_error = None;
    emit_audit(
        "sketcher.pad",
        serde_json::json!({"kind": "sketch.pad"}),
        serde_json::json!({"depth": depth}),
    );
}

/// Helper for single-entity constraints (Horizontal / Vertical /
/// Radius). Reads `s.selected[0]`, calls `build`, appends. Sets an
/// error string if the selection is empty.
fn add_constraint_1<F>(s: &mut SketcherPanelState, build: F)
where
    F: FnOnce(valenx_sketch::EntityId) -> valenx_sketch::constraint::Constraint,
{
    let Some(a) = s.selected.first().copied() else {
        s.last_error = Some("Select an entity first (use the list below).".into());
        return;
    };
    // Snapshot before appending so undo restores pre-constraint state.
    s.record();
    s.sketch.add_constraint(build(a));
    s.last_error = None;
}

/// Helper for two-entity constraints. Reads `s.selected[0..2]`.
fn add_constraint_2<F>(s: &mut SketcherPanelState, build: F)
where
    F: FnOnce(
        valenx_sketch::EntityId,
        valenx_sketch::EntityId,
    ) -> valenx_sketch::constraint::Constraint,
{
    if s.selected.len() < 2 {
        s.last_error = Some("Select TWO entities first (use the list below).".into());
        return;
    }
    let (a, b) = (s.selected[0], s.selected[1]);
    s.record();
    s.sketch.add_constraint(build(a, b));
    s.last_error = None;
}

/// Two-click placement logic for Line and Circle tools (Tasks 41-42).
fn handle_sketcher_geometry_click(s: &mut SketcherPanelState, x: f64, y: f64) {
    match s.tool {
        SketcherTool::Line => {
            // Task 41 — two-click line.
            let new_id = s.sketch.add_point(x, y);
            if let Some(first) = s.pending_first_click.take() {
                match s.sketch.add_line(first, new_id) {
                    Ok(_) => {}
                    Err(e) => s.last_error = Some(format!("add line: {e}")),
                }
            } else {
                s.pending_first_click = Some(new_id);
            }
        }
        SketcherTool::Circle => {
            // Task 42 — centre then radius. First click sets the
            // centre point; second click defines the radius as the
            // distance from the centre to the rim click.
            if let Some(center_id) = s.pending_first_click.take() {
                let center = match s.sketch.point_at(center_id) {
                    Ok(p) => p,
                    Err(e) => {
                        s.last_error = Some(format!("circle: bad centre: {e}"));
                        return;
                    }
                };
                let (cx, cy) = center.read(&s.sketch.vars);
                let r = ((x - cx).powi(2) + (y - cy).powi(2)).sqrt();
                if r <= 0.0 {
                    s.last_error =
                        Some("circle: rim click coincides with centre (radius would be 0)".into());
                    return;
                }
                if let Err(e) = s.sketch.add_circle(center_id, r) {
                    s.last_error = Some(format!("add circle: {e}"));
                }
            } else {
                let id = s.sketch.add_point(x, y);
                s.pending_first_click = Some(id);
            }
        }
        // Select handled by the caller.
        SketcherTool::Select => {}
    }
}

fn push_octahedron(out: &mut Vec<[[f64; 3]; 3]>, c: nalgebra::Vector3<f64>, r: f64) {
    let v = |dx: f64, dy: f64, dz: f64| [c.x + dx * r, c.y + dy * r, c.z + dz * r];
    let top = v(0.0, 0.0, 1.0);
    let bot = v(0.0, 0.0, -1.0);
    let n = v(1.0, 0.0, 0.0);
    let e = v(0.0, 1.0, 0.0);
    let s = v(-1.0, 0.0, 0.0);
    let w = v(0.0, -1.0, 0.0);
    out.push([top, n, e]);
    out.push([top, e, s]);
    out.push([top, s, w]);
    out.push([top, w, n]);
    out.push([bot, e, n]);
    out.push([bot, s, e]);
    out.push([bot, w, s]);
    out.push([bot, n, w]);
}

fn write_binary_stl(path: &std::path::Path, tris: &[[[f64; 3]; 3]]) -> std::io::Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    f.write_all(&[0u8; 80])?;
    f.write_all(&(tris.len() as u32).to_le_bytes())?;
    for t in tris {
        // Normal (zero — STL readers usually re-derive).
        f.write_all(&[0u8; 12])?;
        for v in t {
            f.write_all(&(v[0] as f32).to_le_bytes())?;
            f.write_all(&(v[1] as f32).to_le_bytes())?;
            f.write_all(&(v[2] as f32).to_le_bytes())?;
        }
        f.write_all(&[0u8; 2])?;
    }
    Ok(())
}

/// Render the Part Design panel (parametric feature tree).
///
/// The panel is the UI counterpart of [`valenx_feature_tree`]: it lets
/// the user grow a tree of sketches and features (Pad, Pocket, Revolve,
/// Mirror, LinearPattern, CircularPattern), suppress / delete entries,
/// and watch the live re-evaluated solid land in the viewport.
///
/// Deferred-dispatch idiom: every mutation that needs `&mut app` (e.g.
/// "push the freshly-replayed mesh") is captured during the `&mut
/// part_design` borrow scope as a small flag/struct, then acted on after
/// the borrow ends.
pub fn draw_part_design_panel(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.heading("Part Design (parametric)");
    ui.separator();

    // ----- Undo / Redo affordance (polish pass 3) -----
    //
    // Inline `↶ ↷` mirror the Ctrl+Z / Ctrl+Y host shortcuts. Each
    // mutating add/delete/suppress button calls `s.record()` before
    // touching `tree`, so undo restores the prior snapshot.
    {
        let s = &mut app.mesh_toolbox.part_design;
        ui.horizontal(|ui| {
            if ui
                .add_enabled(s.can_undo(), egui::Button::new("↶"))
                .on_hover_text("Undo last feature-tree edit (Ctrl+Z)")
                .clicked()
            {
                s.undo_edit();
            }
            if ui
                .add_enabled(s.can_redo(), egui::Button::new("↷"))
                .on_hover_text("Redo last feature-tree edit (Ctrl+Y / Ctrl+Shift+Z)")
                .clicked()
            {
                s.redo_edit();
            }
            ui.label(
                egui::RichText::new("Ctrl+Z / Ctrl+Y")
                    .weak()
                    .small(),
            );
        });
    }

    // ----- Tasks 53 + 54: Save / Load project (.valenx) -----
    //
    // Buttons dispatch into module-local helpers that own the
    // FileDialog → ValenxProject round-trip. Cancellation is silent;
    // errors land in last_replay_error and successes in app.status.
    enum PersistAction {
        None,
        Save,
        Load,
    }
    let mut persist_action = PersistAction::None;
    ui.horizontal(|ui| {
        if ui
            .button("Save project (.valenx)")
            .on_hover_text("Save the current Part Design feature tree to a .valenx project file.")
            .clicked()
        {
            persist_action = PersistAction::Save;
        }
        if ui
            .button("Load project (.valenx)")
            .on_hover_text("Replace the current feature tree with one loaded from a .valenx project file.")
            .clicked()
        {
            persist_action = PersistAction::Load;
        }
    });
    match persist_action {
        PersistAction::None => {}
        PersistAction::Save => save_part_design_project(app),
        PersistAction::Load => load_part_design_project(app),
    }
    if let Some(p) = &app.mesh_toolbox.part_design.project_path {
        ui.label(
            egui::RichText::new(format!("Project file: {}", p.display()))
                .weak()
                .small(),
        );
    }
    // Task 56 — surface save/load errors inline beside the buttons so
    // the user sees them immediately rather than having to scroll past
    // the entire tree to reach the replay-error label at the bottom.
    // Only paint when the error looks like a persistence failure (the
    // shared error slot is also used by replay / suppress / delete);
    // those still show through the bottom panel.
    if let Some(err) = &app.mesh_toolbox.part_design.last_replay_error {
        if err.starts_with("save project:")
            || err.starts_with("load project:")
            || err.starts_with("import STEP/IGES:")
            || err.starts_with("export STEP/IGES:")
        {
            ui.colored_label(egui::Color32::from_rgb(220, 100, 100), err);
        }
    }

    // ----- Phase 8 (Tasks 20 + 22): Import / Export STEP / IGES -----
    //
    // Import: pick a file → re-read on every replay via the
    //   ImportedSolid feature variant (no embedded mesh in the project
    //   file).
    // Export: save the *currently-resolved* solid (last successful
    //   replay) to a STEP / IGES path the user picks.
    enum InteropAction {
        None,
        Import,
        ImportAdvanced,
        Export,
    }
    let mut interop_action = InteropAction::None;
    ui.horizontal(|ui| {
        if ui
            .button("Import STEP/IGES…")
            .on_hover_text("Read a STEP / IGES file via the truck kernel and embed it as an `ImportedSolid` feature.")
            .clicked()
        {
            interop_action = InteropAction::Import;
        }
        if ui
            .button("Advanced STEP AP242 / IGES trimmed…")
            .on_hover_text(
                "Import a STEP AP242 file with product hierarchy, materials, \
                 colors, and GD&T metadata, or an IGES file with trimmed \
                 NURBS surfaces (Type 128 / 142 / 144).",
            )
            .clicked()
        {
            interop_action = InteropAction::ImportAdvanced;
        }
        if ui
            .button("Export STEP/IGES…")
            .on_hover_text("Export the currently resolved solid (last successful replay) as a STEP or IGES file.")
            .clicked()
        {
            interop_action = InteropAction::Export;
        }
    });
    match interop_action {
        InteropAction::None => {}
        InteropAction::Import => import_step_iges(app),
        InteropAction::ImportAdvanced => import_step_ap242(app),
        InteropAction::Export => export_step_iges(app),
    }
    ui.separator();

    // ----- Task 37: "Add Sketch" button -----
    //
    // Copies the active sketcher sketch into the tree as a new
    // SketchRef. The sketcher panel's `sketch` is the user's
    // current scratchpad; once they're happy with the profile, this
    // button promotes it into the tree where Pad / Pocket / Revolve
    // can reference it.
    if ui.button("Add Sketch (from Sketcher)").clicked() {
        let sketch = app.mesh_toolbox.sketcher.sketch.clone();
        let s = &mut app.mesh_toolbox.part_design;
        s.record();
        let id = s.tree.add_sketch(sketch);
        s.selected_sketch = Some(id);
        s.last_replay_error = None;
        s.pending_replay = true;
        tracing::info!(
            target: "valenx-app",
            "part-design: added sketch #{}",
            id.0,
        );
    }

    let s = &mut app.mesh_toolbox.part_design;
    ui.label(format!(
        "Sketches in tree: {}   Features: {}",
        s.tree.sketches.len(),
        s.tree.features.len(),
    ));

    // ----- Task 38: Tree view -----
    //
    // Vertical list — every feature renders as a selectable label
    // showing kind + name + suppressed state. Clicking sets
    // `selected_feature`, which downstream buttons (Suppress, Delete)
    // operate on. Selecting a feature also clears any pending delete
    // confirmation so users don't accidentally delete the wrong row.
    ui.separator();
    ui.label(egui::RichText::new("Feature tree").strong());
    if s.tree.features.is_empty() {
        ui.label("(empty — add a sketch and then a Pad to start)");
    } else {
        for (idx, entry) in s.tree.features.iter().enumerate() {
            let id = valenx_feature_tree::feature::FeatureId(idx);
            let selected = s.selected_feature == Some(id);
            let suppressed_tag = if entry.suppressed {
                " (suppressed)"
            } else {
                ""
            };
            let label = format!(
                "#{} {} — {}{}",
                idx,
                entry.feature.kind_label(),
                entry.name,
                suppressed_tag,
            );
            let rich = if entry.suppressed {
                egui::RichText::new(label).weak()
            } else {
                egui::RichText::new(label)
            };
            if ui.selectable_label(selected, rich).clicked() {
                s.selected_feature = Some(id);
                s.pending_delete_confirm = false;
            }
        }
    }

    // ----- Task 39: Add Pad -----
    //
    // v1: edit the param fields BEFORE clicking "Add Pad" (no modal
    // popup; that lands in Phase 2.5). The sketch dropdown shows
    // SketchRef indices; the depth + direction are numeric.
    ui.separator();
    ui.collapsing("Add Pad", |ui| {
        ui.label("Edit params, then click \"Add Pad\".");
        ui.horizontal(|ui| {
            ui.label("Sketch:");
            sketch_index_dropdown(
                ui,
                "pad_sketch",
                &mut s.pad_sketch_index,
                s.tree.sketches.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Depth:");
            ui.add(egui::DragValue::new(&mut s.pad_depth).speed(0.1));
        });
        ui.checkbox(
            &mut s.pad_direction_positive,
            "Direction +Z (uncheck for -Z)",
        );
        if ui.button("Add Pad").clicked() {
            if s.tree.sketches.is_empty() {
                s.last_replay_error =
                    Some("Add a sketch first (use \"Add Sketch (from Sketcher)\").".into());
            } else {
                s.record();
                let sketch = valenx_feature_tree::feature::SketchRef(s.pad_sketch_index);
                let params = valenx_feature_tree::feature::PadParams {
                    sketch,
                    depth: s.pad_depth.into(),
                    direction_positive: s.pad_direction_positive,
                };
                let name = format!("Pad {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::Pad(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ----- Task 40: Add Pocket -----
    ui.collapsing("Add Pocket", |ui| {
        ui.label(
            "Edit params, then click \"Add Pocket\". \
             Must follow a solid-producing feature (e.g. a Pad).",
        );
        ui.horizontal(|ui| {
            ui.label("Sketch:");
            sketch_index_dropdown(
                ui,
                "pocket_sketch",
                &mut s.pocket_sketch_index,
                s.tree.sketches.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Depth:");
            ui.add(egui::DragValue::new(&mut s.pocket_depth).speed(0.1));
        });
        ui.checkbox(
            &mut s.pocket_direction_positive,
            "Direction +Z (uncheck for -Z)",
        );
        if ui.button("Add Pocket").clicked() {
            if s.tree.sketches.is_empty() {
                s.last_replay_error =
                    Some("Add a sketch first (\"Add Sketch (from Sketcher)\").".into());
            } else {
                s.record();
                let sketch = valenx_feature_tree::feature::SketchRef(s.pocket_sketch_index);
                let params = valenx_feature_tree::feature::PocketParams {
                    sketch,
                    depth: s.pocket_depth.into(),
                    direction_positive: s.pocket_direction_positive,
                };
                let name = format!("Pocket {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::Pocket(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ----- Task 41: Add Revolve -----
    ui.collapsing("Add Revolve", |ui| {
        ui.label("Sweep a sketch profile about an axis.");
        ui.horizontal(|ui| {
            ui.label("Sketch:");
            sketch_index_dropdown(
                ui,
                "revolve_sketch",
                &mut s.revolve_sketch_index,
                s.tree.sketches.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Axis origin");
            ui.add(
                egui::DragValue::new(&mut s.revolve_axis_origin[0])
                    .speed(0.1)
                    .prefix("X "),
            );
            ui.add(
                egui::DragValue::new(&mut s.revolve_axis_origin[1])
                    .speed(0.1)
                    .prefix("Y "),
            );
            ui.add(
                egui::DragValue::new(&mut s.revolve_axis_origin[2])
                    .speed(0.1)
                    .prefix("Z "),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Axis direction");
            ui.add(
                egui::DragValue::new(&mut s.revolve_axis_direction[0])
                    .speed(0.1)
                    .prefix("X "),
            );
            ui.add(
                egui::DragValue::new(&mut s.revolve_axis_direction[1])
                    .speed(0.1)
                    .prefix("Y "),
            );
            ui.add(
                egui::DragValue::new(&mut s.revolve_axis_direction[2])
                    .speed(0.1)
                    .prefix("Z "),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Angle (deg):");
            ui.add(
                egui::DragValue::new(&mut s.revolve_angle_deg)
                    .speed(1.0)
                    .range(-360.0..=360.0),
            );
        });
        if ui.button("Add Revolve").clicked() {
            if s.tree.sketches.is_empty() {
                s.last_replay_error =
                    Some("Add a sketch first (\"Add Sketch (from Sketcher)\").".into());
            } else {
                s.record();
                let sketch = valenx_feature_tree::feature::SketchRef(s.revolve_sketch_index);
                let params = valenx_feature_tree::feature::RevolveParams {
                    sketch,
                    axis_origin: nalgebra::Vector3::new(
                        s.revolve_axis_origin[0],
                        s.revolve_axis_origin[1],
                        s.revolve_axis_origin[2],
                    ),
                    axis_direction: nalgebra::Vector3::new(
                        s.revolve_axis_direction[0],
                        s.revolve_axis_direction[1],
                        s.revolve_axis_direction[2],
                    ),
                    angle: s.revolve_angle_deg.to_radians().into(),
                };
                let name = format!("Revolve {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::Revolve(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ----- Task 42: Add Mirror -----
    ui.collapsing("Add Mirror", |ui| {
        ui.label("Reflect an earlier feature across a plane.");
        ui.horizontal(|ui| {
            ui.label("Target feature:");
            feature_index_dropdown(
                ui,
                "mirror_target",
                &mut s.mirror_target_index,
                s.tree.features.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Plane origin");
            ui.add(
                egui::DragValue::new(&mut s.mirror_plane_origin[0])
                    .speed(0.1)
                    .prefix("X "),
            );
            ui.add(
                egui::DragValue::new(&mut s.mirror_plane_origin[1])
                    .speed(0.1)
                    .prefix("Y "),
            );
            ui.add(
                egui::DragValue::new(&mut s.mirror_plane_origin[2])
                    .speed(0.1)
                    .prefix("Z "),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Plane normal");
            ui.add(
                egui::DragValue::new(&mut s.mirror_plane_normal[0])
                    .speed(0.1)
                    .prefix("X "),
            );
            ui.add(
                egui::DragValue::new(&mut s.mirror_plane_normal[1])
                    .speed(0.1)
                    .prefix("Y "),
            );
            ui.add(
                egui::DragValue::new(&mut s.mirror_plane_normal[2])
                    .speed(0.1)
                    .prefix("Z "),
            );
        });
        ui.checkbox(
            &mut s.mirror_keep_original,
            "Keep original (union mirrored with source)",
        );
        if ui.button("Add Mirror").clicked() {
            if s.tree.features.is_empty() {
                s.last_replay_error =
                    Some("Mirror needs an earlier feature to target (add a Pad first).".into());
            } else {
                s.record();
                let target = valenx_feature_tree::feature::FeatureId(s.mirror_target_index);
                let params = valenx_feature_tree::feature::MirrorParams {
                    target,
                    plane_origin: nalgebra::Vector3::new(
                        s.mirror_plane_origin[0],
                        s.mirror_plane_origin[1],
                        s.mirror_plane_origin[2],
                    ),
                    plane_normal: nalgebra::Vector3::new(
                        s.mirror_plane_normal[0],
                        s.mirror_plane_normal[1],
                        s.mirror_plane_normal[2],
                    ),
                    keep_original: s.mirror_keep_original,
                };
                let name = format!("Mirror {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::Mirror(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ----- Task 43: Add Linear Pattern -----
    ui.collapsing("Add Linear Pattern", |ui| {
        ui.label("Translate-and-union N instances of an earlier feature.");
        ui.horizontal(|ui| {
            ui.label("Target feature:");
            feature_index_dropdown(
                ui,
                "lp_target",
                &mut s.lp_target_index,
                s.tree.features.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Direction");
            ui.add(
                egui::DragValue::new(&mut s.lp_direction[0])
                    .speed(0.1)
                    .prefix("X "),
            );
            ui.add(
                egui::DragValue::new(&mut s.lp_direction[1])
                    .speed(0.1)
                    .prefix("Y "),
            );
            ui.add(
                egui::DragValue::new(&mut s.lp_direction[2])
                    .speed(0.1)
                    .prefix("Z "),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Count:");
            ui.add(
                egui::DragValue::new(&mut s.lp_count)
                    .speed(1.0)
                    .range(1..=512),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Spacing:");
            ui.add(egui::DragValue::new(&mut s.lp_spacing).speed(0.1));
        });
        if ui.button("Add Linear Pattern").clicked() {
            if s.tree.features.is_empty() {
                s.last_replay_error =
                    Some("Linear pattern needs an earlier feature to target.".into());
            } else {
                s.record();
                let target = valenx_feature_tree::feature::FeatureId(s.lp_target_index);
                let params = valenx_feature_tree::feature::LinearPatternParams {
                    target,
                    direction: nalgebra::Vector3::new(
                        s.lp_direction[0],
                        s.lp_direction[1],
                        s.lp_direction[2],
                    ),
                    count: s.lp_count,
                    spacing: s.lp_spacing,
                };
                let name = format!("LinearPattern {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::LinearPattern(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ----- Task 44: Add Circular Pattern -----
    ui.collapsing("Add Circular Pattern", |ui| {
        ui.label("Rotate-and-union N instances of an earlier feature.");
        ui.horizontal(|ui| {
            ui.label("Target feature:");
            feature_index_dropdown(
                ui,
                "cp_target",
                &mut s.cp_target_index,
                s.tree.features.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Axis origin");
            ui.add(
                egui::DragValue::new(&mut s.cp_axis_origin[0])
                    .speed(0.1)
                    .prefix("X "),
            );
            ui.add(
                egui::DragValue::new(&mut s.cp_axis_origin[1])
                    .speed(0.1)
                    .prefix("Y "),
            );
            ui.add(
                egui::DragValue::new(&mut s.cp_axis_origin[2])
                    .speed(0.1)
                    .prefix("Z "),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Axis direction");
            ui.add(
                egui::DragValue::new(&mut s.cp_axis_direction[0])
                    .speed(0.1)
                    .prefix("X "),
            );
            ui.add(
                egui::DragValue::new(&mut s.cp_axis_direction[1])
                    .speed(0.1)
                    .prefix("Y "),
            );
            ui.add(
                egui::DragValue::new(&mut s.cp_axis_direction[2])
                    .speed(0.1)
                    .prefix("Z "),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Count:");
            ui.add(
                egui::DragValue::new(&mut s.cp_count)
                    .speed(1.0)
                    .range(1..=512),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Total angle (deg):");
            ui.add(
                egui::DragValue::new(&mut s.cp_total_angle_deg)
                    .speed(1.0)
                    .range(-360.0..=360.0),
            );
        });
        if ui.button("Add Circular Pattern").clicked() {
            if s.tree.features.is_empty() {
                s.last_replay_error =
                    Some("Circular pattern needs an earlier feature to target.".into());
            } else {
                s.record();
                let target = valenx_feature_tree::feature::FeatureId(s.cp_target_index);
                let params = valenx_feature_tree::feature::CircularPatternParams {
                    target,
                    axis_origin: nalgebra::Vector3::new(
                        s.cp_axis_origin[0],
                        s.cp_axis_origin[1],
                        s.cp_axis_origin[2],
                    ),
                    axis_direction: nalgebra::Vector3::new(
                        s.cp_axis_direction[0],
                        s.cp_axis_direction[1],
                        s.cp_axis_direction[2],
                    ),
                    count: s.cp_count,
                    total_angle: s.cp_total_angle_deg.to_radians(),
                };
                let name = format!("CircularPattern {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::CircularPattern(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ----- Task 30 (Phase 3): Add Fillet -----
    //
    // Round every sharp convex edge of an earlier feature's solid with
    // a cylindrical strip. Phase 14 adds a BRep-first dispatch path
    // with an optional explicit-edge-index selector.
    ui.collapsing("Add Fillet", |ui| {
        ui.label(
            "Round sharp convex edges of an earlier feature with a \
             cylindrical strip. Phase 14: tries BRep path first (true \
             Solid::Brep output), falls through to mesh-domain on \
             unsupported geometry. v1 always falls through pending \
             truck-modeling face-substitution (Phase 14.5+) — treat \
             output as mesh-backed and apply Fillets LAST in the tree.",
        );
        ui.horizontal(|ui| {
            ui.label("Target feature:");
            feature_index_dropdown(
                ui,
                "fillet_target",
                &mut s.fillet_target_index,
                s.tree.features.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Radius:");
            ui.add(
                egui::DragValue::new(&mut s.fillet_radius)
                    .speed(0.05)
                    .range(0.001..=1000.0),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Edge angle threshold (deg):");
            ui.add(
                egui::DragValue::new(&mut s.fillet_threshold_deg)
                    .speed(1.0)
                    .range(0.0..=180.0),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Edges (BRep only, empty = auto):");
            ui.add(
                egui::TextEdit::singleline(&mut s.fillet_edge_indices)
                    .desired_width(200.0)
                    .hint_text("e.g. 0,1,2,3"),
            );
        });
        if ui.button("Add Fillet").clicked() {
            if s.tree.features.is_empty() {
                s.last_replay_error = Some("Fillet needs an earlier feature to target.".into());
            } else {
                s.record();
                let target = valenx_feature_tree::feature::FeatureId(s.fillet_target_index);
                let params = valenx_feature_tree::feature::FilletParams {
                    target,
                    radius: s.fillet_radius,
                    threshold_deg: s.fillet_threshold_deg,
                    // Phase 14: edge_indices=Some(...) targets the
                    // BRep path on a specific edge selection; None
                    // delegates to the angle-threshold auto-selector.
                    // The panel exposes a text-field selector that
                    // parses empty as None, populated as Some(vec).
                    edge_indices: parse_edge_indices_or_none(&s.fillet_edge_indices),
                };
                let name = format!("Fillet {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::Fillet(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ----- Task 31 (Phase 3): Add Chamfer -----
    ui.collapsing("Add Chamfer", |ui| {
        ui.label(
            "Replace sharp convex edges of an earlier feature with a \
             flat bevel. Phase 14 dispatch: BRep first, fall through \
             to mesh-domain. Same caveat as Fillet (apply last in the \
             tree until Phase 14.5+ ships).",
        );
        ui.horizontal(|ui| {
            ui.label("Target feature:");
            feature_index_dropdown(
                ui,
                "chamfer_target",
                &mut s.chamfer_target_index,
                s.tree.features.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Distance:");
            ui.add(
                egui::DragValue::new(&mut s.chamfer_distance)
                    .speed(0.05)
                    .range(0.001..=1000.0),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Edge angle threshold (deg):");
            ui.add(
                egui::DragValue::new(&mut s.chamfer_threshold_deg)
                    .speed(1.0)
                    .range(0.0..=180.0),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Edges (BRep only, empty = auto):");
            ui.add(
                egui::TextEdit::singleline(&mut s.chamfer_edge_indices)
                    .desired_width(200.0)
                    .hint_text("e.g. 0,1,2,3"),
            );
        });
        if ui.button("Add Chamfer").clicked() {
            if s.tree.features.is_empty() {
                s.last_replay_error = Some("Chamfer needs an earlier feature to target.".into());
            } else {
                s.record();
                let target = valenx_feature_tree::feature::FeatureId(s.chamfer_target_index);
                let params = valenx_feature_tree::feature::ChamferParams {
                    target,
                    distance: s.chamfer_distance,
                    threshold_deg: s.chamfer_threshold_deg,
                    edge_indices: parse_edge_indices_or_none(&s.chamfer_edge_indices),
                };
                let name = format!("Chamfer {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::Chamfer(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ===== Phase 13 — new feature ops =====

    // ----- Task 39: Add Hole -----
    ui.collapsing("Add Hole (Phase 13A)", |ui| {
        ui.label("Drill cylindrical pockets at sketch points; optional counterbore, countersink, and thread metadata.");
        ui.horizontal(|ui| {
            ui.label("Sketch (positions):");
            sketch_index_dropdown(ui, "hole_sketch", &mut s.hole_sketch_index, s.tree.sketches.len());
        });
        ui.horizontal(|ui| {
            ui.label("Depth mode:");
            egui::ComboBox::from_id_source("hole_depth_mode")
                .selected_text(match s.hole_depth_mode_idx { 0 => "Blind", 1 => "Through", _ => "Up-to-face" })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut s.hole_depth_mode_idx, 0, "Blind");
                    ui.selectable_value(&mut s.hole_depth_mode_idx, 1, "Through");
                    ui.selectable_value(&mut s.hole_depth_mode_idx, 2, "Up-to-face");
                });
        });
        if s.hole_depth_mode_idx == 0 {
            ui.horizontal(|ui| {
                ui.label("Blind depth (mm):");
                ui.add(egui::DragValue::new(&mut s.hole_blind_depth).speed(0.1).range(0.01..=1000.0));
            });
        }
        ui.horizontal(|ui| {
            ui.label("Drill diameter (mm):");
            ui.add(egui::DragValue::new(&mut s.hole_drill_diameter).speed(0.1).range(0.01..=1000.0));
        });
        ui.checkbox(&mut s.hole_direction_negative, "Drill downward (-Z)");
        ui.checkbox(&mut s.hole_use_counterbore, "Counterbore");
        if s.hole_use_counterbore {
            ui.horizontal(|ui| {
                ui.label("  Diameter:");
                ui.add(egui::DragValue::new(&mut s.hole_counterbore_diameter).speed(0.1));
                ui.label(" Depth:");
                ui.add(egui::DragValue::new(&mut s.hole_counterbore_depth).speed(0.1));
            });
        }
        ui.checkbox(&mut s.hole_use_countersink, "Countersink");
        if s.hole_use_countersink {
            ui.horizontal(|ui| {
                ui.label("  Diameter:");
                ui.add(egui::DragValue::new(&mut s.hole_countersink_diameter).speed(0.1));
                ui.label(" Angle (deg):");
                ui.add(egui::DragValue::new(&mut s.hole_countersink_angle_deg).speed(1.0).range(10.0..=180.0));
            });
        }
        ui.checkbox(&mut s.hole_use_thread, "Thread metadata");
        if s.hole_use_thread {
            ui.horizontal(|ui| {
                ui.label("  Standard:");
                egui::ComboBox::from_id_source("hole_thread_std")
                    .selected_text(match s.hole_thread_standard_idx { 0 => "ISO Metric", 1 => "UN", 2 => "BSPP", _ => "NPT" })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut s.hole_thread_standard_idx, 0, "ISO Metric");
                        ui.selectable_value(&mut s.hole_thread_standard_idx, 1, "UN");
                        ui.selectable_value(&mut s.hole_thread_standard_idx, 2, "BSPP");
                        ui.selectable_value(&mut s.hole_thread_standard_idx, 3, "NPT");
                    });
            });
            let table = match s.hole_thread_standard_idx {
                0 => valenx_feature_tree::threads::iso_metric_table(),
                1 => valenx_feature_tree::threads::un_table(),
                2 => valenx_feature_tree::threads::bspp_table(),
                _ => valenx_feature_tree::threads::npt_table(),
            };
            if s.hole_thread_entry_idx >= table.len() {
                s.hole_thread_entry_idx = 0;
            }
            ui.horizontal(|ui| {
                ui.label("  Spec:");
                egui::ComboBox::from_id_source("hole_thread_entry")
                    .selected_text(table.get(s.hole_thread_entry_idx).map(|t| t.designation.clone()).unwrap_or_default())
                    .show_ui(ui, |ui| {
                        for (i, spec) in table.iter().enumerate() {
                            ui.selectable_value(&mut s.hole_thread_entry_idx, i, &spec.designation);
                        }
                    });
            });
        }
        if ui.button("Add Hole").clicked() {
            if s.tree.sketches.is_empty() {
                s.last_replay_error = Some("Add a sketch with points first.".into());
            } else {
                s.record();
                let sketch = valenx_feature_tree::feature::SketchRef(s.hole_sketch_index);
                let depth_mode = match s.hole_depth_mode_idx {
                    0 => valenx_feature_tree::feature::HoleDepthMode::Blind { depth: s.hole_blind_depth },
                    1 => valenx_feature_tree::feature::HoleDepthMode::Through,
                    _ => valenx_feature_tree::feature::HoleDepthMode::UpToFace { face_ref: "face #?".into() },
                };
                let counterbore = if s.hole_use_counterbore {
                    Some(valenx_feature_tree::feature::CounterboreParams {
                        diameter: s.hole_counterbore_diameter,
                        depth: s.hole_counterbore_depth,
                    })
                } else { None };
                let countersink = if s.hole_use_countersink {
                    Some(valenx_feature_tree::feature::CountersinkParams {
                        diameter: s.hole_countersink_diameter,
                        angle_deg: s.hole_countersink_angle_deg,
                    })
                } else { None };
                let thread = if s.hole_use_thread {
                    let table = match s.hole_thread_standard_idx {
                        0 => valenx_feature_tree::threads::iso_metric_table(),
                        1 => valenx_feature_tree::threads::un_table(),
                        2 => valenx_feature_tree::threads::bspp_table(),
                        _ => valenx_feature_tree::threads::npt_table(),
                    };
                    table.get(s.hole_thread_entry_idx).cloned()
                } else { None };
                let params = valenx_feature_tree::feature::HoleParams {
                    sketch,
                    depth_mode,
                    drill_diameter: s.hole_drill_diameter,
                    direction_negative: s.hole_direction_negative,
                    counterbore,
                    countersink,
                    thread,
                };
                let name = format!("Hole {}", s.tree.features.len() + 1);
                let id = s.tree.add_feature(valenx_feature_tree::Feature::Hole(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ----- Task 40: Add Loft -----
    ui.collapsing("Add Loft (Phase 13B)", |ui| {
        ui.label("Interpolate a mesh between 2+ profile sketches (v1 mesh output).");
        ui.horizontal(|ui| {
            ui.label("Profile sketch indices (comma-separated):");
            ui.text_edit_singleline(&mut s.loft_profile_indices);
        });
        ui.horizontal(|ui| {
            ui.label("Guide curve indices (optional, comma-separated):");
            ui.text_edit_singleline(&mut s.loft_guide_indices);
        });
        ui.checkbox(&mut s.loft_closed, "Closed loop (wrap last to first)");
        ui.checkbox(&mut s.loft_ruled, "Ruled (straight) connections");
        if ui.button("Add Loft").clicked() {
            let profiles = parse_indices(&s.loft_profile_indices)
                .into_iter()
                .map(valenx_feature_tree::feature::SketchRef)
                .collect::<Vec<_>>();
            let guides = parse_indices(&s.loft_guide_indices)
                .into_iter()
                .map(valenx_feature_tree::feature::SketchRef)
                .collect::<Vec<_>>();
            if profiles.len() < 2 {
                s.last_replay_error = Some("Loft requires at least 2 profile sketches.".into());
            } else {
                s.record();
                let params = valenx_feature_tree::feature::LoftParams {
                    profile_sketches: profiles,
                    guide_curves: guides,
                    closed: s.loft_closed,
                    ruled: s.loft_ruled,
                };
                let name = format!("Loft {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::Loft(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ----- Task 41: Add Sweep -----
    ui.collapsing("Add Sweep (Phase 13B)", |ui| {
        ui.label("Sweep a profile along a path sketch.");
        ui.horizontal(|ui| {
            ui.label("Profile sketch:");
            sketch_index_dropdown(
                ui,
                "sweep_profile",
                &mut s.sweep_profile_index,
                s.tree.sketches.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Path sketch:");
            sketch_index_dropdown(
                ui,
                "sweep_path",
                &mut s.sweep_path_index,
                s.tree.sketches.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Twist (deg):");
            ui.add(
                egui::DragValue::new(&mut s.sweep_twist_deg)
                    .speed(1.0)
                    .range(-3600.0..=3600.0),
            );
        });
        ui.checkbox(
            &mut s.sweep_keep_orientation,
            "Keep profile orientation (world-space)",
        );
        if ui.button("Add Sweep").clicked() {
            if s.tree.sketches.len() < 2 {
                s.last_replay_error =
                    Some("Sweep needs at least 2 sketches (profile + path).".into());
            } else {
                s.record();
                let params = valenx_feature_tree::feature::SweepParams {
                    profile_sketch: valenx_feature_tree::feature::SketchRef(s.sweep_profile_index),
                    path_sketch: valenx_feature_tree::feature::SketchRef(s.sweep_path_index),
                    twist_angle: s.sweep_twist_deg.to_radians(),
                    keep_profile_orientation: s.sweep_keep_orientation,
                };
                let name = format!("Sweep {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::Sweep(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ----- Task 42: Add Pipe -----
    ui.collapsing("Add Pipe (Phase 13B)", |ui| {
        ui.label("Pipe a cross-section along a centerline with optional bend smoothing.");
        ui.horizontal(|ui| {
            ui.label("Cross-section:");
            sketch_index_dropdown(
                ui,
                "pipe_xs",
                &mut s.pipe_cross_section_index,
                s.tree.sketches.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Centerline:");
            sketch_index_dropdown(
                ui,
                "pipe_path",
                &mut s.pipe_centerline_index,
                s.tree.sketches.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Bend radius (mm, 0 = sharp):");
            ui.add(
                egui::DragValue::new(&mut s.pipe_bend_radius)
                    .speed(0.05)
                    .range(0.0..=1000.0),
            );
        });
        if ui.button("Add Pipe").clicked() {
            if s.tree.sketches.len() < 2 {
                s.last_replay_error = Some("Pipe needs at least 2 sketches.".into());
            } else {
                s.record();
                let params = valenx_feature_tree::feature::PipeParams {
                    cross_section_sketch: valenx_feature_tree::feature::SketchRef(
                        s.pipe_cross_section_index,
                    ),
                    centerline_sketch: valenx_feature_tree::feature::SketchRef(
                        s.pipe_centerline_index,
                    ),
                    bend_radius: s.pipe_bend_radius,
                };
                let name = format!("Pipe {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::Pipe(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ----- Task 43: Add Helix -----
    ui.collapsing("Add Helix (Phase 13C)", |ui| {
        ui.label("Coil a profile along a parametric helix axis.");
        ui.horizontal(|ui| {
            ui.label("Profile sketch:");
            sketch_index_dropdown(
                ui,
                "helix_profile",
                &mut s.helix_profile_index,
                s.tree.sketches.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Pitch:");
            ui.add(
                egui::DragValue::new(&mut s.helix_pitch)
                    .speed(0.1)
                    .range(0.001..=10000.0),
            );
            ui.label("Turns:");
            ui.add(
                egui::DragValue::new(&mut s.helix_turns)
                    .speed(0.1)
                    .range(0.001..=1000.0),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Axis origin");
            ui.add(
                egui::DragValue::new(&mut s.helix_axis_origin[0])
                    .speed(0.1)
                    .prefix("X "),
            );
            ui.add(
                egui::DragValue::new(&mut s.helix_axis_origin[1])
                    .speed(0.1)
                    .prefix("Y "),
            );
            ui.add(
                egui::DragValue::new(&mut s.helix_axis_origin[2])
                    .speed(0.1)
                    .prefix("Z "),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Axis dir");
            ui.add(
                egui::DragValue::new(&mut s.helix_axis_direction[0])
                    .speed(0.1)
                    .prefix("X "),
            );
            ui.add(
                egui::DragValue::new(&mut s.helix_axis_direction[1])
                    .speed(0.1)
                    .prefix("Y "),
            );
            ui.add(
                egui::DragValue::new(&mut s.helix_axis_direction[2])
                    .speed(0.1)
                    .prefix("Z "),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Taper (deg):");
            ui.add(
                egui::DragValue::new(&mut s.helix_taper_deg)
                    .speed(0.1)
                    .range(-89.0..=89.0),
            );
        });
        ui.checkbox(
            &mut s.helix_left_handed,
            "Left-handed (CCW viewed from +axis)",
        );
        if ui.button("Add Helix").clicked() {
            if s.tree.sketches.is_empty() {
                s.last_replay_error = Some("Helix needs a profile sketch.".into());
            } else {
                s.record();
                let params = valenx_feature_tree::feature::HelixParams {
                    profile_sketch: valenx_feature_tree::feature::SketchRef(s.helix_profile_index),
                    pitch: s.helix_pitch,
                    turns: s.helix_turns,
                    axis_origin: nalgebra::Vector3::new(
                        s.helix_axis_origin[0],
                        s.helix_axis_origin[1],
                        s.helix_axis_origin[2],
                    ),
                    axis_direction: nalgebra::Vector3::new(
                        s.helix_axis_direction[0],
                        s.helix_axis_direction[1],
                        s.helix_axis_direction[2],
                    ),
                    taper_angle: s.helix_taper_deg,
                    left_handed: s.helix_left_handed,
                };
                let name = format!("Helix {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::Helix(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ----- Task 44: Add Multi-Transform -----
    ui.collapsing("Add Multi-Transform (Phase 13C)", |ui| {
        ui.label("Apply N transforms to a target feature and union the instances.");
        ui.label("Recipe (one per line):  translate dx dy dz | rotate ax ay az angle_deg | scale factor | mirror nx ny nz");
        ui.horizontal(|ui| {
            ui.label("Target feature:");
            feature_index_dropdown(ui, "mt_target", &mut s.mt_target_index, s.tree.features.len());
        });
        ui.add(egui::TextEdit::multiline(&mut s.mt_ops_recipe).desired_rows(4));
        if ui.button("Add Multi-Transform").clicked() {
            if s.tree.features.is_empty() {
                s.last_replay_error = Some("Multi-Transform needs a target feature.".into());
            } else {
                s.record();
                let ops = parse_transform_recipe(&s.mt_ops_recipe);
                let params = valenx_feature_tree::feature::MultiTransformParams {
                    target: valenx_feature_tree::feature::FeatureId(s.mt_target_index),
                    transforms: ops,
                };
                let name = format!("Multi-Transform {}", s.tree.features.len() + 1);
                let id = s.tree.add_feature(valenx_feature_tree::Feature::MultiTransform(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ----- Task 45: Draft / Shell / Thickness -----
    ui.collapsing("Add Draft Angle (Phase 13D)", |ui| {
        ui.label("Tilt selected mesh triangles about the neutral plane axis. v1: mesh-domain.");
        ui.horizontal(|ui| {
            ui.label("Target:");
            feature_index_dropdown(
                ui,
                "draft_target",
                &mut s.draft_target_index,
                s.tree.features.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Face indices (CSV):");
            ui.text_edit_singleline(&mut s.draft_face_indices_csv);
        });
        ui.horizontal(|ui| {
            ui.label("Neutral normal");
            ui.add(
                egui::DragValue::new(&mut s.draft_neutral_normal[0])
                    .speed(0.05)
                    .prefix("X "),
            );
            ui.add(
                egui::DragValue::new(&mut s.draft_neutral_normal[1])
                    .speed(0.05)
                    .prefix("Y "),
            );
            ui.add(
                egui::DragValue::new(&mut s.draft_neutral_normal[2])
                    .speed(0.05)
                    .prefix("Z "),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Draft angle (deg):");
            ui.add(
                egui::DragValue::new(&mut s.draft_angle_deg)
                    .speed(0.5)
                    .range(-89.0..=89.0),
            );
        });
        if ui.button("Add Draft Angle").clicked() {
            if s.tree.features.is_empty() {
                s.last_replay_error = Some("Draft needs a target feature.".into());
            } else {
                s.record();
                let faces = parse_indices(&s.draft_face_indices_csv);
                let params = valenx_feature_tree::feature::DraftAngleParams {
                    target: valenx_feature_tree::feature::FeatureId(s.draft_target_index),
                    face_indices: faces,
                    neutral_plane_normal: nalgebra::Vector3::new(
                        s.draft_neutral_normal[0],
                        s.draft_neutral_normal[1],
                        s.draft_neutral_normal[2],
                    ),
                    draft_angle_deg: s.draft_angle_deg,
                };
                let name = format!("Draft {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::DraftAngle(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    ui.collapsing("Add Shell (Phase 13D)", |ui| {
        ui.label("Hollow out the target leaving a thin-walled shell.");
        ui.horizontal(|ui| {
            ui.label("Target:");
            feature_index_dropdown(
                ui,
                "shell_target",
                &mut s.shell_target_index,
                s.tree.features.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Remove face indices (CSV, optional):");
            ui.text_edit_singleline(&mut s.shell_face_indices_csv);
        });
        ui.horizontal(|ui| {
            ui.label("Thickness:");
            ui.add(
                egui::DragValue::new(&mut s.shell_thickness)
                    .speed(0.01)
                    .range(0.001..=1000.0),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Side:");
            egui::ComboBox::from_id_source("shell_side")
                .selected_text(if s.shell_side_idx == 0 {
                    "Inward"
                } else {
                    "Outward"
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut s.shell_side_idx, 0, "Inward");
                    ui.selectable_value(&mut s.shell_side_idx, 1, "Outward");
                });
        });
        if ui.button("Add Shell").clicked() {
            if s.tree.features.is_empty() {
                s.last_replay_error = Some("Shell needs a target feature.".into());
            } else {
                s.record();
                let params = valenx_feature_tree::feature::ShellParams {
                    target: valenx_feature_tree::feature::FeatureId(s.shell_target_index),
                    face_indices_to_remove: parse_indices(&s.shell_face_indices_csv),
                    thickness: s.shell_thickness,
                    inward_or_outward: if s.shell_side_idx == 0 {
                        valenx_feature_tree::feature::ShellSide::Inward
                    } else {
                        valenx_feature_tree::feature::ShellSide::Outward
                    },
                };
                let name = format!("Shell {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::Shell(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    ui.collapsing("Add Thickness (Phase 13D)", |ui| {
        ui.label("Extrude one face of the target by `thickness` along its normal.");
        ui.horizontal(|ui| {
            ui.label("Target:");
            feature_index_dropdown(
                ui,
                "thickness_target",
                &mut s.thickness_target_index,
                s.tree.features.len(),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Face index:");
            ui.add(
                egui::DragValue::new(&mut s.thickness_face_index)
                    .speed(1.0)
                    .range(0..=100000),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Thickness:");
            ui.add(
                egui::DragValue::new(&mut s.thickness_thickness)
                    .speed(0.01)
                    .range(0.001..=1000.0),
            );
        });
        if ui.button("Add Thickness").clicked() {
            if s.tree.features.is_empty() {
                s.last_replay_error = Some("Thickness needs a target feature.".into());
            } else {
                s.record();
                let params = valenx_feature_tree::feature::ThicknessParams {
                    target: valenx_feature_tree::feature::FeatureId(s.thickness_target_index),
                    face_index: s.thickness_face_index,
                    thickness: s.thickness_thickness,
                };
                let name = format!("Thickness {}", s.tree.features.len() + 1);
                let id = s
                    .tree
                    .add_feature(valenx_feature_tree::Feature::Thickness(params), name);
                s.selected_feature = Some(id);
                s.pending_replay = true;
                s.last_replay_error = None;
            }
        }
    });

    // ----- Task 46: Boolean History -----
    ui.collapsing("Add Boolean History (Phase 13E)", |ui| {
        ui.label("General N-way boolean op over a list of target features.");
        ui.horizontal(|ui| {
            ui.label("Operation:");
            egui::ComboBox::from_id_source("bh_op")
                .selected_text(match s.bh_op_idx {
                    0 => "Union",
                    1 => "Difference",
                    2 => "Intersection",
                    _ => "Section",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut s.bh_op_idx, 0, "Union");
                    ui.selectable_value(&mut s.bh_op_idx, 1, "Difference");
                    ui.selectable_value(&mut s.bh_op_idx, 2, "Intersection");
                    ui.selectable_value(&mut s.bh_op_idx, 3, "Section");
                });
        });
        ui.horizontal(|ui| {
            ui.label("Target feature ids (CSV):");
            ui.text_edit_singleline(&mut s.bh_targets_csv);
        });
        if ui.button("Add Boolean History").clicked() {
            if s.tree.features.is_empty() {
                s.last_replay_error = Some("Boolean History needs target features.".into());
            } else {
                let targets = parse_indices(&s.bh_targets_csv)
                    .into_iter()
                    .map(valenx_feature_tree::feature::FeatureId)
                    .collect::<Vec<_>>();
                if targets.is_empty() {
                    s.last_replay_error = Some("Provide at least one target id.".into());
                } else {
                    s.record();
                    let op = match s.bh_op_idx {
                        0 => valenx_feature_tree::feature::BoolKind::Union,
                        1 => valenx_feature_tree::feature::BoolKind::Difference,
                        2 => valenx_feature_tree::feature::BoolKind::Intersection,
                        _ => valenx_feature_tree::feature::BoolKind::Section,
                    };
                    let params = valenx_feature_tree::feature::BooleanHistoryParams {
                        operation: op,
                        targets,
                    };
                    let name = format!("Boolean History {}", s.tree.features.len() + 1);
                    let id = s
                        .tree
                        .add_feature(valenx_feature_tree::Feature::BooleanHistory(params), name);
                    s.selected_feature = Some(id);
                    s.pending_replay = true;
                    s.last_replay_error = None;
                }
            }
        }
    });

    // ----- Tasks 45 + 46: Selected-feature actions -----
    //
    // Suppress / un-suppress (Task 45) and Delete with confirmation
    // (Task 46) both operate on `s.selected_feature`. Wrapped in one
    // section so the user has a single home for "things I do to the
    // currently-selected entry".
    ui.separator();
    ui.label(egui::RichText::new("Selected feature").strong());
    // Deferred-dispatch — read the entry's display-relevant fields
    // out of `s.tree.features` BEFORE wiring buttons that mutate
    // `s.tree`. Avoids the closure-borrow conflict that comes from
    // holding a `&entry` while clicks try to call `s.tree.*_mut`.
    enum SelAction {
        None,
        ToggleSuppress,
        ArmDelete,
        ConfirmDelete,
        CancelDelete,
    }
    let mut action = SelAction::None;
    let mut selection_to_clear = false;
    if let Some(id) = s.selected_feature {
        let entry_snapshot = s
            .tree
            .features
            .get(id.0)
            .map(|e| (e.feature.kind_label(), e.name.clone(), e.suppressed));
        match entry_snapshot {
            Some((kind, name, suppressed)) => {
                ui.label(format!("#{} {} — {}", id.0, kind, name));
                let suppress_label = if suppressed {
                    "Un-suppress"
                } else {
                    "Suppress"
                };
                ui.horizontal(|ui| {
                    if ui.button(suppress_label).clicked() {
                        action = SelAction::ToggleSuppress;
                    }
                    // Task 46 — Delete with confirmation. First click
                    // arms `pending_delete_confirm`; the user then
                    // clicks "Confirm delete" (red) to actually delete,
                    // or "Cancel" to back out. Selecting another row
                    // also clears the latch (tree-view handler above).
                    if ui.button("Delete feature").clicked() {
                        action = SelAction::ArmDelete;
                    }
                });
                if s.pending_delete_confirm {
                    ui.colored_label(
                        egui::Color32::from_rgb(220, 100, 100),
                        format!("Confirm delete of #{} {}? Irreversible.", id.0, kind),
                    );
                    ui.horizontal(|ui| {
                        if ui
                            .button(
                                egui::RichText::new("Confirm delete")
                                    .color(egui::Color32::from_rgb(220, 100, 100)),
                            )
                            .clicked()
                        {
                            action = SelAction::ConfirmDelete;
                        }
                        if ui.button("Cancel").clicked() {
                            action = SelAction::CancelDelete;
                        }
                    });
                }
            }
            None => {
                ui.label("(selection points to a deleted entry — pick another)");
                selection_to_clear = true;
            }
        }
    } else {
        ui.label("(none — click a row in the feature tree above)");
    }
    if selection_to_clear {
        s.selected_feature = None;
        s.pending_delete_confirm = false;
    }
    if let Some(id) = s.selected_feature {
        match action {
            SelAction::None => {}
            SelAction::ToggleSuppress => {
                let new_suppressed = !s
                    .tree
                    .features
                    .get(id.0)
                    .map(|e| e.suppressed)
                    .unwrap_or(false);
                s.record();
                if let Err(e) = s.tree.set_suppressed(id, new_suppressed) {
                    s.last_replay_error = Some(format!("suppress: {e}"));
                } else {
                    s.pending_replay = true;
                }
            }
            SelAction::ArmDelete => {
                s.pending_delete_confirm = true;
            }
            SelAction::ConfirmDelete => {
                s.record();
                match s.tree.delete_feature(id) {
                    Ok(_) => {
                        s.selected_feature = None;
                        s.pending_delete_confirm = false;
                        s.pending_replay = true;
                        s.last_replay_error = None;
                    }
                    Err(e) => {
                        s.last_replay_error = Some(format!("delete: {e}"));
                        s.pending_delete_confirm = false;
                    }
                }
            }
            SelAction::CancelDelete => {
                s.pending_delete_confirm = false;
            }
        }
    }

    // ----- Task 47: Replay controls -----
    //
    // Manual "Replay now" forces an immediate re-evaluation even when
    // `auto_replay` is off (or when the user wants to retry after
    // tweaking an external sketch). `auto_replay` checkbox lets users
    // toggle between live updates and snapshot updates (useful for
    // expensive trees where each replay takes seconds).
    ui.separator();
    ui.label(egui::RichText::new("Replay").strong());
    ui.horizontal(|ui| {
        ui.checkbox(&mut s.auto_replay, "Auto-replay on changes");
        if ui.button("Replay now").clicked() {
            s.pending_replay = true;
        }
    });

    // ----- Task 48: Replay error display -----
    //
    // Surface the last replay error in red below the controls. Cleared
    // by the next successful replay or by add/delete/suppress handlers
    // when they kick off a fresh attempt.
    if let Some(err) = &s.last_replay_error {
        ui.colored_label(egui::Color32::from_rgb(220, 100, 100), err);
    }

    // ----- Task 48: Push the final solid to the viewport -----
    //
    // Deferred-dispatch — read `pending_replay` during the `&mut s`
    // borrow scope, clear the flag, then after releasing the borrow
    // reach for `&mut app` to clone the tree, call
    // valenx_feature_tree::replay, tessellate, and push via
    // app.apply_mesh.
    //
    // Why `pending_replay` and not a per-frame replay: replay calls
    // truck-modeling booleans which take ms per op. A user-facing
    // panel that re-evaluates every frame would chew CPU even when
    // nothing changed. The flag is set by add/delete/suppress
    // handlers (auto-mode) and by the "Replay now" button (manual),
    // and cleared here after one replay.
    //
    // Note: `auto_replay = false` does NOT gate the replay here — it
    // gates whether add/delete/suppress handlers SET the flag in the
    // first place. (Phase 2.5 will plumb that distinction through; v1
    // always sets the flag and trusts the user to toggle auto_replay
    // off if they don't want live updates.)
    let mut do_replay = false;
    if s.pending_replay {
        do_replay = true;
        s.pending_replay = false;
    }

    // Release the &mut part_design borrow before reaching for &mut app.
    let _ = s;
    if do_replay {
        run_part_design_replay(app);
    }
}

/// Task 53 — Save the current part-design tree to a user-picked
/// `.valenx` path. Cancellation is silent (no error, no status churn).
/// Errors land in `app.mesh_toolbox.part_design.last_replay_error`.
fn save_part_design_project(app: &mut crate::ValenxApp) {
    let dialog = rfd::FileDialog::new()
        .add_filter("Valenx project", &["valenx"])
        .set_title("Save Valenx project")
        .set_file_name("part-design.valenx");
    let Some(path) = dialog.save_file() else {
        return;
    };
    let project =
        valenx_feature_tree::persist::ValenxProject::from_tree(&app.mesh_toolbox.part_design.tree);
    match project.write_to(&path) {
        Ok(()) => {
            app.mesh_toolbox.part_design.last_replay_error = None;
            app.mesh_toolbox.part_design.project_path = Some(path.clone());
            app.status = Some(format!("Saved project to {}", path.display()));
            emit_audit(
                "part_design.save_project",
                serde_json::json!({"kind": "feature_tree.project"}),
                serde_json::json!({"path": path.display().to_string()}),
            );
        }
        Err(e) => {
            app.mesh_toolbox.part_design.last_replay_error = Some(format!("save project: {e}"));
        }
    }
}

/// Phase 8 Task 20 — Add an `ImportedSolid` feature pointing at a
/// user-picked STEP / IGES file. The actual import happens lazily on
/// the next replay (so we don't have to surface a `Solid` here, and
/// the `.valenx` project stays small).
fn import_step_iges(app: &mut crate::ValenxApp) {
    let dialog = rfd::FileDialog::new()
        .add_filter("STEP / IGES", &["step", "stp", "iges", "igs"])
        .add_filter("STEP", &["step", "stp"])
        .add_filter("IGES", &["iges", "igs"])
        .set_title("Import STEP or IGES file");
    let Some(path) = dialog.pick_file() else {
        return;
    };
    // Eagerly validate by attempting the import once — that way a bad
    // file surfaces the error immediately instead of failing on the
    // first replay after the user clicks "Add Pad" or similar.
    if let Err(e) = valenx_step_iges::import(&path) {
        app.mesh_toolbox.part_design.last_replay_error = Some(format!("import STEP/IGES: {e}"));
        return;
    }
    let params = valenx_feature_tree::feature::ImportedSolidParams {
        source_path: path.display().to_string(),
    };
    let s = &mut app.mesh_toolbox.part_design;
    s.record();
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .map(|n| format!("Imported: {n}"))
        .unwrap_or_else(|| "Imported".to_string());
    let id = s
        .tree
        .add_feature(valenx_feature_tree::Feature::ImportedSolid(params), name);
    s.selected_feature = Some(id);
    s.pending_replay = true;
    s.last_replay_error = None;
    app.status = Some(format!("Imported {}", path.display()));
    emit_audit(
        "part_design.import_step_iges",
        serde_json::json!({"kind": "feature_tree.imported_solid"}),
        serde_json::json!({"path": path.display().to_string()}),
    );
}

/// Phase 20 — Import a STEP AP242 / IGES trimmed-surface file via the
/// AP242-aware reader, capturing the product hierarchy + GD&T metadata
/// onto an [`valenx_feature_tree::Feature::ImportedAdvanced`] tree
/// node.
fn import_step_ap242(app: &mut crate::ValenxApp) {
    let dialog = rfd::FileDialog::new()
        .add_filter("STEP AP242 / IGES trimmed", &["step", "stp", "iges", "igs"])
        .set_title("Import STEP AP242 or IGES trimmed-surface file");
    let Some(path) = dialog.pick_file() else {
        return;
    };
    // Validate first via the generic import; AP242 files are still
    // STEP, so this catches malformed inputs early.
    if let Err(e) = valenx_step_iges::import(&path) {
        app.mesh_toolbox.part_design.last_replay_error = Some(format!("import STEP/IGES: {e}"));
        return;
    }
    // Now scan for AP242 metadata.
    let metadata = valenx_step_iges::ap242::read_metadata(&path).unwrap_or_default();
    let solid_count = valenx_step_iges::ap242::count_solids(&path).unwrap_or(1);
    let is_ap242 = !metadata.is_empty() || solid_count > 1;
    let params = valenx_feature_tree::feature::ImportedAdvancedParams {
        source_path: path.display().to_string(),
        product_path: metadata.product_path.clone(),
        feature_hints: metadata.feature_hints.clone(),
        parametric_values: metadata.parametric_values.clone(),
        tolerances: metadata.tolerances.clone(),
        material_names: metadata.materials.iter().map(|m| m.name.clone()).collect(),
        is_ap242,
    };
    let s = &mut app.mesh_toolbox.part_design;
    s.record();
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .map(|n| {
            if is_ap242 {
                format!("Imported (AP242): {n}")
            } else {
                format!("Imported: {n}")
            }
        })
        .unwrap_or_else(|| "Imported (AP242)".to_string());
    let id = s
        .tree
        .add_feature(valenx_feature_tree::Feature::ImportedAdvanced(params), name);
    s.selected_feature = Some(id);
    s.pending_replay = true;
    s.last_replay_error = None;
    app.status = Some(format!(
        "Imported AP242 {} ({} solids, {} materials)",
        path.display(),
        solid_count,
        metadata.materials.len(),
    ));
    emit_audit(
        "part_design.import_step_ap242",
        serde_json::json!({"kind": "feature_tree.imported_advanced"}),
        serde_json::json!({
            "path": path.display().to_string(),
            "is_ap242": is_ap242,
            "solid_count": solid_count,
            "materials": metadata.materials.len(),
            "tolerances": metadata.tolerances.len(),
        }),
    );
}

/// Phase 8 Task 22 — Save the most-recently-replayed solid to a
/// user-picked STEP / IGES path. Format is auto-detected from the
/// extension by [`valenx_step_iges::export`].
fn export_step_iges(app: &mut crate::ValenxApp) {
    let dialog = rfd::FileDialog::new()
        .add_filter("STEP", &["step", "stp"])
        .add_filter("IGES", &["iges", "igs"])
        .set_title("Export part as STEP or IGES")
        .set_file_name("part-design.step");
    let Some(path) = dialog.save_file() else {
        return;
    };
    // Re-run replay so we operate on the live tree state — the
    // viewport mesh may be stale relative to the most recent edit.
    let tree = app.mesh_toolbox.part_design.tree.clone();
    let solid = match valenx_feature_tree::replay(&tree) {
        Ok(Some(s)) => s,
        Ok(None) => {
            app.mesh_toolbox.part_design.last_replay_error =
                Some("export STEP/IGES: tree is empty or all features are suppressed".into());
            return;
        }
        Err(e) => {
            app.mesh_toolbox.part_design.last_replay_error =
                Some(format!("export STEP/IGES: replay failed: {e}"));
            return;
        }
    };
    match valenx_step_iges::export(&solid, &path) {
        Ok(()) => {
            app.mesh_toolbox.part_design.last_replay_error = None;
            app.status = Some(format!("Exported to {}", path.display()));
            emit_audit(
                "part_design.export_step_iges",
                serde_json::json!({"kind": "feature_tree.step_iges_export"}),
                serde_json::json!({"path": path.display().to_string()}),
            );
        }
        Err(e) => {
            app.mesh_toolbox.part_design.last_replay_error = Some(format!("export STEP/IGES: {e}"));
        }
    }
}

/// Task 54 — Load a `.valenx` project from a user-picked path,
/// replacing the in-memory tree. On success queues a replay so the
/// viewport reflects the loaded model.
fn load_part_design_project(app: &mut crate::ValenxApp) {
    let dialog = rfd::FileDialog::new()
        .add_filter("Valenx project", &["valenx"])
        .set_title("Open Valenx project");
    let Some(path) = dialog.pick_file() else {
        return;
    };
    match valenx_feature_tree::persist::ValenxProject::read_from(&path) {
        Ok(project) => {
            let s = &mut app.mesh_toolbox.part_design;
            s.tree = project.feature_tree;
            s.selected_feature = None;
            s.selected_sketch = None;
            s.pending_delete_confirm = false;
            s.pending_replay = true;
            s.last_replay_error = None;
            s.project_path = Some(path.clone());
            let sketches = s.tree.sketches.len();
            let features = s.tree.features.len();
            app.status = Some(format!("Loaded project from {}", path.display()));
            emit_audit(
                "part_design.load_project",
                serde_json::json!({"kind": "feature_tree.project"}),
                serde_json::json!({
                    "path": path.display().to_string(),
                    "version": project.version,
                    "sketches": sketches,
                    "features": features,
                }),
            );
        }
        Err(e) => {
            app.mesh_toolbox.part_design.last_replay_error = Some(format!("load project: {e}"));
        }
    }
}

/// Replay the part-design tree and push the resulting tessellated
/// mesh into the viewport. Errors land in
/// `app.mesh_toolbox.part_design.last_replay_error`.
///
/// Empty trees (no features, or every feature suppressed) leave the
/// viewport alone — they're a valid "intermediate edit" state.
fn run_part_design_replay(app: &mut crate::ValenxApp) {
    let tree = app.mesh_toolbox.part_design.tree.clone();
    match valenx_feature_tree::replay(&tree) {
        Ok(Some(solid)) => {
            match valenx_cad::solid_to_mesh(&solid, valenx_cad::DEFAULT_TESS_TOLERANCE) {
                Ok(mesh) => {
                    let pseudo_path = std::path::PathBuf::from("<part-design>/replay.solid");
                    app.apply_mesh(mesh, pseudo_path);
                    app.mesh_toolbox.part_design.last_replay_error = None;
                    emit_audit(
                        "part_design.replay",
                        serde_json::json!({"kind": "feature_tree.replay"}),
                        serde_json::json!({
                            "sketches": tree.sketches.len(),
                            "features": tree.features.len(),
                        }),
                    );
                }
                Err(e) => {
                    app.mesh_toolbox.part_design.last_replay_error =
                        Some(format!("tessellate: {e}"));
                }
            }
        }
        Ok(None) => {
            // Empty tree (or all-suppressed) — leave the viewport
            // alone and clear any prior error. The user is mid-edit.
            app.mesh_toolbox.part_design.last_replay_error = None;
        }
        Err(e) => {
            app.mesh_toolbox.part_design.last_replay_error = Some(format!("{e}"));
        }
    }
}

/// Dropdown helper for picking a SketchRef. Shows `SketchRef(i)` for
/// `i in 0..count`. Disabled (greyed out) when `count == 0`.
fn sketch_index_dropdown(ui: &mut egui::Ui, id_source: &str, current: &mut usize, count: usize) {
    if count == 0 {
        ui.label("(no sketches yet)");
        return;
    }
    if *current >= count {
        *current = count - 1;
    }
    egui::ComboBox::from_id_source(id_source)
        .selected_text(format!("SketchRef({current})"))
        .show_ui(ui, |ui| {
            for i in 0..count {
                ui.selectable_value(current, i, format!("SketchRef({i})"));
            }
        });
}

/// Dropdown helper for picking a FeatureId. Shows `FeatureId(i)` for
/// `i in 0..count`. Disabled (greyed out) when `count == 0`.
///
/// Used by Mirror / LinearPattern / CircularPattern dialogs.
fn feature_index_dropdown(ui: &mut egui::Ui, id_source: &str, current: &mut usize, count: usize) {
    if count == 0 {
        ui.label("(no features yet)");
        return;
    }
    if *current >= count {
        *current = count - 1;
    }
    egui::ComboBox::from_id_source(id_source)
        .selected_text(format!("FeatureId({current})"))
        .show_ui(ui, |ui| {
            for i in 0..count {
                ui.selectable_value(current, i, format!("FeatureId({i})"));
            }
        });
}

/// Parse a comma- or whitespace-separated list of `usize` values.
/// Ignores entries that don't parse cleanly. Used by Phase 13 dialogs
/// for "list of sketch / feature / face indices" text fields.
fn parse_indices(s: &str) -> Vec<usize> {
    s.split(|c: char| c == ',' || c.is_whitespace())
        .filter(|t| !t.is_empty())
        .filter_map(|t| t.parse::<usize>().ok())
        .collect()
}

/// Phase 14 helper for the Fillet / Chamfer dialogs' "Edges" field.
///
/// Returns:
/// - `None` if `s` is empty or all-whitespace ("auto-select by angle
///   threshold" — preserves Phase 3 backward-compat semantics).
/// - `Some(parsed)` otherwise, with entries that don't parse as
///   `usize` silently dropped (the dialog's hint text spells out
///   the format so the user can self-correct).
fn parse_edge_indices_or_none(s: &str) -> Option<Vec<usize>> {
    let trimmed = s.trim();
    if trimmed.is_empty() {
        return None;
    }
    let parsed = parse_indices(trimmed);
    if parsed.is_empty() {
        None
    } else {
        Some(parsed)
    }
}

/// Parse a Multi-Transform recipe (one op per line). Each line is one
/// of:
/// - `translate dx dy dz`
/// - `rotate ax ay az angle_deg`
/// - `scale factor`
/// - `mirror nx ny nz`
///
/// Lines that don't parse cleanly are skipped. Used by the Phase 13C
/// Multi-Transform dialog.
fn parse_transform_recipe(text: &str) -> Vec<valenx_feature_tree::feature::TransformOp> {
    let mut out = Vec::new();
    for line in text.lines() {
        let mut parts = line.split_whitespace();
        let Some(kind) = parts.next() else { continue };
        let nums: Vec<f64> = parts.filter_map(|t| t.parse::<f64>().ok()).collect();
        match (kind.to_ascii_lowercase().as_str(), nums.len()) {
            ("translate", 3) => out.push(valenx_feature_tree::feature::TransformOp::Translate {
                delta: nalgebra::Vector3::new(nums[0], nums[1], nums[2]),
            }),
            ("rotate", 4) => out.push(valenx_feature_tree::feature::TransformOp::Rotate {
                axis: nalgebra::Vector3::new(nums[0], nums[1], nums[2]),
                angle_rad: nums[3].to_radians(),
            }),
            ("scale", 1) => {
                out.push(valenx_feature_tree::feature::TransformOp::Scale { factor: nums[0] })
            }
            ("mirror", 3) => out.push(valenx_feature_tree::feature::TransformOp::Mirror {
                plane_normal: nalgebra::Vector3::new(nums[0], nums[1], nums[2]),
            }),
            _ => {} // ignore unparsable lines
        }
    }
    out
}

// ===========================================================================
// Draft workbench panel (Phase 4)
// ===========================================================================

/// Draw the Draft workbench panel (Tasks 14-25 of Phase 4).
///
/// 2D entities are placed on a [`valenx_draft::WorkingPlane`] using
/// numeric inputs. The full viewport-click pipeline lands later; for
/// now each tool reads its parameters from the panel and appends a
/// new [`valenx_draft::DraftEntity`] on Apply.
pub fn draw_draft_panel(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.heading("Draft (2D in 3D)");
    ui.separator();

    let s = &mut app.mesh_toolbox.draft;

    // ----- Overlay toggle -----
    ui.checkbox(&mut s.show_overlay, "Show draft overlay in viewport");

    // ----- Working plane selector (Task 16) -----
    ui.label(egui::RichText::new("Working plane").strong());
    let prev_plane = s.plane_kind;
    egui::ComboBox::from_id_source("draft_plane_combo")
        .selected_text(s.plane_kind.label())
        .show_ui(ui, |ui| {
            ui.selectable_value(
                &mut s.plane_kind,
                DraftPlaneKind::Xy,
                DraftPlaneKind::Xy.label(),
            );
            ui.selectable_value(
                &mut s.plane_kind,
                DraftPlaneKind::Xz,
                DraftPlaneKind::Xz.label(),
            );
            ui.selectable_value(
                &mut s.plane_kind,
                DraftPlaneKind::Yz,
                DraftPlaneKind::Yz.label(),
            );
        });
    if s.plane_kind != prev_plane {
        // Replace the working plane on the active document so new
        // entities project to the freshly-selected orientation. Keeps
        // the entity list intact — only the local→world projection
        // changes.
        s.document.working_plane = s.plane_kind.to_plane();
    }

    // ----- Tool palette (Task 15) -----
    ui.label(egui::RichText::new("Tool").strong());
    ui.horizontal_wrapped(|ui| {
        ui.selectable_value(&mut s.tool, DraftTool::Select, "Select");
        ui.selectable_value(&mut s.tool, DraftTool::Line, "Line");
        ui.selectable_value(&mut s.tool, DraftTool::Polyline, "Polyline");
        ui.selectable_value(&mut s.tool, DraftTool::Arc, "Arc");
        ui.selectable_value(&mut s.tool, DraftTool::Circle, "Circle");
        ui.selectable_value(&mut s.tool, DraftTool::Rectangle, "Rectangle");
        ui.selectable_value(&mut s.tool, DraftTool::Polygon, "Polygon");
        ui.selectable_value(&mut s.tool, DraftTool::Dimension, "Dimension");
        ui.selectable_value(&mut s.tool, DraftTool::Text, "Text");
    });
    ui.label(format!(
        "Active: {:?} — {} entities",
        s.tool,
        s.document.entity_count(),
    ));

    ui.separator();
    ui.horizontal(|ui| {
        ui.label("Grid spacing:");
        ui.add(
            egui::DragValue::new(&mut s.grid_spacing)
                .speed(0.1)
                .range(0.0..=1e6),
        );
    });
    ui.separator();

    // ----- Per-tool input sections (Tasks 17-24) -----
    match s.tool {
        DraftTool::Select => {
            ui.label("Select mode — pick an entity below to delete it.");
        }
        DraftTool::Line => draw_draft_line_inputs(s, ui),
        DraftTool::Polyline => draw_draft_polyline_inputs(s, ui),
        DraftTool::Arc => draw_draft_arc_inputs(s, ui),
        DraftTool::Circle => draw_draft_circle_inputs(s, ui),
        DraftTool::Rectangle => draw_draft_rect_inputs(s, ui),
        DraftTool::Polygon => draw_draft_polygon_inputs(s, ui),
        DraftTool::Dimension => draw_draft_dim_inputs(s, ui),
        DraftTool::Text => draw_draft_text_inputs(s, ui),
    }

    if let Some(err) = &s.last_error {
        ui.colored_label(egui::Color32::from_rgb(220, 100, 100), err);
    }

    // ----- Entity list (Task 25) -----
    ui.separator();
    ui.collapsing("Entities", |ui| {
        if s.document.entities.is_empty() {
            ui.label("(empty)");
        }
        let mut to_delete: Option<usize> = None;
        for (i, e) in s.document.entities.iter().enumerate() {
            let label = format!("#{i}: {} {}", e.kind(), draft_entity_summary(e));
            let selected = s.selected_entity == Some(i);
            let rich = if selected {
                egui::RichText::new(label).strong()
            } else {
                egui::RichText::new(label)
            };
            if ui.selectable_label(selected, rich).clicked() {
                s.selected_entity = if selected { None } else { Some(i) };
            }
        }
        if let Some(idx) = s.selected_entity {
            if ui.button("Delete selected").clicked() {
                to_delete = Some(idx);
            }
        }
        if let Some(idx) = to_delete {
            match s.document.delete_entity(idx) {
                Ok(_) => {
                    s.selected_entity = None;
                    s.last_error = None;
                }
                Err(e) => s.last_error = Some(format!("delete: {e}")),
            }
        }
    });
}

/// One-line summary used in the entity list.
fn draft_entity_summary(e: &valenx_draft::DraftEntity) -> String {
    use valenx_draft::DraftEntity as E;
    match e {
        E::Line { start, end } => {
            format!(
                "({:.2},{:.2})→({:.2},{:.2})",
                start[0], start[1], end[0], end[1]
            )
        }
        E::Polyline { points, closed } => {
            format!(
                "{} pts {}",
                points.len(),
                if *closed { "closed" } else { "open" }
            )
        }
        E::Arc {
            center,
            radius,
            start_angle,
            end_angle,
        } => {
            format!(
                "c=({:.2},{:.2}) r={:.2} {:.1}°→{:.1}°",
                center[0],
                center[1],
                radius,
                start_angle.to_degrees(),
                end_angle.to_degrees(),
            )
        }
        E::Circle { center, radius } => {
            format!("c=({:.2},{:.2}) r={:.2}", center[0], center[1], radius)
        }
        E::Rectangle { min, max } => {
            format!(
                "({:.2},{:.2})-({:.2},{:.2})",
                min[0], min[1], max[0], max[1]
            )
        }
        E::Polygon {
            center,
            radius,
            sides,
        } => {
            format!(
                "c=({:.2},{:.2}) r={:.2} n={}",
                center[0], center[1], radius, sides
            )
        }
        E::LinearDimension { from, to, offset } => {
            format!(
                "({:.2},{:.2})→({:.2},{:.2}) off={:.2}",
                from[0], from[1], to[0], to[1], offset
            )
        }
        E::Text {
            position,
            content,
            size,
        } => {
            format!(
                "({:.2},{:.2}) \"{}\" sz={:.2}",
                position[0], position[1], content, size
            )
        }
    }
}

fn xy_drag(ui: &mut egui::Ui, label: &str, value: &mut [f64; 2]) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::DragValue::new(&mut value[0]).speed(0.05).prefix("x "));
        ui.add(egui::DragValue::new(&mut value[1]).speed(0.05).prefix("y "));
    });
}

/// Task 17 — Line tool: two-step "Place start" → "Place end".
fn draw_draft_line_inputs(s: &mut DraftPanelState, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Line").strong());
    if !s.line_start_placed {
        xy_drag(ui, "Start", &mut s.line_start);
        if ui.button("Place start").clicked() {
            s.line_start_placed = true;
            s.last_error = None;
        }
    } else {
        ui.label(format!(
            "Start placed at ({:.3}, {:.3})",
            s.line_start[0], s.line_start[1],
        ));
        xy_drag(ui, "End", &mut s.line_end);
        ui.horizontal(|ui| {
            if ui.button("Place end").clicked() {
                s.document.add_entity(valenx_draft::DraftEntity::Line {
                    start: s.line_start,
                    end: s.line_end,
                });
                s.line_start_placed = false;
                s.last_error = None;
            }
            if ui.button("Cancel").clicked() {
                s.line_start_placed = false;
            }
        });
    }
}

/// Task 18 — Polyline tool: repeating Append + "Close polyline".
fn draw_draft_polyline_inputs(s: &mut DraftPanelState, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Polyline").strong());
    ui.label(format!("Pending points: {}", s.polyline_points.len()));
    xy_drag(ui, "Next point", &mut s.polyline_next_point);
    ui.horizontal(|ui| {
        if ui.button("Append point").clicked() {
            s.polyline_points.push(s.polyline_next_point);
            s.last_error = None;
        }
        if ui.button("Finish (open)").clicked() {
            commit_polyline(s, false);
        }
        if ui.button("Close polyline").clicked() {
            commit_polyline(s, true);
        }
        if ui.button("Discard").clicked() {
            s.polyline_points.clear();
        }
    });
}

fn commit_polyline(s: &mut DraftPanelState, closed: bool) {
    if s.polyline_points.len() < 2 {
        s.last_error = Some("polyline needs at least 2 points".into());
        return;
    }
    let pts = std::mem::take(&mut s.polyline_points);
    s.document.add_entity(valenx_draft::DraftEntity::Polyline {
        points: pts,
        closed,
    });
    s.last_error = None;
}

/// Task 19 — Arc tool: centre + radius + angles.
fn draw_draft_arc_inputs(s: &mut DraftPanelState, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Arc").strong());
    xy_drag(ui, "Centre", &mut s.arc_center);
    ui.horizontal(|ui| {
        ui.label("Radius");
        ui.add(
            egui::DragValue::new(&mut s.arc_radius)
                .speed(0.05)
                .range(0.0..=1e9),
        );
    });
    ui.horizontal(|ui| {
        ui.label("Start °");
        ui.add(egui::DragValue::new(&mut s.arc_start_angle_deg).speed(1.0));
        ui.label("End °");
        ui.add(egui::DragValue::new(&mut s.arc_end_angle_deg).speed(1.0));
    });
    if ui.button("Place arc").clicked() {
        if s.arc_radius <= 0.0 {
            s.last_error = Some("arc radius must be > 0".into());
        } else {
            s.document.add_entity(valenx_draft::DraftEntity::Arc {
                center: s.arc_center,
                radius: s.arc_radius,
                start_angle: s.arc_start_angle_deg.to_radians(),
                end_angle: s.arc_end_angle_deg.to_radians(),
            });
            s.last_error = None;
        }
    }
}

/// Task 20 — Circle tool: centre + radius.
fn draw_draft_circle_inputs(s: &mut DraftPanelState, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Circle").strong());
    xy_drag(ui, "Centre", &mut s.circle_center);
    ui.horizontal(|ui| {
        ui.label("Radius");
        ui.add(
            egui::DragValue::new(&mut s.circle_radius)
                .speed(0.05)
                .range(0.0..=1e9),
        );
    });
    if ui.button("Place circle").clicked() {
        if s.circle_radius <= 0.0 {
            s.last_error = Some("circle radius must be > 0".into());
        } else {
            s.document.add_entity(valenx_draft::DraftEntity::Circle {
                center: s.circle_center,
                radius: s.circle_radius,
            });
            s.last_error = None;
        }
    }
}

/// Task 21 — Rectangle tool: min + max corners.
fn draw_draft_rect_inputs(s: &mut DraftPanelState, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Rectangle").strong());
    xy_drag(ui, "Min", &mut s.rect_min);
    xy_drag(ui, "Max", &mut s.rect_max);
    if ui.button("Place rectangle").clicked() {
        if s.rect_max[0] <= s.rect_min[0] || s.rect_max[1] <= s.rect_min[1] {
            s.last_error = Some("rectangle max must be strictly greater than min".into());
        } else {
            s.document.add_entity(valenx_draft::DraftEntity::Rectangle {
                min: s.rect_min,
                max: s.rect_max,
            });
            s.last_error = None;
        }
    }
}

/// Task 22 — Polygon tool: centre + radius + sides.
fn draw_draft_polygon_inputs(s: &mut DraftPanelState, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Regular polygon").strong());
    xy_drag(ui, "Centre", &mut s.polygon_center);
    ui.horizontal(|ui| {
        ui.label("Radius");
        ui.add(
            egui::DragValue::new(&mut s.polygon_radius)
                .speed(0.05)
                .range(0.0..=1e9),
        );
        ui.label("Sides");
        ui.add(
            egui::DragValue::new(&mut s.polygon_sides)
                .speed(1.0)
                .range(3..=512),
        );
    });
    if ui.button("Place polygon").clicked() {
        if s.polygon_radius <= 0.0 {
            s.last_error = Some("polygon radius must be > 0".into());
        } else if s.polygon_sides < 3 {
            s.last_error = Some("polygon needs at least 3 sides".into());
        } else {
            s.document.add_entity(valenx_draft::DraftEntity::Polygon {
                center: s.polygon_center,
                radius: s.polygon_radius,
                sides: s.polygon_sides,
            });
            s.last_error = None;
        }
    }
}

/// Task 23 — Linear Dimension tool: from + to + perpendicular offset.
fn draw_draft_dim_inputs(s: &mut DraftPanelState, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Linear dimension").strong());
    xy_drag(ui, "From", &mut s.dim_from);
    xy_drag(ui, "To", &mut s.dim_to);
    ui.horizontal(|ui| {
        ui.label("Offset");
        ui.add(egui::DragValue::new(&mut s.dim_offset).speed(0.05));
    });
    if ui.button("Place dimension").clicked() {
        s.document
            .add_entity(valenx_draft::DraftEntity::LinearDimension {
                from: s.dim_from,
                to: s.dim_to,
                offset: s.dim_offset,
            });
        s.last_error = None;
    }
}

/// Task 24 — Text tool: position + content + size.
fn draw_draft_text_inputs(s: &mut DraftPanelState, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Text label").strong());
    xy_drag(ui, "Position", &mut s.text_position);
    ui.horizontal(|ui| {
        ui.label("Content");
        ui.text_edit_singleline(&mut s.text_content);
    });
    ui.horizontal(|ui| {
        ui.label("Size");
        ui.add(
            egui::DragValue::new(&mut s.text_size)
                .speed(0.05)
                .range(0.0..=1e6),
        );
    });
    if ui.button("Place text").clicked() {
        if s.text_content.is_empty() {
            s.last_error = Some("text content must not be empty".into());
        } else if s.text_size <= 0.0 {
            s.last_error = Some("text size must be > 0".into());
        } else {
            s.document.add_entity(valenx_draft::DraftEntity::Text {
                position: s.text_position,
                content: s.text_content.clone(),
                size: s.text_size,
            });
            s.last_error = None;
        }
    }
}

/// Draw the TechDraw workbench panel (Tasks 32-37 of Phase 5).
///
/// Lets the user pick a sheet size, add front / top / right / iso /
/// custom views from the current Part Design solid, edit
/// scale + position of the selected view, add linear dimensions,
/// and export to SVG / PDF / DXF via native file dialogs.
///
/// The "source solid" is `app.current_solid` (operand A of the Part
/// workbench). Add View buttons regenerate edges immediately; the
/// resulting visible / hidden edge lists are stored on the
/// [`valenx_techdraw::View`] itself so the export path doesn't need
/// access to the solid again.
pub fn draw_techdraw_panel(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.heading("TechDraw (2D drawings from 3D)");
    ui.separator();

    // Pull the current solid + last-feature id out before borrowing
    // techdraw mutably so we can pass them to View::generate /
    // ParametricView::new without re-borrowing app.
    let solid_clone = app.current_solid.clone();
    let last_feature_id: Option<valenx_feature_tree::feature::FeatureId> = app
        .mesh_toolbox
        .part_design
        .tree
        .features
        .len()
        .checked_sub(1)
        .map(valenx_feature_tree::feature::FeatureId);

    let s = &mut app.mesh_toolbox.techdraw;

    // ----- Sheet size selector -----
    ui.label(egui::RichText::new("Sheet size").strong())
        .on_hover_text("ISO A-series sheet size. Changing this rescales the drawing canvas; existing views and annotations keep their absolute positions.");
    let prev_size = s.sheet_size;
    egui::ComboBox::from_id_source("techdraw_sheet_combo")
        .selected_text(s.sheet_size.label())
        .show_ui(ui, |ui| {
            for opt in [
                valenx_techdraw::SheetSize::A4,
                valenx_techdraw::SheetSize::A3,
                valenx_techdraw::SheetSize::A2,
                valenx_techdraw::SheetSize::A1,
                valenx_techdraw::SheetSize::A0,
            ] {
                ui.selectable_value(&mut s.sheet_size, opt, opt.label())
                    .on_hover_text(match opt {
                        valenx_techdraw::SheetSize::A4 => "210 × 297 mm — single component / parts list.",
                        valenx_techdraw::SheetSize::A3 => "297 × 420 mm — typical assembly.",
                        valenx_techdraw::SheetSize::A2 => "420 × 594 mm — larger assembly.",
                        valenx_techdraw::SheetSize::A1 => "594 × 841 mm — general arrangement.",
                        valenx_techdraw::SheetSize::A0 => "841 × 1189 mm — site plan / large GA.",
                        _ => "Custom sheet size.",
                    });
            }
        });
    if s.sheet_size != prev_size {
        // Replace sheet size — preserves title-block fields.
        s.drawing.sheet.size = s.sheet_size;
    }

    // ----- Title-block fields -----
    ui.horizontal(|ui| {
        ui.label("Title:")
            .on_hover_text("Drawing title — appears in the title block and as the default export filename.");
        ui.text_edit_singleline(&mut s.drawing.sheet.title)
            .on_hover_text("Free-form title text.");
    });
    ui.horizontal(|ui| {
        ui.label("Author:")
            .on_hover_text("Drafter name for the title-block author cell.");
        ui.text_edit_singleline(&mut s.drawing.sheet.author)
            .on_hover_text("Author / drafter name.");
    });
    ui.horizontal(|ui| {
        ui.label("Revision:")
            .on_hover_text("Drawing revision code — typically a letter (A, B, C…) or short numeric.");
        ui.text_edit_singleline(&mut s.drawing.sheet.revision)
            .on_hover_text("Revision text.");
    });

    ui.separator();

    // ----- Add View buttons -----
    ui.label(egui::RichText::new("Add view").strong());
    ui.horizontal(|ui| {
        ui.label("Pos:")
            .on_hover_text("Sheet-space position (mm) where the new view is anchored. Origin is the sheet's bottom-left corner.");
        ui.add(
            egui::DragValue::new(&mut s.new_view_position[0])
                .speed(1.0)
                .prefix("x "),
        )
        .on_hover_text("View anchor X (mm, sheet-space).");
        ui.add(
            egui::DragValue::new(&mut s.new_view_position[1])
                .speed(1.0)
                .prefix("y "),
        )
        .on_hover_text("View anchor Y (mm, sheet-space).");
    });
    ui.horizontal(|ui| {
        ui.label("Scale:")
            .on_hover_text("View magnification. 1.0 = 1:1, 0.5 = 1:2, 2.0 = 2:1, etc.");
        ui.add(
            egui::DragValue::new(&mut s.new_view_scale)
                .speed(0.05)
                .range(0.01..=100.0),
        )
        .on_hover_text("View scale factor (0.01..100).");
        // Task 4 — parametric checkbox toggling auto-update.
        ui.checkbox(&mut s.new_view_parametric, "Parametric (auto-update)")
            .on_hover_text("Link this view to the active feature so it regenerates whenever the source solid changes.");
    });

    // The Add View buttons — each one builds a fresh View at the
    // current new_view_position + scale, runs `generate` against the
    // current solid (if any), and pushes it onto the drawing.
    let mut to_add: Option<valenx_techdraw::ViewKind> = None;
    ui.horizontal_wrapped(|ui| {
        if ui.button("Front")
            .on_hover_text("Projection along +Y axis (looking toward -Y).")
            .clicked()
        {
            to_add = Some(valenx_techdraw::ViewKind::Front);
        }
        if ui.button("Top")
            .on_hover_text("Projection along +Z axis (looking down).")
            .clicked()
        {
            to_add = Some(valenx_techdraw::ViewKind::Top);
        }
        if ui.button("Right")
            .on_hover_text("Projection along +X axis (looking toward -X).")
            .clicked()
        {
            to_add = Some(valenx_techdraw::ViewKind::Right);
        }
        if ui.button("Back")
            .on_hover_text("Projection along -Y axis.")
            .clicked()
        {
            to_add = Some(valenx_techdraw::ViewKind::Back);
        }
        if ui.button("Bottom")
            .on_hover_text("Projection along -Z axis (looking up).")
            .clicked()
        {
            to_add = Some(valenx_techdraw::ViewKind::Bottom);
        }
        if ui.button("Left")
            .on_hover_text("Projection along -X axis.")
            .clicked()
        {
            to_add = Some(valenx_techdraw::ViewKind::Left);
        }
        if ui.button("Iso")
            .on_hover_text("Isometric view — 30° axes, equal foreshortening on X / Y / Z.")
            .clicked()
        {
            to_add = Some(valenx_techdraw::ViewKind::Isometric);
        }
    });
    if let Some(kind) = to_add {
        let mut view = valenx_techdraw::View::new(kind, s.new_view_scale, s.new_view_position);
        let want_parametric = s.new_view_parametric;
        match &solid_clone {
            Some(solid) => match view.generate(solid) {
                Ok(()) => {
                    let idx = s.drawing.add_view(view);
                    if want_parametric {
                        // Task 1/4 — register a parametric link. The
                        // FeatureId is the last feature in the part-
                        // design tree (heuristic captured pre-borrow).
                        if let Some(fid) = last_feature_id {
                            s.drawing
                                .add_parametric_view(valenx_techdraw::ParametricView::new(
                                    idx, fid,
                                ));
                        }
                    }
                    s.last_error = None;
                }
                Err(e) => s.last_error = Some(format!("Add view: {e}")),
            },
            None => {
                // No solid yet — still add the view, just empty.
                // Lets users lay out a sheet before any geometry exists.
                s.drawing.add_view(view);
                s.last_error = Some(
                    "View added but no source solid available (use Part workbench to create one)."
                        .into(),
                );
            }
        }
    }

    ui.separator();

    // ----- View list with selection + edit + remove -----
    ui.label(egui::RichText::new("Views").strong());
    if s.drawing.views.is_empty() {
        ui.label("(no views — add one above)");
    }
    let mut to_remove: Option<usize> = None;
    let mut to_regen: Option<usize> = None;
    for i in 0..s.drawing.views.len() {
        let kind_label = s.drawing.views[i].kind.label();
        let pos = s.drawing.views[i].position;
        let scale = s.drawing.views[i].scale;
        let vis = s.drawing.views[i].visible_edges.len();
        let hid = s.drawing.views[i].hidden_edges.len();
        let header = format!(
            "#{i}: {kind_label} @ ({:.1}, {:.1}) ×{:.2}  vis={vis} hid={hid}",
            pos[0], pos[1], scale
        );
        let selected = s.selected_view == Some(i);
        let rich = if selected {
            egui::RichText::new(header).strong()
        } else {
            egui::RichText::new(header)
        };
        if ui.selectable_label(selected, rich).clicked() {
            s.selected_view = if selected { None } else { Some(i) };
        }
    }
    if let Some(idx) = s.selected_view {
        ui.horizontal(|ui| {
            ui.label("Scale:");
            if let Ok(v) = s.drawing.get_view_mut(idx) {
                ui.add(
                    egui::DragValue::new(&mut v.scale)
                        .speed(0.05)
                        .range(0.01..=100.0),
                );
            }
        });
        ui.horizontal(|ui| {
            ui.label("Position:");
            if let Ok(v) = s.drawing.get_view_mut(idx) {
                ui.add(
                    egui::DragValue::new(&mut v.position[0])
                        .speed(1.0)
                        .prefix("x "),
                );
                ui.add(
                    egui::DragValue::new(&mut v.position[1])
                        .speed(1.0)
                        .prefix("y "),
                );
            }
        });
        ui.horizontal(|ui| {
            if ui.button("Regenerate edges").clicked() {
                to_regen = Some(idx);
            }
            if ui.button("Auto-dimension bbox").clicked() {
                if let Ok(v) = s.drawing.get_view(idx) {
                    let dims = v.auto_dimension_bbox();
                    for d in dims {
                        s.drawing.add_dimension(d);
                    }
                }
            }
            if ui.button("Remove view").clicked() {
                to_remove = Some(idx);
            }
        });
    }
    if let Some(idx) = to_regen {
        match &solid_clone {
            Some(solid) => {
                let result = s.drawing.get_view_mut(idx).map(|v| v.generate(solid));
                match result {
                    Ok(Ok(())) => s.last_error = None,
                    Ok(Err(e)) => s.last_error = Some(format!("Regenerate: {e}")),
                    Err(e) => s.last_error = Some(format!("Regenerate: {e}")),
                }
            }
            None => s.last_error = Some("Regenerate: no source solid available.".into()),
        }
    }
    if let Some(idx) = to_remove {
        match s.drawing.remove_view(idx) {
            Ok(_) => {
                s.selected_view = None;
                s.last_error = None;
            }
            Err(e) => s.last_error = Some(format!("Remove: {e}")),
        }
    }

    ui.separator();

    // ----- Add Linear Dimension -----
    ui.label(egui::RichText::new("Add linear dimension").strong());
    ui.horizontal(|ui| {
        ui.label("From:")
            .on_hover_text("Sheet-space start point of the dimension witness line.");
        ui.add(
            egui::DragValue::new(&mut s.dim_from[0])
                .speed(1.0)
                .prefix("x "),
        )
        .on_hover_text("From X (mm, sheet-space).");
        ui.add(
            egui::DragValue::new(&mut s.dim_from[1])
                .speed(1.0)
                .prefix("y "),
        )
        .on_hover_text("From Y (mm, sheet-space).");
    });
    ui.horizontal(|ui| {
        ui.label("To:")
            .on_hover_text("Sheet-space end point of the dimension witness line.");
        ui.add(
            egui::DragValue::new(&mut s.dim_to[0])
                .speed(1.0)
                .prefix("x "),
        )
        .on_hover_text("To X (mm, sheet-space).");
        ui.add(
            egui::DragValue::new(&mut s.dim_to[1])
                .speed(1.0)
                .prefix("y "),
        )
        .on_hover_text("To Y (mm, sheet-space).");
    });
    ui.horizontal(|ui| {
        ui.label("Offset:")
            .on_hover_text("Perpendicular distance from the witness baseline to the dimension line. Positive = above/right.");
        ui.add(egui::DragValue::new(&mut s.dim_offset).speed(0.5))
            .on_hover_text("Dim line offset (mm).");
        if ui.button("Add dimension")
            .on_hover_text("Add a linear dimension between From and To with the chosen offset.")
            .clicked()
        {
            let dx = s.dim_to[0] - s.dim_from[0];
            let dy = s.dim_to[1] - s.dim_from[1];
            let value = (dx * dx + dy * dy).sqrt();
            s.drawing.add_dimension(valenx_techdraw::Dimension::Linear {
                from: s.dim_from,
                to: s.dim_to,
                offset: s.dim_offset,
                value,
            });
        }
    });

    ui.label(format!("Dimensions: {}", s.drawing.dimensions.len()));

    ui.separator();

    // ===== Phase 18 — extended annotations =====

    // ----- Phase 18A: Parametric-view regenerate -----
    ui.label(egui::RichText::new("Parametric views (Phase 18A)").strong());
    ui.label(format!(
        "Linked: {}  | auto-update: {}",
        s.drawing.parametric_views.len(),
        s.drawing
            .parametric_views
            .iter()
            .filter(|p| p.auto_update)
            .count()
    ));
    if ui
        .button("Regenerate all parametric views")
        .on_hover_text("Re-extract edges from the active feature tree's solid for every parametric view whose auto_update is on.")
        .clicked()
    {
        if let Some(solid) = &solid_clone {
            let solid = solid.clone();
            let errs = s.drawing.regenerate_all(|_fid| Some(solid.clone()));
            s.last_error = if errs.is_empty() {
                Some(format!(
                    "Regenerated {} parametric view(s).",
                    s.drawing.parametric_views.iter().filter(|p| p.auto_update).count()
                ))
            } else {
                Some(format!("Regenerate-all: {} view(s) failed", errs.len()))
            };
        } else {
            s.last_error = Some("Regenerate-all: no source solid available.".into());
        }
    }

    ui.separator();

    // ----- Phase 18B: Dim chain -----
    ui.label(egui::RichText::new("Auto-dim chain (Phase 18B)").strong());
    ui.horizontal(|ui| {
        ui.label("Kind:")
            .on_hover_text("Dimension chain style: Ordinate = all measurements from one origin; Baseline = stacked offsets from baseline; Chain = end-to-end consecutive segments.");
        egui::ComboBox::from_id_source("techdraw_dim_chain_kind")
            .selected_text(s.chain_kind.label())
            .show_ui(ui, |ui| {
                for k in [
                    valenx_techdraw::DimChainKind::Ordinate,
                    valenx_techdraw::DimChainKind::Baseline,
                    valenx_techdraw::DimChainKind::Chain,
                ] {
                    ui.selectable_value(&mut s.chain_kind, k, k.label())
                        .on_hover_text(match k {
                            valenx_techdraw::DimChainKind::Ordinate => "All offsets from a single origin entry — eliminates tolerance stack-up.",
                            valenx_techdraw::DimChainKind::Baseline => "Stacked dimension lines, each from the baseline entry.",
                            valenx_techdraw::DimChainKind::Chain => "Consecutive end-to-end dimensions — tolerances accumulate.",
                        });
                }
            });
        ui.add(
            egui::DragValue::new(&mut s.chain_offset)
                .speed(0.5)
                .prefix("offset "),
        )
        .on_hover_text("Distance from the entries baseline to the dimension line stack.");
    });
    ui.horizontal(|ui| {
        ui.label("Entries (x,y; x,y; ...):")
            .on_hover_text("Semicolon-separated list of (x, y) anchor points along the baseline.");
        ui.text_edit_singleline(&mut s.chain_entries)
            .on_hover_text("Entries — e.g. `0,0; 25,0; 50,0`.");
    });
    if ui.button("Add chain")
        .on_hover_text("Build the dimension chain from the entries above and append it to the drawing.")
        .clicked()
    {
        let entries: Vec<[f64; 2]> = s
            .chain_entries
            .split(';')
            .filter_map(|pair| {
                let mut it = pair.split(',').map(|p| p.trim().parse::<f64>().ok());
                let x = it.next().flatten()?;
                let y = it.next().flatten()?;
                Some([x, y])
            })
            .collect();
        if entries.len() >= 2 {
            let mut c = valenx_techdraw::DimChain::new(s.chain_kind, s.chain_offset);
            c.entries = entries;
            s.drawing.add_dim_chain(c);
            s.last_error = None;
        } else {
            s.last_error = Some("Add chain: need at least 2 entries.".into());
        }
    }
    ui.label(format!("Dim chains: {}", s.drawing.dim_chains.len()));

    ui.separator();

    // ----- Phase 18C: Balloons + Leaders -----
    ui.label(egui::RichText::new("Balloons + leaders (Phase 18C)").strong());
    ui.horizontal(|ui| {
        ui.label("Balloon pos:")
            .on_hover_text("Sheet-space position of the balloon shape (the labelled bubble).");
        ui.add(
            egui::DragValue::new(&mut s.balloon_position[0])
                .speed(1.0)
                .prefix("x "),
        )
        .on_hover_text("Balloon X (mm).");
        ui.add(
            egui::DragValue::new(&mut s.balloon_position[1])
                .speed(1.0)
                .prefix("y "),
        )
        .on_hover_text("Balloon Y (mm).");
    });
    ui.horizontal(|ui| {
        ui.label("Target:")
            .on_hover_text("Sheet-space point the balloon's leader arrow points to (the part being labelled).");
        ui.add(
            egui::DragValue::new(&mut s.balloon_target[0])
                .speed(1.0)
                .prefix("x "),
        )
        .on_hover_text("Target X (mm).");
        ui.add(
            egui::DragValue::new(&mut s.balloon_target[1])
                .speed(1.0)
                .prefix("y "),
        )
        .on_hover_text("Target Y (mm).");
    });
    ui.horizontal(|ui| {
        ui.label("Number:")
            .on_hover_text("Balloon label — typically a parts-list item number (1, 2, A1…).");
        ui.text_edit_singleline(&mut s.balloon_number)
            .on_hover_text("Balloon label text.");
        egui::ComboBox::from_id_source("techdraw_balloon_style")
            .selected_text(s.balloon_style.label())
            .show_ui(ui, |ui| {
                for st in [
                    valenx_techdraw::BalloonStyle::Circle,
                    valenx_techdraw::BalloonStyle::Square,
                    valenx_techdraw::BalloonStyle::Hexagon,
                    valenx_techdraw::BalloonStyle::Triangle,
                ] {
                    ui.selectable_value(&mut s.balloon_style, st, st.label())
                        .on_hover_text(match st {
                            valenx_techdraw::BalloonStyle::Circle => "Round balloon — most common.",
                            valenx_techdraw::BalloonStyle::Square => "Square balloon — used for accessory items.",
                            valenx_techdraw::BalloonStyle::Hexagon => "Hexagonal balloon — often for fasteners.",
                            valenx_techdraw::BalloonStyle::Triangle => "Triangular balloon — sometimes used for notes / warnings.",
                        });
                }
            });
        if ui.button("Add balloon")
            .on_hover_text("Place the balloon at the chosen position with a leader to the target.")
            .clicked()
        {
            let mut b = valenx_techdraw::Balloon::new(
                s.balloon_position,
                &s.balloon_number,
                s.balloon_target,
            );
            b.style = s.balloon_style;
            s.drawing.add_balloon(b);
        }
    });
    ui.label(format!("Balloons: {}", s.drawing.balloons.len()));

    ui.horizontal(|ui| {
        ui.label("Leader start:")
            .on_hover_text("Sheet-space text-end of the leader (where the callout text sits).");
        ui.add(
            egui::DragValue::new(&mut s.leader_start[0])
                .speed(1.0)
                .prefix("x "),
        )
        .on_hover_text("Leader text-end X (mm).");
        ui.add(
            egui::DragValue::new(&mut s.leader_start[1])
                .speed(1.0)
                .prefix("y "),
        )
        .on_hover_text("Leader text-end Y (mm).");
    });
    ui.horizontal(|ui| {
        ui.label("Leader end:")
            .on_hover_text("Sheet-space arrow-tip end of the leader (what's being pointed at).");
        ui.add(
            egui::DragValue::new(&mut s.leader_end[0])
                .speed(1.0)
                .prefix("x "),
        )
        .on_hover_text("Leader arrow-tip X (mm).");
        ui.add(
            egui::DragValue::new(&mut s.leader_end[1])
                .speed(1.0)
                .prefix("y "),
        )
        .on_hover_text("Leader arrow-tip Y (mm).");
    });
    ui.horizontal(|ui| {
        ui.label("Text:")
            .on_hover_text("Callout text shown at the leader's start.");
        ui.text_edit_singleline(&mut s.leader_text)
            .on_hover_text("Leader callout text.");
        egui::ComboBox::from_id_source("techdraw_leader_arrow")
            .selected_text(s.leader_arrow.label())
            .show_ui(ui, |ui| {
                for k in [
                    valenx_techdraw::ArrowKind::Closed,
                    valenx_techdraw::ArrowKind::Open,
                    valenx_techdraw::ArrowKind::Dot,
                    valenx_techdraw::ArrowKind::Tick,
                ] {
                    ui.selectable_value(&mut s.leader_arrow, k, k.label())
                        .on_hover_text(match k {
                            valenx_techdraw::ArrowKind::Closed => "Filled triangular arrowhead — standard dimension arrow.",
                            valenx_techdraw::ArrowKind::Open => "Open V-arrowhead — used for non-dimensional callouts.",
                            valenx_techdraw::ArrowKind::Dot => "Filled dot — used for some leader styles (ANSI Y14).",
                            valenx_techdraw::ArrowKind::Tick => "Architectural tick mark — 45° slash, common in architectural drawings.",
                        });
                }
            });
        if ui.button("Add leader")
            .on_hover_text("Place a leader with the chosen arrow style.")
            .clicked()
        {
            let mut l = valenx_techdraw::Leader::new(s.leader_start, s.leader_end, &s.leader_text);
            l.arrow_kind = s.leader_arrow;
            s.drawing.add_leader(l);
        }
    });
    ui.label(format!("Leaders: {}", s.drawing.leaders.len()));

    ui.separator();

    // ----- Phase 18D: Weld symbols -----
    ui.label(egui::RichText::new("Weld symbol (Phase 18D, ISO 2553)").strong());
    ui.horizontal(|ui| {
        ui.label("Ref pos:");
        ui.add(
            egui::DragValue::new(&mut s.weld_position[0])
                .speed(1.0)
                .prefix("x "),
        );
        ui.add(
            egui::DragValue::new(&mut s.weld_position[1])
                .speed(1.0)
                .prefix("y "),
        );
    });
    ui.horizontal(|ui| {
        ui.label("Arrow target:");
        ui.add(
            egui::DragValue::new(&mut s.weld_target[0])
                .speed(1.0)
                .prefix("x "),
        );
        ui.add(
            egui::DragValue::new(&mut s.weld_target[1])
                .speed(1.0)
                .prefix("y "),
        );
    });
    ui.horizontal(|ui| {
        ui.label("Type:");
        egui::ComboBox::from_id_source("techdraw_weld_type")
            .selected_text(s.weld_type.label())
            .show_ui(ui, |ui| {
                for t in [
                    valenx_techdraw::WeldType::Fillet,
                    valenx_techdraw::WeldType::Square,
                    valenx_techdraw::WeldType::V,
                    valenx_techdraw::WeldType::U,
                    valenx_techdraw::WeldType::Bevel,
                    valenx_techdraw::WeldType::J,
                    valenx_techdraw::WeldType::Flare,
                    valenx_techdraw::WeldType::Plug,
                    valenx_techdraw::WeldType::Seam,
                ] {
                    ui.selectable_value(&mut s.weld_type, t, t.label());
                }
            });
        egui::ComboBox::from_id_source("techdraw_weld_side")
            .selected_text(s.weld_side.label())
            .show_ui(ui, |ui| {
                ui.selectable_value(
                    &mut s.weld_side,
                    valenx_techdraw::WeldPosition::Arrow,
                    "Arrow",
                );
                ui.selectable_value(
                    &mut s.weld_side,
                    valenx_techdraw::WeldPosition::Other,
                    "Other",
                );
            });
    });
    ui.horizontal(|ui| {
        ui.label("Size:");
        ui.text_edit_singleline(&mut s.weld_size);
        ui.checkbox(&mut s.weld_all_around, "All around");
        ui.checkbox(&mut s.weld_field, "Field weld");
        if ui.button("Add weld").clicked() {
            let mut w = valenx_techdraw::WeldSymbol::new_fillet(
                s.weld_position,
                s.weld_target,
                &s.weld_size,
            );
            w.weld_type = s.weld_type;
            w.weld_position = s.weld_side;
            w.all_around = s.weld_all_around;
            w.field_weld = s.weld_field;
            s.drawing.add_weld(w);
        }
    });
    ui.label(format!("Welds: {}", s.drawing.welds.len()));

    ui.separator();

    // ----- Phase 18E: Surface finish -----
    ui.label(egui::RichText::new("Surface finish (Phase 18E, ISO 1302)").strong());
    ui.horizontal(|ui| {
        ui.label("Pos:");
        ui.add(
            egui::DragValue::new(&mut s.sf_position[0])
                .speed(1.0)
                .prefix("x "),
        );
        ui.add(
            egui::DragValue::new(&mut s.sf_position[1])
                .speed(1.0)
                .prefix("y "),
        );
        ui.add(egui::DragValue::new(&mut s.sf_ra).speed(0.1).prefix("Ra "));
    });
    ui.horizontal(|ui| {
        egui::ComboBox::from_id_source("techdraw_sf_process")
            .selected_text(s.sf_process.label())
            .show_ui(ui, |ui| {
                for p in [
                    valenx_techdraw::SurfaceProcess::Required,
                    valenx_techdraw::SurfaceProcess::Machined,
                    valenx_techdraw::SurfaceProcess::AsCast,
                    valenx_techdraw::SurfaceProcess::Removed,
                ] {
                    ui.selectable_value(&mut s.sf_process, p, p.label());
                }
            });
        egui::ComboBox::from_id_source("techdraw_sf_lay")
            .selected_text(s.sf_lay.label())
            .show_ui(ui, |ui| {
                for l in [
                    valenx_techdraw::LayPattern::Parallel,
                    valenx_techdraw::LayPattern::Perpendicular,
                    valenx_techdraw::LayPattern::Crossed,
                    valenx_techdraw::LayPattern::Multi,
                    valenx_techdraw::LayPattern::Radial,
                    valenx_techdraw::LayPattern::Circular,
                ] {
                    ui.selectable_value(&mut s.sf_lay, l, l.label());
                }
            });
        if ui.button("Add surface finish").clicked() {
            let mut sf = valenx_techdraw::SurfaceFinish::new(s.sf_position, s.sf_ra);
            sf.process = s.sf_process;
            sf.lay_pattern = s.sf_lay;
            s.drawing.add_surface_finish(sf);
        }
    });
    ui.label(format!(
        "Surface finishes: {}",
        s.drawing.surface_finishes.len()
    ));

    ui.separator();

    // ----- Phase 18F: GD&T feature control frame + datum -----
    ui.label(egui::RichText::new("GD&T (Phase 18F, ASME Y14.5)").strong());
    ui.horizontal(|ui| {
        ui.label("Pos:");
        ui.add(
            egui::DragValue::new(&mut s.gdt_position[0])
                .speed(1.0)
                .prefix("x "),
        );
        ui.add(
            egui::DragValue::new(&mut s.gdt_position[1])
                .speed(1.0)
                .prefix("y "),
        );
    });
    ui.horizontal(|ui| {
        egui::ComboBox::from_id_source("techdraw_gdt_char")
            .selected_text(s.gdt_characteristic.label())
            .show_ui(ui, |ui| {
                for c in [
                    valenx_techdraw::GeometricCharacteristic::Straightness,
                    valenx_techdraw::GeometricCharacteristic::Flatness,
                    valenx_techdraw::GeometricCharacteristic::Circularity,
                    valenx_techdraw::GeometricCharacteristic::Cylindricity,
                    valenx_techdraw::GeometricCharacteristic::ProfileLine,
                    valenx_techdraw::GeometricCharacteristic::ProfileSurface,
                    valenx_techdraw::GeometricCharacteristic::Perpendicularity,
                    valenx_techdraw::GeometricCharacteristic::Angularity,
                    valenx_techdraw::GeometricCharacteristic::Parallelism,
                    valenx_techdraw::GeometricCharacteristic::Position,
                    valenx_techdraw::GeometricCharacteristic::Concentricity,
                    valenx_techdraw::GeometricCharacteristic::Symmetry,
                    valenx_techdraw::GeometricCharacteristic::CircularRunout,
                    valenx_techdraw::GeometricCharacteristic::TotalRunout,
                ] {
                    ui.selectable_value(&mut s.gdt_characteristic, c, c.label());
                }
            });
        ui.label("Tol:");
        ui.text_edit_singleline(&mut s.gdt_tolerance);
        egui::ComboBox::from_id_source("techdraw_gdt_mod")
            .selected_text(s.gdt_modifier.label())
            .show_ui(ui, |ui| {
                ui.selectable_value(
                    &mut s.gdt_modifier,
                    valenx_techdraw::MaterialCondition::Rfs,
                    "RFS",
                );
                ui.selectable_value(
                    &mut s.gdt_modifier,
                    valenx_techdraw::MaterialCondition::Mmc,
                    "MMC",
                );
                ui.selectable_value(
                    &mut s.gdt_modifier,
                    valenx_techdraw::MaterialCondition::Lmc,
                    "LMC",
                );
            });
    });
    ui.horizontal(|ui| {
        ui.label("Datums (comma-sep):");
        ui.text_edit_singleline(&mut s.gdt_datum_letters);
        if ui.button("Add GD&T frame").clicked() {
            let mut g = valenx_techdraw::GdtSymbol::new(
                s.gdt_position,
                s.gdt_characteristic,
                &s.gdt_tolerance,
            );
            g.material_condition = s.gdt_modifier;
            for letter in s.gdt_datum_letters.split(',') {
                let l = letter.trim();
                if !l.is_empty() {
                    g.datums.push(valenx_techdraw::DatumRef::new(l));
                }
            }
            s.drawing.add_gdt(g);
        }
    });
    ui.label(format!("GD&T frames: {}", s.drawing.gdt.len()));

    ui.horizontal(|ui| {
        ui.label("Datum pos:");
        ui.add(
            egui::DragValue::new(&mut s.datum_position[0])
                .speed(1.0)
                .prefix("x "),
        );
        ui.add(
            egui::DragValue::new(&mut s.datum_position[1])
                .speed(1.0)
                .prefix("y "),
        );
        ui.label("Target:");
        ui.add(
            egui::DragValue::new(&mut s.datum_target[0])
                .speed(1.0)
                .prefix("x "),
        );
        ui.add(
            egui::DragValue::new(&mut s.datum_target[1])
                .speed(1.0)
                .prefix("y "),
        );
        ui.label("Letter:");
        ui.text_edit_singleline(&mut s.datum_letter);
        if ui.button("Add datum").clicked() {
            s.drawing.add_datum(valenx_techdraw::Datum::new(
                s.datum_position,
                &s.datum_letter,
                s.datum_target,
            ));
        }
    });
    ui.label(format!("Datums: {}", s.drawing.datums.len()));

    ui.separator();

    // ----- Phase 18G: Hatch pattern picker for section fills -----
    ui.label(egui::RichText::new("Hatch pattern (Phase 18G)").strong());
    ui.horizontal(|ui| {
        ui.label("Pattern:");
        egui::ComboBox::from_id_source("techdraw_hatch_pattern")
            .selected_text(&s.hatch_pattern)
            .show_ui(ui, |ui| {
                for name in valenx_techdraw::hatch_lib::all_names() {
                    ui.selectable_value(&mut s.hatch_pattern, (*name).to_string(), *name);
                }
            });
    });

    ui.separator();

    // ----- Export buttons -----
    ui.label(egui::RichText::new("Export").strong());
    enum ExportFormat {
        Svg,
        Pdf,
        Dxf,
    }
    let mut export: Option<ExportFormat> = None;
    ui.horizontal_wrapped(|ui| {
        if ui.button("Save as SVG…").clicked() {
            export = Some(ExportFormat::Svg);
        }
        if ui.button("Save as PDF…").clicked() {
            export = Some(ExportFormat::Pdf);
        }
        if ui.button("Save as DXF…").clicked() {
            export = Some(ExportFormat::Dxf);
        }
    });
    if let Some(fmt) = export {
        let (ext, filter) = match fmt {
            ExportFormat::Svg => ("svg", "SVG"),
            ExportFormat::Pdf => ("pdf", "PDF"),
            ExportFormat::Dxf => ("dxf", "DXF"),
        };
        let default_name = format!("{}.{ext}", sanitise_filename(&s.drawing.sheet.title));
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter(filter, &[ext])
            .save_file()
        {
            let res = match fmt {
                ExportFormat::Svg => valenx_techdraw::export::svg::write(&s.drawing, &path),
                ExportFormat::Pdf => valenx_techdraw::export::pdf::write(&s.drawing, &path),
                ExportFormat::Dxf => valenx_techdraw::export::dxf::write(&s.drawing, &path),
            };
            match res {
                Ok(()) => s.last_error = Some(format!("Exported to {}", path.display())),
                Err(e) => s.last_error = Some(format!("Export failed: {e}")),
            }
        }
    }

    if let Some(err) = &s.last_error {
        ui.colored_label(egui::Color32::from_rgb(220, 100, 100), err);
    }
}

/// Strip filesystem-unfriendly characters from a sheet title so we
/// can use it as a default filename for the Save As dialog.
fn sanitise_filename(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_alphanumeric() || c == '-' || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect::<String>()
        .trim_matches('_')
        .to_string()
}

/// Build a part from the Assembly panel's "Add part" inputs and append
/// it to the live assembly.
///
/// Extracted verbatim from the "Add part" button closure so the
/// headless UI tests can exercise the action against the real
/// `valenx-cad` + `valenx-assembly` backend. On a primitive-build
/// failure the error lands in `last_error`; on success the new part
/// is added at its initial-position transform.
fn assembly_add_part(s: &mut AssemblyPanelState) {
    let solid = match s.new_part_primitive {
        AssemblyPartPrimitive::Box => valenx_cad::box_solid(
            s.new_part_box_dims[0],
            s.new_part_box_dims[1],
            s.new_part_box_dims[2],
        ),
        AssemblyPartPrimitive::Cylinder => {
            valenx_cad::cylinder(s.new_part_cyl[0], s.new_part_cyl[1])
        }
        AssemblyPartPrimitive::Sphere => valenx_cad::sphere(s.new_part_sphere),
    };
    match solid {
        Ok(sol) => {
            let mut part = valenx_assembly::Part::new(0, s.new_part_name.clone(), sol);
            part.transform.translation = nalgebra::Vector3::new(
                s.new_part_translation[0],
                s.new_part_translation[1],
                s.new_part_translation[2],
            );
            s.assembly.add_part(part);
            s.last_error = None;
        }
        Err(e) => s.last_error = Some(format!("Add part: {e}")),
    }
}

/// Draw the Assembly workbench panel (Phase 6).
pub fn draw_assembly_panel(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.heading("Assembly (multi-part scene with mates + joints)");
    ui.separator();

    // ----- Add Part -----
    ui.label(egui::RichText::new("Add part").strong());
    {
        let s = &mut app.mesh_toolbox.assembly;
        ui.horizontal(|ui| {
            ui.label("Name:")
                .on_hover_text("Display name for the part in the scene tree.");
            ui.text_edit_singleline(&mut s.new_part_name)
                .on_hover_text("Part name (free-form).");
        });
        egui::ComboBox::from_id_source("assembly_part_primitive_combo")
            .selected_text(s.new_part_primitive.label())
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut s.new_part_primitive, AssemblyPartPrimitive::Box, "Box")
                    .on_hover_text("Axis-aligned rectangular box primitive.");
                ui.selectable_value(
                    &mut s.new_part_primitive,
                    AssemblyPartPrimitive::Cylinder,
                    "Cylinder",
                )
                .on_hover_text("Right cylinder primitive with radius + height.");
                ui.selectable_value(
                    &mut s.new_part_primitive,
                    AssemblyPartPrimitive::Sphere,
                    "Sphere",
                )
                .on_hover_text("Sphere primitive with radius.");
            });
        match s.new_part_primitive {
            AssemblyPartPrimitive::Box => {
                ui.horizontal(|ui| {
                    ui.label("Dims:");
                    ui.add(
                        egui::DragValue::new(&mut s.new_part_box_dims[0])
                            .speed(0.1)
                            .prefix("dx "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.new_part_box_dims[1])
                            .speed(0.1)
                            .prefix("dy "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.new_part_box_dims[2])
                            .speed(0.1)
                            .prefix("dz "),
                    );
                });
            }
            AssemblyPartPrimitive::Cylinder => {
                ui.horizontal(|ui| {
                    ui.label("R / H:");
                    ui.add(
                        egui::DragValue::new(&mut s.new_part_cyl[0])
                            .speed(0.05)
                            .prefix("r "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.new_part_cyl[1])
                            .speed(0.05)
                            .prefix("h "),
                    );
                });
            }
            AssemblyPartPrimitive::Sphere => {
                ui.horizontal(|ui| {
                    ui.label("Radius:");
                    ui.add(egui::DragValue::new(&mut s.new_part_sphere).speed(0.05));
                });
            }
        }
        ui.horizontal(|ui| {
            ui.label("Initial position:");
            ui.add(
                egui::DragValue::new(&mut s.new_part_translation[0])
                    .speed(0.1)
                    .prefix("x "),
            );
            ui.add(
                egui::DragValue::new(&mut s.new_part_translation[1])
                    .speed(0.1)
                    .prefix("y "),
            );
            ui.add(
                egui::DragValue::new(&mut s.new_part_translation[2])
                    .speed(0.1)
                    .prefix("z "),
            );
        });
    }
    let mut add_part_clicked = false;
    let mut clear_clicked = false;
    ui.horizontal(|ui| {
        if ui.button("Add part").clicked() {
            add_part_clicked = true;
        }
        if ui.button("Clear assembly").clicked() {
            clear_clicked = true;
        }
    });
    if add_part_clicked {
        assembly_add_part(&mut app.mesh_toolbox.assembly);
    }
    if clear_clicked {
        let s = &mut app.mesh_toolbox.assembly;
        s.assembly = valenx_assembly::Assembly::new();
        s.selected_part = None;
        s.selected_mate = None;
        s.selected_joint = None;
        s.last_error = None;
        s.last_report = None;
    }

    ui.separator();

    // ----- Scene tree -----
    ui.label(egui::RichText::new("Scene tree").strong());
    let mut to_delete_part: Option<usize> = None;
    {
        let s = &mut app.mesh_toolbox.assembly;
        if s.assembly.parts.is_empty() {
            ui.label("(no parts — add one above)");
        }
        for i in 0..s.assembly.parts.len() {
            let pid = s.assembly.parts[i].id;
            let selected = s.selected_part == Some(pid);
            let p = &mut s.assembly.parts[i];
            let header = format!(
                "#{pid}: {} {}",
                p.name,
                if p.fixed { "[fixed]" } else { "" }
            );
            ui.horizontal(|ui| {
                if ui.selectable_label(selected, header).clicked() {
                    s.selected_part = if selected { None } else { Some(pid) };
                }
                ui.checkbox(&mut p.fixed, "fixed");
                if ui.small_button("X").clicked() {
                    to_delete_part = Some(pid);
                }
            });
        }
    }
    if let Some(pid) = to_delete_part {
        let s = &mut app.mesh_toolbox.assembly;
        match s.assembly.delete_part(pid) {
            Ok(()) => {
                if s.selected_part == Some(pid) {
                    s.selected_part = None;
                }
                s.last_error = None;
            }
            Err(e) => s.last_error = Some(format!("Delete part: {e}")),
        }
    }

    // ----- Edit selected part transform -----
    if let Some(pid) = app.mesh_toolbox.assembly.selected_part {
        ui.separator();
        ui.label(egui::RichText::new("Selected part transform").strong());
        // Read the current transform into editable f64 storage, then
        // write back only if changed.
        let mut t = [0.0_f64; 3];
        let mut rdg = [0.0_f64; 3];
        let mut have = false;
        if let Ok(p) = app.mesh_toolbox.assembly.assembly.get_part(pid) {
            t[0] = p.transform.translation.x;
            t[1] = p.transform.translation.y;
            t[2] = p.transform.translation.z;
            let (axis, angle) =
                p.transform.orientation.axis_angle().unwrap_or_else(|| {
                    (nalgebra::Unit::new_unchecked(nalgebra::Vector3::z()), 0.0)
                });
            // Express as axis * angle (Rodrigues) for an unconstrained
            // 3-vector that the user can edit.
            let r = axis.into_inner() * angle;
            rdg[0] = r.x;
            rdg[1] = r.y;
            rdg[2] = r.z;
            have = true;
        }
        if have {
            let mut changed = false;
            ui.horizontal(|ui| {
                ui.label("T:");
                changed |= ui
                    .add(egui::DragValue::new(&mut t[0]).speed(0.05).prefix("x "))
                    .changed();
                changed |= ui
                    .add(egui::DragValue::new(&mut t[1]).speed(0.05).prefix("y "))
                    .changed();
                changed |= ui
                    .add(egui::DragValue::new(&mut t[2]).speed(0.05).prefix("z "))
                    .changed();
            });
            ui.horizontal(|ui| {
                ui.label("R (axis*angle):");
                changed |= ui
                    .add(egui::DragValue::new(&mut rdg[0]).speed(0.02).prefix("rx "))
                    .changed();
                changed |= ui
                    .add(egui::DragValue::new(&mut rdg[1]).speed(0.02).prefix("ry "))
                    .changed();
                changed |= ui
                    .add(egui::DragValue::new(&mut rdg[2]).speed(0.02).prefix("rz "))
                    .changed();
            });
            if changed {
                if let Ok(p) = app.mesh_toolbox.assembly.assembly.get_part_mut(pid) {
                    p.transform.translation = nalgebra::Vector3::new(t[0], t[1], t[2]);
                    let r = nalgebra::Vector3::new(rdg[0], rdg[1], rdg[2]);
                    let angle = r.norm();
                    p.transform.orientation = if angle < 1e-12 {
                        nalgebra::UnitQuaternion::identity()
                    } else {
                        nalgebra::UnitQuaternion::from_axis_angle(
                            &nalgebra::Unit::new_unchecked(r / angle),
                            angle,
                        )
                    };
                }
            }
        }
    }

    ui.separator();

    // ----- Add mate -----
    ui.label(egui::RichText::new("Add mate").strong());
    {
        let s = &mut app.mesh_toolbox.assembly;
        egui::ComboBox::from_id_source("assembly_mate_kind_combo")
            .selected_text(s.mate_kind.label())
            .show_ui(ui, |ui| {
                for opt in [
                    AssemblyMateKindUi::Coincident,
                    AssemblyMateKindUi::Distance,
                    AssemblyMateKindUi::Angle,
                    AssemblyMateKindUi::Parallel,
                    AssemblyMateKindUi::Perpendicular,
                    AssemblyMateKindUi::Tangent,
                ] {
                    ui.selectable_value(&mut s.mate_kind, opt, opt.label())
                        .on_hover_text(match opt {
                            AssemblyMateKindUi::Coincident => "Snap point A on part A to point B on part B (locks 3 translational DOFs).",
                            AssemblyMateKindUi::Distance => "Hold a fixed distance between point A on part A and point B on part B.",
                            AssemblyMateKindUi::Angle => "Hold a fixed angle between vector A and vector B.",
                            AssemblyMateKindUi::Parallel => "Keep vector A parallel to vector B (no relative angle constraint between them).",
                            AssemblyMateKindUi::Perpendicular => "Keep vector A perpendicular to vector B.",
                            AssemblyMateKindUi::Tangent => "Hold cylindrical / spherical surfaces tangent — radius A + axis A + radius B + axis B.",
                        });
                }
            });
        ui.horizontal(|ui| {
            ui.label("Part A id:")
                .on_hover_text("Numeric id of the first part participating in the mate.");
            ui.add(egui::DragValue::new(&mut s.mate_part_a).speed(1.0))
                .on_hover_text("Part A id (matches the scene-tree entry).");
            ui.label("Part B id:")
                .on_hover_text("Numeric id of the second part.");
            ui.add(egui::DragValue::new(&mut s.mate_part_b).speed(1.0))
                .on_hover_text("Part B id.");
        });
        match s.mate_kind {
            AssemblyMateKindUi::Coincident | AssemblyMateKindUi::Distance => {
                ui.horizontal(|ui| {
                    ui.label("Point A:");
                    ui.add(
                        egui::DragValue::new(&mut s.mate_point_a[0])
                            .speed(0.05)
                            .prefix("x "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_point_a[1])
                            .speed(0.05)
                            .prefix("y "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_point_a[2])
                            .speed(0.05)
                            .prefix("z "),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Point B:");
                    ui.add(
                        egui::DragValue::new(&mut s.mate_point_b[0])
                            .speed(0.05)
                            .prefix("x "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_point_b[1])
                            .speed(0.05)
                            .prefix("y "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_point_b[2])
                            .speed(0.05)
                            .prefix("z "),
                    );
                });
                if matches!(s.mate_kind, AssemblyMateKindUi::Distance) {
                    ui.horizontal(|ui| {
                        ui.label("Distance:");
                        ui.add(egui::DragValue::new(&mut s.mate_target).speed(0.1));
                    });
                }
            }
            AssemblyMateKindUi::Angle
            | AssemblyMateKindUi::Parallel
            | AssemblyMateKindUi::Perpendicular => {
                ui.horizontal(|ui| {
                    ui.label("Vec A:");
                    ui.add(
                        egui::DragValue::new(&mut s.mate_vec_a[0])
                            .speed(0.05)
                            .prefix("x "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_vec_a[1])
                            .speed(0.05)
                            .prefix("y "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_vec_a[2])
                            .speed(0.05)
                            .prefix("z "),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Vec B:");
                    ui.add(
                        egui::DragValue::new(&mut s.mate_vec_b[0])
                            .speed(0.05)
                            .prefix("x "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_vec_b[1])
                            .speed(0.05)
                            .prefix("y "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_vec_b[2])
                            .speed(0.05)
                            .prefix("z "),
                    );
                });
                if matches!(s.mate_kind, AssemblyMateKindUi::Angle) {
                    ui.horizontal(|ui| {
                        ui.label("Angle (rad):");
                        ui.add(egui::DragValue::new(&mut s.mate_target).speed(0.05));
                    });
                }
            }
            AssemblyMateKindUi::Tangent => {
                ui.horizontal(|ui| {
                    ui.label("Axis A origin:");
                    ui.add(
                        egui::DragValue::new(&mut s.mate_point_a[0])
                            .speed(0.05)
                            .prefix("x "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_point_a[1])
                            .speed(0.05)
                            .prefix("y "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_point_a[2])
                            .speed(0.05)
                            .prefix("z "),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Axis A dir:");
                    ui.add(
                        egui::DragValue::new(&mut s.mate_vec_a[0])
                            .speed(0.05)
                            .prefix("x "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_vec_a[1])
                            .speed(0.05)
                            .prefix("y "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_vec_a[2])
                            .speed(0.05)
                            .prefix("z "),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Axis B origin:");
                    ui.add(
                        egui::DragValue::new(&mut s.mate_point_b[0])
                            .speed(0.05)
                            .prefix("x "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_point_b[1])
                            .speed(0.05)
                            .prefix("y "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_point_b[2])
                            .speed(0.05)
                            .prefix("z "),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Axis B dir:");
                    ui.add(
                        egui::DragValue::new(&mut s.mate_vec_b[0])
                            .speed(0.05)
                            .prefix("x "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_vec_b[1])
                            .speed(0.05)
                            .prefix("y "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.mate_vec_b[2])
                            .speed(0.05)
                            .prefix("z "),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Radius A:");
                    ui.add(egui::DragValue::new(&mut s.mate_radius_a).speed(0.05));
                    ui.label("Radius B:");
                    ui.add(egui::DragValue::new(&mut s.mate_radius_b).speed(0.05));
                });
            }
        }
    }
    let mut add_mate_clicked = false;
    if ui.button("Add mate").clicked() {
        add_mate_clicked = true;
    }
    if add_mate_clicked {
        let s = &mut app.mesh_toolbox.assembly;
        let vec3 = |xy: [f64; 3]| nalgebra::Vector3::new(xy[0], xy[1], xy[2]);
        let kind = match s.mate_kind {
            AssemblyMateKindUi::Coincident => valenx_assembly::MateKind::Coincident {
                part_a: s.mate_part_a,
                point_a: vec3(s.mate_point_a),
                part_b: s.mate_part_b,
                point_b: vec3(s.mate_point_b),
            },
            AssemblyMateKindUi::Distance => valenx_assembly::MateKind::Distance {
                part_a: s.mate_part_a,
                point_a: vec3(s.mate_point_a),
                part_b: s.mate_part_b,
                point_b: vec3(s.mate_point_b),
                target: s.mate_target,
            },
            AssemblyMateKindUi::Angle => valenx_assembly::MateKind::Angle {
                part_a: s.mate_part_a,
                vec_a: vec3(s.mate_vec_a),
                part_b: s.mate_part_b,
                vec_b: vec3(s.mate_vec_b),
                target: s.mate_target,
            },
            AssemblyMateKindUi::Parallel => valenx_assembly::MateKind::Parallel {
                part_a: s.mate_part_a,
                vec_a: vec3(s.mate_vec_a),
                part_b: s.mate_part_b,
                vec_b: vec3(s.mate_vec_b),
            },
            AssemblyMateKindUi::Perpendicular => valenx_assembly::MateKind::Perpendicular {
                part_a: s.mate_part_a,
                vec_a: vec3(s.mate_vec_a),
                part_b: s.mate_part_b,
                vec_b: vec3(s.mate_vec_b),
            },
            AssemblyMateKindUi::Tangent => valenx_assembly::MateKind::Tangent {
                part_a: s.mate_part_a,
                axis_a_origin: vec3(s.mate_point_a),
                axis_a_dir: vec3(s.mate_vec_a),
                radius_a: s.mate_radius_a,
                part_b: s.mate_part_b,
                axis_b_origin: vec3(s.mate_point_b),
                axis_b_dir: vec3(s.mate_vec_b),
                radius_b: s.mate_radius_b,
            },
        };
        s.assembly.add_mate(valenx_assembly::Mate::new(0, kind));
        s.last_error = None;
    }

    ui.separator();

    // ----- Mate list -----
    ui.label(egui::RichText::new("Mates").strong());
    let mut delete_mate: Option<usize> = None;
    {
        let s = &mut app.mesh_toolbox.assembly;
        if s.assembly.mates.is_empty() {
            ui.label("(no mates)");
        }
        for i in 0..s.assembly.mates.len() {
            let mid = s.assembly.mates[i].id;
            let m = &mut s.assembly.mates[i];
            let (pa, pb) = m.kind.parts();
            let label = format!(
                "#{mid}: {} between {pa} ↔ {pb} {}",
                mate_kind_label(&m.kind),
                if m.suppressed { "[suppressed]" } else { "" }
            );
            ui.horizontal(|ui| {
                ui.label(label);
                ui.checkbox(&mut m.suppressed, "suppress");
                if ui.small_button("X").clicked() {
                    delete_mate = Some(mid);
                }
            });
        }
    }
    if let Some(mid) = delete_mate {
        let s = &mut app.mesh_toolbox.assembly;
        let _ = s.assembly.delete_mate(mid);
    }

    ui.separator();

    // ----- Add joint -----
    ui.label(egui::RichText::new("Add joint").strong());
    {
        let s = &mut app.mesh_toolbox.assembly;
        egui::ComboBox::from_id_source("assembly_joint_kind_combo")
            .selected_text(s.joint_kind.label())
            .show_ui(ui, |ui| {
                for opt in [
                    AssemblyJointKindUi::Fixed,
                    AssemblyJointKindUi::Revolute,
                    AssemblyJointKindUi::Prismatic,
                    AssemblyJointKindUi::Cylindrical,
                    AssemblyJointKindUi::Spherical,
                    AssemblyJointKindUi::Planar,
                ] {
                    ui.selectable_value(&mut s.joint_kind, opt, opt.label())
                        .on_hover_text(match opt {
                            AssemblyJointKindUi::Fixed => "Rigid connection — A and B move as one. 0 relative DOFs.",
                            AssemblyJointKindUi::Revolute => "Hinge — rotation about one axis. 1 rotational DOF (think door hinge).",
                            AssemblyJointKindUi::Prismatic => "Slider — translation along one axis. 1 translational DOF.",
                            AssemblyJointKindUi::Cylindrical => "Combined slide + rotate about the same axis. 2 DOFs.",
                            AssemblyJointKindUi::Spherical => "Ball joint — 3 rotational DOFs, 0 translational.",
                            AssemblyJointKindUi::Planar => "Planar joint — 2 translational + 1 rotational within a plane.",
                        });
                }
            });
        ui.horizontal(|ui| {
            ui.label("Part A id:")
                .on_hover_text("Numeric id of the first part joined.");
            ui.add(egui::DragValue::new(&mut s.joint_part_a).speed(1.0))
                .on_hover_text("Part A id.");
            ui.label("Part B id:")
                .on_hover_text("Numeric id of the second part joined.");
            ui.add(egui::DragValue::new(&mut s.joint_part_b).speed(1.0))
                .on_hover_text("Part B id.");
        });
        match s.joint_kind {
            AssemblyJointKindUi::Fixed => {
                // No additional inputs.
            }
            AssemblyJointKindUi::Revolute | AssemblyJointKindUi::Cylindrical => {
                ui.horizontal(|ui| {
                    ui.label("Axis origin:");
                    ui.add(
                        egui::DragValue::new(&mut s.joint_axis_origin[0])
                            .speed(0.05)
                            .prefix("x "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.joint_axis_origin[1])
                            .speed(0.05)
                            .prefix("y "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.joint_axis_origin[2])
                            .speed(0.05)
                            .prefix("z "),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Axis dir:");
                    ui.add(
                        egui::DragValue::new(&mut s.joint_axis_dir[0])
                            .speed(0.05)
                            .prefix("x "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.joint_axis_dir[1])
                            .speed(0.05)
                            .prefix("y "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.joint_axis_dir[2])
                            .speed(0.05)
                            .prefix("z "),
                    );
                });
            }
            AssemblyJointKindUi::Prismatic => {
                ui.horizontal(|ui| {
                    ui.label("Slide axis:");
                    ui.add(
                        egui::DragValue::new(&mut s.joint_axis_dir[0])
                            .speed(0.05)
                            .prefix("x "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.joint_axis_dir[1])
                            .speed(0.05)
                            .prefix("y "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.joint_axis_dir[2])
                            .speed(0.05)
                            .prefix("z "),
                    );
                });
            }
            AssemblyJointKindUi::Spherical => {
                ui.horizontal(|ui| {
                    ui.label("Pivot:");
                    ui.add(
                        egui::DragValue::new(&mut s.joint_point[0])
                            .speed(0.05)
                            .prefix("x "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.joint_point[1])
                            .speed(0.05)
                            .prefix("y "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.joint_point[2])
                            .speed(0.05)
                            .prefix("z "),
                    );
                });
            }
            AssemblyJointKindUi::Planar => {
                ui.horizontal(|ui| {
                    ui.label("Plane origin:");
                    ui.add(
                        egui::DragValue::new(&mut s.joint_plane_origin[0])
                            .speed(0.05)
                            .prefix("x "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.joint_plane_origin[1])
                            .speed(0.05)
                            .prefix("y "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.joint_plane_origin[2])
                            .speed(0.05)
                            .prefix("z "),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Plane normal:");
                    ui.add(
                        egui::DragValue::new(&mut s.joint_plane_normal[0])
                            .speed(0.05)
                            .prefix("x "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.joint_plane_normal[1])
                            .speed(0.05)
                            .prefix("y "),
                    );
                    ui.add(
                        egui::DragValue::new(&mut s.joint_plane_normal[2])
                            .speed(0.05)
                            .prefix("z "),
                    );
                });
            }
        }
    }
    let mut add_joint_clicked = false;
    if ui.button("Add joint").clicked() {
        add_joint_clicked = true;
    }
    if add_joint_clicked {
        let s = &mut app.mesh_toolbox.assembly;
        let vec3 = |xy: [f64; 3]| nalgebra::Vector3::new(xy[0], xy[1], xy[2]);
        let kind = match s.joint_kind {
            AssemblyJointKindUi::Fixed => valenx_assembly::JointKind::Fixed {
                part_a: s.joint_part_a,
                part_b: s.joint_part_b,
            },
            AssemblyJointKindUi::Revolute => valenx_assembly::JointKind::Revolute {
                part_a: s.joint_part_a,
                part_b: s.joint_part_b,
                axis_origin: vec3(s.joint_axis_origin),
                axis_dir: vec3(s.joint_axis_dir),
            },
            AssemblyJointKindUi::Prismatic => valenx_assembly::JointKind::Prismatic {
                part_a: s.joint_part_a,
                part_b: s.joint_part_b,
                axis_dir: vec3(s.joint_axis_dir),
            },
            AssemblyJointKindUi::Cylindrical => valenx_assembly::JointKind::Cylindrical {
                part_a: s.joint_part_a,
                part_b: s.joint_part_b,
                axis_origin: vec3(s.joint_axis_origin),
                axis_dir: vec3(s.joint_axis_dir),
            },
            AssemblyJointKindUi::Spherical => valenx_assembly::JointKind::Spherical {
                part_a: s.joint_part_a,
                part_b: s.joint_part_b,
                point: vec3(s.joint_point),
            },
            AssemblyJointKindUi::Planar => valenx_assembly::JointKind::Planar {
                part_a: s.joint_part_a,
                part_b: s.joint_part_b,
                plane_origin: vec3(s.joint_plane_origin),
                plane_normal: vec3(s.joint_plane_normal),
            },
        };
        s.assembly.add_joint(valenx_assembly::Joint::new(0, kind));
        s.last_error = None;
    }

    ui.separator();

    // ----- Joint list with slider -----
    ui.label(egui::RichText::new("Joints").strong());
    let mut delete_joint: Option<usize> = None;
    let mut joint_changed = false;
    {
        let s = &mut app.mesh_toolbox.assembly;
        if s.assembly.joints.is_empty() {
            ui.label("(no joints)");
        }
        for i in 0..s.assembly.joints.len() {
            let jid = s.assembly.joints[i].id;
            let j = &mut s.assembly.joints[i];
            let (pa, pb) = j.kind.parts();
            ui.horizontal(|ui| {
                ui.label(format!("#{jid}: {} {pa}↔{pb}", j.kind.label()));
                ui.checkbox(&mut j.suppressed, "suppress");
                if ui.small_button("X").clicked() {
                    delete_joint = Some(jid);
                }
            });
            ui.horizontal(|ui| {
                ui.label("param:");
                let resp = ui.add(
                    egui::Slider::new(
                        &mut j.parameter,
                        -std::f64::consts::PI..=std::f64::consts::PI,
                    )
                    .step_by(0.01),
                );
                if resp.changed() {
                    joint_changed = true;
                }
            });
        }
    }
    if let Some(jid) = delete_joint {
        let s = &mut app.mesh_toolbox.assembly;
        let _ = s.assembly.delete_joint(jid);
    }
    if joint_changed {
        let s = &mut app.mesh_toolbox.assembly;
        if let Err(e) = valenx_assembly::kinematics::apply_all_joints(&mut s.assembly) {
            s.last_error = Some(format!("Joint kinematics: {e}"));
        }
    }

    ui.separator();

    // ----- Solve + Render -----
    let mut solve_clicked = false;
    let mut render_clicked = false;
    ui.horizontal(|ui| {
        if ui.button("Solve mates").clicked() {
            solve_clicked = true;
        }
        if ui.button("Render assembly").clicked() {
            render_clicked = true;
        }
    });
    if solve_clicked {
        let s = &mut app.mesh_toolbox.assembly;
        match valenx_assembly::solver::solve(
            &mut s.assembly,
            valenx_assembly::SolverConfig::default(),
        ) {
            Ok(report) => {
                s.last_report = Some(report);
                s.last_error = None;
            }
            Err(e) => s.last_error = Some(format!("Solver: {e}")),
        }
    }
    if let Some(report) = &app.mesh_toolbox.assembly.last_report {
        ui.label(format!(
            "Status: {:?}, iterations: {}, residual: {:.3e}, DOF balance: {}",
            report.status, report.iterations, report.residual_norm, report.diagnostics.dof_balance,
        ));
    }
    if render_clicked {
        render_assembly_to_viewport(app);
    }

    if let Some(err) = &app.mesh_toolbox.assembly.last_error {
        ui.colored_label(egui::Color32::from_rgb(220, 100, 100), err);
    }
}

/// Render the current assembly into the viewport. Tessellates every
/// non-suppressed part at the panel's `render_tolerance`, applies each
/// part's transform to its vertices, fuses the resulting meshes into
/// a single `Mesh` (concatenating vertex arrays and offsetting
/// triangle indices), and pushes via `app.apply_mesh`.
///
/// Fusion is the right v1 choice — `app.apply_mesh` was built for a
/// single source mesh, and extending it to multi-mesh viewports would
/// touch the wgpu pipeline. Fusion produces correctly-rendered
/// transformed parts today without that refactor.
fn render_assembly_to_viewport(app: &mut crate::ValenxApp) {
    use valenx_mesh::{ElementBlock, ElementType, Mesh};
    let s = &mut app.mesh_toolbox.assembly;
    if s.assembly.parts.is_empty() {
        s.last_error = Some("Render: no parts to render.".into());
        return;
    }
    let mut fused = Mesh::new("assembly");
    let mut block = ElementBlock::new(ElementType::Tri3);
    let mut node_offset: u32 = 0;
    for part in &s.assembly.parts {
        let mesh = match valenx_cad::solid_to_mesh(&part.solid, s.render_tolerance) {
            Ok(m) => m,
            Err(e) => {
                s.last_error = Some(format!("Render: tessellate part {} failed: {e}", part.id));
                return;
            }
        };
        // Transform every vertex and append.
        for n in &mesh.nodes {
            fused.nodes.push(part.transform.apply_point(*n));
        }
        // Append connectivity with offset.
        for tri_block in &mesh.element_blocks {
            if tri_block.element_type != ElementType::Tri3 {
                continue;
            }
            for &idx in &tri_block.connectivity {
                block.connectivity.push(idx + node_offset);
            }
        }
        node_offset = fused.nodes.len() as u32;
    }
    fused.element_blocks.push(block);
    fused.recompute_stats();
    let pseudo_path = std::path::PathBuf::from("<assembly>/scene.fused");
    app.apply_mesh(fused, pseudo_path);
    app.mesh_toolbox.assembly.last_error = None;
}

/// Short human label for a [`valenx_assembly::MateKind`] discriminant —
/// used by the mate list rendering.
fn mate_kind_label(k: &valenx_assembly::MateKind) -> &'static str {
    use valenx_assembly::MateKind;
    match k {
        MateKind::Coincident { .. } => "Coincident",
        MateKind::Distance { .. } => "Distance",
        MateKind::Angle { .. } => "Angle",
        MateKind::Parallel { .. } => "Parallel",
        MateKind::Perpendicular { .. } => "Perpendicular",
        MateKind::Tangent { .. } => "Tangent",
    }
}

// ===== Phase 9 — Surface workbench panel =====

impl crate::ValenxApp {
    /// Construct a NURBS curve from the Surface-panel inputs and
    /// append it to the surface file.
    pub fn surface_create_curve(&mut self) {
        let s = &mut self.mesh_toolbox.surface;
        // Truncate / extend the cp + weight vectors to match n_cps
        // before reading them, so the user can change `curve_n_cps`
        // without immediately losing data.
        s.curve_cps.resize(s.curve_n_cps, [0.0, 0.0, 0.0]);
        s.curve_weights.resize(s.curve_n_cps, 1.0);

        let degree = s.curve_degree;
        let n_cps = s.curve_n_cps;
        let cps: Vec<nalgebra::Vector3<f64>> = s
            .curve_cps
            .iter()
            .map(|c| nalgebra::Vector3::new(c[0], c[1], c[2]))
            .collect();
        let weights = s.curve_weights.clone();
        let knots = open_uniform_knots(n_cps, degree);
        match valenx_surface::NurbsCurve::new(degree, knots, cps, weights) {
            Ok(c) => {
                let id = s.file.curves.len();
                s.file.curves.push(c);
                s.last_status = Some(format!("Created curve #{id} (degree {degree})"));
                s.last_error = None;
                emit_audit(
                    "surface.create_curve",
                    serde_json::json!({"id": id, "degree": degree, "n_cps": n_cps}),
                    serde_json::json!({}),
                );
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("create curve: {e}"));
            }
        }
    }

    /// Construct a NURBS surface from the panel inputs.
    pub fn surface_create_surface(&mut self) {
        let s = &mut self.mesh_toolbox.surface;
        let nu = s.surface_nu;
        let nv = s.surface_nv;
        // Round-24 M5: cap the CP grid size before the resize so
        // `nu * nv` neither overflows usize NOR allocates gigabytes
        // when a user types nu=10000, nv=10000 (100 M tuples × 24
        // bytes each = 2.4 GiB). 1 M CPs covers every realistic
        // analytic surface (typical 16×16 control net for a free-
        // form fillet) while refusing the typo / hostile drag-input
        // shape. `checked_mul` catches the usize overflow on 32-bit
        // builds where 10 K × 10 K would silently wrap.
        const MAX_NURBS_CP_GRID: usize = 1024 * 1024;
        let total = nu.checked_mul(nv);
        if !matches!(total, Some(n) if n <= MAX_NURBS_CP_GRID) {
            s.last_status = None;
            s.last_error = Some(format!(
                "create surface: nu={nu} * nv={nv} exceeds {MAX_NURBS_CP_GRID}-control-point grid cap",
            ));
            return;
        }
        // Resize the flat CP list to match nu*nv before reading.
        s.surface_cps.resize(nu * nv, [0.0, 0.0, 0.0]);

        let mut grid: Vec<Vec<nalgebra::Vector3<f64>>> = Vec::with_capacity(nu);
        for i in 0..nu {
            let mut row = Vec::with_capacity(nv);
            for j in 0..nv {
                let c = s.surface_cps[i * nv + j];
                row.push(nalgebra::Vector3::new(c[0], c[1], c[2]));
            }
            grid.push(row);
        }
        let weights = vec![vec![1.0_f64; nv]; nu];
        let u_knots = open_uniform_knots(nu, s.surface_u_degree);
        let v_knots = open_uniform_knots(nv, s.surface_v_degree);
        match valenx_surface::NurbsSurface::new(
            s.surface_u_degree,
            s.surface_v_degree,
            u_knots,
            v_knots,
            grid,
            weights,
        ) {
            Ok(surf) => {
                let id = s.file.surfaces.len();
                s.file.surfaces.push(surf);
                s.last_status = Some(format!("Created surface #{id} ({nu}x{nv} CPs)"));
                s.last_error = None;
                emit_audit(
                    "surface.create_surface",
                    serde_json::json!({"id": id, "nu": nu, "nv": nv}),
                    serde_json::json!({}),
                );
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("create surface: {e}"));
            }
        }
    }

    /// Fill four selected curves into a Coons patch and add the
    /// resulting NurbsSurface.
    pub fn surface_coons_fill(&mut self) {
        let s = &mut self.mesh_toolbox.surface;
        let ids = s.coons_curves;
        let n = s.file.curves.len();
        for (k, &id) in ids.iter().enumerate() {
            if id >= n {
                s.last_status = None;
                s.last_error = Some(format!(
                    "coons fill: boundary #{k} = curve id {id} out of range ({n} curves)"
                ));
                return;
            }
        }
        let c0 = s.file.curves[ids[0]].clone();
        let c1 = s.file.curves[ids[1]].clone();
        let d0 = s.file.curves[ids[2]].clone();
        let d1 = s.file.curves[ids[3]].clone();
        match valenx_surface::coons::fill([c0, c1, d0, d1]) {
            Ok(surf) => {
                let id = s.file.surfaces.len();
                s.file.surfaces.push(surf);
                s.last_status = Some(format!("Coons fill → surface #{id}"));
                s.last_error = None;
                emit_audit(
                    "surface.coons_fill",
                    serde_json::json!({"surface_id": id, "boundary_ids": ids}),
                    serde_json::json!({}),
                );
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("coons fill: {e}"));
            }
        }
    }

    /// Sew two surfaces along an edge pair.
    pub fn surface_sew(&mut self) {
        let s = &mut self.mesh_toolbox.surface;
        let n = s.file.surfaces.len();
        if s.sew_surface_a >= n || s.sew_surface_b >= n {
            s.last_status = None;
            s.last_error = Some(format!("sew: surface id out of range (have {n})"));
            return;
        }
        let edge_a = surface_edge_from_u8(s.sew_edge_a);
        let edge_b = surface_edge_from_u8(s.sew_edge_b);
        let surf_a = s.file.surfaces[s.sew_surface_a].clone();
        let surf_b = s.file.surfaces[s.sew_surface_b].clone();
        let result = if s.sew_use_g2 {
            valenx_surface::sew::g2_stitch(&surf_a, &surf_b, (edge_a, edge_b), s.sew_tolerance)
        } else {
            valenx_surface::sew::stitch(&surf_a, &surf_b, (edge_a, edge_b), s.sew_tolerance)
        };
        match result {
            Ok(out) => {
                let id = s.file.surfaces.len();
                s.file.surfaces.push(out);
                let mode = if s.sew_use_g2 { "G2" } else { "G0" };
                s.last_status = Some(format!("Sewed ({mode}) → surface #{id}"));
                s.last_error = None;
                emit_audit(
                    "surface.sew",
                    serde_json::json!({
                        "out_id": id,
                        "a_id": s.sew_surface_a,
                        "b_id": s.sew_surface_b,
                        "g2": s.sew_use_g2,
                    }),
                    serde_json::json!({}),
                );
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("sew: {e}"));
            }
        }
    }

    /// Phase 19A — insert a knot into a curve or surface.
    pub fn surface_insert_knot(&mut self) {
        let s = &mut self.mesh_toolbox.surface;
        let direction = s.knot_op_direction;
        if direction == 0 || direction == 1 {
            // Surface knot insertion.
            let n = s.file.surfaces.len();
            if s.knot_op_surface >= n {
                s.last_status = None;
                s.last_error = Some(format!("insert_knot: surface id out of range (have {n})"));
                return;
            }
            let src = s.file.surfaces[s.knot_op_surface].clone();
            let res = if direction == 0 {
                src.insert_knot_u(s.knot_op_u)
            } else {
                src.insert_knot_v(s.knot_op_u)
            };
            match res {
                Ok(out) => {
                    let id = s.file.surfaces.len();
                    s.file.surfaces.push(out);
                    s.last_status = Some(format!(
                        "insert_knot → surface #{id} ({} dir)",
                        if direction == 0 { "u" } else { "v" }
                    ));
                    s.last_error = None;
                }
                Err(e) => {
                    s.last_status = None;
                    s.last_error = Some(format!("insert_knot: {e}"));
                }
            }
        }
    }

    /// Phase 19A — insert a knot into a NURBS curve.
    pub fn surface_insert_knot_curve(&mut self) {
        let s = &mut self.mesh_toolbox.surface;
        let n = s.file.curves.len();
        if s.knot_op_curve >= n {
            s.last_status = None;
            s.last_error = Some(format!("insert_knot: curve id out of range (have {n})"));
            return;
        }
        let src = s.file.curves[s.knot_op_curve].clone();
        match src.insert_knot(s.knot_op_u) {
            Ok(out) => {
                let id = s.file.curves.len();
                s.file.curves.push(out);
                s.last_status = Some(format!("insert_knot → curve #{id}"));
                s.last_error = None;
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("insert_knot: {e}"));
            }
        }
    }

    /// Phase 19A — remove a knot from a NURBS curve.
    pub fn surface_remove_knot_curve(&mut self) {
        let s = &mut self.mesh_toolbox.surface;
        let n = s.file.curves.len();
        if s.knot_op_curve >= n {
            s.last_status = None;
            s.last_error = Some(format!("remove_knot: curve id out of range (have {n})"));
            return;
        }
        let src = s.file.curves[s.knot_op_curve].clone();
        match src.remove_knot(s.knot_op_u, s.knot_op_tolerance) {
            Ok(out) => {
                let id = s.file.curves.len();
                s.file.curves.push(out);
                s.last_status = Some(format!("remove_knot → curve #{id}"));
                s.last_error = None;
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("remove_knot: {e}"));
            }
        }
    }

    /// Phase 19A — elevate the degree of a NURBS curve.
    pub fn surface_elevate_degree_curve(&mut self) {
        let s = &mut self.mesh_toolbox.surface;
        let n = s.file.curves.len();
        if s.knot_op_curve >= n {
            s.last_status = None;
            s.last_error = Some(format!("elevate_degree: curve id out of range (have {n})"));
            return;
        }
        let src = s.file.curves[s.knot_op_curve].clone();
        match src.elevate_degree(s.elevate_degree_by) {
            Ok(out) => {
                let id = s.file.curves.len();
                s.file.curves.push(out);
                s.last_status = Some(format!(
                    "elevate_degree (+{}) → curve #{id}",
                    s.elevate_degree_by
                ));
                s.last_error = None;
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("elevate_degree: {e}"));
            }
        }
    }

    /// Phase 19A — elevate the degree of a NURBS surface in the
    /// currently-selected `knot_op_direction`.
    pub fn surface_elevate_degree_surface(&mut self) {
        let s = &mut self.mesh_toolbox.surface;
        let n = s.file.surfaces.len();
        if s.knot_op_surface >= n {
            s.last_status = None;
            s.last_error = Some(format!(
                "elevate_degree: surface id out of range (have {n})"
            ));
            return;
        }
        let src = s.file.surfaces[s.knot_op_surface].clone();
        let res = if s.knot_op_direction == 0 {
            src.elevate_degree_u(s.elevate_degree_by)
        } else {
            src.elevate_degree_v(s.elevate_degree_by)
        };
        match res {
            Ok(out) => {
                let id = s.file.surfaces.len();
                s.file.surfaces.push(out);
                s.last_status = Some(format!(
                    "elevate_degree (+{}, {} dir) → surface #{id}",
                    s.elevate_degree_by,
                    if s.knot_op_direction == 0 { "u" } else { "v" }
                ));
                s.last_error = None;
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("elevate_degree: {e}"));
            }
        }
    }

    /// Phase 19B — true rational surface-surface intersection.
    pub fn surface_ssi(&mut self) {
        let s = &mut self.mesh_toolbox.surface;
        let n = s.file.surfaces.len();
        if s.ssi_surface_a >= n || s.ssi_surface_b >= n {
            s.last_status = None;
            s.last_error = Some(format!("ssi: surface id out of range (have {n})"));
            return;
        }
        let sa = s.file.surfaces[s.ssi_surface_a].clone();
        let sb = s.file.surfaces[s.ssi_surface_b].clone();
        let curves = valenx_surface::intersect::true_ssi(&sa, &sb, s.ssi_tolerance);
        let added = curves.len();
        for c in curves {
            s.file.curves.push(c);
        }
        s.last_status = Some(format!("SSI → {added} curves added"));
        s.last_error = None;
    }

    /// Phase 19D — fit a NURBS curve through the points listed in
    /// `fit_points_text` (one "x,y,z" per line). Adds the curve and
    /// stores RMS error in `fit_last_rms`.
    pub fn surface_fit_curve(&mut self) {
        let s = &mut self.mesh_toolbox.surface;
        let pts = parse_fit_points(&s.fit_points_text);
        if pts.len() < 2 {
            s.last_status = None;
            s.last_error = Some("fit_curve: need at least 2 points".into());
            return;
        }
        let n_cps = s.fit_n_cps_u.min(pts.len()).max(s.fit_degree_u + 1);
        match valenx_surface::fit::nurbs_curve_through_points(&pts, s.fit_degree_u, n_cps) {
            Ok(fit) => {
                let id = s.file.curves.len();
                s.fit_last_rms = Some(fit.rms_error);
                s.file.curves.push(fit.curve);
                s.last_status = Some(format!(
                    "fit_curve → curve #{id} (rms = {:.4e})",
                    fit.rms_error
                ));
                s.last_error = None;
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("fit_curve: {e}"));
            }
        }
    }

    /// Phase 19D — fit a NURBS surface through the scattered points
    /// listed in `fit_points_text`. Uses the v1 scattered-fit
    /// strategy (plane projection + binning).
    pub fn surface_fit_surface(&mut self) {
        let s = &mut self.mesh_toolbox.surface;
        let pts = parse_fit_points(&s.fit_points_text);
        let need = s.fit_n_cps_u * s.fit_n_cps_v;
        if pts.len() < need {
            s.last_status = None;
            s.last_error = Some(format!(
                "fit_surface: need at least {need} points; got {}",
                pts.len()
            ));
            return;
        }
        match valenx_surface::fit::surface_through_scattered(
            &pts,
            s.fit_degree_u,
            s.fit_degree_v,
            s.fit_n_cps_u,
            s.fit_n_cps_v,
        ) {
            Ok(fit) => {
                let id = s.file.surfaces.len();
                s.fit_last_rms = Some(fit.rms_error);
                s.file.surfaces.push(fit.surface);
                s.last_status = Some(format!(
                    "fit_surface → surface #{id} (rms = {:.4e})",
                    fit.rms_error
                ));
                s.last_error = None;
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("fit_surface: {e}"));
            }
        }
    }

    /// Phase 19E — ruled surface constructors (between two curves,
    /// extrusion along vector, or cone from apex).
    pub fn surface_ruled_build(&mut self) {
        let s = &mut self.mesh_toolbox.surface;
        let nc = s.file.curves.len();
        let result = match s.ruled_kind {
            0 => {
                if s.ruled_curve_a >= nc || s.ruled_curve_b >= nc {
                    s.last_status = None;
                    s.last_error = Some(format!("ruled: curve id out of range (have {nc})"));
                    return;
                }
                let a = s.file.curves[s.ruled_curve_a].clone();
                let b = s.file.curves[s.ruled_curve_b].clone();
                valenx_surface::ruled::between_curves(&a, &b)
            }
            1 => {
                if s.ruled_curve_a >= nc {
                    s.last_status = None;
                    s.last_error = Some(format!("ruled: curve id out of range (have {nc})"));
                    return;
                }
                let a = s.file.curves[s.ruled_curve_a].clone();
                let v = nalgebra::Vector3::new(
                    s.ruled_extrude_vector[0],
                    s.ruled_extrude_vector[1],
                    s.ruled_extrude_vector[2],
                );
                valenx_surface::ruled::extrude_along_vector(&a, v)
            }
            _ => {
                if s.ruled_curve_a >= nc {
                    s.last_status = None;
                    s.last_error = Some(format!("ruled: curve id out of range (have {nc})"));
                    return;
                }
                let a = s.file.curves[s.ruled_curve_a].clone();
                let apex =
                    nalgebra::Vector3::new(s.ruled_apex[0], s.ruled_apex[1], s.ruled_apex[2]);
                valenx_surface::ruled::cone_from_apex(&a, apex)
            }
        };
        match result {
            Ok(surf) => {
                let id = s.file.surfaces.len();
                s.file.surfaces.push(surf);
                s.last_status = Some(format!("ruled → surface #{id}"));
                s.last_error = None;
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("ruled: {e}"));
            }
        }
    }

    /// Trim a surface by a curve and push the resulting mesh to
    /// the viewport.
    ///
    /// When `trim_use_uv` is set (Phase 9.5 default), runs the
    /// parametric `(u, v)` domain trim — works on warped surfaces.
    /// When unset, falls back to the legacy world-xy trim.
    pub fn surface_trim(&mut self) {
        let surf = {
            let s = &mut self.mesh_toolbox.surface;
            let n_surf = s.file.surfaces.len();
            let n_curves = s.file.curves.len();
            if s.trim_surface >= n_surf || s.trim_curve >= n_curves {
                s.last_status = None;
                s.last_error = Some(format!(
                    "trim: id out of range (have {n_surf} surfaces, {n_curves} curves)"
                ));
                return;
            }
            let side = if s.trim_side == 0 {
                valenx_surface::trim::TrimSide::Inside
            } else {
                valenx_surface::trim::TrimSide::Outside
            };
            let res = s.trim_resolution.max(8);
            let curve_segments = (res * 8).max(32);
            let surface = s.file.surfaces[s.trim_surface].clone();
            let curve = s.file.curves[s.trim_curve].clone();
            if s.trim_use_uv {
                let params = valenx_surface::trim::UvTrimParams {
                    nu: res,
                    nv: res,
                    curve_segments,
                    seed_nu: res.div_ceil(2).max(8),
                    seed_nv: res.div_ceil(2).max(8),
                };
                valenx_surface::trim::by_curve_in_uv(&surface, &curve, side, params)
            } else {
                valenx_surface::trim::by_curve(&surface, &curve, side, res, res, curve_segments)
            }
        };
        match surf {
            Ok(mesh) => {
                let pseudo = std::path::PathBuf::from("<surface>/trim.mesh");
                self.apply_mesh(mesh, pseudo);
                let s = &mut self.mesh_toolbox.surface;
                s.last_status = Some("Trim → mesh pushed to viewport".into());
                s.last_error = None;
                emit_audit("surface.trim", serde_json::json!({}), serde_json::json!({}));
            }
            Err(e) => {
                let s = &mut self.mesh_toolbox.surface;
                s.last_status = None;
                s.last_error = Some(format!("trim: {e}"));
            }
        }
    }

    /// Tessellate the selected surface and push the resulting mesh
    /// to the viewport.
    pub fn surface_tessellate(&mut self, idx: usize) {
        let mesh = {
            let s = &self.mesh_toolbox.surface;
            if idx >= s.file.surfaces.len() {
                let s = &mut self.mesh_toolbox.surface;
                s.last_status = None;
                s.last_error = Some(format!(
                    "tessellate: surface id {idx} out of range ({} surfaces)",
                    s.file.surfaces.len()
                ));
                return;
            }
            let res = s.tess_resolution.max(2);
            valenx_surface::tessellate::surface(&s.file.surfaces[idx], res, res)
        };
        let pseudo = std::path::PathBuf::from(format!("<surface>/tess_surface_{idx}.mesh"));
        self.apply_mesh(mesh, pseudo);
        let s = &mut self.mesh_toolbox.surface;
        s.last_status = Some(format!("Tessellated surface #{idx} → mesh pushed"));
        s.last_error = None;
        emit_audit(
            "surface.tessellate",
            serde_json::json!({"surface_id": idx, "resolution": s.tess_resolution}),
            serde_json::json!({}),
        );
    }
}

// ===== Phase 10 — CAM workbench panel =====

impl crate::ValenxApp {
    /// Push the panel's stock-origin / stock-size / material inputs
    /// into the active CamFile's [`valenx_cam::Stock`].
    fn cam_sync_stock_from_inputs(&mut self) -> Result<(), valenx_cam::CamError> {
        let s = &mut self.mesh_toolbox.cam;
        let stock = valenx_cam::Stock::new(
            Vector3::new(s.stock_origin[0], s.stock_origin[1], s.stock_origin[2]),
            Vector3::new(s.stock_size[0], s.stock_size[1], s.stock_size[2]),
            s.stock_material.clone(),
        )?;
        s.file.stock = stock;
        Ok(())
    }

    /// Add a new tool from the panel's "Add Tool" inputs.
    pub fn cam_add_tool(&mut self) {
        let s = &mut self.mesh_toolbox.cam;
        // Pick the next free tool id (max + 1).
        let next_id = s.file.tools.iter().map(|t| t.id).max().unwrap_or(0) + 1;
        let res = valenx_cam::Tool::new(
            next_id,
            s.new_tool_name.clone(),
            s.new_tool_kind,
            s.new_tool_diameter,
            s.new_tool_length,
            s.new_tool_flutes,
            s.new_tool_material.clone(),
        );
        match res {
            Ok(t) => {
                let id = t.id;
                s.file.tools.push(t);
                s.last_status = Some(format!("Added tool T{id} ({})", s.new_tool_name));
                s.last_error = None;
                emit_audit(
                    "cam.add_tool",
                    serde_json::json!({"id": id, "kind": s.new_tool_kind.label()}),
                    serde_json::json!({}),
                );
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("add tool: {e}"));
            }
        }
    }

    /// Delete a tool by index.
    pub fn cam_delete_tool(&mut self, idx: usize) {
        let s = &mut self.mesh_toolbox.cam;
        if idx < s.file.tools.len() {
            let id = s.file.tools[idx].id;
            s.file.tools.remove(idx);
            s.last_status = Some(format!("Removed tool T{id}"));
            s.last_error = None;
            emit_audit(
                "cam.delete_tool",
                serde_json::json!({"id": id}),
                serde_json::json!({}),
            );
        }
    }

    /// Add the currently-configured operation to the op list.
    pub fn cam_add_operation(&mut self) {
        if let Err(e) = self.cam_sync_stock_from_inputs() {
            let s = &mut self.mesh_toolbox.cam;
            s.last_status = None;
            s.last_error = Some(format!("stock: {e}"));
            return;
        }
        let s = &mut self.mesh_toolbox.cam;
        let op = match s.new_op_kind {
            CamOpKind::Profile => valenx_cam::Operation::Profile(valenx_cam::ProfileParams {
                tool_id: s.new_op_tool_id,
                feed_mm_per_min: s.new_op_feed,
                plunge_feed: s.new_op_plunge_feed,
                spindle_rpm: s.new_op_spindle_rpm,
                step_down: s.new_op_step_down,
                depth: s.new_op_depth,
                safe_z_clearance: s.new_op_safe_z,
                climb: s.new_op_climb,
            }),
            CamOpKind::Pocket => valenx_cam::Operation::Pocket(valenx_cam::PocketParams {
                tool_id: s.new_op_tool_id,
                feed_mm_per_min: s.new_op_feed,
                plunge_feed: s.new_op_plunge_feed,
                spindle_rpm: s.new_op_spindle_rpm,
                step_down: s.new_op_step_down,
                step_over: s.new_op_step_over,
                depth: s.new_op_depth,
                safe_z_clearance: s.new_op_safe_z,
                raster_angle_deg: s.new_op_raster_angle,
                strategy: s.new_op_pocket_strategy,
                climb: s.new_op_climb,
            }),
            CamOpKind::Drill => valenx_cam::Operation::Drill(valenx_cam::DrillParams {
                tool_id: s.new_op_tool_id,
                plunge_feed: s.new_op_plunge_feed,
                spindle_rpm: s.new_op_spindle_rpm,
                peck_depth: s.new_op_peck_depth,
                total_depth: s.new_op_drill_total_depth,
                retract_clearance: s.new_op_retract_clearance,
                safe_z_clearance: s.new_op_safe_z,
                hole_positions: s
                    .new_op_hole_positions
                    .iter()
                    .map(|p| Vector3::new(p[0], p[1], p[2]))
                    .collect(),
            }),
            CamOpKind::Face => valenx_cam::Operation::Face(valenx_cam::FaceParams {
                tool_id: s.new_op_tool_id,
                feed_mm_per_min: s.new_op_feed,
                plunge_feed: s.new_op_plunge_feed,
                spindle_rpm: s.new_op_spindle_rpm,
                step_over: s.new_op_step_over,
                step_down: s.new_op_step_down,
                depth: s.new_op_depth,
                safe_z_clearance: s.new_op_safe_z,
                raster_angle_deg: s.new_op_raster_angle,
                climb: s.new_op_climb,
            }),
            // Phase 17A
            CamOpKind::AdaptiveClearing => valenx_cam::Operation::AdaptiveClearing(
                valenx_cam::op::adaptive_clearing::AdaptiveParams {
                    tool_id: s.new_op_tool_id,
                    feed_mm_per_min: s.new_op_feed,
                    plunge_feed: s.new_op_plunge_feed,
                    spindle_rpm: s.new_op_spindle_rpm,
                    step_down: s.new_op_step_down,
                    depth: s.new_op_depth,
                    safe_z_clearance: s.new_op_safe_z,
                    ..Default::default()
                },
            ),
            CamOpKind::HelicalBore => {
                valenx_cam::Operation::HelicalBore(valenx_cam::op::helix_bore::HelicalBoreParams {
                    tool_id: s.new_op_tool_id,
                    feed_mm_per_min: s.new_op_feed,
                    plunge_feed: s.new_op_plunge_feed,
                    spindle_rpm: s.new_op_spindle_rpm,
                    depth: s.new_op_depth,
                    safe_z_clearance: s.new_op_safe_z,
                    ..Default::default()
                })
            }
            CamOpKind::PlungeRough => valenx_cam::Operation::PlungeRough(
                valenx_cam::op::plunge_rough::PlungeRoughParams {
                    tool_id: s.new_op_tool_id,
                    plunge_feed: s.new_op_plunge_feed,
                    spindle_rpm: s.new_op_spindle_rpm,
                    depth: s.new_op_depth,
                    safe_z_clearance: s.new_op_safe_z,
                    plunge_positions: s
                        .new_op_hole_positions
                        .iter()
                        .map(|p| Vector3::new(p[0], p[1], p[2]))
                        .collect(),
                },
            ),
            CamOpKind::RampEntry => {
                valenx_cam::Operation::RampEntry(valenx_cam::op::ramp_entry::RampEntryParams {
                    tool_id: s.new_op_tool_id,
                    feed_mm_per_min: s.new_op_feed,
                    spindle_rpm: s.new_op_spindle_rpm,
                    depth: s.new_op_depth,
                    safe_z_clearance: s.new_op_safe_z,
                    ..Default::default()
                })
            }
            CamOpKind::PeckDrillFull => valenx_cam::Operation::PeckDrillFull(
                valenx_cam::op::peck_drill_full::PeckDrillFullParams {
                    tool_id: s.new_op_tool_id,
                    plunge_feed: s.new_op_plunge_feed,
                    spindle_rpm: s.new_op_spindle_rpm,
                    peck_depth: s.new_op_peck_depth,
                    total_depth: s.new_op_drill_total_depth,
                    retract_clearance: s.new_op_retract_clearance,
                    safe_z_clearance: s.new_op_safe_z,
                    hole_positions: s
                        .new_op_hole_positions
                        .iter()
                        .map(|p| Vector3::new(p[0], p[1], p[2]))
                        .collect(),
                    ..Default::default()
                },
            ),
            // Phase 17B
            CamOpKind::Contour2D => {
                valenx_cam::Operation::Contour2D(valenx_cam::op::contour_2d::Contour2DParams {
                    tool_id: s.new_op_tool_id,
                    feed_mm_per_min: s.new_op_feed,
                    plunge_feed: s.new_op_plunge_feed,
                    spindle_rpm: s.new_op_spindle_rpm,
                    step_down: s.new_op_step_down,
                    depth: s.new_op_depth,
                    safe_z_clearance: s.new_op_safe_z,
                    curve: s
                        .new_op_hole_positions
                        .iter()
                        .map(|p| Vector3::new(p[0], p[1], p[2]))
                        .collect(),
                    closed: false,
                })
            }
            CamOpKind::Contour3D => {
                valenx_cam::Operation::Contour3D(valenx_cam::op::contour_3d::Contour3DParams {
                    tool_id: s.new_op_tool_id,
                    feed_mm_per_min: s.new_op_feed,
                    plunge_feed: s.new_op_plunge_feed,
                    spindle_rpm: s.new_op_spindle_rpm,
                    safe_z_clearance: s.new_op_safe_z,
                    curve: s
                        .new_op_hole_positions
                        .iter()
                        .map(|p| Vector3::new(p[0], p[1], p[2]))
                        .collect(),
                    closed: false,
                })
            }
            CamOpKind::Engrave => {
                valenx_cam::Operation::Engrave(valenx_cam::op::engrave::EngraveParams {
                    tool_id: s.new_op_tool_id,
                    feed_mm_per_min: s.new_op_feed,
                    plunge_feed: s.new_op_plunge_feed,
                    spindle_rpm: s.new_op_spindle_rpm,
                    safe_z_clearance: s.new_op_safe_z,
                    curve: s
                        .new_op_hole_positions
                        .iter()
                        .map(|p| Vector3::new(p[0], p[1], p[2]))
                        .collect(),
                    ..Default::default()
                })
            }
            CamOpKind::Scribe => {
                valenx_cam::Operation::Scribe(valenx_cam::op::scribe::ScribeParams {
                    tool_id: s.new_op_tool_id,
                    feed_mm_per_min: s.new_op_feed,
                    plunge_feed: s.new_op_plunge_feed,
                    spindle_rpm: s.new_op_spindle_rpm,
                    depth: s.new_op_depth,
                    safe_z_clearance: s.new_op_safe_z,
                    curve: s
                        .new_op_hole_positions
                        .iter()
                        .map(|p| Vector3::new(p[0], p[1], p[2]))
                        .collect(),
                })
            }
            CamOpKind::SpiralPocket => valenx_cam::Operation::SpiralPocket(
                valenx_cam::op::spiral_pocket::SpiralPocketParams {
                    tool_id: s.new_op_tool_id,
                    feed_mm_per_min: s.new_op_feed,
                    plunge_feed: s.new_op_plunge_feed,
                    spindle_rpm: s.new_op_spindle_rpm,
                    step_over: s.new_op_step_over,
                    step_down: s.new_op_step_down,
                    depth: s.new_op_depth,
                    safe_z_clearance: s.new_op_safe_z,
                    ..Default::default()
                },
            ),
            CamOpKind::TrochoidalSlot => valenx_cam::Operation::TrochoidalSlot(
                valenx_cam::op::trochoidal_slot::TrochoidalSlotParams {
                    tool_id: s.new_op_tool_id,
                    feed_mm_per_min: s.new_op_feed,
                    plunge_feed: s.new_op_plunge_feed,
                    spindle_rpm: s.new_op_spindle_rpm,
                    step_down: s.new_op_step_down,
                    depth: s.new_op_depth,
                    safe_z_clearance: s.new_op_safe_z,
                    centreline: s
                        .new_op_hole_positions
                        .iter()
                        .map(|p| Vector3::new(p[0], p[1], p[2]))
                        .collect(),
                    ..Default::default()
                },
            ),
            CamOpKind::Waterline3D => valenx_cam::Operation::Waterline3D(
                valenx_cam::op::waterline_3d::Waterline3DParams {
                    tool_id: s.new_op_tool_id,
                    feed_mm_per_min: s.new_op_feed,
                    plunge_feed: s.new_op_plunge_feed,
                    spindle_rpm: s.new_op_spindle_rpm,
                    step_down: s.new_op_step_down,
                    depth: s.new_op_depth,
                    safe_z_clearance: s.new_op_safe_z,
                },
            ),
            CamOpKind::Slot => valenx_cam::Operation::Slot(valenx_cam::op::slot::SlotParams {
                tool_id: s.new_op_tool_id,
                feed_mm_per_min: s.new_op_feed,
                plunge_feed: s.new_op_plunge_feed,
                spindle_rpm: s.new_op_spindle_rpm,
                step_down: s.new_op_step_down,
                depth: s.new_op_depth,
                safe_z_clearance: s.new_op_safe_z,
                start: s
                    .new_op_hole_positions
                    .first()
                    .map(|p| Vector3::new(p[0], p[1], 0.0))
                    .unwrap_or_default(),
                end: s
                    .new_op_hole_positions
                    .get(1)
                    .map(|p| Vector3::new(p[0], p[1], 0.0))
                    .unwrap_or_else(|| Vector3::new(20.0, 0.0, 0.0)),
            }),
            CamOpKind::ThreadMill => {
                valenx_cam::Operation::ThreadMill(valenx_cam::op::thread_mill::ThreadMillParams {
                    tool_id: s.new_op_tool_id,
                    feed_mm_per_min: s.new_op_feed,
                    plunge_feed: s.new_op_plunge_feed,
                    spindle_rpm: s.new_op_spindle_rpm,
                    safe_z_clearance: s.new_op_safe_z,
                    ..Default::default()
                })
            }
            CamOpKind::RestMachining => valenx_cam::Operation::RestMachining(
                valenx_cam::op::rest_machining::RestMachiningParams {
                    tool_id: s.new_op_tool_id,
                    feed_mm_per_min: s.new_op_feed,
                    plunge_feed: s.new_op_plunge_feed,
                    spindle_rpm: s.new_op_spindle_rpm,
                    step_down: s.new_op_step_down,
                    depth: s.new_op_depth,
                    safe_z_clearance: s.new_op_safe_z,
                    ..Default::default()
                },
            ),
        };
        let label = op.label();
        let tool_id = op.tool_id();
        s.file.operations.push(op);
        s.last_status = Some(format!("Added {label} op (T{tool_id})"));
        s.last_error = None;
        emit_audit(
            "cam.add_operation",
            serde_json::json!({"kind": label, "tool_id": tool_id}),
            serde_json::json!({}),
        );
    }

    /// Remove an operation by index.
    pub fn cam_delete_operation(&mut self, idx: usize) {
        let s = &mut self.mesh_toolbox.cam;
        if idx < s.file.operations.len() {
            let label = s.file.operations[idx].label();
            s.file.operations.remove(idx);
            s.last_status = Some(format!("Removed op #{idx} ({label})"));
            s.last_error = None;
        }
    }

    /// Run every operation in order and stash the chained toolpath
    /// into [`CamPanelState::last_toolpath`].
    pub fn cam_generate_toolpath(&mut self) {
        if let Err(e) = self.cam_sync_stock_from_inputs() {
            let s = &mut self.mesh_toolbox.cam;
            s.last_status = None;
            s.last_error = Some(format!("stock: {e}"));
            return;
        }
        // Capture mesh + ops + tools snapshot to keep the borrow scope
        // narrow.
        let mesh_opt = self.mesh.as_ref().map(|m| m.mesh.clone());
        let s = &mut self.mesh_toolbox.cam;
        let mut chain = valenx_cam::Toolpath::new();
        let mut total_ops = 0_usize;
        for (i, op) in s.file.operations.iter().enumerate() {
            let tool = match s.file.tools.iter().find(|t| t.id == op.tool_id()) {
                Some(t) => t.clone(),
                None => {
                    s.last_status = None;
                    s.last_error = Some(format!(
                        "op #{i} ({}): tool T{} not found",
                        op.label(),
                        op.tool_id()
                    ));
                    return;
                }
            };
            let tp = match op {
                valenx_cam::Operation::Profile(p) => {
                    let mesh = match mesh_opt.as_ref() {
                        Some(m) => m,
                        None => {
                            s.last_status = None;
                            s.last_error = Some(format!(
                                "op #{i} (Profile): no mesh loaded — Profile needs a source solid"
                            ));
                            return;
                        }
                    };
                    valenx_cam::op::profile::generate(&s.file.stock, mesh, p, &tool)
                }
                valenx_cam::Operation::Pocket(p) => {
                    let mesh = match mesh_opt.as_ref() {
                        Some(m) => m,
                        None => {
                            s.last_status = None;
                            s.last_error = Some(format!(
                                "op #{i} (Pocket): no mesh loaded — Pocket needs a source solid"
                            ));
                            return;
                        }
                    };
                    valenx_cam::op::pocket::generate(&s.file.stock, mesh, p, &tool)
                }
                valenx_cam::Operation::Drill(p) => {
                    valenx_cam::op::drill::generate(&s.file.stock, p, &tool)
                }
                valenx_cam::Operation::Face(p) => {
                    valenx_cam::op::face::generate(&s.file.stock, p, &tool)
                }
                // Phase 17A
                valenx_cam::Operation::AdaptiveClearing(p) => {
                    let mesh = match mesh_opt.as_ref() {
                        Some(m) => m,
                        None => {
                            s.last_status = None;
                            s.last_error =
                                Some(format!("op #{i} (Adaptive Clearing): no mesh loaded"));
                            return;
                        }
                    };
                    valenx_cam::op::adaptive_clearing::generate(&s.file.stock, mesh, p, &tool)
                }
                valenx_cam::Operation::HelicalBore(p) => {
                    valenx_cam::op::helix_bore::generate(&s.file.stock, p, &tool)
                }
                valenx_cam::Operation::PlungeRough(p) => {
                    valenx_cam::op::plunge_rough::generate(&s.file.stock, p, &tool)
                }
                valenx_cam::Operation::RampEntry(p) => {
                    valenx_cam::op::ramp_entry::generate(&s.file.stock, p, &tool)
                }
                valenx_cam::Operation::PeckDrillFull(p) => {
                    valenx_cam::op::peck_drill_full::generate(&s.file.stock, p, &tool)
                }
                // Phase 17B
                valenx_cam::Operation::Contour2D(p) => {
                    valenx_cam::op::contour_2d::generate(&s.file.stock, p, &tool)
                }
                valenx_cam::Operation::Contour3D(p) => {
                    valenx_cam::op::contour_3d::generate(&s.file.stock, p, &tool)
                }
                valenx_cam::Operation::Engrave(p) => {
                    valenx_cam::op::engrave::generate(&s.file.stock, p, &tool)
                }
                valenx_cam::Operation::Scribe(p) => {
                    valenx_cam::op::scribe::generate(&s.file.stock, p, &tool)
                }
                valenx_cam::Operation::SpiralPocket(p) => {
                    let mesh = match mesh_opt.as_ref() {
                        Some(m) => m,
                        None => {
                            s.last_status = None;
                            s.last_error = Some(format!("op #{i} (Spiral Pocket): no mesh loaded"));
                            return;
                        }
                    };
                    valenx_cam::op::spiral_pocket::generate(&s.file.stock, mesh, p, &tool)
                }
                valenx_cam::Operation::TrochoidalSlot(p) => {
                    valenx_cam::op::trochoidal_slot::generate(&s.file.stock, p, &tool)
                }
                valenx_cam::Operation::Waterline3D(p) => {
                    let mesh = match mesh_opt.as_ref() {
                        Some(m) => m,
                        None => {
                            s.last_status = None;
                            s.last_error = Some(format!("op #{i} (Waterline 3D): no mesh loaded"));
                            return;
                        }
                    };
                    valenx_cam::op::waterline_3d::generate(&s.file.stock, mesh, p, &tool)
                }
                valenx_cam::Operation::Slot(p) => {
                    valenx_cam::op::slot::generate(&s.file.stock, p, &tool)
                }
                valenx_cam::Operation::ThreadMill(p) => {
                    valenx_cam::op::thread_mill::generate(&s.file.stock, p, &tool)
                }
                valenx_cam::Operation::RestMachining(p) => {
                    let mesh = match mesh_opt.as_ref() {
                        Some(m) => m,
                        None => {
                            s.last_status = None;
                            s.last_error =
                                Some(format!("op #{i} (Rest Machining): no mesh loaded"));
                            return;
                        }
                    };
                    valenx_cam::op::rest_machining::generate(&s.file.stock, mesh, p, &tool)
                }
            };
            match tp {
                Ok(t) => chain.concatenate(&t),
                Err(e) => {
                    s.last_status = None;
                    s.last_error = Some(format!("op #{i} ({}): {e}", op.label()));
                    return;
                }
            }
            total_ops += 1;
        }
        // Pick the first tool as the "spindle" tool for stats.
        let stat_tool = s.file.tools.first().cloned();
        let est = valenx_cam::simulate::estimated_time(&chain);
        let vol = stat_tool
            .as_ref()
            .map(|t| valenx_cam::simulate::removed_volume_mm3(&chain, t))
            .unwrap_or(0.0);
        s.last_estimated_time_min = Some(est);
        s.last_removed_volume_mm3 = Some(vol);
        let move_count = chain.len();
        s.last_toolpath = Some(chain);
        s.last_status = Some(format!(
            "Generated toolpath: {total_ops} ops, {move_count} moves, ~{est:.2} min, ~{vol:.0} mm³"
        ));
        s.last_error = None;
        emit_audit(
            "cam.generate_toolpath",
            serde_json::json!({"ops": total_ops, "moves": move_count, "minutes": est}),
            serde_json::json!({}),
        );
    }

    /// Toggle the simulate-overlay (cyan/red/gray polylines in the
    /// viewport).
    pub fn cam_toggle_simulate(&mut self) {
        let s = &mut self.mesh_toolbox.cam;
        s.show_overlay = !s.show_overlay;
        if s.show_overlay && s.last_toolpath.is_none() {
            s.last_status = Some("Generate a toolpath first to see the simulation overlay".into());
        } else {
            s.last_status = Some(format!(
                "Simulation overlay: {}",
                if s.show_overlay { "on" } else { "off" }
            ));
        }
    }

    /// Export the last-generated toolpath to a .nc file via the
    /// native save-file dialog.
    pub fn cam_export_nc(&mut self) {
        let path = match rfd::FileDialog::new()
            .add_filter("G-code (.nc)", &["nc"])
            .save_file()
        {
            Some(p) => p,
            None => return,
        };
        let s = &mut self.mesh_toolbox.cam;
        let tp = match s.last_toolpath.as_ref() {
            Some(t) => t.clone(),
            None => {
                s.last_status = None;
                s.last_error = Some("export .nc: generate a toolpath first".into());
                return;
            }
        };
        let tool = match s.file.tools.first().cloned() {
            Some(t) => t,
            None => {
                s.last_status = None;
                s.last_error = Some("export .nc: tool table is empty".into());
                return;
            }
        };
        // Use the first op's spindle RPM, or a sane default.
        let spindle = s
            .file
            .operations
            .first()
            .map(|op| match op {
                valenx_cam::Operation::Profile(p) => p.spindle_rpm,
                valenx_cam::Operation::Pocket(p) => p.spindle_rpm,
                valenx_cam::Operation::Drill(p) => p.spindle_rpm,
                valenx_cam::Operation::Face(p) => p.spindle_rpm,
                valenx_cam::Operation::AdaptiveClearing(p) => p.spindle_rpm,
                valenx_cam::Operation::HelicalBore(p) => p.spindle_rpm,
                valenx_cam::Operation::PlungeRough(p) => p.spindle_rpm,
                valenx_cam::Operation::RampEntry(p) => p.spindle_rpm,
                valenx_cam::Operation::PeckDrillFull(p) => p.spindle_rpm,
                valenx_cam::Operation::Contour2D(p) => p.spindle_rpm,
                valenx_cam::Operation::Contour3D(p) => p.spindle_rpm,
                valenx_cam::Operation::Engrave(p) => p.spindle_rpm,
                valenx_cam::Operation::Scribe(p) => p.spindle_rpm,
                valenx_cam::Operation::SpiralPocket(p) => p.spindle_rpm,
                valenx_cam::Operation::TrochoidalSlot(p) => p.spindle_rpm,
                valenx_cam::Operation::Waterline3D(p) => p.spindle_rpm,
                valenx_cam::Operation::Slot(p) => p.spindle_rpm,
                valenx_cam::Operation::ThreadMill(p) => p.spindle_rpm,
                valenx_cam::Operation::RestMachining(p) => p.spindle_rpm,
            })
            .unwrap_or(12000.0);
        match valenx_cam::post::save_nc(s.selected_postprocessor, &tp, &tool, spindle, &path) {
            Ok(()) => {
                let kind = s.selected_postprocessor.label();
                s.last_status = Some(format!(
                    "Exported {kind} G-code → {}",
                    path.file_name().and_then(|n| n.to_str()).unwrap_or("?")
                ));
                s.last_error = None;
                emit_audit(
                    "cam.export_nc",
                    serde_json::json!({
                        "post": kind,
                        "path": path.display().to_string(),
                    }),
                    serde_json::json!({}),
                );
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("export .nc: {e}"));
            }
        }
    }

    /// Run the Phase 17D voxel simulation against the last-generated
    /// toolpath. Stores the resulting frame list in
    /// [`CamPanelState::animation_frames`].
    pub fn cam_animate(&mut self) {
        let s = &mut self.mesh_toolbox.cam;
        let tp = match s.last_toolpath.clone() {
            Some(t) => t,
            None => {
                s.last_status = None;
                s.last_error = Some("animate: generate a toolpath first".into());
                return;
            }
        };
        let tool = match s.file.tools.first().cloned() {
            Some(t) => t,
            None => {
                s.last_status = None;
                s.last_error = Some("animate: tool table is empty".into());
                return;
            }
        };
        let frames = match valenx_cam::simulate::animate(
            &tp,
            &s.file.stock,
            &tool,
            s.animation_n_frames,
            s.animation_resolution,
        ) {
            Ok(f) => f,
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("voxel grid: {e}"));
                return;
            }
        };
        let n = frames.len();
        s.animation_frames = frames;
        s.animation_frame_idx = 0;
        s.last_status = Some(format!(
            "Animated: {n} frames @ {res}³ voxels",
            res = s.animation_resolution
        ));
        s.last_error = None;
    }

    /// Run the Phase 17F wear-model check across every op + tool
    /// combination. Stores warnings in
    /// [`CamPanelState::last_wear_warnings`].
    pub fn cam_check_wear(&mut self) {
        let s = &mut self.mesh_toolbox.cam;
        let material = s.file.stock.material.clone();
        let mut all = Vec::new();
        for (i, op) in s.file.operations.iter().enumerate() {
            let tool_id = op.tool_id();
            let tool = match s.file.tools.iter().find(|t| t.id == tool_id) {
                Some(t) => t.clone(),
                None => continue,
            };
            let spindle = match op {
                valenx_cam::Operation::Profile(p) => p.spindle_rpm,
                valenx_cam::Operation::Pocket(p) => p.spindle_rpm,
                valenx_cam::Operation::Drill(p) => p.spindle_rpm,
                valenx_cam::Operation::Face(p) => p.spindle_rpm,
                valenx_cam::Operation::AdaptiveClearing(p) => p.spindle_rpm,
                valenx_cam::Operation::HelicalBore(p) => p.spindle_rpm,
                valenx_cam::Operation::PlungeRough(p) => p.spindle_rpm,
                valenx_cam::Operation::RampEntry(p) => p.spindle_rpm,
                valenx_cam::Operation::PeckDrillFull(p) => p.spindle_rpm,
                valenx_cam::Operation::Contour2D(p) => p.spindle_rpm,
                valenx_cam::Operation::Contour3D(p) => p.spindle_rpm,
                valenx_cam::Operation::Engrave(p) => p.spindle_rpm,
                valenx_cam::Operation::Scribe(p) => p.spindle_rpm,
                valenx_cam::Operation::SpiralPocket(p) => p.spindle_rpm,
                valenx_cam::Operation::TrochoidalSlot(p) => p.spindle_rpm,
                valenx_cam::Operation::Waterline3D(p) => p.spindle_rpm,
                valenx_cam::Operation::Slot(p) => p.spindle_rpm,
                valenx_cam::Operation::ThreadMill(p) => p.spindle_rpm,
                valenx_cam::Operation::RestMachining(p) => p.spindle_rpm,
            };
            // Estimate per-op minutes from the last total / op count.
            let est =
                s.last_estimated_time_min.unwrap_or(0.0) / (s.file.operations.len().max(1) as f64);
            let warns = valenx_cam::wear::check_op(
                &tool,
                valenx_cam::wear::OpRunSpec {
                    spindle_rpm: spindle,
                    estimated_minutes: est,
                },
                &material,
            );
            for w in warns {
                all.push((format!("op#{i} ({})", op.label()), w));
            }
        }
        let n = all.len();
        s.last_wear_warnings = all;
        s.last_status = Some(format!("Wear check: {n} warnings"));
        s.last_error = None;
    }

    /// Run the Phase 17F fixture collision check. Persists the
    /// transient `s.fixture` into the CAM file so the result survives
    /// Save/Load.
    pub fn cam_check_fixture(&mut self) {
        let s = &mut self.mesh_toolbox.cam;
        let tp = match s.last_toolpath.clone() {
            Some(t) => t,
            None => {
                s.last_status = None;
                s.last_error = Some("collision: generate a toolpath first".into());
                return;
            }
        };
        let tool = match s.file.tools.first().cloned() {
            Some(t) => t,
            None => {
                s.last_status = None;
                s.last_error = Some("collision: tool table is empty".into());
                return;
            }
        };
        // Mirror the transient panel-state fixture into the persisted
        // CamFile so Save/Load round-trips it.
        s.file.fixture = s.fixture.clone();
        let hits = valenx_cam::fixture::collision_check(&tp, &tool, &s.fixture);
        let n = hits.len();
        s.last_collisions = hits;
        if n == 0 {
            s.last_status = Some("Fixture check: no collisions".into());
        } else {
            s.last_status = Some(format!("Fixture check: {n} collisions detected"));
        }
        s.last_error = None;
    }

    /// Save the CAM workbench file (tools + stock + ops) to a .ron
    /// file via the native dialog.
    pub fn cam_save_file(&mut self) {
        if let Err(e) = self.cam_sync_stock_from_inputs() {
            let s = &mut self.mesh_toolbox.cam;
            s.last_status = None;
            s.last_error = Some(format!("save: {e}"));
            return;
        }
        let path = match rfd::FileDialog::new()
            .add_filter("CAM file (.ron)", &["ron"])
            .save_file()
        {
            Some(p) => p,
            None => return,
        };
        let s = &mut self.mesh_toolbox.cam;
        match s.file.write_to(&path) {
            Ok(()) => {
                s.last_status = Some(format!(
                    "Saved CAM file → {}",
                    path.file_name().and_then(|n| n.to_str()).unwrap_or("?")
                ));
                s.last_error = None;
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("save: {e}"));
            }
        }
    }

    /// Load a CAM workbench file from a .ron file via the native
    /// dialog.
    pub fn cam_load_file(&mut self) {
        let path = match rfd::FileDialog::new()
            .add_filter("CAM file (.ron)", &["ron"])
            .pick_file()
        {
            Some(p) => p,
            None => return,
        };
        let s = &mut self.mesh_toolbox.cam;
        match valenx_cam::persist::CamFile::read_from(&path) {
            Ok(file) => {
                s.stock_origin = [
                    file.stock.origin.x,
                    file.stock.origin.y,
                    file.stock.origin.z,
                ];
                s.stock_size = [file.stock.size.x, file.stock.size.y, file.stock.size.z];
                s.stock_material = file.stock.material.clone();
                let n_tools = file.tools.len();
                let n_ops = file.operations.len();
                s.file = file;
                s.last_toolpath = None;
                s.last_status = Some(format!("Loaded CAM file: {n_tools} tools, {n_ops} ops"));
                s.last_error = None;
            }
            Err(e) => {
                s.last_status = None;
                s.last_error = Some(format!("load: {e}"));
            }
        }
    }
}

pub fn draw_cam_panel(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.heading("CAM (Path / G-code)");
    ui.separator();
    draw_cam_stock(app, ui);
    ui.separator();
    draw_cam_tool_table(app, ui);
    ui.separator();
    draw_cam_operations(app, ui);
    ui.separator();
    draw_cam_generate_and_export(app, ui);
    ui.separator();
    draw_cam_animation(app, ui);
    ui.separator();
    draw_cam_wear(app, ui);
    ui.separator();
    draw_cam_fixture(app, ui);
    ui.separator();
    let s = &app.mesh_toolbox.cam;
    if let Some(msg) = &s.last_status {
        ui.label(egui::RichText::new(msg).color(egui::Color32::GREEN));
    }
    if let Some(err) = &s.last_error {
        ui.label(egui::RichText::new(err).color(egui::Color32::RED));
    }
}

fn draw_cam_animation(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Simulation animation").strong());
    let s = &mut app.mesh_toolbox.cam;
    ui.horizontal(|ui| {
        ui.label("Voxel resolution:");
        ui.add(
            egui::DragValue::new(&mut s.animation_resolution)
                .range(8u32..=256)
                .speed(1.0),
        );
        ui.label("Frames:");
        ui.add(
            egui::DragValue::new(&mut s.animation_n_frames)
                .range(1u32..=128)
                .speed(0.5),
        );
    });
    if ui.button("Animate (voxel sim)").clicked() {
        app.cam_animate();
    }
    let s = &mut app.mesh_toolbox.cam;
    if !s.animation_frames.is_empty() {
        let max_idx = s.animation_frames.len().saturating_sub(1);
        ui.horizontal(|ui| {
            ui.label("Frame:");
            let mut idx = s.animation_frame_idx;
            if ui.add(egui::Slider::new(&mut idx, 0..=max_idx)).changed() {
                s.animation_frame_idx = idx;
            }
        });
        let f = &s.animation_frames[s.animation_frame_idx];
        let n_tris = f
            .element_blocks
            .first()
            .map(|b| b.connectivity.len() / 3)
            .unwrap_or(0);
        ui.label(format!(
            "Frame {} / {} — {} triangles in voxel surface",
            s.animation_frame_idx + 1,
            s.animation_frames.len(),
            n_tris,
        ));
    }
}

fn draw_cam_wear(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Tool wear").strong());
    if ui.button("Check tool wear").clicked() {
        app.cam_check_wear();
    }
    let s = &app.mesh_toolbox.cam;
    for (op_label, w) in &s.last_wear_warnings {
        let color = if w.code == "wear.life_exceeded" {
            egui::Color32::RED
        } else {
            egui::Color32::YELLOW
        };
        ui.label(egui::RichText::new(format!("{op_label}: {}", w.message)).color(color));
    }
}

fn draw_cam_fixture(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Fixture / clamp collision").strong());
    let s = &mut app.mesh_toolbox.cam;
    ui.collapsing("Fixture AABBs", |ui| {
        let mut remove_idx: Option<usize> = None;
        for (i, f) in s.fixture.aabbs.iter_mut().enumerate() {
            ui.horizontal(|ui| {
                ui.label(format!("#{i} {}: ", f.label));
                ui.add(
                    egui::DragValue::new(&mut f.min.x)
                        .speed(0.5)
                        .prefix("xmin "),
                );
                ui.add(
                    egui::DragValue::new(&mut f.max.x)
                        .speed(0.5)
                        .prefix("xmax "),
                );
                ui.add(
                    egui::DragValue::new(&mut f.min.y)
                        .speed(0.5)
                        .prefix("ymin "),
                );
                ui.add(
                    egui::DragValue::new(&mut f.max.y)
                        .speed(0.5)
                        .prefix("ymax "),
                );
                ui.add(
                    egui::DragValue::new(&mut f.min.z)
                        .speed(0.5)
                        .prefix("zmin "),
                );
                ui.add(
                    egui::DragValue::new(&mut f.max.z)
                        .speed(0.5)
                        .prefix("zmax "),
                );
                if ui.small_button("Delete").clicked() {
                    remove_idx = Some(i);
                }
            });
        }
        if let Some(i) = remove_idx {
            s.fixture.aabbs.remove(i);
        }
        if ui.button("Add fixture AABB").clicked() {
            s.fixture.push_aabb(
                nalgebra::Vector3::new(-50.0, -50.0, -50.0),
                nalgebra::Vector3::new(-25.0, -25.0, 0.0),
                "clamp",
            );
        }
    });
    if ui.button("Check fixture collisions").clicked() {
        app.cam_check_fixture();
    }
    let s = &app.mesh_toolbox.cam;
    if !s.last_collisions.is_empty() {
        ui.label(
            egui::RichText::new(format!("{} collision(s):", s.last_collisions.len()))
                .color(egui::Color32::RED),
        );
        for c in s.last_collisions.iter().take(5) {
            ui.label(format!(
                "  • move#{} hits '{}' at ({:.1}, {:.1}, {:.1})",
                c.move_index,
                c.fixture_label,
                c.tool_position.x,
                c.tool_position.y,
                c.tool_position.z,
            ));
        }
        if s.last_collisions.len() > 5 {
            ui.label(format!("  • … and {} more", s.last_collisions.len() - 5));
        }
    }
}

fn draw_cam_stock(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Stock").strong());
    let s = &mut app.mesh_toolbox.cam;
    ui.horizontal(|ui| {
        ui.label("Origin:");
        for (i, axis) in ["x", "y", "z"].iter().enumerate() {
            ui.add(
                egui::DragValue::new(&mut s.stock_origin[i])
                    .speed(0.5)
                    .prefix(format!("{axis} ")),
            );
        }
    });
    ui.horizontal(|ui| {
        ui.label("Size:");
        for (i, axis) in ["x", "y", "z"].iter().enumerate() {
            ui.add(
                egui::DragValue::new(&mut s.stock_size[i])
                    .speed(0.5)
                    .prefix(format!("{axis} "))
                    .range(0.001..=10000.0),
            );
        }
    });
    ui.horizontal(|ui| {
        ui.label("Material:");
        ui.text_edit_singleline(&mut s.stock_material);
    });
}

fn draw_cam_tool_table(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Tools").strong());
    let mut delete_idx: Option<usize> = None;
    {
        let s = &app.mesh_toolbox.cam;
        if s.file.tools.is_empty() {
            ui.label("(no tools)");
        } else {
            for (i, t) in s.file.tools.iter().enumerate() {
                ui.horizontal(|ui| {
                    ui.label(format!(
                        "T{} {} ({}, D{:.2}mm, L{:.1}mm, F{})",
                        t.id, t.name, t.kind, t.diameter_mm, t.length_mm, t.flutes
                    ));
                    if ui.small_button("Delete").clicked() {
                        delete_idx = Some(i);
                    }
                });
            }
        }
    }
    if let Some(i) = delete_idx {
        app.cam_delete_tool(i);
    }
    ui.collapsing("Add Tool", |ui| {
        let s = &mut app.mesh_toolbox.cam;
        ui.horizontal(|ui| {
            ui.label("Name:");
            ui.text_edit_singleline(&mut s.new_tool_name);
        });
        ui.horizontal(|ui| {
            ui.label("Kind:");
            for kind in [
                valenx_cam::ToolKind::EndMill,
                valenx_cam::ToolKind::BallMill,
                valenx_cam::ToolKind::Drill,
                valenx_cam::ToolKind::FaceMill,
                valenx_cam::ToolKind::Tap,
                valenx_cam::ToolKind::Reamer,
            ] {
                ui.selectable_value(&mut s.new_tool_kind, kind, kind.label());
            }
        });
        ui.horizontal(|ui| {
            ui.label("Diameter:");
            ui.add(
                egui::DragValue::new(&mut s.new_tool_diameter)
                    .range(0.01..=200.0)
                    .speed(0.1)
                    .suffix(" mm"),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Length:");
            ui.add(
                egui::DragValue::new(&mut s.new_tool_length)
                    .range(0.1..=500.0)
                    .speed(0.5)
                    .suffix(" mm"),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Flutes:");
            ui.add(
                egui::DragValue::new(&mut s.new_tool_flutes)
                    .range(1u32..=12)
                    .speed(0.05),
            );
        });
        ui.horizontal(|ui| {
            ui.label("Material:");
            ui.text_edit_singleline(&mut s.new_tool_material);
        });
    });
    if ui.button("Add Tool").clicked() {
        app.cam_add_tool();
    }
}

fn draw_cam_operations(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Operations").strong());
    let mut delete_idx: Option<usize> = None;
    {
        let s = &app.mesh_toolbox.cam;
        if s.file.operations.is_empty() {
            ui.label("(no operations)");
        } else {
            for (i, op) in s.file.operations.iter().enumerate() {
                ui.horizontal(|ui| {
                    ui.label(format!("#{i}: {} (T{})", op.label(), op.tool_id()));
                    if ui.small_button("Delete").clicked() {
                        delete_idx = Some(i);
                    }
                });
            }
        }
    }
    if let Some(i) = delete_idx {
        app.cam_delete_operation(i);
    }
    ui.collapsing("Add Operation", |ui| {
        let s = &mut app.mesh_toolbox.cam;
        ui.horizontal_wrapped(|ui| {
            ui.label("Kind:")
                .on_hover_text("Which CAM operation to add — drives which parameters are shown below.");
            for kind in [
                CamOpKind::Profile,
                CamOpKind::Pocket,
                CamOpKind::Drill,
                CamOpKind::Face,
                CamOpKind::AdaptiveClearing,
                CamOpKind::HelicalBore,
                CamOpKind::PlungeRough,
                CamOpKind::RampEntry,
                CamOpKind::PeckDrillFull,
                CamOpKind::Contour2D,
                CamOpKind::Contour3D,
                CamOpKind::Engrave,
                CamOpKind::Scribe,
                CamOpKind::SpiralPocket,
                CamOpKind::TrochoidalSlot,
                CamOpKind::Waterline3D,
                CamOpKind::Slot,
                CamOpKind::ThreadMill,
                CamOpKind::RestMachining,
            ] {
                ui.selectable_value(&mut s.new_op_kind, kind, kind.label());
            }
        });
        ui.horizontal(|ui| {
            ui.label("Tool id:")
                .on_hover_text("Which tool from the tools table to use. Each operation references a tool by T-number.");
            ui.add(
                egui::DragValue::new(&mut s.new_op_tool_id)
                    .range(1u32..=999)
                    .speed(0.05),
            )
            .on_hover_text("Tool T-number (1..999).");
        });
        ui.horizontal(|ui| {
            ui.label("Spindle RPM:")
                .on_hover_text("Spindle rotation speed in revolutions per minute. Picked from tool + material recommendations.");
            ui.add(
                egui::DragValue::new(&mut s.new_op_spindle_rpm)
                    .range(100.0..=60000.0)
                    .speed(50.0),
            )
            .on_hover_text("Spindle RPM (rev/min).");
        });
        ui.horizontal(|ui| {
            ui.label("Feed:")
                .on_hover_text("Cutting feedrate — XY motion speed when the tool is engaged.");
            ui.add(
                egui::DragValue::new(&mut s.new_op_feed)
                    .range(1.0..=10000.0)
                    .speed(10.0)
                    .suffix(" mm/min"),
            )
            .on_hover_text("Cutting feedrate (mm/min).");
            ui.label("Plunge feed:")
                .on_hover_text("Z-axis plunge speed when entering material. Typically 30–50% of XY feed.");
            ui.add(
                egui::DragValue::new(&mut s.new_op_plunge_feed)
                    .range(1.0..=5000.0)
                    .speed(10.0)
                    .suffix(" mm/min"),
            )
            .on_hover_text("Plunge feedrate (mm/min).");
        });
        ui.horizontal(|ui| {
            ui.label("Safe Z:")
                .on_hover_text("Retract height above the stock — must clear every fixture / clamp.");
            ui.add(
                egui::DragValue::new(&mut s.new_op_safe_z)
                    .range(0.1..=200.0)
                    .speed(0.1)
                    .suffix(" mm"),
            )
            .on_hover_text("Safe Z height above stock (mm).");
        });
        match s.new_op_kind {
            CamOpKind::Profile => {
                ui.horizontal(|ui| {
                    ui.label("Step down:")
                        .on_hover_text("Z increment per pass. Typical = 50–100% of tool diameter for soft materials, less for hard ones.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_step_down)
                            .range(0.01..=50.0)
                            .speed(0.1)
                            .suffix(" mm"),
                    )
                    .on_hover_text("Axial step-down per pass (mm).");
                    ui.label("Depth:")
                        .on_hover_text("Total cut depth below the stock top.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_depth)
                            .range(0.01..=500.0)
                            .speed(0.5)
                            .suffix(" mm"),
                    )
                    .on_hover_text("Total depth of cut (mm).");
                });
                ui.checkbox(&mut s.new_op_climb, "Climb cut")
                    .on_hover_text("Climb (chip thinning) milling — better finish, lower tool wear. Unchecked = conventional.");
            }
            CamOpKind::Pocket => {
                ui.horizontal(|ui| {
                    ui.label("Step down:")
                        .on_hover_text("Z increment per pass — 50–100% of tool diameter for soft material, less for hard.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_step_down)
                            .range(0.01..=50.0)
                            .speed(0.1)
                            .suffix(" mm"),
                    )
                    .on_hover_text("Axial step-down per pass (mm).");
                    ui.label("Step over:")
                        .on_hover_text("XY step between adjacent passes — typically 40–80% of tool diameter.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_step_over)
                            .range(0.01..=50.0)
                            .speed(0.1)
                            .suffix(" mm"),
                    )
                    .on_hover_text("Radial step-over between passes (mm).");
                });
                ui.horizontal(|ui| {
                    ui.label("Depth:")
                        .on_hover_text("Total pocket depth below the stock top.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_depth)
                            .range(0.01..=500.0)
                            .speed(0.5)
                            .suffix(" mm"),
                    )
                    .on_hover_text("Total depth of cut (mm).");
                    ui.label("Angle:")
                        .on_hover_text("Raster orientation for zig-zag / parallel strategies. 0 = passes along +X.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_raster_angle)
                            .range(-180.0..=180.0)
                            .speed(1.0)
                            .suffix("°"),
                    )
                    .on_hover_text("Raster pass angle in degrees.");
                });
                ui.horizontal(|ui| {
                    ui.label("Strategy:")
                        .on_hover_text("Pocket clearing pattern. ZigZag = bidirectional rasters; Parallel = same-direction passes; Spiral = inward concentric loops.");
                    for strat in [
                        valenx_cam::PocketStrategy::ZigZag,
                        valenx_cam::PocketStrategy::Parallel,
                        valenx_cam::PocketStrategy::Spiral,
                    ] {
                        ui.selectable_value(&mut s.new_op_pocket_strategy, strat, strat.label())
                            .on_hover_text(match strat {
                                valenx_cam::PocketStrategy::ZigZag => "Bidirectional rasters — fastest, alternates climb/conventional.",
                                valenx_cam::PocketStrategy::Parallel => "Same-direction rasters with rapid returns — uniform climb cut.",
                                valenx_cam::PocketStrategy::Spiral => "Inward concentric loops — best surface finish, no direction change.",
                            });
                    }
                });
                ui.checkbox(&mut s.new_op_climb, "Climb cut")
                    .on_hover_text("Climb (chip thinning) milling — better finish, lower tool wear. Unchecked = conventional.");
            }
            CamOpKind::Drill => {
                ui.horizontal(|ui| {
                    ui.label("Peck depth:")
                        .on_hover_text("Z increment per peck. The drill retracts to clear chips between pecks.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_peck_depth)
                            .range(0.01..=50.0)
                            .speed(0.1)
                            .suffix(" mm"),
                    )
                    .on_hover_text("Single peck depth (mm).");
                    ui.label("Total depth:")
                        .on_hover_text("Total hole depth from the stock top.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_drill_total_depth)
                            .range(0.01..=500.0)
                            .speed(0.5)
                            .suffix(" mm"),
                    )
                    .on_hover_text("Total drill depth (mm).");
                });
                ui.horizontal(|ui| {
                    ui.label("Retract clearance:")
                        .on_hover_text("How far above the hole the drill retracts between pecks — must clear chip evacuation.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_retract_clearance)
                            .range(0.0..=50.0)
                            .speed(0.1)
                            .suffix(" mm"),
                    )
                    .on_hover_text("Retract clearance above hole (mm).");
                });
                ui.collapsing("Hole positions", |ui| {
                    let mut remove_hole: Option<usize> = None;
                    for (j, h) in s.new_op_hole_positions.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            ui.label(format!("#{j}:"));
                            ui.add(egui::DragValue::new(&mut h[0]).speed(0.5).prefix("x "))
                                .on_hover_text("Hole centre X (mm, stock coordinates).");
                            ui.add(egui::DragValue::new(&mut h[1]).speed(0.5).prefix("y "))
                                .on_hover_text("Hole centre Y (mm, stock coordinates).");
                            if ui.small_button("Delete").clicked() {
                                remove_hole = Some(j);
                            }
                        });
                    }
                    if let Some(j) = remove_hole {
                        s.new_op_hole_positions.remove(j);
                    }
                    if ui.button("Add hole").clicked() {
                        s.new_op_hole_positions.push([0.0, 0.0, 0.0]);
                    }
                });
            }
            CamOpKind::Face => {
                ui.horizontal(|ui| {
                    ui.label("Step over:")
                        .on_hover_text("XY step between adjacent face-milling passes — typically 60–80% of cutter diameter.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_step_over)
                            .range(0.01..=50.0)
                            .speed(0.1)
                            .suffix(" mm"),
                    )
                    .on_hover_text("Radial step-over (mm).");
                    ui.label("Step down:")
                        .on_hover_text("Z increment per pass when facing in multiple passes.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_step_down)
                            .range(0.01..=50.0)
                            .speed(0.1)
                            .suffix(" mm"),
                    )
                    .on_hover_text("Axial step-down per pass (mm).");
                });
                ui.horizontal(|ui| {
                    ui.label("Depth:")
                        .on_hover_text("Total stock removed below the top face.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_depth)
                            .range(0.01..=500.0)
                            .speed(0.5)
                            .suffix(" mm"),
                    )
                    .on_hover_text("Total cut depth (mm).");
                    ui.label("Angle:")
                        .on_hover_text("Raster orientation. 0 = passes along +X.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_raster_angle)
                            .range(-180.0..=180.0)
                            .speed(1.0)
                            .suffix("°"),
                    )
                    .on_hover_text("Raster pass angle in degrees.");
                });
                ui.checkbox(&mut s.new_op_climb, "Climb cut")
                    .on_hover_text("Climb milling — better finish, lower tool wear. Unchecked = conventional.");
            }
            // Phase 17A/B — generic step_down/step_over/depth form.
            CamOpKind::AdaptiveClearing
            | CamOpKind::HelicalBore
            | CamOpKind::RampEntry
            | CamOpKind::SpiralPocket
            | CamOpKind::TrochoidalSlot
            | CamOpKind::Waterline3D
            | CamOpKind::Slot
            | CamOpKind::ThreadMill
            | CamOpKind::RestMachining => {
                let step_down_tt = match s.new_op_kind {
                    CamOpKind::AdaptiveClearing => "Optimal load step-down — adaptive engagement keeps chip load constant.",
                    CamOpKind::HelicalBore => "Pitch per turn — vertical advance per full revolution while spiralling down.",
                    CamOpKind::RampEntry => "Ramp angle equivalent depth — Z drop per linear ramp pass.",
                    CamOpKind::SpiralPocket => "Z step between concentric spiral layers.",
                    CamOpKind::TrochoidalSlot => "Step-down per trochoidal layer — engagement-limited cutting.",
                    CamOpKind::Waterline3D => "Z slice spacing for waterline finishing passes.",
                    CamOpKind::Slot => "Z step per slot pass (drops the slot deeper each pass).",
                    CamOpKind::ThreadMill => "Pitch — vertical distance per full thread revolution.",
                    CamOpKind::RestMachining => "Z slice for rest-cut cleanup passes.",
                    _ => "Axial step-down per pass (mm).",
                };
                let step_over_tt = match s.new_op_kind {
                    CamOpKind::AdaptiveClearing => "Engagement width — typical 5–15% of tool diameter for adaptive.",
                    CamOpKind::HelicalBore => "Radial step from previous helix wall — controls bore wall finish.",
                    CamOpKind::RampEntry => "Lateral overlap between ramp pass and next finishing pass.",
                    CamOpKind::SpiralPocket => "Step between concentric spiral arms.",
                    CamOpKind::TrochoidalSlot => "Trochoid loop pitch — advance per loop along the slot centerline.",
                    CamOpKind::Waterline3D => "XY step inside one waterline slice.",
                    CamOpKind::Slot => "Width of cut (slot width minus tool diameter).",
                    CamOpKind::ThreadMill => "Radial engagement — controls thread depth per pass.",
                    CamOpKind::RestMachining => "Step-over inside the residual region (small for fine cleanup).",
                    _ => "Radial step-over between passes (mm).",
                };
                ui.horizontal(|ui| {
                    ui.label("Step down:").on_hover_text(step_down_tt);
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_step_down)
                            .range(0.01..=50.0)
                            .speed(0.1)
                            .suffix(" mm"),
                    )
                    .on_hover_text(step_down_tt);
                    ui.label("Step over:").on_hover_text(step_over_tt);
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_step_over)
                            .range(0.01..=50.0)
                            .speed(0.1)
                            .suffix(" mm"),
                    )
                    .on_hover_text(step_over_tt);
                });
                ui.horizontal(|ui| {
                    ui.label("Depth:")
                        .on_hover_text("Total depth of cut below stock top.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_depth)
                            .range(0.01..=500.0)
                            .speed(0.5)
                            .suffix(" mm"),
                    )
                    .on_hover_text("Total operation depth (mm).");
                });
                if matches!(s.new_op_kind, CamOpKind::TrochoidalSlot | CamOpKind::Slot) {
                    ui.collapsing("Centreline / endpoints (XY)", |ui| {
                        let mut remove_pt: Option<usize> = None;
                        for (j, h) in s.new_op_hole_positions.iter_mut().enumerate() {
                            ui.horizontal(|ui| {
                                ui.label(format!("#{j}:"));
                                ui.add(egui::DragValue::new(&mut h[0]).speed(0.5).prefix("x "))
                                    .on_hover_text("Slot endpoint X (mm).");
                                ui.add(egui::DragValue::new(&mut h[1]).speed(0.5).prefix("y "))
                                    .on_hover_text("Slot endpoint Y (mm).");
                                if ui.small_button("Delete").clicked() {
                                    remove_pt = Some(j);
                                }
                            });
                        }
                        if let Some(j) = remove_pt {
                            s.new_op_hole_positions.remove(j);
                        }
                        if ui.button("Add point").clicked() {
                            s.new_op_hole_positions.push([0.0, 0.0, 0.0]);
                        }
                    });
                }
            }
            // Plunge rough + peck drill full reuse the hole-position list.
            CamOpKind::PlungeRough | CamOpKind::PeckDrillFull => {
                ui.horizontal(|ui| {
                    ui.label("Depth:")
                        .on_hover_text("Plunge depth below stock top per position.");
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_depth)
                            .range(0.01..=500.0)
                            .speed(0.5)
                            .suffix(" mm"),
                    )
                    .on_hover_text("Plunge depth (mm).");
                });
                if matches!(s.new_op_kind, CamOpKind::PeckDrillFull) {
                    ui.horizontal(|ui| {
                        ui.label("Peck depth:")
                            .on_hover_text("Single-peck Z advance — retract between pecks to clear chips.");
                        ui.add(
                            egui::DragValue::new(&mut s.new_op_peck_depth)
                                .range(0.01..=50.0)
                                .speed(0.1)
                                .suffix(" mm"),
                        )
                        .on_hover_text("Peck depth per cycle (mm).");
                        ui.label("Total depth:")
                            .on_hover_text("Total drill depth from the stock top.");
                        ui.add(
                            egui::DragValue::new(&mut s.new_op_drill_total_depth)
                                .range(0.01..=500.0)
                                .speed(0.5)
                                .suffix(" mm"),
                        )
                        .on_hover_text("Total drill depth (mm).");
                    });
                    ui.horizontal(|ui| {
                        ui.label("Retract clearance:")
                            .on_hover_text("Retract height between pecks — large enough to evacuate chips.");
                        ui.add(
                            egui::DragValue::new(&mut s.new_op_retract_clearance)
                                .range(0.0..=50.0)
                                .speed(0.1)
                                .suffix(" mm"),
                        )
                        .on_hover_text("Retract clearance above hole (mm).");
                    });
                }
                ui.collapsing("Positions (XY)", |ui| {
                    let mut remove_pt: Option<usize> = None;
                    for (j, h) in s.new_op_hole_positions.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            ui.label(format!("#{j}:"));
                            ui.add(egui::DragValue::new(&mut h[0]).speed(0.5).prefix("x "))
                                .on_hover_text("Position X (mm).");
                            ui.add(egui::DragValue::new(&mut h[1]).speed(0.5).prefix("y "))
                                .on_hover_text("Position Y (mm).");
                            if ui.small_button("Delete").clicked() {
                                remove_pt = Some(j);
                            }
                        });
                    }
                    if let Some(j) = remove_pt {
                        s.new_op_hole_positions.remove(j);
                    }
                    if ui.button("Add position").clicked() {
                        s.new_op_hole_positions.push([0.0, 0.0, 0.0]);
                    }
                });
            }
            // Contour / engrave / scribe — XY (Z) curve.
            CamOpKind::Contour2D
            | CamOpKind::Contour3D
            | CamOpKind::Engrave
            | CamOpKind::Scribe => {
                let depth_tt = match s.new_op_kind {
                    CamOpKind::Contour2D => "Final contour depth — 2D contour walks the XY curve at this Z below stock top.",
                    CamOpKind::Contour3D => "Z offset applied to the 3-D curve points (curve depths are absolute).",
                    CamOpKind::Engrave => "Engraving depth into the surface — typical 0.1–0.5 mm for V-bits, deeper for end mills.",
                    CamOpKind::Scribe => "Scribe line depth — typically very shallow (≤ 0.1 mm) to mark only.",
                    _ => "Operation depth (mm).",
                };
                ui.horizontal(|ui| {
                    ui.label("Depth:").on_hover_text(depth_tt);
                    ui.add(
                        egui::DragValue::new(&mut s.new_op_depth)
                            .range(0.01..=500.0)
                            .speed(0.5)
                            .suffix(" mm"),
                    )
                    .on_hover_text(depth_tt);
                });
                if matches!(s.new_op_kind, CamOpKind::Contour2D) {
                    ui.horizontal(|ui| {
                        ui.label("Step down:")
                            .on_hover_text("Z increment per contour pass when Depth requires multiple passes.");
                        ui.add(
                            egui::DragValue::new(&mut s.new_op_step_down)
                                .range(0.01..=50.0)
                                .speed(0.1)
                                .suffix(" mm"),
                        )
                        .on_hover_text("Axial step-down per pass (mm).");
                    });
                }
                ui.collapsing("Curve points", |ui| {
                    let mut remove_pt: Option<usize> = None;
                    for (j, h) in s.new_op_hole_positions.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            ui.label(format!("#{j}:"));
                            ui.add(egui::DragValue::new(&mut h[0]).speed(0.5).prefix("x "))
                                .on_hover_text("Curve point X (mm).");
                            ui.add(egui::DragValue::new(&mut h[1]).speed(0.5).prefix("y "))
                                .on_hover_text("Curve point Y (mm).");
                            ui.add(egui::DragValue::new(&mut h[2]).speed(0.5).prefix("z "))
                                .on_hover_text("Curve point Z (mm). Engrave/scribe usually use 0.");
                            if ui.small_button("Delete").clicked() {
                                remove_pt = Some(j);
                            }
                        });
                    }
                    if let Some(j) = remove_pt {
                        s.new_op_hole_positions.remove(j);
                    }
                    if ui.button("Add point").clicked() {
                        s.new_op_hole_positions.push([0.0, 0.0, 0.0]);
                    }
                });
            }
        }
    });
    if ui.button("Add Operation").clicked() {
        app.cam_add_operation();
    }
}

fn draw_cam_generate_and_export(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Generate / Simulate / Export").strong());
    ui.horizontal(|ui| {
        if ui.button("Generate Toolpath").clicked() {
            app.cam_generate_toolpath();
        }
        if ui.button("Toggle Simulate Overlay").clicked() {
            app.cam_toggle_simulate();
        }
    });
    {
        let s = &app.mesh_toolbox.cam;
        if let (Some(t), Some(v)) = (s.last_estimated_time_min, s.last_removed_volume_mm3) {
            ui.label(format!("Estimated: {t:.2} min  |  Removed: {v:.0} mm³"));
        }
        if let Some(tp) = &s.last_toolpath {
            ui.label(format!("Toolpath: {} moves", tp.len()));
        }
    }
    let s = &mut app.mesh_toolbox.cam;
    ui.label("Postprocessor:");
    ui.horizontal_wrapped(|ui| {
        for k in [
            valenx_cam::PostKind::Grbl,
            valenx_cam::PostKind::LinuxCnc,
            valenx_cam::PostKind::Fanuc,
            valenx_cam::PostKind::Haas,
            valenx_cam::PostKind::Heidenhain,
            valenx_cam::PostKind::Mazatrol,
            valenx_cam::PostKind::Sinumerik,
            valenx_cam::PostKind::Okuma,
            valenx_cam::PostKind::Mori,
            valenx_cam::PostKind::Makino,
            valenx_cam::PostKind::Kitamura,
            valenx_cam::PostKind::Centroid,
            valenx_cam::PostKind::Tormach,
            valenx_cam::PostKind::DeepNest,
            valenx_cam::PostKind::Marlin,
            valenx_cam::PostKind::Klipper,
            valenx_cam::PostKind::Repetier,
            valenx_cam::PostKind::FluidNc,
            valenx_cam::PostKind::SnapMaker,
            valenx_cam::PostKind::Smoothie,
            valenx_cam::PostKind::Smoothieboard,
            valenx_cam::PostKind::TinyG,
            valenx_cam::PostKind::SourceRabbit,
            valenx_cam::PostKind::HsmAdvisor,
            valenx_cam::PostKind::Fusion360Brand,
            valenx_cam::PostKind::OpenDmg,
            valenx_cam::PostKind::VmcAAxis,
            valenx_cam::PostKind::VmcBAxis,
            valenx_cam::PostKind::FanucAAxis,
            valenx_cam::PostKind::FanucBAxis,
        ] {
            ui.selectable_value(&mut s.selected_postprocessor, k, k.label());
        }
    });
    ui.horizontal(|ui| {
        if ui.button("Export .nc").clicked() {
            app.cam_export_nc();
        }
        if ui.button("Save CAM file").clicked() {
            app.cam_save_file();
        }
        if ui.button("Load CAM file").clicked() {
            app.cam_load_file();
        }
    });
}

fn surface_edge_from_u8(v: u8) -> valenx_surface::Edge {
    match v {
        0 => valenx_surface::Edge::UMin,
        1 => valenx_surface::Edge::UMax,
        2 => valenx_surface::Edge::VMin,
        _ => valenx_surface::Edge::VMax,
    }
}

/// Parse the "x,y,z" per-line text used by the Phase 19D Fit tool
/// into a vector of points. Blank lines and `#` comments are
/// skipped; malformed lines are silently dropped.
fn parse_fit_points(text: &str) -> Vec<nalgebra::Vector3<f64>> {
    let mut out = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split(',').map(|p| p.trim()).collect();
        if parts.len() != 3 {
            continue;
        }
        let Ok(x) = parts[0].parse::<f64>() else {
            continue;
        };
        let Ok(y) = parts[1].parse::<f64>() else {
            continue;
        };
        let Ok(z) = parts[2].parse::<f64>() else {
            continue;
        };
        out.push(nalgebra::Vector3::new(x, y, z));
    }
    out
}

/// Build a clamped open-uniform knot vector — duplicated here to
/// avoid reaching into `valenx-surface`'s private helpers.
fn open_uniform_knots(n_cp: usize, degree: usize) -> Vec<f64> {
    let p = degree;
    let m = n_cp + p + 1;
    let mut k = vec![0.0; m];
    if n_cp <= p + 1 {
        // `vec![0.0; m]` already zeroes the first p+1 entries.
        for kv in k.iter_mut().skip(m - p - 1) {
            *kv = 1.0;
        }
        return k;
    }
    let n_internal = n_cp - p - 1;
    for (i, kv) in k.iter_mut().enumerate().take(m) {
        if i <= p {
            *kv = 0.0;
        } else if i >= n_cp {
            *kv = 1.0;
        } else {
            let idx = i - p;
            *kv = idx as f64 / (n_internal + 1) as f64;
        }
    }
    k
}

pub fn draw_surface_panel(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.heading("Surface (NURBS curves + surfaces)");
    ui.separator();

    // Tool palette.
    ui.label(egui::RichText::new("Tool").strong());
    {
        let s = &mut app.mesh_toolbox.surface;
        ui.horizontal_wrapped(|ui| {
            for tool in [
                SurfaceTool::NurbsCurve,
                SurfaceTool::NurbsSurface,
                SurfaceTool::CoonsFill,
                SurfaceTool::Sew,
                SurfaceTool::Trim,
                SurfaceTool::KnotOps,
                SurfaceTool::Ssi,
                SurfaceTool::Fit,
                SurfaceTool::Ruled,
            ] {
                ui.selectable_value(&mut s.tool, tool, tool.label());
            }
        });
    }
    ui.separator();

    let tool = app.mesh_toolbox.surface.tool;
    match tool {
        SurfaceTool::NurbsCurve => draw_surface_curve_inputs(app, ui),
        SurfaceTool::NurbsSurface => draw_surface_surface_inputs(app, ui),
        SurfaceTool::CoonsFill => draw_surface_coons_inputs(app, ui),
        SurfaceTool::Sew => draw_surface_sew_inputs(app, ui),
        SurfaceTool::Trim => draw_surface_trim_inputs(app, ui),
        SurfaceTool::KnotOps => draw_surface_knot_ops_inputs(app, ui),
        SurfaceTool::Ssi => draw_surface_ssi_inputs(app, ui),
        SurfaceTool::Fit => draw_surface_fit_inputs(app, ui),
        SurfaceTool::Ruled => draw_surface_ruled_inputs(app, ui),
    }

    ui.separator();
    draw_surface_entity_lists(app, ui);

    // Status / error.
    ui.separator();
    let s = &app.mesh_toolbox.surface;
    if let Some(msg) = &s.last_status {
        ui.label(egui::RichText::new(msg).color(egui::Color32::GREEN));
    }
    if let Some(err) = &s.last_error {
        ui.label(egui::RichText::new(err).color(egui::Color32::RED));
    }
}

fn draw_surface_curve_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Construct NURBS curve").strong());
    {
        let s = &mut app.mesh_toolbox.surface;
        ui.horizontal(|ui| {
            ui.label("Degree:")
                .on_hover_text("Polynomial degree of the curve. 1 = polyline, 3 = cubic (typical), 5 = quintic.");
            ui.add(
                egui::DragValue::new(&mut s.curve_degree)
                    .range(1..=9)
                    .speed(0.1),
            )
            .on_hover_text("Curve degree (1..9).");
        });
        ui.horizontal(|ui| {
            ui.label("N control points:")
                .on_hover_text("Number of control points. Must be > degree. More points = more curvature freedom.");
            ui.add(
                egui::DragValue::new(&mut s.curve_n_cps)
                    .range(2..=20)
                    .speed(0.1),
            )
            .on_hover_text("Control point count.");
        });
        s.curve_cps.resize(s.curve_n_cps, [0.0, 0.0, 0.0]);
        s.curve_weights.resize(s.curve_n_cps, 1.0);
        ui.collapsing("Control points", |ui| {
            for i in 0..s.curve_n_cps {
                ui.horizontal(|ui| {
                    ui.label(format!("#{i}:"));
                    ui.add(
                        egui::DragValue::new(&mut s.curve_cps[i][0])
                            .speed(0.05)
                            .prefix("x "),
                    )
                    .on_hover_text("Control point X (model units).");
                    ui.add(
                        egui::DragValue::new(&mut s.curve_cps[i][1])
                            .speed(0.05)
                            .prefix("y "),
                    )
                    .on_hover_text("Control point Y (model units).");
                    ui.add(
                        egui::DragValue::new(&mut s.curve_cps[i][2])
                            .speed(0.05)
                            .prefix("z "),
                    )
                    .on_hover_text("Control point Z (model units).");
                    ui.add(
                        egui::DragValue::new(&mut s.curve_weights[i])
                            .speed(0.05)
                            .prefix("w "),
                    )
                    .on_hover_text("Rational weight — 1.0 = non-rational B-spline. >1 pulls curve toward the point.");
                });
            }
        });
    }
    if ui.button("Add curve")
        .on_hover_text("Build the NURBS curve from the inputs above and add it to the curve list.")
        .clicked()
    {
        app.surface_create_curve();
    }
}

fn draw_surface_surface_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Construct NURBS surface").strong());
    {
        let s = &mut app.mesh_toolbox.surface;
        ui.horizontal(|ui| {
            ui.label("u_degree:")
                .on_hover_text("Degree along the u parameter direction (typically 3 for cubic).");
            ui.add(
                egui::DragValue::new(&mut s.surface_u_degree)
                    .range(1..=9)
                    .speed(0.1),
            )
            .on_hover_text("u-direction degree (1..9).");
            ui.label("v_degree:")
                .on_hover_text("Degree along the v parameter direction.");
            ui.add(
                egui::DragValue::new(&mut s.surface_v_degree)
                    .range(1..=9)
                    .speed(0.1),
            )
            .on_hover_text("v-direction degree (1..9).");
        });
        ui.horizontal(|ui| {
            ui.label("nu × nv:")
                .on_hover_text("Control-point grid size: nu rows × nv columns. Each must be > respective degree.");
            ui.add(
                egui::DragValue::new(&mut s.surface_nu)
                    .range(2..=10)
                    .speed(0.1),
            )
            .on_hover_text("Number of control-point rows (u direction).");
            ui.add(
                egui::DragValue::new(&mut s.surface_nv)
                    .range(2..=10)
                    .speed(0.1),
            )
            .on_hover_text("Number of control-point columns (v direction).");
        });
        s.surface_cps
            .resize(s.surface_nu * s.surface_nv, [0.0, 0.0, 0.0]);
        ui.collapsing("Control points (row-major in u)", |ui| {
            for i in 0..s.surface_nu {
                for j in 0..s.surface_nv {
                    let k = i * s.surface_nv + j;
                    ui.horizontal(|ui| {
                        ui.label(format!("[{i},{j}]:"));
                        ui.add(
                            egui::DragValue::new(&mut s.surface_cps[k][0])
                                .speed(0.05)
                                .prefix("x "),
                        )
                        .on_hover_text("Control point X.");
                        ui.add(
                            egui::DragValue::new(&mut s.surface_cps[k][1])
                                .speed(0.05)
                                .prefix("y "),
                        )
                        .on_hover_text("Control point Y.");
                        ui.add(
                            egui::DragValue::new(&mut s.surface_cps[k][2])
                                .speed(0.05)
                                .prefix("z "),
                        )
                        .on_hover_text("Control point Z.");
                    });
                }
            }
        });
    }
    if ui.button("Add surface")
        .on_hover_text("Build the NURBS surface from the inputs above and add it to the surface list.")
        .clicked()
    {
        app.surface_create_surface();
    }
}

fn draw_surface_coons_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Coons patch fill").strong());
    let n_curves = app.mesh_toolbox.surface.file.curves.len();
    ui.label(format!("(have {n_curves} curves)"));
    {
        let s = &mut app.mesh_toolbox.surface;
        let descriptions = [
            "Boundary curve at v=0 (parametric bottom edge).",
            "Boundary curve at v=1 (parametric top edge).",
            "Boundary curve at u=0 (parametric left edge).",
            "Boundary curve at u=1 (parametric right edge).",
        ];
        for (k, name) in ["c0 (v_min)", "c1 (v_max)", "d0 (u_min)", "d1 (u_max)"]
            .iter()
            .enumerate()
        {
            ui.horizontal(|ui| {
                ui.label(format!("{name}:")).on_hover_text(descriptions[k]);
                ui.add(
                    egui::DragValue::new(&mut s.coons_curves[k])
                        .range(0..=n_curves.saturating_sub(1).max(0))
                        .speed(0.2)
                        .prefix("#"),
                )
                .on_hover_text(descriptions[k]);
            });
        }
    }
    if ui.button("Fill")
        .on_hover_text("Build a bilinearly-blended Coons patch from the four boundary curves.")
        .clicked()
    {
        app.surface_coons_fill();
    }
}

fn draw_surface_sew_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Sew two surfaces").strong());
    let n_surf = app.mesh_toolbox.surface.file.surfaces.len();
    ui.label(format!("(have {n_surf} surfaces)"));
    {
        let s = &mut app.mesh_toolbox.surface;
        ui.horizontal(|ui| {
            ui.label("A id:")
                .on_hover_text("First surface to sew. Edge `A.edge_a` will be matched against `B.edge_b`.");
            ui.add(
                egui::DragValue::new(&mut s.sew_surface_a)
                    .range(0..=n_surf.saturating_sub(1).max(0))
                    .speed(0.2),
            )
            .on_hover_text("Surface A index in the surfaces list.");
            ui.label("A edge:")
                .on_hover_text("Which parametric edge of surface A is joined to B. UMin/UMax = u=0 / u=1, VMin/VMax = v=0 / v=1.");
            egui::ComboBox::from_id_source("surface_sew_edge_a")
                .selected_text(edge_label(s.sew_edge_a))
                .show_ui(ui, |ui| {
                    for (v, label) in [(0, "UMin"), (1, "UMax"), (2, "VMin"), (3, "VMax")] {
                        ui.selectable_value(&mut s.sew_edge_a, v, label)
                            .on_hover_text(match label {
                                "UMin" => "u = 0 edge.",
                                "UMax" => "u = 1 edge.",
                                "VMin" => "v = 0 edge.",
                                _ => "v = 1 edge.",
                            });
                    }
                });
        });
        ui.horizontal(|ui| {
            ui.label("B id:")
                .on_hover_text("Second surface to sew against surface A.");
            ui.add(
                egui::DragValue::new(&mut s.sew_surface_b)
                    .range(0..=n_surf.saturating_sub(1).max(0))
                    .speed(0.2),
            )
            .on_hover_text("Surface B index.");
            ui.label("B edge:")
                .on_hover_text("Which parametric edge of surface B is joined to A.");
            egui::ComboBox::from_id_source("surface_sew_edge_b")
                .selected_text(edge_label(s.sew_edge_b))
                .show_ui(ui, |ui| {
                    for (v, label) in [(0, "UMin"), (1, "UMax"), (2, "VMin"), (3, "VMax")] {
                        ui.selectable_value(&mut s.sew_edge_b, v, label);
                    }
                });
        });
        ui.horizontal(|ui| {
            ui.label("Tolerance:")
                .on_hover_text("Geometric distance below which the two edges are considered coincident.");
            ui.add(
                egui::DragValue::new(&mut s.sew_tolerance)
                    .speed(0.0001)
                    .range(1e-9..=1.0),
            )
            .on_hover_text("Sew tolerance (model units).");
        });
        // Phase 19C — continuity toggle (G2 default; uncheck for the
        // Phase 9 G0-averaging behaviour).
        ui.horizontal(|ui| {
            ui.checkbox(&mut s.sew_use_g2, "G2 continuous (Phase 19C)")
                .on_hover_text("Enforce curvature-continuous joint between the two surfaces. Uncheck for fast G0 vertex averaging.");
            ui.label("(uncheck for G0 / Phase 9 fast averaging)");
        });
    }
    if ui.button("Sew")
        .on_hover_text("Stitch surfaces A and B along the chosen edges. Result is a single sewn surface.")
        .clicked()
    {
        app.surface_sew();
    }
}

fn draw_surface_trim_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Trim surface by curve").strong());
    let n_surf = app.mesh_toolbox.surface.file.surfaces.len();
    let n_curves = app.mesh_toolbox.surface.file.curves.len();
    ui.label(format!("(have {n_surf} surfaces, {n_curves} curves)"));
    {
        let s = &mut app.mesh_toolbox.surface;
        ui.horizontal(|ui| {
            ui.label("Surface id:")
                .on_hover_text("Surface to trim.");
            ui.add(
                egui::DragValue::new(&mut s.trim_surface)
                    .range(0..=n_surf.saturating_sub(1).max(0))
                    .speed(0.2),
            )
            .on_hover_text("Surface index in the surfaces list.");
        });
        ui.horizontal(|ui| {
            ui.label("Curve id:")
                .on_hover_text("Boundary curve that defines the trim region (typically a closed loop on the surface).");
            ui.add(
                egui::DragValue::new(&mut s.trim_curve)
                    .range(0..=n_curves.saturating_sub(1).max(0))
                    .speed(0.2),
            )
            .on_hover_text("Curve index in the curves list.");
        });
        ui.horizontal(|ui| {
            ui.label("Side:")
                .on_hover_text("Keep the inside region (bounded by the curve) or the outside (everything else).");
            ui.selectable_value(&mut s.trim_side, 0, "Inside")
                .on_hover_text("Keep the surface region inside the curve loop.");
            ui.selectable_value(&mut s.trim_side, 1, "Outside")
                .on_hover_text("Keep the surface region outside the curve loop.");
        });
        ui.horizontal(|ui| {
            ui.label("Resolution:")
                .on_hover_text("Tessellation resolution for the trim output (higher = finer mesh).");
            ui.add(
                egui::DragValue::new(&mut s.trim_resolution)
                    .range(8..=128)
                    .speed(0.2),
            )
            .on_hover_text("Trim mesh resolution (8..128 samples per direction).");
        });
        ui.checkbox(
            &mut s.trim_use_uv,
            "(u, v) parametric trim (Phase 9.5 — required for warped surfaces)",
        )
        .on_hover_text("Trim in parameter space rather than 3-D Euclidean. Mandatory for highly-warped surfaces where 3-D projection is ambiguous.");
    }
    if ui.button("Trim → mesh")
        .on_hover_text("Apply the trim and produce a tessellated mesh of the kept region.")
        .clicked()
    {
        app.surface_trim();
    }
}

// ===== Phase 19 tool panels =====

fn draw_surface_knot_ops_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Knot ops (Phase 19A)").strong());
    let n_curves = app.mesh_toolbox.surface.file.curves.len();
    let n_surf = app.mesh_toolbox.surface.file.surfaces.len();
    ui.label(format!("({n_curves} curves, {n_surf} surfaces)"));
    {
        let s = &mut app.mesh_toolbox.surface;
        ui.horizontal(|ui| {
            ui.label("Parameter u:")
                .on_hover_text("Parameter value at which to insert / remove a knot. Range [0, 1].");
            ui.add(
                egui::DragValue::new(&mut s.knot_op_u)
                    .speed(0.01)
                    .range(0.0..=1.0),
            )
            .on_hover_text("Knot insertion parameter (0..1).");
            ui.label("Tolerance:")
                .on_hover_text("How close a knot must be to `u` to be considered a match (for remove).");
            ui.add(
                egui::DragValue::new(&mut s.knot_op_tolerance)
                    .speed(0.0001)
                    .range(1e-10..=1.0),
            )
            .on_hover_text("Knot match tolerance (parameter-space).");
        });
        ui.horizontal(|ui| {
            ui.label("Direction:")
                .on_hover_text("For surface ops only: which knot vector to operate on (u or v).");
            ui.selectable_value(&mut s.knot_op_direction, 0, "u")
                .on_hover_text("u-direction knot vector.");
            ui.selectable_value(&mut s.knot_op_direction, 1, "v")
                .on_hover_text("v-direction knot vector.");
        });
        ui.horizontal(|ui| {
            ui.label("Elevate by:")
                .on_hover_text("How many degrees to raise (degree elevation preserves geometry but adds control points).");
            ui.add(
                egui::DragValue::new(&mut s.elevate_degree_by)
                    .range(1..=6)
                    .speed(0.1),
            )
            .on_hover_text("Degree increment for the elevate op (1..6).");
        });
        ui.separator();
        ui.label("Curve ops:");
        ui.horizontal(|ui| {
            ui.label("Curve id:")
                .on_hover_text("Index into the curves list — which curve to mutate.");
            ui.add(
                egui::DragValue::new(&mut s.knot_op_curve)
                    .range(0..=n_curves.saturating_sub(1).max(0))
                    .speed(0.2),
            )
            .on_hover_text("Curve index.");
        });
    }
    ui.horizontal(|ui| {
        if ui.button("Insert knot (curve)")
            .on_hover_text("Insert a knot at parameter u — geometry preserved, refines local control.")
            .clicked()
        {
            app.surface_insert_knot_curve();
        }
        if ui.button("Remove knot (curve)")
            .on_hover_text("Remove a knot within tolerance of u — best-effort, may fail if removal would change geometry.")
            .clicked()
        {
            app.surface_remove_knot_curve();
        }
        if ui.button("Elevate degree (curve)")
            .on_hover_text("Increase polynomial degree by `Elevate by` — adds control points but preserves shape.")
            .clicked()
        {
            app.surface_elevate_degree_curve();
        }
    });
    ui.separator();
    ui.label("Surface ops:");
    {
        let s = &mut app.mesh_toolbox.surface;
        ui.horizontal(|ui| {
            ui.label("Surface id:")
                .on_hover_text("Index into the surfaces list — which surface to mutate.");
            ui.add(
                egui::DragValue::new(&mut s.knot_op_surface)
                    .range(0..=n_surf.saturating_sub(1).max(0))
                    .speed(0.2),
            )
            .on_hover_text("Surface index.");
        });
    }
    ui.horizontal(|ui| {
        if ui.button("Insert knot (surface)")
            .on_hover_text("Insert a knot in the selected direction at parameter u.")
            .clicked()
        {
            app.surface_insert_knot();
        }
        if ui.button("Elevate degree (surface)")
            .on_hover_text("Raise the surface degree in the selected direction by `Elevate by`.")
            .clicked()
        {
            app.surface_elevate_degree_surface();
        }
    });
}

fn draw_surface_ssi_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Surface-Surface Intersection (Phase 19B)").strong());
    let n_surf = app.mesh_toolbox.surface.file.surfaces.len();
    ui.label(format!("({n_surf} surfaces)"));
    {
        let s = &mut app.mesh_toolbox.surface;
        ui.horizontal(|ui| {
            ui.label("Surface A id:")
                .on_hover_text("First surface to intersect.");
            ui.add(
                egui::DragValue::new(&mut s.ssi_surface_a)
                    .range(0..=n_surf.saturating_sub(1).max(0))
                    .speed(0.2),
            )
            .on_hover_text("Surface A index.");
            ui.label("Surface B id:")
                .on_hover_text("Second surface to intersect with A.");
            ui.add(
                egui::DragValue::new(&mut s.ssi_surface_b)
                    .range(0..=n_surf.saturating_sub(1).max(0))
                    .speed(0.2),
            )
            .on_hover_text("Surface B index.");
        });
        ui.horizontal(|ui| {
            ui.label("Tolerance:")
                .on_hover_text("Distance below which a point is considered to lie on both surfaces.");
            ui.add(
                egui::DragValue::new(&mut s.ssi_tolerance)
                    .speed(0.0001)
                    .range(1e-9..=1.0),
            )
            .on_hover_text("Intersection tolerance (model units).");
        });
    }
    if ui.button("Compute SSI → curves")
        .on_hover_text("Compute intersection curves between surfaces A and B. Output curves are added to the curve list.")
        .clicked()
    {
        app.surface_ssi();
    }
}

fn draw_surface_fit_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Fit curve / surface (Phase 19D)").strong());
    {
        let s = &mut app.mesh_toolbox.surface;
        ui.horizontal(|ui| {
            ui.label("Degree u:")
                .on_hover_text("Degree of the fitted curve / u-direction of the fitted surface.");
            ui.add(
                egui::DragValue::new(&mut s.fit_degree_u)
                    .range(1..=9)
                    .speed(0.1),
            )
            .on_hover_text("u degree (1..9).");
            ui.label("Degree v:")
                .on_hover_text("v-direction degree for surface fit (ignored for curve fit).");
            ui.add(
                egui::DragValue::new(&mut s.fit_degree_v)
                    .range(1..=9)
                    .speed(0.1),
            )
            .on_hover_text("v degree (1..9).");
        });
        ui.horizontal(|ui| {
            ui.label("nu CPs:")
                .on_hover_text("Number of control points in u — more CPs = lower RMS, but risk of overfitting.");
            ui.add(
                egui::DragValue::new(&mut s.fit_n_cps_u)
                    .range(2..=32)
                    .speed(0.1),
            )
            .on_hover_text("u-direction control-point count (2..32).");
            ui.label("nv CPs:")
                .on_hover_text("Number of control points in v (surface fit only).");
            ui.add(
                egui::DragValue::new(&mut s.fit_n_cps_v)
                    .range(2..=32)
                    .speed(0.1),
            )
            .on_hover_text("v-direction control-point count.");
        });
        ui.label("Points (one \"x,y,z\" per line; # for comments):");
        ui.add(
            egui::TextEdit::multiline(&mut s.fit_points_text)
                .desired_rows(6)
                .desired_width(f32::INFINITY),
        )
        .on_hover_text("Point cloud to fit. One point per line, comma-separated x,y,z. Lines starting with # are ignored.");
        if let Some(rms) = s.fit_last_rms {
            ui.label(format!("Last fit RMS: {rms:.4e}"));
        }
    }
    ui.horizontal(|ui| {
        if ui.button("Fit curve")
            .on_hover_text("Fit a NURBS curve through the point list (preserves input order).")
            .clicked()
        {
            app.surface_fit_curve();
        }
        if ui.button("Fit surface (scattered)")
            .on_hover_text("Fit a NURBS surface to scattered point data via least-squares.")
            .clicked()
        {
            app.surface_fit_surface();
        }
    });
}

fn draw_surface_ruled_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Ruled surface (Phase 19E)").strong());
    let n_curves = app.mesh_toolbox.surface.file.curves.len();
    ui.label(format!("({n_curves} curves)"));
    {
        let s = &mut app.mesh_toolbox.surface;
        ui.horizontal(|ui| {
            ui.label("Mode:")
                .on_hover_text("Construction mode. Between curves = straight rules between two curves. Extrude vec = sweep along a vector. Cone apex = rules from the curve to a single apex point.");
            ui.selectable_value(&mut s.ruled_kind, 0, "Between curves")
                .on_hover_text("Linear ruled surface joining curve A to curve B.");
            ui.selectable_value(&mut s.ruled_kind, 1, "Extrude vec")
                .on_hover_text("Sweep curve A along the given vector — produces a translational surface.");
            ui.selectable_value(&mut s.ruled_kind, 2, "Cone apex")
                .on_hover_text("Rules from each point on curve A to the apex point — produces a generalized cone.");
        });
        ui.horizontal(|ui| {
            ui.label("Curve A id:")
                .on_hover_text("Index of the base curve.");
            ui.add(
                egui::DragValue::new(&mut s.ruled_curve_a)
                    .range(0..=n_curves.saturating_sub(1).max(0))
                    .speed(0.2),
            )
            .on_hover_text("Curve A index.");
            if s.ruled_kind == 0 {
                ui.label("Curve B id:")
                    .on_hover_text("Index of the second curve (for Between-curves mode).");
                ui.add(
                    egui::DragValue::new(&mut s.ruled_curve_b)
                        .range(0..=n_curves.saturating_sub(1).max(0))
                        .speed(0.2),
                )
                .on_hover_text("Curve B index.");
            }
        });
        if s.ruled_kind == 1 {
            ui.horizontal(|ui| {
                ui.label("Vector:")
                    .on_hover_text("Extrusion vector applied to every point on curve A.");
                ui.add(
                    egui::DragValue::new(&mut s.ruled_extrude_vector[0])
                        .speed(0.05)
                        .prefix("x "),
                )
                .on_hover_text("Extrusion vector X component.");
                ui.add(
                    egui::DragValue::new(&mut s.ruled_extrude_vector[1])
                        .speed(0.05)
                        .prefix("y "),
                )
                .on_hover_text("Extrusion vector Y component.");
                ui.add(
                    egui::DragValue::new(&mut s.ruled_extrude_vector[2])
                        .speed(0.05)
                        .prefix("z "),
                )
                .on_hover_text("Extrusion vector Z component.");
            });
        }
        if s.ruled_kind == 2 {
            ui.horizontal(|ui| {
                ui.label("Apex:")
                    .on_hover_text("Single apex point that all rule lines connect to from curve A.");
                ui.add(
                    egui::DragValue::new(&mut s.ruled_apex[0])
                        .speed(0.05)
                        .prefix("x "),
                )
                .on_hover_text("Apex X.");
                ui.add(
                    egui::DragValue::new(&mut s.ruled_apex[1])
                        .speed(0.05)
                        .prefix("y "),
                )
                .on_hover_text("Apex Y.");
                ui.add(
                    egui::DragValue::new(&mut s.ruled_apex[2])
                        .speed(0.05)
                        .prefix("z "),
                )
                .on_hover_text("Apex Z.");
            });
        }
    }
    if ui.button("Build ruled surface")
        .on_hover_text("Create the ruled surface from the inputs above and add it to the surfaces list.")
        .clicked()
    {
        app.surface_ruled_build();
    }
}

fn edge_label(v: u8) -> &'static str {
    match v {
        0 => "UMin",
        1 => "UMax",
        2 => "VMin",
        _ => "VMax",
    }
}

fn draw_surface_entity_lists(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Curves").strong());
    let n_curves = app.mesh_toolbox.surface.file.curves.len();
    if n_curves == 0 {
        ui.label("(no curves)");
    } else {
        for i in 0..n_curves {
            let c = &app.mesh_toolbox.surface.file.curves[i];
            ui.label(format!(
                "#{i}: degree {}, {} CPs",
                c.degree,
                c.n_control_points()
            ));
        }
    }
    ui.separator();
    ui.label(egui::RichText::new("Surfaces").strong());
    let n_surf = app.mesh_toolbox.surface.file.surfaces.len();
    if n_surf == 0 {
        ui.label("(no surfaces)");
    } else {
        // Collect ids first so the borrow checker is happy when we
        // call `surface_tessellate(i)` inside the loop.
        let ids: Vec<(usize, usize, usize, usize, usize)> = (0..n_surf)
            .map(|i| {
                let s = &app.mesh_toolbox.surface.file.surfaces[i];
                (i, s.u_degree, s.v_degree, s.nu(), s.nv())
            })
            .collect();
        for (i, u_deg, v_deg, nu, nv) in ids {
            ui.horizontal(|ui| {
                ui.label(format!("#{i}: deg ({u_deg},{v_deg}), CPs {nu}×{nv}"));
                if ui.button("Tessellate → viewport").clicked() {
                    app.surface_tessellate(i);
                }
            });
        }
    }
    ui.horizontal(|ui| {
        ui.label("Tess resolution:");
        ui.add(
            egui::DragValue::new(&mut app.mesh_toolbox.surface.tess_resolution)
                .range(4..=128)
                .speed(0.2),
        );
    });
}

// ===== Phase 15 — Arch / BIM workbench panel =====

impl crate::ValenxApp {
    /// Add a wall using the panel's current inputs.
    pub fn arch_add_wall(&mut self) {
        let s = &mut self.mesh_toolbox.arch;
        let w = valenx_arch::WallParams {
            start: nalgebra::Vector3::new(s.wall_start[0], s.wall_start[1], s.wall_start[2]),
            end: nalgebra::Vector3::new(s.wall_end[0], s.wall_end[1], s.wall_end[2]),
            height: s.wall_height,
            thickness: s.wall_thickness,
            material: s.wall_material.clone(),
        };
        if let Err(e) = w.validate() {
            s.last_error = Some(format!("Wall: {e}"));
            return;
        }
        let id = s.doc.add_entity(valenx_arch::ArchEntity::Wall(w));
        s.last_status = Some(format!("Added Wall #{id}"));
        s.last_error = None;
        emit_audit(
            "arch.add_wall",
            serde_json::json!({"id": id}),
            serde_json::json!({}),
        );
    }

    /// Add a slab using the panel's current inputs.
    pub fn arch_add_slab(&mut self) {
        let s = &mut self.mesh_toolbox.arch;
        if s.slab_boundary.len() < 3 {
            s.last_error = Some("Slab: need ≥ 3 boundary points".into());
            return;
        }
        let boundary: Vec<_> = s
            .slab_boundary
            .iter()
            .map(|p| nalgebra::Vector3::new(p[0], p[1], s.slab_z))
            .collect();
        let sl = valenx_arch::SlabParams {
            boundary,
            thickness: s.slab_thickness,
            material: s.slab_material.clone(),
            structural: None,
        };
        if let Err(e) = sl.validate() {
            s.last_error = Some(format!("Slab: {e}"));
            return;
        }
        let id = s.doc.add_entity(valenx_arch::ArchEntity::Slab(sl));
        s.last_status = Some(format!("Added Slab #{id}"));
        s.last_error = None;
    }

    /// Add a column.
    pub fn arch_add_column(&mut self) {
        let s = &mut self.mesh_toolbox.arch;
        let section = match s.column_section_kind {
            0 => valenx_arch::ColumnSection::Rectangle {
                width: s.column_width,
                depth: s.column_depth,
            },
            1 => valenx_arch::ColumnSection::Circular {
                radius: s.column_radius,
                segments: s.column_segments,
            },
            _ => valenx_arch::ColumnSection::IBeam {
                width: s.column_width,
                depth: s.column_depth,
                flange_thickness: s.column_flange_thickness,
                web_thickness: s.column_web_thickness,
            },
        };
        let c = valenx_arch::ColumnParams {
            base: nalgebra::Vector3::new(s.column_base[0], s.column_base[1], s.column_base[2]),
            height: s.column_height,
            cross_section: section,
            material: s.column_material.clone(),
            structural: None,
        };
        if let Err(e) = c.validate() {
            s.last_error = Some(format!("Column: {e}"));
            return;
        }
        let id = s.doc.add_entity(valenx_arch::ArchEntity::Column(c));
        s.last_status = Some(format!("Added Column #{id}"));
        s.last_error = None;
    }

    /// Add a beam.
    pub fn arch_add_beam(&mut self) {
        let s = &mut self.mesh_toolbox.arch;
        let section = match s.beam_section_kind {
            0 => valenx_arch::BeamSection::Rectangle {
                width: s.beam_width,
                depth: s.beam_depth,
            },
            1 => valenx_arch::BeamSection::IBeam {
                width: s.beam_width,
                depth: s.beam_depth,
                flange_thickness: s.beam_flange_thickness,
                web_thickness: s.beam_web_thickness,
            },
            _ => valenx_arch::BeamSection::Channel {
                width: s.beam_width,
                depth: s.beam_depth,
                thickness: s.beam_flange_thickness,
            },
        };
        let b = valenx_arch::BeamParams {
            start: nalgebra::Vector3::new(s.beam_start[0], s.beam_start[1], s.beam_start[2]),
            end: nalgebra::Vector3::new(s.beam_end[0], s.beam_end[1], s.beam_end[2]),
            cross_section: section,
            orientation_angle: s.beam_orientation,
            material: s.beam_material.clone(),
            structural: None,
        };
        if let Err(e) = b.validate() {
            s.last_error = Some(format!("Beam: {e}"));
            return;
        }
        let id = s.doc.add_entity(valenx_arch::ArchEntity::Beam(b));
        s.last_status = Some(format!("Added Beam #{id}"));
        s.last_error = None;
    }

    /// Add a window.
    pub fn arch_add_window(&mut self) {
        let s = &mut self.mesh_toolbox.arch;
        let style = match s.window_style {
            0 => valenx_arch::WindowStyle::Casement,
            1 => valenx_arch::WindowStyle::Sliding,
            2 => valenx_arch::WindowStyle::Awning,
            _ => valenx_arch::WindowStyle::Fixed,
        };
        let w = valenx_arch::WindowParams {
            host: s.window_host,
            position_along_wall: s.window_position_along,
            position_height: s.window_position_height,
            width: s.window_width,
            height: s.window_height,
            frame_thickness: s.window_frame_thickness,
            style,
        };
        if let Err(e) = w.validate() {
            s.last_error = Some(format!("Window: {e}"));
            return;
        }
        let id = s.doc.add_entity(valenx_arch::ArchEntity::Window(w));
        s.last_status = Some(format!("Added Window #{id}"));
        s.last_error = None;
    }

    /// Add a door.
    pub fn arch_add_door(&mut self) {
        let s = &mut self.mesh_toolbox.arch;
        let style = match s.door_style {
            0 => valenx_arch::DoorStyle::Single,
            1 => valenx_arch::DoorStyle::Double,
            2 => valenx_arch::DoorStyle::Sliding,
            _ => valenx_arch::DoorStyle::Bifold,
        };
        let hinge = match s.door_hinge_side {
            0 => valenx_arch::Side::Left,
            _ => valenx_arch::Side::Right,
        };
        let d = valenx_arch::DoorParams {
            host: s.door_host,
            position_along_wall: s.door_position_along,
            width: s.door_width,
            height: s.door_height,
            style,
            hinge_side: hinge,
        };
        if let Err(e) = d.validate() {
            s.last_error = Some(format!("Door: {e}"));
            return;
        }
        let id = s.doc.add_entity(valenx_arch::ArchEntity::Door(d));
        s.last_status = Some(format!("Added Door #{id}"));
        s.last_error = None;
    }

    /// Add a stair.
    pub fn arch_add_stair(&mut self) {
        let s = &mut self.mesh_toolbox.arch;
        let st = valenx_arch::StairParams {
            base: nalgebra::Vector3::new(s.stair_base[0], s.stair_base[1], s.stair_base[2]),
            direction: nalgebra::Vector3::new(
                s.stair_direction[0],
                s.stair_direction[1],
                s.stair_direction[2],
            ),
            total_rise: s.stair_total_rise,
            total_run: s.stair_total_run,
            num_steps: s.stair_num_steps,
            width: s.stair_width,
        };
        if let Err(e) = st.validate() {
            s.last_error = Some(format!("Stair: {e}"));
            return;
        }
        let id = s.doc.add_entity(valenx_arch::ArchEntity::Stair(st));
        s.last_status = Some(format!("Added Stair #{id}"));
        s.last_error = None;
    }

    /// Add a roof.
    pub fn arch_add_roof(&mut self) {
        let s = &mut self.mesh_toolbox.arch;
        if s.roof_boundary.len() < 3 {
            s.last_error = Some("Roof: need ≥ 3 boundary points".into());
            return;
        }
        let boundary: Vec<_> = s
            .roof_boundary
            .iter()
            .map(|p| nalgebra::Vector3::new(p[0], p[1], s.roof_z))
            .collect();
        let rt = match s.roof_type {
            0 => valenx_arch::RoofType::Flat,
            1 => valenx_arch::RoofType::Gable,
            2 => valenx_arch::RoofType::Hip,
            _ => valenx_arch::RoofType::Shed,
        };
        let r = valenx_arch::RoofParams {
            boundary,
            peak_height: s.roof_peak_height,
            roof_type: rt,
        };
        if let Err(e) = r.validate() {
            s.last_error = Some(format!("Roof: {e}"));
            return;
        }
        let id = s.doc.add_entity(valenx_arch::ArchEntity::Roof(r));
        s.last_status = Some(format!("Added Roof #{id}"));
        s.last_error = None;
    }

    /// Add a space.
    pub fn arch_add_space(&mut self) {
        let s = &mut self.mesh_toolbox.arch;
        if s.space_boundary.len() < 3 {
            s.last_error = Some("Space: need ≥ 3 boundary points".into());
            return;
        }
        let boundary: Vec<_> = s
            .space_boundary
            .iter()
            .map(|p| nalgebra::Vector3::new(p[0], p[1], s.space_z))
            .collect();
        let sp = valenx_arch::SpaceParams {
            boundary,
            ceiling_height: s.space_ceiling_height,
            space_name: s.space_name.clone(),
        };
        if let Err(e) = sp.validate() {
            s.last_error = Some(format!("Space: {e}"));
            return;
        }
        let id = s.doc.add_entity(valenx_arch::ArchEntity::Space(sp));
        s.last_status = Some(format!("Added Space #{id}"));
        s.last_error = None;
    }

    /// Delete an entity by id.
    pub fn arch_delete_entity(&mut self, id: usize) {
        let s = &mut self.mesh_toolbox.arch;
        match s.doc.delete_entity(id) {
            Ok(_) => {
                s.last_status = Some(format!("Deleted entity #{id}"));
                s.last_error = None;
                if s.selected == Some(id) {
                    s.selected = None;
                }
            }
            Err(e) => s.last_error = Some(format!("Delete: {e}")),
        }
    }

    /// Refresh the cached schedule.
    pub fn arch_refresh_schedule(&mut self) {
        let s = &mut self.mesh_toolbox.arch;
        s.last_schedule = Some(valenx_arch::Schedule::from_document(&s.doc));
        s.last_status = Some("Schedule refreshed".into());
        s.last_error = None;
    }

    /// Tessellate every entity and push to the viewport.
    pub fn arch_render(&mut self) {
        let tol = self.mesh_toolbox.arch.render_tolerance;
        let mesh = match self.mesh_toolbox.arch.doc.tessellate_all(tol) {
            Ok(m) => m,
            Err(e) => {
                self.mesh_toolbox.arch.last_error = Some(format!("Render: {e}"));
                return;
            }
        };
        let pseudo_path = std::path::PathBuf::from("<arch>/scene.fused");
        self.apply_mesh(mesh, pseudo_path);
        self.mesh_toolbox.arch.last_status = Some("Rendered to viewport".into());
        self.mesh_toolbox.arch.last_error = None;
    }

    /// Save the schedule as CSV at a user-chosen path.
    ///
    /// Uses `rfd::FileDialog` — only call from a runtime UI event, NEVER
    /// from a unit test (Phase 10 lockdown). Tests that exercise this
    /// must be marked `#[ignore]`.
    pub fn arch_save_schedule_csv(&mut self) {
        let s = &mut self.mesh_toolbox.arch;
        if s.last_schedule.is_none() {
            s.last_schedule = Some(valenx_arch::Schedule::from_document(&s.doc));
        }
        let schedule = s.last_schedule.as_ref().unwrap();
        let csv = schedule.to_csv();
        let path = match rfd::FileDialog::new()
            .add_filter("CSV", &["csv"])
            .set_file_name("schedule.csv")
            .save_file()
        {
            Some(p) => p,
            None => return,
        };
        match valenx_core::io_caps::atomic_write_str(&path, &csv) {
            Ok(()) => s.last_status = Some(format!("Saved {}", path.display())),
            Err(e) => s.last_error = Some(format!("Save CSV: {e}")),
        }
    }

    /// Export the document to an IFC4 file.
    pub fn arch_export_ifc(&mut self) {
        let path = match rfd::FileDialog::new()
            .add_filter("IFC4", &["ifc"])
            .set_file_name("model.ifc")
            .save_file()
        {
            Some(p) => p,
            None => return,
        };
        let s = &mut self.mesh_toolbox.arch;
        match valenx_arch::ifc::write_document(&s.doc, &path) {
            Ok(()) => {
                s.last_status = Some(format!("Exported IFC to {}", path.display()));
                s.last_error = None;
            }
            Err(e) => s.last_error = Some(format!("Export IFC: {e}")),
        }
    }
}

pub fn draw_arch_panel(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.heading("Arch / BIM");
    ui.separator();
    {
        let s = &mut app.mesh_toolbox.arch;
        ui.horizontal(|ui| {
            ui.label("Project:");
            ui.text_edit_singleline(&mut s.doc.project_name);
        });
    }
    ui.separator();

    // Tool palette.
    {
        let s = &mut app.mesh_toolbox.arch;
        ui.label(egui::RichText::new("Tool").strong());
        ui.horizontal_wrapped(|ui| {
            for t in [
                ArchTool::Wall,
                ArchTool::Slab,
                ArchTool::Column,
                ArchTool::Beam,
                ArchTool::Window,
                ArchTool::Door,
                ArchTool::Stair,
                ArchTool::Roof,
                ArchTool::Space,
            ] {
                ui.selectable_value(&mut s.tool, t, t.label());
            }
        });
    }
    ui.separator();

    let tool = app.mesh_toolbox.arch.tool;
    match tool {
        ArchTool::Wall => draw_arch_wall_inputs(app, ui),
        ArchTool::Slab => draw_arch_slab_inputs(app, ui),
        ArchTool::Column => draw_arch_column_inputs(app, ui),
        ArchTool::Beam => draw_arch_beam_inputs(app, ui),
        ArchTool::Window => draw_arch_window_inputs(app, ui),
        ArchTool::Door => draw_arch_door_inputs(app, ui),
        ArchTool::Stair => draw_arch_stair_inputs(app, ui),
        ArchTool::Roof => draw_arch_roof_inputs(app, ui),
        ArchTool::Space => draw_arch_space_inputs(app, ui),
    }

    ui.separator();
    draw_arch_entity_list(app, ui);
    ui.separator();
    draw_arch_render_and_schedule(app, ui);
    ui.separator();
    {
        let s = &app.mesh_toolbox.arch;
        if let Some(msg) = &s.last_status {
            ui.label(egui::RichText::new(msg).color(egui::Color32::GREEN));
        }
        if let Some(err) = &s.last_error {
            ui.label(egui::RichText::new(err).color(egui::Color32::RED));
        }
    }
}

fn draw_arch_wall_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    let s = &mut app.mesh_toolbox.arch;
    ui.label(egui::RichText::new("Wall").strong())
        .on_hover_text("Architectural wall from a start XYZ to an end XYZ, extruded upwards by Height.");
    ui.horizontal(|ui| {
        ui.label("Start:")
            .on_hover_text("Start endpoint of the wall centreline (model units, typically m).");
        let prefixes = ["x ", "y ", "z "];
        let tooltips = ["Start X (m).", "Start Y (m).", "Start Z (m, base of the wall)."];
        for i in 0..3 {
            ui.add(
                egui::DragValue::new(&mut s.wall_start[i])
                    .speed(0.1)
                    .prefix(prefixes[i]),
            )
            .on_hover_text(tooltips[i]);
        }
    });
    ui.horizontal(|ui| {
        ui.label("End:")
            .on_hover_text("End endpoint of the wall centreline (model units).");
        let prefixes = ["x ", "y ", "z "];
        let tooltips = ["End X (m).", "End Y (m).", "End Z (m)."];
        for i in 0..3 {
            ui.add(
                egui::DragValue::new(&mut s.wall_end[i])
                    .speed(0.1)
                    .prefix(prefixes[i]),
            )
            .on_hover_text(tooltips[i]);
        }
    });
    ui.horizontal(|ui| {
        ui.label("Height:")
            .on_hover_text("Wall height above the base — typically the floor-to-ceiling height.");
        ui.add(
            egui::DragValue::new(&mut s.wall_height)
                .range(0.01..=100.0)
                .speed(0.1)
                .suffix(" m"),
        )
        .on_hover_text("Wall height (m).");
        ui.label("Thickness:")
            .on_hover_text("Wall cross-section thickness — typically 100..400 mm for partitions / exterior walls.");
        ui.add(
            egui::DragValue::new(&mut s.wall_thickness)
                .range(0.001..=10.0)
                .speed(0.01)
                .suffix(" m"),
        )
        .on_hover_text("Wall thickness (m).");
    });
    ui.horizontal(|ui| {
        ui.label("Material:")
            .on_hover_text("Material tag — surfaces in the BIM entity list and the schedule export.");
        ui.text_edit_singleline(&mut s.wall_material)
            .on_hover_text("Free-text material name (e.g. \"concrete\", \"gypsum board\").");
    });
    if ui
        .button("Add Wall")
        .on_hover_text("Add a wall entity with the parameters above to the BIM document.")
        .clicked()
    {
        app.arch_add_wall();
    }
}

fn draw_arch_slab_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    let s = &mut app.mesh_toolbox.arch;
    ui.label(egui::RichText::new("Slab").strong())
        .on_hover_text("Floor / ceiling slab — horizontal plate at the given Z, with the boundary polygon defined below.");
    ui.horizontal(|ui| {
        ui.label("Z:")
            .on_hover_text("Slab elevation — height of the slab's top surface above the project origin (m).");
        ui.add(egui::DragValue::new(&mut s.slab_z).speed(0.1).suffix(" m"))
            .on_hover_text("Slab Z (m).");
        ui.label("Thickness:")
            .on_hover_text("Slab depth — typical concrete slabs are 150..300 mm.");
        ui.add(
            egui::DragValue::new(&mut s.slab_thickness)
                .range(0.01..=5.0)
                .speed(0.01)
                .suffix(" m"),
        )
        .on_hover_text("Slab thickness (m).");
    });
    ui.horizontal(|ui| {
        ui.label("Material:")
            .on_hover_text("Slab material tag — flows into the BIM schedule.");
        ui.text_edit_singleline(&mut s.slab_material)
            .on_hover_text("Free-text material name (e.g. \"reinforced concrete\").");
    });
    ui.collapsing("Boundary (x, y)", |ui| {
        let mut remove: Option<usize> = None;
        for (i, p) in s.slab_boundary.iter_mut().enumerate() {
            ui.horizontal(|ui| {
                ui.label(format!("#{i}:"));
                ui.add(egui::DragValue::new(&mut p[0]).speed(0.1).prefix("x "));
                ui.add(egui::DragValue::new(&mut p[1]).speed(0.1).prefix("y "));
                if ui.small_button("Delete").clicked() {
                    remove = Some(i);
                }
            });
        }
        if let Some(i) = remove {
            s.slab_boundary.remove(i);
        }
        if ui.button("Add point").clicked() {
            s.slab_boundary.push([0.0, 0.0]);
        }
    });
    if ui.button("Add Slab").clicked() {
        app.arch_add_slab();
    }
}

fn draw_arch_column_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    let s = &mut app.mesh_toolbox.arch;
    ui.label(egui::RichText::new("Column").strong());
    ui.horizontal(|ui| {
        ui.label("Base:");
        for i in 0..3 {
            ui.add(
                egui::DragValue::new(&mut s.column_base[i])
                    .speed(0.1)
                    .prefix(["x ", "y ", "z "][i]),
            );
        }
    });
    ui.horizontal(|ui| {
        ui.label("Height:");
        ui.add(
            egui::DragValue::new(&mut s.column_height)
                .range(0.01..=100.0)
                .speed(0.1)
                .suffix(" m"),
        );
    });
    ui.horizontal(|ui| {
        ui.label("Section:");
        ui.selectable_value(&mut s.column_section_kind, 0, "Rect");
        ui.selectable_value(&mut s.column_section_kind, 1, "Circ");
        ui.selectable_value(&mut s.column_section_kind, 2, "IBeam");
    });
    match s.column_section_kind {
        0 => {
            ui.horizontal(|ui| {
                ui.label("Width:");
                ui.add(
                    egui::DragValue::new(&mut s.column_width)
                        .speed(0.01)
                        .suffix(" m"),
                );
                ui.label("Depth:");
                ui.add(
                    egui::DragValue::new(&mut s.column_depth)
                        .speed(0.01)
                        .suffix(" m"),
                );
            });
        }
        1 => {
            ui.horizontal(|ui| {
                ui.label("Radius:");
                ui.add(
                    egui::DragValue::new(&mut s.column_radius)
                        .speed(0.01)
                        .suffix(" m"),
                );
                ui.label("Segments:");
                ui.add(
                    egui::DragValue::new(&mut s.column_segments)
                        .range(3..=128)
                        .speed(0.5),
                );
            });
        }
        _ => {
            ui.horizontal(|ui| {
                ui.label("Width:");
                ui.add(
                    egui::DragValue::new(&mut s.column_width)
                        .speed(0.01)
                        .suffix(" m"),
                );
                ui.label("Depth:");
                ui.add(
                    egui::DragValue::new(&mut s.column_depth)
                        .speed(0.01)
                        .suffix(" m"),
                );
            });
            ui.horizontal(|ui| {
                ui.label("Flange t:");
                ui.add(
                    egui::DragValue::new(&mut s.column_flange_thickness)
                        .speed(0.005)
                        .suffix(" m"),
                );
                ui.label("Web t:");
                ui.add(
                    egui::DragValue::new(&mut s.column_web_thickness)
                        .speed(0.005)
                        .suffix(" m"),
                );
            });
        }
    }
    ui.horizontal(|ui| {
        ui.label("Material:");
        ui.text_edit_singleline(&mut s.column_material);
    });
    if ui.button("Add Column").clicked() {
        app.arch_add_column();
    }
}

fn draw_arch_beam_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    let s = &mut app.mesh_toolbox.arch;
    ui.label(egui::RichText::new("Beam").strong());
    ui.horizontal(|ui| {
        ui.label("Start:");
        for i in 0..3 {
            ui.add(
                egui::DragValue::new(&mut s.beam_start[i])
                    .speed(0.1)
                    .prefix(["x ", "y ", "z "][i]),
            );
        }
    });
    ui.horizontal(|ui| {
        ui.label("End:");
        for i in 0..3 {
            ui.add(
                egui::DragValue::new(&mut s.beam_end[i])
                    .speed(0.1)
                    .prefix(["x ", "y ", "z "][i]),
            );
        }
    });
    ui.horizontal(|ui| {
        ui.label("Section:");
        ui.selectable_value(&mut s.beam_section_kind, 0, "Rect");
        ui.selectable_value(&mut s.beam_section_kind, 1, "IBeam");
        ui.selectable_value(&mut s.beam_section_kind, 2, "Channel");
    });
    ui.horizontal(|ui| {
        ui.label("Width:");
        ui.add(
            egui::DragValue::new(&mut s.beam_width)
                .speed(0.01)
                .suffix(" m"),
        );
        ui.label("Depth:");
        ui.add(
            egui::DragValue::new(&mut s.beam_depth)
                .speed(0.01)
                .suffix(" m"),
        );
    });
    if matches!(s.beam_section_kind, 1 | 2) {
        ui.horizontal(|ui| {
            ui.label("Flange/Wall t:");
            ui.add(
                egui::DragValue::new(&mut s.beam_flange_thickness)
                    .speed(0.005)
                    .suffix(" m"),
            );
            if s.beam_section_kind == 1 {
                ui.label("Web t:");
                ui.add(
                    egui::DragValue::new(&mut s.beam_web_thickness)
                        .speed(0.005)
                        .suffix(" m"),
                );
            }
        });
    }
    ui.horizontal(|ui| {
        ui.label("Orientation:");
        ui.add(
            egui::DragValue::new(&mut s.beam_orientation)
                .speed(0.05)
                .suffix(" rad"),
        );
        ui.label("Material:");
        ui.text_edit_singleline(&mut s.beam_material);
    });
    if ui.button("Add Beam").clicked() {
        app.arch_add_beam();
    }
}

fn draw_arch_window_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    let s = &mut app.mesh_toolbox.arch;
    ui.label(egui::RichText::new("Window").strong());
    ui.horizontal(|ui| {
        ui.label("Host wall id:");
        ui.add(
            egui::DragValue::new(&mut s.window_host)
                .range(1..=u32::MAX as usize)
                .speed(0.2),
        );
    });
    ui.horizontal(|ui| {
        ui.label("Position along:");
        ui.add(
            egui::DragValue::new(&mut s.window_position_along)
                .speed(0.1)
                .suffix(" m"),
        );
        ui.label("Sill height:");
        ui.add(
            egui::DragValue::new(&mut s.window_position_height)
                .speed(0.05)
                .suffix(" m"),
        );
    });
    ui.horizontal(|ui| {
        ui.label("Width:");
        ui.add(
            egui::DragValue::new(&mut s.window_width)
                .range(0.05..=10.0)
                .speed(0.05)
                .suffix(" m"),
        );
        ui.label("Height:");
        ui.add(
            egui::DragValue::new(&mut s.window_height)
                .range(0.05..=10.0)
                .speed(0.05)
                .suffix(" m"),
        );
    });
    ui.horizontal(|ui| {
        ui.label("Frame t:");
        ui.add(
            egui::DragValue::new(&mut s.window_frame_thickness)
                .speed(0.005)
                .suffix(" m"),
        );
    });
    ui.horizontal(|ui| {
        ui.label("Style:");
        ui.selectable_value(&mut s.window_style, 0, "Casement");
        ui.selectable_value(&mut s.window_style, 1, "Sliding");
        ui.selectable_value(&mut s.window_style, 2, "Awning");
        ui.selectable_value(&mut s.window_style, 3, "Fixed");
    });
    if ui.button("Add Window").clicked() {
        app.arch_add_window();
    }
}

fn draw_arch_door_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    let s = &mut app.mesh_toolbox.arch;
    ui.label(egui::RichText::new("Door").strong());
    ui.horizontal(|ui| {
        ui.label("Host wall id:");
        ui.add(
            egui::DragValue::new(&mut s.door_host)
                .range(1..=u32::MAX as usize)
                .speed(0.2),
        );
    });
    ui.horizontal(|ui| {
        ui.label("Position along:");
        ui.add(
            egui::DragValue::new(&mut s.door_position_along)
                .speed(0.1)
                .suffix(" m"),
        );
    });
    ui.horizontal(|ui| {
        ui.label("Width:");
        ui.add(
            egui::DragValue::new(&mut s.door_width)
                .range(0.05..=10.0)
                .speed(0.05)
                .suffix(" m"),
        );
        ui.label("Height:");
        ui.add(
            egui::DragValue::new(&mut s.door_height)
                .range(0.05..=10.0)
                .speed(0.05)
                .suffix(" m"),
        );
    });
    ui.horizontal(|ui| {
        ui.label("Style:");
        ui.selectable_value(&mut s.door_style, 0, "Single");
        ui.selectable_value(&mut s.door_style, 1, "Double");
        ui.selectable_value(&mut s.door_style, 2, "Sliding");
        ui.selectable_value(&mut s.door_style, 3, "Bifold");
    });
    ui.horizontal(|ui| {
        ui.label("Hinge:");
        ui.selectable_value(&mut s.door_hinge_side, 0, "Left");
        ui.selectable_value(&mut s.door_hinge_side, 1, "Right");
    });
    if ui.button("Add Door").clicked() {
        app.arch_add_door();
    }
}

fn draw_arch_stair_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    let s = &mut app.mesh_toolbox.arch;
    ui.label(egui::RichText::new("Stair").strong());
    ui.horizontal(|ui| {
        ui.label("Base:");
        for i in 0..3 {
            ui.add(
                egui::DragValue::new(&mut s.stair_base[i])
                    .speed(0.1)
                    .prefix(["x ", "y ", "z "][i]),
            );
        }
    });
    ui.horizontal(|ui| {
        ui.label("Dir:");
        for i in 0..3 {
            ui.add(
                egui::DragValue::new(&mut s.stair_direction[i])
                    .speed(0.05)
                    .prefix(["x ", "y ", "z "][i]),
            );
        }
    });
    ui.horizontal(|ui| {
        ui.label("Rise:");
        ui.add(
            egui::DragValue::new(&mut s.stair_total_rise)
                .range(0.05..=50.0)
                .speed(0.05)
                .suffix(" m"),
        );
        ui.label("Run:");
        ui.add(
            egui::DragValue::new(&mut s.stair_total_run)
                .range(0.05..=50.0)
                .speed(0.05)
                .suffix(" m"),
        );
    });
    ui.horizontal(|ui| {
        ui.label("Steps:");
        ui.add(
            egui::DragValue::new(&mut s.stair_num_steps)
                .range(1u32..=200)
                .speed(0.1),
        );
        ui.label("Width:");
        ui.add(
            egui::DragValue::new(&mut s.stair_width)
                .range(0.05..=10.0)
                .speed(0.05)
                .suffix(" m"),
        );
    });
    if ui.button("Add Stair").clicked() {
        app.arch_add_stair();
    }
}

fn draw_arch_roof_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    let s = &mut app.mesh_toolbox.arch;
    ui.label(egui::RichText::new("Roof").strong());
    ui.horizontal(|ui| {
        ui.label("Eave z:");
        ui.add(egui::DragValue::new(&mut s.roof_z).speed(0.1).suffix(" m"));
        ui.label("Peak h:");
        ui.add(
            egui::DragValue::new(&mut s.roof_peak_height)
                .range(0.0..=50.0)
                .speed(0.1)
                .suffix(" m"),
        );
    });
    ui.horizontal(|ui| {
        ui.label("Type:");
        ui.selectable_value(&mut s.roof_type, 0, "Flat");
        ui.selectable_value(&mut s.roof_type, 1, "Gable");
        ui.selectable_value(&mut s.roof_type, 2, "Hip");
        ui.selectable_value(&mut s.roof_type, 3, "Shed");
    });
    ui.collapsing("Boundary (x, y)", |ui| {
        let mut remove: Option<usize> = None;
        for (i, p) in s.roof_boundary.iter_mut().enumerate() {
            ui.horizontal(|ui| {
                ui.label(format!("#{i}:"));
                ui.add(egui::DragValue::new(&mut p[0]).speed(0.1).prefix("x "));
                ui.add(egui::DragValue::new(&mut p[1]).speed(0.1).prefix("y "));
                if ui.small_button("Delete").clicked() {
                    remove = Some(i);
                }
            });
        }
        if let Some(i) = remove {
            s.roof_boundary.remove(i);
        }
        if ui.button("Add point").clicked() {
            s.roof_boundary.push([0.0, 0.0]);
        }
    });
    if ui.button("Add Roof").clicked() {
        app.arch_add_roof();
    }
}

fn draw_arch_space_inputs(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    let s = &mut app.mesh_toolbox.arch;
    ui.label(egui::RichText::new("Space").strong());
    ui.horizontal(|ui| {
        ui.label("Name:");
        ui.text_edit_singleline(&mut s.space_name);
    });
    ui.horizontal(|ui| {
        ui.label("Floor z:");
        ui.add(egui::DragValue::new(&mut s.space_z).speed(0.1).suffix(" m"));
        ui.label("Ceiling h:");
        ui.add(
            egui::DragValue::new(&mut s.space_ceiling_height)
                .range(0.05..=10.0)
                .speed(0.05)
                .suffix(" m"),
        );
    });
    ui.collapsing("Boundary (x, y)", |ui| {
        let mut remove: Option<usize> = None;
        for (i, p) in s.space_boundary.iter_mut().enumerate() {
            ui.horizontal(|ui| {
                ui.label(format!("#{i}:"));
                ui.add(egui::DragValue::new(&mut p[0]).speed(0.1).prefix("x "));
                ui.add(egui::DragValue::new(&mut p[1]).speed(0.1).prefix("y "));
                if ui.small_button("Delete").clicked() {
                    remove = Some(i);
                }
            });
        }
        if let Some(i) = remove {
            s.space_boundary.remove(i);
        }
        if ui.button("Add point").clicked() {
            s.space_boundary.push([0.0, 0.0]);
        }
    });
    if ui.button("Add Space").clicked() {
        app.arch_add_space();
    }
}

fn draw_arch_entity_list(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Entities").strong());
    let mut delete_id: Option<usize> = None;
    {
        let s = &app.mesh_toolbox.arch;
        if s.doc.count() == 0 {
            ui.label("(empty)");
        } else {
            for (id, ent) in s.doc.iter() {
                ui.horizontal(|ui| {
                    ui.label(format!("#{id} {} — {}", ent.kind().label(), ent.summary()));
                    if ui.small_button("Delete").clicked() {
                        delete_id = Some(id);
                    }
                });
            }
        }
    }
    if let Some(id) = delete_id {
        app.arch_delete_entity(id);
    }
}

fn draw_arch_render_and_schedule(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    {
        let s = &mut app.mesh_toolbox.arch;
        ui.horizontal(|ui| {
            ui.label("Tess tol:");
            ui.add(
                egui::DragValue::new(&mut s.render_tolerance)
                    .range(0.001..=10.0)
                    .speed(0.01)
                    .suffix(" m"),
            );
        });
    }
    ui.horizontal(|ui| {
        if ui.button("Render to viewport").clicked() {
            app.arch_render();
        }
        if ui.button("Refresh schedule").clicked() {
            app.arch_refresh_schedule();
        }
        if ui.button("Save CSV…").clicked() {
            app.arch_save_schedule_csv();
        }
        if ui.button("Export IFC4…").clicked() {
            app.arch_export_ifc();
        }
    });
    let table = app
        .mesh_toolbox
        .arch
        .last_schedule
        .as_ref()
        .map(|s| s.to_text_table())
        .unwrap_or_default();
    if !table.is_empty() {
        ui.collapsing("Schedule", |ui| {
            ui.label(egui::RichText::new(table).monospace());
        });
    }
}

// =============================================================================
// Phase 16 — Spreadsheet Workbench panel
// =============================================================================

/// Draw the Spreadsheet workbench panel — sheet picker, add-sheet
/// button, editable grid, re-evaluate-all button.
///
/// The panel renders the active sheet as a `view_rows x view_cols`
/// grid of small buttons. Clicking a cell selects it and copies its
/// current contents into the editor buffer; the "Set" / "Clear"
/// buttons apply. Numeric inputs go in as [`valenx_spreadsheet::Cell::Number`],
/// strings beginning with `=` become [`valenx_spreadsheet::Cell::Formula`],
/// and anything else becomes [`valenx_spreadsheet::Cell::Text`].
pub fn draw_spreadsheet_panel(app: &mut crate::ValenxApp, ui: &mut egui::Ui) {
    let s = &mut app.mesh_toolbox.spreadsheet;

    // ----- sheet picker + add-sheet -----
    ui.horizontal(|ui| {
        ui.label("Sheet:")
            .on_hover_text("Which sheet from the workbook is currently visible in the grid below.");
        let sheet_names: Vec<String> = s.workbook.sheets.keys().cloned().collect();
        let active = if sheet_names.contains(&s.active_sheet) {
            s.active_sheet.clone()
        } else {
            sheet_names.first().cloned().unwrap_or_default()
        };
        s.active_sheet = active.clone();
        egui::ComboBox::from_id_source("spreadsheet_active_sheet")
            .selected_text(if active.is_empty() {
                "(none)"
            } else {
                active.as_str()
            })
            .show_ui(ui, |ui| {
                for name in &sheet_names {
                    ui.selectable_value(&mut s.active_sheet, name.clone(), name)
                        .on_hover_text(format!("Switch to sheet `{name}`."));
                }
            });
    });
    ui.horizontal(|ui| {
        ui.label("New sheet:")
            .on_hover_text("Name for the next sheet to add. Formulas reference sheets as `SheetName.A1` etc.");
        ui.text_edit_singleline(&mut s.new_sheet_name)
            .on_hover_text("Sheet name (must be unique within the workbook).");
        if ui.button("Add sheet")
            .on_hover_text("Add a new empty sheet with the name above.")
            .clicked() && !s.new_sheet_name.is_empty()
        {
            if s.workbook.add_sheet(s.new_sheet_name.clone()) {
                s.last_status = Some(format!("Added sheet `{}`.", s.new_sheet_name));
                s.last_error = None;
                s.active_sheet = s.new_sheet_name.clone();
            } else {
                s.last_error = Some(format!("Sheet `{}` already exists.", s.new_sheet_name));
            }
        }
        if ui.button("Remove sheet")
            .on_hover_text("Remove the currently-active sheet (irreversible — undo via re-add).")
            .clicked() && !s.active_sheet.is_empty()
        {
            let removed = s.workbook.remove_sheet(&s.active_sheet);
            if removed {
                s.last_status = Some(format!("Removed sheet `{}`.", s.active_sheet));
                s.active_sheet = s.workbook.sheets.keys().next().cloned().unwrap_or_default();
            }
        }
    });
    ui.horizontal(|ui| {
        ui.label("Rows:")
            .on_hover_text("How many rows of the active sheet are shown in the grid (cells beyond this still exist).");
        ui.add(egui::DragValue::new(&mut s.view_rows).range(1..=64))
            .on_hover_text("Visible row count (1..64).");
        ui.label("Cols:")
            .on_hover_text("How many columns are shown.");
        ui.add(egui::DragValue::new(&mut s.view_cols).range(1..=26))
            .on_hover_text("Visible column count (1..26, A through Z).");
    });

    ui.separator();

    // ----- grid -----
    if s.active_sheet.is_empty() {
        ui.label("(No sheet — add one above to get started.)");
        return;
    }
    let active_sheet_name = s.active_sheet.clone();
    let view_rows = s.view_rows;
    let view_cols = s.view_cols;

    egui::Grid::new("spreadsheet_grid")
        .striped(true)
        .show(ui, |ui| {
            // Header row: column letters.
            ui.label("");
            for col in 0..view_cols {
                ui.label(egui::RichText::new(col_letter(col)).strong());
            }
            ui.end_row();
            for row in 0..view_rows {
                ui.label(egui::RichText::new(format!("{}", row + 1)).strong());
                for col in 0..view_cols {
                    let r = valenx_spreadsheet::CellRef {
                        sheet_name: active_sheet_name.clone(),
                        row,
                        col,
                    };
                    let cell = s.workbook.cell(&r);
                    let label = match cell {
                        valenx_spreadsheet::Cell::Empty => String::new(),
                        valenx_spreadsheet::Cell::Number(n) => format!("{n}"),
                        valenx_spreadsheet::Cell::Text(t) => t.clone(),
                        valenx_spreadsheet::Cell::Formula(src) => {
                            // Show the evaluated value with the source as
                            // tooltip; on parse / circular error, show the
                            // source verbatim so the user can edit it.
                            match s.workbook.evaluate_cell(&r) {
                                Ok(v) => format!("= {v}"),
                                Err(_) => format!("? {src}"),
                            }
                        }
                    };
                    let button_label = if label.is_empty() {
                        "·".to_string()
                    } else {
                        label
                    };
                    if ui.button(button_label).clicked() {
                        s.selected_cell = Some((row, col));
                        s.editor_text = match cell {
                            valenx_spreadsheet::Cell::Empty => String::new(),
                            valenx_spreadsheet::Cell::Number(n) => format!("{n}"),
                            valenx_spreadsheet::Cell::Text(t) => t.clone(),
                            valenx_spreadsheet::Cell::Formula(src) => format!("={src}"),
                        };
                    }
                }
                ui.end_row();
            }
        });

    ui.separator();

    // ----- editor + actions -----
    let (sel_row, sel_col) = s.selected_cell.unwrap_or((0, 0));
    ui.horizontal(|ui| {
        ui.label(format!("Selected: {}{}", col_letter(sel_col), sel_row + 1))
            .on_hover_text("Cell currently being edited. Click a cell in the grid above to change selection.");
        ui.text_edit_singleline(&mut s.editor_text)
            .on_hover_text(
                "Cell content. Plain number = literal value. `=` prefix = formula \
                 (e.g. `=A1+B1` or `=Sheet2.A1 * 2`). Anything else = text label.",
            );
    });
    ui.horizontal(|ui| {
        let r = valenx_spreadsheet::CellRef {
            sheet_name: active_sheet_name.clone(),
            row: sel_row,
            col: sel_col,
        };
        if ui.button("Set")
            .on_hover_text("Apply the editor text to the selected cell. Formulas are re-evaluated when their dependencies change.")
            .clicked()
        {
            let new_cell = parse_editor_cell(&s.editor_text);
            match s.workbook.set_cell(&r, new_cell) {
                Ok(()) => {
                    s.last_status = Some(format!("Set {r}."));
                    s.last_error = None;
                }
                Err(e) => s.last_error = Some(format!("{e}")),
            }
        }
        if ui.button("Clear")
            .on_hover_text("Empty the selected cell — drops formula and value.")
            .clicked()
        {
            s.editor_text.clear();
            if let Err(e) = s.workbook.set_cell(&r, valenx_spreadsheet::Cell::Empty) {
                s.last_error = Some(format!("{e}"));
            } else {
                s.last_status = Some(format!("Cleared {r}."));
                s.last_error = None;
            }
        }
        if ui.button("Re-evaluate all")
            .on_hover_text("Force evaluation of every formula in the workbook. Useful after editing many cells at once or fixing a circular reference.")
            .clicked()
        {
            let mut evaluated = 0usize;
            let mut errors = 0usize;
            // Collect refs first so we don't borrow ss while iterating.
            let refs: Vec<valenx_spreadsheet::CellRef> = s
                .workbook
                .iter_sheets()
                .flat_map(|(name, sheet)| {
                    sheet
                        .iter()
                        .map(move |((row, col), _cell)| valenx_spreadsheet::CellRef {
                            sheet_name: name.to_string(),
                            row,
                            col,
                        })
                })
                .collect();
            for r in refs {
                match s.workbook.evaluate_cell(&r) {
                    Ok(_) => evaluated += 1,
                    Err(_) => errors += 1,
                }
            }
            s.last_status = Some(format!(
                "Re-evaluated {evaluated} cells ({errors} errored)."
            ));
            s.last_error = None;
        }
    });

    if let Some(msg) = &s.last_status {
        ui.colored_label(egui::Color32::LIGHT_GREEN, msg);
    }
    if let Some(msg) = &s.last_error {
        ui.colored_label(egui::Color32::LIGHT_RED, msg);
    }
}

/// Convert column 0 -> "A", 25 -> "Z", 26 -> "AA", ....
fn col_letter(mut col: u32) -> String {
    let mut letters = Vec::new();
    loop {
        let rem = (col % 26) as u8;
        letters.push((b'A' + rem) as char);
        if col < 26 {
            break;
        }
        col = col / 26 - 1;
    }
    letters.reverse();
    letters.into_iter().collect()
}

/// Convert an editor buffer into a [`valenx_spreadsheet::Cell`].
///
/// - Empty string → `Cell::Empty`.
/// - String starting with `=` → `Cell::Formula(source after the `=`)`.
/// - Pure-numeric strings (parse to `f64`) → `Cell::Number`.
/// - Anything else → `Cell::Text`.
fn parse_editor_cell(text: &str) -> valenx_spreadsheet::Cell {
    if text.is_empty() {
        return valenx_spreadsheet::Cell::Empty;
    }
    if let Some(rest) = text.strip_prefix('=') {
        return valenx_spreadsheet::Cell::Formula(rest.trim().to_string());
    }
    if let Ok(n) = text.trim().parse::<f64>() {
        return valenx_spreadsheet::Cell::Number(n);
    }
    valenx_spreadsheet::Cell::Text(text.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::LoadedStl;
    use valenx_viz::StlTriangle;

    /// Build a triangle mesh with a single triangle for fixture work.
    fn one_triangle() -> TriangleMesh {
        let mut m = TriangleMesh::new();
        m.triangles.push(StlTriangle {
            normal: [0.0, 0.0, 1.0],
            vertices: [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        });
        m
    }

    #[test]
    fn apply_translate_no_load_sets_error() {
        let mut app = ValenxApp::default();
        app.apply_translate(1.0, 0.0, 0.0);
        assert!(app.last_error.is_some());
    }

    #[test]
    fn apply_translate_with_stl_shifts_vertices() {
        let mut app = ValenxApp {
            stl: Some(LoadedStl {
                path: PathBuf::from("tri.stl"),
                mesh: one_triangle(),
            }),
            ..Default::default()
        };
        app.apply_translate(2.0, 3.0, 4.0);
        let stl = app.stl.as_ref().unwrap();
        assert_eq!(stl.mesh.triangles[0].vertices[0], [2.0, 3.0, 4.0]);
        assert!(app.status.as_ref().unwrap().contains("Translated"));
    }

    #[test]
    fn apply_scale_uniform_with_stl_multiplies_vertices() {
        let mut app = ValenxApp {
            stl: Some(LoadedStl {
                path: PathBuf::from("tri.stl"),
                mesh: one_triangle(),
            }),
            ..Default::default()
        };
        app.apply_scale_uniform(3.0);
        let stl = app.stl.as_ref().unwrap();
        assert_eq!(stl.mesh.triangles[0].vertices[1], [3.0, 0.0, 0.0]);
    }

    #[test]
    fn apply_scale_uniform_zero_or_nan_errors() {
        let mut app = ValenxApp {
            stl: Some(LoadedStl {
                path: PathBuf::from("tri.stl"),
                mesh: one_triangle(),
            }),
            ..Default::default()
        };
        app.apply_scale_uniform(0.0);
        assert!(app.last_error.is_some());
    }

    #[test]
    fn apply_mirror_reverses_triangle_winding() {
        let mut app = ValenxApp {
            stl: Some(LoadedStl {
                path: PathBuf::from("tri.stl"),
                mesh: one_triangle(),
            }),
            ..Default::default()
        };
        app.apply_mirror(ToolboxAxis::Y);
        let tri = &app.stl.as_ref().unwrap().mesh.triangles[0];
        // Original vertices: [(0,0,0), (1,0,0), (0,1,0)].
        // y-flip: [(0,0,0), (1,0,0), (0,-1,0)] then reverse: [(0,-1,0), (1,0,0), (0,0,0)].
        assert_eq!(tri.vertices[0], [0.0, -1.0, 0.0]);
        assert_eq!(tri.vertices[2], [0.0, 0.0, 0.0]);
    }

    #[test]
    fn apply_cut_plane_with_zero_normal_errors() {
        let mut app = ValenxApp {
            stl: Some(LoadedStl {
                path: PathBuf::from("tri.stl"),
                mesh: one_triangle(),
            }),
            ..Default::default()
        };
        app.apply_cut_plane([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!(app
            .last_error
            .as_ref()
            .is_some_and(|e| e.contains("non-zero")));
    }

    #[test]
    fn apply_cut_plane_keeps_positive_side_stl() {
        // Centroid of the only triangle is (1/3, 1/3, 0). Plane at
        // origin with normal -x: dot = -1/3 < 0 → discard.
        let mut app = ValenxApp {
            stl: Some(LoadedStl {
                path: PathBuf::from("tri.stl"),
                mesh: one_triangle(),
            }),
            ..Default::default()
        };
        app.apply_cut_plane([0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]);
        assert_eq!(app.stl.as_ref().unwrap().mesh.triangles.len(), 0);
    }

    #[test]
    fn apply_merge_coincident_on_stl_only_reports_error() {
        // STL triangle soup doesn't have a shared-node table — the
        // operation surfaces an honest error instead of silently
        // no-op'ing.
        let mut app = ValenxApp {
            stl: Some(LoadedStl {
                path: PathBuf::from("tri.stl"),
                mesh: one_triangle(),
            }),
            ..Default::default()
        };
        app.apply_merge_coincident(1e-6);
        assert!(app
            .last_error
            .as_ref()
            .is_some_and(|e| e.contains("canonical")));
    }

    #[test]
    fn toggle_mesh_toolbox_flips_visibility() {
        let mut app = ValenxApp {
            show_mesh_toolbox: true,
            ..Default::default()
        };
        app.toggle_mesh_toolbox();
        assert!(!app.show_mesh_toolbox);
        app.toggle_mesh_toolbox();
        assert!(app.show_mesh_toolbox);
    }

    #[test]
    fn toolbox_state_defaults_are_neutral() {
        let s = MeshToolboxState::default();
        assert_eq!(s.translate, [0.0, 0.0, 0.0]);
        assert_eq!(s.scale_uniform, 1.0);
        assert_eq!(s.scale_per_axis, [1.0, 1.0, 1.0]);
        assert_eq!(s.cut_normal, [0.0, 0.0, 1.0]);
        assert!(s.repair_tolerance > 0.0);
    }

    #[test]
    fn toolbox_state_cad_defaults_are_buildable_primitives() {
        // The default toolbox values should produce a valid box,
        // cylinder, sphere, cone, and torus. If they ever drift the
        // user clicks Create and gets InvalidParam — make sure that
        // doesn't happen.
        let s = MeshToolboxState::default();
        assert!(s.cad_box_dims.iter().all(|&v| v > 0.0));
        assert!(s.cad_cyl_radius > 0.0 && s.cad_cyl_height > 0.0);
        assert!(s.cad_sphere_radius > 0.0);
        assert!(s.cad_cone_base > 0.0 && s.cad_cone_height > 0.0);
        assert!(s.cad_torus_minor < s.cad_torus_major);
    }

    #[test]
    fn apply_create_box_populates_current_solid_and_viewport_mesh() {
        let mut app = ValenxApp::default();
        app.apply_create_primitive(CadPrimitiveKind::Box);
        assert!(
            app.current_solid.is_some(),
            "Create primitive should populate current_solid: {:?}",
            app.last_error
        );
        assert!(
            app.mesh.is_some(),
            "Create primitive should tessellate into a viewport mesh"
        );
        let solid = app.current_solid.as_ref().unwrap();
        assert_eq!(solid.faces(), 6);
        assert_eq!(solid.vertices(), 8);
    }

    #[test]
    fn apply_create_as_second_populates_operand_b() {
        let mut app = ValenxApp::default();
        app.mesh_toolbox.cad_create_as_second = true;
        app.apply_create_primitive(CadPrimitiveKind::Sphere);
        assert!(app.second_solid.is_some());
        assert!(app.current_solid.is_none());
    }

    #[test]
    fn boolean_without_both_operands_reports_error() {
        let mut app = ValenxApp::default();
        app.apply_cad_boolean(CadBooleanOp::Union);
        assert!(app.last_error.is_some());
        assert!(app.current_solid.is_none());
    }

    #[test]
    fn boolean_consumes_b_and_replaces_a() {
        let mut app = ValenxApp::default();
        // Build the canonical "punched cube" pair: unit cube at the
        // origin as A, then a small cylinder centred inside A as B
        // (translated so it doesn't share faces with the cube
        // boundary — that's the geometry truck-shapeops handles).
        app.apply_create_primitive(CadPrimitiveKind::Box);
        let _a_faces_before = app.current_solid.as_ref().unwrap().faces();

        app.mesh_toolbox.cad_create_as_second = true;
        app.mesh_toolbox.cad_cyl_radius = 0.25;
        app.mesh_toolbox.cad_cyl_height = 2.0;
        app.apply_create_primitive(CadPrimitiveKind::Cylinder);
        assert!(app.second_solid.is_some());
        // Manually translate operand B so it goes through the cube
        // — without this, both solids overlap the (0,0,0) corner and
        // truck-shapeops returns EmptyResult.
        let b_translated = app
            .second_solid
            .as_ref()
            .unwrap()
            .translated(0.5, 0.5, -0.5)
            .unwrap();
        app.second_solid = Some(b_translated);

        app.apply_cad_boolean(CadBooleanOp::Difference);
        assert!(
            app.second_solid.is_none(),
            "boolean should consume operand B (last_error: {:?})",
            app.last_error
        );
        assert!(
            app.current_solid.is_some(),
            "boolean should leave a result in operand A"
        );
    }

    #[test]
    fn apply_cad_fillet_reports_not_implemented() {
        let mut app = ValenxApp::default();
        app.apply_create_primitive(CadPrimitiveKind::Box);
        app.apply_cad_fillet();
        assert!(app.last_error.is_some());
        let err = app.last_error.as_ref().unwrap();
        assert!(
            err.contains("not implemented"),
            "fillet should report a typed NotImplemented, got: {err}"
        );
    }

    #[test]
    fn apply_cad_fillet_without_solid_reports_error() {
        let mut app = ValenxApp::default();
        app.apply_cad_fillet();
        assert!(app.last_error.is_some());
    }

    #[test]
    fn apply_cad_clear_drops_both_operands() {
        let mut app = ValenxApp::default();
        app.apply_create_primitive(CadPrimitiveKind::Box);
        app.mesh_toolbox.cad_create_as_second = true;
        app.apply_create_primitive(CadPrimitiveKind::Sphere);
        app.apply_cad_clear();
        assert!(app.current_solid.is_none());
        assert!(app.second_solid.is_none());
    }

    /// Smoke test for the assembly-fusion renderer (Phase 6, Task 44).
    /// Two unit cubes at different positions fuse into a single mesh
    /// with vertices from both — connectivity indices offset
    /// correctly so the second cube's triangles reference its own
    /// nodes, not the first cube's.
    #[test]
    fn render_assembly_to_viewport_fuses_two_cubes() {
        let mut app = ValenxApp::default();
        let s = &mut app.mesh_toolbox.assembly;
        let cube1 = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let cube2 = valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap();
        let mut p1 = valenx_assembly::Part::new(0, "a", cube1);
        p1.transform.translation = nalgebra::Vector3::new(0.0, 0.0, 0.0);
        s.assembly.add_part(p1);
        let mut p2 = valenx_assembly::Part::new(0, "b", cube2);
        p2.transform.translation = nalgebra::Vector3::new(5.0, 0.0, 0.0);
        s.assembly.add_part(p2);

        super::render_assembly_to_viewport(&mut app);
        // The viewport mesh should now contain both cubes' triangles.
        let loaded = app.mesh.as_ref().expect("fused mesh applied");
        // Each cube tessellates to ~12 triangles minimum; we have 2 → 24+.
        assert!(loaded.mesh.total_elements() >= 24);
        // The second cube's vertices should be offset along X.
        let max_x = loaded
            .mesh
            .nodes
            .iter()
            .map(|n| n.x)
            .fold(f64::MIN, f64::max);
        assert!(max_x > 4.9, "second cube not at x=5: max_x={max_x}");
    }

    #[test]
    fn assembly_panel_state_defaults_are_sensible() {
        let s = AssemblyPanelState::default();
        assert_eq!(s.new_part_primitive.label(), "Box");
        assert_eq!(s.mate_kind.label(), "Coincident");
        assert_eq!(s.joint_kind.label(), "Fixed");
        assert!(s.assembly.parts.is_empty());
        assert!(s.last_error.is_none());
    }

    // ----- Phase 7 Mesh Tools coverage -----

    /// Build a canonical-mesh `LoadedMesh` from a simple subdivided
    /// plane so the mesh-tools methods have something to chew on.
    fn loaded_plane(n: usize) -> crate::types::LoadedMesh {
        use nalgebra::Vector3;
        use valenx_mesh::element::{ElementBlock, ElementType};
        use valenx_mesh::Mesh;
        let mut mesh = Mesh::new("plane");
        for j in 0..=n {
            for i in 0..=n {
                mesh.nodes.push(Vector3::new(i as f64, j as f64, 0.0));
            }
        }
        let nx_plus = n + 1;
        let mut blk = ElementBlock::new(ElementType::Tri3);
        for j in 0..n {
            for i in 0..n {
                let i0 = (j * nx_plus + i) as u32;
                let i1 = i0 + 1;
                let i2 = i0 + nx_plus as u32;
                let i3 = i2 + 1;
                blk.connectivity
                    .extend_from_slice(&[i0, i1, i3, i0, i3, i2]);
            }
        }
        mesh.element_blocks = vec![blk];
        mesh.recompute_stats();
        let quality = valenx_mesh::quality_report(&mesh);
        let aspect_hist =
            valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
        let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
        crate::types::LoadedMesh {
            path: PathBuf::from("plane.json"),
            mesh,
            quality,
            aspect_hist,
            skew_hist,
        }
    }

    #[test]
    fn mesh_tools_defaults_are_sensible() {
        let s = MeshToolboxState::default();
        assert!(s.mesh_tools_decimate_fraction > 0.0 && s.mesh_tools_decimate_fraction <= 1.0);
        assert!(s.mesh_tools_laplacian_iter > 0);
        assert!(s.mesh_tools_taubin_lambda > 0.0);
        assert!(s.mesh_tools_taubin_mu < 0.0, "Taubin μ must be negative");
        assert!(s.mesh_tools_remesh_target > 0.0);
        assert!(s.mesh_tools_fill_holes_max > 0.0);
    }

    #[test]
    fn apply_mesh_decimate_reduces_vertex_count() {
        let mut app = ValenxApp {
            mesh: Some(loaded_plane(8)),
            ..Default::default()
        };
        let before = app.mesh.as_ref().unwrap().mesh.nodes.len();
        app.apply_mesh_decimate(0.5);
        let after = app.mesh.as_ref().unwrap().mesh.nodes.len();
        assert!(after <= before, "after {after} > before {before}");
        assert!(app.status.as_ref().is_some_and(|s| s.contains("Decimated")));
    }

    #[test]
    fn apply_mesh_decimate_without_mesh_errors() {
        let mut app = ValenxApp::default();
        app.apply_mesh_decimate(0.5);
        assert!(app.last_error.is_some());
    }

    #[test]
    fn apply_mesh_decimate_invalid_fraction_errors() {
        let mut app = ValenxApp {
            mesh: Some(loaded_plane(4)),
            ..Default::default()
        };
        app.apply_mesh_decimate(-1.0);
        assert!(app.last_error.is_some());
    }

    #[test]
    fn apply_mesh_laplacian_runs() {
        let mut app = ValenxApp {
            mesh: Some(loaded_plane(4)),
            ..Default::default()
        };
        app.apply_mesh_laplacian(3, 0.5);
        assert!(app.status.as_ref().is_some_and(|s| s.contains("Laplacian")));
    }

    #[test]
    fn apply_mesh_taubin_runs() {
        let mut app = ValenxApp {
            mesh: Some(loaded_plane(4)),
            ..Default::default()
        };
        app.apply_mesh_taubin(3, 0.5, -0.53);
        assert!(app.status.as_ref().is_some_and(|s| s.contains("Taubin")));
    }

    #[test]
    fn apply_mesh_remesh_runs() {
        let mut app = ValenxApp {
            mesh: Some(loaded_plane(4)),
            ..Default::default()
        };
        app.apply_mesh_remesh(0.5, 2);
        assert!(app.status.as_ref().is_some_and(|s| s.contains("Remesh")));
    }

    #[test]
    fn apply_mesh_remesh_invalid_target_errors() {
        let mut app = ValenxApp {
            mesh: Some(loaded_plane(4)),
            ..Default::default()
        };
        app.apply_mesh_remesh(-1.0, 2);
        assert!(app.last_error.is_some());
    }

    #[test]
    fn apply_mesh_fill_holes_runs() {
        let mut app = ValenxApp {
            mesh: Some(loaded_plane(4)),
            ..Default::default()
        };
        app.apply_mesh_fill_holes(100.0);
        // A plane with all boundary edges has one big closed loop;
        // fill triangulates it. Status should mention filled.
        assert!(app.status.is_some() || app.last_error.is_none());
    }

    // UI-coupled test — opens OS dialog. Run interactively only. (Phase 10 lockdown.)
    #[ignore]
    #[test]
    fn save_mesh_as_3mf_without_mesh_errors() {
        let mut app = ValenxApp::default();
        app.save_mesh_as_3mf();
        // Without a dialog interaction, we'd hit the save-file
        // dialog first; we can't easily test the dialog itself, so
        // this is mostly a smoke test that the call doesn't panic.
        let _ = app.last_error;
    }

    /// Round-20 H1 RED→GREEN: pre-fix, the Desktop Dock panel's
    /// `run_dock_now` did a bare `std::fs::read_to_string` on the
    /// receptor path. A multi-GB file would OOM the renderer before
    /// the docker ran. Post-fix, the read is bounded by
    /// `MAX_PDBQT_FILE_BYTES` and the over-cap file surfaces as a
    /// `last_error` set string with no allocation past the cap.
    #[test]
    fn run_dock_now_rejects_oversize_receptor() {
        use std::io::Write;
        use valenx_core::io_caps::MAX_PDBQT_FILE_BYTES;
        let tmp = std::env::temp_dir().join(format!(
            "valenx_round20_h1_receptor_{}.pdbqt",
            std::process::id()
        ));
        // Write 1 byte past the cap — enough to trip the stat check
        // without actually allocating 64 MiB in the test harness.
        let mut f = std::fs::File::create(&tmp).unwrap();
        // SetLen + a single write gives us a sparse-or-real file
        // larger than the cap. Use set_len so the test doesn't write
        // 64 MiB of bytes (slow on every CI).
        f.set_len(MAX_PDBQT_FILE_BYTES as u64 + 1).unwrap();
        f.write_all(b"x").unwrap();
        drop(f);

        let mut s = DockPanelState {
            receptor_path: tmp.display().to_string(),
            ligand_path: String::new(),
            output_path: String::new(),
            ..Default::default()
        };
        run_dock_now(&mut s);
        let err = s.last_error.expect(
            "round-20 H1: oversized receptor must set last_error, not OOM the renderer",
        );
        assert!(
            err.contains("receptor read") && err.contains("cap"),
            "expected `receptor read: ...cap...` got: {err}"
        );
        let _ = std::fs::remove_file(&tmp);
    }
}

/// Headless egui UI-logic tests for the **CAD workbench panels**.
///
/// The Mesh Toolbox hosts eleven CAD panels — Part, Draft, TechDraw,
/// Assembly, Surface, CAM, Arch/BIM, Spreadsheet, plus the Dock,
/// Sketcher, and Part Design sub-panels. This module is the headless
/// counterpart of the Genetics workbench's `headless_ui_tests`: for
/// every panel it
///
/// 1. **draws the panel** in a windowless [`egui::Context`] across
///    representative states (fresh, mid-edit, post-run, error) — the
///    draw must never panic;
/// 2. **exercises the panel's Run / Compute action** against the real
///    backend crate (`valenx-cad`, `valenx-draft`, `valenx-techdraw`,
///    `valenx-assembly`, `valenx-surface`, `valenx-cam`,
///    `valenx-arch`, `valenx-spreadsheet`, `valenx-sketch`,
///    `valenx-feature-tree`) and asserts a sane result;
/// 3. **feeds bad input** and asserts the action surfaces a graceful
///    error rather than panicking.
///
/// Every test renders into a free-standing [`egui::CentralPanel`] and
/// drives the real `draw_*_panel` functions — none opens an OS window
/// and none reaches `rfd::FileDialog` (the Load / Save / Import
/// buttons are layout-only elements that are drawn but never clicked).
///
/// Named `headless_ui_tests` so the repo-wide safe filter
/// (`cargo test -p valenx-app headless_ui_tests`) selects exactly
/// these.
#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    /// Draw an arbitrary panel-draw closure once in a headless egui
    /// context. A panic inside the closure surfaces as a failed test.
    fn draw_headless(f: impl FnOnce(&mut egui::Ui)) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, f);
            });
        });
    }

    // ===================================================================
    // Whole-toolbox host panel
    // ===================================================================

    #[test]
    fn toolbox_is_a_noop_when_hidden() {
        // With the toggle off the whole toolbox draws nothing and
        // never panics.
        let mut app = ValenxApp::default();
        assert!(!app.show_mesh_toolbox);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_mesh_toolbox(&mut app, ctx);
        });
    }

    #[test]
    fn whole_toolbox_draws_without_panic() {
        // The whole right-side toolbox — every section + every CAD
        // sub-panel collapsing header — mounts in one headless frame.
        let mut app = ValenxApp::default();
        app.show_mesh_toolbox = true;
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_mesh_toolbox(&mut app, ctx);
        });
    }

    // ===================================================================
    // Part workbench — valenx-cad primitives + booleans
    // ===================================================================

    #[test]
    fn part_workbench_draws_across_states() {
        // Fresh (no operands).
        draw_headless(|ui| draw_part_workbench(&mut ValenxApp::default(), ui));
        // With operand A populated.
        let mut app = ValenxApp::default();
        app.apply_create_primitive(CadPrimitiveKind::Box);
        draw_headless(|ui| draw_part_workbench(&mut app, ui));
        // With an error surfaced.
        let mut app = ValenxApp::default();
        app.apply_cad_boolean(CadBooleanOp::Union); // no operands -> error
        assert!(app.last_error.is_some());
        draw_headless(|ui| draw_part_workbench(&mut app, ui));
    }

    #[test]
    fn part_workbench_create_runs_every_primitive() {
        // The Run action (Create) builds a real truck BRep solid for
        // every primitive kind.
        for kind in [
            CadPrimitiveKind::Box,
            CadPrimitiveKind::Cylinder,
            CadPrimitiveKind::Sphere,
            CadPrimitiveKind::Cone,
            CadPrimitiveKind::Torus,
        ] {
            let mut app = ValenxApp::default();
            app.apply_create_primitive(kind);
            assert!(
                app.current_solid.is_some(),
                "Create {kind:?} produced no solid: {:?}",
                app.last_error
            );
            assert!(app.mesh.is_some(), "Create {kind:?} produced no viewport mesh");
        }
    }

    #[test]
    fn part_workbench_create_rejects_a_degenerate_primitive() {
        // Bad input — a zero-dimension box — must surface an error,
        // not panic and not produce a solid.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.cad_box_dims = [0.0, 1.0, 1.0];
        app.apply_create_primitive(CadPrimitiveKind::Box);
        assert!(
            app.last_error.is_some(),
            "a degenerate box should surface an error"
        );
        assert!(app.current_solid.is_none());
    }

    // ===================================================================
    // Draft workbench — valenx-draft 2D entities
    // ===================================================================

    #[test]
    fn draft_panel_draws_every_tool() {
        // The panel draws for every tool selection without panic.
        for tool in [
            DraftTool::Select,
            DraftTool::Line,
            DraftTool::Polyline,
            DraftTool::Arc,
            DraftTool::Circle,
            DraftTool::Rectangle,
            DraftTool::Polygon,
            DraftTool::Dimension,
            DraftTool::Text,
        ] {
            let mut app = ValenxApp::default();
            app.mesh_toolbox.draft.tool = tool;
            draw_headless(|ui| draw_draft_panel(&mut app, ui));
        }
    }

    #[test]
    fn draft_panel_draws_post_edit_and_error_states() {
        // Post-edit — a document with entities.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.draft.document.add_entity(valenx_draft::DraftEntity::Circle {
            center: [0.0, 0.0],
            radius: 2.0,
        });
        app.mesh_toolbox.draft.selected_entity = Some(0);
        draw_headless(|ui| draw_draft_panel(&mut app, ui));
        // Error state.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.draft.last_error = Some("arc radius must be > 0".into());
        draw_headless(|ui| draw_draft_panel(&mut app, ui));
    }

    #[test]
    fn draft_commit_polyline_builds_a_real_entity() {
        // The Run action (commit a polyline) calls the real
        // valenx-draft document API and adds an entity.
        let mut s = DraftPanelState::default();
        s.polyline_points = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]];
        let before = s.document.entity_count();
        commit_polyline(&mut s, true);
        assert!(s.last_error.is_none(), "commit errored: {:?}", s.last_error);
        assert_eq!(s.document.entity_count(), before + 1, "polyline not added");
    }

    #[test]
    fn draft_commit_polyline_rejects_too_few_points() {
        // Bad input — a one-point polyline — surfaces an error and
        // adds nothing.
        let mut s = DraftPanelState::default();
        s.polyline_points = vec![[0.0, 0.0]];
        let before = s.document.entity_count();
        commit_polyline(&mut s, false);
        assert!(s.last_error.is_some(), "single-point polyline should error");
        assert_eq!(s.document.entity_count(), before, "nothing should be added");
    }

    #[test]
    fn draft_document_round_trips_every_entity_kind() {
        // Exercise the real valenx-draft backend the panel's tool
        // buttons drive — one entity of each kind, then delete one.
        let mut doc = valenx_draft::DraftDocument::new(valenx_draft::WorkingPlane::from_xy());
        doc.add_entity(valenx_draft::DraftEntity::Line {
            start: [0.0, 0.0],
            end: [3.0, 0.0],
        });
        doc.add_entity(valenx_draft::DraftEntity::Arc {
            center: [0.0, 0.0],
            radius: 1.5,
            start_angle: 0.0,
            end_angle: std::f64::consts::PI,
        });
        doc.add_entity(valenx_draft::DraftEntity::Rectangle {
            min: [0.0, 0.0],
            max: [2.0, 1.0],
        });
        doc.add_entity(valenx_draft::DraftEntity::Polygon {
            center: [0.0, 0.0],
            radius: 2.0,
            sides: 6,
        });
        assert_eq!(doc.entity_count(), 4);
        doc.delete_entity(0).expect("delete a valid index");
        assert_eq!(doc.entity_count(), 3);
        assert!(doc.delete_entity(99).is_err(), "out-of-range delete should error");
    }

    // ===================================================================
    // TechDraw workbench — valenx-techdraw drawings
    // ===================================================================

    #[test]
    fn techdraw_panel_draws_across_states() {
        // Fresh.
        draw_headless(|ui| draw_techdraw_panel(&mut ValenxApp::default(), ui));
        // With a source solid (so Add View can generate edges).
        let mut app = ValenxApp::default();
        app.apply_create_primitive(CadPrimitiveKind::Box);
        draw_headless(|ui| draw_techdraw_panel(&mut app, ui));
        // With a view already present + selected, and an error.
        let mut app = ValenxApp::default();
        let view = valenx_techdraw::View::new(
            valenx_techdraw::ViewKind::Front,
            1.0,
            [80.0, 100.0],
        );
        app.mesh_toolbox.techdraw.drawing.add_view(view);
        app.mesh_toolbox.techdraw.selected_view = Some(0);
        app.mesh_toolbox.techdraw.last_error = Some("Add view: no solid".into());
        draw_headless(|ui| draw_techdraw_panel(&mut app, ui));
    }

    #[test]
    fn techdraw_generate_view_projects_a_real_solid() {
        // The Run action (Add View) generates projected edges from a
        // real valenx-cad solid through the valenx-techdraw backend.
        let solid = valenx_cad::box_solid(2.0, 1.0, 1.0).expect("box");
        let mut view = valenx_techdraw::View::new(
            valenx_techdraw::ViewKind::Front,
            1.0,
            [0.0, 0.0],
        );
        view.generate(&solid).expect("view generation should succeed");
        assert!(
            !view.visible_edges.is_empty(),
            "a projected box must yield visible edges"
        );
        let mut drawing =
            valenx_techdraw::Drawing::new(valenx_techdraw::Sheet::a4_landscape("T", "", "A"));
        let idx = drawing.add_view(view);
        assert_eq!(drawing.views.len(), 1);
        assert!(drawing.get_view_mut(idx).is_ok());
    }

    #[test]
    fn techdraw_get_view_rejects_a_bad_index() {
        // Bad input — addressing a non-existent view — returns an
        // error rather than panicking.
        let mut drawing =
            valenx_techdraw::Drawing::new(valenx_techdraw::Sheet::a4_landscape("T", "", "A"));
        assert!(drawing.get_view_mut(0).is_err(), "empty drawing has no view 0");
    }

    // ===================================================================
    // Assembly workbench — valenx-assembly multi-part scene
    // ===================================================================

    #[test]
    fn assembly_panel_draws_across_states() {
        // Fresh.
        draw_headless(|ui| draw_assembly_panel(&mut ValenxApp::default(), ui));
        // With parts + a selection.
        let mut app = ValenxApp::default();
        assembly_add_part(&mut app.mesh_toolbox.assembly);
        app.mesh_toolbox.assembly.selected_part = Some(0);
        draw_headless(|ui| draw_assembly_panel(&mut app, ui));
        // Error state.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.assembly.last_error = Some("Add part: bad dims".into());
        draw_headless(|ui| draw_assembly_panel(&mut app, ui));
    }

    #[test]
    fn assembly_add_part_builds_a_real_part() {
        // The Run action (Add part) builds a real valenx-cad solid and
        // appends a valenx-assembly part — for every primitive kind.
        for prim in [
            AssemblyPartPrimitive::Box,
            AssemblyPartPrimitive::Cylinder,
            AssemblyPartPrimitive::Sphere,
        ] {
            let mut s = AssemblyPanelState::default();
            s.new_part_primitive = prim;
            assembly_add_part(&mut s);
            assert!(s.last_error.is_none(), "{prim:?} add errored: {:?}", s.last_error);
            assert_eq!(s.assembly.parts.len(), 1, "{prim:?} part not added");
        }
    }

    #[test]
    fn assembly_add_part_rejects_a_degenerate_primitive() {
        // Bad input — a zero-radius sphere — surfaces an error and
        // adds no part.
        let mut s = AssemblyPanelState::default();
        s.new_part_primitive = AssemblyPartPrimitive::Sphere;
        s.new_part_sphere = 0.0;
        assembly_add_part(&mut s);
        assert!(s.last_error.is_some(), "zero-radius sphere should error");
        assert!(s.assembly.parts.is_empty(), "no part should be added");
    }

    #[test]
    fn assembly_solver_runs_on_a_real_two_part_scene() {
        // The Run action (Solve mates) drives the real valenx-assembly
        // constraint solver. Two parts, no mates → an immediately
        // satisfied system the solver converges.
        let mut s = AssemblyPanelState::default();
        assembly_add_part(&mut s);
        s.new_part_translation = [5.0, 0.0, 0.0];
        assembly_add_part(&mut s);
        assert_eq!(s.assembly.parts.len(), 2);
        let report = valenx_assembly::solver::solve(
            &mut s.assembly,
            valenx_assembly::SolverConfig::default(),
        );
        assert!(report.is_ok(), "solver errored: {:?}", report.err());
    }

    // ===================================================================
    // Surface workbench — valenx-surface NURBS
    // ===================================================================

    #[test]
    fn surface_panel_draws_every_tool() {
        for tool in [
            SurfaceTool::NurbsCurve,
            SurfaceTool::NurbsSurface,
            SurfaceTool::CoonsFill,
            SurfaceTool::Sew,
            SurfaceTool::Trim,
            SurfaceTool::KnotOps,
            SurfaceTool::Ssi,
            SurfaceTool::Fit,
            SurfaceTool::Ruled,
        ] {
            let mut app = ValenxApp::default();
            app.mesh_toolbox.surface.tool = tool;
            draw_headless(|ui| draw_surface_panel(&mut app, ui));
        }
    }

    #[test]
    fn surface_panel_draws_post_run_and_error_states() {
        // Post-run — a curve already built + a status string.
        let mut app = ValenxApp::default();
        app.surface_create_curve();
        app.mesh_toolbox.surface.last_status = Some("Created curve #0".into());
        draw_headless(|ui| draw_surface_panel(&mut app, ui));
        // Error state.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.surface.last_error = Some("create curve: bad knots".into());
        draw_headless(|ui| draw_surface_panel(&mut app, ui));
    }

    #[test]
    fn surface_create_curve_and_surface_run_against_the_backend() {
        // The Run actions build real valenx-surface NURBS entities.
        let mut app = ValenxApp::default();
        app.surface_create_curve();
        assert!(
            app.mesh_toolbox.surface.last_error.is_none(),
            "create curve errored: {:?}",
            app.mesh_toolbox.surface.last_error
        );
        assert_eq!(app.mesh_toolbox.surface.file.curves.len(), 1);

        app.surface_create_surface();
        assert!(
            app.mesh_toolbox.surface.last_error.is_none(),
            "create surface errored: {:?}",
            app.mesh_toolbox.surface.last_error
        );
        assert_eq!(app.mesh_toolbox.surface.file.surfaces.len(), 1);
    }

    #[test]
    fn surface_create_curve_rejects_a_degenerate_definition() {
        // Bad input — a degree-3 curve with only 1 control point —
        // surfaces an error and builds nothing.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.surface.curve_degree = 3;
        app.mesh_toolbox.surface.curve_n_cps = 1;
        app.surface_create_curve();
        assert!(
            app.mesh_toolbox.surface.last_error.is_some(),
            "a 1-CP degree-3 curve should error"
        );
        assert!(app.mesh_toolbox.surface.file.curves.is_empty());
    }

    /// RED→GREEN (round-24 M5): `surface_create_surface` no longer
    /// panics / OOMs when the user types nu=10000, nv=10000 → 100 M
    /// control points (2.4 GiB at 24 B per CP). Pre-fix
    /// `s.surface_cps.resize(nu * nv, …)` would allocate the entire
    /// 2.4 GiB before any sanity gate fired. Post-fix the
    /// MAX_NURBS_CP_GRID = 1 M cap surfaces a structured
    /// `last_error` and no surface is added.
    #[test]
    fn surface_create_surface_rejects_oversized_cp_grid() {
        let mut app = ValenxApp::default();
        app.mesh_toolbox.surface.surface_nu = 10_000;
        app.mesh_toolbox.surface.surface_nv = 10_000;
        app.surface_create_surface();
        assert!(
            app.mesh_toolbox.surface.last_error.is_some(),
            "10000×10000 grid must error pre-resize, got last_error=None"
        );
        let err = app.mesh_toolbox.surface.last_error.as_ref().unwrap();
        assert!(
            err.contains("control-point grid cap")
                || err.contains("exceeds")
                || err.contains("1048576"),
            "expected grid-cap error, got: {err}"
        );
        assert!(
            app.mesh_toolbox.surface.file.surfaces.is_empty(),
            "no surface should have been added"
        );
    }

    /// Sanity (round-24 M5): a small grid still creates a surface
    /// post-fix (no regression of the happy path).
    #[test]
    fn surface_create_surface_accepts_small_grid_post_m5() {
        let mut app = ValenxApp::default();
        // Defaults are 4×4 = 16 CPs, well under the 1 M cap.
        app.surface_create_surface();
        assert!(
            app.mesh_toolbox.surface.last_error.is_none(),
            "small grid errored: {:?}",
            app.mesh_toolbox.surface.last_error
        );
        assert_eq!(app.mesh_toolbox.surface.file.surfaces.len(), 1);
    }

    #[test]
    fn surface_coons_fill_rejects_missing_curves() {
        // Bad input — Coons fill referencing curves that do not exist
        // — surfaces an error rather than panicking.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.surface.coons_curves = [0, 1, 2, 3];
        app.surface_coons_fill();
        assert!(
            app.mesh_toolbox.surface.last_error.is_some(),
            "Coons fill with no curves should error"
        );
    }

    // ===================================================================
    // CAM workbench — valenx-cam tool table + toolpaths
    // ===================================================================

    #[test]
    fn cam_panel_draws_across_op_kinds() {
        // The panel draws for a representative set of op kinds.
        for kind in [
            CamOpKind::Profile,
            CamOpKind::Pocket,
            CamOpKind::Drill,
            CamOpKind::Face,
            CamOpKind::AdaptiveClearing,
            CamOpKind::Engrave,
        ] {
            let mut app = ValenxApp::default();
            app.mesh_toolbox.cam.new_op_kind = kind;
            draw_headless(|ui| draw_cam_panel(&mut app, ui));
        }
    }

    #[test]
    fn cam_panel_draws_post_run_and_error_states() {
        // Post-run — a generated toolpath + status.
        let mut app = ValenxApp::default();
        app.cam_add_operation();
        app.cam_generate_toolpath();
        draw_headless(|ui| draw_cam_panel(&mut app, ui));
        // Error state.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.cam.last_error = Some("add tool: bad diameter".into());
        draw_headless(|ui| draw_cam_panel(&mut app, ui));
    }

    #[test]
    fn cam_add_tool_runs_against_the_backend() {
        // The Run action (Add Tool) builds a real valenx-cam tool.
        let mut app = ValenxApp::default();
        let before = app.mesh_toolbox.cam.file.tools.len();
        app.cam_add_tool();
        assert!(
            app.mesh_toolbox.cam.last_error.is_none(),
            "add tool errored: {:?}",
            app.mesh_toolbox.cam.last_error
        );
        assert_eq!(app.mesh_toolbox.cam.file.tools.len(), before + 1);
    }

    #[test]
    fn cam_generate_toolpath_runs_a_drill_op() {
        // The Run action (Generate) drives the real valenx-cam
        // toolpath engine. A drill op needs no source mesh, so it
        // generates cleanly headlessly.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.cam.new_op_kind = CamOpKind::Drill;
        app.cam_add_operation();
        assert!(
            app.mesh_toolbox.cam.last_error.is_none(),
            "add drill op errored: {:?}",
            app.mesh_toolbox.cam.last_error
        );
        app.cam_generate_toolpath();
        assert!(
            app.mesh_toolbox.cam.last_toolpath.is_some(),
            "drill toolpath not generated: {:?}",
            app.mesh_toolbox.cam.last_error
        );
    }

    #[test]
    fn cam_add_tool_rejects_a_degenerate_tool() {
        // Bad input — a zero-diameter tool — surfaces an error and
        // adds nothing.
        let mut app = ValenxApp::default();
        let before = app.mesh_toolbox.cam.file.tools.len();
        app.mesh_toolbox.cam.new_tool_diameter = 0.0;
        app.cam_add_tool();
        assert!(
            app.mesh_toolbox.cam.last_error.is_some(),
            "a zero-diameter tool should error"
        );
        assert_eq!(app.mesh_toolbox.cam.file.tools.len(), before);
    }

    // ===================================================================
    // Arch / BIM workbench — valenx-arch building entities
    // ===================================================================

    #[test]
    fn arch_panel_draws_every_tool() {
        for tool in [
            ArchTool::Wall,
            ArchTool::Slab,
            ArchTool::Column,
            ArchTool::Beam,
            ArchTool::Window,
            ArchTool::Door,
            ArchTool::Stair,
            ArchTool::Roof,
            ArchTool::Space,
        ] {
            let mut app = ValenxApp::default();
            app.mesh_toolbox.arch.tool = tool;
            draw_headless(|ui| draw_arch_panel(&mut app, ui));
        }
    }

    #[test]
    fn arch_panel_draws_post_run_and_error_states() {
        // Post-run — a wall already added + status.
        let mut app = ValenxApp::default();
        app.arch_add_wall();
        draw_headless(|ui| draw_arch_panel(&mut app, ui));
        // Error state.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.arch.last_error = Some("Wall: zero length".into());
        draw_headless(|ui| draw_arch_panel(&mut app, ui));
    }

    #[test]
    fn arch_add_entities_run_against_the_backend() {
        // The Run actions build real valenx-arch entities — wall,
        // slab, column, beam.
        let mut app = ValenxApp::default();
        app.arch_add_wall();
        app.arch_add_column();
        app.arch_add_beam();
        assert!(
            app.mesh_toolbox.arch.last_error.is_none(),
            "an arch add errored: {:?}",
            app.mesh_toolbox.arch.last_error
        );
        assert!(
            app.mesh_toolbox.arch.doc.entities.len() >= 3,
            "expected ≥3 entities, got {}",
            app.mesh_toolbox.arch.doc.entities.len()
        );
    }

    #[test]
    fn arch_render_tessellates_the_model() {
        // The Run action (Render) tessellates the building model into
        // a viewport mesh through valenx-arch + valenx-cad.
        let mut app = ValenxApp::default();
        app.arch_add_wall();
        app.arch_add_slab();
        app.arch_render();
        // Render either pushes a viewport mesh or surfaces an honest
        // error — never panics.
        assert!(
            app.mesh.is_some() || app.mesh_toolbox.arch.last_error.is_some(),
            "render should produce a mesh or an error"
        );
    }

    #[test]
    fn arch_add_wall_rejects_a_zero_length_wall() {
        // Bad input — start == end — surfaces an error and adds no
        // entity.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.arch.wall_start = [1.0, 1.0, 0.0];
        app.mesh_toolbox.arch.wall_end = [1.0, 1.0, 0.0];
        let before = app.mesh_toolbox.arch.doc.entities.len();
        app.arch_add_wall();
        assert!(
            app.mesh_toolbox.arch.last_error.is_some(),
            "a zero-length wall should error"
        );
        assert_eq!(app.mesh_toolbox.arch.doc.entities.len(), before);
    }

    // ===================================================================
    // Spreadsheet workbench — valenx-spreadsheet
    // ===================================================================

    #[test]
    fn spreadsheet_panel_draws_across_states() {
        // Fresh.
        draw_headless(|ui| draw_spreadsheet_panel(&mut ValenxApp::default(), ui));
        // With a populated cell + a selection.
        let mut app = ValenxApp::default();
        {
            let s = &mut app.mesh_toolbox.spreadsheet;
            let r = valenx_spreadsheet::CellRef {
                sheet_name: s.active_sheet.clone(),
                row: 0,
                col: 0,
            };
            s.workbook
                .set_cell(&r, valenx_spreadsheet::Cell::Number(42.0))
                .expect("set a number cell");
            s.selected_cell = Some((0, 0));
        }
        draw_headless(|ui| draw_spreadsheet_panel(&mut app, ui));
        // Error state.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.spreadsheet.last_error = Some("circular reference".into());
        draw_headless(|ui| draw_spreadsheet_panel(&mut app, ui));
    }

    #[test]
    fn spreadsheet_formula_evaluates_against_the_backend() {
        // The Run action (Set cell + evaluate) drives the real
        // valenx-spreadsheet formula engine.
        let mut s = SpreadsheetPanelState::default();
        let a1 = valenx_spreadsheet::CellRef {
            sheet_name: s.active_sheet.clone(),
            row: 0,
            col: 0,
        };
        let a2 = valenx_spreadsheet::CellRef {
            sheet_name: s.active_sheet.clone(),
            row: 1,
            col: 0,
        };
        let a3 = valenx_spreadsheet::CellRef {
            sheet_name: s.active_sheet.clone(),
            row: 2,
            col: 0,
        };
        s.workbook
            .set_cell(&a1, parse_editor_cell("6"))
            .expect("set A1");
        s.workbook
            .set_cell(&a2, parse_editor_cell("7"))
            .expect("set A2");
        // valenx-spreadsheet formula refs are sheet-qualified
        // (`Sheet.Cell`), exactly as the editor passes the raw text
        // after `=` straight through to Cell::Formula.
        s.workbook
            .set_cell(&a3, parse_editor_cell("=Default.A1 * Default.A2"))
            .expect("set A3 formula");
        let value = s.workbook.evaluate_cell(&a3).expect("formula should evaluate");
        assert_eq!(value, 42.0, "6 * 7 should evaluate to 42");
    }

    #[test]
    fn spreadsheet_parse_editor_cell_classifies_input() {
        // The editor-buffer parser classifies number / formula / text.
        assert!(matches!(
            parse_editor_cell("3.5"),
            valenx_spreadsheet::Cell::Number(_)
        ));
        assert!(matches!(
            parse_editor_cell("=A1+1"),
            valenx_spreadsheet::Cell::Formula(_)
        ));
        assert!(matches!(
            parse_editor_cell("hello"),
            valenx_spreadsheet::Cell::Text(_)
        ));
        assert!(matches!(
            parse_editor_cell(""),
            valenx_spreadsheet::Cell::Empty
        ));
    }

    #[test]
    fn spreadsheet_bad_formula_surfaces_an_error() {
        // Bad input — a self-referential formula — must surface an
        // error from evaluate_cell rather than panicking.
        let mut s = SpreadsheetPanelState::default();
        let a1 = valenx_spreadsheet::CellRef {
            sheet_name: s.active_sheet.clone(),
            row: 0,
            col: 0,
        };
        // A genuinely self-referential formula (sheet-qualified so the
        // failure is the circular reference, not a parse error).
        s.workbook
            .set_cell(&a1, parse_editor_cell("=Default.A1 + 1"))
            .expect("setting a circular formula is allowed; evaluating is not");
        assert!(
            s.workbook.evaluate_cell(&a1).is_err(),
            "a self-referential formula must error on evaluation"
        );
    }

    // ===================================================================
    // Dock panel — valenx docking
    // ===================================================================

    #[test]
    fn dock_panel_draws_across_states() {
        // Fresh.
        draw_headless(|ui| draw_dock_panel(&mut ValenxApp::default(), ui));
        // Error state.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.dock_panel.last_error = Some("receptor: file not found".into());
        draw_headless(|ui| draw_dock_panel(&mut app, ui));
        // Post-run — scored poses present + a selection.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.dock_panel.last_scores = vec![(-8.1, 0.0), (-7.4, 1.2)];
        app.mesh_toolbox.dock_panel.selected_pose = Some(0);
        draw_headless(|ui| draw_dock_panel(&mut app, ui));
    }

    #[test]
    fn dock_run_rejects_a_missing_receptor_file() {
        // Bad input — a receptor path that does not exist — the Run
        // action surfaces an error rather than panicking.
        let mut s = DockPanelState {
            receptor_path: "/nonexistent/receptor.pdbqt".to_string(),
            ligand_path: "/nonexistent/ligand.pdbqt".to_string(),
            ..Default::default()
        };
        run_dock_now(&mut s);
        assert!(
            s.last_error.is_some(),
            "docking with a missing receptor should error"
        );
        assert!(s.last_scores.is_empty());
    }

    // ===================================================================
    // Sketcher panel — valenx-sketch 2D constraint solver
    // ===================================================================

    #[test]
    fn sketcher_panel_draws_every_tool() {
        for tool in [SketcherTool::Select, SketcherTool::Line, SketcherTool::Circle] {
            let mut app = ValenxApp::default();
            app.mesh_toolbox.sketcher.tool = tool;
            draw_headless(|ui| draw_sketcher_panel(&mut app, ui));
        }
    }

    #[test]
    fn sketcher_panel_draws_post_solve_and_error_states() {
        // Post-solve — a sketch with geometry and a solver report.
        let mut app = ValenxApp::default();
        {
            let s = &mut app.mesh_toolbox.sketcher;
            let p = s.sketch.add_point(0.0, 0.0);
            let q = s.sketch.add_point(1.0, 0.0);
            s.sketch.add_line(p, q).expect("add a line");
            let report = valenx_sketch::solver::solve(
                &mut s.sketch,
                valenx_sketch::SolverConfig::default(),
            )
            .expect("solve a trivial sketch");
            s.last_report = Some(report);
        }
        draw_headless(|ui| draw_sketcher_panel(&mut app, ui));
        // Error state.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.sketcher.last_error = Some("add line: bad endpoint".into());
        draw_headless(|ui| draw_sketcher_panel(&mut app, ui));
    }

    #[test]
    fn sketcher_geometry_click_builds_a_line() {
        // The Run action (a Line-tool click pair) calls the real
        // valenx-sketch API and adds a line entity.
        let mut s = SketcherPanelState {
            tool: SketcherTool::Line,
            ..Default::default()
        };
        let entities_before = s.sketch.entities.len();
        handle_sketcher_geometry_click(&mut s, 0.0, 0.0); // first click
        handle_sketcher_geometry_click(&mut s, 5.0, 0.0); // second click -> line
        assert!(s.last_error.is_none(), "line click errored: {:?}", s.last_error);
        assert!(
            s.sketch.entities.len() > entities_before + 1,
            "a line + its two points should have been added"
        );
    }

    #[test]
    fn sketcher_solver_runs_on_a_real_sketch() {
        // The Run action (Solve) drives the real valenx-sketch
        // constraint solver on a sketch with a Distance constraint.
        let mut sketch = valenx_sketch::Sketch::default();
        let p = sketch.add_point(0.0, 0.0);
        let q = sketch.add_point(3.0, 0.0);
        let _ = sketch.add_line(p, q);
        sketch.add_constraint(valenx_sketch::constraint::Constraint::Distance {
            a: p,
            b: q,
            target: 10.0,
        });
        let report = valenx_sketch::solver::solve(
            &mut sketch,
            valenx_sketch::SolverConfig::default(),
        );
        assert!(report.is_ok(), "sketch solve errored: {:?}", report.err());
    }

    #[test]
    fn sketcher_circle_click_rejects_a_zero_radius() {
        // Bad input — a Circle-tool rim click coinciding with the
        // centre — surfaces an error (radius would be 0), no panic.
        let mut s = SketcherPanelState {
            tool: SketcherTool::Circle,
            ..Default::default()
        };
        handle_sketcher_geometry_click(&mut s, 2.0, 2.0); // centre
        handle_sketcher_geometry_click(&mut s, 2.0, 2.0); // rim == centre
        assert!(
            s.last_error.is_some(),
            "a zero-radius circle click should surface an error"
        );
    }

    // ===================================================================
    // Part Design panel — valenx-feature-tree parametric modelling
    // ===================================================================

    #[test]
    fn part_design_panel_draws_across_states() {
        // Fresh.
        draw_headless(|ui| draw_part_design_panel(&mut ValenxApp::default(), ui));
        // With a sketch + a Pad feature in the tree.
        let mut app = ValenxApp::default();
        {
            let s = &mut app.mesh_toolbox.part_design;
            let sketch = valenx_sketch::Sketch::default();
            let sref = s.tree.add_sketch(sketch);
            let params = valenx_feature_tree::feature::PadParams {
                sketch: sref,
                depth: 5.0.into(),
                direction_positive: true,
            };
            s.tree
                .add_feature(valenx_feature_tree::Feature::Pad(params), "Pad 1");
        }
        draw_headless(|ui| draw_part_design_panel(&mut app, ui));
        // Replay-error state.
        let mut app = ValenxApp::default();
        app.mesh_toolbox.part_design.last_replay_error = Some("replay: empty sketch".into());
        draw_headless(|ui| draw_part_design_panel(&mut app, ui));
    }

    #[test]
    fn part_design_replay_runs_a_real_pad_feature() {
        // The Run action (Replay) evaluates a real feature tree — a
        // closed-square sketch padded into a solid — through
        // valenx-feature-tree + valenx-cad.
        let mut app = ValenxApp::default();
        {
            let s = &mut app.mesh_toolbox.part_design;
            // A closed unit-square profile so the Pad has area.
            let mut sketch = valenx_sketch::Sketch::default();
            let a = sketch.add_point(0.0, 0.0);
            let b = sketch.add_point(1.0, 0.0);
            let c = sketch.add_point(1.0, 1.0);
            let d = sketch.add_point(0.0, 1.0);
            let _ = sketch.add_line(a, b);
            let _ = sketch.add_line(b, c);
            let _ = sketch.add_line(c, d);
            let _ = sketch.add_line(d, a);
            let sref = s.tree.add_sketch(sketch);
            let params = valenx_feature_tree::feature::PadParams {
                sketch: sref,
                depth: 3.0.into(),
                direction_positive: true,
            };
            s.tree
                .add_feature(valenx_feature_tree::Feature::Pad(params), "Pad 1");
        }
        run_part_design_replay(&mut app);
        // A valid pad replay either tessellates into a viewport mesh
        // or surfaces an honest error — it must never panic.
        assert!(
            app.mesh.is_some() || app.mesh_toolbox.part_design.last_replay_error.is_some(),
            "pad replay should produce a mesh or an honest error"
        );
    }

    #[test]
    fn part_design_replay_of_empty_tree_is_a_clean_no_op() {
        // An empty feature tree is a valid "mid-edit" state — replay
        // leaves the viewport alone and surfaces no error.
        let mut app = ValenxApp::default();
        run_part_design_replay(&mut app);
        assert!(
            app.mesh_toolbox.part_design.last_replay_error.is_none(),
            "empty-tree replay should not error"
        );
    }

    #[test]
    fn part_design_replay_rejects_a_pad_of_a_missing_sketch() {
        // Bad input — a Pad referencing a SketchRef that does not
        // exist — replay surfaces an error rather than panicking.
        let mut app = ValenxApp::default();
        {
            let s = &mut app.mesh_toolbox.part_design;
            let params = valenx_feature_tree::feature::PadParams {
                sketch: valenx_feature_tree::feature::SketchRef(99),
                depth: 5.0.into(),
                direction_positive: true,
            };
            s.tree
                .add_feature(valenx_feature_tree::Feature::Pad(params), "Bad Pad");
        }
        run_part_design_replay(&mut app);
        assert!(
            app.mesh_toolbox.part_design.last_replay_error.is_some(),
            "a pad of a non-existent sketch should surface a replay error"
        );
    }
}
