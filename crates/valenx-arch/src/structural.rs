//! Structural-analysis integration.
//!
//! Beams, columns, and slabs in an [`crate::ArchDocument`] carry an
//! optional [`StructuralMember`] payload — material grade + an applied
//! load — that lets [`export_structural_model`] convert a building's
//! frame into a `valenx-fem`-compatible problem definition: nodes
//! (from element end points), beam elements (one per structural beam
//! / column), materials (steel / concrete grade), loads (applied
//! forces and moments at element ends), and BCs (fixed / pinned
//! supports at ground level).
//!
//! The output is the [`StructuralModel`] struct — a freestanding
//! representation that callers can pipe into
//! `valenx_fem::solve_beam_static` / `valenx_fem::solve_beam_modal`
//! by translating the model's fields to the FEM solver's vectors.
//!
//! ## Honest scope
//!
//! - The element library is the **prismatic Timoshenko beam** that
//!   `valenx-fem` already ships — every beam / column maps to one
//!   beam element; an end-to-end portal-frame analysis is end-to-end
//!   real.
//! - **Slabs** are exported as a placeholder shell-element-like
//!   metadata row (`slab_count`); the v1 FEM solver does not assemble
//!   shells, and we honestly report rather than fabricate.
//! - Loads default to **self-weight as a node force** at each
//!   element's end when [`StructuralMember::self_weight_load`] is
//!   true; otherwise the user supplies the per-member applied force
//!   and moment.
//! - Supports come from one of three sources: an explicit
//!   `pinned`/`clamped` flag on a [`StructuralMember`]; an end node
//!   whose Z lies at or below the document's
//!   [`StructuralModelOptions::support_z`] (auto-ground); or no
//!   support at all (returns the model with `supports` empty, which
//!   the caller can validate).
//!
//! ## Example
//!
//! ```
//! use nalgebra::Vector3;
//! use valenx_arch::{ArchDocument, ArchEntity, BeamParams, BeamSection,
//!     ColumnParams, ColumnSection,
//!     structural::{export_structural_model, StructuralMember,
//!         StructuralMaterial, SupportKind, StructuralModelOptions}};
//!
//! // A 2-column / 1-beam portal frame.
//! let mut doc = ArchDocument::new("Portal");
//! // Left column.
//! doc.add_entity(ArchEntity::Column(ColumnParams {
//!     base: Vector3::new(0.0, 0.0, 0.0),
//!     height: 3.0,
//!     cross_section: ColumnSection::Rectangle { width: 0.3, depth: 0.3 },
//!     material: "Steel".into(),
//!     structural: Some(StructuralMember {
//!         material: StructuralMaterial::SteelS355,
//!         support: SupportKind::Clamped,
//!         applied_force: [0.0; 3],
//!         applied_moment: [0.0; 3],
//!         self_weight_load: false,
//!     }),
//! }));
//! // Right column.
//! doc.add_entity(ArchEntity::Column(ColumnParams {
//!     base: Vector3::new(5.0, 0.0, 0.0),
//!     height: 3.0,
//!     cross_section: ColumnSection::Rectangle { width: 0.3, depth: 0.3 },
//!     material: "Steel".into(),
//!     structural: Some(StructuralMember {
//!         material: StructuralMaterial::SteelS355,
//!         support: SupportKind::Clamped,
//!         applied_force: [0.0; 3],
//!         applied_moment: [0.0; 3],
//!         self_weight_load: false,
//!     }),
//! }));
//! // Crowning beam.
//! doc.add_entity(ArchEntity::Beam(BeamParams {
//!     start: Vector3::new(0.0, 0.0, 3.0),
//!     end: Vector3::new(5.0, 0.0, 3.0),
//!     cross_section: BeamSection::IBeam {
//!         width: 0.2, depth: 0.4,
//!         flange_thickness: 0.02, web_thickness: 0.01,
//!     },
//!     orientation_angle: 0.0,
//!     material: "Steel".into(),
//!     structural: Some(StructuralMember {
//!         material: StructuralMaterial::SteelS355,
//!         support: SupportKind::Free,
//!         applied_force: [0.0, 0.0, -10_000.0],
//!         applied_moment: [0.0; 3],
//!         self_weight_load: false,
//!     }),
//! }));
//!
//! let model = export_structural_model(&doc, &StructuralModelOptions::default()).unwrap();
//! assert_eq!(model.elements.len(), 3); // 2 columns + 1 beam
//! assert!(model.nodes.len() >= 4); // 4 unique end points (joined at midspan)
//! assert!(!model.supports.is_empty()); // 2 clamped column bases
//! assert!(!model.loads.is_empty()); // crown point load
//! ```

use std::collections::HashMap;

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::document::ArchDocument;
use crate::entity::ArchEntity;
use crate::error::ArchError;

/// Structural material grade — a curated set of common BIM materials
/// with characteristic mechanical properties.
///
/// Values are taken from Eurocode 3 (steel), Eurocode 2 (concrete),
/// and Eurocode 5 (timber) representative grades. Used by
/// [`StructuralMember`] to map a BIM-level material descriptor onto
/// the elastic constants the FEM beam solver consumes.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub enum StructuralMaterial {
    /// Structural steel — EN 10025-2 S235 (`fy = 235 MPa`).
    SteelS235,
    /// Structural steel — EN 10025-2 S355 (`fy = 355 MPa`, the most
    /// common European structural grade).
    SteelS355,
    /// Reinforced concrete — EN 1992 C25/30 (cylinder strength
    /// `fck = 25 MPa`).
    ConcreteC25,
    /// Reinforced concrete — EN 1992 C30/37.
    ConcreteC30,
    /// Glulam timber — EN 14080 GL24h (characteristic bending
    /// strength `fmk = 24 MPa`).
    TimberGL24,
}

impl StructuralMaterial {
    /// Short human label used in the schedule / IFC `Material` field.
    pub fn label(self) -> &'static str {
        match self {
            StructuralMaterial::SteelS235 => "Steel S235",
            StructuralMaterial::SteelS355 => "Steel S355",
            StructuralMaterial::ConcreteC25 => "Concrete C25/30",
            StructuralMaterial::ConcreteC30 => "Concrete C30/37",
            StructuralMaterial::TimberGL24 => "Timber GL24h",
        }
    }

    /// Young's modulus E in Pa.
    pub fn youngs_modulus(self) -> f64 {
        match self {
            StructuralMaterial::SteelS235 | StructuralMaterial::SteelS355 => 210.0e9,
            StructuralMaterial::ConcreteC25 => 31.0e9,
            StructuralMaterial::ConcreteC30 => 33.0e9,
            StructuralMaterial::TimberGL24 => 11.5e9,
        }
    }

    /// Poisson's ratio (dimensionless).
    pub fn poisson_ratio(self) -> f64 {
        match self {
            StructuralMaterial::SteelS235 | StructuralMaterial::SteelS355 => 0.30,
            StructuralMaterial::ConcreteC25 | StructuralMaterial::ConcreteC30 => 0.20,
            StructuralMaterial::TimberGL24 => 0.42,
        }
    }

    /// Density in kg/m³.
    pub fn density(self) -> f64 {
        match self {
            StructuralMaterial::SteelS235 | StructuralMaterial::SteelS355 => 7850.0,
            StructuralMaterial::ConcreteC25 | StructuralMaterial::ConcreteC30 => 2400.0,
            StructuralMaterial::TimberGL24 => 470.0,
        }
    }

    /// Characteristic yield / compressive strength in Pa — used by
    /// the IFC Pset_*Common `LoadBearing` / `MaterialGrade`
    /// attribution, not by the FEM solve.
    pub fn characteristic_strength(self) -> f64 {
        match self {
            StructuralMaterial::SteelS235 => 235.0e6,
            StructuralMaterial::SteelS355 => 355.0e6,
            StructuralMaterial::ConcreteC25 => 25.0e6,
            StructuralMaterial::ConcreteC30 => 30.0e6,
            StructuralMaterial::TimberGL24 => 24.0e6,
        }
    }
}

/// Support type at the structural member's "ground" end.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum SupportKind {
    /// No support — the member is free (typical for a beam whose
    /// columns provide the supports).
    Free,
    /// Pinned (translations fixed, rotations free).
    Pinned,
    /// Clamped / fully fixed (all 6 DOFs zero).
    Clamped,
}

impl SupportKind {
    /// Per-DOF fix mask `[ux,uy,uz, θx,θy,θz]` — `Some(0.0)` =
    /// fixed-to-zero, `None` = free.
    pub fn dof_mask(self) -> [Option<f64>; 6] {
        match self {
            SupportKind::Free => [None; 6],
            SupportKind::Pinned => [
                Some(0.0),
                Some(0.0),
                Some(0.0),
                None,
                None,
                None,
            ],
            SupportKind::Clamped => [Some(0.0); 6],
        }
    }
}

/// Structural attributes attached to a beam, column, or slab.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct StructuralMember {
    /// Material grade.
    pub material: StructuralMaterial,
    /// Support at the member's "ground" end (`base` for a column,
    /// `start` for a beam).
    pub support: SupportKind,
    /// Concentrated applied force `[Fx, Fy, Fz]` in newtons at the
    /// member's "tip" end (`base + height·ẑ` for a column, `end` for
    /// a beam).
    pub applied_force: [f64; 3],
    /// Concentrated applied moment `[Mx, My, Mz]` in newton-metres
    /// at the member's tip end.
    pub applied_moment: [f64; 3],
    /// When `true`, add an extra downward `[0, 0, -ρ·A·L·g]` node
    /// force at the tip end (a coarse self-weight lumping).
    pub self_weight_load: bool,
}

impl Default for StructuralMember {
    /// A free, unloaded steel-S355 default — useful for "make this a
    /// structural member but I'll set the rest later" UI calls.
    fn default() -> Self {
        Self {
            material: StructuralMaterial::SteelS355,
            support: SupportKind::Free,
            applied_force: [0.0; 3],
            applied_moment: [0.0; 3],
            self_weight_load: false,
        }
    }
}

/// One node of the exported structural model.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct StructuralNode {
    /// Position in world coordinates (metres).
    pub position: Vector3<f64>,
}

/// Cross-section properties for one structural element — the
/// inputs `valenx_fem::beam::BeamSection` needs.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct StructuralSection {
    /// Cross-section area `A` in m².
    pub area: f64,
    /// Second moment about the local-y axis `Iy` in m⁴.
    pub iy: f64,
    /// Second moment about the local-z axis `Iz` in m⁴.
    pub iz: f64,
    /// St-Venant torsion constant `J` in m⁴.
    pub j: f64,
}

/// One beam element of the exported model — a 2-node 3D beam.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct StructuralElement {
    /// Start-node index in [`StructuralModel::nodes`].
    pub start_node: usize,
    /// End-node index in [`StructuralModel::nodes`].
    pub end_node: usize,
    /// Cross-section properties.
    pub section: StructuralSection,
    /// Material grade.
    pub material: StructuralMaterial,
    /// Which arch entity id this element came from (back-ref into
    /// the [`crate::ArchDocument`]).
    pub source_entity_id: usize,
}

/// A single nodal support — one element of [`StructuralModel::supports`].
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct StructuralSupport {
    /// Node index.
    pub node: usize,
    /// Per-DOF fix mask `[ux,uy,uz, θx,θy,θz]`.
    pub fixed: [Option<f64>; 6],
}

/// A single nodal load — one element of [`StructuralModel::loads`].
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct StructuralLoad {
    /// Node index.
    pub node: usize,
    /// Force `[Fx,Fy,Fz]` in newtons.
    pub force: [f64; 3],
    /// Moment `[Mx,My,Mz]` in newton-metres.
    pub moment: [f64; 3],
}

/// The exported structural model — directly translatable to the
/// `valenx_fem::beam` solver's input vectors.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct StructuralModel {
    /// Unique nodes — element end points deduplicated by position.
    pub nodes: Vec<StructuralNode>,
    /// Beam elements (one per structural beam / column).
    pub elements: Vec<StructuralElement>,
    /// Nodal supports.
    pub supports: Vec<StructuralSupport>,
    /// Nodal loads.
    pub loads: Vec<StructuralLoad>,
    /// Material grades referenced by [`StructuralElement::material`],
    /// deduplicated for the IFC writer / structural-analysis input
    /// deck.
    pub materials: Vec<StructuralMaterial>,
    /// Count of slabs in the source document carrying structural
    /// attributes — emitted as a reporting field because the v1 FEM
    /// solver does not assemble shell elements (so we don't fabricate
    /// shell elements; we honestly carry the metadata).
    pub slab_count: usize,
}

impl StructuralModel {
    /// Total degree-of-freedom count (6 × nodes for a 3D beam model).
    pub fn dof_count(&self) -> usize {
        6 * self.nodes.len()
    }

    /// Count of constrained DOFs across all supports.
    pub fn constrained_dof_count(&self) -> usize {
        self.supports
            .iter()
            .map(|s| s.fixed.iter().filter(|f| f.is_some()).count())
            .sum()
    }
}

/// Tunable knobs for [`export_structural_model`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StructuralModelOptions {
    /// Position-merge tolerance in metres — two element end points
    /// closer than this are treated as a shared node.
    pub join_tolerance: f64,
    /// Optional auto-ground support: any element end point whose
    /// `z ≤ support_z` gets an implicit clamped support if the
    /// member's own `support` field is `Free`.
    pub support_z: Option<f64>,
    /// Acceleration due to gravity in m/s² — used by
    /// [`StructuralMember::self_weight_load`].
    pub gravity: f64,
}

impl Default for StructuralModelOptions {
    fn default() -> Self {
        Self {
            join_tolerance: 1.0e-6,
            support_z: None,
            gravity: 9.81,
        }
    }
}

/// Build a [`StructuralModel`] from every structural-tagged beam /
/// column / slab in `doc`.
///
/// Slabs contribute to [`StructuralModel::slab_count`] but do not
/// produce elements (see the honest-scope note in the module header).
/// Returns [`ArchError::BadDimension`] when an element's geometry is
/// degenerate.
pub fn export_structural_model(
    doc: &ArchDocument,
    opts: &StructuralModelOptions,
) -> Result<StructuralModel, ArchError> {
    let mut model = StructuralModel::default();
    let mut node_map: NodeMap = NodeMap::new(opts.join_tolerance);
    let mut materials: HashMap<StructuralMaterial, ()> = HashMap::new();

    for (eid, entity) in &doc.entities {
        match entity {
            ArchEntity::Beam(beam) => {
                if let Some(s) = &beam.structural {
                    add_beam_element(&mut model, &mut node_map, beam, s, *eid, opts)?;
                    materials.insert(s.material, ());
                }
            }
            ArchEntity::Column(col) => {
                if let Some(s) = &col.structural {
                    add_column_element(&mut model, &mut node_map, col, s, *eid, opts)?;
                    materials.insert(s.material, ());
                }
            }
            ArchEntity::Slab(slab) => {
                if slab.structural.is_some() {
                    model.slab_count += 1;
                    if let Some(s) = &slab.structural {
                        materials.insert(s.material, ());
                    }
                }
            }
            _ => {}
        }
    }

    // Stable order: sort materials by their string label so the
    // output is deterministic.
    let mut mats: Vec<StructuralMaterial> = materials.into_keys().collect();
    mats.sort_by_key(|m| m.label());
    model.materials = mats;

    Ok(model)
}

/// Translate a [`crate::beam::BeamSection`] into a
/// [`StructuralSection`] suitable for FEM.
fn section_from_beam(cs: &crate::beam::BeamSection) -> StructuralSection {
    match cs {
        crate::beam::BeamSection::Rectangle { width, depth } => {
            let a = width * depth;
            // Iz governs x-y plane bending (about local z) — depends on
            // the y-extent (here = width); standard `h·w³/12`.
            let iz = depth * width.powi(3) / 12.0;
            let iy = width * depth.powi(3) / 12.0;
            // St-Venant J for a thin rectangle (Roark).
            let (long, short) = if width >= depth {
                (*width, *depth)
            } else {
                (*depth, *width)
            };
            let a_h = long / 2.0;
            let b_h = short / 2.0;
            let ratio = b_h / a_h;
            let j = a_h * b_h.powi(3) * (16.0 / 3.0 - 3.36 * ratio * (1.0 - ratio.powi(4) / 12.0));
            StructuralSection {
                area: a,
                iy,
                iz,
                j,
            }
        }
        crate::beam::BeamSection::IBeam {
            width,
            depth,
            flange_thickness,
            web_thickness,
        } => {
            // True I-section properties (web depth = depth - 2 t_f).
            let inner_d = (depth - 2.0 * flange_thickness).max(0.0);
            let a = 2.0 * width * flange_thickness + inner_d * web_thickness;
            // Iy about the local y axis (along the web — strong axis
            // for an upright I). Parallel-axis sum of the two flanges
            // (each at z ≈ ±(d/2 − t_f/2)) and the web.
            let flange_z = depth * 0.5 - flange_thickness * 0.5;
            let iy_flange =
                width * flange_thickness.powi(3) / 12.0 + width * flange_thickness * flange_z.powi(2);
            let iy_web = web_thickness * inner_d.powi(3) / 12.0;
            let iy = 2.0 * iy_flange + iy_web;
            // Iz about the weak axis (perpendicular to the web).
            let iz_flange = flange_thickness * width.powi(3) / 12.0;
            let iz_web = inner_d * web_thickness.powi(3) / 12.0;
            let iz = 2.0 * iz_flange + iz_web;
            // J ≈ Σ b·t³/3 over rectangles (open thin-walled).
            let j = (2.0 * width * flange_thickness.powi(3) + inner_d * web_thickness.powi(3)) / 3.0;
            StructuralSection {
                area: a.max(1.0e-12),
                iy: iy.max(1.0e-12),
                iz: iz.max(1.0e-12),
                j: j.max(1.0e-12),
            }
        }
        crate::beam::BeamSection::Channel {
            width,
            depth,
            thickness,
        } => {
            // Treat as three rectangles: web (h × t) + 2 flanges
            // ((w - t) × t).
            let inner_w = (width - thickness).max(0.0);
            let a = depth * thickness + 2.0 * inner_w * thickness;
            // Iy about the symmetry axis through the web.
            let iy_web = thickness * depth.powi(3) / 12.0;
            let flange_z = depth * 0.5 - thickness * 0.5;
            let iy_flange =
                inner_w * thickness.powi(3) / 12.0 + inner_w * thickness * flange_z.powi(2);
            let iy = iy_web + 2.0 * iy_flange;
            let iz_web = depth * thickness.powi(3) / 12.0;
            let iz = iz_web + 2.0 * (thickness * inner_w.powi(3) / 12.0);
            let j = (depth * thickness.powi(3) + 2.0 * inner_w * thickness.powi(3)) / 3.0;
            StructuralSection {
                area: a.max(1.0e-12),
                iy: iy.max(1.0e-12),
                iz: iz.max(1.0e-12),
                j: j.max(1.0e-12),
            }
        }
    }
}

/// Translate a [`crate::column::ColumnSection`] into a
/// [`StructuralSection`].
fn section_from_column(cs: &crate::column::ColumnSection) -> StructuralSection {
    match cs {
        crate::column::ColumnSection::Rectangle { width, depth } => {
            let a = width * depth;
            let iz = depth * width.powi(3) / 12.0;
            let iy = width * depth.powi(3) / 12.0;
            let (long, short) = if width >= depth {
                (*width, *depth)
            } else {
                (*depth, *width)
            };
            let a_h = long / 2.0;
            let b_h = short / 2.0;
            let ratio = b_h / a_h;
            let j = a_h * b_h.powi(3) * (16.0 / 3.0 - 3.36 * ratio * (1.0 - ratio.powi(4) / 12.0));
            StructuralSection {
                area: a,
                iy,
                iz,
                j,
            }
        }
        crate::column::ColumnSection::Circular { radius, .. } => {
            let a = std::f64::consts::PI * radius * radius;
            let i = std::f64::consts::PI * radius.powi(4) / 4.0;
            StructuralSection {
                area: a,
                iy: i,
                iz: i,
                j: 2.0 * i,
            }
        }
        crate::column::ColumnSection::IBeam {
            width,
            depth,
            flange_thickness,
            web_thickness,
        } => {
            let inner_d = (depth - 2.0 * flange_thickness).max(0.0);
            let a = 2.0 * width * flange_thickness + inner_d * web_thickness;
            let flange_z = depth * 0.5 - flange_thickness * 0.5;
            let iy_flange =
                width * flange_thickness.powi(3) / 12.0 + width * flange_thickness * flange_z.powi(2);
            let iy_web = web_thickness * inner_d.powi(3) / 12.0;
            let iy = 2.0 * iy_flange + iy_web;
            let iz_flange = flange_thickness * width.powi(3) / 12.0;
            let iz_web = inner_d * web_thickness.powi(3) / 12.0;
            let iz = 2.0 * iz_flange + iz_web;
            let j = (2.0 * width * flange_thickness.powi(3) + inner_d * web_thickness.powi(3)) / 3.0;
            StructuralSection {
                area: a.max(1.0e-12),
                iy: iy.max(1.0e-12),
                iz: iz.max(1.0e-12),
                j: j.max(1.0e-12),
            }
        }
    }
}

/// Merge structural-element end points by position so two members
/// joining at a corner share one node.
struct NodeMap {
    tol: f64,
}

impl NodeMap {
    fn new(tol: f64) -> Self {
        Self { tol: tol.max(0.0) }
    }

    /// Insert `p` (or return the existing index if a node within
    /// `self.tol` already exists).
    fn intern(&mut self, model: &mut StructuralModel, p: Vector3<f64>) -> usize {
        for (i, n) in model.nodes.iter().enumerate() {
            if (n.position - p).norm() <= self.tol {
                return i;
            }
        }
        model.nodes.push(StructuralNode { position: p });
        model.nodes.len() - 1
    }
}

fn add_beam_element(
    model: &mut StructuralModel,
    nodes: &mut NodeMap,
    beam: &crate::beam::BeamParams,
    s: &StructuralMember,
    source_id: usize,
    opts: &StructuralModelOptions,
) -> Result<(), ArchError> {
    beam.validate()?;
    let n_start = nodes.intern(model, beam.start);
    let n_end = nodes.intern(model, beam.end);
    let section = section_from_beam(&beam.cross_section);

    model.elements.push(StructuralElement {
        start_node: n_start,
        end_node: n_end,
        section,
        material: s.material,
        source_entity_id: source_id,
    });

    apply_support(model, n_start, s);
    apply_loads(model, n_end, s, section, beam.length(), opts);
    apply_ground_support(model, n_start, opts);
    apply_ground_support(model, n_end, opts);

    Ok(())
}

fn add_column_element(
    model: &mut StructuralModel,
    nodes: &mut NodeMap,
    col: &crate::column::ColumnParams,
    s: &StructuralMember,
    source_id: usize,
    opts: &StructuralModelOptions,
) -> Result<(), ArchError> {
    col.validate()?;
    let bottom = col.base;
    let top = col.base + Vector3::new(0.0, 0.0, col.height);
    let n_start = nodes.intern(model, bottom);
    let n_end = nodes.intern(model, top);
    let section = section_from_column(&col.cross_section);

    model.elements.push(StructuralElement {
        start_node: n_start,
        end_node: n_end,
        section,
        material: s.material,
        source_entity_id: source_id,
    });

    apply_support(model, n_start, s);
    apply_loads(model, n_end, s, section, col.height, opts);
    apply_ground_support(model, n_start, opts);
    apply_ground_support(model, n_end, opts);

    Ok(())
}

/// Apply the explicit `StructuralMember::support` to `node`. Skips
/// `SupportKind::Free`. Merges with an existing support at the same
/// node by OR-ing fixed DOFs (a clamp wins over a pin).
fn apply_support(model: &mut StructuralModel, node: usize, s: &StructuralMember) {
    let mask = s.support.dof_mask();
    if mask.iter().all(Option::is_none) {
        return;
    }
    if let Some(existing) = model.supports.iter_mut().find(|sp| sp.node == node) {
        for (i, dof) in existing.fixed.iter_mut().enumerate() {
            if dof.is_none() {
                *dof = mask[i];
            }
        }
    } else {
        model.supports.push(StructuralSupport { node, fixed: mask });
    }
}

/// Auto-ground: if `opts.support_z` is set and `node.z ≤ support_z`
/// and no support already exists at this node, clamp it.
fn apply_ground_support(model: &mut StructuralModel, node: usize, opts: &StructuralModelOptions) {
    let Some(z_g) = opts.support_z else {
        return;
    };
    let z = model.nodes[node].position.z;
    if z > z_g {
        return;
    }
    if model.supports.iter().any(|sp| sp.node == node) {
        return;
    }
    model.supports.push(StructuralSupport {
        node,
        fixed: SupportKind::Clamped.dof_mask(),
    });
}

/// Apply the member's `applied_force` + `applied_moment` plus
/// `self_weight_load` (if enabled) to `node`. A repeat node sums
/// loads, mirroring the FEM solver's nodal-force superposition.
fn apply_loads(
    model: &mut StructuralModel,
    node: usize,
    s: &StructuralMember,
    section: StructuralSection,
    length: f64,
    opts: &StructuralModelOptions,
) {
    let mut force = s.applied_force;
    let moment = s.applied_moment;
    if s.self_weight_load {
        let weight = s.material.density() * section.area * length * opts.gravity;
        force[2] -= weight;
    }
    if force == [0.0; 3] && moment == [0.0; 3] {
        return;
    }
    if let Some(existing) = model.loads.iter_mut().find(|l| l.node == node) {
        for i in 0..3 {
            existing.force[i] += force[i];
            existing.moment[i] += moment[i];
        }
    } else {
        model.loads.push(StructuralLoad {
            node,
            force,
            moment,
        });
    }
}

/// LRFD factored load combination per ASCE 7 §2.3 strength design: `1.2·D + 1.6·L`,
/// combining a dead (permanent) load `dead_load` with a live (transient) load `live_load`
/// into the ultimate-strength design demand. Units pass through unchanged (N or Pa) and the
/// sign of each load is preserved (uplift / negative permitted). Returns `0` for non-finite input.
pub fn lrfd_factored_load(dead_load: f64, live_load: f64) -> f64 {
    if !dead_load.is_finite() || !live_load.is_finite() {
        return 0.0;
    }
    1.2 * dead_load + 1.6 * live_load
}

/// ASD load combination per ASCE 7 §2.4 basic combination 1: `D + L`, combining a dead
/// (permanent) load `dead_load` with a live (transient) load `live_load` into the
/// allowable-stress-design service demand. Unfactored (unity factors), so the result is
/// always ≤ the LRFD [`lrfd_factored_load`] for the same non-negative inputs. Returns `0`
/// for non-finite input; the sign of each load is preserved.
pub fn asd_load_combination(dead_load: f64, live_load: f64) -> f64 {
    if !dead_load.is_finite() || !live_load.is_finite() {
        return 0.0;
    }
    dead_load + live_load
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::beam::{BeamParams, BeamSection};
    use crate::column::{ColumnParams, ColumnSection};
    use crate::entity::ArchEntity;
    use crate::slab::SlabParams;

    fn portal_frame() -> ArchDocument {
        let mut d = ArchDocument::new("Portal");
        // Left column at x=0.
        d.add_entity(ArchEntity::Column(ColumnParams {
            base: Vector3::new(0.0, 0.0, 0.0),
            height: 3.0,
            cross_section: ColumnSection::Rectangle {
                width: 0.3,
                depth: 0.3,
            },
            material: "Steel".into(),
            structural: Some(StructuralMember {
                material: StructuralMaterial::SteelS355,
                support: SupportKind::Clamped,
                applied_force: [0.0; 3],
                applied_moment: [0.0; 3],
                self_weight_load: false,
            }),
        }));
        // Right column at x=5.
        d.add_entity(ArchEntity::Column(ColumnParams {
            base: Vector3::new(5.0, 0.0, 0.0),
            height: 3.0,
            cross_section: ColumnSection::Rectangle {
                width: 0.3,
                depth: 0.3,
            },
            material: "Steel".into(),
            structural: Some(StructuralMember {
                material: StructuralMaterial::SteelS355,
                support: SupportKind::Clamped,
                applied_force: [0.0; 3],
                applied_moment: [0.0; 3],
                self_weight_load: false,
            }),
        }));
        // Beam spanning the columns at z=3.
        d.add_entity(ArchEntity::Beam(BeamParams {
            start: Vector3::new(0.0, 0.0, 3.0),
            end: Vector3::new(5.0, 0.0, 3.0),
            cross_section: BeamSection::IBeam {
                width: 0.2,
                depth: 0.4,
                flange_thickness: 0.02,
                web_thickness: 0.01,
            },
            orientation_angle: 0.0,
            material: "Steel".into(),
            structural: Some(StructuralMember {
                material: StructuralMaterial::SteelS355,
                support: SupportKind::Free,
                applied_force: [0.0, 0.0, -10_000.0],
                applied_moment: [0.0; 3],
                self_weight_load: false,
            }),
        }));
        d
    }

    #[test]
    fn portal_frame_three_elements_four_nodes() {
        let doc = portal_frame();
        let m = export_structural_model(&doc, &StructuralModelOptions::default()).unwrap();
        // 2 columns + 1 beam.
        assert_eq!(m.elements.len(), 3);
        // 4 unique end points: (0,0,0), (5,0,0), (0,0,3), (5,0,3).
        assert_eq!(m.nodes.len(), 4);
        // The 2 column bases are clamped, so 12 DOFs constrained.
        assert_eq!(m.supports.len(), 2);
        assert_eq!(m.constrained_dof_count(), 12);
        assert_eq!(m.dof_count(), 24);
        // One node load (crown point load at one of the top corners).
        assert_eq!(m.loads.len(), 1);
        // Materials deduplicated to one entry.
        assert_eq!(m.materials.len(), 1);
        assert_eq!(m.materials[0], StructuralMaterial::SteelS355);
    }

    #[test]
    fn portal_frame_top_nodes_are_joined() {
        let doc = portal_frame();
        let m = export_structural_model(&doc, &StructuralModelOptions::default()).unwrap();
        // The beam ends at (0,0,3) and (5,0,3); each column's top is
        // the same point. So the column-beam joint must dedupe.
        let positions: Vec<Vector3<f64>> = m.nodes.iter().map(|n| n.position).collect();
        let crown_count = positions
            .iter()
            .filter(|p| (p.z - 3.0).abs() < 1e-9)
            .count();
        assert_eq!(crown_count, 2, "got crown positions {positions:?}");
    }

    #[test]
    fn ignores_non_structural_members() {
        let mut d = ArchDocument::new("Test");
        d.add_entity(ArchEntity::Beam(BeamParams {
            start: Vector3::zeros(),
            end: Vector3::new(2.0, 0.0, 0.0),
            cross_section: BeamSection::Rectangle {
                width: 0.1,
                depth: 0.1,
            },
            orientation_angle: 0.0,
            material: "Wood".into(),
            structural: None, // not structural — should be skipped.
        }));
        let m = export_structural_model(&d, &StructuralModelOptions::default()).unwrap();
        assert!(m.elements.is_empty());
        assert!(m.nodes.is_empty());
    }

    #[test]
    fn auto_ground_support_pulls_in_columns_without_explicit_support() {
        let mut d = ArchDocument::new("Auto");
        d.add_entity(ArchEntity::Column(ColumnParams {
            base: Vector3::new(0.0, 0.0, 0.0),
            height: 3.0,
            cross_section: ColumnSection::Circular {
                radius: 0.15,
                segments: 12,
            },
            material: "Steel".into(),
            structural: Some(StructuralMember {
                support: SupportKind::Free, // no explicit support.
                ..StructuralMember::default()
            }),
        }));
        let opts = StructuralModelOptions {
            support_z: Some(0.001),
            ..StructuralModelOptions::default()
        };
        let m = export_structural_model(&d, &opts).unwrap();
        // The base at z=0 picks up an auto-ground clamp.
        assert_eq!(m.supports.len(), 1);
        assert_eq!(m.constrained_dof_count(), 6);
    }

    #[test]
    fn self_weight_adds_negative_z_force() {
        let mut d = ArchDocument::new("SelfWeight");
        d.add_entity(ArchEntity::Beam(BeamParams {
            start: Vector3::new(0.0, 0.0, 3.0),
            end: Vector3::new(2.0, 0.0, 3.0),
            cross_section: BeamSection::Rectangle {
                width: 0.1,
                depth: 0.1,
            },
            orientation_angle: 0.0,
            material: "Steel".into(),
            structural: Some(StructuralMember {
                material: StructuralMaterial::SteelS355,
                support: SupportKind::Free,
                applied_force: [0.0; 3],
                applied_moment: [0.0; 3],
                self_weight_load: true,
            }),
        }));
        let m = export_structural_model(&d, &StructuralModelOptions::default()).unwrap();
        assert_eq!(m.loads.len(), 1);
        let l = m.loads[0];
        assert!(l.force[2] < 0.0, "got self-weight load {l:?}");
        // ρ·A·L·g = 7850 · 0.01 · 2 · 9.81 ≈ 1540 N (downward).
        let expected = 7850.0 * 0.01 * 2.0 * 9.81;
        assert!(
            (l.force[2].abs() - expected).abs() / expected < 1e-6,
            "load {} not close to {}",
            l.force[2],
            -expected
        );
    }

    #[test]
    fn slab_count_includes_structural_slabs_without_emitting_elements() {
        let mut d = ArchDocument::new("Slabs");
        d.add_entity(ArchEntity::Slab(SlabParams {
            boundary: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(2.0, 0.0, 0.0),
                Vector3::new(2.0, 2.0, 0.0),
                Vector3::new(0.0, 2.0, 0.0),
            ],
            thickness: 0.2,
            material: "Concrete".into(),
            structural: Some(StructuralMember {
                material: StructuralMaterial::ConcreteC25,
                support: SupportKind::Pinned,
                applied_force: [0.0; 3],
                applied_moment: [0.0; 3],
                self_weight_load: false,
            }),
        }));
        let m = export_structural_model(&d, &StructuralModelOptions::default()).unwrap();
        assert_eq!(m.slab_count, 1);
        assert!(m.elements.is_empty());
        assert!(m.materials.contains(&StructuralMaterial::ConcreteC25));
    }

    #[test]
    fn material_grade_constants_are_finite_positive() {
        for grade in [
            StructuralMaterial::SteelS235,
            StructuralMaterial::SteelS355,
            StructuralMaterial::ConcreteC25,
            StructuralMaterial::ConcreteC30,
            StructuralMaterial::TimberGL24,
        ] {
            assert!(grade.youngs_modulus() > 0.0);
            assert!(grade.density() > 0.0);
            assert!(grade.characteristic_strength() > 0.0);
            assert!(grade.poisson_ratio() > -1.0 && grade.poisson_ratio() < 0.5);
            assert!(!grade.label().is_empty());
        }
    }

    #[test]
    fn section_from_rectangle_matches_handbook() {
        let cs = crate::beam::BeamSection::Rectangle {
            width: 0.2,
            depth: 0.3,
        };
        let s = section_from_beam(&cs);
        assert!((s.area - 0.06).abs() < 1e-9);
        // Iy = w·h³/12 = 0.2·0.027/12 = 4.5e-4
        assert!((s.iy - 0.2 * 0.027 / 12.0).abs() < 1e-12);
        // Iz = h·w³/12 = 0.3·0.008/12 = 2.0e-4
        assert!((s.iz - 0.3 * 0.008 / 12.0).abs() < 1e-12);
        assert!(s.j > 0.0);
    }

    #[test]
    fn lrfd_factored_load_basic() {
        // ASCE 7 §2.3 basic strength combo: 1.2·D + 1.6·L.
        assert!((lrfd_factored_load(10.0, 5.0) - 20.0).abs() < 1e-9); // 12 + 8
        assert!((lrfd_factored_load(10.0, 0.0) - 12.0).abs() < 1e-9); // dead only
        assert!((lrfd_factored_load(0.0, 5.0) - 8.0).abs() < 1e-9); // live only
        // Factors are NOT swapped: 1.2·10 + 1.6·5 = 20 ≠ 1.6·10 + 1.2·5 = 22.
        assert!((lrfd_factored_load(10.0, 5.0) - (1.6 * 10.0 + 1.2 * 5.0)).abs() > 1.0);
        // Sign preserved (uplift); non-finite → 0.
        assert!((lrfd_factored_load(-5.0, 3.0) - (-1.2)).abs() < 1e-9);
        assert_eq!(lrfd_factored_load(f64::NAN, 5.0), 0.0);
        assert_eq!(lrfd_factored_load(10.0, f64::INFINITY), 0.0);
    }

    #[test]
    fn asd_load_combination_basic() {
        // ASCE 7 §2.4 ASD basic combo: D + L (unity factors).
        assert!((asd_load_combination(10.0, 5.0) - 15.0).abs() < 1e-9);
        assert!((asd_load_combination(10.0, 0.0) - 10.0).abs() < 1e-9); // dead only
        assert!((asd_load_combination(0.0, 5.0) - 5.0).abs() < 1e-9); // live only
        // Non-tautological: ASD < LRFD for the same non-negative inputs (15 < 20).
        assert!(asd_load_combination(10.0, 5.0) < lrfd_factored_load(10.0, 5.0));
        // Sign preserved (uplift); non-finite → 0.
        assert!((asd_load_combination(-5.0, 3.0) - (-2.0)).abs() < 1e-9);
        assert_eq!(asd_load_combination(f64::NAN, 5.0), 0.0);
        assert_eq!(asd_load_combination(10.0, f64::INFINITY), 0.0);
    }
}
