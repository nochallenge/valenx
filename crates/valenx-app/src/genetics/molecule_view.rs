//! Molecular 3D-viewport integration for the Genetics workbench.
//!
//! The thirteen Genetics panels render structures and molecules as
//! text and tables. This module bridges three of the bio-crate data
//! models —
//!
//! - [`valenx_biostruct::structure::Structure`] (a macromolecular
//!   PDB / mmCIF structure),
//! - [`valenx_cheminf::Molecule`] (a small-molecule graph carrying a
//!   3-D conformer),
//! - a [`valenx_md::System`] frame (an MD snapshot),
//!
//! — into a [`valenx_viz::TriangleMesh`] the app's wgpu 3-D viewport
//! already knows how to draw. Two representations are produced:
//!
//! - [`ball_and_stick`] — a small sphere at every atom (radius scaled
//!   from the element van-der-Waals radius) plus a cylinder for every
//!   bond,
//! - [`spacefill`] — full van-der-Waals spheres, no bonds (the CPK /
//!   space-filling model).
//!
//! A panel pushes the resulting mesh into the viewport with
//! [`show_molecule`], which sets `ValenxApp::stl` exactly as the STL
//! importer does so the viewport renders it like any other triangle
//! soup.
//!
//! ## Honest scope
//!
//! Ball-and-stick + spacefill *geometry* is the shipped deliverable:
//! correctly sized element spheres and bond cylinders, framed in the
//! existing wgpu 3-D viewport. The standard CPK element palette is
//! computed by [`element_color`] and is what an STL/colour-aware
//! consumer would use, but the app's current shaded viewport renders
//! every triangle with one material (the [`valenx_viz::TriangleMesh`]
//! the viewport consumes carries no per-triangle colour, and adding a
//! colour channel would mean changing the shared CAD rendering
//! pipeline) — so on screen the molecule appears as a single-material
//! shaded ball-and-stick model. Per-atom colour in the viewport, and
//! protein cartoon / ribbon rendering (a spline through the Cα trace,
//! sheet arrows, helix tubes), are documented follow-ons.

use std::path::PathBuf;

use valenx_viz::stl::{StlTriangle, TriangleMesh};

use crate::types::LoadedStl;
use crate::ValenxApp;

/// One atom reduced to what the mesh builder needs: a position
/// (ångström) and an element symbol.
#[derive(Clone, Debug, PartialEq)]
pub struct ViewAtom {
    /// Cartesian position in ångström.
    pub pos: [f32; 3],
    /// Upper-cased element symbol (`"C"`, `"N"`, `"FE"`, …). An empty
    /// or unrecognised symbol falls back to a neutral grey carbon-like
    /// atom.
    pub element: String,
}

impl ViewAtom {
    /// A view atom at `pos` of element `element`.
    pub fn new(pos: [f32; 3], element: impl Into<String>) -> Self {
        ViewAtom {
            pos,
            element: element.into(),
        }
    }
}

/// A molecule reduced to atoms + bonds, ready for meshing. Bonds are
/// index pairs into `atoms`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ViewMolecule {
    /// Every atom.
    pub atoms: Vec<ViewAtom>,
    /// Bonds as `(atom_a, atom_b)` index pairs.
    pub bonds: Vec<(usize, usize)>,
}

impl ViewMolecule {
    /// An empty molecule.
    pub fn new() -> Self {
        ViewMolecule::default()
    }

    /// Whether the molecule has no atoms.
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Build a [`ViewMolecule`] from a [`valenx_biostruct`] structure.
    ///
    /// Every atom of the first model becomes a [`ViewAtom`]; bonds are
    /// inferred by the covalent-radius distance rule
    /// ([`detect_bonds`]) since a PDB structure carries no explicit
    /// heavy-atom connectivity in this crate's model.
    pub fn from_biostruct(s: &valenx_biostruct::structure::Structure) -> Self {
        let mut atoms = Vec::new();
        for chain in &s.first_model().chains {
            for res in &chain.residues {
                for atom in &res.atoms {
                    atoms.push(ViewAtom::new(
                        [
                            atom.coord.x as f32,
                            atom.coord.y as f32,
                            atom.coord.z as f32,
                        ],
                        atom.element.trim().to_ascii_uppercase(),
                    ));
                }
            }
        }
        let bonds = detect_bonds(&atoms);
        ViewMolecule { atoms, bonds }
    }

    /// Build a [`ViewMolecule`] from a [`valenx_cheminf`] molecule that
    /// already carries a 3-D conformer (`coords_3d == true`).
    ///
    /// Returns `None` if the molecule has no atoms or its coordinates
    /// are absent / only a 2-D depiction — the caller should generate a
    /// conformer (`valenx_cheminf::coords::embed_3d`) first.
    pub fn from_cheminf(mol: &valenx_cheminf::Molecule) -> Option<Self> {
        if mol.atoms.is_empty() || !mol.coords_3d || mol.coords.len() != mol.atoms.len() {
            return None;
        }
        let atoms = mol
            .atoms
            .iter()
            .zip(&mol.coords)
            .map(|(a, c)| {
                ViewAtom::new(
                    [c[0] as f32, c[1] as f32, c[2] as f32],
                    a.symbol().to_ascii_uppercase(),
                )
            })
            .collect();
        // The cheminf graph carries explicit bonds — use them directly.
        let bonds = mol.bonds.iter().map(|b| (b.a, b.b)).collect();
        Some(ViewMolecule { atoms, bonds })
    }

    /// Build a [`ViewMolecule`] from a [`valenx_md`] system frame.
    ///
    /// MD positions are in nanometres; they are converted to ångström
    /// here so the mesh shares the biostruct / cheminf length scale.
    /// Bonds come from the system topology if it has any, else they are
    /// inferred by [`detect_bonds`].
    pub fn from_md_system(system: &valenx_md::System) -> Self {
        const ANGSTROM_PER_NM: f32 = 10.0;
        let atoms: Vec<ViewAtom> = system
            .topology
            .atoms
            .iter()
            .zip(&system.positions)
            .map(|(a, p)| {
                let element = if a.element.is_empty() {
                    a.type_name.clone()
                } else {
                    a.element.clone()
                };
                ViewAtom::new(
                    [
                        p.x as f32 * ANGSTROM_PER_NM,
                        p.y as f32 * ANGSTROM_PER_NM,
                        p.z as f32 * ANGSTROM_PER_NM,
                    ],
                    element.trim().to_ascii_uppercase(),
                )
            })
            .collect();
        let bonds = if system.topology.bonds.is_empty() {
            detect_bonds(&atoms)
        } else {
            system.topology.bonds.iter().map(|b| (b.i, b.j)).collect()
        };
        ViewMolecule { atoms, bonds }
    }
}

/// Covalent radius (ångström) for an element symbol — Cordero 2008
/// values for the common biomolecular elements, a carbon-ish `0.75`
/// fallback. Used by [`detect_bonds`].
pub fn covalent_radius(element: &str) -> f32 {
    match element.trim().to_ascii_uppercase().as_str() {
        "H" | "D" => 0.31,
        "B" => 0.84,
        "C" => 0.76,
        "N" => 0.71,
        "O" => 0.66,
        "F" => 0.57,
        "NA" => 1.66,
        "MG" => 1.41,
        "SI" => 1.11,
        "P" => 1.07,
        "S" => 1.05,
        "CL" => 1.02,
        "K" => 2.03,
        "CA" => 1.76,
        "FE" => 1.32,
        "ZN" => 1.22,
        "BR" => 1.20,
        "I" => 1.39,
        _ => 0.75,
    }
}

/// Van-der-Waals radius (ångström) for an element symbol — Bondi 1964
/// values for the common elements, a `1.7` carbon-ish fallback. Drives
/// the sphere size in both representations.
pub fn vdw_radius(element: &str) -> f32 {
    match element.trim().to_ascii_uppercase().as_str() {
        "H" | "D" => 1.20,
        "B" => 1.92,
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "F" => 1.47,
        "NA" => 2.27,
        "MG" => 1.73,
        "SI" => 2.10,
        "P" => 1.80,
        "S" => 1.80,
        "CL" => 1.75,
        "K" => 2.75,
        "CA" => 2.31,
        "FE" => 2.04,
        "ZN" => 2.10,
        "BR" => 1.85,
        "I" => 1.98,
        _ => 1.70,
    }
}

/// CPK element colour as linear-ish `[r, g, b]` in `0.0..=1.0`. Maps
/// the common biomolecular elements to the conventional palette;
/// unrecognised elements get a soft pink (the CPK "other" colour).
pub fn element_color(element: &str) -> [f32; 3] {
    match element.trim().to_ascii_uppercase().as_str() {
        "H" | "D" => [0.95, 0.95, 0.95],
        "C" => [0.30, 0.30, 0.32],
        "N" => [0.19, 0.31, 0.97],
        "O" => [0.94, 0.15, 0.10],
        "F" | "CL" => [0.22, 0.86, 0.22],
        "BR" => [0.59, 0.21, 0.14],
        "I" => [0.39, 0.16, 0.55],
        "S" => [0.94, 0.78, 0.20],
        "P" => [0.96, 0.55, 0.16],
        "NA" => [0.55, 0.36, 0.92],
        "MG" => [0.13, 0.54, 0.13],
        "K" => [0.50, 0.22, 0.78],
        "CA" => [0.24, 0.78, 0.78],
        "FE" => [0.88, 0.40, 0.05],
        "ZN" => [0.49, 0.50, 0.69],
        "B" => [1.0, 0.71, 0.71],
        _ => [1.0, 0.41, 0.71],
    }
}

/// Detect bonds from interatomic distances.
///
/// Two atoms are bonded when their separation is within
/// `(r_cov_a + r_cov_b + TOLERANCE)` — the standard covalent-radius
/// rule used by every molecular viewer to add connectivity to a
/// bond-free structure. A `0.45 Å` slack absorbs coordinate noise and
/// the radius table's approximations; an absolute `0.4 Å` floor drops
/// pathological near-coincident atoms.
///
/// `O(n²)`; for the structure sizes a desktop panel meshes (a few
/// thousand atoms) that is comfortably fast.
pub fn detect_bonds(atoms: &[ViewAtom]) -> Vec<(usize, usize)> {
    const TOLERANCE: f32 = 0.45;
    const MIN_DIST: f32 = 0.4;
    let mut bonds = Vec::new();
    for i in 0..atoms.len() {
        let ri = covalent_radius(&atoms[i].element);
        for j in (i + 1)..atoms.len() {
            let rj = covalent_radius(&atoms[j].element);
            let d = distance(atoms[i].pos, atoms[j].pos);
            let max = ri + rj + TOLERANCE;
            if d > MIN_DIST && d <= max {
                bonds.push((i, j));
            }
        }
    }
    bonds
}

/// Euclidean distance between two points.
fn distance(a: [f32; 3], b: [f32; 3]) -> f32 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// How finely spheres / cylinders are tessellated. `8` longitude
/// segments keeps an all-atom structure's triangle count manageable
/// while still reading as round at viewport scale.
const SEGMENTS: usize = 8;

/// Build a **ball-and-stick** triangle mesh of a molecule.
///
/// Each atom is a sphere of radius `vdw_radius(element) * ball_scale`
/// — `ball_scale` around `0.25` gives the classic ball-and-stick look
/// (small balls, visible sticks). Each bond is a cylinder of radius
/// `stick_radius` split at its midpoint so each half takes its atom's
/// colour. An empty molecule yields an empty mesh.
pub fn ball_and_stick(mol: &ViewMolecule, ball_scale: f32, stick_radius: f32) -> TriangleMesh {
    let mut tris: Vec<StlTriangle> = Vec::new();
    for atom in &mol.atoms {
        let r = (vdw_radius(&atom.element) * ball_scale).max(0.05);
        push_sphere(&mut tris, atom.pos, r);
    }
    for &(a, b) in &mol.bonds {
        let (Some(atom_a), Some(atom_b)) = (mol.atoms.get(a), mol.atoms.get(b)) else {
            continue;
        };
        // Split each bond at its midpoint — half-bonds are kept as
        // distinct cylinders so a future per-atom-colour renderer can
        // tint each half by its end atom's element.
        let mid = midpoint(atom_a.pos, atom_b.pos);
        push_cylinder(&mut tris, atom_a.pos, mid, stick_radius);
        push_cylinder(&mut tris, mid, atom_b.pos, stick_radius);
    }
    TriangleMesh {
        format: None,
        name: Some("genetics-ball-and-stick".to_string()),
        triangles: tris,
    }
}

/// Build a **ball-and-stick** triangle mesh of a molecule *with a paired
/// per-triangle CPK colour* — the colour-aware sibling of [`ball_and_stick`].
///
/// Returns the same geometry as [`ball_and_stick`] (identical triangle order)
/// plus a `Vec<[f32; 3]>` carrying **one colour per triangle**, in lockstep with
/// `mesh.triangles`: every triangle of an atom's sphere takes that atom's
/// [`element_color`], and each half of a midpoint-split bond cylinder takes its
/// own end atom's colour. The colour for triangle `k` is `colors[k]`, so
/// `colors.len() == mesh.triangles.len()`.
///
/// This is what lets a colour-aware consumer (the Workbench+Agent product tile,
/// which *does* render per-vertex colour) tint the molecule by element: promote
/// the mesh to a `valenx_mesh::Mesh` one-`Tri3`-per-triangle in order
/// (`products_registry::mesh_from_triangle_soup`) and expand these per-triangle
/// colours to triangle-major per-vertex colours
/// (`products_registry::per_triangle_to_vertex_colors`). The colours are
/// recovered by snapshotting the triangle count before/after each `push_sphere`
/// / `push_cylinder`, so the mapping stays correct regardless of the
/// sphere/cylinder tessellation density.
pub fn ball_and_stick_colored(
    mol: &ViewMolecule,
    ball_scale: f32,
    stick_radius: f32,
) -> (TriangleMesh, Vec<[f32; 3]>) {
    let mut tris: Vec<StlTriangle> = Vec::new();
    let mut colors: Vec<[f32; 3]> = Vec::new();
    // Push `tris` for one primitive, then tag every triangle it added with
    // `color` — independent of how many triangles the tessellator emitted.
    let mut tag = |tris: &mut Vec<StlTriangle>, before: usize, color: [f32; 3]| {
        for _ in before..tris.len() {
            colors.push(color);
        }
    };
    for atom in &mol.atoms {
        let r = (vdw_radius(&atom.element) * ball_scale).max(0.05);
        let before = tris.len();
        push_sphere(&mut tris, atom.pos, r);
        tag(&mut tris, before, element_color(&atom.element));
    }
    for &(a, b) in &mol.bonds {
        let (Some(atom_a), Some(atom_b)) = (mol.atoms.get(a), mol.atoms.get(b)) else {
            continue;
        };
        let mid = midpoint(atom_a.pos, atom_b.pos);
        // Each half-cylinder takes its near atom's colour.
        let before = tris.len();
        push_cylinder(&mut tris, atom_a.pos, mid, stick_radius);
        tag(&mut tris, before, element_color(&atom_a.element));
        let before = tris.len();
        push_cylinder(&mut tris, mid, atom_b.pos, stick_radius);
        tag(&mut tris, before, element_color(&atom_b.element));
    }
    let mesh = TriangleMesh {
        format: None,
        name: Some("genetics-ball-and-stick".to_string()),
        triangles: tris,
    };
    (mesh, colors)
}

/// Build a **spacefill** (CPK / van-der-Waals) triangle mesh.
///
/// Every atom is a full van-der-Waals sphere; bonds are not drawn (the
/// spheres overlap to show connectivity). An empty molecule yields an
/// empty mesh.
pub fn spacefill(mol: &ViewMolecule) -> TriangleMesh {
    let mut tris: Vec<StlTriangle> = Vec::new();
    for atom in &mol.atoms {
        let r = vdw_radius(&atom.element).max(0.05);
        push_sphere(&mut tris, atom.pos, r);
    }
    TriangleMesh {
        format: None,
        name: Some("genetics-spacefill".to_string()),
        triangles: tris,
    }
}

/// Midpoint of two points.
fn midpoint(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        0.5 * (a[0] + b[0]),
        0.5 * (a[1] + b[1]),
        0.5 * (a[2] + b[2]),
    ]
}

/// Append a UV sphere centred at `center` with radius `r` to `tris`.
///
/// A latitude/longitude tessellation with `SEGMENTS` longitude bands
/// and `SEGMENTS` latitude rings. Pole caps are triangles; the body is
/// quads split into two triangles. Each triangle's normal is computed
/// from its winding so STL exports of the mesh are well-formed.
fn push_sphere(tris: &mut Vec<StlTriangle>, center: [f32; 3], r: f32) {
    let lat_steps = SEGMENTS;
    let lon_steps = SEGMENTS * 2;
    // Pre-compute the (lat, lon) vertex grid.
    let vert = |lat: usize, lon: usize| -> [f32; 3] {
        let theta = std::f32::consts::PI * lat as f32 / lat_steps as f32;
        let phi = 2.0 * std::f32::consts::PI * lon as f32 / lon_steps as f32;
        [
            center[0] + r * theta.sin() * phi.cos(),
            center[1] + r * theta.cos(),
            center[2] + r * theta.sin() * phi.sin(),
        ]
    };
    for lat in 0..lat_steps {
        for lon in 0..lon_steps {
            let v00 = vert(lat, lon);
            let v01 = vert(lat, lon + 1);
            let v10 = vert(lat + 1, lon);
            let v11 = vert(lat + 1, lon + 1);
            if lat == 0 {
                // Top cap — one triangle per longitude.
                push_tri(tris, v00, v10, v11);
            } else if lat == lat_steps - 1 {
                // Bottom cap.
                push_tri(tris, v00, v10, v01);
            } else {
                // Body quad → two triangles.
                push_tri(tris, v00, v10, v11);
                push_tri(tris, v00, v11, v01);
            }
        }
    }
}

/// Append a cylinder from `a` to `b` of radius `r` to `tris`.
///
/// The cylinder is a `SEGMENTS`-sided prism — a side wall of quads
/// (two triangles each); the flat end caps are omitted because in a
/// ball-and-stick model the atom spheres always cover them.
fn push_cylinder(tris: &mut Vec<StlTriangle>, a: [f32; 3], b: [f32; 3], r: f32) {
    let axis = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    if len < 1e-6 {
        return;
    }
    let dir = [axis[0] / len, axis[1] / len, axis[2] / len];
    // Two unit vectors spanning the plane perpendicular to `dir`.
    let (u, v) = perpendicular_basis(dir);
    let ring = |center: [f32; 3], seg: usize| -> [f32; 3] {
        let ang = 2.0 * std::f32::consts::PI * seg as f32 / SEGMENTS as f32;
        let (c, s) = (ang.cos(), ang.sin());
        [
            center[0] + r * (c * u[0] + s * v[0]),
            center[1] + r * (c * u[1] + s * v[1]),
            center[2] + r * (c * u[2] + s * v[2]),
        ]
    };
    for seg in 0..SEGMENTS {
        let a0 = ring(a, seg);
        let a1 = ring(a, seg + 1);
        let b0 = ring(b, seg);
        let b1 = ring(b, seg + 1);
        push_tri(tris, a0, b0, b1);
        push_tri(tris, a0, b1, a1);
    }
}

/// An orthonormal pair perpendicular to `dir` (assumed unit length).
fn perpendicular_basis(dir: [f32; 3]) -> ([f32; 3], [f32; 3]) {
    // Pick the world axis least aligned with `dir` as a seed so the
    // cross product is well-conditioned.
    let seed = if dir[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = normalize(cross(dir, seed));
    let v = normalize(cross(dir, u));
    (u, v)
}

/// Cross product.
fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Normalise a vector; returns `[0, 0, 1]` for a (near-)zero input.
fn normalize(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-9 {
        [0.0, 0.0, 1.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

/// Append one triangle (with a winding-derived normal) to `tris`.
fn push_tri(tris: &mut Vec<StlTriangle>, v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) {
    let tri = StlTriangle {
        normal: [0.0, 0.0, 0.0],
        vertices: [v0, v1, v2],
    };
    let normal = tri.computed_normal();
    tris.push(StlTriangle {
        normal,
        vertices: [v0, v1, v2],
    });
}

/// Push a built molecular mesh into the app's wgpu 3-D viewport.
///
/// Sets `ValenxApp::stl` to the triangle soup — the same field the
/// STL importer fills — and clears any canonical `ValenxApp::mesh`
/// so the molecule is the geometry on screen (the viewport overlays a
/// loaded mesh on top of the STL otherwise). The camera is framed onto
/// the new geometry. `source_label` becomes the viewport header text.
///
/// Returns an error string if the mesh is empty (nothing to show).
pub fn show_molecule(
    app: &mut ValenxApp,
    mesh: TriangleMesh,
    source_label: &str,
) -> Result<usize, String> {
    if mesh.triangles.is_empty() {
        return Err("nothing to display — the molecule has no atoms".to_string());
    }
    let count = mesh.triangle_count();
    app.mesh = None;
    app.stl = Some(LoadedStl::new(
        PathBuf::from(format!("<genetics>/{source_label}")),
        mesh,
    ));
    app.frame_current_stl();
    Ok(count)
}

/// Push a built molecular mesh into the viewport **carrying per-vertex colours**
/// so the shaded central viewport tints each atom instead of rendering one
/// material. Same viewport plumbing as [`show_molecule`] — sets `ValenxApp::stl`
/// and clears `ValenxApp::mesh`, then frames the camera — but attaches a
/// per-vertex colour buffer (`LoadedStl::colors`) that the viewport's shaded
/// wgpu path uploads via
/// [`crate::wgpu_renderer::triangles_to_vertices_colored`].
///
/// `per_tri_colors` is **one colour per triangle**, in lockstep with
/// `mesh.triangles` (the shape [`crate::molviz::build_mesh_colored`] returns);
/// it is expanded here to the triangle-major per-vertex layout the renderer
/// emits (three copies per triangle, via
/// `crate::products_registry::per_triangle_to_vertex_colors`) so it lines up
/// 1:1 with the surface vertices — exactly the per-triangle→per-vertex bridge
/// the Workbench+Agent molecule tile (`molecule_product`) uses. A length
/// mismatch can never half-colour the mesh: the viewport length-guards the
/// buffer and falls back to neutral metal.
///
/// Returns an error string (and pushes nothing) if the mesh is empty.
pub fn show_molecule_colored(
    app: &mut ValenxApp,
    mesh: TriangleMesh,
    per_tri_colors: &[[f32; 3]],
    source_label: &str,
) -> Result<usize, String> {
    if mesh.triangles.is_empty() {
        return Err("nothing to display — the molecule has no atoms".to_string());
    }
    let count = mesh.triangle_count();
    let vertex_colors = crate::products_registry::per_triangle_to_vertex_colors(per_tri_colors);
    app.mesh = None;
    app.stl = Some(LoadedStl::with_colors(
        PathBuf::from(format!("<genetics>/{source_label}")),
        mesh,
        vertex_colors,
    ));
    app.frame_current_stl();
    Ok(count)
}

/// A small canonical demo molecule — a single water (H₂O) with realistic
/// geometry (O at the origin, two O–H bonds at ~104.5°, 0.96 Å), bonds detected
/// by the covalent-radius rule. Used by `molecule_product` so the
/// agent-bridge molecule tile renders a real, correctly-coloured structure with
/// no external data.
fn demo_molecule() -> ViewMolecule {
    let mut mol = ViewMolecule {
        atoms: vec![
            ViewAtom::new([0.000, 0.000, 0.1173], "O"),
            ViewAtom::new([0.000, 0.7572, -0.4692], "H"),
            ViewAtom::new([0.000, -0.7572, -0.4692], "H"),
        ],
        bonds: Vec::new(),
    };
    mol.bonds = detect_bonds(&mol.atoms);
    mol
}

/// Count atoms by element and emit a Hill-ish formula string (C, then H, then
/// the rest alphabetically) for a [`ViewMolecule`] — a compact readout row.
fn molecule_formula(mol: &ViewMolecule) -> String {
    use std::collections::BTreeMap;
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for a in &mol.atoms {
        *counts
            .entry(a.element.trim().to_ascii_uppercase())
            .or_default() += 1;
    }
    let mut order: Vec<String> = Vec::new();
    if counts.contains_key("C") {
        order.push("C".to_string());
    }
    if counts.contains_key("H") {
        order.push("H".to_string());
    }
    for k in counts.keys() {
        if k != "C" && k != "H" {
            order.push(k.clone());
        }
    }
    let mut s = String::new();
    for el in order {
        let n = counts[&el];
        // Title-case the symbol (e.g. "FE" → "Fe") for display.
        let mut sym = String::new();
        for (i, ch) in el.chars().enumerate() {
            if i == 0 {
                sym.extend(ch.to_uppercase());
            } else {
                sym.extend(ch.to_lowercase());
            }
        }
        s.push_str(&sym);
        if n > 1 {
            s.push_str(&n.to_string());
        }
    }
    s
}

/// The agent-bridge product for the molecule view (`show_3d{kind="molecule"}`).
///
/// Builds the [`demo_molecule`] (a water H₂O), meshes it as a **colour-aware**
/// ball-and-stick model via [`ball_and_stick_colored`], promotes the triangle
/// soup to a `Tri3` [`valenx_mesh::Mesh`]
/// ([`crate::products_registry::mesh_from_triangle_soup`]), and expands the
/// per-triangle CPK colours to the triangle-major per-vertex `vertex_colors` the
/// tile renderer paints
/// (`crate::products_registry::per_triangle_to_vertex_colors`) — so the
/// molecule renders coloured by element (O red, H white) rather than flat metal.
/// Pure and app-state-free. The readout reports the formula and atom/bond
/// counts.
pub(crate) fn molecule_product() -> crate::WorkspaceProduct {
    let mol = demo_molecule();
    let (soup, per_tri_colors) = ball_and_stick_colored(&mol, 0.25, 0.15);
    let mesh = crate::products_registry::mesh_from_triangle_soup(&soup, "valenx-molecule");
    let vertex_colors = crate::products_registry::per_triangle_to_vertex_colors(&per_tri_colors);
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<molecule>/ball-and-stick");
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    let lines = vec![
        format!("molecule: {} (water)", molecule_formula(&mol)),
        format!("{} atoms · {} bonds", mol.atoms.len(), mol.bonds.len()),
        "ball-and-stick, coloured by CPK element palette".to_string(),
    ];
    crate::WorkspaceProduct {
        title: "Molecule (ball-and-stick)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: Some(vertex_colors),
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

    /// A water-like 3-atom cluster: O at the origin, two H at ~0.96 Å.
    fn water() -> ViewMolecule {
        ViewMolecule {
            atoms: vec![
                ViewAtom::new([0.0, 0.0, 0.0], "O"),
                ViewAtom::new([0.96, 0.0, 0.0], "H"),
                ViewAtom::new([-0.24, 0.93, 0.0], "H"),
            ],
            bonds: vec![],
        }
    }

    #[test]
    fn detect_bonds_finds_the_two_oh_bonds() {
        let atoms = water().atoms;
        let bonds = detect_bonds(&atoms);
        // O-H1 and O-H2 are within range; H-H is not.
        assert_eq!(bonds.len(), 2);
        assert!(bonds.contains(&(0, 1)));
        assert!(bonds.contains(&(0, 2)));
        assert!(!bonds.contains(&(1, 2)));
    }

    #[test]
    fn detect_bonds_skips_far_apart_atoms() {
        let atoms = vec![
            ViewAtom::new([0.0, 0.0, 0.0], "C"),
            ViewAtom::new([10.0, 0.0, 0.0], "C"),
        ];
        assert!(detect_bonds(&atoms).is_empty());
    }

    #[test]
    fn detect_bonds_rejects_coincident_atoms() {
        // Two atoms on top of each other must not bond (MIN_DIST floor).
        let atoms = vec![
            ViewAtom::new([0.0, 0.0, 0.0], "C"),
            ViewAtom::new([0.01, 0.0, 0.0], "C"),
        ];
        assert!(detect_bonds(&atoms).is_empty());
    }

    #[test]
    fn ball_and_stick_meshes_atoms_and_bonds() {
        let mut mol = water();
        mol.bonds = detect_bonds(&mol.atoms);
        let mesh = ball_and_stick(&mol, 0.25, 0.15);
        // 3 spheres + 2 bonds (each split in two) → non-trivial soup.
        assert!(mesh.triangle_count() > 100);
        // Every triangle must have three finite vertices.
        for t in &mesh.triangles {
            for v in &t.vertices {
                assert!(v.iter().all(|c| c.is_finite()));
            }
        }
    }

    #[test]
    fn spacefill_has_no_bonds_only_spheres() {
        let mol = water();
        let spheres_only = spacefill(&mol);
        // One sphere per atom. The same molecule with bonds drawn must
        // have strictly more triangles.
        let mut bonded = mol.clone();
        bonded.bonds = detect_bonds(&bonded.atoms);
        let with_sticks = ball_and_stick(&bonded, 0.25, 0.15);
        assert!(with_sticks.triangle_count() > spheres_only.triangle_count());
    }

    #[test]
    fn empty_molecule_yields_empty_mesh() {
        let empty = ViewMolecule::new();
        assert!(ball_and_stick(&empty, 0.25, 0.15).triangles.is_empty());
        assert!(spacefill(&empty).triangles.is_empty());
    }

    #[test]
    fn ball_and_stick_bounding_box_encloses_all_atoms() {
        let mut mol = water();
        mol.bonds = detect_bonds(&mol.atoms);
        let mesh = ball_and_stick(&mol, 0.25, 0.15);
        let (min, max) = mesh.bounding_box().expect("non-empty");
        // Every atom centre must lie inside the mesh AABB.
        for atom in &mol.atoms {
            for k in 0..3 {
                assert!(atom.pos[k] >= min[k] - 1e-3);
                assert!(atom.pos[k] <= max[k] + 1e-3);
            }
        }
    }

    #[test]
    fn show_molecule_rejects_empty_mesh() {
        let mut app = ValenxApp::default();
        let empty = TriangleMesh::new();
        assert!(super::show_molecule(&mut app, empty, "test").is_err());
    }

    #[test]
    fn show_molecule_sets_the_viewport_stl() {
        let mut app = ValenxApp::default();
        let mut mol = water();
        mol.bonds = detect_bonds(&mol.atoms);
        let mesh = ball_and_stick(&mol, 0.25, 0.15);
        let n = super::show_molecule(&mut app, mesh, "water.bs").expect("non-empty mesh");
        assert!(n > 0);
        assert!(app.stl.is_some());
        // The canonical mesh slot is cleared so the molecule is what
        // the viewport draws.
        assert!(app.mesh.is_none());
        // The plain (uncoloured) path attaches no per-vertex colour buffer.
        assert!(app.stl.as_ref().unwrap().colors.is_none());
    }

    #[test]
    fn show_molecule_colored_attaches_per_vertex_colors() {
        // The colour-aware path attaches a per-vertex colour buffer of exactly
        // 3 × triangle count (one per surface vertex), expanded from the
        // per-triangle colours the colour-aware builder returns.
        let mut app = ValenxApp::default();
        let mut mol = water();
        mol.bonds = detect_bonds(&mol.atoms);
        let (mesh, per_tri) = ball_and_stick_colored(&mol, 0.25, 0.15);
        let tri_count = mesh.triangle_count();
        assert_eq!(per_tri.len(), tri_count, "one colour per triangle");
        let n =
            super::show_molecule_colored(&mut app, mesh, &per_tri, "water.bs").expect("non-empty");
        assert!(n > 0);
        let stl = app.stl.as_ref().expect("viewport STL set");
        let colors = stl.colors.as_ref().expect("per-vertex colours attached");
        assert_eq!(
            colors.len(),
            tri_count * 3,
            "per-vertex colours = 3 × triangle count (triangles_to_vertices order)"
        );
        assert!(app.mesh.is_none());
    }

    #[test]
    fn show_molecule_colored_rejects_empty_mesh() {
        let mut app = ValenxApp::default();
        let empty = TriangleMesh::new();
        assert!(super::show_molecule_colored(&mut app, empty, &[], "test").is_err());
    }

    #[test]
    fn covalent_and_vdw_radii_have_sane_values() {
        // Hydrogen is the smallest; carbon mid-range; the fallback is
        // carbon-ish for an unknown element.
        assert!(covalent_radius("H") < covalent_radius("C"));
        assert!(vdw_radius("H") < vdw_radius("C"));
        assert!(vdw_radius("Xx") > 1.0);
        assert!(covalent_radius("") > 0.0);
    }

    #[test]
    fn element_color_distinguishes_common_elements() {
        // Carbon, nitrogen and oxygen must not collide.
        let c = element_color("C");
        let n = element_color("N");
        let o = element_color("O");
        assert_ne!(c, n);
        assert_ne!(n, o);
        assert_ne!(c, o);
    }

    #[test]
    fn from_md_system_converts_nm_to_angstrom() {
        use nalgebra::Vector3;
        use valenx_md::system::{Atom, System, Topology};

        let mut top = Topology::new();
        top.push_atom(Atom::new("C", 12.011, 0.0).unwrap().with_element("C"));
        top.push_atom(Atom::new("O", 15.999, 0.0).unwrap().with_element("O"));
        // 0.12 nm apart — a C=O-ish bond.
        let pos = vec![Vector3::zeros(), Vector3::new(0.12, 0.0, 0.0)];
        let system = System::new(top, pos).unwrap();
        let view = ViewMolecule::from_md_system(&system);
        assert_eq!(view.atoms.len(), 2);
        // 0.12 nm → 1.2 Å.
        assert!((view.atoms[1].pos[0] - 1.2).abs() < 1e-4);
        // Distance-rule bond detection should join the C and O.
        assert_eq!(view.bonds.len(), 1);
    }
}
