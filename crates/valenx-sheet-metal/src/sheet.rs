//! Sheet — outline polygon + thickness + material + bend recipe list.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::bend::{Bend, Flange};
use crate::error::SheetMetalError;

/// Sheet material — currently a textual tag. The structural property
/// table lives downstream (e.g. valenx-fem's FemMaterial); this
/// crate's surface stays declarative.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SheetMaterial {
    /// Display name (e.g. "Aluminium 6061-T6", "Steel AISI 1045").
    pub name: String,
}

impl SheetMaterial {
    /// Construct a material.
    pub fn new(name: impl Into<String>) -> Self {
        Self { name: name.into() }
    }
}

/// A parametric sheet metal plate.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Sheet {
    /// Outline polygon (closed loop, `[u, v]` coords — last vertex
    /// is NOT a duplicate of the first).
    pub outline: Vec<[f64; 2]>,
    /// Plate thickness (sheet units).
    pub thickness: f64,
    /// Material tag.
    pub material: SheetMaterial,
    /// K-factor (neutral-axis fraction, typical 0.33–0.5 for soft
    /// aluminium / mild steel; 0.44 is a sane default).
    pub k_factor: f64,
    /// Bend recipe (applied in order by [`Sheet::to_solid`]).
    pub bends: Vec<Bend>,
    /// Flange recipe.
    pub flanges: Vec<Flange>,
    /// Cutout list (each is a 2D polygon subtracted from the
    /// outline).
    pub cutouts: Vec<Vec<[f64; 2]>>,
}

impl Default for Sheet {
    fn default() -> Self {
        // Default 100 x 100 plate, 1.5 mm thickness, aluminium k=0.44.
        Self {
            outline: vec![[0.0, 0.0], [100.0, 0.0], [100.0, 100.0], [0.0, 100.0]],
            thickness: 1.5,
            material: SheetMaterial::new("Aluminium 6061-T6"),
            k_factor: 0.44,
            bends: Vec::new(),
            flanges: Vec::new(),
            cutouts: Vec::new(),
        }
    }
}

impl Sheet {
    /// Construct with a rectangular outline `width × height`.
    pub fn rectangle(width: f64, height: f64, thickness: f64) -> Result<Self, SheetMetalError> {
        if width <= 0.0 || height <= 0.0 || thickness <= 0.0 {
            return Err(SheetMetalError::BadParameter {
                name: "rectangle dims",
                reason: format!("must be > 0, got {width}x{height} t={thickness}"),
            });
        }
        Ok(Self {
            outline: vec![[0.0, 0.0], [width, 0.0], [width, height], [0.0, height]],
            thickness,
            ..Default::default()
        })
    }

    /// Number of outline edges (one per polygon segment, with the
    /// loop closing back to vertex 0).
    pub fn edge_count(&self) -> usize {
        self.outline.len()
    }

    /// Append a bend.
    pub fn add_bend(self, bend: Bend) -> Self {
        let mut s = self;
        s.bends.push(bend);
        s
    }

    /// Append a flange.
    pub fn add_flange(self, flange: Flange) -> Result<Self, SheetMetalError> {
        if flange.edge_id >= self.edge_count() {
            return Err(SheetMetalError::BadEdge {
                got: flange.edge_id,
                n: self.edge_count(),
            });
        }
        let mut s = self;
        s.flanges.push(flange);
        Ok(s)
    }

    /// Append a cutout polygon (caller's responsibility to ensure
    /// the polygon lies inside the outline — v1 doesn't enforce).
    pub fn cutout(self, polygon: Vec<[f64; 2]>) -> Result<Self, SheetMetalError> {
        if polygon.len() < 3 {
            return Err(SheetMetalError::BadPolygon(
                "cutout has fewer than 3 vertices".into(),
            ));
        }
        let mut s = self;
        s.cutouts.push(polygon);
        Ok(s)
    }

    /// Tessellate the sheet into a thick plate, then approximate each
    /// bend by rotating the downstream node set around the bend axis.
    /// Returns a [`Solid::Mesh`] suitable for the viewport.
    ///
    /// **v1 limitation:** does not subdivide along the bend line, so
    /// the bend approximation is faceted — Phase 34.5 introduces a
    /// proper bend region. Bends and flanges share the same
    /// approximation kernel for visual purposes; for accurate fabrication
    /// data, call [`Sheet::unfold`].
    pub fn to_solid(&self) -> Result<Solid, SheetMetalError> {
        let mesh = build_thick_plate_mesh(&self.outline, self.thickness)?;
        let mut deformed = mesh;
        for b in &self.bends {
            apply_bend_inplace(&mut deformed, b);
        }
        Ok(Solid::from_mesh(deformed))
    }

    /// Flatten the sheet to a 2D pattern, accounting for k-factor
    /// bend allowance. Currently returns the outline plus per-bend
    /// allowance strips appended to the right edge of the outline
    /// (caller's responsibility to nest into a stock blank).
    pub fn unfold(&self) -> Result<Vec<[f64; 2]>, SheetMetalError> {
        if self.outline.len() < 3 {
            return Err(SheetMetalError::BadPolygon(
                "outline has fewer than 3 vertices".into(),
            ));
        }
        let mut flat = self.outline.clone();
        // Append bend-allowance strips to the right of the outline so
        // the nest stays connected. Each strip is a thin rectangle of
        // width `bend_allowance(thickness, k_factor)` × `bend_line_length`.
        let bbox_w = flat.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max);
        let mut cursor_x = bbox_w + 5.0;
        for b in &self.bends {
            let line_len = {
                let dx = b.end[0] - b.start[0];
                let dy = b.end[1] - b.start[1];
                (dx * dx + dy * dy).sqrt()
            };
            let ba = b.bend_allowance(self.thickness, self.k_factor);
            // Add a rectangle at (cursor_x, 0) of size ba × line_len.
            flat.push([cursor_x, 0.0]);
            flat.push([cursor_x + ba, 0.0]);
            flat.push([cursor_x + ba, line_len]);
            flat.push([cursor_x, line_len]);
            cursor_x += ba + 3.0;
        }
        Ok(flat)
    }
}

// ---------- internal helpers ----------

fn build_thick_plate_mesh(outline: &[[f64; 2]], thickness: f64) -> Result<Mesh, SheetMetalError> {
    if outline.len() < 3 {
        return Err(SheetMetalError::BadPolygon(
            "outline has fewer than 3 vertices".into(),
        ));
    }
    if thickness <= 0.0 {
        return Err(SheetMetalError::BadParameter {
            name: "thickness",
            reason: format!("must be > 0, got {thickness}"),
        });
    }
    let n = outline.len();
    let mut mesh = Mesh::new("sheet_plate");
    // Bottom + top vertices.
    for v in outline {
        mesh.nodes.push(Vector3::new(v[0], v[1], 0.0));
    }
    for v in outline {
        mesh.nodes.push(Vector3::new(v[0], v[1], thickness));
    }
    let mut block = ElementBlock::new(ElementType::Tri3);
    // Top + bottom face fans (triangulate around vertex 0). This
    // works for convex outlines; concave outlines need
    // valenx-meshpart::triangulate_polygon (Phase 34.5 follow-up).
    for i in 1..(n - 1) {
        // Bottom: 0, i, i+1.
        block.connectivity.push(0);
        block.connectivity.push(i as u32);
        block.connectivity.push((i + 1) as u32);
        // Top: n + 0, n + i + 1, n + i (reverse winding).
        let base = n as u32;
        block.connectivity.push(base);
        block.connectivity.push(base + (i + 1) as u32);
        block.connectivity.push(base + i as u32);
    }
    // Side walls: one quad per outline edge.
    for i in 0..n {
        let j = (i + 1) % n;
        let bi = i as u32;
        let bj = j as u32;
        let ti = (n + i) as u32;
        let tj = (n + j) as u32;
        block.connectivity.extend_from_slice(&[bi, bj, tj]);
        block.connectivity.extend_from_slice(&[bi, tj, ti]);
    }
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(mesh)
}

fn apply_bend_inplace(mesh: &mut Mesh, bend: &Bend) {
    // Determine which side of the bend line each node is on (in the
    // XY plane), rotate downstream nodes around the bend axis (the
    // bend line lifted to 3D along z = 0).
    let s = Vector3::new(bend.start[0], bend.start[1], 0.0);
    let e = Vector3::new(bend.end[0], bend.end[1], 0.0);
    let dir = (e - s).try_normalize(1e-12).unwrap_or_else(Vector3::x);
    let (sin_a, cos_a) = bend.angle_rad.sin_cos();
    for node in &mut mesh.nodes {
        // Vector from bend line origin to node.
        let rel = *node - s;
        // Project rel onto bend direction → parallel component.
        let parallel = dir * rel.dot(&dir);
        let perp = rel - parallel; // perpendicular component in 3D
                                   // Determine "downstream" side: cross product with world +Z
                                   // gives the in-plane normal to the bend line.
        let in_plane_normal = dir.cross(&Vector3::z());
        let side = perp.dot(&in_plane_normal);
        if side <= 0.0 {
            continue;
        }
        // Rotate perp around dir by angle_rad (Rodrigues').
        let perp_rot =
            perp * cos_a + dir.cross(&perp) * sin_a + dir * (dir.dot(&perp)) * (1.0 - cos_a);
        *node = s + parallel + perp_rot;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rectangle_validates() {
        let r = Sheet::rectangle(-1.0, 100.0, 1.5);
        assert!(matches!(r, Err(SheetMetalError::BadParameter { .. })));
        let s = Sheet::rectangle(100.0, 50.0, 1.5).unwrap();
        assert_eq!(s.edge_count(), 4);
    }

    #[test]
    fn bend_added() {
        let s = Sheet::rectangle(100.0, 50.0, 1.5).unwrap();
        let s2 = s.add_bend(Bend::new([50.0, 0.0], [50.0, 50.0], 1.0, 1.0));
        assert_eq!(s2.bends.len(), 1);
    }

    #[test]
    fn flange_validates_edge() {
        let s = Sheet::rectangle(100.0, 50.0, 1.5).unwrap();
        let bad = s.clone().add_flange(Flange::new(99, 10.0, 1.0));
        assert!(matches!(bad, Err(SheetMetalError::BadEdge { .. })));
        let ok = s.add_flange(Flange::new(0, 10.0, 1.0)).unwrap();
        assert_eq!(ok.flanges.len(), 1);
    }

    #[test]
    fn cutout_polygon_min_3_verts() {
        let s = Sheet::rectangle(100.0, 50.0, 1.5).unwrap();
        let bad = s.clone().cutout(vec![[10.0, 10.0], [20.0, 10.0]]);
        assert!(matches!(bad, Err(SheetMetalError::BadPolygon(_))));
        let ok = s
            .cutout(vec![[10.0, 10.0], [20.0, 10.0], [15.0, 20.0]])
            .unwrap();
        assert_eq!(ok.cutouts.len(), 1);
    }

    #[test]
    fn to_solid_basic() {
        let s = Sheet::rectangle(100.0, 50.0, 1.5).unwrap();
        let solid = s.to_solid().unwrap();
        // Expect Mesh-backed solid.
        match solid {
            Solid::Mesh(m) => assert!(!m.nodes.is_empty()),
            _ => panic!("expected mesh-backed"),
        }
    }

    #[test]
    fn unfold_returns_outline_plus_bend_strips() {
        let s = Sheet::rectangle(100.0, 50.0, 1.5)
            .unwrap()
            .add_bend(Bend::new([50.0, 0.0], [50.0, 50.0], 1.5, 1.0));
        let flat = s.unfold().unwrap();
        // Outline + 4 strip vertices.
        assert!(flat.len() >= 8);
    }
}
