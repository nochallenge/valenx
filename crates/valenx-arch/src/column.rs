//! Column entity — a vertical prismatic column extruded from a
//! cross-section profile.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::error::ArchError;

/// Round-4 DoS hardening: upper bound on the polygon segment count
/// for circular columns. 4096 leaves plenty of headroom for the most
/// ornamental column anyone would realistically draw while keeping
/// the per-column allocation bounded.
pub const MAX_COLUMN_SEGMENTS: u32 = 4096;

/// Standard cross-sections a column can carry. Real structural design
/// requires many more (W-shapes, channels, angles, hollow sections);
/// v1 covers the three a residential/light-commercial floor plan
/// usually needs.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum ColumnSection {
    /// Rectangular column — width × depth.
    Rectangle {
        /// Width along +X (in world units).
        width: f64,
        /// Depth along +Y (in world units).
        depth: f64,
    },
    /// Circular column.
    Circular {
        /// Outer radius (in world units).
        radius: f64,
        /// Polygon segment count for tessellation. 12 is a sensible
        /// v1 default — high enough that the column looks round at
        /// normal viewport zooms.
        segments: u32,
    },
    /// I-beam / wide-flange column. Doubly-symmetric profile with two
    /// horizontal flanges connected by a vertical web. Outer
    /// bounding box is `width × depth`. v1 emits the I as the outer
    /// bounding box for tessellation simplicity — Phase 15.5 will
    /// emit the true I-shape via a polygon extrude. The label still
    /// reports the proper kind for the schedule.
    IBeam {
        /// Flange width (along +X).
        width: f64,
        /// Profile depth / overall height (along +Y).
        depth: f64,
        /// Flange thickness — stored for IFC + schedule.
        flange_thickness: f64,
        /// Web thickness — stored for IFC + schedule.
        web_thickness: f64,
    },
}

impl ColumnSection {
    /// Short human label for UI lists / schedules.
    pub fn label(&self) -> &'static str {
        match self {
            ColumnSection::Rectangle { .. } => "Rect",
            ColumnSection::Circular { .. } => "Circ",
            ColumnSection::IBeam { .. } => "IBeam",
        }
    }

    /// Cross-section area (m² assuming inputs are in metres).
    pub fn area(&self) -> f64 {
        match self {
            ColumnSection::Rectangle { width, depth } => width * depth,
            ColumnSection::Circular { radius, .. } => std::f64::consts::PI * radius * radius,
            ColumnSection::IBeam {
                width,
                depth,
                flange_thickness,
                web_thickness,
            } => {
                // Flange contribution + (web depth) × web thickness.
                let inner_depth = (depth - 2.0 * flange_thickness).max(0.0);
                2.0 * width * flange_thickness + inner_depth * web_thickness
            }
        }
    }
}

/// Parameters describing a column.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ColumnParams {
    /// Base point (centre of the section), in world space. The
    /// column extrudes from `base` along +Z by [`Self::height`].
    pub base: Vector3<f64>,
    /// Column height (along +Z).
    pub height: f64,
    /// Cross-section profile.
    pub cross_section: ColumnSection,
    /// Material descriptor.
    pub material: String,
    /// Optional structural attributes — material grade + axial load.
    /// `None` for non-structural columns. Consumed by
    /// [`crate::structural::export_structural_model`].
    #[serde(default)]
    pub structural: Option<crate::structural::StructuralMember>,
}

impl ColumnParams {
    /// Validate dimensions.
    pub fn validate(&self) -> Result<(), ArchError> {
        if !self.height.is_finite() || self.height <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "height",
                reason: format!("must be > 0 (got {})", self.height),
            });
        }
        match &self.cross_section {
            ColumnSection::Rectangle { width, depth } => {
                if *width <= 0.0 || *depth <= 0.0 {
                    return Err(ArchError::BadDimension {
                        name: "rect",
                        reason: format!("width/depth must be > 0 (got {width}, {depth})"),
                    });
                }
            }
            ColumnSection::Circular { radius, segments } => {
                if *radius <= 0.0 {
                    return Err(ArchError::BadDimension {
                        name: "radius",
                        reason: format!("must be > 0 (got {radius})"),
                    });
                }
                if *segments < 3 {
                    return Err(ArchError::BadDimension {
                        name: "segments",
                        reason: format!("must be ≥ 3 (got {segments})"),
                    });
                }
                // Round-4 DoS hardening: a malicious BCF can advertise
                // `segments = u32::MAX` which then drives
                // `Vec::with_capacity(*segments as usize)` to OOM the
                // host. Cap at 4096 — production columns are smooth at
                // ~32 segments, complex ornamental columns at ~256, and
                // 4096 leaves plenty of headroom for the worst legitimate
                // case while keeping the allocation under 200 KB.
                if *segments > MAX_COLUMN_SEGMENTS {
                    return Err(ArchError::BadDimension {
                        name: "segments",
                        reason: format!(
                            "must be ≤ {MAX_COLUMN_SEGMENTS} (got {segments}) — \
                             DoS guard against pathological circular columns"
                        ),
                    });
                }
            }
            ColumnSection::IBeam {
                width,
                depth,
                flange_thickness,
                web_thickness,
            } => {
                if *width <= 0.0
                    || *depth <= 0.0
                    || *flange_thickness <= 0.0
                    || *web_thickness <= 0.0
                {
                    return Err(ArchError::BadDimension {
                        name: "ibeam",
                        reason: format!(
                            "dims must be > 0 (got w={width}, d={depth}, ft={flange_thickness}, wt={web_thickness})"
                        ),
                    });
                }
            }
        }
        Ok(())
    }

    /// Generate the cross-section's outline polygon in the X-Y plane
    /// centred on the origin. The column extrudes this polygon along
    /// +Z then translates by [`Self::base`].
    fn profile_xy(&self) -> Vec<(f64, f64)> {
        match &self.cross_section {
            ColumnSection::Rectangle { width, depth } => {
                let hw = width * 0.5;
                let hd = depth * 0.5;
                vec![(-hw, -hd), (hw, -hd), (hw, hd), (-hw, hd)]
            }
            ColumnSection::Circular { radius, segments } => {
                let mut v = Vec::with_capacity(*segments as usize);
                for k in 0..*segments {
                    let a = (k as f64) * 2.0 * std::f64::consts::PI / (*segments as f64);
                    v.push((radius * a.cos(), radius * a.sin()));
                }
                v
            }
            // I-Beam: emit the outer bounding rectangle for tessellation.
            // v1 simplification documented above.
            ColumnSection::IBeam { width, depth, .. } => {
                let hw = width * 0.5;
                let hd = depth * 0.5;
                vec![(-hw, -hd), (hw, -hd), (hw, hd), (-hw, hd)]
            }
        }
    }

    /// Tessellate the column into a [`valenx_mesh::Mesh`].
    ///
    /// Same fan-triangulation strategy as [`crate::SlabParams`] but
    /// the profile is convex and centred on the origin, then
    /// translated.
    pub fn tessellate_mesh(&self) -> Result<Mesh, ArchError> {
        self.validate()?;
        let profile = self.profile_xy();
        let n = profile.len();
        let z0 = self.base.z;
        let z1 = z0 + self.height;
        let (cx, cy) = (self.base.x, self.base.y);

        let mut mesh = Mesh::new("column");
        for (x, y) in &profile {
            mesh.nodes.push(Vector3::new(cx + x, cy + y, z0));
        }
        for (x, y) in &profile {
            mesh.nodes.push(Vector3::new(cx + x, cy + y, z1));
        }
        let mut block = ElementBlock::new(ElementType::Tri3);
        for i in 1..(n - 1) {
            block
                .connectivity
                .extend_from_slice(&[0_u32, (i + 1) as u32, i as u32]);
        }
        for i in 1..(n - 1) {
            block
                .connectivity
                .extend_from_slice(&[n as u32, (n + i) as u32, (n + i + 1) as u32]);
        }
        for i in 0..n {
            let j = (i + 1) % n;
            let b0 = i as u32;
            let b1 = j as u32;
            let t0 = (n + i) as u32;
            let t1 = (n + j) as u32;
            block.connectivity.extend_from_slice(&[b0, b1, t1]);
            block.connectivity.extend_from_slice(&[b0, t1, t0]);
        }
        mesh.element_blocks.push(block);
        mesh.recompute_stats();
        Ok(mesh)
    }

    /// Wrap [`Self::tessellate_mesh`] in a mesh-backed
    /// [`valenx_cad::Solid`].
    pub fn tessellate(&self) -> Result<valenx_cad::Solid, ArchError> {
        let m = self.tessellate_mesh()?;
        Ok(valenx_cad::Solid::from_mesh(m))
    }

    /// Hint points for [`crate::ArchDocument::bbox`].
    pub fn bbox_hint_points(&self) -> impl Iterator<Item = Vector3<f64>> + '_ {
        let profile = self.profile_xy();
        let mut v = Vec::with_capacity(profile.len() * 2);
        let z0 = self.base.z;
        let z1 = z0 + self.height;
        for (x, y) in profile {
            v.push(Vector3::new(self.base.x + x, self.base.y + y, z0));
            v.push(Vector3::new(self.base.x + x, self.base.y + y, z1));
        }
        v.into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rect_validation() {
        let c = ColumnParams {
            base: Vector3::zeros(),
            height: 3.0,
            cross_section: ColumnSection::Rectangle {
                width: 0.3,
                depth: 0.3,
            },
            material: "Steel".into(),
            structural: None,
        };
        c.validate().unwrap();
        assert!(c.cross_section.area() > 0.0);
        assert_eq!(c.cross_section.label(), "Rect");
    }

    #[test]
    fn circular_validation_segments() {
        let mut c = ColumnParams {
            base: Vector3::zeros(),
            height: 3.0,
            cross_section: ColumnSection::Circular {
                radius: 0.2,
                segments: 12,
            },
            material: "Concrete".into(),
            structural: None,
        };
        c.validate().unwrap();
        if let ColumnSection::Circular { segments, .. } = &mut c.cross_section {
            *segments = 2;
        }
        assert!(matches!(
            c.validate(),
            Err(ArchError::BadDimension {
                name: "segments",
                ..
            })
        ));
    }

    #[test]
    fn rejects_negative_height() {
        let c = ColumnParams {
            base: Vector3::zeros(),
            height: -1.0,
            cross_section: ColumnSection::Rectangle {
                width: 0.1,
                depth: 0.1,
            },
            material: "x".into(),
            structural: None,
        };
        assert!(matches!(
            c.validate(),
            Err(ArchError::BadDimension { name: "height", .. })
        ));
    }

    #[test]
    fn ibeam_area_uses_flange_and_web() {
        let s = ColumnSection::IBeam {
            width: 0.2,
            depth: 0.3,
            flange_thickness: 0.02,
            web_thickness: 0.01,
        };
        // A = 2 * 0.2 * 0.02 + (0.3 - 0.04) * 0.01 = 0.008 + 0.0026 = 0.0106
        assert!((s.area() - 0.0106).abs() < 1e-9);
    }

    #[test]
    fn rect_column_tessellation_topology() {
        let c = ColumnParams {
            base: Vector3::new(1.0, 2.0, 3.0),
            height: 3.0,
            cross_section: ColumnSection::Rectangle {
                width: 0.4,
                depth: 0.4,
            },
            material: "Steel".into(),
            structural: None,
        };
        let m = c.tessellate_mesh().unwrap();
        // 4-sided prism: 8 nodes, 12 tris (2 top + 2 bot + 4 sides * 2).
        assert_eq!(m.nodes.len(), 8);
        assert_eq!(m.total_elements(), 12);
    }
}
