//! Beam entity — a horizontal (or sloped) prismatic member extruded
//! along an axis from `start` to `end`.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::error::ArchError;

/// Standard cross-sections a beam can carry. Mirrors
/// [`crate::ColumnSection`] but lives in its own enum because
/// extending one without the other later (e.g. adding a T-shape only
/// to beams) is cleaner.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum BeamSection {
    /// Rectangular beam.
    Rectangle {
        /// Width perpendicular to the beam axis.
        width: f64,
        /// Depth perpendicular to width.
        depth: f64,
    },
    /// I-beam — see [`crate::ColumnSection::IBeam`] for the v1
    /// simplification (tessellated as the outer bounding rectangle).
    IBeam {
        /// Flange width.
        width: f64,
        /// Profile depth.
        depth: f64,
        /// Flange thickness.
        flange_thickness: f64,
        /// Web thickness.
        web_thickness: f64,
    },
    /// Channel (C-shape) — open back of the C faces the +X direction
    /// in the beam's local frame.
    Channel {
        /// Overall width.
        width: f64,
        /// Overall depth.
        depth: f64,
        /// Thickness of the back wall + flanges.
        thickness: f64,
    },
}

impl BeamSection {
    /// Short human label.
    pub fn label(&self) -> &'static str {
        match self {
            BeamSection::Rectangle { .. } => "Rect",
            BeamSection::IBeam { .. } => "IBeam",
            BeamSection::Channel { .. } => "Channel",
        }
    }

    /// Outer bounding rectangle (width, depth) used for the v1
    /// tessellation.
    pub fn outer_box(&self) -> (f64, f64) {
        match self {
            BeamSection::Rectangle { width, depth } => (*width, *depth),
            BeamSection::IBeam { width, depth, .. } => (*width, *depth),
            BeamSection::Channel { width, depth, .. } => (*width, *depth),
        }
    }
}

/// Parameters describing a beam.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BeamParams {
    /// Start of the centre line, in world space.
    pub start: Vector3<f64>,
    /// End of the centre line, in world space.
    pub end: Vector3<f64>,
    /// Cross-section profile.
    pub cross_section: BeamSection,
    /// Rotation of the cross-section about the beam axis, in radians.
    /// 0.0 keeps the cross-section's "depth" axis vertical (+Z) when
    /// the beam is horizontal.
    pub orientation_angle: f64,
    /// Material descriptor.
    pub material: String,
    /// Optional structural attributes — material grade + applied
    /// distributed load. `None` for non-structural / aesthetic beams.
    /// Consumed by [`crate::structural::export_structural_model`].
    #[serde(default)]
    pub structural: Option<crate::structural::StructuralMember>,
}

impl BeamParams {
    /// Validate dimensions.
    pub fn validate(&self) -> Result<(), ArchError> {
        if self.length() < 1e-9 {
            return Err(ArchError::BadDimension {
                name: "length",
                reason: format!("start and end coincide (length {})", self.length()),
            });
        }
        let (w, d) = self.cross_section.outer_box();
        if w <= 0.0 || d <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "cross_section",
                reason: format!("width/depth must be > 0 (got {w}, {d})"),
            });
        }
        Ok(())
    }

    /// Beam length (Euclidean distance from `start` to `end`).
    pub fn length(&self) -> f64 {
        (self.end - self.start).norm()
    }

    /// Local frame: (axis, side, up). `axis` is the unit vector from
    /// `start` to `end`. `up` is chosen so the cross-section's depth
    /// axis points along +Z when possible (when beam is vertical,
    /// `up` defaults to +Y). `side` = axis × up rotated by
    /// `orientation_angle`.
    fn local_frame(&self) -> (Vector3<f64>, Vector3<f64>, Vector3<f64>) {
        // A zero-length beam (start == end) would make normalize() divide by a
        // zero norm and yield NaN corners; fall back to a finite default frame.
        let axis = (self.end - self.start)
            .try_normalize(1e-12)
            .unwrap_or_else(|| Vector3::new(1.0, 0.0, 0.0));
        let world_up = Vector3::new(0.0, 0.0, 1.0);
        let mut up = if (axis.dot(&world_up)).abs() > 0.99 {
            // Beam ~parallel to Z, fall back to +Y.
            Vector3::new(0.0, 1.0, 0.0)
        } else {
            world_up
        };
        let side = axis
            .cross(&up)
            .try_normalize(1e-12)
            .unwrap_or_else(|| Vector3::new(0.0, 1.0, 0.0));
        up = side.cross(&axis).try_normalize(1e-12).unwrap_or(up);
        // Apply orientation_angle rotation in the (side, up) plane.
        let (sa, ca) = self.orientation_angle.sin_cos();
        let side_rot = side * ca + up * sa;
        let up_rot = -side * sa + up * ca;
        (axis, side_rot, up_rot)
    }

    /// Compute the 8 corner positions of the beam's outer bounding
    /// box. Used for tessellation + bbox.
    pub fn corners(&self) -> [Vector3<f64>; 8] {
        let (axis, side, up) = self.local_frame();
        let (w, d) = self.cross_section.outer_box();
        let len = self.length();
        let hw = side * (w * 0.5);
        let hd = up * (d * 0.5);
        let s = self.start;
        let e = s + axis * len;
        [
            s - hw - hd,
            s + hw - hd,
            s + hw + hd,
            s - hw + hd,
            e - hw - hd,
            e + hw - hd,
            e + hw + hd,
            e - hw + hd,
        ]
    }

    /// Tessellate to a [`valenx_mesh::Mesh`] — same 6-face box layout
    /// as [`crate::WallParams::tessellate_mesh`] but oriented by the
    /// local frame.
    pub fn tessellate_mesh(&self) -> Result<Mesh, ArchError> {
        self.validate()?;
        let c = self.corners();
        let mut mesh = Mesh::new("beam");
        for v in &c {
            mesh.nodes.push(*v);
        }
        let mut block = ElementBlock::new(ElementType::Tri3);
        let quads: [[u32; 4]; 6] = [
            [0, 1, 2, 3],
            [4, 7, 6, 5],
            [1, 5, 6, 2],
            [0, 3, 7, 4],
            [0, 4, 5, 1],
            [3, 2, 6, 7],
        ];
        for q in quads {
            block.connectivity.extend_from_slice(&[q[0], q[1], q[2]]);
            block.connectivity.extend_from_slice(&[q[0], q[2], q[3]]);
        }
        mesh.element_blocks.push(block);
        mesh.recompute_stats();
        Ok(mesh)
    }

    /// Wrap [`Self::tessellate_mesh`] in a mesh-backed solid.
    pub fn tessellate(&self) -> Result<valenx_cad::Solid, ArchError> {
        let m = self.tessellate_mesh()?;
        Ok(valenx_cad::Solid::from_mesh(m))
    }

    /// Hint points for [`crate::ArchDocument::bbox`].
    pub fn bbox_hint_points(&self) -> impl Iterator<Item = Vector3<f64>> + '_ {
        let c = self.corners();
        c.into_iter().collect::<Vec<_>>().into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_beam() -> BeamParams {
        BeamParams {
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
            structural: None,
        }
    }

    #[test]
    fn length_5m() {
        let b = sample_beam();
        assert!((b.length() - 5.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_degenerate_length() {
        let mut b = sample_beam();
        b.end = b.start;
        assert!(matches!(
            b.validate(),
            Err(ArchError::BadDimension { name: "length", .. })
        ));
    }

    #[test]
    fn rejects_bad_cross_section() {
        let b = BeamParams {
            start: Vector3::zeros(),
            end: Vector3::new(1.0, 0.0, 0.0),
            cross_section: BeamSection::Rectangle {
                width: 0.0,
                depth: 0.1,
            },
            orientation_angle: 0.0,
            material: "x".into(),
            structural: None,
        };
        assert!(matches!(
            b.validate(),
            Err(ArchError::BadDimension {
                name: "cross_section",
                ..
            })
        ));
    }

    #[test]
    fn tessellation_has_12_triangles() {
        let b = sample_beam();
        let m = b.tessellate_mesh().unwrap();
        assert_eq!(m.nodes.len(), 8);
        assert_eq!(m.total_elements(), 12);
    }

    #[test]
    fn vertical_beam_uses_y_fallback() {
        // A beam along +Z should not panic on the local-frame
        // computation.
        let b = BeamParams {
            start: Vector3::zeros(),
            end: Vector3::new(0.0, 0.0, 3.0),
            cross_section: BeamSection::Rectangle {
                width: 0.1,
                depth: 0.1,
            },
            orientation_angle: 0.0,
            material: "x".into(),
            structural: None,
        };
        let m = b.tessellate_mesh().unwrap();
        assert_eq!(m.total_elements(), 12);
    }
}
