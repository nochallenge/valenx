//! Mechanical / electrical / plumbing (MEP) entities.
//!
//! Production BIM tools (Revit MEP, ArchiCAD MEP Modeler, FreeCAD's
//! Arch MEP add-on) carry first-class entities for ducting, piping,
//! electrical conduit, and discipline equipment placements. This
//! module ships representative entities for the four discipline
//! categories plus an `MepEquipment` generic placement:
//!
//! - [`DuctSegmentParams`] — HVAC supply / return / exhaust ducts;
//!   round, rectangular, or oval cross-section.
//! - [`PipeSegmentParams`] — plumbing / process pipes; round
//!   cross-section, carries fluid name and operating pressure.
//! - [`CableSegmentParams`] — electrical cables; gauge + voltage
//!   class.
//! - [`ConduitSegmentParams`] — electrical conduit ducts; nominal
//!   diameter + cable trays. The IFC4 distinction between conduit
//!   and cable matters because conduit is the carrier (physical
//!   tube) and cable is the conductor.
//! - [`MepEquipmentParams`] — a generic placement for AHUs,
//!   electrical panels, pumps, valves, sprinkler heads, etc., with
//!   an [`EquipmentKind`] discriminant.
//!
//! Each entity carries a parametric path (`start`/`end` for the four
//! segments) or anchor + bounding-box (for equipment) and tessellates
//! to a swept profile / box. The IFC4 writer emits the matching
//! `IfcDuctSegment` / `IfcPipeSegment` / `IfcCableSegment` /
//! `IfcConduitSegment` / `IfcDistributionElement` entity per the
//! IFC4 schema's `IfcSharedBldgServiceElements` and
//! `IfcDistributionFlowElements` packages.
//!
//! ## Honest scope
//!
//! Single-segment cylindrical / rectangular swept solids — fittings
//! (elbows, tees, transitions) are represented by separate segments
//! abutting at a node; a true MEP system carries fitting libraries
//! and connector ports. That's a follow-up. The schedule integration
//! aggregates length and (where applicable) flow-area for each
//! discipline.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::error::ArchError;

/// Duct cross-section shape (HVAC).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum DuctShape {
    /// Round duct — `diameter` in metres.
    Round {
        /// Outer diameter (m).
        diameter: f64,
    },
    /// Rectangular duct — `width × height` in metres.
    Rectangular {
        /// Width (m, horizontal in the local frame).
        width: f64,
        /// Height (m, vertical in the local frame).
        height: f64,
    },
    /// Oval (flat-oval) duct — `width × height` with rounded
    /// short-axis ends. v1 tessellates as the bounding rectangle and
    /// records the shape kind for IFC.
    Oval {
        /// Major axis (m).
        width: f64,
        /// Minor axis (m).
        height: f64,
    },
}

impl DuctShape {
    /// Short human label.
    pub fn label(self) -> &'static str {
        match self {
            DuctShape::Round { .. } => "Round",
            DuctShape::Rectangular { .. } => "Rect",
            DuctShape::Oval { .. } => "Oval",
        }
    }

    /// Bounding-box width × height for the v1 box tessellation.
    pub fn outer_box(self) -> (f64, f64) {
        match self {
            DuctShape::Round { diameter } => (diameter, diameter),
            DuctShape::Rectangular { width, height } => (width, height),
            DuctShape::Oval { width, height } => (width, height),
        }
    }

    /// Cross-section flow area in m².
    pub fn flow_area(self) -> f64 {
        match self {
            DuctShape::Round { diameter } => std::f64::consts::PI * (diameter * 0.5).powi(2),
            DuctShape::Rectangular { width, height } => width * height,
            // Oval area ≈ π·a·b for an ellipse with semi-axes a, b.
            DuctShape::Oval { width, height } => {
                std::f64::consts::PI * (width * 0.5) * (height * 0.5)
            }
        }
    }
}

/// HVAC flow direction along a segment — informational, used by IFC's
/// `IfcDuctSegmentType.PredefinedType` and by the schedule.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum FlowDirection {
    /// Air / fluid flows from `start` toward `end`.
    SourceToSink,
    /// Bidirectional / branch (no fixed flow direction).
    Bidirectional,
    /// Return — flows from `end` toward `start` in the larger
    /// system's topology (the segment is drawn `start → end` for
    /// authoring convenience).
    Return,
}

impl FlowDirection {
    /// Short label.
    pub fn label(self) -> &'static str {
        match self {
            FlowDirection::SourceToSink => "Supply",
            FlowDirection::Bidirectional => "Bidir",
            FlowDirection::Return => "Return",
        }
    }
}

/// Discipline kind discriminator for [`MepEquipmentParams`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum EquipmentKind {
    /// Air-handling unit (HVAC).
    AirHandlingUnit,
    /// Variable air volume terminal (HVAC).
    VavBox,
    /// Pump (HVAC / plumbing).
    Pump,
    /// Valve (HVAC / plumbing).
    Valve,
    /// Sprinkler head (fire-protection / plumbing).
    SprinklerHead,
    /// Electrical panel / distribution board.
    ElectricalPanel,
    /// Light fitting (electrical).
    LightFitting,
}

impl EquipmentKind {
    /// Short human label.
    pub fn label(self) -> &'static str {
        match self {
            EquipmentKind::AirHandlingUnit => "AHU",
            EquipmentKind::VavBox => "VAV",
            EquipmentKind::Pump => "Pump",
            EquipmentKind::Valve => "Valve",
            EquipmentKind::SprinklerHead => "Sprinkler",
            EquipmentKind::ElectricalPanel => "Panel",
            EquipmentKind::LightFitting => "Light",
        }
    }

    /// IFC4 entity type the equipment maps to in
    /// `crate::ifc::writer::write_mep_equipment`.
    pub fn ifc_entity(self) -> &'static str {
        match self {
            EquipmentKind::AirHandlingUnit => "IFCAIRTERMINALBOX",
            EquipmentKind::VavBox => "IFCAIRTERMINALBOX",
            EquipmentKind::Pump => "IFCPUMP",
            EquipmentKind::Valve => "IFCVALVE",
            EquipmentKind::SprinklerHead => "IFCFIRESUPPRESSIONTERMINAL",
            EquipmentKind::ElectricalPanel => "IFCELECTRICDISTRIBUTIONBOARD",
            EquipmentKind::LightFitting => "IFCLIGHTFIXTURE",
        }
    }
}

/// HVAC duct segment.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DuctSegmentParams {
    /// Start of the centre line, in world space.
    pub start: Vector3<f64>,
    /// End of the centre line.
    pub end: Vector3<f64>,
    /// Cross-section shape + size.
    pub shape: DuctShape,
    /// Wall material (e.g. "Galv steel", "Stainless steel", "PVC").
    pub material: String,
    /// Air flow direction along the segment.
    pub flow_direction: FlowDirection,
}

impl DuctSegmentParams {
    /// Validate the duct geometry.
    pub fn validate(&self) -> Result<(), ArchError> {
        if (self.end - self.start).norm() < 1e-9 {
            return Err(ArchError::BadDimension {
                name: "length",
                reason: format!(
                    "duct length must be > 0 (got {})",
                    (self.end - self.start).norm()
                ),
            });
        }
        let (w, h) = self.shape.outer_box();
        if w <= 0.0 || h <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "shape",
                reason: format!("dimensions must be > 0 (got {w}x{h})"),
            });
        }
        Ok(())
    }

    /// Path length in metres.
    pub fn length(&self) -> f64 {
        (self.end - self.start).norm()
    }

    /// Tessellate as a swept box along the centreline.
    pub fn tessellate_mesh(&self) -> Result<Mesh, ArchError> {
        self.validate()?;
        let (w, h) = self.shape.outer_box();
        swept_box_mesh("duct", self.start, self.end, w, h, 0.0)
    }

    /// Hint points for [`crate::ArchDocument::bbox`].
    pub fn bbox_hint_points(&self) -> impl Iterator<Item = Vector3<f64>> + '_ {
        let (w, h) = self.shape.outer_box();
        swept_box_corners(self.start, self.end, w, h).into_iter()
    }
}

/// Plumbing / process pipe segment.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PipeSegmentParams {
    /// Start of the centre line.
    pub start: Vector3<f64>,
    /// End of the centre line.
    pub end: Vector3<f64>,
    /// Nominal outer diameter (m).
    pub diameter: f64,
    /// Wall material ("Copper", "PEX", "Cast iron", ...).
    pub material: String,
    /// Working fluid name ("Cold water", "Hot water return", ...).
    pub fluid: String,
    /// Operating pressure in Pa — informational, used by IFC's
    /// `Pset_PipeSegmentTypeCommon.Pressure`.
    pub operating_pressure: f64,
}

impl PipeSegmentParams {
    /// Validate the pipe geometry.
    pub fn validate(&self) -> Result<(), ArchError> {
        if (self.end - self.start).norm() < 1e-9 {
            return Err(ArchError::BadDimension {
                name: "length",
                reason: format!(
                    "pipe length must be > 0 (got {})",
                    (self.end - self.start).norm()
                ),
            });
        }
        if !self.diameter.is_finite() || self.diameter <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "diameter",
                reason: format!("must be > 0 (got {})", self.diameter),
            });
        }
        Ok(())
    }

    /// Path length in metres.
    pub fn length(&self) -> f64 {
        (self.end - self.start).norm()
    }

    /// Flow area in m².
    pub fn flow_area(&self) -> f64 {
        std::f64::consts::PI * (self.diameter * 0.5).powi(2)
    }

    /// Tessellate as a swept square box (the bounding box of the
    /// round cross-section).
    pub fn tessellate_mesh(&self) -> Result<Mesh, ArchError> {
        self.validate()?;
        swept_box_mesh(
            "pipe",
            self.start,
            self.end,
            self.diameter,
            self.diameter,
            0.0,
        )
    }

    /// Hint points for [`crate::ArchDocument::bbox`].
    pub fn bbox_hint_points(&self) -> impl Iterator<Item = Vector3<f64>> + '_ {
        swept_box_corners(self.start, self.end, self.diameter, self.diameter).into_iter()
    }
}

/// Electrical cable run (the conductor).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CableSegmentParams {
    /// Start of the run.
    pub start: Vector3<f64>,
    /// End of the run.
    pub end: Vector3<f64>,
    /// Outer diameter of the cable bundle (m).
    pub diameter: f64,
    /// Conductor cross-section in mm² (AWG-equivalent area). Stored
    /// for the IFC `Pset_CableSegmentTypeCommon.CrossSectionalArea`.
    pub conductor_csa_mm2: f64,
    /// Nominal voltage class in volts ("400", "230", "12", ...).
    pub voltage: f64,
    /// Insulation / sheath material ("PVC", "XLPE", "Rubber", ...).
    pub material: String,
}

impl CableSegmentParams {
    /// Validate the cable geometry.
    pub fn validate(&self) -> Result<(), ArchError> {
        if (self.end - self.start).norm() < 1e-9 {
            return Err(ArchError::BadDimension {
                name: "length",
                reason: format!(
                    "cable length must be > 0 (got {})",
                    (self.end - self.start).norm()
                ),
            });
        }
        if !self.diameter.is_finite() || self.diameter <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "diameter",
                reason: format!("must be > 0 (got {})", self.diameter),
            });
        }
        if !self.conductor_csa_mm2.is_finite() || self.conductor_csa_mm2 <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "conductor_csa_mm2",
                reason: format!("must be > 0 (got {})", self.conductor_csa_mm2),
            });
        }
        Ok(())
    }

    /// Path length in metres.
    pub fn length(&self) -> f64 {
        (self.end - self.start).norm()
    }

    /// Tessellate as a swept square box.
    pub fn tessellate_mesh(&self) -> Result<Mesh, ArchError> {
        self.validate()?;
        swept_box_mesh(
            "cable",
            self.start,
            self.end,
            self.diameter,
            self.diameter,
            0.0,
        )
    }

    /// Hint points for [`crate::ArchDocument::bbox`].
    pub fn bbox_hint_points(&self) -> impl Iterator<Item = Vector3<f64>> + '_ {
        swept_box_corners(self.start, self.end, self.diameter, self.diameter).into_iter()
    }
}

/// Electrical conduit (the protective tube that cables run inside).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ConduitSegmentParams {
    /// Start of the conduit.
    pub start: Vector3<f64>,
    /// End of the conduit.
    pub end: Vector3<f64>,
    /// Outer diameter (m).
    pub outer_diameter: f64,
    /// Inner diameter (m) — used by the schedule's free-area
    /// reporting and the IFC `Pset_ConduitSegmentTypeCommon`.
    pub inner_diameter: f64,
    /// Material ("EMT", "RMC", "PVC", "Flex steel", ...).
    pub material: String,
}

impl ConduitSegmentParams {
    /// Validate the conduit geometry.
    pub fn validate(&self) -> Result<(), ArchError> {
        if (self.end - self.start).norm() < 1e-9 {
            return Err(ArchError::BadDimension {
                name: "length",
                reason: format!(
                    "conduit length must be > 0 (got {})",
                    (self.end - self.start).norm()
                ),
            });
        }
        if !self.outer_diameter.is_finite() || self.outer_diameter <= 0.0 {
            return Err(ArchError::BadDimension {
                name: "outer_diameter",
                reason: format!("must be > 0 (got {})", self.outer_diameter),
            });
        }
        if !self.inner_diameter.is_finite()
            || self.inner_diameter <= 0.0
            || self.inner_diameter >= self.outer_diameter
        {
            return Err(ArchError::BadDimension {
                name: "inner_diameter",
                reason: format!(
                    "must be > 0 and < outer (got inner={}, outer={})",
                    self.inner_diameter, self.outer_diameter
                ),
            });
        }
        Ok(())
    }

    /// Path length in metres.
    pub fn length(&self) -> f64 {
        (self.end - self.start).norm()
    }

    /// Free cross-section area (m²) — the conductor-fill calculation
    /// downstream tools care about.
    pub fn free_area(&self) -> f64 {
        std::f64::consts::PI * (self.inner_diameter * 0.5).powi(2)
    }

    /// Tessellate as a swept square box.
    pub fn tessellate_mesh(&self) -> Result<Mesh, ArchError> {
        self.validate()?;
        swept_box_mesh(
            "conduit",
            self.start,
            self.end,
            self.outer_diameter,
            self.outer_diameter,
            0.0,
        )
    }

    /// Hint points for [`crate::ArchDocument::bbox`].
    pub fn bbox_hint_points(&self) -> impl Iterator<Item = Vector3<f64>> + '_ {
        swept_box_corners(
            self.start,
            self.end,
            self.outer_diameter,
            self.outer_diameter,
        )
        .into_iter()
    }
}

/// Generic MEP equipment placement — an axis-aligned box anchored at
/// `position` with `size` extents, classified by [`EquipmentKind`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MepEquipmentParams {
    /// Lower-back-left corner of the equipment's bounding box.
    pub position: Vector3<f64>,
    /// `(width, depth, height)` extents in metres.
    pub size: [f64; 3],
    /// Discipline-kind discriminator (drives the IFC entity choice).
    pub kind: EquipmentKind,
    /// User tag / asset id ("AHU-101", "P-203", "VAV-3-12", ...).
    pub tag: String,
    /// Free-form property string used for additional Pset fields.
    pub description: String,
}

impl MepEquipmentParams {
    /// Validate the equipment geometry.
    pub fn validate(&self) -> Result<(), ArchError> {
        for (i, name) in ["width", "depth", "height"].iter().enumerate() {
            let s = self.size[i];
            if !s.is_finite() || s <= 0.0 {
                return Err(ArchError::BadDimension {
                    name: match *name {
                        "width" => "size.width",
                        "depth" => "size.depth",
                        "height" => "size.height",
                        _ => "size",
                    },
                    reason: format!("must be > 0 (got {s})"),
                });
            }
        }
        Ok(())
    }

    /// Bounding-box max corner (`position + size`).
    pub fn max_corner(&self) -> Vector3<f64> {
        self.position + Vector3::new(self.size[0], self.size[1], self.size[2])
    }

    /// Tessellate as the equipment's axis-aligned bounding box.
    pub fn tessellate_mesh(&self) -> Result<Mesh, ArchError> {
        self.validate()?;
        let p0 = self.position;
        let p1 = self.max_corner();
        let mut mesh = Mesh::new("mep_equipment");
        let corners = [
            Vector3::new(p0.x, p0.y, p0.z),
            Vector3::new(p1.x, p0.y, p0.z),
            Vector3::new(p1.x, p1.y, p0.z),
            Vector3::new(p0.x, p1.y, p0.z),
            Vector3::new(p0.x, p0.y, p1.z),
            Vector3::new(p1.x, p0.y, p1.z),
            Vector3::new(p1.x, p1.y, p1.z),
            Vector3::new(p0.x, p1.y, p1.z),
        ];
        for c in &corners {
            mesh.nodes.push(*c);
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

    /// Hint points for [`crate::ArchDocument::bbox`].
    pub fn bbox_hint_points(&self) -> impl Iterator<Item = Vector3<f64>> + '_ {
        let p0 = self.position;
        let p1 = self.max_corner();
        vec![p0, p1].into_iter()
    }
}

/// Build a swept-box mesh between `start` and `end` with the given
/// `width` / `height` cross-section. Used by every MEP segment that
/// renders as an axis-aligned bounding box around its centreline.
///
/// The box is oriented so its long axis lies along `start → end`,
/// with the local-y axis (width) horizontal and local-z (height)
/// vertical when the segment isn't itself vertical.
fn swept_box_mesh(
    name: &str,
    start: Vector3<f64>,
    end: Vector3<f64>,
    width: f64,
    height: f64,
    _roll: f64,
) -> Result<Mesh, ArchError> {
    let corners = swept_box_corners(start, end, width, height);
    let mut mesh = Mesh::new(name);
    for c in &corners {
        mesh.nodes.push(*c);
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

/// 8 corners of a box swept from `start` to `end` with `width` (local
/// y) and `height` (local z) cross-section. Ordering matches
/// [`crate::beam::BeamParams::corners`] so the swept-box mesh
/// connectivity stays in lockstep.
fn swept_box_corners(
    start: Vector3<f64>,
    end: Vector3<f64>,
    width: f64,
    height: f64,
) -> [Vector3<f64>; 8] {
    let axis = end - start;
    let len = axis.norm();
    if len < 1e-12 {
        return [start; 8];
    }
    let axis = axis / len;
    let world_up = Vector3::new(0.0, 0.0, 1.0);
    let mut up = if (axis.dot(&world_up)).abs() > 0.99 {
        Vector3::new(0.0, 1.0, 0.0)
    } else {
        world_up
    };
    let side = axis.cross(&up).normalize();
    up = side.cross(&axis).normalize();
    let hw = side * (width * 0.5);
    let hd = up * (height * 0.5);
    let s = start;
    let e = end;
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

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_duct() -> DuctSegmentParams {
        DuctSegmentParams {
            start: Vector3::zeros(),
            end: Vector3::new(5.0, 0.0, 0.0),
            shape: DuctShape::Round { diameter: 0.3 },
            material: "Galv steel".into(),
            flow_direction: FlowDirection::SourceToSink,
        }
    }

    #[test]
    fn duct_validation_and_length() {
        let d = sample_duct();
        d.validate().unwrap();
        assert!((d.length() - 5.0).abs() < 1e-9);
        assert!((d.shape.flow_area() - std::f64::consts::PI * 0.0225).abs() < 1e-9);
    }

    #[test]
    fn duct_rejects_zero_diameter() {
        let mut d = sample_duct();
        d.shape = DuctShape::Round { diameter: 0.0 };
        assert!(matches!(
            d.validate(),
            Err(ArchError::BadDimension { name: "shape", .. })
        ));
    }

    #[test]
    fn duct_tessellates_to_12_tris() {
        let d = sample_duct();
        let m = d.tessellate_mesh().unwrap();
        assert_eq!(m.nodes.len(), 8);
        assert_eq!(m.total_elements(), 12);
    }

    #[test]
    fn pipe_validation_and_area() {
        let p = PipeSegmentParams {
            start: Vector3::zeros(),
            end: Vector3::new(3.0, 0.0, 0.0),
            diameter: 0.05,
            material: "Copper".into(),
            fluid: "Cold water".into(),
            operating_pressure: 4.0e5,
        };
        p.validate().unwrap();
        assert!((p.length() - 3.0).abs() < 1e-9);
        // π · (0.025)² ≈ 1.963e-3
        assert!((p.flow_area() - std::f64::consts::PI * 0.025 * 0.025).abs() < 1e-12);
    }

    #[test]
    fn pipe_rejects_zero_length() {
        let p = PipeSegmentParams {
            start: Vector3::zeros(),
            end: Vector3::zeros(),
            diameter: 0.05,
            material: "Copper".into(),
            fluid: "Water".into(),
            operating_pressure: 1e5,
        };
        assert!(matches!(
            p.validate(),
            Err(ArchError::BadDimension { name: "length", .. })
        ));
    }

    #[test]
    fn cable_validation() {
        let c = CableSegmentParams {
            start: Vector3::zeros(),
            end: Vector3::new(10.0, 0.0, 0.0),
            diameter: 0.012,
            conductor_csa_mm2: 4.0,
            voltage: 230.0,
            material: "PVC".into(),
        };
        c.validate().unwrap();
        let m = c.tessellate_mesh().unwrap();
        assert_eq!(m.total_elements(), 12);
    }

    #[test]
    fn conduit_inner_must_be_smaller_than_outer() {
        let c = ConduitSegmentParams {
            start: Vector3::zeros(),
            end: Vector3::new(5.0, 0.0, 0.0),
            outer_diameter: 0.025,
            inner_diameter: 0.030, // larger than outer.
            material: "EMT".into(),
        };
        assert!(matches!(
            c.validate(),
            Err(ArchError::BadDimension {
                name: "inner_diameter",
                ..
            })
        ));
    }

    #[test]
    fn conduit_free_area_uses_inner() {
        let c = ConduitSegmentParams {
            start: Vector3::zeros(),
            end: Vector3::new(5.0, 0.0, 0.0),
            outer_diameter: 0.030,
            inner_diameter: 0.025,
            material: "EMT".into(),
        };
        c.validate().unwrap();
        // π · (0.0125)² ≈ 4.909e-4
        let expected = std::f64::consts::PI * 0.0125 * 0.0125;
        assert!((c.free_area() - expected).abs() < 1e-12);
    }

    #[test]
    fn equipment_validation_and_corners() {
        let e = MepEquipmentParams {
            position: Vector3::new(1.0, 2.0, 3.0),
            size: [1.5, 0.8, 1.2],
            kind: EquipmentKind::AirHandlingUnit,
            tag: "AHU-101".into(),
            description: "Roof-top air handler".into(),
        };
        e.validate().unwrap();
        let max = e.max_corner();
        assert!((max - Vector3::new(2.5, 2.8, 4.2)).norm() < 1e-9);
        let m = e.tessellate_mesh().unwrap();
        assert_eq!(m.nodes.len(), 8);
        assert_eq!(m.total_elements(), 12);
    }

    #[test]
    fn equipment_rejects_bad_size() {
        let e = MepEquipmentParams {
            position: Vector3::zeros(),
            size: [-1.0, 1.0, 1.0],
            kind: EquipmentKind::Pump,
            tag: "P-1".into(),
            description: "".into(),
        };
        assert!(matches!(e.validate(), Err(ArchError::BadDimension { .. })));
    }

    #[test]
    fn equipment_kind_ifc_mapping_is_stable() {
        for k in [
            EquipmentKind::AirHandlingUnit,
            EquipmentKind::VavBox,
            EquipmentKind::Pump,
            EquipmentKind::Valve,
            EquipmentKind::SprinklerHead,
            EquipmentKind::ElectricalPanel,
            EquipmentKind::LightFitting,
        ] {
            assert!(!k.label().is_empty());
            assert!(k.ifc_entity().starts_with("IFC"));
        }
    }

    #[test]
    fn flow_direction_label() {
        assert_eq!(FlowDirection::SourceToSink.label(), "Supply");
        assert_eq!(FlowDirection::Bidirectional.label(), "Bidir");
        assert_eq!(FlowDirection::Return.label(), "Return");
    }

    #[test]
    fn duct_shape_outer_box_and_area() {
        let r = DuctShape::Rectangular {
            width: 0.4,
            height: 0.3,
        };
        assert_eq!(r.outer_box(), (0.4, 0.3));
        assert!((r.flow_area() - 0.12).abs() < 1e-12);
        let o = DuctShape::Oval {
            width: 0.5,
            height: 0.3,
        };
        // ellipse area π·a·b = π · 0.25 · 0.15
        assert!((o.flow_area() - std::f64::consts::PI * 0.25 * 0.15).abs() < 1e-12);
    }
}
