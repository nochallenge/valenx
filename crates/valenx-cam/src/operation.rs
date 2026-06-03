//! The four v1 operation types + per-op params.
//!
//! Every operation borrows `tool_id` to look up the [`crate::Tool`] in
//! the host's tool table; carries feeds and a step-down / step-over;
//! plus op-specific fields (e.g. drill `hole_positions`).
//!
//! The actual generation logic lives in [`crate::op`]; this module
//! only defines the data shapes so the persistence layer
//! ([`crate::persist::CamFile`]) can round-trip them.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Strategy for filling a pocket (or face) area.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub enum PocketStrategy {
    /// Back-and-forth raster lines spaced by `step_over`.
    #[default]
    ZigZag,
    /// Parallel one-way lines (raster but lift between passes).
    Parallel,
    /// Concentric inward offset polygons.
    Spiral,
}

impl PocketStrategy {
    /// Short panel label.
    pub fn label(self) -> &'static str {
        match self {
            PocketStrategy::ZigZag => "ZigZag",
            PocketStrategy::Parallel => "Parallel",
            PocketStrategy::Spiral => "Spiral",
        }
    }
}

/// Profile-op params. Cuts the boundary of the source mesh's
/// cross-section at each Z step-down.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ProfileParams {
    /// Tool id (looked up in the host's [`crate::Tool`] table).
    pub tool_id: u32,
    /// Cutting feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min) — used for Z-decreasing moves.
    pub plunge_feed: f64,
    /// Spindle RPM (referenced by the postprocessor `M3` header).
    pub spindle_rpm: f64,
    /// Z step-down per pass (mm, must be > 0).
    pub step_down: f64,
    /// Total cut depth below `stock.top_z()` (mm, must be > 0).
    pub depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
    /// `true` = climb cut (CCW for outside profile), `false` =
    /// conventional cut (CW). Reverses the polygon winding.
    pub climb: bool,
}

impl Default for ProfileParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 600.0,
            plunge_feed: 200.0,
            spindle_rpm: 12000.0,
            step_down: 1.0,
            depth: 5.0,
            safe_z_clearance: 5.0,
            climb: true,
        }
    }
}

/// Pocket-op params. Hollows out the cross-section interior at each
/// Z step-down, using [`PocketStrategy`] to fill the area.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PocketParams {
    /// Tool id.
    pub tool_id: u32,
    /// Cutting feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// Z step-down per pass (mm).
    pub step_down: f64,
    /// XY step-over between raster lines or spiral rings (mm).
    /// Must be ≤ `tool.diameter * 0.5` for safe engagement.
    pub step_over: f64,
    /// Total cut depth below `stock.top_z()` (mm).
    pub depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
    /// Raster angle in degrees (only used by ZigZag / Parallel).
    pub raster_angle_deg: f64,
    /// Fill strategy.
    pub strategy: PocketStrategy,
    /// Climb-vs-conventional (raster direction reversal).
    pub climb: bool,
}

impl Default for PocketParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 600.0,
            plunge_feed: 200.0,
            spindle_rpm: 12000.0,
            step_down: 1.0,
            step_over: 2.0,
            depth: 5.0,
            safe_z_clearance: 5.0,
            raster_angle_deg: 0.0,
            strategy: PocketStrategy::ZigZag,
            climb: true,
        }
    }
}

/// Drill-op params. Pecks vertical holes at each XY position.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DrillParams {
    /// Tool id (should reference a `ToolKind::Drill`).
    pub tool_id: u32,
    /// Plunge feed (mm/min) — also used for the cut feed.
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// Per-peck depth (mm). Must be > 0.
    pub peck_depth: f64,
    /// Total drill depth from stock top (mm). Must be > 0.
    pub total_depth: f64,
    /// Retract height above stock top between pecks (mm).
    pub retract_clearance: f64,
    /// Safe-Z clearance above stock top for rapid moves (mm).
    pub safe_z_clearance: f64,
    /// XY positions of the holes (Z is ignored; stock top is used).
    pub hole_positions: Vec<Vector3<f64>>,
}

impl Default for DrillParams {
    fn default() -> Self {
        Self {
            tool_id: 2,
            plunge_feed: 100.0,
            spindle_rpm: 1500.0,
            peck_depth: 1.0,
            total_depth: 5.0,
            retract_clearance: 1.0,
            safe_z_clearance: 5.0,
            hole_positions: Vec::new(),
        }
    }
}

/// Face-op params. Levels the top of the stock by parallel raster.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FaceParams {
    /// Tool id (face-mill or large end-mill).
    pub tool_id: u32,
    /// Cutting feed (mm/min).
    pub feed_mm_per_min: f64,
    /// Plunge feed (mm/min).
    pub plunge_feed: f64,
    /// Spindle RPM.
    pub spindle_rpm: f64,
    /// XY step-over between passes (mm).
    pub step_over: f64,
    /// Z step-down per pass (mm).
    pub step_down: f64,
    /// Total face depth below `stock.top_z()` (mm). Often 0.5–1.0 mm.
    pub depth: f64,
    /// Safe-Z clearance above stock top (mm).
    pub safe_z_clearance: f64,
    /// Raster angle in degrees.
    pub raster_angle_deg: f64,
    /// Climb-vs-conventional.
    pub climb: bool,
}

impl Default for FaceParams {
    fn default() -> Self {
        Self {
            tool_id: 1,
            feed_mm_per_min: 800.0,
            plunge_feed: 200.0,
            spindle_rpm: 12000.0,
            step_over: 4.0,
            step_down: 0.5,
            depth: 0.5,
            safe_z_clearance: 5.0,
            raster_angle_deg: 0.0,
            climb: true,
        }
    }
}

/// One CAM operation. The host stores a `Vec<Operation>` and runs
/// them in order to produce a chained [`crate::Toolpath`].
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Operation {
    /// Boundary profile.
    Profile(ProfileParams),
    /// Area pocket.
    Pocket(PocketParams),
    /// Drill cycle.
    Drill(DrillParams),
    /// Face / surface mill.
    Face(FaceParams),
    // Phase 17A — Adaptive clearing + new entry primitives.
    /// High-MRR trochoidal pocketing.
    AdaptiveClearing(crate::op::adaptive_clearing::AdaptiveParams),
    /// Helical bore.
    HelicalBore(crate::op::helix_bore::HelicalBoreParams),
    /// Plunge rough.
    PlungeRough(crate::op::plunge_rough::PlungeRoughParams),
    /// Ramp entry.
    RampEntry(crate::op::ramp_entry::RampEntryParams),
    /// Peck-drill full (dwell + chip retract).
    PeckDrillFull(crate::op::peck_drill_full::PeckDrillFullParams),
    // Phase 17B — More 2D + 3D ops.
    /// Contour-2D (XY curve at constant Z).
    Contour2D(crate::op::contour_2d::Contour2DParams),
    /// Contour-3D (3D curve at its actual Z).
    Contour3D(crate::op::contour_3d::Contour3DParams),
    /// V-bit engrave.
    Engrave(crate::op::engrave::EngraveParams),
    /// Shallow scribe.
    Scribe(crate::op::scribe::ScribeParams),
    /// Archimedean-spiral pocket.
    SpiralPocket(crate::op::spiral_pocket::SpiralPocketParams),
    /// Trochoidal slot.
    TrochoidalSlot(crate::op::trochoidal_slot::TrochoidalSlotParams),
    /// 3D waterline.
    Waterline3D(crate::op::waterline_3d::Waterline3DParams),
    /// Straight slot.
    Slot(crate::op::slot::SlotParams),
    /// Thread mill.
    ThreadMill(crate::op::thread_mill::ThreadMillParams),
    /// Rest machining.
    RestMachining(crate::op::rest_machining::RestMachiningParams),
}

impl Operation {
    /// Short label for the panel + audit log.
    pub fn label(&self) -> &'static str {
        match self {
            Operation::Profile(_) => "Profile",
            Operation::Pocket(_) => "Pocket",
            Operation::Drill(_) => "Drill",
            Operation::Face(_) => "Face",
            Operation::AdaptiveClearing(_) => "Adaptive Clearing",
            Operation::HelicalBore(_) => "Helical Bore",
            Operation::PlungeRough(_) => "Plunge Rough",
            Operation::RampEntry(_) => "Ramp Entry",
            Operation::PeckDrillFull(_) => "Peck Drill Full",
            Operation::Contour2D(_) => "Contour 2D",
            Operation::Contour3D(_) => "Contour 3D",
            Operation::Engrave(_) => "Engrave",
            Operation::Scribe(_) => "Scribe",
            Operation::SpiralPocket(_) => "Spiral Pocket",
            Operation::TrochoidalSlot(_) => "Trochoidal Slot",
            Operation::Waterline3D(_) => "Waterline 3D",
            Operation::Slot(_) => "Slot",
            Operation::ThreadMill(_) => "Thread Mill",
            Operation::RestMachining(_) => "Rest Machining",
        }
    }

    /// Which tool id this op references.
    pub fn tool_id(&self) -> u32 {
        match self {
            Operation::Profile(p) => p.tool_id,
            Operation::Pocket(p) => p.tool_id,
            Operation::Drill(p) => p.tool_id,
            Operation::Face(p) => p.tool_id,
            Operation::AdaptiveClearing(p) => p.tool_id,
            Operation::HelicalBore(p) => p.tool_id,
            Operation::PlungeRough(p) => p.tool_id,
            Operation::RampEntry(p) => p.tool_id,
            Operation::PeckDrillFull(p) => p.tool_id,
            Operation::Contour2D(p) => p.tool_id,
            Operation::Contour3D(p) => p.tool_id,
            Operation::Engrave(p) => p.tool_id,
            Operation::Scribe(p) => p.tool_id,
            Operation::SpiralPocket(p) => p.tool_id,
            Operation::TrochoidalSlot(p) => p.tool_id,
            Operation::Waterline3D(p) => p.tool_id,
            Operation::Slot(p) => p.tool_id,
            Operation::ThreadMill(p) => p.tool_id,
            Operation::RestMachining(p) => p.tool_id,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labels() {
        assert_eq!(
            Operation::Profile(ProfileParams::default()).label(),
            "Profile"
        );
        assert_eq!(Operation::Pocket(PocketParams::default()).label(), "Pocket");
        assert_eq!(Operation::Drill(DrillParams::default()).label(), "Drill");
        assert_eq!(Operation::Face(FaceParams::default()).label(), "Face");
    }

    #[test]
    fn tool_ids() {
        assert_eq!(Operation::Profile(ProfileParams::default()).tool_id(), 1);
        assert_eq!(Operation::Drill(DrillParams::default()).tool_id(), 2);
    }

    #[test]
    fn strategy_label() {
        assert_eq!(PocketStrategy::ZigZag.label(), "ZigZag");
        assert_eq!(PocketStrategy::Spiral.label(), "Spiral");
    }
}
