//! Top-level [`ArchEntity`] enum and kind discriminant.

use serde::{Deserialize, Serialize};

use crate::beam::BeamParams;
use crate::column::ColumnParams;
use crate::door::DoorParams;
use crate::mep::{
    CableSegmentParams, ConduitSegmentParams, DuctSegmentParams, MepEquipmentParams,
    PipeSegmentParams,
};
use crate::roof::RoofParams;
use crate::slab::SlabParams;
use crate::space::SpaceParams;
use crate::stair::StairParams;
use crate::wall::WallParams;
use crate::window::WindowParams;

/// A single element of an [`super::ArchDocument`].
///
/// All variants carry a strongly-typed parameter struct (see each
/// `Params` definition for fields + units). Geometry is computed
/// on-demand by [`crate::wall::WallParams::tessellate`] and friends —
/// the document only stores the parametric description.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum ArchEntity {
    /// A wall segment between two points.
    Wall(WallParams),
    /// A flat slab (floor / ceiling) defined by a closed boundary.
    Slab(SlabParams),
    /// A vertical column.
    Column(ColumnParams),
    /// A beam between two points.
    Beam(BeamParams),
    /// A window in a host wall.
    Window(WindowParams),
    /// A door in a host wall.
    Door(DoorParams),
    /// A staircase.
    Stair(StairParams),
    /// A roof over a footprint.
    Roof(RoofParams),
    /// A named space / room enclosing a footprint.
    Space(SpaceParams),
    /// HVAC duct segment ([`crate::DuctSegmentParams`]).
    DuctSegment(DuctSegmentParams),
    /// Plumbing / process pipe segment ([`crate::PipeSegmentParams`]).
    PipeSegment(PipeSegmentParams),
    /// Electrical cable run ([`crate::CableSegmentParams`]).
    CableSegment(CableSegmentParams),
    /// Electrical conduit ([`crate::ConduitSegmentParams`]).
    ConduitSegment(ConduitSegmentParams),
    /// MEP equipment placement ([`crate::MepEquipmentParams`]).
    MepEquipment(MepEquipmentParams),
}

/// Discriminant for [`ArchEntity`] without the heavy param payload.
///
/// Used by [`crate::Schedule`] (BOM grouping) and the UI panel for
/// concise list rendering.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd, Serialize, Deserialize)]
pub enum ArchEntityKind {
    /// [`ArchEntity::Wall`].
    Wall,
    /// [`ArchEntity::Slab`].
    Slab,
    /// [`ArchEntity::Column`].
    Column,
    /// [`ArchEntity::Beam`].
    Beam,
    /// [`ArchEntity::Window`].
    Window,
    /// [`ArchEntity::Door`].
    Door,
    /// [`ArchEntity::Stair`].
    Stair,
    /// [`ArchEntity::Roof`].
    Roof,
    /// [`ArchEntity::Space`].
    Space,
    /// [`ArchEntity::DuctSegment`].
    DuctSegment,
    /// [`ArchEntity::PipeSegment`].
    PipeSegment,
    /// [`ArchEntity::CableSegment`].
    CableSegment,
    /// [`ArchEntity::ConduitSegment`].
    ConduitSegment,
    /// [`ArchEntity::MepEquipment`].
    MepEquipment,
}

impl ArchEntityKind {
    /// Stable short label used in lists / schedules / IFC entity
    /// names.
    pub fn label(self) -> &'static str {
        match self {
            ArchEntityKind::Wall => "Wall",
            ArchEntityKind::Slab => "Slab",
            ArchEntityKind::Column => "Column",
            ArchEntityKind::Beam => "Beam",
            ArchEntityKind::Window => "Window",
            ArchEntityKind::Door => "Door",
            ArchEntityKind::Stair => "Stair",
            ArchEntityKind::Roof => "Roof",
            ArchEntityKind::Space => "Space",
            ArchEntityKind::DuctSegment => "Duct",
            ArchEntityKind::PipeSegment => "Pipe",
            ArchEntityKind::CableSegment => "Cable",
            ArchEntityKind::ConduitSegment => "Conduit",
            ArchEntityKind::MepEquipment => "Equipment",
        }
    }
}

impl ArchEntity {
    /// Tessellate this entity in the context of a document, honoring
    /// opening cuts for hosted windows / doors on walls.
    ///
    /// `walls` is a lookup `host_id → &WallParams` so window / door
    /// entities can find their host. Walls aggregate their hosted
    /// openings by scanning every entity in the doc — but to avoid
    /// re-walking the doc for every wall here, the public entry
    /// [`crate::ArchDocument::tessellate_all`] pre-builds the wall
    /// map and the per-wall opening lists.
    ///
    /// Windows and doors render their visible leaf / pane (not the
    /// void); the void only matters when computing the host wall's
    /// pierced geometry, which `ArchDocument::tessellate_all` handles
    /// before invoking this method.
    pub fn tessellate_in_doc(
        &self,
        _tolerance: f64,
        walls: &std::collections::HashMap<usize, &crate::wall::WallParams>,
    ) -> Result<valenx_mesh::Mesh, String> {
        match self {
            ArchEntity::Wall(w) => w.tessellate_mesh().map_err(|e| e.to_string()),
            ArchEntity::Slab(s) => s.tessellate_mesh().map_err(|e| e.to_string()),
            ArchEntity::Column(c) => c.tessellate_mesh().map_err(|e| e.to_string()),
            ArchEntity::Beam(b) => b.tessellate_mesh().map_err(|e| e.to_string()),
            ArchEntity::Window(w) => {
                if let Some(host) = walls.get(&w.host) {
                    w.tessellate_in_wall(host).map_err(|e| e.to_string())
                } else {
                    Err(format!("window references unknown host wall {}", w.host))
                }
            }
            ArchEntity::Door(d) => {
                if let Some(host) = walls.get(&d.host) {
                    d.tessellate_in_wall(host).map_err(|e| e.to_string())
                } else {
                    Err(format!("door references unknown host wall {}", d.host))
                }
            }
            ArchEntity::Stair(s) => s.tessellate_mesh().map_err(|e| e.to_string()),
            ArchEntity::Roof(r) => r.tessellate_mesh().map_err(|e| e.to_string()),
            ArchEntity::Space(s) => s.tessellate_mesh().map_err(|e| e.to_string()),
            ArchEntity::DuctSegment(d) => d.tessellate_mesh().map_err(|e| e.to_string()),
            ArchEntity::PipeSegment(p) => p.tessellate_mesh().map_err(|e| e.to_string()),
            ArchEntity::CableSegment(c) => c.tessellate_mesh().map_err(|e| e.to_string()),
            ArchEntity::ConduitSegment(c) => c.tessellate_mesh().map_err(|e| e.to_string()),
            ArchEntity::MepEquipment(e) => e.tessellate_mesh().map_err(|e| e.to_string()),
        }
    }

    /// Iterator over points that bound this entity. Used by
    /// [`crate::ArchDocument::bbox`]. Each implementor materialises
    /// into a Vec to keep the trait object-safe and the borrow
    /// lifetimes simple.
    pub fn bbox_hint_points(&self) -> Vec<nalgebra::Vector3<f64>> {
        match self {
            ArchEntity::Wall(w) => w.bbox_hint_points().collect(),
            ArchEntity::Slab(s) => s.bbox_hint_points().collect(),
            ArchEntity::Column(c) => c.bbox_hint_points().collect(),
            ArchEntity::Beam(b) => b.bbox_hint_points().collect(),
            ArchEntity::Window(_) | ArchEntity::Door(_) => Vec::new(),
            ArchEntity::Stair(s) => s.bbox_hint_points().collect(),
            ArchEntity::Roof(r) => r.bbox_hint_points().collect(),
            ArchEntity::Space(s) => s.bbox_hint_points().collect(),
            ArchEntity::DuctSegment(d) => d.bbox_hint_points().collect(),
            ArchEntity::PipeSegment(p) => p.bbox_hint_points().collect(),
            ArchEntity::CableSegment(c) => c.bbox_hint_points().collect(),
            ArchEntity::ConduitSegment(c) => c.bbox_hint_points().collect(),
            ArchEntity::MepEquipment(e) => e.bbox_hint_points().collect(),
        }
    }

    /// Discriminant for this entity (BOM grouping, UI lists).
    pub fn kind(&self) -> ArchEntityKind {
        match self {
            ArchEntity::Wall(_) => ArchEntityKind::Wall,
            ArchEntity::Slab(_) => ArchEntityKind::Slab,
            ArchEntity::Column(_) => ArchEntityKind::Column,
            ArchEntity::Beam(_) => ArchEntityKind::Beam,
            ArchEntity::Window(_) => ArchEntityKind::Window,
            ArchEntity::Door(_) => ArchEntityKind::Door,
            ArchEntity::Stair(_) => ArchEntityKind::Stair,
            ArchEntity::Roof(_) => ArchEntityKind::Roof,
            ArchEntity::Space(_) => ArchEntityKind::Space,
            ArchEntity::DuctSegment(_) => ArchEntityKind::DuctSegment,
            ArchEntity::PipeSegment(_) => ArchEntityKind::PipeSegment,
            ArchEntity::CableSegment(_) => ArchEntityKind::CableSegment,
            ArchEntity::ConduitSegment(_) => ArchEntityKind::ConduitSegment,
            ArchEntity::MepEquipment(_) => ArchEntityKind::MepEquipment,
        }
    }

    /// One-line "id-less" summary for the entity list (e.g.
    /// `"Wall 5.00 m × 2.70 m (Concrete)"`).
    pub fn summary(&self) -> String {
        match self {
            ArchEntity::Wall(w) => format!(
                "Wall {:.2} m × {:.2} m ({})",
                (w.end - w.start).norm(),
                w.height,
                w.material
            ),
            ArchEntity::Slab(s) => format!(
                "Slab ({} pts, thk {:.2} m, {})",
                s.boundary.len(),
                s.thickness,
                s.material
            ),
            ArchEntity::Column(c) => format!(
                "Column h={:.2} m ({}, {})",
                c.height,
                c.cross_section.label(),
                c.material
            ),
            ArchEntity::Beam(b) => format!(
                "Beam {:.2} m ({}, {})",
                (b.end - b.start).norm(),
                b.cross_section.label(),
                b.material
            ),
            ArchEntity::Window(w) => format!(
                "Window {:.2}×{:.2} m on wall {} ({})",
                w.width,
                w.height,
                w.host,
                w.style.label()
            ),
            ArchEntity::Door(d) => format!(
                "Door {:.2}×{:.2} m on wall {} ({})",
                d.width,
                d.height,
                d.host,
                d.style.label()
            ),
            ArchEntity::Stair(s) => format!(
                "Stair {} steps, rise {:.2} m × run {:.2} m",
                s.num_steps, s.total_rise, s.total_run
            ),
            ArchEntity::Roof(r) => format!(
                "Roof ({} pts, {}, peak {:.2} m)",
                r.boundary.len(),
                r.roof_type.label(),
                r.peak_height
            ),
            ArchEntity::Space(s) => format!(
                "Space \"{}\" ({} pts, h {:.2} m)",
                s.space_name,
                s.boundary.len(),
                s.ceiling_height
            ),
            ArchEntity::DuctSegment(d) => format!(
                "Duct {} {:.2} m ({}, {})",
                d.shape.label(),
                d.length(),
                d.flow_direction.label(),
                d.material
            ),
            ArchEntity::PipeSegment(p) => format!(
                "Pipe ø{:.3} m × {:.2} m ({}, {})",
                p.diameter,
                p.length(),
                p.fluid,
                p.material
            ),
            ArchEntity::CableSegment(c) => format!(
                "Cable {:.1} mm² × {:.2} m @ {:.0} V ({})",
                c.conductor_csa_mm2,
                c.length(),
                c.voltage,
                c.material
            ),
            ArchEntity::ConduitSegment(c) => format!(
                "Conduit ø{:.3}/{:.3} m × {:.2} m ({})",
                c.outer_diameter,
                c.inner_diameter,
                c.length(),
                c.material
            ),
            ArchEntity::MepEquipment(e) => format!(
                "Equipment {} \"{}\" ({:.2}×{:.2}×{:.2} m)",
                e.kind.label(),
                e.tag,
                e.size[0],
                e.size[1],
                e.size[2]
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn kind_and_summary_for_each_variant() {
        let w = ArchEntity::Wall(WallParams {
            start: Vector3::zeros(),
            end: Vector3::new(3.0, 0.0, 0.0),
            height: 2.5,
            thickness: 0.2,
            material: "Brick".into(),
        });
        assert_eq!(w.kind(), ArchEntityKind::Wall);
        assert!(w.summary().contains("Wall"));
        assert!(w.summary().contains("Brick"));
        assert_eq!(ArchEntityKind::Wall.label(), "Wall");
    }
}
