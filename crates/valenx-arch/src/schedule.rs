//! Schedule (Bill-of-Materials-like) reporting over an
//! [`crate::ArchDocument`].
//!
//! Walks every entity in the document, groups by
//! [`crate::ArchEntityKind`], and accumulates a quantity per group:
//! - linear metres for walls + beams (sum of lengths)
//! - square metres for slabs + roofs (sum of footprint area)
//! - count for columns / windows / doors / stairs / spaces
//!
//! Output can be rendered to a CSV or a plain-text aligned table.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::document::ArchDocument;
use crate::entity::{ArchEntity, ArchEntityKind};

/// One row of the schedule.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct ScheduleEntry {
    /// Total count of entities of this kind.
    pub count: u32,
    /// Aggregated linear quantity (metres) — sum of wall/beam lengths.
    /// Zero for kinds where this metric doesn't apply.
    pub linear_m: f64,
    /// Aggregated area quantity (m²) — sum of slab/roof footprint
    /// area. Zero for kinds where this metric doesn't apply.
    pub area_m2: f64,
    /// Aggregated volume (m³) — currently used for spaces only.
    pub volume_m3: f64,
}

/// Full schedule for an [`ArchDocument`].
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Schedule {
    /// Per-kind aggregated entries, sorted by kind.
    pub entries: BTreeMap<ArchEntityKind, ScheduleEntry>,
}

impl Schedule {
    /// Build a schedule from a document.
    pub fn from_document(doc: &ArchDocument) -> Self {
        let mut s = Schedule::default();
        for (_id, ent) in &doc.entities {
            let kind = ent.kind();
            let row = s.entries.entry(kind).or_default();
            row.count += 1;
            match ent {
                ArchEntity::Wall(w) => row.linear_m += w.length(),
                ArchEntity::Beam(b) => row.linear_m += b.length(),
                ArchEntity::Slab(sl) => row.area_m2 += sl.area_m2(),
                ArchEntity::Roof(r) => {
                    // Use the boundary footprint area (same shoelace).
                    let n = r.boundary.len();
                    let mut a = 0.0;
                    for i in 0..n {
                        let j = (i + 1) % n;
                        a += r.boundary[i].x * r.boundary[j].y;
                        a -= r.boundary[j].x * r.boundary[i].y;
                    }
                    row.area_m2 += (a * 0.5).abs();
                }
                ArchEntity::Space(sp) => {
                    row.area_m2 += sp.floor_area();
                    row.volume_m3 += sp.volume();
                }
                ArchEntity::DuctSegment(d) => row.linear_m += d.length(),
                ArchEntity::PipeSegment(p) => row.linear_m += p.length(),
                ArchEntity::CableSegment(c) => row.linear_m += c.length(),
                ArchEntity::ConduitSegment(c) => row.linear_m += c.length(),
                ArchEntity::MepEquipment(e) => {
                    row.volume_m3 += e.size[0] * e.size[1] * e.size[2];
                }
                ArchEntity::Column(_)
                | ArchEntity::Window(_)
                | ArchEntity::Door(_)
                | ArchEntity::Stair(_) => {
                    // count-only.
                }
            }
        }
        s
    }

    /// Render to a UTF-8 CSV string. First row is the header.
    /// Columns: `kind, count, linear_m, area_m2, volume_m3`.
    pub fn to_csv(&self) -> String {
        let mut out = String::from("kind,count,linear_m,area_m2,volume_m3\n");
        for (k, e) in &self.entries {
            out.push_str(&format!(
                "{},{},{:.3},{:.3},{:.3}\n",
                k.label(),
                e.count,
                e.linear_m,
                e.area_m2,
                e.volume_m3
            ));
        }
        out
    }

    /// Render to a plain-text aligned table.
    pub fn to_text_table(&self) -> String {
        let mut out = String::new();
        out.push_str(&format!(
            "{:<10} {:>6} {:>10} {:>10} {:>10}\n",
            "kind", "count", "lin (m)", "area (m²)", "vol (m³)"
        ));
        for (k, e) in &self.entries {
            out.push_str(&format!(
                "{:<10} {:>6} {:>10.2} {:>10.2} {:>10.2}\n",
                k.label(),
                e.count,
                e.linear_m,
                e.area_m2,
                e.volume_m3
            ));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::door::{DoorParams, DoorStyle, Side};
    use crate::wall::WallParams;
    use nalgebra::Vector3;

    fn sample_doc() -> ArchDocument {
        let mut d = ArchDocument::new("house");
        // 4 walls of length 3.0 each.
        for (start, end) in [
            (Vector3::new(0.0, 0.0, 0.0), Vector3::new(3.0, 0.0, 0.0)),
            (Vector3::new(3.0, 0.0, 0.0), Vector3::new(3.0, 3.0, 0.0)),
            (Vector3::new(3.0, 3.0, 0.0), Vector3::new(0.0, 3.0, 0.0)),
            (Vector3::new(0.0, 3.0, 0.0), Vector3::new(0.0, 0.0, 0.0)),
        ] {
            d.add_entity(ArchEntity::Wall(WallParams {
                start,
                end,
                height: 2.5,
                thickness: 0.2,
                material: "Brick".into(),
            }));
        }
        // 1 door on wall id 1.
        d.add_entity(ArchEntity::Door(DoorParams {
            host: 1,
            position_along_wall: 1.5,
            width: 0.9,
            height: 2.1,
            style: DoorStyle::Single,
            hinge_side: Side::Left,
        }));
        d
    }

    #[test]
    fn four_walls_one_door_groups_correctly() {
        let s = Schedule::from_document(&sample_doc());
        let wall_row = s.entries.get(&ArchEntityKind::Wall).unwrap();
        assert_eq!(wall_row.count, 4);
        assert!((wall_row.linear_m - 12.0).abs() < 1e-9);
        let door_row = s.entries.get(&ArchEntityKind::Door).unwrap();
        assert_eq!(door_row.count, 1);
        // Door rows don't contribute linear/area.
        assert!(door_row.linear_m.abs() < 1e-9);
    }

    #[test]
    fn csv_has_header_and_one_row_per_kind() {
        let s = Schedule::from_document(&sample_doc());
        let csv = s.to_csv();
        let lines: Vec<_> = csv.lines().collect();
        // header + 2 kinds (Wall, Door).
        assert_eq!(lines.len(), 3);
        assert!(lines[0].starts_with("kind,"));
        assert!(lines.iter().any(|l| l.starts_with("Wall,")));
        assert!(lines.iter().any(|l| l.starts_with("Door,")));
    }

    #[test]
    fn text_table_includes_labels() {
        let s = Schedule::from_document(&sample_doc());
        let t = s.to_text_table();
        assert!(t.contains("Wall"));
        assert!(t.contains("Door"));
        assert!(t.contains("count"));
    }
}
