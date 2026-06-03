//! Top-level [`ArchDocument`] — owns the entity table and exposes
//! `add` / `delete` / `get` / `count`.
//!
//! Stable per-entity ids are issued by a monotonic counter
//! (`next_id`); deleting an entity does NOT renumber the remaining
//! entries, so callers can keep an id long-term without re-fetching
//! it from a sliding index. Insertion order is preserved.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::entity::ArchEntity;
use crate::error::ArchError;

/// An architectural project — list of `(id, entity)` pairs plus
/// `next_id` counter and a project name.
///
/// Add / delete / lookup are O(n) (linear scan) because the entity
/// count is small in practice (a residential floor is ~30 walls / 10
/// doors / 30 windows / 4 slabs, well under any threshold where a
/// HashMap would beat a Vec). Trade-off favours predictable iteration
/// order for the UI list.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ArchDocument {
    /// Project / building name used by IFC's `IfcProject.Name`.
    pub project_name: String,
    /// All entities in insertion order with their stable id.
    pub entities: Vec<(usize, ArchEntity)>,
    /// Next id to issue. Always > the max id in `entities`.
    pub next_id: usize,
}

impl ArchDocument {
    /// Create an empty document with the given project name.
    pub fn new(project_name: impl Into<String>) -> Self {
        Self {
            project_name: project_name.into(),
            entities: Vec::new(),
            next_id: 1,
        }
    }

    /// Append `entity` and return its stable id.
    pub fn add_entity(&mut self, entity: ArchEntity) -> usize {
        let id = self.next_id;
        self.next_id += 1;
        self.entities.push((id, entity));
        id
    }

    /// Remove the entity with id `id`. Returns the removed entity on
    /// success, or [`ArchError::UnknownEntity`] when no entity carries
    /// that id.
    pub fn delete_entity(&mut self, id: usize) -> Result<ArchEntity, ArchError> {
        let pos = self
            .entities
            .iter()
            .position(|(e_id, _)| *e_id == id)
            .ok_or(ArchError::UnknownEntity(id))?;
        let (_, ent) = self.entities.remove(pos);
        Ok(ent)
    }

    /// Borrow the entity with id `id`. Returns `None` when no entity
    /// carries that id (lets the UI render a "not found" placeholder
    /// without a Result panic).
    pub fn get_entity(&self, id: usize) -> Option<&ArchEntity> {
        self.entities
            .iter()
            .find_map(|(e_id, ent)| (*e_id == id).then_some(ent))
    }

    /// Total entity count.
    pub fn count(&self) -> usize {
        self.entities.len()
    }

    /// Iterate `(id, &entity)` pairs in insertion order.
    pub fn iter(&self) -> impl Iterator<Item = (usize, &ArchEntity)> + '_ {
        self.entities.iter().map(|(id, e)| (*id, e))
    }

    /// Axis-aligned bounding box across every entity's geometry hints.
    ///
    /// Returns `None` when the document is empty. The bbox is
    /// computed from the parametric inputs (no tessellation): walls
    /// expand by half-thickness in their perpendicular direction,
    /// slabs use the boundary footprint and thickness, etc. This is
    /// fast (no truck calls) and good enough to size the viewport.
    pub fn bbox(&self) -> Option<(Vector3<f64>, Vector3<f64>)> {
        if self.entities.is_empty() {
            return None;
        }
        let mut min = Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut max = Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        let mut accept = |p: Vector3<f64>| {
            for i in 0..3 {
                if p[i] < min[i] {
                    min[i] = p[i];
                }
                if p[i] > max[i] {
                    max[i] = p[i];
                }
            }
        };
        for (_, ent) in &self.entities {
            for p in ent.bbox_hint_points() {
                accept(p);
            }
        }
        Some((min, max))
    }

    /// Fuse every entity's tessellated triangles into a single
    /// viewport-ready [`valenx_mesh::Mesh`]. Empty document → empty
    /// mesh (zero nodes, zero elements).
    ///
    /// Each entity is tessellated with `tolerance` chord-error in
    /// model units. The fused mesh has one `Tri3` element block and a
    /// concatenated node array with per-entity offsets applied to the
    /// triangle indices (same pattern as Phase 6
    /// `render_assembly_to_viewport`).
    pub fn tessellate_all(&self, tolerance: f64) -> Result<valenx_mesh::Mesh, ArchError> {
        use valenx_mesh::{ElementBlock, ElementType, Mesh};

        let mut fused = Mesh::new("arch");
        let mut block = ElementBlock::new(ElementType::Tri3);
        let mut node_offset: u32 = 0;

        // Build a quick wall lookup so opening tessellation can find
        // its host wall for placement.
        let walls: std::collections::HashMap<usize, &crate::wall::WallParams> = self
            .entities
            .iter()
            .filter_map(|(id, e)| match e {
                ArchEntity::Wall(w) => Some((*id, w)),
                _ => None,
            })
            .collect();

        // Pre-bucket openings (windows + doors) by host wall id so
        // we can call `wall_with_openings` once per wall.
        let mut wall_windows: std::collections::HashMap<usize, Vec<&crate::window::WindowParams>> =
            std::collections::HashMap::new();
        let mut wall_doors: std::collections::HashMap<usize, Vec<&crate::door::DoorParams>> =
            std::collections::HashMap::new();
        for (_id, ent) in &self.entities {
            match ent {
                ArchEntity::Window(w) => {
                    wall_windows.entry(w.host).or_default().push(w);
                }
                ArchEntity::Door(d) => {
                    wall_doors.entry(d.host).or_default().push(d);
                }
                _ => {}
            }
        }

        for (id, ent) in &self.entities {
            let mesh = match ent {
                ArchEntity::Wall(w) => {
                    let windows = wall_windows.get(id).map(|v| v.as_slice()).unwrap_or(&[]);
                    let doors = wall_doors.get(id).map(|v| v.as_slice()).unwrap_or(&[]);
                    crate::opening::wall_with_openings(w, windows, doors)
                        .map_err(|e| ArchError::Cad(format!("wall {id}: {e}")))?
                }
                _ => ent
                    .tessellate_in_doc(tolerance, &walls)
                    .map_err(|e| ArchError::Cad(format!("entity {id}: {e}")))?,
            };
            for n in &mesh.nodes {
                fused.nodes.push(*n);
            }
            for tri_block in &mesh.element_blocks {
                if tri_block.element_type != valenx_mesh::ElementType::Tri3 {
                    continue;
                }
                for &i in &tri_block.connectivity {
                    block.connectivity.push(i + node_offset);
                }
            }
            node_offset = fused.nodes.len() as u32;
        }
        fused.element_blocks.push(block);
        fused.recompute_stats();
        Ok(fused)
    }
}

impl Default for ArchDocument {
    /// Empty document named `"Untitled"`.
    fn default() -> Self {
        Self::new("Untitled")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::wall::WallParams;

    fn sample_wall() -> WallParams {
        WallParams {
            start: Vector3::new(0.0, 0.0, 0.0),
            end: Vector3::new(2.0, 0.0, 0.0),
            height: 2.5,
            thickness: 0.2,
            material: "Concrete".into(),
        }
    }

    #[test]
    fn empty_document_has_zero_entities() {
        let d = ArchDocument::new("Test");
        assert_eq!(d.count(), 0);
        assert_eq!(d.project_name, "Test");
        assert!(d.bbox().is_none());
    }

    #[test]
    fn add_then_delete_uses_stable_id() {
        let mut d = ArchDocument::new("p");
        let a = d.add_entity(ArchEntity::Wall(sample_wall()));
        let b = d.add_entity(ArchEntity::Wall(sample_wall()));
        let c = d.add_entity(ArchEntity::Wall(sample_wall()));
        assert_eq!((a, b, c), (1, 2, 3));
        assert_eq!(d.count(), 3);
        d.delete_entity(b).unwrap();
        assert_eq!(d.count(), 2);
        // id `a` still resolves after deleting `b`.
        assert!(d.get_entity(a).is_some());
        assert!(d.get_entity(b).is_none());
        assert!(d.get_entity(c).is_some());
    }

    #[test]
    fn delete_unknown_id_errors() {
        let mut d = ArchDocument::new("p");
        match d.delete_entity(99).unwrap_err() {
            ArchError::UnknownEntity(99) => {}
            other => panic!("wrong variant: {other:?}"),
        }
    }

    #[test]
    fn bbox_walls_only() {
        let mut d = ArchDocument::new("p");
        d.add_entity(ArchEntity::Wall(sample_wall()));
        let bb = d.bbox().unwrap();
        // Min < max on every axis, and height should be at least 2.5.
        assert!(bb.1.x >= 2.0);
        assert!(bb.1.z >= 2.5);
    }
}
