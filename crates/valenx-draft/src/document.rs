//! Top-level draft document.
//!
//! A [`DraftDocument`] owns a [`WorkingPlane`] and a flat list of
//! [`DraftEntity`] values. Insertion order is preserved and used as
//! the entity's stable index for selection / deletion. There is no
//! parametric tree — each entity carries its own coordinates.

use serde::{Deserialize, Serialize};

use crate::entity::DraftEntity;
use crate::error::DraftError;
use crate::plane::WorkingPlane;

/// A 2D drawing on a [`WorkingPlane`] in 3D space.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DraftDocument {
    /// The plane every entity is expressed in.
    pub working_plane: WorkingPlane,
    /// All entities, in insertion order. Index is the entity's
    /// stable identifier for UI selection and `delete_entity`.
    pub entities: Vec<DraftEntity>,
}

impl DraftDocument {
    /// Create an empty document on the given plane.
    pub fn new(working_plane: WorkingPlane) -> Self {
        Self {
            working_plane,
            entities: Vec::new(),
        }
    }

    /// Append `entity` and return its index in `entities`.
    pub fn add_entity(&mut self, entity: DraftEntity) -> usize {
        let idx = self.entities.len();
        self.entities.push(entity);
        idx
    }

    /// Remove the entity at `idx`. Indices past `idx` shift down by
    /// one — callers that store entity ids long-term must refresh.
    ///
    /// Returns [`DraftError::UnknownEntity`] when `idx` is out of
    /// range so the caller can surface a friendly message instead of
    /// panicking.
    pub fn delete_entity(&mut self, idx: usize) -> Result<DraftEntity, DraftError> {
        if idx >= self.entities.len() {
            return Err(DraftError::UnknownEntity(idx));
        }
        Ok(self.entities.remove(idx))
    }

    /// Number of entities currently in the document.
    pub fn entity_count(&self) -> usize {
        self.entities.len()
    }
}

impl Default for DraftDocument {
    /// A document on the world XY plane with no entities.
    fn default() -> Self {
        Self::new(WorkingPlane::from_xy())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_document_has_zero_entities() {
        let d = DraftDocument::new(WorkingPlane::from_xy());
        assert_eq!(d.entity_count(), 0);
    }

    #[test]
    fn add_then_delete_round_trips_entity_count() {
        let mut d = DraftDocument::new(WorkingPlane::from_xy());
        let idx_a = d.add_entity(DraftEntity::Line {
            start: [0.0, 0.0],
            end: [1.0, 0.0],
        });
        let idx_b = d.add_entity(DraftEntity::Circle {
            center: [0.0, 0.0],
            radius: 1.0,
        });
        assert_eq!(idx_a, 0);
        assert_eq!(idx_b, 1);
        assert_eq!(d.entity_count(), 2);
        let removed = d.delete_entity(0).unwrap();
        assert_eq!(
            removed,
            DraftEntity::Line {
                start: [0.0, 0.0],
                end: [1.0, 0.0]
            }
        );
        assert_eq!(d.entity_count(), 1);
    }

    #[test]
    fn delete_out_of_range_returns_unknown_entity() {
        let mut d = DraftDocument::new(WorkingPlane::from_xy());
        match d.delete_entity(7).unwrap_err() {
            DraftError::UnknownEntity(7) => {}
            other => panic!("wrong variant: {other:?}"),
        }
    }
}
