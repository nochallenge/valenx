//! The [`Assembly`] container — owns parts, mates, and joints.

use crate::error::AssemblyError;
use crate::joint::Joint;
use crate::mate::Mate;
use crate::part::Part;

/// Multi-part assembly scene.
///
/// The container is intentionally minimal — collections are public
/// `Vec`s indexed by stable id (assigned at insert time and reused
/// by [`Mate`] / [`Joint`] references). `add_*` returns the id; the
/// caller is responsible for stashing it.
///
/// `delete_part` *also* removes any mate / joint that referenced the
/// deleted part — leaving dangling references would leak into
/// solver lookups and produce [`AssemblyError::UnknownPart`].
#[derive(Clone, Debug, Default)]
pub struct Assembly {
    /// All parts. `parts[i].id` is the stable id; `i` is not the id
    /// (parts get inserted in id order so they happen to align, but
    /// the solver looks up by id, not index — `delete_part` leaves a
    /// hole rather than re-indexing).
    pub parts: Vec<Part>,
    /// All mates. See [`Mate`].
    pub mates: Vec<Mate>,
    /// All joints. See [`Joint`].
    pub joints: Vec<Joint>,

    /// Next part id to hand out — monotonically increasing.
    next_part_id: usize,
    /// Next mate id.
    next_mate_id: usize,
    /// Next joint id.
    next_joint_id: usize,
}

impl Assembly {
    /// Build an empty assembly.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a part. Ignores the caller's `part.id` field and replaces
    /// it with a freshly allocated stable id. Returns the id.
    pub fn add_part(&mut self, mut part: Part) -> usize {
        let id = self.next_part_id;
        self.next_part_id += 1;
        part.id = id;
        self.parts.push(part);
        id
    }

    /// Add a mate. Same id-rewrite semantics as [`Self::add_part`].
    pub fn add_mate(&mut self, mut mate: Mate) -> usize {
        let id = self.next_mate_id;
        self.next_mate_id += 1;
        mate.id = id;
        self.mates.push(mate);
        id
    }

    /// Add a joint. Same id-rewrite semantics as [`Self::add_part`].
    pub fn add_joint(&mut self, mut joint: Joint) -> usize {
        let id = self.next_joint_id;
        self.next_joint_id += 1;
        joint.id = id;
        self.joints.push(joint);
        id
    }

    /// Delete a part by id. Also removes any mate / joint that
    /// referenced it (leaving them in would surface as
    /// [`AssemblyError::UnknownPart`] inside the solver).
    pub fn delete_part(&mut self, id: usize) -> Result<(), AssemblyError> {
        let pos = self
            .parts
            .iter()
            .position(|p| p.id == id)
            .ok_or(AssemblyError::UnknownPart(id))?;
        self.parts.remove(pos);
        self.mates.retain(|m| {
            let (a, b) = m.kind.parts();
            a != id && b != id
        });
        self.joints.retain(|j| {
            let (a, b) = j.kind.parts();
            a != id && b != id
        });
        Ok(())
    }

    /// Delete a mate by id.
    pub fn delete_mate(&mut self, id: usize) -> Result<(), AssemblyError> {
        let pos = self
            .mates
            .iter()
            .position(|m| m.id == id)
            .ok_or(AssemblyError::UnknownMate(id))?;
        self.mates.remove(pos);
        Ok(())
    }

    /// Delete a joint by id.
    pub fn delete_joint(&mut self, id: usize) -> Result<(), AssemblyError> {
        let pos = self
            .joints
            .iter()
            .position(|j| j.id == id)
            .ok_or(AssemblyError::UnknownJoint(id))?;
        self.joints.remove(pos);
        Ok(())
    }

    /// Borrow a part by id.
    pub fn get_part(&self, id: usize) -> Result<&Part, AssemblyError> {
        self.parts
            .iter()
            .find(|p| p.id == id)
            .ok_or(AssemblyError::UnknownPart(id))
    }

    /// Borrow a part mutably by id.
    pub fn get_part_mut(&mut self, id: usize) -> Result<&mut Part, AssemblyError> {
        self.parts
            .iter_mut()
            .find(|p| p.id == id)
            .ok_or(AssemblyError::UnknownPart(id))
    }

    /// Borrow a mate by id.
    pub fn get_mate(&self, id: usize) -> Result<&Mate, AssemblyError> {
        self.mates
            .iter()
            .find(|m| m.id == id)
            .ok_or(AssemblyError::UnknownMate(id))
    }

    /// Borrow a mate mutably by id.
    pub fn get_mate_mut(&mut self, id: usize) -> Result<&mut Mate, AssemblyError> {
        self.mates
            .iter_mut()
            .find(|m| m.id == id)
            .ok_or(AssemblyError::UnknownMate(id))
    }

    /// Borrow a joint by id.
    pub fn get_joint(&self, id: usize) -> Result<&Joint, AssemblyError> {
        self.joints
            .iter()
            .find(|j| j.id == id)
            .ok_or(AssemblyError::UnknownJoint(id))
    }

    /// Borrow a joint mutably by id.
    pub fn get_joint_mut(&mut self, id: usize) -> Result<&mut Joint, AssemblyError> {
        self.joints
            .iter_mut()
            .find(|j| j.id == id)
            .ok_or(AssemblyError::UnknownJoint(id))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::joint::JointKind;
    use crate::mate::MateKind;
    use crate::part::PartTransform;
    use nalgebra::Vector3;

    fn unit_cube_part(name: &str) -> Part {
        Part {
            id: 0,
            name: name.into(),
            solid: valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap(),
            transform: PartTransform::identity(),
            fixed: false,
        }
    }

    #[test]
    fn add_part_returns_unique_ids() {
        let mut a = Assembly::new();
        let id0 = a.add_part(unit_cube_part("a"));
        let id1 = a.add_part(unit_cube_part("b"));
        assert_eq!(id0, 0);
        assert_eq!(id1, 1);
        assert_eq!(a.parts.len(), 2);
    }

    #[test]
    fn delete_part_cascades_to_mates_and_joints() {
        let mut a = Assembly::new();
        let id0 = a.add_part(unit_cube_part("a"));
        let id1 = a.add_part(unit_cube_part("b"));
        a.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: id0,
                point_a: Vector3::zeros(),
                part_b: id1,
                point_b: Vector3::zeros(),
            },
        ));
        a.add_joint(Joint::new(
            0,
            JointKind::Revolute {
                part_a: id0,
                part_b: id1,
                axis_origin: Vector3::zeros(),
                axis_dir: Vector3::z(),
            },
        ));
        assert_eq!(a.mates.len(), 1);
        assert_eq!(a.joints.len(), 1);
        a.delete_part(id0).unwrap();
        // Mate + joint both referenced part 0 → both dropped.
        assert!(a.mates.is_empty());
        assert!(a.joints.is_empty());
    }

    #[test]
    fn delete_unknown_part_errors() {
        let mut a = Assembly::new();
        assert!(matches!(
            a.delete_part(42),
            Err(AssemblyError::UnknownPart(42))
        ));
    }

    #[test]
    fn get_part_returns_borrow_by_id() {
        let mut a = Assembly::new();
        let id = a.add_part(unit_cube_part("a"));
        assert_eq!(a.get_part(id).unwrap().name, "a");
        assert!(matches!(
            a.get_part(99),
            Err(AssemblyError::UnknownPart(99))
        ));
    }

    #[test]
    fn add_joint_returns_unique_ids() {
        let mut a = Assembly::new();
        let id0 = a.add_part(unit_cube_part("a"));
        let id1 = a.add_part(unit_cube_part("b"));
        let j0 = a.add_joint(Joint::new(
            0,
            JointKind::Fixed {
                part_a: id0,
                part_b: id1,
            },
        ));
        let j1 = a.add_joint(Joint::new(
            0,
            JointKind::Fixed {
                part_a: id0,
                part_b: id1,
            },
        ));
        assert_eq!(j0, 0);
        assert_eq!(j1, 1);
    }
}
