//! RON-based persistence for assemblies.
//!
//! Stub for Phase 6A. Phase 6D wires the full envelope, with the
//! caveat that [`crate::Part::solid`] does **not** round-trip (truck
//! BRep isn't serializable) — callers re-attach geometry after
//! loading.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::assembly::Assembly;
use crate::error::AssemblyError;
use crate::joint::Joint;
use crate::mate::Mate;
use crate::part::PartTransform;

/// On-disk envelope wrapping an assembly with format-version metadata.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AssemblyFile {
    /// Format version — bumped when on-disk schema changes.
    pub version: u32,
    /// Persistent payload — parts (without their solid geometry),
    /// mates, joints.
    pub parts: Vec<PartShell>,
    /// All mates.
    pub mates: Vec<Mate>,
    /// All joints.
    pub joints: Vec<Joint>,
}

/// Persistent shape of a [`crate::Part`] — the BRep solid is *not*
/// serialized (truck doesn't support it). Callers re-attach geometry
/// after [`AssemblyFile::into_assembly`] using the saved name as a
/// lookup key.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PartShell {
    /// Stable id.
    pub id: usize,
    /// Display name — used as the lookup key when re-attaching solids.
    pub name: String,
    /// Current world-space pose.
    pub transform: PartTransform,
    /// "Don't move me" flag.
    pub fixed: bool,
}

impl AssemblyFile {
    /// Current on-disk format version.
    pub const VERSION: u32 = 1;

    /// Wrap an assembly. The `solid` field of each [`crate::Part`] is
    /// dropped — see the module doc for the rationale.
    pub fn from_assembly(a: &Assembly) -> Self {
        Self {
            version: Self::VERSION,
            parts: a
                .parts
                .iter()
                .map(|p| PartShell {
                    id: p.id,
                    name: p.name.clone(),
                    transform: p.transform.clone(),
                    fixed: p.fixed,
                })
                .collect(),
            mates: a.mates.clone(),
            joints: a.joints.clone(),
        }
    }

    /// Serialize to a RON string.
    pub fn to_ron(&self) -> Result<String, AssemblyError> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
            .map_err(|e| AssemblyError::Ron(e.to_string()))
    }

    /// Write to a file. Round-28 H2: routed through the canonical
    /// `valenx_core::io_caps::atomic_write_str` (sidecar O_NOFOLLOW,
    /// fsync-before-rename, parent-dir fsync on Unix). Pre-fix this
    /// was a bare `std::fs::write` which silently followed leaf
    /// symlinks and was non-atomic.
    pub fn write_to(&self, path: &Path) -> Result<(), AssemblyError> {
        let ron = self.to_ron()?;
        valenx_core::io_caps::atomic_write_str(path, &ron)?;
        Ok(())
    }

    /// Parse from a RON string.
    pub fn from_ron(s: &str) -> Result<Self, AssemblyError> {
        ron::from_str(s).map_err(|e| AssemblyError::Ron(e.to_string()))
    }

    /// Read from a file.
    ///
    /// R29 D: bounded at [`valenx_core::io_caps::MAX_DOC_FILE_BYTES`]
    /// (16 MiB) via the canonical helper, replacing the per-crate
    /// private `read_capped_to_string` duplicate.
    pub fn read_from(path: &Path) -> Result<Self, AssemblyError> {
        let s = valenx_core::io_caps::read_capped_to_string(
            path,
            valenx_core::io_caps::MAX_DOC_FILE_BYTES,
        )?;
        Self::from_ron(&s)
    }

    /// Hydrate a partial [`Assembly`] — mates + joints + part shells
    /// (without solids). The caller supplies a `solid_for` closure
    /// that maps a part name → solid; if the closure returns `None`,
    /// the part is dropped on the floor (silently, so that loading a
    /// project where some referenced parts are missing surfaces as
    /// "missing parts in scene tree" rather than a hard error).
    pub fn into_assembly<F>(self, mut solid_for: F) -> Assembly
    where
        F: FnMut(&str) -> Option<valenx_cad::Solid>,
    {
        let mut a = Assembly::new();
        // Re-allocate ids inside `add_part` — but mates / joints reference
        // the *old* ids. We use a name→new-id map to remap.
        let mut id_remap: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();
        for shell in self.parts {
            let Some(solid) = solid_for(&shell.name) else {
                continue;
            };
            let part = crate::part::Part {
                id: 0, // overwritten by add_part
                name: shell.name,
                solid,
                transform: shell.transform,
                fixed: shell.fixed,
            };
            let new_id = a.add_part(part);
            id_remap.insert(shell.id, new_id);
        }
        for mut m in self.mates {
            // Remap part ids inside the mate kind.
            if let Some((new_a, new_b)) = remap_pair(m.kind.parts(), &id_remap) {
                rewrite_mate_parts(&mut m, new_a, new_b);
                a.add_mate(m);
            }
        }
        for mut j in self.joints {
            if let Some((new_a, new_b)) = remap_pair(j.kind.parts(), &id_remap) {
                rewrite_joint_parts(&mut j, new_a, new_b);
                a.add_joint(j);
            }
        }
        a
    }
}

fn remap_pair(
    (a, b): (usize, usize),
    map: &std::collections::HashMap<usize, usize>,
) -> Option<(usize, usize)> {
    Some((*map.get(&a)?, *map.get(&b)?))
}

fn rewrite_mate_parts(m: &mut Mate, new_a: usize, new_b: usize) {
    use crate::mate::MateKind;
    match &mut m.kind {
        MateKind::Coincident { part_a, part_b, .. }
        | MateKind::Distance { part_a, part_b, .. }
        | MateKind::Angle { part_a, part_b, .. }
        | MateKind::Parallel { part_a, part_b, .. }
        | MateKind::Perpendicular { part_a, part_b, .. }
        | MateKind::Tangent { part_a, part_b, .. } => {
            *part_a = new_a;
            *part_b = new_b;
        }
    }
}

fn rewrite_joint_parts(j: &mut Joint, new_a: usize, new_b: usize) {
    use crate::joint::JointKind;
    match &mut j.kind {
        JointKind::Fixed { part_a, part_b }
        | JointKind::Revolute { part_a, part_b, .. }
        | JointKind::Prismatic { part_a, part_b, .. }
        | JointKind::Cylindrical { part_a, part_b, .. }
        | JointKind::Spherical { part_a, part_b, .. }
        | JointKind::Planar { part_a, part_b, .. } => {
            *part_a = new_a;
            *part_b = new_b;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::joint::JointKind;
    use crate::mate::MateKind;
    use crate::part::Part;
    use nalgebra::Vector3;

    fn unit_cube(name: &str) -> Part {
        Part::new(0, name, valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap())
    }

    #[test]
    fn empty_assembly_round_trips() {
        let a = Assembly::new();
        let ron = AssemblyFile::from_assembly(&a).to_ron().unwrap();
        assert!(ron.contains("version: 1"));
        let parsed = AssemblyFile::from_ron(&ron).unwrap();
        assert_eq!(parsed.version, 1);
        assert!(parsed.parts.is_empty());
        assert!(parsed.mates.is_empty());
        assert!(parsed.joints.is_empty());
    }

    /// Task 28 — round-trip an assembly with parts, mates, and joints.
    /// Solids are dropped (truck BRep isn't serializable) and re-attached
    /// via the `solid_for` closure during load.
    #[test]
    fn round_trip_with_parts_mates_joints() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("b"));
        a.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: id_a,
                point_a: Vector3::new(0.0, 0.0, 0.0),
                part_b: id_b,
                point_b: Vector3::new(0.5, 0.5, 0.5),
            },
        ));
        a.add_joint(Joint::new(
            0,
            JointKind::Revolute {
                part_a: id_a,
                part_b: id_b,
                axis_origin: Vector3::new(1.0, 0.0, 0.0),
                axis_dir: Vector3::z(),
            },
        ));

        let ron = AssemblyFile::from_assembly(&a).to_ron().unwrap();
        let parsed = AssemblyFile::from_ron(&ron).unwrap();
        assert_eq!(parsed.parts.len(), 2);
        assert_eq!(parsed.mates.len(), 1);
        assert_eq!(parsed.joints.len(), 1);
        assert!(parsed.parts[0].fixed);
        assert_eq!(parsed.parts[0].name, "a");
        assert_eq!(parsed.parts[1].name, "b");

        // Re-hydrate into an Assembly using a solid-lookup closure.
        let rebuilt = parsed.into_assembly(|name| match name {
            "a" | "b" => Some(valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap()),
            _ => None,
        });
        assert_eq!(rebuilt.parts.len(), 2);
        assert_eq!(rebuilt.mates.len(), 1);
        assert_eq!(rebuilt.joints.len(), 1);
        // Confirm mate's part references were remapped to the new ids.
        let (ma, mb) = rebuilt.mates[0].kind.parts();
        assert_eq!(ma, rebuilt.parts[0].id);
        assert_eq!(mb, rebuilt.parts[1].id);
    }

    #[test]
    fn write_and_read_round_trip_via_file() {
        let mut a = Assembly::new();
        a.add_part(unit_cube("solo"));
        let f = AssemblyFile::from_assembly(&a);
        let tmp = std::env::temp_dir().join("valenx_assembly_test.ron");
        f.write_to(&tmp).unwrap();
        let read = AssemblyFile::read_from(&tmp).unwrap();
        assert_eq!(read.parts.len(), 1);
        assert_eq!(read.parts[0].name, "solo");
        let _ = std::fs::remove_file(&tmp);
    }

    /// Round-12 M1: a file larger than MAX_DOC_FILE_BYTES must be
    /// rejected at the read-cap layer rather than slurped into
    /// memory. Assembly is the biggest workbench document type, so
    /// this test exercises the most demanding case in the persist.rs
    /// family swept in round-12.
    #[test]
    fn read_from_rejects_17mib_file() {
        let tmp = std::env::temp_dir().join("valenx_assembly_oversize_test.ron");
        // 17 MiB of zeros — over the 16 MiB cap.
        let oversize = valenx_core::io_caps::MAX_DOC_FILE_BYTES + 1024 * 1024;
        std::fs::write(&tmp, vec![0u8; oversize]).unwrap();
        let err = AssemblyFile::read_from(&tmp).expect_err("must reject oversize file");
        match err {
            AssemblyError::Io(io_err) => {
                assert_eq!(io_err.kind(), std::io::ErrorKind::InvalidData);
            }
            other => panic!("expected Io(InvalidData), got: {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn into_assembly_drops_parts_without_solids() {
        let mut a = Assembly::new();
        a.add_part(unit_cube("present"));
        a.add_part(unit_cube("missing"));
        let f = AssemblyFile::from_assembly(&a);
        let rebuilt = f.into_assembly(|name| {
            if name == "present" {
                Some(valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap())
            } else {
                None
            }
        });
        assert_eq!(rebuilt.parts.len(), 1);
        assert_eq!(rebuilt.parts[0].name, "present");
    }
}
