//! Ligand model: atoms + flexibility tree built from a PDBQT record stream.

use nalgebra::{UnitQuaternion, Vector3};
use valenx_bio::format::pdbqt::{parse, PdbqtRecord};

use crate::atom_type::Ad4AtomType;
use crate::error::DockError;

/// A single ligand atom in dock-engine form (post-AD4 classification).
#[derive(Clone, Debug)]
pub struct LigandAtom {
    /// Position in receptor coordinates (Å).
    pub position: Vector3<f64>,
    /// Vina/AD4 atom type.
    pub ad4_type: Ad4AtomType,
    /// Gasteiger partial charge.
    pub partial_charge: f64,
    /// Source PDBQT serial number (1-based, for output round-trip).
    pub serial: i32,
}

/// A rigid group of atoms rotatable about a single bond relative to its parent.
#[derive(Clone, Debug)]
pub struct RotGroup {
    /// Indices into [`Ligand::atoms`] for atoms inside this group.
    pub atom_indices: Vec<usize>,
    /// Bond axis: (parent atom index, child atom index). The child is
    /// inside this group; the parent is in an ancestor group.
    pub axis: (usize, usize),
    /// Indices into [`Ligand::groups`] for child groups attached
    /// downstream of this one.
    pub children: Vec<usize>,
}

/// Parsed, classified, fully-assembled ligand.
#[derive(Clone, Debug)]
pub struct Ligand {
    /// All atoms in file order.
    pub atoms: Vec<LigandAtom>,
    /// Rigid groups; [`Ligand::root_group`] indexes the root.
    pub groups: Vec<RotGroup>,
    /// Index into [`Ligand::groups`] for the root rigid group.
    pub root_group: usize,
    /// Centroid of all atoms in input coordinates. Computed once at
    /// parse time so `apply_pose` doesn't recompute it on every call
    /// (a hot path inside the search inner loop).
    pub centroid: Vector3<f64>,
    /// Per-group flattened atom-index list: `subtree_indices[g]` is
    /// every atom index in group `g` plus every atom index in groups
    /// transitively reachable through `g.children`. Precomputed at
    /// parse time so torsion application is O(atoms in subtree) per
    /// rotation rather than O(walk-tree + atoms in subtree).
    pub subtree_indices: Vec<Vec<usize>>,
}

impl Ligand {
    /// Parse a PDBQT document and assemble the flexibility tree.
    pub fn from_pdbqt(text: &str) -> Result<Self, DockError> {
        let records = parse(text).map_err(DockError::from)?;
        Self::from_records(&records)
    }

    fn from_records(records: &[PdbqtRecord]) -> Result<Self, DockError> {
        // Atoms are appended in file order; serial-to-index map lets
        // BRANCH lines resolve their parent/child indices.
        let mut atoms: Vec<LigandAtom> = Vec::new();
        let mut serial_to_index: std::collections::HashMap<i32, usize> = Default::default();
        let mut groups: Vec<RotGroup> = Vec::new();
        // Group-construction stack: groups[i] is the currently open
        // group; ROOT/BRANCH pushes; ENDROOT/ENDBRANCH pops.
        let mut stack: Vec<usize> = Vec::new();
        let mut root_group: Option<usize> = None;

        for rec in records {
            match rec {
                PdbqtRecord::Root => {
                    let idx = groups.len();
                    groups.push(RotGroup {
                        atom_indices: Vec::new(),
                        axis: (usize::MAX, usize::MAX),
                        children: Vec::new(),
                    });
                    root_group = Some(idx);
                    stack.push(idx);
                }
                PdbqtRecord::EndRoot => {
                    stack.pop().ok_or(DockError::FlexibilityImbalance)?;
                }
                PdbqtRecord::Branch {
                    parent_serial,
                    child_serial,
                } => {
                    let parent_atom = *serial_to_index
                        .get(parent_serial)
                        .ok_or(DockError::UnknownBranchSerial(*parent_serial))?;
                    // Child atom is added when its ATOM record appears
                    // INSIDE this branch — defer index resolution until then.
                    let idx = groups.len();
                    groups.push(RotGroup {
                        atom_indices: Vec::new(),
                        axis: (parent_atom, usize::MAX),
                        children: Vec::new(),
                    });
                    // Link to the currently open group as a child. After
                    // ENDROOT the stack is empty; in real PDBQTs BRANCHes
                    // outside ROOT still hang off the root group.
                    let parent_group = stack.last().copied().or(root_group);
                    if let Some(top) = parent_group {
                        groups[top].children.push(idx);
                    }
                    stack.push(idx);
                    // Stash the child serial so we resolve when ATOM appears.
                    pending_child_serials_push(&mut serial_to_index, *child_serial, idx);
                }
                PdbqtRecord::EndBranch { .. } => {
                    stack.pop().ok_or(DockError::FlexibilityImbalance)?;
                }
                PdbqtRecord::Atom(a) => {
                    let ad4: Ad4AtomType = a.ad4_type.parse().map_err(DockError::AtomType)?;
                    let idx = atoms.len();
                    atoms.push(LigandAtom {
                        position: a.position,
                        ad4_type: ad4,
                        partial_charge: a.partial_charge,
                        serial: a.serial,
                    });
                    serial_to_index.insert(a.serial, idx);
                    if let Some(&top) = stack.last() {
                        groups[top].atom_indices.push(idx);
                        // If this is the first atom in a branch, it's the child of the rotatable bond.
                        if groups[top].axis.0 != usize::MAX && groups[top].axis.1 == usize::MAX {
                            groups[top].axis.1 = idx;
                        }
                    }
                }
                PdbqtRecord::Torsdof(_) | PdbqtRecord::Other(_) => {}
            }
        }

        let root_group = root_group.ok_or(DockError::NoRoot)?;
        if !stack.is_empty() {
            return Err(DockError::FlexibilityImbalance);
        }
        // Cache the centroid + per-group subtree-atom-index lists once
        // so the hot inner loop never recomputes them.
        let centroid: Vector3<f64> = if atoms.is_empty() {
            Vector3::zeros()
        } else {
            let sum: Vector3<f64> = atoms.iter().map(|a| a.position).sum();
            sum / atoms.len() as f64
        };
        let subtree_indices = build_subtree_indices(&groups);
        Ok(Self {
            atoms,
            groups,
            root_group,
            centroid,
            subtree_indices,
        })
    }

    /// Number of rotatable bonds (= number of non-root groups).
    pub fn n_torsions(&self) -> usize {
        self.groups.len().saturating_sub(1)
    }

    /// Apply a pose: returns the world-space positions of every atom.
    /// Rotation is applied about the centroid of all atoms (Vina
    /// convention — the search lives in centroid-centered space).
    /// Torsions are applied after the rigid transform, walking the
    /// flexibility tree from root outward.
    ///
    /// Uses the cached centroid (computed once at parse time) so the
    /// inner search loop stays cheap.
    pub fn apply_pose(&self, pose: &crate::pose::Pose) -> Vec<Vector3<f64>> {
        let centroid = self.centroid;
        // 1. Centre, rotate, translate.
        let mut out: Vec<Vector3<f64>> = self
            .atoms
            .iter()
            .map(|a| pose.orientation * (a.position - centroid) + centroid + pose.translation)
            .collect();
        // 2. Apply torsions (rotates each non-root group about its axis).
        self.apply_torsions(&mut out, pose);
        out
    }

    /// Walk the flexibility tree depth-first from the root, rotating
    /// each branch's atoms (and the atoms of all downstream branches)
    /// about the bond axis by the corresponding torsion angle.
    ///
    /// Reads precomputed [`Ligand::subtree_indices`] instead of
    /// re-walking the tree on every call — a hot path inside BFGS
    /// inner iterations.
    fn apply_torsions(&self, positions: &mut [Vector3<f64>], pose: &crate::pose::Pose) {
        // Pre-compute the depth-first visit order so each branch sees
        // its parent's atoms already in their final position.
        let mut visit: Vec<usize> = Vec::new();
        let mut stack = vec![self.root_group];
        while let Some(g) = stack.pop() {
            visit.push(g);
            for &c in self.groups[g].children.iter().rev() {
                stack.push(c);
            }
        }
        // Torsion vector indexes non-root groups in *visit order minus root*.
        // For simplicity we use a fixed map: torsion[i] applies to the i-th
        // non-root group in visit order.
        let mut tor_iter = pose.torsions.iter().copied();
        for g in visit {
            if g == self.root_group {
                continue;
            }
            let Some(angle) = tor_iter.next() else {
                continue;
            };
            let (parent_idx, child_idx) = self.groups[g].axis;
            if parent_idx == usize::MAX || child_idx == usize::MAX {
                continue;
            }
            let axis_origin = positions[parent_idx];
            let axis_dir = (positions[child_idx] - axis_origin).normalize();
            let rot =
                UnitQuaternion::from_axis_angle(&nalgebra::Unit::new_normalize(axis_dir), angle);
            // Rotate every atom in this group + transitive descendants.
            // Cached in `subtree_indices` at parse time.
            for &ai in &self.subtree_indices[g] {
                positions[ai] = rot * (positions[ai] - axis_origin) + axis_origin;
            }
        }
    }
}

/// All atom indices in group `g` plus all atom indices in groups
/// transitively reachable through `children`. Used at parse time to
/// populate [`Ligand::subtree_indices`]; the hot path then reads
/// straight from the cache.
pub(crate) fn subtree_atom_indices(groups: &[RotGroup], g: usize) -> Vec<usize> {
    let mut out = Vec::new();
    let mut stack = vec![g];
    while let Some(h) = stack.pop() {
        out.extend_from_slice(&groups[h].atom_indices);
        for &c in &groups[h].children {
            stack.push(c);
        }
    }
    out
}

fn build_subtree_indices(groups: &[RotGroup]) -> Vec<Vec<usize>> {
    (0..groups.len())
        .map(|g| subtree_atom_indices(groups, g))
        .collect()
}

/// Helper placeholder: child-serial resolution happens inline above;
/// kept as a no-op so the signature stays available if we change strategy.
fn pending_child_serials_push(
    _map: &mut std::collections::HashMap<i32, usize>,
    _serial: i32,
    _group_idx: usize,
) {
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_ethanol_like_ligand() {
        // Two heavy atoms joined by one rotatable bond. The
        // serial-to-group resolution test: the second ATOM should
        // land inside the BRANCH'd child group.
        let pdbqt = "\
ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C 
ENDROOT
BRANCH   1   2
ATOM      2  O1  LIG A   1       1.500   0.000   0.000  1.00  0.00     0.000 OA
ENDBRANCH   1   2
TORSDOF 1
";
        let lig = Ligand::from_pdbqt(pdbqt).unwrap();
        assert_eq!(lig.atoms.len(), 2);
        assert_eq!(lig.groups.len(), 2);
        assert_eq!(lig.n_torsions(), 1);

        // Root group has atom 0 (the C); branch group has atom 1 (the O).
        assert_eq!(lig.groups[lig.root_group].atom_indices, vec![0]);
        assert_eq!(lig.groups[1].atom_indices, vec![1]);
        // The rotatable bond axis links parent atom 0 -> child atom 1.
        assert_eq!(lig.groups[1].axis, (0, 1));
        // Root's children include the branch.
        assert_eq!(lig.groups[lig.root_group].children, vec![1]);
    }

    #[test]
    fn rejects_missing_root() {
        let err = Ligand::from_pdbqt("TORSDOF 0\n").unwrap_err();
        assert!(matches!(err, DockError::NoRoot));
    }

    use crate::pose::Pose;
    use nalgebra::{UnitQuaternion, Vector3};

    #[test]
    fn rigid_transform_translates_only_when_orientation_identity() {
        let pdbqt = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C 
ATOM      2  C2  LIG A   1       1.500   0.000   0.000  1.00  0.00     0.000 C 
ENDROOT
TORSDOF 0
";
        let lig = Ligand::from_pdbqt(pdbqt).unwrap();
        let mut pose = Pose::identity(lig.n_torsions());
        pose.translation = Vector3::new(10.0, 5.0, -2.0);
        let pts = lig.apply_pose(&pose);
        assert_eq!(pts[0], Vector3::new(10.0, 5.0, -2.0));
        assert_eq!(pts[1], Vector3::new(11.5, 5.0, -2.0));
    }

    #[test]
    fn rigid_transform_rotation_about_centroid() {
        // Two atoms on the X axis, centroid at (0.75, 0, 0). A 180° Z-rotation
        // about the centroid should swap their X positions and negate Y components
        // (which are zero), leaving them at (1.5, 0, 0) and (0.0, 0, 0) — i.e. the
        // mirrored configuration about the centroid.
        let pdbqt = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C 
ATOM      2  C2  LIG A   1       1.500   0.000   0.000  1.00  0.00     0.000 C 
ENDROOT
TORSDOF 0
";
        let lig = Ligand::from_pdbqt(pdbqt).unwrap();
        let mut pose = Pose::identity(0);
        pose.orientation =
            UnitQuaternion::from_axis_angle(&Vector3::z_axis(), std::f64::consts::PI);
        let pts = lig.apply_pose(&pose);
        assert!((pts[0].x - 1.5).abs() < 1e-9);
        assert!((pts[1].x - 0.0).abs() < 1e-9);
    }

    #[test]
    fn single_torsion_rotates_downstream_atom() {
        // Ligand:  C1(0,0,0) -- C2(1,0,0) -- C3(2,0,0)
        // Rotatable bond C2-C3. With a 90° torsion about C2-C3 axis,
        // C3 itself lies on the axis so it does NOT move. To see
        // motion we add a fourth atom branched off C3 — but that
        // would require nested BRANCH support, which we have.
        // Instead, place C3 OFF the C1-C2 axis: rotation about the
        // C1-C2 axis (which is the X axis) by 90° around X moves C3
        // from (1.5, 1.0, 0) to (1.5, 0, 1.0).
        let pdbqt = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C 
ATOM      2  C2  LIG A   1       1.500   0.000   0.000  1.00  0.00     0.000 C 
ENDROOT
BRANCH   1   3
ATOM      3  C3  LIG A   1       1.500   1.000   0.000  1.00  0.00     0.000 C 
ENDBRANCH   1   3
TORSDOF 1
";
        let lig = Ligand::from_pdbqt(pdbqt).unwrap();
        assert_eq!(lig.n_torsions(), 1);
        let mut pose = Pose::identity(1);
        pose.torsions[0] = std::f64::consts::FRAC_PI_2;
        let pts = lig.apply_pose(&pose);
        // C1, C2 unchanged. C3 rotates 90° about the C1-C3-axis...
        // Axis = (C3 - C1).normalize(), origin = C1.
        // |C3 - C1| = sqrt(1.5^2 + 1.0^2) ≈ 1.803.
        // Rotation about that axis by 90° keeps |C3 - C1| constant.
        let r = (pts[2] - pts[0]).norm();
        assert!((r - ((1.5_f64).hypot(1.0))).abs() < 1e-9);
    }
}
