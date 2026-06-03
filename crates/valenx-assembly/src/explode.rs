//! Exploded-view layout — auto-arrange the parts of an assembly along a
//! direction vector for documentation / animation purposes.
//!
//! ## Algorithm
//!
//! 1. **Topological ordering** — build a graph from `Mate`s (treating
//!    each mate as an undirected edge between the two referenced
//!    parts). The "root" of the explode is the union of all *fixed*
//!    parts; un-fixed parts get a depth equal to their shortest path
//!    (in mates) to any fixed part. Parts with no path stay at depth 0
//!    (they're treated as standalone roots).
//!
//!    Why mates and not joints: in SolidWorks-class exploded views,
//!    the explode hierarchy follows the assembly tree (which mates
//!    define) rather than the kinematic hierarchy (which joints
//!    define). A user who's wired a Revolute joint typically doesn't
//!    want the exploded view to follow it.
//!
//! 2. **Per-part offset** — every part is translated by
//!    `direction.normalized() * depth * spacing` where `spacing` is
//!    configurable. Roots (depth 0, including fixed parts) stay put;
//!    the rest move outward in proportion to their depth.
//!
//! 3. **Animation steps** — [`linear_explode_steps`] returns the
//!    per-part offsets, sorted by depth, so a UI loop can interpolate
//!    `t ∈ [0, 1]` linearly: at `t = 0` parts are at their solved
//!    poses; at `t = 1` they're at the fully-exploded layout.
//!
//! ## Honest scope
//!
//! - **Spacing is uniform.** A "real" exploded view often uses
//!   non-uniform spacing scaled by the parts' size so a tiny screw
//!   doesn't crowd a large casing — that's a follow-up that requires
//!   per-part bounding-box queries; the uniform spacing is honest
//!   and produces useable layouts for assemblies whose parts are of
//!   similar scale.
//! - **Direction is a single vector.** SolidWorks supports an
//!   "explode axis per step" UI; that's also a follow-up.
//! - **Parts not connected to any fixed part get depth 0** (treated
//!   as roots) rather than being implicitly rooted at the centroid
//!   — that matches the "rigid bodies float separately" semantics
//!   of an assembly without any anchors.

use std::collections::VecDeque;

use nalgebra::Vector3;

use crate::assembly::Assembly;
use crate::error::AssemblyError;
use crate::part::PartTransform;

/// Tunable knobs for [`auto_explode`].
#[derive(Clone, Debug)]
pub struct ExplodeConfig {
    /// Distance between adjacent depths along the explode direction.
    /// Default `2.0` — generally larger than the typical part size so
    /// parts visibly separate.
    pub spacing: f64,
}

impl Default for ExplodeConfig {
    fn default() -> Self {
        Self { spacing: 2.0 }
    }
}

/// One step in an animated explode sequence — the source part id and the
/// offset to apply to its solved pose to reach the exploded layout.
#[derive(Clone, Debug, PartialEq)]
pub struct ExplodeStep {
    /// Source part id.
    pub part_id: usize,
    /// Depth in the mate graph (0 = root / fixed, increasing outward).
    pub depth: usize,
    /// World-space offset to apply on top of the solved pose to reach
    /// the exploded position.
    pub offset: Vector3<f64>,
}

/// Output of [`auto_explode`] — the exploded poses keyed by part id,
/// plus the per-part metadata.
#[derive(Clone, Debug)]
pub struct ExplodedAssembly {
    /// Per-part `(part_id, exploded_transform)`. Same ordering as the
    /// source assembly's `parts` field.
    pub poses: Vec<(usize, PartTransform)>,
    /// Per-part metadata sorted by depth → animation step order.
    pub steps: Vec<ExplodeStep>,
}

/// Auto-arrange the parts of an assembly along `direction` into an
/// exploded layout.
///
/// Returns [`AssemblyError::BadParameter`] if `direction` is zero or
/// non-finite.
pub fn auto_explode(
    a: &Assembly,
    direction: Vector3<f64>,
    cfg: ExplodeConfig,
) -> Result<ExplodedAssembly, AssemblyError> {
    let dir_len = direction.norm();
    if !dir_len.is_finite() || dir_len < 1e-12 {
        return Err(AssemblyError::BadParameter {
            name: "direction",
            reason: "must be a non-zero finite vector".into(),
        });
    }
    let unit = direction / dir_len;

    let depths = compute_depths(a);

    let mut poses = Vec::with_capacity(a.parts.len());
    let mut steps = Vec::with_capacity(a.parts.len());
    for part in &a.parts {
        let depth = *depths.get(&part.id).unwrap_or(&0);
        let offset = unit * (depth as f64) * cfg.spacing;
        let exploded = PartTransform {
            translation: part.transform.translation + offset,
            orientation: part.transform.orientation,
        };
        poses.push((part.id, exploded));
        steps.push(ExplodeStep {
            part_id: part.id,
            depth,
            offset,
        });
    }
    steps.sort_by_key(|s| (s.depth, s.part_id));
    Ok(ExplodedAssembly { poses, steps })
}

/// Return only the per-part offsets, in animation step order. The
/// step-0 entries are the roots (zero offset) and entries grow
/// outward by depth.
pub fn linear_explode_steps(
    a: &Assembly,
    direction: Vector3<f64>,
    cfg: ExplodeConfig,
) -> Result<Vec<ExplodeStep>, AssemblyError> {
    Ok(auto_explode(a, direction, cfg)?.steps)
}

/// Breadth-first depth assignment from the fixed parts. Edges come from
/// mates; parts unreachable from any fixed part get depth 0.
fn compute_depths(a: &Assembly) -> std::collections::HashMap<usize, usize> {
    let mut depths = std::collections::HashMap::new();
    let mut q: VecDeque<usize> = VecDeque::new();

    // Roots: every fixed part starts at depth 0.
    let mut any_fixed = false;
    for p in &a.parts {
        if p.fixed {
            depths.insert(p.id, 0usize);
            q.push_back(p.id);
            any_fixed = true;
        }
    }

    // If no fixed parts, treat the first part (in storage order) as the
    // implicit root so the layout still has a meaningful depth gradient
    // for a free-floating assembly. (Without this, every part gets
    // depth 0 and nothing explodes.)
    if !any_fixed {
        if let Some(first) = a.parts.first() {
            depths.insert(first.id, 0usize);
            q.push_back(first.id);
        }
    }

    // Build undirected adjacency from un-suppressed mates.
    let mut adj: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for m in &a.mates {
        if m.suppressed {
            continue;
        }
        let (pa, pb) = m.kind.parts();
        adj.entry(pa).or_default().push(pb);
        adj.entry(pb).or_default().push(pa);
    }

    // BFS outward, assigning the shortest-path depth.
    while let Some(id) = q.pop_front() {
        // Invariant: every id pushed to `q` was also inserted into
        // `depths`; the `else` branch defensively yields 0 rather than
        // panicking, mirroring the no-fixed-root fallback above.
        let d = depths.get(&id).copied().unwrap_or(0);
        if let Some(nbrs) = adj.get(&id) {
            for &n in nbrs {
                // Only assign if not already visited (BFS first-touch =
                // shortest path).
                if let std::collections::hash_map::Entry::Vacant(slot) = depths.entry(n) {
                    slot.insert(d + 1);
                    q.push_back(n);
                }
            }
        }
    }

    // Parts not in `depths` (no path to any root) get depth 0 — they're
    // treated as their own roots.
    for p in &a.parts {
        depths.entry(p.id).or_insert(0);
    }
    depths
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mate::{Mate, MateKind};
    use crate::part::Part;
    use nalgebra::Vector3;

    fn unit_cube(name: &str) -> Part {
        Part::new(0, name, valenx_cad::box_solid(1.0, 1.0, 1.0).unwrap())
    }

    /// Task spec — a 4-part chain a → b → c → d with a fixed at depth 0.
    /// Exploded along +Z should produce monotonically increasing Z
    /// offsets ordered by mate hierarchy.
    #[test]
    fn four_part_chain_explodes_along_z() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("b"));
        let id_c = a.add_part(unit_cube("c"));
        let id_d = a.add_part(unit_cube("d"));
        for (pa, pb) in [(id_a, id_b), (id_b, id_c), (id_c, id_d)] {
            a.add_mate(Mate::new(
                0,
                MateKind::Coincident {
                    part_a: pa,
                    point_a: Vector3::zeros(),
                    part_b: pb,
                    point_b: Vector3::zeros(),
                },
            ));
        }
        let cfg = ExplodeConfig { spacing: 3.0 };
        let exp = auto_explode(&a, Vector3::z(), cfg).unwrap();

        // Steps sorted by depth: a (0), b (1), c (2), d (3).
        let z = |i: usize| exp.steps[i].offset.z;
        assert_eq!(exp.steps[0].part_id, id_a);
        assert_eq!(exp.steps[0].depth, 0);
        assert!(z(0).abs() < 1e-12, "root should not move");
        assert_eq!(exp.steps[1].part_id, id_b);
        assert_eq!(exp.steps[1].depth, 1);
        assert!((z(1) - 3.0).abs() < 1e-12, "b at 3.0");
        assert_eq!(exp.steps[2].part_id, id_c);
        assert_eq!(exp.steps[2].depth, 2);
        assert!((z(2) - 6.0).abs() < 1e-12, "c at 6.0");
        assert_eq!(exp.steps[3].part_id, id_d);
        assert_eq!(exp.steps[3].depth, 3);
        assert!((z(3) - 9.0).abs() < 1e-12, "d at 9.0");

        // Z offsets strictly increase by mate hierarchy.
        for w in exp.steps.windows(2) {
            assert!(w[0].offset.z <= w[1].offset.z);
        }
    }

    /// Multiple fixed parts (multiple roots) get depth 0; their immediate
    /// children get depth 1.
    #[test]
    fn multiple_fixed_parts_are_all_depth_zero() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("root1");
        p0.fixed = true;
        let id_r1 = a.add_part(p0);
        let mut p1 = unit_cube("root2");
        p1.fixed = true;
        let id_r2 = a.add_part(p1);
        let id_c = a.add_part(unit_cube("child"));
        a.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: id_r1,
                point_a: Vector3::zeros(),
                part_b: id_c,
                point_b: Vector3::zeros(),
            },
        ));
        let exp = auto_explode(&a, Vector3::x(), ExplodeConfig::default()).unwrap();
        let depth = |id: usize| exp.steps.iter().find(|s| s.part_id == id).unwrap().depth;
        assert_eq!(depth(id_r1), 0);
        assert_eq!(depth(id_r2), 0);
        assert_eq!(depth(id_c), 1);
    }

    /// A part with no path to any fixed root stays at depth 0 (treated
    /// as a standalone root).
    #[test]
    fn disconnected_part_is_treated_as_root() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("a");
        p0.fixed = true;
        a.add_part(p0);
        let id_lonely = a.add_part(unit_cube("lonely"));
        let exp = auto_explode(&a, Vector3::y(), ExplodeConfig::default()).unwrap();
        let lonely_depth = exp
            .steps
            .iter()
            .find(|s| s.part_id == id_lonely)
            .unwrap()
            .depth;
        assert_eq!(lonely_depth, 0, "disconnected part should be a root");
    }

    /// Zero-vector direction is rejected.
    #[test]
    fn zero_direction_errors() {
        let mut a = Assembly::new();
        a.add_part(unit_cube("p"));
        let err = auto_explode(&a, Vector3::zeros(), ExplodeConfig::default()).unwrap_err();
        assert_eq!(err.code(), "assembly.bad_parameter");
    }

    /// Spacing scales linearly. Doubling spacing doubles every non-root
    /// offset.
    #[test]
    fn spacing_scales_offsets_linearly() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("root");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("b"));
        a.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
            },
        ));
        let dir = Vector3::z();
        let off_at = |spacing: f64| {
            let cfg = ExplodeConfig { spacing };
            auto_explode(&a, dir, cfg)
                .unwrap()
                .steps
                .iter()
                .find(|s| s.part_id == id_b)
                .unwrap()
                .offset
                .z
        };
        let o1 = off_at(1.0);
        let o2 = off_at(2.0);
        assert!((o2 - 2.0 * o1).abs() < 1e-12, "{o2} != 2 * {o1}");
    }

    /// Exploded poses preserve orientation; only translation moves.
    #[test]
    fn explode_preserves_orientation() {
        use nalgebra::UnitQuaternion;
        let mut a = Assembly::new();
        let mut p0 = unit_cube("root");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let mut p1 = unit_cube("b");
        // Give b a non-trivial orientation.
        p1.transform.orientation = UnitQuaternion::from_axis_angle(&Vector3::y_axis(), 1.2);
        let id_b = a.add_part(p1);
        a.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
            },
        ));
        let exp = auto_explode(&a, Vector3::z(), ExplodeConfig::default()).unwrap();
        let (_, xform_b) = exp.poses.iter().find(|(id, _)| *id == id_b).unwrap();
        let original = UnitQuaternion::from_axis_angle(&Vector3::y_axis(), 1.2);
        let dot = xform_b.orientation.coords.dot(&original.coords).abs();
        assert!((dot - 1.0).abs() < 1e-12, "orientation drifted");
    }

    /// `linear_explode_steps` is the same as `auto_explode(...).steps`.
    #[test]
    fn linear_steps_match_auto_explode() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("root");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("b"));
        a.add_mate(Mate::new(
            0,
            MateKind::Coincident {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
            },
        ));
        let dir = Vector3::x();
        let cfg = ExplodeConfig::default();
        let steps_direct = linear_explode_steps(&a, dir, cfg.clone()).unwrap();
        let steps_via_full = auto_explode(&a, dir, cfg).unwrap().steps;
        assert_eq!(steps_direct, steps_via_full);
    }

    /// Suppressed mates don't contribute to the depth graph — a depth-1
    /// part with its only mate suppressed becomes depth 0.
    #[test]
    fn suppressed_mates_dont_shape_depth() {
        let mut a = Assembly::new();
        let mut p0 = unit_cube("root");
        p0.fixed = true;
        let id_a = a.add_part(p0);
        let id_b = a.add_part(unit_cube("b"));
        let mut m = Mate::new(
            0,
            MateKind::Coincident {
                part_a: id_a,
                point_a: Vector3::zeros(),
                part_b: id_b,
                point_b: Vector3::zeros(),
            },
        );
        m.suppressed = true;
        a.add_mate(m);
        let exp = auto_explode(&a, Vector3::z(), ExplodeConfig::default()).unwrap();
        let b_depth = exp
            .steps
            .iter()
            .find(|s| s.part_id == id_b)
            .unwrap()
            .depth;
        assert_eq!(b_depth, 0, "b should not be reachable through suppressed");
    }
}
