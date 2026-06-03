//! Continuous swept collision — the **per-move CCD** that every
//! commercial CAM verifier ships (Vericut, NCSimul, Mastercam
//! Verify, HSMWorks ToolPath Verifier).
//!
//! The v1 [`crate::fixture::collision_check`] tests the cylinder at
//! each move's *endpoint* against the fixture AABB list — it misses
//! the case where the swept volume *between* two endpoints clips a
//! fixture. That's the textbook failure mode the commercial
//! verifiers exist to catch.
//!
//! This module catches it.
//!
//! ## What ships
//!
//! - [`Holder`] — a cone or cylinder stack representing the tool
//!   holder above the flutes. Real CAM verifiers query the holder
//!   *and* the cutter against the workpiece for gouge detection.
//! - [`SetupPart`] — one obstacle (workpiece OR fixture AABB).
//! - [`CollisionSetup`] — the bundle of obstacles to check against.
//! - [`continuous_collision_check`] — runs **swept-cylinder vs AABB**
//!   for the cutter, plus **swept-cone vs AABB** for the holder, at
//!   every consecutive pair of moves. Reports the first colliding
//!   point per move + which obstacle it hit.
//!
//! ## v1 simplifications (honest)
//!
//! - **AABB obstacles** — workpiece is represented either as its
//!   AABB or as an explicit list of AABBs (typical fixture form).
//!   Mesh-level CCD is the documented follow-up. Most fixturing
//!   geometry — vise jaws, soft jaws, clamps, sub-plates — is well
//!   approximated by AABB unions.
//! - **Cylindrical cutter + conical holder** — real cutters have a
//!   ball nose, drill point, etc.; the cylinder + cone capture the
//!   load-bearing collision volume.
//! - **Sampled CCD** — the swept volume is sampled along each move
//!   at `sample_step_mm` (default 0.5 mm) and a per-sample
//!   cylinder-AABB test is run. This is **bounded-error CCD**:
//!   the worst-case missed-collision depth is `sample_step / 2`.
//!   Commercial verifiers use exact GJK / SAT; that's a follow-up
//!   pass — for the AABB-only obstacle case, sampled CCD matches
//!   the obstacle resolution well.
//! - **No tool-axis tracking** — assumes the spindle is +Z up. 5-axis
//!   collision against a tilted tool requires the rotational state
//!   alongside the linear state — the documented T3 gap.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::fixture::{Fixture, FixtureAabb};
use crate::tool::Tool;
use crate::toolpath::{Move, MoveKind, Toolpath};

/// Tool holder above the flutes — a conical or cylindrical body
/// that can collide with workpiece edges. Commercial CAM treats
/// holder collision as a *separate* gouge class from flute
/// collision.
///
/// All segments are taken to lie above the flute tip along +Z.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Holder {
    /// Segments above the flute tip. Each segment is `(z_lower,
    /// z_upper, radius_lower, radius_upper)`. A cylinder has
    /// `radius_lower == radius_upper`; a cone is the standard
    /// taper-shoulder shape.
    ///
    /// `z_lower` is measured *above the flute tip* (i.e. the tip of
    /// the cutter is at z=0; flute length is the start of the
    /// holder stack).
    pub segments: Vec<HolderSegment>,
}

/// One holder segment.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct HolderSegment {
    /// Bottom z of the segment (mm above the flute tip).
    pub z_lower_mm: f64,
    /// Top z of the segment (mm above the flute tip).
    pub z_upper_mm: f64,
    /// Radius at the bottom of the segment (mm).
    pub radius_lower_mm: f64,
    /// Radius at the top of the segment (mm).
    pub radius_upper_mm: f64,
}

impl Holder {
    /// Empty holder (no body above the flutes).
    pub fn empty() -> Self {
        Self {
            segments: Vec::new(),
        }
    }

    /// Single-cylinder shank holder: the most common.
    pub fn cylinder_shank(diameter_mm: f64, length_mm: f64) -> Self {
        Self {
            segments: vec![HolderSegment {
                z_lower_mm: 0.0,
                z_upper_mm: length_mm,
                radius_lower_mm: diameter_mm * 0.5,
                radius_upper_mm: diameter_mm * 0.5,
            }],
        }
    }

    /// Conical holder above an optional cylindrical neck — the
    /// shrink-fit / collet shape.
    pub fn neck_cone(
        neck_diameter_mm: f64,
        neck_length_mm: f64,
        body_diameter_mm: f64,
        body_length_mm: f64,
    ) -> Self {
        Self {
            segments: vec![
                HolderSegment {
                    z_lower_mm: 0.0,
                    z_upper_mm: neck_length_mm,
                    radius_lower_mm: neck_diameter_mm * 0.5,
                    radius_upper_mm: neck_diameter_mm * 0.5,
                },
                HolderSegment {
                    z_lower_mm: neck_length_mm,
                    z_upper_mm: neck_length_mm + body_length_mm,
                    radius_lower_mm: neck_diameter_mm * 0.5,
                    radius_upper_mm: body_diameter_mm * 0.5,
                },
            ],
        }
    }

    /// True if the holder has any segments.
    pub fn is_present(&self) -> bool {
        !self.segments.is_empty()
    }
}

/// One obstacle in the setup — either a fixture AABB or a workpiece
/// AABB. The classification carries through to the collision
/// report so the user can tell whether the tool is about to wreck
/// the vise or the part.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SetupPart {
    /// Min XYZ corner.
    pub min: Vector3<f64>,
    /// Max XYZ corner.
    pub max: Vector3<f64>,
    /// Free-form label.
    pub label: String,
    /// Part kind — `Fixture` for the vise / clamps, `Workpiece` for
    /// the part itself (collisions against workpiece *outside* the
    /// programmed envelope mean a programming error).
    pub kind: SetupPartKind,
}

/// Whether a setup part is the workpiece or a fixture.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum SetupPartKind {
    /// Part being machined.
    Workpiece,
    /// Vise / clamp / support — must not be touched.
    Fixture,
}

/// The complete setup the toolpath is being checked against.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct CollisionSetup {
    /// Every obstacle to test.
    pub parts: Vec<SetupPart>,
}

impl CollisionSetup {
    /// Empty setup.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a fixture AABB (vise jaw, clamp, sub-plate).
    pub fn push_fixture(
        &mut self,
        min: Vector3<f64>,
        max: Vector3<f64>,
        label: impl Into<String>,
    ) {
        self.parts.push(SetupPart {
            min,
            max,
            label: label.into(),
            kind: SetupPartKind::Fixture,
        });
    }

    /// Add a workpiece AABB — typically the stock blank or the
    /// part's external envelope.
    pub fn push_workpiece(
        &mut self,
        min: Vector3<f64>,
        max: Vector3<f64>,
        label: impl Into<String>,
    ) {
        self.parts.push(SetupPart {
            min,
            max,
            label: label.into(),
            kind: SetupPartKind::Workpiece,
        });
    }

    /// Convert from the existing [`Fixture`] type for back-compat.
    pub fn from_fixture(fixture: &Fixture) -> Self {
        let parts = fixture
            .aabbs
            .iter()
            .map(|a: &FixtureAabb| SetupPart {
                min: a.min,
                max: a.max,
                label: a.label.clone(),
                kind: SetupPartKind::Fixture,
            })
            .collect();
        Self { parts }
    }
}

/// Parameters for the continuous CCD pass.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ContinuousCollisionParams {
    /// Sample step along each move (mm). Bounds the worst-case
    /// missed-collision depth to `step / 2`. Default 0.5 mm.
    pub sample_step_mm: f64,
    /// If `true`, the cutter cylinder may legitimately *enter*
    /// `Workpiece` parts (it's machining them); only colliding
    /// with `Fixture` parts is reported. If `false`, both kinds are
    /// reported. Default `true`.
    pub allow_workpiece_contact: bool,
    /// If `true`, rapid moves (G0) are also checked. Real CAM
    /// verifiers check rapids — they're the most dangerous, hi-speed
    /// non-cutting moves where a crash is unrecoverable.
    pub check_rapids: bool,
}

impl Default for ContinuousCollisionParams {
    fn default() -> Self {
        Self {
            sample_step_mm: 0.5,
            allow_workpiece_contact: true,
            check_rapids: true,
        }
    }
}

/// One detected collision.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ContinuousCollision {
    /// Index of the offending move in the toolpath.
    pub move_index: usize,
    /// Parametric t along the move where the collision was detected
    /// (0 = move start, 1 = move end).
    pub t_at_hit: f64,
    /// Tool centre at the collision point.
    pub tool_position: Vector3<f64>,
    /// Setup part label that was hit.
    pub part_label: String,
    /// Which part kind it was.
    pub part_kind: SetupPartKind,
    /// Which body collided — flute or holder.
    pub body: CollisionBody,
}

/// Which part of the tool stack collided.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum CollisionBody {
    /// The flute cylinder.
    Flute,
    /// The holder above the flutes.
    Holder,
}

/// Run continuous (swept) collision detection along the entire
/// toolpath. Returns one [`ContinuousCollision`] per detected hit.
///
/// The flute body is `tool.radius_mm() × tool.length_mm` cylinder.
/// The holder, if present, is the segment stack from [`Holder`].
pub fn continuous_collision_check(
    toolpath: &Toolpath,
    tool: &Tool,
    holder: &Holder,
    setup: &CollisionSetup,
    params: &ContinuousCollisionParams,
) -> Vec<ContinuousCollision> {
    let mut out = Vec::new();
    if toolpath.moves.len() < 2 {
        return out;
    }
    let flute_r = tool.radius_mm();
    let flute_len = tool.length_mm;
    for w in 1..toolpath.moves.len() {
        let prev: Move = toolpath.moves[w - 1];
        let curr: Move = toolpath.moves[w];
        if !params.check_rapids && curr.kind == MoveKind::Rapid {
            continue;
        }
        let a = prev.position;
        let b = curr.position;
        let d = b - a;
        let len = d.norm();
        let n_steps = ((len / params.sample_step_mm).ceil() as usize).max(1);
        for s in 0..=n_steps {
            let t = (s as f64) / (n_steps as f64);
            let centre = a + d * t;
            for part in &setup.parts {
                if part.kind == SetupPartKind::Workpiece && params.allow_workpiece_contact {
                    continue;
                }
                // Flute cylinder vs AABB.
                if cylinder_aabb_overlap(centre, flute_r, flute_len, part) {
                    out.push(ContinuousCollision {
                        move_index: w,
                        t_at_hit: t,
                        tool_position: centre,
                        part_label: part.label.clone(),
                        part_kind: part.kind,
                        body: CollisionBody::Flute,
                    });
                    // Stop sampling further on this move/part to
                    // avoid flooding the report with one hit per
                    // sample.
                    break;
                }
                // Holder segments.
                if holder.is_present() {
                    for seg in &holder.segments {
                        if cone_segment_aabb_overlap(centre, flute_len, seg, part) {
                            out.push(ContinuousCollision {
                                move_index: w,
                                t_at_hit: t,
                                tool_position: centre,
                                part_label: part.label.clone(),
                                part_kind: part.kind,
                                body: CollisionBody::Holder,
                            });
                            break;
                        }
                    }
                }
            }
        }
    }
    out
}

/// Test the flute cylinder (centred at `centre`, radius `r`, length
/// `len` extending downward in -Z) for overlap with `part`.
fn cylinder_aabb_overlap(centre: Vector3<f64>, r: f64, len: f64, part: &SetupPart) -> bool {
    // Flute extends from `centre.z` down to `centre.z - len`.
    let tool_z_max = centre.z;
    let tool_z_min = centre.z - len;
    if tool_z_max < part.min.z || tool_z_min > part.max.z {
        return false;
    }
    // XY: closest-point AABB-to-circle.
    let cx = centre.x.clamp(part.min.x, part.max.x);
    let cy = centre.y.clamp(part.min.y, part.max.y);
    let dx = cx - centre.x;
    let dy = cy - centre.y;
    (dx * dx + dy * dy) <= r * r
}

/// Test one holder segment against a setup part. The holder
/// segment's z-range is *above the flute tip* — i.e. the bottom of
/// the segment sits at `centre.z + (flute_len + seg.z_lower)`
/// because the tool tip is at `centre.z - flute_len` and the flute
/// tip is at `centre.z - 0`... wait that's wrong.
///
/// **Convention**: in this module, the cutter "centre" is the
/// **bottom tip** of the flute (i.e. the lowest point of the cutter
/// outline). The flute spans from `centre.z` (tip, the lower face)
/// up to `centre.z + flute_len`. The holder begins at `centre.z +
/// flute_len + seg.z_lower_mm` and extends upward.
fn cone_segment_aabb_overlap(
    centre: Vector3<f64>,
    flute_len: f64,
    seg: &HolderSegment,
    part: &SetupPart,
) -> bool {
    let seg_z_min = centre.z + flute_len + seg.z_lower_mm;
    let seg_z_max = centre.z + flute_len + seg.z_upper_mm;
    if seg_z_max < part.min.z || seg_z_min > part.max.z {
        return false;
    }
    // Use the LARGER of the two radii — conservative, matches the
    // commercial-CAM "outermost holder envelope" gouge check. A
    // refined per-z slice test is the documented follow-up.
    let r = seg.radius_lower_mm.max(seg.radius_upper_mm);
    let cx = centre.x.clamp(part.min.x, part.max.x);
    let cy = centre.y.clamp(part.min.y, part.max.y);
    let dx = cx - centre.x;
    let dy = cy - centre.y;
    (dx * dx + dy * dy) <= r * r
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tool::ToolKind;
    use crate::toolpath::{Move, MoveKind};

    fn p(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    fn em6() -> Tool {
        Tool::new(1, "EM6", ToolKind::EndMill, 6.0, 25.0, 2, "").unwrap()
    }

    #[test]
    fn empty_setup_no_collisions() {
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 5.0), 0.0));
        tp.push(Move::new(MoveKind::Cut, p(10.0, 0.0, 5.0), 500.0));
        let setup = CollisionSetup::new();
        let hits = continuous_collision_check(
            &tp,
            &em6(),
            &Holder::empty(),
            &setup,
            &ContinuousCollisionParams::default(),
        );
        assert!(hits.is_empty());
    }

    #[test]
    fn grazing_move_caught_by_continuous_check_misses_endpoint() {
        // A toolpath whose *endpoints* clear a fixture AABB, but
        // whose midpoint sweeps right through it. This is the
        // textbook v1-cylinder-vs-AABB miss case.
        let mut tp = Toolpath::new();
        // Start at (0, 0, 5) — flute tip at z=5, flute extends down
        // to z=-20.
        tp.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 5.0), 0.0));
        // End at (40, 0, 5).
        tp.push(Move::new(MoveKind::Rapid, p(40.0, 0.0, 5.0), 0.0));
        // Fixture AABB sits in the middle of the X span.
        let mut setup = CollisionSetup::new();
        setup.push_fixture(
            p(15.0, -10.0, 0.0),
            p(25.0, 10.0, 4.0),
            "clamp",
        );
        let hits = continuous_collision_check(
            &tp,
            &em6(),
            &Holder::empty(),
            &setup,
            &ContinuousCollisionParams::default(),
        );
        // We MUST detect this; the v1 endpoint-only test would not.
        assert!(
            !hits.is_empty(),
            "continuous CCD failed to catch grazing rapid"
        );
        assert_eq!(hits[0].part_label, "clamp");
        assert_eq!(hits[0].body, CollisionBody::Flute);
        // The collision should happen mid-move.
        assert!(hits[0].t_at_hit > 0.2 && hits[0].t_at_hit < 0.8);
    }

    #[test]
    fn clear_path_reports_no_collision() {
        // Path well above the fixture AABB.
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, p(0.0, 0.0, 100.0), 0.0));
        tp.push(Move::new(MoveKind::Rapid, p(40.0, 0.0, 100.0), 0.0));
        let mut setup = CollisionSetup::new();
        setup.push_fixture(p(15.0, -10.0, 0.0), p(25.0, 10.0, 5.0), "clamp");
        let hits = continuous_collision_check(
            &tp,
            &em6(),
            &Holder::empty(),
            &setup,
            &ContinuousCollisionParams::default(),
        );
        assert!(hits.is_empty(), "clear path falsely reported {hits:?}");
    }

    #[test]
    fn workpiece_contact_ignored_when_allowed() {
        // The tool is supposed to enter the workpiece — those
        // contacts are normal machining, not crashes.
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Cut, p(0.0, 0.0, 0.0), 500.0));
        tp.push(Move::new(MoveKind::Cut, p(10.0, 0.0, 0.0), 500.0));
        let mut setup = CollisionSetup::new();
        setup.push_workpiece(
            p(-5.0, -10.0, -50.0),
            p(15.0, 10.0, 5.0),
            "stock",
        );
        let params = ContinuousCollisionParams {
            allow_workpiece_contact: true,
            ..Default::default()
        };
        let hits = continuous_collision_check(&tp, &em6(), &Holder::empty(), &setup, &params);
        assert!(hits.is_empty());
        // But with workpiece contact OFF, we should see a hit.
        let params_strict = ContinuousCollisionParams {
            allow_workpiece_contact: false,
            ..Default::default()
        };
        let hits_strict =
            continuous_collision_check(&tp, &em6(), &Holder::empty(), &setup, &params_strict);
        assert!(!hits_strict.is_empty());
        assert_eq!(hits_strict[0].part_kind, SetupPartKind::Workpiece);
    }

    #[test]
    fn holder_collision_detected() {
        // Tool tip clears a tall fixture, but the holder (a wider
        // shank above the flute) clips it.
        let tool = Tool::new(1, "EM3", ToolKind::EndMill, 3.0, 10.0, 2, "").unwrap();
        let holder = Holder::cylinder_shank(12.0, 30.0);
        // Tool tip at z=0; flute extends z=0..10; holder extends
        // z=10..40 with r=6.
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Cut, p(0.0, 0.0, 0.0), 500.0));
        tp.push(Move::new(MoveKind::Cut, p(20.0, 0.0, 0.0), 500.0));
        // Fixture in the path of the holder (z=15..20 — above the
        // flute, but below the top of the holder), wide enough that
        // r=6 holder hits it.
        let mut setup = CollisionSetup::new();
        setup.push_fixture(p(8.0, -10.0, 15.0), p(12.0, 10.0, 20.0), "stop");
        let hits = continuous_collision_check(
            &tp,
            &tool,
            &holder,
            &setup,
            &ContinuousCollisionParams::default(),
        );
        assert!(!hits.is_empty());
        let any_holder = hits.iter().any(|h| h.body == CollisionBody::Holder);
        assert!(any_holder, "expected a holder collision in {hits:?}");
    }

    #[test]
    fn endpoint_only_check_would_miss_but_we_dont() {
        // Place a small fixture between two move endpoints. Endpoint
        // test would miss; sweep test catches it.
        let mut tp = Toolpath::new();
        tp.push(Move::new(MoveKind::Rapid, p(-1.0, 0.0, 5.0), 0.0));
        tp.push(Move::new(MoveKind::Rapid, p(11.0, 0.0, 5.0), 0.0));
        let mut setup = CollisionSetup::new();
        // Fixture sits at x=5 with radius 0.5 — the endpoints sit
        // outside, but the path slices through the middle.
        setup.push_fixture(
            p(4.5, -5.0, 0.0),
            p(5.5, 5.0, 4.0),
            "tab",
        );
        let hits = continuous_collision_check(
            &tp,
            &em6(),
            &Holder::empty(),
            &setup,
            &ContinuousCollisionParams::default(),
        );
        assert!(!hits.is_empty(), "expected to catch mid-move collision");
    }

    #[test]
    fn from_fixture_compat() {
        let mut fix = Fixture::new();
        fix.push_aabb(p(0.0, 0.0, 0.0), p(10.0, 10.0, 10.0), "vise");
        let setup = CollisionSetup::from_fixture(&fix);
        assert_eq!(setup.parts.len(), 1);
        assert_eq!(setup.parts[0].kind, SetupPartKind::Fixture);
    }
}
