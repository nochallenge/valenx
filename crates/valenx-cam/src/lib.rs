//! # valenx-cam
//!
//! CAM (Computer-Aided Manufacturing) workbench for Valenx (Phase 10
//! of the FreeCAD-parity roadmap).
//!
//! ## What's in here
//!
//! - [`Tool`] / [`ToolKind`] — tool library (end-mill, ball-mill,
//!   drill, face-mill, tap, reamer) with validation.
//! - [`Stock`] — rectangular stock block (origin + extents).
//! - [`Toolpath`] / [`Move`] / [`MoveKind`] — the canonical
//!   sequence-of-moves output every operation emits.
//! - [`Operation`] enum — `Profile`, `Pocket`, `Drill`, `Face`, each
//!   with its own params struct.
//! - [`offset`] / [`raster`] — 2D polygon primitives used by Pocket
//!   and Face operations.
//! - [`op`] — per-operation `generate` functions that take a stock +
//!   mesh + params + tool and produce a [`Toolpath`].
//! - [`post`] — postprocessor trait + `Grbl` / `LinuxCnc` / `Fanuc`
//!   implementations that turn a Toolpath into G-code.
//! - [`simulate`] — group consecutive moves into polylines for
//!   rendering + estimate cycle time / removed volume.
//! - [`persist`] — RON envelope `CamFile` that round-trips tools +
//!   stock + operations (the regenerated toolpath is not persisted).
//!
//! ## v1 scope honesty
//!
//! - **Rectangular stock only** — no shaped stock or wrap-around.
//! - **4 op types** — Profile, Pocket, Drill, Face. No adaptive
//!   clearing, no 4+ axis, no automatic tool-change sequencing,
//!   no multi-setup.
//! - **Geometric algorithms only** — operations work on the
//!   tessellated mesh of the source solid (via Phase 0's
//!   [`valenx_mesh::cut::intersect_plane_triangles`]).
//! - **Polygon offset** (`offset::polygon`) — per-vertex bisector
//!   offset. Works on convex / mostly-convex polygons; fails on
//!   highly-concave geometry or self-intersection.
//! - **Toolpath simulation** — draws coloured polylines in the
//!   viewport via egui's `Painter`. No material-removal animation.
//!
//! ## Commercial-depth modules (post-v1)
//!
//! - [`engagement`] — XY occupancy-grid engagement-angle analysis;
//!   the load proxy modern HSM (Mastercam Dynamic / HSMWorks /
//!   Fusion 360 Adaptive) bounds.
//! - [`op::adaptive_constant_engagement`] — toolpath generator that
//!   keeps the engagement angle under a configurable bound by
//!   inserting trochoidal roll-overs at corners (the real HSM
//!   adaptive-clearing algorithm).
//! - [`arcfit`] — G2/G3 circular-arc fitting via least-squares
//!   Kåsa fit + chord-error tolerance. Reduces G-code line count
//!   on rounded paths by 50-95 %, lets the controller's
//!   centripetal-acceleration lookahead plan correctly.
//! - [`feedrate`] — three-pass feedrate optimization
//!   (centripetal bound + sharp-corner bound + lookahead
//!   backward decel ramp). Production-CAM-class.
//! - [`collision`] — continuous swept-volume collision detection
//!   of cutter + holder against workpiece + fixture AABBs along
//!   the entire toolpath. Catches grazing mid-move collisions
//!   the [`fixture`] endpoint-only check misses.

#![forbid(unsafe_code)]
#![warn(missing_docs)]
// We deliberately use `!(x > 0.0)` in validators so NaN is treated
// as invalid (the equivalent `x <= 0.0` returns `false` for NaN
// and would let bad parameters through). Suppress the clippy lint
// for the whole crate.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

pub mod arcfit;
pub mod axis;
pub mod collision;
pub mod engagement;
pub mod error;
pub mod feedrate;
pub mod fixture;
pub mod gcode;
pub mod offset;
pub mod op;
pub mod operation;
pub mod persist;
pub mod post;
pub mod raster;
pub mod setup;
pub mod simulate;
pub mod stock;
pub mod tool;
pub mod toolpath;
pub mod voxel;
pub mod wear;

pub use arcfit::{fit_arcs, ArcDir, ArcFitParams, ArcFitReport};
pub use collision::{
    continuous_collision_check, CollisionBody, CollisionSetup, ContinuousCollision,
    ContinuousCollisionParams, Holder, HolderSegment, SetupPart, SetupPartKind,
};
pub use engagement::{engagement_along, engagement_at, EngagementSample, StockGrid};
pub use error::CamError;
pub use feedrate::{FeedrateParams, FeedrateReport};
pub use gcode::{to_gcode, to_gcode_checked, to_gcode_with, GcodeOptions};
pub use op::adaptive_constant_engagement::{
    AdaptiveConstantEngagementParams, AdaptiveEngagementReport,
};
pub use operation::{
    DrillParams, FaceParams, Operation, PocketParams, PocketStrategy, ProfileParams,
};
pub use post::{PostKind, Postprocessor};
pub use stock::Stock;
pub use tool::{Tool, ToolKind};
pub use toolpath::{Move, MoveKind, Toolpath};
