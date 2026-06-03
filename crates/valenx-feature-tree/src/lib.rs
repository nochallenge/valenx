//! # valenx-feature-tree
//!
//! Parametric feature history for Valenx. A [`FeatureTree`] holds an
//! ordered list of [`Feature`]s; [`replay::replay`] walks the tree
//! top-to-bottom and produces a final [`valenx_cad::Solid`].
//!
//! Phase 2 of the FreeCAD-parity roadmap.
//!
//! ## End-to-end example
//!
//! Build a sketch, drop it into the tree, pad it, pocket a hole out of
//! the result, replay the tree, and tessellate the final solid for the
//! viewport:
//!
//! ```
//! use valenx_feature_tree::{Feature, FeatureTree, replay};
//! use valenx_feature_tree::feature::{PadParams, PocketParams};
//!
//! // 1. Build a sketch — a 4×4 square on the XY plane centered at origin.
//! let mut profile = valenx_sketch::Sketch::new();
//! let a = profile.add_point(-2.0, -2.0);
//! let b = profile.add_point(2.0, -2.0);
//! let c = profile.add_point(2.0, 2.0);
//! let d = profile.add_point(-2.0, 2.0);
//! profile.add_line(a, b).unwrap();
//! profile.add_line(b, c).unwrap();
//! profile.add_line(c, d).unwrap();
//! profile.add_line(d, a).unwrap();
//!
//! // 2. Drop it into a fresh tree and Pad it to a 2-unit-thick block.
//! let mut tree = FeatureTree::new();
//! let sketch_ref = tree.add_sketch(profile);
//! tree.add_feature(
//!     Feature::Pad(PadParams {
//!         sketch: sketch_ref,
//!         depth: 2.0.into(),
//!         direction_positive: true,
//!     }),
//!     "Base Pad",
//! );
//!
//! // 3. Add a smaller square to pocket a hole out of the top.
//! let mut hole = valenx_sketch::Sketch::new();
//! let p1 = hole.add_point(-1.0, -1.0);
//! let p2 = hole.add_point(1.0, -1.0);
//! let p3 = hole.add_point(1.0, 1.0);
//! let p4 = hole.add_point(-1.0, 1.0);
//! hole.add_line(p1, p2).unwrap();
//! hole.add_line(p2, p3).unwrap();
//! hole.add_line(p3, p4).unwrap();
//! hole.add_line(p4, p1).unwrap();
//! let hole_ref = tree.add_sketch(hole);
//! tree.add_feature(
//!     Feature::Pocket(PocketParams {
//!         sketch: hole_ref,
//!         depth: 1.0.into(),
//!         direction_positive: true,
//!     }),
//!     "Central Pocket",
//! );
//!
//! // 4. Replay the tree to get the final solid, then tessellate.
//! let solid = replay(&tree).unwrap().expect("non-empty tree → Some(solid)");
//! let mesh = valenx_cad::solid_to_mesh(&solid, 0.25).unwrap();
//! assert!(mesh.total_elements() > 0);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod feature;
pub mod ops;
pub mod persist;
pub mod replay;
pub mod threads;
pub mod tree;

pub use error::FeatureError;
pub use feature::Feature;
pub use replay::{replay, replay_with_spreadsheet};
pub use tree::FeatureTree;
