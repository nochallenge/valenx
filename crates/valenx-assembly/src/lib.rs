//! # valenx-assembly
//!
//! Assembly workbench — multi-part scene with rigid-body placement,
//! geometric mate constraints, joint types (Fixed / Revolute /
//! Prismatic / Cylindrical / Spherical / Planar), 3D constraint
//! solver, and kinematic motion preview.
//!
//! This is Phase 6 of the FreeCAD-parity roadmap.
//!
//! # Example
//!
//! Two unit cubes — one fixed at origin, one floating — connected by a
//! single Coincident mate that pins a vertex of the moving cube to the
//! fixed cube's origin:
//!
//! ```
//! use nalgebra::Vector3;
//! use valenx_assembly::{Assembly, Part, PartTransform, Mate, MateKind};
//! use valenx_assembly::solver::{solve, SolverConfig};
//! use valenx_cad::box_solid;
//!
//! let mut a = Assembly::new();
//! let cube_a = box_solid(1.0, 1.0, 1.0).unwrap();
//! let cube_b = box_solid(1.0, 1.0, 1.0).unwrap();
//!
//! // Anchor part — won't move during solve.
//! let mut fixed = Part::new(0, "Fixed", cube_a);
//! fixed.fixed = true;
//! let id_a = a.add_part(fixed);
//!
//! // Moving part — initial transform far from origin.
//! let mut moving = Part::new(1, "Moving", cube_b);
//! moving.transform.translation = Vector3::new(3.0, 4.0, 0.0);
//! let id_b = a.add_part(moving);
//!
//! // Coincident mate: pin vertex (0,0,0) of moving to (0,0,0) of fixed.
//! a.add_mate(Mate::new(0, MateKind::Coincident {
//!     part_a: id_a, point_a: Vector3::new(0.0, 0.0, 0.0),
//!     part_b: id_b, point_b: Vector3::new(0.0, 0.0, 0.0),
//! }));
//!
//! let report = solve(&mut a, SolverConfig::default()).unwrap();
//! assert!(report.residual_norm < 1e-6);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod assembly;
pub mod diagnostics;
pub mod drag;
pub mod error;
pub mod explode;
pub mod interference;
pub mod joint;
pub mod kinematics;
pub mod mate;
pub mod part;
pub mod persist;
pub mod solver;
pub mod urdf;

pub use assembly::Assembly;
pub use diagnostics::{diagnose, ConstraintState, DiagnosticsConfig, DiagnosticsReport};
pub use drag::{drag_part, DragOutcome};
pub use error::AssemblyError;
pub use explode::{
    auto_explode, linear_explode_steps, ExplodeConfig, ExplodeStep, ExplodedAssembly,
};
pub use interference::{detect_interference, Interference, InterferenceConfig};
pub use joint::{Joint, JointKind};
pub use mate::{Mate, MateKind};
pub use part::{Part, PartTransform};
pub use persist::AssemblyFile;
pub use solver::{SolverConfig, SolverDiagnostics, SolverReport, SolverStatus};
pub use urdf::{assembly_to_mesh, demo_hand_urdf, import_urdf, UrdfError, UrdfRobot};
