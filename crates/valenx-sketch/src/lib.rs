//! # valenx-sketch
//!
//! 2D parametric sketcher for Valenx. Primitives (points, lines,
//! circles, arcs) live in a [`Sketch`]; geometric constraints
//! (coincident, horizontal, distance, etc.) are added between them;
//! a Newton-Raphson solver with Levenberg-Marquardt damping drives
//! the variables to a configuration that satisfies all constraints
//! simultaneously.
//!
//! Sketches can be extruded into 3-D solids via [`extrude::extrude`]
//! using the `valenx-cad` BRep kernel.
//!
//! This is Phase 1 of the FreeCAD-parity roadmap.
//!
//! # Example
//!
//! Sketch a rough quadrilateral, lock it into a unit square with
//! [horizontal](constraint::Constraint::Horizontal),
//! [vertical](constraint::Constraint::Vertical), and
//! [distance](constraint::Constraint::Distance) constraints, run the
//! solver, then extrude into a 3-D solid:
//!
//! ```
//! use valenx_sketch::{Sketch, SolverConfig};
//! use valenx_sketch::constraint::Constraint;
//! use valenx_sketch::solver::{solve, SolverStatus};
//!
//! // 1. Build the sketch — four points roughly arranged in a square.
//! let mut sketch = Sketch::new();
//! let a = sketch.add_point(0.0, 0.0);
//! let b = sketch.add_point(0.9, 0.1); // intentionally off
//! let c = sketch.add_point(1.0, 1.0);
//! let d = sketch.add_point(0.0, 1.0);
//!
//! // 2. Connect with four line segments forming a closed loop.
//! let ab = sketch.add_line(a, b).unwrap();
//! let bc = sketch.add_line(b, c).unwrap();
//! let cd = sketch.add_line(c, d).unwrap();
//! let da = sketch.add_line(d, a).unwrap();
//!
//! // 3. Apply geometric constraints — make it a unit square.
//! sketch.add_constraint(Constraint::Horizontal(ab));
//! sketch.add_constraint(Constraint::Vertical(bc));
//! sketch.add_constraint(Constraint::Horizontal(cd));
//! sketch.add_constraint(Constraint::Vertical(da));
//! sketch.add_constraint(Constraint::Distance { a, b, target: 1.0 });
//! sketch.add_constraint(Constraint::Distance { a: b, b: c, target: 1.0 });
//!
//! // 4. Solve. The solver drives the residuals to zero.
//! let report = solve(&mut sketch, SolverConfig::default()).unwrap();
//! assert_eq!(report.status, SolverStatus::Converged);
//! assert!(report.residual_norm < 1e-6);
//!
//! // 5. Extrude the closed profile +Z by 0.5 to produce a 3-D solid.
//! let _solid = sketch.extrude(0.5).unwrap();
//! ```
//!
//! ## Persistence
//!
//! Sketches serialize to RON via [`persist::SketchFile::write_to`] /
//! [`persist::SketchFile::read_from`].
//!
//! ## Diagnostics
//!
//! [`solver::SolverReport::diagnostics`] reports residual count vs DOF
//! count so under- and over-constrained sketches are caught instead of
//! silently producing garbage geometry. See
//! [`solver::SolverDiagnostics`] for field-level meaning.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod constraint;
pub mod construction;
pub mod error;
pub mod external_geom;
pub mod extrude;
pub mod geom;
pub mod geom_bspline;
pub mod geom_ellipse;
pub mod ops;
pub mod persist;
pub mod sketch;
pub mod solver;

pub use error::SketchError;
pub use geom::{Arc2, Circle2, EntityId, Line2, Point2};
pub use geom_bspline::BSpline2;
pub use geom_ellipse::{Ellipse2, EllipticalArc2};
pub use sketch::Sketch;
pub use solver::{SolverConfig, SolverDiagnostics, SolverReport, SolverStatus};
