//! Constraint-based modelling — features 22 and 23.
//!
//! Flux-balance analysis and its relatives — the COBRA-toolbox family
//! of methods. Unlike the ODE / SSA layers these do not simulate
//! dynamics; they solve a linear program over the steady-state flux
//! cone of a metabolic network.
//!
//! - [`simplex`] — a from-scratch two-phase primal-simplex LP solver
//!   (no external LP dependency).
//! - [`flux`] — [`FbaProblem`]: flux-balance analysis (feature 22),
//!   flux variability analysis and parsimonious FBA (feature 23).

pub mod flux;
pub mod simplex;

pub use flux::{FbaProblem, FbaSolution, FluxRange};
pub use simplex::{solve_lp, ConstraintSense, LinearProgram, LpSolution, LpStatus};
