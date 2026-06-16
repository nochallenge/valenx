//! # valenx-truss
//!
//! Planar **pin-jointed truss statics by the method of joints**: define
//! the joints, supports, two-force members and point loads of a 2-D
//! truss, then solve one linear system for the axial force in every
//! member (tension positive) and the reaction at every support.
//!
//! ## What
//!
//! A self-contained, in-process statics solver for the classic statically
//! determinate planar truss — the structure made of straight bars pinned
//! at their ends, loaded only at the joints. You assemble a [`Truss`]
//! from [`Node`]s (each optionally a [`Support`] and/or a point [`Load`])
//! and [`Member`]s, call [`solve`], and read the [`TrussSolution`]:
//!
//! - the axial force in each member, sign-conventioned **tension
//!   positive** (so a negative number is a member in compression), and
//! - the support reactions, each resolved to global `(fx, fy)`.
//!
//! Before solving you can classify the structure by the Maxwell counting
//! rule: [`Truss::static_determinacy_number`] returns the signed
//! `(m + r) - 2N`, and [`Truss::determinacy`] sorts its sign into the
//! [`Determinacy`] cases — mechanism (deficient), determinate, or
//! indeterminate (redundant).
//!
//! It is the in-process counterpart to a hand worked-example or a
//! teaching FE tool: no meshing, no external solver, runs anywhere the
//! crate compiles.
//!
//! ```
//! use valenx_truss::{Load, Member, Node, Support, Truss, solve};
//!
//! // The canonical loaded triangle:
//! //          C (2,3)  ← 60 kN downward
//! //         / \
//! //  A(0,0)*---*B(4,0)
//! //  pin      roller (slides horizontally → vertical reaction only)
//! let mut truss = Truss::new();
//! let a = truss.add_node(Node::new(0.0, 0.0)?.with_support(Support::Pin));
//! let b = truss.add_node(
//!     Node::new(4.0, 0.0)?.with_support(Support::horizontal_roller()),
//! );
//! let c = truss.add_node(Node::new(2.0, 3.0)?.with_load(Load::down(60.0)?));
//! let ab = truss.add_member(Member::new(a, b))?;
//! let _bc = truss.add_member(Member::new(b, c))?;
//! let _ca = truss.add_member(Member::new(c, a))?;
//!
//! let sol = solve(&truss)?;
//!
//! // Bottom chord AB is the analytic value P/3 = 20, in tension.
//! assert!((sol.force(ab) - 20.0).abs() < 1e-9);
//! assert!(sol.member_forces[ab].is_tension());
//! // The two inclined struts carry −10·√13 each (compression).
//! assert!((sol.force(_bc) - (-10.0 * 13.0_f64.sqrt())).abs() < 1e-9);
//! // Global equilibrium holds: ΣFx = ΣFy = 0.
//! let (rx, ry) = sol.global_residual(&truss);
//! assert!(rx.abs() < 1e-9 && ry.abs() < 1e-9);
//! # Ok::<(), valenx_truss::TrussError>(())
//! ```
//!
//! ## Model
//!
//! Every joint of a truss in equilibrium satisfies `ΣFx = 0` and
//! `ΣFy = 0`. Stacking those two equations over all `N` joints gives a
//! linear system in the unknown member forces `S₀ … S_{m-1}` and support
//! reactions `R₀ … R_{r-1}`:
//!
//! ```text
//!   A · x = b ,   x = [S₀ … S_{m-1}, R₀ … R_{r-1}]ᵀ
//! ```
//!
//! where each member of tension `S` joining nodes `a` and `b` (unit axis
//! `û = (b − a)/‖b − a‖`) contributes `+S·û` to joint `a`'s balance and
//! `−S·û` to joint `b`'s, a pin contributes its `(Rx, Ry)`, a roller its
//! single reaction along the slide-normal, and `b` is the negated applied
//! load at each joint. For a statically determinate, stable truss
//! `m + r = 2N`, so `A` is square and invertible and the system has the
//! unique solution returned by [`solve`] via [`nalgebra`]'s pivoted LU.
//! The assembly and the linear algebra live in [`solver`]; the validated
//! geometry / topology in [`model`]; the solved quantities in [`result`].
//!
//! Determinacy is enforced in two layers: the counting condition
//! `m + r = 2N` ([`Truss::is_count_determinate`], refined by the signed
//! [`Truss::static_determinacy_number`] and its [`Determinacy`]
//! classification) up front, and rank (catching count-balanced
//! *mechanisms*) by the invertibility test on the assembled matrix — a
//! singular system returns [`TrussError::Singular`] rather than a garbage
//! force state.
//!
//! ## Honest scope
//!
//! This is research / educational grade. It is the genuine textbook
//! method of joints — the member forces, the tension/compression signs,
//! and global equilibrium are all exact for the model it implements, and
//! they are checked against analytic hand-solutions in the test-suite
//! (the loaded triangle above, a square-with-diagonal, a cantilever, and
//! known zero-force members). The model it implements is the standard
//! *idealisation*, and that idealisation is the limit of its validity:
//!
//! - **Statically determinate trusses only.** The method of joints
//!   solves `m + r = 2N`; statically *indeterminate* (redundant) trusses
//!   need member stiffnesses and the displacement/stiffness method, which
//!   this crate does not implement (it reports
//!   [`TrussError::NotDeterminate`]).
//! - **Pin joints, two-force members, joint loads only.** Bars carry
//!   pure axial force; there are no bending moments, no rigid (welded)
//!   connections, and loads applied along a member's span are not
//!   modelled — put loads at joints.
//! - **Rigid-body, small-displacement, linear statics.** Equilibrium is
//!   written on the *undeformed* geometry: no member axial flexibility,
//!   no joint displacements, no second-order / large-displacement (P-Δ)
//!   effects, and therefore no buckling check on the compression members
//!   this solver happily reports.
//! - **No self-weight, dynamics, or temperature.** Members are
//!   weightless; there is no inertia, no time history, and no thermal
//!   strain.
//!
//! It is **not** a clinical/medical tool and **not** a production or
//! code-checked structural-engineering tool: it does not size members,
//! apply load factors, or check any building code. For a real structure,
//! use qualified engineering software and a licensed engineer; for
//! flexible / indeterminate / 3-D frames, reach for Valenx's `valenx-fem`
//! finite-element crate.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod model;
pub mod result;
pub mod solver;

pub use error::{ErrorCategory, TrussError};
pub use model::{Determinacy, Load, Member, Node, Support, Truss};
pub use result::{AxialState, MemberForce, Reaction, TrussSolution};
pub use solver::solve;
