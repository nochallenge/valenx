//! # valenx-diffusion
//!
//! One-dimensional **Fickian diffusion** — the two Fick laws, an
//! explicit finite-difference solver for the transient problem, and the
//! closed-form point-source solution, in pure Rust with no heavy or
//! platform dependencies.
//!
//! ## What
//!
//! Three small, well-tested pieces of the classical diffusion theory:
//!
//! - **Fick's first law** — the flux/gradient relation `J = -D dC/dx`
//!   ([`first_law_flux`], [`flux_central`]). Matter moves down the
//!   gradient.
//! - **Fick's second law** — the diffusion (heat) equation
//!   `dC/dt = D d2C/dx2`, solved transiently by an explicit
//!   forward-time centred-space (FTCS) scheme on a uniform grid
//!   ([`Grid`], [`Field`], [`Field::step`]), with fixed-concentration
//!   ([`Boundary::Dirichlet`]) and zero-flux ([`Boundary::NoFlux`])
//!   boundaries.
//! - **The instantaneous point source** — the analytic Gaussian
//!   (heat-kernel) solution
//!   `C(x, t) = M / sqrt(4 pi D t) * exp(-x^2 / (4 D t))`
//!   ([`gaussian_point_source`]), whose spreading variance is the exact
//!   `2 D t` ([`gaussian_variance`]), and the inverse time-to-spread
//!   `t = var / (2 D)` ([`time_to_reach_variance`] /
//!   [`time_to_reach_std`]).
//!
//! Plus the **steady state** between two fixed walls — the linear
//! gradient that solves `d2C/dx2 = 0` and the uniform flux it carries
//! ([`steady_profile`], [`steady_gradient`], [`steady_flux`]).
//!
//! ## Model
//!
//! The transient solver discretises Fick's second law as
//!
//! ```text
//!   C_i^{k+1} = C_i^k + r ( C_{i+1}^k - 2 C_i^k + C_{i-1}^k ),
//!   r = D dt / dx^2,
//! ```
//!
//! the explicit Euler / centred-Laplacian stencil. It is conditionally
//! stable: the diffusion number `r` must not exceed `1/2`, i.e. the
//! time step is clamped to
//!
//! ```text
//!   dt <= dx^2 / (2 D)
//! ```
//!
//! ([`max_stable_dt`], [`is_stable_dt`], [`Grid::stable_dt`]). A step
//! above the limit is rejected with [`DiffusionError::Unstable`] rather
//! than allowed to blow up. Zero-flux walls use a ghost node mirroring
//! the first interior node, so a closed domain conserves total mass to
//! round-off; Dirichlet walls pin their end node each step.
//!
//! These are the textbook results: the FTCS stability bound, the
//! one-dimensional heat kernel / Green's function, the linear
//! steady-state Laplace solution, and the `var = 2 D t` spreading law.
//!
//! ```
//! use valenx_diffusion::{Boundary, Field, Grid};
//!
//! // A unit spike in a closed (no-flux) box spreads but conserves mass.
//! let grid = Grid::new(101, 0.1).unwrap();
//! let mut field =
//!     Field::point_source(grid, 50, Boundary::NoFlux, Boundary::NoFlux).unwrap();
//! let mass_before = field.total_mass();
//!
//! let dt = grid.stable_dt(1.0, 0.4).unwrap(); // 40% of the FTCS limit
//! field.advance(1.0, dt, 500).unwrap();
//!
//! let mass_after = field.total_mass();
//! assert!((mass_after - mass_before).abs() < 1e-9); // mass conserved
//! ```
//!
//! ## Honest scope
//!
//! Research / educational grade. Everything here is a **textbook
//! closed-form or well-established numerical model**: the one-
//! dimensional Fick laws, the explicit FTCS scheme and its standard
//! stability limit, the heat-kernel point source, and the linear
//! steady state. It is **not** a clinical, medical, or production
//! engineering tool, and it is deliberately narrow:
//!
//! - strictly **1-D**, single species, with a **constant scalar**
//!   diffusion coefficient — no concentration- or position-dependent
//!   `D`, no anisotropy, no coupling between species;
//! - **no advection, no reaction, no sources/sinks** beyond the initial
//!   condition (pure diffusion only);
//! - the transient solver is **first-order explicit** — accurate when
//!   well resolved and inside the stability bound, but with no implicit
//!   / Crank-Nicolson option and no adaptive time stepping;
//! - only **uniform grids** and only Dirichlet / zero-flux boundaries
//!   (no Robin / flux-prescribed walls).
//!
//! For multi-dimensional, reactive, advective, or implicit transport,
//! dispatch to a full PDE / CFD solver.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod analytic;
pub mod error;
pub mod explicit;
pub mod steady;

pub use analytic::{
    first_law_flux, flux_central, gaussian_point_source, gaussian_std, gaussian_variance,
    time_to_reach_std, time_to_reach_variance,
};
pub use error::{DiffusionError, ErrorCategory, Result};
pub use explicit::{is_stable_dt, max_stable_dt, Boundary, Field, Grid};
pub use steady::{steady_flux, steady_gradient, steady_profile, steady_value};
