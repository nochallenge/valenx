//! # valenx-cfd-native
//!
//! A native, in-process **2-D laminar computational-fluid-dynamics
//! solver** — steady, incompressible Navier-Stokes by the SIMPLE
//! algorithm on a staggered finite-volume grid.
//!
//! ## What this is
//!
//! A genuine, self-contained laminar-flow solver. It discretises the
//! steady incompressible Navier-Stokes equations
//!
//! ```text
//!   ∇·u = 0
//!   (u·∇)u = −(1/ρ)∇p + ν∇²u
//! ```
//!
//! by the finite-volume method on a Harlow-Welch **staggered grid**
//! ([`Grid`]) — pressure at cell centres, the two velocity components
//! at the cell faces, the layout that defeats the checkerboard-pressure
//! mode. The convective flux uses the unconditionally-stable hybrid
//! scheme; the pressure-velocity coupling is **SIMPLE** (Patankar &
//! Spalding 1972): a momentum predictor, a pressure-correction Poisson
//! solve ([`crate::linsolve`]), and an under-relaxed correction,
//! iterated to a divergence-free steady field.
//!
//! It is the in-process counterpart to dispatching to the OpenFOAM /
//! SU2 subprocess adapters: no external binary, no case directory, runs
//! anywhere the crate compiles.
//!
//! ## What you do with it
//!
//! 1. Build a [`Grid`] over the rectangular domain.
//! 2. Pick a [`Fluid`] (density + kinematic viscosity).
//! 3. Set the four-sided [`Boundaries`] — the
//!    [`Boundaries::lid_driven_cavity`] and [`Boundaries::channel_flow`]
//!    constructors cover the two canonical cases.
//! 4. Call [`solve_simple`]; read the velocity / pressure fields off
//!    the returned [`FlowSolution`].
//!
//! ```
//! use valenx_cfd_native::{Boundaries, Fluid, Grid, SimpleControls, solve_simple};
//!
//! // A 20×20 lid-driven cavity at a modest Reynolds number.
//! let grid = Grid::new(20, 20, 1.0, 1.0);
//! let fluid = Fluid::new(1.0, 0.1);            // Re = U·L/ν = 10
//! let bcs = Boundaries::lid_driven_cavity(1.0); // lid slides at 1 m/s
//! let solution = solve_simple(&grid, &fluid, &bcs, &SimpleControls::default());
//!
//! // The fluid under the lid is dragged along with it.
//! assert!(solution.u_at_cell(10, 18) > 0.0);
//! ```
//!
//! ## Turbulence and unsteady flow
//!
//! Beyond the steady laminar core, three further subsystems ship —
//! enough to take the solver into the turbulent and the time-dependent
//! regimes:
//!
//! - **Turbulence — the standard k-ε and Menter k-ω SST models**
//!   ([`turbulence`]). The two canonical RANS closures. **k-ε**
//!   (Launder & Spalding 1974): two transport equations for the
//!   turbulent kinetic energy `k` and its dissipation `ε`, eddy
//!   viscosity `μ_t = ρ·C_μ·k²/ε`. **k-ω SST** (Menter 1994): two
//!   transport equations for `k` and the specific dissipation `ω`,
//!   the `F₁` / `F₂` blending functions that switch the model from
//!   k-ω near a wall to k-ε in the free stream, and the SST-limited
//!   eddy viscosity `μ_t = ρ·a₁·k/max(a₁·ω, S·F₂)`. Both use
//!   equilibrium log-law wall functions. The SIMPLE driver consumes
//!   either through the [`EffectiveViscosity`] selector, so a
//!   turbulent run is the laminar call with one extra argument.
//! - **Transient (unsteady) flow** ([`transient`]). The
//!   time-derivative term put back into the momentum equation and an
//!   implicit (backward-Euler) **time-marching loop** — a
//!   transient-SIMPLE solver. A flow started from rest relaxes onto
//!   the steady solution; the field is divergence-free at every step.
//! - **Geometric multigrid for the pressure-Poisson solve**
//!   ([`multigrid`]). A V-cycle on the staggered grid with a
//!   weighted-Jacobi smoother, full-weighting restriction,
//!   bilinear-interpolation prolongation, and a Galerkin coarse
//!   operator — the convergence-rate-per-cycle is essentially
//!   grid-independent (SOR alone scales poorly on fine 2-D grids).
//!   Wired into the SIMPLE driver through [`PressurePoissonSolver`]
//!   alongside the original SOR.
//!
//! ## Honest scope — a real v1, not OpenFOAM
//!
//! Every algorithm here is the genuine article and the solver produces
//! the physically correct flow — the lid-driven-cavity recirculation,
//! the developing-channel parabolic profile, the flatter turbulent
//! profile, and the transient relaxation onto steady state are all
//! verified in the per-module tests. It is deliberately a **v1**:
//!
//! - **2-D only.** The staggered grid and the discretisation are
//!   two-dimensional; 3-D is the same method with a third momentum
//!   equation and is a substantial but well-understood extension.
//! - **Turbulence: the standard k-ε and Menter k-ω SST models.** No
//!   realizable / RNG k-ε, no Reynolds-stress model, no LES / DES.
//!   The laminar solver remains the right tool while the Reynolds
//!   number stays modest.
//! - **Transient: first-order implicit Euler, transient SIMPLE.** A
//!   second-order (Crank-Nicolson / BDF2) scheme or a non-iterative
//!   PISO loop are accuracy / efficiency follow-ups.
//! - **Structured uniform rectangular grid**, single-phase, constant
//!   properties, hybrid (not higher-order TVD) convection scheme.
//!
//! None of those omissions affects the correctness of the flows it
//! does solve; each is a documented, well-understood extension. For
//! 3-D / complex-geometry / multiphase work, use Valenx's OpenFOAM or
//! SU2 subprocess adapters.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod benchmark;
pub mod error;
pub mod grid;
pub mod linsolve;
pub mod multigrid;
pub mod solver;
pub mod transient;
pub mod turbulence;

pub use benchmark::{
    backward_facing_step_reattachment, compare_to_ghia_cavity, poiseuille_centerline_check,
    sample_u, sample_v, BackwardStepResult, GhiaError, PoiseuilleError, GHIA_U_RE_100,
    GHIA_U_RE_1000, GHIA_U_RE_400, GHIA_V_RE_100, GHIA_V_RE_1000, GHIA_V_RE_400, GHIA_X, GHIA_Y,
};
pub use error::CfdError;
pub use grid::{Field, Grid};
pub use linsolve::{poisson_residual, solve_sor, PoissonCoeffs, SorResult};
pub use multigrid::{
    coarsen_coefficients, prolong_bilinear, restrict_full_weighting, solve_multigrid,
    solve_pressure_poisson, v_cycle, weighted_jacobi_sweep, MultigridControls, MultigridResult,
    PressurePoissonSolver,
};
pub use solver::{
    solve_simple, solve_simple_with, Boundaries, EffectiveViscosity, FlowSolution, Fluid, SideBc,
    SimpleControls, TurbulenceSnapshot,
};
pub use transient::{solve_transient, TransientControls, TransientSolution, TransientStep};
pub use turbulence::{
    advance_k_epsilon, advance_k_omega_sst, f1_blend, f2_blend, solve_turbulent_channel,
    solve_turbulent_channel_sst, strain_rate_magnitude, strain_rate_squared, wall_distance_field,
    ChannelProfile, KEpsilonModel, SstField, SstModel, SstSet, TurbulenceField, WallFunction,
    WallMask,
};
