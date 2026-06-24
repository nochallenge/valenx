//! # valenx-fluids — particle-based fluid simulation (SPH family)
//!
//! An **in-house, dependency-light** library for *particle-based* fluid
//! simulation. Fluid is represented as a cloud of Lagrangian particles; forces
//! are evaluated by **smoothed-particle hydrodynamics (SPH)** kernel sums over
//! nearby particles, found with a **uniform-grid spatial hash**. Two solvers are
//! provided:
//!
//! * [`SphSolver`] — weakly-compressible SPH after **Müller, Charypar & Gross
//!   2003** (*"Particle-Based Fluid Simulation for Interactive Applications"*):
//!   a poly6 density kernel, a spiky pressure-gradient kernel, a
//!   viscosity-Laplacian kernel; a smoothed density estimate; pressure from a
//!   Tait or linearised-ideal-gas [equation of state](EquationOfState);
//!   momentum-conserving pressure + viscosity + gravity forces; semi-implicit
//!   (symplectic) Euler integration; and a reflect/clamp [`BoxBoundary`].
//! * [`PciSphSolver`] — **predictive-corrective incompressible SPH**
//!   (**Solenthaler & Pajarola 2009**): instead of a stiff EOS it iteratively
//!   *solves* for the pressure that cancels the predicted density error,
//!   enforcing near-incompressibility at a larger stable time step. This is the
//!   one method implemented "beyond plain SPH".
//!
//! ## Clean-room provenance
//!
//! This is a **clean-room reimplementation from the published methods** (the
//! 2003 and 2009 papers above), *algorithmically inspired* by the structure of
//! Doyub Kim's `fluid-engine-dev` ("Jet", MIT) — a graphics/real-time fluid
//! engine. No source was copied; the kernels, EOS, and PCISPH δ-factor are coded
//! directly from the equations in the papers.
//!
//! ## ⚠ Honesty / scope — graphics-grade, NOT validated CFD
//!
//! **This is graphics / real-time-grade interactive fluid. It is NOT validated
//! against analytic computational fluid dynamics or experiment, and must not be
//! used for engineering CFD.** The methods here trade physical fidelity for
//! speed and stability the way a real-time graphics fluid does:
//!
//! * **Weakly compressible, not incompressible** — even PCISPH only drives the
//!   density error *small*, it does not solve a true pressure-Poisson problem.
//! * **No validated turbulence, surface tension, or thermal coupling**;
//!   viscosity is the simple Müller Laplacian model, not a calibrated stress
//!   tensor.
//! * **The boundary is a positional reflect/clamp box** — it contributes no
//!   boundary-particle pressure to the density sum, so wall behaviour (no-slip,
//!   stand-off layer, pressure at the wall) is **not** physical.
//! * **Semi-implicit Euler** is first-order; it has a known `O(dt)` integration
//!   bias (made explicit in the free-fall benchmark below) and the usual
//!   sound-speed CFL stability limit on the time step.
//!
//! What *is* pinned in the tests are the **achievable analytic checks** — the
//! kernel mathematics and the conservation laws the discretisation must obey
//! exactly — not a comparison to a Navier–Stokes reference solution.
//!
//! ## Pinned benchmarks (fail-loud, never fudged)
//!
//! 1. **Free-fall** — a single particle with no neighbours under gravity
//!    integrates as `x = ½ g t²` (within the semi-implicit-Euler `O(dt)` bias,
//!    which is itself checked in closed form). *(`sph::tests`)*
//! 2. **Kernel normalisation & support** — the poly6 kernel integrates to `≈ 1`
//!    over its support by spherical-shell quadrature, and is *exactly* `0`
//!    beyond `h`. *(`kernels::tests`)*
//! 3. **Momentum conservation** — a closed clump with gravity off conserves
//!    total linear momentum to round-off over several steps (the symmetric
//!    `p/ρ²` pressure form gives equal-and-opposite pair forces). *(`sph::tests`)*
//! 4. **Newton's third law** — two particles compressed below rest spacing are
//!    pushed apart *symmetrically*: the centre of mass and net momentum are
//!    unchanged. Checked for both [`SphSolver`] and [`PciSphSolver`].
//! 5. **Hydrostatic gradient** — a settled column is denser (higher pressure)
//!    at the bottom than the top. *(`sph::tests`)*
//!
//! ## Fail-loud configuration & guarded divides
//!
//! Every constructor validates its parameters and returns
//! [`Result<_, FluidError>`] (never a panic or a silent `NaN`): a non-positive
//! smoothing length `h`, particle mass, rest density, sound speed, or time step;
//! a degenerate box; a non-finite gravity vector or particle position are all
//! caught up front. In the solver every divide is guarded: density is floored at
//! the lone-particle self-density `m·W_poly6(0) > 0` before any `1/ρ` or `1/ρ²`,
//! the spiky gradient returns the zero vector at `r = 0` (radial direction
//! undefined — no divide by `r`), and every kernel only evaluates `r ∈ [0, h]`.
//!
//! ## Dependencies & purity
//!
//! Only [`nalgebra`](https://docs.rs/nalgebra) (vectors) and
//! [`thiserror`](https://docs.rs/thiserror) (the error type), both already in
//! the workspace — **no new external dependency**. Pure Rust, `#![forbid(unsafe_code)]`,
//! no `rand` (the solvers are deterministic).
//!
//! ## Example
//!
//! ```
//! use nalgebra::Vector3;
//! use valenx_fluids::{SphSolver, SphConfig, ParticleSystem, Particle, BoxBoundary};
//!
//! // A weakly-compressible water-like solver at smoothing length h = 0.1 m.
//! let mut cfg = SphConfig::water(0.1).unwrap();
//! cfg.gravity = Vector3::new(0.0, 0.0, -9.806_65);
//! let mut solver = SphSolver::new(cfg).unwrap();
//!
//! // A small block of particles dropped into a box.
//! let spacing = 0.05;
//! let mut system = ParticleSystem::new();
//! for ix in 0..4 {
//!     for iy in 0..4 {
//!         for iz in 0..4 {
//!             let p = Vector3::new(ix as f64, iy as f64, iz as f64) * spacing
//!                 + Vector3::new(0.1, 0.1, 0.3);
//!             system.push(Particle::at(p)).unwrap();
//!         }
//!     }
//! }
//! let tank = BoxBoundary::new(
//!     Vector3::new(0.0, 0.0, 0.0),
//!     Vector3::new(0.5, 0.5, 0.5),
//!     0.3, // wall restitution
//! ).unwrap();
//!
//! // Step the simulation; the block falls, compresses, and stays in the tank.
//! for _ in 0..50 {
//!     solver.step(&mut system, 1e-3, Some(&tank)).unwrap();
//! }
//! assert!(system.particles().iter().all(|p| tank.contains(p.position)));
//! // Densities are populated and physical (≥ the lone-particle floor).
//! assert!(system.particles().iter().all(|p| p.density >= solver.density_floor()));
//! ```

#![forbid(unsafe_code)]

pub mod boundary;
pub mod eos;
pub mod grid;
pub mod kernels;
pub mod particle;
pub mod pcisph;
pub mod sph;

mod error;

pub use boundary::BoxBoundary;
pub use eos::EquationOfState;
pub use error::FluidError;
pub use grid::SpatialHash;
pub use kernels::SmoothingKernels;
pub use particle::{Particle, ParticleSystem};
pub use pcisph::{PciSphConfig, PciSphSolver};
pub use sph::{SphConfig, SphSolver};
