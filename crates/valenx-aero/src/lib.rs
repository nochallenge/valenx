//! # valenx-aero
//!
//! A native **3-D external-aerodynamics CFD engine** — a virtual wind
//! tunnel for cars, wings, aircraft and arbitrary bodies.
//!
//! ## What this is
//!
//! Drop a 3-D triangle-mesh body into a virtual wind tunnel, set the
//! wind, run, and get the engineering answer: the **drag, lift and
//! moment coefficients**, the surface **pressure-coefficient** field,
//! and the 3-D velocity / pressure flow field plus a wake survey.
//!
//! The engine solves the steady incompressible Reynolds-averaged
//! Navier-Stokes equations
//!
//! ```text
//!   ∇·u = 0
//!   (u·∇)u = −(1/ρ)∇p + ∇·((ν + ν_t)∇u)
//! ```
//!
//! by the finite-volume method on a 3-D Cartesian Harlow-Welch
//! **staggered (MAC) grid** ([`grid`]) — pressure at cell centres, the
//! three velocity components at the cell faces, the layout that
//! defeats the checkerboard-pressure mode. The pressure-velocity
//! coupling is **SIMPLE**; the pressure-correction Poisson system is
//! solved with a **geometric-multigrid V-cycle** ([`poisson`]) because
//! a 3-D Poisson system converges far too slowly under plain SOR.
//!
//! The body is handled by the **immersed-boundary method**
//! ([`immersed`]): an arbitrary triangle mesh is *voxelized* into the
//! Cartesian grid (robust ray-cast inside/outside classification), and
//! a no-slip wall is forced on the body's cut cells. This is what lets
//! a car or a wing be dropped into the flow with **no body-fitted
//! mesher** in the loop. Turbulence — the flow is firmly turbulent at
//! road / flight Reynolds numbers — is closed with a **k-ε** or
//! **k-ω SST** RANS model ([`turbulence`]).
//!
//! ## What you do with it
//!
//! ```no_run
//! use valenx_aero::{geometry, run_windtunnel, AeroRequest, TurbulenceModel};
//!
//! // A simple bluff body — a 4×2×1.5 m box "car".
//! let body = geometry::box_body(
//!     nalgebra::Vector3::new(0.0, 0.0, 0.0),
//!     nalgebra::Vector3::new(4.0, 2.0, 1.5),
//! );
//! // Run it at 30 m/s in sea-level air.
//! let request = AeroRequest::new(30.0).with_turbulence(TurbulenceModel::KEpsilon);
//! let result = run_windtunnel(&body, &request).expect("valid case");
//! // The drag coefficient is the headline number.
//! let _cd = result.coefficients.cd;
//! ```
//!
//! (The example is marked `no_run`: a full wind-tunnel solve at the
//! default resolution takes seconds — it is type-checked here and
//! exercised end-to-end by the crate's test suite.)
//!
//! ## Honest scope — a real v1, not Fluent / STAR-CCM+
//!
//! Every algorithm here is the genuine article and the engine produces
//! a physically meaningful external-aero result — the immersed-
//! boundary tests verify a sphere / cube drag coefficient lands in the
//! right ballpark for the Reynolds regime, an empty tunnel preserves
//! the free-stream, and the multigrid recovers a known Poisson
//! solution. It is deliberately a **v1**, **not** ANSYS Fluent /
//! Siemens STAR-CCM+ parity:
//!
//! - **Immersed-boundary, not body-fitted meshing.** The body is
//!   voxelized into a uniform Cartesian grid; the wall inside a cut
//!   cell is the true sub-cell clipped surface ([`cutcell`]), and a
//!   law-of-the-wall **near-wall model** ([`wallmodel`]) reconstructs
//!   the turbulent boundary-layer profile so the wall shear and the
//!   surface forces are accurate without the grid resolving the
//!   boundary layer. There is still no body-fitted unstructured mesher
//!   and no near-wall prism layer — a production wind-tunnel run uses
//!   a body-fitted near-wall mesh.
//! - **RANS only.** The turbulence models are steady k-ε and k-ω SST
//!   with a high-Reynolds wall function. There is **no DES / LES** —
//!   the engine resolves the *mean* flow, not the unsteady eddies.
//! - **Incompressible.** Valid up to roughly Mach 0.3; a
//!   Prandtl-Glauert compressibility *correction* ([`compressible`])
//!   extends the coefficient estimate into the subsonic range, but
//!   there is no compressible / transonic solver.
//! - **Uniform Cartesian grid**, single-phase, constant properties,
//!   the hybrid (not a higher-order TVD) convection scheme.
//!
//! None of those omissions makes the result meaningless — the trends,
//! the pressure field, the wake and the coefficient *ballpark* are
//! all real. Each is a documented, well-understood extension. For
//! production-accuracy automotive / aerospace CFD, dispatch to the
//! OpenFOAM or SU2 subprocess adapters.
//!
//! ## Module map
//!
//! - [`grid`] — the 3-D Cartesian staggered (MAC) grid + field
//!   storage.
//! - [`poisson`] — the pressure Poisson solver: SOR + geometric
//!   multigrid.
//! - [`geometry`] — the body triangle mesh + extraction from
//!   `valenx-mesh` / `valenx-cad`.
//! - [`immersed`] — the immersed-boundary voxelizer.
//! - [`wind`] — the free-stream wind specification + air properties.
//! - [`domain`] — the virtual-wind-tunnel domain builder + boundary
//!   conditions.
//! - [`turbulence`] — k-ε and k-ω SST RANS turbulence models.
//! - [`solver`] — the 3-D incompressible Navier-Stokes SIMPLE solver.
//! - [`forces`] — surface-force integration + the aerodynamic
//!   coefficients.
//! - [`postprocess`] — wake survey, streamlines, vorticity,
//!   Q-criterion, cut planes.
//! - [`transient`] — the unsteady solver for wake shedding.
//! - [`compressible`] — the Prandtl-Glauert compressibility
//!   correction.
//! - [`wallmodel`] — the near-wall model: a law-of-the-wall
//!   reconstruction of the turbulent boundary-layer profile so the
//!   wall shear, the momentum sink and the surface forces are accurate
//!   without resolving the boundary layer with the grid.
//! - [`refine`] — local grid refinement near the body and in the
//!   wake.
//! - [`benchmark`] — the published-reference validation suite (sphere
//!   drag curve, flat-plate skin friction, NACA-airfoil lift / drag).
//! - [`presets`] — automotive / aircraft case presets.
//! - [`sweep`] — the angle-of-attack sweep → a lift / drag polar.
//! - [`report`] — the human-readable run summary.
//! - [`api`] — the typed request / response surface for an external
//!   agent.
//! - [`error`] — the crate error taxonomy.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod api;
pub mod benchmark;
pub mod compressible;
pub mod cutcell;
pub mod domain;
pub mod error;
pub mod forces;
pub mod geometry;
pub mod grid;
pub mod immersed;
pub mod poisson;
pub mod postprocess;
pub mod presets;
pub mod refine;
pub mod report;
pub mod solver;
pub mod sweep;
pub mod transient;
pub mod turbulence;
pub mod wallmodel;
pub mod wind;

pub use api::{run_windtunnel, AeroRequest, AeroResult};
pub use benchmark::{
    run_flat_plate, run_naca_airfoil, run_sphere_drag, AirfoilResult,
    FlatPlateResult, SphereDragPoint,
};
pub use domain::{BoundaryConditions, FaceBc, TunnelSizing, WindTunnel};
pub use error::{AeroError, ErrorCategory};
pub use forces::{
    coefficients, integrate_forces, integrate_forces_with, surface_field,
    surface_stats, AeroCoefficients, AeroForces, SurfacePoint, SurfaceStats,
    WindFrame,
};
pub use cutcell::{CellGeometry, CutFace, WallMethod};
pub use geometry::{
    box_body, naca4_half_thickness, naca_wing, sphere_body, Aabb, TriMesh, Triangle,
};
pub use grid::{Field3, Grid3};
pub use immersed::{voxelize, voxelize_with, CellTag, ImmersedBody};
pub use poisson::{solve_multigrid, solve_sor, PoissonStencil};
pub use postprocess::{
    q_criterion, slice_field, trace_streamline, vorticity, wake_survey, FieldSlice,
    SliceAxis, Streamline, VorticityField, WakeSurvey,
};
pub use report::AeroReport;
pub use solver::{solve_steady, BodyMotion, FlowField, SolverControls};
pub use sweep::{aoa_sweep, PolarCurve, PolarPoint};
pub use transient::{solve_transient, TransientControls, TransientHistory};
pub use turbulence::{TurbulenceModel, TurbulenceState};
pub use wallmodel::{
    friction_velocity, wall_effective_viscosity, wall_shear_stress, y_plus,
};
pub use wind::{Air, Wind};
