//! The 3-D incompressible Navier-Stokes solver — the SIMPLE algorithm
//! on the Cartesian staggered grid, with immersed-boundary wall
//! forcing and an optional RANS turbulence model.
//!
//! # The equations
//!
//! Steady incompressible RANS over the 3-D tunnel domain:
//!
//! ```text
//!   ∇·u = 0                                              (continuity)
//!   (u·∇)u = −(1/ρ)∇p + ∇·((ν + ν_t)∇u)                  (momentum)
//! ```
//!
//! discretised by the finite-volume method on the staggered grid of
//! [`crate::grid`]. Each velocity component's momentum equation is
//! integrated over its own staggered control volume; the convective
//! flux uses the unconditionally-stable **hybrid (Spalding) scheme**,
//! the diffusive flux is central, and the effective viscosity adds the
//! eddy viscosity `ν_t` from the turbulence model
//! ([`crate::turbulence`]) to the molecular `ν`.
//!
//! # SIMPLE in 3-D
//!
//! Identical in structure to the 2-D SIMPLE of `valenx-cfd-native`,
//! with a third momentum equation and a seven-point pressure-
//! correction Poisson stencil:
//!
//! 1. with the current pressure, solve the three momentum equations
//!    for a provisional velocity `u*` (satisfies momentum, not
//!    continuity);
//! 2. assemble the pressure-correction Poisson equation whose source
//!    is each cell's mass imbalance and solve it with multigrid
//!    ([`crate::poisson`]);
//! 3. correct the pressure and the three velocity components so the
//!    corrected velocity is divergence-free;
//! 4. advance the turbulence fields, re-impose the boundary and
//!    immersed-boundary conditions, under-relax, and repeat.
//!
//! # The immersed boundary
//!
//! A cell tagged `Solid` by [`crate::immersed`] is excluded from the
//! fluid solve — its velocity is frozen at the body velocity (zero for
//! a static body, the wheel-tangential velocity for a rotating wheel).
//! A face between a fluid and a solid cell is a wall: the velocity
//! there is driven to the wall velocity (direct-forcing IBM).
//!
//! Under the default **cut-cell** method ([`crate::cutcell`]) a `Cut`
//! cell carries genuine sub-cell geometry — a fluid volume fraction,
//! apertured Cartesian face areas, and a true clipped wall face. The
//! finite-volume discretisation uses it: the **continuity / pressure-
//! correction** equation scales every Cartesian face flux by the
//! face's open aperture (a closed-off face carries no mass), and the
//! **momentum** equation adds the clipped wall face as a no-slip
//! shear-drag term that pulls the cut cell's tangential velocity toward
//! the wall. A genuinely tiny cut cell is merged into the solid at
//! voxelisation time and the kept cut cells are floored
//! ([`crate::immersed::SMALL_CELL_FLOOR`]) so the segregated solver
//! stays stable.
//!
//! # Honest scope
//!
//! A real working 3-D incompressible solver. It is a v1: the
//! convection scheme is hybrid (not a higher-order TVD scheme), the
//! pressure-velocity coupling is SIMPLE (not the faster SIMPLEC /
//! PISO), the turbulence wall treatment is a high-Re wall function, and
//! the cut-cell method runs on a uniform Cartesian background grid with
//! no near-wall prism layer. Each is a documented, well-understood
//! extension.

use crate::cutcell::WallMethod;
use crate::domain::{BoundaryConditions, FaceBc, WindTunnel};
use crate::grid::{Field3, Grid3};
use crate::immersed::{stabilized_volume, CellTag, ImmersedBody};
use crate::poisson::{solve_multigrid_anchored, PoissonStencil};
use crate::turbulence::{advance_turbulence, TurbulenceModel, TurbulenceState};

/// Solver controls — under-relaxation, iteration limits, convergence.
#[derive(Clone, Copy, Debug)]
pub struct SolverControls {
    /// Momentum under-relaxation factor `αu` (`0 < αu ≤ 1`).
    pub relax_u: f64,
    /// Pressure under-relaxation factor `αp` (`0 < αp ≤ 1`).
    pub relax_p: f64,
    /// Turbulence under-relaxation / pseudo-time step.
    pub relax_turb: f64,
    /// Maximum number of SIMPLE outer iterations.
    pub max_iterations: usize,
    /// Convergence tolerance on the scaled mass-imbalance residual.
    pub tolerance: f64,
    /// Maximum multigrid V-cycles per pressure solve.
    pub mg_cycles: usize,
    /// The turbulence model to run.
    pub turbulence: TurbulenceModel,
    /// Whether the near-wall model ([`crate::wallmodel`]) drives the
    /// cut-cell wall-drag term.
    ///
    /// With it `true` (the default) the no-slip wall-drag a cut cell
    /// exerts on the momentum equation is reconstructed from the
    /// turbulent law of the wall — the near-wall momentum sink, hence
    /// the separation point and the pressure drag, are physically
    /// correct on a grid that does not resolve the boundary layer.
    /// Set it `false` for the legacy crude-linear-gradient wall drag
    /// (`μ·u/d`) — kept for an explicit before/after comparison; the
    /// near-wall model is the accurate path and the wind-tunnel default.
    pub near_wall_model: bool,
}

impl Default for SolverControls {
    /// Conservative, broadly-stable defaults: `αu = 0.5`, `αp = 0.2`,
    /// up to 400 outer iterations to a `1e-4` residual, k-ω SST, the
    /// near-wall model on.
    fn default() -> Self {
        SolverControls {
            relax_u: 0.5,
            relax_p: 0.2,
            relax_turb: 0.3,
            max_iterations: 400,
            tolerance: 1e-4,
            mg_cycles: 12,
            turbulence: TurbulenceModel::KOmegaSST,
            near_wall_model: true,
        }
    }
}

/// The solved (or iteration-capped) 3-D flow field.
#[derive(Clone, Debug)]
pub struct FlowField {
    /// The grid the field lives on.
    pub grid: Grid3,
    /// Staggered `u`-velocity, `(nx+1)·ny·nz`.
    pub u: Field3,
    /// Staggered `v`-velocity, `nx·(ny+1)·nz`.
    pub v: Field3,
    /// Staggered `w`-velocity, `nx·ny·(nz+1)`.
    pub w: Field3,
    /// Cell-centred pressure, `nx·ny·nz` (gauge: pinned at the outlet).
    pub pressure: Field3,
    /// The turbulence state.
    pub turbulence: TurbulenceState,
    /// Number of SIMPLE outer iterations performed.
    pub iterations: usize,
    /// Final scaled mass-imbalance residual.
    pub residual: f64,
    /// Per-iteration residual history.
    pub residual_history: Vec<f64>,
    /// True if the residual reached [`SolverControls::tolerance`].
    pub converged: bool,
}

impl FlowField {
    /// Cell-centred `u`-velocity of pressure cell `(i, j, k)` — the
    /// average of the two bracketing `u`-faces.
    pub fn u_at_cell(&self, i: usize, j: usize, k: usize) -> f64 {
        0.5 * (self.u.at(i, j, k) + self.u.at(i + 1, j, k))
    }

    /// Cell-centred `v`-velocity of pressure cell `(i, j, k)`.
    pub fn v_at_cell(&self, i: usize, j: usize, k: usize) -> f64 {
        0.5 * (self.v.at(i, j, k) + self.v.at(i, j + 1, k))
    }

    /// Cell-centred `w`-velocity of pressure cell `(i, j, k)`.
    pub fn w_at_cell(&self, i: usize, j: usize, k: usize) -> f64 {
        0.5 * (self.w.at(i, j, k) + self.w.at(i, j, k + 1))
    }

    /// Velocity magnitude at the centre of pressure cell `(i, j, k)`.
    pub fn speed_at_cell(&self, i: usize, j: usize, k: usize) -> f64 {
        let u = self.u_at_cell(i, j, k);
        let v = self.v_at_cell(i, j, k);
        let w = self.w_at_cell(i, j, k);
        (u * u + v * v + w * w).sqrt()
    }

    /// All three cell-centred velocity components as a `Field3` triple
    /// — convenient for post-processing (vorticity, turbulence
    /// production, streamlines all want cell-centred velocity).
    pub fn cell_centred_velocity(&self) -> (Field3, Field3, Field3) {
        let g = &self.grid;
        let mut uc = g.scalar_field();
        let mut vc = g.scalar_field();
        let mut wc = g.scalar_field();
        for k in 0..g.nz {
            for j in 0..g.ny {
                for i in 0..g.nx {
                    uc.set(i, j, k, self.u_at_cell(i, j, k));
                    vc.set(i, j, k, self.v_at_cell(i, j, k));
                    wc.set(i, j, k, self.w_at_cell(i, j, k));
                }
            }
        }
        (uc, vc, wc)
    }
}

/// Per-cell body velocity override — used so a rotating wheel can spin
/// the solid cells in its region. The default (`None`) means a static
/// body (zero velocity).
#[derive(Clone, Debug, Default)]
pub struct BodyMotion {
    /// Optional per-cell `(u, v, w)` velocity for solid cells. Empty
    /// = a static body.
    pub solid_velocity: Vec<[f64; 3]>,
}

impl BodyMotion {
    /// A static body — every solid cell at rest.
    pub fn static_body() -> BodyMotion {
        BodyMotion::default()
    }

    /// The body velocity at a solid cell index, or zero if no motion
    /// is specified.
    #[inline]
    fn at(&self, idx: usize) -> [f64; 3] {
        if idx < self.solid_velocity.len() {
            self.solid_velocity[idx]
        } else {
            [0.0, 0.0, 0.0]
        }
    }
}

/// Solve the steady 3-D flow over the body in the wind tunnel.
///
/// `tunnel` carries the grid, the voxelized body, the wind and the
/// boundary set; `controls` tunes the iteration; `motion` optionally
/// spins solid cells (a rotating wheel). Returns the [`FlowField`] —
/// the velocity / pressure / turbulence fields plus the iteration and
/// residual diagnostics. The solver always returns a field, even if it
/// hits the iteration cap before the tolerance.
pub fn solve_steady(
    tunnel: &WindTunnel,
    controls: &SolverControls,
    motion: &BodyMotion,
) -> FlowField {
    let grid = tunnel.grid;
    let rho = tunnel.wind.air.density;
    let mu_lam = tunnel.wind.air.dynamic_viscosity;
    let free = tunnel.free_stream();

    // A boolean solid mask, indexed like a cell-centred field.
    let solid: Vec<bool> = tunnel
        .body
        .tags
        .iter()
        .map(|&t| t == CellTag::Solid)
        .collect();

    // Fields.
    let mut u = grid.u_field();
    let mut v = grid.v_field();
    let mut w = grid.w_field();
    let mut p = grid.scalar_field();

    // Initialise the interior with the free-stream — a good guess that
    // speeds convergence and keeps the early iterations stable.
    u.fill(free.x);
    v.fill(free.y);
    w.fill(free.z);

    // Turbulence state.
    let length_scale = 0.07 * tunnel.reference_length.max(1e-6);
    let mut turb = match controls.turbulence {
        TurbulenceModel::Laminar => TurbulenceState::laminar(&grid),
        model => TurbulenceState::uniform(
            &grid,
            model,
            tunnel.wind.inlet_tke().max(1e-6),
            match model {
                TurbulenceModel::KEpsilon => tunnel.wind.inlet_epsilon(length_scale).max(1e-8),
                _ => tunnel.wind.inlet_omega(length_scale).max(1e-3),
            },
            rho,
        ),
    };

    // Per-face momentum diagonals — reused for the pressure correction.
    let mut apu = grid.u_field();
    let mut apv = grid.v_field();
    let mut apw = grid.w_field();

    apply_boundary_velocities(&mut u, &mut v, &mut w, &grid, &tunnel.bc, free);
    apply_immersed_velocities(&mut u, &mut v, &mut w, &grid, &solid, motion);

    let mut residual = f64::INFINITY;
    let mut residual_history = Vec::new();
    let mut converged = false;
    let mut iterations = 0;

    // Residual normalisation — a characteristic mass flow.
    let ref_speed = tunnel.wind.speed.max(1e-6);
    let cell_face = grid.dy() * grid.dz();
    let mass_scale = (rho * ref_speed * cell_face).max(1e-30);

    // The near-wall model reconstructs a steep boundary-layer profile
    // in the wall-adjacent cells, which stiffens the segregated
    // pressure-velocity coupling; it runs at a more conservative
    // under-relaxation so the SIMPLE outer iteration stays stable.
    let (relax_u, relax_p) = if controls.near_wall_model {
        (controls.relax_u.min(0.3), controls.relax_p.min(0.1))
    } else {
        (controls.relax_u, controls.relax_p)
    };

    for outer in 0..controls.max_iterations.max(1) {
        iterations = outer + 1;

        // --- (1) momentum predictor ---
        solve_momentum_u(
            &mut u,
            &v,
            &w,
            &p,
            &mut apu,
            &grid,
            rho,
            mu_lam,
            &turb,
            &solid,
            &tunnel.body,
            &tunnel.bc,
            free,
            relax_u,
            controls.near_wall_model,
        );
        solve_momentum_v(
            &u,
            &mut v,
            &w,
            &p,
            &mut apv,
            &grid,
            rho,
            mu_lam,
            &turb,
            &solid,
            &tunnel.body,
            &tunnel.bc,
            free,
            relax_u,
            controls.near_wall_model,
        );
        solve_momentum_w(
            &u,
            &v,
            &mut w,
            &p,
            &mut apw,
            &grid,
            rho,
            mu_lam,
            &turb,
            &solid,
            &tunnel.body,
            &tunnel.bc,
            free,
            relax_u,
            controls.near_wall_model,
        );
        apply_boundary_velocities(&mut u, &mut v, &mut w, &grid, &tunnel.bc, free);
        apply_immersed_velocities(&mut u, &mut v, &mut w, &grid, &solid, motion);

        // --- (2) pressure-correction Poisson equation ---
        let mut stencil = PoissonStencil::zeros(grid.nx, grid.ny, grid.nz);
        let mass_imbalance = assemble_pressure_correction(
            &mut stencil,
            &u,
            &v,
            &w,
            &apu,
            &apv,
            &apw,
            &grid,
            rho,
            &solid,
            &tunnel.body,
            &tunnel.bc,
        );

        // --- (3) solve for p' ---
        // The outlet anchor is a genuine Dirichlet pin and is threaded
        // through every multigrid level, so the operator is non-singular
        // on the coarse grids too (a purely-singular coarse Neumann
        // operator cannot consistently solve the restricted residual,
        // which carries a net boundary mass imbalance).
        let mut pcorr = grid.scalar_field();
        solve_multigrid_anchored(
            &stencil,
            &mut pcorr,
            controls.tolerance * 1e-2,
            controls.mg_cycles,
            2,
            2,
            false,
            outlet_anchor(&grid, &tunnel.bc),
        );

        // --- (4) correct pressure + velocity ---
        correct_pressure(&mut p, &pcorr, &solid, relax_p);
        correct_velocity(
            &mut u,
            &mut v,
            &mut w,
            &pcorr,
            &apu,
            &apv,
            &apw,
            &grid,
            &solid,
            &tunnel.body,
        );
        apply_boundary_velocities(&mut u, &mut v, &mut w, &grid, &tunnel.bc, free);
        apply_immersed_velocities(&mut u, &mut v, &mut w, &grid, &solid, motion);

        // --- (5) advance turbulence ---
        if controls.turbulence != TurbulenceModel::Laminar {
            let (uc, vc, wc) = cell_centre_velocity(&u, &v, &w, &grid);
            let wall_d = wall_distance_field(&tunnel.body);
            // With the near-wall model the cut cells take the
            // wall-function equilibrium turbulence; otherwise the
            // turbulence transport runs free everywhere.
            let wall = if controls.near_wall_model {
                Some(&tunnel.body)
            } else {
                None
            };
            advance_turbulence(
                &mut turb,
                &uc,
                &vc,
                &wc,
                &wall_d,
                &solid,
                &grid,
                rho,
                mu_lam,
                controls.relax_turb * 1e-2,
                wall,
            );
        }

        residual = mass_imbalance / mass_scale;
        residual_history.push(residual);
        if residual.is_finite() && residual <= controls.tolerance {
            converged = true;
            break;
        }
        if !residual.is_finite() || residual > 1e14 {
            break;
        }
    }

    FlowField {
        grid,
        u,
        v,
        w,
        pressure: p,
        turbulence: turb,
        iterations,
        residual,
        residual_history,
        converged,
    }
}

/// Build a cell-centred wall-distance field from the immersed body.
fn wall_distance_field(body: &crate::immersed::ImmersedBody) -> Field3 {
    let g = body.grid;
    let mut f = g.scalar_field();
    for k in 0..g.nz {
        for j in 0..g.ny {
            for i in 0..g.nx {
                let d = body.wall_distance(i, j, k);
                f.set(i, j, k, if d.is_finite() { d } else { 1e9 });
            }
        }
    }
    f
}

/// Average the staggered velocity faces onto the pressure-cell centres.
fn cell_centre_velocity(
    u: &Field3,
    v: &Field3,
    w: &Field3,
    grid: &Grid3,
) -> (Field3, Field3, Field3) {
    let mut uc = grid.scalar_field();
    let mut vc = grid.scalar_field();
    let mut wc = grid.scalar_field();
    for k in 0..grid.nz {
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                uc.set(i, j, k, 0.5 * (u.at(i, j, k) + u.at(i + 1, j, k)));
                vc.set(i, j, k, 0.5 * (v.at(i, j, k) + v.at(i, j + 1, k)));
                wc.set(i, j, k, 0.5 * (w.at(i, j, k) + w.at(i, j, k + 1)));
            }
        }
    }
    (uc, vc, wc)
}

/// The hybrid-scheme convection-diffusion neighbour coefficient.
#[inline]
fn hybrid_coeff(d: f64, f: f64, upwind_sign: f64) -> f64 {
    (d - 0.5 * f.abs()).max(0.0) + (upwind_sign * f).max(0.0)
}

/// The effective (molecular + eddy) kinematic viscosity at a cell.
#[inline]
fn nu_eff(turb: &TurbulenceState, mu_lam: f64, rho: f64, i: usize, j: usize, k: usize) -> f64 {
    (mu_lam + turb.mu_t.at(i, j, k)) / rho
}

/// The open (fluid) fraction of one Cartesian face of pressure cell
/// `(i, j, k)` — `1.0` everywhere under the staircased method, or the
/// cut-cell aperture under the cut-cell method. `face` is `0..6` in the
/// order `[-x, +x, -y, +y, -z, +z]`. A fully-fluid cell always reports
/// `1.0`; a solid cell reports `0.0`.
#[inline]
fn face_aperture(body: &ImmersedBody, i: usize, j: usize, k: usize, face: usize) -> f64 {
    match body.tag(i, j, k) {
        CellTag::Fluid => 1.0,
        CellTag::Solid => 0.0,
        CellTag::Cut => {
            if body.method == WallMethod::CutCell {
                body.geometry(i, j, k).face_aperture[face]
            } else {
                1.0
            }
        }
    }
}

/// The shared open aperture of the Cartesian face *between* two
/// pressure cells. The two cells each see this face — `cell_a`'s `+`
/// face and `cell_b`'s `−` face on the same axis — so the consistent
/// shared aperture is the *minimum* of the two (the flow is throttled
/// by whichever side is more closed off). Returns `1.0` under the
/// staircased method.
#[allow(clippy::too_many_arguments)]
#[inline]
fn shared_aperture(
    body: &ImmersedBody,
    ai: usize,
    aj: usize,
    ak: usize,
    a_face: usize,
    bi: usize,
    bj: usize,
    bk: usize,
    b_face: usize,
) -> f64 {
    if body.method != WallMethod::CutCell {
        return 1.0;
    }
    face_aperture(body, ai, aj, ak, a_face).min(face_aperture(body, bi, bj, bk, b_face))
}

/// The smallest fluid volume fraction among the cells bordering a
/// velocity face — `1.0` if neither bordering cell is a cut cell. Used
/// to size the small-cell `aP` floor.
#[inline]
fn min_border_volume_fraction(
    body: &ImmersedBody,
    ai: usize,
    aj: usize,
    ak: usize,
    bi: usize,
    bj: usize,
    bk: usize,
) -> f64 {
    if body.method != WallMethod::CutCell {
        return 1.0;
    }
    let fa = if body.tag(ai, aj, ak) == CellTag::Cut {
        body.volume_fraction(ai, aj, ak)
    } else {
        1.0
    };
    let fb = if body.tag(bi, bj, bk) == CellTag::Cut {
        body.volume_fraction(bi, bj, bk)
    } else {
        1.0
    };
    fa.min(fb)
}

/// Small-cell stabilisation of a momentum diagonal `aP`.
///
/// If a velocity face borders a (kept) cut cell, its momentum diagonal
/// is floored at a volume-fraction-scaled share of the wall-normal
/// diffusive conductance, so a small cut cell cannot produce a
/// near-zero `aP` and hence a runaway velocity correction
/// `d = A/aP` in the pressure-correction step. A face with no cut-cell
/// neighbour is returned unchanged — the floor must not alter the
/// hybrid-scheme `aP` of an ordinary fluid face. This is the
/// flux-redistribution-style floor that pairs with the cell-merging of
/// genuinely tiny cells at voxelisation time.
#[allow(clippy::too_many_arguments)]
#[inline]
fn small_cell_floor_ap(
    body: &ImmersedBody,
    ai: usize,
    aj: usize,
    ak: usize,
    bi: usize,
    bj: usize,
    bk: usize,
    ap: f64,
    diff_conductance: f64,
) -> f64 {
    if body.method != WallMethod::CutCell {
        return ap;
    }
    let border_cut = body.tag(ai, aj, ak) == CellTag::Cut || body.tag(bi, bj, bk) == CellTag::Cut;
    if !border_cut {
        return ap;
    }
    let theta = min_border_volume_fraction(body, ai, aj, ak, bi, bj, bk);
    let floor = stabilized_volume(theta) * diff_conductance;
    ap.max(floor)
}

/// The wall-normal distance from a cut cell's centre to its clipped
/// wall face — the boundary-layer sampling height the near-wall model
/// uses.
///
/// The cut face carries an exact centroid and outward normal; the
/// signed distance of the cell centre off that wall plane is the proper
/// wall-normal sampling height (better than a blanket half-cell). It is
/// floored at a small fraction of the cell so a cell whose centroid
/// almost touches the wall does not produce a singular `1/y`.
#[inline]
fn cut_wall_distance(body: &ImmersedBody, i: usize, j: usize, k: usize) -> f64 {
    let grid = body.grid;
    let cut = body.cut_face(i, j, k);
    let (cx, cy, cz) = grid.cell_centre(i, j, k);
    let centre = nalgebra::Vector3::new(cx, cy, cz);
    // Distance of the cell centre off the wall plane along the outward
    // normal — the wall points from solid into fluid, so a fluid cell
    // centre sits a positive distance out.
    let d = (centre - cut.centroid).dot(&cut.normal).abs();
    let cell = (grid.dx() * grid.dx() + grid.dy() * grid.dy() + grid.dz() * grid.dz()).sqrt()
        / 3.0_f64.sqrt();
    // Floor at 25 % of a cell, cap at one cell — the centroid-based
    // distance is the BL sampling height, kept in a sane band. The
    // `is_finite` guard rejects a degenerate cut face (a NaN distance
    // would otherwise pass straight through `clamp`).
    if d.is_finite() && d > 0.0 {
        d.clamp(0.25 * cell, cell)
    } else {
        0.5 * cell
    }
}

/// The cut-cell wall-drag conductance a clipped wall face adds to the
/// momentum balance of a velocity face that borders cut cell
/// `(i, j, k)`.
///
/// The clipped wall carries no-slip and removes tangential momentum at
/// the rate `μ_w · A_wall,tan / y`, where `A_wall,tan` is the wall
/// area's projection that opposes the velocity component being solved,
/// `y` is the wall-normal distance from the cell centre, and `μ_w` is
/// the near-wall momentum-transport coefficient.
///
/// With `use_wall_model` set (the default) `μ_w` is the **wall-model
/// effective viscosity** ([`crate::wallmodel`]) — the turbulent
/// transport the law of the wall implies for the cell's tangential
/// speed `u_tan_cell`. That makes the near-wall momentum sink match the
/// turbulent boundary layer instead of a crude laminar gradient, which
/// is what corrects the separation point and the pressure drag. With it
/// clear, `μ_w` is the plain molecular+eddy viscosity `mu_eff` — the
/// legacy crude-linear-gradient wall drag, kept for comparison.
///
/// Returned as a conductance (units of `aP`) so the caller adds it to
/// `aP` (driving the face velocity toward the zero wall velocity).
/// `axis` is the velocity component, 0/1/2 for u/v/w.
///
/// The conductance is always strictly positive for a cut cell that has
/// a wall (the wall-model effective viscosity is floored at the
/// molecular viscosity), so it keeps a near-fully-closed cut cell's
/// momentum diagonal `aP` away from zero — which the segregated solver
/// relies on (a zero `aP` face is skipped and its pressure-correction
/// coefficient would be left undefined). A non-finite result is
/// rejected.
#[allow(clippy::too_many_arguments)]
#[inline]
fn cut_wall_drag(
    body: &ImmersedBody,
    i: usize,
    j: usize,
    k: usize,
    axis: usize,
    mu_lam: f64,
    rho: f64,
    u_tan_cell: f64,
    mu_eff: f64,
    use_wall_model: bool,
) -> f64 {
    if body.method != WallMethod::CutCell || body.tag(i, j, k) != CellTag::Cut {
        return 0.0;
    }
    let cut = body.cut_face(i, j, k);
    if !cut.has_wall() {
        return 0.0;
    }
    // The tangential wall area for this velocity component: the wall
    // patch area times the sine between the wall normal and the
    // component axis — `A·(1 − n_axis²)` is the area "seen" by a
    // tangential shear on that component.
    let n_axis = cut.normal[axis];
    let tan_factor = (1.0 - n_axis * n_axis).max(0.0);
    let area_tan = cut.area * tan_factor;
    if area_tan <= 0.0 {
        return 0.0;
    }
    let drag = if use_wall_model {
        // The law-of-the-wall wall shear, linearised as a momentum
        // sink. The wall-model effective viscosity `μ_w = τ_w·y/u_t`
        // over the wall-normal sampling height `y` gives a conductance
        // `μ_w·A_tan/y = τ_w·A_tan/u_t`, so the implicit term extracts
        // exactly the law-of-the-wall wall-shear force `τ_w·A_tan` —
        // and the conductance is the *physically small* fraction of the
        // diagonal a real wall shear is (it does NOT drive the
        // first-cell velocity to zero; on a grid that does not resolve
        // the boundary layer the first cell sits at the boundary-layer
        // edge, not the wall). `μ_w` is floored at the molecular
        // viscosity inside `wall_effective_viscosity`.
        let y = cut_wall_distance(body, i, j, k);
        let nu = (mu_lam / rho).max(1e-12);
        let mu_w = crate::wallmodel::wall_effective_viscosity(rho, u_tan_cell, y, nu, mu_lam);
        mu_w * area_tan / y.max(1e-9)
    } else {
        // Legacy crude-linear-gradient wall drag — molecular + eddy
        // viscosity over a half-cell wall distance, the pre-near-wall-
        // model treatment, kept verbatim for the before/after compare.
        let grid = body.grid;
        let d = 0.5
            * (grid.dx() * grid.dx() + grid.dy() * grid.dy() + grid.dz() * grid.dz()).sqrt()
            / 3.0_f64.sqrt();
        mu_eff * area_tan / d.max(1e-9)
    };
    if drag.is_finite() && drag > 0.0 {
        drag
    } else {
        0.0
    }
}

/// The magnitude of cut cell `(i, j, k)`'s cell-centred velocity
/// projected onto its wall tangent plane — the tangential speed `u_t`
/// the near-wall model samples.
///
/// The cell-centred velocity is the average of the cell's bracketing
/// staggered faces; its component along the cut-face outward normal is
/// removed, leaving the wall-parallel speed that drives the wall shear.
/// A non-cut cell (or one with no wall) returns `0`.
#[inline]
fn cut_cell_tangential_speed(
    body: &ImmersedBody,
    u: &Field3,
    v: &Field3,
    w: &Field3,
    i: usize,
    j: usize,
    k: usize,
) -> f64 {
    if body.tag(i, j, k) != CellTag::Cut {
        return 0.0;
    }
    let cut = body.cut_face(i, j, k);
    if !cut.has_wall() {
        return 0.0;
    }
    let vel = nalgebra::Vector3::new(
        0.5 * (u.at(i, j, k) + u.at(i + 1, j, k)),
        0.5 * (v.at(i, j, k) + v.at(i, j + 1, k)),
        0.5 * (w.at(i, j, k) + w.at(i, j, k + 1)),
    );
    let tan = vel - vel.dot(&cut.normal) * cut.normal;
    tan.norm()
}

/// Solve the discrete `u`-momentum equation with a couple of
/// Gauss-Seidel sweeps. The `u` control volume is staggered in `x`.
///
/// Under the cut-cell method the convective / diffusive neighbour
/// coefficients of a face bordering a cut cell are throttled by the
/// shared open aperture of the relevant Cartesian faces, and the
/// clipped wall faces of the two bracketing cut cells add a no-slip
/// shear-drag term to `aP` — the cut face carrying the wall BC.
#[allow(clippy::too_many_arguments)]
fn solve_momentum_u(
    u: &mut Field3,
    v: &Field3,
    w: &Field3,
    p: &Field3,
    apu: &mut Field3,
    grid: &Grid3,
    rho: f64,
    mu_lam: f64,
    turb: &TurbulenceState,
    solid: &[bool],
    body: &ImmersedBody,
    bc: &BoundaryConditions,
    free: nalgebra::Vector3<f64>,
    relax: f64,
    wall_model: bool,
) {
    let (nx, ny, nz) = (grid.nx, grid.ny, grid.nz);
    let (dx, dy, dz) = (grid.dx(), grid.dy(), grid.dz());
    let is_solid = |i: usize, j: usize, k: usize| solid[i + nx * (j + ny * k)];

    for _sweep in 0..2 {
        for k in 0..nz {
            for j in 0..ny {
                for i in 1..nx {
                    // The u-face (i,j,k) sits between cells i-1 and i.
                    // If either bracketing cell is solid, this is a
                    // wall face — drive it to zero (no-slip), the
                    // immersed-boundary forcing.
                    if is_solid(i - 1, j, k) || is_solid(i, j, k) {
                        u.set(i, j, k, 0.0);
                        apu.set(i, j, k, 1.0);
                        continue;
                    }
                    // Effective viscosity averaged from the two cells.
                    let nu = 0.5
                        * (nu_eff(turb, mu_lam, rho, i - 1, j, k)
                            + nu_eff(turb, mu_lam, rho, i, j, k));
                    let dx_diff = nu * dy * dz / dx;
                    let dy_diff = nu * dx * dz / dy;
                    let dz_diff = nu * dx * dy / dz;

                    // Convective mass flows on the six faces.
                    let fe = rho * dy * dz * 0.5 * (u.at(i, j, k) + u.at(i + 1, j, k));
                    let fw = rho * dy * dz * 0.5 * (u.at(i - 1, j, k) + u.at(i, j, k));
                    let fn_ = rho * dx * dz * 0.5 * (v.at(i - 1, j + 1, k) + v.at(i, j + 1, k));
                    let fs = rho * dx * dz * 0.5 * (v.at(i - 1, j, k) + v.at(i, j, k));
                    let ft = rho * dx * dy * 0.5 * (w.at(i - 1, j, k + 1) + w.at(i, j, k + 1));
                    let fb = rho * dx * dy * 0.5 * (w.at(i - 1, j, k) + w.at(i, j, k));

                    // Cut-cell aperture throttling. The e-face of the
                    // u-CV is throttled by cell i's +x aperture, the
                    // w-face by cell (i-1)'s −x aperture; the four
                    // tangential faces by the averaged aperture of the
                    // two bracketing cells on that axis.
                    let ap_e = face_aperture(body, i, j, k, 1);
                    let ap_w = face_aperture(body, i - 1, j, k, 0);
                    let ap_n = 0.5
                        * (face_aperture(body, i - 1, j, k, 3) + face_aperture(body, i, j, k, 3));
                    let ap_s = 0.5
                        * (face_aperture(body, i - 1, j, k, 2) + face_aperture(body, i, j, k, 2));
                    let ap_t = 0.5
                        * (face_aperture(body, i - 1, j, k, 5) + face_aperture(body, i, j, k, 5));
                    let ap_b = 0.5
                        * (face_aperture(body, i - 1, j, k, 4) + face_aperture(body, i, j, k, 4));

                    let ae = ap_e * hybrid_coeff(dx_diff, fe, -1.0);
                    let aw = ap_w * hybrid_coeff(dx_diff, fw, 1.0);

                    let mut a_p = ae + aw;
                    let mut su = 0.0;
                    let (mut an, mut as_, mut at, mut ab) = (0.0, 0.0, 0.0, 0.0);

                    // y faces.
                    handle_tangential(
                        j,
                        ny,
                        dy_diff,
                        fn_,
                        fs,
                        bc.y_min,
                        bc.y_max,
                        free.x,
                        &mut a_p,
                        &mut su,
                        &mut an,
                        &mut as_,
                        nu * dx * dz,
                    );
                    // z faces.
                    handle_tangential(
                        k,
                        nz,
                        dz_diff,
                        ft,
                        fb,
                        bc.z_min,
                        bc.z_max,
                        free.x,
                        &mut a_p,
                        &mut su,
                        &mut at,
                        &mut ab,
                        nu * dx * dy,
                    );

                    // Throttle the interior tangential coefficients by
                    // their apertures (boundary faces are untouched —
                    // a domain face has no body aperture).
                    if j + 1 < ny {
                        a_p -= an;
                        an *= ap_n;
                        a_p += an;
                    }
                    if j > 0 {
                        a_p -= as_;
                        as_ *= ap_s;
                        a_p += as_;
                    }
                    if k + 1 < nz {
                        a_p -= at;
                        at *= ap_t;
                        a_p += at;
                    }
                    if k > 0 {
                        a_p -= ab;
                        ab *= ap_b;
                        a_p += ab;
                    }

                    // Cut-cell wall drag — the clipped wall faces of the
                    // two bracketing cut cells pull the face velocity
                    // toward the (zero) no-slip wall velocity, the
                    // momentum sink reconstructed from the near-wall
                    // model's law-of-the-wall profile.
                    let mu_e = mu_lam + 0.5 * (turb.mu_t.at(i - 1, j, k) + turb.mu_t.at(i, j, k));
                    let ut_wm = cut_cell_tangential_speed(body, u, v, w, i - 1, j, k);
                    a_p +=
                        cut_wall_drag(body, i - 1, j, k, 0, mu_lam, rho, ut_wm, mu_e, wall_model);
                    let ut_we = cut_cell_tangential_speed(body, u, v, w, i, j, k);
                    a_p += cut_wall_drag(body, i, j, k, 0, mu_lam, rho, ut_we, mu_e, wall_model);

                    if a_p.abs() < 1e-30 {
                        continue;
                    }
                    // Pressure gradient — exact cell-pressure
                    // difference on the staggered CV. The face area is
                    // apertured by the shared open fraction of the two
                    // bracketing cells so the predictor's pressure
                    // force is consistent with the apertured continuity
                    // equation and the apertured velocity correction.
                    let q_pg = shared_aperture(body, i - 1, j, k, 1, i, j, k, 0);
                    let dp = (p.at(i - 1, j, k) - p.at(i, j, k)) * q_pg * dy * dz;
                    let mut nb = su + dp;
                    nb += ae * u.at(i + 1, j, k);
                    nb += aw * u.at(i - 1, j, k);
                    if j + 1 < ny {
                        nb += an * u.at(i, j + 1, k);
                    }
                    if j > 0 {
                        nb += as_ * u.at(i, j - 1, k);
                    }
                    if k + 1 < nz {
                        nb += at * u.at(i, j, k + 1);
                    }
                    if k > 0 {
                        nb += ab * u.at(i, j, k - 1);
                    }
                    let u_new = nb / a_p;
                    let u_old = u.at(i, j, k);
                    u.set(i, j, k, u_old + relax * (u_new - u_old));
                    apu.set(
                        i,
                        j,
                        k,
                        small_cell_floor_ap(body, i - 1, j, k, i, j, k, a_p / relax, dx_diff),
                    );
                }
            }
        }
    }
}

/// Solve the discrete `v`-momentum equation. The `v` control volume is
/// staggered in `y`. The cut-cell aperture throttling and wall-drag are
/// applied as in [`solve_momentum_u`].
#[allow(clippy::too_many_arguments)]
fn solve_momentum_v(
    u: &Field3,
    v: &mut Field3,
    w: &Field3,
    p: &Field3,
    apv: &mut Field3,
    grid: &Grid3,
    rho: f64,
    mu_lam: f64,
    turb: &TurbulenceState,
    solid: &[bool],
    body: &ImmersedBody,
    bc: &BoundaryConditions,
    free: nalgebra::Vector3<f64>,
    relax: f64,
    wall_model: bool,
) {
    let (nx, ny, nz) = (grid.nx, grid.ny, grid.nz);
    let (dx, dy, dz) = (grid.dx(), grid.dy(), grid.dz());
    let is_solid = |i: usize, j: usize, k: usize| solid[i + nx * (j + ny * k)];

    for _sweep in 0..2 {
        for k in 0..nz {
            for j in 1..ny {
                for i in 0..nx {
                    if is_solid(i, j - 1, k) || is_solid(i, j, k) {
                        v.set(i, j, k, 0.0);
                        apv.set(i, j, k, 1.0);
                        continue;
                    }
                    let nu = 0.5
                        * (nu_eff(turb, mu_lam, rho, i, j - 1, k)
                            + nu_eff(turb, mu_lam, rho, i, j, k));
                    let dx_diff = nu * dy * dz / dx;
                    let dy_diff = nu * dx * dz / dy;
                    let dz_diff = nu * dx * dy / dz;

                    let fn_ = rho * dx * dz * 0.5 * (v.at(i, j, k) + v.at(i, j + 1, k));
                    let fs = rho * dx * dz * 0.5 * (v.at(i, j - 1, k) + v.at(i, j, k));
                    let fe = rho * dy * dz * 0.5 * (u.at(i + 1, j - 1, k) + u.at(i + 1, j, k));
                    let fw = rho * dy * dz * 0.5 * (u.at(i, j - 1, k) + u.at(i, j, k));
                    let ft = rho * dx * dy * 0.5 * (w.at(i, j - 1, k + 1) + w.at(i, j, k + 1));
                    let fb = rho * dx * dy * 0.5 * (w.at(i, j - 1, k) + w.at(i, j, k));

                    // Cut-cell apertures — n/s faces from the two
                    // bracketing cells' ±y faces, the tangential faces
                    // averaged on their axis.
                    let ap_n = face_aperture(body, i, j, k, 3);
                    let ap_s = face_aperture(body, i, j - 1, k, 2);
                    let ap_e = 0.5
                        * (face_aperture(body, i, j - 1, k, 1) + face_aperture(body, i, j, k, 1));
                    let ap_w = 0.5
                        * (face_aperture(body, i, j - 1, k, 0) + face_aperture(body, i, j, k, 0));
                    let ap_t = 0.5
                        * (face_aperture(body, i, j - 1, k, 5) + face_aperture(body, i, j, k, 5));
                    let ap_b = 0.5
                        * (face_aperture(body, i, j - 1, k, 4) + face_aperture(body, i, j, k, 4));

                    let an = ap_n * hybrid_coeff(dy_diff, fn_, -1.0);
                    let as_ = ap_s * hybrid_coeff(dy_diff, fs, 1.0);

                    let mut a_p = an + as_;
                    let mut sv = 0.0;
                    let (mut ae, mut aw, mut at, mut ab) = (0.0, 0.0, 0.0, 0.0);

                    handle_tangential(
                        i,
                        nx,
                        dx_diff,
                        fe,
                        fw,
                        bc.x_min,
                        bc.x_max,
                        free.y,
                        &mut a_p,
                        &mut sv,
                        &mut ae,
                        &mut aw,
                        nu * dy * dz,
                    );
                    handle_tangential(
                        k,
                        nz,
                        dz_diff,
                        ft,
                        fb,
                        bc.z_min,
                        bc.z_max,
                        free.y,
                        &mut a_p,
                        &mut sv,
                        &mut at,
                        &mut ab,
                        nu * dx * dy,
                    );

                    if i + 1 < nx {
                        a_p -= ae;
                        ae *= ap_e;
                        a_p += ae;
                    }
                    if i > 0 {
                        a_p -= aw;
                        aw *= ap_w;
                        a_p += aw;
                    }
                    if k + 1 < nz {
                        a_p -= at;
                        at *= ap_t;
                        a_p += at;
                    }
                    if k > 0 {
                        a_p -= ab;
                        ab *= ap_b;
                        a_p += ab;
                    }

                    let mu_e = mu_lam + 0.5 * (turb.mu_t.at(i, j - 1, k) + turb.mu_t.at(i, j, k));
                    let ut_ws = cut_cell_tangential_speed(body, u, v, w, i, j - 1, k);
                    a_p +=
                        cut_wall_drag(body, i, j - 1, k, 1, mu_lam, rho, ut_ws, mu_e, wall_model);
                    let ut_wn = cut_cell_tangential_speed(body, u, v, w, i, j, k);
                    a_p += cut_wall_drag(body, i, j, k, 1, mu_lam, rho, ut_wn, mu_e, wall_model);

                    if a_p.abs() < 1e-30 {
                        continue;
                    }
                    let q_pg = shared_aperture(body, i, j - 1, k, 3, i, j, k, 2);
                    let dp = (p.at(i, j - 1, k) - p.at(i, j, k)) * q_pg * dx * dz;
                    let mut nb = sv + dp;
                    nb += an * v.at(i, j + 1, k);
                    nb += as_ * v.at(i, j - 1, k);
                    if i + 1 < nx {
                        nb += ae * v.at(i + 1, j, k);
                    }
                    if i > 0 {
                        nb += aw * v.at(i - 1, j, k);
                    }
                    if k + 1 < nz {
                        nb += at * v.at(i, j, k + 1);
                    }
                    if k > 0 {
                        nb += ab * v.at(i, j, k - 1);
                    }
                    let v_new = nb / a_p;
                    let v_old = v.at(i, j, k);
                    v.set(i, j, k, v_old + relax * (v_new - v_old));
                    apv.set(
                        i,
                        j,
                        k,
                        small_cell_floor_ap(body, i, j - 1, k, i, j, k, a_p / relax, dy_diff),
                    );
                }
            }
        }
    }
}

/// Solve the discrete `w`-momentum equation. The `w` control volume is
/// staggered in `z`. The cut-cell aperture throttling and wall-drag are
/// applied as in [`solve_momentum_u`].
#[allow(clippy::too_many_arguments)]
fn solve_momentum_w(
    u: &Field3,
    v: &Field3,
    w: &mut Field3,
    p: &Field3,
    apw: &mut Field3,
    grid: &Grid3,
    rho: f64,
    mu_lam: f64,
    turb: &TurbulenceState,
    solid: &[bool],
    body: &ImmersedBody,
    bc: &BoundaryConditions,
    free: nalgebra::Vector3<f64>,
    relax: f64,
    wall_model: bool,
) {
    let (nx, ny, nz) = (grid.nx, grid.ny, grid.nz);
    let (dx, dy, dz) = (grid.dx(), grid.dy(), grid.dz());
    let is_solid = |i: usize, j: usize, k: usize| solid[i + nx * (j + ny * k)];

    for _sweep in 0..2 {
        for k in 1..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if is_solid(i, j, k - 1) || is_solid(i, j, k) {
                        w.set(i, j, k, 0.0);
                        apw.set(i, j, k, 1.0);
                        continue;
                    }
                    let nu = 0.5
                        * (nu_eff(turb, mu_lam, rho, i, j, k - 1)
                            + nu_eff(turb, mu_lam, rho, i, j, k));
                    let dx_diff = nu * dy * dz / dx;
                    let dy_diff = nu * dx * dz / dy;
                    let dz_diff = nu * dx * dy / dz;

                    let ft = rho * dx * dy * 0.5 * (w.at(i, j, k) + w.at(i, j, k + 1));
                    let fb = rho * dx * dy * 0.5 * (w.at(i, j, k - 1) + w.at(i, j, k));
                    let fe = rho * dy * dz * 0.5 * (u.at(i + 1, j, k - 1) + u.at(i + 1, j, k));
                    let fw = rho * dy * dz * 0.5 * (u.at(i, j, k - 1) + u.at(i, j, k));
                    let fn_ = rho * dx * dz * 0.5 * (v.at(i, j + 1, k - 1) + v.at(i, j + 1, k));
                    let fs = rho * dx * dz * 0.5 * (v.at(i, j, k - 1) + v.at(i, j, k));

                    // Cut-cell apertures — t/b faces from the two
                    // bracketing cells' ±z faces, the tangential faces
                    // averaged on their axis.
                    let ap_t = face_aperture(body, i, j, k, 5);
                    let ap_b = face_aperture(body, i, j, k - 1, 4);
                    let ap_e = 0.5
                        * (face_aperture(body, i, j, k - 1, 1) + face_aperture(body, i, j, k, 1));
                    let ap_w = 0.5
                        * (face_aperture(body, i, j, k - 1, 0) + face_aperture(body, i, j, k, 0));
                    let ap_n = 0.5
                        * (face_aperture(body, i, j, k - 1, 3) + face_aperture(body, i, j, k, 3));
                    let ap_s = 0.5
                        * (face_aperture(body, i, j, k - 1, 2) + face_aperture(body, i, j, k, 2));

                    let at = ap_t * hybrid_coeff(dz_diff, ft, -1.0);
                    let ab = ap_b * hybrid_coeff(dz_diff, fb, 1.0);

                    let mut a_p = at + ab;
                    let mut sw = 0.0;
                    let (mut ae, mut aw, mut an, mut as_) = (0.0, 0.0, 0.0, 0.0);

                    handle_tangential(
                        i,
                        nx,
                        dx_diff,
                        fe,
                        fw,
                        bc.x_min,
                        bc.x_max,
                        free.z,
                        &mut a_p,
                        &mut sw,
                        &mut ae,
                        &mut aw,
                        nu * dy * dz,
                    );
                    handle_tangential(
                        j,
                        ny,
                        dy_diff,
                        fn_,
                        fs,
                        bc.y_min,
                        bc.y_max,
                        free.z,
                        &mut a_p,
                        &mut sw,
                        &mut an,
                        &mut as_,
                        nu * dx * dz,
                    );

                    if i + 1 < nx {
                        a_p -= ae;
                        ae *= ap_e;
                        a_p += ae;
                    }
                    if i > 0 {
                        a_p -= aw;
                        aw *= ap_w;
                        a_p += aw;
                    }
                    if j + 1 < ny {
                        a_p -= an;
                        an *= ap_n;
                        a_p += an;
                    }
                    if j > 0 {
                        a_p -= as_;
                        as_ *= ap_s;
                        a_p += as_;
                    }

                    let mu_e = mu_lam + 0.5 * (turb.mu_t.at(i, j, k - 1) + turb.mu_t.at(i, j, k));
                    let ut_wb = cut_cell_tangential_speed(body, u, v, w, i, j, k - 1);
                    a_p +=
                        cut_wall_drag(body, i, j, k - 1, 2, mu_lam, rho, ut_wb, mu_e, wall_model);
                    let ut_wt = cut_cell_tangential_speed(body, u, v, w, i, j, k);
                    a_p += cut_wall_drag(body, i, j, k, 2, mu_lam, rho, ut_wt, mu_e, wall_model);

                    if a_p.abs() < 1e-30 {
                        continue;
                    }
                    let q_pg = shared_aperture(body, i, j, k - 1, 5, i, j, k, 4);
                    let dp = (p.at(i, j, k - 1) - p.at(i, j, k)) * q_pg * dx * dy;
                    let mut nb = sw + dp;
                    nb += at * w.at(i, j, k + 1);
                    nb += ab * w.at(i, j, k - 1);
                    if i + 1 < nx {
                        nb += ae * w.at(i + 1, j, k);
                    }
                    if i > 0 {
                        nb += aw * w.at(i - 1, j, k);
                    }
                    if j + 1 < ny {
                        nb += an * w.at(i, j + 1, k);
                    }
                    if j > 0 {
                        nb += as_ * w.at(i, j - 1, k);
                    }
                    let w_new = nb / a_p;
                    let w_old = w.at(i, j, k);
                    w.set(i, j, k, w_old + relax * (w_new - w_old));
                    apw.set(
                        i,
                        j,
                        k,
                        small_cell_floor_ap(body, i, j, k - 1, i, j, k, a_p / relax, dz_diff),
                    );
                }
            }
        }
    }
}

/// Handle one pair of opposed tangential faces of a momentum control
/// volume — pick the regular hybrid neighbour coefficient when the
/// face is interior, or apply the domain-face boundary treatment
/// (slip = zero shear, wall / inlet = wall-shear term) at a boundary.
#[allow(clippy::too_many_arguments)]
fn handle_tangential(
    idx: usize,
    n: usize,
    diff: f64,
    f_plus: f64,
    f_minus: f64,
    bc_min: FaceBc,
    bc_max: FaceBc,
    free_tangential: f64,
    a_p: &mut f64,
    src: &mut f64,
    a_plus: &mut f64,
    a_minus: &mut f64,
    wall_face_area: f64,
) {
    // The "plus" side (idx == n-1 is the domain max face).
    if idx == n - 1 {
        match bc_max {
            FaceBc::Slip | FaceBc::Outlet => { /* zero shear — no term */ }
            FaceBc::Inlet => {
                // The inlet carries the free-stream tangential value.
                let dwall = diff * 2.0; // half-cell wall distance
                let _ = wall_face_area;
                *a_p += dwall;
                *src += dwall * free_tangential;
            }
            FaceBc::MovingWall(speed) => {
                let dwall = diff * 2.0;
                *a_p += dwall;
                *src += dwall * speed;
            }
        }
    } else {
        *a_plus = hybrid_coeff(diff, f_plus, -1.0);
        *a_p += *a_plus;
    }
    // The "minus" side (idx == 0 is the domain min face).
    if idx == 0 {
        match bc_min {
            FaceBc::Slip | FaceBc::Outlet => {}
            FaceBc::Inlet => {
                let dwall = diff * 2.0;
                *a_p += dwall;
                *src += dwall * free_tangential;
            }
            FaceBc::MovingWall(speed) => {
                let dwall = diff * 2.0;
                *a_p += dwall;
                *src += dwall * speed;
            }
        }
    } else {
        *a_minus = hybrid_coeff(diff, f_minus, 1.0);
        *a_p += *a_minus;
    }
}

/// Stamp the prescribed boundary velocities onto the staggered fields.
///
/// On the staggered grid the *normal* velocity of an inlet / wall lies
/// exactly on the boundary face, so it is set directly here; the
/// tangential component has no boundary face and is handled inside the
/// momentum sweep. A slip face sets the normal velocity to zero (no
/// penetration); an outlet copies the last interior face (zero
/// gradient / convective).
fn apply_boundary_velocities(
    u: &mut Field3,
    v: &mut Field3,
    w: &mut Field3,
    grid: &Grid3,
    bc: &BoundaryConditions,
    free: nalgebra::Vector3<f64>,
) {
    let (nx, ny, nz) = (grid.nx, grid.ny, grid.nz);
    // x faces — u at i = 0 and i = nx.
    for k in 0..nz {
        for j in 0..ny {
            match bc.x_min {
                FaceBc::Inlet => u.set(0, j, k, free.x),
                FaceBc::Slip => u.set(0, j, k, 0.0),
                FaceBc::MovingWall(_) => u.set(0, j, k, 0.0),
                FaceBc::Outlet => {
                    let interior = u.at(1, j, k);
                    u.set(0, j, k, interior);
                }
            }
            match bc.x_max {
                FaceBc::Inlet => u.set(nx, j, k, free.x),
                FaceBc::Slip => u.set(nx, j, k, 0.0),
                FaceBc::MovingWall(_) => u.set(nx, j, k, 0.0),
                FaceBc::Outlet => {
                    let interior = u.at(nx - 1, j, k);
                    u.set(nx, j, k, interior);
                }
            }
        }
    }
    // y faces — v at j = 0 and j = ny.
    for k in 0..nz {
        for i in 0..nx {
            match bc.y_min {
                FaceBc::Inlet => v.set(i, 0, k, free.y),
                FaceBc::Slip | FaceBc::MovingWall(_) => v.set(i, 0, k, 0.0),
                FaceBc::Outlet => {
                    let interior = v.at(i, 1, k);
                    v.set(i, 0, k, interior);
                }
            }
            match bc.y_max {
                FaceBc::Inlet => v.set(i, ny, k, free.y),
                FaceBc::Slip | FaceBc::MovingWall(_) => v.set(i, ny, k, 0.0),
                FaceBc::Outlet => {
                    let interior = v.at(i, ny - 1, k);
                    v.set(i, ny, k, interior);
                }
            }
        }
    }
    // z faces — w at k = 0 and k = nz.
    for j in 0..ny {
        for i in 0..nx {
            match bc.z_min {
                FaceBc::Inlet => w.set(i, j, 0, free.z),
                FaceBc::Slip | FaceBc::MovingWall(_) => w.set(i, j, 0, 0.0),
                FaceBc::Outlet => {
                    let interior = w.at(i, j, 1);
                    w.set(i, j, 0, interior);
                }
            }
            match bc.z_max {
                FaceBc::Inlet => w.set(i, j, nz, free.z),
                FaceBc::Slip | FaceBc::MovingWall(_) => w.set(i, j, nz, 0.0),
                FaceBc::Outlet => {
                    let interior = w.at(i, j, nz - 1);
                    w.set(i, j, nz, interior);
                }
            }
        }
    }
}

/// Freeze the velocity faces of solid cells at the body velocity (the
/// immersed-boundary no-slip / moving-wall forcing).
fn apply_immersed_velocities(
    u: &mut Field3,
    v: &mut Field3,
    w: &mut Field3,
    grid: &Grid3,
    solid: &[bool],
    motion: &BodyMotion,
) {
    let (nx, ny, nz) = (grid.nx, grid.ny, grid.nz);
    let is_solid = |i: usize, j: usize, k: usize| solid[i + nx * (j + ny * k)];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                if !is_solid(i, j, k) {
                    continue;
                }
                let vel = motion.at(i + nx * (j + ny * k));
                // Every face of a solid cell carries the body velocity.
                u.set(i, j, k, vel[0]);
                u.set(i + 1, j, k, vel[0]);
                v.set(i, j, k, vel[1]);
                v.set(i, j + 1, k, vel[1]);
                w.set(i, j, k, vel[2]);
                w.set(i, j, k + 1, vel[2]);
            }
        }
    }
}

/// Assemble the pressure-correction Poisson stencil and return the L2
/// norm of the per-cell mass imbalance (the convergence residual).
///
/// This is the discrete continuity equation. Under the cut-cell method
/// every Cartesian face flux is scaled by the face's open aperture — a
/// fully-closed face passes no mass and contributes no `p'` coupling —
/// and the clipped wall face, being a no-slip wall (zero normal
/// velocity), contributes nothing, exactly as a closed face does. Under
/// the staircased method every aperture is `1` and this reduces to the
/// plain seven-point Poisson stencil.
#[allow(clippy::too_many_arguments)]
fn assemble_pressure_correction(
    stencil: &mut PoissonStencil,
    u: &Field3,
    v: &Field3,
    w: &Field3,
    apu: &Field3,
    apv: &Field3,
    apw: &Field3,
    grid: &Grid3,
    rho: f64,
    solid: &[bool],
    body: &ImmersedBody,
    bc: &BoundaryConditions,
) -> f64 {
    let (nx, ny, nz) = (grid.nx, grid.ny, grid.nz);
    let (dx, dy, dz) = (grid.dx(), grid.dy(), grid.dz());
    let ax = dy * dz;
    let ay = dx * dz;
    let az = dx * dy;
    let is_solid = |i: usize, j: usize, k: usize| solid[i + nx * (j + ny * k)];

    let mut sum_sq = 0.0;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                if is_solid(i, j, k) {
                    // A solid cell has no pressure correction.
                    stencil.ap.set(i, j, k, 1.0);
                    stencil.b.set(i, j, k, 0.0);
                    continue;
                }
                // The six Cartesian-face open apertures of this cell —
                // the shared aperture with each neighbour. `1.0` for a
                // staircased run or a fully-open face.
                let qe = if i + 1 < nx {
                    shared_aperture(body, i, j, k, 1, i + 1, j, k, 0)
                } else {
                    1.0
                };
                let qw = if i > 0 {
                    shared_aperture(body, i, j, k, 0, i - 1, j, k, 1)
                } else {
                    1.0
                };
                let qn = if j + 1 < ny {
                    shared_aperture(body, i, j, k, 3, i, j + 1, k, 2)
                } else {
                    1.0
                };
                let qs = if j > 0 {
                    shared_aperture(body, i, j, k, 2, i, j - 1, k, 3)
                } else {
                    1.0
                };
                let qt = if k + 1 < nz {
                    shared_aperture(body, i, j, k, 5, i, j, k + 1, 4)
                } else {
                    1.0
                };
                let qb = if k > 0 {
                    shared_aperture(body, i, j, k, 4, i, j, k - 1, 5)
                } else {
                    1.0
                };

                // Velocity-correction coefficient on a face:
                // ρ·(aperture·A)²/aP. A face touching a solid cell or a
                // domain boundary gets zero (homogeneous Neumann on p').
                let ae = if i + 1 < nx && !is_solid(i + 1, j, k) {
                    rho * (qe * ax) * (qe * ax) / apu.at(i + 1, j, k).max(1e-30)
                } else {
                    0.0
                };
                let aw = if i > 0 && !is_solid(i - 1, j, k) {
                    rho * (qw * ax) * (qw * ax) / apu.at(i, j, k).max(1e-30)
                } else {
                    0.0
                };
                let an = if j + 1 < ny && !is_solid(i, j + 1, k) {
                    rho * (qn * ay) * (qn * ay) / apv.at(i, j + 1, k).max(1e-30)
                } else {
                    0.0
                };
                let as_ = if j > 0 && !is_solid(i, j - 1, k) {
                    rho * (qs * ay) * (qs * ay) / apv.at(i, j, k).max(1e-30)
                } else {
                    0.0
                };
                let at = if k + 1 < nz && !is_solid(i, j, k + 1) {
                    rho * (qt * az) * (qt * az) / apw.at(i, j, k + 1).max(1e-30)
                } else {
                    0.0
                };
                let ab = if k > 0 && !is_solid(i, j, k - 1) {
                    rho * (qb * az) * (qb * az) / apw.at(i, j, k).max(1e-30)
                } else {
                    0.0
                };

                // The cell's net mass outflow with the provisional
                // velocity, each face flux scaled by its open aperture.
                // The pressure-correction source is its negative —
                // continuity wants the imbalance → 0.
                let mass_out = rho * ax * (qe * u.at(i + 1, j, k) - qw * u.at(i, j, k))
                    + rho * ay * (qn * v.at(i, j + 1, k) - qs * v.at(i, j, k))
                    + rho * az * (qt * w.at(i, j, k + 1) - qb * w.at(i, j, k));

                stencil.ae.set(i, j, k, ae);
                stencil.aw.set(i, j, k, aw);
                stencil.an.set(i, j, k, an);
                stencil.as_.set(i, j, k, as_);
                stencil.at.set(i, j, k, at);
                stencil.ab.set(i, j, k, ab);
                stencil.ap.set(i, j, k, ae + aw + an + as_ + at + ab);
                stencil.b.set(i, j, k, -mass_out);

                sum_sq += mass_out * mass_out;
            }
        }
    }

    // Pin one outlet cell's correction to zero — the natural pressure
    // reference for an open external-aero domain.
    //
    // The pin is a clean identity row: `ap = 1`, every off-diagonal and
    // the source `0`, so the cell's value is driven to exactly `0`.
    // It must NOT be a `1e30`-diagonal "big number" row — a `1e30`
    // entry makes the cell's residual `1e30·p'` blow up the global
    // residual norm and the multigrid restriction the instant `p'` is
    // even slightly non-zero, which sends the V-cycle to infinity.
    let anchor = outlet_anchor(grid, bc);
    if let Some((ai, aj, ak)) = anchor {
        stencil.ap.set(ai, aj, ak, 1.0);
        stencil.ae.set(ai, aj, ak, 0.0);
        stencil.aw.set(ai, aj, ak, 0.0);
        stencil.an.set(ai, aj, ak, 0.0);
        stencil.as_.set(ai, aj, ak, 0.0);
        stencil.at.set(ai, aj, ak, 0.0);
        stencil.ab.set(ai, aj, ak, 0.0);
        stencil.b.set(ai, aj, ak, 0.0);
    }

    let n = (nx * ny * nz) as f64;
    if n > 0.0 {
        (sum_sq / n).sqrt()
    } else {
        0.0
    }
}

/// Pick a cell on an outlet face to anchor the pressure gauge.
fn outlet_anchor(grid: &Grid3, bc: &BoundaryConditions) -> Option<(usize, usize, usize)> {
    let (nx, ny, nz) = (grid.nx, grid.ny, grid.nz);
    if bc.x_max == FaceBc::Outlet {
        Some((nx - 1, ny / 2, nz / 2))
    } else if bc.x_min == FaceBc::Outlet {
        Some((0, ny / 2, nz / 2))
    } else if bc.y_max == FaceBc::Outlet {
        Some((nx / 2, ny - 1, nz / 2))
    } else if bc.y_min == FaceBc::Outlet {
        Some((nx / 2, 0, nz / 2))
    } else if bc.z_max == FaceBc::Outlet {
        Some((nx / 2, ny / 2, nz - 1))
    } else if bc.z_min == FaceBc::Outlet {
        Some((nx / 2, ny / 2, 0))
    } else {
        None
    }
}

/// Apply the under-relaxed pressure correction `p ← p + αp·p'` to the
/// fluid cells.
fn correct_pressure(p: &mut Field3, pcorr: &Field3, solid: &[bool], relax_p: f64) {
    for (idx, is_solid) in solid.iter().enumerate().take(p.data.len()) {
        if !is_solid {
            p.data[idx] += relax_p * pcorr.data[idx];
        }
    }
}

/// Apply the velocity correction so the corrected field is
/// divergence-free. Each interior fluid face is corrected by
/// `u' = d·(p'_P − p'_E)` with `d = aperture·A/aP`; a wall face is left
/// at the no-slip value. The aperture factor keeps the correction
/// consistent with the apertured continuity equation that produced the
/// pressure correction.
#[allow(clippy::too_many_arguments)]
fn correct_velocity(
    u: &mut Field3,
    v: &mut Field3,
    w: &mut Field3,
    pcorr: &Field3,
    apu: &Field3,
    apv: &Field3,
    apw: &Field3,
    grid: &Grid3,
    solid: &[bool],
    body: &ImmersedBody,
) {
    let (nx, ny, nz) = (grid.nx, grid.ny, grid.nz);
    let (dx, dy, dz) = (grid.dx(), grid.dy(), grid.dz());
    let ax = dy * dz;
    let ay = dx * dz;
    let az = dx * dy;
    let is_solid = |i: usize, j: usize, k: usize| solid[i + nx * (j + ny * k)];

    // Interior u-faces.
    for k in 0..nz {
        for j in 0..ny {
            for i in 1..nx {
                if is_solid(i - 1, j, k) || is_solid(i, j, k) {
                    continue;
                }
                let q = shared_aperture(body, i - 1, j, k, 1, i, j, k, 0);
                let d = q * ax / apu.at(i, j, k).max(1e-30);
                let du = d * (pcorr.at(i - 1, j, k) - pcorr.at(i, j, k));
                u.add(i, j, k, du);
            }
        }
    }
    // Interior v-faces.
    for k in 0..nz {
        for j in 1..ny {
            for i in 0..nx {
                if is_solid(i, j - 1, k) || is_solid(i, j, k) {
                    continue;
                }
                let q = shared_aperture(body, i, j - 1, k, 3, i, j, k, 2);
                let d = q * ay / apv.at(i, j, k).max(1e-30);
                let dv = d * (pcorr.at(i, j - 1, k) - pcorr.at(i, j, k));
                v.add(i, j, k, dv);
            }
        }
    }
    // Interior w-faces.
    for k in 1..nz {
        for j in 0..ny {
            for i in 0..nx {
                if is_solid(i, j, k - 1) || is_solid(i, j, k) {
                    continue;
                }
                let q = shared_aperture(body, i, j, k - 1, 5, i, j, k, 4);
                let d = q * az / apw.at(i, j, k).max(1e-30);
                let dw = d * (pcorr.at(i, j, k - 1) - pcorr.at(i, j, k));
                w.add(i, j, k, dw);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::WindTunnel;
    use crate::geometry::{box_body, sphere_body};
    use crate::wind::Wind;
    use nalgebra::Vector3;

    #[test]
    fn empty_tunnel_preserves_the_free_stream() {
        // The defining sanity test: with no body, the solver must
        // leave the uniform free-stream untouched — every cell stays
        // at U∞ along x.
        let body = box_body(Vector3::zeros(), Vector3::new(0.2, 0.2, 0.2));
        let wind = Wind::straight(10.0).unwrap();
        let mut tunnel = WindTunnel::build(&body, wind).unwrap();
        // Strip the body out — make every cell fluid.
        for t in tunnel.body.tags.iter_mut() {
            *t = CellTag::Fluid;
        }
        for d in tunnel.body.wall_distance.iter_mut() {
            *d = 1e9;
        }
        let controls = SolverControls {
            max_iterations: 20,
            turbulence: TurbulenceModel::Laminar,
            ..SolverControls::default()
        };
        let flow = solve_steady(&tunnel, &controls, &BodyMotion::static_body());
        // Sample a cell in the middle of the domain.
        let (i, j, k) = (tunnel.grid.nx / 2, tunnel.grid.ny / 2, tunnel.grid.nz / 2);
        let u = flow.u_at_cell(i, j, k);
        assert!(
            (u - 10.0).abs() < 0.5,
            "empty tunnel should keep the free-stream, u = {u}"
        );
        // The transverse components stay near zero.
        assert!(flow.v_at_cell(i, j, k).abs() < 0.5);
        assert!(flow.w_at_cell(i, j, k).abs() < 0.5);
    }

    #[test]
    fn flow_blocks_inside_a_solid_body() {
        // Inside the body the velocity must be the (zero) body velocity
        // — the immersed boundary forcing.
        let body = box_body(Vector3::new(-1.0, -1.0, -1.0), Vector3::new(1.0, 1.0, 1.0));
        let wind = Wind::straight(15.0).unwrap();
        let tunnel = WindTunnel::build(&body, wind).unwrap();
        let controls = SolverControls {
            max_iterations: 30,
            turbulence: TurbulenceModel::Laminar,
            ..SolverControls::default()
        };
        let flow = solve_steady(&tunnel, &controls, &BodyMotion::static_body());
        // Find a solid cell and confirm its velocity is ~0.
        let mut checked = false;
        for k in 0..tunnel.grid.nz {
            for j in 0..tunnel.grid.ny {
                for i in 0..tunnel.grid.nx {
                    if tunnel.body.is_solid(i, j, k) {
                        assert!(
                            flow.speed_at_cell(i, j, k) < 0.5,
                            "solid cell should have ~zero velocity"
                        );
                        checked = true;
                    }
                }
            }
        }
        assert!(checked, "the box should have produced solid cells");
    }

    #[test]
    fn flow_accelerates_around_a_bluff_body() {
        // Mass conservation: the flow squeezed past the sides of a
        // bluff body must speed up somewhere relative to the
        // free-stream, and slow down (stagnate) right in front. A
        // coarse grid keeps the steady SIMPLE solve fast — this is a
        // qualitative mass-conservation check.
        let body = sphere_body(Vector3::zeros(), 0.5, 20, 40);
        let wind = Wind::straight(20.0).unwrap();
        let tunnel = WindTunnel::build_with(
            &body,
            wind,
            crate::domain::BoundaryConditions::external_aero(),
            crate::domain::TunnelSizing {
                cells_across_body: 4,
                max_cells: 40_000,
                ..crate::domain::TunnelSizing::default()
            },
        )
        .unwrap();
        let controls = SolverControls {
            max_iterations: 60,
            turbulence: TurbulenceModel::KEpsilon,
            ..SolverControls::default()
        };
        let flow = solve_steady(&tunnel, &controls, &BodyMotion::static_body());
        // The peak speed in the domain should exceed the free-stream
        // (flow speeds up squeezing past the body).
        let mut peak = 0.0f64;
        for k in 0..tunnel.grid.nz {
            for j in 0..tunnel.grid.ny {
                for i in 0..tunnel.grid.nx {
                    if tunnel.body.is_fluid_like(i, j, k) {
                        peak = peak.max(flow.speed_at_cell(i, j, k));
                    }
                }
            }
        }
        assert!(
            peak > 20.0,
            "flow should accelerate past the body, peak {peak} vs U∞ 20"
        );
        // The field must remain finite (no blow-up).
        assert!(flow.residual.is_finite());
    }

    #[test]
    fn residual_history_is_recorded() {
        let body = box_body(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let tunnel = WindTunnel::build(&body, Wind::straight(10.0).unwrap()).unwrap();
        let controls = SolverControls {
            max_iterations: 15,
            turbulence: TurbulenceModel::Laminar,
            ..SolverControls::default()
        };
        let flow = solve_steady(&tunnel, &controls, &BodyMotion::static_body());
        assert_eq!(flow.residual_history.len(), flow.iterations);
        assert!(flow.residual_history.iter().all(|r| r.is_finite()));
    }

    #[test]
    fn solver_does_not_blow_up_on_a_car_like_body() {
        // A box car at a road-traffic speed with k-ω SST must finish
        // with a finite residual — the robustness clamp + the
        // turbulence cap keep it stable. A coarse grid keeps the run
        // fast; the test only checks the field stays finite.
        let body = box_body(Vector3::zeros(), Vector3::new(4.0, 1.8, 1.4));
        let tunnel = WindTunnel::build_with(
            &body,
            Wind::straight(30.0).unwrap(),
            BoundaryConditions::automotive(30.0),
            crate::domain::TunnelSizing {
                cells_across_body: 4,
                max_cells: 40_000,
                ..crate::domain::TunnelSizing::default()
            },
        )
        .unwrap();
        let controls = SolverControls {
            max_iterations: 60,
            turbulence: TurbulenceModel::KOmegaSST,
            ..SolverControls::default()
        };
        let flow = solve_steady(&tunnel, &controls, &BodyMotion::static_body());
        assert!(flow.residual.is_finite(), "car run blew up");
        assert!(flow.turbulence.mean_eddy_viscosity().is_finite());
    }
}
