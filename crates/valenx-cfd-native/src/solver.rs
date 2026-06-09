//! The SIMPLE pressure-velocity-coupling solver for steady,
//! incompressible, laminar 2-D flow.
//!
//! # The equations
//!
//! Steady incompressible Navier-Stokes on the 2-D domain:
//!
//! ```text
//!   ∇·u = 0                                    (continuity)
//!   (u·∇)u = −(1/ρ)∇p + ν∇²u                   (momentum)
//! ```
//!
//! discretised by the **finite-volume method** on the staggered grid
//! of [`crate::grid`]: each velocity component's momentum equation is
//! integrated over its own staggered control volume, the convective
//! flux uses the **hybrid (Spalding) scheme** — central differencing
//! where the cell Péclet number is low, upwind where it is high, which
//! is unconditionally stable — and the diffusive flux is central.
//!
//! # SIMPLE
//!
//! The pressure does not have its own transport equation; SIMPLE
//! (Patankar & Spalding 1972, *Semi-Implicit Method for
//! Pressure-Linked Equations*) recovers it by iteration:
//!
//! 1. with the current pressure `p*`, solve the momentum equations for
//!    a provisional velocity `u*` — it satisfies momentum but **not**
//!    continuity;
//! 2. assemble the **pressure-correction** Poisson equation whose
//!    source is each cell's mass imbalance and solve it for `p'`
//!    ([`crate::linsolve`]);
//! 3. correct the pressure (`p = p* + αp·p'`) and the velocities
//!    (`u = u* + d·∂p'/∂x`) so the corrected velocity *is*
//!    divergence-free;
//! 4. **under-relax** the momentum update (factor `αu`) for stability
//!    and repeat until the mass imbalance vanishes.
//!
//! At convergence the velocity satisfies both momentum and continuity —
//! it is the steady laminar solution.
//!
//! # Honest scope
//!
//! This is a **real working laminar solver** — the lid-driven-cavity
//! and developing-channel-flow tests in this module verify it produces
//! the physically correct circulation and the parabolic profile. It is
//! deliberately a v1, **not OpenFOAM**: 2-D only, structured uniform
//! rectangular grid only, **laminar** (no turbulence model — valid
//! while the Reynolds number stays modest), steady-state (no transient
//! time-stepping), single-phase, constant properties, and the
//! hybrid convection scheme rather than a higher-order TVD scheme.
//! Each of those is a documented, well-understood extension.

use crate::grid::{Field, Grid};
use crate::linsolve::PoissonCoeffs;
use crate::multigrid::{solve_pressure_poisson, PressurePoissonSolver};
use crate::turbulence::{
    advance_k_epsilon, advance_k_omega_sst, wall_distance_field,
    KEpsilonModel, SstField, SstModel, TurbulenceField, WallFunction, WallMask,
};

/// Fluid properties — constant density and kinematic viscosity.
#[derive(Clone, Copy, Debug)]
pub struct Fluid {
    /// Density `ρ` (kg·m⁻³).
    pub density: f64,
    /// Kinematic viscosity `ν` (m²·s⁻¹).
    pub viscosity: f64,
}

impl Fluid {
    /// Build a fluid; both properties must be positive.
    pub fn new(density: f64, viscosity: f64) -> Fluid {
        Fluid {
            density,
            viscosity,
        }
    }
}

impl Default for Fluid {
    /// A unit fluid (`ρ = 1`, `ν = 1`) — convenient for non-dimensional
    /// runs where the Reynolds number is set through the geometry /
    /// boundary velocity instead.
    fn default() -> Self {
        Fluid {
            density: 1.0,
            viscosity: 1.0,
        }
    }
}

/// The boundary condition on one side of the rectangular domain.
///
/// The two velocity components are prescribed independently so a side
/// can be a moving wall (lid-driven cavity), a uniform inlet, a
/// stationary no-slip wall, or a zero-gradient outlet.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SideBc {
    /// A solid wall moving with the given tangential / normal velocity
    /// `(u, v)`. A stationary no-slip wall is `Wall { u: 0, v: 0 }`; a
    /// lid is a wall with a non-zero tangential component.
    Wall {
        /// Prescribed `u` (x) velocity on the wall.
        u: f64,
        /// Prescribed `v` (y) velocity on the wall.
        v: f64,
    },
    /// A velocity inlet — the flow enters with the prescribed
    /// `(u, v)`.
    Inlet {
        /// Inflow `u` (x) velocity.
        u: f64,
        /// Inflow `v` (y) velocity.
        v: f64,
    },
    /// A zero-gradient outflow — the velocity normal gradient is set to
    /// zero (fully-developed-flow assumption) and the cell pressure is
    /// pinned so the singular system has a reference.
    Outlet,
}

/// The four-sided boundary specification of the rectangular domain.
#[derive(Clone, Copy, Debug)]
pub struct Boundaries {
    /// Left edge (`x = 0`).
    pub west: SideBc,
    /// Right edge (`x = lx`).
    pub east: SideBc,
    /// Bottom edge (`y = 0`).
    pub south: SideBc,
    /// Top edge (`y = ly`).
    pub north: SideBc,
}

impl Boundaries {
    /// The lid-driven-cavity boundary set: stationary no-slip walls on
    /// three sides, the top wall sliding at speed `lid_speed` in `+x`.
    pub fn lid_driven_cavity(lid_speed: f64) -> Boundaries {
        let wall = SideBc::Wall { u: 0.0, v: 0.0 };
        Boundaries {
            west: wall,
            east: wall,
            south: wall,
            north: SideBc::Wall {
                u: lid_speed,
                v: 0.0,
            },
        }
    }

    /// The plane-channel-flow boundary set: a uniform `inlet_speed`
    /// inlet on the west edge, a zero-gradient outlet on the east edge,
    /// and stationary no-slip walls top and bottom.
    pub fn channel_flow(inlet_speed: f64) -> Boundaries {
        let wall = SideBc::Wall { u: 0.0, v: 0.0 };
        Boundaries {
            west: SideBc::Inlet {
                u: inlet_speed,
                v: 0.0,
            },
            east: SideBc::Outlet,
            south: wall,
            north: wall,
        }
    }
}

/// Solver controls — under-relaxation factors and convergence limits.
#[derive(Clone, Copy, Debug)]
pub struct SimpleControls {
    /// Momentum under-relaxation factor `αu` (`0 < αu ≤ 1`). SIMPLE
    /// needs `αu ≈ 0.5–0.7` for stability — the momentum update is
    /// damped to keep the nonlinear iteration from diverging.
    pub relax_u: f64,
    /// Pressure under-relaxation factor `αp` (`0 < αp ≤ 1`). Typically
    /// `αp ≈ 0.2–0.3` — SIMPLE over-predicts the pressure correction,
    /// so only a fraction is applied.
    pub relax_p: f64,
    /// Maximum number of SIMPLE outer iterations.
    pub max_iterations: usize,
    /// Convergence tolerance on the scaled mass-imbalance residual —
    /// the run stops when it drops below this.
    pub tolerance: f64,
    /// Over-relaxation factor for the inner SOR pressure-correction
    /// solve (used when `pressure_solver = SOR`).
    pub sor_omega: f64,
    /// Maximum inner SOR sweeps per outer iteration (used when
    /// `pressure_solver = SOR`).
    pub sor_iterations: usize,
    /// **Which pressure-correction Poisson solver to use.** The
    /// historical default is [`PressurePoissonSolver::Sor`]; the
    /// production choice on fine 2-D grids is
    /// [`PressurePoissonSolver::Multigrid`] (whose V-cycle is
    /// grid-independent — SOR alone scales poorly above ~64²).
    pub pressure_solver: PressurePoissonSolver,
}

impl Default for SimpleControls {
    /// Conservative, broadly-stable defaults: `αu = 0.6`, `αp = 0.25`,
    /// up to 4000 outer iterations to a `1e-5` residual, SOR for the
    /// pressure correction (the historical default).
    fn default() -> Self {
        SimpleControls {
            relax_u: 0.6,
            relax_p: 0.25,
            max_iterations: 4000,
            tolerance: 1e-5,
            sor_omega: 1.7,
            sor_iterations: 40,
            pressure_solver: PressurePoissonSolver::Sor,
        }
    }
}

/// Which **effective viscosity** the SIMPLE momentum solve uses —
/// the laminar molecular `μ`, or `μ_eff = μ + μ_t` produced by one of
/// the two-equation RANS closures.
///
/// The two-equation models live in [`crate::turbulence`]; this enum
/// is the agent-facing selector the SIMPLE driver consults each outer
/// iteration. Choosing a turbulent variant tells [`solve_simple_with`]
/// to (1) carry a cell-centred eddy-viscosity field alongside the flow,
/// (2) advance the model's transport equations between SIMPLE outer
/// iterations, and (3) add `μ_t` to the molecular `μ` in every
/// momentum face diffusion coefficient.
#[derive(Clone, Copy, Debug)]
pub enum EffectiveViscosity {
    /// **Laminar** — the momentum solve uses the molecular viscosity
    /// alone. Equivalent to the historical [`solve_simple`] path.
    Laminar,
    /// **Standard k-ε** (Launder & Spalding 1974). Two transport
    /// equations for `k` and `ε`, eddy viscosity
    /// `ν_t = C_μ·k²/ε`, log-law wall functions on the four sides
    /// flagged in `walls`.
    KEpsilon {
        /// The five model constants.
        model: KEpsilonModel,
        /// The wall-function constants.
        wall: WallFunction,
        /// Which sides of the domain are no-slip walls — wall
        /// functions apply only there.
        walls: WallMask,
    },
    /// **Menter k-ω SST** (Menter 1994). Two transport equations for
    /// `k` and `ω`, the SST-limited eddy viscosity
    /// `ν_t = a₁·k / max(a₁·ω, S·F₂)`, the `F₁` / `F₂` blending
    /// functions built from a wall-distance field, log-law wall
    /// functions for the wall cells, the production limiter and the
    /// `(1−F₁)` cross-diffusion term.
    Sst {
        /// The Menter SST constants and the two coefficient sets.
        model: SstModel,
        /// The wall-function constants.
        wall: WallFunction,
        /// Which sides of the domain are no-slip walls.
        walls: WallMask,
    },
}

impl Default for EffectiveViscosity {
    /// The historical default — a laminar solve.
    fn default() -> Self {
        EffectiveViscosity::Laminar
    }
}

/// The converged (or iteration-capped) flow field.
#[derive(Clone, Debug)]
pub struct FlowSolution {
    /// The grid the solution lives on.
    pub grid: Grid,
    /// Staggered `u`-velocity, `(nx+1) × ny`.
    pub u: Field,
    /// Staggered `v`-velocity, `nx × (ny+1)`.
    pub v: Field,
    /// Cell-centred pressure, `nx × ny` (gauge: zero-mean unless an
    /// outlet pinned it).
    pub pressure: Field,
    /// Number of SIMPLE outer iterations performed.
    pub iterations: usize,
    /// Final scaled mass-imbalance residual.
    pub residual: f64,
    /// True if the residual reached [`SimpleControls::tolerance`].
    pub converged: bool,
}

impl FlowSolution {
    /// Cell-centred `u`-velocity of pressure cell `(i, j)` — the
    /// average of the two `u`-faces bracketing the cell. Useful for
    /// post-processing / plotting on the pressure grid.
    pub fn u_at_cell(&self, i: usize, j: usize) -> f64 {
        0.5 * (self.u.at(i, j) + self.u.at(i + 1, j))
    }

    /// Cell-centred `v`-velocity of pressure cell `(i, j)` — the
    /// average of the two `v`-faces bracketing the cell.
    pub fn v_at_cell(&self, i: usize, j: usize) -> f64 {
        0.5 * (self.v.at(i, j) + self.v.at(i, j + 1))
    }

    /// Velocity magnitude at the centre of pressure cell `(i, j)`.
    pub fn speed_at_cell(&self, i: usize, j: usize) -> f64 {
        let u = self.u_at_cell(i, j);
        let v = self.v_at_cell(i, j);
        (u * u + v * v).sqrt()
    }

    /// The **kinetic-energy density** `½·ρ·|u|²` (J/m³ = Pa) at the centre of cell
    /// `(i, j)`, for fluid density `density` `ρ` (kg/m³) — the per-cell scalar field
    /// that [`FlowSolution::mean_kinetic_energy_density`] area-averages over the grid,
    /// and the local form of the freestream dynamic pressure `½ρU²`. Built from the
    /// cell-centred [`FlowSolution::speed_at_cell`], so it is defined on every cell.
    pub fn kinetic_energy_density_at_cell(&self, i: usize, j: usize, density: f64) -> f64 {
        let speed = self.speed_at_cell(i, j);
        0.5 * density * speed * speed
    }

    /// The **total (stagnation) pressure** `p₀ = p + ½·ρ·|u|²` (Pa) at the centre of
    /// cell `(i, j)`, for fluid density `density` `ρ` (kg/m³) — the sum of the static
    /// pressure and the dynamic head, i.e. the pressure the flow would reach if brought
    /// reversibly to rest (Bernoulli's constant). It is the per-cell field that
    /// [`FlowSolution::total_pressure_range`] takes the max–min of, formed from the
    /// static [`FlowSolution::pressure`] and the dynamic
    /// [`FlowSolution::kinetic_energy_density_at_cell`]. In a steady inviscid flow it is
    /// constant along a streamline; its decline marks irreversible (viscous) loss.
    pub fn total_pressure_at_cell(&self, i: usize, j: usize, density: f64) -> f64 {
        self.pressure.at(i, j) + self.kinetic_energy_density_at_cell(i, j, density)
    }

    /// The **pressure coefficient** `Cp = (p − p_ref) / (½·ρ·U_ref²)` at the centre of
    /// cell `(i, j)` — the static pressure normalised by the reference dynamic head,
    /// the dimensionless pressure field of aerodynamic post-processing.
    /// `reference_pressure` `p_ref` (Pa) is the freestream static pressure, `density`
    /// `ρ` (kg/m³) and `reference_speed` `U_ref` (m/s) set the reference dynamic
    /// pressure `½ρU_ref²`.
    ///
    /// `Cp = 1` at a stagnation point (the flow brought fully to rest), `0` in the
    /// undisturbed freestream, and negative where the flow accelerates past `U_ref`
    /// (suction). For an incompressible inviscid flow it reduces to `Cp = 1 −
    /// (V/U_ref)²` by Bernoulli. Returns `0` if the reference dynamic pressure is not
    /// positive (`density` or `reference_speed` non-positive).
    pub fn pressure_coefficient_at_cell(
        &self,
        i: usize,
        j: usize,
        reference_pressure: f64,
        density: f64,
        reference_speed: f64,
    ) -> f64 {
        let q_ref = 0.5 * density * reference_speed * reference_speed;
        if q_ref <= 0.0 {
            return 0.0;
        }
        (self.pressure.at(i, j) - reference_pressure) / q_ref
    }

    /// Static-pressure range `Δp = p_max − p_min` (Pa) over the field — the
    /// total pressure variation. Gauge-independent (a difference), so it is
    /// meaningful despite the arbitrary absolute level: for channel flow it is
    /// essentially the streamwise pressure drop driving the throughflow; for
    /// the lid-driven cavity it is the swing that sustains the recirculation.
    /// Returns `0` for an empty grid.
    pub fn pressure_range(&self) -> f64 {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                let p = self.pressure.at(i, j);
                lo = lo.min(p);
                hi = hi.max(p);
            }
        }
        if hi >= lo {
            hi - lo
        } else {
            0.0
        }
    }

    /// The **area-averaged static pressure** `⟨p⟩ = (1/N)·Σ p_ij` (Pa) over the `N = nx·ny`
    /// cells — the mean of the static-pressure field, the pressure analogue of
    /// [`FlowSolution::mean_speed`]. Unlike the gauge-independent
    /// [`FlowSolution::pressure_range`] it carries the absolute level, so it shifts
    /// one-for-one with the (arbitrary) datum; it is the reference level the range swings
    /// about, and the spatial mean a net-force balance integrates. Returns `0` for an empty
    /// grid.
    pub fn mean_pressure(&self) -> f64 {
        let n = self.grid.nx * self.grid.ny;
        if n == 0 {
            return 0.0;
        }
        let mut sum = 0.0;
        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                sum += self.pressure.at(i, j);
            }
        }
        sum / n as f64
    }

    /// The **minimum static pressure** `min p_ij` (Pa) over the `nx·ny` cells — the
    /// suction-peak value, the lowest pressure anywhere in the field (where cavitation first
    /// onsets and the strongest suction load acts). It is the lower companion to the
    /// gauge-independent [`FlowSolution::pressure_range`] (`Δp = max − min`) and bounds the
    /// area-averaged [`FlowSolution::mean_pressure`] from below (`min ≤ ⟨p⟩`). Returns `0`
    /// for an empty grid.
    pub fn min_pressure(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny == 0 {
            return 0.0;
        }
        let mut lo = f64::INFINITY;
        for j in 0..ny {
            for i in 0..nx {
                lo = lo.min(self.pressure.at(i, j));
            }
        }
        lo
    }

    /// The **maximum static pressure** `max p_ij` (Pa) over the `nx·ny` cells — the
    /// stagnation-peak value, the highest pressure anywhere in the field (where the flow
    /// stagnates against the body). It is the upper companion to [`FlowSolution::min_pressure`];
    /// their difference is exactly the gauge-independent [`FlowSolution::pressure_range`]
    /// (`max − min = Δp`), and it bounds the area-averaged [`FlowSolution::mean_pressure`] from
    /// above (`⟨p⟩ ≤ max`). Returns `0` for an empty grid.
    pub fn max_pressure(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny == 0 {
            return 0.0;
        }
        let mut hi = f64::NEG_INFINITY;
        for j in 0..ny {
            for i in 0..nx {
                hi = hi.max(self.pressure.at(i, j));
            }
        }
        hi
    }

    /// The **spatial standard deviation of the static-pressure field**
    /// `σ_p = √(⟨(p − ⟨p⟩)²⟩)` (Pa) — the root-mean-square fluctuation of the cell
    /// pressures about their area-average [`FlowSolution::mean_pressure`], a measure of
    /// how non-uniform the pressure field is (the static-pressure analogue of the
    /// [`FlowSolution::speed_std_dev`] and the acoustic/turbulent `p_rms`). It is `0` for
    /// a perfectly uniform field and **gauge-independent** — shifting every pressure by a
    /// constant datum leaves it unchanged, like the [`FlowSolution::pressure_range`].
    /// Returns `0` for an empty grid.
    pub fn pressure_std_dev(&self) -> f64 {
        let n = self.grid.nx * self.grid.ny;
        if n == 0 {
            return 0.0;
        }
        let mean = self.mean_pressure();
        let mut sum_sq = 0.0;
        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                let d = self.pressure.at(i, j) - mean;
                sum_sq += d * d;
            }
        }
        (sum_sq / n as f64).sqrt()
    }

    /// The **root-mean-square static pressure** `p_rms = √(⟨p²⟩)` (Pa) over the `nx·ny`
    /// cells — the quadratic-mean pressure magnitude, the static-pressure analogue of the
    /// [`FlowSolution::rms_speed`]. By the variance identity it satisfies
    /// `p_rms² = ⟨p⟩² + σ_p²`, so it combines the area-averaged
    /// [`FlowSolution::mean_pressure`] with the fluctuation [`FlowSolution::pressure_std_dev`];
    /// it is `≥ |⟨p⟩|`, equal only for a uniform field, and unlike the signed mean it is
    /// always non-negative. Returns `0` for an empty grid.
    pub fn rms_pressure(&self) -> f64 {
        let n = self.grid.nx * self.grid.ny;
        if n == 0 {
            return 0.0;
        }
        let mut sum = 0.0;
        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                let p = self.pressure.at(i, j);
                sum += p * p;
            }
        }
        (sum / n as f64).sqrt()
    }

    /// The **total (stagnation) pressure range** `Δp₀ = p₀_max − p₀_min` (Pa) over
    /// the field, where `p₀ = p + ½ρ|u|²` is the total pressure (static plus the
    /// dynamic head). For an ideal inviscid flow Bernoulli's theorem keeps `p₀`
    /// constant along a streamline, so this variation is a direct read on the
    /// **total-pressure loss** — the irreversible (viscous) degradation a real
    /// flow suffers, the swing the static [`FlowSolution::pressure_range`] alone
    /// cannot see (static pressure trades *reversibly* against the dynamic head
    /// `½ρ|u|²`, but their sum only falls). `density` is the fluid density
    /// (kg/m³). Gauge-independent (a difference). Returns `0` for an empty grid.
    pub fn total_pressure_range(&self, density: f64) -> f64 {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                let speed = self.speed_at_cell(i, j);
                let p0 = self.pressure.at(i, j) + 0.5 * density * speed * speed;
                lo = lo.min(p0);
                hi = hi.max(p0);
            }
        }
        if hi >= lo {
            hi - lo
        } else {
            0.0
        }
    }

    /// The **area-averaged total (stagnation) pressure** `⟨p₀⟩ = ⟨p⟩ + ½ρ⟨|u|²⟩` (Pa)
    /// over the `nx·ny` cells — the mean of the per-cell total pressure
    /// `p₀ = p + ½ρ|u|²` ([`FlowSolution::total_pressure_at_cell`]), equal to the static
    /// mean [`FlowSolution::mean_pressure`] plus the mean dynamic head, which is exactly
    /// the area-averaged [`FlowSolution::mean_kinetic_energy_density`]. Where
    /// [`FlowSolution::mean_pressure`] tracks the field's net static level and
    /// [`FlowSolution::total_pressure_range`] its loss *spread*, this is the mean
    /// **mechanical-energy level**: in an ideal inviscid flow Bernoulli holds `p₀`
    /// constant, so a fall in `⟨p₀⟩` is the signature of irreversible (viscous) loss. It
    /// is bounded below by `⟨p⟩` (the dynamic head is non-negative), with equality only
    /// for a fluid at rest. `density` is the fluid density (kg/m³). Returns `0` for an
    /// empty grid.
    pub fn mean_total_pressure(&self, density: f64) -> f64 {
        self.mean_pressure() + self.mean_kinetic_energy_density(density)
    }

    /// Area-averaged velocity magnitude over the cell grid (m/s) — the typical
    /// flow speed, as opposed to the peak that `speed_at_cell` reveals. Returns
    /// `0` for an empty grid.
    pub fn mean_speed(&self) -> f64 {
        let n = self.grid.nx * self.grid.ny;
        if n == 0 {
            return 0.0;
        }
        let mut sum = 0.0;
        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                sum += self.speed_at_cell(i, j);
            }
        }
        sum / n as f64
    }

    /// The **root-mean-square speed** `u_rms = √(⟨|u|²⟩)` (m/s) — the square root of the
    /// cell-averaged squared speed. It is the velocity scale that carries the flow's
    /// kinetic energy: `½·ρ·u_rms²` equals the mean kinetic-energy density
    /// [`FlowSolution::mean_kinetic_energy_density`]. By the Cauchy–Schwarz inequality
    /// it never falls below the arithmetic-mean speed [`FlowSolution::mean_speed`] (the
    /// two are equal only for a perfectly uniform field) and never exceeds
    /// [`FlowSolution::max_speed`]. Returns `0` for an empty grid.
    pub fn rms_speed(&self) -> f64 {
        let n = self.grid.nx * self.grid.ny;
        if n == 0 {
            return 0.0;
        }
        let mut sum = 0.0;
        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                let speed = self.speed_at_cell(i, j);
                sum += speed * speed;
            }
        }
        (sum / n as f64).sqrt()
    }

    /// The **spatial standard deviation of the speed field** `σ_u = √(⟨|u|²⟩ − ⟨|u|⟩²)`
    /// (m/s) — the root-mean-square spread of the cell speeds about their area-average,
    /// a measure of how non-uniform the flow is (and the numerator of a
    /// turbulence-intensity estimate `σ_u / ⟨|u|⟩`). By the variance identity it is built
    /// from the [`FlowSolution::rms_speed`] and [`FlowSolution::mean_speed`] reducers; it
    /// is `0` for a perfectly uniform field and never exceeds the RMS speed. Returns `0`
    /// for an empty grid.
    pub fn speed_std_dev(&self) -> f64 {
        let variance = self.rms_speed().powi(2) - self.mean_speed().powi(2);
        variance.max(0.0).sqrt()
    }

    /// The **coefficient of variation of the speed field** `σ_u / ⟨|u|⟩` (dimensionless)
    /// — the spatial standard deviation [`FlowSolution::speed_std_dev`] normalised by the
    /// mean speed [`FlowSolution::mean_speed`], a scale-free measure of how non-uniform
    /// the flow is (the resolved-field analogue of a turbulence-intensity level, distinct
    /// from the modelled RANS intensity of [`crate::turbulence`]). It is `0` for a
    /// perfectly uniform field. Returns `0` for a quiescent or empty field (non-positive
    /// mean speed), so it never yields `NaN`.
    pub fn speed_coefficient_of_variation(&self) -> f64 {
        let mean = self.mean_speed();
        if mean <= 0.0 {
            return 0.0;
        }
        self.speed_std_dev() / mean
    }

    /// Peak velocity magnitude `max √(u²+v²)` (m/s) over the cell grid — the
    /// fastest the flow moves anywhere, the peak counterpart to the area-averaged
    /// [`FlowSolution::mean_speed`]. It is the speed that sets the convective CFL
    /// limit `U_max·Δt/Δx`, the maximum dynamic pressure `½ρ·U_max²`, and (in the
    /// lid-driven cavity) it cannot exceed the lid speed that drives the flow.
    /// Returns `0` for an empty grid.
    pub fn max_speed(&self) -> f64 {
        let mut peak = 0.0_f64;
        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                peak = peak.max(self.speed_at_cell(i, j));
            }
        }
        peak
    }

    /// The **minimum speed** `min √(u²+v²)` (m/s) over the cell grid — the slowest the flow
    /// moves anywhere, marking stagnation points and recirculation cores. It is the
    /// lower-bound counterpart to [`FlowSolution::max_speed`]; together they bracket the
    /// area-averaged [`FlowSolution::mean_speed`] (`min ≤ ⟨|u|⟩ ≤ max`). Returns `0` for an
    /// empty grid.
    pub fn min_speed(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny == 0 {
            return 0.0;
        }
        let mut slowest = f64::INFINITY;
        for j in 0..ny {
            for i in 0..nx {
                slowest = slowest.min(self.speed_at_cell(i, j));
            }
        }
        slowest
    }

    /// The **speed range** `Δ|u| = max√(u²+v²) − min√(u²+v²)` (m/s) over the cell grid —
    /// the spread between the fastest and slowest cell, the velocity analogue of the
    /// [`FlowSolution::pressure_range`]. It is a gauge-independent measure of how
    /// non-uniform the flow is: `0` for plug / uniform flow, growing as the field
    /// develops shear, recirculation, or stagnation. The peak-to-trough companion to the
    /// individual [`FlowSolution::max_speed`] and [`FlowSolution::min_speed`] (whose
    /// difference it equals). Returns `0` for an empty grid.
    pub fn speed_range(&self) -> f64 {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                let s = self.speed_at_cell(i, j);
                lo = lo.min(s);
                hi = hi.max(s);
            }
        }
        if hi >= lo {
            hi - lo
        } else {
            0.0
        }
    }

    /// The cell-centre location `(x, y)` (m) of the **fastest cell** — the *where*
    /// of the peak speed, the locational companion to [`FlowSolution::max_speed`]
    /// (which gives the magnitude). It marks the convective hot-spot that sets the
    /// CFL limit and, in a constriction, the lowest-pressure point by Bernoulli.
    /// Like [`FlowSolution::min_pressure_location`] it scans every cell-centred
    /// speed `√(u²+v²)` and returns the centre of the argmax cell. `None` for an
    /// empty grid.
    pub fn max_speed_location(&self) -> Option<(f64, f64)> {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny == 0 {
            return None;
        }
        let (dx, dy) = (self.grid.dx(), self.grid.dy());
        let mut fastest = f64::NEG_INFINITY;
        let mut loc = (0.0, 0.0);
        for j in 0..ny {
            for i in 0..nx {
                let speed = self.speed_at_cell(i, j);
                if speed > fastest {
                    fastest = speed;
                    loc = ((i as f64 + 0.5) * dx, (j as f64 + 0.5) * dy);
                }
            }
        }
        Some(loc)
    }

    /// The **convective Courant–Friedrichs–Lewy (CFL) number** `C = U_max·Δt /
    /// min(Δx, Δy)` for an explicit time step `dt` (s) — how many cells the fastest
    /// parcel crosses in one step, the stability gauge of an explicit advection
    /// scheme (which needs `C ≤ 1`; at `C = 1` a parcel traverses exactly one cell per
    /// step). It is built from [`FlowSolution::max_speed`] (whose doc notes it "sets
    /// the convective CFL limit") and the finer grid spacing. Returns `0` for a
    /// non-positive grid spacing or a non-finite / negative `dt`.
    pub fn convective_cfl_number(&self, dt: f64) -> f64 {
        let h = self.grid.dx().min(self.grid.dy());
        if h <= 0.0 || !dt.is_finite() || dt < 0.0 {
            return 0.0;
        }
        self.max_speed() * dt / h
    }

    /// The **diffusive (von Neumann) stability number** `d = ν·Δt / min(Δx, Δy)²` for
    /// kinematic viscosity `kinematic_viscosity` `ν` (m²/s) and explicit time step `dt`
    /// (s) — the viscous-diffusion counterpart of the [`convective_cfl_number`]. An
    /// explicit diffusion update is stable only for `d ≤ 1/4` in 2-D, so together with
    /// the convective CFL it bounds the stable time step. Their ratio is the cell
    /// Reynolds number `CFL/d = U_max·h/ν`. Returns `0` for a non-positive grid
    /// spacing, or non-finite / negative `dt` or `ν`.
    pub fn diffusive_number(&self, kinematic_viscosity: f64, dt: f64) -> f64 {
        let h = self.grid.dx().min(self.grid.dy());
        if h <= 0.0
            || !dt.is_finite()
            || dt < 0.0
            || !kinematic_viscosity.is_finite()
            || kinematic_viscosity < 0.0
        {
            return 0.0;
        }
        kinematic_viscosity * dt / (h * h)
    }

    /// The **cell (grid) Reynolds number** `Re_h = U_max·h / ν` (dimensionless), with
    /// `h = min(Δx, Δy)` and `U_max` the peak speed [`max_speed`] — the ratio of
    /// convective to diffusive transport across one grid cell, for kinematic viscosity
    /// `kinematic_viscosity` `ν` (m²/s). It is the [`convective_cfl_number`] to
    /// [`diffusive_number`] ratio (`Re_h = CFL/d`, independent of the time step), and a
    /// central-difference convection scheme stays free of spurious node-to-node
    /// oscillations only while `Re_h ≤ 2`. Returns `0` for a non-positive grid spacing
    /// or non-finite / non-positive `ν`.
    pub fn cell_reynolds_number(&self, kinematic_viscosity: f64) -> f64 {
        let h = self.grid.dx().min(self.grid.dy());
        if h <= 0.0 || !kinematic_viscosity.is_finite() || kinematic_viscosity <= 0.0 {
            return 0.0;
        }
        self.max_speed() * h / kinematic_viscosity
    }

    /// The **domain (flow) Reynolds number** `Re_L = U·Lₓ / ν` (dimensionless) — the peak
    /// speed `U` ([`FlowSolution::max_speed`]) times the streamwise domain length
    /// `Lₓ = nx·dx`, divided by the kinematic viscosity `kinematic_viscosity` `ν` (m²/s).
    /// Unlike the [`FlowSolution::cell_reynolds_number`] (a per-cell numerical-resolution
    /// measure on the grid spacing `h`), this is the *physical* flow-regime classifier:
    /// `Re_L ≲ 5×10⁵` is laminar over a flat plate, above it transitions to turbulence.
    /// Returns `0` for non-physical viscosity (`ν ≤ 0` or non-finite).
    pub fn domain_reynolds_number(&self, kinematic_viscosity: f64) -> f64 {
        if !kinematic_viscosity.is_finite() || kinematic_viscosity <= 0.0 {
            return 0.0;
        }
        let lx = (self.grid.nx as f64) * self.grid.dx();
        self.max_speed() * lx / kinematic_viscosity
    }

    /// The **flow-through (residence) time** `τ = Lₓ / U` (s) — the streamwise domain
    /// length `Lₓ = nx·dx` divided by the peak speed `U` ([`FlowSolution::max_speed`]), the
    /// convective time scale over which the flow crosses the domain. A simulation is
    /// typically advanced (and time-averaged) over several flow-through times to reach a
    /// statistically stationary state, so it sets the natural convergence / averaging
    /// window — distinct from the per-step numerical [`FlowSolution::convective_cfl_number`].
    /// Returns `0` for a quiescent field (`U ≤ 0`), where the residence time is unbounded.
    pub fn flow_through_time(&self) -> f64 {
        let u = self.max_speed();
        if u <= 0.0 {
            return 0.0;
        }
        let lx = (self.grid.nx as f64) * self.grid.dx();
        lx / u
    }

    /// Peak vorticity magnitude `|∂v/∂x − ∂u/∂y|` (1/s) over the interior cells,
    /// from central differences of the cell-centred velocity — the strongest
    /// local rotation in the flow (the lid-driven cavity's central vortex, a
    /// channel's wall shear layers). Returns `0` for a grid too small to take an
    /// interior central difference (`nx < 3` or `ny < 3`).
    pub fn max_vorticity(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx < 3 || ny < 3 {
            return 0.0;
        }
        let two_dx = 2.0 * self.grid.dx();
        let two_dy = 2.0 * self.grid.dy();
        let mut peak = 0.0_f64;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let dv_dx = (self.v_at_cell(i + 1, j) - self.v_at_cell(i - 1, j)) / two_dx;
                let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / two_dy;
                peak = peak.max((dv_dx - du_dy).abs());
            }
        }
        peak
    }

    /// The cell-centre location `(x, y)` (m) of the peak vorticity magnitude —
    /// the vortex core (e.g. the centre of the lid-driven cavity's primary
    /// vortex). Uses the same interior central-difference stencil as
    /// [`FlowSolution::max_vorticity`]. `None` for a grid too small to take an
    /// interior difference (`nx < 3` or `ny < 3`).
    pub fn peak_vorticity_location(&self) -> Option<(f64, f64)> {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx < 3 || ny < 3 {
            return None;
        }
        let (dx, dy) = (self.grid.dx(), self.grid.dy());
        let (two_dx, two_dy) = (2.0 * dx, 2.0 * dy);
        let mut peak = -1.0_f64;
        let mut loc = (0.0, 0.0);
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let dv_dx = (self.v_at_cell(i + 1, j) - self.v_at_cell(i - 1, j)) / two_dx;
                let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / two_dy;
                let w = (dv_dx - du_dy).abs();
                if w > peak {
                    peak = w;
                    loc = ((i as f64 + 0.5) * dx, (j as f64 + 0.5) * dy);
                }
            }
        }
        Some(loc)
    }

    /// The cell-centre location `(x, y)` (m) of the **lowest static pressure** —
    /// the *suction peak*: the cavitation-inception site and the low-pressure
    /// core of a vortex (the heart of the lid-driven cavity's primary vortex, or
    /// the suction side of a bluff-body wake). Where
    /// [`FlowSolution::pressure_range`] gives the *magnitude* of the pressure
    /// swing, this gives the *position* of its minimum — the locational reading
    /// none of the field-magnitude diagnostics expose. Unlike the vorticity-core
    /// locator it needs no interior stencil — pressure is stored cell-centred
    /// everywhere — so it scans every cell. `None` for an empty grid.
    pub fn min_pressure_location(&self) -> Option<(f64, f64)> {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny == 0 {
            return None;
        }
        let (dx, dy) = (self.grid.dx(), self.grid.dy());
        let mut lowest = f64::INFINITY;
        let mut loc = (0.0, 0.0);
        for j in 0..ny {
            for i in 0..nx {
                let p = self.pressure.at(i, j);
                if p < lowest {
                    lowest = p;
                    loc = ((i as f64 + 0.5) * dx, (j as f64 + 0.5) * dy);
                }
            }
        }
        Some(loc)
    }

    /// The cell-centre location `(x, y)` (m) of the **highest static pressure** —
    /// the *stagnation point*: where the flow is brought to rest and the Bernoulli
    /// pressure peaks (the windward face of a bluff body, or the downstream corner
    /// where a wall turns the flow). It is the high-pressure companion to
    /// [`FlowSolution::min_pressure_location`] (the suction peak / vortex core):
    /// together they bracket where the pressure field does its work, and where
    /// [`FlowSolution::pressure_range`] gives the *magnitude* of the swing this
    /// gives the *position* of its maximum. Like the minimum locator it needs no
    /// interior stencil — pressure is stored cell-centred everywhere — so it scans
    /// every cell, and like it the result is gauge-invariant (adding a constant to
    /// the whole field leaves the argmax fixed). `None` for an empty grid.
    pub fn max_pressure_location(&self) -> Option<(f64, f64)> {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny == 0 {
            return None;
        }
        let (dx, dy) = (self.grid.dx(), self.grid.dy());
        let mut highest = f64::NEG_INFINITY;
        let mut loc = (0.0, 0.0);
        for j in 0..ny {
            for i in 0..nx {
                let p = self.pressure.at(i, j);
                if p > highest {
                    highest = p;
                    loc = ((i as f64 + 0.5) * dx, (j as f64 + 0.5) * dy);
                }
            }
        }
        Some(loc)
    }

    /// Net volumetric flow rate through the **inlet** (left, `x = 0`) face, per
    /// unit depth (m²/s in 2-D): `Σ_j u(0, j)·dy`, the discrete `∫ u dy` along
    /// the inlet. Positive for fluid entering the domain; ~0 for an enclosed
    /// case (the lid-driven cavity).
    pub fn inlet_flow_rate(&self) -> f64 {
        let dy = self.grid.dy();
        (0..self.grid.ny).map(|j| self.u.at(0, j) * dy).sum()
    }

    /// Net volumetric flow rate through the **outlet** (right, `x = lx`) face,
    /// per unit depth (m²/s in 2-D): `Σ_j u(nx, j)·dy`. For an incompressible
    /// steady solution it should equal [`FlowSolution::inlet_flow_rate`]; the
    /// gap is the global continuity residual.
    pub fn outlet_flow_rate(&self) -> f64 {
        let dy = self.grid.dy();
        (0..self.grid.ny)
            .map(|j| self.u.at(self.grid.nx, j) * dy)
            .sum()
    }

    /// Relative mass-continuity error `|Q_out − Q_in| / |Q_in|` between the
    /// outlet and inlet volumetric flow rates — a global check that the
    /// incompressible solver conserves mass across the domain (a well-converged
    /// through-flow drives this toward zero). Returns `0` when there is no net
    /// inflow (an enclosed case), where the relative measure is undefined.
    pub fn continuity_error(&self) -> f64 {
        let q_in = self.inlet_flow_rate();
        if q_in.abs() < 1e-12 {
            0.0
        } else {
            (self.outlet_flow_rate() - q_in).abs() / q_in.abs()
        }
    }

    /// The bulk (mean-throughflow) velocity `U_bulk = Q_in / H` (m/s) — the inlet
    /// volumetric flow rate per unit channel height, the reference speed behind
    /// the bulk Reynolds number. For fully-developed plane Poiseuille flow it is
    /// exactly ⅔ of the centreline peak. Returns `0` for a zero-height grid.
    pub fn bulk_velocity(&self) -> f64 {
        if self.grid.ly > 0.0 {
            self.inlet_flow_rate() / self.grid.ly
        } else {
            0.0
        }
    }

    /// The **displacement thickness** `δ* = ∫₀^Ly (1 − u/u_ref) dy` (m) at the
    /// outlet cross-section (`x = lx`) — the distance the wall would have to be
    /// displaced into an inviscid stream of edge speed `u_ref` to impose the same
    /// mass-flow deficit the viscous boundary layer does. `u_ref` is the
    /// edge / free-stream (or channel-centreline) reference speed, typically the
    /// [`FlowSolution::bulk_velocity`] or a known inlet speed. By construction it
    /// obeys the mass-deficit identity `δ* = Ly − Q_out/u_ref` against the
    /// [`FlowSolution::outlet_flow_rate`] `Q_out`: the channel height minus the
    /// height an inviscid `u_ref` stream would need to carry the actual outflow.
    /// Plug flow (uniform `u = u_ref`) gives `0`; the more the near-wall fluid
    /// lags, the larger `δ*`. Returns `0` for a non-finite or non-positive `u_ref`
    /// or an empty grid.
    pub fn displacement_thickness(&self, u_ref: f64) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny == 0 || !u_ref.is_finite() || u_ref <= 0.0 {
            return 0.0;
        }
        let dy = self.grid.dy();
        (0..ny).map(|j| (1.0 - self.u.at(nx, j) / u_ref) * dy).sum()
    }

    /// The **momentum thickness** `θ = ∫₀^Ly (u/u_ref)(1 − u/u_ref) dy` (m) at the
    /// outlet cross-section (`x = lx`) — the height of an inviscid stream of edge
    /// speed `u_ref` whose momentum flux equals the momentum *deficit* the viscous
    /// boundary layer carries. `u_ref` is the edge / free-stream (or channel-
    /// centreline) reference speed (the [`FlowSolution::bulk_velocity`] or a known
    /// inlet speed), matching [`FlowSolution::displacement_thickness`]. It is the
    /// momentum-deficit companion to that mass-deficit `δ*`, and the two define the
    /// shape factor `H = δ*/θ` (≈ 2.6 for a laminar Blasius layer, falling toward
    /// ~1.3–1.4 as a profile fills out turbulently). Since `δ* − θ = ∫(1 − u/u_ref)²
    /// dy ≥ 0`, always `θ ≤ δ*`. Plug flow (uniform `u = u_ref`) gives `0`; the more
    /// the near-wall fluid lags, the larger `θ`. Returns `0` for a non-finite or
    /// non-positive `u_ref` or an empty grid.
    pub fn momentum_thickness(&self, u_ref: f64) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny == 0 || !u_ref.is_finite() || u_ref <= 0.0 {
            return 0.0;
        }
        let dy = self.grid.dy();
        (0..ny)
            .map(|j| {
                let r = self.u.at(nx, j) / u_ref;
                r * (1.0 - r) * dy
            })
            .sum()
    }

    /// The **energy thickness** `δ_E = ∫₀^Ly (u/u_ref)(1 − (u/u_ref)²) dy` (m) at the
    /// outlet cross-section (`x = lx`) — the height of an inviscid stream of edge speed
    /// `u_ref` whose *kinetic-energy* flux equals the energy deficit the viscous boundary
    /// layer carries. It is the third of the boundary-layer integral thicknesses,
    /// completing the family with the mass-deficit [`FlowSolution::displacement_thickness`]
    /// `δ*` and the momentum-deficit [`FlowSolution::momentum_thickness`] `θ`; their ratio
    /// `δ_E/θ` is the energy shape factor `H*` (≈ 1.5 for a uniform-defect or linear
    /// profile, ≈ 1.7 for a laminar Blasius layer). Since `δ_E − θ = ∫(u/u_ref)²(1 −
    /// u/u_ref) dy ≥ 0` for a profile bounded by `u_ref`, always `δ_E ≥ θ`. Plug flow
    /// (uniform `u = u_ref`) carries no deficit and gives `0`. Returns `0` for a
    /// non-finite or non-positive `u_ref` or an empty grid.
    pub fn energy_thickness(&self, u_ref: f64) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny == 0 || !u_ref.is_finite() || u_ref <= 0.0 {
            return 0.0;
        }
        let dy = self.grid.dy();
        (0..ny)
            .map(|j| {
                let r = self.u.at(nx, j) / u_ref;
                r * (1.0 - r * r) * dy
            })
            .sum()
    }

    /// The boundary-layer **shape factor** `H = δ*/θ` (dimensionless) — the ratio of
    /// the [`FlowSolution::displacement_thickness`] to the
    /// [`FlowSolution::momentum_thickness`] at the outlet, evaluated against edge
    /// speed `u_ref`. It is the single number that classifies how *full* the velocity
    /// profile is: `H ≈ 2.6` for a laminar Blasius layer, falling toward `~1.3–1.4`
    /// as the profile fills out turbulently, and rising toward `~3.5–4` as the
    /// near-wall flow stalls — the classic precursor of **separation**. Because
    /// `δ* ≥ θ` always, `H ≥ 1`. It completes the boundary-layer integral-thickness
    /// family with `δ*` and `θ`. Returns `0` when the momentum thickness is `0` (plug
    /// flow, or a non-finite / non-positive `u_ref` or empty grid, where both
    /// thicknesses vanish).
    pub fn shape_factor(&self, u_ref: f64) -> f64 {
        let theta = self.momentum_thickness(u_ref);
        if theta.abs() < 1e-12 {
            return 0.0;
        }
        self.displacement_thickness(u_ref) / theta
    }

    /// The **energy shape factor** `H* = δ_E/θ` (dimensionless) — the ratio of the
    /// boundary-layer [`FlowSolution::energy_thickness`] `δ_E` to the
    /// [`FlowSolution::momentum_thickness`] `θ` at the outlet, against edge speed `u_ref`.
    /// It is the *kinetic-energy* analog of the mass/momentum [`FlowSolution::shape_factor`]
    /// `H = δ*/θ`, and the quantity advanced by boundary-layer integral methods (e.g. Head's
    /// entrainment closure). It measures how full the profile is in an energy sense: `H* ≈
    /// 1.5` for a uniform-defect or linear profile, rising toward `≈ 1.6–1.7` for a laminar
    /// Blasius layer. Because `δ_E ≥ θ` always (`δ_E − θ = ∫(u/u_ref)²(1 − u/u_ref) dy ≥ 0`),
    /// `H* ≥ 1`. Returns `0` when the momentum thickness is `0` (plug flow, or a non-finite
    /// / non-positive `u_ref` or empty grid, where both thicknesses vanish).
    pub fn energy_shape_factor(&self, u_ref: f64) -> f64 {
        let theta = self.momentum_thickness(u_ref);
        if theta.abs() < 1e-12 {
            return 0.0;
        }
        self.energy_thickness(u_ref) / theta
    }

    /// The **momentum-thickness Reynolds number** `Re_θ = u_ref·θ / ν` (dimensionless) —
    /// the Reynolds number formed on the boundary-layer [`FlowSolution::momentum_thickness`]
    /// `θ` and the edge speed `u_ref`, with kinematic viscosity `kinematic_viscosity` `ν`
    /// (m²/s). It is the standard *transition* parameter of a boundary layer: a flat-plate
    /// layer trips from laminar to turbulent near `Re_θ ≈ 350–400` and is fully turbulent
    /// by `~ 1000`, so it is the single number that decides a station's regime. `u_ref` is
    /// the same edge / free-stream reference speed the thickness integrals use (the
    /// [`FlowSolution::bulk_velocity`] or a known inlet speed). Distinct from the
    /// geometry-based [`FlowSolution::domain_reynolds_number`] — this is built on the
    /// boundary layer's own momentum thickness. Returns `0` for non-positive or non-finite
    /// `u_ref` or viscosity, or whenever the momentum thickness vanishes (plug flow).
    pub fn momentum_thickness_reynolds_number(&self, u_ref: f64, kinematic_viscosity: f64) -> f64 {
        if !u_ref.is_finite()
            || u_ref <= 0.0
            || !kinematic_viscosity.is_finite()
            || kinematic_viscosity <= 0.0
        {
            return 0.0;
        }
        u_ref * self.momentum_thickness(u_ref) / kinematic_viscosity
    }

    /// Area-averaged kinetic-energy density `⟨½ρ|u|²⟩` (J/m³ = Pa) over the cell
    /// grid — the mean specific kinetic energy of the flow, in the same units as
    /// (and directly comparable to) the freestream dynamic pressure `½ρU²`.
    /// `density` is the fluid density (kg/m³). Returns `0` for an empty grid.
    pub fn mean_kinetic_energy_density(&self, density: f64) -> f64 {
        let n = self.grid.nx * self.grid.ny;
        if n == 0 {
            return 0.0;
        }
        let mut sum = 0.0;
        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                let speed = self.speed_at_cell(i, j);
                sum += 0.5 * density * speed * speed;
            }
        }
        sum / n as f64
    }

    /// The **total kinetic energy** of the flow `∫ ½ρ|u|² dA = ⟨KE density⟩ · A` (J per
    /// unit depth in 2-D) — the extensive companion to the area-averaged
    /// [`FlowSolution::mean_kinetic_energy_density`], obtained by multiplying it by the
    /// domain area `A = Lx·Ly`. With density `density` `ρ` (kg/m³) it is the energy budget
    /// that decays under viscous dissipation; equivalently `½·ρ·u_rms²·A`. Returns `0` for
    /// an empty grid.
    pub fn total_kinetic_energy(&self, density: f64) -> f64 {
        let area = (self.grid.nx as f64) * self.grid.dx() * (self.grid.ny as f64) * self.grid.dy();
        self.mean_kinetic_energy_density(density) * area
    }

    /// The **peak dynamic pressure** `q = ½·ρ·U_max²` (Pa) — half the density `density`
    /// `ρ` (kg/m³) times the square of the peak speed `U_max` ([`FlowSolution::max_speed`]).
    /// It is the reference pressure aerodynamic loads are normalised against — lift and
    /// drag scale as `q·S·C`, and the pressure coefficient is `Cp = (p − p∞)/q` — and it is
    /// the kinetic-energy density evaluated at the fastest point. It therefore bounds the
    /// area-averaged [`FlowSolution::mean_kinetic_energy_density`] from above (equal only
    /// when the flow is uniform), since `U_max² ≥ ⟨|u|²⟩`.
    pub fn dynamic_pressure(&self, density: f64) -> f64 {
        let u_max = self.max_speed();
        0.5 * density * u_max * u_max
    }

    /// The **Euler number** `Eu = Δp / (ρ·U²)` (dimensionless) — the field's static
    /// pressure variation [`FlowSolution::pressure_range`] `Δp = p_max − p_min`
    /// normalised by the dynamic head `ρ·U_max²` (twice the peak
    /// [`FlowSolution::dynamic_pressure`] `q = ½ρU²`). It measures how large the
    /// pressure forces are relative to the inertial ones: the global, single-scalar
    /// companion to the *local* per-cell pressure coefficient
    /// [`FlowSolution::pressure_coefficient_at_cell`], and the natural figure of merit
    /// for internal (duct / pump / cavity) flow where the pressure *drop* — not a
    /// free-stream reference — is what matters.
    ///
    /// `Eu = pressure_range / (2·dynamic_pressure)`. Independent of the pressure datum
    /// (it uses the gauge-independent range), it scales as `1/ρ` and as `1/U²`. Returns
    /// `0` for a non-positive or non-finite density, or when the dynamic head vanishes
    /// (a quiescent field, `U_max = 0`).
    pub fn euler_number(&self, density: f64) -> f64 {
        if !density.is_finite() || density <= 0.0 {
            return 0.0;
        }
        let q = self.dynamic_pressure(density);
        if q <= 0.0 {
            return 0.0;
        }
        self.pressure_range() / (2.0 * q)
    }

    /// The mean wall-normal velocity gradient `⟨|∂u/∂y|⟩` at the bottom wall
    /// (1/s) — estimated one-sidedly from the no-slip wall (`u = 0`) to the first
    /// cell centre, `2·u_cell(i, 0)/dy`, averaged over the streamwise cells.
    /// Multiplying by the dynamic viscosity `μ = ρν` gives the wall shear stress
    /// `τ_w` — the skin friction the wall exerts on the flow (and that the
    /// streamwise pressure drop must overcome in a channel). Returns `0` for an
    /// empty grid.
    pub fn bottom_wall_shear_rate(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny == 0 {
            return 0.0;
        }
        let dy = self.grid.dy();
        let mut sum = 0.0;
        for i in 0..nx {
            sum += (2.0 * self.u_at_cell(i, 0) / dy).abs();
        }
        sum / nx as f64
    }

    /// The mean wall-normal velocity gradient `⟨|∂u/∂y|⟩` at the **top** wall
    /// (`y = ly`) (1/s) — the one-sided estimate `2·u_cell(i, ny−1)/dy` from the
    /// no-slip wall to the first cell centre, averaged over the streamwise cells:
    /// the companion to [`FlowSolution::bottom_wall_shear_rate`]. Multiplying by
    /// the dynamic viscosity `μ = ρν` gives the top-wall shear stress `τ_w`. In the
    /// **lid-driven cavity** the top wall is the moving lid, so this is the drag
    /// the lid feels — the skin friction through which it pumps power into the
    /// flow; in a channel it is the second, symmetric wall. Returns `0` for an
    /// empty grid.
    pub fn top_wall_shear_rate(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny == 0 {
            return 0.0;
        }
        let dy = self.grid.dy();
        let mut sum = 0.0;
        for i in 0..nx {
            sum += (2.0 * self.u_at_cell(i, ny - 1) / dy).abs();
        }
        sum / nx as f64
    }

    /// The **friction velocity** `u_τ = √(ν · ⟨|∂u/∂y|⟩_wall)` (m/s) at the bottom
    /// wall — the characteristic velocity scale of wall turbulence, formed from the
    /// kinematic viscosity `kinematic_viscosity` `ν` (m²/s) and the bottom-wall shear
    /// rate [`FlowSolution::bottom_wall_shear_rate`] (equivalently `u_τ = √(τ_w/ρ)`).
    /// It is the velocity that non-dimensionalises the near-wall region: the wall
    /// coordinate `y⁺ = u_τ·y/ν` and the law-of-the-wall velocity `u⁺ = u/u_τ` are
    /// both built on it, and it sets the inner-layer scaling of any wall-bounded
    /// turbulent flow. Unlike the shear *rate* (1/s) it is a *velocity* (m/s).
    ///
    /// This is the wall-**resolved** `u_τ`, read straight from the computed wall
    /// gradient (the DNS / wall-resolved route); the wall-function turbulence model
    /// (`turbulence::WallFunction::friction_velocity`) gives the complementary
    /// wall-**modelled** `u_τ`, solved from a near-wall velocity sample via the log
    /// law when the near-wall region is left unresolved. Returns `0` for a
    /// non-positive or non-finite viscosity, or an empty grid / zero wall shear.
    pub fn friction_velocity(&self, kinematic_viscosity: f64) -> f64 {
        if !kinematic_viscosity.is_finite() || kinematic_viscosity <= 0.0 {
            return 0.0;
        }
        (kinematic_viscosity * self.bottom_wall_shear_rate()).sqrt()
    }

    /// The **wall shear stress** `τ_w = μ · ⟨|∂u/∂y|⟩_wall` (Pa) at the bottom wall —
    /// the dynamic viscosity `dynamic_viscosity` `μ` (Pa·s) times the bottom-wall shear
    /// rate [`FlowSolution::bottom_wall_shear_rate`] `∂u/∂y` (1/s). This is the
    /// *dimensional* skin-friction traction the wall exerts on the flow (force per unit
    /// wall area) — the quantity the streamwise pressure drop must overcome to push
    /// fluid through a channel.
    ///
    /// It is distinct from the shear *rate* (1/s) it is built on, and is the
    /// stress-domain partner of the [`FlowSolution::friction_velocity`] `u_τ`: since
    /// `u_τ = √(ν · ∂u/∂y)` with `ν = μ/ρ`, the two satisfy `τ_w = ρ · u_τ²` exactly.
    /// Non-dimensionalising `τ_w` by the dynamic pressure `½ρU²` gives the
    /// skin-friction coefficient `C_f`. Returns `0` for a non-positive or non-finite
    /// viscosity, or an empty grid / zero wall shear.
    pub fn bottom_wall_shear_stress(&self, dynamic_viscosity: f64) -> f64 {
        if !dynamic_viscosity.is_finite() || dynamic_viscosity <= 0.0 {
            return 0.0;
        }
        dynamic_viscosity * self.bottom_wall_shear_rate()
    }

    /// The **skin-friction coefficient** `C_f = τ_w / (½ρU²)` (dimensionless) at the
    /// bottom wall — the wall shear stress [`FlowSolution::bottom_wall_shear_stress`]
    /// `τ_w` normalised by the dynamic pressure `½ρU²` built from the density `density`
    /// `ρ` (kg/m³) and the reference speed `U` (the bulk through-flow
    /// [`FlowSolution::bulk_velocity`]). It is the dimensionless headline of the
    /// wall-friction family: where `τ_w` (Pa) and `u_τ` (m/s) carry units, `C_f`
    /// collapses the wall drag onto the universal scale every boundary-layer
    /// skin-friction correlation (`C_f = 0.664/√Re_x` laminar, the `1/7`-power law
    /// turbulent) is plotted on.
    ///
    /// Equivalently `C_f = 2·(u_τ/U)²` through the friction velocity
    /// [`FlowSolution::friction_velocity`] (since `τ_w = ρ·u_τ²`). Returns `0` for a
    /// non-positive or non-finite density, or when the reference dynamic pressure
    /// vanishes (a quiescent field, `U = 0`); the viscosity guard is inherited from
    /// `bottom_wall_shear_stress`.
    pub fn skin_friction_coefficient(&self, dynamic_viscosity: f64, density: f64) -> f64 {
        if !density.is_finite() || density <= 0.0 {
            return 0.0;
        }
        let u_ref = self.bulk_velocity();
        let dyn_pressure = 0.5 * density * u_ref * u_ref;
        if dyn_pressure <= 0.0 {
            return 0.0;
        }
        self.bottom_wall_shear_stress(dynamic_viscosity) / dyn_pressure
    }

    /// The mean wall-tangential velocity gradient `⟨|∂v/∂x|⟩` at the **left** wall
    /// (`x = 0`) (1/s) — the one-sided estimate `2·v_cell(0, j)/dx` from the no-slip
    /// wall to the first cell centre, averaged over the wall-normal (vertical)
    /// cells: the *vertical-wall* companion to the horizontal-wall
    /// [`FlowSolution::bottom_wall_shear_rate`] / [`FlowSolution::top_wall_shear_rate`].
    /// Where those measure `∂u/∂y` on the top and bottom, this measures the
    /// tangential velocity `v` shearing past the side wall (`∂v/∂x`); multiplying by
    /// the dynamic viscosity `μ = ρν` gives the side-wall shear stress `τ_w`. In the
    /// **lid-driven cavity** the side walls carry the vertical drag of the
    /// recirculating vortex. Returns `0` for an empty grid.
    pub fn left_wall_shear_rate(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny == 0 {
            return 0.0;
        }
        let dx = self.grid.dx();
        let mut sum = 0.0;
        for j in 0..ny {
            sum += (2.0 * self.v_at_cell(0, j) / dx).abs();
        }
        sum / ny as f64
    }

    /// The mean wall-tangential velocity gradient `⟨|∂v/∂x|⟩` at the **right** wall
    /// (`x = lx`) (1/s) — the one-sided estimate `2·v_cell(nx−1, j)/dx` from the
    /// no-slip wall to the last cell centre, averaged over the wall-normal (vertical)
    /// cells: the companion to [`FlowSolution::left_wall_shear_rate`] on the opposite
    /// side, completing the four-wall skin-friction set with the horizontal-wall
    /// [`FlowSolution::bottom_wall_shear_rate`] / [`FlowSolution::top_wall_shear_rate`].
    /// Multiplying by the dynamic viscosity `μ = ρν` gives the right-wall shear
    /// stress `τ_w`. Returns `0` for an empty grid.
    pub fn right_wall_shear_rate(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny == 0 {
            return 0.0;
        }
        let dx = self.grid.dx();
        let mut sum = 0.0;
        for j in 0..ny {
            sum += (2.0 * self.v_at_cell(nx - 1, j) / dx).abs();
        }
        sum / ny as f64
    }

    /// The **vorticity** `ω = ∂v/∂x − ∂u/∂y` (1/s) at the centre of cell `(i, j)` —
    /// the local signed rate of fluid rotation, the per-cell field that
    /// [`FlowSolution::circulation`] (`∫ω dA`) and [`FlowSolution::enstrophy`]
    /// (`½∫ω² dA`) integrate and that [`FlowSolution::max_vorticity`] takes the peak
    /// of. It complements the [`FlowSolution::u_at_cell`] / [`FlowSolution::v_at_cell`]
    /// / [`FlowSolution::speed_at_cell`] accessors, exposing the rotation a user would
    /// otherwise have to difference by hand. Computed with the same interior
    /// central-difference stencil as those integrals; **boundary cells** (the first or
    /// last row/column, where the centred stencil has no neighbour) return `0`,
    /// matching the interior-only domain of `circulation`/`enstrophy`. Positive `ω` is
    /// anticlockwise swirl; pure translation or a uniform field gives `0`, and
    /// solid-body rotation at rate `Ω` gives `2Ω`.
    pub fn vorticity_at_cell(&self, i: usize, j: usize) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if i == 0 || j == 0 || i + 1 >= nx || j + 1 >= ny {
            return 0.0;
        }
        let (dx, dy) = (self.grid.dx(), self.grid.dy());
        let dv_dx = (self.v_at_cell(i + 1, j) - self.v_at_cell(i - 1, j)) / (2.0 * dx);
        let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / (2.0 * dy);
        dv_dx - du_dy
    }

    /// The **velocity divergence** `∇·u = ∂u/∂x + ∂v/∂y` (1/s) at the centre of cell
    /// `(i, j)` — the per-cell field that [`FlowSolution::max_divergence`] takes the
    /// peak `|∇·u|` of, and the local mass-continuity residual the pressure-projection
    /// step drives to zero. It is the **divergence** companion to the **curl**
    /// [`FlowSolution::vorticity_at_cell`]: together `∇·u` and `ω = ∇×u` are the two
    /// fundamental first-order differential invariants of a 2-D velocity field, and a
    /// flow can carry either without the other — incompressible **solid-body rotation**
    /// (`u = −Ω·y`, `v = Ω·x`) has `∇·u = 0` yet `ω = 2Ω`, while a **radial source**
    /// (`u = a·x`, `v = a·y`) has `∇·u = 2a` yet `ω = 0`. Formed from the MAC face
    /// velocities exactly as the projection measures continuity —
    /// `(u_{i+1,j} − u_{i,j})/dx + (v_{i,j+1} − v_{i,j})/dy` — so unlike the centred
    /// vorticity stencil it is defined on **every** cell, boundaries included. A
    /// well-converged incompressible solution reads ≈ 0 everywhere; out-of-range
    /// indices return `0`.
    pub fn divergence_at_cell(&self, i: usize, j: usize) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if i >= nx || j >= ny {
            return 0.0;
        }
        let (dx, dy) = (self.grid.dx(), self.grid.dy());
        let du_dx = (self.u.at(i + 1, j) - self.u.at(i, j)) / dx;
        let dv_dy = (self.v.at(i, j + 1) - self.v.at(i, j)) / dy;
        du_dx + dv_dy
    }

    /// The **Q-criterion** `Q = ½(‖Ω‖² − ‖S‖²)` (1/s²) at the centre of cell
    /// `(i, j)` — the per-cell vortex-identification scalar (Hunt, Wray & Moin 1988)
    /// that [`FlowSolution::max_q_criterion`] takes the peak of. It completes the
    /// per-cell first-order field-accessor trio with the curl
    /// [`FlowSolution::vorticity_at_cell`] (`ω = ∂v/∂x − ∂u/∂y`) and the divergence
    /// [`FlowSolution::divergence_at_cell`] (`∇·u`): `Q` weighs the rotation tensor
    /// `Ω` against the strain-rate tensor `S`, so it is **positive only where the
    /// fluid genuinely swirls**. Solid-body rotation at rate `Ω` gives `Q = Ω²`
    /// (`= ω²/4`, since `ω = 2Ω`); a pure shear layer — equal rotation and strain —
    /// gives `Q = 0` despite its non-zero vorticity; and a pure straining flow gives
    /// `Q < 0`. Unlike the clamped `max_q_criterion`, this returns the **raw** `Q`
    /// (which may be negative). Computed with the same interior central-difference
    /// stencil as the other accessors; **boundary cells** (first/last row or column)
    /// return `0`.
    pub fn q_criterion_at_cell(&self, i: usize, j: usize) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if i == 0 || j == 0 || i + 1 >= nx || j + 1 >= ny {
            return 0.0;
        }
        let (dx, dy) = (self.grid.dx(), self.grid.dy());
        let du_dx = (self.u_at_cell(i + 1, j) - self.u_at_cell(i - 1, j)) / (2.0 * dx);
        let dv_dy = (self.v_at_cell(i, j + 1) - self.v_at_cell(i, j - 1)) / (2.0 * dy);
        let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / (2.0 * dy);
        let dv_dx = (self.v_at_cell(i + 1, j) - self.v_at_cell(i - 1, j)) / (2.0 * dx);
        let omega = dv_dx - du_dy;
        let s_off = 0.5 * (du_dy + dv_dx);
        let omega_norm_sq = 0.5 * omega * omega;
        let s_norm_sq = du_dx * du_dx + dv_dy * dv_dy + 2.0 * s_off * s_off;
        0.5 * (omega_norm_sq - s_norm_sq)
    }

    /// The **strain rate** `√(2·S:S)` (1/s) at the centre of cell `(i, j)` — the
    /// local *rate of deformation*, the symmetric (rate-of-strain `S`) part of the
    /// velocity gradient, completing the per-cell decomposition with the
    /// antisymmetric [`FlowSolution::vorticity_at_cell`] (spin `ω`), the trace
    /// [`FlowSolution::divergence_at_cell`] (dilatation `∇·u`), and their balance
    /// [`FlowSolution::q_criterion_at_cell`] (`Q`). It is the per-cell field that
    /// [`FlowSolution::max_strain_rate`] takes the peak of, with
    /// `S:S = u_x² + v_y² + ½(u_y+v_x)²`. A **solid-body rotation** (pure spin) reads
    /// `0` here, a **pure strain** reads `0` on the vorticity, and the three are
    /// locked by `Q = (ω² − strain²)/4`. Computed with the same interior
    /// central-difference stencil as the other accessors; **boundary cells**
    /// (first/last row or column) return `0`.
    pub fn strain_rate_at_cell(&self, i: usize, j: usize) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if i == 0 || j == 0 || i + 1 >= nx || j + 1 >= ny {
            return 0.0;
        }
        let (dx, dy) = (self.grid.dx(), self.grid.dy());
        let du_dx = (self.u_at_cell(i + 1, j) - self.u_at_cell(i - 1, j)) / (2.0 * dx);
        let dv_dy = (self.v_at_cell(i, j + 1) - self.v_at_cell(i, j - 1)) / (2.0 * dy);
        let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / (2.0 * dy);
        let dv_dx = (self.v_at_cell(i + 1, j) - self.v_at_cell(i - 1, j)) / (2.0 * dx);
        let s_off = 0.5 * (du_dy + dv_dx);
        let s_dd = du_dx * du_dx + dv_dy * dv_dy + 2.0 * s_off * s_off;
        (2.0 * s_dd).sqrt()
    }

    /// The total circulation `Γ = ∫ ω dA` (m²/s) over the interior cells — the
    /// signed net rotation of the flow (Kelvin's circulation theorem; by
    /// Kutta–Joukowski it ties to lift on an immersed body). Uses the same
    /// interior central-difference vorticity stencil as
    /// [`FlowSolution::max_vorticity`], summed and weighted by the cell area.
    /// Unlike the always-positive peak vorticity it carries a sign (the sense of
    /// the net swirl). Returns `0` for a grid too small for an interior
    /// difference (`nx < 3` or `ny < 3`).
    pub fn circulation(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx < 3 || ny < 3 {
            return 0.0;
        }
        let (dx, dy) = (self.grid.dx(), self.grid.dy());
        let (two_dx, two_dy) = (2.0 * dx, 2.0 * dy);
        let mut sum = 0.0;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let dv_dx = (self.v_at_cell(i + 1, j) - self.v_at_cell(i - 1, j)) / two_dx;
                let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / two_dy;
                sum += dv_dx - du_dy;
            }
        }
        sum * dx * dy
    }

    /// The **area-averaged interior vorticity** `⟨ω⟩ = (1/A)∫∫ω dA` (1/s) over the
    /// interior cells — the mean of the central-difference vorticity `ω = ∂v/∂x − ∂u/∂y`,
    /// the *net* solid-body rotation rate of the field. By Stokes' theorem it equals the
    /// [`FlowSolution::circulation`] per unit interior area (`⟨ω⟩ = Γ/A`); unlike the
    /// strictly-positive [`FlowSolution::enstrophy`] it is *signed*, so equal-and-opposite
    /// swirl cancels and a symmetric vortex pair averages to `0`. Uses the same interior
    /// stencil as [`FlowSolution::max_vorticity`]. Returns `0` for a grid too small for an
    /// interior difference (`nx < 3` or `ny < 3`).
    pub fn mean_vorticity(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx < 3 || ny < 3 {
            return 0.0;
        }
        let (two_dx, two_dy) = (2.0 * self.grid.dx(), 2.0 * self.grid.dy());
        let mut sum = 0.0;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let dv_dx = (self.v_at_cell(i + 1, j) - self.v_at_cell(i - 1, j)) / two_dx;
                let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / two_dy;
                sum += dv_dx - du_dy;
            }
        }
        sum / ((nx - 2) * (ny - 2)) as f64
    }

    /// The **root-mean-square interior vorticity** `ω_rms = √(⟨ω²⟩)` (1/s) over the
    /// interior cells — the magnitude scale of the central-difference vorticity, completing
    /// the peak / mean / rms set ([`FlowSolution::max_vorticity`],
    /// [`FlowSolution::mean_vorticity`]). It relates to the [`FlowSolution::enstrophy`] by
    /// `ω_rms = √(2E/A)` (enstrophy is `½∫ω² dA`), and it is always `≥ |⟨ω⟩|` (Cauchy–
    /// Schwarz), with equality only for a spatially uniform vorticity. Uses the same
    /// interior stencil as [`FlowSolution::mean_vorticity`]; returns `0` for a grid too
    /// small for an interior difference (`nx < 3` or `ny < 3`).
    pub fn rms_vorticity(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx < 3 || ny < 3 {
            return 0.0;
        }
        let (two_dx, two_dy) = (2.0 * self.grid.dx(), 2.0 * self.grid.dy());
        let mut sum_sq = 0.0;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let dv_dx = (self.v_at_cell(i + 1, j) - self.v_at_cell(i - 1, j)) / two_dx;
                let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / two_dy;
                let omega = dv_dx - du_dy;
                sum_sq += omega * omega;
            }
        }
        (sum_sq / ((nx - 2) * (ny - 2)) as f64).sqrt()
    }

    /// The **standard deviation of the interior vorticity** `σ_ω = √(⟨(ω − ⟨ω⟩)²⟩)` (1/s)
    /// over the interior cells — the spatial fluctuation of the central-difference vorticity
    /// about its area-average, completing the mean / rms / std triad
    /// ([`FlowSolution::mean_vorticity`], [`FlowSolution::rms_vorticity`]). By the variance
    /// identity it satisfies `σ_ω² = ω_rms² − ⟨ω⟩²`; it is `0` exactly for a spatially
    /// uniform vorticity and grows as the swirl becomes more non-uniform. Uses the same
    /// interior stencil as [`FlowSolution::mean_vorticity`] (deviation-sum form); returns `0`
    /// for a grid too small for an interior difference (`nx < 3` or `ny < 3`).
    pub fn vorticity_std_dev(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx < 3 || ny < 3 {
            return 0.0;
        }
        let mean = self.mean_vorticity();
        let (two_dx, two_dy) = (2.0 * self.grid.dx(), 2.0 * self.grid.dy());
        let mut sum_sq = 0.0;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let dv_dx = (self.v_at_cell(i + 1, j) - self.v_at_cell(i - 1, j)) / two_dx;
                let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / two_dy;
                let d = (dv_dx - du_dy) - mean;
                sum_sq += d * d;
            }
        }
        (sum_sq / ((nx - 2) * (ny - 2)) as f64).sqrt()
    }

    /// The total **enstrophy** `E = ½∫ω² dA` (m²/s²) over the interior cells — the
    /// integrated *squared* vorticity, a strictly non-negative measure of the
    /// flow's total rotational intensity. Where the circulation `Γ = ∫ω dA` is a
    /// signed sum in which equal-and-opposite swirl cancels, and the peak vorticity
    /// is a single point, enstrophy accumulates rotation everywhere it occurs — the
    /// quantity whose downscale cascade governs 2-D turbulence and that scales the
    /// viscous dissipation in incompressible flow. Uses the same interior
    /// central-difference vorticity stencil as [`FlowSolution::max_vorticity`],
    /// summed with `½ω²` weighted by the cell area. Returns `0` for a grid too
    /// small for an interior difference (`nx < 3` or `ny < 3`).
    pub fn enstrophy(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx < 3 || ny < 3 {
            return 0.0;
        }
        let (dx, dy) = (self.grid.dx(), self.grid.dy());
        let (two_dx, two_dy) = (2.0 * dx, 2.0 * dy);
        let mut sum = 0.0;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let dv_dx = (self.v_at_cell(i + 1, j) - self.v_at_cell(i - 1, j)) / two_dx;
                let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / two_dy;
                let omega = dv_dx - du_dy;
                sum += omega * omega;
            }
        }
        0.5 * sum * dx * dy
    }

    /// The peak **Q-criterion** `Q = ½(‖Ω‖² − ‖S‖²)` (1/s²) over the interior
    /// cells — the standard scalar that *identifies* a vortex as a region where
    /// rotation (the spin tensor `Ω`) outweighs strain (the rate-of-strain tensor
    /// `S`). Where the vorticity is also large in a straight shear layer, `Q > 0`
    /// only where the fluid genuinely swirls: a **pure shear** flow has equal
    /// rotation and strain and gives `Q = 0`, while **solid-body rotation** at
    /// rate `Ω` gives `Q = Ω²`. This is the workhorse vortex-identification
    /// criterion (Hunt, Wray & Moin 1988). Uses the same interior
    /// central-difference stencil as [`FlowSolution::max_vorticity`], now also
    /// forming the normal strains `∂u/∂x` and `∂v/∂y`. Returns `0` for a grid too
    /// small for an interior difference (`nx < 3` or `ny < 3`), or when no cell is
    /// rotation-dominated (`Q ≤ 0` everywhere — no vortex).
    pub fn max_q_criterion(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx < 3 || ny < 3 {
            return 0.0;
        }
        let two_dx = 2.0 * self.grid.dx();
        let two_dy = 2.0 * self.grid.dy();
        let mut peak = 0.0_f64;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let du_dx = (self.u_at_cell(i + 1, j) - self.u_at_cell(i - 1, j)) / two_dx;
                let dv_dy = (self.v_at_cell(i, j + 1) - self.v_at_cell(i, j - 1)) / two_dy;
                let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / two_dy;
                let dv_dx = (self.v_at_cell(i + 1, j) - self.v_at_cell(i - 1, j)) / two_dx;
                // ‖Ω‖² = ½ω² with ω = ∂v/∂x − ∂u/∂y; ‖S‖² from the symmetric part.
                let omega = dv_dx - du_dy;
                let s_off = 0.5 * (du_dy + dv_dx);
                let omega_norm_sq = 0.5 * omega * omega;
                let s_norm_sq = du_dx * du_dx + dv_dy * dv_dy + 2.0 * s_off * s_off;
                let q = 0.5 * (omega_norm_sq - s_norm_sq);
                peak = peak.max(q);
            }
        }
        peak
    }

    /// The total **viscous dissipation rate** `Φ = 2μ·∫(S:S) dA` (W·m⁻¹, i.e. per
    /// unit depth in 2-D) over the interior cells — the rate at which the flow
    /// *irreversibly* converts mechanical energy into heat, where `S` is the
    /// rate-of-strain tensor and `S:S = S_ij S_ij`. `dynamic_viscosity` is the
    /// fluid's dynamic viscosity `μ = ρ·ν` (Pa·s). Unlike the rotation diagnostics
    /// it is built from the *full* strain rate, not just the vorticity, so a pure
    /// straining flow — which has zero vorticity and zero enstrophy — still
    /// dissipates: `Φ` is strictly non-negative and vanishes only for rigid
    /// motion. Uses the same interior central-difference stencil as
    /// [`FlowSolution::max_q_criterion`], with `S:S = u_x² + v_y² + ½(u_y+v_x)²`.
    /// Returns `0` for a grid too small for an interior difference (`nx < 3` or
    /// `ny < 3`).
    pub fn viscous_dissipation(&self, dynamic_viscosity: f64) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx < 3 || ny < 3 {
            return 0.0;
        }
        let (dx, dy) = (self.grid.dx(), self.grid.dy());
        let (two_dx, two_dy) = (2.0 * dx, 2.0 * dy);
        let mut sum = 0.0;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let du_dx = (self.u_at_cell(i + 1, j) - self.u_at_cell(i - 1, j)) / two_dx;
                let dv_dy = (self.v_at_cell(i, j + 1) - self.v_at_cell(i, j - 1)) / two_dy;
                let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / two_dy;
                let dv_dx = (self.v_at_cell(i + 1, j) - self.v_at_cell(i - 1, j)) / two_dx;
                let s_off = 0.5 * (du_dy + dv_dx);
                // S:S = u_x² + v_y² + 2·(½(u_y+v_x))².
                sum += du_dx * du_dx + dv_dy * dv_dy + 2.0 * s_off * s_off;
            }
        }
        2.0 * dynamic_viscosity * sum * dx * dy
    }

    /// The total **palinstrophy** `P = ½∫|∇ω|² dA` (s⁻²) over the doubly-interior
    /// cells — the integrated squared *gradient* of vorticity. It completes the
    /// 2-D-turbulence cascade trio (kinetic energy → enstrophy → palinstrophy):
    /// where the enstrophy `½∫ω²` measures how much rotation there is, the
    /// palinstrophy measures how *finely structured* it is, and it sets the rate
    /// at which viscosity destroys enstrophy. Because it sees the gradient, a flow
    /// with strong but *uniform* vorticity (a plain shear layer) has large
    /// enstrophy yet zero palinstrophy. Forms `ω = ∂v/∂x − ∂u/∂y` at each interior
    /// cell with the same central-difference stencil as
    /// [`FlowSolution::max_vorticity`], then central-differences `ω` itself for
    /// `∇ω`. Returns `0` for a grid too small for a doubly-interior difference
    /// (`nx < 5` or `ny < 5`).
    pub fn palinstrophy(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx < 5 || ny < 5 {
            return 0.0;
        }
        let (dx, dy) = (self.grid.dx(), self.grid.dy());
        let (two_dx, two_dy) = (2.0 * dx, 2.0 * dy);
        // Vorticity at an interior cell (1 ≤ i ≤ nx−2, 1 ≤ j ≤ ny−2).
        let vort = |i: usize, j: usize| -> f64 {
            let dv_dx = (self.v_at_cell(i + 1, j) - self.v_at_cell(i - 1, j)) / two_dx;
            let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / two_dy;
            dv_dx - du_dy
        };
        let mut sum = 0.0;
        for j in 2..ny - 2 {
            for i in 2..nx - 2 {
                let dw_dx = (vort(i + 1, j) - vort(i - 1, j)) / two_dx;
                let dw_dy = (vort(i, j + 1) - vort(i, j - 1)) / two_dy;
                sum += dw_dx * dw_dx + dw_dy * dw_dy;
            }
        }
        0.5 * sum * dx * dy
    }

    /// The peak **strain rate** `√(2·S:S)` (1/s) over the interior cells — the
    /// strongest local *rate of deformation*, the symmetric-tensor companion to
    /// the [`FlowSolution::max_vorticity`] (rotation). The velocity gradient
    /// splits into a symmetric rate-of-strain `S` and an antisymmetric spin: the
    /// vorticity sees only the spin, this sees only `S = √(2·S_ij S_ij)`, so a
    /// solid-body rotation (all spin, no deformation) reads zero here while a
    /// pure strain (all deformation, no spin) reads zero on the vorticity. It is
    /// the rate that stretches fluid elements and sets mixing, non-Newtonian
    /// viscosity and cavitation onset. Uses the same interior central-difference
    /// stencil as [`FlowSolution::max_q_criterion`], with
    /// `S:S = u_x² + v_y² + ½(u_y+v_x)²`. Returns `0` for a grid too small for an
    /// interior difference (`nx < 3` or `ny < 3`).
    pub fn max_strain_rate(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx < 3 || ny < 3 {
            return 0.0;
        }
        let two_dx = 2.0 * self.grid.dx();
        let two_dy = 2.0 * self.grid.dy();
        let mut peak = 0.0_f64;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let du_dx = (self.u_at_cell(i + 1, j) - self.u_at_cell(i - 1, j)) / two_dx;
                let dv_dy = (self.v_at_cell(i, j + 1) - self.v_at_cell(i, j - 1)) / two_dy;
                let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / two_dy;
                let dv_dx = (self.v_at_cell(i + 1, j) - self.v_at_cell(i - 1, j)) / two_dx;
                let s_off = 0.5 * (du_dy + dv_dx);
                // S:S = u_x² + v_y² + 2·(½(u_y+v_x))²; the strain rate is √(2·S:S).
                let s_dd = du_dx * du_dx + dv_dy * dv_dy + 2.0 * s_off * s_off;
                peak = peak.max((2.0 * s_dd).sqrt());
            }
        }
        peak
    }

    /// The **area-averaged rate-of-strain magnitude** `⟨|S|⟩` (1/s) over the interior cells —
    /// the deformation analog of [`FlowSolution::mean_vorticity`]: where the vorticity mean
    /// captures the antisymmetric (rotation) part of `∇u`, this captures the symmetric
    /// (deformation/shear) part, `|S| = √(2 SᵢⱼSᵢⱼ)`. It is the companion to the peak
    /// [`FlowSolution::max_strain_rate`], summing the per-cell
    /// [`FlowSolution::strain_rate_at_cell`]. Returns `0` for a grid too small for an interior
    /// central difference (`nx < 3` or `ny < 3`).
    pub fn mean_strain_rate(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx < 3 || ny < 3 {
            return 0.0;
        }
        let mut sum = 0.0;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                sum += self.strain_rate_at_cell(i, j);
            }
        }
        sum / ((nx - 2) * (ny - 2)) as f64
    }

    /// The **root-mean-square rate-of-strain magnitude** `√(⟨|S|²⟩)` (1/s) over the interior
    /// cells — the L2 deformation rate. Since `|S|² = 2 SᵢⱼSᵢⱼ`, this is the quantity that
    /// sets the viscous dissipation `ε = 2ν⟨SᵢⱼSᵢⱼ⟩ = ν⟨|S|²⟩`, so it measures how strongly
    /// the flow is being strained (and thus dissipating energy). It completes the strain-rate
    /// family with [`FlowSolution::mean_strain_rate`] and the peak
    /// [`FlowSolution::max_strain_rate`], mirroring [`FlowSolution::rms_vorticity`]; `0` for a
    /// grid too small for an interior central difference (`nx < 3` or `ny < 3`).
    pub fn rms_strain_rate(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx < 3 || ny < 3 {
            return 0.0;
        }
        let mut sum_sq = 0.0;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let sr = self.strain_rate_at_cell(i, j);
                sum_sq += sr * sr;
            }
        }
        (sum_sq / ((nx - 2) * (ny - 2)) as f64).sqrt()
    }

    /// The **mean viscous dissipation rate** `ε = ν·⟨|S|²⟩` (m²/s³, specific — per unit mass)
    /// for a fluid of kinematic viscosity `kinematic_viscosity` `ν` (m²/s) — the rate at which
    /// the flow irreversibly converts kinetic energy to heat by viscosity, the quantity that
    /// closes the turbulent energy budget (the Kolmogorov cascade rate). Since
    /// `|S|² = 2 SᵢⱼSᵢⱼ`, `ε = 2ν⟨SᵢⱼSᵢⱼ⟩ = ν·⟨|S|²⟩ = ν·(rms_strain_rate)²`; it is the
    /// strain-squared companion to the [`FlowSolution::rms_strain_rate`]. A pure rotation
    /// strains nothing and dissipates nothing; only deformation dissipates. Returns `0` for
    /// zero viscosity or a grid too small for an interior strain-rate stencil.
    pub fn mean_dissipation_rate(&self, kinematic_viscosity: f64) -> f64 {
        kinematic_viscosity * self.rms_strain_rate().powi(2)
    }

    /// The **Kolmogorov length scale** `η = (ν³/ε)^(1/4)` (m) — the size of the smallest
    /// turbulent eddies, where viscous dissipation balances the energy cascade and gradients
    /// are smoothed out. It is set by the kinematic viscosity `kinematic_viscosity` `ν`
    /// (m²/s) and the mean dissipation rate `ε` ([`FlowSolution::mean_dissipation_rate`]):
    /// stronger dissipation (or lower viscosity) drives `η` smaller. Returns `0` when there is
    /// no dissipation (`ε ≤ 0` — a strainless or quiescent flow), where the scale is
    /// undefined.
    pub fn kolmogorov_length_scale(&self, kinematic_viscosity: f64) -> f64 {
        let eps = self.mean_dissipation_rate(kinematic_viscosity);
        if eps <= 0.0 {
            return 0.0;
        }
        (kinematic_viscosity.powi(3) / eps).powf(0.25)
    }

    /// The **Kolmogorov time scale** `τ_η = (ν/ε)^(1/2)` (s) — the turnover time (lifetime) of
    /// the smallest dissipative eddies, set by the kinematic viscosity `kinematic_viscosity`
    /// `ν` (m²/s) and the mean dissipation rate `ε` ([`FlowSolution::mean_dissipation_rate`]).
    /// It is the temporal companion to the [`FlowSolution::kolmogorov_length_scale`] `η`: the
    /// two satisfy the unit Kolmogorov-scale Reynolds number `η²/(τ_η·ν) = 1` by construction.
    /// Returns `0` when there is no dissipation (`ε ≤ 0`), where the scale is undefined.
    pub fn kolmogorov_time_scale(&self, kinematic_viscosity: f64) -> f64 {
        let eps = self.mean_dissipation_rate(kinematic_viscosity);
        if eps <= 0.0 {
            return 0.0;
        }
        (kinematic_viscosity / eps).sqrt()
    }

    /// The **Kolmogorov velocity scale** `u_η = (ν·ε)^(1/4)` (m/s) — the characteristic
    /// velocity of the smallest dissipative eddies, set by the kinematic viscosity
    /// `kinematic_viscosity` `ν` (m²/s) and the mean dissipation rate `ε`
    /// ([`FlowSolution::mean_dissipation_rate`]). It completes the Kolmogorov microscale triad
    /// with the [`FlowSolution::kolmogorov_length_scale`] `η` and
    /// [`FlowSolution::kolmogorov_time_scale`] `τ_η`: the three are self-consistent
    /// (`u_η = η/τ_η`) and define the unit Kolmogorov-scale Reynolds number `u_η·η/ν = 1`.
    /// Returns `0` when there is no dissipation (`ε ≤ 0`).
    pub fn kolmogorov_velocity_scale(&self, kinematic_viscosity: f64) -> f64 {
        let eps = self.mean_dissipation_rate(kinematic_viscosity);
        if eps <= 0.0 {
            return 0.0;
        }
        (kinematic_viscosity * eps).powf(0.25)
    }

    /// The **Taylor microscale** `λ = u'·√(15ν/ε)` (m) — the intermediate turbulence
    /// length lying between the energy-containing integral scale and the dissipative
    /// [`FlowSolution::kolmogorov_length_scale`] `η`. It is the length over which the
    /// velocity field stays correlated to the strain that dissipates it, defined (for
    /// isotropic turbulence) by `ε = 15·ν·u'²/λ²`. Here `u'` is the rms velocity
    /// ([`FlowSolution::rms_speed`]), `ε` the mean dissipation rate
    /// ([`FlowSolution::mean_dissipation_rate`]), and `kinematic_viscosity` `ν` (m²/s).
    /// This is the standard *isotropic estimate* — a research-grade diagnostic, not a
    /// resolved two-point correlation; unlike the Kolmogorov scales it carries the
    /// large-eddy velocity `u'`, so it is not a purely viscous scale. Returns `0` when
    /// there is no dissipation (`ε ≤ 0` — a strainless or quiescent flow), where the
    /// scale is undefined.
    pub fn taylor_microscale(&self, kinematic_viscosity: f64) -> f64 {
        let eps = self.mean_dissipation_rate(kinematic_viscosity);
        if eps <= 0.0 {
            return 0.0;
        }
        self.rms_speed() * (15.0 * kinematic_viscosity / eps).sqrt()
    }

    /// The **Taylor-microscale Reynolds number** `Re_λ = u'·λ / ν` (dimensionless) — the
    /// Reynolds number formed on the rms velocity `u'` ([`FlowSolution::rms_speed`]) and
    /// the [`FlowSolution::taylor_microscale`] `λ`, with kinematic viscosity
    /// `kinematic_viscosity` `ν` (m²/s). It is the standard *intensity* measure of
    /// turbulence — a large `Re_λ` means a wide separation between the energy-containing
    /// and dissipative scales (vigorous, multi-scale turbulence), a small `Re_λ` a flow
    /// the viscosity keeps close to laminar. Distinct from the geometry-based
    /// [`FlowSolution::domain_reynolds_number`] (built on the domain length and peak
    /// speed): this one is built on the *internal* turbulence scales. Like the Taylor
    /// microscale it is the isotropic estimate (research-grade). Returns `0` for
    /// non-positive or non-finite viscosity, or when there is no dissipation (a strainless
    /// or quiescent flow, where `λ = 0`).
    pub fn taylor_reynolds_number(&self, kinematic_viscosity: f64) -> f64 {
        if !kinematic_viscosity.is_finite() || kinematic_viscosity <= 0.0 {
            return 0.0;
        }
        self.rms_speed() * self.taylor_microscale(kinematic_viscosity) / kinematic_viscosity
    }

    /// The **integral (energy-containing) length scale** `L = u'³ / ε` (m) — the standard
    /// scaling estimate of the *largest* turbulence scale, the size of the energy-bearing
    /// eddies that hold most of the kinetic energy. It tops the length-scale hierarchy
    /// above the [`FlowSolution::taylor_microscale`] `λ` and the dissipative
    /// [`FlowSolution::kolmogorov_length_scale`] `η`. Here `u'` is the rms velocity
    /// ([`FlowSolution::rms_speed`]) and `ε` the mean dissipation rate
    /// ([`FlowSolution::mean_dissipation_rate`]), via the Taylor energy-cascade scaling
    /// `ε ~ u'³/L`. This is the *standard scaling estimate* with the order-one cascade
    /// constant absorbed (research-grade) — not a resolved two-point-correlation integral.
    /// Returns `0` when there is no dissipation (`ε ≤ 0` — a strainless or quiescent flow),
    /// where the scale is undefined.
    pub fn integral_length_scale(&self, kinematic_viscosity: f64) -> f64 {
        let eps = self.mean_dissipation_rate(kinematic_viscosity);
        if eps <= 0.0 {
            return 0.0;
        }
        self.rms_speed().powi(3) / eps
    }

    /// The **integral (eddy-turnover) time scale** `T = u'² / ε` (s) — the standard
    /// scaling estimate of the *largest* turbulence time, the turnover time of the
    /// energy-containing eddies (equivalently `T = L/u'`, the
    /// [`FlowSolution::integral_length_scale`] over the rms velocity). It is the rate at
    /// which the large eddies hand their energy to the cascade — the `k/ε` time of RANS
    /// closure (up to the ½ in `k`) — and it caps the time hierarchy above the dissipative
    /// [`FlowSolution::kolmogorov_time_scale`] `τ_η`, parallel to `L` above `η`. Here `u'`
    /// is the rms velocity ([`FlowSolution::rms_speed`]) and `ε` the mean dissipation rate
    /// ([`FlowSolution::mean_dissipation_rate`]). This is the *standard scaling estimate*
    /// with the order-one constant absorbed (research-grade). Returns `0` when there is no
    /// dissipation (`ε ≤ 0` — a strainless or quiescent flow), where the scale is undefined.
    pub fn integral_time_scale(&self, kinematic_viscosity: f64) -> f64 {
        let eps = self.mean_dissipation_rate(kinematic_viscosity);
        if eps <= 0.0 {
            return 0.0;
        }
        self.rms_speed().powi(2) / eps
    }

    /// The **Taylor-to-Kolmogorov scale ratio** `λ/η` (dimensionless) — the ratio of the
    /// intermediate [`FlowSolution::taylor_microscale`] `λ` to the dissipative
    /// [`FlowSolution::kolmogorov_length_scale`] `η`, a measure of the *width of the
    /// inertial subrange*: the larger the separation between the energy-carrying and
    /// dissipative scales, the broader the range of eddies over which the turbulent cascade
    /// operates. It grows with Reynolds number (for a fixed flow, `∝ 1/√ν`), and in
    /// isotropic turbulence satisfies `λ/η = 15^(1/4)·√(Re_λ)` with the Taylor-microscale
    /// Reynolds number [`FlowSolution::taylor_reynolds_number`]. Returns `0` when there is
    /// no dissipation (a strainless or quiescent flow, where both scales vanish).
    pub fn taylor_to_kolmogorov_ratio(&self, kinematic_viscosity: f64) -> f64 {
        let eta = self.kolmogorov_length_scale(kinematic_viscosity);
        if eta <= 0.0 {
            return 0.0;
        }
        self.taylor_microscale(kinematic_viscosity) / eta
    }

    /// The **integral-to-Kolmogorov scale ratio** `L/η` (dimensionless) — the ratio of the
    /// energy-containing [`FlowSolution::integral_length_scale`] `L` to the dissipative
    /// [`FlowSolution::kolmogorov_length_scale`] `η`, the *total* span of the turbulent
    /// cascade from the largest eddies to the smallest. It is the companion to the
    /// [`FlowSolution::taylor_to_kolmogorov_ratio`] `λ/η` and the widest of the
    /// scale-separation measures: it grows strongly with Reynolds number (for a fixed flow,
    /// `∝ 1/ν^(3/2)`), and in isotropic turbulence satisfies `L/η = Re_λ^(3/2) / 15^(3/4)`
    /// with the Taylor-microscale Reynolds number [`FlowSolution::taylor_reynolds_number`].
    /// Returns `0` when there is no dissipation (a strainless or quiescent flow, where both
    /// scales vanish).
    pub fn integral_to_kolmogorov_ratio(&self, kinematic_viscosity: f64) -> f64 {
        let eta = self.kolmogorov_length_scale(kinematic_viscosity);
        if eta <= 0.0 {
            return 0.0;
        }
        self.integral_length_scale(kinematic_viscosity) / eta
    }

    /// The **integral-to-Taylor scale ratio** `L/λ` (dimensionless) — the ratio of the
    /// energy-containing [`FlowSolution::integral_length_scale`] `L` to the intermediate
    /// [`FlowSolution::taylor_microscale`] `λ`. It completes the length-scale-ratio trio
    /// with `λ/η` ([`FlowSolution::taylor_to_kolmogorov_ratio`]) and `L/η`
    /// ([`FlowSolution::integral_to_kolmogorov_ratio`]). It grows *linearly* with Reynolds
    /// number (for a fixed flow, `∝ 1/ν`), and in isotropic turbulence is exactly
    /// `L/λ = Re_λ / 15` with the Taylor-microscale Reynolds number
    /// [`FlowSolution::taylor_reynolds_number`]. Returns `0` when there is no dissipation
    /// (a strainless or quiescent flow, where both scales vanish).
    pub fn integral_to_taylor_ratio(&self, kinematic_viscosity: f64) -> f64 {
        let lambda = self.taylor_microscale(kinematic_viscosity);
        if lambda <= 0.0 {
            return 0.0;
        }
        self.integral_length_scale(kinematic_viscosity) / lambda
    }

    /// The **integral-to-Kolmogorov time-scale ratio** `T/τ_η` (dimensionless) — the ratio
    /// of the large-eddy turnover time [`FlowSolution::integral_time_scale`] `T` to the
    /// dissipative [`FlowSolution::kolmogorov_time_scale`] `τ_η`, the time-domain analog of
    /// the length ratio `L/η` ([`FlowSolution::integral_to_kolmogorov_ratio`]). It measures
    /// how many Kolmogorov turnovers fit inside one large-eddy turnover; it grows with
    /// Reynolds number (for a fixed flow, `∝ 1/ν`), and in isotropic turbulence is exactly
    /// `T/τ_η = Re_λ / √15` with the Taylor-microscale Reynolds number
    /// [`FlowSolution::taylor_reynolds_number`]. Returns `0` when there is no dissipation (a
    /// strainless or quiescent flow, where both time scales vanish).
    pub fn integral_to_kolmogorov_time_ratio(&self, kinematic_viscosity: f64) -> f64 {
        let tau_eta = self.kolmogorov_time_scale(kinematic_viscosity);
        if tau_eta <= 0.0 {
            return 0.0;
        }
        self.integral_time_scale(kinematic_viscosity) / tau_eta
    }

    /// The peak **local mass-continuity residual** `max |∇·u|` (1/s) — the
    /// largest pointwise velocity divergence `∂u/∂x + ∂v/∂y` over the cells,
    /// formed straight from the MAC face velocities the way the projection step
    /// measures it: `(u_{i+1,j} − u_{i,j})/dx + (v_{i,j+1} − v_{i,j})/dy`. For an
    /// incompressible solve this is exactly the field the pressure-Poisson
    /// equation drives to zero, so its maximum is the honest *local* convergence
    /// check — a well-converged flow reads ≈ 0 in every cell, and a large value
    /// pinpoints a cell where continuity is still violated. This is distinct
    /// from [`FlowSolution::continuity_error`], which is the *global*
    /// inlet-to-outlet mass-balance ratio rather than the worst single cell.
    /// Returns `0` for an empty grid. Note the discrete divergence uses the
    /// cell's own bracketing faces directly, so it is defined on every cell (no
    /// interior-stencil restriction).
    pub fn max_divergence(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        let (dx, dy) = (self.grid.dx(), self.grid.dy());
        let mut worst = 0.0_f64;
        for j in 0..ny {
            for i in 0..nx {
                let du_dx = (self.u.at(i + 1, j) - self.u.at(i, j)) / dx;
                let dv_dy = (self.v.at(i, j + 1) - self.v.at(i, j)) / dy;
                worst = worst.max((du_dx + dv_dy).abs());
            }
        }
        worst
    }

    /// The **area-averaged divergence** `⟨∇·u⟩ = (1/N)·Σ ∇·u` (1/s) over all `nx·ny` cells —
    /// the mean continuity residual. For an incompressible flow `∇·u = 0` pointwise, so this
    /// is `≈ 0` for a well-converged solution and grows as the discrete continuity constraint
    /// is violated: it is the signed-average companion to the peak
    /// [`FlowSolution::max_divergence`], summing the per-cell
    /// [`FlowSolution::divergence_at_cell`]. Returns `0` for an empty grid.
    pub fn mean_divergence(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        let n = nx * ny;
        if n == 0 {
            return 0.0;
        }
        let mut sum = 0.0;
        for j in 0..ny {
            for i in 0..nx {
                sum += self.divergence_at_cell(i, j);
            }
        }
        sum / n as f64
    }

    /// The **root-mean-square divergence** `√(⟨(∇·u)²⟩)` (1/s) over all `nx·ny` cells — the
    /// L2 continuity residual. Unlike the signed [`FlowSolution::mean_divergence`] (whose
    /// positive and negative cells cancel), the RMS sums squares, so it never cancels and is
    /// the standard convergence metric for how well the discrete continuity constraint
    /// `∇·u = 0` is satisfied. The square-summing companion to the signed mean and the peak
    /// [`FlowSolution::max_divergence`]; `0` for an empty grid.
    pub fn rms_divergence(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        let n = nx * ny;
        if n == 0 {
            return 0.0;
        }
        let mut sum_sq = 0.0;
        for j in 0..ny {
            for i in 0..nx {
                let d = self.divergence_at_cell(i, j);
                sum_sq += d * d;
            }
        }
        (sum_sq / n as f64).sqrt()
    }

    /// The **vorticity (shear-layer) thickness** `δ_ω = (u_max − u_min)/max|∂u/∂y|`
    /// (m) — the canonical mixing-layer length scale: the streamwise-velocity span
    /// across the flow divided by the steepest wall-normal velocity gradient. Where
    /// the other field diagnostics report a magnitude, rate, location, or area, this
    /// is a *length* — the characteristic distance over which the shear is spread,
    /// the scale that sets a free-shear layer's growth and instability. `u_max` and
    /// `u_min` are taken over the cell-centred streamwise velocity, and `max|∂u/∂y|`
    /// from a central difference over the interior rows. Returns `0` for a grid too
    /// small for an interior row (`ny < 3`), an empty grid, or a uniform
    /// (shear-free) flow, where no shear layer is defined.
    pub fn vorticity_thickness(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny < 3 {
            return 0.0;
        }
        let dy = self.grid.dy();
        let mut u_min = f64::INFINITY;
        let mut u_max = f64::NEG_INFINITY;
        for j in 0..ny {
            for i in 0..nx {
                let u = self.u_at_cell(i, j);
                u_min = u_min.min(u);
                u_max = u_max.max(u);
            }
        }
        let mut max_du_dy = 0.0_f64;
        for j in 1..ny - 1 {
            for i in 0..nx {
                let du_dy = (self.u_at_cell(i, j + 1) - self.u_at_cell(i, j - 1)) / (2.0 * dy);
                max_du_dy = max_du_dy.max(du_dy.abs());
            }
        }
        if max_du_dy <= 0.0 {
            return 0.0;
        }
        (u_max - u_min) / max_du_dy
    }

    /// The **reverse-flow (recirculation) area fraction** — the share of cells
    /// whose streamwise (`x`) velocity `u` is negative, i.e. running *back*
    /// against the main throughflow. This is the quantitative separation /
    /// recirculation diagnostic: a fully attached forward flow reads `0`, while
    /// a recirculation zone — the lid-driven cavity's returning lower vortex, or
    /// a separation bubble behind a step or bluff body — reads the share of the
    /// domain flowing backward. Where every other field measure here reads a
    /// *magnitude* (speed, vorticity, strain, pressure), this reads the *sign*
    /// of the streamwise velocity, the flow-reversal topology those miss.
    /// Returns `0` for an empty grid.
    pub fn reverse_flow_fraction(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        let n = nx * ny;
        if n == 0 {
            return 0.0;
        }
        let mut reversed = 0usize;
        for j in 0..ny {
            for i in 0..nx {
                if self.u_at_cell(i, j) < 0.0 {
                    reversed += 1;
                }
            }
        }
        reversed as f64 / n as f64
    }

    /// The **stream-function span** `max ψ − min ψ` (m²/s) over the domain — the
    /// total volumetric throughflow carried between the extreme streamlines (for
    /// an enclosed cavity, the strength of its recirculating vortex). The stream
    /// function `ψ`, defined by `u = ∂ψ/∂y` and `v = −∂ψ/∂x`, is the foundational
    /// 2-D incompressible field whose iso-contours *are* the streamlines — the
    /// single scalar that encodes the whole flow's path structure, which none of
    /// the velocity / vorticity / pressure diagnostics expose directly.
    ///
    /// `ψ` is reconstructed at the `(nx+1)×(ny+1)` grid corners by integrating
    /// the staggered face velocities from a `ψ = 0` corner: up a column
    /// `ψ_{i,j+1} = ψ_{i,j} + u_{i,j}·dy`, along the base row
    /// `ψ_{i+1,0} = ψ_{i,0} − v_{i,0}·dx`. For a divergence-free field this is
    /// path-independent; the converged incompressible solution is divergence-free
    /// to the solver tolerance, so the span is well defined. Returns `0` for an
    /// empty grid.
    pub fn stream_function_range(&self) -> f64 {
        let (nx, ny) = (self.grid.nx, self.grid.ny);
        if nx == 0 || ny == 0 {
            return 0.0;
        }
        let (dx, dy) = (self.grid.dx(), self.grid.dy());
        // ψ at the (nx+1)×(ny+1) corners, indexed [column i][row j], ψ[0][0] = 0.
        let mut psi = vec![vec![0.0_f64; ny + 1]; nx + 1];
        // Base row (j = 0): integrate the v faces, ψ_{i+1,0} = ψ_{i,0} − v·dx.
        for i in 0..nx {
            psi[i + 1][0] = psi[i][0] - self.v.at(i, 0) * dx;
        }
        // Each column: integrate the u faces upward, ψ_{i,j+1} = ψ_{i,j} + u·dy.
        for (i, col) in psi.iter_mut().enumerate() {
            for j in 0..ny {
                col[j + 1] = col[j] + self.u.at(i, j) * dy;
            }
        }
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for &p in psi.iter().flatten() {
            lo = lo.min(p);
            hi = hi.max(p);
        }
        hi - lo
    }
}

/// Solve a steady laminar flow with the SIMPLE algorithm.
///
/// `grid` is the staggered finite-volume grid; `fluid` supplies the
/// constant properties; `bcs` sets the four-sided boundary; `controls`
/// tunes the under-relaxation and the convergence criteria.
///
/// Returns the [`FlowSolution`] — the converged velocity and pressure
/// fields plus the iteration / residual diagnostics. The solver always
/// returns a field, even if it hits the iteration cap before the
/// tolerance ([`FlowSolution::converged`] reports which).
pub fn solve_simple(
    grid: &Grid,
    fluid: &Fluid,
    bcs: &Boundaries,
    controls: &SimpleControls,
) -> FlowSolution {
    let (sol, _turb) = solve_simple_with(
        grid,
        fluid,
        bcs,
        controls,
        &EffectiveViscosity::Laminar,
    );
    sol
}

/// The optional turbulence-field bundle [`solve_simple_with`] returns
/// alongside the flow solution — `None` for a laminar run, `Some(...)`
/// when a two-equation model produced the eddy-viscosity field used in
/// the momentum solve.
#[derive(Clone, Debug)]
pub enum TurbulenceSnapshot {
    /// The converged k-ε `k`, `ε`, `ν_t` fields.
    KEpsilon(TurbulenceField),
    /// The converged k-ω SST `k`, `ω`, `ν_t` fields.
    Sst(SstField),
}

/// Solve a steady incompressible flow with the SIMPLE algorithm and
/// the requested **effective-viscosity** model.
///
/// `EffectiveViscosity::Laminar` is the historical laminar path
/// (equivalent to [`solve_simple`]). `EffectiveViscosity::KEpsilon`
/// and `EffectiveViscosity::Sst` couple a two-equation RANS closure
/// into the momentum solve: at every outer iteration the turbulence
/// transport equations are advanced, the cell-centred eddy viscosity
/// `ν_t` is refreshed, and the momentum-equation face diffusion uses
/// `ν_eff = ν + ν_t` (averaged from the two bracketing cells onto the
/// staggered face). The pressure-correction step is unchanged.
///
/// Returns the [`FlowSolution`] and, for a turbulent run, the
/// `TurbulenceSnapshot` carrying the converged turbulence fields.
pub fn solve_simple_with(
    grid: &Grid,
    fluid: &Fluid,
    bcs: &Boundaries,
    controls: &SimpleControls,
    visc: &EffectiveViscosity,
) -> (FlowSolution, Option<TurbulenceSnapshot>) {
    let nx = grid.nx;
    let ny = grid.ny;
    // `dy` is needed for the mass-flux residual scale below; `dx` is
    // only used inside the per-equation kernels, which read it from
    // `grid` themselves.
    let dy = grid.dy();
    let rho = fluid.density;
    let nu = fluid.viscosity;

    // Field storage.
    let mut u = grid.u_field(); // (nx+1) × ny
    let mut v = grid.v_field(); // nx × (ny+1)
    let mut p = grid.pressure_field(); // nx × ny
    // Apply the boundary velocities to the field once up front.
    apply_velocity_bcs(&mut u, &mut v, bcs, nx, ny);

    // Per-face momentum-equation diagonal coefficients — reused to
    // build the pressure-correction stencil.
    let mut apu = Field::zeros(nx + 1, ny);
    let mut apv = Field::zeros(nx, ny + 1);

    // Turbulence state — initialised lazily for the turbulent paths.
    let mut turb_ke: Option<TurbulenceField> = None;
    let mut turb_sst: Option<SstField> = None;
    let mut wall_dist: Option<Field> = None;
    let ref_velocity = reference_velocity(bcs).max(1e-6);
    match visc {
        EffectiveViscosity::Laminar => {}
        EffectiveViscosity::KEpsilon { model, .. } => {
            turb_ke =
                Some(TurbulenceField::initialise(grid, ref_velocity, 0.05, model));
        }
        EffectiveViscosity::Sst { model, walls, .. } => {
            turb_sst = Some(SstField::initialise(grid, ref_velocity, 0.05, model));
            wall_dist = Some(wall_distance_field(grid, *walls));
        }
    }

    let mut residual = f64::INFINITY;
    let mut converged = false;
    let mut iterations = 0;

    // A normalising scale for the mass-imbalance residual so the
    // tolerance is geometry-independent (a characteristic mass flow).
    let mass_scale = (rho * ref_velocity * dy).max(1e-30);

    for outer in 0..controls.max_iterations.max(1) {
        iterations = outer + 1;

        // The cell-centred eddy viscosity for this outer iteration —
        // None means laminar (the kernels use ν alone).
        let nu_t_ref: Option<&Field> = match (&turb_ke, &turb_sst) {
            (Some(t), _) => Some(&t.nu_t),
            (_, Some(s)) => Some(&s.nu_t),
            _ => None,
        };

        // --- (1) momentum predictor: solve u*, v* with p frozen ---
        solve_u_momentum(
            &mut u, &v, &p, &mut apu, grid, rho, nu, nu_t_ref, bcs,
            controls.relax_u,
        );
        solve_v_momentum(
            &u, &mut v, &p, &mut apv, grid, rho, nu, nu_t_ref, bcs,
            controls.relax_u,
        );
        // Re-impose the prescribed boundary velocities (the sweep only
        // updates interior faces; this keeps the BC faces exact).
        apply_velocity_bcs(&mut u, &mut v, bcs, nx, ny);

        // --- (2) pressure-correction Poisson equation ---
        let mut coeffs = PoissonCoeffs::zeros(nx, ny);
        let mass_imbalance =
            assemble_pressure_correction(&mut coeffs, &u, &v, &apu, &apv, grid, rho, bcs);

        // --- (3) solve for p' ---
        // A closed domain (no outlet) gives a fully-singular all-Neumann
        // pressure-correction system → pin the gauge by zero-mean. An
        // open domain has an outlet cell anchored to p' = 0 inside the
        // assembly, so the system is already non-singular.
        let pin_mean = !has_outlet(bcs);
        let mut pcorr = Field::zeros(nx, ny);
        let _ = solve_pressure_poisson(
            &controls.pressure_solver,
            &coeffs,
            &mut pcorr,
            controls.sor_omega,
            controls.sor_iterations,
            controls.tolerance * 1e-3,
            pin_mean,
        );

        // --- (4) correct pressure and velocity ---
        correct_pressure(&mut p, &pcorr, controls.relax_p);
        correct_velocity(&mut u, &mut v, &pcorr, &apu, &apv, grid, bcs);
        apply_velocity_bcs(&mut u, &mut v, bcs, nx, ny);

        // --- (5) advance the turbulence transport ---
        match visc {
            EffectiveViscosity::Laminar => {}
            EffectiveViscosity::KEpsilon { model, wall, walls } => {
                if let Some(turb) = turb_ke.as_mut() {
                    advance_k_epsilon(
                        turb, &u, &v, grid, rho, nu, model, wall,
                        walls.south, walls.north, 0.5,
                    );
                }
            }
            EffectiveViscosity::Sst { model, wall, walls } => {
                if let (Some(state), Some(wd)) =
                    (turb_sst.as_mut(), wall_dist.as_ref())
                {
                    advance_k_omega_sst(
                        state, &u, &v, grid, rho, nu, model, wall, wd,
                        walls.south, walls.north, 0.5,
                    );
                }
            }
        }

        // Convergence: the L2 norm of the cell mass imbalance, scaled.
        residual = mass_imbalance / mass_scale;
        if residual.is_finite() && residual <= controls.tolerance {
            converged = true;
            break;
        }
        // Divergence guard — a blown-up field is reported, not looped.
        if !residual.is_finite() || residual > 1e12 {
            break;
        }
    }

    let snapshot = match (turb_ke, turb_sst) {
        (Some(t), _) => Some(TurbulenceSnapshot::KEpsilon(t)),
        (_, Some(s)) => Some(TurbulenceSnapshot::Sst(s)),
        _ => None,
    };

    (
        FlowSolution {
            grid: *grid,
            u,
            v,
            pressure: p,
            iterations,
            residual,
            converged,
        },
        snapshot,
    )
}

/// A characteristic velocity for the residual normalisation — the
/// largest boundary speed in the problem.
fn reference_velocity(bcs: &Boundaries) -> f64 {
    let mag = |s: &SideBc| -> f64 {
        match s {
            SideBc::Wall { u, v } | SideBc::Inlet { u, v } => (u * u + v * v).sqrt(),
            SideBc::Outlet => 0.0,
        }
    };
    mag(&bcs.west)
        .max(mag(&bcs.east))
        .max(mag(&bcs.south))
        .max(mag(&bcs.north))
}

/// Stamp the prescribed boundary velocities onto the velocity fields.
///
/// On the staggered grid the *normal* velocity of a wall / inlet lives
/// exactly on the boundary face, so it is set directly. The *tangential*
/// velocity has no face on the boundary; the no-slip / moving-wall
/// condition is imposed inside the momentum sweep via a wall-shear term,
/// so this routine only fixes the normal-component faces here.
fn apply_velocity_bcs(u: &mut Field, v: &mut Field, bcs: &Boundaries, nx: usize, ny: usize) {
    // West / East edges carry u-faces at i = 0 and i = nx.
    for j in 0..ny {
        match bcs.west {
            SideBc::Wall { u: uw, .. } => u.set(0, j, uw),
            SideBc::Inlet { u: uw, .. } => u.set(0, j, uw),
            SideBc::Outlet => { /* set later from the interior */ }
        }
        match bcs.east {
            SideBc::Wall { u: ue, .. } => u.set(nx, j, ue),
            SideBc::Inlet { u: ue, .. } => u.set(nx, j, ue),
            SideBc::Outlet => {
                // Zero-gradient: the outlet face copies the last
                // interior u-face.
                let interior = u.at(nx - 1, j);
                u.set(nx, j, interior);
            }
        }
    }
    // South / North edges carry v-faces at j = 0 and j = ny.
    for i in 0..nx {
        match bcs.south {
            SideBc::Wall { v: vs, .. } => v.set(i, 0, vs),
            SideBc::Inlet { v: vs, .. } => v.set(i, 0, vs),
            SideBc::Outlet => {}
        }
        match bcs.north {
            SideBc::Wall { v: vn, .. } => v.set(i, ny, vn),
            SideBc::Inlet { v: vn, .. } => v.set(i, ny, vn),
            SideBc::Outlet => {
                let interior = v.at(i, ny - 1);
                v.set(i, ny, interior);
            }
        }
    }
}

/// The hybrid-scheme convection-diffusion coefficient for one face.
///
/// Combines the diffusive conductance `D = ν·A/δ` and the convective
/// mass flow `F` into a single neighbour coefficient by Spalding's
/// hybrid rule `a = max(|F|/2 + D − |F|/2 , D − |F|/2 , 0)` written in
/// the compact `max(D − |F|/2, 0) + max(±F, 0)` upwind form. This is
/// unconditionally stable for any cell Péclet number.
#[inline]
fn hybrid_coeff(d: f64, f: f64, upwind_sign: f64) -> f64 {
    // upwind_sign = +1 for the coefficient that receives inflow from
    // the positive-side neighbour, −1 for the other.
    (d - 0.5 * f.abs()).max(0.0) + (upwind_sign * f).max(0.0)
}

/// Solve the discrete `u`-momentum equation with a few Gauss-Seidel
/// sweeps, writing the provisional `u*` back into `u` and the diagonal
/// coefficients into `apu`.
///
/// The `u` control volume is staggered in `x`: it is centred on the
/// `u`-face and straddles the two pressure cells either side. The
/// convective mass flows on its four faces are interpolated from the
/// surrounding `u` / `v` faces; the pressure gradient driving it is the
/// exact difference of the two bracketing cell pressures (the staggered
/// grid's payoff).
#[allow(clippy::too_many_arguments)]
fn solve_u_momentum(
    u: &mut Field,
    v: &Field,
    p: &Field,
    apu: &mut Field,
    grid: &Grid,
    rho: f64,
    nu: f64,
    nu_t: Option<&Field>,
    bcs: &Boundaries,
    relax: f64,
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx();
    let dy = grid.dy();

    // Effective viscosity at the u-CV's four faces. With a turbulence
    // model the cell-centred ν_t is averaged onto each face; without
    // one, ν_eff = ν everywhere.
    //
    // The u-CV centre sits on a vertical face at (i, j) bracketed by
    // pressure cells (i-1, j) and (i, j). Its **east** face is at
    // pressure cell (i, j) — ν_eff_e = ν + ν_t(i, j). Its **west**
    // face is at pressure cell (i-1, j) — ν_eff_w = ν + ν_t(i-1, j).
    // Its **north** face is the corner straddling cells
    // (i-1, j), (i, j), (i-1, j+1), (i, j+1) — ν_eff is the average of
    // the four (or of the two on a boundary). Likewise for **south**.
    let nut_cell = |i: usize, j: usize| -> f64 {
        nu_t.map(|f| f.at(i, j).max(0.0)).unwrap_or(0.0)
    };

    for _sweep in 0..2 {
        for j in 0..ny {
            for i in 1..nx {
                // Convective mass flows on the four faces of the
                // u(i,j) control volume.
                let fe = rho * dy * 0.5 * (u.at(i, j) + u.at(i + 1, j));
                let fw = rho * dy * 0.5 * (u.at(i - 1, j) + u.at(i, j));
                let fn_ = rho * dx * 0.5 * (v.at(i - 1, j + 1) + v.at(i, j + 1));
                let fs = rho * dx * 0.5 * (v.at(i - 1, j) + v.at(i, j));

                // Face viscosities — ν + ν_t averaged to the face.
                let nu_e = nu + nut_cell(i, j);
                let nu_w = nu + nut_cell(i - 1, j);
                // North-face corner average.
                let mut nu_n = nu + 0.5 * (nut_cell(i - 1, j) + nut_cell(i, j));
                if j + 1 < ny {
                    nu_n = nu
                        + 0.25
                            * (nut_cell(i - 1, j)
                                + nut_cell(i, j)
                                + nut_cell(i - 1, j + 1)
                                + nut_cell(i, j + 1));
                }
                let mut nu_s = nu + 0.5 * (nut_cell(i - 1, j) + nut_cell(i, j));
                if j > 0 {
                    nu_s = nu
                        + 0.25
                            * (nut_cell(i - 1, j)
                                + nut_cell(i, j)
                                + nut_cell(i - 1, j - 1)
                                + nut_cell(i, j - 1));
                }
                let d_e = nu_e * dy / dx;
                let d_w = nu_w * dy / dx;
                let d_n = nu_n * dx / dy;
                let d_s = nu_s * dx / dy;

                let ae = hybrid_coeff(d_e, fe, -1.0);
                let aw = hybrid_coeff(d_w, fw, 1.0);

                let mut a_p = ae + aw;
                let mut su = 0.0;
                let mut an = 0.0;
                let mut as_ = 0.0;
                // North side.
                if j == ny - 1 {
                    let wall_u = match bcs.north {
                        SideBc::Wall { u, .. } | SideBc::Inlet { u, .. } => u,
                        SideBc::Outlet => u.at(i, j),
                    };
                    let dwall_n = nu_n * dx / (0.5 * dy);
                    a_p += dwall_n;
                    su += dwall_n * wall_u;
                } else {
                    an = hybrid_coeff(d_n, fn_, -1.0);
                    a_p += an;
                }
                // South side.
                if j == 0 {
                    let wall_u = match bcs.south {
                        SideBc::Wall { u, .. } | SideBc::Inlet { u, .. } => u,
                        SideBc::Outlet => u.at(i, j),
                    };
                    let dwall_s = nu_s * dx / (0.5 * dy);
                    a_p += dwall_s;
                    su += dwall_s * wall_u;
                } else {
                    as_ = hybrid_coeff(d_s, fs, 1.0);
                    a_p += as_;
                }

                if a_p.abs() < 1e-30 {
                    continue;
                }
                let dp = (p.at(i - 1, j) - p.at(i, j)) * dy;
                let mut nb = su + dp;
                nb += ae * u.at(i + 1, j);
                nb += aw * u.at(i - 1, j);
                if j + 1 < ny {
                    nb += an * u.at(i, j + 1);
                }
                if j > 0 {
                    nb += as_ * u.at(i, j - 1);
                }
                let u_new = nb / a_p;
                let u_old = u.at(i, j);
                u.set(i, j, u_old + relax * (u_new - u_old));
                apu.set(i, j, a_p / relax);
            }
        }
    }
}

/// Solve the discrete `v`-momentum equation — the `y`-direction mirror
/// of [`solve_u_momentum`]. The `v` control volume is staggered in `y`.
#[allow(clippy::too_many_arguments)]
fn solve_v_momentum(
    u: &Field,
    v: &mut Field,
    p: &Field,
    apv: &mut Field,
    grid: &Grid,
    rho: f64,
    nu: f64,
    nu_t: Option<&Field>,
    bcs: &Boundaries,
    relax: f64,
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx();
    let dy = grid.dy();

    let nut_cell = |i: usize, j: usize| -> f64 {
        nu_t.map(|f| f.at(i, j).max(0.0)).unwrap_or(0.0)
    };

    for _sweep in 0..2 {
        for j in 1..ny {
            for i in 0..nx {
                // Mass flows on the four faces of the v(i,j) CV.
                let fn_ = rho * dx * 0.5 * (v.at(i, j) + v.at(i, j + 1));
                let fs = rho * dx * 0.5 * (v.at(i, j - 1) + v.at(i, j));
                let fe = rho * dy * 0.5 * (u.at(i + 1, j - 1) + u.at(i + 1, j));
                let fw = rho * dy * 0.5 * (u.at(i, j - 1) + u.at(i, j));

                // The v-CV is centred on a horizontal face at (i, j),
                // bracketed by pressure cells (i, j-1) and (i, j).
                let nu_n = nu + nut_cell(i, j);
                let nu_s = nu + nut_cell(i, j - 1);
                // East corner: average of the four (or two on boundary).
                let mut nu_e = nu + 0.5 * (nut_cell(i, j - 1) + nut_cell(i, j));
                if i + 1 < nx {
                    nu_e = nu
                        + 0.25
                            * (nut_cell(i, j - 1)
                                + nut_cell(i, j)
                                + nut_cell(i + 1, j - 1)
                                + nut_cell(i + 1, j));
                }
                let mut nu_w = nu + 0.5 * (nut_cell(i, j - 1) + nut_cell(i, j));
                if i > 0 {
                    nu_w = nu
                        + 0.25
                            * (nut_cell(i, j - 1)
                                + nut_cell(i, j)
                                + nut_cell(i - 1, j - 1)
                                + nut_cell(i - 1, j));
                }
                let d_n = nu_n * dx / dy;
                let d_s = nu_s * dx / dy;
                let d_e = nu_e * dy / dx;
                let d_w = nu_w * dy / dx;

                let an = hybrid_coeff(d_n, fn_, -1.0);
                let as_ = hybrid_coeff(d_s, fs, 1.0);

                let mut a_p = an + as_;
                let mut sv = 0.0;
                let mut ae = 0.0;
                let mut aw = 0.0;
                // East side.
                if i == nx - 1 {
                    let wall_v = match bcs.east {
                        SideBc::Wall { v, .. } | SideBc::Inlet { v, .. } => v,
                        SideBc::Outlet => v.at(i, j),
                    };
                    let dwall_e = nu_e * dy / (0.5 * dx);
                    a_p += dwall_e;
                    sv += dwall_e * wall_v;
                } else {
                    ae = hybrid_coeff(d_e, fe, -1.0);
                    a_p += ae;
                }
                // West side.
                if i == 0 {
                    let wall_v = match bcs.west {
                        SideBc::Wall { v, .. } | SideBc::Inlet { v, .. } => v,
                        SideBc::Outlet => v.at(i, j),
                    };
                    let dwall_w = nu_w * dy / (0.5 * dx);
                    a_p += dwall_w;
                    sv += dwall_w * wall_v;
                } else {
                    aw = hybrid_coeff(d_w, fw, 1.0);
                    a_p += aw;
                }

                if a_p.abs() < 1e-30 {
                    continue;
                }
                let dp = (p.at(i, j - 1) - p.at(i, j)) * dx;
                let mut nb = sv + dp;
                nb += an * v.at(i, j + 1);
                nb += as_ * v.at(i, j - 1);
                if i + 1 < nx {
                    nb += ae * v.at(i + 1, j);
                }
                if i > 0 {
                    nb += aw * v.at(i - 1, j);
                }
                let v_new = nb / a_p;
                let v_old = v.at(i, j);
                v.set(i, j, v_old + relax * (v_new - v_old));
                apv.set(i, j, a_p / relax);
            }
        }
    }
}

/// Assemble the pressure-correction Poisson equation and return the L2
/// norm of the per-cell mass imbalance (the convergence residual).
///
/// SIMPLE's pressure correction relates a cell's `p'` to its neighbours
/// through the momentum-equation diagonals: a face velocity correction
/// is `u' = d·(p'_P − p'_E)` with `d = A/aP`. Substituting that into
/// the discrete continuity equation gives the five-point Poisson
/// stencil whose source `b` is exactly the provisional velocity's mass
/// imbalance. Driving `b → 0` is what makes the corrected velocity
/// divergence-free.
#[allow(clippy::too_many_arguments)]
fn assemble_pressure_correction(
    coeffs: &mut PoissonCoeffs,
    u: &Field,
    v: &Field,
    apu: &Field,
    apv: &Field,
    grid: &Grid,
    rho: f64,
    bcs: &Boundaries,
) -> f64 {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx();
    let dy = grid.dy();

    // The velocity-correction coefficient on a face: d = ρ·A²/aP. The
    // pressure-correction neighbour coefficient is then ρ·A·d.
    let mut sum_sq = 0.0;
    for j in 0..ny {
        for i in 0..nx {
            // East face (between cell i and i+1) — interior only.
            let ae = if i + 1 < nx {
                let ap_face = apu.at(i + 1, j).max(1e-30);
                rho * dy * dy / ap_face
            } else {
                0.0
            };
            // West face.
            let aw = if i > 0 {
                let ap_face = apu.at(i, j).max(1e-30);
                rho * dy * dy / ap_face
            } else {
                0.0
            };
            // North face.
            let an = if j + 1 < ny {
                let ap_face = apv.at(i, j + 1).max(1e-30);
                rho * dx * dx / ap_face
            } else {
                0.0
            };
            // South face.
            let as_ = if j > 0 {
                let ap_face = apv.at(i, j).max(1e-30);
                rho * dx * dx / ap_face
            } else {
                0.0
            };

            // The cell's net outflow of mass with the provisional
            // velocity: (ρ·u_e − ρ·u_w)·dy + (ρ·v_n − ρ·v_s)·dx. The
            // pressure-correction *source* is the negative of this
            // imbalance (continuity wants it driven to zero).
            let mass_out = rho * dy * (u.at(i + 1, j) - u.at(i, j))
                + rho * dx * (v.at(i, j + 1) - v.at(i, j));
            let b = -mass_out;

            coeffs.ae.set(i, j, ae);
            coeffs.aw.set(i, j, aw);
            coeffs.an.set(i, j, an);
            coeffs.as_.set(i, j, as_);
            coeffs.ap.set(i, j, ae + aw + an + as_);
            coeffs.b.set(i, j, b);

            sum_sq += mass_out * mass_out;
        }
    }

    // If an outlet exists, pin one outlet-boundary cell's correction to
    // zero so the otherwise-singular system has a unique (non-singular)
    // solution — the outlet is the natural pressure reference. The SOR
    // solver's zero-mean projection handles the all-Neumann (closed)
    // case instead; for an open domain it is told `pin_mean = false`.
    if has_outlet(bcs) {
        // Pick a cell that sits on the outlet boundary so `p' = 0` is
        // imposed where an outlet pressure reference physically belongs.
        let anchor = if matches!(bcs.east, SideBc::Outlet) {
            (nx - 1, ny / 2)
        } else if matches!(bcs.west, SideBc::Outlet) {
            (0, ny / 2)
        } else if matches!(bcs.north, SideBc::Outlet) {
            (nx / 2, ny - 1)
        } else {
            (nx / 2, 0) // south outlet
        };
        // A huge diagonal + a zero source forces p' ≈ 0 at the anchor.
        let big = 1e30;
        coeffs.ap.set(anchor.0, anchor.1, big);
        coeffs.ae.set(anchor.0, anchor.1, 0.0);
        coeffs.aw.set(anchor.0, anchor.1, 0.0);
        coeffs.an.set(anchor.0, anchor.1, 0.0);
        coeffs.as_.set(anchor.0, anchor.1, 0.0);
        coeffs.b.set(anchor.0, anchor.1, 0.0);
    }

    let n = (nx * ny) as f64;
    if n > 0.0 {
        (sum_sq / n).sqrt()
    } else {
        0.0
    }
}

/// True if any side of the domain is an outlet.
fn has_outlet(bcs: &Boundaries) -> bool {
    matches!(bcs.west, SideBc::Outlet)
        || matches!(bcs.east, SideBc::Outlet)
        || matches!(bcs.south, SideBc::Outlet)
        || matches!(bcs.north, SideBc::Outlet)
}

/// Apply the under-relaxed pressure correction: `p ← p + αp·p'`.
fn correct_pressure(p: &mut Field, pcorr: &Field, relax_p: f64) {
    for k in 0..p.data.len() {
        p.data[k] += relax_p * pcorr.data[k];
    }
}

/// Apply the velocity correction so the corrected field is
/// divergence-free.
///
/// Each interior face velocity is corrected by `u' = d·(p'_P − p'_E)`
/// with the velocity-correction coefficient `d = A/aP` — the face area
/// over the momentum-equation diagonal. The momentum equation's
/// pressure-force term is `Δp·A` (the `aP` already absorbs the density
/// through the convective mass fluxes), so `d` carries **no** explicit
/// `ρ`. The Poisson stencil in [`assemble_pressure_correction`] then
/// uses the *mass*-flux form `aE = ρ·A·d`; that shared `d` is exactly
/// what makes the corrected velocity satisfy the discrete continuity
/// equation.
fn correct_velocity(
    u: &mut Field,
    v: &mut Field,
    pcorr: &Field,
    apu: &Field,
    apv: &Field,
    grid: &Grid,
    _bcs: &Boundaries,
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let dy = grid.dy();
    let dx = grid.dx();

    // Interior u-faces: i = 1 .. nx-1. The velocity-correction
    // coefficient is d = A/aP = dy/apu (no ρ — see the doc comment).
    for j in 0..ny {
        for i in 1..nx {
            let ap = apu.at(i, j).max(1e-30);
            let d = dy / ap;
            // u' = d·(p'_W − p'_E) where W = cell i-1, E = cell i.
            let du = d * (pcorr.at(i - 1, j) - pcorr.at(i, j));
            u.add(i, j, du);
        }
    }
    // Interior v-faces: j = 1 .. ny-1.
    for j in 1..ny {
        for i in 0..nx {
            let ap = apv.at(i, j).max(1e-30);
            let d = dx / ap;
            let dv = d * (pcorr.at(i, j - 1) - pcorr.at(i, j));
            v.add(i, j, dv);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lid_driven_cavity_develops_a_circulation() {
        // The canonical CFD verification case. A square cavity with a
        // sliding top lid drives a single large recirculating vortex:
        // the lid drags the top fluid in +x, it turns down the right
        // wall, returns in −x along the floor, and rises up the left
        // wall. The defining signature is that the flow *reverses
        // direction* between the top and the bottom — the floor fluid
        // moves opposite to the lid.
        let grid = Grid::new(24, 24, 1.0, 1.0);
        let fluid = Fluid::new(1.0, 0.05); // Re = U·L/ν = 1·1/0.05 = 20
        let bcs = Boundaries::lid_driven_cavity(1.0);
        let controls = SimpleControls {
            max_iterations: 3000,
            tolerance: 1e-5,
            ..SimpleControls::default()
        };
        let sol = solve_simple(&grid, &fluid, &bcs, &controls);
        assert!(
            sol.converged,
            "lid-driven cavity should converge (residual {})",
            sol.residual
        );

        // Sample u on the vertical centreline. Near the lid u must be
        // positive (dragged along by the lid); near the floor it must
        // be negative (the return flow of the vortex).
        let mid_i = grid.nx / 2;
        let u_near_lid = sol.u_at_cell(mid_i, grid.ny - 2);
        let u_near_floor = sol.u_at_cell(mid_i, 1);
        assert!(
            u_near_lid > 0.05,
            "fluid under the lid should move with it, u = {u_near_lid}"
        );
        assert!(
            u_near_floor < -0.01,
            "fluid near the floor should return upstream, u = {u_near_floor}"
        );
    }

    #[test]
    fn lid_driven_cavity_velocity_is_bounded_by_the_lid_speed() {
        // Energy sanity: nothing in the cavity can move faster than the
        // lid that drives it — the lid is the only momentum source.
        let grid = Grid::new(20, 20, 1.0, 1.0);
        let fluid = Fluid::new(1.0, 0.1);
        let lid = 1.0;
        let bcs = Boundaries::lid_driven_cavity(lid);
        let sol = solve_simple(&grid, &fluid, &bcs, &SimpleControls::default());
        let max_speed = sol.max_speed();
        assert!(
            max_speed <= lid * 1.15,
            "interior speed {max_speed} should not exceed the lid speed {lid}"
        );
    }

    #[test]
    fn channel_flow_develops_a_parabolic_profile() {
        // Plane Poiseuille flow: fluid enters a long channel with a
        // flat (uniform) velocity profile; viscosity drags it to zero
        // at the walls, and far downstream the profile becomes the
        // analytic parabola — fast in the middle, zero at the walls.
        // A long, thin channel ensures the flow is fully developed
        // before the outlet.
        let grid = Grid::new(60, 16, 6.0, 1.0);
        let fluid = Fluid::new(1.0, 0.05);
        let inlet = 1.0;
        let bcs = Boundaries::channel_flow(inlet);
        let controls = SimpleControls {
            max_iterations: 4000,
            tolerance: 1e-5,
            ..SimpleControls::default()
        };
        let sol = solve_simple(&grid, &fluid, &bcs, &controls);
        assert!(
            sol.converged,
            "channel flow should converge (residual {})",
            sol.residual
        );

        // Sample the downstream u-profile across the channel height.
        let out_i = grid.nx - 3; // a column near the outlet
        let centre_j = grid.ny / 2;
        let u_centre = sol.u_at_cell(out_i, centre_j);
        let u_near_wall = sol.u_at_cell(out_i, 0);
        // The profile must be fastest at the centre and slow at the
        // wall — the parabola's defining shape.
        assert!(
            u_centre > u_near_wall + 0.1,
            "developed profile must peak in the centre: centre {u_centre}, wall {u_near_wall}"
        );
        // Mass conservation: the developed centreline speed of a
        // parabolic profile is 1.5× the mean (the inlet) speed.
        assert!(
            u_centre > 1.2 * inlet && u_centre < 1.8 * inlet,
            "developed centreline u {u_centre} should be ~1.5× the inlet {inlet}"
        );
        // The profile must be roughly symmetric about the channel
        // centreline.
        let u_low = sol.u_at_cell(out_i, 2);
        let u_high = sol.u_at_cell(out_i, grid.ny - 3);
        assert!(
            (u_low - u_high).abs() < 0.2 * u_centre,
            "profile should be symmetric: low {u_low}, high {u_high}"
        );
    }

    #[test]
    fn channel_flow_has_a_streamwise_pressure_drop() {
        // Pushing a viscous flow down a channel needs a pressure drop, so the
        // static-pressure range (gauge-independent max − min) must be strictly
        // positive and finite.
        let grid = Grid::new(60, 16, 6.0, 1.0);
        let fluid = Fluid::new(1.0, 0.05);
        let bcs = Boundaries::channel_flow(1.0);
        let controls = SimpleControls {
            max_iterations: 4000,
            tolerance: 1e-5,
            ..SimpleControls::default()
        };
        let sol = solve_simple(&grid, &fluid, &bcs, &controls);
        let dp = sol.pressure_range();
        assert!(dp.is_finite() && dp > 0.0, "channel flow needs a pressure drop, got {dp}");
        // It matches an independent max − min scan of the cell pressures.
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let p = sol.pressure.at(i, j);
                lo = lo.min(p);
                hi = hi.max(p);
            }
        }
        assert!((dp - (hi - lo)).abs() < 1e-12, "Δp {dp} vs max−min {}", hi - lo);
    }

    #[test]
    fn channel_flow_mean_speed_is_positive_and_below_the_peak() {
        let grid = Grid::new(60, 16, 6.0, 1.0);
        let fluid = Fluid::new(1.0, 0.05);
        let bcs = Boundaries::channel_flow(1.0);
        let controls = SimpleControls {
            max_iterations: 4000,
            tolerance: 1e-5,
            ..SimpleControls::default()
        };
        let sol = solve_simple(&grid, &fluid, &bcs, &controls);
        let mean = sol.mean_speed();
        let peak = sol.max_speed();
        // A driven flow has a positive mean that cannot exceed the peak.
        assert!(mean > 0.0, "driven flow should have a positive mean speed");
        assert!(mean <= peak + 1e-9, "mean {mean} cannot exceed peak {peak}");
    }

    #[test]
    fn flow_through_time_is_the_convective_residence_time() {
        // Grid 5×4 over 5×8 → Lₓ = nx·dx = 5·1 = 5; uniform u = 2 → max_speed = 2.
        let grid = Grid::new(5, 4, 5.0, 8.0);
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, 2.0);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // Worked: τ = Lₓ/U = 5/2 = 2.5 s.
        assert!((sol.flow_through_time() - 2.5).abs() <= 1e-9 * 2.5, "τ = Lₓ/U = 2.5");

        // Threads max_speed: τ·U = Lₓ.
        let lx = grid.nx as f64 * grid.dx();
        assert!(
            (sol.flow_through_time() * sol.max_speed() - lx).abs() <= 1e-9 * lx,
            "τ·U = Lₓ"
        );

        // Threads domain_reynolds_number: τ = Re_L·ν / U² (since Re_L·ν = U·Lₓ).
        assert!(
            (sol.flow_through_time()
                - sol.domain_reynolds_number(0.01) * 0.01 / (sol.max_speed() * sol.max_speed()))
            .abs()
                <= 1e-9 * sol.flow_through_time(),
            "τ = Re_L·ν/U²"
        );

        // Inverse in U: a faster field gives half the flow-through time.
        let mut u4 = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u4.set(fi, j, 4.0);
            }
        }
        let fast = FlowSolution {
            grid,
            u: u4,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (fast.flow_through_time() - 0.5 * sol.flow_through_time()).abs()
                <= 1e-9 * sol.flow_through_time(),
            "inverse in U"
        );

        // 0 sentinel for a quiescent (all-zero velocity) field.
        let quiescent = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(quiescent.flow_through_time(), 0.0, "quiescent → 0");
    }

    #[test]
    fn domain_reynolds_number_is_the_flow_regime_classifier() {
        // Grid 5×4 over 5×8 → Lₓ = nx·dx = 5·1 = 5; uniform u = 2 → max_speed = 2.
        let grid = Grid::new(5, 4, 5.0, 8.0);
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, 2.0);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // Worked: Re_L = U·Lₓ/ν = 2·5/0.01 = 1000.
        assert!(
            (sol.domain_reynolds_number(0.01) - 1000.0).abs() <= 1e-9 * 1000.0,
            "Re_L = U·Lₓ/ν"
        );

        // Threads max_speed: Re_L = max_speed · Lₓ / ν.
        let lx = grid.nx as f64 * grid.dx();
        assert!(
            (sol.domain_reynolds_number(0.01) - sol.max_speed() * lx / 0.01).abs()
                <= 1e-12 * sol.domain_reynolds_number(0.01),
            "Re_L via max_speed"
        );

        // Threads cell_reynolds_number: Re_L / Re_cell = Lₓ / h (h = min cell size).
        let h = grid.dx().min(grid.dy());
        assert!(
            (sol.domain_reynolds_number(0.01) - sol.cell_reynolds_number(0.01) * (lx / h)).abs()
                <= 1e-9 * sol.domain_reynolds_number(0.01),
            "Re_L = Re_cell · Lₓ/h"
        );

        // Inverse in ν.
        assert!(
            (sol.domain_reynolds_number(0.02) - 0.5 * sol.domain_reynolds_number(0.01)).abs()
                <= 1e-9 * sol.domain_reynolds_number(0.01),
            "inverse in ν"
        );

        // Linear in U: doubling the speed doubles Re_L.
        let mut u4 = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u4.set(fi, j, 4.0);
            }
        }
        let fast = FlowSolution {
            grid,
            u: u4,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (fast.domain_reynolds_number(0.01) - 2.0 * sol.domain_reynolds_number(0.01)).abs()
                <= 1e-9 * fast.domain_reynolds_number(0.01),
            "linear in U"
        );

        // 0 sentinel for non-physical viscosity.
        assert_eq!(sol.domain_reynolds_number(0.0), 0.0);
        assert_eq!(sol.domain_reynolds_number(-1.0), 0.0);
        assert_eq!(sol.domain_reynolds_number(f64::NAN), 0.0);
    }

    #[test]
    fn cell_reynolds_number_is_the_cfl_to_diffusion_ratio() {
        // Grid 5×4 over 5×8 → dx = 1.0, dy = 2.0, h = min = 1.0; uniform u = 3.
        let grid = Grid::new(5, 4, 5.0, 8.0);
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, 3.0);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let h = grid.dx().min(grid.dy());
        let nu = 0.1;

        // Definition + worked value: Re_h = U_max·h/ν = 3·1/0.1 = 30.
        assert!(
            (sol.cell_reynolds_number(nu) - sol.max_speed() * h / nu).abs() < 1e-12,
            "Re_h = U·h/ν"
        );
        assert!((sol.cell_reynolds_number(0.1) - 30.0).abs() < 1e-12);

        // Threads both stability numbers: Re_h = CFL(dt)/d(ν,dt), independent of dt.
        for &dt in &[0.01, 0.05, 0.2] {
            let ratio = sol.convective_cfl_number(dt) / sol.diffusive_number(nu, dt);
            assert!((sol.cell_reynolds_number(nu) - ratio).abs() <= 1e-9 * ratio, "Re_h = CFL/d");
        }

        // Linear in 1/ν: doubling ν halves Re_h.
        assert!(
            (sol.cell_reynolds_number(0.2) - 0.5 * sol.cell_reynolds_number(0.1)).abs() < 1e-9,
            "Re_h ∝ 1/ν"
        );

        // Sentinel: zero / negative / non-finite viscosity → 0.
        assert_eq!(sol.cell_reynolds_number(0.0), 0.0);
        assert_eq!(sol.cell_reynolds_number(-0.1), 0.0);
        assert_eq!(sol.cell_reynolds_number(f64::NAN), 0.0);
    }

    #[test]
    fn diffusive_number_pairs_with_the_convective_cfl() {
        // Grid 5×4 over 5×8 → dx = 1.0, dy = 2.0, h = min = 1.0; uniform u = 3.
        let grid = Grid::new(5, 4, 5.0, 8.0);
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, 3.0);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let h = grid.dx().min(grid.dy());
        let nu = 0.1;
        let dt = 0.05;

        // Definition + worked value: d = ν·Δt/h² = 0.1·0.05/1 = 0.005.
        assert!((sol.diffusive_number(nu, dt) - nu * dt / (h * h)).abs() < 1e-12, "d = ν·Δt/h²");
        assert!((sol.diffusive_number(0.1, 0.05) - 0.005).abs() < 1e-12);

        // Threads convective_cfl_number via the cell Reynolds number: CFL/d = U_max·h/ν.
        let re_cell = sol.max_speed() * h / nu;
        assert!(
            (sol.convective_cfl_number(dt) - sol.diffusive_number(nu, dt) * re_cell).abs() < 1e-12,
            "CFL = d · Re_cell"
        );

        // Linear in ν; the d = 1/4 stability limit at Δt = h²/(4ν).
        assert!(
            (sol.diffusive_number(0.2, dt) - 2.0 * sol.diffusive_number(0.1, dt)).abs() < 1e-12,
            "linear in ν"
        );
        assert!(
            (sol.diffusive_number(nu, h * h / (4.0 * nu)) - 0.25).abs() < 1e-12,
            "d = 1/4 at Δt = h²/4ν"
        );

        // Sentinel: zero / negative / non-finite inputs → 0.
        assert_eq!(sol.diffusive_number(nu, 0.0), 0.0);
        assert_eq!(sol.diffusive_number(-0.1, dt), 0.0);
        assert_eq!(sol.diffusive_number(nu, f64::NAN), 0.0);
    }

    #[test]
    fn convective_cfl_number_threads_max_speed() {
        // Grid 5×4 over 5×8 → dx = 1.0, dy = 2.0, so h = min = 1.0.
        let grid = Grid::new(5, 4, 5.0, 8.0);
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, 3.0); // uniform u = 3 → max_speed = 3
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let h = grid.dx().min(grid.dy());

        // Threads max_speed: C = U_max·Δt / min(Δx, Δy).
        assert!(
            (sol.convective_cfl_number(0.1) - sol.max_speed() * 0.1 / h).abs() < 1e-12,
            "C = U_max·Δt/h"
        );
        // Worked: max_speed = 3, Δt = 0.1, h = 1 → C = 0.3.
        assert!((sol.convective_cfl_number(0.1) - 0.3).abs() < 1e-12, "C(0.1) = 0.3");
        // Linear in Δt.
        assert!(
            (sol.convective_cfl_number(0.2) - 2.0 * sol.convective_cfl_number(0.1)).abs() < 1e-12,
            "linear in Δt"
        );
        // The unit-CFL stability limit: Δt = h/U_max → C = 1.
        assert!(
            (sol.convective_cfl_number(h / sol.max_speed()) - 1.0).abs() < 1e-12,
            "C = 1 at Δt = h/U_max"
        );
        // Sentinel: zero / negative / non-finite Δt → 0.
        assert_eq!(sol.convective_cfl_number(0.0), 0.0);
        assert_eq!(sol.convective_cfl_number(-0.1), 0.0);
        assert_eq!(sol.convective_cfl_number(f64::NAN), 0.0);
    }

    #[test]
    fn max_speed_is_the_peak_velocity_magnitude() {
        // A uniform field u = 3, v = 4 gives every cell the speed √(3²+4²) = 5,
        // so the peak equals the mean equals 5.
        let grid = Grid::new(3, 2, 3.0, 2.0);
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                u.set(i, j, 3.0);
            }
        }
        let mut v = grid.v_field();
        for j in 0..=grid.ny {
            for i in 0..grid.nx {
                v.set(i, j, 4.0);
            }
        }
        let uniform = FlowSolution {
            grid,
            u,
            v,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((uniform.max_speed() - 5.0).abs() < 1e-12, "uniform peak speed = 5");
        assert!(
            (uniform.max_speed() - uniform.mean_speed()).abs() < 1e-12,
            "a uniform field has peak == mean"
        );

        // A non-uniform field whose cell speeds are {5, 10, 5}: v = 0 and the
        // middle cell's two shared u-faces both set to 10 (u_at_cell averages the
        // two faces), so cell 1 reaches 10 while its neighbours sit at 5. The peak
        // picks out the fast cell (10), strictly above the mean (20/3 ≈ 6.67).
        let g2 = Grid::new(3, 1, 3.0, 1.0);
        let mut u2 = g2.u_field();
        u2.set(1, 0, 10.0); // right face of cell 0 / left face of cell 1
        u2.set(2, 0, 10.0); // right face of cell 1 / left face of cell 2
        let peaked = FlowSolution {
            grid: g2,
            u: u2,
            v: g2.v_field(),
            pressure: g2.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (peaked.max_speed() - 10.0).abs() < 1e-12,
            "peak = 10, got {}",
            peaked.max_speed()
        );
        assert!(
            peaked.max_speed() > peaked.mean_speed(),
            "the peak strictly exceeds the mean for a non-uniform field"
        );
    }

    #[test]
    fn min_speed_is_the_slowest_cell_and_bounds_the_mean() {
        // Uniform field u = 3 → zero spread: min = mean = max.
        let grid = Grid::new(5, 4, 5.0, 8.0);
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, 3.0);
            }
        }
        let uniform = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((uniform.min_speed() - 3.0).abs() <= 1e-9 * 3.0, "uniform → min = 3");
        assert!((uniform.min_speed() - uniform.max_speed()).abs() <= 1e-9, "uniform: min = max");
        assert!((uniform.min_speed() - uniform.mean_speed()).abs() <= 1e-9, "uniform: min = mean");

        // Non-uniform field (streamwise gradient).
        let mut ug = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                ug.set(fi, j, fi as f64);
            }
        }
        let varied = FlowSolution {
            grid,
            u: ug,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // Statistical ordering min ≤ mean ≤ max, strict for a varying field (threads both).
        assert!(varied.min_speed() <= varied.mean_speed(), "min ≤ mean");
        assert!(varied.mean_speed() <= varied.max_speed(), "mean ≤ max");
        assert!(varied.min_speed() < varied.max_speed(), "min < max for a varying field");

        // It is the actual minimum: ≤ every cell speed, and non-negative.
        let m = varied.min_speed();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                assert!(m <= varied.speed_at_cell(i, j) + 1e-12, "≤ every cell");
            }
        }
        assert!(m >= 0.0, "non-negative");
    }

    #[test]
    fn speed_range_is_the_velocity_peak_to_trough() {
        // (a) WORKED: a 3×1 grid (v = 0), the middle cell's two shared u-faces set to
        // 10 → cell speeds {5, 10, 5}, so the speed range = 10 − 5 = 5.
        let g = Grid::new(3, 1, 3.0, 1.0);
        let mut u = g.u_field();
        u.set(1, 0, 10.0);
        u.set(2, 0, 10.0);
        let sol = FlowSolution {
            grid: g,
            u,
            v: g.v_field(),
            pressure: g.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let sr = sol.speed_range();
        assert!((sr - 5.0).abs() <= 1e-12 * 5.0, "speed_range = 10 − 5 = 5, got {sr}");

        // (b) THREAD max_speed + min_speed (non-tautological — the body is one loop,
        // these are separate loops): speed_range == max_speed − min_speed.
        assert!(
            (sr - (sol.max_speed() - sol.min_speed())).abs() <= 1e-12 * sr,
            "speed_range = max_speed − min_speed"
        );

        // (c) UNIFORM: a uniform u = 4 field → all cell speeds 4 → range 0.
        let grid = Grid::new(5, 4, 5.0, 8.0);
        let mut u4 = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u4.set(fi, j, 4.0);
            }
        }
        let uniform = FlowSolution {
            grid,
            u: u4,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(uniform.speed_range(), 0.0, "uniform field → range 0");

        // (d) NON-NEGATIVE + BOUND: range ≥ 0 always, and ≤ max_speed (min_speed ≥ 0).
        assert!(sol.speed_range() >= 0.0, "non-negative");
        assert!(sol.speed_range() <= sol.max_speed(), "range ≤ max_speed");

        // (e) MONOTONE SPREAD: a narrower field has a smaller range. Middle faces set
        // to 6 give cell speeds {3, 6, 3} (range 3) vs the {5, 10, 5} field's range 5.
        let mut u2 = g.u_field();
        u2.set(1, 0, 6.0);
        u2.set(2, 0, 6.0);
        let narrow = FlowSolution {
            grid: g,
            u: u2,
            v: g.v_field(),
            pressure: g.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(narrow.speed_range() < sol.speed_range(), "narrower spread → smaller range");
    }

    #[test]
    fn max_speed_location_finds_the_fastest_cell() {
        // A 3×1 grid (dx = dy = 1), v = 0, the middle cell's two shared u-faces
        // set to 10 → cell speeds {5, 10, 5}; the fastest is cell (1, 0), whose
        // centre is ((1+0.5)·dx, (0+0.5)·dy) = (1.5, 0.5).
        let grid = Grid::new(3, 1, 3.0, 1.0);
        let mut u = grid.u_field();
        u.set(1, 0, 10.0);
        u.set(2, 0, 10.0);
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let (x, y) = sol.max_speed_location().expect("a non-empty grid has a fastest cell");
        assert!(
            (x - 1.5).abs() < 1e-12 && (y - 0.5).abs() < 1e-12,
            "fastest cell centre (1.5, 0.5), got ({x}, {y})"
        );
        // The located cell's speed is indeed the peak (ties to max_speed).
        assert!((sol.speed_at_cell(1, 0) - sol.max_speed()).abs() < 1e-12, "cell (1,0) is the peak");

        // A 2×2 grid (dx = dy = 1) with the peak planted at cell (1, 1) → centre
        // (1.5, 1.5) — confirms the locator picks the right cell in 2-D.
        let g2 = Grid::new(2, 2, 2.0, 2.0);
        let mut u2 = g2.u_field();
        u2.set(1, 1, 8.0);
        u2.set(2, 1, 8.0);
        let sol2 = FlowSolution {
            grid: g2,
            u: u2,
            v: g2.v_field(),
            pressure: g2.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let (x2, y2) = sol2.max_speed_location().unwrap();
        assert!(
            (x2 - 1.5).abs() < 1e-12 && (y2 - 1.5).abs() < 1e-12,
            "fastest cell (1,1) centre (1.5, 1.5), got ({x2}, {y2})"
        );
    }

    #[test]
    fn lid_driven_cavity_has_positive_vorticity_that_scales_with_the_lid() {
        let grid = Grid::new(24, 24, 1.0, 1.0);
        let fluid = Fluid::new(1.0, 0.05);
        let controls = SimpleControls {
            max_iterations: 2000,
            tolerance: 1e-5,
            ..SimpleControls::default()
        };
        let slow = solve_simple(&grid, &fluid, &Boundaries::lid_driven_cavity(1.0), &controls);
        let fast = solve_simple(&grid, &fluid, &Boundaries::lid_driven_cavity(2.0), &controls);
        let w_slow = slow.max_vorticity();
        let w_fast = fast.max_vorticity();
        // The lid shears the fluid into a rotating vortex: vorticity is positive
        // and finite, and driving the lid harder spins it up.
        assert!(w_slow > 0.0 && w_slow.is_finite(), "cavity vorticity {w_slow}");
        assert!(w_fast > w_slow, "a faster lid raises vorticity: {w_slow} → {w_fast}");
    }

    #[test]
    fn peak_vorticity_location_finds_the_strongest_shear_cell() {
        // A 5×5 unit grid (dx = dy = 1). With v = 0 and u varying only with
        // height as U(j) = j², the vorticity ω = −∂u/∂y grows with j, so the
        // strongest shear is at the top interior row.
        let grid = Grid::new(5, 5, 5.0, 5.0);
        let mut u = grid.u_field(); // 6 × 5
        for j in 0..grid.ny {
            let val = (j * j) as f64;
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let (x, y) = sol.peak_vorticity_location().unwrap();
        // Peak central-difference shear is at interior row j = 3 → y = 3.5;
        // it appears first at the leftmost interior column i = 1 → x = 1.5.
        assert!((y - 3.5).abs() < 1e-12, "peak vort y {y}");
        assert!((x - 1.5).abs() < 1e-12, "peak vort x {x}");
        // A grid too small for an interior difference → None.
        let tg = Grid::new(2, 2, 1.0, 1.0);
        let tiny = FlowSolution {
            grid: tg,
            u: tg.u_field(),
            v: tg.v_field(),
            pressure: tg.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(tiny.peak_vorticity_location().is_none());
    }

    #[test]
    fn enstrophy_of_a_uniform_shear_matches_the_analytic_value() {
        // A pure horizontal shear u(y) = γ·y, v = 0 has a *constant* vorticity
        // ω = ∂v/∂x − ∂u/∂y = −γ everywhere, so its enstrophy is exactly
        // ½·γ²·(area the interior stencil covers) = ½·γ²·(nx−2)(ny−2)·dx·dy.
        let gamma = 2.0_f64;
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let mut u = grid.u_field(); // 6 × 5
        for j in 0..grid.ny {
            // u_at_cell(i, j) = γ·y_cell = γ·(j + 0.5)·dy → set every u-face on the row.
            let val = gamma * (j as f64 + 0.5) * grid.dy();
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let interior = ((grid.nx - 2) * (grid.ny - 2)) as f64;
        let expected = 0.5 * gamma * gamma * interior * grid.dx() * grid.dy();
        let e = sol.enstrophy();
        assert!((e - expected).abs() < 1e-9, "enstrophy {e} vs analytic {expected}");
        // Strictly positive for a rotational field (and never negative — it is ½∑ω²).
        assert!(e > 0.0);

        // A grid too small for an interior central difference → 0.
        let tg = Grid::new(2, 2, 1.0, 1.0);
        let tiny = FlowSolution {
            grid: tg,
            u: tg.u_field(),
            v: tg.v_field(),
            pressure: tg.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(tiny.enstrophy(), 0.0);
    }

    #[test]
    fn mean_vorticity_is_the_area_averaged_interior_vorticity() {
        // (a) ANALYTIC: a pure horizontal shear u(y) = γ·y, v = 0 has a constant
        // vorticity ω = ∂v/∂x − ∂u/∂y = −γ everywhere, so its area-average is −γ.
        let gamma = 2.0_f64;
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * grid.dy();
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((sol.mean_vorticity() - (-gamma)).abs() < 1e-9, "⟨ω⟩ = −γ = −2");

        // (b) STOKES cross-check threading circulation (non-tautological — a different
        // code path): the area-average equals the circulation per unit interior area.
        let interior_area = (grid.nx - 2) as f64 * (grid.ny - 2) as f64 * grid.dx() * grid.dy();
        assert!(
            (sol.mean_vorticity() - sol.circulation() / interior_area).abs() <= 1e-9 * gamma,
            "⟨ω⟩ = Γ/A"
        );

        // (c) IRROTATIONAL: a uniform (zero) field has zero mean vorticity.
        let uni = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(uni.mean_vorticity(), 0.0, "uniform → ⟨ω⟩ = 0");

        // (d) SMALL GRID: too small for an interior central difference → 0.
        let tg = Grid::new(2, 2, 1.0, 1.0);
        let tiny = FlowSolution {
            grid: tg,
            u: tg.u_field(),
            v: tg.v_field(),
            pressure: tg.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(tiny.mean_vorticity(), 0.0);
    }

    #[test]
    fn rms_vorticity_is_the_root_mean_square_interior_vorticity() {
        // A pure horizontal shear u(y) = γ·y has constant vorticity ω = −γ everywhere, so
        // its rms is exactly |ω| = γ.
        let gamma = 2.0_f64;
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * grid.dy();
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) ANALYTIC: ω_rms = |ω| = γ for the constant-vorticity field.
        assert!((sol.rms_vorticity() - gamma).abs() < 1e-9, "ω_rms = γ = 2");

        // (b) THREAD enstrophy (non-tautological — a different code path): ω_rms = √(2E/A).
        let interior_area = (grid.nx - 2) as f64 * (grid.ny - 2) as f64 * grid.dx() * grid.dy();
        assert!(
            (sol.rms_vorticity() - (2.0 * sol.enstrophy() / interior_area).sqrt()).abs()
                <= 1e-9 * sol.rms_vorticity(),
            "ω_rms = √(2E/A)"
        );

        // (c) RMS ≥ |MEAN| (Cauchy–Schwarz), with equality for this uniform field.
        assert!(sol.rms_vorticity() >= sol.mean_vorticity().abs() - 1e-12, "rms ≥ |mean|");
        assert!(
            (sol.rms_vorticity() - sol.mean_vorticity().abs()).abs() < 1e-9,
            "uniform ω → rms = |mean|"
        );

        // (d) SMALL GRID: too small for an interior central difference → 0.
        let tg = Grid::new(2, 2, 1.0, 1.0);
        let tiny = FlowSolution {
            grid: tg,
            u: tg.u_field(),
            v: tg.v_field(),
            pressure: tg.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(tiny.rms_vorticity(), 0.0);
    }

    #[test]
    fn vorticity_std_dev_completes_the_variance_identity() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1

        // (a) UNIFORM shear u(y) = γ·y has constant vorticity ω = −γ, so its spatial
        // standard deviation is exactly 0.
        let gamma = 2.0_f64;
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * grid.dy();
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let uniform = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(uniform.vorticity_std_dev() < 1e-9, "constant ω → σ = 0");

        // (b) QUADRATIC shear u(y) = a·y² has ω = −2a·y, varying in space, so σ > 0 — and the
        // variance identity σ² = ω_rms² − ⟨ω⟩² threads mean_vorticity + rms_vorticity (three
        // separate loops, non-tautological).
        let a = 0.5_f64;
        let mut u2 = grid.u_field();
        for j in 0..grid.ny {
            let val = a * ((j as f64 + 0.5) * grid.dy()).powi(2);
            for i in 0..=grid.nx {
                u2.set(i, j, val);
            }
        }
        let quad = FlowSolution {
            grid,
            u: u2,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let sigma = quad.vorticity_std_dev();
        assert!(sigma > 0.0, "varying ω → σ > 0");
        let identity = (quad.rms_vorticity().powi(2) - quad.mean_vorticity().powi(2)).sqrt();
        assert!((sigma - identity).abs() <= 1e-9 * sigma, "σ = √(ω_rms² − ⟨ω⟩²)");

        // (c) SMALL GRID: too small for an interior central difference → 0.
        let tg = Grid::new(2, 2, 1.0, 1.0);
        let tiny = FlowSolution {
            grid: tg,
            u: tg.u_field(),
            v: tg.v_field(),
            pressure: tg.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(tiny.vorticity_std_dev(), 0.0);
    }

    #[test]
    fn max_q_criterion_tells_a_vortex_from_a_shear_layer() {
        // The whole point of Q: vorticity alone cannot tell a vortex from a
        // straight shear layer — Q can. Solid-body rotation u=−Ωy, v=Ωx has
        // Q=Ω² everywhere; a pure shear u=γy, v=0 (with the *same* nonzero
        // vorticity and enstrophy) has Q=0 — it is not a vortex.
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // --- solid-body rotation at rate Ω ---
        let omega = 2.0_f64;
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy; // u_at_cell = −Ω·y
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let mut v = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx; // v_at_cell = Ω·x
            for j in 0..=grid.ny {
                v.set(i, j, val);
            }
        }
        let rot = FlowSolution {
            grid,
            u,
            v,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let q_rot = rot.max_q_criterion();
        assert!((q_rot - omega * omega).abs() < 1e-9, "rotation Q {q_rot} vs Ω² {}", omega * omega);

        // --- pure shear u = γ·y, v = 0 ---
        let gamma = 3.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(), // v ≡ 0
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        // Q rejects the shear (Q=0) even though its enstrophy is clearly non-zero.
        assert!(shear.max_q_criterion() < 1e-9, "shear Q {}", shear.max_q_criterion());
        assert!(shear.enstrophy() > 0.0, "the shear still carries enstrophy");

        // A grid too small for an interior central difference → 0.
        let tg = Grid::new(2, 2, 1.0, 1.0);
        let tiny = FlowSolution {
            grid: tg,
            u: tg.u_field(),
            v: tg.v_field(),
            pressure: tg.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(tiny.max_q_criterion(), 0.0);
    }

    #[test]
    fn viscous_dissipation_is_nonzero_for_a_vorticity_free_straining_flow() {
        // The signature that separates dissipation from the rotation diagnostics:
        // a pure (irrotational) straining flow u = ε·x, v = −ε·y has ZERO
        // vorticity and zero enstrophy, yet still dissipates — S:S = 2ε², so the
        // integral is 2ε²·(interior area).
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());
        let eps = 2.0_f64;
        let mu = 0.5_f64;
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            // u at the x-face of index i sits at x = i·dx → u = ε·i·dx.
            for i in 0..=grid.nx {
                u.set(i, j, eps * i as f64 * dx);
            }
        }
        let mut v = grid.v_field();
        for j in 0..=grid.ny {
            // v at the y-face of index j sits at y = j·dy → v = −ε·j·dy.
            for i in 0..grid.nx {
                v.set(i, j, -eps * j as f64 * dy);
            }
        }
        let strain = FlowSolution {
            grid,
            u,
            v,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            strain.enstrophy() < 1e-9,
            "straining flow is irrotational: enstrophy {}",
            strain.enstrophy()
        );
        let interior = ((grid.nx - 2) * (grid.ny - 2)) as f64;
        // Φ = 2μ·(S:S = 2ε²)·area.
        let expected = 2.0 * mu * 2.0 * eps * eps * interior * dx * dy;
        let phi = strain.viscous_dissipation(mu);
        assert!((phi - expected).abs() < 1e-9, "dissipation {phi} vs analytic {expected}");
        assert!(phi > 0.0, "a straining flow dissipates despite zero enstrophy");

        // A uniform shear u = γ·y gives S:S = ½γ², so Φ = 2μ·½γ²·area = μγ²·area.
        let gamma = 3.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let expected_shear = mu * gamma * gamma * interior * dx * dy;
        assert!(
            (shear.viscous_dissipation(mu) - expected_shear).abs() < 1e-9,
            "shear dissipation {}",
            shear.viscous_dissipation(mu)
        );

        // A grid too small for an interior central difference → 0.
        let tg = Grid::new(2, 2, 1.0, 1.0);
        let tiny = FlowSolution {
            grid: tg,
            u: tg.u_field(),
            v: tg.v_field(),
            pressure: tg.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(tiny.viscous_dissipation(mu), 0.0);
    }

    #[test]
    fn palinstrophy_of_a_quadratic_shear_matches_the_analytic_value() {
        // A quadratic shear u = ½γy², v = 0 has ω = −∂u/∂y = −γy (linear), so
        // ∇ω = (0, −γ) is constant and the palinstrophy is exactly
        // ½γ²·(area the doubly-interior stencil covers) = ½γ²·(nx−4)(ny−4)·dx·dy.
        let gamma = 2.0_f64;
        let grid = Grid::new(7, 7, 7.0, 7.0); // dx = dy = 1
        let dy = grid.dy();
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            // u_at_cell(i, j) = ½γ·y_cell² with y_cell = (j + 0.5)·dy.
            let y = (j as f64 + 0.5) * dy;
            let val = 0.5 * gamma * y * y;
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let interior = ((grid.nx - 4) * (grid.ny - 4)) as f64;
        let expected = 0.5 * gamma * gamma * interior * grid.dx() * grid.dy();
        let p = sol.palinstrophy();
        assert!((p - expected).abs() < 1e-9, "palinstrophy {p} vs analytic {expected}");
        assert!(p > 0.0, "a graded shear has positive palinstrophy");

        // The discriminator vs enstrophy: a *uniform* shear u = γy has constant
        // vorticity ω = −γ, so ∇ω = 0 and the palinstrophy vanishes — even though
        // the enstrophy is large. Palinstrophy sees the gradient, not the magnitude.
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let uniform = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(uniform.palinstrophy() < 1e-9, "uniform shear → ∇ω=0 → 0, got {}", uniform.palinstrophy());
        assert!(uniform.enstrophy() > 0.0, "but the uniform shear still has enstrophy");

        // A grid too small for a doubly-interior difference → 0.
        let tg = Grid::new(4, 4, 1.0, 1.0);
        let tiny = FlowSolution {
            grid: tg,
            u: tg.u_field(),
            v: tg.v_field(),
            pressure: tg.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(tiny.palinstrophy(), 0.0);
    }

    #[test]
    fn max_strain_rate_is_deformation_not_rotation() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // A uniform shear u = γy, v = 0: the strain rate is √(2·½γ²) = γ.
        let gamma = 3.0_f64;
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((shear.max_strain_rate() - gamma).abs() < 1e-9, "shear strain {}", shear.max_strain_rate());

        // Solid-body rotation u=−Ωy, v=Ωx is *pure spin*: zero strain, even though
        // its vorticity is 2Ω. This is the whole point — strain ≠ rotation.
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(rotation.max_strain_rate() < 1e-9, "rotation has no strain: {}", rotation.max_strain_rate());
        assert!(
            (rotation.max_vorticity() - 2.0 * omega).abs() < 1e-9,
            "but it does rotate (ω = 2Ω)"
        );

        // A grid too small for an interior central difference → 0.
        let tg = Grid::new(2, 2, 1.0, 1.0);
        let tiny = FlowSolution {
            grid: tg,
            u: tg.u_field(),
            v: tg.v_field(),
            pressure: tg.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(tiny.max_strain_rate(), 0.0);
    }

    #[test]
    fn mean_strain_rate_is_the_interior_average_deformation() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // (a) SOLID-BODY ROTATION (u = −Ωy, v = Ωx) is pure spin: zero deformation → mean = 0.
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(rotation.mean_strain_rate().abs() < 1e-9, "rotation → 0 strain");

        // (b) PURE SHEAR u(y) = γy, v = 0: |S| = γ uniformly over the interior, so the mean
        // equals γ AND equals the peak max_strain_rate (constant field; threads max).
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((shear.mean_strain_rate() - gamma).abs() <= 1e-9, "uniform shear → ⟨|S|⟩ = γ");
        assert!(
            (shear.mean_strain_rate() - shear.max_strain_rate()).abs()
                <= 1e-9 * shear.max_strain_rate(),
            "constant strain: mean = peak"
        );

        // (c) NON-UNIFORM quadratic shear u(y) = a·y² (|S| = 2a·y varies by row): the strain
        // varies, so the mean lies strictly between 0 and the peak.
        let a = 0.5_f64;
        let mut uq = grid.u_field();
        for j in 0..grid.ny {
            let yc = (j as f64 + 0.5) * dy;
            let val = a * yc * yc;
            for i in 0..=grid.nx {
                uq.set(i, j, val);
            }
        }
        let quad = FlowSolution {
            grid,
            u: uq,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            quad.mean_strain_rate() > 0.0 && quad.mean_strain_rate() < quad.max_strain_rate(),
            "varying strain: 0 < mean < peak"
        );

        // (d) A grid too small for an interior central difference → 0.
        let tg = Grid::new(2, 2, 1.0, 1.0);
        let tiny = FlowSolution {
            grid: tg,
            u: tg.u_field(),
            v: tg.v_field(),
            pressure: tg.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(tiny.mean_strain_rate(), 0.0);
    }

    #[test]
    fn rms_strain_rate_is_the_l2_deformation_rate() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // (a) PURE SHEAR u(y) = γy, v = 0: |S| = γ uniformly → rms = γ. For a pure shear the
        // strain-rate magnitude equals the vorticity magnitude (|S| = |ω| = γ), so
        // rms_strain_rate == rms_vorticity (threads rms_vorticity via the ∇u decomposition).
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((shear.rms_strain_rate() - gamma).abs() <= 1e-9, "uniform shear → rms = γ");
        assert!(
            (shear.rms_strain_rate() - shear.rms_vorticity()).abs() <= 1e-9 * shear.rms_vorticity(),
            "pure shear: |S| = |ω| ⟹ rms_strain_rate = rms_vorticity"
        );

        // (b) SOLID-BODY ROTATION (u = −Ωy, v = Ωx) is pure spin: rms strain = 0, but the rms
        // vorticity is 2Ω ≠ 0 (strain ≠ rotation).
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(rotation.rms_strain_rate().abs() < 1e-9, "rotation → 0 strain");
        assert!(rotation.rms_vorticity() > 1e-9, "but it does rotate");

        // (c) NON-UNIFORM quadratic shear u(y) = a·y²: the rms is ≥ the mean (Jensen) and
        // strictly positive.
        let a = 0.5_f64;
        let mut uq = grid.u_field();
        for j in 0..grid.ny {
            let yc = (j as f64 + 0.5) * dy;
            let val = a * yc * yc;
            for i in 0..=grid.nx {
                uq.set(i, j, val);
            }
        }
        let quad = FlowSolution {
            grid,
            u: uq,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            quad.rms_strain_rate() >= quad.mean_strain_rate() && quad.rms_strain_rate() > 0.0,
            "rms ≥ mean > 0"
        );

        // (d) A grid too small for an interior central difference → 0.
        let tg = Grid::new(2, 2, 1.0, 1.0);
        let tiny = FlowSolution {
            grid: tg,
            u: tg.u_field(),
            v: tg.v_field(),
            pressure: tg.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(tiny.rms_strain_rate(), 0.0);
    }

    #[test]
    fn mean_dissipation_rate_is_viscosity_times_strain_squared() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // (a) WORKED: a pure shear u(y) = γy has |S| = γ uniformly, so with ν = 0.1 the
        // dissipation is ε = ν·γ² = 0.1·4 = 0.4.
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((shear.mean_dissipation_rate(0.1) - 0.4).abs() <= 1e-9, "ε = ν·γ² = 0.4");

        // (b) THREAD rms_vorticity (#406) (non-tautological): for a pure shear |S| = |ω| = γ,
        // so the dissipation also equals ν·rms_vorticity².
        assert!(
            (shear.mean_dissipation_rate(0.1) - 0.1 * shear.rms_vorticity().powi(2)).abs()
                <= 1e-9 * shear.mean_dissipation_rate(0.1),
            "ε = ν·|ω|² for a pure shear"
        );

        // (c) SOLID-BODY ROTATION (u = −Ωy, v = Ωx): pure spin strains nothing, so it
        // dissipates nothing even though it rotates.
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(rotation.mean_dissipation_rate(0.1).abs() < 1e-9, "rotation → no dissipation");

        // (d) Proportional to viscosity, and zero at ν = 0.
        assert!(
            (shear.mean_dissipation_rate(0.2) - 2.0 * shear.mean_dissipation_rate(0.1)).abs()
                <= 1e-9 * shear.mean_dissipation_rate(0.2),
            "ε ∝ ν"
        );
        assert_eq!(shear.mean_dissipation_rate(0.0), 0.0);
    }

    #[test]
    fn kolmogorov_length_scale_is_nu3_over_eps_quarter() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // Pure shear u(y) = γy: |S| = γ uniformly, so ε = ν·γ² and η = (ν³/ε)^(1/4) = √(ν/γ).
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) WORKED (closed form, non-tautological — field-based body vs analytic): with
        // ν = 0.1, η = √(ν/γ) = √(0.1/2) = √0.05 ≈ 0.223607.
        assert!(
            (shear.kolmogorov_length_scale(0.1) - (0.1_f64 / 2.0).sqrt()).abs() <= 1e-9,
            "η = √(ν/γ) for a shear"
        );

        // (b) THREAD mean_dissipation_rate (#431): by definition η⁴·ε = ν³.
        let eps = shear.mean_dissipation_rate(0.1);
        assert!(
            (shear.kolmogorov_length_scale(0.1).powi(4) * eps - 0.1_f64.powi(3)).abs()
                <= 1e-9 * 0.1_f64.powi(3),
            "η⁴·ε = ν³"
        );

        // (c) SCALING with viscosity: for a shear η = √(ν/γ) ∝ √ν, so η(0.4) = 2·η(0.1).
        assert!(
            (shear.kolmogorov_length_scale(0.4) - 2.0 * shear.kolmogorov_length_scale(0.1)).abs()
                <= 1e-9 * shear.kolmogorov_length_scale(0.4),
            "η ∝ √ν"
        );

        // (d) NO DISSIPATION → 0 sentinel: a solid-body rotation strains nothing (ε = 0), so
        // the Kolmogorov scale is undefined.
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(rotation.kolmogorov_length_scale(0.1), 0.0);
    }

    #[test]
    fn kolmogorov_time_scale_is_nu_over_eps_root() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // Pure shear u(y) = γy: |S| = γ uniformly, so ε = ν·γ² and τ_η = (ν/ε)^(1/2) = 1/γ.
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) WORKED (closed form, ν-independent): for a shear ε = ν·γ² so τ_η = 1/γ = 0.5,
        // the same at any viscosity.
        assert!((shear.kolmogorov_time_scale(0.1) - 1.0 / 2.0).abs() <= 1e-9, "τ_η = 1/γ = 0.5");
        assert!(
            (shear.kolmogorov_time_scale(0.4) - shear.kolmogorov_time_scale(0.1)).abs() <= 1e-9,
            "τ_η is ν-independent for a fixed-strain shear"
        );

        // (b) THREAD kolmogorov_length_scale (#434) (non-tautological): the Kolmogorov-scale
        // Reynolds number η²/(τ_η·ν) is exactly 1 by construction.
        let nu = 0.1_f64;
        let eta = shear.kolmogorov_length_scale(nu);
        let tau = shear.kolmogorov_time_scale(nu);
        assert!((eta * eta / (tau * nu) - 1.0).abs() <= 1e-9, "η²/(τ_η·ν) = 1");

        // (c) NO DISSIPATION → 0 sentinel: a solid-body rotation strains nothing (ε = 0).
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(rotation.kolmogorov_time_scale(0.1), 0.0);
    }

    #[test]
    fn kolmogorov_velocity_scale_completes_the_microscale_triad() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // Pure shear u(y) = γy: |S| = γ uniformly, so ε = ν·γ² and u_η = (ν·ε)^(1/4) = √(ν·γ).
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) WORKED (closed form): for a shear u_η = √(ν·γ); with ν = 0.1, √(0.1·2) = √0.2 ≈
        // 0.447214.
        assert!(
            (shear.kolmogorov_velocity_scale(0.1) - (0.1_f64 * 2.0).sqrt()).abs() <= 1e-9,
            "u_η = √(ν·γ) for a shear"
        );

        // (b) THREAD kolmogorov_length_scale (#434) + kolmogorov_time_scale (#437)
        // (non-tautological): the three scales are self-consistent — u_η = η/τ_η.
        let nu = 0.1_f64;
        assert!(
            (shear.kolmogorov_velocity_scale(nu)
                - shear.kolmogorov_length_scale(nu) / shear.kolmogorov_time_scale(nu))
            .abs()
                <= 1e-9 * shear.kolmogorov_velocity_scale(nu),
            "u_η = η/τ_η"
        );

        // (c) UNIT REYNOLDS: by construction the Kolmogorov-scale Reynolds number u_η·η/ν = 1.
        assert!(
            (shear.kolmogorov_velocity_scale(nu) * shear.kolmogorov_length_scale(nu) / nu - 1.0)
                .abs()
                <= 1e-9,
            "u_η·η/ν = 1"
        );

        // (d) NO DISSIPATION → 0 sentinel: a solid-body rotation strains nothing (ε = 0).
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(rotation.kolmogorov_velocity_scale(0.1), 0.0);
    }

    #[test]
    fn taylor_microscale_is_the_intermediate_turbulence_length() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // Pure shear u(y) = γy: |S| = γ uniformly, so ε = ν·γ², and the Taylor
        // microscale λ = u'·√(15ν/ε) = u_rms·√15/γ (independent of ν).
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) WORKED + THREAD rms_speed: for the shear ε = νγ², so λ = u_rms·√15/γ.
        let s = shear.taylor_microscale(0.1);
        assert!(
            (s - shear.rms_speed() * 15.0_f64.sqrt() / gamma).abs() <= 1e-9 * s,
            "λ = u_rms·√15/γ for a pure shear"
        );

        // (b) ν-INVARIANCE (non-obvious): ε ∝ ν for a pure shear, so the ν cancels and
        // λ is identical at very different viscosities.
        assert!(
            (shear.taylor_microscale(0.1) - shear.taylor_microscale(0.9)).abs()
                <= 1e-9 * shear.taylor_microscale(0.1),
            "λ independent of ν for a pure shear"
        );

        // (c) DEFINING IDENTITY (non-tautological round-trip): ε = 15ν·u'²/λ², i.e.
        // λ²·ε = 15ν·u_rms².
        let nu = 0.1_f64;
        let lam = shear.taylor_microscale(nu);
        let rhs = 15.0 * nu * shear.rms_speed().powi(2);
        assert!(
            (lam * lam * shear.mean_dissipation_rate(nu) - rhs).abs() <= 1e-9 * rhs,
            "λ²·ε = 15ν·u'²"
        );

        // (d) NO DISSIPATION → 0 sentinel: a solid-body rotation strains nothing (ε = 0).
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(rotation.taylor_microscale(0.1), 0.0);
    }

    #[test]
    fn taylor_to_kolmogorov_ratio_is_the_inertial_range_width() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // Pure shear u(y) = γy → ε = ν·γ²; λ = u_rms·√15/γ (ν-independent) and η =
        // (ν³/ε)^(1/4) = √(ν/γ), so λ/η = u_rms·√(15/(γν)) ∝ 1/√ν.
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) THREAD delegation: λ/η = taylor_microscale / kolmogorov_length_scale.
        let nu = 0.1_f64;
        let r = shear.taylor_to_kolmogorov_ratio(nu);
        assert!(
            (r - shear.taylor_microscale(nu) / shear.kolmogorov_length_scale(nu)).abs() <= 1e-9 * r,
            "λ/η = λ ÷ η"
        );

        // (b) IDENTITY (non-tautological, threads taylor_reynolds_number): in isotropic
        // turbulence λ/η = 15^(1/4)·√(Re_λ).
        let rhs = 15.0_f64.powf(0.25) * shear.taylor_reynolds_number(nu).sqrt();
        assert!((r - rhs).abs() <= 1e-9 * r, "λ/η = 15^(1/4)·√Re_λ");

        // (c) 1/√ν SCALING: the shear velocity field is fixed, so λ/η ∝ 1/√ν — quartering
        // ν doubles the ratio.
        assert!(
            (shear.taylor_to_kolmogorov_ratio(0.025) - 2.0 * shear.taylor_to_kolmogorov_ratio(0.1))
                .abs()
                <= 1e-9 * shear.taylor_to_kolmogorov_ratio(0.025),
            "λ/η ∝ 1/√ν"
        );

        // (d) NO DISSIPATION → 0: a solid-body rotation strains nothing (η = 0).
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(rotation.taylor_to_kolmogorov_ratio(0.1), 0.0);
    }

    #[test]
    fn integral_to_kolmogorov_ratio_is_the_total_cascade_span() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // Pure shear u(y) = γy → ε = ν·γ²; L = u_rms³/ε, η = (ν³/ε)^(1/4), so L/η ∝
        // 1/ν^(3/2) and L/η = Re_λ^(3/2)/15^(3/4).
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) THREAD delegation: L/η = integral_length_scale / kolmogorov_length_scale.
        let nu = 0.1_f64;
        let r = shear.integral_to_kolmogorov_ratio(nu);
        assert!(
            (r - shear.integral_length_scale(nu) / shear.kolmogorov_length_scale(nu)).abs()
                <= 1e-9 * r,
            "L/η = L ÷ η"
        );

        // (b) IDENTITY (non-tautological, threads taylor_reynolds_number): in isotropic
        // turbulence L/η = Re_λ^(3/2) / 15^(3/4).
        let rhs = shear.taylor_reynolds_number(nu).powf(1.5) / 15.0_f64.powf(0.75);
        assert!((r - rhs).abs() <= 1e-9 * r, "L/η = Re_λ^(3/2)/15^(3/4)");

        // (c) 1/ν^(3/2) SCALING: the shear velocity field is fixed, so L/η ∝ ν^(−3/2) —
        // quartering ν multiplies the ratio by 8 (4^(3/2)).
        assert!(
            (shear.integral_to_kolmogorov_ratio(0.025)
                - 8.0 * shear.integral_to_kolmogorov_ratio(0.1))
            .abs()
                <= 1e-9 * shear.integral_to_kolmogorov_ratio(0.025),
            "L/η ∝ 1/ν^(3/2)"
        );

        // (d) NO DISSIPATION → 0: a solid-body rotation strains nothing (η = 0).
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(rotation.integral_to_kolmogorov_ratio(0.1), 0.0);
    }

    #[test]
    fn integral_to_taylor_ratio_completes_the_scale_ratio_trio() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // Pure shear u(y) = γy → ε = ν·γ²; L = u_rms³/ε, λ = u_rms·√(15ν/ε), so L/λ ∝ 1/ν
        // and L/λ = Re_λ/15 exactly.
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) THREAD delegation: L/λ = integral_length_scale / taylor_microscale.
        let nu = 0.1_f64;
        let r = shear.integral_to_taylor_ratio(nu);
        assert!(
            (r - shear.integral_length_scale(nu) / shear.taylor_microscale(nu)).abs() <= 1e-9 * r,
            "L/λ = L ÷ λ"
        );

        // (b) IDENTITY (non-tautological, threads taylor_reynolds_number): in isotropic
        // turbulence L/λ = Re_λ / 15 exactly.
        let rhs = shear.taylor_reynolds_number(nu) / 15.0;
        assert!((r - rhs).abs() <= 1e-9 * r, "L/λ = Re_λ/15");

        // (c) 1/ν SCALING: the shear velocity field is fixed, so L/λ ∝ 1/ν — halving ν
        // doubles the ratio.
        assert!(
            (shear.integral_to_taylor_ratio(0.05) - 2.0 * shear.integral_to_taylor_ratio(0.1))
                .abs()
                <= 1e-9 * shear.integral_to_taylor_ratio(0.05),
            "L/λ ∝ 1/ν"
        );

        // (d) NO DISSIPATION → 0: a solid-body rotation strains nothing (λ = 0).
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(rotation.integral_to_taylor_ratio(0.1), 0.0);
    }

    #[test]
    fn integral_to_kolmogorov_time_ratio_is_the_time_scale_separation() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // Pure shear u(y) = γy → ε = ν·γ²; T = u_rms²/ε, τ_η = √(ν/ε), so T/τ_η ∝ 1/ν and
        // T/τ_η = Re_λ/√15 exactly.
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) THREAD delegation: T/τ_η = integral_time_scale / kolmogorov_time_scale.
        let nu = 0.1_f64;
        let r = shear.integral_to_kolmogorov_time_ratio(nu);
        assert!(
            (r - shear.integral_time_scale(nu) / shear.kolmogorov_time_scale(nu)).abs() <= 1e-9 * r,
            "T/τ_η = T ÷ τ_η"
        );

        // (b) IDENTITY (non-tautological, threads taylor_reynolds_number): in isotropic
        // turbulence T/τ_η = Re_λ / √15.
        let rhs = shear.taylor_reynolds_number(nu) / 15.0_f64.sqrt();
        assert!((r - rhs).abs() <= 1e-9 * r, "T/τ_η = Re_λ/√15");

        // (c) 1/ν SCALING: the shear velocity field is fixed, so T/τ_η ∝ 1/ν — halving ν
        // doubles the ratio.
        assert!(
            (shear.integral_to_kolmogorov_time_ratio(0.05)
                - 2.0 * shear.integral_to_kolmogorov_time_ratio(0.1))
            .abs()
                <= 1e-9 * shear.integral_to_kolmogorov_time_ratio(0.05),
            "T/τ_η ∝ 1/ν"
        );

        // (d) NO DISSIPATION → 0: a solid-body rotation strains nothing (τ_η = 0).
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(rotation.integral_to_kolmogorov_time_ratio(0.1), 0.0);
    }

    #[test]
    fn taylor_reynolds_number_is_the_turbulence_intensity() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // Pure shear u(y) = γy: ε = ν·γ², λ = u_rms·√15/γ, so Re_λ = u'·λ/ν =
        // u_rms²·√15/(γ·ν).
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) WORKED + THREAD rms_speed (non-tautological — uses rms_speed², not the
        // taylor_microscale path): Re_λ = u_rms²·√15/(γ·ν).
        let re = shear.taylor_reynolds_number(0.1);
        assert!(
            (re - shear.rms_speed().powi(2) * 15.0_f64.sqrt() / (gamma * 0.1)).abs() <= 1e-9 * re,
            "Re_λ = u_rms²·√15/(γν) for a pure shear"
        );

        // (b) DEFINING IDENTITY (thread rms_speed + mean_dissipation_rate): Re_λ² =
        // 15·u'⁴/(ν·ε).
        let nu = 0.1_f64;
        let re = shear.taylor_reynolds_number(nu);
        let rhs = 15.0 * shear.rms_speed().powi(4) / (nu * shear.mean_dissipation_rate(nu));
        assert!((re * re - rhs).abs() <= 1e-9 * rhs, "Re_λ² = 15·u'⁴/(νε)");

        // (c) 1/ν SCALING: the shear velocity field is fixed, so Re_λ ∝ 1/ν — halving ν
        // doubles Re_λ.
        assert!(
            (shear.taylor_reynolds_number(0.05) - 2.0 * shear.taylor_reynolds_number(0.1)).abs()
                <= 1e-9 * shear.taylor_reynolds_number(0.05),
            "Re_λ ∝ 1/ν"
        );

        // (d) NO DISSIPATION → 0: a solid-body rotation strains nothing (λ = 0).
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(rotation.taylor_reynolds_number(0.1), 0.0);
    }

    #[test]
    fn integral_length_scale_tops_the_turbulence_hierarchy() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // Pure shear u(y) = γy: ε = ν·γ², so the integral scale L = u'³/ε = u_rms³/(νγ²).
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) WORKED + THREAD rms_speed (substitutes the analytic ε = νγ² = 0.4):
        // L = u_rms³/(νγ²).
        let l = shear.integral_length_scale(0.1);
        assert!(
            (l - shear.rms_speed().powi(3) / (0.1 * 4.0)).abs() <= 1e-9 * l,
            "L = u_rms³/(νγ²) for a pure shear"
        );

        // (b) THREAD taylor_microscale (non-tautological EXACT identity): since L = u'³/ε
        // and λ² = 15ν·u'²/ε, we have 15ν·L = λ²·u'.
        let nu = 0.1_f64;
        let l = shear.integral_length_scale(nu);
        let rhs = shear.taylor_microscale(nu).powi(2) * shear.rms_speed();
        assert!((15.0 * nu * l - rhs).abs() <= 1e-9 * rhs, "15ν·L = λ²·u'");

        // (c) 1/ν SCALING: the shear velocity field is fixed, so L ∝ 1/ν — halving ν
        // doubles L.
        assert!(
            (shear.integral_length_scale(0.05) - 2.0 * shear.integral_length_scale(0.1)).abs()
                <= 1e-9 * shear.integral_length_scale(0.05),
            "L ∝ 1/ν"
        );

        // (d) NO DISSIPATION → 0: a solid-body rotation strains nothing (ε = 0).
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(rotation.integral_length_scale(0.1), 0.0);
    }

    #[test]
    fn integral_time_scale_caps_the_time_hierarchy() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1
        let (dx, dy) = (grid.dx(), grid.dy());

        // Pure shear u(y) = γy: ε = ν·γ², so the integral time T = u'²/ε = u_rms²/(νγ²).
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) WORKED + THREAD rms_speed (substitutes the analytic ε = νγ² = 0.4):
        // T = u_rms²/(νγ²).
        let t = shear.integral_time_scale(0.1);
        assert!(
            (t - shear.rms_speed().powi(2) / (0.1 * 4.0)).abs() <= 1e-9 * t,
            "T = u_rms²/(νγ²) for a pure shear"
        );

        // (b) THREAD integral_length_scale (non-tautological EXACT identity): T = L/u'.
        let nu = 0.1_f64;
        let t = shear.integral_time_scale(nu);
        assert!(
            (t - shear.integral_length_scale(nu) / shear.rms_speed()).abs() <= 1e-9 * t,
            "T = L/u'"
        );

        // (c) 1/ν SCALING: the shear velocity field is fixed, so T ∝ 1/ν — halving ν
        // doubles T.
        assert!(
            (shear.integral_time_scale(0.05) - 2.0 * shear.integral_time_scale(0.1)).abs()
                <= 1e-9 * shear.integral_time_scale(0.05),
            "T ∝ 1/ν"
        );

        // (d) NO DISSIPATION → 0: a solid-body rotation strains nothing (ε = 0).
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            let val = -omega * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                ur.set(i, j, val);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            let val = omega * (i as f64 + 0.5) * dx;
            for j in 0..=grid.ny {
                vr.set(i, j, val);
            }
        }
        let rotation = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(rotation.integral_time_scale(0.1), 0.0);
    }

    #[test]
    fn channel_flow_conserves_mass_inlet_to_outlet() {
        // The integrated mass flux across the inlet must equal the flux
        // across the outlet — continuity, integrated over the domain.
        let grid = Grid::new(40, 12, 4.0, 1.0);
        let fluid = Fluid::new(1.0, 0.08);
        let bcs = Boundaries::channel_flow(1.0);
        let sol = solve_simple(&grid, &fluid, &bcs, &SimpleControls::default());

        let dy = grid.dy();
        // Inlet flux: sum of u over the west boundary faces × dy.
        let inlet_flux: f64 = (0..grid.ny).map(|j| sol.u.at(0, j) * dy).sum();
        // Outlet flux: sum of u over the east boundary faces × dy.
        let outlet_flux: f64 = (0..grid.ny).map(|j| sol.u.at(grid.nx, j) * dy).sum();
        let rel = (inlet_flux - outlet_flux).abs() / inlet_flux.abs().max(1e-12);
        assert!(
            rel < 0.02,
            "mass must be conserved: inlet {inlet_flux}, outlet {outlet_flux} (rel {rel})"
        );
    }

    #[test]
    fn quiescent_fluid_with_no_driving_stays_at_rest() {
        // With every wall stationary and no inlet, nothing drives the
        // flow — the solution must be the trivial zero-velocity field.
        let grid = Grid::new(16, 16, 1.0, 1.0);
        let fluid = Fluid::new(1.0, 0.1);
        let bcs = Boundaries::lid_driven_cavity(0.0); // lid not moving
        let sol = solve_simple(&grid, &fluid, &bcs, &SimpleControls::default());
        assert!(sol.u.abs_max() < 1e-8, "still fluid should not move (u)");
        assert!(sol.v.abs_max() < 1e-8, "still fluid should not move (v)");
    }

    #[test]
    fn higher_reynolds_number_pushes_the_vortex_centre_downstream() {
        // A physical trend check: as Re rises, the lid-driven cavity's
        // vortex centre migrates toward the geometric centre and the
        // flow becomes more inertia-dominated. We check the simpler
        // robust consequence — a higher-Re run still converges and
        // still shows the reversed return flow.
        let grid = Grid::new(20, 20, 1.0, 1.0);
        let bcs = Boundaries::lid_driven_cavity(1.0);
        for nu in [0.1, 0.02] {
            let fluid = Fluid::new(1.0, nu);
            let controls = SimpleControls {
                max_iterations: 5000,
                relax_u: 0.5,
                relax_p: 0.2,
                ..SimpleControls::default()
            };
            let sol = solve_simple(&grid, &fluid, &bcs, &controls);
            assert!(
                sol.converged,
                "Re = {} run should converge (residual {})",
                1.0 / nu,
                sol.residual
            );
            // The return flow along the floor is always present.
            let u_floor = sol.u_at_cell(grid.nx / 2, 1);
            assert!(u_floor < 0.0, "return flow expected at Re = {}", 1.0 / nu);
        }
    }

    #[test]
    fn solution_helpers_average_staggered_faces() {
        // u_at_cell / v_at_cell average the two bracketing faces.
        let grid = Grid::new(2, 2, 1.0, 1.0);
        let mut u = grid.u_field();
        u.set(0, 0, 2.0);
        u.set(1, 0, 4.0);
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        // Cell (0,0) u = average of faces 0 and 1 = 3.
        assert!((sol.u_at_cell(0, 0) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn flow_rate_integrates_the_inlet_and_uniform_flow_conserves_mass() {
        // A 4×3 grid over a 4×1 channel (dy = 1/3). Every u-face = 2 (plug flow).
        let grid = Grid::new(4, 3, 4.0, 1.0);
        let mut u = grid.u_field(); // (nx+1) × ny = 5 × 3
        for i in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(i, j, 2.0);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        // Q = U·height = 2·1 = 2 per unit depth, identical at inlet and outlet.
        assert!((sol.inlet_flow_rate() - 2.0).abs() < 1e-12, "Q_in {}", sol.inlet_flow_rate());
        assert!((sol.outlet_flow_rate() - 2.0).abs() < 1e-12, "Q_out {}", sol.outlet_flow_rate());
        // Uniform throughflow conserves mass exactly.
        assert!(sol.continuity_error() < 1e-12, "continuity {}", sol.continuity_error());
    }

    #[test]
    fn displacement_thickness_matches_the_couette_profile_and_the_deficit_identity() {
        // A 4×6 grid over a 4×3 channel (dy = 0.5) carrying a linear Couette
        // profile u(y) = U·y/Ly, U = 2 (zero at the bottom wall, U at the top).
        let (u_top, ly) = (2.0_f64, 3.0_f64);
        let grid = Grid::new(4, 6, 4.0, ly);
        let dy = grid.dy();
        let mut u = grid.u_field(); // (nx+1) × ny
        for i in 0..=grid.nx {
            for j in 0..grid.ny {
                let y = (j as f64 + 0.5) * dy; // cell-centre height
                u.set(i, j, u_top * y / ly);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        // δ* = ∫(1 − u/U) dy = ∫(1 − y/Ly) dy = Ly/2 for a Couette profile; the
        // midpoint rule integrates the linear integrand exactly.
        let dstar = sol.displacement_thickness(u_top);
        assert!((dstar - ly / 2.0).abs() < 1e-12, "δ* = Ly/2 = 1.5, got {dstar}");
        // STRONG non-tautological cross-check: the mass-deficit identity
        // δ* = Ly − Q_out/U against the independently-implemented outlet_flow_rate.
        let identity = ly - sol.outlet_flow_rate() / u_top;
        assert!((dstar - identity).abs() < 1e-12, "δ* = Ly − Q_out/U: {dstar} vs {identity}");
        // Bounded by the channel height (0 ≤ δ* ≤ Ly for u ∈ [0, U]).
        assert!((0.0..=ly + 1e-12).contains(&dstar));

        // Plug flow (uniform u = U) carries no deficit → δ* = 0.
        let mut plug_u = grid.u_field();
        for i in 0..=grid.nx {
            for j in 0..grid.ny {
                plug_u.set(i, j, u_top);
            }
        }
        let plug = FlowSolution {
            grid,
            u: plug_u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(plug.displacement_thickness(u_top).abs() < 1e-12, "plug flow → δ* = 0");

        // A non-physical reference speed is undefined → 0.
        assert_eq!(sol.displacement_thickness(0.0), 0.0); // u_ref = 0
        assert_eq!(sol.displacement_thickness(-1.0), 0.0); // u_ref < 0
        assert_eq!(sol.displacement_thickness(f64::NAN), 0.0); // non-finite
    }

    #[test]
    fn momentum_thickness_matches_the_definitions_and_bounds_displacement() {
        // EXACT anchor: a uniform half-speed profile u = U/2 over a 4×8 channel of
        // height Ly = 3, referenced to u_ref = U. The momentum integrand
        // (u/U)(1−u/U) = ½·½ = ¼ is CONSTANT, so the midpoint sum is exact:
        // θ = ¼·Ly and δ* = ½·Ly with no discretisation error.
        let (u_full, ly) = (2.0_f64, 3.0_f64);
        let grid = Grid::new(4, 8, 4.0, ly);
        let mut u = grid.u_field();
        for i in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(i, j, 0.5 * u_full);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let theta = sol.momentum_thickness(u_full);
        assert!((theta - ly / 4.0).abs() < 1e-12, "θ = Ly/4 = 0.75, got {theta}");
        let dstar = sol.displacement_thickness(u_full);
        assert!((dstar - ly / 2.0).abs() < 1e-12, "δ* = Ly/2 = 1.5, got {dstar}");
        // θ ≤ δ* always (δ* − θ = ∫(1−u/U)² dy ≥ 0) — threads displacement_thickness.
        assert!(theta <= dstar + 1e-12, "θ ≤ δ*: {theta} vs {dstar}");

        // Physical Couette profile u(y) = U·y/Ly on a fine grid → θ ≈ Ly/6. The
        // momentum integrand is quadratic in y, so the midpoint rule carries an
        // O(dy²) error; a fine grid drives it well below 1e-3 relative.
        let cgrid = Grid::new(4, 200, 4.0, ly);
        let dy = cgrid.dy();
        let mut cu = cgrid.u_field();
        for i in 0..=cgrid.nx {
            for j in 0..cgrid.ny {
                let y = (j as f64 + 0.5) * dy;
                cu.set(i, j, u_full * y / ly);
            }
        }
        let couette = FlowSolution {
            grid: cgrid,
            u: cu,
            v: cgrid.v_field(),
            pressure: cgrid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let theta_c = couette.momentum_thickness(u_full);
        assert!((theta_c - ly / 6.0).abs() / (ly / 6.0) < 1e-3, "Couette θ ≈ Ly/6, got {theta_c}");
        // Couette δ* = Ly/2 (exact, linear integrand) ⇒ shape factor H = δ*/θ ≈ 3.
        let dstar_c = couette.displacement_thickness(u_full);
        let h = dstar_c / theta_c;
        assert!((h - 3.0).abs() < 5e-3, "Couette H = δ*/θ ≈ 3, got {h}");
        assert!(theta_c <= dstar_c, "θ ≤ δ* for Couette");

        // Plug flow (uniform u = U) carries no momentum deficit → θ = 0.
        let mut plug_u = cgrid.u_field();
        for i in 0..=cgrid.nx {
            for j in 0..cgrid.ny {
                plug_u.set(i, j, u_full);
            }
        }
        let plug = FlowSolution {
            grid: cgrid,
            u: plug_u,
            v: cgrid.v_field(),
            pressure: cgrid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(plug.momentum_thickness(u_full).abs() < 1e-12, "plug flow → θ = 0");

        // A non-physical reference speed is undefined → 0.
        assert_eq!(sol.momentum_thickness(0.0), 0.0); // u_ref = 0
        assert_eq!(sol.momentum_thickness(-1.0), 0.0); // u_ref < 0
        assert_eq!(sol.momentum_thickness(f64::NAN), 0.0); // non-finite
    }

    #[test]
    fn energy_thickness_completes_the_integral_thickness_family() {
        // EXACT anchor: a uniform half-speed profile u = U/2 over a 4×8 channel of height
        // Ly = 3, referenced to u_ref = U. The energy integrand (u/U)(1−(u/U)²) =
        // ½·(1−¼) = 0.375 is CONSTANT, so the midpoint sum is exact: δ_E = 0.375·Ly.
        let (u_full, ly) = (2.0_f64, 3.0_f64);
        let grid = Grid::new(4, 8, 4.0, ly);
        let mut u = grid.u_field();
        for i in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(i, j, 0.5 * u_full);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) WORKED EXACT: δ_E = 0.375·Ly = 1.125.
        let de = sol.energy_thickness(u_full);
        assert!((de - 0.375 * ly).abs() < 1e-12, "δ_E = 0.375·Ly = 1.125, got {de}");

        // (b) THREAD momentum_thickness (non-tautological): for this profile the energy
        // shape factor δ_E/θ = 0.375/0.25 = 1.5.
        assert!(
            (de - 1.5 * sol.momentum_thickness(u_full)).abs() < 1e-12,
            "δ_E = 1.5·θ for the half-speed profile"
        );

        // (c) PLUG FLOW → 0: a uniform u = U carries no kinetic-energy deficit.
        let mut plug_u = grid.u_field();
        for i in 0..=grid.nx {
            for j in 0..grid.ny {
                plug_u.set(i, j, u_full);
            }
        }
        let plug = FlowSolution {
            grid,
            u: plug_u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(plug.energy_thickness(u_full).abs() < 1e-12, "plug flow → δ_E = 0");

        // (d) A non-physical reference speed is undefined → 0.
        assert_eq!(sol.energy_thickness(0.0), 0.0); // u_ref = 0
        assert_eq!(sol.energy_thickness(-1.0), 0.0); // u_ref < 0
        assert_eq!(sol.energy_thickness(f64::NAN), 0.0); // non-finite
    }

    #[test]
    fn energy_shape_factor_is_the_energy_analog_of_the_shape_factor() {
        // EXACT anchor: a uniform half-speed profile u = U/2 over a 4×8 channel of height
        // Ly = 3, referenced to u_ref = U. Both integrands are constant, so δ_E = 0.375·Ly
        // and θ = 0.25·Ly exactly → H* = δ_E/θ = 1.5.
        let (u_full, ly) = (2.0_f64, 3.0_f64);
        let grid = Grid::new(4, 8, 4.0, ly);
        let mut u = grid.u_field();
        for i in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(i, j, 0.5 * u_full);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) WORKED EXACT: H* = δ_E/θ = 0.375/0.25 = 1.5.
        let hstar = sol.energy_shape_factor(u_full);
        assert!((hstar - 1.5).abs() < 1e-12, "H* = 1.5, got {hstar}");

        // (b) THREAD energy_thickness + momentum_thickness (delegation consistency).
        assert!(
            (hstar - sol.energy_thickness(u_full) / sol.momentum_thickness(u_full)).abs() < 1e-12,
            "H* = δ_E/θ"
        );

        // (c) BOUND: H* ≥ 1 always, since δ_E ≥ θ (δ_E − θ = ∫(u/U)²(1−u/U) dy ≥ 0).
        assert!(hstar >= 1.0, "H* ≥ 1, got {hstar}");

        // (d) PLUG FLOW → 0 (θ = 0): a uniform u = U has no momentum deficit.
        let mut plug_u = grid.u_field();
        for i in 0..=grid.nx {
            for j in 0..grid.ny {
                plug_u.set(i, j, u_full);
            }
        }
        let plug = FlowSolution {
            grid,
            u: plug_u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(plug.energy_shape_factor(u_full), 0.0); // plug flow → θ = 0 → 0
        assert_eq!(sol.energy_shape_factor(0.0), 0.0); // u_ref = 0
        assert_eq!(sol.energy_shape_factor(f64::NAN), 0.0); // non-finite
    }

    #[test]
    fn momentum_thickness_reynolds_number_is_the_transition_parameter() {
        // EXACT anchor: a uniform half-speed profile u = U/2 over a 4×8 channel of height
        // Ly = 3, referenced to u_ref = U. The integrand (u/U)(1−u/U) = ¼ is constant, so
        // θ = Ly/4 = 0.75 exactly (no discretisation error).
        let (u_full, ly) = (2.0_f64, 3.0_f64);
        let grid = Grid::new(4, 8, 4.0, ly);
        let mut u = grid.u_field();
        for i in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(i, j, 0.5 * u_full);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // (a) WORKED (exact, analytic θ = Ly/4, independent of the impl): Re_θ =
        // u_ref·θ/ν = 2·0.75/0.1 = 15.
        let nu = 0.1_f64;
        let re_theta = sol.momentum_thickness_reynolds_number(u_full, nu);
        assert!(
            (re_theta - u_full * (ly / 4.0) / nu).abs() <= 1e-9 * re_theta,
            "Re_θ = u_ref·(Ly/4)/ν = 15"
        );

        // (b) THREAD momentum_thickness (delegation consistency): Re_θ = u_ref·θ/ν.
        assert!(
            (re_theta - u_full * sol.momentum_thickness(u_full) / nu).abs() <= 1e-9 * re_theta,
            "Re_θ = u_ref·θ/ν"
        );

        // (c) 1/ν SCALING (θ and u_ref fixed, ν only in the denominator): halve ν →
        // double Re_θ. (Note Re_θ is NOT ∝ u_ref, since θ itself depends on u_ref.)
        assert!(
            (sol.momentum_thickness_reynolds_number(u_full, nu / 2.0) - 2.0 * re_theta).abs()
                <= 1e-9 * 2.0 * re_theta,
            "Re_θ ∝ 1/ν"
        );

        // (d) GUARDS: non-physical viscosity or reference speed → 0.
        assert_eq!(sol.momentum_thickness_reynolds_number(u_full, 0.0), 0.0);
        assert_eq!(sol.momentum_thickness_reynolds_number(u_full, -1.0), 0.0);
        assert_eq!(sol.momentum_thickness_reynolds_number(u_full, f64::NAN), 0.0);
        assert_eq!(sol.momentum_thickness_reynolds_number(0.0, nu), 0.0);
        assert_eq!(sol.momentum_thickness_reynolds_number(f64::NAN, nu), 0.0);
    }

    #[test]
    fn shape_factor_completes_the_boundary_layer_family() {
        // EXACT anchor: a uniform half-speed profile u = U/2 (u_ref = U) → δ* = Ly/2,
        // θ = Ly/4 (both exact, constant integrands) → H = δ*/θ = 2 exactly.
        let (u_full, ly) = (2.0_f64, 3.0_f64);
        let grid = Grid::new(4, 8, 4.0, ly);
        let mut u = grid.u_field();
        for i in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(i, j, 0.5 * u_full);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let h = sol.shape_factor(u_full);
        assert!((h - 2.0).abs() < 1e-12, "uniform-half H = 2 exactly, got {h}");
        assert!(h >= 1.0, "H ≥ 1 (δ* ≥ θ)");
        // Sanity tie H·θ == δ* (catches a ratio inversion), threading #232/#238.
        assert!(
            (h * sol.momentum_thickness(u_full) - sol.displacement_thickness(u_full)).abs() < 1e-12,
            "H·θ = δ*"
        );

        // Physical Couette profile u = U·y/Ly on a fine grid → H ≈ 3 (the linear-
        // profile shape factor); θ's quadratic integrand carries the only error.
        let cgrid = Grid::new(4, 200, 4.0, ly);
        let dy = cgrid.dy();
        let mut cu = cgrid.u_field();
        for i in 0..=cgrid.nx {
            for j in 0..cgrid.ny {
                cu.set(i, j, u_full * (j as f64 + 0.5) * dy / ly);
            }
        }
        let couette = FlowSolution {
            grid: cgrid,
            u: cu,
            v: cgrid.v_field(),
            pressure: cgrid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let hc = couette.shape_factor(u_full);
        assert!((hc - 3.0).abs() < 5e-3, "Couette H ≈ 3, got {hc}");
        assert!(hc >= 1.0, "H ≥ 1 for Couette");

        // Plug flow (θ = 0) and non-physical u_ref → 0 (shape factor undefined).
        let mut plug_u = cgrid.u_field();
        for i in 0..=cgrid.nx {
            for j in 0..cgrid.ny {
                plug_u.set(i, j, u_full);
            }
        }
        let plug = FlowSolution {
            grid: cgrid,
            u: plug_u,
            v: cgrid.v_field(),
            pressure: cgrid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(plug.shape_factor(u_full), 0.0, "plug flow → H = 0 (θ = 0)");
        assert_eq!(sol.shape_factor(0.0), 0.0); // u_ref = 0
        assert_eq!(sol.shape_factor(f64::NAN), 0.0); // non-finite
    }

    #[test]
    fn continuity_error_is_relative_and_safe_when_enclosed() {
        // Inlet U = 1 over height 1 → Q_in = 1; outlet faces 10% higher.
        let grid = Grid::new(3, 4, 3.0, 1.0);
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            u.set(0, j, 1.0); // inlet face
            u.set(grid.nx, j, 1.1); // outlet face, +10%
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((sol.inlet_flow_rate() - 1.0).abs() < 1e-12);
        assert!(
            (sol.continuity_error() - 0.1).abs() < 1e-12,
            "relative err {}",
            sol.continuity_error()
        );
        // An enclosed case (no inflow) → 0, never a divide-by-zero.
        let calm = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(calm.continuity_error(), 0.0);
    }

    #[test]
    fn max_pressure_is_the_stagnation_peak() {
        let grid = Grid::new(5, 4, 5.0, 8.0);

        // Uniform p = 5 → max = mean = 5.
        let mut p = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                p.set(i, j, 5.0);
            }
        }
        let uniform = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: p,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((uniform.max_pressure() - 5.0).abs() <= 1e-9, "uniform → max = 5");
        assert!(
            (uniform.max_pressure() - uniform.mean_pressure()).abs() <= 1e-9,
            "uniform: max = mean"
        );

        // Linear ramp p = i → max = nx−1; threads pressure_range + min_pressure.
        let mut pr = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                pr.set(i, j, i as f64);
            }
        }
        let ramp = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: pr,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((ramp.max_pressure() - (grid.nx - 1) as f64).abs() <= 1e-9, "ramp → max = nx−1");

        // Exact span: max − min = range (cross-checks three independent loops).
        assert!(
            (ramp.max_pressure() - ramp.min_pressure() - ramp.pressure_range()).abs() <= 1e-9,
            "max − min = range"
        );
        // Ordering min ≤ ⟨p⟩ ≤ max.
        assert!(
            ramp.min_pressure() <= ramp.mean_pressure()
                && ramp.mean_pressure() <= ramp.max_pressure(),
            "min ≤ ⟨p⟩ ≤ max"
        );

        // Negative pressures: NEG_INFINITY init, not a 0-ceiling (p = i − 10 → max = (nx−1) − 10).
        let mut pn = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                pn.set(i, j, i as f64 - 10.0);
            }
        }
        let neg = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: pn,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (neg.max_pressure() - ((grid.nx - 1) as f64 - 10.0)).abs() <= 1e-9,
            "negative max = (nx−1) − 10"
        );
    }

    #[test]
    fn pressure_std_dev_is_the_rms_pressure_fluctuation() {
        // Builder seeding a per-cell pressure field (velocity left zero).
        let build = |nx: usize, ny: usize, vals: &[f64]| {
            let grid = Grid::new(nx, ny, nx as f64, ny as f64);
            let mut p = grid.pressure_field();
            for j in 0..grid.ny {
                for i in 0..grid.nx {
                    p.set(i, j, vals[j * grid.nx + i]);
                }
            }
            FlowSolution {
                grid,
                u: grid.u_field(),
                v: grid.v_field(),
                pressure: p,
                iterations: 0,
                residual: 0.0,
                converged: true,
            }
        };

        // (a) WORKED: a 2×1 field {3, 7} → mean 5, deviations ±2 → σ_p = 2.
        let sol = build(2, 1, &[3.0, 7.0]);
        assert!((sol.pressure_std_dev() - 2.0).abs() <= 1e-12 * 2.0, "σ_p = 2");

        // (b) VARIANCE IDENTITY (non-tautological — a different formula): σ_p² = ⟨p²⟩ − ⟨p⟩².
        let mut sum_p2 = 0.0;
        for j in 0..sol.grid.ny {
            for i in 0..sol.grid.nx {
                let p = sol.pressure.at(i, j);
                sum_p2 += p * p;
            }
        }
        let mean_p2 = sum_p2 / (sol.grid.nx * sol.grid.ny) as f64;
        assert!(
            (sol.pressure_std_dev().powi(2) - (mean_p2 - sol.mean_pressure().powi(2))).abs()
                <= 1e-9 * mean_p2,
            "σ_p² = ⟨p²⟩ − ⟨p⟩²"
        );

        // (c) UNIFORM: a constant field has no spread.
        let uniform = build(3, 2, &[5.0; 6]);
        assert_eq!(uniform.pressure_std_dev(), 0.0, "uniform → σ_p = 0");

        // (d) SHIFT-INVARIANCE (gauge-independent): {103,107} = {3,7}+100 → same σ_p.
        let shifted = build(2, 1, &[103.0, 107.0]);
        assert!(
            (shifted.pressure_std_dev() - sol.pressure_std_dev()).abs() <= 1e-12,
            "σ_p unchanged by a constant datum shift"
        );

        // (e) BOUND: non-negative, and ≤ the pressure range (here σ_p = range/2).
        assert!(sol.pressure_std_dev() >= 0.0, "non-negative");
        assert!(sol.pressure_std_dev() <= sol.pressure_range(), "σ_p ≤ range");
    }

    #[test]
    fn rms_pressure_is_the_quadratic_mean_pressure() {
        // Builder seeding a per-cell pressure field (velocity left zero).
        let build = |nx: usize, ny: usize, vals: &[f64]| {
            let grid = Grid::new(nx, ny, nx as f64, ny as f64);
            let mut p = grid.pressure_field();
            for j in 0..grid.ny {
                for i in 0..grid.nx {
                    p.set(i, j, vals[j * grid.nx + i]);
                }
            }
            FlowSolution {
                grid,
                u: grid.u_field(),
                v: grid.v_field(),
                pressure: p,
                iterations: 0,
                residual: 0.0,
                converged: true,
            }
        };

        // (a) WORKED: a 2×1 field {3, 7} → ⟨p²⟩ = (9+49)/2 = 29 → p_rms = √29 ≈ 5.385.
        let sol = build(2, 1, &[3.0, 7.0]);
        assert!(
            (sol.rms_pressure() - 29.0_f64.sqrt()).abs() <= 1e-12 * sol.rms_pressure(),
            "p_rms = √29"
        );

        // (b) VARIANCE IDENTITY (non-tautological): p_rms² = ⟨p⟩² + σ_p² (threads
        // mean_pressure + pressure_std_dev).
        assert!(
            (sol.rms_pressure().powi(2)
                - (sol.mean_pressure().powi(2) + sol.pressure_std_dev().powi(2)))
            .abs()
                <= 1e-9 * sol.rms_pressure().powi(2),
            "p_rms² = ⟨p⟩² + σ_p²"
        );

        // (c) UNIFORM: a constant p = 5 field → p_rms = |mean| = 5.
        let uniform = build(3, 2, &[5.0; 6]);
        assert!((uniform.rms_pressure() - 5.0).abs() <= 1e-12 * 5.0, "uniform → p_rms = |mean|");

        // (d) BOUND: p_rms ≥ |⟨p⟩| and p_rms ≥ σ_p (both follow from p_rms² = ⟨p⟩² + σ_p²).
        assert!(sol.rms_pressure() >= sol.mean_pressure().abs(), "p_rms ≥ |mean|");
        assert!(sol.rms_pressure() >= sol.pressure_std_dev(), "p_rms ≥ σ_p");

        // (e) ZERO-MEAN field {−3, 3}: p_rms = σ_p = 3 (mean = 0), always non-negative.
        let zero_mean = build(2, 1, &[-3.0, 3.0]);
        assert!((zero_mean.rms_pressure() - 3.0).abs() <= 1e-12 * 3.0, "{{-3,3}} → p_rms = 3");
        assert!(zero_mean.rms_pressure() >= 0.0, "non-negative");
    }

    #[test]
    fn min_pressure_is_the_suction_peak() {
        let grid = Grid::new(5, 4, 5.0, 8.0);

        // Uniform p = 5 → min = mean = 5, range = 0.
        let mut p = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                p.set(i, j, 5.0);
            }
        }
        let uniform = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: p,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((uniform.min_pressure() - 5.0).abs() <= 1e-9, "uniform → min = 5");
        assert!(
            (uniform.min_pressure() - uniform.mean_pressure()).abs() <= 1e-9,
            "uniform: min = mean"
        );
        assert!(uniform.pressure_range().abs() <= 1e-9, "uniform: range = 0");

        // Linear ramp p = i → min = 0; threads pressure_range + mean_pressure.
        let mut pr = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                pr.set(i, j, i as f64);
            }
        }
        let ramp = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: pr,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((ramp.min_pressure() - 0.0).abs() <= 1e-9, "ramp → min = 0");

        // min + range = max (= nx−1 for the ramp); ordering min ≤ ⟨p⟩ ≤ max.
        let max_p = ramp.min_pressure() + ramp.pressure_range();
        assert!((max_p - (grid.nx - 1) as f64).abs() <= 1e-9, "min + range = max = nx−1");
        assert!(
            ramp.min_pressure() <= ramp.mean_pressure() && ramp.mean_pressure() <= max_p,
            "min ≤ ⟨p⟩ ≤ max"
        );

        // Negative pressures: INFINITY-init, not a 0-floor (p = i − 10 → min = −10).
        let mut pn = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                pn.set(i, j, i as f64 - 10.0);
            }
        }
        let neg = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: pn,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((neg.min_pressure() - (-10.0)).abs() <= 1e-9, "negative min = −10");
    }

    #[test]
    fn mean_pressure_is_the_area_average_of_the_pressure_field() {
        let grid = Grid::new(5, 4, 5.0, 8.0);

        // Uniform pressure p = 5 → ⟨p⟩ = 5.
        let mut p = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                p.set(i, j, 5.0);
            }
        }
        let uniform = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: p,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((uniform.mean_pressure() - 5.0).abs() <= 1e-9, "uniform → ⟨p⟩ = 5");

        // Linear ramp p = i → ⟨p⟩ = (nx−1)/2 (mean of 0..nx, independent of j).
        let mut pr = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                pr.set(i, j, i as f64);
            }
        }
        let ramp = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: pr,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let expected = (grid.nx - 1) as f64 / 2.0;
        assert!((ramp.mean_pressure() - expected).abs() <= 1e-9 * expected, "ramp → (nx−1)/2");

        // Bounded by the field extremes (the ramp spans 0 .. nx−1).
        assert!(
            0.0 <= ramp.mean_pressure() && ramp.mean_pressure() <= (grid.nx - 1) as f64,
            "min ≤ ⟨p⟩ ≤ max"
        );

        // Linearity: shifting the whole field by c raises the mean by exactly c.
        let mut ps = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                ps.set(i, j, i as f64 + 10.0);
            }
        }
        let shifted = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: ps,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (shifted.mean_pressure() - (ramp.mean_pressure() + 10.0)).abs() <= 1e-9,
            "shift by c → mean + c"
        );
    }

    #[test]
    fn pressure_coefficient_at_cell_follows_bernoulli() {
        let grid = Grid::new(5, 5, 5.0, 5.0);
        let rho = 1.225_f64;
        let u_ref = 10.0_f64;
        let p_ref = 101_325.0_f64;
        let gamma = 0.6_f64;
        let dy = grid.dy();

        // A shear field u = γ·y plus a Bernoulli-consistent static pressure
        // p = p_ref + ½ρ(U_ref² − V²), so Cp should reduce to 1 − (V/U_ref)².
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, gamma * (j as f64 + 0.5) * dy);
            }
        }
        let mut p = grid.pressure_field();
        for j in 0..grid.ny {
            let v = gamma * (j as f64 + 0.5) * dy;
            for i in 0..grid.nx {
                p.set(i, j, p_ref + 0.5 * rho * (u_ref * u_ref - v * v));
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: p,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        let q_ref = 0.5 * rho * u_ref * u_ref;
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let cp = sol.pressure_coefficient_at_cell(i, j, p_ref, rho, u_ref);
                // Threads speed_at_cell via the incompressible Bernoulli relation.
                let speed = sol.speed_at_cell(i, j);
                assert!(
                    (cp - (1.0 - (speed / u_ref).powi(2))).abs() < 1e-9,
                    "Cp = 1 − (V/U_ref)² at ({i},{j})"
                );
                // Definition: Cp = (p − p_ref)/q_ref.
                assert!(
                    (cp - (sol.pressure.at(i, j) - p_ref) / q_ref).abs() < 1e-9,
                    "Cp = (p − p_ref)/(½ρU_ref²)"
                );
            }
        }
        // Non-positive reference dynamic pressure → 0 sentinel.
        assert_eq!(sol.pressure_coefficient_at_cell(2, 2, p_ref, rho, 0.0), 0.0);
        assert_eq!(sol.pressure_coefficient_at_cell(2, 2, p_ref, 0.0, u_ref), 0.0);
    }

    #[test]
    fn total_pressure_at_cell_threads_ke_and_the_range() {
        let grid = Grid::new(5, 5, 5.0, 5.0);
        let rho = 1.225_f64;
        let gamma = 2.0_f64;
        let dy = grid.dy();

        // A shear velocity field u = γ·y plus an arbitrary varying static pressure.
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, gamma * (j as f64 + 0.5) * dy);
            }
        }
        let mut p = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                p.set(i, j, 100.0 + (i as f64) * 3.0);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: p,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // Threads pressure + kinetic_energy_density_at_cell: p₀ = p + ½ρ|u|².
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let expect = sol.pressure.at(i, j) + sol.kinetic_energy_density_at_cell(i, j, rho);
                assert!(
                    (sol.total_pressure_at_cell(i, j, rho) - expect).abs() < 1e-12,
                    "p₀ = p + ½ρ|u|² at ({i},{j})"
                );
            }
        }

        // STRONG cross-check threading total_pressure_range: the spread of the per-cell
        // total pressure reproduces the independently-implemented range reducer.
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let p0 = sol.total_pressure_at_cell(i, j, rho);
                lo = lo.min(p0);
                hi = hi.max(p0);
            }
        }
        assert!(
            (hi - lo - sol.total_pressure_range(rho)).abs() < 1e-9,
            "max − min p₀ = total_pressure_range"
        );

        // Bernoulli: a field whose static pressure compensates the local dynamic head
        // (p = p_ref − ½ρ|u|²) has a spatially uniform total pressure.
        let p_ref = 500.0_f64;
        let mut ub = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                ub.set(fi, j, gamma * (j as f64 + 0.5) * dy);
            }
        }
        let mut pb = grid.pressure_field();
        for j in 0..grid.ny {
            let speed = gamma * (j as f64 + 0.5) * dy; // cell-centre speed of the shear
            for i in 0..grid.nx {
                pb.set(i, j, p_ref - 0.5 * rho * speed * speed);
            }
        }
        let bern = FlowSolution {
            grid,
            u: ub,
            v: grid.v_field(),
            pressure: pb,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                assert!(
                    (bern.total_pressure_at_cell(i, j, rho) - p_ref).abs() < 1e-9,
                    "Bernoulli: uniform total pressure at ({i},{j})"
                );
            }
        }
    }

    #[test]
    fn kinetic_energy_density_at_cell_threads_speed_and_the_mean() {
        let grid = Grid::new(5, 5, 5.0, 5.0);
        let rho = 1.225_f64; // air

        // A non-uniform field: linear shear u = γ·y, v = 0.
        let gamma = 2.0_f64;
        let dy = grid.dy();
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, gamma * (j as f64 + 0.5) * dy);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // Threads speed_at_cell: KE density = ½·ρ·speed² at every cell.
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let s = sol.speed_at_cell(i, j);
                assert!(
                    (sol.kinetic_energy_density_at_cell(i, j, rho) - 0.5 * rho * s * s).abs() < 1e-12,
                    "KE density = ½ρ·speed² at ({i},{j})"
                );
            }
        }

        // STRONG cross-check threading mean_kinetic_energy_density: the cell-average
        // of the per-cell KE density reproduces the independently-implemented mean.
        let mut sum = 0.0;
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                sum += sol.kinetic_energy_density_at_cell(i, j, rho);
            }
        }
        let mean = sum / (grid.nx * grid.ny) as f64;
        assert!(
            (mean - sol.mean_kinetic_energy_density(rho)).abs() / mean < 1e-12,
            "⟨KE density⟩ = mean_kinetic_energy_density"
        );

        // Uniform flow u = U → KE density = ½ρU² (the dynamic pressure) everywhere.
        let u_inf = 4.0_f64;
        let mut uu = grid.u_field();
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                uu.set(i, j, u_inf);
            }
        }
        let uni = FlowSolution {
            grid,
            u: uu,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let q_dyn = 0.5 * rho * u_inf * u_inf;
        assert!(
            (uni.kinetic_energy_density_at_cell(2, 2, rho) - q_dyn).abs() < 1e-12,
            "uniform flow → ½ρU²"
        );

        // Linear in density.
        assert!(
            (uni.kinetic_energy_density_at_cell(2, 2, 2.0 * rho)
                - 2.0 * uni.kinetic_energy_density_at_cell(2, 2, rho))
            .abs()
                < 1e-9,
            "linear in ρ"
        );
    }

    #[test]
    fn speed_coefficient_of_variation_normalises_the_spread() {
        // Uniform field u = 3 → zero dispersion.
        let grid = Grid::new(5, 4, 5.0, 8.0);
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, 3.0);
            }
        }
        let uniform = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(uniform.speed_coefficient_of_variation().abs() < 1e-12, "uniform → CV = 0");

        // Quiescent field (all-zero velocity) → 0, not NaN (mean ≤ 0 guard).
        let quiescent = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(quiescent.speed_coefficient_of_variation(), 0.0, "zero field → 0, no NaN");

        // Non-uniform field (streamwise gradient).
        let mut ug = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                ug.set(fi, j, fi as f64);
            }
        }
        let varied = FlowSolution {
            grid,
            u: ug,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // Non-tautological: CV = σ_direct / mean, with σ from the deviation-sum form.
        let mean = varied.mean_speed();
        let mut acc = 0.0;
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let d = varied.speed_at_cell(i, j) - mean;
                acc += d * d;
            }
        }
        let sigma_direct = (acc / (grid.nx * grid.ny) as f64).sqrt();
        assert!(
            (varied.speed_coefficient_of_variation() - sigma_direct / mean).abs()
                <= 1e-9 * (sigma_direct / mean),
            "CV = σ / ⟨|u|⟩"
        );

        // Recovers the spread: CV · mean = σ; and CV > 0 for a varying field.
        assert!(
            (varied.speed_coefficient_of_variation() * varied.mean_speed()
                - varied.speed_std_dev())
            .abs()
                <= 1e-9 * varied.speed_std_dev(),
            "CV · mean = σ"
        );
        assert!(varied.speed_coefficient_of_variation() > 0.0, "varying field → CV > 0");
    }

    #[test]
    fn speed_std_dev_matches_the_deviation_sum() {
        // Uniform field u = 3 → zero spread.
        let grid = Grid::new(5, 4, 5.0, 8.0);
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, 3.0);
            }
        }
        let uniform = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(uniform.speed_std_dev().abs() < 1e-12, "uniform → σ = 0");
        assert!(uniform.speed_std_dev() <= uniform.rms_speed() + 1e-12, "σ ≤ rms");

        // Non-uniform field (streamwise gradient).
        let mut ug = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                ug.set(fi, j, fi as f64);
            }
        }
        let varied = FlowSolution {
            grid,
            u: ug,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // Non-tautological: the rms²−mean² identity equals the mean-of-squared-deviations
        // form, computed here over the same cells by a different code path.
        let mean = varied.mean_speed();
        let mut acc = 0.0;
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let d = varied.speed_at_cell(i, j) - mean;
                acc += d * d;
            }
        }
        let sigma_direct = (acc / (grid.nx * grid.ny) as f64).sqrt();
        assert!(
            (varied.speed_std_dev() - sigma_direct).abs() <= 1e-9 * sigma_direct,
            "σ = √⟨(u − ⟨u⟩)²⟩"
        );

        // The spread is positive for a varying field and bounded by the RMS speed.
        assert!(varied.speed_std_dev() > 0.0, "varying field → σ > 0");
        assert!(varied.speed_std_dev() <= varied.rms_speed(), "σ ≤ rms");
    }

    #[test]
    fn rms_speed_is_the_kinetic_energy_velocity_scale() {
        // Uniform field u = 3 → rms = mean = max = 3.
        let grid = Grid::new(5, 4, 5.0, 8.0);
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, 3.0);
            }
        }
        let uniform = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((uniform.rms_speed() - 3.0).abs() < 1e-12, "uniform rms = 3");
        assert!((uniform.rms_speed() - uniform.mean_speed()).abs() < 1e-12, "uniform: rms = mean");
        assert!((uniform.rms_speed() - uniform.max_speed()).abs() < 1e-12, "uniform: rms = max");

        // Non-uniform field (streamwise gradient) → rms strictly exceeds the mean.
        let mut ug = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                ug.set(fi, j, fi as f64);
            }
        }
        let varied = FlowSolution {
            grid,
            u: ug,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(varied.rms_speed() > varied.mean_speed(), "RMS > mean for a varying field");
        assert!(varied.rms_speed() <= varied.max_speed() + 1e-12, "rms ≤ max");

        // Threads mean_kinetic_energy_density exactly: ½·ρ·u_rms² = ⟨½ρ|u|²⟩.
        let rho = 1.225;
        for sol in [&uniform, &varied] {
            let from_rms = 0.5 * rho * sol.rms_speed().powi(2);
            let mean_ke = sol.mean_kinetic_energy_density(rho);
            assert!((from_rms - mean_ke).abs() <= 1e-9 * mean_ke, "½ρ·u_rms² = ⟨KE⟩");
        }
    }

    #[test]
    fn total_kinetic_energy_is_the_extensive_energy() {
        // Grid 5×4 over 5×8 → domain area = 5·8 = 40; uniform u = 3; ρ = 1.225.
        let grid = Grid::new(5, 4, 5.0, 8.0);
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, 3.0);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let rho = 1.225;
        let area = 40.0;

        // Worked: total KE = ½ρU²·A = 0.5·1.225·9·40 = 220.5 J.
        assert!(
            (sol.total_kinetic_energy(rho) - 0.5 * rho * 9.0 * area).abs() <= 1e-9 * 220.5,
            "½ρU²·A = 220.5"
        );

        // Threads mean_kinetic_energy_density: total = mean·area.
        assert!(
            (sol.total_kinetic_energy(rho) - sol.mean_kinetic_energy_density(rho) * area).abs()
                <= 1e-12 * sol.total_kinetic_energy(rho),
            "total = ⟨KE⟩·A"
        );

        // Threads rms_speed (#310): total = ½ρ·u_rms²·A.
        assert!(
            (sol.total_kinetic_energy(rho) - 0.5 * rho * sol.rms_speed().powi(2) * area).abs()
                <= 1e-9 * sol.total_kinetic_energy(rho),
            "total = ½ρ·u_rms²·A"
        );

        // Cell-sum cross-check: Σ per-cell KE density · cell area (a different path).
        let mut s = 0.0;
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                s += sol.kinetic_energy_density_at_cell(i, j, rho);
            }
        }
        let from_cells = s * grid.dx() * grid.dy();
        assert!(
            (sol.total_kinetic_energy(rho) - from_cells).abs() <= 1e-9 * sol.total_kinetic_energy(rho),
            "Σ KE_cell·cell_area = total"
        );

        // Linear in density.
        assert!(
            (sol.total_kinetic_energy(2.0 * rho) - 2.0 * sol.total_kinetic_energy(rho)).abs()
                <= 1e-9 * sol.total_kinetic_energy(rho),
            "linear in ρ"
        );
    }

    #[test]
    fn dynamic_pressure_is_the_peak_half_rho_u_squared() {
        // Uniform field u = 3, ρ = 1.225 → q = ½ρU² = 0.5·1.225·9 = 5.5125 Pa.
        let grid = Grid::new(5, 4, 5.0, 8.0);
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, 3.0);
            }
        }
        let uniform = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let rho = 1.225;
        assert!(
            (uniform.dynamic_pressure(rho) - 5.5125).abs() <= 1e-9 * 5.5125,
            "q = ½ρU² = 5.5125"
        );

        // Threads max_speed: q = ½ρ·U_max².
        assert!(
            (uniform.dynamic_pressure(rho) - 0.5 * rho * uniform.max_speed() * uniform.max_speed())
                .abs()
                <= 1e-12 * uniform.dynamic_pressure(rho),
            "q = ½ρ·U_max²"
        );

        // For a uniform field the peak dynamic pressure equals the mean KE density.
        assert!(
            (uniform.dynamic_pressure(rho) - uniform.mean_kinetic_energy_density(rho)).abs()
                <= 1e-9 * uniform.dynamic_pressure(rho),
            "uniform → q = ⟨KE⟩"
        );

        // Linear in density; a 2× faster field gives 4× the q.
        assert!(
            (uniform.dynamic_pressure(2.0 * rho) - 2.0 * uniform.dynamic_pressure(rho)).abs()
                <= 1e-12 * uniform.dynamic_pressure(rho),
            "linear in ρ"
        );
        let mut u6 = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u6.set(fi, j, 6.0);
            }
        }
        let fast = FlowSolution {
            grid,
            u: u6,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (fast.dynamic_pressure(rho) - 4.0 * uniform.dynamic_pressure(rho)).abs()
                <= 1e-9 * fast.dynamic_pressure(rho),
            "quadratic in U"
        );

        // A NON-uniform field (streamwise gradient): q strictly exceeds the mean KE.
        let mut ug = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                ug.set(fi, j, fi as f64);
            }
        }
        let varied = FlowSolution {
            grid,
            u: ug,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            varied.dynamic_pressure(rho) > varied.mean_kinetic_energy_density(rho),
            "non-uniform → q > ⟨KE⟩ (U_max² > ⟨u²⟩)"
        );
    }

    #[test]
    fn euler_number_normalises_the_pressure_variation() {
        // A 4×4 grid with a uniform velocity u = 2 (so U_max = 2) and a linear pressure
        // ramp p = i (so the pressure range Δp = nx − 1 = 3). With ρ = 1.25 the Euler
        // number is Eu = Δp/(ρU²) = 3/(1.25·4) = 0.6.
        let grid = Grid::new(4, 4, 4.0, 4.0);
        let u_speed = 2.0_f64;
        let mut u = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, u_speed);
            }
        }
        let mut p = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                p.set(i, j, i as f64);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: p,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let rho = 1.25_f64;
        let eu = sol.euler_number(rho);
        // Worked value Eu = 0.6, and delegation Eu = Δp/(2q).
        assert!((eu - 0.6).abs() < 1e-9, "Eu = Δp/(ρU²) = 0.6, got {eu}");
        assert!(
            (eu - sol.pressure_range() / (2.0 * sol.dynamic_pressure(rho))).abs() <= 1e-9 * eu,
            "Eu = pressure_range/(2·dynamic_pressure)"
        );
        // NON-TAUTOLOGICAL thread through the independent max_speed: Eu = Δp/(ρ·U_max²).
        let umax = sol.max_speed();
        assert!(
            (eu - sol.pressure_range() / (rho * umax * umax)).abs() <= 1e-9 * eu,
            "Eu = Δp/(ρU_max²)"
        );
        // Scales as 1/ρ: doubling the density halves Eu.
        assert!(
            (sol.euler_number(2.0 * rho) - 0.5 * eu).abs() <= 1e-9 * eu,
            "Eu ∝ 1/ρ"
        );
        // Guards: non-physical density → 0.
        assert_eq!(sol.euler_number(0.0), 0.0);
        assert_eq!(sol.euler_number(-1.0), 0.0);
        // A quiescent field (U_max = 0 ⇒ q = 0) → 0, even though a pressure range exists.
        let mut p_still = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                p_still.set(i, j, i as f64);
            }
        }
        let still = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: p_still,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(still.pressure_range() > 0.0, "still field retains a pressure range");
        assert_eq!(still.euler_number(rho), 0.0, "U_max = 0 ⇒ q = 0 ⇒ Eu = 0");
    }

    #[test]
    fn mean_kinetic_energy_density_matches_dynamic_pressure_for_uniform_flow() {
        // A uniform u = 3, v = 0 field over a 4×2 grid.
        let grid = Grid::new(4, 2, 4.0, 2.0);
        let mut u = grid.u_field(); // 5 × 2
        for i in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(i, j, 3.0);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let rho = 1.2;
        // Every cell has speed 3 → ⟨½ρ|u|²⟩ = ½·1.2·9 = 5.4 Pa, exactly the
        // freestream dynamic pressure ½ρU².
        let ke = sol.mean_kinetic_energy_density(rho);
        assert!((ke - 0.5 * rho * 9.0).abs() < 1e-12, "mean KE density {ke}");
        // Linear in density: doubling ρ doubles the energy.
        assert!(
            (sol.mean_kinetic_energy_density(2.0 * rho) - 2.0 * ke).abs() < 1e-12
        );
    }

    #[test]
    fn total_pressure_range_is_the_stagnation_pressure_variation() {
        let density = 1.2;
        // (a) A stationary fluid (u = v = 0) with a varying pressure field: the
        // dynamic head is zero everywhere, so p0 = p and Δp0 equals the static
        // range — ties total_pressure_range to pressure_range (non-tautological).
        let grid = Grid::new(3, 2, 3.0, 2.0);
        let mut p = grid.pressure_field();
        p.set(0, 0, 100.0);
        p.set(2, 1, -20.0);
        let stationary = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: p,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (stationary.total_pressure_range(density) - stationary.pressure_range()).abs() < 1e-12,
            "stationary fluid: Δp0 = static Δp"
        );

        // (b) A uniform-pressure field with the {5,10,5}-speed pattern (v = 0, the
        // middle cell's two shared u-faces = 10): p0 = p + ½ρ·speed² with p uniform,
        // so Δp0 = ½·ρ·(10² − 5²) = ½·1.2·75 = 45 Pa, while the STATIC range is 0.
        let g2 = Grid::new(3, 1, 3.0, 1.0);
        let mut u2 = g2.u_field();
        u2.set(1, 0, 10.0);
        u2.set(2, 0, 10.0);
        let speedy = FlowSolution {
            grid: g2,
            u: u2,
            v: g2.v_field(),
            pressure: g2.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let expected = 0.5 * density * (10.0 * 10.0 - 5.0 * 5.0);
        let dp0 = speedy.total_pressure_range(density);
        assert!((dp0 - expected).abs() / expected < 1e-9, "Δp0 = ½ρ·75 = {expected}, got {dp0}");
        assert_eq!(speedy.pressure_range(), 0.0, "static range is zero — the loss is all dynamic");

        // (c) A fully uniform field (u = 3, v = 4, p = 0): p0 constant → Δp0 = 0.
        let g3 = Grid::new(2, 2, 2.0, 2.0);
        let mut u3 = g3.u_field();
        for j in 0..g3.ny {
            for i in 0..=g3.nx {
                u3.set(i, j, 3.0);
            }
        }
        let mut v3 = g3.v_field();
        for j in 0..=g3.ny {
            for i in 0..g3.nx {
                v3.set(i, j, 4.0);
            }
        }
        let uniform = FlowSolution {
            grid: g3,
            u: u3,
            v: v3,
            pressure: g3.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(uniform.total_pressure_range(density), 0.0, "uniform field → Δp0 = 0");
    }

    #[test]
    fn mean_total_pressure_is_static_plus_dynamic_head() {
        // Builder seeding BOTH a uniform velocity and a uniform pressure field.
        let build = |nx: usize, ny: usize, lx: f64, ly: f64, uval: f64, pval: f64| {
            let grid = Grid::new(nx, ny, lx, ly);
            let mut u = grid.u_field();
            for fi in 0..=grid.nx {
                for j in 0..grid.ny {
                    u.set(fi, j, uval);
                }
            }
            let mut p = grid.pressure_field();
            for j in 0..grid.ny {
                for i in 0..grid.nx {
                    p.set(i, j, pval);
                }
            }
            FlowSolution {
                grid,
                u,
                v: grid.v_field(),
                pressure: p,
                iterations: 0,
                residual: 0.0,
                converged: true,
            }
        };
        let rho = 1.225;

        // (a) WORKED UNIFORM (independent anchor): u = 3, p = 5, ρ = 1.225. Every
        // cell has speed 3, so ⟨p0⟩ = 5 + ½·1.225·3² = 5 + 5.5125 = 10.5125 Pa,
        // computed DIRECTLY from the worked numbers (not via the two helper fns).
        let sol = build(5, 4, 5.0, 8.0, 3.0, 5.0);
        let expected = 10.5125;
        let p0 = sol.mean_total_pressure(rho);
        assert!(
            (p0 - expected).abs() / expected < 1e-9,
            "⟨p0⟩ = 5 + ½·1.225·9 = {expected}, got {p0}"
        );

        // (b) DECOMPOSITION: ⟨p0⟩ == ⟨p⟩ + ⟨½ρ|u|²⟩ (the static + dynamic split).
        assert!(
            (sol.mean_total_pressure(rho)
                - (sol.mean_pressure() + sol.mean_kinetic_energy_density(rho)))
            .abs()
                < 1e-12,
            "decomposition into mean_pressure + mean_kinetic_energy_density"
        );

        // (c) DIRECT FIELD AVERAGE (independent path): the area-average of the
        // per-cell total_pressure_at_cell must equal mean_total_pressure — a
        // non-tautological cross-check through a different function family.
        let (nx, ny) = (sol.grid.nx, sol.grid.ny);
        let mut sum = 0.0;
        for j in 0..ny {
            for i in 0..nx {
                sum += sol.total_pressure_at_cell(i, j, rho);
            }
        }
        let direct = sum / (nx * ny) as f64;
        assert!(
            (sol.mean_total_pressure(rho) - direct).abs() < 1e-9,
            "⟨p0⟩ = (1/N)·Σ total_pressure_at_cell = {direct}"
        );

        // (d) ZERO-VELOCITY: a quiescent field has no dynamic head ⇒ ⟨p0⟩ = ⟨p⟩.
        let still = build(4, 3, 4.0, 3.0, 0.0, 7.0);
        assert!(
            (still.mean_total_pressure(rho) - still.mean_pressure()).abs() < 1e-12,
            "u = 0 ⇒ ⟨p0⟩ = ⟨p⟩"
        );

        // (e) BOUND + MONOTONICITY: ⟨p0⟩ ≥ ⟨p⟩ always (dynamic head ≥ 0), and a
        // faster field carries a higher ⟨p0⟩ at the same static pressure.
        assert!(sol.mean_total_pressure(rho) >= sol.mean_pressure());
        let fast = build(5, 4, 5.0, 8.0, 6.0, 5.0);
        assert!(
            fast.mean_total_pressure(rho) > sol.mean_total_pressure(rho),
            "u = 6 field ⇒ higher ⟨p0⟩ than u = 3"
        );
    }

    #[test]
    fn bulk_velocity_is_the_flow_rate_per_height() {
        // Uniform plug flow u = 2 over a 4×3 grid of height 1 → U_bulk = Q/H = 2.
        let grid = Grid::new(4, 3, 4.0, 1.0);
        let mut u = grid.u_field();
        for i in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(i, j, 2.0);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((sol.bulk_velocity() - 2.0).abs() < 1e-12, "U_bulk {}", sol.bulk_velocity());
        // By definition U_bulk · H = Q_in.
        assert!((sol.bulk_velocity() * grid.ly - sol.inlet_flow_rate()).abs() < 1e-12);
    }

    #[test]
    fn bottom_wall_shear_rate_recovers_a_linear_shear() {
        // A 4×4 unit grid (dy = 1) with a linear shear u(y) = γ·y, γ = 3. The
        // one-sided wall gradient 2·u_cell(i,0)/dy recovers γ exactly.
        let grid = Grid::new(4, 4, 4.0, 4.0);
        let gamma = 3.0;
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5); // u_at_cell(i,j) = γ·(j+0.5)·dy
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        // ⟨2·u_cell(i,0)/dy⟩ = 2·(γ·0.5)/1 = γ.
        assert!(
            (sol.bottom_wall_shear_rate() - gamma).abs() < 1e-9,
            "wall shear rate {}",
            sol.bottom_wall_shear_rate()
        );
    }

    #[test]
    fn friction_velocity_is_the_wall_turbulence_scale() {
        // A 4×4 unit grid with a linear shear u(y) = γ·y, γ = 3 → bottom-wall shear
        // rate = γ (per the bottom_wall_shear_rate test).
        let grid = Grid::new(4, 4, 4.0, 4.0);
        let gamma = 3.0;
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5);
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let nu = 1.0e-3;
        // u_τ = √(ν·(∂u/∂y)_wall) = √(ν·γ) — the worked value for the imposed shear γ.
        let u_tau = sol.friction_velocity(nu);
        assert!((u_tau - (nu * gamma).sqrt()).abs() < 1e-12, "u_τ = √(νγ), got {u_tau}");
        // The defining relation u_τ² = ν·(wall shear), against the field's own wall
        // shear (independent path through bottom_wall_shear_rate).
        assert!(
            (u_tau * u_tau - nu * sol.bottom_wall_shear_rate()).abs() < 1e-12,
            "u_τ² = ν·τ_w/ρ"
        );
        // It is a velocity that scales as √ν: quadrupling ν doubles u_τ.
        assert!((sol.friction_velocity(4.0 * nu) - 2.0 * u_tau).abs() < 1e-12, "u_τ ∝ √ν");
        // Non-physical viscosity → 0.
        assert_eq!(sol.friction_velocity(0.0), 0.0);
        assert_eq!(sol.friction_velocity(-1.0e-3), 0.0);
        // A quiescent field (no wall shear) → u_τ = 0.
        let still = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(still.friction_velocity(nu), 0.0, "no wall shear → u_τ = 0");
    }

    #[test]
    fn bottom_wall_shear_stress_is_the_dimensional_wall_traction() {
        // A 4×4 unit grid with a linear shear u(y) = γ·y, γ = 3 → bottom-wall shear
        // rate = γ (per the bottom_wall_shear_rate test).
        let grid = Grid::new(4, 4, 4.0, 4.0);
        let gamma = 3.0_f64;
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5);
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        // τ_w = μ·(∂u/∂y)_wall = μ·γ — the worked value for the imposed shear γ.
        let mu = 2.0_f64;
        let tau = sol.bottom_wall_shear_stress(mu);
        assert!((tau - mu * gamma).abs() < 1e-9, "τ_w = μγ, got {tau}");
        // Delegation: τ_w is exactly μ times the field's own wall shear rate.
        assert!(
            (tau - mu * sol.bottom_wall_shear_rate()).abs() <= 1e-9 * tau,
            "τ_w = μ·(wall shear rate)"
        );
        // NON-TAUTOLOGICAL thread through friction_velocity: τ_w = ρ·u_τ² with ν = μ/ρ,
        // since u_τ = √(ν·rate) ⇒ ρ·u_τ² = ρ·(μ/ρ)·rate = μ·rate = τ_w. This routes the
        // stress through the sqrt-and-square of the independent u_τ path.
        let rho = 1.25_f64;
        let nu = mu / rho;
        assert!(
            (tau - rho * sol.friction_velocity(nu).powi(2)).abs() <= 1e-9 * tau,
            "τ_w = ρ·u_τ²"
        );
        // Linear in μ: doubling the viscosity doubles the stress.
        assert!(
            (sol.bottom_wall_shear_stress(4.0) - 2.0 * sol.bottom_wall_shear_stress(2.0)).abs()
                <= 1e-9 * sol.bottom_wall_shear_stress(4.0),
            "τ_w ∝ μ"
        );
        // Non-physical viscosity → 0.
        assert_eq!(sol.bottom_wall_shear_stress(0.0), 0.0);
        assert_eq!(sol.bottom_wall_shear_stress(-1.0), 0.0);
        // A quiescent field (no wall shear) → τ_w = 0.
        let still = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(still.bottom_wall_shear_stress(mu), 0.0, "no wall shear → τ_w = 0");
    }

    #[test]
    fn skin_friction_coefficient_normalises_the_wall_drag() {
        // A 4×4 unit grid with a linear shear u(y) = γ·y, γ = 3 → bottom-wall shear
        // rate = γ and a non-zero bulk velocity, so C_f is well-defined.
        let grid = Grid::new(4, 4, 4.0, 4.0);
        let gamma = 3.0_f64;
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5);
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let mu = 2.0_f64;
        let rho = 1.25_f64;
        let u_ref = sol.bulk_velocity();
        assert!(u_ref > 0.0, "builder must have non-zero bulk velocity");
        let cf = sol.skin_friction_coefficient(mu, rho);
        // Delegation: C_f = τ_w / (½ρU²), against the field's own τ_w and U.
        assert!(
            (cf - sol.bottom_wall_shear_stress(mu) / (0.5 * rho * u_ref * u_ref)).abs() <= 1e-9 * cf,
            "C_f = τ_w/(½ρU²)"
        );
        // NON-TAUTOLOGICAL thread C_f = 2·(u_τ/U)² through the independent
        // friction_velocity (τ_w = ρ·u_τ² ⇒ C_f = ρu_τ²/(½ρU²) = 2u_τ²/U²).
        let nu = mu / rho;
        assert!(
            (cf - 2.0 * sol.friction_velocity(nu).powi(2) / (u_ref * u_ref)).abs() <= 1e-9 * cf,
            "C_f = 2(u_τ/U)²"
        );
        // Linear in μ at fixed ρ and field: doubling μ doubles C_f.
        assert!(
            (sol.skin_friction_coefficient(4.0, rho)
                - 2.0 * sol.skin_friction_coefficient(2.0, rho))
            .abs()
                <= 1e-9 * sol.skin_friction_coefficient(4.0, rho),
            "C_f ∝ μ"
        );
        // Guards: non-physical density → 0.
        assert_eq!(sol.skin_friction_coefficient(mu, 0.0), 0.0);
        assert_eq!(sol.skin_friction_coefficient(mu, -1.0), 0.0);
        // A quiescent field (U = 0 ⇒ no dynamic pressure) → C_f = 0.
        let still = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(still.skin_friction_coefficient(mu, rho), 0.0, "no flow → C_f = 0");
    }

    #[test]
    fn top_wall_shear_rate_recovers_a_linear_shear() {
        // A 4×4 unit grid (dy = 1, ly = 4) with a linear shear that is ZERO at the
        // TOP wall: u(y) = γ·(ly − y), γ = 3 — no-slip at the top, exactly as the
        // bottom test's u = γ·y is no-slip at the bottom. The one-sided top-wall
        // gradient 2·u_cell(i, ny−1)/dy then recovers γ.
        let grid = Grid::new(4, 4, 4.0, 4.0);
        let gamma = 3.0;
        let dy = grid.dy();
        let ly = grid.ly;
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (ly - (j as f64 + 0.5) * dy); // u_at_cell(i,j) = γ·(ly − y_cell)
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        // ⟨2·u_cell(i, ny−1)/dy⟩ = 2·(γ·0.5·dy)/dy = γ.
        assert!(
            (sol.top_wall_shear_rate() - gamma).abs() < 1e-9,
            "top wall shear rate {}",
            sol.top_wall_shear_rate()
        );
        // This profile is no-slip at the TOP only, so the top estimate recovers γ
        // while the bottom (where u ≠ 0) does not — confirming the two read
        // different boundaries.
        assert!(
            sol.top_wall_shear_rate() < sol.bottom_wall_shear_rate(),
            "the top wall is the no-slip one for this profile"
        );
    }

    #[test]
    fn left_wall_shear_rate_measures_the_side_wall_gradient() {
        // A 4×4 unit grid (dx = 1, lx = 4) with a vertical-velocity field linear in x
        // that is ZERO at the LEFT wall (x = 0): v(x) = γ·x, γ = 3 — no-slip at the
        // left, the vertical-wall analogue of the bottom test's u = γ·y. The one-sided
        // left-wall gradient 2·v_cell(0, j)/dx then recovers γ.
        let grid = Grid::new(4, 4, 4.0, 4.0);
        let gamma = 3.0;
        let dx = grid.dx();
        let mut v = grid.v_field();
        for i in 0..grid.nx {
            let val = gamma * (i as f64 + 0.5) * dx; // v_at_cell(i,j) = γ·x_cell
            for j in 0..=grid.ny {
                v.set(i, j, val);
            }
        }
        let sol = FlowSolution {
            grid,
            u: grid.u_field(),
            v,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        // ⟨2·v_cell(0, j)/dx⟩ = 2·(γ·0.5·dx)/dx = γ.
        assert!(
            (sol.left_wall_shear_rate() - gamma).abs() < 1e-9,
            "left wall shear rate {}",
            sol.left_wall_shear_rate()
        );
        // Distinct axis from the horizontal walls: this field has no u, so the
        // top/bottom (∂u/∂y) shear is zero while the left wall (∂v/∂x) reads γ —
        // confirming the new method measures a genuinely different quantity.
        assert!(
            sol.bottom_wall_shear_rate().abs() < 1e-12,
            "a pure-v field has no horizontal-wall shear"
        );
        assert!(sol.top_wall_shear_rate().abs() < 1e-12);
    }

    #[test]
    fn right_wall_shear_rate_measures_the_side_wall_gradient() {
        // A 4×4 unit grid (dx = 1, lx = 4) with a vertical-velocity field linear in x
        // that is ZERO at the RIGHT wall (x = lx): v(x) = γ·(lx − x), γ = 3 — no-slip
        // at the right, mirroring the left test's v = γ·x. The one-sided right-wall
        // gradient 2·v_cell(nx−1, j)/dx then recovers γ.
        let grid = Grid::new(4, 4, 4.0, 4.0);
        let gamma = 3.0;
        let (dx, lx) = (grid.dx(), grid.lx);
        let mut v = grid.v_field();
        for i in 0..grid.nx {
            let val = gamma * (lx - (i as f64 + 0.5) * dx); // v_at_cell(i,j) = γ·(lx − x_cell)
            for j in 0..=grid.ny {
                v.set(i, j, val);
            }
        }
        let sol = FlowSolution {
            grid,
            u: grid.u_field(),
            v,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        // ⟨2·v_cell(nx−1, j)/dx⟩ = 2·(γ·0.5·dx)/dx = γ.
        assert!(
            (sol.right_wall_shear_rate() - gamma).abs() < 1e-9,
            "right wall shear rate {}",
            sol.right_wall_shear_rate()
        );
        // The field is largest near the LEFT wall and zero at the right, so the left
        // wall reads a much larger gradient — confirming the two read DIFFERENT walls
        // (this method the i = nx−1 column, not i = 0).
        assert!(
            sol.left_wall_shear_rate() > sol.right_wall_shear_rate(),
            "asymmetric field: left {} > right {}",
            sol.left_wall_shear_rate(),
            sol.right_wall_shear_rate()
        );
    }

    #[test]
    fn circulation_integrates_vorticity_over_the_interior() {
        // A 5×5 unit grid (dx=dy=1), v=0, linear shear u(y)=γ·y (γ=2) → ω = −γ.
        // Γ = ∫ω dA = −γ·(interior area) = −2·9 = −18, carrying the sign of ω.
        let grid = Grid::new(5, 5, 5.0, 5.0);
        let gamma = 2.0;
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5); // u_at_cell(i,j) = γ·(j+0.5)
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((sol.circulation() - (-18.0)).abs() < 1e-9, "Γ {}", sol.circulation());
        // A grid too small for an interior difference → 0.
        let tg = Grid::new(2, 2, 1.0, 1.0);
        let tiny = FlowSolution {
            grid: tg,
            u: tg.u_field(),
            v: tg.v_field(),
            pressure: tg.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(tiny.circulation(), 0.0);
    }

    #[test]
    fn vorticity_at_cell_integrates_to_circulation_and_enstrophy() {
        // Solid-body rotation u = −Ω·y, v = Ω·x → ω = ∂v/∂x − ∂u/∂y = 2Ω everywhere.
        let omega = 2.0_f64;
        let grid = Grid::new(6, 6, 6.0, 6.0);
        let (dx, dy) = (grid.dx(), grid.dy());
        let mut u = grid.u_field(); // (nx+1) × ny
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                u.set(fi, j, -omega * (j as f64 + 0.5) * dy); // u_at_cell(i,j) = −Ω·y_j
            }
        }
        let mut v = grid.v_field(); // nx × (ny+1)
        for i in 0..grid.nx {
            for fj in 0..=grid.ny {
                v.set(i, fj, omega * (i as f64 + 0.5) * dx); // v_at_cell(i,j) = Ω·x_i
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        // Each interior cell sees ω = 2Ω = 4.
        assert!((sol.vorticity_at_cell(3, 3) - 2.0 * omega).abs() < 1e-9, "interior ω = 2Ω");
        // Boundary cells (centred stencil unavailable) → 0.
        assert_eq!(sol.vorticity_at_cell(0, 3), 0.0);
        assert_eq!(sol.vorticity_at_cell(3, 0), 0.0);
        assert_eq!(sol.vorticity_at_cell(grid.nx - 1, 3), 0.0);
        // STRONG cross-checks: integrating the per-cell field reproduces the
        // independently-implemented circulation (∫ω dA) and enstrophy (½∫ω² dA).
        let mut gamma = 0.0;
        let mut ens = 0.0;
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let w = sol.vorticity_at_cell(i, j);
                gamma += w * dx * dy;
                ens += w * w;
            }
        }
        ens *= 0.5 * dx * dy;
        assert!((gamma - sol.circulation()).abs() / gamma.abs() < 1e-9, "Σω·dA = circulation");
        assert!((ens - sol.enstrophy()).abs() / ens < 1e-9, "½Σω²·dA = enstrophy");
        // Solid-body rotation: Γ = 2Ω · interior area.
        let interior_area = (grid.nx - 2) as f64 * (grid.ny - 2) as f64 * dx * dy;
        assert!((gamma - 2.0 * omega * interior_area).abs() / gamma < 1e-9, "Γ = 2Ω·A_interior");
    }

    #[test]
    fn strain_rate_at_cell_completes_the_velocity_gradient_decomposition() {
        // 6×6 unit grid (dx = dy = 1).
        let grid = Grid::new(6, 6, 6.0, 6.0);
        let (dx, dy) = (grid.dx(), grid.dy());

        // (1) Solid-body rotation u = −Ω·y, v = +Ω·x — pure spin → strain rate 0.
        let rate = 2.0_f64;
        let mut ur = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                ur.set(fi, j, -rate * (j as f64 + 0.5) * dy);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            for fj in 0..=grid.ny {
                vr.set(i, fj, rate * (i as f64 + 0.5) * dx);
            }
        }
        let rot = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(rot.strain_rate_at_cell(3, 3).abs() < 1e-9, "solid-body rotation has zero strain");

        // (2) Pure straining u = a·x, v = −a·y — strain rate 2a, zero spin.
        let a = 1.5_f64;
        let mut ue = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                ue.set(fi, j, a * (fi as f64) * dx);
            }
        }
        let mut ve = grid.v_field();
        for i in 0..grid.nx {
            for fj in 0..=grid.ny {
                ve.set(i, fj, -a * (fj as f64) * dy);
            }
        }
        let strain = FlowSolution {
            grid,
            u: ue,
            v: ve,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((strain.strain_rate_at_cell(3, 3) - 2.0 * a).abs() < 1e-9, "pure strain rate = 2a");
        // STRONG cross-check threading the independent max_strain_rate reducer.
        assert!((strain.max_strain_rate() - 2.0 * a).abs() < 1e-9, "max strain rate = 2a");

        // (3) Pure shear u = γ·y, v = 0 — strain rate |γ| (equals the vorticity
        // magnitude, which is exactly why a shear layer has Q = 0).
        let gamma = 3.0_f64;
        let mut us = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                us.set(fi, j, gamma * (j as f64 + 0.5) * dy);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((shear.strain_rate_at_cell(3, 3) - gamma).abs() < 1e-9, "pure shear strain = γ");

        // STRONG three-way identity tying all per-cell accessors together: the
        // Q-criterion is Q = (ω² − strain²)/4 (since Q = ½(‖Ω‖² − ‖S‖²), ‖Ω‖² = ½ω²,
        // ‖S‖² = strain²/2). Holds on every field — rotation, strain, shear.
        for sol in [&rot, &strain, &shear] {
            let q = sol.q_criterion_at_cell(3, 3);
            let w = sol.vorticity_at_cell(3, 3);
            let s = sol.strain_rate_at_cell(3, 3);
            assert!((q - (w * w - s * s) / 4.0).abs() < 1e-9, "Q = (ω² − strain²)/4");
        }

        // Boundary cells have no centred stencil → 0.
        assert_eq!(rot.strain_rate_at_cell(0, 3), 0.0);
    }

    #[test]
    fn q_criterion_at_cell_identifies_rotation_over_strain() {
        // 6×6 unit grid (dx = dy = 1).
        let grid = Grid::new(6, 6, 6.0, 6.0);
        let (dx, dy) = (grid.dx(), grid.dy());

        // (1) Solid-body rotation u = −Ω·y, v = +Ω·x — pure rotation (strain S = 0):
        // Q = Ω² > 0 everywhere, and Q = ω²/4 (ω = vorticity_at_cell = 2Ω).
        let rate = 2.0_f64;
        let mut ur = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                ur.set(fi, j, -rate * (j as f64 + 0.5) * dy);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            for fj in 0..=grid.ny {
                vr.set(i, fj, rate * (i as f64 + 0.5) * dx);
            }
        }
        let rot = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let q = rot.q_criterion_at_cell(3, 3);
        assert!((q - rate * rate).abs() < 1e-9, "solid-body rotation Q = Ω², got {q}");
        let w = rot.vorticity_at_cell(3, 3);
        assert!((q - w * w / 4.0).abs() < 1e-9, "pure rotation Q = ω²/4");
        assert!(q > 0.0, "rotation-dominated → Q > 0");
        // Boundary cells have no centred stencil → 0.
        assert_eq!(rot.q_criterion_at_cell(0, 3), 0.0);
        assert_eq!(rot.q_criterion_at_cell(3, grid.ny - 1), 0.0);
        // STRONG cross-check threading the independent max_q_criterion reducer.
        assert!((rot.max_q_criterion() - rate * rate).abs() < 1e-9, "max Q = Ω²");

        // (2) Pure shear u = γ·y, v = 0 — equal rotation and strain → Q = 0, even
        // though the vorticity is non-zero (a shear layer is not a vortex).
        let gamma = 3.0_f64;
        let mut us = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                us.set(fi, j, gamma * (j as f64 + 0.5) * dy);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(shear.q_criterion_at_cell(3, 3).abs() < 1e-9, "pure shear → Q = 0");
        assert!(shear.vorticity_at_cell(3, 3).abs() > 0.1, "...yet shear has vorticity");

        // (3) Pure straining u = a·x, v = −a·y (incompressible, irrotational) → Q = −a²
        // < 0 (strain-dominated, not a vortex), and max_q_criterion clamps to 0.
        let a = 1.5_f64;
        let mut ue = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                ue.set(fi, j, a * (fi as f64) * dx);
            }
        }
        let mut ve = grid.v_field();
        for i in 0..grid.nx {
            for fj in 0..=grid.ny {
                ve.set(i, fj, -a * (fj as f64) * dy);
            }
        }
        let strain = FlowSolution {
            grid,
            u: ue,
            v: ve,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let qs = strain.q_criterion_at_cell(3, 3);
        assert!((qs + a * a).abs() < 1e-9, "pure strain Q = −a², got {qs}");
        assert!(strain.vorticity_at_cell(3, 3).abs() < 1e-9, "pure strain is irrotational");
        assert_eq!(strain.max_q_criterion(), 0.0, "no vortex → max Q clamps to 0");
    }

    #[test]
    fn divergence_at_cell_is_the_per_cell_continuity_residual() {
        // 5×5 unit grid (dx = dy = 1).
        let grid = Grid::new(5, 5, 5.0, 5.0);
        let (dx, dy) = (grid.dx(), grid.dy());

        // Incompressible solid-body rotation u = −Ω·y, v = +Ω·x: ∂u/∂x = ∂v/∂y = 0, so
        // it is divergence-free in EVERY cell — yet it genuinely spins
        // (vorticity_at_cell = 2Ω). Divergence (∇·u) and curl (ω) are independent.
        let omega = 2.0_f64;
        let mut ur = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                ur.set(fi, j, -omega * (j as f64 + 0.5) * dy);
            }
        }
        let mut vr = grid.v_field();
        for i in 0..grid.nx {
            for fj in 0..=grid.ny {
                vr.set(i, fj, omega * (i as f64 + 0.5) * dx);
            }
        }
        let rot = FlowSolution {
            grid,
            u: ur,
            v: vr,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                assert!(
                    rot.divergence_at_cell(i, j).abs() < 1e-12,
                    "solid-body rotation is divergence-free at ({i},{j})"
                );
            }
        }
        assert!(
            (rot.vorticity_at_cell(2, 2) - 2.0 * omega).abs() < 1e-9,
            "...but it carries vorticity 2Ω (curl is not divergence)"
        );

        // Radial source u = a·x, v = a·y: ∂u/∂x = ∂v/∂y = a → ∇·u = 2a in every cell,
        // and it is irrotational (vorticity 0) — the exact dual of the rotation case.
        let a = 3.0;
        let mut us = grid.u_field();
        for fi in 0..=grid.nx {
            for j in 0..grid.ny {
                us.set(fi, j, a * (fi as f64) * dx); // u-face at x = fi·dx
            }
        }
        let mut vs = grid.v_field();
        for i in 0..grid.nx {
            for fj in 0..=grid.ny {
                vs.set(i, fj, a * (fj as f64) * dy); // v-face at y = fj·dy
            }
        }
        let src = FlowSolution {
            grid,
            u: us,
            v: vs,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (src.divergence_at_cell(2, 2) - 2.0 * a).abs() < 1e-9,
            "radial source ∇·u = 2a"
        );
        assert!(
            src.vorticity_at_cell(2, 2).abs() < 1e-9,
            "radial source is irrotational"
        );

        // STRONG cross-check (independent implementation): the peak |∇·u| over the
        // cells reproduces max_divergence, which forms the MAC divergence on its own.
        let mut worst = 0.0_f64;
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                worst = worst.max(src.divergence_at_cell(i, j).abs());
            }
        }
        assert!(
            (worst - src.max_divergence()).abs() < 1e-12,
            "max|∇·u| (per-cell accessor) = max_divergence (reducer)"
        );
        assert!(
            (src.max_divergence() - 2.0 * a).abs() < 1e-9,
            "peak source divergence = 2a"
        );

        // Out-of-range indices return 0 (the accessor is defined on valid cells only).
        assert_eq!(src.divergence_at_cell(grid.nx, 0), 0.0);
        assert_eq!(src.divergence_at_cell(0, grid.ny), 0.0);
    }

    #[test]
    fn max_divergence_is_the_peak_local_continuity_residual() {
        // A 5×5 unit grid (dx = dy = 1).
        let grid = Grid::new(5, 5, 5.0, 5.0);

        // Uniform flow u = U, v = 0 is divergence-free → max|∇·u| = 0.
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                u.set(i, j, 4.0);
            }
        }
        let uniform = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            uniform.max_divergence() < 1e-12,
            "uniform flow is divergence-free: {}",
            uniform.max_divergence()
        );

        // A pure horizontal shear u(y) = γ·y, v = 0 is also divergence-free
        // (∂u/∂x = 0) — divergence is not vorticity.
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = 2.0 * (j as f64 + 0.5);
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            shear.max_divergence() < 1e-12,
            "shear is divergence-free: {}",
            shear.max_divergence()
        );

        // A pure x-expansion u(x) = a·x, v = 0 → ∂u/∂x = a in every cell, so the
        // peak local residual is exactly a.
        let a = 3.0;
        let mut ue = grid.u_field();
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                ue.set(i, j, a * (i as f64) * grid.dx()); // u-face at x = i·dx
            }
        }
        let expanding = FlowSolution {
            grid,
            u: ue,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (expanding.max_divergence() - a).abs() < 1e-9,
            "∇·u = a = {a}, got {}",
            expanding.max_divergence()
        );

        // |∇·u| is sign-blind: a contraction reads the same magnitude.
        let mut uc = grid.u_field();
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                uc.set(i, j, -a * (i as f64) * grid.dx());
            }
        }
        let contracting = FlowSolution {
            grid,
            u: uc,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (contracting.max_divergence() - a).abs() < 1e-9,
            "|∇·u| is sign-blind"
        );
    }

    #[test]
    fn mean_divergence_is_the_area_averaged_continuity_residual() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1

        // (a) SOLENOIDAL: a pure horizontal shear u(y) = γ·y, v = 0 is divergence-free
        // (∂u/∂x = 0 along each row, ∂v/∂y = 0), so the mean residual is 0.
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * grid.dy();
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(shear.mean_divergence().abs() < 1e-9, "solenoidal → ⟨∇·u⟩ = 0");

        // (b) UNIFORM EXPANSION: u(x) = c·x, v = 0 → ∂u/∂x = c in every cell, so ⟨∇·u⟩ = c;
        // and (c) for this CONSTANT-divergence field the mean equals the peak max_divergence
        // (threaded, non-tautological).
        let c = 0.3_f64;
        let mut ue = grid.u_field();
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                ue.set(i, j, c * (i as f64) * grid.dx());
            }
        }
        let expanding = FlowSolution {
            grid,
            u: ue,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((expanding.mean_divergence() - c).abs() <= 1e-9 * c, "⟨∇·u⟩ = c");
        assert!(
            (expanding.mean_divergence() - expanding.max_divergence()).abs()
                <= 1e-9 * expanding.max_divergence(),
            "constant divergence: mean = peak"
        );

        // (d) ZERO-FLOW: a quiescent field has zero divergence.
        let zero = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(zero.mean_divergence(), 0.0);
    }

    #[test]
    fn rms_divergence_is_the_l2_continuity_residual() {
        let grid = Grid::new(5, 5, 5.0, 5.0); // dx = dy = 1

        // (a) SOLENOIDAL: a pure shear u(y) = γ·y, v = 0 is divergence-free → rms = 0.
        let gamma = 2.0_f64;
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * grid.dy();
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(shear.rms_divergence().abs() < 1e-9, "solenoidal → rms = 0");

        // (b) UNIFORM EXPANSION: u(x) = c·x, v = 0 → ∇·u = c everywhere → rms = c; and for
        // this CONSTANT-divergence field rms == mean == max (threads both, non-tautological).
        let c = 0.3_f64;
        let mut ue = grid.u_field();
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                ue.set(i, j, c * (i as f64) * grid.dx());
            }
        }
        let expanding = FlowSolution {
            grid,
            u: ue,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((expanding.rms_divergence() - c).abs() <= 1e-9 * c, "rms = c");
        assert!(
            (expanding.rms_divergence() - expanding.mean_divergence()).abs() <= 1e-9 * c,
            "constant field: rms = mean"
        );
        assert!(
            (expanding.rms_divergence() - expanding.max_divergence()).abs() <= 1e-9 * c,
            "constant field: rms = peak"
        );

        // (c) SIGN-VARYING "tent": u rises then falls in i (faces c·[0,1,2,3,2,1]) so per row
        // ∇·u = [c,c,c,−c,−c]. The signed mean cancels to c/5 but the rms = c — the RMS catches
        // the local violations the mean misses.
        let mut ut = grid.u_field();
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                let f = if i <= 3 { i } else { 6 - i };
                ut.set(i, j, c * (f as f64) * grid.dx());
            }
        }
        let tent = FlowSolution {
            grid,
            u: ut,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!((tent.rms_divergence() - c).abs() <= 1e-9 * c, "tent rms = c");
        assert!(
            (tent.mean_divergence() - c / 5.0).abs() <= 1e-9 * (c / 5.0),
            "tent mean = c/5 (cancellation)"
        );
        assert!(
            tent.rms_divergence() > tent.mean_divergence().abs(),
            "rms > |mean| when divergence changes sign"
        );

        // (d) ZERO-FLOW: a quiescent field has zero divergence.
        let zero = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(zero.rms_divergence(), 0.0);
    }

    #[test]
    fn vorticity_thickness_measures_the_shear_layer() {
        // A 5×5 unit grid (dy = 1, ly = 5) with a linear shear u(y) = γ·y, γ = 2, v = 0.
        let grid = Grid::new(5, 5, 5.0, 5.0);
        let gamma = 2.0;
        let (dy, ly) = (grid.dy(), grid.ly);
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy; // u_at_cell(i,j) = γ·y_cell
            for i in 0..=grid.nx {
                u.set(i, j, val);
            }
        }
        let sol = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        // The cell-centre u spans γ·(ly−dy) and the central-difference ∂u/∂y is
        // exactly γ, so δ_ω = γ(ly−dy)/γ = ly − dy (the discrete full-domain shear
        // resolves the span between the outermost cell centres).
        assert!(
            (sol.vorticity_thickness() - (ly - dy)).abs() < 1e-9,
            "δ_ω = ly−dy = {}, got {}",
            ly - dy,
            sol.vorticity_thickness()
        );
        // δ_ω is a LENGTH — independent of the shear magnitude: doubling γ scales the
        // velocity span and the gradient equally, leaving the thickness unchanged.
        let mut u2 = grid.u_field();
        for j in 0..grid.ny {
            let val = 2.0 * gamma * (j as f64 + 0.5) * dy;
            for i in 0..=grid.nx {
                u2.set(i, j, val);
            }
        }
        let sol2 = FlowSolution {
            grid,
            u: u2,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (sol2.vorticity_thickness() - sol.vorticity_thickness()).abs() < 1e-9,
            "δ_ω is independent of the shear rate γ"
        );
        // A uniform flow has no shear layer → 0.
        let mut uc = grid.u_field();
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                uc.set(i, j, 3.0);
            }
        }
        let uniform = FlowSolution {
            grid,
            u: uc,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(uniform.vorticity_thickness(), 0.0, "uniform flow → no shear layer");
    }

    #[test]
    fn reverse_flow_fraction_measures_the_recirculating_share() {
        let grid = Grid::new(4, 4, 4.0, 4.0);

        // A uniform forward flow u = +U has no reverse flow → fraction 0.
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                u.set(i, j, 2.0);
            }
        }
        let forward = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(forward.reverse_flow_fraction(), 0.0);

        // Every cell running backward (u = −U) → fraction 1.
        let mut ur = grid.u_field();
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                ur.set(i, j, -2.0);
            }
        }
        let reverse = FlowSolution {
            grid,
            u: ur,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(reverse.reverse_flow_fraction(), 1.0);

        // Half the rows flow back (j=0,1 → u<0) and half forward (j=2,3 → u>0):
        // 8 of the 16 cells are reversed → exactly 0.5.
        let mut uh = grid.u_field();
        for j in 0..grid.ny {
            let val = if j < 2 { -1.0 } else { 1.0 };
            for i in 0..=grid.nx {
                uh.set(i, j, val);
            }
        }
        let split = FlowSolution {
            grid,
            u: uh,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (split.reverse_flow_fraction() - 0.5).abs() < 1e-12,
            "half-and-half → 0.5, got {}",
            split.reverse_flow_fraction()
        );
    }

    #[test]
    fn stream_function_range_recovers_the_throughflow_and_vortex_strength() {
        let grid = Grid::new(5, 4, 5.0, 8.0); // dx = 1, dy = 2, Lx = 5, Ly = 8
        let (lx, ly) = (5.0_f64, 8.0_f64);

        // (1) Uniform forward flow u = U, v = 0: horizontal streamlines, the span
        // is the throughflow U·Ly.
        let big_u = 3.0;
        let mut u = grid.u_field();
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                u.set(i, j, big_u);
            }
        }
        let uflow = FlowSolution {
            grid,
            u,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (uflow.stream_function_range() - big_u * ly).abs() < 1e-9,
            "u-flow span = U·Ly, got {}",
            uflow.stream_function_range()
        );

        // (2) Uniform upward flow u = 0, v = V: span = V·Lx. This exercises the
        // v / base-row term — a sign error there passes (1) but fails here.
        let big_v = 2.0;
        let mut v = grid.v_field();
        for j in 0..=grid.ny {
            for i in 0..grid.nx {
                v.set(i, j, big_v);
            }
        }
        let vflow = FlowSolution {
            grid,
            u: grid.u_field(),
            v,
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (vflow.stream_function_range() - big_v * lx).abs() < 1e-9,
            "v-flow span = V·Lx, got {}",
            vflow.stream_function_range()
        );

        // (3) Horizontal shear u(y) = γ·y, v = 0: span = ∫₀^Ly γy dy = γ·Ly²/2.
        let gamma = 0.5;
        let dy = grid.dy();
        let mut us = grid.u_field();
        for j in 0..grid.ny {
            let val = gamma * (j as f64 + 0.5) * dy; // u-face value = γ·y_cell
            for i in 0..=grid.nx {
                us.set(i, j, val);
            }
        }
        let shear = FlowSolution {
            grid,
            u: us,
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert!(
            (shear.stream_function_range() - gamma * ly * ly / 2.0).abs() < 1e-9,
            "shear span = γ·Ly²/2, got {}",
            shear.stream_function_range()
        );

        // (4) A quiescent field has no spanning streamlines → 0.
        let calm = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: grid.pressure_field(),
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        assert_eq!(calm.stream_function_range(), 0.0);
    }

    #[test]
    fn min_pressure_location_finds_the_suction_peak() {
        let grid = Grid::new(5, 4, 5.0, 8.0); // dx = 1, dy = 2
        let (dx, dy) = (grid.dx(), grid.dy());

        // Plant the pressure minimum at cell (3, 1): a uniform 100 Pa field with a dip.
        let mut p = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                p.set(i, j, 100.0);
            }
        }
        p.set(3, 1, -50.0);
        let sol = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: p,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let (x, y) = sol.min_pressure_location().expect("non-empty grid");
        assert!(
            (x - (3.0 + 0.5) * dx).abs() < 1e-9 && (y - (1.0 + 0.5) * dy).abs() < 1e-9,
            "suction peak at cell (3,1) centre, got ({x}, {y})"
        );

        // Move the minimum to cell (0, 3) → the located point follows.
        let mut p2 = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                p2.set(i, j, 10.0);
            }
        }
        p2.set(0, 3, -1.0);
        let sol2 = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: p2,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };
        let (x2, y2) = sol2.min_pressure_location().unwrap();
        assert!(
            (x2 - 0.5 * dx).abs() < 1e-9 && (y2 - (3.0 + 0.5) * dy).abs() < 1e-9,
            "suction peak at cell (0,3) centre, got ({x2}, {y2})"
        );
    }

    #[test]
    fn max_pressure_location_finds_the_stagnation_point() {
        let grid = Grid::new(5, 4, 5.0, 8.0); // dx = 1, dy = 2
        let (dx, dy) = (grid.dx(), grid.dy());

        // Uniform 100 Pa field with a unique spike at cell (3, 1) — the stagnation
        // point — and a distinct dip at (0, 3) so the max and min sit apart.
        let mut p = grid.pressure_field();
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                p.set(i, j, 100.0);
            }
        }
        p.set(3, 1, 250.0); // the stagnation spike (unique maximum)
        p.set(0, 3, -50.0); // a distinct minimum elsewhere
        let sol = FlowSolution {
            grid,
            u: grid.u_field(),
            v: grid.v_field(),
            pressure: p,
            iterations: 0,
            residual: 0.0,
            converged: true,
        };

        // The located point is the centre of the peak cell (3, 1).
        let (x, y) = sol.max_pressure_location().expect("non-empty grid");
        assert!(
            (x - (3.0 + 0.5) * dx).abs() < 1e-9 && (y - (1.0 + 0.5) * dy).abs() < 1e-9,
            "stagnation point at cell (3,1) centre, got ({x}, {y})"
        );

        // Non-tautological re-scan: the located cell's pressure is ≥ every cell's.
        let peak = sol.pressure.at(3, 1);
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                assert!(
                    sol.pressure.at(i, j) <= peak,
                    "cell ({i},{j}) exceeds the located peak"
                );
            }
        }

        // The stagnation point (max) is a distinct cell from the suction peak (min).
        assert_ne!(
            sol.max_pressure_location(),
            sol.min_pressure_location(),
            "max and min pressure sit at distinct cells for a non-constant field"
        );
        // (The `None`-on-empty-grid branch mirrors min_pressure_location but is
        // unreachable through the public API — `Grid::new` panics on zero cells.)
    }

    // ----- EffectiveViscosity / SST in the SIMPLE driver ------------

    #[test]
    fn simple_with_laminar_matches_solve_simple() {
        // The Laminar effective-viscosity path is the historical solver.
        let grid = Grid::new(16, 16, 1.0, 1.0);
        let fluid = Fluid::new(1.0, 0.1);
        let bcs = Boundaries::lid_driven_cavity(1.0);
        let ctrl = SimpleControls::default();
        let baseline = solve_simple(&grid, &fluid, &bcs, &ctrl);
        let (sol, snap) = solve_simple_with(
            &grid, &fluid, &bcs, &ctrl, &EffectiveViscosity::Laminar,
        );
        assert!(snap.is_none(), "laminar path returns no turbulence snapshot");
        for k in 0..baseline.u.data.len() {
            assert!(
                (baseline.u.data[k] - sol.u.data[k]).abs() < 1e-12,
                "laminar path must match solve_simple exactly"
            );
        }
    }

    #[test]
    fn simple_with_sst_returns_a_field_with_positive_eddy_viscosity() {
        // Run a turbulent lid-driven cavity with k-ω SST and confirm
        // the eddy viscosity is everywhere finite and somewhere
        // strictly positive (the model must have produced *some* μ_t
        // from the lid-driven shear).
        use crate::turbulence::{SstModel, WallFunction, WallMask};
        let grid = Grid::new(20, 20, 1.0, 1.0);
        let fluid = Fluid::new(1.0, 1e-3); // Re = 1000
        let bcs = Boundaries::lid_driven_cavity(1.0);
        let ctrl = SimpleControls {
            max_iterations: 800,
            tolerance: 1e-4,
            relax_u: 0.5,
            relax_p: 0.2,
            ..SimpleControls::default()
        };
        let visc = EffectiveViscosity::Sst {
            model: SstModel::standard(),
            wall: WallFunction::standard(),
            walls: WallMask::closed(),
        };
        let (sol, snap) = solve_simple_with(&grid, &fluid, &bcs, &ctrl, &visc);
        let snap = snap.expect("SST path returns a snapshot");
        let nu_t = match &snap {
            crate::solver::TurbulenceSnapshot::Sst(s) => &s.nu_t,
            _ => panic!("expected SST snapshot"),
        };
        // Every value is finite ≥ 0.
        let mut max_nu_t = 0.0_f64;
        for &v in &nu_t.data {
            assert!(v.is_finite() && v >= 0.0, "ν_t must be finite ≥ 0: {v}");
            max_nu_t = max_nu_t.max(v);
        }
        // The lid-driven cavity has strong shear under the lid — the
        // model should produce a positive eddy viscosity somewhere.
        assert!(
            max_nu_t > 1e-8,
            "SST should produce a positive ν_t in a sheared flow"
        );
        // The mean flow is still bounded by the lid speed (sanity).
        let mut max_speed = 0.0_f64;
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                max_speed = max_speed.max(sol.speed_at_cell(i, j));
            }
        }
        assert!(
            max_speed <= 1.25,
            "interior speed {max_speed} should not blow up past the lid"
        );
    }

    #[test]
    fn simple_with_multigrid_reaches_the_same_steady_state_as_sor() {
        // The pressure-Poisson solver is an internal choice; the
        // converged SIMPLE solution must be the *same* whether SOR or
        // multigrid produced each p' field.
        use crate::multigrid::{MultigridControls, PressurePoissonSolver};
        let grid = Grid::new(32, 32, 1.0, 1.0);
        let fluid = Fluid::new(1.0, 0.05);
        let bcs = Boundaries::lid_driven_cavity(1.0);
        let sor_ctrl = SimpleControls {
            max_iterations: 4000,
            tolerance: 1e-6,
            ..SimpleControls::default()
        };
        let mg_ctrl = SimpleControls {
            pressure_solver: PressurePoissonSolver::Multigrid(
                MultigridControls {
                    max_cycles: 4,
                    tolerance: 1e-8,
                    ..MultigridControls::default()
                },
            ),
            ..sor_ctrl
        };
        let (sor_sol, _) = solve_simple_with(
            &grid, &fluid, &bcs, &sor_ctrl, &EffectiveViscosity::Laminar,
        );
        let (mg_sol, _) = solve_simple_with(
            &grid, &fluid, &bcs, &mg_ctrl, &EffectiveViscosity::Laminar,
        );
        assert!(sor_sol.converged, "SOR baseline must converge");
        assert!(mg_sol.converged, "multigrid path must converge");
        // The two velocity fields agree to a small absolute tolerance —
        // both have reached the same fixed point of the SIMPLE iteration.
        let mut max_diff = 0.0_f64;
        for k in 0..sor_sol.u.data.len() {
            max_diff =
                max_diff.max((sor_sol.u.data[k] - mg_sol.u.data[k]).abs());
        }
        assert!(
            max_diff < 5e-3,
            "SOR vs multigrid u differ by {max_diff} (should be ~1e-3)"
        );
    }

    #[test]
    fn simple_with_kepsilon_runs_and_keeps_k_positive() {
        // The k-ε coupling path: a channel-flow run with the standard
        // k-ε model. The flow still develops a parabolic-ish profile
        // (this is a coarse, high-ν setup so the k-ε turbulence stays
        // small) and the k field stays non-negative.
        use crate::turbulence::{KEpsilonModel, WallFunction, WallMask};
        let grid = Grid::new(30, 10, 3.0, 1.0);
        let fluid = Fluid::new(1.0, 0.05);
        let bcs = Boundaries::channel_flow(1.0);
        let ctrl = SimpleControls {
            max_iterations: 1500,
            tolerance: 1e-4,
            relax_u: 0.5,
            relax_p: 0.2,
            ..SimpleControls::default()
        };
        let visc = EffectiveViscosity::KEpsilon {
            model: KEpsilonModel::standard(),
            wall: WallFunction::standard(),
            walls: WallMask::channel(),
        };
        let (sol, snap) = solve_simple_with(&grid, &fluid, &bcs, &ctrl, &visc);
        let snap = snap.expect("k-ε path returns a snapshot");
        match snap {
            crate::solver::TurbulenceSnapshot::KEpsilon(t) => {
                for &k in &t.k.data {
                    assert!(k >= 0.0 && k.is_finite(), "k must stay ≥ 0: {k}");
                }
                for &eps in &t.epsilon.data {
                    assert!(eps > 0.0 && eps.is_finite(), "ε must stay > 0: {eps}");
                }
            }
            _ => panic!("expected k-ε snapshot"),
        }
        // The flow still peaks at the centre (slower-relaxed but
        // physically the same parabola family).
        let out_i = grid.nx - 3;
        let centre = sol.u_at_cell(out_i, grid.ny / 2);
        let near_wall = sol.u_at_cell(out_i, 0);
        assert!(
            centre > near_wall,
            "centreline must exceed wall velocity"
        );
    }
}
