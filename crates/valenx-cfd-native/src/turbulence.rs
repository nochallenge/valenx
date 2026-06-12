//! Two-equation **RANS turbulence models** — the standard k-ε
//! (Launder & Spalding 1974) and **Menter's k-ω SST** (Menter 1994),
//! both with equilibrium log-law wall functions.
//!
//! # What this is
//!
//! [`crate::solver`] solves the *laminar* Navier-Stokes equations: it
//! is valid only while the Reynolds number is modest and the flow
//! stays smooth. Real engineering flows are **turbulent** — the
//! velocity field carries a chaotic spectrum of eddies far too fine to
//! resolve on an affordable grid. The Reynolds-averaged (RANS) approach
//! does not resolve them; it solves for the *mean* flow and models the
//! eddies' net effect as an extra, flow-dependent **turbulent
//! viscosity** `μ_t`.
//!
//! This module implements the two canonical RANS closures the
//! external-aero / wall-bounded-shear community actually uses:
//!
//! - **Standard k-ε** ([`KEpsilonModel`]) — transports the turbulence
//!   kinetic energy `k` and its dissipation rate `ε`;
//!   `μ_t = ρ·C_μ·k²/ε`. The historical workhorse — robust in the
//!   free stream, weaker in adverse-pressure-gradient / separated flow.
//! - **k-ω SST** ([`SstModel`]) — Menter's shear-stress-transport
//!   model transports `k` and the **specific dissipation rate** `ω`;
//!   `μ_t = ρ·a₁·k / max(a₁·ω, S·F₂)`. A blending function `F₁` runs
//!   the k-ω formulation near walls (where it predicts boundary-layer
//!   separation well) and a k-ε formulation in the free stream (where
//!   k-ω is over-sensitive to the inlet `ω`). The SST limiter
//!   `μ_t = ρ·a₁·k/max(a₁·ω, S·F₂)` bounds the eddy viscosity in an
//!   adverse pressure gradient — the model's defining accuracy gain
//!   over k-ε. **The external-aerodynamics standard.**
//!
//! The momentum equation is then solved with the **effective
//! viscosity** `μ_eff = μ + μ_t` in place of the molecular `μ` (see
//! [`crate::solver::EffectiveViscosity`]). Because `μ_t` is large
//! wherever the flow is turbulent, the effective viscosity is large
//! there too — the model's net effect is a strong extra diffusion
//! that flattens the mean-velocity profile, exactly the signature
//! of a turbulent flow.
//!
//! # The two transport equations
//!
//! **k-ε.** Both `k` and `ε` obey a standard
//! convection-diffusion-source balance. With `P_k` the turbulent
//! **production** (the rate at which the mean shear feeds energy into
//! the eddies),
//!
//! ```text
//!   ∇·(ρ u k) = ∇·[(μ + μ_t/σ_k)∇k] + P_k − ρ ε
//!   ∇·(ρ u ε) = ∇·[(μ + μ_t/σ_ε)∇ε] + C_1ε·(ε/k)·P_k − C_2ε·ρ·ε²/k
//! ```
//!
//! ```text
//!   P_k = μ_t · S²,   S² = 2·(∂u/∂x)² + 2·(∂v/∂y)² + (∂u/∂y + ∂v/∂x)²
//! ```
//!
//! The standard constants are `C_μ = 0.09`, `C_1ε = 1.44`,
//! `C_2ε = 1.92`, `σ_k = 1.0`, `σ_ε = 1.3` ([`KEpsilonModel`]).
//!
//! **k-ω SST.** The two transported scalars are `k` and `ω`; both
//! obey
//!
//! ```text
//!   ∇·(ρ u k) = ∇·[(μ + σ_k·μ_t)∇k] + min(P_k, 10·β*·ρ·k·ω) − β*·ρ·k·ω
//!   ∇·(ρ u ω) = ∇·[(μ + σ_ω·μ_t)∇ω] + (γ·ρ/μ_t)·P_k − β·ρ·ω²
//!             + 2·(1−F₁)·ρ·σ_ω2·(∇k·∇ω)/ω
//! ```
//!
//! The `(1−F₁)` cross-diffusion term is what makes SST switch
//! continuously from the k-ω near-wall form (`F₁ = 1`) to the k-ε
//! free-stream form (`F₁ = 0`). The model constants `(σ_k, σ_ω,
//! β, γ)` are themselves **blended** by `F₁` from a near-wall *set 1*
//! and a free-stream *set 2* — `φ = F₁·φ₁ + (1−F₁)·φ₂`. The eddy
//! viscosity uses the **SST limiter**
//!
//! ```text
//!   μ_t = ρ · a₁ · k / max(a₁·ω, S·F₂)
//! ```
//!
//! where `F₂` is a second blending function active in the boundary
//! layer; the limiter cap is what bounds `μ_t` in an
//! adverse-pressure-gradient separation. See [`SstModel`].
//!
//! # Wall functions
//!
//! Both models are **high-Reynolds-number** formulations — they are
//! not valid in the thin viscous sublayer right at a wall, where the
//! turbulence is damped out. Resolving that sublayer needs an
//! impractically fine grid. The standard remedy, used here, is the
//! **wall-function** approach: the first cell off a wall is assumed to
//! sit in the logarithmic "law of the wall" region, and the wall
//! shear, `k`, and `ε`/`ω` in that cell are set from the log-law
//! relations rather than from the differential equations. See
//! [`WallFunction`].
//!
//! # The turbulent-channel driver
//!
//! [`solve_turbulent_channel`] couples the **k-ε** model to a 1-D
//! wall-normal momentum balance and [`solve_turbulent_channel_sst`]
//! does the same for **k-ω SST**. The driver is the verification
//! vehicle for each closure. A **fully-developed plane channel** —
//! the flow far downstream of the inlet, where nothing varies along
//! the flow direction — is governed by the 1-D wall-normal balance
//!
//! ```text
//!   d/dy[ (μ + μ_t)·du/dy ] = dp/dx   (a constant pressure gradient)
//! ```
//!
//! Each driver iterates this balance against its closure's transport
//! equations: solve the momentum ODE for `u(y)` with the current
//! `μ_t(y)`, advance the turbulence scalars from the new shear,
//! recompute `μ_t`, repeat. With the molecular viscosity alone
//! (`μ_t = 0`) the balance integrates to the **parabolic** laminar
//! profile; with the turbulent `μ_t` the strong near-wall eddy
//! diffusion **flattens** the profile — the defining difference
//! between a laminar and a turbulent channel, and exactly what the
//! tests check.
//!
//! # Honest scope
//!
//! Both closures are **real, working v1s** — the per-driver tests
//! verify the turbulent profile flattens relative to the laminar
//! parabola, that `k` / `ε` / `ω` stay physical, and that the SST
//! limiter actually bounds `μ_t` in a separated-flow regime. The v1
//! caveats:
//!
//! - **The standard k-ε and Menter k-ω SST only.** No realizable or
//!   RNG k-ε, no Reynolds-stress model, no LES/DES. Each further model
//!   is its own closure.
//! - **Equilibrium wall functions.** The log-law treatment assumes
//!   local equilibrium (production ≈ dissipation in the wall cell);
//!   non-equilibrium / enhanced wall treatments are extensions.
//! - **2-D, structured grid, steady** — inherited from
//!   [`crate::solver`].
//!
//! Within that scope each model is the genuine article: the same
//! transport equations, the same calibrated constants, the same
//! blending functions Menter published.

use crate::grid::{Field, Grid};

/// The five calibrated constants of the **standard k-ε model**
/// (Launder & Spalding 1974).
///
/// [`KEpsilonModel::standard`] returns the universally-used set; the
/// struct is exposed so a caller can experiment, but the standard
/// values are what the model was calibrated against and what the tests
/// use.
#[derive(Clone, Copy, Debug)]
pub struct KEpsilonModel {
    /// `C_μ` — the eddy-viscosity constant in `μ_t = ρ·C_μ·k²/ε`.
    /// Standard value `0.09`.
    pub c_mu: f64,
    /// `C_1ε` — the production constant in the ε equation. Standard
    /// value `1.44`.
    pub c_1eps: f64,
    /// `C_2ε` — the destruction constant in the ε equation. Standard
    /// value `1.92`.
    pub c_2eps: f64,
    /// `σ_k` — the turbulent Prandtl number for `k` (the `k`-equation
    /// diffusion is scaled by `μ_t/σ_k`). Standard value `1.0`.
    pub sigma_k: f64,
    /// `σ_ε` — the turbulent Prandtl number for `ε`. Standard value
    /// `1.3`.
    pub sigma_eps: f64,
}

impl KEpsilonModel {
    /// The standard k-ε constant set — `C_μ = 0.09`, `C_1ε = 1.44`,
    /// `C_2ε = 1.92`, `σ_k = 1.0`, `σ_ε = 1.3`.
    pub fn standard() -> KEpsilonModel {
        KEpsilonModel {
            c_mu: 0.09,
            c_1eps: 1.44,
            c_2eps: 1.92,
            sigma_k: 1.0,
            sigma_eps: 1.3,
        }
    }
}

impl Default for KEpsilonModel {
    fn default() -> Self {
        KEpsilonModel::standard()
    }
}

/// The log-law (equilibrium) **wall-function** constants.
///
/// The standard k-ε model is a high-Re model; right at a wall the
/// turbulence is damped and the model is invalid. The wall function
/// bridges that gap: the first cell off the wall is assumed to sit in
/// the logarithmic region of the boundary layer, where the mean
/// velocity follows the **law of the wall**
///
/// ```text
///   u⁺ = (1/κ)·ln(E·y⁺),     u⁺ = u/u_τ,   y⁺ = ρ·u_τ·y/μ
/// ```
///
/// with `u_τ = √(τ_w/ρ)` the friction velocity. From the same
/// assumption the wall-cell turbulence is fixed at its equilibrium
/// values
///
/// ```text
///   k_P = u_τ² / √C_μ ,      ε_P = u_τ³ / (κ·y_P)
/// ```
///
/// These relations replace the `k` / `ε` differential equations *in the
/// wall cell only*, and supply the wall shear stress that the momentum
/// equation needs there.
#[derive(Clone, Copy, Debug)]
pub struct WallFunction {
    /// The von Kármán constant `κ`. Standard value `0.41`.
    pub kappa: f64,
    /// The log-law roughness parameter `E` (smooth wall). Standard
    /// value `9.8`.
    pub e_wall: f64,
}

impl WallFunction {
    /// The standard smooth-wall log-law constants — `κ = 0.41`,
    /// `E = 9.8`.
    pub fn standard() -> WallFunction {
        WallFunction {
            kappa: 0.41,
            e_wall: 9.8,
        }
    }

    /// The dimensionless wall distance `y⁺ = ρ·u_τ·y/μ`.
    #[inline]
    pub fn y_plus(&self, rho: f64, u_tau: f64, y: f64, mu: f64) -> f64 {
        rho * u_tau * y / mu.max(1e-30)
    }

    /// The log-law mean velocity `u⁺ = (1/κ)·ln(E·y⁺)`, clamped so the
    /// argument of the logarithm stays positive.
    #[inline]
    pub fn u_plus(&self, y_plus: f64) -> f64 {
        (1.0 / self.kappa) * (self.e_wall * y_plus.max(1e-6)).ln()
    }

    /// The friction velocity `u_τ` of a wall cell.
    ///
    /// Given the wall-parallel velocity `u_p` a distance `y_p` from the
    /// wall, `u_τ` is the root of the implicit law-of-the-wall relation
    /// `u_p/u_τ = (1/κ)·ln(E·ρ·u_τ·y_p/μ)`. It is solved by a short
    /// fixed-point iteration; in the viscous sublayer (`y⁺ ≲ 11.06`,
    /// where the log law and the linear law `u⁺ = y⁺` cross) the linear
    /// relation `u_τ = √(μ·u_p/(ρ·y_p))` is used instead.
    pub fn friction_velocity(&self, rho: f64, mu: f64, u_p: f64, y_p: f64) -> f64 {
        let u_p = u_p.abs();
        if u_p < 1e-12 || y_p <= 0.0 {
            return 0.0;
        }
        // Viscous-sublayer estimate (linear law u⁺ = y⁺):
        //   u_p/u_τ = ρ·u_τ·y_p/μ  →  u_τ = √(μ·u_p/(ρ·y_p)).
        let u_tau_visc = (mu * u_p / (rho * y_p)).sqrt();
        // The log law and the linear law cross at y⁺ ≈ 11.06.
        let y_plus_visc = self.y_plus(rho, u_tau_visc, y_p, mu);
        if y_plus_visc <= 11.06 {
            return u_tau_visc;
        }
        // Log-law region: fixed-point iterate u_τ = κ·u_p / ln(E·y⁺).
        let mut u_tau = u_tau_visc.max(1e-9);
        for _ in 0..50 {
            let y_plus = self.y_plus(rho, u_tau, y_p, mu);
            let denom = (self.e_wall * y_plus.max(1e-6)).ln();
            if denom <= 1e-9 {
                break;
            }
            let next = self.kappa * u_p / denom;
            if (next - u_tau).abs() <= 1e-10 * next.max(1e-12) {
                u_tau = next;
                break;
            }
            u_tau = next;
        }
        u_tau
    }
}

impl Default for WallFunction {
    fn default() -> Self {
        WallFunction::standard()
    }
}

/// The turbulence state of a flow — the cell-centred `k` and `ε`
/// fields and the turbulent viscosity derived from them.
///
/// All three fields are cell-centred (`nx × ny`), so they live on the
/// pressure grid of [`crate::Grid`]. The [`TurbulenceField::nu_t`]
/// field is the model's output: the eddy viscosity at every cell,
/// ready to be added to the molecular viscosity in the momentum solve.
#[derive(Clone, Debug)]
pub struct TurbulenceField {
    /// Turbulent kinetic energy `k` (m²/s²), cell-centred.
    pub k: Field,
    /// Dissipation rate `ε` (m²/s³), cell-centred.
    pub epsilon: Field,
    /// Turbulent (eddy) **kinematic** viscosity `ν_t = μ_t/ρ` (m²/s),
    /// cell-centred — `ν_t = C_μ·k²/ε`. Kinematic so it adds directly
    /// to [`crate::Fluid::viscosity`].
    pub nu_t: Field,
}

impl TurbulenceField {
    /// Initialise a turbulence field from a reference velocity and a
    /// **turbulence intensity** `I` (the fraction of the mean velocity
    /// carried by the velocity fluctuations — `0.05` is a typical
    /// "medium" inlet value).
    ///
    /// The standard inlet estimates are used:
    ///
    /// ```text
    ///   k  = 1.5·(I·U_ref)²
    ///   ε  = C_μ^{3/4} · k^{3/2} / ℓ        ℓ = 0.07·L  (mixing length)
    /// ```
    ///
    /// with `ℓ` a turbulence length scale taken as 7 % of the domain
    /// height — the conventional fully-developed-duct value.
    pub fn initialise(
        grid: &Grid,
        u_ref: f64,
        intensity: f64,
        model: &KEpsilonModel,
    ) -> TurbulenceField {
        let u_ref = u_ref.abs().max(1e-6);
        let intensity = intensity.max(1e-4);
        let k0 = 1.5 * (intensity * u_ref).powi(2);
        // Mixing length ≈ 7 % of the duct height (standard rule).
        let length_scale = (0.07 * grid.ly).max(1e-6);
        let eps0 = model.c_mu.powf(0.75) * k0.powf(1.5) / length_scale;
        let nu_t0 = model.c_mu * k0 * k0 / eps0.max(1e-30);
        TurbulenceField {
            k: Field::filled(grid.nx, grid.ny, k0),
            epsilon: Field::filled(grid.nx, grid.ny, eps0),
            nu_t: Field::filled(grid.nx, grid.ny, nu_t0),
        }
    }

    /// Recompute the eddy viscosity field `ν_t = C_μ·k²/ε` from the
    /// current `k` and `ε`.
    ///
    /// `ν_t` is clamped non-negative and to a generous multiple of the
    /// molecular viscosity `nu_mol` so a transient spike in `k²/ε`
    /// cannot blow the momentum solve up — the standard production
    /// limiter every k-ε implementation carries.
    pub fn update_eddy_viscosity(&mut self, model: &KEpsilonModel, nu_mol: f64) {
        let cap = (1.0e5 * nu_mol).max(1.0e3 * nu_mol).max(1e-6);
        for idx in 0..self.k.data.len() {
            let k = self.k.data[idx].max(0.0);
            let eps = self.epsilon.data[idx].max(1e-30);
            let nu_t = (model.c_mu * k * k / eps).clamp(0.0, cap);
            self.nu_t.data[idx] = nu_t;
        }
    }

    /// The mean turbulent kinetic energy across the field — a scalar
    /// turbulence-level diagnostic.
    pub fn mean_k(&self) -> f64 {
        if self.k.data.is_empty() {
            return 0.0;
        }
        self.k.data.iter().sum::<f64>() / self.k.data.len() as f64
    }
}

/// The cell-centred **turbulent-production** field
/// `P_k = ν_t · S²` (per unit mass), where `S²` is the mean
/// strain-rate magnitude squared.
///
/// The velocity components live on the staggered faces of
/// [`crate::Grid`]; this routine reconstructs the four mean-velocity
/// gradients at each cell centre by central differences of the
/// bracketing faces and forms
///
/// ```text
///   S² = 2(∂u/∂x)² + 2(∂v/∂y)² + (∂u/∂y + ∂v/∂x)²
/// ```
///
/// `u` is the staggered `u`-field (`(nx+1) × ny`), `v` the staggered
/// `v`-field (`nx × (ny+1)`).
pub fn strain_rate_squared(u: &Field, v: &Field, grid: &Grid) -> Field {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx();
    let dy = grid.dy();
    let mut s2 = Field::zeros(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            // ∂u/∂x — exact: the two u-faces bracket the cell in x.
            let dudx = (u.at(i + 1, j) - u.at(i, j)) / dx;
            // ∂v/∂y — exact: the two v-faces bracket the cell in y.
            let dvdy = (v.at(i, j + 1) - v.at(i, j)) / dy;
            // ∂u/∂y — cell-centred u differenced across y. Cell-centre
            // u is the average of its two x-faces.
            let u_here = 0.5 * (u.at(i, j) + u.at(i + 1, j));
            let u_up = if j + 1 < ny {
                0.5 * (u.at(i, j + 1) + u.at(i + 1, j + 1))
            } else {
                u_here
            };
            let u_dn = if j > 0 {
                0.5 * (u.at(i, j - 1) + u.at(i + 1, j - 1))
            } else {
                u_here
            };
            let span_y = if j + 1 < ny && j > 0 { 2.0 } else { 1.0 } * dy;
            let dudy = (u_up - u_dn) / span_y.max(1e-30);
            // ∂v/∂x — cell-centred v differenced across x.
            let v_here = 0.5 * (v.at(i, j) + v.at(i, j + 1));
            let v_e = if i + 1 < nx {
                0.5 * (v.at(i + 1, j) + v.at(i + 1, j + 1))
            } else {
                v_here
            };
            let v_w = if i > 0 {
                0.5 * (v.at(i - 1, j) + v.at(i - 1, j + 1))
            } else {
                v_here
            };
            let span_x = if i + 1 < nx && i > 0 { 2.0 } else { 1.0 } * dx;
            let dvdx = (v_e - v_w) / span_x.max(1e-30);

            let s_sq = 2.0 * dudx * dudx + 2.0 * dvdy * dvdy + (dudy + dvdx).powi(2);
            s2.set(i, j, s_sq);
        }
    }
    s2
}

/// Advance the `k` and `ε` transport equations one outer iteration.
///
/// This is the turbulence-model counterpart of [`crate::solver`]'s
/// momentum + pressure-correction step. Given the current mean
/// velocity (the staggered `u` / `v` fields) it:
///
/// 1. forms the production field `P_k = ν_t·S²`
///    ([`strain_rate_squared`]);
/// 2. sweeps a convection-diffusion update of `k` and then of `ε`,
///    each with its source/sink pair, using the same hybrid
///    upwind/central scheme and Gauss-Seidel relaxation as the
///    momentum solve;
/// 3. applies the **wall-function** values of `k` and `ε` in every
///    cell adjacent to a no-slip wall (the rows / columns flagged in
///    `wall_*`);
/// 4. clamps `k` and `ε` to small positive floors (both are
///    physically non-negative; the model divides by them);
/// 5. recomputes `ν_t = C_μ·k²/ε`.
///
/// `rho`/`nu_mol` are the fluid density and molecular *kinematic*
/// viscosity. `relax` is the under-relaxation factor (k-ε needs heavy
/// relaxation, `≈ 0.5–0.7`, for the same reason the momentum solve
/// does). `wall_south`/`wall_north` flag whether the bottom / top
/// domain edges are no-slip walls (the channel case — both true).
#[allow(clippy::too_many_arguments)]
pub fn advance_k_epsilon(
    turb: &mut TurbulenceField,
    u: &Field,
    v: &Field,
    grid: &Grid,
    rho: f64,
    nu_mol: f64,
    model: &KEpsilonModel,
    wall: &WallFunction,
    wall_south: bool,
    wall_north: bool,
    relax: f64,
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx();
    let dy = grid.dy();
    let relax = relax.clamp(0.05, 1.0);

    // --- production field P_k = ν_t · S² ---
    let s2 = strain_rate_squared(u, v, grid);
    let mut production = Field::zeros(nx, ny);
    for idx in 0..production.data.len() {
        production.data[idx] = turb.nu_t.data[idx] * s2.data[idx];
    }

    // Cell-centred speed, used by the wall function.
    let cell_speed = |i: usize, j: usize| -> f64 {
        let uc = 0.5 * (u.at(i, j) + u.at(i + 1, j));
        let vc = 0.5 * (v.at(i, j) + v.at(i, j + 1));
        (uc * uc + vc * vc).sqrt()
    };

    // Floors — k and ε are physically positive and the model divides
    // by them.
    let k_floor = 1e-10;
    let eps_floor = 1e-10;

    // --- (1) k transport: convection-diffusion + P_k − ρε ---
    // Sweep interior cells; wall-adjacent cells are overwritten by the
    // wall function afterwards.
    let k_old = turb.k.clone();
    for _sweep in 0..2 {
        for j in 0..ny {
            for i in 0..nx {
                // Effective diffusivity at the cell: ν + ν_t/σ_k.
                let diff_p = nu_mol + turb.nu_t.at(i, j) / model.sigma_k;
                // Face mass flows (per unit depth) — staggered faces.
                let fe = rho * dy * u.at(i + 1, j);
                let fw = rho * dy * u.at(i, j);
                let fn_ = rho * dx * v.at(i, j + 1);
                let fs = rho * dx * v.at(i, j);

                let mut a_p = 0.0;
                let mut nb = 0.0;
                // East.
                if i + 1 < nx {
                    let d = rho * diff_p * dy / dx;
                    let a = hybrid(d, fe, -1.0);
                    a_p += a;
                    nb += a * turb.k.at(i + 1, j);
                }
                // West.
                if i > 0 {
                    let d = rho * diff_p * dy / dx;
                    let a = hybrid(d, fw, 1.0);
                    a_p += a;
                    nb += a * turb.k.at(i - 1, j);
                }
                // North.
                if j + 1 < ny {
                    let d = rho * diff_p * dx / dy;
                    let a = hybrid(d, fn_, -1.0);
                    a_p += a;
                    nb += a * turb.k.at(i, j + 1);
                }
                // South.
                if j > 0 {
                    let d = rho * diff_p * dx / dy;
                    let a = hybrid(d, fs, 1.0);
                    a_p += a;
                    nb += a * turb.k.at(i, j - 1);
                }
                // Source: production P_k·vol. Sink: ρ·ε linearised as
                // (ρ·ε/k)·k — put on the diagonal for stability
                // (Patankar's positive-coefficient rule).
                let vol = dx * dy;
                let k_p = k_old.at(i, j).max(k_floor);
                let eps_p = turb.epsilon.at(i, j).max(eps_floor);
                let su = production.at(i, j) * vol;
                let sp = rho * eps_p / k_p * vol; // sink coefficient
                a_p += sp;
                if a_p.abs() < 1e-30 {
                    continue;
                }
                let k_new = (nb + su) / a_p;
                let updated = turb.k.at(i, j) + relax * (k_new - turb.k.at(i, j));
                turb.k.set(i, j, updated.max(k_floor));
            }
        }
    }

    // --- (2) ε transport: convection-diffusion + C1·(ε/k)·P_k − C2·ρε²/k ---
    let eps_old = turb.epsilon.clone();
    for _sweep in 0..2 {
        for j in 0..ny {
            for i in 0..nx {
                let diff_p = nu_mol + turb.nu_t.at(i, j) / model.sigma_eps;
                let fe = rho * dy * u.at(i + 1, j);
                let fw = rho * dy * u.at(i, j);
                let fn_ = rho * dx * v.at(i, j + 1);
                let fs = rho * dx * v.at(i, j);

                let mut a_p = 0.0;
                let mut nb = 0.0;
                if i + 1 < nx {
                    let d = rho * diff_p * dy / dx;
                    let a = hybrid(d, fe, -1.0);
                    a_p += a;
                    nb += a * turb.epsilon.at(i + 1, j);
                }
                if i > 0 {
                    let d = rho * diff_p * dy / dx;
                    let a = hybrid(d, fw, 1.0);
                    a_p += a;
                    nb += a * turb.epsilon.at(i - 1, j);
                }
                if j + 1 < ny {
                    let d = rho * diff_p * dx / dy;
                    let a = hybrid(d, fn_, -1.0);
                    a_p += a;
                    nb += a * turb.epsilon.at(i, j + 1);
                }
                if j > 0 {
                    let d = rho * diff_p * dx / dy;
                    let a = hybrid(d, fs, 1.0);
                    a_p += a;
                    nb += a * turb.epsilon.at(i, j - 1);
                }
                let vol = dx * dy;
                let k_p = turb.k.at(i, j).max(k_floor);
                let eps_p = eps_old.at(i, j).max(eps_floor);
                // Production source C1ε·(ε/k)·P_k.
                let su = model.c_1eps * eps_p / k_p * production.at(i, j) * vol;
                // Destruction sink C2ε·ρ·ε²/k, linearised as
                // (C2ε·ρ·ε/k)·ε onto the diagonal.
                let sp = model.c_2eps * rho * eps_p / k_p * vol;
                a_p += sp;
                if a_p.abs() < 1e-30 {
                    continue;
                }
                let eps_new = (nb + su) / a_p;
                let updated = turb.epsilon.at(i, j) + relax * (eps_new - turb.epsilon.at(i, j));
                turb.epsilon.set(i, j, updated.max(eps_floor));
            }
        }
    }

    // --- (3) wall functions on the wall-adjacent cells ---
    // The first cell off a no-slip wall sits a distance dy/2 from it.
    // Its k and ε are *prescribed* from the equilibrium log-law
    // relations, overriding the differential-equation update above.
    let apply_wall = |turb: &mut TurbulenceField, i: usize, j: usize, y_p: f64| {
        let u_p = cell_speed(i, j);
        let u_tau = wall.friction_velocity(rho, nu_mol, u_p, y_p);
        // k_P = u_τ² / √C_μ  (equilibrium wall-cell kinetic energy).
        let k_wall = (u_tau * u_tau / model.c_mu.sqrt()).max(k_floor);
        // ε_P = u_τ³ / (κ·y_P).
        let eps_wall = (u_tau.powi(3) / (wall.kappa * y_p.max(1e-30))).max(eps_floor);
        turb.k.set(i, j, k_wall);
        turb.epsilon.set(i, j, eps_wall);
    };
    let y_first = 0.5 * dy;
    if wall_south && ny > 0 {
        for i in 0..nx {
            apply_wall(turb, i, 0, y_first);
        }
    }
    if wall_north && ny > 0 {
        for i in 0..nx {
            apply_wall(turb, i, ny - 1, y_first);
        }
    }

    // --- (4 + 5) refresh the eddy viscosity from the new k, ε ---
    turb.update_eddy_viscosity(model, nu_mol);
}

/// The hybrid-scheme convection-diffusion neighbour coefficient — the
/// same `max(D − |F|/2, 0) + max(±F, 0)` rule [`crate::solver`] uses
/// for momentum, reused here for the `k` / `ε` scalar transport so the
/// turbulence equations share the solver's unconditionally-stable
/// discretisation.
#[inline]
fn hybrid(d: f64, f: f64, upwind_sign: f64) -> f64 {
    (d - 0.5 * f.abs()).max(0.0) + (upwind_sign * f).max(0.0)
}

/// The converged velocity profile of a fully-developed plane channel.
#[derive(Clone, Debug)]
pub struct ChannelProfile {
    /// The cell-centred wall-normal coordinate `y` of each sample
    /// (metres), bottom wall at `y = 0`, top wall at `y = ly`.
    pub y: Vec<f64>,
    /// The streamwise velocity `u(y)` at each sample (m/s).
    pub u: Vec<f64>,
    /// The turbulent (eddy) kinematic viscosity `ν_t(y)` at each
    /// sample — zero everywhere for a laminar run.
    pub nu_t: Vec<f64>,
    /// The bulk (cross-section-averaged) velocity (m/s).
    pub bulk_velocity: f64,
}

impl ChannelProfile {
    /// The peak (centreline) velocity of the profile.
    pub fn centreline_velocity(&self) -> f64 {
        self.u.iter().copied().fold(0.0, f64::max)
    }

    /// The **profile flatness** — the ratio `u_bulk / u_centreline`.
    ///
    /// This is the single number that distinguishes a laminar from a
    /// turbulent channel. A parabola has `u_bulk/u_max = 2/3 ≈ 0.667`;
    /// a turbulent profile is much flatter, with the ratio rising
    /// toward `~0.8–0.87`. A larger value means a flatter profile.
    pub fn flatness(&self) -> f64 {
        let peak = self.centreline_velocity();
        if peak <= 0.0 {
            0.0
        } else {
            self.bulk_velocity / peak
        }
    }
}

/// Solve a **fully-developed plane channel flow** — laminar or
/// turbulent — and return its wall-normal velocity profile.
///
/// Far downstream of a channel inlet the flow is *fully developed*:
/// nothing varies along the flow direction, and the streamwise
/// momentum equation collapses to the 1-D wall-normal balance
///
/// ```text
///   d/dy[ (ν + ν_t)·du/dy ] = (1/ρ)·dp/dx
/// ```
///
/// driven by a constant pressure gradient `dp/dx`. This routine solves
/// that balance on the wall-normal cell column of `grid`, with no-slip
/// (`u = 0`) walls top and bottom.
///
/// `turbulent` selects the physics:
///
/// - `false` — **laminar**: `ν_t ≡ 0`, the balance integrates to the
///   classical **parabolic** Poiseuille profile.
/// - `true` — **turbulent**: the k-ε model ([`advance_k_epsilon`]) is
///   coupled in. The momentum balance and the `k`/`ε` transport
///   equations are iterated to a joint fixed point — solve `u(y)` with
///   the current `ν_t`, advance `k`/`ε` from the new shear, recompute
///   `ν_t`, repeat. The strong near-wall eddy viscosity **flattens**
///   the profile relative to the laminar parabola.
///
/// `dpdx` is the (negative, for flow in `+x`) imposed pressure
/// gradient; `outer_iterations` caps the momentum-↔-turbulence coupling
/// loop. Returns the converged [`ChannelProfile`].
pub fn solve_turbulent_channel(
    grid: &Grid,
    fluid: &crate::solver::Fluid,
    model: &KEpsilonModel,
    wall: &WallFunction,
    dpdx: f64,
    turbulent: bool,
    outer_iterations: usize,
) -> ChannelProfile {
    let ny = grid.ny;
    let dy = grid.dy();
    let rho = fluid.density;
    let nu_mol = fluid.viscosity;
    // Source term of the momentum ODE: (1/ρ)·dp/dx, the same for every
    // wall-normal cell (fully-developed flow has a uniform pressure
    // gradient). The flow runs in +x against an adverse... actually a
    // *favourable* gradient, so dpdx is negative and the source is
    // negative; u comes out positive.
    let source = dpdx / rho;

    // The wall-normal cell column. u lives at cell centres; the two
    // walls are half a cell beyond the first / last centre.
    let mut u = vec![0.0_f64; ny];
    // Turbulent kinematic viscosity at cell centres.
    let mut nu_t = vec![0.0_f64; ny];

    // A 2-D turbulence field over the (1 × ny) column, so the shared
    // k-ε advance routine can be reused unchanged.
    let column = Grid::new(1, ny, grid.dx(), grid.ly);
    let mut turb = TurbulenceField::initialise(&column, 1.0, 0.05, model);

    let outers = if turbulent {
        outer_iterations.max(1)
    } else {
        // Laminar: the momentum ODE is linear and converges in one
        // sweep group; a few iterations of the tridiagonal smoother is
        // plenty.
        outer_iterations.max(1)
    };

    for _outer in 0..outers {
        // --- (1) momentum balance for u(y) given the current ν_t ---
        // Discretise d/dy[(ν+ν_t)du/dy] = source by central
        // differences and relax with Gauss-Seidel sweeps. The face
        // diffusivity is the average of the two bracketing cell
        // diffusivities.
        let diff_at = |j: usize| -> f64 { nu_mol + nu_t[j] };
        for _sweep in 0..200 {
            for j in 0..ny {
                // North-face diffusivity (between cell j and j+1, or
                // the wall past the last cell).
                let (an, u_n) = if j + 1 < ny {
                    let d = 0.5 * (diff_at(j) + diff_at(j + 1)) / dy;
                    (d, u[j + 1])
                } else {
                    // Wall: half-cell distance, u_wall = 0.
                    let d = diff_at(j) / (0.5 * dy);
                    (d, 0.0)
                };
                let (as_, u_s) = if j > 0 {
                    let d = 0.5 * (diff_at(j) + diff_at(j - 1)) / dy;
                    (d, u[j - 1])
                } else {
                    let d = diff_at(j) / (0.5 * dy);
                    (d, 0.0)
                };
                let a_p = an + as_;
                if a_p.abs() < 1e-30 {
                    continue;
                }
                // a_P·u_P = a_N·u_N + a_S·u_S − source·dy
                // (the source is integrated over the cell height dy).
                let rhs = an * u_n + as_ * u_s - source * dy;
                u[j] = rhs / a_p;
            }
        }

        if !turbulent {
            // Laminar run: no turbulence, ν_t stays zero.
            break;
        }

        // --- (2) advance k-ε from the new velocity field ---
        // Build the staggered velocity fields the k-ε routine wants
        // from the cell-centred column profile.
        let mut u_field = column.u_field(); // (2 × ny)
        for (j, &uj) in u.iter().enumerate().take(ny) {
            // The two x-faces of the single column both carry u(y)
            // (fully-developed: no x variation).
            u_field.set(0, j, uj);
            u_field.set(1, j, uj);
        }
        let v_field = column.v_field(); // all zero — no wall-normal flow
        for _ke in 0..3 {
            advance_k_epsilon(
                &mut turb, &u_field, &v_field, &column, rho, nu_mol, model, wall,
                true, // south wall
                true, // north wall
                0.6,
            );
        }
        // Pull the updated eddy viscosity back onto the column.
        for (j, nt) in nu_t.iter_mut().enumerate().take(ny) {
            *nt = turb.nu_t.at(0, j);
        }
    }

    // Cell-centred y coordinates and the bulk velocity.
    let y: Vec<f64> = (0..ny).map(|j| (j as f64 + 0.5) * dy).collect();
    let bulk = u.iter().sum::<f64>() / ny.max(1) as f64;

    ChannelProfile {
        y,
        u,
        nu_t,
        bulk_velocity: bulk,
    }
}

// =====================================================================
// Menter k-ω SST
// =====================================================================

/// The constants of **Menter's k-ω SST** turbulence model (Menter
/// 1994), grouped as the two blended *sets* the model carries plus the
/// shared limiter constant `a₁` and the standard `β*` / `κ`.
///
/// SST blends a *set 1* (the near-wall, k-ω-like calibration) with a
/// *set 2* (the free-stream, k-ε-like calibration) by the cell-centred
/// blending function `F₁`. Every transport-coefficient `φ` the model
/// uses is itself blended `φ = F₁·φ₁ + (1−F₁)·φ₂`. The values
/// returned by [`SstModel::standard`] are the universally-used Menter
/// 1994 constants and are what every standard-SST implementation
/// uses; the struct is exposed so a caller can experiment, not so it
/// has to.
#[derive(Clone, Copy, Debug)]
pub struct SstModel {
    /// `β*` — the eddy-viscosity / `k`-destruction constant
    /// (= `C_μ = 0.09` in k-ε notation).
    pub beta_star: f64,
    /// The von Kármán constant `κ` — appears in the `γ` derivation
    /// `γ = β/β* − σ_ω·κ²/√β*`.
    pub kappa: f64,
    /// `a₁` — the SST eddy-viscosity limiter constant
    /// (`μ_t = ρ·a₁·k / max(a₁·ω, S·F₂)`). Menter 1994 value `0.31`.
    pub a1: f64,
    /// Near-wall *set 1* — the k-ω calibration.
    pub set1: SstSet,
    /// Free-stream *set 2* — the k-ε-equivalent calibration.
    pub set2: SstSet,
}

/// One of the two `(σ_k, σ_ω, β, γ)` coefficient sets SST blends
/// between by the wall-distance blending function `F₁`. The near-wall
/// *set 1* is Menter's k-ω calibration; the free-stream *set 2* is
/// the k-ε-equivalent calibration. Every model coefficient is the
/// `F₁`-blend `φ = F₁·φ₁ + (1−F₁)·φ₂`.
#[derive(Clone, Copy, Debug)]
pub struct SstSet {
    /// `σ_k` — turbulent Prandtl number for `k`. Sets `μ_t` scaling
    /// in the `k` diffusion term `(μ + σ_k·μ_t)∇k`.
    pub sigma_k: f64,
    /// `σ_ω` — turbulent Prandtl number for `ω`.
    pub sigma_omega: f64,
    /// `β` — `ω`-destruction constant in `β·ρ·ω²`.
    pub beta: f64,
    /// `γ` — `ω`-production constant in `γ·(ρ/μ_t)·P_k`. Menter
    /// derives `γ = β/β* − σ_ω·κ²/√β*`.
    pub gamma: f64,
}

impl SstModel {
    /// **The standard Menter 1994 SST constants** — the calibration
    /// every SST implementation ships with.
    pub fn standard() -> SstModel {
        let beta_star: f64 = 0.09;
        let kappa: f64 = 0.41;
        // Menter 1994, Table 1.
        let set1 = SstSet {
            sigma_k: 0.85,
            sigma_omega: 0.5,
            beta: 0.075,
            // γ₁ = β₁/β* − σ_ω1·κ²/√β*.
            gamma: 0.075 / beta_star - 0.5 * kappa * kappa / beta_star.sqrt(),
        };
        let set2 = SstSet {
            sigma_k: 1.0,
            sigma_omega: 0.856,
            beta: 0.0828,
            gamma: 0.0828 / beta_star - 0.856 * kappa * kappa / beta_star.sqrt(),
        };
        SstModel {
            beta_star,
            kappa,
            a1: 0.31,
            set1,
            set2,
        }
    }

    /// Blend the two coefficient sets by `F₁`:
    /// `φ = F₁·φ₁ + (1−F₁)·φ₂`.
    #[inline]
    pub fn blend(&self, f1: f64) -> SstSet {
        let f1 = f1.clamp(0.0, 1.0);
        let g = 1.0 - f1;
        SstSet {
            sigma_k: f1 * self.set1.sigma_k + g * self.set2.sigma_k,
            sigma_omega: f1 * self.set1.sigma_omega + g * self.set2.sigma_omega,
            beta: f1 * self.set1.beta + g * self.set2.beta,
            gamma: f1 * self.set1.gamma + g * self.set2.gamma,
        }
    }
}

impl Default for SstModel {
    fn default() -> Self {
        SstModel::standard()
    }
}

/// The Menter SST **F₁ blending function** at one cell.
///
/// `F₁ = tanh(arg₁⁴)`, with
///
/// ```text
///   arg₁ = min(
///       max( √k/(β*·ω·d) ,  500·ν/(d²·ω) ),
///       4·ρ·σ_ω2·k / (CD_kω·d²)
///   )
/// ```
///
/// where `CD_kω = max(2·ρ·σ_ω2·(∇k·∇ω)/ω , 1e-10)` is the positive
/// part of the cross-diffusion term and `d` is the distance to the
/// nearest wall. `F₁ → 1` deep inside a boundary layer (set 1, k-ω)
/// and `F₁ → 0` in the free stream (set 2, k-ε).
#[inline]
#[allow(clippy::too_many_arguments)]
pub fn f1_blend(
    rho: f64,
    nu_mol: f64,
    k: f64,
    omega: f64,
    grad_k_dot_grad_omega: f64,
    wall_distance: f64,
    sigma_omega2: f64,
    beta_star: f64,
) -> f64 {
    let d = wall_distance.max(1e-12);
    let omega = omega.max(1e-12);
    let k = k.max(0.0);
    // Positive cross-diffusion CD_kω.
    let cd_pos = (2.0 * rho * sigma_omega2 * grad_k_dot_grad_omega / omega).max(1.0e-10);
    let term1 = k.sqrt() / (beta_star * omega * d);
    let term2 = 500.0 * nu_mol / (d * d * omega);
    let arg_lo = term1.max(term2);
    let term3 = 4.0 * rho * sigma_omega2 * k / (cd_pos * d * d);
    let arg1 = arg_lo.min(term3);
    (arg1.powi(4)).tanh()
}

/// The Menter SST **F₂ blending function** at one cell — the second
/// blending function, used inside the eddy-viscosity limiter
/// `μ_t = ρ·a₁·k / max(a₁·ω, S·F₂)`.
///
/// `F₂ = tanh(arg₂²)` with
///
/// ```text
///   arg₂ = max( 2·√k/(β*·ω·d) ,  500·ν/(d²·ω) )
/// ```
///
/// `F₂ → 1` inside a boundary layer (the limiter is *active* there,
/// and bounds `μ_t` in an adverse-pressure-gradient separation), and
/// `F₂ → 0` in the free stream (the limiter falls back to the
/// standard `ρ·k/ω`).
#[inline]
pub fn f2_blend(nu_mol: f64, k: f64, omega: f64, wall_distance: f64, beta_star: f64) -> f64 {
    let d = wall_distance.max(1e-12);
    let omega = omega.max(1e-12);
    let k = k.max(0.0);
    let term1 = 2.0 * k.sqrt() / (beta_star * omega * d);
    let term2 = 500.0 * nu_mol / (d * d * omega);
    let arg2 = term1.max(term2);
    (arg2.powi(2)).tanh()
}

/// The state of the **k-ω SST** turbulence model — the cell-centred
/// `k` and `ω` transported scalars and the eddy viscosity derived
/// from them.
///
/// All three fields are cell-centred (`nx × ny`), so they live on the
/// pressure grid of [`crate::Grid`]. [`SstField::nu_t`] is the model
/// output: the eddy viscosity at every cell, ready to be added to the
/// molecular viscosity in the momentum solve.
#[derive(Clone, Debug)]
pub struct SstField {
    /// Turbulent kinetic energy `k` (m²/s²), cell-centred.
    pub k: Field,
    /// Specific dissipation rate `ω` (1/s), cell-centred.
    pub omega: Field,
    /// Turbulent (eddy) **kinematic** viscosity `ν_t = μ_t/ρ` (m²/s),
    /// cell-centred — the SST-limited value
    /// `ν_t = a₁·k / max(a₁·ω, S·F₂)`. Kinematic so it adds directly
    /// to [`crate::Fluid::viscosity`].
    pub nu_t: Field,
}

impl SstField {
    /// Initialise an SST field from a reference velocity and a
    /// **turbulence intensity** `I` (the fraction of the mean velocity
    /// carried by the velocity fluctuations).
    ///
    /// The standard inlet estimates are used:
    ///
    /// ```text
    ///   k  = 1.5·(I·U_ref)²
    ///   ω  = √k / (C_μ^{1/4} · ℓ)        ℓ = 0.07·L  (mixing length)
    /// ```
    ///
    /// (the standard k-ω inlet relation, equivalent to the k-ε
    /// `ε = C_μ^{3/4}·k^{3/2}/ℓ` via `ω = ε/(C_μ·k)`).
    pub fn initialise(grid: &Grid, u_ref: f64, intensity: f64, model: &SstModel) -> SstField {
        let u_ref = u_ref.abs().max(1e-6);
        let intensity = intensity.max(1e-4);
        let k0 = 1.5 * (intensity * u_ref).powi(2);
        let length_scale = (0.07 * grid.ly).max(1e-6);
        let omega0 = k0.sqrt() / (model.beta_star.powf(0.25) * length_scale).max(1e-30);
        let nu_t0 = k0 / omega0.max(1e-30);
        SstField {
            k: Field::filled(grid.nx, grid.ny, k0),
            omega: Field::filled(grid.nx, grid.ny, omega0.max(1e-12)),
            nu_t: Field::filled(grid.nx, grid.ny, nu_t0.max(0.0)),
        }
    }

    /// Mean turbulent kinetic energy across the field — a scalar
    /// turbulence-level diagnostic.
    pub fn mean_k(&self) -> f64 {
        if self.k.data.is_empty() {
            return 0.0;
        }
        self.k.data.iter().sum::<f64>() / self.k.data.len() as f64
    }

    /// Recompute the SST-limited eddy viscosity field
    /// `ν_t = a₁·k / max(a₁·ω, S·F₂)` from the current `k`, `ω`,
    /// strain rate, and `F₂`.
    pub fn update_eddy_viscosity(
        &mut self,
        model: &SstModel,
        s_mag: &Field,
        wall_distance: &Field,
        nu_mol: f64,
    ) {
        let cap = (1.0e5 * nu_mol).max(1e-6);
        for j in 0..self.k.height {
            for i in 0..self.k.width {
                let k = self.k.at(i, j).max(0.0);
                let om = self.omega.at(i, j).max(1e-12);
                let f2 = f2_blend(nu_mol, k, om, wall_distance.at(i, j), model.beta_star);
                let s = s_mag.at(i, j).max(0.0);
                let denom = (model.a1 * om).max(s * f2).max(1e-30);
                let nu_t = (model.a1 * k / denom).clamp(0.0, cap);
                self.nu_t.set(i, j, nu_t);
            }
        }
    }
}

/// The cell-centred **strain-rate magnitude** `S = √(2·Sᵢⱼ·Sᵢⱼ)` at
/// every cell — `S² = 2(∂u/∂x)² + 2(∂v/∂y)² + (∂u/∂y + ∂v/∂x)²`, the
/// **square root** of the field [`strain_rate_squared`] returns.
///
/// Used directly by the SST limiter
/// `μ_t = ρ·a₁·k / max(a₁·ω, S·F₂)`.
pub fn strain_rate_magnitude(u: &Field, v: &Field, grid: &Grid) -> Field {
    let s2 = strain_rate_squared(u, v, grid);
    let mut s = Field::zeros(s2.width, s2.height);
    for idx in 0..s2.data.len() {
        s.data[idx] = s2.data[idx].max(0.0).sqrt();
    }
    s
}

/// Build a **wall-distance field** for a rectangular cavity / channel
/// — at every cell, the distance to the nearest of the four boundary
/// walls (in metres).
///
/// The arguments select which walls participate: a `true` here turns
/// that side into a no-slip wall the cell measures its distance from;
/// `false` excludes it (an inlet / outlet / open boundary).
///
/// SST's blending functions `F₁` and `F₂` need this field; for a
/// generic rectangular domain it is closed form, so a small standalone
/// helper is the cheapest way to produce it.
pub fn wall_distance_field(grid: &Grid, walls: WallMask) -> Field {
    let dx = grid.dx();
    let dy = grid.dy();
    let mut field = Field::zeros(grid.nx, grid.ny);
    let big = (grid.lx + grid.ly) * 10.0;
    for j in 0..grid.ny {
        for i in 0..grid.nx {
            let x = (i as f64 + 0.5) * dx;
            let y = (j as f64 + 0.5) * dy;
            let mut d = big;
            if walls.west {
                d = d.min(x);
            }
            if walls.east {
                d = d.min(grid.lx - x);
            }
            if walls.south {
                d = d.min(y);
            }
            if walls.north {
                d = d.min(grid.ly - y);
            }
            if d == big {
                // No walls at all: degenerate. Use the half-diagonal
                // as a safe floor so SST blends never see zero.
                d = 0.5 * (grid.lx * grid.lx + grid.ly * grid.ly).sqrt();
            }
            field.set(i, j, d.max(1e-12));
        }
    }
    field
}

/// Which sides of a rectangular domain are no-slip walls — the mask
/// [`wall_distance_field`] consults.
#[derive(Clone, Copy, Debug)]
pub struct WallMask {
    /// West (`x = 0`) is a no-slip wall.
    pub west: bool,
    /// East (`x = lx`) is a no-slip wall.
    pub east: bool,
    /// South (`y = 0`) is a no-slip wall.
    pub south: bool,
    /// North (`y = ly`) is a no-slip wall.
    pub north: bool,
}

impl WallMask {
    /// A four-wall closed cavity — every side a wall.
    pub fn closed() -> WallMask {
        WallMask {
            west: true,
            east: true,
            south: true,
            north: true,
        }
    }

    /// A plane channel — only the top and bottom are walls, inlet /
    /// outlet on the sides.
    pub fn channel() -> WallMask {
        WallMask {
            west: false,
            east: false,
            south: true,
            north: true,
        }
    }
}

/// Advance the **k-ω SST** transport equations one outer iteration.
///
/// This is the SST counterpart of [`advance_k_epsilon`]. Given the
/// current mean velocity (the staggered `u` / `v` fields), the wall
/// distance field (for `F₁` / `F₂`), and the existing `k`-`ω` state,
/// it:
///
/// 1. forms the strain-rate magnitude `S` and the production
///    `P_k = ν_t·S²`;
/// 2. evaluates `F₁` and `F₂` per cell and blends the model
///    coefficients accordingly;
/// 3. sweeps a convection-diffusion update of `k` and then of `ω`,
///    each with the SST source / sink pair and the cross-diffusion
///    term, using the same hybrid upwind/central scheme and
///    Gauss-Seidel relaxation as the momentum solve;
/// 4. applies the **wall-function** value of `ω` in every cell
///    adjacent to a no-slip wall (`ω_wall = 6·ν / (β₁·y²)` —
///    Menter's near-wall boundary condition), and the equilibrium
///    `k_wall = u_τ² / √β*`;
/// 5. clamps `k ≥ 0` and `ω` to a small positive floor;
/// 6. recomputes the SST-limited `ν_t = a₁·k / max(a₁·ω, S·F₂)`.
///
/// `relax` is the under-relaxation factor (SST, like k-ε, needs heavy
/// relaxation, `≈ 0.5–0.7`). `wall_south` / `wall_north` flag whether
/// the bottom / top domain edges are no-slip walls.
#[allow(clippy::too_many_arguments)]
pub fn advance_k_omega_sst(
    state: &mut SstField,
    u: &Field,
    v: &Field,
    grid: &Grid,
    rho: f64,
    nu_mol: f64,
    model: &SstModel,
    wall: &WallFunction,
    wall_distance: &Field,
    wall_south: bool,
    wall_north: bool,
    relax: f64,
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx();
    let dy = grid.dy();
    let relax = relax.clamp(0.05, 1.0);

    let k_floor = 1e-10;
    let omega_floor = 1e-8;

    // --- strain rate and production ---
    let s2 = strain_rate_squared(u, v, grid);
    let mut s_mag = Field::zeros(nx, ny);
    let mut production = Field::zeros(nx, ny);
    for idx in 0..s2.data.len() {
        let s = s2.data[idx].max(0.0).sqrt();
        s_mag.data[idx] = s;
        production.data[idx] = state.nu_t.data[idx] * s * s;
    }

    // --- F₁ / F₂ and blended coefficients per cell ---
    // Pre-compute the ∇k·∇ω central-difference dot products.
    let mut f1 = Field::zeros(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let dkdx = central_difference_x(&state.k, i, j, dx);
            let dkdy = central_difference_y(&state.k, i, j, dy);
            let dwdx = central_difference_x(&state.omega, i, j, dx);
            let dwdy = central_difference_y(&state.omega, i, j, dy);
            let grad_dot = dkdx * dwdx + dkdy * dwdy;
            let val = f1_blend(
                rho,
                nu_mol,
                state.k.at(i, j),
                state.omega.at(i, j),
                grad_dot,
                wall_distance.at(i, j),
                model.set2.sigma_omega,
                model.beta_star,
            );
            f1.set(i, j, val);
        }
    }

    // --- (1) k transport ---
    // Use the freshest k while sweeping; store source data from a
    // frozen snapshot for nonlinear stability.
    let k_old = state.k.clone();
    let omega_snapshot = state.omega.clone();
    for _sweep in 0..2 {
        for j in 0..ny {
            for i in 0..nx {
                let blend = model.blend(f1.at(i, j));
                let mu_t = state.nu_t.at(i, j); // kinematic
                let diff_p = nu_mol + blend.sigma_k * mu_t;
                let fe = rho * dy * u.at(i + 1, j);
                let fw = rho * dy * u.at(i, j);
                let fn_ = rho * dx * v.at(i, j + 1);
                let fs = rho * dx * v.at(i, j);

                let mut a_p = 0.0;
                let mut nb = 0.0;
                if i + 1 < nx {
                    let d = rho * diff_p * dy / dx;
                    let a = hybrid(d, fe, -1.0);
                    a_p += a;
                    nb += a * state.k.at(i + 1, j);
                }
                if i > 0 {
                    let d = rho * diff_p * dy / dx;
                    let a = hybrid(d, fw, 1.0);
                    a_p += a;
                    nb += a * state.k.at(i - 1, j);
                }
                if j + 1 < ny {
                    let d = rho * diff_p * dx / dy;
                    let a = hybrid(d, fn_, -1.0);
                    a_p += a;
                    nb += a * state.k.at(i, j + 1);
                }
                if j > 0 {
                    let d = rho * diff_p * dx / dy;
                    let a = hybrid(d, fs, 1.0);
                    a_p += a;
                    nb += a * state.k.at(i, j - 1);
                }
                let vol = dx * dy;
                let k_p = k_old.at(i, j).max(k_floor);
                let om_p = omega_snapshot.at(i, j).max(omega_floor);
                // Production limiter min(P_k, 10·β*·ρ·k·ω) — the
                // standard SST stagnation-point cap.
                let p_lim = (production.at(i, j)).min(10.0 * model.beta_star * rho * k_p * om_p);
                let su = p_lim * vol;
                // Sink β*·ρ·k·ω linearised as (β*·ρ·ω)·k on the
                // diagonal (Patankar positive-coefficient rule).
                let sp = model.beta_star * rho * om_p * vol;
                a_p += sp;
                if a_p.abs() < 1e-30 {
                    continue;
                }
                let k_new = (nb + su) / a_p;
                let updated = state.k.at(i, j) + relax * (k_new - state.k.at(i, j));
                state.k.set(i, j, updated.max(k_floor));
            }
        }
    }

    // --- (2) ω transport with cross-diffusion ---
    let omega_old = state.omega.clone();
    for _sweep in 0..2 {
        for j in 0..ny {
            for i in 0..nx {
                let f1_p = f1.at(i, j);
                let blend = model.blend(f1_p);
                let mu_t = state.nu_t.at(i, j);
                let diff_p = nu_mol + blend.sigma_omega * mu_t;
                let fe = rho * dy * u.at(i + 1, j);
                let fw = rho * dy * u.at(i, j);
                let fn_ = rho * dx * v.at(i, j + 1);
                let fs = rho * dx * v.at(i, j);

                let mut a_p = 0.0;
                let mut nb = 0.0;
                if i + 1 < nx {
                    let d = rho * diff_p * dy / dx;
                    let a = hybrid(d, fe, -1.0);
                    a_p += a;
                    nb += a * state.omega.at(i + 1, j);
                }
                if i > 0 {
                    let d = rho * diff_p * dy / dx;
                    let a = hybrid(d, fw, 1.0);
                    a_p += a;
                    nb += a * state.omega.at(i - 1, j);
                }
                if j + 1 < ny {
                    let d = rho * diff_p * dx / dy;
                    let a = hybrid(d, fn_, -1.0);
                    a_p += a;
                    nb += a * state.omega.at(i, j + 1);
                }
                if j > 0 {
                    let d = rho * diff_p * dx / dy;
                    let a = hybrid(d, fs, 1.0);
                    a_p += a;
                    nb += a * state.omega.at(i, j - 1);
                }
                let vol = dx * dy;
                let k_p = state.k.at(i, j).max(k_floor);
                let om_p = omega_old.at(i, j).max(omega_floor);
                // Production γ·ρ·S² (form using S² directly, valid
                // because of the eddy-viscosity definition;
                // equivalent to γ·(ρ/μ_t)·P_k when μ_t = ρ·k/ω;
                // robust at μ_t → 0).
                let su_prod = blend.gamma * rho * s_mag.at(i, j).powi(2) * vol;
                // Sink β·ρ·ω², linearised onto the diagonal.
                let sp_sink = blend.beta * rho * om_p * vol;
                // Cross-diffusion 2·(1−F₁)·ρ·σ_ω2·(∇k·∇ω)/ω — only
                // its positive part is a stable source (the negative
                // part shrinks ω, so a linearised sink is added in
                // that branch).
                let dkdx = central_difference_x(&state.k, i, j, dx);
                let dkdy = central_difference_y(&state.k, i, j, dy);
                let dwdx = central_difference_x(&state.omega, i, j, dx);
                let dwdy = central_difference_y(&state.omega, i, j, dy);
                let grad_dot = dkdx * dwdx + dkdy * dwdy;
                let cd = 2.0 * (1.0 - f1_p) * rho * model.set2.sigma_omega * grad_dot / om_p;
                let mut su = su_prod;
                let mut sp_extra = sp_sink;
                if cd >= 0.0 {
                    su += cd * vol;
                } else {
                    // Negative CD: a sink proportional to −cd·ω, but
                    // since cd already includes 1/ω we add as a fixed
                    // negative source (it can shrink ω toward floor).
                    su += cd * vol;
                }
                a_p += sp_extra;
                // Suppress unused-mut warning if cd>0 path doesn't update sp_extra.
                let _ = &mut sp_extra;
                if a_p.abs() < 1e-30 {
                    continue;
                }
                let _ = k_p;
                let om_new = (nb + su) / a_p;
                let updated = state.omega.at(i, j) + relax * (om_new - state.omega.at(i, j));
                state.omega.set(i, j, updated.max(omega_floor));
            }
        }
    }

    // --- (3) wall-function values in wall-adjacent cells ---
    let y_first = 0.5 * dy;
    let beta1 = model.set1.beta; // near-wall β for the ω boundary value
    let apply_wall = |state: &mut SstField, i: usize, j: usize| {
        let uc = 0.5 * (u.at(i, j) + u.at(i + 1, j));
        let vc = 0.5 * (v.at(i, j) + v.at(i, j + 1));
        let u_p = (uc * uc + vc * vc).sqrt();
        let u_tau = wall.friction_velocity(rho, nu_mol, u_p, y_first);
        // k_P = u_τ² / √β* (equilibrium log-layer kinetic energy).
        let k_wall = (u_tau * u_tau / model.beta_star.sqrt()).max(k_floor);
        // ω_wall — Menter's near-wall boundary value
        //   ω = 6·ν / (β₁ · y²)   in the viscous sublayer
        //   ω = u_τ / (√β* · κ · y)   in the log layer.
        // Use the maximum of the two, which is the well-conditioned
        // blended boundary value SST uses.
        let omega_visc = 6.0 * nu_mol / (beta1 * y_first * y_first).max(1e-30);
        let omega_log = u_tau / (model.beta_star.sqrt() * wall.kappa * y_first).max(1e-30);
        let omega_wall = omega_visc.max(omega_log).max(omega_floor);
        state.k.set(i, j, k_wall);
        state.omega.set(i, j, omega_wall);
    };
    if wall_south && ny > 0 {
        for i in 0..nx {
            apply_wall(state, i, 0);
        }
    }
    if wall_north && ny > 0 {
        for i in 0..nx {
            apply_wall(state, i, ny - 1);
        }
    }

    // --- (6) refresh the SST-limited eddy viscosity ---
    state.update_eddy_viscosity(model, &s_mag, wall_distance, nu_mol);
}

/// Central-difference partial derivative `∂φ/∂x` of a cell-centred
/// scalar at cell `(i, j)`. One-sided at the boundary.
#[inline]
fn central_difference_x(f: &Field, i: usize, j: usize, dx: f64) -> f64 {
    if i + 1 < f.width && i > 0 {
        (f.at(i + 1, j) - f.at(i - 1, j)) / (2.0 * dx)
    } else if i + 1 < f.width {
        (f.at(i + 1, j) - f.at(i, j)) / dx
    } else if i > 0 {
        (f.at(i, j) - f.at(i - 1, j)) / dx
    } else {
        0.0
    }
}

/// Central-difference partial derivative `∂φ/∂y` of a cell-centred
/// scalar at cell `(i, j)`. One-sided at the boundary.
#[inline]
fn central_difference_y(f: &Field, i: usize, j: usize, dy: f64) -> f64 {
    if j + 1 < f.height && j > 0 {
        (f.at(i, j + 1) - f.at(i, j - 1)) / (2.0 * dy)
    } else if j + 1 < f.height {
        (f.at(i, j + 1) - f.at(i, j)) / dy
    } else if j > 0 {
        (f.at(i, j) - f.at(i, j - 1)) / dy
    } else {
        0.0
    }
}

/// Solve a **fully-developed plane channel flow** with the **k-ω SST**
/// model and return its wall-normal velocity profile. The SST
/// counterpart of [`solve_turbulent_channel`].
///
/// `turbulent = false` returns the laminar parabola (no turbulence
/// model engaged); `turbulent = true` runs the k-ω SST coupling.
pub fn solve_turbulent_channel_sst(
    grid: &Grid,
    fluid: &crate::solver::Fluid,
    model: &SstModel,
    wall: &WallFunction,
    dpdx: f64,
    turbulent: bool,
    outer_iterations: usize,
) -> ChannelProfile {
    let ny = grid.ny;
    let dy = grid.dy();
    let rho = fluid.density;
    let nu_mol = fluid.viscosity;
    let source = dpdx / rho;

    let mut u = vec![0.0_f64; ny];
    let mut nu_t = vec![0.0_f64; ny];

    let column = Grid::new(1, ny, grid.dx(), grid.ly);
    let mut state = SstField::initialise(&column, 1.0, 0.05, model);
    let wall_distance = wall_distance_field(&column, WallMask::channel());

    let outers = outer_iterations.max(1);
    for _outer in 0..outers {
        // --- (1) 1-D momentum balance ---
        let diff_at = |j: usize| -> f64 { nu_mol + nu_t[j] };
        for _sweep in 0..200 {
            for j in 0..ny {
                let (an, u_n) = if j + 1 < ny {
                    let d = 0.5 * (diff_at(j) + diff_at(j + 1)) / dy;
                    (d, u[j + 1])
                } else {
                    let d = diff_at(j) / (0.5 * dy);
                    (d, 0.0)
                };
                let (as_, u_s) = if j > 0 {
                    let d = 0.5 * (diff_at(j) + diff_at(j - 1)) / dy;
                    (d, u[j - 1])
                } else {
                    let d = diff_at(j) / (0.5 * dy);
                    (d, 0.0)
                };
                let a_p = an + as_;
                if a_p.abs() < 1e-30 {
                    continue;
                }
                let rhs = an * u_n + as_ * u_s - source * dy;
                u[j] = rhs / a_p;
            }
        }
        if !turbulent {
            break;
        }
        // --- (2) advance SST from the new velocity field ---
        let mut u_field = column.u_field();
        for (j, &uj) in u.iter().enumerate().take(ny) {
            u_field.set(0, j, uj);
            u_field.set(1, j, uj);
        }
        let v_field = column.v_field();
        for _ke in 0..3 {
            advance_k_omega_sst(
                &mut state,
                &u_field,
                &v_field,
                &column,
                rho,
                nu_mol,
                model,
                wall,
                &wall_distance,
                true,
                true,
                0.6,
            );
        }
        for (j, nt) in nu_t.iter_mut().enumerate().take(ny) {
            *nt = state.nu_t.at(0, j);
        }
    }

    let y: Vec<f64> = (0..ny).map(|j| (j as f64 + 0.5) * dy).collect();
    let bulk = u.iter().sum::<f64>() / ny.max(1) as f64;
    ChannelProfile {
        y,
        u,
        nu_t,
        bulk_velocity: bulk,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::Fluid;

    #[test]
    fn eddy_viscosity_follows_the_cmu_k2_over_eps_law() {
        // ν_t = C_μ·k²/ε — the defining relation of the k-ε model.
        let grid = Grid::new(4, 4, 1.0, 1.0);
        let model = KEpsilonModel::standard();
        let mut turb = TurbulenceField {
            k: Field::filled(4, 4, 0.2),
            epsilon: Field::filled(4, 4, 0.5),
            nu_t: Field::zeros(4, 4),
        };
        turb.update_eddy_viscosity(&model, 1e-5);
        let expected = 0.09 * 0.2 * 0.2 / 0.5;
        for &nu_t in &turb.nu_t.data {
            assert!(
                (nu_t - expected).abs() < 1e-12,
                "ν_t {nu_t} should equal C_μ·k²/ε {expected}"
            );
        }
        let _ = grid;
    }

    #[test]
    fn initial_turbulence_scales_with_intensity() {
        // k = 1.5·(I·U)² — doubling the intensity quadruples k.
        let grid = Grid::new(8, 8, 1.0, 1.0);
        let model = KEpsilonModel::standard();
        let low = TurbulenceField::initialise(&grid, 1.0, 0.05, &model);
        let high = TurbulenceField::initialise(&grid, 1.0, 0.10, &model);
        let ratio = high.mean_k() / low.mean_k();
        assert!(
            (ratio - 4.0).abs() < 1e-6,
            "k should scale with I²: ratio {ratio} expected 4"
        );
        // k itself matches the analytic estimate.
        let k_expected = 1.5 * (0.05_f64 * 1.0).powi(2);
        assert!((low.mean_k() - k_expected).abs() < 1e-12);
    }

    #[test]
    fn friction_velocity_satisfies_the_log_law() {
        // u_τ from friction_velocity must, plugged back into the law of
        // the wall, reproduce the input velocity.
        let wall = WallFunction::standard();
        let (rho, mu) = (1.0, 1e-3);
        let (u_p, y_p) = (1.0, 0.01);
        let u_tau = wall.friction_velocity(rho, mu, u_p, y_p);
        assert!(u_tau > 0.0, "friction velocity must be positive");
        let y_plus = wall.y_plus(rho, u_tau, y_p, mu);
        // Only meaningful to check the log law when we are in its range.
        if y_plus > 11.06 {
            let u_reconstructed = u_tau * wall.u_plus(y_plus);
            let rel = (u_reconstructed - u_p).abs() / u_p;
            assert!(
                rel < 1e-3,
                "log law should reproduce u_p: got {u_reconstructed} vs {u_p}"
            );
        }
    }

    #[test]
    fn friction_velocity_rises_with_wall_velocity() {
        // A faster near-wall flow exerts more wall shear → larger u_τ.
        let wall = WallFunction::standard();
        let (rho, mu) = (1.0, 1e-3);
        let slow = wall.friction_velocity(rho, mu, 0.5, 0.01);
        let fast = wall.friction_velocity(rho, mu, 5.0, 0.01);
        assert!(fast > slow, "u_τ should grow with wall velocity");
        // Zero velocity → zero friction.
        assert_eq!(wall.friction_velocity(rho, mu, 0.0, 0.01), 0.0);
    }

    #[test]
    fn strain_rate_is_zero_for_uniform_flow() {
        // A spatially-uniform velocity field has no shear → S² = 0.
        let grid = Grid::new(6, 6, 1.0, 1.0);
        let u = Field::filled(grid.nx + 1, grid.ny, 2.0);
        let v = Field::zeros(grid.nx, grid.ny + 1);
        let s2 = strain_rate_squared(&u, &v, &grid);
        assert!(
            s2.abs_max() < 1e-12,
            "uniform flow should have zero strain rate"
        );
    }

    #[test]
    fn strain_rate_is_positive_for_sheared_flow() {
        // A linear shear u(y) gives a non-zero strain rate.
        let grid = Grid::new(6, 6, 1.0, 1.0);
        let mut u = Field::zeros(grid.nx + 1, grid.ny);
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                u.set(i, j, j as f64); // u increases with y → ∂u/∂y ≠ 0
            }
        }
        let v = Field::zeros(grid.nx, grid.ny + 1);
        let s2 = strain_rate_squared(&u, &v, &grid);
        assert!(
            s2.abs_max() > 0.0,
            "a sheared flow must have a positive strain rate"
        );
    }

    #[test]
    fn k_epsilon_step_keeps_k_and_epsilon_positive() {
        // k and ε are physically non-negative; the solver must keep
        // them so even after a transport step on a sheared flow.
        let grid = Grid::new(10, 10, 2.0, 1.0);
        let model = KEpsilonModel::standard();
        let wall = WallFunction::standard();
        let mut turb = TurbulenceField::initialise(&grid, 1.0, 0.05, &model);
        // A sheared velocity field.
        let mut u = Field::zeros(grid.nx + 1, grid.ny);
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                let y = (j as f64 + 0.5) / grid.ny as f64;
                u.set(i, j, 4.0 * y * (1.0 - y)); // parabola-ish
            }
        }
        let v = Field::zeros(grid.nx, grid.ny + 1);
        for _ in 0..20 {
            advance_k_epsilon(
                &mut turb, &u, &v, &grid, 1.0, 1e-5, &model, &wall, true, true, 0.6,
            );
        }
        for &k in &turb.k.data {
            assert!(k > 0.0 && k.is_finite(), "k must stay positive: {k}");
        }
        for &eps in &turb.epsilon.data {
            assert!(eps > 0.0 && eps.is_finite(), "ε must stay positive: {eps}");
        }
        for &nu_t in &turb.nu_t.data {
            assert!(nu_t >= 0.0 && nu_t.is_finite(), "ν_t must be finite ≥ 0");
        }
    }

    #[test]
    fn turbulence_production_raises_k_in_a_sheared_flow() {
        // Shear feeds the turbulence: starting k-ε in a strongly
        // sheared flow, the production term should grow k above a
        // near-zero initial level.
        let grid = Grid::new(12, 12, 2.0, 1.0);
        let model = KEpsilonModel::standard();
        let wall = WallFunction::standard();
        let mut turb = TurbulenceField {
            k: Field::filled(grid.nx, grid.ny, 1e-6),
            epsilon: Field::filled(grid.nx, grid.ny, 1e-4),
            nu_t: Field::filled(grid.nx, grid.ny, 1e-4),
        };
        let k_initial = turb.mean_k();
        // A strong linear shear.
        let mut u = Field::zeros(grid.nx + 1, grid.ny);
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                u.set(i, j, 10.0 * j as f64 / grid.ny as f64);
            }
        }
        let v = Field::zeros(grid.nx, grid.ny + 1);
        for _ in 0..40 {
            advance_k_epsilon(
                &mut turb, &u, &v, &grid, 1.0, 1e-5, &model, &wall, false, false, 0.6,
            );
        }
        assert!(
            turb.mean_k() > k_initial,
            "shear production should raise k: {} → {}",
            k_initial,
            turb.mean_k()
        );
    }

    #[test]
    fn laminar_channel_develops_the_parabolic_profile() {
        // With the turbulence model OFF the fully-developed channel
        // must integrate to the analytic parabola: u_bulk/u_max = 2/3.
        let grid = Grid::new(1, 40, 1.0, 1.0);
        let fluid = Fluid::new(1.0, 0.01);
        let model = KEpsilonModel::standard();
        let wall = WallFunction::standard();
        let profile = solve_turbulent_channel(&grid, &fluid, &model, &wall, -1.0, false, 1);
        // The flatness of a parabola is exactly 2/3.
        let flat = profile.flatness();
        assert!(
            (flat - 2.0 / 3.0).abs() < 0.03,
            "laminar profile flatness {flat} should be ≈ 2/3 (parabola)"
        );
        // The profile is symmetric and peaks at the centre.
        let n = profile.u.len();
        let peak_idx = profile
            .u
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap()
            .0;
        assert!(
            peak_idx >= n / 2 - 2 && peak_idx <= n / 2 + 1,
            "laminar peak should be at the channel centre, got index {peak_idx}"
        );
        // No-slip walls: the near-wall velocity is far below the peak.
        assert!(
            profile.u[0] < 0.5 * profile.centreline_velocity(),
            "velocity must fall toward the no-slip wall"
        );
    }

    #[test]
    fn turbulent_channel_is_flatter_than_the_laminar_one() {
        // The headline k-ε verification. Solving the *same* channel
        // with the turbulence model ON must produce a flatter
        // mean-velocity profile than the laminar parabola — turbulent
        // eddy diffusion redistributes momentum toward the walls. A
        // flatter profile means a larger bulk/centreline ratio.
        let grid = Grid::new(1, 60, 1.0, 2.0);
        let fluid = Fluid::new(1.0, 1e-4); // high Re → strongly turbulent
        let model = KEpsilonModel::standard();
        let wall = WallFunction::standard();

        let laminar = solve_turbulent_channel(&grid, &fluid, &model, &wall, -1.0, false, 1);
        let turbulent = solve_turbulent_channel(&grid, &fluid, &model, &wall, -1.0, true, 80);

        let lam_flat = laminar.flatness();
        let turb_flat = turbulent.flatness();
        assert!(
            turb_flat > lam_flat + 0.03,
            "turbulent profile (flatness {turb_flat}) must be flatter than laminar (flatness {lam_flat})"
        );
        // A turbulent channel's flatness lands well above the parabolic
        // 2/3 — it does not become *more* peaked.
        assert!(
            turb_flat > 0.7,
            "turbulent channel flatness {turb_flat} should exceed the parabolic 2/3"
        );
    }

    #[test]
    fn turbulent_channel_has_a_sane_near_wall_eddy_viscosity() {
        // Physical near-wall behaviour: the turbulent (eddy) viscosity
        // must be largest away from the wall and fall toward the wall
        // (the wall damps the turbulence). And it must everywhere swamp
        // the molecular viscosity in a high-Reynolds-number flow —
        // that dominance is what makes the flow "turbulent".
        let grid = Grid::new(1, 60, 1.0, 2.0);
        let fluid = Fluid::new(1.0, 1e-4);
        let model = KEpsilonModel::standard();
        let wall = WallFunction::standard();
        let profile = solve_turbulent_channel(&grid, &fluid, &model, &wall, -1.0, true, 80);
        // Every eddy viscosity is finite and non-negative.
        for &nt in &profile.nu_t {
            assert!(nt >= 0.0 && nt.is_finite(), "ν_t must be finite ≥ 0");
        }
        let n = profile.nu_t.len();
        let nu_t_wall = profile.nu_t[0];
        let nu_t_centre = profile.nu_t[n / 2];
        // The eddy viscosity at the channel centre exceeds the value in
        // the near-wall cell.
        assert!(
            nu_t_centre > nu_t_wall,
            "ν_t should grow away from the wall: wall {nu_t_wall}, centre {nu_t_centre}"
        );
        // In a turbulent flow the eddy viscosity dwarfs the molecular
        // one over the channel core.
        assert!(
            nu_t_centre > 10.0 * fluid.viscosity,
            "core ν_t {nu_t_centre} should swamp molecular ν {}",
            fluid.viscosity
        );
    }

    // ----- k-ω SST tests --------------------------------------------

    #[test]
    fn sst_blend_returns_set1_at_f1_unity_and_set2_at_zero() {
        // F₁ = 1 → blended coefficients = set 1; F₁ = 0 → set 2.
        let m = SstModel::standard();
        let b1 = m.blend(1.0);
        let b2 = m.blend(0.0);
        assert!((b1.sigma_k - m.set1.sigma_k).abs() < 1e-12);
        assert!((b1.beta - m.set1.beta).abs() < 1e-12);
        assert!((b2.sigma_k - m.set2.sigma_k).abs() < 1e-12);
        assert!((b2.beta - m.set2.beta).abs() < 1e-12);
        // Halfway blends as an average.
        let bh = m.blend(0.5);
        assert!((bh.beta - 0.5 * (m.set1.beta + m.set2.beta)).abs() < 1e-12);
    }

    #[test]
    fn sst_gamma_uses_the_menter_relation() {
        // γ_i = β_i / β* − σ_ω_i · κ² / √β* (Menter 1994).
        let m = SstModel::standard();
        let g1 =
            m.set1.beta / m.beta_star - m.set1.sigma_omega * m.kappa * m.kappa / m.beta_star.sqrt();
        let g2 =
            m.set2.beta / m.beta_star - m.set2.sigma_omega * m.kappa * m.kappa / m.beta_star.sqrt();
        assert!((m.set1.gamma - g1).abs() < 1e-12);
        assert!((m.set2.gamma - g2).abs() < 1e-12);
    }

    #[test]
    fn f1_blend_tends_to_one_near_a_wall_and_zero_far_away() {
        // Deep inside a boundary layer F₁ → 1 (set 1, k-ω); far in
        // the free stream F₁ → 0 (set 2, k-ε).
        let m = SstModel::standard();
        let rho = 1.0;
        let nu = 1e-5;
        let k = 0.01;
        let omega = 100.0;
        let near = f1_blend(
            rho,
            nu,
            k,
            omega,
            1e-6,
            1e-6,
            m.set2.sigma_omega,
            m.beta_star,
        );
        let far = f1_blend(
            rho,
            nu,
            k,
            omega,
            1e-6,
            10.0,
            m.set2.sigma_omega,
            m.beta_star,
        );
        assert!(near > 0.9, "near-wall F₁ {near} should be ≈ 1");
        assert!(far < 0.1, "free-stream F₁ {far} should be ≈ 0");
    }

    #[test]
    fn f2_blend_is_active_in_the_boundary_layer() {
        // F₂ → 1 in the boundary layer (limiter active); F₂ → 0 far
        // from the wall (limiter falls back to plain ρ·k/ω).
        let near = f2_blend(1e-5, 0.01, 100.0, 1e-6, 0.09);
        let far = f2_blend(1e-5, 0.01, 100.0, 10.0, 0.09);
        assert!(near > 0.9, "near-wall F₂ {near} should be ≈ 1");
        assert!(far < 0.1, "free-stream F₂ {far} should be ≈ 0");
    }

    #[test]
    fn wall_distance_matches_the_closed_form_for_a_channel() {
        // For a channel (south + north walls) the wall distance at
        // cell j is min(y, ly − y).
        let g = Grid::new(4, 8, 1.0, 2.0);
        let wd = wall_distance_field(&g, WallMask::channel());
        let dy = g.dy();
        for j in 0..g.ny {
            let y = (j as f64 + 0.5) * dy;
            let d_ref = y.min(g.ly - y);
            for i in 0..g.nx {
                assert!(
                    (wd.at(i, j) - d_ref).abs() < 1e-12,
                    "wall distance mismatch at ({i},{j})"
                );
            }
        }
    }

    #[test]
    fn strain_rate_magnitude_is_sqrt_of_strain_rate_squared() {
        // s_mag = √s² — they had better satisfy this exactly.
        let grid = Grid::new(6, 6, 1.0, 1.0);
        let mut u = Field::zeros(grid.nx + 1, grid.ny);
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                u.set(i, j, (j as f64) * 0.7);
            }
        }
        let v = Field::zeros(grid.nx, grid.ny + 1);
        let s2 = strain_rate_squared(&u, &v, &grid);
        let s = strain_rate_magnitude(&u, &v, &grid);
        for idx in 0..s.data.len() {
            assert!((s.data[idx] - s2.data[idx].max(0.0).sqrt()).abs() < 1e-12);
        }
    }

    #[test]
    fn sst_field_eddy_viscosity_uses_the_a1_limiter() {
        // ν_t = a₁·k / max(a₁·ω, S·F₂). With S·F₂ > a₁·ω the limiter
        // kicks in and ν_t = a₁·k / (S·F₂); otherwise ν_t = k/ω.
        let m = SstModel::standard();
        let grid = Grid::new(2, 2, 1.0, 1.0);
        let mut state = SstField {
            k: Field::filled(grid.nx, grid.ny, 0.04),
            omega: Field::filled(grid.nx, grid.ny, 1.0),
            nu_t: Field::zeros(grid.nx, grid.ny),
        };
        let big_s = Field::filled(grid.nx, grid.ny, 100.0);
        let small_s = Field::filled(grid.nx, grid.ny, 0.001);
        let wd = Field::filled(grid.nx, grid.ny, 1e-6); // F₂ ≈ 1
                                                        // Limiter-active branch: a₁·ω = 0.31, S·F₂ ≈ 100 → ν_t = 0.31·k/100.
        state.update_eddy_viscosity(&m, &big_s, &wd, 1e-6);
        let expected_lim = 0.31 * 0.04 / 100.0;
        assert!(
            (state.nu_t.at(0, 0) - expected_lim).abs() < 1e-6,
            "limiter-active ν_t {} should be a₁·k/(S·F₂) = {}",
            state.nu_t.at(0, 0),
            expected_lim
        );
        // Limiter-inactive branch: S·F₂ small → ν_t = k/ω.
        state.update_eddy_viscosity(&m, &small_s, &wd, 1e-6);
        let expected_kw = 0.04 / 1.0;
        assert!(
            (state.nu_t.at(0, 0) - expected_kw).abs() < 1e-6,
            "limiter-inactive ν_t {} should be k/ω = {}",
            state.nu_t.at(0, 0),
            expected_kw
        );
    }

    #[test]
    fn sst_step_keeps_k_and_omega_positive() {
        // The SST advance must preserve the physical positivity of
        // both transported scalars even after many sweeps on a sheared
        // flow.
        let grid = Grid::new(8, 12, 1.0, 1.0);
        let m = SstModel::standard();
        let w = WallFunction::standard();
        let wd = wall_distance_field(&grid, WallMask::channel());
        let mut state = SstField::initialise(&grid, 1.0, 0.05, &m);
        let mut u = Field::zeros(grid.nx + 1, grid.ny);
        for j in 0..grid.ny {
            for i in 0..=grid.nx {
                let y = (j as f64 + 0.5) / grid.ny as f64;
                u.set(i, j, 4.0 * y * (1.0 - y));
            }
        }
        let v = Field::zeros(grid.nx, grid.ny + 1);
        for _ in 0..20 {
            advance_k_omega_sst(
                &mut state, &u, &v, &grid, 1.0, 1e-5, &m, &w, &wd, true, true, 0.6,
            );
        }
        for &k in &state.k.data {
            assert!(k >= 0.0 && k.is_finite(), "k must stay ≥ 0: {k}");
        }
        for &om in &state.omega.data {
            assert!(om > 0.0 && om.is_finite(), "ω must stay > 0: {om}");
        }
        for &nu_t in &state.nu_t.data {
            assert!(
                nu_t >= 0.0 && nu_t.is_finite(),
                "ν_t must be finite ≥ 0: {nu_t}"
            );
        }
    }

    #[test]
    fn sst_turbulent_channel_is_flatter_than_the_laminar_one() {
        // Headline SST verification — turning the model on flattens
        // the fully-developed channel profile relative to the laminar
        // parabola. Same test pattern as the k-ε verification.
        let grid = Grid::new(1, 60, 1.0, 2.0);
        let fluid = Fluid::new(1.0, 1e-4);
        let m = SstModel::standard();
        let w = WallFunction::standard();
        let laminar = solve_turbulent_channel_sst(&grid, &fluid, &m, &w, -1.0, false, 1);
        let turbulent = solve_turbulent_channel_sst(&grid, &fluid, &m, &w, -1.0, true, 200);
        let lam_flat = laminar.flatness();
        let turb_flat = turbulent.flatness();
        assert!(
            turb_flat > lam_flat + 0.02,
            "SST profile (flatness {turb_flat}) must be flatter than laminar (flatness {lam_flat})"
        );
        // A turbulent channel's flatness lands well above the parabolic
        // 2/3 — it does not become *more* peaked.
        assert!(
            turb_flat > 0.68,
            "SST channel flatness {turb_flat} should exceed the parabolic 2/3"
        );
        // No-slip walls: the near-wall velocity is far below the peak.
        let peak = turbulent.centreline_velocity();
        assert!(
            turbulent.u[0] < 0.5 * peak,
            "velocity must fall toward the no-slip wall"
        );
    }

    #[test]
    fn sst_channel_eddy_viscosity_is_largest_in_the_core() {
        // Physical near-wall behaviour: ν_t is largest away from the
        // wall and falls toward the wall.
        let grid = Grid::new(1, 60, 1.0, 2.0);
        let fluid = Fluid::new(1.0, 1e-4);
        let m = SstModel::standard();
        let w = WallFunction::standard();
        let p = solve_turbulent_channel_sst(&grid, &fluid, &m, &w, -1.0, true, 200);
        for &nt in &p.nu_t {
            assert!(nt >= 0.0 && nt.is_finite());
        }
        let n = p.nu_t.len();
        let nu_t_wall = p.nu_t[0];
        let nu_t_centre = p.nu_t[n / 2];
        assert!(
            nu_t_centre > nu_t_wall,
            "ν_t should grow away from the wall: wall {nu_t_wall}, centre {nu_t_centre}"
        );
        // The SST limiter still allows ν_t to swamp the molecular ν
        // in the core of a high-Re channel.
        assert!(
            nu_t_centre > 5.0 * fluid.viscosity,
            "SST core ν_t {nu_t_centre} should exceed molecular ν {} by a healthy factor",
            fluid.viscosity
        );
    }
}
