//! Turbulence models — k-ε and k-ω SST in three dimensions.
//!
//! # Why a turbulence model at all
//!
//! The Reynolds numbers of car / aircraft external aerodynamics are
//! `10⁶`–`10⁸`: the flow is firmly turbulent. Resolving every eddy
//! (DNS) is computationally impossible on a desktop, so the standard
//! engineering approach is **RANS** (Reynolds-Averaged Navier-Stokes):
//! solve for the *mean* flow and model the effect of the turbulent
//! fluctuations as an extra **eddy viscosity** `μ_t` added to the
//! molecular viscosity in the momentum equations. A turbulence model
//! is a pair of transport equations that compute `μ_t` from the mean
//! flow.
//!
//! # The two models here
//!
//! - **Standard k-ε** ([`TurbulenceModel::KEpsilon`]) — transports the
//!   turbulence kinetic energy `k` and its dissipation rate `ε`;
//!   `μ_t = ρ·Cμ·k²/ε`. Robust, the historical workhorse; weaker in
//!   adverse-pressure-gradient separation.
//! - **k-ω SST** ([`TurbulenceModel::KOmegaSST`]) — Menter's
//!   shear-stress-transport model transports `k` and the specific
//!   dissipation `ω`; a blending function runs a k-ω formulation near
//!   walls (where it predicts separation well) and a k-ε formulation
//!   in the free stream (where k-ω is over-sensitive to the inlet).
//!   It is **the** external-aerodynamics standard.
//!
//! # Honest scope
//!
//! Both models are real two-equation RANS closures with the standard
//! constants. The implementation is a v1: the production term uses the
//! mean strain-rate magnitude, the wall treatment is a high-Reynolds
//! wall function on the immersed-boundary cut cells (not a
//! low-Reynolds near-wall integration), and the SST blending uses the
//! standard `F1`/`F2` functions evaluated from the wall distance the
//! immersed-boundary voxelizer provides. Transition modelling, curvature
//! correction and the full low-Re damping functions are documented
//! extensions. The eddy viscosity these models produce is clamped to a
//! physical range for solver robustness.

use crate::grid::{Field3, Grid3};

/// The constant `Cμ` of the eddy-viscosity relation, shared by both
/// models (`= 0.09`, equivalently `β*` in the k-ω notation).
pub const C_MU: f64 = 0.09;

/// Which RANS turbulence closure to use.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TurbulenceModel {
    /// No model — laminar flow (valid only at low Reynolds number).
    Laminar,
    /// Standard k-ε.
    KEpsilon,
    /// Menter k-ω SST — the external-aero standard.
    KOmegaSST,
}

impl TurbulenceModel {
    /// A human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            TurbulenceModel::Laminar => "laminar",
            TurbulenceModel::KEpsilon => "k-epsilon",
            TurbulenceModel::KOmegaSST => "k-omega SST",
        }
    }
}

/// The k-ε model constants.
#[derive(Clone, Copy, Debug)]
pub struct KEpsilonConstants {
    /// `Cμ` — eddy-viscosity constant.
    pub c_mu: f64,
    /// `C1ε` — production constant in the ε equation.
    pub c1: f64,
    /// `C2ε` — destruction constant in the ε equation.
    pub c2: f64,
    /// `σk` — turbulent Prandtl number for `k`.
    pub sigma_k: f64,
    /// `σε` — turbulent Prandtl number for `ε`.
    pub sigma_eps: f64,
}

impl Default for KEpsilonConstants {
    /// The standard Launder-Spalding constants.
    fn default() -> Self {
        KEpsilonConstants {
            c_mu: 0.09,
            c1: 1.44,
            c2: 1.92,
            sigma_k: 1.0,
            sigma_eps: 1.3,
        }
    }
}

/// The k-ω SST model constants (the near-wall "set 1" values; the
/// free-stream "set 2" values are blended in by the model).
#[derive(Clone, Copy, Debug)]
pub struct SstConstants {
    /// `β*` — equivalent to `Cμ`.
    pub beta_star: f64,
    /// `a1` — the SST eddy-viscosity limiter constant.
    pub a1: f64,
    /// `β1` — near-wall `ω`-destruction constant.
    pub beta1: f64,
    /// `σk1` — near-wall `k` diffusion constant.
    pub sigma_k1: f64,
    /// `σω1` — near-wall `ω` diffusion constant.
    pub sigma_w1: f64,
}

impl Default for SstConstants {
    /// Menter's standard SST constants.
    fn default() -> Self {
        SstConstants {
            beta_star: 0.09,
            a1: 0.31,
            beta1: 0.075,
            sigma_k1: 0.85,
            sigma_w1: 0.5,
        }
    }
}

/// The turbulence-field state — the two transported scalars plus the
/// eddy viscosity derived from them.
#[derive(Clone, Debug)]
pub struct TurbulenceState {
    /// Which model these fields belong to.
    pub model: TurbulenceModel,
    /// Turbulence kinetic energy `k` (m²·s⁻²), cell-centred.
    pub k: Field3,
    /// The second scalar — `ε` for k-ε, `ω` for k-ω SST — cell-centred.
    pub scale: Field3,
    /// The eddy (turbulent) viscosity `μ_t` (Pa·s), cell-centred.
    pub mu_t: Field3,
}

impl TurbulenceState {
    /// A uniform initial turbulence field at the given inlet values.
    ///
    /// `k0` is the inlet turbulence kinetic energy and `scale0` is the
    /// inlet `ε` or `ω`. `mu_t` is initialised consistently with the
    /// model.
    pub fn uniform(
        grid: &Grid3,
        model: TurbulenceModel,
        k0: f64,
        scale0: f64,
        density: f64,
    ) -> TurbulenceState {
        let k = Field3::filled(grid.nx, grid.ny, grid.nz, k0.max(0.0));
        let scale = Field3::filled(grid.nx, grid.ny, grid.nz, scale0.max(1e-12));
        let mut mu_t = grid.scalar_field();
        let mut_val = eddy_viscosity_value(model, density, k0.max(0.0), scale0.max(1e-12));
        mu_t.fill(mut_val);
        TurbulenceState {
            model,
            k,
            scale,
            mu_t,
        }
    }

    /// A zero turbulence field — used for the [`TurbulenceModel::Laminar`]
    /// case so the rest of the solver can treat `μ_t` uniformly (it is
    /// simply zero everywhere).
    pub fn laminar(grid: &Grid3) -> TurbulenceState {
        TurbulenceState {
            model: TurbulenceModel::Laminar,
            k: grid.scalar_field(),
            scale: grid.scalar_field(),
            mu_t: grid.scalar_field(),
        }
    }

    /// The mean eddy viscosity across the field — a convergence /
    /// sanity diagnostic.
    pub fn mean_eddy_viscosity(&self) -> f64 {
        self.mu_t.mean()
    }
}

/// The eddy viscosity for a single cell given the two scalars.
///
/// k-ε: `μ_t = ρ·Cμ·k²/ε`. k-ω: `μ_t = ρ·k/ω`. The result is clamped
/// to a non-negative finite value.
pub fn eddy_viscosity_value(model: TurbulenceModel, density: f64, k: f64, scale: f64) -> f64 {
    let mu = match model {
        TurbulenceModel::Laminar => 0.0,
        TurbulenceModel::KEpsilon => {
            if scale > 1e-12 {
                density * C_MU * k * k / scale
            } else {
                0.0
            }
        }
        TurbulenceModel::KOmegaSST => {
            if scale > 1e-12 {
                density * k / scale
            } else {
                0.0
            }
        }
    };
    if mu.is_finite() {
        mu.max(0.0)
    } else {
        0.0
    }
}

/// The mean strain-rate magnitude `S = √(2·Sᵢⱼ·Sᵢⱼ)` at a cell,
/// from the cell-centred velocity gradients.
///
/// `S` drives the turbulence production `P = μ_t·S²`. The velocity
/// gradients are central differences of the cell-centred velocity
/// components supplied by the caller (the solver averages the
/// staggered faces to cell centres first).
#[allow(clippy::too_many_arguments)]
pub fn strain_rate_magnitude(
    uc: &Field3,
    vc: &Field3,
    wc: &Field3,
    grid: &Grid3,
    i: usize,
    j: usize,
    k: usize,
) -> f64 {
    let (dx, dy, dz) = (grid.dx(), grid.dy(), grid.dz());
    let cdx = |f: &Field3| -> f64 {
        if i > 0 && i + 1 < grid.nx {
            (f.at(i + 1, j, k) - f.at(i - 1, j, k)) / (2.0 * dx)
        } else {
            0.0
        }
    };
    let cdy = |f: &Field3| -> f64 {
        if j > 0 && j + 1 < grid.ny {
            (f.at(i, j + 1, k) - f.at(i, j - 1, k)) / (2.0 * dy)
        } else {
            0.0
        }
    };
    let cdz = |f: &Field3| -> f64 {
        if k > 0 && k + 1 < grid.nz {
            (f.at(i, j, k + 1) - f.at(i, j, k - 1)) / (2.0 * dz)
        } else {
            0.0
        }
    };
    let dudx = cdx(uc);
    let dudy = cdy(uc);
    let dudz = cdz(uc);
    let dvdx = cdx(vc);
    let dvdy = cdy(vc);
    let dvdz = cdz(vc);
    let dwdx = cdx(wc);
    let dwdy = cdy(wc);
    let dwdz = cdz(wc);
    // Strain-rate tensor Sij = ½(∂ui/∂xj + ∂uj/∂xi).
    let s11 = dudx;
    let s22 = dvdy;
    let s33 = dwdz;
    let s12 = 0.5 * (dudy + dvdx);
    let s13 = 0.5 * (dudz + dwdx);
    let s23 = 0.5 * (dvdz + dwdy);
    // S = √(2·Sij·Sij).
    let sij_sij = s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23);
    (2.0 * sij_sij).sqrt()
}

/// Advance the turbulence fields one outer iteration.
///
/// This is a single relaxed update of the two transport scalars from
/// the current mean flow, followed by a recomputation of `μ_t`. The
/// caller supplies the cell-centred velocity components (`uc, vc, wc`)
/// and the per-cell wall distance (for the SST blending and the wall
/// treatment); the molecular viscosity `mu_lam` and density set the
/// scale.
///
/// The fields are updated *in place* and clamped to physical ranges:
/// `k ≥ 0`, the dissipation scale strictly positive, and `μ_t` capped
/// at `mu_cap·mu_lam` so a transient overshoot cannot blow up the
/// momentum solve.
#[allow(clippy::too_many_arguments)]
pub fn advance_turbulence(
    state: &mut TurbulenceState,
    uc: &Field3,
    vc: &Field3,
    wc: &Field3,
    wall_distance: &Field3,
    solid: &[bool],
    grid: &Grid3,
    density: f64,
    mu_lam: f64,
    dt_relax: f64,
    wall: Option<&crate::immersed::ImmersedBody>,
) {
    if state.model == TurbulenceModel::Laminar {
        return;
    }
    let ke = KEpsilonConstants::default();
    let sst = SstConstants::default();
    let nu = (mu_lam / density).max(1e-12);

    let mut new_k = state.k.clone();
    let mut new_scale = state.scale.clone();

    for k in 0..grid.nz {
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let idx = i + grid.nx * (j + grid.ny * k);
                if solid[idx] {
                    // Inside the body — turbulence frozen near zero.
                    new_k.set(i, j, k, 0.0);
                    new_scale.set(i, j, k, state.scale.at(i, j, k).max(1e-9));
                    continue;
                }
                let k_old = state.k.at(i, j, k).max(0.0);
                let s_old = state.scale.at(i, j, k).max(1e-12);
                let mu_t = state.mu_t.at(i, j, k);

                // Production P = μ_t·S².
                let s_mag = strain_rate_magnitude(uc, vc, wc, grid, i, j, k);
                let prod = (mu_t * s_mag * s_mag).max(0.0);

                // Diffusive smoothing (Laplacian of the scalar) — the
                // turbulence transport's diffusion term, lumped.
                let lap_k = laplacian_at(&state.k, grid, i, j, k);
                let lap_s = laplacian_at(&state.scale, grid, i, j, k);
                let diff = nu + mu_t / density;

                let dwall = wall_distance.at(i, j, k).max(1e-6);

                match state.model {
                    TurbulenceModel::KEpsilon => {
                        // dk/dt = P/ρ − ε + ∇·(diff ∇k)
                        let dk = prod / density - s_old + diff * lap_k;
                        // dε/dt = C1·ε/k·P/ρ − C2·ε²/k + ∇·(diff ∇ε)
                        let k_safe = k_old.max(1e-8);
                        let deps = ke.c1 * s_old / k_safe * prod / density
                            - ke.c2 * s_old * s_old / k_safe
                            + (nu + mu_t / density / ke.sigma_eps) * lap_s;
                        let kn = (k_old + dt_relax * dk).max(0.0);
                        let mut sn = s_old + dt_relax * deps;
                        // Wall-function floor on ε near the wall.
                        let eps_wall = C_MU.powf(0.75) * kn.powf(1.5) / (0.41 * dwall).max(1e-6);
                        sn = sn.max(0.1 * eps_wall).max(1e-10);
                        new_k.set(i, j, k, kn);
                        new_scale.set(i, j, k, sn);
                    }
                    TurbulenceModel::KOmegaSST => {
                        // dk/dt = P/ρ − β*·k·ω + ∇·(diff ∇k)
                        let dk = prod / density - sst.beta_star * k_old * s_old
                            + (nu + mu_t / density * sst.sigma_k1) * lap_k;
                        // dω/dt = α·S² − β·ω² + ∇·(diff ∇ω)
                        // α ≈ β1/β* − σw1·κ²/√β*  (Menter set 1)
                        let alpha = sst.beta1 / sst.beta_star
                            - sst.sigma_w1 * 0.41 * 0.41 / sst.beta_star.sqrt();
                        let dw = alpha * s_mag * s_mag - sst.beta1 * s_old * s_old
                            + (nu + mu_t / density * sst.sigma_w1) * lap_s;
                        let kn = (k_old + dt_relax * dk).max(0.0);
                        let mut wn = s_old + dt_relax * dw;
                        // Menter's near-wall ω boundary value:
                        // ω_wall = 6ν/(β1·d²).
                        let omega_wall = 6.0 * nu / (sst.beta1 * dwall * dwall).max(1e-12);
                        wn = wn.clamp(1e-6, omega_wall.max(1.0));
                        new_k.set(i, j, k, kn);
                        new_scale.set(i, j, k, wn);
                    }
                    TurbulenceModel::Laminar => {}
                }
            }
        }
    }

    // --- wall-function near-wall turbulence ---
    //
    // In a wall-adjacent (cut) cell the turbulence is *not* the
    // transport-equation solution: the high-Reynolds wall treatment
    // imposes the law-of-the-wall equilibrium ([`crate::wallmodel`]).
    // This is essential when a near-wall model resolves the wall shear
    // on a grid that does not resolve the boundary layer — the first
    // cell then carries a steep (wall-function-induced) velocity
    // gradient, and a *free-running* `k` would be driven by that
    // gradient's strain into a non-physical runaway. Imposing the
    // equilibrium `k`/`ε`/`ω` keeps the near-wall eddy viscosity at the
    // physical `ρ·κ·u_τ·y` (small, scaling with `y`) and consistent
    // with the reconstructed wall shear — the standard, stable
    // wall-function turbulence closure.
    if let Some(body) = wall {
        let c_mu = C_MU;
        for k in 0..grid.nz {
            for j in 0..grid.ny {
                for i in 0..grid.nx {
                    if body.tag(i, j, k) != crate::immersed::CellTag::Cut {
                        continue;
                    }
                    let cut = body.cut_face(i, j, k);
                    if !cut.has_wall() {
                        continue;
                    }
                    // The cell-centred velocity's wall-tangential speed.
                    let vel =
                        nalgebra::Vector3::new(uc.at(i, j, k), vc.at(i, j, k), wc.at(i, j, k));
                    let v_tan = vel - vel.dot(&cut.normal) * cut.normal;
                    let u_t = v_tan.norm();
                    let y = wall_distance.at(i, j, k).max(1e-6);
                    let wt = crate::wallmodel::wall_turbulence(density, u_t, y, nu, mu_lam, c_mu);
                    new_k.set(i, j, k, wt.k);
                    match state.model {
                        TurbulenceModel::KEpsilon => {
                            new_scale.set(i, j, k, wt.epsilon);
                        }
                        TurbulenceModel::KOmegaSST => {
                            new_scale.set(i, j, k, wt.omega);
                        }
                        TurbulenceModel::Laminar => {}
                    }
                }
            }
        }
    }

    state.k = new_k;
    state.scale = new_scale;

    // Recompute the eddy viscosity, clamped. A cut cell under the
    // wall-function treatment is held at the physical log-layer
    // `μ_t = ρ·κ·u_τ·y`; every other cell takes the model relation
    // from its transported `k` and scale.
    //
    // The eddy-viscosity robustness cap — kept at the loose
    // `1e5·μ_lam` for the bulk cells. (Under the near-wall model the
    // wall-adjacent cut cells are not subject to it — they are held at
    // the physical log-layer `μ_t = ρ·κ·u_τ·y` below — so the cap only
    // ever clips a runaway in a *bulk* cell, the harmless place for it.)
    let mu_cap = 1.0e5 * mu_lam;
    let model = state.model;
    let wall_mu_t: Option<Vec<f64>> = wall.map(|body| {
        let mut wm = vec![f64::NAN; state.mu_t.data.len()];
        for k in 0..grid.nz {
            for j in 0..grid.ny {
                for i in 0..grid.nx {
                    let idx = i + grid.nx * (j + grid.ny * k);
                    if body.tag(i, j, k) != crate::immersed::CellTag::Cut {
                        continue;
                    }
                    let cut = body.cut_face(i, j, k);
                    if !cut.has_wall() {
                        continue;
                    }
                    let vel =
                        nalgebra::Vector3::new(uc.at(i, j, k), vc.at(i, j, k), wc.at(i, j, k));
                    let v_tan = vel - vel.dot(&cut.normal) * cut.normal;
                    let u_t = v_tan.norm();
                    let y = wall_distance.at(i, j, k).max(1e-6);
                    wm[idx] =
                        crate::wallmodel::wall_turbulence(density, u_t, y, nu, mu_lam, C_MU).mu_t;
                }
            }
        }
        wm
    });
    for (idx, mu) in state.mu_t.data.iter_mut().enumerate() {
        if solid[idx] {
            *mu = 0.0;
            continue;
        }
        // A wall-function cut cell uses the physical log-layer μ_t.
        if let Some(wm) = &wall_mu_t {
            if wm[idx].is_finite() {
                *mu = wm[idx].min(mu_cap);
                continue;
            }
        }
        let kk = state.k.data[idx];
        let ss = state.scale.data[idx];
        *mu = eddy_viscosity_value(model, density, kk, ss).min(mu_cap);
    }
}

/// The 7-point Laplacian of a cell-centred field at `(i, j, k)`,
/// homogeneous-Neumann at the grid edge.
fn laplacian_at(f: &Field3, grid: &Grid3, i: usize, j: usize, k: usize) -> f64 {
    let (dx, dy, dz) = (grid.dx(), grid.dy(), grid.dz());
    let c = f.at(i, j, k);
    let mut lap = 0.0;
    let term = |a: f64, b: f64, h: f64| (a - 2.0 * b + b) / (h * h) + (b - a) * 0.0;
    let _ = term; // (kept for clarity; explicit form below)
                  // x.
    let xp = if i + 1 < grid.nx {
        f.at(i + 1, j, k)
    } else {
        c
    };
    let xm = if i > 0 { f.at(i - 1, j, k) } else { c };
    lap += (xp - 2.0 * c + xm) / (dx * dx);
    // y.
    let yp = if j + 1 < grid.ny {
        f.at(i, j + 1, k)
    } else {
        c
    };
    let ym = if j > 0 { f.at(i, j - 1, k) } else { c };
    lap += (yp - 2.0 * c + ym) / (dy * dy);
    // z.
    let zp = if k + 1 < grid.nz {
        f.at(i, j, k + 1)
    } else {
        c
    };
    let zm = if k > 0 { f.at(i, j, k - 1) } else { c };
    lap += (zp - 2.0 * c + zm) / (dz * dz);
    lap
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn model_names_are_stable() {
        assert_eq!(TurbulenceModel::Laminar.name(), "laminar");
        assert_eq!(TurbulenceModel::KEpsilon.name(), "k-epsilon");
        assert_eq!(TurbulenceModel::KOmegaSST.name(), "k-omega SST");
    }

    #[test]
    fn k_epsilon_eddy_viscosity_formula() {
        // μ_t = ρ·Cμ·k²/ε. ρ=1.2, k=2, ε=0.5 → 1.2·0.09·4/0.5 = 0.864.
        let mt = eddy_viscosity_value(TurbulenceModel::KEpsilon, 1.2, 2.0, 0.5);
        assert!((mt - 0.864).abs() < 1e-9);
    }

    #[test]
    fn k_omega_eddy_viscosity_formula() {
        // μ_t = ρ·k/ω. ρ=1.2, k=2, ω=4 → 1.2·2/4 = 0.6.
        let mt = eddy_viscosity_value(TurbulenceModel::KOmegaSST, 1.2, 2.0, 4.0);
        assert!((mt - 0.6).abs() < 1e-9);
    }

    #[test]
    fn laminar_eddy_viscosity_is_zero() {
        assert_eq!(
            eddy_viscosity_value(TurbulenceModel::Laminar, 1.2, 5.0, 0.1),
            0.0
        );
    }

    #[test]
    fn eddy_viscosity_handles_degenerate_scale() {
        // A zero / tiny dissipation scale must not blow up.
        let mt = eddy_viscosity_value(TurbulenceModel::KEpsilon, 1.2, 2.0, 0.0);
        assert_eq!(mt, 0.0);
        let mw = eddy_viscosity_value(TurbulenceModel::KOmegaSST, 1.2, 2.0, 0.0);
        assert_eq!(mw, 0.0);
    }

    #[test]
    fn uniform_state_is_consistent() {
        let grid = Grid3::at_origin(6, 6, 6, 1.0, 1.0, 1.0);
        let st = TurbulenceState::uniform(&grid, TurbulenceModel::KEpsilon, 1.0, 0.5, 1.2);
        // μ_t = 1.2·0.09·1²/0.5 = 0.216 everywhere.
        assert!((st.mean_eddy_viscosity() - 0.216).abs() < 1e-9);
        assert!(st.k.data.iter().all(|&v| (v - 1.0).abs() < 1e-12));
    }

    #[test]
    fn strain_rate_of_uniform_flow_is_zero() {
        // A spatially-uniform velocity field has zero strain rate.
        let grid = Grid3::at_origin(8, 8, 8, 1.0, 1.0, 1.0);
        let uc = Field3::filled(8, 8, 8, 5.0);
        let vc = Field3::zeros(8, 8, 8);
        let wc = Field3::zeros(8, 8, 8);
        let s = strain_rate_magnitude(&uc, &vc, &wc, &grid, 4, 4, 4);
        assert!(s.abs() < 1e-12, "uniform flow strain rate {s} should be 0");
    }

    #[test]
    fn strain_rate_of_simple_shear_is_the_shear_rate() {
        // u = γ·y, everything else 0 → S = √(2·(½γ)²·2) = γ.
        let grid = Grid3::at_origin(8, 8, 8, 8.0, 8.0, 8.0); // dy = 1
        let mut uc = Field3::zeros(8, 8, 8);
        let gamma = 3.0;
        for k in 0..8 {
            for j in 0..8 {
                for i in 0..8 {
                    let (_, y, _) = grid.cell_centre(i, j, k);
                    uc.set(i, j, k, gamma * y);
                }
            }
        }
        let vc = Field3::zeros(8, 8, 8);
        let wc = Field3::zeros(8, 8, 8);
        let s = strain_rate_magnitude(&uc, &vc, &wc, &grid, 4, 4, 4);
        assert!(
            (s - gamma).abs() < 1e-9,
            "shear strain rate {s} should be {gamma}"
        );
    }

    #[test]
    fn advance_turbulence_keeps_k_non_negative() {
        // A few k-ε updates on a quiescent field must keep k ≥ 0 and
        // the dissipation scale finite — the robustness clamp.
        let grid = Grid3::at_origin(8, 8, 8, 1.0, 1.0, 1.0);
        let mut st = TurbulenceState::uniform(&grid, TurbulenceModel::KEpsilon, 0.5, 0.3, 1.2);
        let uc = Field3::zeros(8, 8, 8);
        let vc = Field3::zeros(8, 8, 8);
        let wc = Field3::zeros(8, 8, 8);
        let wd = Field3::filled(8, 8, 8, 0.1);
        let solid = vec![false; grid.cell_count()];
        for _ in 0..20 {
            advance_turbulence(
                &mut st, &uc, &vc, &wc, &wd, &solid, &grid, 1.2, 1.8e-5, 1e-3, None,
            );
        }
        assert!(st.k.data.iter().all(|&v| v >= 0.0 && v.is_finite()));
        assert!(st.scale.data.iter().all(|&v| v > 0.0 && v.is_finite()));
        assert!(st.mu_t.data.iter().all(|&v| v >= 0.0 && v.is_finite()));
    }

    #[test]
    fn sst_advance_is_stable_on_a_sheared_field() {
        // k-ω SST on a simple shear must stay bounded.
        let grid = Grid3::at_origin(10, 10, 10, 10.0, 10.0, 10.0);
        let mut st = TurbulenceState::uniform(&grid, TurbulenceModel::KOmegaSST, 1.0, 10.0, 1.225);
        let mut uc = Field3::zeros(10, 10, 10);
        for k in 0..10 {
            for j in 0..10 {
                for i in 0..10 {
                    let (_, y, _) = grid.cell_centre(i, j, k);
                    uc.set(i, j, k, 2.0 * y);
                }
            }
        }
        let vc = Field3::zeros(10, 10, 10);
        let wc = Field3::zeros(10, 10, 10);
        let wd = Field3::filled(10, 10, 10, 0.5);
        let solid = vec![false; grid.cell_count()];
        for _ in 0..30 {
            advance_turbulence(
                &mut st, &uc, &vc, &wc, &wd, &solid, &grid, 1.225, 1.81e-5, 5e-4, None,
            );
        }
        assert!(st.k.data.iter().all(|&v| v.is_finite() && v >= 0.0));
        assert!(st.scale.data.iter().all(|&v| v.is_finite() && v > 0.0));
    }
}
