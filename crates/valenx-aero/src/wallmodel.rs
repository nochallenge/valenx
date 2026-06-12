//! The near-wall model — a boundary-layer-resolving wall treatment for
//! the immersed-boundary cut cells.
//!
//! # Why a wall model at all
//!
//! The accuracy ceiling of an immersed-boundary CFD on a *uniform*
//! Cartesian background grid is the near-wall region. At a car / wing
//! Reynolds number `10⁶`–`10⁸` the turbulent boundary layer is
//! microscopically thin — for a `1 m` body at `Re = 10⁶` the layer is
//! `δ ≈ 0.072·L·Re⁻¹ᐟ⁵ ≈ 4.5 mm`. A practical Cartesian wind-tunnel
//! grid resolves the body with tens of cells, so the *first fluid cell*
//! sits one whole cell — often several boundary-layer thicknesses —
//! from the wall. Treating the wall shear as a plain linear gradient
//! `τ_w = μ·u₁/y₁` over that distance is then badly wrong: the velocity
//! does **not** vary linearly across a turbulent boundary layer, it
//! follows the **law of the wall**. A linear gradient under-resolves
//! the near-wall momentum loss, the boundary layer comes out too thick,
//! it separates too early, and the integrated **pressure drag** is
//! over-predicted. That is the documented reason a sphere's `Cd` stayed
//! above the textbook `≈ 0.5` even with the cut-cell wall geometry.
//!
//! # What this module does
//!
//! It reconstructs the turbulent boundary-layer velocity profile and
//! recovers the wall shear stress `τ_w` *self-consistently* from the
//! first-cell velocity, using **Spalding's law of the wall** — a single
//! smooth implicit relation valid continuously across the viscous
//! sublayer, the buffer layer and the logarithmic layer:
//!
//! ```text
//!   y⁺ = u⁺ + e^(−κB)·[ e^(κu⁺) − 1 − κu⁺ − (κu⁺)²/2 − (κu⁺)³/6 ]
//! ```
//!
//! with `y⁺ = y·u_τ/ν`, `u⁺ = u_t/u_τ`, `κ = 0.41`, `B = 5.2`. Given
//! the tangential speed `u_t` at the first-cell wall distance `y`, this
//! is one nonlinear equation in the **friction velocity** `u_τ`; a
//! Newton iteration solves it, and the wall shear is `τ_w = ρ·u_τ²`.
//!
//! Unlike a plain log-law wall function — which is only valid for a
//! first cell that happens to land in the log layer (`30 ≲ y⁺ ≲ 300`)
//! and gives nonsense for a cell in the viscous sublayer or the buffer
//! layer — the Spalding blend is correct **wherever the first cell
//! lands**. On a uniform Cartesian grid the first-cell `y⁺` is not
//! controllable (it varies cell to cell and case to case), so a blended
//! wall model, not a bare log law, is what the geometry demands.
//!
//! The recovered `τ_w` is used two ways:
//!
//! 1. as a **wall-function effective viscosity** `μ_w = τ_w·y/u_t` fed
//!    into the momentum equation's no-slip wall-drag term — so the
//!    near-wall momentum sink, hence the separation point and the
//!    pressure drag, are physically correct;
//! 2. as the **wall shear stress** in the surface-force integration —
//!    so the friction drag and the reported skin-friction coefficient
//!    `Cf` reflect the real turbulent profile, not a linear guess.
//!
//! # Honest scope
//!
//! A real, standard high-Reynolds **wall-function** treatment — the
//! Spalding all-`y⁺` law of the wall with the Reichardt-consistent
//! constants. It reconstructs the *equilibrium* turbulent boundary-layer
//! profile, which is the right model for an attached or mildly
//! adverse-pressure-gradient boundary layer. It is **not** a
//! non-equilibrium / pressure-gradient-sensitised wall model, and it is
//! still not a body-fitted near-wall **prism layer** that resolves the
//! sublayer directly — that remains the documented Tier-3 residue. What
//! the wall model *does* close is the crude-linear-gradient error: the
//! near-wall shear and the momentum loss are now reconstructed from the
//! turbulent law of the wall, which measurably moves the surface forces
//! toward the published references (see the benchmark suite).

/// The von Kármán constant `κ` of the law of the wall.
pub const KAPPA: f64 = 0.41;

/// The additive constant `B` of the logarithmic law of the wall for a
/// smooth wall (`u⁺ = (1/κ)·ln y⁺ + B`). The paired smooth-wall value
/// for `κ = 0.41` is `B ≈ 5.2`.
pub const B_LOGLAW: f64 = 5.2;

/// The lower edge of the logarithmic layer in wall units — below this
/// `y⁺` the buffer / viscous sublayer dominates.
pub const Y_PLUS_LOG_LOWER: f64 = 11.0;

/// Spalding's law of the wall: the wall-normal coordinate `y⁺` as the
/// implicit function of the velocity `u⁺`.
///
/// ```text
///   y⁺(u⁺) = u⁺ + e^(−κB)·[ e^(κu⁺) − 1 − κu⁺ − (κu⁺)²/2 − (κu⁺)³/6 ]
/// ```
///
/// This single relation reproduces `u⁺ = y⁺` in the viscous sublayer
/// (the bracket vanishes to fourth order in `κu⁺`) and `u⁺ = (1/κ)·ln
/// y⁺ + B` in the log layer (the `e^(κu⁺)` term dominates), with the
/// buffer layer smoothly blended between — so it is valid for any first
/// cell, sublayer to log layer.
#[inline]
pub fn spalding_y_plus(u_plus: f64) -> f64 {
    let ku = KAPPA * u_plus;
    let bracket = ku.exp() - 1.0 - ku - 0.5 * ku * ku - ku * ku * ku / 6.0;
    u_plus + (-KAPPA * B_LOGLAW).exp() * bracket
}

/// The derivative `d y⁺ / d u⁺` of [`spalding_y_plus`] — used by the
/// Newton solve for the friction velocity.
#[inline]
fn spalding_dy_plus_du_plus(u_plus: f64) -> f64 {
    let ku = KAPPA * u_plus;
    // d/du⁺ of the bracket: κ·[ e^(κu⁺) − 1 − κu⁺ − (κu⁺)²/2 ].
    let dbracket = KAPPA * (ku.exp() - 1.0 - ku - 0.5 * ku * ku);
    1.0 + (-KAPPA * B_LOGLAW).exp() * dbracket
}

/// The friction velocity `u_τ` recovered from a near-wall sample by
/// solving Spalding's law of the wall.
///
/// `u_t` is the wall-tangential speed at wall-normal distance `y`;
/// `nu` is the (kinematic) viscosity. Returns the `u_τ ≥ 0` for which
/// the sample `(u_t, y)` lies on the Spalding profile. The relation
/// `y⁺ = y·u_τ/ν` and `u⁺ = u_t/u_τ` turns the law into one nonlinear
/// equation in `u_τ`; a damped Newton iteration (started from the
/// log-law estimate) solves it robustly.
///
/// Degenerate inputs — a non-positive `u_t`, `y` or `nu` — return `0`.
pub fn friction_velocity(u_t: f64, y: f64, nu: f64) -> f64 {
    let u_t = u_t.abs();
    if !(u_t > 0.0 && y > 0.0 && nu > 0.0 && u_t.is_finite() && y.is_finite()) {
        return 0.0;
    }
    // Residual in u_τ: F(u_τ) = spalding_y_plus(u_t/u_τ) − y·u_τ/ν = 0.
    //
    // Newton on u_τ. dF/du_τ = y⁺'(u⁺)·(−u_t/u_τ²) − y/ν.
    // Both terms of dF are negative, so F is strictly decreasing — the
    // root is unique and Newton converges from any positive start.
    //
    // Initial guess: the pure-log-law u_τ if the cell looks like it is
    // in the log layer, else a viscous-sublayer estimate. A robust
    // bracket then guards the iteration.
    let nu_y = nu / y;
    // Viscous-sublayer guess: u⁺ = y⁺ ⇒ u_τ = √(u_t·ν/y).
    let visc_guess = (u_t * nu_y).sqrt();
    // Log-layer guess: solve u_t/u_τ = (1/κ)ln(y u_τ/ν)+B by a couple
    // of fixed-point steps.
    let mut log_guess = visc_guess.max(1e-12);
    for _ in 0..6 {
        let y_plus = (y * log_guess / nu).max(1.0001);
        let u_plus = (1.0 / KAPPA) * y_plus.ln() + B_LOGLAW;
        if u_plus > 0.0 {
            log_guess = u_t / u_plus;
        }
    }
    let mut u_tau = if u_t * y / nu > 30.0 {
        log_guess.max(visc_guess * 0.1)
    } else {
        visc_guess
    }
    .max(1e-12);

    for _ in 0..25 {
        let u_plus = u_t / u_tau;
        let f = spalding_y_plus(u_plus) - y * u_tau / nu;
        // dF/du_τ.
        let df = spalding_dy_plus_du_plus(u_plus) * (-u_t / (u_tau * u_tau)) - y / nu;
        if !df.is_finite() || df.abs() < 1e-300 {
            break;
        }
        let mut step = f / df;
        // Damp so a Newton step never overshoots past zero or explodes
        // — the Spalding exponential is stiff for a small u_τ.
        let max_step = 0.5 * u_tau;
        if step > max_step {
            step = max_step;
        } else if step < -max_step {
            step = -max_step;
        }
        let next = (u_tau - step).max(1e-13);
        if (next - u_tau).abs() <= 1e-12 * next.max(1e-9) {
            u_tau = next;
            break;
        }
        u_tau = next;
    }
    if u_tau.is_finite() {
        u_tau.max(0.0)
    } else {
        0.0
    }
}

/// The wall shear stress `τ_w` from a near-wall velocity sample, using
/// the reconstructed turbulent boundary-layer profile.
///
/// `rho` is the density, `u_t` the wall-tangential speed at wall-normal
/// distance `y`, `nu` the kinematic viscosity. `τ_w = ρ·u_τ²` with the
/// friction velocity `u_τ` from [`friction_velocity`]. This is the
/// physically-correct wall shear — the law-of-the-wall reconstruction,
/// not the crude linear gradient `μ·u_t/y` which under-resolves the
/// turbulent profile.
pub fn wall_shear_stress(rho: f64, u_t: f64, y: f64, nu: f64) -> f64 {
    let u_tau = friction_velocity(u_t, y, nu);
    rho * u_tau * u_tau
}

/// The dimensionless wall distance `y⁺` of a near-wall sample.
///
/// `y⁺ = y·u_τ/ν` with `u_τ` from [`friction_velocity`]. A `y⁺` in the
/// rough band `30–300` means the first cell lands in the log layer, the
/// classic wall-function-friendly placement; the Spalding model used
/// here stays valid outside that band too.
pub fn y_plus(u_t: f64, y: f64, nu: f64) -> f64 {
    let u_tau = friction_velocity(u_t, y, nu);
    y * u_tau / nu.max(1e-30)
}

/// The **wall-model effective viscosity** at a near-wall sample — the
/// turbulent momentum-transport coefficient the wall shear implies.
///
/// The wall shear stress recovered from the law of the wall is
/// `τ_w = μ_w · u_t / y`, so the effective viscosity that, applied as a
/// plain gradient over the first cell, reproduces the *correct*
/// turbulent wall shear is
///
/// ```text
///   μ_w = τ_w · y / u_t = ρ · u_τ² · y / u_t .
/// ```
///
/// This is what the momentum solver's no-slip wall-drag term consumes:
/// substituting `μ_w` for the plain molecular+eddy viscosity makes the
/// near-wall momentum sink match the turbulent boundary layer instead
/// of a laminar linear gradient. It is floored at the molecular
/// viscosity `mu_lam` (the wall can never transport *less* momentum
/// than molecular diffusion) and is finite for a vanishing `u_t`.
pub fn wall_effective_viscosity(rho: f64, u_t: f64, y: f64, nu: f64, mu_lam: f64) -> f64 {
    let u_t = u_t.abs();
    if !(u_t > 1e-12 && y > 0.0 && u_t.is_finite()) {
        return mu_lam.max(0.0);
    }
    let tau_w = wall_shear_stress(rho, u_t, y, nu);
    let mu_w = tau_w * y / u_t;
    if mu_w.is_finite() {
        mu_w.max(mu_lam.max(0.0))
    } else {
        mu_lam.max(0.0)
    }
}

/// The local turbulence kinetic energy `k` implied by the wall shear,
/// for the near-wall production / boundary value of the turbulence
/// model — the equilibrium relation `k = u_τ² / √(Cμ)`.
///
/// In an equilibrium turbulent boundary layer the turbulence energy and
/// the wall shear are tied: `−u'v' ≈ u_τ²` and `k ≈ u_τ²/√(Cμ)`. A
/// high-Reynolds wall treatment imposes this `k` (rather than solving
/// the `k` transport equation down to the wall), which keeps the eddy
/// viscosity consistent with the reconstructed wall shear.
pub fn wall_tke(u_tau: f64, c_mu: f64) -> f64 {
    if c_mu > 0.0 {
        u_tau * u_tau / c_mu.sqrt()
    } else {
        0.0
    }
}

/// The wall-function equilibrium turbulence state of a wall-adjacent
/// cell — the values a high-Reynolds wall treatment *imposes* there
/// instead of integrating the turbulence transport equations down to
/// the wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallTurbulence {
    /// The friction velocity `u_τ` (m·s⁻¹).
    pub u_tau: f64,
    /// The turbulence kinetic energy `k = u_τ²/√Cμ` (m²·s⁻²).
    pub k: f64,
    /// The eddy viscosity `μ_t = ρ·κ·u_τ·y` (Pa·s) — the log-layer
    /// mixing-length value.
    pub mu_t: f64,
    /// The dissipation rate `ε = u_τ³/(κ·y)` (m²·s⁻³).
    pub epsilon: f64,
    /// The specific dissipation `ω = u_τ/(√Cμ·κ·y)` (s⁻¹).
    pub omega: f64,
}

/// The wall-function equilibrium turbulence of a wall-adjacent cell.
///
/// In the logarithmic layer the turbulence is in local equilibrium —
/// production balances dissipation — and every quantity follows from
/// the friction velocity `u_τ` and the wall distance `y`:
///
/// ```text
///   k   = u_τ² / √Cμ                 (the equilibrium TKE)
///   μ_t = ρ · κ · u_τ · y            (the log-layer mixing length)
///   ε   = u_τ³ / (κ · y)             (production = dissipation)
///   ω   = u_τ / (√Cμ · κ · y) = ε/(Cμ·k)
/// ```
///
/// A high-Reynolds wall treatment **imposes** these in the wall-adjacent
/// cell rather than integrating the `k`/`ε`/`ω` transport equations to
/// the wall — the standard wall-function turbulence closure. Imposing
/// them keeps the near-wall eddy viscosity *physical* (`μ_t ∝ y`, small
/// at the wall) and consistent with the reconstructed wall shear, which
/// is what a near-wall model on a boundary-layer-under-resolving grid
/// needs: without it the free-running `k`/`μ_t` near a wall is driven by
/// the steep (wall-function-induced) grid velocity gradient into a
/// non-physical runaway.
///
/// `rho`, `u_t`, `y`, `nu`, `mu_lam` are as in [`friction_velocity`] /
/// [`wall_effective_viscosity`]; `c_mu` is the model constant `Cμ`.
pub fn wall_turbulence(
    rho: f64,
    u_t: f64,
    y: f64,
    nu: f64,
    mu_lam: f64,
    c_mu: f64,
) -> WallTurbulence {
    let u_tau = friction_velocity(u_t, y, nu);
    let sqrt_cmu = c_mu.max(1e-12).sqrt();
    let k = u_tau * u_tau / sqrt_cmu;
    let y_safe = y.max(1e-9);
    // The log-layer mixing-length eddy viscosity, floored at molecular.
    let mu_t = (rho * KAPPA * u_tau * y_safe).max(mu_lam.max(0.0));
    let epsilon = u_tau * u_tau * u_tau / (KAPPA * y_safe);
    // ω = ε/(Cμ·k); equivalently u_τ/(√Cμ·κ·y).
    let omega = if k > 1e-12 {
        epsilon / (c_mu.max(1e-12) * k)
    } else {
        u_tau / (sqrt_cmu * KAPPA * y_safe)
    };
    WallTurbulence {
        u_tau,
        k: if k.is_finite() { k.max(0.0) } else { 0.0 },
        mu_t: if mu_t.is_finite() {
            mu_t
        } else {
            mu_lam.max(0.0)
        },
        epsilon: if epsilon.is_finite() {
            epsilon.max(1e-12)
        } else {
            1e-12
        },
        omega: if omega.is_finite() {
            omega.max(1e-6)
        } else {
            1e-6
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spalding_reduces_to_the_linear_law_in_the_viscous_sublayer() {
        // For a small u⁺ the Spalding bracket vanishes to 4th order, so
        // y⁺ ≈ u⁺ — the viscous sublayer u⁺ = y⁺.
        for u_plus in [0.1, 0.5, 1.0, 2.0] {
            let y_plus = spalding_y_plus(u_plus);
            assert!(
                (y_plus - u_plus).abs() < 0.15 * u_plus.max(0.2),
                "sublayer: y+({u_plus}) = {y_plus} should ≈ u+"
            );
        }
    }

    #[test]
    fn spalding_reduces_to_the_log_law_in_the_log_layer() {
        // For a large u⁺ the Spalding relation must reproduce the log
        // law u⁺ = (1/κ)·ln y⁺ + B.
        for u_plus in [16.0, 20.0, 24.0] {
            let y_plus = spalding_y_plus(u_plus);
            let log_u_plus = (1.0 / KAPPA) * y_plus.ln() + B_LOGLAW;
            assert!(
                (log_u_plus - u_plus).abs() < 0.25,
                "log layer: u+ = {u_plus}, log-law gives {log_u_plus}"
            );
        }
    }

    #[test]
    fn friction_velocity_inverts_the_spalding_profile() {
        // Pick a u_τ, place a sample exactly on the Spalding profile,
        // and confirm friction_velocity recovers the u_τ.
        let nu = 1.5e-5;
        let y = 0.01;
        for u_tau_true in [0.05, 0.2, 0.8, 2.0] {
            // y⁺ and u⁺ on the profile for this u_τ.
            let y_plus = y * u_tau_true / nu;
            // Invert Spalding numerically for u⁺(y⁺): bisection.
            let mut lo = 0.0;
            let mut hi = 60.0;
            for _ in 0..100 {
                let mid = 0.5 * (lo + hi);
                if spalding_y_plus(mid) < y_plus {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let u_plus = 0.5 * (lo + hi);
            let u_t = u_plus * u_tau_true;
            let recovered = friction_velocity(u_t, y, nu);
            let rel = (recovered - u_tau_true).abs() / u_tau_true;
            assert!(
                rel < 1e-3,
                "u_τ recovery: true {u_tau_true}, got {recovered} (rel {rel})"
            );
        }
    }

    #[test]
    fn wall_shear_stress_exceeds_the_linear_estimate_for_a_turbulent_cell() {
        // The whole point: for a first cell well outside the viscous
        // sublayer (a large y⁺), the law-of-the-wall shear is much
        // larger than the crude linear gradient μ·u_t/y, because the
        // turbulent profile is far steeper at the wall than a straight
        // line from the first cell would suggest.
        let rho = 1.225;
        let nu = 1.5e-5;
        let mu = rho * nu;
        let u_t = 20.0;
        let y = 0.02; // a coarse-grid first-cell distance
        let tau_wall_model = wall_shear_stress(rho, u_t, y, nu);
        let tau_linear = mu * u_t / y;
        assert!(
            tau_wall_model > 5.0 * tau_linear,
            "law-of-the-wall τ_w {tau_wall_model} should greatly exceed \
             the linear estimate {tau_linear}"
        );
        // And it must be a finite positive stress.
        assert!(tau_wall_model.is_finite() && tau_wall_model > 0.0);
    }

    #[test]
    fn wall_effective_viscosity_is_floored_at_molecular() {
        // A vanishing tangential velocity (a stagnation cell) must fall
        // back to the molecular viscosity, never below it.
        let rho = 1.225;
        let nu = 1.5e-5;
        let mu = rho * nu;
        let mu_w0 = wall_effective_viscosity(rho, 0.0, 0.01, nu, mu);
        assert!((mu_w0 - mu).abs() < 1e-12);
        // A moving cell gives an effective viscosity at or above
        // molecular (turbulence only ever adds transport).
        let mu_w = wall_effective_viscosity(rho, 15.0, 0.02, nu, mu);
        assert!(mu_w >= mu, "wall μ_w {mu_w} must be ≥ molecular {mu}");
    }

    #[test]
    fn y_plus_grows_with_wall_distance_and_speed() {
        let nu = 1.5e-5;
        // Faster flow → larger y⁺.
        let slow = y_plus(5.0, 0.01, nu);
        let fast = y_plus(40.0, 0.01, nu);
        assert!(fast > slow, "y+ should grow with speed");
        // Cell further from the wall → larger y⁺.
        let near = y_plus(20.0, 0.005, nu);
        let far = y_plus(20.0, 0.04, nu);
        assert!(far > near, "y+ should grow with wall distance");
        assert!(slow.is_finite() && fast.is_finite());
    }

    #[test]
    fn degenerate_inputs_are_handled() {
        let nu = 1.5e-5;
        // Zero / negative velocity, distance or viscosity → zero u_τ.
        assert_eq!(friction_velocity(0.0, 0.01, nu), 0.0);
        assert_eq!(friction_velocity(10.0, 0.0, nu), 0.0);
        assert_eq!(friction_velocity(10.0, 0.01, 0.0), 0.0);
        assert_eq!(friction_velocity(f64::NAN, 0.01, nu), 0.0);
        // A negative u_t is treated by magnitude (the shear sign is the
        // caller's; the model needs the speed).
        assert!(friction_velocity(-10.0, 0.01, nu) > 0.0);
    }

    #[test]
    fn wall_tke_follows_the_equilibrium_relation() {
        // k = u_τ²/√(Cμ). With Cμ = 0.09, √Cμ = 0.3, so k = u_τ²/0.3.
        let k = wall_tke(0.6, 0.09);
        assert!((k - 0.36 / 0.3).abs() < 1e-9);
        assert_eq!(wall_tke(0.5, 0.0), 0.0);
    }

    #[test]
    fn wall_turbulence_is_a_consistent_equilibrium_state() {
        // The wall-function equilibrium turbulence: k, μ_t, ε, ω all
        // derived from one friction velocity, mutually consistent.
        let rho = 1.225;
        let nu = 1.5e-5;
        let mu = rho * nu;
        let wt = wall_turbulence(rho, 20.0, 0.02, nu, mu, 0.09);
        // Every quantity is finite and physical.
        assert!(wt.u_tau > 0.0 && wt.u_tau.is_finite());
        assert!(wt.k > 0.0 && wt.k.is_finite());
        assert!(wt.mu_t >= mu && wt.mu_t.is_finite());
        assert!(wt.epsilon > 0.0 && wt.epsilon.is_finite());
        assert!(wt.omega > 0.0 && wt.omega.is_finite());
        // k = u_τ²/√Cμ.
        assert!((wt.k - wt.u_tau * wt.u_tau / 0.09_f64.sqrt()).abs() < 1e-9);
        // The eddy-viscosity identity μ_t = ρ·k/ω must hold (k-ω) — the
        // wall ω is defined so it does.
        let mu_t_komega = rho * wt.k / wt.omega;
        assert!(
            (mu_t_komega - wt.mu_t).abs() < 1e-6 * wt.mu_t.max(1e-9),
            "μ_t {} vs ρk/ω {}",
            wt.mu_t,
            mu_t_komega
        );
        // μ_t = ρ·κ·u_τ·y — the log-layer mixing length.
        assert!(
            (wt.mu_t - rho * KAPPA * wt.u_tau * 0.02).abs() < 1e-9 * wt.mu_t,
            "μ_t should be the ρκu_τy mixing-length value"
        );
    }

    #[test]
    fn wall_turbulence_eddy_viscosity_is_physical_not_runaway() {
        // The key property: the wall-function eddy viscosity is the
        // physical ρκu_τy — small near the wall, scaling with y — not a
        // runaway. For a typical near-wall cell it is O(0.01–0.1) Pa·s,
        // far below a runaway 1e5·μ_lam.
        let rho = 1.225;
        let nu = 1.5e-5;
        let mu = rho * nu;
        let wt = wall_turbulence(rho, 20.0, 0.02, nu, mu, 0.09);
        assert!(
            wt.mu_t < 1.0,
            "wall-function μ_t {} should be a modest physical value",
            wt.mu_t
        );
        // It grows with wall distance (mixing length ∝ y).
        let near = wall_turbulence(rho, 20.0, 0.005, nu, mu, 0.09);
        let far = wall_turbulence(rho, 20.0, 0.05, nu, mu, 0.09);
        assert!(far.mu_t > near.mu_t, "μ_t should grow with y");
    }

    #[test]
    fn friction_velocity_scales_linearly_with_speed_in_the_log_layer() {
        // Deep in the log layer u_τ is very nearly proportional to u_t
        // (the log term varies slowly), so doubling the speed roughly
        // doubles u_τ — a physical sanity check.
        let nu = 1.5e-5;
        let y = 0.03;
        let utau1 = friction_velocity(20.0, y, nu);
        let utau2 = friction_velocity(40.0, y, nu);
        let ratio = utau2 / utau1;
        assert!(
            (1.7..=2.3).contains(&ratio),
            "u_τ should scale ≈ linearly with speed, ratio {ratio}"
        );
    }
}
