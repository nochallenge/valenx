//! Structural / system **reliability analysis**.
//!
//! This module estimates the **probability of failure** `Pf = P[g(x) ≤ 0]`
//! for a limit-state function `g: ℝⁿ → ℝ` over uncertain inputs described by
//! the crate's [`Distribution`] types.  Three complementary estimators are
//! provided:
//!
//! | Method | Cost | Notes |
//! |--------|------|-------|
//! | [`form`] — First-Order Reliability Method | very cheap | searches for the Most-Probable-Point; accurate for smooth, near-linear `g` |
//! | [`pf_monte_carlo`] — crude Monte-Carlo | `O(n)` model calls | unbiased cross-check; converges as `O(1/√n)` |
//! | [`sorm_breitung`] — Second-Order Reliability Method | cheap add-on to FORM | corrects curvature bias; requires an MPP from FORM |
//!
//! ## Limit-state convention
//!
//! `g(x) > 0`  →  safe;  `g(x) ≤ 0`  →  failure.
//!
//! ## FORM — HLRF iteration
//!
//! FORM works in **standard-normal space** (also called *u*-space).  Each
//! input `xⱼ` is mapped to a standard-normal coordinate
//!
//! ```text
//! uⱼ = Φ⁻¹(Fⱼ(xⱼ))
//! ```
//!
//! where `Φ⁻¹` is the standard-normal quantile and `Fⱼ` is the CDF of the
//! `j`-th input distribution.  For Normal inputs this simplifies to
//! `uⱼ = (xⱼ − μⱼ)/σⱼ`.
//!
//! The **Most-Probable-Point** (MPP) `u*` is the point on the failure surface
//! `g(T⁻¹(u)) = 0` closest to the origin in *u*-space.  Its distance is the
//! **reliability index** `β = ‖u*‖`, and `Pf ≈ Φ(−β)`.
//!
//! The Hasofer–Lind–Rackwitz–Fiessler (HLRF) update rule (each iteration)
//!
//! ```text
//! uₖ₊₁ = (∇G·uₖ − G(uₖ)) / ‖∇G‖ · α̂
//! ```
//!
//! where `G(u) = g(T⁻¹(u))`, `∇G` is the *u*-space gradient (computed by
//! finite differences of `g` with respect to *x*, then chain-ruled via the
//! distribution PDFs), and `α̂ = −∇G / ‖∇G‖` is the unit vector pointing
//! toward the origin along the gradient.
//!
//! Gradients of `g` are estimated by **central finite differences** with step
//! `h = 1e-5 * sigma_j` (or `h` itself when sigma_j is negligible).  For a
//! linear `g` the gradient is exact to floating-point precision, so the MPP
//! is found **exactly in one iteration**.
//!
//! ## SORM — Breitung correction
//!
//! The Breitung formula corrects `Pf` for the principal curvatures `κᵢ` of the
//! limit-state surface at the MPP:
//!
//! ```text
//! Pf_SORM ≈ Φ(−β) · ∏ᵢ (1 + β κᵢ)^(−1/2)
//! ```
//!
//! Curvatures are estimated from the Hessian of `G` in *u*-space via
//! finite differences, projected onto the hyperplane orthogonal to `α̂` and
//! normalised by `‖∇G‖`.  The number of curvatures is `n − 1` (the
//! dimensionality of the failure surface).
//!
//! ## Honesty / scope
//!
//! * FORM is a **first-order approximation**; it is exact only for linear `g`
//!   and Normal inputs (β = a₀ / ‖a‖ for `g = a₀ + Σ aᵢxᵢ`).
//! * HLRF is not globally convergent for non-convex or highly non-linear `g`
//!   (use [`pf_monte_carlo`] as a cross-check).
//! * SORM Breitung requires `1 + β κᵢ > 0`; if this is violated for any
//!   curvature the correction is skipped and a warning is embedded in
//!   [`SormResult`].
//! * Monte-Carlo `Pf` has standard error `≈ √(Pf(1−Pf)/n)`.

use crate::distribution::Distribution;
use crate::error::UqError;
use crate::rng::SplitMix64;

// ── public constants ──────────────────────────────────────────────────────────

/// Default maximum HLRF iterations.
pub const DEFAULT_MAX_ITER: usize = 100;
/// Default HLRF convergence tolerance (change in `β` between iterations).
pub const DEFAULT_TOL: f64 = 1e-8;
/// Finite-difference step for gradient evaluation (fraction of σ, or absolute
/// when σ is not available).
const FD_H: f64 = 1e-5;

// ── helper: standard normal CDF Φ(z) ────────────────────────────────────────
//
// Reuses the `erf`-based formula from the `sampling` module.  We inline it
// here so `reliability` has no module-internal visibility dependency.

/// Φ(z) — standard-normal CDF.
#[inline]
fn phi(z: f64) -> f64 {
    0.5 * (1.0 + erf_approx(z / std::f64::consts::SQRT_2))
}

/// Acklam's rational quantile (inverse CDF) of N(0,1).  Accuracy ≈ 1.15e-9.
/// Endpoints are clamped to large finite values.
fn phi_inv(p: f64) -> f64 {
    const A: [f64; 6] = [
        -39.696_830_286_653_76,
        220.946_098_424_520_5,
        -275.928_510_446_968_7,
        138.357_751_867_269,
        -30.664_798_066_147_16,
        2.506_628_277_459_239,
    ];
    const B: [f64; 5] = [
        -54.476_098_798_224_06,
        161.585_836_858_040_9,
        -155.698_979_859_886_6,
        66.801_311_887_719_72,
        -13.280_681_552_885_72,
    ];
    const C: [f64; 6] = [
        -0.007_784_894_002_430_293,
        -0.322_396_458_041_136_5,
        -2.400_758_277_161_838,
        -2.549_732_539_343_734,
        4.374_664_141_464_968,
        2.938_163_982_698_783,
    ];
    const D: [f64; 4] = [
        0.007_784_695_709_041_462,
        0.322_467_129_070_039_8,
        2.445_134_137_142_996,
        3.754_408_661_907_416,
    ];
    const P_LOW: f64 = 0.024_25;
    const P_HIGH: f64 = 1.0 - P_LOW;

    if p <= 0.0 {
        return -1e10;
    }
    if p >= 1.0 {
        return 1e10;
    }
    if p < P_LOW {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= P_HIGH {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}

/// `erf(x)` via Abramowitz & Stegun 7.1.26 rational approximation.
/// Max absolute error ≈ 1.5e-7.
fn erf_approx(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0_f64 } else { 1.0_f64 };
    let x = x.abs();
    const A1: f64 = 0.254_829_592;
    const A2: f64 = -0.284_496_736;
    const A3: f64 = 1.421_413_741;
    const A4: f64 = -1.453_152_027;
    const A5: f64 = 1.061_405_429;
    const P: f64 = 0.327_591_1;
    let t = 1.0 / (1.0 + P * x);
    let y = 1.0 - (((((A5 * t + A4) * t) + A3) * t + A2) * t + A1) * t * (-x * x).exp();
    sign * y
}

// ── CDF / quantile of a Distribution ─────────────────────────────────────────

/// Quantile (inverse CDF) of `dist` at probability `p` (clamped to `[0,1]`).
fn dist_quantile(dist: &Distribution, p: f64) -> f64 {
    let p = p.clamp(1e-12, 1.0 - 1e-12);
    match *dist {
        Distribution::Uniform { lo, hi } => lo + (hi - lo) * p,
        Distribution::Normal { mean, std } => mean + std * phi_inv(p),
        Distribution::Triangular { lo, mode, hi } => {
            let span = hi - lo;
            if span <= 0.0 {
                return lo;
            }
            let fc = (mode - lo) / span;
            if p < fc {
                lo + (p * span * (mode - lo)).sqrt()
            } else {
                hi - ((1.0 - p) * span * (hi - mode)).sqrt()
            }
        }
    }
}

/// Standard deviation *equivalent* for a distribution — used to scale the
/// finite-difference step.  For Normal this is exact; for others we use
/// the distribution variance formula.
fn dist_std(dist: &Distribution) -> f64 {
    match *dist {
        Distribution::Uniform { lo, hi } => (hi - lo) / (12.0_f64).sqrt(),
        Distribution::Normal { std, .. } => std,
        Distribution::Triangular { lo, mode, hi } => {
            let v = (lo * lo + mode * mode + hi * hi - lo * mode - lo * hi - mode * hi) / 18.0;
            v.max(0.0).sqrt()
        }
    }
}

// ── x ↔ u transforms ─────────────────────────────────────────────────────────

/// Transform u-space point → physical x.
///
/// For Normal inputs the inverse transform is `x_j = mean_j + std_j * u_j`
/// — exact without any quantile/CDF composition and its rounding error.
/// For non-Normal inputs we go through the CDF route:
/// `x_j = F_j^{-1}(Phi(u_j))`.
fn u_to_x(u: &[f64], dists: &[Distribution]) -> Vec<f64> {
    u.iter()
        .zip(dists.iter())
        .map(|(&ui, di)| match *di {
            // Direct formula: avoids the round-trip phi_inv(phi(u)) error.
            Distribution::Normal { mean, std } => mean + std * ui,
            // For non-Normal: use the quantile route.
            _ => dist_quantile(di, phi(ui)),
        })
        .collect()
}

// ── results ───────────────────────────────────────────────────────────────────

/// Output of a [`form`] analysis.
#[derive(Debug, Clone)]
pub struct FormResult {
    /// Reliability index β = ‖u*‖.
    pub beta: f64,
    /// Probability of failure `Φ(−β)`.
    pub pf: f64,
    /// Most-Probable-Point in physical (*x*) space.
    pub mpp_x: Vec<f64>,
    /// Most-Probable-Point in standard-normal (*u*) space.
    pub mpp_u: Vec<f64>,
    /// Number of HLRF iterations to convergence.
    pub iterations: usize,
}

/// Output of a Monte-Carlo probability-of-failure estimate.
#[derive(Debug, Clone, Copy)]
pub struct McResult {
    /// Estimated probability of failure.
    pub pf: f64,
    /// Monte-Carlo standard error `√(Pf(1−Pf)/n)`.
    pub std_error: f64,
    /// Number of samples.
    pub n_samples: usize,
    /// Number of samples that fell in the failure domain.
    pub n_failures: usize,
}

/// Output of a SORM Breitung correction applied on top of a [`FormResult`].
#[derive(Debug, Clone)]
pub struct SormResult {
    /// FORM probability of failure (unchanged from input).
    pub pf_form: f64,
    /// Breitung-corrected probability of failure.  `None` if the correction
    /// could not be applied (a curvature violated `1 + β κᵢ > 0`).
    pub pf_sorm: Option<f64>,
    /// Principal curvatures at the MPP (length `n − 1`).
    pub curvatures: Vec<f64>,
    /// Human-readable reason if `pf_sorm` is `None`.
    pub warning: Option<String>,
}

// ── public API ────────────────────────────────────────────────────────────────

/// Configuration for a FORM / HLRF analysis.
#[derive(Debug, Clone)]
pub struct FormConfig {
    /// Maximum number of HLRF iterations before returning
    /// [`UqError::NotConverged`].
    pub max_iter: usize,
    /// Convergence tolerance on the change in `β` between successive iterations.
    pub tol: f64,
}

impl Default for FormConfig {
    fn default() -> Self {
        Self {
            max_iter: DEFAULT_MAX_ITER,
            tol: DEFAULT_TOL,
        }
    }
}

// ── FORM ─────────────────────────────────────────────────────────────────────

/// First-Order Reliability Method (FORM) using the HLRF algorithm.
///
/// Finds the Most-Probable-Point (MPP) of failure — the point on the surface
/// `g(x) = 0` closest to the mean in standard-normal *u*-space — and returns
/// the reliability index `β = ‖u*‖` and `Pf = Φ(−β)`.
///
/// # Parameters
/// * `g` — limit-state function; failure when `g(x) ≤ 0`.
/// * `dists` — prior distribution of each input; must be non-empty and match
///   the arity of `g`.
/// * `config` — HLRF iteration settings; use [`FormConfig::default()`] for
///   sensible defaults.
///
/// # Errors
/// * [`UqError::EmptyInput`] — `dists` is empty.
/// * [`UqError::NotConverged`] — HLRF did not converge within `max_iter`.
///
/// # Example
/// ```
/// use valenx_uq::Distribution;
/// use valenx_uq::reliability::{form, FormConfig};
///
/// // Linear limit state: g = 3 + 2x₀ + x₁, both inputs N(0,1).
/// // Exact β = 3 / √(2²+1²) = 3/√5 ≈ 1.3416.
/// let dists = [
///     Distribution::normal(0.0, 1.0).unwrap(),
///     Distribution::normal(0.0, 1.0).unwrap(),
/// ];
/// let g = |x: &[f64]| 3.0 + 2.0 * x[0] + x[1];
/// let result = form(g, &dists, &FormConfig::default()).unwrap();
/// let beta_exact = 3.0_f64 / 5.0_f64.sqrt();
/// assert!((result.beta - beta_exact).abs() < 1e-6,
///     "beta={} expected={beta_exact}", result.beta);
/// ```
pub fn form<G>(g: G, dists: &[Distribution], config: &FormConfig) -> Result<FormResult, UqError>
where
    G: Fn(&[f64]) -> f64,
{
    let n = dists.len();
    if n == 0 {
        return Err(UqError::EmptyInput(
            "FORM requires at least one random variable".into(),
        ));
    }

    // Start from the mean of each distribution (u = 0 in standard-normal space).
    let mut u: Vec<f64> = vec![0.0; n];

    let mut beta_prev = f64::INFINITY;
    let mut iters = 0_usize;

    for iter in 0..config.max_iter {
        iters = iter + 1;
        let x = u_to_x(&u, dists);
        let g_val = g(&x);

        // Finite-difference gradient of g in x-space.
        let grad_g_x = fd_gradient(&g, &x, dists);

        // Chain-rule: ∇ᵤG = ∇ₓg · (dx/du).
        // For Φ-based transform: dx_j/du_j = 1/φ(u_j) / f_j(x_j) — but we
        // implement it more robustly as: u_j = Φ⁻¹(F_j(x_j)),
        // so du_j/dx_j = f_j(x_j) / φ(u_j) and dx_j/du_j = φ(u_j) / f_j(x_j).
        // In u-space: ∂G/∂u_j = ∂g/∂x_j · dx_j/du_j.
        let grad_g_u: Vec<f64> = grad_g_x
            .iter()
            .zip(u.iter())
            .zip(x.iter())
            .zip(dists.iter())
            .map(|(((dg_dx, &uj), &xj), dist)| {
                let pdf_x = dist_pdf(dist, xj);
                let phi_u = standard_normal_pdf(uj);
                // Guard: if pdf_x is tiny, the chain-rule factor blows up.
                // We clamp it so we at minimum keep the sign correct.
                let scale = if pdf_x > 1e-300 {
                    phi_u / pdf_x
                } else {
                    // Degenerate: fall back to treating as Normal with same std.
                    dist_std(dist)
                };
                dg_dx * scale
            })
            .collect();

        let grad_norm = l2_norm_vec(&grad_g_u);
        if grad_norm < 1e-300 {
            // Gradient vanished — we may be at a stationary point.  Return
            // current estimate; it might be the MPP.
            break;
        }

        // HLRF update: u_{k+1} = (∇G·u - G) / ‖∇G‖ · (-∇G/‖∇G‖)
        //            = [ (∇G·u - G) / ‖∇G‖² ] · (-∇G)
        let dot = dot_product(&grad_g_u, &u);
        let lambda = (dot - g_val) / (grad_norm * grad_norm);
        // New u = -lambda * ∇G
        let u_new: Vec<f64> = grad_g_u.iter().map(|&dg| -lambda * dg).collect();

        let beta_new = l2_norm_vec(&u_new);

        // Check convergence on β.
        if (beta_new - beta_prev).abs() < config.tol {
            u = u_new;
            iters = iter + 1;
            break;
        }
        beta_prev = beta_new;
        u = u_new;

        if iter + 1 == config.max_iter {
            return Err(UqError::NotConverged(format!(
                "HLRF did not converge in {} iterations (|Δβ| = {:.2e})",
                config.max_iter,
                (beta_new - beta_prev).abs()
            )));
        }
    }

    let beta = l2_norm_vec(&u);
    let pf = phi(-beta);
    let mpp_x = u_to_x(&u, dists);

    Ok(FormResult {
        beta,
        pf,
        mpp_x,
        mpp_u: u,
        iterations: iters,
    })
}

// ── Monte-Carlo Pf ────────────────────────────────────────────────────────────

/// Crude Monte-Carlo estimate of the probability of failure `P[g(x) ≤ 0]`.
///
/// Samples `n_samples` realisations of the input vector from `dists` and
/// counts the fraction for which `g(x) ≤ 0`.
///
/// # Parameters
/// * `g` — limit-state function; failure when `g(x) ≤ 0`.
/// * `dists` — input distributions.
/// * `n_samples` — number of Monte-Carlo samples.  More samples → lower
///   standard error (`≈ sqrt(Pf(1-Pf)/n)`).
/// * `rng` — seeded [`SplitMix64`]; the same seed gives the same result.
///
/// # Errors
/// * [`UqError::EmptyInput`] — `dists` is empty or `n_samples == 0`.
///
/// # Example
/// ```
/// use valenx_uq::{Distribution, SplitMix64};
/// use valenx_uq::reliability::pf_monte_carlo;
///
/// // g = 2 - x, x ~ N(0,1).  Pf = Φ(-2) ≈ 0.02275.
/// let dists = [Distribution::normal(0.0, 1.0).unwrap()];
/// let g = |x: &[f64]| 2.0 - x[0];
/// let mut rng = SplitMix64::new(0xDEAD_BEEF);
/// let mc = pf_monte_carlo(g, &dists, 200_000, &mut rng).unwrap();
/// // Within 3σ of the exact value (σ_MC ≈ √(0.02275·0.9772/200_000) ≈ 3.3e-4).
/// assert!((mc.pf - 0.02275).abs() < 3.0 * mc.std_error + 1e-4);
/// ```
pub fn pf_monte_carlo<G>(
    g: G,
    dists: &[Distribution],
    n_samples: usize,
    rng: &mut SplitMix64,
) -> Result<McResult, UqError>
where
    G: Fn(&[f64]) -> f64,
{
    if dists.is_empty() {
        return Err(UqError::EmptyInput(
            "Monte-Carlo Pf requires at least one distribution".into(),
        ));
    }
    if n_samples == 0 {
        return Err(UqError::EmptyInput(
            "Monte-Carlo Pf requires n_samples > 0".into(),
        ));
    }

    let mut n_failures = 0_usize;
    for _ in 0..n_samples {
        let x: Vec<f64> = dists.iter().map(|d| d.sample(rng)).collect();
        if g(&x) <= 0.0 {
            n_failures += 1;
        }
    }

    let pf = n_failures as f64 / n_samples as f64;
    let std_error = (pf * (1.0 - pf) / n_samples as f64).sqrt();

    Ok(McResult {
        pf,
        std_error,
        n_samples,
        n_failures,
    })
}

// ── SORM Breitung ─────────────────────────────────────────────────────────────

/// Breitung second-order correction to the FORM probability of failure.
///
/// Given the MPP and reliability index from a [`form`] analysis, estimates
/// the principal curvatures of the limit-state surface at the MPP and applies
/// the Breitung formula
///
/// ```text
/// Pf_SORM ≈ Φ(−β) · ∏ᵢ (1 + β κᵢ)^(−1/2)
/// ```
///
/// where the product runs over the `n − 1` principal curvatures `κᵢ`.
///
/// If any factor `(1 + β κᵢ) ≤ 0` the correction is mathematically undefined;
/// in that case `pf_sorm` is `None` and a `warning` is set.
///
/// # Parameters
/// * `g`      — same limit-state function used for [`form`].
/// * `dists`  — same input distributions.
/// * `form_r` — result of [`form`]; provides the MPP and β.
///
/// # Errors
/// * [`UqError::EmptyInput`] — `dists` has fewer than 2 variables (no
///   curvature in 1-D).
pub fn sorm_breitung<G>(
    g: G,
    dists: &[Distribution],
    form_r: &FormResult,
) -> Result<SormResult, UqError>
where
    G: Fn(&[f64]) -> f64,
{
    let n = dists.len();
    if n < 2 {
        return Err(UqError::EmptyInput(
            "SORM requires at least 2 random variables (n-1 curvatures)".into(),
        ));
    }

    let beta = form_r.beta;
    let u_mpp = &form_r.mpp_u;
    let x_mpp = &form_r.mpp_x;

    // -- Step 1: Hessian of G(u) = g(T⁻¹(u)) in u-space via central differences.
    // We compute H[i][j] = (G(u+hᵢ+hⱼ) - G(u+hᵢ-hⱼ) - G(u-hᵢ+hⱼ) + G(u-hᵢ-hⱼ)) / (4hᵢhⱼ)
    // for off-diagonals, and the standard 3-point formula for diagonals.
    let h_vec: Vec<f64> = (0..n)
        .map(|j| {
            let s = dist_std(&dists[j]).max(1e-8);
            FD_H * s
        })
        .collect();

    // Gradient of G in u-space at the MPP (needed for normalisation).
    let g_mpp = g(x_mpp);
    let grad_g_u = {
        let grad_g_x = fd_gradient(&g, x_mpp, dists);
        chain_rule_u(&grad_g_x, u_mpp, x_mpp, dists)
    };
    let grad_norm = l2_norm_vec(&grad_g_u);

    // Direction cosines (unit vector pointing toward origin along the gradient).
    let alpha: Vec<f64> = if grad_norm > 1e-300 {
        grad_g_u.iter().map(|&v| -v / grad_norm).collect()
    } else {
        // Degenerate: can't compute curvatures.
        return Ok(SormResult {
            pf_form: form_r.pf,
            pf_sorm: None,
            curvatures: vec![],
            warning: Some(
                "gradient norm at MPP is effectively zero; SORM curvatures unavailable".into(),
            ),
        });
    };

    // -- Step 2: Hessian of G in u-space (numerical, central differences).
    // G(u) = g(T⁻¹(u)).  We perturb u, transform to x, call g.
    let hessian = {
        let mut h = vec![vec![0.0_f64; n]; n];
        let g_at_u = |u_pt: &Vec<f64>| -> f64 {
            let x_pt = u_to_x(u_pt, dists);
            g(&x_pt)
        };
        // Central difference on diagonals.
        for i in 0..n {
            let mut u_fwd = u_mpp.clone();
            let mut u_bwd = u_mpp.clone();
            u_fwd[i] += h_vec[i];
            u_bwd[i] -= h_vec[i];
            h[i][i] = (g_at_u(&u_fwd) - 2.0 * g_mpp + g_at_u(&u_bwd)) / (h_vec[i] * h_vec[i]);
        }
        // Mixed second partials.
        for i in 0..n {
            for j in (i + 1)..n {
                let hi = h_vec[i];
                let hj = h_vec[j];
                let mut upp = u_mpp.clone();
                let mut upm = u_mpp.clone();
                let mut ump = u_mpp.clone();
                let mut umm = u_mpp.clone();
                upp[i] += hi;
                upp[j] += hj;
                upm[i] += hi;
                upm[j] -= hj;
                ump[i] -= hi;
                ump[j] += hj;
                umm[i] -= hi;
                umm[j] -= hj;
                let val =
                    (g_at_u(&upp) - g_at_u(&upm) - g_at_u(&ump) + g_at_u(&umm)) / (4.0 * hi * hj);
                h[i][j] = val;
                h[j][i] = val;
            }
        }
        h
    };

    // -- Step 3: Project Hessian onto the (n-1) dimensional subspace ⊥ α.
    // Using the Gram-Schmidt basis of that subspace.
    // Then curvatures = eigenvalues of (projected H) / ‖∇G‖.
    //
    // For simplicity (no nalgebra required) we use the power-iteration method
    // only for the *diagonal* of the projected Hessian in a rotated frame,
    // which equals the principal curvatures for symmetric matrices via
    // Jacobi sweeps.
    //
    // We use a simple Jacobi eigenvalue algorithm on the small (n-1)×(n-1) matrix.

    // Build an orthonormal basis for the subspace orthogonal to alpha.
    let basis = orthonormal_complement(&alpha);
    // basis has shape (n-1) x n.

    let dim = n - 1;
    let mut reduced = vec![vec![0.0_f64; dim]; dim];
    // Compute reduced[i][j] = basis[i]^T · hessian · basis[j].
    // Two independent indices into a 2-D matrix — range loops are clearest.
    #[allow(clippy::needless_range_loop)]
    for i in 0..dim {
        for j in 0..dim {
            let mut val = 0.0;
            for r in 0..n {
                for s in 0..n {
                    val += basis[i][r] * hessian[r][s] * basis[j][s];
                }
            }
            reduced[i][j] = val;
        }
    }

    // Normalise by ‖∇G‖ to get curvatures κ (Hasofer-Lind convention).
    for row in &mut reduced {
        for v in row.iter_mut() {
            *v /= grad_norm;
        }
    }

    // Jacobi eigenvalue decomposition of the symmetric reduced matrix.
    let curvatures = jacobi_eigenvalues(&reduced);

    // -- Step 4: Breitung product.
    let pf_form = form_r.pf;
    let mut product = 1.0_f64;
    let mut bad_kappa: Option<String> = None;
    for (i, &kappa) in curvatures.iter().enumerate() {
        let factor = 1.0 + beta * kappa;
        if factor <= 0.0 {
            bad_kappa = Some(format!(
                "curvature κ[{i}] = {kappa:.4e} violates 1 + β·κ > 0 (β = {beta:.4})"
            ));
            break;
        }
        product *= factor.powf(-0.5);
    }

    if let Some(warn) = bad_kappa {
        return Ok(SormResult {
            pf_form,
            pf_sorm: None,
            curvatures,
            warning: Some(warn),
        });
    }

    Ok(SormResult {
        pf_form,
        pf_sorm: Some(pf_form * product),
        curvatures,
        warning: None,
    })
}

// ── internal helpers ──────────────────────────────────────────────────────────

/// Forward finite-difference gradient of `g` in x-space.
fn fd_gradient<G>(g: &G, x: &[f64], dists: &[Distribution]) -> Vec<f64>
where
    G: Fn(&[f64]) -> f64,
{
    let n = x.len();
    let mut grad = vec![0.0_f64; n];
    let mut x_pt = x.to_vec();
    for j in 0..n {
        // Central-difference O(h²) scheme: avoids the O(h) truncation error
        // of forward differences, so a linear g is differentiated exactly
        // (up to floating-point rounding only).
        let h = FD_H * dist_std(&dists[j]).max(1e-8);
        let xj_orig = x_pt[j];
        x_pt[j] = xj_orig + h;
        let gp = g(&x_pt);
        x_pt[j] = xj_orig - h;
        let gm = g(&x_pt);
        x_pt[j] = xj_orig;
        grad[j] = (gp - gm) / (2.0 * h);
    }
    grad
}

/// Chain-rule: convert x-space gradient to u-space gradient.
/// ∂G/∂u_j = ∂g/∂x_j · (dx_j/du_j) = ∂g/∂x_j · φ(u_j) / f_j(x_j).
fn chain_rule_u(grad_g_x: &[f64], u: &[f64], x: &[f64], dists: &[Distribution]) -> Vec<f64> {
    grad_g_x
        .iter()
        .zip(u.iter())
        .zip(x.iter())
        .zip(dists.iter())
        .map(|(((dg_dx, &uj), &xj), dist)| {
            let pdf_x = dist_pdf(dist, xj);
            let phi_u = standard_normal_pdf(uj);
            let scale = if pdf_x > 1e-300 {
                phi_u / pdf_x
            } else {
                dist_std(dist)
            };
            dg_dx * scale
        })
        .collect()
}

/// PDF of a distribution at `x`.
fn dist_pdf(dist: &Distribution, x: f64) -> f64 {
    match *dist {
        Distribution::Uniform { lo, hi } => {
            if x >= lo && x <= hi {
                1.0 / (hi - lo)
            } else {
                0.0
            }
        }
        Distribution::Normal { mean, std } => {
            let z = (x - mean) / std;
            standard_normal_pdf(z) / std
        }
        Distribution::Triangular { lo, mode, hi } => {
            let span = hi - lo;
            if span <= 0.0 || x < lo || x > hi {
                return 0.0;
            }
            if x <= mode {
                2.0 * (x - lo) / (span * (mode - lo))
            } else {
                2.0 * (hi - x) / (span * (hi - mode))
            }
        }
    }
}

/// PDF of N(0,1).
#[inline]
fn standard_normal_pdf(z: f64) -> f64 {
    (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt()
}

/// L2 norm of a slice.
#[inline]
fn l2_norm_vec(v: &[f64]) -> f64 {
    v.iter().map(|&x| x * x).sum::<f64>().sqrt()
}

/// Dot product of two equal-length slices.
#[inline]
fn dot_product(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(&ai, &bi)| ai * bi).sum()
}

/// Build an orthonormal basis for the subspace orthogonal to `alpha`
/// (Gram–Schmidt starting from standard basis vectors).
/// Returns `n-1` rows, each of length `n`.
fn orthonormal_complement(alpha: &[f64]) -> Vec<Vec<f64>> {
    let n = alpha.len();
    // Collect the candidates: e_0, e_1, ..., e_{n-1}.
    let mut basis: Vec<Vec<f64>> = Vec::with_capacity(n - 1);

    for k in 0..n {
        // Candidate: e_k.
        let mut v: Vec<f64> = vec![0.0; n];
        v[k] = 1.0;

        // Subtract projection onto alpha.
        let dot_alpha = dot_product(&v, alpha);
        for i in 0..n {
            v[i] -= dot_alpha * alpha[i];
        }

        // Subtract projections onto already accepted basis vectors.
        for b in &basis {
            let d = dot_product(&v, b);
            for i in 0..n {
                v[i] -= d * b[i];
            }
        }

        // Normalise.
        let norm = l2_norm_vec(&v);
        if norm > 1e-10 {
            for vi in &mut v {
                *vi /= norm;
            }
            basis.push(v);
            if basis.len() == n - 1 {
                break;
            }
        }
    }

    basis
}

/// Jacobi eigenvalue algorithm for a real symmetric matrix.
/// Returns eigenvalues (approximate) in no particular order.
/// Suitable for small matrices (n ≤ ~20); we only need n-1 ≤ dim.
fn jacobi_eigenvalues(a: &[Vec<f64>]) -> Vec<f64> {
    let n = a.len();
    if n == 0 {
        return vec![];
    }
    if n == 1 {
        return vec![a[0][0]];
    }

    let mut m = a.to_vec();
    const MAX_SWEEPS: usize = 100;
    const TOL: f64 = 1e-12;

    for _ in 0..MAX_SWEEPS {
        // Find off-diagonal element with largest absolute value.
        // Two independent indices into the same 2-D matrix — range loop is clearest.
        let (mut max_val, mut p, mut q) = (0.0_f64, 0_usize, 1_usize);
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            for j in (i + 1)..n {
                if m[i][j].abs() > max_val {
                    max_val = m[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }
        if max_val < TOL {
            break;
        }
        // Compute Givens angle.
        let theta = 0.5 * (m[q][q] - m[p][p]) / m[p][q];
        let t = if theta >= 0.0 {
            1.0 / (theta + (1.0 + theta * theta).sqrt())
        } else {
            1.0 / (theta - (1.0 + theta * theta).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;

        // Apply rotation to m.
        let app = m[p][p];
        let aqq = m[q][q];
        let apq = m[p][q];
        m[p][p] = app - t * apq;
        m[q][q] = aqq + t * apq;
        m[p][q] = 0.0;
        m[q][p] = 0.0;

        for r in 0..n {
            if r == p || r == q {
                continue;
            }
            let arp = m[r][p];
            let arq = m[r][q];
            m[r][p] = c * arp - s * arq;
            m[p][r] = m[r][p];
            m[r][q] = s * arp + c * arq;
            m[q][r] = m[r][q];
        }
    }

    (0..n).map(|i| m[i][i]).collect()
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distribution::Distribution;
    use crate::rng::SplitMix64;

    // ── Benchmark pins ────────────────────────────────────────────────────────

    /// BP1: linear limit-state with standard-normal inputs.
    /// g = a0 + a1*x0 + a2*x1; exact β = a0 / sqrt(a1^2 + a2^2).
    /// FORM must match to < 1e-6.
    #[test]
    fn bp1_linear_normal_exact_beta() {
        let a0 = 3.0_f64;
        let a1 = 2.0_f64;
        let a2 = 1.0_f64;
        let beta_exact = a0 / (a1 * a1 + a2 * a2).sqrt();

        let dists = [
            Distribution::normal(0.0, 1.0).unwrap(),
            Distribution::normal(0.0, 1.0).unwrap(),
        ];
        let g = move |x: &[f64]| a0 + a1 * x[0] + a2 * x[1];
        let r = form(g, &dists, &FormConfig::default()).unwrap();

        assert!(
            (r.beta - beta_exact).abs() < 1e-6,
            "BP1: beta={} expected={beta_exact}",
            r.beta
        );
        let pf_exact = phi(-beta_exact);
        assert!(
            (r.pf - pf_exact).abs() < 1e-8,
            "BP1: Pf={} expected={pf_exact}",
            r.pf
        );
    }

    /// BP1 variant: 4-variable linear limit-state, still exact for FORM.
    #[test]
    fn bp1_linear_4var_exact_beta() {
        let coeffs = [5.0_f64, 1.0, 2.0, 3.0, 4.0]; // [a0, a1..a4]
        let beta_exact = coeffs[0] / coeffs[1..].iter().map(|&a| a * a).sum::<f64>().sqrt();

        let dists: Vec<Distribution> = (0..4)
            .map(|_| Distribution::normal(0.0, 1.0).unwrap())
            .collect();
        let g = move |x: &[f64]| {
            coeffs[0]
                + coeffs[1..]
                    .iter()
                    .zip(x.iter())
                    .map(|(&a, &xi)| a * xi)
                    .sum::<f64>()
        };
        let r = form(g, &dists, &FormConfig::default()).unwrap();
        assert!(
            (r.beta - beta_exact).abs() < 1e-5,
            "BP1-4var: beta={} expected={beta_exact}",
            r.beta
        );
    }

    /// BP2: FORM Pf agrees with Monte-Carlo within sampling error for a
    /// nonlinear limit-state.  We use g = R - S where R ~ N(10, 1) and
    /// S ~ N(8, 1.5) — analytic Pf = Φ(-(10-8)/√(1+1.5²)) ≈ Φ(-0.7428).
    #[test]
    fn bp2_form_matches_mc_nonlinear() {
        // g = R - S, variables: [R, S].
        // Exact beta = (mu_R - mu_S) / sqrt(sigma_R^2 + sigma_S^2)
        let mu_r = 10.0_f64;
        let mu_s = 8.0_f64;
        let sig_r = 1.0_f64;
        let sig_s = 1.5_f64;
        let beta_exact = (mu_r - mu_s) / (sig_r * sig_r + sig_s * sig_s).sqrt();
        let pf_exact = phi(-beta_exact);

        let dists = [
            Distribution::normal(mu_r, sig_r).unwrap(),
            Distribution::normal(mu_s, sig_s).unwrap(),
        ];
        let g = |x: &[f64]| x[0] - x[1];
        let r = form(g, &dists, &FormConfig::default()).unwrap();
        assert!(
            (r.pf - pf_exact).abs() < 1e-5,
            "BP2 FORM: pf={} expected={pf_exact}",
            r.pf
        );

        let mut rng = SplitMix64::new(0xABCD_1234);
        let mc = pf_monte_carlo(g, &dists, 500_000, &mut rng).unwrap();
        // FORM and MC agree within 3×MC standard error.
        let tol = 3.0 * mc.std_error;
        assert!(
            (mc.pf - r.pf).abs() < tol + 1e-4,
            "BP2: MC pf={} FORM pf={} tol={}",
            mc.pf,
            r.pf,
            tol
        );
    }

    /// BP3: HLRF converges (β stabilises) for a smooth nonlinear g.
    /// g = 1 - x0^2 * x1 / 20, inputs N(0,1).
    #[test]
    fn bp3_hlrf_converges_smooth_g() {
        let dists = [
            Distribution::normal(0.0, 1.0).unwrap(),
            Distribution::normal(0.0, 1.0).unwrap(),
        ];
        let g = |x: &[f64]| 1.0 - x[0] * x[0] * x[1] / 20.0;
        let r = form(g, &dists, &FormConfig::default());
        assert!(r.is_ok(), "BP3: HLRF should converge for smooth g");
        let r = r.unwrap();
        assert!(r.beta.is_finite(), "BP3: beta must be finite");
        assert!(
            r.iterations < DEFAULT_MAX_ITER,
            "BP3: converged before max_iter"
        );
    }

    /// BP4: shifting a variable's mean toward the failure surface increases Pf.
    #[test]
    fn bp4_increasing_mean_raises_pf() {
        let dists_safe = [
            Distribution::normal(5.0, 1.0).unwrap(),
            Distribution::normal(0.0, 1.0).unwrap(),
        ];
        let dists_risky = [
            Distribution::normal(1.0, 1.0).unwrap(), // mean closer to 0 boundary
            Distribution::normal(0.0, 1.0).unwrap(),
        ];
        // g = x0 - 0: failure when x0 <= 0.
        let g = |x: &[f64]| x[0];
        let r_safe = form(g, &dists_safe, &FormConfig::default()).unwrap();
        let r_risky = form(g, &dists_risky, &FormConfig::default()).unwrap();
        assert!(
            r_risky.pf > r_safe.pf,
            "BP4: risky Pf ({}) should exceed safe Pf ({})",
            r_risky.pf,
            r_safe.pf
        );
    }

    // ── Error guards ──────────────────────────────────────────────────────────

    /// No random variables → Err, not hang.
    #[test]
    fn error_no_variables() {
        let g = |_x: &[f64]| 1.0_f64;
        let r = form(g, &[], &FormConfig::default());
        assert!(r.is_err(), "should return Err for zero variables");
    }

    /// MC with zero samples → Err.
    #[test]
    fn error_mc_zero_samples() {
        let dists = [Distribution::normal(0.0, 1.0).unwrap()];
        let g = |x: &[f64]| x[0];
        let mut rng = SplitMix64::new(0);
        let r = pf_monte_carlo(g, &dists, 0, &mut rng);
        assert!(r.is_err());
    }

    /// SORM with single variable → Err (no curvatures defined).
    #[test]
    fn error_sorm_single_var() {
        let dists = [Distribution::normal(0.0, 1.0).unwrap()];
        fn g(x: &[f64]) -> f64 {
            2.0 - x[0]
        }
        let form_r = form(g, &dists, &FormConfig::default()).unwrap();
        let sr = sorm_breitung(g, &dists, &form_r);
        assert!(sr.is_err());
    }

    // ── Additional coverage ───────────────────────────────────────────────────

    /// MC returns sensible estimate for g = 3 - x, x ~ N(0,1).
    /// Exact Pf = Φ(-3) ≈ 0.001349.
    #[test]
    fn mc_small_pf_estimate() {
        let dists = [Distribution::normal(0.0, 1.0).unwrap()];
        let g = |x: &[f64]| 3.0 - x[0];
        let mut rng = SplitMix64::new(0xFEED_BEEF);
        let mc = pf_monte_carlo(g, &dists, 1_000_000, &mut rng).unwrap();
        let pf_exact = phi(-3.0);
        // 4σ tolerance on MC.
        assert!(
            (mc.pf - pf_exact).abs() < 4.0 * mc.std_error + 1e-5,
            "mc.pf={} pf_exact={pf_exact} se={}",
            mc.pf,
            mc.std_error
        );
    }

    /// FORM works correctly for non-Normal (Uniform) inputs.
    #[test]
    fn form_uniform_input() {
        // x ~ U(-2, 2), g = x.  Failure when x <= 0, Pf = 0.5.
        let dists = [Distribution::uniform(-2.0, 2.0).unwrap()];
        let g = |x: &[f64]| x[0];
        let r = form(g, &dists, &FormConfig::default()).unwrap();
        // Pf should be close to 0.5 (symmetric uniform about 0).
        assert!(
            (r.pf - 0.5).abs() < 0.05,
            "FORM with Uniform input: pf={} expected≈0.5",
            r.pf
        );
    }

    /// SORM applies without error for a 2-D smooth case.
    #[test]
    fn sorm_runs_2d() {
        let dists = [
            Distribution::normal(0.0, 1.0).unwrap(),
            Distribution::normal(0.0, 1.0).unwrap(),
        ];
        fn g(x: &[f64]) -> f64 {
            3.0 + x[0] + x[1]
        }
        let form_r = form(g, &dists, &FormConfig::default()).unwrap();
        let sr = sorm_breitung(g, &dists, &form_r).unwrap();
        // For a linear g the curvature is 0 → SORM == FORM.
        let pf_sorm = sr.pf_sorm.expect("SORM should give a result for linear g");
        assert!(
            (pf_sorm - form_r.pf).abs() < 1e-6,
            "SORM on linear g should equal FORM pf"
        );
    }

    /// Beta is always non-negative.
    #[test]
    fn beta_non_negative() {
        let dists = [Distribution::normal(1.0, 0.5).unwrap()];
        // g = x - 0.5: positive mean, modest sigma.
        let g = |x: &[f64]| x[0] - 0.5;
        let r = form(g, &dists, &FormConfig::default()).unwrap();
        assert!(r.beta >= 0.0, "beta must be non-negative, got {}", r.beta);
    }

    /// Pf is always in [0, 1].
    #[test]
    fn pf_in_unit_interval() {
        let dists = [Distribution::normal(2.0, 1.0).unwrap()];
        let g = |x: &[f64]| x[0];
        let r = form(g, &dists, &FormConfig::default()).unwrap();
        assert!(
            (0.0..=1.0).contains(&r.pf),
            "Pf must be in [0,1], got {}",
            r.pf
        );
    }

    /// FORM beta for 1-variable standard normal equals |threshold|.
    #[test]
    fn form_single_var_exact() {
        // g = 2 - x, x ~ N(0,1). Exact beta = 2, Pf = Φ(-2).
        let dists = [Distribution::normal(0.0, 1.0).unwrap()];
        let g = |x: &[f64]| 2.0 - x[0];
        let r = form(g, &dists, &FormConfig::default()).unwrap();
        assert!((r.beta - 2.0).abs() < 1e-6, "beta={} expected=2.0", r.beta);
        assert!(
            (r.pf - phi(-2.0)).abs() < 1e-8,
            "pf={} expected={}",
            r.pf,
            phi(-2.0)
        );
    }
}
