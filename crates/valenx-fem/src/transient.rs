//! **Blast / shock survivability** — transient structural *response* of a
//! linear system to a fast transient (blast, shock, impact) load, by
//! Newmark-β implicit time integration on **supplied** mass, damping, and
//! stiffness matrices.
//!
//! ## Defensive scope — survivability / protection only
//!
//! This module answers one question: **given a structure (its `M`, `C`, `K`)
//! and a transient load `F(t)`, how does the structure respond, and does it
//! survive?** It is the structural-dynamics core of *protective* /
//! survivability design — the same physics as civil blast-resistant design
//! (ASCE 59, UFC 3-340-02 single-degree-of-freedom analysis) and automotive
//! crash-pulse response. It models how a structure *takes* a load: the
//! displacement, velocity, and acceleration histories, the peak response, and
//! the dynamic amplification.
//!
//! It is deliberately **not** a weapon, warhead, penetration, fragmentation,
//! or lethality model. The blast *source* is represented only by the standard
//! idealized far-field overpressure-time curve ([`FriedlanderPulse`]) used as a
//! protective design load; nothing here models how a charge produces that
//! pressure, nor how a structure is defeated. The output is a survivability
//! *screen*: "what is the peak deflection / acceleration this protective
//! structure sees, and how does it compare to its static capacity?"
//!
//! ## What this is, and how it relates to [`crate::dynamics`]
//!
//! [`crate::dynamics::solve_transient_dynamics`] is the **full-FEM** transient
//! solver: it takes a tetrahedral [`valenx_mesh::Mesh`] plus a material,
//! assembles the consistent mass and elastic stiffness itself, and marches a
//! **constant-in-time** load. That is the right tool when you have a meshed
//! solid and a step load.
//!
//! This module is the complementary **reduced / supplied-matrix** survivability
//! analyzer. It does **not** assemble anything from a mesh. The caller supplies
//! the system matrices directly:
//!
//! - a **single-degree-of-freedom (SDOF)** equivalent — the workhorse of
//!   blast-resistant design, where a wall/panel/frame is reduced to one
//!   effective mass `m`, stiffness `k`, and damping `c` (transformation
//!   factors are the engineer's input, not computed here); or
//! - a small **multi-degree-of-freedom (MDOF)** model assembled by the caller
//!   (for instance the reduced free-DOF `K`/`M` recovered elsewhere in this
//!   crate, or a hand-built lumped-mass chain).
//!
//! and a **time-varying** load `F(t)` (e.g. the Friedlander blast pulse). The
//! payoff is twofold: (1) the load may vary arbitrarily in time, which the
//! mesh-based solver's constant load cannot express, and (2) it works on any
//! matrices, so an SDOF blast screen needs no mesh at all. We say this plainly:
//! **this is a supplied-`M`/`C`/`K` SDOF/MDOF analyzer, not a full-FEM
//! assembler.**
//!
//! ## Newmark-β time integration
//!
//! The semidiscrete equation of motion
//!
//! ```text
//!   M·ü + C·u̇ + K·u = F(t)
//! ```
//!
//! is marched by the **Newmark-β** method (Newmark 1959). Over one step `Δt`,
//!
//! ```text
//!   u_{n+1} = u_n + Δt·u̇_n + Δt²·[(½−β)·ü_n + β·ü_{n+1}]
//!   u̇_{n+1} = u̇_n + Δt·[(1−γ)·ü_n + γ·ü_{n+1}]
//! ```
//!
//! Substituting into the equation of motion at `t_{n+1}` gives a linear system
//! in `u_{n+1}` with the **effective stiffness**
//!
//! ```text
//!   K_eff = K + (γ/(β·Δt))·C + (1/(β·Δt²))·M
//! ```
//!
//! and an effective load built from the previous step's state **plus the new
//! external load `F(t_{n+1})`** — this is what lets the load vary in time.
//! `K_eff` is constant for a linear system with a fixed step, so it is
//! factorised **once** (LU) and only back-substituted each step.
//!
//! With the **average-acceleration** parameters `β = 1/4, γ = 1/2`
//! ([`NewmarkBeta::average_acceleration`], the default) the scheme is
//! **unconditionally stable** and **second-order accurate**, and conserves
//! energy for an undamped system — `Δt` is an *accuracy* choice, not a
//! stability bound. For a blast/shock pulse, resolve the rise and the
//! structural period: `Δt ≤ min(t_rise, T)/20` is a safe rule of thumb.
//!
//! ## Honest scope — linear, small-deflection
//!
//! This is a **linear-elastic, small-deflection** survivability *screen*:
//!
//! - `K`, `C`, `M` are **constant** — no plasticity (no resistance-function
//!   elastic-plastic clamp), no geometric (large-deflection) stiffening, no
//!   contact, no fracture. A real protective element loaded into the plastic
//!   range absorbs far more energy than the linear model predicts, so a linear
//!   peak deflection is **conservative for elastic response and unconservative
//!   as a ductility/energy estimate** past yield. For the elastic-plastic SDOF
//!   resistance function and pressure-impulse (P-I) diagrams of full
//!   blast-resistant practice, a nonlinear resistance model is the follow-up.
//! - It is therefore a *preliminary-design / screening* tool, **not** a
//!   nonlinear hydrocode (LS-DYNA, Autodyn) and not a substitute for one.
//! - The Friedlander curve is the standard *idealized* free-field pulse; real
//!   reflected-pressure loading, clearing, and confinement are not modelled.
//!
//! Each limitation is a well-understood extension; none affects the
//! correctness of the linear transient response computed here, which is
//! verified against closed-form SDOF solutions in the [tests](self).

use std::f64::consts::PI;

use nalgebra::{DMatrix, DVector};
use thiserror::Error;

/// Errors from a transient survivability solve. Every degenerate or
/// non-physical input fails **loud** with one of these — the integrator never
/// panics and never returns a silently-corrupt `Ok`.
#[derive(Debug, Error, PartialEq)]
pub enum SurvivabilityError {
    /// A supplied system matrix was empty (zero rows) — there is nothing to
    /// integrate.
    #[error("system matrices are empty (0 DOFs); supply at least a 1×1 SDOF system")]
    EmptySystem,
    /// The three system matrices (`M`, `C`, `K`) are not all the same square
    /// `n × n` shape.
    #[error(
        "M, C, K must be the same square shape: got M {m_rows}×{m_cols}, \
         C {c_rows}×{c_cols}, K {k_rows}×{k_cols}"
    )]
    ShapeMismatch {
        /// Rows of `M`.
        m_rows: usize,
        /// Columns of `M`.
        m_cols: usize,
        /// Rows of `C`.
        c_rows: usize,
        /// Columns of `C`.
        c_cols: usize,
        /// Rows of `K`.
        k_rows: usize,
        /// Columns of `K`.
        k_cols: usize,
    },
    /// The mass matrix is not positive-definite (a zero or negative effective
    /// mass). Dynamics needs `M` SPD — the initial acceleration solve
    /// `M·ü₀ = …` is otherwise singular or yields a non-physical sign.
    #[error("mass matrix is not positive-definite (zero/negative effective mass): {0}")]
    NonPositiveMass(String),
    /// The stiffness matrix is not positive-definite (a zero or negative
    /// effective stiffness — an unrestrained / mechanism system). A
    /// survivability screen needs a restrained, stable structure.
    #[error("stiffness matrix is not positive-definite (unrestrained or non-positive stiffness)")]
    NonPositiveStiffness,
    /// The effective Newmark stiffness `K_eff` could not be factorised — the
    /// supplied matrices do not form a solvable linear system.
    #[error("Newmark effective-stiffness factorisation failed (system is singular)")]
    SolveFailed,
    /// A time-step, Newmark parameter, or step count was non-physical.
    #[error("invalid parameter: {0}")]
    InvalidParameter(String),
    /// A supplied matrix entry, initial-state value, or evaluated load was not
    /// finite (`NaN` or `±∞`). A non-finite value flows through every step and
    /// silently corrupts the whole history, so it is rejected up front.
    #[error("non-finite value in {what}")]
    NonFinite {
        /// Where the offending value was found (matrix / load / initial state).
        what: String,
    },
}

/// Upper bound on the number of integration steps a single
/// [`solve_transient_response`] call will accept, mirroring
/// [`crate::dynamics::MAX_TIME_STEPS`]. Guards against a runaway or hostile
/// `n_steps` OOM-ing the host inside `Vec::with_capacity`.
pub const MAX_TIME_STEPS: usize = 1_000_000;

/// The two Newmark-β integration parameters `β` and `γ`.
///
/// [`NewmarkBeta::average_acceleration`] (`β = 1/4, γ = 1/2`) is the default
/// and the recommended choice — it is the only *unconditionally stable* member
/// of the family, which is what a blast/shock screen wants.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct NewmarkBeta {
    /// Newmark `β` — how the acceleration is assumed to vary across a step.
    pub beta: f64,
    /// Newmark `γ` — the algorithmic damping. `1/2` gives none (second-order,
    /// energy-conserving); `γ > 1/2` damps spurious high modes at first order.
    pub gamma: f64,
}

impl NewmarkBeta {
    /// The **average-acceleration** ("trapezoidal", constant-average-
    /// acceleration) scheme — `β = 1/4, γ = 1/2`. Unconditionally stable,
    /// second-order accurate, energy-conserving for an undamped system. The
    /// standard choice for survivability analysis.
    pub fn average_acceleration() -> NewmarkBeta {
        NewmarkBeta {
            beta: 0.25,
            gamma: 0.5,
        }
    }

    /// The **linear-acceleration** scheme — `β = 1/6, γ = 1/2`. Marginally more
    /// accurate but only *conditionally* stable; offered for completeness.
    pub fn linear_acceleration() -> NewmarkBeta {
        NewmarkBeta {
            beta: 1.0 / 6.0,
            gamma: 0.5,
        }
    }
}

impl Default for NewmarkBeta {
    fn default() -> Self {
        NewmarkBeta::average_acceleration()
    }
}

/// The standard idealized **Friedlander** decaying-exponential blast
/// overpressure pulse — the canonical free-field protective design load.
///
/// After the shock front arrives, the overpressure jumps to its peak `p₀` and
/// then decays through the **positive phase** of duration `t_d`, crossing zero
/// at `t = t_d` and going slightly negative (the suction phase) before
/// returning to ambient. The Friedlander form captures this with a single
/// dimensionless **waveform (decay) parameter** `b`:
///
/// ```text
///   p(τ) = p₀ · (1 − τ/t_d) · exp(−b·τ/t_d),   τ = t − t_arrival ≥ 0,  τ ≤ ?
/// ```
///
/// where `τ` is time measured from the arrival of the front. For `τ < 0`
/// (before arrival) the overpressure is zero. The `(1 − τ/t_d)` factor makes
/// the curve cross zero exactly at the end of the positive phase `τ = t_d`;
/// past that the same expression is negative (the suction phase) and decays
/// back to zero.
///
/// This is the overpressure curve specified for blast-resistant design in
/// UFC 3-340-02 / ASCE references. It is an *idealization of the load*, not a
/// model of the explosive source.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FriedlanderPulse {
    /// Peak overpressure `p₀` at the shock front (pressure units, e.g. Pa, or
    /// force if you have already lumped the loaded area in). Must be finite.
    pub peak_overpressure: f64,
    /// Positive-phase duration `t_d` (seconds) — the time from arrival to the
    /// first zero crossing. Must be finite and strictly positive.
    pub positive_duration: f64,
    /// Dimensionless waveform / decay parameter `b` (`> 0`). Larger `b` means a
    /// faster exponential fall-off within the positive phase. A common default
    /// for far-field air blast is `b ≈ 1`.
    pub decay: f64,
    /// Time of arrival `t_arrival` (seconds) of the shock front. The pulse is
    /// zero before this. Defaults to `0`.
    pub arrival_time: f64,
}

impl FriedlanderPulse {
    /// Build a Friedlander pulse arriving at `t = 0`.
    ///
    /// # Errors
    ///
    /// [`SurvivabilityError::InvalidParameter`] if `peak` is not finite,
    /// `positive_duration` is not finite-and-positive, or `decay` is not
    /// finite-and-positive.
    pub fn new(
        peak_overpressure: f64,
        positive_duration: f64,
        decay: f64,
    ) -> Result<FriedlanderPulse, SurvivabilityError> {
        FriedlanderPulse {
            peak_overpressure,
            positive_duration,
            decay,
            arrival_time: 0.0,
        }
        .validated()
    }

    /// Set the shock-front arrival time (seconds), returning the updated pulse.
    pub fn with_arrival(mut self, arrival_time: f64) -> FriedlanderPulse {
        self.arrival_time = arrival_time;
        self
    }

    /// Validate the parameters, returning `self` if they are physical.
    fn validated(self) -> Result<FriedlanderPulse, SurvivabilityError> {
        if !self.peak_overpressure.is_finite() {
            return Err(SurvivabilityError::InvalidParameter(format!(
                "peak overpressure must be finite, got {}",
                self.peak_overpressure
            )));
        }
        if !(self.positive_duration.is_finite() && self.positive_duration > 0.0) {
            return Err(SurvivabilityError::InvalidParameter(format!(
                "positive-phase duration must be finite and > 0, got {}",
                self.positive_duration
            )));
        }
        if !(self.decay.is_finite() && self.decay > 0.0) {
            return Err(SurvivabilityError::InvalidParameter(format!(
                "decay parameter b must be finite and > 0, got {}",
                self.decay
            )));
        }
        if !self.arrival_time.is_finite() {
            return Err(SurvivabilityError::InvalidParameter(format!(
                "arrival time must be finite, got {}",
                self.arrival_time
            )));
        }
        Ok(self)
    }

    /// The overpressure at absolute time `t` (seconds). Zero before the front
    /// arrives; the Friedlander decay afterwards (including the negative
    /// suction phase past the positive duration).
    pub fn pressure_at(&self, t: f64) -> f64 {
        let tau = t - self.arrival_time;
        if tau < 0.0 {
            return 0.0;
        }
        let x = tau / self.positive_duration;
        self.peak_overpressure * (1.0 - x) * (-self.decay * x).exp()
    }

    /// The **positive-phase specific impulse** `i = ∫₀^{t_d} p(τ) dτ` — the
    /// area under the overpressure curve over the positive phase only (arrival
    /// to first zero crossing). This is the load's impulse, the second key
    /// blast parameter alongside the peak.
    ///
    /// Closed form of `∫₀^{t_d} p₀ (1 − τ/t_d) e^{−bτ/t_d} dτ`, substituting
    /// `x = τ/t_d`:
    ///
    /// ```text
    ///   i = p₀ · t_d · [ 1/b − (1 − e^{−b})/b² ]
    /// ```
    ///
    /// (For small `b` this tends to the triangular-pulse limit `p₀·t_d/2`.)
    pub fn positive_impulse(&self) -> f64 {
        let b = self.decay;
        let p0 = self.peak_overpressure;
        let td = self.positive_duration;
        p0 * td * (1.0 / b - (1.0 - (-b).exp()) / (b * b))
    }
}

/// Build a Rayleigh proportional damping matrix `C = α·M + β·K` from supplied
/// mass and stiffness matrices.
///
/// A convenience for callers who want the standard proportional damping form
/// but already hold dense `M`/`K`. (For *fitting* `α`, `β` to modal damping
/// targets, see [`crate::rayleigh::RayleighDamping`].) Returns an error if the
/// shapes disagree.
pub fn rayleigh_damping_matrix(
    mass: &DMatrix<f64>,
    stiffness: &DMatrix<f64>,
    alpha: f64,
    beta: f64,
) -> Result<DMatrix<f64>, SurvivabilityError> {
    if mass.shape() != stiffness.shape() {
        return Err(SurvivabilityError::ShapeMismatch {
            m_rows: mass.nrows(),
            m_cols: mass.ncols(),
            c_rows: mass.nrows(),
            c_cols: mass.ncols(),
            k_rows: stiffness.nrows(),
            k_cols: stiffness.ncols(),
        });
    }
    Ok(mass * alpha + stiffness * beta)
}

/// Controls for a Newmark-β transient survivability run.
#[derive(Copy, Clone, Debug)]
pub struct TransientControls {
    /// The time step `Δt` (seconds). With average-acceleration Newmark this is
    /// an *accuracy* choice; resolve the load rise and the structural period
    /// (`Δt ≤ min(t_rise, T)/20`).
    pub dt: f64,
    /// How many steps to march. The run covers `[0, n_steps·Δt]`.
    pub n_steps: usize,
    /// The Newmark parameters.
    pub newmark: NewmarkBeta,
}

impl Default for TransientControls {
    fn default() -> Self {
        TransientControls {
            dt: 1.0e-4,
            n_steps: 2000,
            newmark: NewmarkBeta::average_acceleration(),
        }
    }
}

/// The full transient response history of a survivability solve.
///
/// `time[i]`, `displacement[i]`, `velocity[i]`, and `acceleration[i]` are the
/// state at step `i`. Index `0` is the **initial state** at `t = 0`; there are
/// `n_steps + 1` recorded instants in total (the initial state plus one per
/// step).
#[derive(Clone, Debug)]
pub struct TransientResponse {
    /// Physical time `t` (seconds) at each recorded instant.
    pub time: Vec<f64>,
    /// Per-DOF displacement vector at each instant.
    pub displacement: Vec<DVector<f64>>,
    /// Per-DOF velocity vector at each instant.
    pub velocity: Vec<DVector<f64>>,
    /// Per-DOF acceleration vector at each instant.
    pub acceleration: Vec<DVector<f64>>,
}

impl TransientResponse {
    /// The displacement time history of one DOF — `(time, displacement)` pairs.
    /// Empty if `dof` is out of range.
    pub fn dof_history(&self, dof: usize) -> Vec<(f64, f64)> {
        self.time
            .iter()
            .zip(&self.displacement)
            .filter_map(|(&t, u)| u.get(dof).map(|&d| (t, d)))
            .collect()
    }

    /// The **peak absolute displacement** of one DOF over the whole history —
    /// the headline survivability readout (the maximum deflection the
    /// structure reaches). `0.0` if `dof` is out of range or the history is
    /// empty.
    pub fn peak_abs_displacement(&self, dof: usize) -> f64 {
        self.displacement
            .iter()
            .filter_map(|u| u.get(dof))
            .map(|d| d.abs())
            .fold(0.0, f64::max)
    }

    /// The **peak absolute acceleration** of one DOF over the whole history —
    /// the readout that drives inertial / shock damage and occupant g-loading.
    /// `0.0` if `dof` is out of range or the history is empty.
    pub fn peak_abs_acceleration(&self, dof: usize) -> f64 {
        self.acceleration
            .iter()
            .filter_map(|a| a.get(dof))
            .map(|a| a.abs())
            .fold(0.0, f64::max)
    }

    /// The **dynamic amplification factor (DAF)** at one DOF: the ratio of the
    /// peak *dynamic* displacement to the supplied static reference deflection
    /// `static_deflection`. A DAF of `1` means the dynamic peak equals the
    /// static deflection; for a suddenly-applied (step) load on an undamped
    /// SDOF the DAF is exactly `2` (the structure overshoots to twice static).
    ///
    /// The caller supplies the static reference (typically `F_peak / k`, the
    /// deflection the *peak* load would cause if applied slowly). Returns
    /// `None` if `static_deflection` is zero or non-finite (no meaningful
    /// ratio).
    pub fn dynamic_amplification_factor(&self, dof: usize, static_deflection: f64) -> Option<f64> {
        if !static_deflection.is_finite() || static_deflection == 0.0 {
            return None;
        }
        Some(self.peak_abs_displacement(dof) / static_deflection.abs())
    }
}

/// Solve the linear transient survivability problem
/// `M·ü + C·u̇ + K·u = F(t)` for the full displacement / velocity /
/// acceleration history, by Newmark-β implicit time integration on the
/// **supplied** dense matrices.
///
/// This is the supplied-matrix SDOF/MDOF analyzer described in the [module
/// docs](self) — it assembles nothing from a mesh. `load(t)` is a closure
/// returning the external force vector (length `n`, the DOF count) at absolute
/// time `t`; use [`FriedlanderPulse`] to drive it with a blast overpressure, or
/// any other transient. `u0` / `v0` are the initial displacement / velocity
/// (each length `n`; pass freshly-zeroed vectors to start from rest).
///
/// Returns the [`TransientResponse`] — `n_steps + 1` recorded instants
/// including the initial state.
///
/// # Method
///
/// 1. Validate shapes / finiteness; check `M` and `K` are positive-definite
///    (Cholesky).
/// 2. Recover the initial acceleration from `M·ü₀ = F(0) − C·u̇₀ − K·u₀`.
/// 3. Form the constant Newmark effective stiffness `K_eff = K + a₀·M + a₁·C`
///    and factorise it **once** (LU).
/// 4. Each step: evaluate `F(t_{n+1})`, build the effective load from the
///    previous state, back-substitute for `u_{n+1}`, and update velocity and
///    acceleration from the Newmark relations.
///
/// # Errors
///
/// See [`SurvivabilityError`] — empty/mismatched matrices, a non-finite entry,
/// a non-positive-definite `M` or `K`, a bad time step / step count, or a
/// factorisation failure. A `load(t)` that returns the wrong length or a
/// non-finite force is rejected on the first evaluation.
pub fn solve_transient_response<F>(
    mass: &DMatrix<f64>,
    damping: &DMatrix<f64>,
    stiffness: &DMatrix<f64>,
    u0: &DVector<f64>,
    v0: &DVector<f64>,
    mut load: F,
    controls: &TransientControls,
) -> Result<TransientResponse, SurvivabilityError>
where
    F: FnMut(f64) -> DVector<f64>,
{
    // --- step-count guard (mirror dynamics.rs DoS cap) ---
    if controls.n_steps > MAX_TIME_STEPS {
        return Err(SurvivabilityError::InvalidParameter(format!(
            "n_steps {} exceeds the {MAX_TIME_STEPS}-step cap (lower the step count or coarsen Δt)",
            controls.n_steps
        )));
    }

    let n = mass.nrows();
    if n == 0 || mass.ncols() == 0 {
        return Err(SurvivabilityError::EmptySystem);
    }

    // --- shape agreement: M, C, K all n×n ---
    if mass.nrows() != mass.ncols()
        || damping.nrows() != n
        || damping.ncols() != n
        || stiffness.nrows() != n
        || stiffness.ncols() != n
    {
        return Err(SurvivabilityError::ShapeMismatch {
            m_rows: mass.nrows(),
            m_cols: mass.ncols(),
            c_rows: damping.nrows(),
            c_cols: damping.ncols(),
            k_rows: stiffness.nrows(),
            k_cols: stiffness.ncols(),
        });
    }
    if u0.len() != n || v0.len() != n {
        return Err(SurvivabilityError::ShapeMismatch {
            m_rows: n,
            m_cols: n,
            c_rows: u0.len(),
            c_cols: v0.len(),
            k_rows: n,
            k_cols: n,
        });
    }

    // --- finiteness of every supplied entry (loud, not silent corruption) ---
    ensure_all_finite(mass, "mass matrix M")?;
    ensure_all_finite(damping, "damping matrix C")?;
    ensure_all_finite(stiffness, "stiffness matrix K")?;
    ensure_vec_finite(u0, "initial displacement u0")?;
    ensure_vec_finite(v0, "initial velocity v0")?;

    // --- time-step / Newmark validity ---
    if !(controls.dt.is_finite() && controls.dt > 0.0) {
        return Err(SurvivabilityError::InvalidParameter(format!(
            "time step dt must be finite and > 0, got {}",
            controls.dt
        )));
    }
    let beta = controls.newmark.beta;
    let gamma = controls.newmark.gamma;
    if !(beta.is_finite() && beta > 0.0) {
        return Err(SurvivabilityError::InvalidParameter(format!(
            "Newmark beta must be finite and > 0 (beta = 0 is explicit central difference, \
             which this implicit integrator does not handle), got {beta}"
        )));
    }
    if !gamma.is_finite() {
        return Err(SurvivabilityError::InvalidParameter(format!(
            "Newmark gamma must be finite, got {gamma}"
        )));
    }

    // --- M positive-definite (Cholesky succeeds) ---
    // A zero/negative effective mass makes the dynamics non-physical and the
    // initial-acceleration solve singular. nalgebra's Cholesky returns None
    // for a non-SPD matrix.
    let m_chol = mass
        .clone()
        .cholesky()
        .ok_or_else(|| SurvivabilityError::NonPositiveMass("Cholesky of M failed".to_string()))?;

    // --- K positive-definite (restrained, stable structure) ---
    if stiffness.clone().cholesky().is_none() {
        return Err(SurvivabilityError::NonPositiveStiffness);
    }

    // --- initial acceleration: M·a0 = F(0) − C·v0 − K·u0 ---
    let f0 = eval_load(&mut load, 0.0, n)?;
    let rhs0 = &f0 - damping * v0 - stiffness * u0;
    let a0_vec = m_chol.solve(&rhs0);

    // --- Newmark integration constants ---
    let dt = controls.dt;
    let c0 = 1.0 / (beta * dt * dt);
    let c1 = gamma / (beta * dt);
    let c2 = 1.0 / (beta * dt);
    let c3 = 1.0 / (2.0 * beta) - 1.0;
    let c4 = gamma / beta - 1.0;
    let c5 = dt * (gamma / (2.0 * beta) - 1.0);
    let c6 = dt * (1.0 - gamma);
    let c7 = dt * gamma;

    // --- constant effective stiffness K_eff = K + c0·M + c1·C, factorised once ---
    let k_eff = stiffness + mass * c0 + damping * c1;
    let k_eff_lu = k_eff.lu();

    // Pre-size the history (initial state + one per step).
    let cap = controls.n_steps + 1;
    let mut time = Vec::with_capacity(cap);
    let mut disp = Vec::with_capacity(cap);
    let mut vel = Vec::with_capacity(cap);
    let mut acc = Vec::with_capacity(cap);

    let mut u = u0.clone();
    let mut v = v0.clone();
    let mut a = a0_vec;

    time.push(0.0);
    disp.push(u.clone());
    vel.push(v.clone());
    acc.push(a.clone());

    let mut t = 0.0;
    for _ in 0..controls.n_steps {
        t += dt;
        let f_next = eval_load(&mut load, t, n)?;

        // f_eff = F(t_{n+1}) + M·(c0·u + c2·v + c3·a) + C·(c1·u + c4·v + c5·a)
        let m_term = mass * (&u * c0 + &v * c2 + &a * c3);
        let c_term = damping * (&u * c1 + &v * c4 + &a * c5);
        let f_eff = &f_next + m_term + c_term;

        let u_next = k_eff_lu
            .solve(&f_eff)
            .ok_or(SurvivabilityError::SolveFailed)?;

        // a_{n+1} = c0·(u_{n+1} − u_n) − c2·v_n − c3·a_n
        let a_next = (&u_next - &u) * c0 - &v * c2 - &a * c3;
        // v_{n+1} = v_n + c6·a_n + c7·a_{n+1}
        let v_next = &v + &a * c6 + &a_next * c7;

        u = u_next;
        v = v_next;
        a = a_next;

        time.push(t);
        disp.push(u.clone());
        vel.push(v.clone());
        acc.push(a.clone());
    }

    Ok(TransientResponse {
        time,
        displacement: disp,
        velocity: vel,
        acceleration: acc,
    })
}

/// Convenience SDOF survivability screen: build the `1×1` `M`/`C`/`K` from the
/// scalar effective mass, damping, and stiffness, drive it with a
/// [`FriedlanderPulse`] (the pressure is taken as the scalar force on the one
/// DOF — the caller has already lumped the loaded area into `peak_overpressure`
/// if needed), march it, and return the response. Starts from rest.
///
/// This is the textbook single-degree-of-freedom blast analysis: a wall/panel
/// reduced to one effective mass and stiffness, hit by the Friedlander pulse.
///
/// # Errors
///
/// [`SurvivabilityError::NonPositiveMass`] / [`SurvivabilityError::NonPositiveStiffness`]
/// for a non-positive `mass` / `stiffness`; otherwise as
/// [`solve_transient_response`].
pub fn solve_sdof_blast(
    mass: f64,
    damping: f64,
    stiffness: f64,
    pulse: &FriedlanderPulse,
    controls: &TransientControls,
) -> Result<TransientResponse, SurvivabilityError> {
    let m = DMatrix::from_element(1, 1, mass);
    let c = DMatrix::from_element(1, 1, damping);
    let k = DMatrix::from_element(1, 1, stiffness);
    let u0 = DVector::zeros(1);
    let v0 = DVector::zeros(1);
    let pulse = *pulse;
    solve_transient_response(
        &m,
        &c,
        &k,
        &u0,
        &v0,
        |t| DVector::from_element(1, pulse.pressure_at(t)),
        controls,
    )
}

/// The undamped natural angular frequency `ω = √(k/m)` of an SDOF system —
/// a small helper for sizing the time step and the analytic period.
///
/// Returns `None` if `mass` is not finite-and-positive (the frequency is then
/// undefined).
pub fn sdof_natural_frequency(mass: f64, stiffness: f64) -> Option<f64> {
    // Mass must be finite-and-positive; stiffness must be finite-and-non-negative
    // (a negative stiffness is an unstable mechanism, not an oscillator).
    if !(mass.is_finite() && mass > 0.0) {
        return None;
    }
    if !(stiffness.is_finite() && stiffness >= 0.0) {
        return None;
    }
    Some((stiffness / mass).sqrt())
}

/// The undamped natural period `T = 2π/ω = 2π·√(m/k)` of an SDOF system.
/// `None` if `stiffness` is zero (infinite period) or the frequency is
/// undefined (see [`sdof_natural_frequency`]).
pub fn sdof_natural_period(mass: f64, stiffness: f64) -> Option<f64> {
    match sdof_natural_frequency(mass, stiffness) {
        Some(w) if w > 0.0 => Some(2.0 * PI / w),
        _ => None,
    }
}

// --- internal helpers ---------------------------------------------------------

/// Evaluate the load closure at time `t`, validating the returned length and
/// finiteness. A bad load is the most likely caller error, so it fails loud.
fn eval_load<F>(load: &mut F, t: f64, n: usize) -> Result<DVector<f64>, SurvivabilityError>
where
    F: FnMut(f64) -> DVector<f64>,
{
    let f = load(t);
    if f.len() != n {
        return Err(SurvivabilityError::ShapeMismatch {
            m_rows: n,
            m_cols: n,
            c_rows: f.len(),
            c_cols: 1,
            k_rows: n,
            k_cols: n,
        });
    }
    for (i, &val) in f.iter().enumerate() {
        if !val.is_finite() {
            return Err(SurvivabilityError::NonFinite {
                what: format!("evaluated load F(t={t}) component {i}"),
            });
        }
    }
    Ok(f)
}

/// Reject a matrix carrying any non-finite entry.
fn ensure_all_finite(m: &DMatrix<f64>, what: &str) -> Result<(), SurvivabilityError> {
    if m.iter().any(|v| !v.is_finite()) {
        return Err(SurvivabilityError::NonFinite {
            what: what.to_string(),
        });
    }
    Ok(())
}

/// Reject a vector carrying any non-finite entry.
fn ensure_vec_finite(v: &DVector<f64>, what: &str) -> Result<(), SurvivabilityError> {
    if v.iter().any(|x| !x.is_finite()) {
        return Err(SurvivabilityError::NonFinite {
            what: what.to_string(),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Default Newmark is the unconditionally-stable average-acceleration
    /// (β=1/4, γ=1/2) scheme.
    #[test]
    fn newmark_defaults_to_average_acceleration() {
        let p = NewmarkBeta::default();
        assert!((p.beta - 0.25).abs() < 1e-12);
        assert!((p.gamma - 0.5).abs() < 1e-12);
        assert_eq!(p, NewmarkBeta::average_acceleration());
    }

    /// **Test 1 — SDOF undamped free vibration vs analytic `x(t)=x0·cos(ωt)`.**
    /// Release an undamped oscillator from `x0` at rest with no external load;
    /// the Newmark history must track the closed-form cosine over several
    /// periods.
    #[test]
    fn sdof_free_vibration_matches_cosine() {
        let m = 2.0_f64;
        let k = 800.0_f64;
        let omega = (k / m).sqrt(); // 20 rad/s
        let period = 2.0 * PI / omega;
        let x0 = 0.01;

        let mass = DMatrix::from_element(1, 1, m);
        let damping = DMatrix::zeros(1, 1);
        let stiffness = DMatrix::from_element(1, 1, k);
        let u0 = DVector::from_element(1, x0);
        let v0 = DVector::zeros(1);

        // ~200 steps/period over 5 periods — accuracy choice, not stability.
        let dt = period / 200.0;
        let n_steps = 1000;
        let resp = solve_transient_response(
            &mass,
            &damping,
            &stiffness,
            &u0,
            &v0,
            |_t| DVector::zeros(1),
            &TransientControls {
                dt,
                n_steps,
                newmark: NewmarkBeta::average_acceleration(),
            },
        )
        .unwrap();

        assert_eq!(resp.time.len(), n_steps + 1);
        // Compare against the analytic cosine at every recorded instant.
        let mut max_err = 0.0_f64;
        for (&t, u) in resp.time.iter().zip(&resp.displacement) {
            let analytic = x0 * (omega * t).cos();
            max_err = max_err.max((u[0] - analytic).abs());
        }
        // Average-acceleration Newmark at 200 steps/period tracks the cosine
        // tightly (period error ~1e-4 relative); tolerance on absolute
        // displacement is a small fraction of the amplitude.
        assert!(
            max_err < 0.02 * x0,
            "free-vibration history should match x0·cos(ωt); max abs error {max_err} vs x0 {x0}"
        );
        // Energy-conserving: the peak stays at the release amplitude.
        let peak = resp.peak_abs_displacement(0);
        assert!(
            (peak - x0).abs() < 0.02 * x0,
            "undamped peak {peak} should stay at the release amplitude {x0}"
        );
    }

    /// **Test 2 — SDOF undamped step load → DAF = 2.** A suddenly-applied
    /// constant load on an undamped oscillator overshoots to exactly twice the
    /// static deflection (`u_static = F/k`); the dynamic amplification factor
    /// is 2.
    #[test]
    fn sdof_step_load_dynamic_amplification_is_two() {
        let m = 1.0_f64;
        let k = 1000.0_f64;
        let omega = (k / m).sqrt();
        let period = 2.0 * PI / omega;
        let f = 50.0;
        let u_static = f / k;

        let mass = DMatrix::from_element(1, 1, m);
        let damping = DMatrix::zeros(1, 1);
        let stiffness = DMatrix::from_element(1, 1, k);
        let u0 = DVector::zeros(1);
        let v0 = DVector::zeros(1);

        let dt = period / 400.0;
        let n_steps = 1000; // ≥ 2 periods, enough to capture the first peak
        let resp = solve_transient_response(
            &mass,
            &damping,
            &stiffness,
            &u0,
            &v0,
            |_t| DVector::from_element(1, f), // step load applied at t=0, held
            &TransientControls {
                dt,
                n_steps,
                newmark: NewmarkBeta::average_acceleration(),
            },
        )
        .unwrap();

        let daf = resp
            .dynamic_amplification_factor(0, u_static)
            .expect("nonzero static deflection");
        assert!(
            (daf - 2.0).abs() < 0.01,
            "undamped step-load DAF should be 2.0, got {daf}"
        );
    }

    /// **Test 3 — SDOF damped response decays as `e^{−ζω t}`.** With damping
    /// ratio `ζ`, the successive peaks of the free vibration fall on the
    /// analytic envelope `x0·e^{−ζω t}`. We compare an early peak and a late
    /// peak to the envelope evaluated at their times.
    #[test]
    fn sdof_damped_envelope_matches_exponential() {
        let m = 1.0_f64;
        let k = 400.0_f64;
        let omega = (k / m).sqrt(); // 20 rad/s
        let zeta = 0.05;
        // c = 2ζ√(km) = 2ζωm
        let c = 2.0 * zeta * omega * m;
        let period = 2.0 * PI / omega;
        let x0 = 0.02;

        let mass = DMatrix::from_element(1, 1, m);
        let damping = DMatrix::from_element(1, 1, c);
        let stiffness = DMatrix::from_element(1, 1, k);
        let u0 = DVector::from_element(1, x0);
        let v0 = DVector::zeros(1);

        let dt = period / 400.0;
        let n_steps = 400 * 8; // 8 periods
        let resp = solve_transient_response(
            &mass,
            &damping,
            &stiffness,
            &u0,
            &v0,
            |_t| DVector::zeros(1),
            &TransientControls {
                dt,
                n_steps,
                newmark: NewmarkBeta::average_acceleration(),
            },
        )
        .unwrap();

        // Collect local maxima (upward peaks) of the displacement history.
        let h = resp.dof_history(0);
        let mut peaks: Vec<(f64, f64)> = Vec::new();
        for w in h.windows(3) {
            if w[1].1 > w[0].1 && w[1].1 >= w[2].1 && w[1].1 > 0.0 {
                peaks.push((w[1].0, w[1].1));
            }
        }
        assert!(
            peaks.len() >= 4,
            "damped oscillation should show several peaks, got {}",
            peaks.len()
        );

        // Each peak must sit on the envelope x0·e^{−ζω t} within tolerance.
        for &(t, amp) in &peaks {
            let envelope = x0 * (-zeta * omega * t).exp();
            let rel = (amp - envelope).abs() / envelope;
            assert!(
                rel < 0.05,
                "peak at t={t} amp={amp} should match envelope {envelope} (rel {rel})"
            );
        }

        // Sanity: a late peak is markedly smaller than an early one (it decays).
        let first = peaks.first().unwrap().1;
        let last = peaks.last().unwrap().1;
        assert!(
            last < 0.7 * first,
            "envelope must decay: first {first}, last {last}"
        );

        // And the log-decrement between the first two peaks recovers ζ.
        // δ = ln(x_i / x_{i+1}) ≈ 2πζ/√(1−ζ²) ≈ 2πζ for small ζ.
        let delta = (peaks[0].1 / peaks[1].1).ln();
        let zeta_est = delta / (2.0 * PI);
        assert!(
            (zeta_est - zeta).abs() < 0.02,
            "log-decrement should recover ζ≈{zeta}, got {zeta_est}"
        );
    }

    /// **Test 4 — Friedlander pulse impulse equals the analytic `∫F dt`.**
    /// Numerically integrate `pressure_at` over the positive phase with a fine
    /// trapezoidal rule and compare to the closed-form `positive_impulse`.
    /// Also checks the peak and the zero crossing.
    #[test]
    fn friedlander_impulse_matches_analytic() {
        let p0 = 5.0e4; // 50 kPa peak overpressure
        let td = 8.0e-3; // 8 ms positive phase
        let b = 1.3;
        let pulse = FriedlanderPulse::new(p0, td, b).unwrap();

        // Peak is at t=0 (arrival) and equals p0.
        assert!((pulse.pressure_at(0.0) - p0).abs() < 1e-6 * p0);
        // Zero crossing at exactly t = td.
        assert!(pulse.pressure_at(td).abs() < 1e-9 * p0);
        // Zero before arrival.
        assert_eq!(pulse.pressure_at(-1.0e-3), 0.0);

        // Fine trapezoidal integral over the positive phase [0, td].
        let n = 200_000;
        let h = td / n as f64;
        let mut integral = 0.0;
        for i in 0..n {
            let t0 = i as f64 * h;
            let t1 = (i + 1) as f64 * h;
            integral += 0.5 * (pulse.pressure_at(t0) + pulse.pressure_at(t1)) * h;
        }
        let analytic = pulse.positive_impulse();
        let rel = (integral - analytic).abs() / analytic.abs();
        assert!(
            rel < 1e-4,
            "numeric impulse {integral} should match analytic {analytic} (rel {rel})"
        );

        // Cross-check the closed form against the small-b triangular limit:
        // as b→0, i → p0·td/2.
        let tri = FriedlanderPulse::new(p0, td, 1.0e-6).unwrap();
        let tri_i = tri.positive_impulse();
        assert!(
            (tri_i - p0 * td / 2.0).abs() / (p0 * td / 2.0) < 1e-3,
            "small-b impulse {tri_i} should approach triangular p0·td/2 {}",
            p0 * td / 2.0
        );
    }

    /// The SDOF blast convenience driver runs end-to-end and produces a bounded,
    /// physically-sensible peak response to a Friedlander pulse.
    #[test]
    fn sdof_blast_screen_runs_and_is_bounded() {
        let m = 500.0_f64; // effective mass
        let k = 2.0e7_f64; // effective stiffness
        let omega = (k / m).sqrt();
        let period = 2.0 * PI / omega;
        let zeta = 0.02;
        let c = 2.0 * zeta * omega * m;

        let pulse = FriedlanderPulse::new(3.0e4, period * 0.5, 1.0).unwrap();
        let resp = solve_sdof_blast(
            m,
            c,
            k,
            &pulse,
            &TransientControls {
                dt: period / 200.0,
                n_steps: 1200,
                newmark: NewmarkBeta::average_acceleration(),
            },
        )
        .unwrap();

        let peak = resp.peak_abs_displacement(0);
        // A finite, positive peak; the structure responds and then (lightly
        // damped) rings down.
        assert!(
            peak.is_finite() && peak > 0.0,
            "blast peak should be finite-positive, got {peak}"
        );

        // DAF against the peak-load static deflection is a sensible O(1) number
        // for a pulse comparable to the period (impulsive-to-quasi-static
        // regime), and certainly below the step-load ceiling of 2.
        let u_static = pulse.peak_overpressure / k;
        let daf = resp.dynamic_amplification_factor(0, u_static).unwrap();
        assert!(
            daf > 0.0 && daf < 2.0,
            "pulse DAF should be in (0,2), got {daf}"
        );

        // Acceleration readout is finite.
        assert!(resp.peak_abs_acceleration(0).is_finite());
    }

    /// A 2-DOF (MDOF) lumped-mass system integrates without error — the
    /// analyzer is not SDOF-only.
    #[test]
    fn mdof_two_dof_runs() {
        // Fixed–m–m chain: K = [[2k,−k],[−k,2k]], M = diag(m,m).
        let m = 1.0;
        let k = 1000.0;
        let mass = DMatrix::from_diagonal(&DVector::from_vec(vec![m, m]));
        let stiffness = DMatrix::from_row_slice(2, 2, &[2.0 * k, -k, -k, 2.0 * k]);
        let damping = rayleigh_damping_matrix(&mass, &stiffness, 0.5, 1.0e-4).unwrap();
        let u0 = DVector::from_vec(vec![0.005, 0.0]);
        let v0 = DVector::zeros(2);

        let resp = solve_transient_response(
            &mass,
            &damping,
            &stiffness,
            &u0,
            &v0,
            |_t| DVector::zeros(2),
            &TransientControls {
                dt: 1.0e-4,
                n_steps: 500,
                newmark: NewmarkBeta::average_acceleration(),
            },
        )
        .unwrap();
        assert_eq!(resp.displacement.len(), 501);
        // Damped → DOF 0 amplitude decays below its start.
        assert!(resp.peak_abs_displacement(0) <= 0.005 + 1e-9);
        // Both DOFs have a finite history.
        assert!(resp.peak_abs_displacement(1).is_finite());
    }

    /// `rayleigh_damping_matrix` returns `C = αM + βK` and rejects mismatched
    /// shapes.
    #[test]
    fn rayleigh_matrix_builds_and_checks_shape() {
        let mass = DMatrix::identity(2, 2);
        let stiffness = DMatrix::from_row_slice(2, 2, &[4.0, -1.0, -1.0, 4.0]);
        let c = rayleigh_damping_matrix(&mass, &stiffness, 2.0, 3.0).unwrap();
        // C = 2·I + 3·K
        assert!((c[(0, 0)] - (2.0 + 3.0 * 4.0)).abs() < 1e-12);
        assert!((c[(0, 1)] - (-3.0)).abs() < 1e-12);

        let bad = DMatrix::identity(3, 3);
        assert!(matches!(
            rayleigh_damping_matrix(&mass, &bad, 1.0, 1.0),
            Err(SurvivabilityError::ShapeMismatch { .. })
        ));
    }

    /// SDOF natural-frequency / period helpers.
    #[test]
    fn sdof_frequency_and_period_helpers() {
        let w = sdof_natural_frequency(2.0, 800.0).unwrap();
        assert!((w - 20.0).abs() < 1e-9);
        let t = sdof_natural_period(2.0, 800.0).unwrap();
        assert!((t - 2.0 * PI / 20.0).abs() < 1e-9);
        // Degenerate: zero stiffness → no finite period; bad mass → None.
        assert!(sdof_natural_period(2.0, 0.0).is_none());
        assert!(sdof_natural_frequency(0.0, 800.0).is_none());
        assert!(sdof_natural_frequency(-1.0, 800.0).is_none());
    }

    // ---- Test 5: degenerate inputs fail loud (no panic) ----

    /// Zero / negative effective mass → `NonPositiveMass`.
    #[test]
    fn rejects_non_positive_mass() {
        let stiffness = DMatrix::from_element(1, 1, 100.0);
        let damping = DMatrix::zeros(1, 1);
        let u0 = DVector::zeros(1);
        let v0 = DVector::zeros(1);

        for bad_mass in [0.0, -5.0] {
            let mass = DMatrix::from_element(1, 1, bad_mass);
            let err = solve_transient_response(
                &mass,
                &damping,
                &stiffness,
                &u0,
                &v0,
                |_t| DVector::zeros(1),
                &TransientControls::default(),
            )
            .unwrap_err();
            assert!(
                matches!(err, SurvivabilityError::NonPositiveMass(_)),
                "mass {bad_mass} should be rejected, got {err:?}"
            );
        }
    }

    /// Non-positive-definite (zero / negative) stiffness → `NonPositiveStiffness`.
    #[test]
    fn rejects_non_positive_stiffness() {
        let mass = DMatrix::from_element(1, 1, 1.0);
        let damping = DMatrix::zeros(1, 1);
        let u0 = DVector::zeros(1);
        let v0 = DVector::zeros(1);

        for bad_k in [0.0, -100.0] {
            let stiffness = DMatrix::from_element(1, 1, bad_k);
            let err = solve_transient_response(
                &mass,
                &damping,
                &stiffness,
                &u0,
                &v0,
                |_t| DVector::zeros(1),
                &TransientControls::default(),
            )
            .unwrap_err();
            assert!(
                matches!(err, SurvivabilityError::NonPositiveStiffness),
                "stiffness {bad_k} should be rejected, got {err:?}"
            );
        }
    }

    /// Bad time step (zero, negative, NaN) and bad Newmark β → `InvalidParameter`.
    #[test]
    fn rejects_bad_dt_and_beta() {
        let mass = DMatrix::from_element(1, 1, 1.0);
        let damping = DMatrix::zeros(1, 1);
        let stiffness = DMatrix::from_element(1, 1, 100.0);
        let u0 = DVector::zeros(1);
        let v0 = DVector::zeros(1);

        for bad_dt in [0.0, -1.0e-3, f64::NAN, f64::INFINITY] {
            let err = solve_transient_response(
                &mass,
                &damping,
                &stiffness,
                &u0,
                &v0,
                |_t| DVector::zeros(1),
                &TransientControls {
                    dt: bad_dt,
                    n_steps: 10,
                    newmark: NewmarkBeta::average_acceleration(),
                },
            )
            .unwrap_err();
            assert!(
                matches!(err, SurvivabilityError::InvalidParameter(_)),
                "dt {bad_dt} should be rejected"
            );
        }

        // Newmark beta = 0 (explicit central difference, unsupported here).
        let err = solve_transient_response(
            &mass,
            &damping,
            &stiffness,
            &u0,
            &v0,
            |_t| DVector::zeros(1),
            &TransientControls {
                dt: 1.0e-3,
                n_steps: 10,
                newmark: NewmarkBeta {
                    beta: 0.0,
                    gamma: 0.5,
                },
            },
        )
        .unwrap_err();
        assert!(matches!(err, SurvivabilityError::InvalidParameter(_)));
    }

    /// Empty system and shape mismatch → the corresponding errors.
    #[test]
    fn rejects_empty_and_mismatched_shapes() {
        let empty = DMatrix::<f64>::zeros(0, 0);
        let v_empty = DVector::<f64>::zeros(0);
        let err = solve_transient_response(
            &empty,
            &empty,
            &empty,
            &v_empty,
            &v_empty,
            |_t| DVector::zeros(0),
            &TransientControls::default(),
        )
        .unwrap_err();
        assert_eq!(err, SurvivabilityError::EmptySystem);

        // K is 2×2 but M is 1×1.
        let m = DMatrix::from_element(1, 1, 1.0);
        let c = DMatrix::zeros(1, 1);
        let k = DMatrix::identity(2, 2);
        let u0 = DVector::zeros(1);
        let v0 = DVector::zeros(1);
        let err = solve_transient_response(
            &m,
            &c,
            &k,
            &u0,
            &v0,
            |_t| DVector::zeros(1),
            &TransientControls::default(),
        )
        .unwrap_err();
        assert!(matches!(err, SurvivabilityError::ShapeMismatch { .. }));
    }

    /// A non-finite matrix entry and a non-finite evaluated load both fail loud.
    #[test]
    fn rejects_non_finite_inputs() {
        let damping = DMatrix::zeros(1, 1);
        let stiffness = DMatrix::from_element(1, 1, 100.0);
        let u0 = DVector::zeros(1);
        let v0 = DVector::zeros(1);

        // NaN in the mass matrix.
        let mass_nan = DMatrix::from_element(1, 1, f64::NAN);
        let err = solve_transient_response(
            &mass_nan,
            &damping,
            &stiffness,
            &u0,
            &v0,
            |_t| DVector::zeros(1),
            &TransientControls::default(),
        )
        .unwrap_err();
        assert!(matches!(err, SurvivabilityError::NonFinite { .. }));

        // A load closure that returns a NaN force.
        let mass = DMatrix::from_element(1, 1, 1.0);
        let err = solve_transient_response(
            &mass,
            &damping,
            &stiffness,
            &u0,
            &v0,
            |_t| DVector::from_element(1, f64::NAN),
            &TransientControls {
                dt: 1.0e-3,
                n_steps: 5,
                newmark: NewmarkBeta::average_acceleration(),
            },
        )
        .unwrap_err();
        assert!(matches!(err, SurvivabilityError::NonFinite { .. }));

        // A load closure returning the wrong length is a shape mismatch.
        let err = solve_transient_response(
            &mass,
            &damping,
            &stiffness,
            &u0,
            &v0,
            |_t| DVector::zeros(3),
            &TransientControls {
                dt: 1.0e-3,
                n_steps: 5,
                newmark: NewmarkBeta::average_acceleration(),
            },
        )
        .unwrap_err();
        assert!(matches!(err, SurvivabilityError::ShapeMismatch { .. }));
    }

    /// Friedlander constructor rejects non-physical parameters.
    #[test]
    fn friedlander_rejects_bad_params() {
        assert!(matches!(
            FriedlanderPulse::new(f64::NAN, 1.0e-3, 1.0),
            Err(SurvivabilityError::InvalidParameter(_))
        ));
        assert!(matches!(
            FriedlanderPulse::new(1.0e4, 0.0, 1.0),
            Err(SurvivabilityError::InvalidParameter(_))
        ));
        assert!(matches!(
            FriedlanderPulse::new(1.0e4, -1.0e-3, 1.0),
            Err(SurvivabilityError::InvalidParameter(_))
        ));
        assert!(matches!(
            FriedlanderPulse::new(1.0e4, 1.0e-3, 0.0),
            Err(SurvivabilityError::InvalidParameter(_))
        ));
        assert!(matches!(
            FriedlanderPulse::new(1.0e4, 1.0e-3, -1.0),
            Err(SurvivabilityError::InvalidParameter(_))
        ));
    }

    /// The n_steps DoS cap is enforced.
    #[test]
    fn rejects_n_steps_past_cap() {
        let mass = DMatrix::from_element(1, 1, 1.0);
        let damping = DMatrix::zeros(1, 1);
        let stiffness = DMatrix::from_element(1, 1, 100.0);
        let u0 = DVector::zeros(1);
        let v0 = DVector::zeros(1);
        let err = solve_transient_response(
            &mass,
            &damping,
            &stiffness,
            &u0,
            &v0,
            |_t| DVector::zeros(1),
            &TransientControls {
                dt: 1.0e-3,
                n_steps: MAX_TIME_STEPS + 1,
                newmark: NewmarkBeta::average_acceleration(),
            },
        )
        .unwrap_err();
        assert!(matches!(err, SurvivabilityError::InvalidParameter(_)));
    }
}
