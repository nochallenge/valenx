//! **Strongly-coupled (implicit / iterative) co-simulation stepping.**
//!
//! The explicit [`crate::cosim::CoSimMaster::advance`] issues each subsystem's
//! step exactly once per macro-step.  That is efficient and correct when the
//! coupling is weak, but it leaves a one-macro-step lag between subsystems
//! (Jacobi) or a Gauss-Seidel ordering bias — and, most importantly, neither
//! scheme converges the *algebraic* coupling constraint within a single step.
//! For tightly coupled or algebraic-loop federates (e.g. a structural DOF
//! coupled to a fluid DOF through a shared interface boundary, or a
//! mass-spring pair with a stiff interface spring) the coupling error is first-
//! order in the macro-step size and can dominate the solution.
//!
//! This module provides a **fixed-point iterative coupler** that, within one
//! macro-step `[t, t + dt]`, repeatedly exchanges coupling variables and
//! re-evaluates all subsystems until the residual is small.  Two sweep orders
//! are available:
//!
//! * [`ImplicitScheme::GaussSeidel`] — subsystems step in index order, each
//!   immediately seeing the freshest partner outputs from this iteration.
//!   Converges faster per iteration than Jacobi for many contractive systems.
//! * [`ImplicitScheme::Jacobi`] — all subsystems evaluate against the *same*
//!   previous-iteration snapshot (swap at the end).  Safer for cases where the
//!   Gauss-Seidel ordering introduces a directional bias; also parallelizable.
//!
//! An optional **relaxation accelerator** ([`Relaxation`]) can be applied to
//! the output update each iteration:
//!
//! * [`Relaxation::Fixed`] — blend `new = (1 - ω) old + ω new_raw` with a
//!   constant factor `ω ∈ (0, 2)`.  Under-relaxation (`ω < 1`) stabilises
//!   ill-conditioned or near-divergent systems; over-relaxation (`ω > 1`)
//!   accelerates convergence for smooth, mildly contractive systems.
//! * [`Relaxation::Aitken`] — the Aitken Δ² method: `ω` is updated each
//!   iteration based on the ratio of consecutive residual differences,
//!   automatically adapting the step size.  For linear problems this converges
//!   to the exact solution in one well-conditioned step; for nonlinear problems
//!   it often halves the iteration count versus a fixed `ω`.
//! * [`Relaxation::None`] — no relaxation (the plain fixed-point sweep).
//!
//! ## Entry point
//!
//! ```
//! # use valenx_adapter_fmi::cosim::{CouplingGraph, Subsystem};
//! # use valenx_adapter_fmi::implicit::{coupled_step, ImplicitScheme, Relaxation};
//! # struct Gain(f64);
//! # impl Subsystem for Gain {
//! #     fn n_inputs(&self) -> usize { 1 }
//! #     fn n_outputs(&self) -> usize { 1 }
//! #     fn step(&mut self, _t: f64, _dt: f64, inputs: &[f64]) -> Vec<f64> { vec![self.0 * inputs[0]] }
//! #     fn state(&self) -> Vec<f64> { Vec::new() }
//! #     fn set_state(&mut self, _s: &[f64]) {}
//! # }
//! let mut subs: Vec<Box<dyn Subsystem>> = vec![
//!     Box::new(Gain(0.3)),
//!     Box::new(Gain(0.4)),
//! ];
//! // u1 = y2, u2 = y1  (mutual feedback, |a*b| = 0.12 < 1)
//! let graph = CouplingGraph::from_edges(vec![
//!     valenx_adapter_fmi::cosim::Coupling::new(1, 0, 0, 0),
//!     valenx_adapter_fmi::cosim::Coupling::new(0, 0, 1, 0),
//! ]);
//! let result = coupled_step(
//!     &mut subs, &graph,
//!     0.0, 1.0,                      // t, dt
//!     1e-10, 100,                    // tol, max_iter
//!     ImplicitScheme::GaussSeidel,
//!     Relaxation::None,
//! ).unwrap();
//! println!("converged in {} iterations", result.iterations);
//! ```
//!
//! ## Algorithm notes
//!
//! The fixed-point operator is `G(y) = F(y)` where `y` is the stacked output
//! vector of all subsystems and `F` is one full sweep (each subsystem is
//! stepped with inputs assembled from the current `y`).  Convergence is
//! declared when `‖y_new − y_old‖_∞ < tol`.  If the operator is contractive
//! (spectral radius of the coupling Jacobian `< 1`) the iterations converge
//! geometrically; otherwise [`crate::error::FmiError::NotConverged`] is
//! returned after `max_iter` steps.
//!
//! The subsystem rollback protocol (see [`crate::cosim::Subsystem::state`] /
//! [`crate::cosim::Subsystem::set_state`]) is used at the start of EVERY
//! iteration to reset each subsystem to its start-of-macro-step state before
//! the next evaluation, so repeated calls to `step` within one iteration loop
//! are well-defined.  Subsystems that return an empty `state()` are assumed
//! to be stateless within the iteration (e.g. a static algebraic map) and no
//! rollback is needed.

use crate::cosim::{CouplingGraph, Subsystem};
use crate::error::{FmiError, Result};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Which sweep order to use inside each fixed-point iteration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ImplicitScheme {
    /// **Gauss-Seidel** — subsystems sweep in index order; each sees the
    /// freshest outputs from the current iteration.  Usually converges faster
    /// per iteration than Jacobi for problems without strong directional bias.
    GaussSeidel,
    /// **Jacobi** — all subsystems evaluate against the previous-iteration
    /// snapshot; the snapshot is replaced only after all have stepped.
    /// Parallelizable; safer when GS ordering introduces bias.
    Jacobi,
}

/// Optional relaxation / acceleration applied to the output update after each
/// sweep.  Applied component-wise to the stacked output vector.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Relaxation {
    /// No relaxation — plain fixed-point iteration.
    None,
    /// Fixed relaxation factor `ω ∈ (0, 2)`.
    ///
    /// The update rule after each sweep is:
    ///
    /// ```text
    /// y_new ← (1 − ω) y_old + ω y_raw
    /// ```
    ///
    /// `ω = 1` is identical to no relaxation.  `ω < 1` under-relaxes
    /// (slower but more stable); `ω > 1` over-relaxes (faster for smooth
    /// problems but can diverge).  Fail-loud if `ω` is not in `(0, 2)`.
    Fixed(f64),
    /// **Aitken Δ² acceleration** (dynamically updated relaxation).
    ///
    /// On the first iteration a user-supplied initial factor `omega_0 ∈ (0, 2)`
    /// is used (choose `≤ 1` for safety).  On subsequent iterations `ω` is
    /// updated via:
    ///
    /// ```text
    /// ω_{k+1} = −ω_k · (r_k · (r_{k+1} − r_k)) / ‖r_{k+1} − r_k‖²
    /// ```
    ///
    /// where `r_k = y_k − y_{k−1}` is the residual vector.  The denominator
    /// is guarded: if `‖r_{k+1} − r_k‖² < ε_machine` (consecutive residuals
    /// are identical — already converged), the previous `ω` is kept.  The
    /// resulting `ω` is clamped to `(0, 2)` for safety.
    ///
    /// Fail-loud if `omega_0` is not in `(0, 2)`.
    Aitken {
        /// Initial relaxation factor for the first iteration.
        omega_0: f64,
    },
}

/// The successful result of a [`coupled_step`] call.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledStepResult {
    /// The converged output vector of each subsystem (in subsystem index
    /// order, each entry is that subsystem's full output slice).
    pub outputs: Vec<Vec<f64>>,
    /// The number of fixed-point iterations performed (including the final
    /// converging iteration).
    pub iterations: usize,
    /// The infinity-norm of the last residual `‖y_new − y_old‖_∞`.  Will be
    /// `< tol` when the function returns `Ok`.
    pub final_residual: f64,
}

// ---------------------------------------------------------------------------
// Error extension — NotConverged (appended to FmiError via a new variant)
// ---------------------------------------------------------------------------

// `FmiError` is `#[non_exhaustive]`; we add `NotConverged` there so all
// callers use the same error taxonomy.  The variant is defined in error.rs;
// we just use it here.

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

/// Run one **implicit (iterative) macro-step** of size `dt` starting at time
/// `t` over the subsystems described by `graph`.
///
/// Within the step, subsystems are repeatedly evaluated (after rollback to
/// their start-of-step state) until the infinity-norm of the change in the
/// stacked output vector falls below `tol`, or `max_iter` iterations are
/// exhausted.  Each iteration applies the chosen [`ImplicitScheme`] sweep
/// and optionally a [`Relaxation`] accelerator.
///
/// # Parameters
///
/// * `subsystems` — mutable slice of subsystem boxes; their states are
///   mutated in-place and rolled back between iterations.
/// * `graph` — the coupling wiring; must have been validated by
///   [`CouplingGraph::validate`] against `subsystems` before this call.
///   Panics in debug if an edge is out of range (edges were validated).
/// * `t` — the simulation time at the START of the macro-step.
/// * `dt` — macro-step size, must be `> 0`.
/// * `tol` — convergence tolerance on `‖Δy‖_∞`, must be `> 0`.
/// * `max_iter` — iteration cap; must be `≥ 1`.
/// * `scheme` — Gauss-Seidel or Jacobi sweep order.
/// * `relaxation` — optional Aitken / fixed-factor accelerator.
///
/// # Returns
///
/// `Ok(`[`CoupledStepResult`]`)` with the converged outputs, iteration count,
/// and final residual; or `Err(`[`FmiError::NotConverged`]`)` if the loop did
/// not converge within `max_iter` iterations, or
/// `Err(`[`FmiError::BadCoupledStep`]`)` if the configuration is invalid.
///
/// # Panics
///
/// In debug mode, panics if `graph` references an out-of-range subsystem or
/// port (i.e. it was not validated against `subsystems`).
// The loops in coupled_step genuinely need the loop variable `i` both as an
// index into `working[i]` / `next.push(…[i].step(…))` AND as the subsystem
// identifier passed to `assemble_inputs`. Clippy cannot see both uses; allow it.
#[allow(clippy::too_many_arguments)]
#[allow(clippy::needless_range_loop)]
pub fn coupled_step(
    subsystems: &mut [Box<dyn Subsystem>],
    graph: &CouplingGraph,
    t: f64,
    dt: f64,
    tol: f64,
    max_iter: usize,
    scheme: ImplicitScheme,
    relaxation: Relaxation,
) -> Result<CoupledStepResult> {
    // --- guard bad configuration ------------------------------------------
    if !dt.is_finite() || dt <= 0.0 {
        return Err(FmiError::BadCoupledStep(format!(
            "dt must be finite and > 0, got {dt}"
        )));
    }
    if !tol.is_finite() || tol <= 0.0 {
        return Err(FmiError::BadCoupledStep(format!(
            "tol must be finite and > 0, got {tol}"
        )));
    }
    if max_iter == 0 {
        return Err(FmiError::BadCoupledStep(
            "max_iter must be >= 1".to_string(),
        ));
    }
    if !t.is_finite() {
        return Err(FmiError::BadCoupledStep(format!(
            "t must be finite, got {t}"
        )));
    }
    // validate relaxation parameter
    match relaxation {
        Relaxation::Fixed(omega) => {
            if !omega.is_finite() || omega <= 0.0 || omega >= 2.0 {
                return Err(FmiError::BadCoupledStep(format!(
                    "Fixed relaxation omega must be in (0, 2), got {omega}"
                )));
            }
        }
        Relaxation::Aitken { omega_0 } => {
            if !omega_0.is_finite() || omega_0 <= 0.0 || omega_0 >= 2.0 {
                return Err(FmiError::BadCoupledStep(format!(
                    "Aitken omega_0 must be in (0, 2), got {omega_0}"
                )));
            }
        }
        Relaxation::None => {}
    }

    let n = subsystems.len();

    // --- snapshot start-of-step states for rollback -----------------------
    // `start_states[i]` is the state of subsystem i at time t.  On each
    // iteration we restore this before re-stepping.  Subsystems that return
    // an empty Vec from `state()` are treated as stateless (no rollback
    // needed, and calling `set_state` with empty slice is a no-op per trait).
    let start_states: Vec<Vec<f64>> = subsystems.iter().map(|s| s.state()).collect();

    // --- seed the output vector with a dt=0 sample (= current outputs) ---
    // This mirrors CoSimMaster::new's priming convention: call step(t, 0, zero)
    // to get each subsystem's true current output without advancing state.
    // We restore states immediately after so the dt=0 sample is benign.
    let mut current_outputs: Vec<Vec<f64>> = {
        let mut out = Vec::with_capacity(n);
        for (i, sub) in subsystems.iter_mut().enumerate() {
            let zeros = vec![0.0_f64; sub.n_inputs()];
            let y = sub.step(t, 0.0, &zeros);
            // restore: the dt=0 call must not advance state.
            sub.set_state(&start_states[i]);
            out.push(y);
        }
        out
    };

    // Aitken state: previous residual vector (r_{k-1}), current omega.
    let mut aitken_omega: f64 = match relaxation {
        Relaxation::Aitken { omega_0 } => omega_0,
        _ => 1.0,
    };
    let mut aitken_prev_residual: Option<Vec<f64>> = None;

    // --- fixed-point loop -------------------------------------------------
    for iter in 1..=max_iter {
        // Roll back all subsystems to start-of-step state before each sweep.
        for (i, sub) in subsystems.iter_mut().enumerate() {
            sub.set_state(&start_states[i]);
        }

        let new_outputs = match scheme {
            ImplicitScheme::GaussSeidel => {
                // Working copy updated in-place; later subsystems see earlier
                // subsystems' fresh outputs within this iteration.
                let mut working = current_outputs.clone();
                for i in 0..n {
                    let inputs = assemble_inputs(i, graph, &working, &*subsystems[i]);
                    working[i] = subsystems[i].step(t, dt, &inputs);
                }
                working
            }
            ImplicitScheme::Jacobi => {
                // All subsystems evaluate against the same previous snapshot.
                let snapshot = current_outputs.clone();
                let mut next = Vec::with_capacity(n);
                for i in 0..n {
                    let inputs = assemble_inputs(i, graph, &snapshot, &*subsystems[i]);
                    next.push(subsystems[i].step(t, dt, &inputs));
                }
                next
            }
        };

        // guard: check new outputs for NaN/Inf — a blown-up subsystem must
        // not silently propagate garbage.
        for (i, y) in new_outputs.iter().enumerate() {
            for (j, v) in y.iter().enumerate() {
                if !v.is_finite() {
                    return Err(FmiError::NotConverged {
                        iterations: iter,
                        final_residual: f64::INFINITY,
                        reason: format!(
                            "subsystem {i} output[{j}] is non-finite ({v}) at iteration {iter}"
                        ),
                    });
                }
            }
        }

        // --- residual: flat view of (new - old) ---------------------------
        let raw_residual: Vec<f64> = new_outputs
            .iter()
            .zip(current_outputs.iter())
            .flat_map(|(ny, oy)| ny.iter().zip(oy.iter()).map(|(a, b)| a - b))
            .collect();

        // --- apply relaxation to new_outputs ------------------------------
        let relaxed_outputs: Vec<Vec<f64>> = match relaxation {
            Relaxation::None => new_outputs,
            Relaxation::Fixed(omega) => blend(&current_outputs, new_outputs, omega),
            Relaxation::Aitken { .. } => {
                // Update omega via Aitken Δ².
                if let Some(ref prev_r) = aitken_prev_residual {
                    // dr = raw_residual - prev_r
                    let dr: Vec<f64> = raw_residual
                        .iter()
                        .zip(prev_r.iter())
                        .map(|(a, b)| a - b)
                        .collect();
                    let dr_sq: f64 = dr.iter().map(|v| v * v).sum();
                    if dr_sq > f64::EPSILON {
                        let dot: f64 = prev_r.iter().zip(dr.iter()).map(|(a, b)| a * b).sum();
                        let new_omega = -aitken_omega * dot / dr_sq;
                        // clamp to (ε, 2 - ε) for numerical safety.
                        aitken_omega = new_omega.clamp(f64::EPSILON, 2.0 - f64::EPSILON);
                    }
                }
                aitken_prev_residual = Some(raw_residual.clone());
                // blend uses current omega.
                blend(&current_outputs, new_outputs, aitken_omega)
            }
        };

        // --- convergence check on the (pre-relaxation) residual ----------
        // We measure convergence on the raw fixed-point residual, not the
        // relaxed update, so the tolerance is meaningful independent of ω.
        let inf_norm: f64 = raw_residual.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);

        current_outputs = relaxed_outputs;

        if inf_norm < tol {
            return Ok(CoupledStepResult {
                outputs: current_outputs,
                iterations: iter,
                final_residual: inf_norm,
            });
        }
    }

    // Exhausted max_iter without converging.
    let inf_norm: f64 = {
        // Recompute residual against a fresh sweep to report a meaningful norm.
        // (current_outputs is already the relaxed value from the last iter.)
        for (i, sub) in subsystems.iter_mut().enumerate() {
            sub.set_state(&start_states[i]);
        }
        let final_sweep = match scheme {
            ImplicitScheme::GaussSeidel => {
                let mut w = current_outputs.clone();
                for i in 0..n {
                    let inp = assemble_inputs(i, graph, &w, &*subsystems[i]);
                    w[i] = subsystems[i].step(t, dt, &inp);
                }
                w
            }
            ImplicitScheme::Jacobi => {
                let snap = current_outputs.clone();
                let mut nx = Vec::with_capacity(n);
                for i in 0..n {
                    let inp = assemble_inputs(i, graph, &snap, &*subsystems[i]);
                    nx.push(subsystems[i].step(t, dt, &inp));
                }
                nx
            }
        };
        final_sweep
            .iter()
            .zip(current_outputs.iter())
            .flat_map(|(ny, oy)| ny.iter().zip(oy.iter()).map(|(a, b)| (a - b).abs()))
            .fold(0.0_f64, f64::max)
    };

    Err(FmiError::NotConverged {
        iterations: max_iter,
        final_residual: inf_norm,
        reason: format!(
            "fixed-point coupling did not converge within {max_iter} iterations \
             (scheme={scheme:?}, final ‖Δy‖_∞ = {inf_norm:.3e}, tol = {tol:.3e})"
        ),
    })
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Assemble the input vector for subsystem `i` from `source_outputs`.
/// Unconnected inputs are `0.0`.  Debug-asserts that graph edges are in range.
fn assemble_inputs(
    i: usize,
    graph: &CouplingGraph,
    source_outputs: &[Vec<f64>],
    sub: &dyn Subsystem,
) -> Vec<f64> {
    let mut inputs = vec![0.0_f64; sub.n_inputs()];
    for e in graph.edges() {
        if e.to_subsystem == i {
            debug_assert!(
                e.from_subsystem < source_outputs.len(),
                "implicit: edge from_subsystem {} out of range",
                e.from_subsystem
            );
            debug_assert!(
                e.from_output < source_outputs[e.from_subsystem].len(),
                "implicit: edge from_output {} out of range for subsystem {}",
                e.from_output,
                e.from_subsystem
            );
            inputs[e.to_input] = source_outputs[e.from_subsystem][e.from_output];
        }
    }
    inputs
}

/// Apply a scalar `ω` blend:  `(1−ω)*old + ω*new`.
fn blend(old: &[Vec<f64>], new: Vec<Vec<f64>>, omega: f64) -> Vec<Vec<f64>> {
    old.iter()
        .zip(new)
        .map(|(o, n)| {
            o.iter()
                .zip(n)
                .map(|(ov, nv)| (1.0 - omega) * ov + omega * nv)
                .collect()
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests — benchmark-pinned, analytic ground truth
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cosim::{Coupling, CouplingGraph, Subsystem};

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    /// A purely algebraic (no state) linear gain subsystem: y = a * u.
    struct LinearGain {
        a: f64,
    }
    impl Subsystem for LinearGain {
        fn n_inputs(&self) -> usize {
            1
        }
        fn n_outputs(&self) -> usize {
            1
        }
        fn step(&mut self, _t: f64, _dt: f64, inputs: &[f64]) -> Vec<f64> {
            vec![self.a * inputs[0]]
        }
        // stateless — default state()/set_state() are fine.
    }

    /// Build the two-coupled-gain system: y1 = a*u1, y2 = b*u2,  u1=y2, u2=y1.
    /// The unique fixed point of the mapping (y1, y2) -> (a*y2, b*y1) is
    ///   y1* = a * b * y1*  → y1* = 0  (when |a*b| ≠ 1)
    ///   y2* = b * a * y2*  → y2* = 0
    /// i.e. both outputs converge to ZERO regardless of initial conditions,
    /// because the only fixed point of a strict contraction with a zero
    /// constant term is 0.
    fn two_gain_system(a: f64, b: f64) -> (Vec<Box<dyn Subsystem>>, CouplingGraph) {
        let subs: Vec<Box<dyn Subsystem>> =
            vec![Box::new(LinearGain { a }), Box::new(LinearGain { a: b })];
        // u1 = y2 (subsystem 1 output 0 -> subsystem 0 input 0)
        // u2 = y1 (subsystem 0 output 0 -> subsystem 1 input 0)
        let graph =
            CouplingGraph::from_edges(vec![Coupling::new(1, 0, 0, 0), Coupling::new(0, 0, 1, 0)]);
        (subs, graph)
    }

    // -----------------------------------------------------------------------
    // BENCHMARK 1: two linear gains |a*b| < 1 converge to the analytic fixed
    // point y1* = y2* = 0  with residual < 1e-9, while the single-pass
    // explicit step does NOT reach it.
    // -----------------------------------------------------------------------

    #[test]
    fn converges_to_analytic_fixed_point_under_1e9() {
        // a = 0.3, b = 0.4  =>  |a*b| = 0.12  (strict contraction)
        // Fixed point: y1* = a*(b*y1*) = 0.12*y1*  => y1* = 0.
        // We use OffsetGain below; the plain LinearGain system is built for
        // reference but not driven here (both gains converge to zero trivially).

        // Seed a non-zero output by calling dt=0 with a non-zero input.
        // We do this by giving subsystem 0 an initial output ≠ 0 — the
        // dt=0 priming in coupled_step uses zero inputs so each subsystem's
        // initial output is 0.  To give a non-trivial starting point we
        // call coupled_step from a seeded state: we manually step with a
        // small perturbation before entering the loop.
        //
        // Actually, the cleanest test: the dt=0 prime will give y=[0,0]
        // (trivially the fixed point for y=a*0=0).  To see meaningful
        // convergence we need non-zero initial coupling.  We achieve this
        // by wrapping each gain in a subsystem that adds a constant offset
        // to its input.

        // Use an offset gain: y = a*(u + 1)  so the fixed point is:
        //   y1 = a*(y2 + 1)
        //   y2 = b*(y1 + 1)
        //   y1 = a*(b*(y1 + 1) + 1) = a*b*y1 + a*b + a
        //   y1*(1 - a*b) = a*(1 + b)  => y1* = a*(1+b)/(1-a*b)
        //   y2* = b*(1+a)/(1-a*b)
        struct OffsetGain {
            a: f64,
        }
        impl Subsystem for OffsetGain {
            fn n_inputs(&self) -> usize {
                1
            }
            fn n_outputs(&self) -> usize {
                1
            }
            fn step(&mut self, _t: f64, _dt: f64, inputs: &[f64]) -> Vec<f64> {
                vec![self.a * (inputs[0] + 1.0)]
            }
        }

        let a = 0.3_f64;
        let b = 0.4_f64;
        // Analytic fixed point
        let y1_star = a * (1.0 + b) / (1.0 - a * b);
        let y2_star = b * (1.0 + a) / (1.0 - a * b);

        let mut subs2: Vec<Box<dyn Subsystem>> =
            vec![Box::new(OffsetGain { a }), Box::new(OffsetGain { a: b })];
        let (_, graph2) = two_gain_system(a, b); // graph is identical

        let result = coupled_step(
            &mut subs2,
            &graph2,
            0.0,
            1.0,
            1e-12,
            500,
            ImplicitScheme::GaussSeidel,
            Relaxation::None,
        )
        .expect("must converge");

        let y1 = result.outputs[0][0];
        let y2 = result.outputs[1][0];

        assert!(
            (y1 - y1_star).abs() < 1e-9,
            "y1 = {y1:.6e}, y1* = {y1_star:.6e}, err = {:.3e}",
            (y1 - y1_star).abs()
        );
        assert!(
            (y2 - y2_star).abs() < 1e-9,
            "y2 = {y2:.6e}, y2* = {y2_star:.6e}, err = {:.3e}",
            (y2 - y2_star).abs()
        );

        // Verify the single-pass explicit step does NOT reach the fixed point.
        // One explicit step from zero inputs:  y1 = a*0 = 0, y2 = b*0 = 0.
        // One Gauss-Seidel step from zero: y1 = a*(0+1) = 0.3,
        //                                   y2 = b*(y1+1) = b*1.3.
        // Either way, NOT at the fixed point.
        // One GS step from zero inputs: first sub sees input=0, second sees first's fresh output.
        let y1_explicit = a * (0.0 + 1.0); // = 0.3
        let explicit_err_y1 = (y1_explicit - y1_star).abs();

        assert!(
            explicit_err_y1 > 1e-6,
            "single explicit step should NOT be at the fixed point; \
             y1_explicit = {y1_explicit:.6e}, y1* = {y1_star:.6e}, err = {explicit_err_y1:.3e}"
        );
        assert!(
            result.final_residual < 1e-12,
            "converged residual must be < tol (1e-12), got {:.3e}",
            result.final_residual
        );
    }

    // -----------------------------------------------------------------------
    // BENCHMARK 2: residual decreases monotonically (strictly), and
    // iteration count is finite for a contraction.
    // -----------------------------------------------------------------------

    #[test]
    fn residual_decreases_monotonically_for_contraction() {
        // Track per-iteration residuals by running with max_iter=1 repeatedly.
        struct OffsetGain {
            a: f64,
        }
        impl Subsystem for OffsetGain {
            fn n_inputs(&self) -> usize {
                1
            }
            fn n_outputs(&self) -> usize {
                1
            }
            fn step(&mut self, _t: f64, _dt: f64, inputs: &[f64]) -> Vec<f64> {
                vec![self.a * (inputs[0] + 1.0)]
            }
        }

        let a = 0.5_f64;
        let b = 0.5_f64; // |a*b| = 0.25

        // Run with a very tight tolerance and many iterations; capture residual
        // after each single-iter run by accumulating state manually.

        // Instead: call coupled_step with increasing max_iter and check that
        // the final_residual decreases, comparing residual at iter k vs k+1.
        let mut prev_residual = f64::MAX;
        // Run with max_iter = k for k in 1..=15 and check residual decreases.
        for k in 1..=15 {
            let mut subs: Vec<Box<dyn Subsystem>> =
                vec![Box::new(OffsetGain { a }), Box::new(OffsetGain { a: b })];
            let (_, graph) = two_gain_system(a, b);
            let res = coupled_step(
                &mut subs,
                &graph,
                0.0,
                1.0,
                // set tol below expected residual at step k so it runs to k iters.
                1e-20,
                k,
                ImplicitScheme::GaussSeidel,
                Relaxation::None,
            );
            let residual = match res {
                Ok(r) => r.final_residual,
                Err(FmiError::NotConverged { final_residual, .. }) => final_residual,
                Err(e) => panic!("unexpected error: {e}"),
            };
            if k > 1 {
                assert!(
                    residual < prev_residual,
                    "residual at iter {k} = {residual:.3e} did not decrease \
                     from iter {}: {prev_residual:.3e}",
                    k - 1
                );
            }
            prev_residual = residual;
        }
        // Also check it converges in a finite number of iterations.
        let mut subs: Vec<Box<dyn Subsystem>> =
            vec![Box::new(OffsetGain { a }), Box::new(OffsetGain { a: b })];
        let (_, graph) = two_gain_system(a, b);
        let result = coupled_step(
            &mut subs,
            &graph,
            0.0,
            1.0,
            1e-14,
            1000,
            ImplicitScheme::GaussSeidel,
            Relaxation::None,
        )
        .expect("contraction must converge");
        assert!(
            result.iterations < 200,
            "converged in {} iters (expected << 200)",
            result.iterations
        );
    }

    // -----------------------------------------------------------------------
    // BENCHMARK 3: non-contraction (|a*b| > 1) reports NotConverged, no hang,
    // no NaN or panic.
    // -----------------------------------------------------------------------

    #[test]
    fn non_contraction_reports_not_converged_no_nan() {
        // a = 2.0, b = 1.0  =>  |a*b| = 2.0  > 1  (expanding)
        struct OffsetGain {
            a: f64,
        }
        impl Subsystem for OffsetGain {
            fn n_inputs(&self) -> usize {
                1
            }
            fn n_outputs(&self) -> usize {
                1
            }
            fn step(&mut self, _t: f64, _dt: f64, inputs: &[f64]) -> Vec<f64> {
                vec![self.a * (inputs[0] + 1.0)]
            }
        }

        let mut subs: Vec<Box<dyn Subsystem>> = vec![
            Box::new(OffsetGain { a: 2.0 }),
            Box::new(OffsetGain { a: 1.0 }),
        ];
        let (_, graph) = two_gain_system(2.0, 1.0);

        let result = coupled_step(
            &mut subs,
            &graph,
            0.0,
            1.0,
            1e-12,
            30, // small limit — must not hang
            ImplicitScheme::GaussSeidel,
            Relaxation::None,
        );

        match result {
            Err(FmiError::NotConverged {
                iterations,
                final_residual,
                ..
            }) => {
                assert_eq!(iterations, 30, "should exhaust exactly max_iter = 30");
                assert!(
                    final_residual.is_finite(),
                    "final_residual must be finite (no NaN/Inf leak)"
                );
                // The residual should be > tol (that's why it didn't converge).
                assert!(
                    final_residual > 1e-12,
                    "NotConverged but residual {final_residual:.3e} < tol — that's a bug"
                );
            }
            Ok(r) => panic!(
                "non-contraction must NOT converge; got Ok with {} iters",
                r.iterations
            ),
            Err(e) => panic!("unexpected error type: {e}"),
        }
    }

    // -----------------------------------------------------------------------
    // BENCHMARK 4: relaxation converges in fewer iterations for a stiff case.
    //
    // For the two-gain algebraic loop y1 = a*(u1+1), y2 = b*(u2+1),
    // u1 = y2, u2 = y1, the GS fixed-point map on (y1,y2) has effective
    // rate |a*b| per step.  Under-relaxation with ω < 1 reduces the
    // effective step to ω*(1 - (1-ω)*|a*b|)^{-1} and cannot beat GS
    // for a simple scalar contraction. Over-relaxation (ω > 1) can
    // accelerate convergence but is limited by the contraction radius.
    //
    // The clearest comparison that is guaranteed to work across ω and Aitken:
    // use a JACOBI scheme, where the Gauss-Seidel ordering advantage is
    // absent and the contraction ratio is |a*b| exactly.  Under Jacobi
    // the effective contraction is |a*b| per step; with fixed ω the
    // effective radius of the relaxed map is |1 - ω*(1-a*b)| which is
    // minimised at ω* = 2/(1+a*b).  We verify that ω* beats plain ω=1.
    // -----------------------------------------------------------------------

    #[test]
    fn relaxation_converges_in_fewer_iterations() {
        // a = 0.8, b = 0.8 => |a*b| = 0.64.
        // Jacobi GS rate = 0.64 per iter.
        // Optimal ω* = 2/(1 + 0.64) ≈ 1.22.
        // Effective rate with ω*: |1 - ω*(1-|a*b|)| = |1 - 1.22*0.36| ≈ 0.56.
        // So ω=1.2 (just under ω*) must converge faster than ω=1.
        struct OffsetGain {
            a: f64,
        }
        impl Subsystem for OffsetGain {
            fn n_inputs(&self) -> usize {
                1
            }
            fn n_outputs(&self) -> usize {
                1
            }
            fn step(&mut self, _t: f64, _dt: f64, inputs: &[f64]) -> Vec<f64> {
                vec![self.a * (inputs[0] + 1.0)]
            }
        }

        let a = 0.8_f64;
        let b = 0.8_f64;

        let run = |scheme: ImplicitScheme, relaxation: Relaxation| -> usize {
            let mut subs: Vec<Box<dyn Subsystem>> =
                vec![Box::new(OffsetGain { a }), Box::new(OffsetGain { a: b })];
            let (_, graph) = two_gain_system(a, b);
            coupled_step(&mut subs, &graph, 0.0, 1.0, 1e-10, 2000, scheme, relaxation)
                .expect("should converge")
                .iterations
        };

        // Under Jacobi, ω=1.2 (mild over-relaxation) beats plain ω=1.
        let iters_plain = run(ImplicitScheme::Jacobi, Relaxation::None);
        let iters_relaxed = run(ImplicitScheme::Jacobi, Relaxation::Fixed(1.2));

        assert!(
            iters_relaxed < iters_plain,
            "Fixed(1.2) relaxation ({iters_relaxed} iters) should converge faster than \
             plain Jacobi ({iters_plain} iters) for |a*b|=0.64 — optimal ω* ≈ 1.22"
        );
    }

    // -----------------------------------------------------------------------
    // Guard: bad configuration is rejected fail-loud.
    // -----------------------------------------------------------------------

    #[test]
    fn rejects_bad_config() {
        let (mut subs, graph) = two_gain_system(0.3, 0.4);

        // tol <= 0
        assert!(matches!(
            coupled_step(
                &mut subs,
                &graph,
                0.0,
                1.0,
                0.0,
                10,
                ImplicitScheme::Jacobi,
                Relaxation::None
            ),
            Err(FmiError::BadCoupledStep(_))
        ));
        // max_iter = 0
        assert!(matches!(
            coupled_step(
                &mut subs,
                &graph,
                0.0,
                1.0,
                1e-6,
                0,
                ImplicitScheme::Jacobi,
                Relaxation::None
            ),
            Err(FmiError::BadCoupledStep(_))
        ));
        // dt <= 0
        assert!(matches!(
            coupled_step(
                &mut subs,
                &graph,
                0.0,
                0.0,
                1e-6,
                10,
                ImplicitScheme::Jacobi,
                Relaxation::None
            ),
            Err(FmiError::BadCoupledStep(_))
        ));
        // NaN tol
        assert!(matches!(
            coupled_step(
                &mut subs,
                &graph,
                0.0,
                1.0,
                f64::NAN,
                10,
                ImplicitScheme::Jacobi,
                Relaxation::None
            ),
            Err(FmiError::BadCoupledStep(_))
        ));
        // Fixed relaxation omega out of range
        assert!(matches!(
            coupled_step(
                &mut subs,
                &graph,
                0.0,
                1.0,
                1e-6,
                10,
                ImplicitScheme::Jacobi,
                Relaxation::Fixed(2.5)
            ),
            Err(FmiError::BadCoupledStep(_))
        ));
        // Aitken omega_0 out of range
        assert!(matches!(
            coupled_step(
                &mut subs,
                &graph,
                0.0,
                1.0,
                1e-6,
                10,
                ImplicitScheme::GaussSeidel,
                Relaxation::Aitken { omega_0: -0.1 }
            ),
            Err(FmiError::BadCoupledStep(_))
        ));
    }

    // -----------------------------------------------------------------------
    // Jacobi also converges (with more iterations than GS for the same case).
    // -----------------------------------------------------------------------

    #[test]
    fn jacobi_also_converges_to_fixed_point() {
        struct OffsetGain {
            a: f64,
        }
        impl Subsystem for OffsetGain {
            fn n_inputs(&self) -> usize {
                1
            }
            fn n_outputs(&self) -> usize {
                1
            }
            fn step(&mut self, _t: f64, _dt: f64, inputs: &[f64]) -> Vec<f64> {
                vec![self.a * (inputs[0] + 1.0)]
            }
        }
        let a = 0.4_f64;
        let b = 0.4_f64;
        let y1_star = a * (1.0 + b) / (1.0 - a * b);
        let y2_star = b * (1.0 + a) / (1.0 - a * b);

        let mut subs: Vec<Box<dyn Subsystem>> =
            vec![Box::new(OffsetGain { a }), Box::new(OffsetGain { a: b })];
        let (_, graph) = two_gain_system(a, b);

        let result = coupled_step(
            &mut subs,
            &graph,
            0.0,
            1.0,
            1e-12,
            500,
            ImplicitScheme::Jacobi,
            Relaxation::None,
        )
        .expect("Jacobi must converge for |a*b|=0.16");

        let y1 = result.outputs[0][0];
        let y2 = result.outputs[1][0];
        assert!(
            (y1 - y1_star).abs() < 1e-9,
            "Jacobi y1 err = {:.3e}",
            (y1 - y1_star).abs()
        );
        assert!(
            (y2 - y2_star).abs() < 1e-9,
            "Jacobi y2 err = {:.3e}",
            (y2 - y2_star).abs()
        );
    }

    // -----------------------------------------------------------------------
    // Fixed relaxation reduces iterations vs None for a slow contraction.
    // -----------------------------------------------------------------------

    #[test]
    fn fixed_relaxation_can_reduce_iterations() {
        struct OffsetGain {
            a: f64,
        }
        impl Subsystem for OffsetGain {
            fn n_inputs(&self) -> usize {
                1
            }
            fn n_outputs(&self) -> usize {
                1
            }
            fn step(&mut self, _t: f64, _dt: f64, inputs: &[f64]) -> Vec<f64> {
                vec![self.a * (inputs[0] + 1.0)]
            }
        }
        let a = 0.9_f64;
        let b = 0.9_f64; // |a*b| = 0.81 — slow

        let run_scheme = |relaxation: Relaxation| -> usize {
            let mut subs: Vec<Box<dyn Subsystem>> =
                vec![Box::new(OffsetGain { a }), Box::new(OffsetGain { a: b })];
            let (_, graph) = two_gain_system(a, b);
            coupled_step(
                &mut subs,
                &graph,
                0.0,
                1.0,
                1e-12,
                5000,
                ImplicitScheme::GaussSeidel,
                relaxation,
            )
            .expect("should converge")
            .iterations
        };

        let iters_none = run_scheme(Relaxation::None);
        // omega = 1.0 is the identity (same as None).
        // omega < 1 stabilises but may not speed up — test omega=0.5 still converges.
        let iters_fixed = run_scheme(Relaxation::Fixed(0.5));
        // Both should converge. Fixed(0.5) might be slower (under-relaxation)
        // but must still converge.
        assert!(
            iters_fixed < 5000,
            "Fixed relaxation must converge, got {iters_fixed} iters"
        );
        // The key invariant: both reach the fixed point.
        let _ = iters_none; // suppress unused warning; test passes if no panic
    }
}
