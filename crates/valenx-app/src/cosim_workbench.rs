//! The right-side **Co-Simulation Workbench** panel — a native front-end over
//! the in-house [`valenx_adapter_fmi`] crate (Valenx's FMI / HELICS-style
//! co-simulation coupling engine).
//!
//! Co-simulation couples several independently-integrated subsystems that
//! exchange scalar signals only at discrete *macro-steps*. A real deployment
//! couples heterogeneous solvers (an FMU exported from one tool, a controller
//! from another, an HLA/AFSIM federate). Loading an external binary FMU needs
//! a `.so`/`.dll` that cannot exist in the headless CI environment, so this
//! workbench instead drives the **real**
//! [`valenx_adapter_fmi::cosim::CoSimMaster`] / [`valenx_adapter_fmi::implicit`]
//! coordinator on a fully-native, fully-transparent **demo benchmark**: a
//! classic *two coupled mass-spring-dampers* problem.
//!
//! ```text
//! wall |--k1,c1--[ m1 ]--kc,cc--[ m2 ]
//!
//! m1 x1'' = -k1 x1 - c1 x1' + kc (x2 - x1) + cc (x2' - x1')
//! m2 x2'' =                   -kc (x2 - x1) - cc (x2' - x1')
//! ```
//!
//! The system is split into two [`valenx_adapter_fmi::cosim::Subsystem`]s:
//! subsystem **A** owns `(x1, v1)` and reads `(x2, v2)`; subsystem **B** owns
//! `(x2, v2)` and reads `(x1, v1)`. Each exposes `[position, velocity]` as its
//! two outputs and reads the partner's `[position, velocity]` as its two
//! inputs. The [`valenx_adapter_fmi::cosim::CouplingGraph`] wires
//! `A.out -> B.in` and `B.out -> A.in`. Each subsystem integrates itself with
//! its own fine RK4 substeps over a macro-step, holding the partner state
//! constant across the step (standard loosely-coupled co-simulation).
//!
//! The user sets the coupling parameters — macro-step `H`, number of steps,
//! the explicit [`valenx_adapter_fmi::cosim::Scheme`] (Jacobi vs
//! Gauss-Seidel), and (optionally) **implicit coupling** with its fixed-point
//! tolerance + iteration cap — clicks **Run**, and valenx-adapter-fmi then
//! advances the coupled system one macro-step at a time:
//!
//! * **Explicit** ([`valenx_adapter_fmi::cosim::CoSimMaster::advance`]): each
//!   subsystem steps exactly once per macro-step against the other's
//!   previous-step (Jacobi) or freshest-this-step (Gauss-Seidel) outputs.
//! * **Implicit** ([`valenx_adapter_fmi::implicit::coupled_step`]): within each
//!   macro-step the two subsystems are re-evaluated (rolled back between
//!   sweeps) until the coupling residual `‖Δy‖_∞ < tol`, recording the
//!   iteration count per step.
//!
//! At every macro-step the workbench records the four exchanged interface
//! signals `(x1, v1, x2, v2)`, the per-step coupling-iteration count (implicit
//! only), and a **coupling error** — the max position deviation of the co-sim
//! state from a *monolithic* reference (the full 4-state coupled ODE integrated
//! directly with RK4 at the co-sim's fine resolution, with no coupling lag).
//! The numbers are exactly what the coordinator returns — the workbench invents
//! none of them.
//!
//! Two readouts close the panel: a painter **time-series plot** of the four
//! exchanged signals over the horizon, plus a stats grid (final coupling
//! residual, total coupling iterations, final coupling error vs monolithic, and
//! a stability indicator).
//!
//! Mirrors the other workbenches (`uq_workbench`, `photogrammetry_workbench`):
//! a [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_cosim_workbench`], toggled from the View menu and
//! openable by the agent bridge under the workbench id `"cosim"` (aliases
//! `"co-simulation"` / `"fmi"`; see [`crate::project_tabs::TabKind`]). Every
//! numeric control is `.labelled_by` an accessible caption so the panel is
//! AI-drivable by name.
//!
//! Honesty: valenx-adapter-fmi's default path is the in-house native
//! `Subsystem` master — co-simulation *import only* (not model-exchange), and a
//! real binary FMU sits behind an off-by-default cargo feature this workbench
//! does not enable. The demo couples two textbook linear oscillators; it is a
//! research / educational illustration of co-sim coupling schemes, not a
//! certified multi-physics product. For a linear coupled system at small `H`
//! the co-sim solution tracks the monolithic/analytic solution within the
//! coupling tolerance, implicit Gauss-Seidel coupling converges (residual ->
//! 0, iteration count bounded), and a degenerate input (`H <= 0`, or zero
//! steps) surfaces an in-panel error — **not** a panic. The tests pin all
//! three.

use eframe::egui;
use valenx_adapter_fmi::cosim::{CoSimMaster, Coupling, CouplingGraph, Scheme, Subsystem};
use valenx_adapter_fmi::implicit::{coupled_step, ImplicitScheme, Relaxation};

use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Physical constants of the two-mass coupled spring-damper benchmark.
//
// These mirror the in-crate `spring_damper_validation` benchmark so the demo
// exercises the exact physical model the adapter is pinned against. The
// interface is the coupling spring/damper (kc, cc) between the masses.
// ---------------------------------------------------------------------------

/// Mass of body 1 (kg).
const M1: f64 = 1.0;
/// Mass of body 2 (kg).
const M2: f64 = 1.5;
/// Wall spring stiffness on mass 1 (N/m).
const K1: f64 = 30.0;
/// Wall damper on mass 1 (N·s/m).
const C1: f64 = 0.4;
/// Coupling spring stiffness between the two masses (N/m).
const KC: f64 = 20.0;
/// Coupling damper between the two masses (N·s/m).
const CC: f64 = 0.3;

/// Initial position of mass 1 (m).
const X1_0: f64 = 1.0;
/// Initial position of mass 2 (m).
const X2_0: f64 = -0.5;

/// Fine RK4 substeps taken inside each subsystem per macro-step. The
/// monolithic reference uses the same fine step (`H / SUBSTEPS`) so the RK4
/// truncation error cancels in the comparison and what remains is the pure
/// inter-subsystem coupling error.
const SUBSTEPS: usize = 10;

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// The explicit-coupling sweep order the master uses when implicit coupling is
/// off. A thin UI mirror of [`valenx_adapter_fmi::cosim::Scheme`] (so the type
/// can derive the traits the egui `selectable_value` widget needs and so the
/// label text lives here).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CouplingScheme {
    /// Subsystems step in index order; each sees the freshest partner output
    /// produced earlier within the same macro-step (less coupling lag).
    GaussSeidel,
    /// Every subsystem steps against the previous macro-step's outputs (a
    /// consistent one-step-old picture; parallelizable, one-step lag).
    Jacobi,
}

impl CouplingScheme {
    /// Human-readable label for the combo box / status line.
    fn label(self) -> &'static str {
        match self {
            CouplingScheme::GaussSeidel => "Gauss-Seidel (sequential)",
            CouplingScheme::Jacobi => "Jacobi (parallel)",
        }
    }

    /// Map to the explicit master's scheme.
    fn to_scheme(self) -> Scheme {
        match self {
            CouplingScheme::GaussSeidel => Scheme::GaussSeidel,
            CouplingScheme::Jacobi => Scheme::Jacobi,
        }
    }

    /// Map to the implicit coupler's sweep order.
    fn to_implicit(self) -> ImplicitScheme {
        match self {
            CouplingScheme::GaussSeidel => ImplicitScheme::GaussSeidel,
            CouplingScheme::Jacobi => ImplicitScheme::Jacobi,
        }
    }
}

/// Editable co-simulation inputs shown in the workbench.
#[derive(Clone, Copy, Debug)]
pub struct CosimParams {
    /// Macro-step size `H` (s): how often the two subsystems exchange signals.
    /// The coupling error of an explicit scheme is first-order in `H`. Must be
    /// `> 0`.
    pub macro_step: f64,
    /// Number of macro-steps to advance (the horizon is `num_steps * H`). Must
    /// be `>= 1`.
    pub num_steps: usize,
    /// Which explicit coupling sweep order to use (also selects the implicit
    /// sweep order when implicit coupling is on).
    pub scheme: CouplingScheme,
    /// Whether to use **implicit** (iterative, strongly-coupled) coupling: each
    /// macro-step is iterated to a converged coupling residual instead of a
    /// single explicit exchange.
    pub implicit: bool,
    /// Fixed-point convergence tolerance on `‖Δy‖_∞` for implicit coupling.
    /// Must be `> 0`.
    pub tol: f64,
    /// Maximum fixed-point iterations per macro-step (implicit coupling). Must
    /// be `>= 1`.
    pub max_iters: usize,
}

impl Default for CosimParams {
    fn default() -> Self {
        Self {
            // A modest, well-conditioned default: a small macro-step (so both
            // schemes track the monolithic reference closely) over a 1 s
            // horizon, explicit Gauss-Seidel coupling.
            macro_step: 2.0e-3,
            num_steps: 500,
            scheme: CouplingScheme::GaussSeidel,
            implicit: false,
            tol: 1.0e-10,
            max_iters: 50,
        }
    }
}

// ---------------------------------------------------------------------------
// Simulation result
// ---------------------------------------------------------------------------

/// One recorded macro-step: the simulation time and the four exchanged
/// interface signals at the end of the step.
#[derive(Clone, Copy, Debug)]
pub struct StepSample {
    /// Simulation time at the end of this macro-step (s).
    pub t: f64,
    /// Position of mass 1 (m).
    pub x1: f64,
    /// Velocity of mass 1 (m/s).
    pub v1: f64,
    /// Position of mass 2 (m).
    pub x2: f64,
    /// Velocity of mass 2 (m/s).
    pub v2: f64,
}

/// The full recorded co-simulation history + summary metrics. Everything here
/// is read straight from the coordinator's outputs.
#[derive(Default, Clone)]
pub struct CosimResult {
    /// Per-macro-step exchanged-signal history (one entry per step).
    pub samples: Vec<StepSample>,
    /// Per-macro-step coupling-iteration count (implicit coupling only; empty
    /// for explicit, where each step is exactly one exchange).
    pub iters_per_step: Vec<usize>,
    /// Total coupling iterations summed over the horizon (implicit only).
    pub total_iterations: usize,
    /// The final macro-step's coupling residual `‖Δy‖_∞` (implicit only; `0`
    /// for explicit, which performs no iteration).
    pub final_residual: f64,
    /// Largest coupling residual seen over any macro-step (implicit only).
    pub max_residual: f64,
    /// The final coupling error: the max absolute position deviation of the
    /// co-sim end state `(x1, x2)` from the monolithic reference integrated to
    /// the same horizon.
    pub final_coupling_error: f64,
    /// Whether the trajectory stayed bounded (no non-finite or runaway state)
    /// — a simple stability indicator for the readout.
    pub stable: bool,
    /// Whether implicit coupling was used for this run (affects which readouts
    /// are meaningful).
    pub implicit: bool,
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the co-simulation workbench.
#[derive(Default)]
pub struct CosimWorkbenchState {
    /// User-editable parameters.
    pub params: CosimParams,
    /// Last successful result (populated after a successful run).
    pub result: Option<CosimResult>,
    /// Status / error line shown below the controls.
    pub status: String,
}

// ---------------------------------------------------------------------------
// The two coupled subsystems (implement the real `Subsystem` trait)
// ---------------------------------------------------------------------------

/// Mass-1 subsystem. State `(x1, v1)`; inputs `(x2, v2)`; outputs `(x1, v1)`.
struct MassA {
    x: f64,
    v: f64,
}

impl MassA {
    /// Acceleration of mass 1 given its state and the partner state.
    fn accel(x1: f64, v1: f64, x2: f64, v2: f64) -> f64 {
        (-K1 * x1 - C1 * v1 + KC * (x2 - x1) + CC * (v2 - v1)) / M1
    }
}

impl Subsystem for MassA {
    fn n_inputs(&self) -> usize {
        2
    }
    fn n_outputs(&self) -> usize {
        2
    }
    fn step(&mut self, _t: f64, dt: f64, inputs: &[f64]) -> Vec<f64> {
        // Partner state held constant across the macro-step (standard
        // loosely-coupled co-simulation); dt == 0 only samples the outputs.
        let (x2, v2) = (inputs[0], inputs[1]);
        if dt > 0.0 {
            let h = dt / SUBSTEPS as f64;
            for _ in 0..SUBSTEPS {
                rk4_1dof(&mut self.x, &mut self.v, h, |x, v| {
                    MassA::accel(x, v, x2, v2)
                });
            }
        }
        vec![self.x, self.v]
    }
    fn state(&self) -> Vec<f64> {
        vec![self.x, self.v]
    }
    fn set_state(&mut self, s: &[f64]) {
        self.x = s[0];
        self.v = s[1];
    }
}

/// Mass-2 subsystem. State `(x2, v2)`; inputs `(x1, v1)`; outputs `(x2, v2)`.
struct MassB {
    x: f64,
    v: f64,
}

impl MassB {
    /// Acceleration of mass 2 given its state and the partner state.
    fn accel(x2: f64, v2: f64, x1: f64, v1: f64) -> f64 {
        (-KC * (x2 - x1) - CC * (v2 - v1)) / M2
    }
}

impl Subsystem for MassB {
    fn n_inputs(&self) -> usize {
        2
    }
    fn n_outputs(&self) -> usize {
        2
    }
    fn step(&mut self, _t: f64, dt: f64, inputs: &[f64]) -> Vec<f64> {
        let (x1, v1) = (inputs[0], inputs[1]);
        if dt > 0.0 {
            let h = dt / SUBSTEPS as f64;
            for _ in 0..SUBSTEPS {
                rk4_1dof(&mut self.x, &mut self.v, h, |x, v| {
                    MassB::accel(x, v, x1, v1)
                });
            }
        }
        vec![self.x, self.v]
    }
    fn state(&self) -> Vec<f64> {
        vec![self.x, self.v]
    }
    fn set_state(&mut self, s: &[f64]) {
        self.x = s[0];
        self.v = s[1];
    }
}

/// One classic RK4 step of a single second-order DOF `x'' = a(x, v)`.
fn rk4_1dof(x: &mut f64, v: &mut f64, h: f64, accel: impl Fn(f64, f64) -> f64) {
    let (x0, v0) = (*x, *v);

    let a1 = accel(x0, v0);
    let (k1x, k1v) = (v0, a1);

    let a2 = accel(x0 + 0.5 * h * k1x, v0 + 0.5 * h * k1v);
    let (k2x, k2v) = (v0 + 0.5 * h * k1v, a2);

    let a3 = accel(x0 + 0.5 * h * k2x, v0 + 0.5 * h * k2v);
    let (k3x, k3v) = (v0 + 0.5 * h * k2v, a3);

    let a4 = accel(x0 + h * k3x, v0 + h * k3v);
    let (k4x, k4v) = (v0 + h * k3v, a4);

    *x = x0 + (h / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
    *v = v0 + (h / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
}

/// The coupling graph wiring `A.out -> B.in` and `B.out -> A.in`:
/// A outputs `(x1, v1)` into B's inputs, B outputs `(x2, v2)` into A's inputs.
fn coupling_graph() -> CouplingGraph {
    CouplingGraph::from_edges(vec![
        Coupling::new(0, 0, 1, 0), // A.x1 -> B.in[0]
        Coupling::new(0, 1, 1, 1), // A.v1 -> B.in[1]
        Coupling::new(1, 0, 0, 0), // B.x2 -> A.in[0]
        Coupling::new(1, 1, 0, 1), // B.v2 -> A.in[1]
    ])
}

/// Build the two fresh subsystems at their initial state.
fn build_subsystems() -> Vec<Box<dyn Subsystem>> {
    vec![
        Box::new(MassA { x: X1_0, v: 0.0 }),
        Box::new(MassB { x: X2_0, v: 0.0 }),
    ]
}

// ---------------------------------------------------------------------------
// Monolithic reference (no coupling lag) — the analytic-grade ground truth
// ---------------------------------------------------------------------------

/// Monolithic reference: integrate the full 4-state coupled system
/// `[x1, v1, x2, v2]` directly with RK4 at the co-sim's fine resolution
/// (`macro_step / SUBSTEPS`), with no inter-subsystem coupling lag. Returns the
/// state at `t = num_steps * macro_step`.
fn monolithic_reference(macro_step: f64, num_steps: usize) -> [f64; 4] {
    let mut s = [X1_0, 0.0, X2_0, 0.0];
    let h = macro_step / SUBSTEPS as f64;
    let n = num_steps * SUBSTEPS;

    let deriv = |s: &[f64; 4]| -> [f64; 4] {
        let (x1, v1, x2, v2) = (s[0], s[1], s[2], s[3]);
        let a1 = (-K1 * x1 - C1 * v1 + KC * (x2 - x1) + CC * (v2 - v1)) / M1;
        let a2 = (-KC * (x2 - x1) - CC * (v2 - v1)) / M2;
        [v1, a1, v2, a2]
    };
    let add = |a: &[f64; 4], b: &[f64; 4]| [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]];
    let scale = |a: &[f64; 4], k: f64| [a[0] * k, a[1] * k, a[2] * k, a[3] * k];

    for _ in 0..n {
        let k1 = deriv(&s);
        let s2 = add(&s, &scale(&k1, 0.5 * h));
        let k2 = deriv(&s2);
        let s3 = add(&s, &scale(&k2, 0.5 * h));
        let k3 = deriv(&s3);
        let s4 = add(&s, &scale(&k3, h));
        let k4 = deriv(&s4);
        for i in 0..4 {
            s[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
    s
}

// ---------------------------------------------------------------------------
// Run the co-simulation through the REAL valenx-adapter-fmi coordinator
// ---------------------------------------------------------------------------

impl CosimWorkbenchState {
    /// Step the demo co-simulation through the **real** valenx-adapter-fmi
    /// coordinator and collect the per-step exchanged-signal histories, the
    /// per-step coupling-iteration counts (implicit), and the coupling-error /
    /// stability metrics, fail-loud.
    ///
    /// Every failure path returns an `Err(String)` — never a panic, never an
    /// invented number. Degenerate inputs (`macro_step <= 0`, `num_steps == 0`,
    /// or — for implicit coupling — `tol <= 0` / `max_iters == 0`) are rejected
    /// up front, and a non-converging implicit step surfaces the coordinator's
    /// own [`valenx_adapter_fmi::error::FmiError::NotConverged`] message
    /// verbatim.
    pub fn run(&self) -> Result<CosimResult, String> {
        let p = &self.params;

        // --- Degenerate-input guards (fail loud, no panic) ------------------
        if !p.macro_step.is_finite() || p.macro_step <= 0.0 {
            return Err(format!(
                "macro-step H must be finite and > 0 (got {})",
                p.macro_step
            ));
        }
        if p.num_steps == 0 {
            return Err("number of macro-steps must be >= 1 (got 0)".to_string());
        }
        if p.implicit {
            if !p.tol.is_finite() || p.tol <= 0.0 {
                return Err(format!(
                    "implicit coupling tolerance must be finite and > 0 (got {})",
                    p.tol
                ));
            }
            if p.max_iters == 0 {
                return Err("implicit coupling max-iters must be >= 1 (got 0)".to_string());
            }
        }

        if p.implicit {
            self.run_implicit(p)
        } else {
            self.run_explicit(p)
        }
    }

    /// Explicit (single-exchange-per-step) coupling via the real
    /// [`CoSimMaster`].
    fn run_explicit(&self, p: &CosimParams) -> Result<CosimResult, String> {
        let mut master =
            CoSimMaster::new(build_subsystems(), coupling_graph(), p.scheme.to_scheme())
                .map_err(|e| format!("co-sim master rejected the coupling graph: {e}"))?;

        let mut samples = Vec::with_capacity(p.num_steps);
        let mut stable = true;
        for _ in 0..p.num_steps {
            master.advance(p.macro_step);
            let a = master
                .outputs_of(0)
                .ok_or_else(|| "co-sim master lost subsystem A outputs".to_string())?;
            let b = master
                .outputs_of(1)
                .ok_or_else(|| "co-sim master lost subsystem B outputs".to_string())?;
            let sample = StepSample {
                t: master.time(),
                x1: a[0],
                v1: a[1],
                x2: b[0],
                v2: b[1],
            };
            if !sample_is_finite(&sample) {
                stable = false;
                break;
            }
            samples.push(sample);
        }

        let final_coupling_error = self.coupling_error(&samples, p, stable);

        Ok(CosimResult {
            samples,
            iters_per_step: Vec::new(),
            total_iterations: 0,
            final_residual: 0.0,
            max_residual: 0.0,
            final_coupling_error,
            stable,
            implicit: false,
        })
    }

    /// Implicit (strongly-coupled, iterated-to-convergence) coupling via the
    /// real [`coupled_step`]. We drive the macro-step loop ourselves, carrying
    /// each subsystem's converged state forward between steps, so we can record
    /// the per-step iteration count + residual.
    fn run_implicit(&self, p: &CosimParams) -> Result<CosimResult, String> {
        let mut subs = build_subsystems();
        let graph = coupling_graph();
        // Validate the graph once up front (fail-loud), mirroring what the
        // explicit master's constructor does — `coupled_step` only
        // debug-asserts edge ranges.
        graph
            .validate(&subs)
            .map_err(|e| format!("co-sim coupling graph is invalid: {e}"))?;

        let mut samples = Vec::with_capacity(p.num_steps);
        let mut iters_per_step = Vec::with_capacity(p.num_steps);
        let mut total_iterations = 0usize;
        let mut final_residual = 0.0;
        let mut max_residual: f64 = 0.0;
        let mut stable = true;
        let mut t = 0.0;

        for step in 0..p.num_steps {
            let res = coupled_step(
                &mut subs,
                &graph,
                t,
                p.macro_step,
                p.tol,
                p.max_iters,
                p.scheme.to_implicit(),
                Relaxation::None,
            )
            .map_err(|e| {
                format!(
                    "implicit coupling failed at macro-step {} (t = {:.4} s): {e}",
                    step + 1,
                    t
                )
            })?;

            // `coupled_step` leaves each subsystem rolled back to its
            // start-of-step state, so commit the converged outputs as the new
            // continuous state before advancing time.
            for (sub, out) in subs.iter_mut().zip(res.outputs.iter()) {
                sub.set_state(out);
            }

            t += p.macro_step;
            total_iterations += res.iterations;
            final_residual = res.final_residual;
            max_residual = max_residual.max(res.final_residual);

            let sample = StepSample {
                t,
                x1: res.outputs[0][0],
                v1: res.outputs[0][1],
                x2: res.outputs[1][0],
                v2: res.outputs[1][1],
            };
            if !sample_is_finite(&sample) {
                stable = false;
                break;
            }
            samples.push(sample);
            iters_per_step.push(res.iterations);
        }

        let final_coupling_error = self.coupling_error(&samples, p, stable);

        Ok(CosimResult {
            samples,
            iters_per_step,
            total_iterations,
            final_residual,
            max_residual,
            final_coupling_error,
            stable,
            implicit: true,
        })
    }

    /// The final coupling error: the max absolute position deviation of the
    /// co-sim end state `(x1, x2)` from the monolithic reference integrated to
    /// the same number of steps the co-sim actually completed (so a truncated /
    /// blown-up run compares against the matching horizon). Returns `+inf` if
    /// the run produced no finite samples.
    fn coupling_error(&self, samples: &[StepSample], p: &CosimParams, stable: bool) -> f64 {
        let Some(last) = samples.last() else {
            return f64::INFINITY;
        };
        if !stable {
            return f64::INFINITY;
        }
        let completed = samples.len();
        let reference = monolithic_reference(p.macro_step, completed);
        (last.x1 - reference[0])
            .abs()
            .max((last.x2 - reference[2]).abs())
    }
}

/// Whether all four signals in a sample are finite.
fn sample_is_finite(s: &StepSample) -> bool {
    s.x1.is_finite() && s.v1.is_finite() && s.x2.is_finite() && s.v2.is_finite()
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the co-simulation workbench. A no-op unless toggled on via
/// View → Co-Simulation.
///
/// Mirrors [`crate::photogrammetry_workbench::draw_photogrammetry_workbench`].
pub fn draw_cosim_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_cosim_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_cosim_workbench",
        "Co-Simulation (FMI / HELICS)",
        cosim_workbench_body,
    );
    if close {
        app.show_cosim_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn cosim_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "Co-simulation \u{2014} two coupled mass-spring-dampers, split into two subsystems \
             that exchange position/velocity at each macro-step, advanced through the REAL \
             in-house valenx-adapter-fmi coordinator (Jacobi / Gauss-Seidel explicit coupling, \
             or a strongly-coupled fixed-point implicit coupler). The exchanged signals + \
             coupling residual + error vs a monolithic reference are exactly what the \
             coordinator returns. [research / educational \u{2014} native co-sim master; FMI \
             import only]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let mut do_run = false;

    {
        let s = &mut app.cosim;
        let p = &mut s.params;

        ui.label(egui::RichText::new("Coupling").strong());
        egui::Grid::new("cosim_coupling_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("macro-step H (s)");
                ui.add(
                    egui::DragValue::new(&mut p.macro_step)
                        .speed(1.0e-4)
                        .range(0.0..=1.0)
                        .max_decimals(6),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Macro-step size H: how often the two subsystems exchange signals. The \
                     coupling error of an explicit scheme is first-order in H, so a smaller H \
                     tracks the monolithic reference more closely. Must be > 0.",
                );
                ui.end_row();

                let lbl = ui.label("macro-steps");
                ui.add(
                    egui::DragValue::new(&mut p.num_steps)
                        .speed(10)
                        .range(0..=20000),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Number of macro-steps to advance. The horizon is num_steps * H. Must \
                     be >= 1.",
                );
                ui.end_row();

                let lbl = ui.label("coupling scheme");
                egui::ComboBox::from_id_source("cosim_scheme_combo")
                    .selected_text(p.scheme.label())
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut p.scheme,
                            CouplingScheme::GaussSeidel,
                            CouplingScheme::GaussSeidel.label(),
                        );
                        ui.selectable_value(
                            &mut p.scheme,
                            CouplingScheme::Jacobi,
                            CouplingScheme::Jacobi.label(),
                        );
                    })
                    .response
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Explicit coupling sweep order (also selects the sweep order used \
                         inside each implicit iteration). Gauss-Seidel carries less coupling \
                         lag than Jacobi at the same macro-step.",
                    );
                ui.end_row();

                let lbl = ui.label("implicit coupling");
                ui.checkbox(&mut p.implicit, "")
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "When on, each macro-step is iterated (subsystems rolled back between \
                         sweeps) to a converged coupling residual instead of a single explicit \
                         exchange. Strongly-coupled / algebraic-loop systems need this.",
                    );
                ui.end_row();
            });

        // Implicit-only controls. Shown always (greyed when implicit is off) so
        // the form layout is stable and the controls keep stable accessible
        // names for the a11y tree (mirrors uq_workbench's enabled-ui pattern).
        ui.add_space(4.0);
        ui.label(egui::RichText::new("Implicit fixed-point").strong());
        egui::Grid::new("cosim_implicit_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let implicit = p.implicit;
                ui.add_enabled_ui(implicit, |ui| {
                    let lbl = ui.label("tolerance");
                    ui.add(
                        egui::DragValue::new(&mut p.tol)
                            .speed(1.0e-11)
                            .range(0.0..=1.0)
                            .max_decimals(12),
                    )
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Fixed-point convergence tolerance on the infinity-norm of the change \
                         in the exchanged signals. The step converges when ‖Δy‖_∞ < tol. \
                         Must be > 0.",
                    );
                });
                ui.end_row();

                ui.add_enabled_ui(implicit, |ui| {
                    let lbl = ui.label("max iters / step");
                    ui.add(
                        egui::DragValue::new(&mut p.max_iters)
                            .speed(1)
                            .range(0..=1000),
                    )
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Maximum fixed-point iterations per macro-step. If the coupling does \
                         not converge within this many iterations the run fails loud with the \
                         coordinator's NotConverged error. Must be >= 1.",
                    );
                });
                ui.end_row();
            });

        ui.add_space(6.0);
        ui.horizontal(|ui| {
            if ui
                .button(egui::RichText::new("Run").strong())
                .on_hover_text(
                    "Advance the coupled two-mass system through the valenx-adapter-fmi \
                     coordinator and report the exchanged-signal history, coupling residual / \
                     iteration counts, and the coupling error vs the monolithic reference.",
                )
                .clicked()
            {
                do_run = true;
            }
        });
    }

    // --- Execute (outside the params borrow) --------------------------------
    if do_run {
        run_and_store(app);
    }

    // --- Status line ---------------------------------------------------------
    let s = &app.cosim;
    if !s.status.is_empty() {
        ui.add_space(6.0);
        let color = if s.status.starts_with('\u{26A0}') {
            egui::Color32::from_rgb(220, 120, 60)
        } else {
            egui::Color32::from_rgb(90, 180, 110)
        };
        ui.label(egui::RichText::new(&s.status).color(color).strong());
    }

    // --- Visualisation -------------------------------------------------------
    ui.add_space(6.0);
    ui.separator();
    draw_cosim_viz(s, ui);
}

/// Run the co-simulation and fold the result (or error) into the workbench
/// status. Factored out so the Run button (and tests) can share it.
pub(crate) fn run_and_store(app: &mut ValenxApp) {
    let s = &mut app.cosim;
    match s.run() {
        Ok(res) => {
            let stab = if res.stable { "stable" } else { "UNSTABLE" };
            if res.implicit {
                s.status = format!(
                    "\u{2714} {} steps \u{00B7} implicit {} \u{00B7} {} total iters \u{00B7} \
                     final residual {:.2e} \u{00B7} coupling err {:.3e} \u{00B7} {}",
                    res.samples.len(),
                    s.params.scheme.label(),
                    res.total_iterations,
                    res.final_residual,
                    res.final_coupling_error,
                    stab,
                );
            } else {
                s.status = format!(
                    "\u{2714} {} steps \u{00B7} explicit {} \u{00B7} coupling err {:.3e} vs \
                     monolithic \u{00B7} {}",
                    res.samples.len(),
                    s.params.scheme.label(),
                    res.final_coupling_error,
                    stab,
                );
            }
            s.result = Some(res);
        }
        Err(e) => {
            s.status = format!("\u{26A0} {e}");
            s.result = None;
        }
    }
}

// ---------------------------------------------------------------------------
// 2-D visualisation (painter time-series of the four exchanged signals)
// ---------------------------------------------------------------------------

fn draw_cosim_viz(s: &CosimWorkbenchState, ui: &mut egui::Ui) {
    let Some(res) = &s.result else {
        ui.label(
            egui::RichText::new(
                "press \"Run\" to advance the coupled co-simulation and plot the exchanged \
                 interface signals over the horizon",
            )
            .weak(),
        );
        return;
    };

    ui.label(egui::RichText::new("Exchanged interface signals").strong());
    ui.label(
        egui::RichText::new(
            "cyan = x1 \u{00B7} sky = v1 \u{00B7} amber = x2 \u{00B7} red = v2 \u{00B7} \
             horizontal axis = simulation time",
        )
        .weak()
        .small(),
    );

    draw_timeseries(res, ui);

    // Readouts grid below the plot.
    ui.add_space(6.0);
    egui::Grid::new("cosim_stats")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(v);
                ui.end_row();
            };
            row(
                ui,
                "macro-steps completed",
                format!("{}", res.samples.len()),
            );
            row(
                ui,
                "coupling",
                if res.implicit {
                    format!("implicit ({})", s.params.scheme.label())
                } else {
                    format!("explicit ({})", s.params.scheme.label())
                },
            );
            if res.implicit {
                row(
                    ui,
                    "total coupling iterations",
                    format!("{}", res.total_iterations),
                );
                row(
                    ui,
                    "final / max residual",
                    format!("{:.3e} / {:.3e}", res.final_residual, res.max_residual),
                );
            }
            row(
                ui,
                "coupling error vs monolithic",
                format!("{:.4e}", res.final_coupling_error),
            );
            row(
                ui,
                "stability",
                if res.stable {
                    "stable (bounded)".to_string()
                } else {
                    "UNSTABLE (diverged)".to_string()
                },
            );
        });
}

/// Draw the time-series plot of the four exchanged signals with the egui
/// painter (auto-scaled to fit the data bounds). Self-contained — no plotting
/// dependency.
fn draw_timeseries(res: &CosimResult, ui: &mut egui::Ui) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(460.0, 220.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(14, 22, 34));

    if res.samples.len() < 2 {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "need >= 2 macro-steps to plot",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    // Data bounds: t spans [0, t_last]; the value axis spans min/max over all
    // four signals across every sample. Guard against a degenerate (zero-
    // extent) range so a flat trace still draws.
    let t_lo = 0.0_f64;
    let t_hi = res.samples.last().map(|s| s.t).unwrap_or(1.0).max(1e-9);
    let mut v_lo = f64::INFINITY;
    let mut v_hi = f64::NEG_INFINITY;
    for s in &res.samples {
        for val in [s.x1, s.v1, s.x2, s.v2] {
            if val.is_finite() {
                v_lo = v_lo.min(val);
                v_hi = v_hi.max(val);
            }
        }
    }
    if !(v_lo.is_finite() && v_hi.is_finite()) {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "no finite signal to plot",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }
    if (v_hi - v_lo).abs() < 1e-9 {
        v_lo -= 0.5;
        v_hi += 0.5;
    }

    let margin = 14.0_f32;
    let inner = rect.shrink(margin);
    let t_span = (t_hi - t_lo) as f32;
    let v_span = (v_hi - v_lo) as f32;

    // Map (t, value) to a painter position (value axis flipped so +value is up).
    let to_screen = |t: f64, val: f64| -> egui::Pos2 {
        let ft = ((t - t_lo) as f32 / t_span).clamp(0.0, 1.0);
        let fv = ((val - v_lo) as f32 / v_span).clamp(0.0, 1.0);
        egui::pos2(
            inner.left() + ft * inner.width(),
            inner.bottom() - fv * inner.height(),
        )
    };

    // Zero baseline (if 0 is within the value range).
    if v_lo < 0.0 && v_hi > 0.0 {
        let y0 = to_screen(t_lo, 0.0).y;
        painter.line_segment(
            [egui::pos2(inner.left(), y0), egui::pos2(inner.right(), y0)],
            egui::Stroke::new(1.0, egui::Color32::from_gray(60)),
        );
    }

    // One polyline per exchanged signal.
    let polyline = |accessor: &dyn Fn(&StepSample) -> f64, color: egui::Color32| {
        let pts: Vec<egui::Pos2> = res
            .samples
            .iter()
            .map(|s| to_screen(s.t, accessor(s)))
            .collect();
        painter.add(egui::Shape::line(pts, egui::Stroke::new(1.4, color)));
    };
    polyline(&|s| s.x1, egui::Color32::from_rgb(70, 200, 210));
    polyline(&|s| s.v1, egui::Color32::from_rgb(120, 170, 230));
    polyline(&|s| s.x2, egui::Color32::from_rgb(230, 180, 70));
    polyline(&|s| s.v2, egui::Color32::from_rgb(230, 110, 90));
}

// ---------------------------------------------------------------------------
// Tests (unit + headless_ui_tests, mirroring photogrammetry_workbench)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_run_succeeds_and_is_populated() {
        let s = CosimWorkbenchState::default();
        let res = s.run().expect("default co-sim run should succeed");
        assert_eq!(
            res.samples.len(),
            s.params.num_steps,
            "every requested macro-step is recorded"
        );
        assert!(res.stable, "the default linear system stays bounded");
        assert!(
            res.final_coupling_error.is_finite(),
            "coupling error must be finite for a stable run"
        );
        // Explicit run: no iteration recorded.
        assert!(res.iters_per_step.is_empty());
        assert_eq!(res.total_iterations, 0);
    }

    #[test]
    fn cosim_tracks_monolithic_reference_within_tol_pin() {
        // PIN (analytic): for the linear coupled system the explicit
        // Gauss-Seidel co-sim end state tracks the monolithic / analytic-grade
        // reference (the full 4-state ODE integrated at the co-sim's fine
        // resolution with no coupling lag), and — because the explicit coupling
        // error is FIRST-ORDER in the macro-step H — that error shrinks with H.
        // We pin both: (a) at a small H the error is below a tight tolerance,
        // and (b) halving H roughly halves the error (first-order convergence).
        let err_at = |h: f64, n: usize| -> f64 {
            let mut s = CosimWorkbenchState::default();
            s.params.scheme = CouplingScheme::GaussSeidel;
            s.params.macro_step = h;
            s.params.num_steps = n; // same 2 s horizon for both
            s.params.implicit = false;
            s.run()
                .expect("explicit GS run should succeed")
                .final_coupling_error
        };

        // (a) Small macro-step tracks the monolithic reference tightly.
        let err_fine = err_at(2.0e-4, 10_000); // T = 2 s
        assert!(
            err_fine < 1.0e-3,
            "Gauss-Seidel coupling error {err_fine:.3e} at H=2e-4 must be < 1e-3 vs monolithic"
        );

        // (b) First-order convergence: doubling H (same horizon) grows the
        // error by roughly 2x (allow a generous band — it is provably O(H)).
        let err_coarse = err_at(4.0e-4, 5_000); // T = 2 s, H doubled
        let ratio = err_coarse / err_fine;
        assert!(
            (1.5..=3.0).contains(&ratio),
            "doubling H should ~double the coupling error (first-order); ratio = {ratio:.2} \
             (fine {err_fine:.3e}, coarse {err_coarse:.3e})"
        );
    }

    #[test]
    fn jacobi_coupling_error_exceeds_gauss_seidel() {
        // Documented expectation (matches the adapter's own benchmark, which
        // samples the T = 2 s horizon): Jacobi carries a one-macro-step coupling
        // lag, so its end-state error against the monolithic reference exceeds
        // Gauss-Seidel's at the same macro-step. (The single-sampled error
        // oscillates with the dynamics, so this is pinned at the same 2 s
        // horizon the adapter validates at, where the ordering is robust.)
        let run = |scheme: CouplingScheme| -> f64 {
            let mut s = CosimWorkbenchState::default();
            s.params.scheme = scheme;
            s.params.macro_step = 1.0e-3;
            s.params.num_steps = 2000; // T = 2 s
            s.params.implicit = false;
            s.run().expect("run should succeed").final_coupling_error
        };
        let err_gs = run(CouplingScheme::GaussSeidel);
        let err_jac = run(CouplingScheme::Jacobi);
        assert!(
            err_jac > err_gs,
            "Jacobi error {err_jac:.3e} should exceed Gauss-Seidel error {err_gs:.3e} \
             (Jacobi has a one-macro-step coupling lag)"
        );
    }

    #[test]
    fn implicit_gauss_seidel_converges_with_bounded_iterations_pin() {
        // PIN: implicit Gauss-Seidel coupling converges every macro-step
        // (residual -> 0, below tol) with a bounded per-step iteration count,
        // and tracks the monolithic reference at least as well as the explicit
        // scheme.
        let mut s = CosimWorkbenchState::default();
        s.params.scheme = CouplingScheme::GaussSeidel;
        s.params.macro_step = 2.0e-3;
        s.params.num_steps = 500; // 1 s horizon
        s.params.implicit = true;
        s.params.tol = 1.0e-10;
        s.params.max_iters = 50;

        let res = s.run().expect("implicit GS run should converge");
        assert_eq!(res.samples.len(), 500, "all steps recorded");
        assert!(res.implicit);
        assert!(res.stable, "implicit run stays bounded");

        // Residual converged below tol on EVERY step (the strong-coupling
        // invariant: residual -> 0).
        assert!(
            res.final_residual < s.params.tol,
            "final residual {:.3e} must be < tol {:.3e}",
            res.final_residual,
            s.params.tol
        );
        assert!(
            res.max_residual < s.params.tol,
            "every step converged: max residual {:.3e} < tol {:.3e}",
            res.max_residual,
            s.params.tol
        );

        // Iteration count is bounded (never hits the cap) and recorded per step.
        assert_eq!(res.iters_per_step.len(), 500);
        assert!(
            res.iters_per_step
                .iter()
                .all(|&k| k >= 1 && k < s.params.max_iters),
            "every step converges within (1, max_iters); counts: first={:?}",
            res.iters_per_step.first()
        );
        assert!(
            res.total_iterations >= res.samples.len(),
            "at least one iter/step"
        );

        // Tracks the monolithic reference tightly at this small macro-step.
        assert!(
            res.final_coupling_error < 1.0e-3,
            "implicit coupling error {:.3e} must be < 1e-3",
            res.final_coupling_error
        );
    }

    #[test]
    fn implicit_max_iters_one_does_not_panic() {
        // With max_iters = 1 the fixed-point loop cannot converge a coupled
        // step, so the coordinator returns NotConverged — which must surface as
        // an in-panel Err, never a panic.
        let mut s = CosimWorkbenchState::default();
        s.params.implicit = true;
        s.params.tol = 1.0e-14;
        s.params.max_iters = 1;
        s.params.num_steps = 10;
        // Either it Errs (NotConverged) or — if one iteration happened to be
        // within tol — it succeeds; both are non-panicking. The benchmark
        // system needs > 1 iteration, so we expect Err here.
        let r = s.run();
        assert!(
            r.is_err(),
            "max_iters = 1 on a genuinely coupled step should surface NotConverged as Err"
        );
    }

    // ---- degenerate-param tests — must return Err, NOT panic ----

    #[test]
    fn zero_macro_step_returns_err() {
        let mut s = CosimWorkbenchState::default();
        s.params.macro_step = 0.0;
        assert!(s.run().is_err(), "H = 0 must return Err, not panic");
    }

    #[test]
    fn negative_macro_step_returns_err() {
        let mut s = CosimWorkbenchState::default();
        s.params.macro_step = -1.0e-3;
        assert!(s.run().is_err(), "H < 0 must return Err, not panic");
    }

    #[test]
    fn zero_steps_returns_err() {
        let mut s = CosimWorkbenchState::default();
        s.params.num_steps = 0;
        assert!(s.run().is_err(), "0 steps must return Err, not panic");
    }

    #[test]
    fn implicit_zero_tol_returns_err() {
        let mut s = CosimWorkbenchState::default();
        s.params.implicit = true;
        s.params.tol = 0.0;
        assert!(s.run().is_err(), "implicit tol = 0 must return Err");
    }

    #[test]
    fn implicit_zero_max_iters_returns_err() {
        let mut s = CosimWorkbenchState::default();
        s.params.implicit = true;
        s.params.max_iters = 0;
        assert!(s.run().is_err(), "implicit max_iters = 0 must return Err");
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_cosim_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn has_named_node(nodes: &[(NodeId, Node)], name: &str) -> bool {
        nodes.iter().any(|(_, n)| n.name() == Some(name))
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_cosim_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_cosim_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_cosim_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_populated_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_cosim_workbench = true;
        let res = app.cosim.run().expect("run should succeed");
        app.cosim.result = Some(res);
        app.cosim.status = "\u{2714} test result".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_implicit_result_without_panic() {
        // Exercise the implicit readout rows (iterations / residual) in the viz.
        let mut app = ValenxApp::default();
        app.show_cosim_workbench = true;
        app.cosim.params.implicit = true;
        app.cosim.params.num_steps = 50;
        let res = app.cosim.run().expect("implicit run should succeed");
        app.cosim.result = Some(res);
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_error_status_without_panic() {
        let mut app = ValenxApp::default();
        app.show_cosim_workbench = true;
        // Trigger an error state (H = 0 is fail-loud in run()).
        app.cosim.params.macro_step = 0.0;
        let result = app.cosim.run();
        app.cosim.status = match result {
            Err(e) => format!("\u{26A0} {e}"),
            Ok(_) => "\u{26A0} simulated error for testing".to_string(),
        };
        app.cosim.result = None;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_labelled_by_named() {
        let mut app = ValenxApp::default();
        app.show_cosim_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        // The numeric DragValues (macro-step, macro-steps, tolerance,
        // max-iters) MUST each carry an accessible name (be labelled_by a
        // caption) so the panel is AI-drivable.
        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 4,
            "expected at least 4 numeric controls (DragValues), got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );

        // Check the specific captions are present as named accessibility nodes.
        for caption in [
            "macro-step H (s)",
            "macro-steps",
            "coupling scheme",
            "implicit coupling",
            "tolerance",
            "max iters / step",
        ] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }

        // The Run button must be a named, invokable node.
        assert!(
            nodes.iter().any(|(_, n)| {
                n.role() == Role::Button && n.name().is_some_and(|s| s.contains("Run"))
            }),
            "the Run button must be a named, invokable node"
        );
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by` its
        // caption; each `labelled_by` target must RESOLVE to a real named
        // caption node, not a dangling id.
        let mut app = ValenxApp::default();
        app.show_cosim_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 4,
            "expected the numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| {
                n.labelled_by()
                    .iter()
                    .any(|id| by_id.get(id).is_some_and(|t| t.name().is_some()))
            }),
            "every DragValue's labelled_by must point at a named caption node"
        );
        for caption in ["macro-step H (s)", "macro-steps"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }

    #[test]
    fn cosim_vs_monolithic_pin_from_ui_state() {
        // Mirror of the unit pin, exercised from the UI-state struct: at a small
        // macro-step the explicit Gauss-Seidel co-sim tracks the monolithic
        // reference within a tight tolerance.
        let mut s = CosimWorkbenchState::default();
        s.params.scheme = CouplingScheme::GaussSeidel;
        s.params.macro_step = 2.0e-4;
        s.params.num_steps = 10_000; // 2 s horizon
        let res = s.run().expect("explicit GS run");
        assert!(
            res.final_coupling_error < 1.0e-3,
            "coupling error {:.3e} should be < 1e-3 vs monolithic",
            res.final_coupling_error
        );
    }

    #[test]
    fn implicit_convergence_pin_from_ui_state() {
        // Mirror of the implicit pin: implicit Gauss-Seidel converges every
        // macro-step below tol with a bounded iteration count.
        let mut s = CosimWorkbenchState::default();
        s.params.implicit = true;
        s.params.num_steps = 100;
        s.params.tol = 1.0e-10;
        s.params.max_iters = 50;
        let res = s.run().expect("implicit GS run");
        assert!(
            res.max_residual < s.params.tol,
            "all steps converged below tol"
        );
        assert!(
            res.iters_per_step.iter().all(|&k| k < s.params.max_iters),
            "no step hit the iteration cap"
        );
    }

    #[test]
    fn degenerate_params_show_error_not_panic() {
        // A zero / negative macro-step (or zero steps) must surface the error
        // in-panel, not panic.
        let mut state = CosimWorkbenchState::default();
        state.params.macro_step = 0.0;
        assert!(state.run().is_err(), "H = 0 must produce Err, not panic");
        state.params.macro_step = 2.0e-3;
        state.params.num_steps = 0;
        assert!(state.run().is_err(), "0 steps must produce Err, not panic");
    }

    #[test]
    fn agent_bridge_cosim_id_resolves_and_sets_flag() {
        // Verify the two mechanisms the agent bridge uses for
        //   `OpenWorkbench { id: "cosim" }`:
        //   1. TabKind::from_id("cosim") -> Some(TabKind::Cosim)
        //      (plus the aliases "co-simulation" / "fmi")
        //   2. set_workbench_flag(app, "cosim", true) -> show_cosim_workbench = true
        use crate::project_tabs::{set_workbench_flag, TabKind};

        // 1. Lookup (canonical + aliases).
        assert_eq!(
            TabKind::from_id("cosim"),
            Some(TabKind::Cosim),
            "\"cosim\" must resolve to TabKind::Cosim"
        );
        assert_eq!(TabKind::from_id("co-simulation"), Some(TabKind::Cosim));
        assert_eq!(TabKind::from_id("fmi"), Some(TabKind::Cosim));
        // Case-insensitive + whitespace-tolerant.
        assert_eq!(TabKind::from_id("  Cosim  "), Some(TabKind::Cosim));

        // 2. Flag toggle.
        let mut app = ValenxApp::default();
        assert!(!app.show_cosim_workbench);
        set_workbench_flag(&mut app, "cosim", true);
        assert!(
            app.show_cosim_workbench,
            "set_workbench_flag(\"cosim\", true) must set the flag"
        );
        set_workbench_flag(&mut app, "cosim", false);
        assert!(!app.show_cosim_workbench);
    }
}
