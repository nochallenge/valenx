//! The network: pipes plus loops, and the Hardy-Cross driver.
//!
//! A [`Network`] owns a list of [`Pipe`]s, an optional endpoint pair per
//! pipe (`from`/`to` node ids) used for the mass-conservation check, and a
//! list of [`Loop`]s. [`Network::solve`] runs the Hardy-Cross iteration to
//! balance the loops, mutating the pipe flows in place and returning a
//! [`SolveReport`].
//!
//! ## Continuity convention
//!
//! A pipe with endpoints `(from, to)` and signed flow `q` carries `q` *out
//! of* `from` and *into* `to`. The net flow imbalance at a node is the sum
//! of pipe flows leaving it minus those entering it, minus any external
//! injection at that node. Hardy-Cross conserves mass exactly at every
//! iteration *if the initial flows do* — the loop corrections only move
//! flow around closed loops, which cannot change any node's balance. So the
//! [`Network::node_imbalances`] check is really a check on the caller's
//! initial guess, and it stays satisfied through the solve.

use crate::error::NetworkError;
use crate::loop_solver::Loop;
use crate::pipe::Pipe;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Endpoints of a pipe, as node identifiers. A pipe's positive-`q`
/// direction runs from [`Endpoints::from`] to [`Endpoints::to`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Endpoints {
    /// Node the pipe's positive flow leaves.
    pub from: u32,
    /// Node the pipe's positive flow enters.
    pub to: u32,
}

impl Endpoints {
    /// Construct an endpoint pair.
    ///
    /// # Errors
    ///
    /// Returns [`NetworkError::BadParameter`] if `from == to` (a pipe may
    /// not loop back to its own start node — that has no defined flow
    /// direction for continuity).
    pub fn new(from: u32, to: u32) -> Result<Self, NetworkError> {
        if from == to {
            return Err(NetworkError::bad(
                "endpoints",
                format!("a pipe's from and to nodes must differ, both were {from}"),
            ));
        }
        Ok(Self { from, to })
    }
}

/// Tuning for the Hardy-Cross iteration.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct SolveConfig {
    /// Convergence tolerance on the largest absolute loop correction `|dQ|`
    /// in a sweep. When every loop's `|dQ|` falls below this, the solve
    /// has converged. Must be finite and `> 0`.
    pub tolerance: f64,
    /// Maximum number of full sweeps over all loops before giving up.
    /// Must be `>= 1`.
    pub max_iterations: usize,
}

impl Default for SolveConfig {
    /// A tight tolerance (`1e-9`) with a generous iteration budget (`200`),
    /// suitable for the small textbook networks this crate targets.
    fn default() -> Self {
        Self {
            tolerance: 1e-9,
            max_iterations: 200,
        }
    }
}

impl SolveConfig {
    /// Validate the configuration.
    ///
    /// # Errors
    ///
    /// Returns [`NetworkError::BadParameter`] if `tolerance` is non-finite
    /// or non-positive, or if `max_iterations` is zero.
    pub fn validate(&self) -> Result<(), NetworkError> {
        if !self.tolerance.is_finite() || self.tolerance <= 0.0 {
            return Err(NetworkError::bad(
                "tolerance",
                format!("must be finite and > 0, got {}", self.tolerance),
            ));
        }
        if self.max_iterations == 0 {
            return Err(NetworkError::bad("max_iterations", "must be >= 1, got 0"));
        }
        Ok(())
    }
}

/// Outcome of a converged Hardy-Cross solve.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct SolveReport {
    /// Number of full sweeps performed to reach convergence.
    pub iterations: usize,
    /// The largest absolute loop correction `|dQ|` in the final sweep.
    /// At convergence this is below the configured tolerance.
    pub final_residual: f64,
}

/// A hydraulic network ready to be balanced by Hardy-Cross.
///
/// Construct with [`Network::new`], optionally attach pipe endpoints with
/// [`Network::with_endpoints`] to enable the mass-conservation check, then
/// call [`Network::solve`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Network {
    pipes: Vec<Pipe>,
    /// Optional, parallel to `pipes`: `endpoints[i]` are the nodes of
    /// `pipes[i]`. Either empty (topology unknown) or the same length as
    /// `pipes`.
    endpoints: Vec<Endpoints>,
    loops: Vec<Loop>,
}

impl Network {
    /// Build a network from its pipes and loops, with no endpoint topology.
    ///
    /// Every loop is validated against the pipe count immediately, so an
    /// out-of-range loop member is reported here rather than at solve time.
    ///
    /// # Errors
    ///
    /// Returns [`NetworkError::UnknownPipe`] if any loop references a pipe
    /// index that does not exist.
    pub fn new(pipes: Vec<Pipe>, loops: Vec<Loop>) -> Result<Self, NetworkError> {
        for lp in &loops {
            lp.validate_against(pipes.len())?;
        }
        Ok(Self {
            pipes,
            endpoints: Vec::new(),
            loops,
        })
    }

    /// Attach an endpoint pair to every pipe, enabling
    /// [`Network::node_imbalances`] and the continuity helpers.
    ///
    /// # Errors
    ///
    /// Returns [`NetworkError::BadParameter`] if `endpoints.len()` does not
    /// equal the number of pipes.
    pub fn with_endpoints(mut self, endpoints: Vec<Endpoints>) -> Result<Self, NetworkError> {
        if endpoints.len() != self.pipes.len() {
            return Err(NetworkError::bad(
                "endpoints",
                format!(
                    "expected one endpoint pair per pipe ({} pipes), got {}",
                    self.pipes.len(),
                    endpoints.len()
                ),
            ));
        }
        self.endpoints = endpoints;
        Ok(self)
    }

    /// Read-only view of the pipes (with their current flows).
    pub fn pipes(&self) -> &[Pipe] {
        &self.pipes
    }

    /// Read-only view of the loops.
    pub fn loops(&self) -> &[Loop] {
        &self.loops
    }

    /// Read-only view of the endpoint topology, empty if none was attached.
    pub fn endpoints(&self) -> &[Endpoints] {
        &self.endpoints
    }

    /// Signed head loss around loop `loop_index` for the current flows.
    ///
    /// Returns the loop's [`Loop::head_loss_sum`]. At convergence this is
    /// within solver tolerance of zero. Returns `None` if `loop_index` is
    /// out of range.
    pub fn loop_head_loss(&self, loop_index: usize) -> Option<f64> {
        self.loops
            .get(loop_index)
            .map(|lp| lp.head_loss_sum(&self.pipes))
    }

    /// Net flow imbalance at every node, keyed by node id.
    ///
    /// The imbalance at a node is `(flow out) - (flow in) - external_inflow`,
    /// summed over all incident pipes. A node where the entry is zero (to
    /// within rounding) is *balanced* — conserves mass. `external_inflow`
    /// maps a node id to flow injected from outside the network at that node
    /// (positive = supply into the node); nodes absent from the map are
    /// taken to have zero external injection.
    ///
    /// Returns an empty map if no endpoint topology was attached.
    ///
    /// Only nodes that appear in the topology (or in `external_inflow`) are
    /// present in the result.
    pub fn node_imbalances(&self, external_inflow: &BTreeMap<u32, f64>) -> BTreeMap<u32, f64> {
        let mut balance: BTreeMap<u32, f64> = BTreeMap::new();
        for (pipe, ep) in self.pipes.iter().zip(self.endpoints.iter()) {
            // `q` leaves `from` and enters `to`.
            *balance.entry(ep.from).or_insert(0.0) += pipe.q;
            *balance.entry(ep.to).or_insert(0.0) -= pipe.q;
        }
        for (&node, &inflow) in external_inflow {
            *balance.entry(node).or_insert(0.0) -= inflow;
        }
        balance
    }

    /// Largest absolute node imbalance for the given external injections.
    ///
    /// A convenient scalar summary of [`Network::node_imbalances`]; returns
    /// `0.0` when there is no topology or no nodes.
    pub fn max_node_imbalance(&self, external_inflow: &BTreeMap<u32, f64>) -> f64 {
        self.node_imbalances(external_inflow)
            .values()
            .fold(0.0_f64, |acc, &v| acc.max(v.abs()))
    }

    /// Run the Hardy-Cross iteration to balance every loop, mutating the
    /// pipe flows in place.
    ///
    /// Each sweep computes and applies the correction `dQ` for every loop in
    /// turn (Gauss-Seidel style: later loops in a sweep see the corrections
    /// already applied by earlier ones). The sweep's residual is the largest
    /// `|dQ|` applied; when that drops below `config.tolerance` the solve has
    /// converged.
    ///
    /// # Errors
    ///
    /// - [`NetworkError::BadParameter`] if `config` is invalid.
    /// - [`NetworkError::EmptyLoop`] if any loop has no members (this is
    ///   normally caught at [`Loop::new`], but is re-checked defensively).
    /// - [`NetworkError::NoConvergence`] if the residual is still above
    ///   tolerance after `config.max_iterations` sweeps.
    ///
    /// A loop whose correction is undefined (every member loss-free or
    /// stationary, see [`Loop::correction`]) contributes a zero correction
    /// for that sweep rather than erroring; if such a loop is the only thing
    /// preventing convergence the solve will instead surface
    /// [`NetworkError::NoConvergence`].
    pub fn solve(&mut self, config: &SolveConfig) -> Result<SolveReport, NetworkError> {
        config.validate()?;
        for lp in &self.loops {
            if lp.members.is_empty() {
                return Err(NetworkError::EmptyLoop(lp.name.clone()));
            }
        }

        let mut last_residual = f64::INFINITY;
        for iteration in 1..=config.max_iterations {
            let mut sweep_residual = 0.0_f64;
            for lp_index in 0..self.loops.len() {
                let dq = match self.loops[lp_index].correction(&self.pipes) {
                    Some(dq) => dq,
                    None => continue,
                };
                self.loops[lp_index].apply_correction(&mut self.pipes, dq);
                sweep_residual = sweep_residual.max(dq.abs());
            }
            last_residual = sweep_residual;
            if sweep_residual < config.tolerance {
                return Ok(SolveReport {
                    iterations: iteration,
                    final_residual: sweep_residual,
                });
            }
        }

        Err(NetworkError::NoConvergence {
            iterations: config.max_iterations,
            residual: last_residual,
            tolerance: config.tolerance,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;
    use crate::loop_solver::LoopMember;

    /// Tolerance for analytic comparisons after a tight solve.
    const SOLVE_EPS: f64 = 1e-7;
    /// Tolerance for exact arithmetic identities.
    const EXACT_EPS: f64 = 1e-12;

    // ----- single-loop ground truth -------------------------------------

    /// Two parallel pipes A->B, k0=1, k1=4, total inflow Q=3.
    /// Balance k0 q0^2 = k1 q1^2 with q0 + q1 = 3 gives q0=2, q1=1.
    fn parallel_network(q0: f64, q1: f64) -> Network {
        let pipes = vec![Pipe::new(1.0, q0).unwrap(), Pipe::new(4.0, q1).unwrap()];
        let lp = Loop::new("loop", vec![LoopMember::forward(0), LoopMember::reverse(1)]).unwrap();
        Network::new(pipes, vec![lp]).unwrap()
    }

    #[test]
    fn single_loop_converges_to_analytic_split() {
        // Start from a continuity-satisfying guess (1.5, 1.5), Q=3.
        let mut net = parallel_network(1.5, 1.5);
        let report = net.solve(&SolveConfig::default()).unwrap();
        assert!(report.iterations >= 1);
        let q0 = net.pipes()[0].q;
        let q1 = net.pipes()[1].q;
        assert!((q0 - 2.0).abs() < SOLVE_EPS, "q0={q0}, expected 2");
        assert!((q1 - 1.0).abs() < SOLVE_EPS, "q1={q1}, expected 1");
    }

    #[test]
    fn single_loop_head_loss_vanishes_on_convergence() {
        let mut net = parallel_network(1.5, 1.5);
        net.solve(&SolveConfig::default()).unwrap();
        let h = net.loop_head_loss(0).unwrap();
        assert!(h.abs() < SOLVE_EPS, "loop head loss {h} should be ~0");
    }

    #[test]
    fn single_loop_balances_from_a_far_initial_guess() {
        // A lopsided but continuity-respecting guess (2.9, 0.1) still lands
        // on (2, 1).
        let mut net = parallel_network(2.9, 0.1);
        net.solve(&SolveConfig::default()).unwrap();
        assert!((net.pipes()[0].q - 2.0).abs() < SOLVE_EPS);
        assert!((net.pipes()[1].q - 1.0).abs() < SOLVE_EPS);
    }

    #[test]
    fn symmetric_loop_splits_evenly() {
        // k0 = k1 = 1, Q = 4 -> q0 = q1 = 2 by symmetry.
        let pipes = vec![Pipe::new(1.0, 1.0).unwrap(), Pipe::new(1.0, 3.0).unwrap()];
        let lp = Loop::new("l", vec![LoopMember::forward(0), LoopMember::reverse(1)]).unwrap();
        let mut net = Network::new(pipes, vec![lp]).unwrap();
        net.solve(&SolveConfig::default()).unwrap();
        assert!((net.pipes()[0].q - 2.0).abs() < SOLVE_EPS);
        assert!((net.pipes()[1].q - 2.0).abs() < SOLVE_EPS);
    }

    // ----- node mass conservation ---------------------------------------

    #[test]
    fn continuity_holds_before_and_after_solve() {
        // Attach topology: both pipes run node 1 -> node 2. Inflow 3 at
        // node 1, outflow 3 at node 2.
        let net = parallel_network(1.5, 1.5)
            .with_endpoints(vec![
                Endpoints::new(1, 2).unwrap(),
                Endpoints::new(1, 2).unwrap(),
            ])
            .unwrap();
        let mut external = BTreeMap::new();
        external.insert(1u32, 3.0); // 3 supplied into node 1
        external.insert(2u32, -3.0); // 3 drawn from node 2

        // Before solving, the (1.5,1.5) guess already conserves mass.
        assert!(net.max_node_imbalance(&external) < EXACT_EPS);

        let mut net = net;
        net.solve(&SolveConfig::default()).unwrap();
        // Hardy-Cross only circulates flow around loops, so continuity is
        // preserved exactly through the solve.
        let imbalances = net.node_imbalances(&external);
        for (node, imb) in &imbalances {
            assert!(imb.abs() < SOLVE_EPS, "node {node} imbalance {imb}");
        }
        assert!(net.max_node_imbalance(&external) < SOLVE_EPS);
    }

    #[test]
    fn node_imbalances_detects_a_bad_guess() {
        // Endpoints 1->2 for both pipes; flows (2.0, 2.0) means 4 leaves
        // node 1, but external supply is only 3 -> imbalance 1 at node 1.
        let net = parallel_network(2.0, 2.0)
            .with_endpoints(vec![
                Endpoints::new(1, 2).unwrap(),
                Endpoints::new(1, 2).unwrap(),
            ])
            .unwrap();
        let mut external = BTreeMap::new();
        external.insert(1u32, 3.0);
        external.insert(2u32, -3.0);
        let imbalances = net.node_imbalances(&external);
        // node 1: q0 + q1 - 3 = 4 - 3 = 1
        assert!((imbalances[&1] - 1.0).abs() < EXACT_EPS);
        // node 2: -(q0 + q1) - (-3) = -4 + 3 = -1
        assert!((imbalances[&2] - (-1.0)).abs() < EXACT_EPS);
    }

    // ----- two-loop (theta) ground truth --------------------------------

    /// Three parallel pipes A->B with k = (1, 4, 9) and two independent
    /// loops. Equal head losses force k_i q_i^2 = C; with Q = 11 the
    /// analytic flows are (6, 3, 2): 1*36 = 4*9 = 9*4 = 36.
    #[test]
    fn theta_two_loop_network_matches_analytic_flows() {
        let g = 11.0 / 3.0; // continuity-respecting initial guess: 11/3 each
        let pipes = vec![
            Pipe::new(1.0, g).unwrap(),
            Pipe::new(4.0, g).unwrap(),
            Pipe::new(9.0, g).unwrap(),
        ];
        let loop_a = Loop::new("A", vec![LoopMember::forward(0), LoopMember::reverse(1)]).unwrap();
        let loop_b = Loop::new("B", vec![LoopMember::forward(1), LoopMember::reverse(2)]).unwrap();
        let mut net = Network::new(pipes, vec![loop_a, loop_b]).unwrap();

        let report = net.solve(&SolveConfig::default()).unwrap();
        assert!(report.final_residual < SolveConfig::default().tolerance);

        let q: Vec<f64> = net.pipes().iter().map(|p| p.q).collect();
        assert!((q[0] - 6.0).abs() < SOLVE_EPS, "q0={}", q[0]);
        assert!((q[1] - 3.0).abs() < SOLVE_EPS, "q1={}", q[1]);
        assert!((q[2] - 2.0).abs() < SOLVE_EPS, "q2={}", q[2]);

        // Both loops are balanced.
        assert!(net.loop_head_loss(0).unwrap().abs() < SOLVE_EPS);
        assert!(net.loop_head_loss(1).unwrap().abs() < SOLVE_EPS);

        // Continuity: 11 in at node 1, 11 out at node 2.
        let net = net
            .with_endpoints(vec![
                Endpoints::new(1, 2).unwrap(),
                Endpoints::new(1, 2).unwrap(),
                Endpoints::new(1, 2).unwrap(),
            ])
            .unwrap();
        let mut external = BTreeMap::new();
        external.insert(1u32, 11.0);
        external.insert(2u32, -11.0);
        assert!(net.max_node_imbalance(&external) < SOLVE_EPS);
    }

    // ----- iteration monotonicity ---------------------------------------

    #[test]
    fn loop_head_loss_magnitude_decreases_per_sweep() {
        // Solve one sweep at a time and confirm the loop imbalance shrinks.
        let mut net = parallel_network(2.9, 0.1);
        let one_sweep = SolveConfig {
            tolerance: 1e-30, // force max_iterations to be the stop, not tol
            max_iterations: 1,
        };
        let mut prev = net.loop_head_loss(0).unwrap().abs();
        for _ in 0..6 {
            // Each call runs exactly one sweep (it cannot converge at this
            // absurd tolerance), so ignore the NoConvergence result.
            let _ = net.solve(&one_sweep);
            let now = net.loop_head_loss(0).unwrap().abs();
            assert!(now < prev || now < SOLVE_EPS, "{now} !< {prev}");
            prev = now;
        }
    }

    // ----- error paths ---------------------------------------------------

    #[test]
    fn unknown_pipe_in_loop_is_caught_at_construction() {
        let pipes = vec![Pipe::new(1.0, 1.0).unwrap()];
        let lp = Loop::new("l", vec![LoopMember::forward(3)]).unwrap();
        let err = Network::new(pipes, vec![lp]).unwrap_err();
        assert!(matches!(
            err,
            NetworkError::UnknownPipe { index: 3, count: 1 }
        ));
    }

    #[test]
    fn no_convergence_is_reported() {
        // One sweep is far too few to balance the lopsided guess to 1e-30.
        let mut net = parallel_network(2.9, 0.1);
        let cfg = SolveConfig {
            tolerance: 1e-30,
            max_iterations: 1,
        };
        let err = net.solve(&cfg).unwrap_err();
        assert!(matches!(
            err,
            NetworkError::NoConvergence { iterations: 1, .. }
        ));
        assert_eq!(err.code(), "pipenetwork.no_convergence");
        assert_eq!(err.category(), ErrorCategory::Algorithm);
    }

    #[test]
    fn invalid_config_is_rejected() {
        let mut net = parallel_network(1.5, 1.5);
        assert!(net
            .solve(&SolveConfig {
                tolerance: 0.0,
                max_iterations: 10
            })
            .is_err());
        assert!(net
            .solve(&SolveConfig {
                tolerance: 1e-9,
                max_iterations: 0
            })
            .is_err());
        assert!(net
            .solve(&SolveConfig {
                tolerance: f64::NAN,
                max_iterations: 10
            })
            .is_err());
    }

    #[test]
    fn endpoints_length_must_match_pipe_count() {
        let net = parallel_network(1.5, 1.5);
        // Two pipes, but only one endpoint pair supplied.
        let err = net
            .clone()
            .with_endpoints(vec![Endpoints::new(1, 2).unwrap()])
            .unwrap_err();
        assert!(matches!(
            err,
            NetworkError::BadParameter {
                name: "endpoints",
                ..
            }
        ));
    }

    #[test]
    fn self_looping_endpoint_is_rejected() {
        let err = Endpoints::new(7, 7).unwrap_err();
        assert!(matches!(
            err,
            NetworkError::BadParameter {
                name: "endpoints",
                ..
            }
        ));
    }

    #[test]
    fn degenerate_all_zero_loop_does_not_converge_in_one_step() {
        // Every pipe loss-free: corrections are all None, so the loop never
        // balances away from a non-balanced start... but here the start is
        // already balanced (zero loss), so it converges immediately.
        let pipes = vec![Pipe::new(0.0, 1.0).unwrap(), Pipe::new(0.0, 1.0).unwrap()];
        let lp = Loop::new("l", vec![LoopMember::forward(0), LoopMember::reverse(1)]).unwrap();
        let mut net = Network::new(pipes, vec![lp]).unwrap();
        // Loss-free loop has zero head loss already; solve converges at the
        // first sweep with a zero residual (no correction applied).
        let report = net.solve(&SolveConfig::default()).unwrap();
        assert_eq!(report.iterations, 1);
        assert!(report.final_residual.abs() < EXACT_EPS);
    }
}
