//! Native in-process co-simulation master.
//!
//! This is the **house-style native path** (per `AGENTS.md`: native-first,
//! fail-loud). A co-simulation couples several independently-integrated
//! [`Subsystem`]s that exchange scalar signals at discrete *macro-steps*.
//! Each subsystem advances itself over one macro-step of size `dt` given
//! its current input vector, and returns its current output vector; the
//! [`CoSimMaster`] routes outputs to inputs through a [`CouplingGraph`]
//! using either a [`Scheme::Jacobi`] or [`Scheme::GaussSeidel`] update.
//!
//! The trait is deliberately tiny so that an HLA / AFSIM federate bridge,
//! or a binary FMU (see [`crate::fmi`]), can implement it later without
//! changing the master.

use crate::error::{FmiError, Result};

/// One in-process coupling unit.
///
/// A subsystem owns its internal state and integrator. `step` advances
/// that internal state from time `t` to `t + dt` using the supplied
/// `inputs` (held constant across the macro-step, as is standard for
/// loosely-coupled co-simulation), and returns the subsystem's output
/// signals **sampled at the end of the step**.
///
/// Contract:
/// - `inputs.len()` equals [`Subsystem::n_inputs`].
/// - the returned vector's length equals [`Subsystem::n_outputs`].
///
/// `step` takes `&mut self` because advancing typically mutates internal
/// state; a master may call it more than once per macro-step when an
/// iterative (implicit) coupling scheme is layered on top, so it must be
/// safe to re-enter — see [`Subsystem::set_state`] / [`Subsystem::state`]
/// if your subsystem needs to roll back between iterations.
pub trait Subsystem {
    /// Number of scalar input channels this subsystem reads.
    fn n_inputs(&self) -> usize;

    /// Number of scalar output channels this subsystem writes.
    fn n_outputs(&self) -> usize;

    /// Advance internal state over `[t, t + dt]` with `inputs` held
    /// constant, returning the output vector sampled at `t + dt`.
    fn step(&mut self, t: f64, dt: f64, inputs: &[f64]) -> Vec<f64>;

    /// Snapshot the internal continuous state (for rollback during an
    /// iterative coupling sweep). Default: no rollback support.
    fn state(&self) -> Vec<f64> {
        Vec::new()
    }

    /// Restore a snapshot previously returned by [`Subsystem::state`].
    /// Default: no-op (matches the default empty `state`).
    fn set_state(&mut self, _state: &[f64]) {}
}

/// Which subsystem output feeds which subsystem input.
///
/// Identifies a single scalar channel by `(subsystem index, port index)`
/// on each end. The output of `(from_subsystem, from_output)` is copied
/// into the input slot `(to_subsystem, to_input)` of the consumer before
/// it steps.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Coupling {
    /// Index of the producing subsystem.
    pub from_subsystem: usize,
    /// Output port index on the producing subsystem.
    pub from_output: usize,
    /// Index of the consuming subsystem.
    pub to_subsystem: usize,
    /// Input port index on the consuming subsystem.
    pub to_input: usize,
}

impl Coupling {
    /// Construct a coupling edge `(from_subsystem.from_output) ->
    /// (to_subsystem.to_input)`.
    pub fn new(
        from_subsystem: usize,
        from_output: usize,
        to_subsystem: usize,
        to_input: usize,
    ) -> Self {
        Self {
            from_subsystem,
            from_output,
            to_subsystem,
            to_input,
        }
    }
}

/// The full set of coupling edges over a list of subsystems.
#[derive(Clone, Debug, Default)]
pub struct CouplingGraph {
    edges: Vec<Coupling>,
}

impl CouplingGraph {
    /// An empty graph (no signals exchanged).
    pub fn new() -> Self {
        Self { edges: Vec::new() }
    }

    /// Build a graph directly from a list of edges.
    pub fn from_edges(edges: Vec<Coupling>) -> Self {
        Self { edges }
    }

    /// Append one coupling edge.
    pub fn add(&mut self, edge: Coupling) -> &mut Self {
        self.edges.push(edge);
        self
    }

    /// The edges, in declaration order.
    pub fn edges(&self) -> &[Coupling] {
        &self.edges
    }

    /// Validate every edge against the port counts of `subsystems`.
    ///
    /// Fail-loud: returns the FIRST problem found rather than silently
    /// dropping a bad edge. Per edge: (1) both subsystem indices are in
    /// range, and (2) the output port exists on the producer and the input
    /// port exists on the consumer. Across edges: (3) no two edges drive
    /// the same `(consumer, input)` slot.
    pub fn validate(&self, subsystems: &[Box<dyn Subsystem>]) -> Result<()> {
        let n = subsystems.len();
        let mut driven: Vec<(usize, usize)> = Vec::new();
        for (edge_index, e) in self.edges.iter().enumerate() {
            // 1. subsystem indices in range.
            if e.from_subsystem >= n {
                return Err(FmiError::SubsystemIndexOutOfRange {
                    edge_index,
                    bad_index: e.from_subsystem,
                    n_subsystems: n,
                });
            }
            if e.to_subsystem >= n {
                return Err(FmiError::SubsystemIndexOutOfRange {
                    edge_index,
                    bad_index: e.to_subsystem,
                    n_subsystems: n,
                });
            }

            // 2. port indices in range.
            let n_out = subsystems[e.from_subsystem].n_outputs();
            if e.from_output >= n_out {
                return Err(FmiError::PortIndexOutOfRange {
                    edge_index,
                    side: "output",
                    subsystem: e.from_subsystem,
                    bad_port: e.from_output,
                    n_ports: n_out,
                });
            }
            let n_in = subsystems[e.to_subsystem].n_inputs();
            if e.to_input >= n_in {
                return Err(FmiError::PortIndexOutOfRange {
                    edge_index,
                    side: "input",
                    subsystem: e.to_subsystem,
                    bad_port: e.to_input,
                    n_ports: n_in,
                });
            }

            // 3. each input fed by at most one edge.
            let slot = (e.to_subsystem, e.to_input);
            if driven.contains(&slot) {
                return Err(FmiError::DuplicateInputSource {
                    subsystem: e.to_subsystem,
                    port: e.to_input,
                });
            }
            driven.push(slot);
        }
        Ok(())
    }
}

/// Macro-step coupling scheme.
///
/// The two classic explicit (non-iterative) co-simulation update orders:
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scheme {
    /// **Jacobi (parallel).** Every subsystem steps using the input
    /// values assembled from the *previous* macro-step's outputs; all
    /// subsystems see a consistent, one-macro-step-old picture. Outputs
    /// are committed only after every subsystem has stepped. Parallelizable
    /// but introduces a one-step coupling lag.
    Jacobi,
    /// **Gauss-Seidel (sequential).** Subsystems step in index order, and
    /// each one sees the *latest* outputs produced earlier within the same
    /// macro-step. Tighter coupling (less lag) at the cost of a fixed
    /// sequential order.
    GaussSeidel,
}

/// The native co-simulation master: subsystems + coupling graph + scheme.
///
/// Construct with [`CoSimMaster::new`] (which validates the graph),
/// then drive the coupled system one macro-step at a time with
/// [`CoSimMaster::advance`].
pub struct CoSimMaster {
    subsystems: Vec<Box<dyn Subsystem>>,
    graph: CouplingGraph,
    scheme: Scheme,
    /// Latest committed output vector per subsystem.
    outputs: Vec<Vec<f64>>,
    /// Current simulation time (advances by `dt` per macro-step).
    t: f64,
}

// Manual `Debug`: `Box<dyn Subsystem>` is not `Debug`, so summarise the
// master by its shape (subsystem count, edges, scheme, time) instead of
// requiring every `Subsystem` impl to also be `Debug`. This lets call
// sites use `Result::unwrap_err` etc. in tests.
impl std::fmt::Debug for CoSimMaster {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CoSimMaster")
            .field("n_subsystems", &self.subsystems.len())
            .field("n_edges", &self.graph.edges().len())
            .field("scheme", &self.scheme)
            .field("t", &self.t)
            .finish()
    }
}

impl CoSimMaster {
    /// Assemble a master from subsystems, a coupling graph, and a scheme.
    ///
    /// Validates the graph up front (fail-loud): a dangling subsystem or
    /// port reference, or a doubly-driven input, returns `Err` here rather
    /// than panicking later inside [`CoSimMaster::advance`].
    ///
    /// The output channels are then **primed with the subsystems' true
    /// initial outputs** by calling [`Subsystem::step`] once with
    /// `dt = 0.0` (a zero-length step must not advance state — it only
    /// samples the current outputs). This matters for accuracy: the very
    /// first macro-step must route each subsystem's *actual* initial
    /// output (e.g. a mass's initial position), not a placeholder zero,
    /// otherwise both schemes carry a spurious first-step transient that
    /// pollutes a benchmark comparison. A subsystem with state-dependent
    /// outputs that cannot honour a `dt == 0` sample may override
    /// [`Subsystem::state`]/[`Subsystem::set_state`] to roll back instead.
    pub fn new(
        mut subsystems: Vec<Box<dyn Subsystem>>,
        graph: CouplingGraph,
        scheme: Scheme,
    ) -> Result<Self> {
        graph.validate(&subsystems)?;
        // Prime each subsystem's output channel. There are no committed
        // outputs yet, so feed zero inputs for the dt=0 sample — a
        // well-behaved subsystem ignores inputs when dt is zero (it does
        // not integrate), reporting its initial outputs.
        let outputs: Vec<Vec<f64>> = subsystems
            .iter_mut()
            .map(|s| {
                let zero_inputs = vec![0.0; s.n_inputs()];
                let out = s.step(0.0, 0.0, &zero_inputs);
                debug_assert_eq!(out.len(), s.n_outputs(), "subsystem output arity mismatch");
                out
            })
            .collect();
        Ok(Self {
            subsystems,
            graph,
            scheme,
            outputs,
            t: 0.0,
        })
    }

    /// Current simulation time.
    pub fn time(&self) -> f64 {
        self.t
    }

    /// The latest committed output vector of subsystem `i`, or `None` if
    /// `i` is out of range.
    pub fn outputs_of(&self, i: usize) -> Option<&[f64]> {
        self.outputs.get(i).map(|v| v.as_slice())
    }

    /// Number of registered subsystems.
    pub fn len(&self) -> usize {
        self.subsystems.len()
    }

    /// Whether no subsystems are registered.
    pub fn is_empty(&self) -> bool {
        self.subsystems.is_empty()
    }

    /// Assemble the input vector for subsystem `i` by gathering the
    /// driving outputs from `source_outputs`. Unconnected inputs are 0.0.
    fn assemble_inputs(&self, i: usize, source_outputs: &[Vec<f64>]) -> Vec<f64> {
        let n_in = self.subsystems[i].n_inputs();
        let mut inputs = vec![0.0; n_in];
        for e in self.graph.edges() {
            if e.to_subsystem == i {
                // Bounds were validated in `new`; indexing is safe.
                inputs[e.to_input] = source_outputs[e.from_subsystem][e.from_output];
            }
        }
        inputs
    }

    /// Run exactly one macro-step of size `dt`, advancing simulation time
    /// from `t` to `t + dt` and committing each subsystem's new outputs.
    ///
    /// - [`Scheme::Jacobi`]: inputs for every subsystem are gathered from
    ///   the outputs committed at the END of the previous macro-step, so
    ///   all subsystems step against the same snapshot.
    /// - [`Scheme::GaussSeidel`]: subsystems step in index order, and each
    ///   gathers inputs from a working copy that already contains the
    ///   outputs of the lower-indexed subsystems stepped earlier in THIS
    ///   macro-step.
    pub fn advance(&mut self, dt: f64) {
        let t = self.t;
        match self.scheme {
            Scheme::Jacobi => {
                // Freeze the previous snapshot; every subsystem reads it.
                let snapshot = self.outputs.clone();
                let mut new_outputs: Vec<Vec<f64>> = Vec::with_capacity(self.subsystems.len());
                for i in 0..self.subsystems.len() {
                    let inputs = self.assemble_inputs(i, &snapshot);
                    let out = self.subsystems[i].step(t, dt, &inputs);
                    new_outputs.push(out);
                }
                self.outputs = new_outputs;
            }
            Scheme::GaussSeidel => {
                // Working copy is updated in place as we go, so a later
                // subsystem sees earlier subsystems' fresh outputs.
                let mut working = self.outputs.clone();
                for i in 0..self.subsystems.len() {
                    let inputs = self.assemble_inputs(i, &working);
                    let out = self.subsystems[i].step(t, dt, &inputs);
                    working[i] = out;
                }
                self.outputs = working;
            }
        }
        self.t += dt;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A trivial passthrough: y = u (one in, one out), with no state.
    struct Passthrough;
    impl Subsystem for Passthrough {
        fn n_inputs(&self) -> usize {
            1
        }
        fn n_outputs(&self) -> usize {
            1
        }
        fn step(&mut self, _t: f64, _dt: f64, inputs: &[f64]) -> Vec<f64> {
            vec![inputs[0]]
        }
    }

    /// A constant source: zero inputs, one output that holds `value`.
    struct Source {
        value: f64,
    }
    impl Subsystem for Source {
        fn n_inputs(&self) -> usize {
            0
        }
        fn n_outputs(&self) -> usize {
            1
        }
        fn step(&mut self, _t: f64, _dt: f64, _inputs: &[f64]) -> Vec<f64> {
            vec![self.value]
        }
    }

    #[test]
    fn empty_graph_validates_and_is_inert() {
        let subs: Vec<Box<dyn Subsystem>> = vec![Box::new(Passthrough), Box::new(Passthrough)];
        let mut m = CoSimMaster::new(subs, CouplingGraph::new(), Scheme::GaussSeidel).unwrap();
        m.advance(0.1);
        // Nothing connected → all outputs read zero inputs → 0.0.
        assert_eq!(m.outputs_of(0), Some(&[0.0][..]));
        assert!((m.time() - 0.1).abs() < 1e-12);
    }

    #[test]
    fn gauss_seidel_propagates_within_one_macro_step() {
        // source(3) -> A.in ; A.out -> B.in. Gauss-Seidel: in ONE step
        // B should already see A's fresh output (= 3), not a lagged 0.
        let subs: Vec<Box<dyn Subsystem>> = vec![
            Box::new(Source { value: 3.0 }),
            Box::new(Passthrough),
            Box::new(Passthrough),
        ];
        let graph = CouplingGraph::from_edges(vec![
            Coupling::new(0, 0, 1, 0), // source -> A
            Coupling::new(1, 0, 2, 0), // A -> B
        ]);
        let mut m = CoSimMaster::new(subs, graph, Scheme::GaussSeidel).unwrap();
        m.advance(1.0);
        assert_eq!(m.outputs_of(2), Some(&[3.0][..]));
    }

    #[test]
    fn jacobi_lags_one_macro_step_per_link() {
        // source -> A -> B chain under Jacobi. The master primes outputs at
        // construction with a dt=0 sample, so the constant source already
        // reports 3 and the two passthroughs report 0. Thereafter each
        // Jacobi macro-step reads the PREVIOUS snapshot, so the signal
        // advances one link per step:
        //   primed:  source=3, A=0, B=0
        //   step 1:  A sees source=3 -> 3 ; B sees A(old)=0 -> 0
        //   step 2:  B sees A=3 -> 3
        let subs: Vec<Box<dyn Subsystem>> = vec![
            Box::new(Source { value: 3.0 }),
            Box::new(Passthrough),
            Box::new(Passthrough),
        ];
        let graph =
            CouplingGraph::from_edges(vec![Coupling::new(0, 0, 1, 0), Coupling::new(1, 0, 2, 0)]);
        let mut m = CoSimMaster::new(subs, graph, Scheme::Jacobi).unwrap();
        // Primed state: source reports 3, passthroughs report 0.
        assert_eq!(m.outputs_of(0), Some(&[3.0][..]), "primed: source = 3");
        assert_eq!(m.outputs_of(1), Some(&[0.0][..]), "primed: A = 0");
        m.advance(1.0);
        assert_eq!(m.outputs_of(1), Some(&[3.0][..]), "step 1: A now driven");
        assert_eq!(m.outputs_of(2), Some(&[0.0][..]), "step 1: B still lagging");
        m.advance(1.0);
        assert_eq!(
            m.outputs_of(2),
            Some(&[3.0][..]),
            "step 2: signal reaches B"
        );
    }

    #[test]
    fn validate_rejects_dangling_subsystem_index() {
        let subs: Vec<Box<dyn Subsystem>> = vec![Box::new(Passthrough)];
        // Edge references subsystem 5, which doesn't exist.
        let graph = CouplingGraph::from_edges(vec![Coupling::new(0, 0, 5, 0)]);
        let err = CoSimMaster::new(subs, graph, Scheme::Jacobi).unwrap_err();
        assert!(matches!(
            err,
            FmiError::SubsystemIndexOutOfRange { bad_index: 5, .. }
        ));
    }

    #[test]
    fn validate_rejects_out_of_range_port() {
        let subs: Vec<Box<dyn Subsystem>> = vec![Box::new(Passthrough), Box::new(Passthrough)];
        // Passthrough has only output port 0; reference port 2.
        let graph = CouplingGraph::from_edges(vec![Coupling::new(0, 2, 1, 0)]);
        let err = CoSimMaster::new(subs, graph, Scheme::Jacobi).unwrap_err();
        assert!(matches!(
            err,
            FmiError::PortIndexOutOfRange {
                side: "output",
                bad_port: 2,
                ..
            }
        ));
    }

    #[test]
    fn validate_rejects_double_driven_input() {
        let subs: Vec<Box<dyn Subsystem>> = vec![
            Box::new(Source { value: 1.0 }),
            Box::new(Source { value: 2.0 }),
            Box::new(Passthrough),
        ];
        // Both sources drive (subsystem 2, input 0).
        let graph =
            CouplingGraph::from_edges(vec![Coupling::new(0, 0, 2, 0), Coupling::new(1, 0, 2, 0)]);
        let err = CoSimMaster::new(subs, graph, Scheme::GaussSeidel).unwrap_err();
        assert!(matches!(
            err,
            FmiError::DuplicateInputSource {
                subsystem: 2,
                port: 0
            }
        ));
    }
}
