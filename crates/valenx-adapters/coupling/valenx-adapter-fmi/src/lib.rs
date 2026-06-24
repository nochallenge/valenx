//! # valenx-adapter-fmi
//!
//! A **native co-simulation master** for Valenx, plus the connective layer
//! to bring outside simulation interfaces in as coupling units.
//!
//! Per `AGENTS.md` (native-first, fail-loud, benchmark-pinned tests), the
//! primary path is fully in-house: a [`cosim::Subsystem`] is an in-process
//! Rust coupling unit, and [`cosim::CoSimMaster`] advances a graph of them
//! one macro-step at a time under a [`cosim::Scheme`] (Jacobi or
//! Gauss-Seidel). Nothing here requires an external binary to run.
//!
//! Four connective pieces sit on top of that core:
//!
//! * [`cosim`] — the native co-simulation master (the house-style path; an
//!   HLA / AFSIM federate, or the binary FMU below, can implement
//!   `Subsystem` too).
//! * [`federation`] — a **HELICS-style co-simulation federation** with
//!   distributed (dependency-aware) **time coordination**: a
//!   [`federation::Broker`] grants each [`federation::Federate`] the largest
//!   time consistent with every federate's lookahead, generalizing the
//!   fixed-step master into coordinated multi-rate advancement, with named
//!   pub/sub value exchange and timed message endpoints. Pure-Rust, in-process,
//!   reimplemented from the published HELICS algorithm (BSD-3 — credited in the
//!   module). A co-sim FMU rides it via [`federation::SubsystemFederate`].
//! * [`fmi`] — **FMI 2.0 / 3.0 co-simulation import**: a hand-rolled,
//!   dependency-free parser for an FMU's `modelDescription.xml` that pulls
//!   out the model name and the scalar interface variables (name, value
//!   reference, causality). This is co-simulation *import only* — not a
//!   model-exchange importer.
//! * [`dis`] — an **IEEE 1278.1 DIS Entity State PDU** codec with the
//!   bit-exact, big-endian (network-order) standard byte layout. Entity
//!   State PDU only — not the full PDU family.
//!
//! Loading an actual compiled binary FMU (`.so` / `.dll` / `.dylib`) is the
//! optional `binary` module, gated behind the off-by-default `binary-fmu`
//! cargo feature; the default build neither pulls in `libloading` nor
//! requires a real FMU. (That module — and its intra-doc links — only
//! exist when the feature is enabled, so it is referred to here by name.)
//!
//! ## Features
//!
//! * `binary-fmu` (off by default) — enable the `binary` module's
//!   `BinaryFmu`, which `dlopen`s a real co-simulation FMU shared library.
//!   Untested in CI (no FMU binary available); fail-loud on any
//!   unimplemented step.
//!
//! ## Fail-loud
//!
//! Every fallible operation returns [`error::FmiError`]. A dangling
//! coupling reference, a malformed `modelDescription.xml`, or a truncated
//! DIS PDU is an `Err` with a precise message — never a silent default or
//! a plausible-but-wrong number.

// The default build is 100% safe code and forbids `unsafe`. The optional
// binary-FMU bridge must `dlopen` a native library (inherently unsafe), so
// when — and only when — that feature is enabled we relax to `deny` (which,
// unlike `forbid`, can be locally overridden by the one audited `unsafe fn`
// that loads the library). `forbid` cannot be downgraded by a later
// `allow`, hence the cfg split rather than a blanket `forbid`.
#![cfg_attr(not(feature = "binary-fmu"), forbid(unsafe_code))]
#![cfg_attr(feature = "binary-fmu", deny(unsafe_code))]

pub mod cosim;
pub mod dis;
pub mod error;
pub mod federation;
pub mod fmi;

#[cfg(feature = "binary-fmu")]
pub mod binary;

pub use cosim::{CoSimMaster, Coupling, CouplingGraph, Scheme, Subsystem};
pub use dis::EntityStatePdu;
pub use error::{FmiError, Result};
pub use federation::{
    Broker, Federate, FederateBehavior, FederateId, FederationError, GrantContext, GrantRecord,
    Message, SubsystemFederate, Time, TimePolicy, Value,
};
pub use fmi::{Causality, ModelDescription, ScalarVariable};

#[cfg(test)]
mod spring_damper_validation {
    //! Benchmark-pinned validation (test requirement #1): a two-mass
    //! spring-damper, split into two [`Subsystem`]s, co-simulated and
    //! checked against the monolithic 2-DOF ODE integrated directly.
    //!
    //! Physical model (interface = the coupling spring/damper between the
    //! masses):
    //!
    //! ```text
    //! wall |--k1,c1--[ m1 ]--kc,cc--[ m2 ]
    //!
    //! m1 x1'' = -k1 x1 - c1 x1' + kc (x2 - x1) + cc (x2' - x1')
    //! m2 x2'' =                   -kc (x2 - x1) - cc (x2' - x1')
    //! ```
    //!
    //! Co-sim split: subsystem A owns (x1, v1) and reads (x2, v2);
    //! subsystem B owns (x2, v2) and reads (x1, v1). Each exposes
    //! [position, velocity] as its 2 outputs and reads the partner's
    //! [position, velocity] as its 2 inputs. The coupling graph wires
    //! A.out -> B.in and B.out -> A.in.

    use crate::cosim::{CoSimMaster, Coupling, CouplingGraph, Scheme, Subsystem};

    // Shared physical constants.
    const M1: f64 = 1.0;
    const M2: f64 = 1.5;
    const K1: f64 = 30.0; // wall spring on mass 1
    const C1: f64 = 0.4; // wall damper on mass 1
    const KC: f64 = 20.0; // coupling spring
    const CC: f64 = 0.3; // coupling damper

    // Co-sim macro-step. The explicit (non-iterative) coupling error of
    // both schemes is first-order in the macro-step; at 2e-4 the
    // Gauss-Seidel position error is ~5e-4 (comfortably < 1e-3) while
    // Jacobi's one-step lag makes it ~10% larger — see the two tests.
    const DT_MACRO: f64 = 2.0e-4;
    // Fine RK4 substeps PER macro-step. The monolithic reference uses the
    // same fine step (DT_MACRO / SUBSTEPS) so RK4 truncation cancels in the
    // comparison and what remains is the pure inter-subsystem coupling
    // error.
    const SUBSTEPS: usize = 10;
    const T_END: f64 = 2.0;

    /// Mass 1 subsystem. State (x1, v1); inputs (x2, v2); outputs (x1, v1).
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
            // Partner state held constant across the macro-step.
            let (x2, v2) = (inputs[0], inputs[1]);
            let h = dt / SUBSTEPS as f64;
            for _ in 0..SUBSTEPS {
                rk4_1dof(&mut self.x, &mut self.v, h, |x, v| {
                    MassA::accel(x, v, x2, v2)
                });
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

    /// Mass 2 subsystem. State (x2, v2); inputs (x1, v1); outputs (x2, v2).
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
            let h = dt / SUBSTEPS as f64;
            for _ in 0..SUBSTEPS {
                rk4_1dof(&mut self.x, &mut self.v, h, |x, v| {
                    MassB::accel(x, v, x1, v1)
                });
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

    /// Monolithic reference: integrate the full 4-state coupled system
    /// `[x1, v1, x2, v2]` directly with RK4 at the co-sim's fine
    /// resolution (macro-step / SUBSTEPS), no coupling lag.
    fn monolithic_reference(t_end: f64) -> [f64; 4] {
        let mut s = [1.0_f64, 0.0, -0.5, 0.0]; // initial state
        let h = DT_MACRO / SUBSTEPS as f64;
        let n = (t_end / h).round() as usize;

        let deriv = |s: &[f64; 4]| -> [f64; 4] {
            let (x1, v1, x2, v2) = (s[0], s[1], s[2], s[3]);
            let a1 = (-K1 * x1 - C1 * v1 + KC * (x2 - x1) + CC * (v2 - v1)) / M1;
            let a2 = (-KC * (x2 - x1) - CC * (v2 - v1)) / M2;
            [v1, a1, v2, a2]
        };

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

    fn add(a: &[f64; 4], b: &[f64; 4]) -> [f64; 4] {
        [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]]
    }
    fn scale(a: &[f64; 4], s: f64) -> [f64; 4] {
        [a[0] * s, a[1] * s, a[2] * s, a[3] * s]
    }

    /// Co-simulate to `t_end` and return the worst-case absolute error in
    /// (x1, x2) against the monolithic reference, sampled at the end.
    fn cosim_position_error(scheme: Scheme) -> f64 {
        let subs: Vec<Box<dyn Subsystem>> = vec![
            Box::new(MassA { x: 1.0, v: 0.0 }),
            Box::new(MassB { x: -0.5, v: 0.0 }),
        ];
        // A outputs (x1, v1) -> B inputs; B outputs (x2, v2) -> A inputs.
        let graph = CouplingGraph::from_edges(vec![
            Coupling::new(0, 0, 1, 0), // A.x1 -> B.in[0]
            Coupling::new(0, 1, 1, 1), // A.v1 -> B.in[1]
            Coupling::new(1, 0, 0, 0), // B.x2 -> A.in[0]
            Coupling::new(1, 1, 0, 1), // B.v2 -> A.in[1]
        ]);
        let mut master = CoSimMaster::new(subs, graph, scheme).expect("valid graph");

        let n = (T_END / DT_MACRO).round() as usize;
        for _ in 0..n {
            master.advance(DT_MACRO);
        }

        let reference = monolithic_reference(T_END);
        let x1 = master.outputs_of(0).unwrap()[0];
        let x2 = master.outputs_of(1).unwrap()[0];
        (x1 - reference[0]).abs().max((x2 - reference[2]).abs())
    }

    #[test]
    fn gauss_seidel_matches_monolithic_reference_under_1e_3() {
        let err = cosim_position_error(Scheme::GaussSeidel);
        assert!(
            err < 1.0e-3,
            "Gauss-Seidel co-sim error {err:.3e} must be < 1e-3 vs monolithic ref"
        );
    }

    #[test]
    fn jacobi_coupling_error_is_larger_than_gauss_seidel() {
        // Documented expectation: Jacobi carries a one-macro-step coupling
        // lag, so its error against the monolithic reference is strictly
        // larger than Gauss-Seidel's at the same macro-step.
        let err_gs = cosim_position_error(Scheme::GaussSeidel);
        let err_jac = cosim_position_error(Scheme::Jacobi);
        assert!(
            err_jac > err_gs,
            "Jacobi error {err_jac:.3e} should exceed Gauss-Seidel error {err_gs:.3e} \
             (Jacobi has a one-macro-step coupling lag)"
        );
    }
}
