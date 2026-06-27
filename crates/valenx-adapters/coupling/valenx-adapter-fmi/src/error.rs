//! Fail-loud error taxonomy for the co-simulation master.
//!
//! Every fallible path in this crate returns [`FmiError`] rather than a
//! plausible-but-wrong default. Per the Valenx non-negotiables, a wrong
//! number is worse than a crash: if a coupling graph is malformed, a
//! `modelDescription.xml` cannot be parsed, or a DIS PDU buffer is the
//! wrong length, we return an `Err` describing exactly what was wrong —
//! we never silently clamp, pad, or fabricate.

use thiserror::Error;

/// Errors raised by the native co-simulation master, the FMI importer,
/// and the DIS Entity State PDU codec.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum FmiError {
    /// A coupling edge referenced a subsystem index that does not exist.
    #[error(
        "coupling edge {edge_index}: subsystem index {bad_index} is out of \
         range (only {n_subsystems} subsystem(s) registered)"
    )]
    SubsystemIndexOutOfRange {
        /// Position of the offending edge in the coupling graph.
        edge_index: usize,
        /// The out-of-range subsystem index that was referenced.
        bad_index: usize,
        /// Number of subsystems actually registered on the master.
        n_subsystems: usize,
    },

    /// A coupling edge referenced an output or input port (a scalar
    /// channel) that the named subsystem does not expose.
    #[error(
        "coupling edge {edge_index}: {side} port {bad_port} is out of range \
         for subsystem {subsystem} (which has {n_ports} {side} port(s))"
    )]
    PortIndexOutOfRange {
        /// Position of the offending edge in the coupling graph.
        edge_index: usize,
        /// Which side of the edge was bad: `"output"` or `"input"`.
        side: &'static str,
        /// The subsystem index whose port count was exceeded.
        subsystem: usize,
        /// The out-of-range port index that was referenced.
        bad_port: usize,
        /// Number of ports of that side the subsystem declares.
        n_ports: usize,
    },

    /// Two coupling edges drive the SAME (subsystem, input) port. A
    /// co-simulation input may be fed by at most one source, otherwise
    /// the value the subsystem sees is ambiguous.
    #[error(
        "coupling input (subsystem {subsystem}, port {port}) is driven by \
         more than one edge — an input may have at most one source"
    )]
    DuplicateInputSource {
        /// Subsystem index whose input is over-driven.
        subsystem: usize,
        /// Input port index that is over-driven.
        port: usize,
    },

    /// `modelDescription.xml` could not be parsed (malformed XML, missing
    /// required attribute, or unexpected structure).
    #[error("modelDescription.xml parse error: {0}")]
    ModelDescriptionParse(String),

    /// A DIS PDU byte buffer was shorter than the fixed standard layout
    /// requires (so a field would read past the end of the slice).
    #[error(
        "DIS PDU buffer too short: need at least {needed} bytes for the \
         {what}, got {got}"
    )]
    PduTooShort {
        /// Number of bytes the layout requires.
        needed: usize,
        /// Number of bytes actually supplied.
        got: usize,
        /// Which structure was being decoded when the shortfall hit.
        what: &'static str,
    },

    /// A DIS PDU header declared a PDU type other than Entity State (1).
    /// This crate's codec only handles the Entity State PDU.
    #[error(
        "unsupported DIS PDU type {got}: this codec only decodes the Entity \
         State PDU (type {expected})"
    )]
    UnsupportedPduType {
        /// PDU type byte found in the header.
        got: u8,
        /// The only PDU type this codec supports.
        expected: u8,
    },

    /// A binary FMU could not be loaded (only reachable with the
    /// `binary-fmu` feature enabled).
    #[error("binary FMU load error: {0}")]
    BinaryFmu(String),

    /// A fault in the HELICS-style co-simulation **federation** / distributed
    /// time-coordination layer (see [`crate::federation`]): an unknown or
    /// duplicate federate / publication / subscription / endpoint name, a
    /// negative period or lookahead, a non-monotone time request, a message
    /// sent into the past, and so on. The machine-readable
    /// [`crate::federation::FederationError`] `code` is the supported way to
    /// match a specific failure; `message` is the human-readable detail.
    #[error("federation error [{code:?}]: {message}")]
    Federation {
        /// Stable, machine-readable classification of the fault.
        code: crate::federation::FederationError,
        /// Human-readable detail of exactly what went wrong.
        message: String,
    },

    /// The implicit (iterative) coupling solver did not converge within the
    /// allowed number of iterations. Returned by
    /// [`crate::implicit::coupled_step`] when the fixed-point loop exhausts
    /// `max_iter` without the infinity-norm of the residual falling below
    /// `tol`. Also returned when a subsystem produces a non-finite output
    /// (NaN / Inf) mid-iteration.
    #[error(
        "implicit coupling did not converge after {iterations} iterations \
         (‖Δy‖_∞ = {final_residual:.3e}): {reason}"
    )]
    NotConverged {
        /// Number of iterations performed before giving up.
        iterations: usize,
        /// Infinity-norm of the residual at the last iteration.
        final_residual: f64,
        /// Human-readable reason (scheme, tolerance, etc.).
        reason: String,
    },

    /// A configuration error in [`crate::implicit::coupled_step`]: invalid
    /// `tol`, `max_iter`, `dt`, or relaxation parameter.
    #[error("coupled_step configuration error: {0}")]
    BadCoupledStep(String),
}

/// Convenience result alias for this crate.
pub type Result<T> = std::result::Result<T, FmiError>;
