//! Truss-statics error taxonomy.
//!
//! Every fallible constructor and the solver funnel their failure modes
//! through [`TrussError`]. The variants split cleanly into *modelling*
//! mistakes the caller can fix (a bad node index, two coincident nodes,
//! a non-finite load) and the *physics* outcome that the assembled truss
//! is not a single statically-determinate, stable structure
//! ([`TrussError::Singular`] / [`TrussError::NotDeterminate`]).

use thiserror::Error;

/// Errors raised while building or solving a [`crate::Truss`].
#[derive(Debug, Error, Clone, PartialEq)]
pub enum TrussError {
    /// A coordinate, load, or direction component was not finite
    /// (NaN or ±∞), which would poison the assembled linear system.
    #[error("non-finite value for `{what}`: {value}")]
    NonFinite {
        /// Name of the offending quantity (e.g. `"node.x"`).
        what: &'static str,
        /// The bad value, formatted for the message.
        value: f64,
    },

    /// A member or load referenced a node index that does not exist.
    #[error("node index {index} out of range (truss has {count} nodes)")]
    NodeIndexOutOfRange {
        /// The out-of-range index that was requested.
        index: usize,
        /// Number of nodes currently in the truss.
        count: usize,
    },

    /// A member's two end nodes are the same node, or are coincident in
    /// space, so the member has zero length and no defined direction.
    #[error("member {member} is degenerate (zero length between nodes {a} and {b})")]
    DegenerateMember {
        /// Index of the offending member.
        member: usize,
        /// First end-node index.
        a: usize,
        /// Second end-node index.
        b: usize,
    },

    /// A roller-support sliding direction (the line the joint is free to
    /// move along) was given as a zero-length vector, so the reaction
    /// normal is undefined.
    #[error("support on node {node} has a zero-length roller direction")]
    ZeroRollerDirection {
        /// Index of the node carrying the bad support.
        node: usize,
    },

    /// The number of unknowns (member axial forces + reaction
    /// components) does not equal the number of joint-equilibrium
    /// equations `2·N`, so the truss is not statically determinate and
    /// the method of joints cannot solve it in closed form.
    #[error(
        "truss is not statically determinate: {unknowns} unknowns \
         (m={members} member forces + r={reactions} reactions) \
         vs {equations} equilibrium equations (2·{nodes} nodes); \
         need unknowns == equations"
    )]
    NotDeterminate {
        /// Total unknowns `m + r`.
        unknowns: usize,
        /// Number of member axial-force unknowns `m`.
        members: usize,
        /// Number of reaction-component unknowns `r`.
        reactions: usize,
        /// Number of equilibrium equations `2·N`.
        equations: usize,
        /// Number of nodes `N`.
        nodes: usize,
    },

    /// The equilibrium system is square but singular: the geometry is a
    /// mechanism (insufficient or badly-aligned constraints) or contains
    /// redundant/parallel reactions, so no unique force state exists.
    #[error(
        "equilibrium system is singular (rank-deficient): the truss is a \
         mechanism or is improperly constrained — check supports and geometry"
    )]
    Singular,

    /// A truss with no nodes (or no members) was handed to the solver.
    #[error("empty truss: {what}")]
    Empty {
        /// Which part was empty (`"no nodes"` / `"no members"`).
        what: &'static str,
    },
}

/// Coarse category for routing / metrics, mirroring the sibling crates.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Caller-supplied modelling input is malformed.
    Input,
    /// The assembled structure is not a solvable determinate system
    /// (a property of the model as a whole, not one bad field).
    Algorithm,
}

impl TrussError {
    /// Stable kebab-cased identifier, handy for logging / telemetry.
    pub fn code(&self) -> &'static str {
        match self {
            TrussError::NonFinite { .. } => "truss.non_finite",
            TrussError::NodeIndexOutOfRange { .. } => "truss.node_index_out_of_range",
            TrussError::DegenerateMember { .. } => "truss.degenerate_member",
            TrussError::ZeroRollerDirection { .. } => "truss.zero_roller_direction",
            TrussError::NotDeterminate { .. } => "truss.not_determinate",
            TrussError::Singular => "truss.singular",
            TrussError::Empty { .. } => "truss.empty",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            TrussError::NonFinite { .. }
            | TrussError::NodeIndexOutOfRange { .. }
            | TrussError::DegenerateMember { .. }
            | TrussError::ZeroRollerDirection { .. }
            | TrussError::Empty { .. } => ErrorCategory::Input,
            TrussError::NotDeterminate { .. } | TrussError::Singular => ErrorCategory::Algorithm,
        }
    }
}
