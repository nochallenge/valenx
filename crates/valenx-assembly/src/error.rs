//! Assembly error taxonomy.

use thiserror::Error;

/// Errors raised by assembly construction, mutation, persistence, or
/// solving.
#[derive(Debug, Error)]
pub enum AssemblyError {
    /// Tried to look up a part by id that doesn't exist.
    #[error("part {0} not found in assembly")]
    UnknownPart(usize),

    /// Tried to look up a mate by id that doesn't exist.
    #[error("mate {0} not found in assembly")]
    UnknownMate(usize),

    /// Tried to look up a joint by id that doesn't exist.
    #[error("joint {0} not found in assembly")]
    UnknownJoint(usize),

    /// Solver ran out of iterations without satisfying all the mates.
    #[error("mate system is infeasible (residual {residual_norm} after {iterations} iterations)")]
    MateInfeasible {
        /// Iterations performed before giving up.
        iterations: usize,
        /// Final L2 norm of the residual vector (above tolerance).
        residual_norm: f64,
    },

    /// Caller passed something the assembly rejected — negative
    /// distance target, zero-norm axis vector, mismatched part ids
    /// across a mate, etc. Always carries the parameter name and a
    /// human-readable reason.
    #[error("bad parameter '{name}': {reason}")]
    BadParameter {
        /// Parameter name (e.g. `"axis_dir"`, `"target"`).
        name: &'static str,
        /// Short explanation of why the value is bad.
        reason: String,
    },

    /// IO error wrapping std::io (persistence path).
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON parse / serialize error (persistence path).
    #[error("ron: {0}")]
    Ron(String),
}

impl AssemblyError {
    /// Stable kebab-cased identifier; never changes across versions.
    /// Same pattern as `SketchError::code()`.
    pub fn code(&self) -> &'static str {
        match self {
            AssemblyError::UnknownPart(_) => "assembly.unknown_part",
            AssemblyError::UnknownMate(_) => "assembly.unknown_mate",
            AssemblyError::UnknownJoint(_) => "assembly.unknown_joint",
            AssemblyError::MateInfeasible { .. } => "assembly.mate_infeasible",
            AssemblyError::BadParameter { .. } => "assembly.bad_parameter",
            AssemblyError::Io(_) => "assembly.io",
            AssemblyError::Ron(_) => "assembly.ron",
        }
    }

    /// Broad category bucket for telemetry — `"input"` for caller
    /// mistakes, `"solver"` for numerics, `"io"` for persistence,
    /// `"lookup"` for missing entity ids.
    pub fn category(&self) -> &'static str {
        match self {
            AssemblyError::UnknownPart(_)
            | AssemblyError::UnknownMate(_)
            | AssemblyError::UnknownJoint(_) => "lookup",
            AssemblyError::MateInfeasible { .. } => "solver",
            AssemblyError::BadParameter { .. } => "input",
            AssemblyError::Io(_) | AssemblyError::Ron(_) => "io",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code_and_category() {
        let cases: Vec<(AssemblyError, &'static str, &'static str)> = vec![
            (
                AssemblyError::UnknownPart(3),
                "assembly.unknown_part",
                "lookup",
            ),
            (
                AssemblyError::UnknownMate(0),
                "assembly.unknown_mate",
                "lookup",
            ),
            (
                AssemblyError::UnknownJoint(7),
                "assembly.unknown_joint",
                "lookup",
            ),
            (
                AssemblyError::MateInfeasible {
                    iterations: 50,
                    residual_norm: 2.0,
                },
                "assembly.mate_infeasible",
                "solver",
            ),
            (
                AssemblyError::BadParameter {
                    name: "axis_dir",
                    reason: "must be non-zero".into(),
                },
                "assembly.bad_parameter",
                "input",
            ),
            (
                AssemblyError::Ron("syntax error".into()),
                "assembly.ron",
                "io",
            ),
        ];
        for (err, code, cat) in cases {
            assert_eq!(err.code(), code, "wrong code for {err:?}");
            assert_eq!(err.category(), cat, "wrong category for {err:?}");
        }
    }

    #[test]
    fn display_includes_key_fields() {
        let err = AssemblyError::MateInfeasible {
            iterations: 100,
            residual_norm: 1.5,
        };
        let s = err.to_string();
        assert!(s.contains("100"));
        assert!(s.contains("1.5"));

        let err = AssemblyError::BadParameter {
            name: "target",
            reason: "distance must be positive".into(),
        };
        let s = err.to_string();
        assert!(s.contains("target"));
        assert!(s.contains("positive"));
    }
}
