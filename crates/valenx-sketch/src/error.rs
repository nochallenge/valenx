//! Sketcher error taxonomy.

use thiserror::Error;

/// Errors raised by sketch construction, mutation, or solving.
#[derive(Debug, Error)]
pub enum SketchError {
    /// Tried to look up an entity that doesn't exist in this sketch.
    #[error("entity {0} not found in sketch")]
    UnknownEntity(usize),

    /// Tried to apply a constraint to entities of the wrong kind
    /// (e.g. Tangent with two points).
    #[error("constraint type mismatch: {0}")]
    ConstraintTypeMismatch(String),

    /// Constraint references an entity that doesn't exist.
    #[error("constraint references unknown entity {0}")]
    ConstraintUnknownEntity(usize),

    /// Solver hit max iterations without converging.
    #[error("solver did not converge in {iterations} iterations (final residual {residual_norm})")]
    SolverDiverged {
        /// Total Newton iterations attempted before giving up.
        iterations: usize,
        /// Final L2 norm of the residual vector.
        residual_norm: f64,
    },

    /// Two or more constraints contradict each other (no solution exists).
    #[error("constraint system is inconsistent (residual {residual_norm} after {iterations} iterations)")]
    SolverInconsistent {
        /// Iterations performed before declaring the system inconsistent.
        iterations: usize,
        /// Final L2 norm of the residual vector (above tolerance).
        residual_norm: f64,
    },

    /// A loaded/deserialized sketch carries an entity whose variable
    /// handle indexes past the end of `vars` (corrupt, hand-edited, or
    /// version-skewed `.valenx`). Caught by [`crate::Sketch::validate`]
    /// at load so it can never panic during feature-tree replay. R33 H1.
    #[error("entity {entity} field `{field}` references variable index {var} but the sketch has only {len} variables")]
    CorruptHandle {
        /// 0-based index of the offending entity in `Sketch::entities`.
        entity: usize,
        /// Name of the handle field (e.g. `"x_var"`, `"radius_var"`).
        field: &'static str,
        /// The out-of-range variable index that was found.
        var: usize,
        /// Number of variables actually present in `Sketch::vars`.
        len: usize,
    },

    /// IO error wrapping std::io.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON parse error.
    #[error("ron: {0}")]
    Ron(String),
}

/// Stable code an LLM or scripting layer can branch on. Same pattern
/// as `valenx-dock`'s `DockError::code()`.
impl SketchError {
    /// Stable kebab-cased identifier; never changes across versions.
    pub fn code(&self) -> &'static str {
        match self {
            SketchError::UnknownEntity(_) => "sketch.unknown_entity",
            SketchError::ConstraintTypeMismatch(_) => "sketch.constraint_type_mismatch",
            SketchError::ConstraintUnknownEntity(_) => "sketch.constraint_unknown_entity",
            SketchError::SolverDiverged { .. } => "sketch.solver_diverged",
            SketchError::SolverInconsistent { .. } => "sketch.solver_inconsistent",
            SketchError::CorruptHandle { .. } => "sketch.corrupt_handle",
            SketchError::Io(_) => "sketch.io",
            SketchError::Ron(_) => "sketch.ron",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code() {
        let cases = vec![
            (SketchError::UnknownEntity(7), "sketch.unknown_entity"),
            (
                SketchError::ConstraintTypeMismatch("tangent on 2 points".into()),
                "sketch.constraint_type_mismatch",
            ),
            (
                SketchError::ConstraintUnknownEntity(3),
                "sketch.constraint_unknown_entity",
            ),
            (
                SketchError::SolverDiverged {
                    iterations: 100,
                    residual_norm: 1.5,
                },
                "sketch.solver_diverged",
            ),
            (
                SketchError::SolverInconsistent {
                    iterations: 50,
                    residual_norm: 2.3,
                },
                "sketch.solver_inconsistent",
            ),
            (
                SketchError::CorruptHandle {
                    entity: 1,
                    field: "x_var",
                    var: 999,
                    len: 4,
                },
                "sketch.corrupt_handle",
            ),
        ];
        for (err, expected) in cases {
            assert_eq!(err.code(), expected, "wrong code for {err:?}");
        }
    }
}
