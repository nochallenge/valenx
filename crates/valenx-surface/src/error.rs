//! Surface workbench error taxonomy.
//!
//! Same stable `code()` / `category()` pattern as `valenx-sketch`,
//! `valenx-cad`, etc. — every variant maps to a kebab-cased string
//! identifier that LLM / scripting layers can branch on without
//! parsing the human-readable message.

use thiserror::Error;

/// Errors raised by NURBS construction, evaluation, fill, sew,
/// extend, intersect, or trim operations.
#[derive(Debug, Error)]
pub enum SurfaceError {
    /// Knot vector failed validation — wrong length, not
    /// non-decreasing, or wrong multiplicity at the endpoints.
    #[error("bad knot vector: {reason}")]
    BadKnotVector {
        /// Human-readable reason — included in the displayed message.
        reason: String,
    },

    /// A negative or unreasonably large NURBS degree was requested.
    /// Practical degree is in [1, 9]; this guards against accidental
    /// overflow in basis-function recursion.
    #[error("bad NURBS degree: {0}")]
    BadDegree(usize),

    /// The requested parameter is outside the valid `[knots[degree],
    /// knots[n_cp]]` range.
    #[error("evaluation parameter {u} is out of the valid knot-vector range")]
    EvaluationOutOfRange {
        /// The offending parameter.
        u: f64,
    },

    /// Surface-surface intersection (or similar geometric op) failed
    /// to converge or produced no useful result.
    #[error("intersection failed: {0}")]
    IntersectionFailed(String),

    /// A differential-geometry quantity (normal, fundamental form,
    /// curvature) is undefined at the requested point because the
    /// parameterisation is locally singular — the tangents are parallel
    /// or vanish (`|Sᵤ × Sᵥ| ≈ 0`, e.g. a pole), or the metric
    /// determinant `EG − F²` is non-positive. The fail-loud counterpart
    /// of the silent `0.0` returned by the infallible curvature methods.
    #[error("degenerate geometry: {reason}")]
    DegenerateGeometry {
        /// Human-readable reason — included in the displayed message.
        reason: String,
    },

    /// Two boundary curves passed to `coons::fill` don't share an
    /// endpoint within tolerance.
    #[error("coons fill: {0}")]
    BadBoundary(String),

    /// Mismatched parameterisation between two surfaces passed to
    /// `sew::stitch` — shared edge control points don't align.
    #[error("sew: {0}")]
    SewMismatch(String),

    /// IO error wrapping std::io.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON serialise / parse error.
    #[error("ron: {0}")]
    Ron(String),
}

impl SurfaceError {
    /// Stable kebab-cased identifier; never changes across versions.
    ///
    /// Mirrors `SketchError::code()` / `CadError::code()` so the same
    /// scripting layer can branch on errors from any workbench.
    pub fn code(&self) -> &'static str {
        match self {
            SurfaceError::BadKnotVector { .. } => "surface.bad_knot_vector",
            SurfaceError::BadDegree(_) => "surface.bad_degree",
            SurfaceError::EvaluationOutOfRange { .. } => "surface.evaluation_out_of_range",
            SurfaceError::IntersectionFailed(_) => "surface.intersection_failed",
            SurfaceError::DegenerateGeometry { .. } => "surface.degenerate_geometry",
            SurfaceError::BadBoundary(_) => "surface.bad_boundary",
            SurfaceError::SewMismatch(_) => "surface.sew_mismatch",
            SurfaceError::Io(_) => "surface.io",
            SurfaceError::Ron(_) => "surface.ron",
        }
    }

    /// Coarse classification — `"input"`, `"geometry"`, or `"io"`.
    /// Useful when the caller wants to retry on IO failures but
    /// surface input/geometry failures directly to the user.
    pub fn category(&self) -> &'static str {
        match self {
            SurfaceError::BadKnotVector { .. }
            | SurfaceError::BadDegree(_)
            | SurfaceError::EvaluationOutOfRange { .. }
            | SurfaceError::BadBoundary(_) => "input",
            SurfaceError::IntersectionFailed(_)
            | SurfaceError::SewMismatch(_)
            | SurfaceError::DegenerateGeometry { .. } => "geometry",
            SurfaceError::Io(_) | SurfaceError::Ron(_) => "io",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code() {
        let cases: Vec<(SurfaceError, &str, &str)> = vec![
            (
                SurfaceError::BadKnotVector {
                    reason: "len mismatch".into(),
                },
                "surface.bad_knot_vector",
                "input",
            ),
            (SurfaceError::BadDegree(99), "surface.bad_degree", "input"),
            (
                SurfaceError::EvaluationOutOfRange { u: 2.0 },
                "surface.evaluation_out_of_range",
                "input",
            ),
            (
                SurfaceError::IntersectionFailed("no convergence".into()),
                "surface.intersection_failed",
                "geometry",
            ),
            (
                SurfaceError::BadBoundary("corners don't match".into()),
                "surface.bad_boundary",
                "input",
            ),
            (
                SurfaceError::SewMismatch("edges differ".into()),
                "surface.sew_mismatch",
                "geometry",
            ),
            (
                SurfaceError::DegenerateGeometry {
                    reason: "EG − F² ≤ 0".into(),
                },
                "surface.degenerate_geometry",
                "geometry",
            ),
            (SurfaceError::Ron("malformed".into()), "surface.ron", "io"),
        ];
        for (err, expected_code, expected_cat) in cases {
            assert_eq!(err.code(), expected_code, "wrong code for {err:?}");
            assert_eq!(err.category(), expected_cat, "wrong category for {err:?}");
        }
    }

    #[test]
    fn io_error_wraps_std() {
        let e: SurfaceError = std::io::Error::other("disk gone").into();
        assert_eq!(e.code(), "surface.io");
        assert_eq!(e.category(), "io");
    }
}
