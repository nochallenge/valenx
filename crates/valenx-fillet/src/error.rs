//! Fillet / chamfer error taxonomy.
//!
//! Same shape as `valenx-feature-tree`'s `FeatureError`: stable
//! `code()` strings an LLM can branch on plus a coarse
//! [`ErrorCategory`] for UI / escalation routing.

use thiserror::Error;

/// Errors raised by fillet / chamfer application.
#[derive(Debug, Error)]
pub enum FilletError {
    /// Input mesh had no triangle elements at all.
    #[error("input mesh has no triangles")]
    EmptyMesh,

    /// A filletable edge had effectively zero length (its two
    /// endpoints were coincident); cannot offset along a zero vector.
    #[error("degenerate edge from vertex {from} to {to}")]
    DegenerateEdge {
        /// Index of the edge's first endpoint in the source mesh.
        from: usize,
        /// Index of the edge's second endpoint in the source mesh.
        to: usize,
    },

    /// Requested radius is larger than the available edge length —
    /// the strip would self-intersect.
    #[error("radius {radius} too large for edge of length {edge_length}")]
    RadiusTooLarge {
        /// The radius the caller asked for.
        radius: f64,
        /// Length of the shortest edge that violated the bound.
        edge_length: f64,
    },

    /// Generic mesh-side failure (degenerate triangle, missing
    /// element block, etc.). String-wrapped — `valenx-mesh` does not
    /// yet expose a typed error enum.
    #[error("mesh: {0}")]
    Mesh(String),

    /// User-supplied parameter was out of range or nonsense.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter (e.g. `"radius"`, `"angle"`).
        name: &'static str,
        /// Human-readable explanation.
        reason: String,
    },
}

/// Coarse category an LLM / UI can branch on to decide who to escalate
/// the error to.
///
/// Mirrors `valenx-feature-tree::ErrorCategory`.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User supplied bad geometry (fix the source mesh / model).
    Input,
    /// User-tunable knob out of range (fix the radius / angle).
    Config,
    /// Transient or environmental failure (retry may help).
    Runtime,
    /// Bug in valenx-fillet or its mesh layer (file a report).
    Internal,
}

impl FilletError {
    /// Stable kebab-cased identifier; never changes across versions.
    pub fn code(&self) -> &'static str {
        match self {
            FilletError::EmptyMesh => "fillet.empty_mesh",
            FilletError::DegenerateEdge { .. } => "fillet.degenerate_edge",
            FilletError::RadiusTooLarge { .. } => "fillet.radius_too_large",
            FilletError::Mesh(_) => "fillet.mesh",
            FilletError::BadParameter { .. } => "fillet.bad_parameter",
        }
    }

    /// High-level classification for LLM / UI routing.
    pub fn category(&self) -> ErrorCategory {
        match self {
            FilletError::EmptyMesh | FilletError::DegenerateEdge { .. } => ErrorCategory::Input,
            FilletError::RadiusTooLarge { .. } | FilletError::BadParameter { .. } => {
                ErrorCategory::Config
            }
            FilletError::Mesh(_) => ErrorCategory::Runtime,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code_and_category() {
        let cases: Vec<(FilletError, &'static str, ErrorCategory)> = vec![
            (
                FilletError::EmptyMesh,
                "fillet.empty_mesh",
                ErrorCategory::Input,
            ),
            (
                FilletError::DegenerateEdge { from: 1, to: 2 },
                "fillet.degenerate_edge",
                ErrorCategory::Input,
            ),
            (
                FilletError::RadiusTooLarge {
                    radius: 1.0,
                    edge_length: 0.5,
                },
                "fillet.radius_too_large",
                ErrorCategory::Config,
            ),
            (
                FilletError::Mesh("triangle missing nodes".into()),
                "fillet.mesh",
                ErrorCategory::Runtime,
            ),
            (
                FilletError::BadParameter {
                    name: "radius",
                    reason: "must be positive".into(),
                },
                "fillet.bad_parameter",
                ErrorCategory::Config,
            ),
        ];
        for (err, expected_code, expected_cat) in cases {
            assert_eq!(err.code(), expected_code, "wrong code for {err:?}");
            assert_eq!(err.category(), expected_cat, "wrong category for {err:?}");
        }
    }
}
