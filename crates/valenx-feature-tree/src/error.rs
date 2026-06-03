//! Feature-tree error taxonomy.

use thiserror::Error;

use valenx_cad::CadError;
use valenx_fillet::FilletError;
use valenx_sketch::SketchError;

/// Errors raised by feature-tree construction, mutation, or replay.
#[derive(Debug, Error)]
pub enum FeatureError {
    /// Tried to look up a sketch by index that doesn't exist in the tree.
    #[error("sketch {0} not found in feature tree")]
    UnknownSketch(usize),

    /// Tried to look up a feature by index that doesn't exist in the tree.
    #[error("feature {0} not found in feature tree")]
    UnknownFeature(usize),

    /// Wrapped sketcher error (constraint / solver / persistence).
    #[error("sketch: {0}")]
    SketchError(#[from] SketchError),

    /// Wrapped CAD-kernel error (boolean / primitive / tessellation).
    #[error("cad: {0}")]
    CadError(#[from] CadError),

    /// Wrapped fillet/chamfer error (mesh-domain fillet pipeline).
    #[error("fillet: {0}")]
    FilletError(#[from] FilletError),

    /// Sketch produced a closed profile with too few edges to extrude
    /// (i.e. no closed loop, or a degenerate triangle).
    #[error("sketch has no extrudable profile")]
    EmptyProfile,

    /// User-supplied parameter to a feature was out of range or nonsense.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter (e.g. `"depth"`, `"angle"`).
        name: &'static str,
        /// Human-readable explanation of what was wrong.
        reason: String,
    },

    /// IO error wrapping std::io.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// RON serialization / deserialization failure.
    #[error("ron: {0}")]
    Ron(String),
}

/// Coarse category an LLM can use to decide who to escalate to.
///
/// Matches the pattern used by `valenx-dock`'s `ErrorCategory`.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User supplied bad input data (fix the sketch / tree / file).
    Input,
    /// User-tunable knob out of range (fix the parameters).
    Config,
    /// Transient runtime failure (retry may help).
    Runtime,
    /// Bug in valenx-feature-tree or an underlying kernel (file a report).
    Internal,
}

impl FeatureError {
    /// Stable kebab-cased identifier; never changes across versions.
    /// New variants get new codes.
    pub fn code(&self) -> &'static str {
        match self {
            FeatureError::UnknownSketch(_) => "feature_tree.unknown_sketch",
            FeatureError::UnknownFeature(_) => "feature_tree.unknown_feature",
            FeatureError::SketchError(_) => "feature_tree.sketch_error",
            FeatureError::CadError(_) => "feature_tree.cad_error",
            FeatureError::FilletError(_) => "feature_tree.fillet_error",
            FeatureError::EmptyProfile => "feature_tree.empty_profile",
            FeatureError::BadParameter { .. } => "feature_tree.bad_parameter",
            FeatureError::Io(_) => "feature_tree.io",
            FeatureError::Ron(_) => "feature_tree.ron",
        }
    }

    /// High-level classification for LLM / UI routing.
    pub fn category(&self) -> ErrorCategory {
        match self {
            FeatureError::UnknownSketch(_)
            | FeatureError::UnknownFeature(_)
            | FeatureError::EmptyProfile
            | FeatureError::SketchError(_) => ErrorCategory::Input,
            FeatureError::BadParameter { .. } => ErrorCategory::Config,
            FeatureError::Io(_) | FeatureError::Ron(_) => ErrorCategory::Runtime,
            FeatureError::CadError(_) | FeatureError::FilletError(_) => ErrorCategory::Internal,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code_and_category() {
        let cases: Vec<(FeatureError, &'static str, ErrorCategory)> = vec![
            (
                FeatureError::UnknownSketch(2),
                "feature_tree.unknown_sketch",
                ErrorCategory::Input,
            ),
            (
                FeatureError::UnknownFeature(5),
                "feature_tree.unknown_feature",
                ErrorCategory::Input,
            ),
            (
                FeatureError::EmptyProfile,
                "feature_tree.empty_profile",
                ErrorCategory::Input,
            ),
            (
                FeatureError::BadParameter {
                    name: "depth",
                    reason: "must be nonzero".into(),
                },
                "feature_tree.bad_parameter",
                ErrorCategory::Config,
            ),
            (
                FeatureError::Ron("parse fail".into()),
                "feature_tree.ron",
                ErrorCategory::Runtime,
            ),
            (
                FeatureError::CadError(CadError::EmptyResult),
                "feature_tree.cad_error",
                ErrorCategory::Internal,
            ),
            (
                FeatureError::FilletError(FilletError::EmptyMesh),
                "feature_tree.fillet_error",
                ErrorCategory::Internal,
            ),
        ];
        for (err, expected_code, expected_cat) in cases {
            assert_eq!(err.code(), expected_code, "wrong code for {err:?}");
            assert_eq!(err.category(), expected_cat, "wrong category for {err:?}");
        }
    }
}
