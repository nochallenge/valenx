//! Replay evaluator for [`crate::feature::Feature::ImportedAdvanced`]
//! (Phase 20 — STEP AP242 + IGES trimmed surface).
//!
//! Re-reads the source file at replay time, like
//! [`crate::ops::imported_solid::evaluate`], but uses the AP242-aware
//! import path that recovers product hierarchy + GD&T metadata too.
//! The metadata stays on the [`crate::feature::ImportedAdvancedParams`]
//! struct so dependent features and the UI can query it without
//! re-parsing.

use std::path::Path;

use valenx_cad::Solid;

use crate::feature::ImportedAdvancedParams;
use crate::FeatureError;

/// Replay an [`crate::feature::Feature::ImportedAdvanced`] feature.
///
/// v1 dispatches the geometry through [`valenx_step_iges::import`] (the
/// same path the basic [`crate::ops::imported_solid::evaluate`] takes).
/// AP242 metadata captured at original import is held on the
/// `params` struct and re-presented to the UI without re-parsing.
///
/// # Errors
///
/// - [`FeatureError::Io`] if the source file is missing / unreadable.
/// - [`FeatureError::BadParameter`] for unsupported formats or parse
///   failures.
pub fn evaluate(params: &ImportedAdvancedParams) -> Result<Solid, FeatureError> {
    let path = Path::new(&params.source_path);
    valenx_step_iges::import(path).map_err(|e| {
        if let valenx_step_iges::StepIgesError::Io(io) = e {
            FeatureError::Io(io)
        } else {
            FeatureError::BadParameter {
                name: "source_path",
                reason: format!("AP242 import failed: {e}"),
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_file_yields_io_error() {
        let params = ImportedAdvancedParams {
            source_path: "C:\\definitely\\does\\not\\exist\\ap242.step".to_string(),
            ..Default::default()
        };
        let err = evaluate(&params).unwrap_err();
        assert!(matches!(err, FeatureError::Io(_)), "got {err:?}");
    }

    #[test]
    fn unsupported_extension_yields_bad_parameter() {
        let tmp = std::env::temp_dir().join("valenx_imported_advanced_test.obj");
        std::fs::write(&tmp, "fake obj contents").unwrap();
        let params = ImportedAdvancedParams {
            source_path: tmp.to_string_lossy().to_string(),
            ..Default::default()
        };
        let err = evaluate(&params).unwrap_err();
        let _ = std::fs::remove_file(&tmp);
        assert!(
            matches!(err, FeatureError::BadParameter { .. }),
            "got {err:?}",
        );
    }

    #[test]
    fn params_default_is_empty() {
        let p = ImportedAdvancedParams::default();
        assert!(p.source_path.is_empty());
        assert!(p.product_path.is_empty());
        assert!(!p.is_ap242);
    }
}
