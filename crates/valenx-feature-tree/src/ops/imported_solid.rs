//! Replay evaluator for [`crate::feature::Feature::ImportedSolid`].
//!
//! The variant holds only a `source_path`; we re-read the STEP / IGES
//! file via [`valenx_step_iges::import`] every replay. That keeps the
//! `.valenx` project envelope small (no embedded triangle meshes) at
//! the cost of repeating the parse on every dependent edit. v2 could
//! cache the imported solid in-memory keyed by `(path, mtime)`.

use std::path::Path;

use valenx_cad::Solid;

use crate::feature::ImportedSolidParams;
use crate::FeatureError;

/// Replay an [`crate::feature::Feature::ImportedSolid`] feature.
///
/// # Errors
///
/// - [`FeatureError::Io`] if the source file is missing / unreadable.
/// - [`FeatureError::BadParameter`] if the import returns any other
///   error (unsupported format, malformed STEP, etc.).
pub fn evaluate(params: &ImportedSolidParams) -> Result<Solid, FeatureError> {
    let path = Path::new(&params.source_path);
    valenx_step_iges::import(path).map_err(|e| {
        // Pre-wrap any IO failure as FeatureError::Io so callers see
        // a consistent "file went missing" surface.
        if let valenx_step_iges::StepIgesError::Io(io) = e {
            FeatureError::Io(io)
        } else {
            FeatureError::BadParameter {
                name: "source_path",
                reason: format!("import failed: {e}"),
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_file_yields_io_error() {
        let params = ImportedSolidParams {
            source_path: "C:\\definitely\\does\\not\\exist\\nope.step".to_string(),
        };
        let err = evaluate(&params).unwrap_err();
        assert!(matches!(err, FeatureError::Io(_)), "got {err:?}");
    }

    #[test]
    fn unsupported_extension_yields_bad_parameter() {
        // Use a temp file with an .obj extension — extension dispatch
        // returns Unsupported, which maps to BadParameter.
        let tmp = std::env::temp_dir().join("valenx_imported_solid_test.obj");
        std::fs::write(&tmp, "fake obj contents").unwrap();
        let params = ImportedSolidParams {
            source_path: tmp.to_string_lossy().to_string(),
        };
        let err = evaluate(&params).unwrap_err();
        let _ = std::fs::remove_file(&tmp);
        assert!(
            matches!(err, FeatureError::BadParameter { .. }),
            "got {err:?}"
        );
    }
}
