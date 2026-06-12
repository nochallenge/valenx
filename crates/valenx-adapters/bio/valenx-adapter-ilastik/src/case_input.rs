//! `[bio.ilastik]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "ilastik.classify"
//!
//! [bio.ilastik]
//! ilastik_app     = "/opt/ilastik-1.4.0/run_ilastik.sh"
//! project         = "trained_classifier.ilp"
//! input_images    = ["raw_001.tif", "raw_002.tif"]
//! output_basename = "predictions"
//! workflow        = "Pixel Classification"
//! extra_args      = []
//! ```
//!
//! Ilastik (Hamprecht et al, GPL-3.0) is the canonical interactive
//! machine-learning tool for bioimage pixel / object classification.
//! Users train a classifier in the GUI on a handful of labelled
//! examples; production then re-runs that trained `.ilp` project
//! headlessly across new image batches.
//!
//! Headless invocation:
//! `<ilastik_app> --headless --project=<project>
//!  --output_filename_format=<output_basename>_{nickname}.h5
//!  <input_images...>`. The literal `{nickname}` token is Ilastik's
//! per-input substitution placeholder (replaced at run time with
//! each input file's stem) and must reach Ilastik unmodified.
//!
//! Ilastik ships as a portable launcher
//! (`run_ilastik.sh` on Linux / macOS, `ilastik.exe` on Windows)
//! whose absolute path is install-site-specific, so the case-input
//! schema carries the launcher path verbatim
//! (`[bio.ilastik].ilastik_app`).

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct IlastikInput {
    /// Absolute path to the Ilastik launcher
    /// (`run_ilastik.sh` on Linux / macOS or `ilastik.exe` on
    /// Windows). Ilastik installs into a versioned vendor directory
    /// (`/opt/ilastik-1.4.0/`, `C:\Program Files\ilastik-1.4.0\`)
    /// that isn't predictable from PATH alone.
    pub ilastik_app: PathBuf,
    /// Path to the trained Ilastik project file (`.ilp`) authored
    /// in the Ilastik GUI. The project embeds the workflow type,
    /// the trained classifier, and the feature-selection config.
    pub project: PathBuf,
    /// Input images to classify. Required to be non-empty —
    /// Ilastik consumes one or more raw image files (TIFF, HDF5,
    /// PNG, etc.) and emits a probability map / segmentation per
    /// input.
    pub input_images: Vec<PathBuf>,
    /// Filename stem used to build Ilastik's
    /// `--output_filename_format=<basename>_{nickname}.h5` argument
    /// and to filter `collect()`-time output artefacts.
    pub output_basename: String,
    /// Workflow name passed to Ilastik. Stored on the input for
    /// schema completeness — the workflow type is actually fixed
    /// by the trained `.ilp` project itself, so we don't pass it
    /// on the command line. Defaults to `"Pixel Classification"`.
    pub workflow: String,
    /// Additional CLI arguments appended to the Ilastik command —
    /// useful for `--export_source=Probabilities|Simple Segmentation`
    /// or `--pipeline_result_drange` overrides.
    pub extra_args: Vec<String>,
}

impl IlastikInput {
    pub fn from_case_dir(case_dir: &std::path::Path) -> Result<Self, AdapterError> {
        let case_toml = case_dir.join("case.toml");
        let text = valenx_core::io_caps::read_capped_to_string(
            &case_toml,
            valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
        )
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("read {}: {e}", case_toml.display())))?;
        let parsed: toml::Value = toml::from_str(&text).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("parse {}: {e}", case_toml.display()))
        })?;
        let block = parsed
            .get("bio")
            .and_then(|v| v.get("ilastik"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.ilastik] section",
                    case_toml.display()
                ))
            })?;

        let ilastik_app = block
            .get("ilastik_app")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.ilastik].ilastik_app required (path to \
                     run_ilastik.sh or ilastik.exe)"
                ))
            })?;
        if ilastik_app.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.ilastik].ilastik_app must not be empty"
            )));
        }

        let project = block
            .get("project")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.ilastik].project required (path to .ilp)"
                ))
            })?;
        if project.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.ilastik].project must not be empty"
            )));
        }

        let input_images: Vec<PathBuf> = match block.get("input_images") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.ilastik].input_images must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.ilastik].input_images entries must be strings"
                        ))
                    })?;
                    if s.is_empty() {
                        return Err(AdapterError::Other(anyhow::anyhow!(
                            "[bio.ilastik].input_images entries must not be empty"
                        )));
                    }
                    out.push(PathBuf::from(s));
                }
                out
            }
            None => {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.ilastik].input_images required (array of >= 1 \
                     image paths)"
                )));
            }
        };
        if input_images.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.ilastik].input_images must contain >= 1 entry"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.ilastik].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.ilastik].output_basename must not be empty"
            )));
        }

        // `workflow` is captured for schema completeness but isn't
        // passed on the command line — the trained `.ilp` project
        // itself fixes the workflow type. Default to the canonical
        // workflow name so the case.toml stays self-documenting.
        let workflow = block
            .get("workflow")
            .and_then(|v| v.as_str())
            .unwrap_or("Pixel Classification")
            .to_string();
        if workflow.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.ilastik].workflow must not be empty when set"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.ilastik].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.ilastik].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            ilastik_app: PathBuf::from(ilastik_app),
            project: PathBuf::from(project),
            input_images,
            output_basename: output_basename.to_string(),
            workflow,
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_case_with_defaults() {
        // The canonical layout: launcher + project + a couple of
        // input images. `workflow` defaults to Pixel Classification
        // and `extra_args` to empty — these are the most common
        // settings for the trained-classifier inference path.
        let d = tempdir("ilastik-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "ilastik.classify"

[bio.ilastik]
ilastik_app     = "/opt/ilastik-1.4.0/run_ilastik.sh"
project         = "trained_classifier.ilp"
input_images    = ["raw_001.tif", "raw_002.tif"]
output_basename = "predictions"
"#,
        )
        .unwrap();
        let input = IlastikInput::from_case_dir(&d).unwrap();
        assert_eq!(
            input.ilastik_app,
            PathBuf::from("/opt/ilastik-1.4.0/run_ilastik.sh")
        );
        assert_eq!(input.project, PathBuf::from("trained_classifier.ilp"));
        assert_eq!(
            input.input_images,
            vec![PathBuf::from("raw_001.tif"), PathBuf::from("raw_002.tif"),]
        );
        assert_eq!(input.output_basename, "predictions");
        assert_eq!(input.workflow, "Pixel Classification");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_input_images() {
        // Ilastik requires at least one input image to classify —
        // an empty array is a configuration error. Reject up front
        // so the user catches it at validation time rather than
        // after a confusing Ilastik usage error at run time.
        let d = tempdir("ilastik-noimgs");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "ilastik.classify"

[bio.ilastik]
ilastik_app     = "/opt/ilastik-1.4.0/run_ilastik.sh"
project         = "trained_classifier.ilp"
input_images    = []
output_basename = "predictions"
"#,
        )
        .unwrap();
        let err = IlastikInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("input_images"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        // Without `[bio.ilastik]` we have nothing to work with —
        // surface the missing-section error up front rather than
        // letting prepare() fail with a cryptic key lookup.
        let d = tempdir("ilastik-nosec");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = IlastikInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.ilastik]"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
