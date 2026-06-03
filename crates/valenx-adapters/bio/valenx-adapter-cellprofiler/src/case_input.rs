//! `[bio.cellprofiler]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "cellprofiler.segment"
//!
//! [bio.cellprofiler]
//! pipeline        = "segmentation.cppipe"
//! input_dir       = "raw_images"
//! output_basename = "results"
//! python          = "python3"
//! extra_args      = []
//! ```
//!
//! CellProfiler is the Broad Institute's pipeline-driven cell
//! segmentation + measurement suite. Pipelines are authored in the
//! GUI and saved as `.cppipe` (text) or `.cpproj` (binary project)
//! files; the headless invocation
//! (`cellprofiler -c -r -p <pipeline> -i <input_dir> -o <output>`)
//! consumes that pipeline and a directory of input images, then
//! writes per-image segmentation masks + per-object measurements
//! into the output directory.
//!
//! `python` is kept around for the fallback path
//! (`<python> -m cellprofiler ...`) when the standalone
//! `cellprofiler` shim isn't on PATH but the package was installed
//! into a Python environment (typical conda layout).

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct CellProfilerInput {
    /// Path to the pipeline file (`.cppipe` text pipeline or
    /// `.cpproj` binary project) authored in the CellProfiler GUI.
    pub pipeline: PathBuf,
    /// Directory containing the input images CellProfiler should
    /// process. CellProfiler walks this directory looking for image
    /// files matching the pipeline's `Images` module configuration.
    pub input_dir: PathBuf,
    /// Output directory name (resolved relative to the workdir).
    /// CellProfiler writes per-image segmentation masks +
    /// `.csv` measurement tables underneath.
    pub output_basename: String,
    /// Python interpreter to use for the `<python> -m cellprofiler`
    /// fallback when the standalone `cellprofiler` shim isn't on
    /// PATH. Defaults to `"python3"`.
    pub python: String,
    /// Additional CLI arguments appended to the CellProfiler
    /// command — useful for module-specific overrides like
    /// `--first-image-set` / `--last-image-set` or
    /// `--data-file=<groups.csv>`.
    pub extra_args: Vec<String>,
}

impl CellProfilerInput {
    pub fn from_case_dir(case_dir: &std::path::Path) -> Result<Self, AdapterError> {
        let case_toml = case_dir.join("case.toml");
        let text = valenx_core::io_caps::read_capped_to_string(
            &case_toml,
            valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
        )
        .map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("read {}: {e}", case_toml.display()))
        })?;
        let parsed: toml::Value = toml::from_str(&text).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("parse {}: {e}", case_toml.display()))
        })?;
        let block = parsed
            .get("bio")
            .and_then(|v| v.get("cellprofiler"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.cellprofiler] section",
                    case_toml.display()
                ))
            })?;

        let pipeline = block
            .get("pipeline")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.cellprofiler].pipeline required (path to .cppipe / .cpproj)"
                ))
            })?;
        if pipeline.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cellprofiler].pipeline must not be empty"
            )));
        }

        let input_dir = block
            .get("input_dir")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.cellprofiler].input_dir required"))
            })?;
        if input_dir.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cellprofiler].input_dir must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.cellprofiler].output_basename required"
                ))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cellprofiler].output_basename must not be empty"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.cellprofiler].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.cellprofiler].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            pipeline: PathBuf::from(pipeline),
            input_dir: PathBuf::from(input_dir),
            output_basename: output_basename.to_string(),
            python,
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
        // A bare pipeline + input_dir + output_basename — the
        // python knob defaults to `python3` and extra_args to
        // empty. This is the canonical "run a saved pipeline
        // against a folder of images" setup.
        let d = tempdir("cellprofiler-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cellprofiler.segment"

[bio.cellprofiler]
pipeline        = "segmentation.cppipe"
input_dir       = "raw_images"
output_basename = "results"
"#,
        )
        .unwrap();
        let input = CellProfilerInput::from_case_dir(&d).unwrap();
        assert_eq!(input.pipeline, PathBuf::from("segmentation.cppipe"));
        assert_eq!(input.input_dir, PathBuf::from("raw_images"));
        assert_eq!(input.output_basename, "results");
        assert_eq!(input.python, "python3");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_python_and_extra_args_overrides() {
        // Pinned conda interpreter + image-set range pass-through.
        // CellProfiler's `--first-image-set` / `--last-image-set`
        // pair is the standard way to slice a large batch into
        // parallel chunks; the adapter just passes them through.
        let d = tempdir("cellprofiler-over");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cellprofiler.segment"

[bio.cellprofiler]
pipeline        = "pipelines/seg.cppipe"
input_dir       = "/data/plate1/images"
output_basename = "results"
python          = "/opt/conda/envs/cellprofiler/bin/python"
extra_args      = ["--first-image-set", "1", "--last-image-set", "100"]
"#,
        )
        .unwrap();
        let input = CellProfilerInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/cellprofiler/bin/python");
        assert_eq!(
            input.extra_args,
            vec![
                "--first-image-set".to_string(),
                "1".to_string(),
                "--last-image-set".to_string(),
                "100".to_string(),
            ]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        // Without `[bio.cellprofiler]` we have nothing to work
        // with — surface the missing-section error up front rather
        // than letting prepare() fail with a cryptic key lookup.
        let d = tempdir("cellprofiler-nosec");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = CellProfilerInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.cellprofiler]"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
