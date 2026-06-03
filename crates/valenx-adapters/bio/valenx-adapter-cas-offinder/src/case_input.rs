//! `[bio.cas_offinder]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "cas-offinder.search"
//!
//! [bio.cas_offinder]
//! input      = "query.in"
//! output     = "hits.tsv"
//! backend    = "C"               # "C" (CPU) | "G" (GPU) | "A" (Auto)
//! extra_args = []                # optional, defaults to []
//! ```
//!
//! Cas-OFFinder is a CLI off-target scanner: given a list of guide
//! sequences + PAM patterns + mismatch budget in a plain-text input
//! file, it walks a reference genome and reports every position whose
//! sequence matches one of the guides within the configured Hamming
//! distance. The CLI shape is fixed:
//!
//! ```sh
//! cas-offinder <input> {C|G|A} <output>
//! ```
//!
//! The `backend` letter selects between OpenCL device classes — `C`
//! pins to CPU, `G` to GPU, `A` lets Cas-OFFinder pick the fastest
//! available device at startup.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct CasOffinderInput {
    pub input: PathBuf,
    pub output: PathBuf,
    pub backend: String,
    pub extra_args: Vec<String>,
}

/// Recognised `backend` values. Cas-OFFinder's CLI accepts exactly one
/// of these single-letter codes — `C` for CPU, `G` for GPU, `A` for
/// "auto, pick the fastest device available at runtime".
pub const BACKENDS: &[&str] = &["C", "G", "A"];

impl CasOffinderInput {
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
        // Try both spellings — `bio.cas_offinder` (TOML-friendly,
        // underscore) and `bio.cas-offinder` (matches the binary
        // name). TOML allows hyphens in bare keys but most users
        // prefer the underscore form for IDE / lint friendliness.
        let block = parsed
            .get("bio")
            .and_then(|v| v.get("cas_offinder").or_else(|| v.get("cas-offinder")))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.cas_offinder] section",
                    case_toml.display()
                ))
            })?;

        let input_str = block.get("input").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!("[bio.cas_offinder].input required"))
        })?;
        if input_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cas_offinder].input must not be empty"
            )));
        }

        let output_str = block
            .get("output")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.cas_offinder].output required"))
            })?;
        if output_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cas_offinder].output must not be empty"
            )));
        }

        let backend = block
            .get("backend")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.cas_offinder].backend required (one of {BACKENDS:?})"
                ))
            })?;
        if !BACKENDS.contains(&backend) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cas_offinder].backend `{backend}` not recognised; \
                 must be one of {BACKENDS:?}"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.cas_offinder].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.cas_offinder].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            input: PathBuf::from(input_str),
            output: PathBuf::from(output_str),
            backend: backend.to_string(),
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_with_cpu_backend() {
        let d = tempdir("cas-offinder-cpu");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cas-offinder.search"

[bio.cas_offinder]
input   = "query.in"
output  = "hits.tsv"
backend = "C"
"#,
        )
        .unwrap();
        let input = CasOffinderInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("query.in"));
        assert_eq!(input.output, PathBuf::from("hits.tsv"));
        assert_eq!(input.backend, "C");
        // Default: no extra CLI flags appended.
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_gpu_backend() {
        // GPU mode is the typical production setting — Cas-OFFinder's
        // OpenCL path scales linearly with available compute units.
        let d = tempdir("cas-offinder-gpu");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cas-offinder.search"

[bio.cas_offinder]
input      = "query.in"
output     = "hits.tsv"
backend    = "G"
extra_args = ["--quiet"]
"#,
        )
        .unwrap();
        let input = CasOffinderInput::from_case_dir(&d).unwrap();
        assert_eq!(input.backend, "G");
        assert_eq!(input.extra_args, vec!["--quiet".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_auto_backend() {
        // `A` lets Cas-OFFinder pick the fastest available device at
        // startup. Useful for portable case definitions that should
        // run on either GPU-equipped workstations or CPU-only HPC
        // queues without per-host config.
        let d = tempdir("cas-offinder-auto");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cas-offinder.search"

[bio.cas_offinder]
input   = "query.in"
output  = "hits.tsv"
backend = "A"
"#,
        )
        .unwrap();
        let input = CasOffinderInput::from_case_dir(&d).unwrap();
        assert_eq!(input.backend, "A");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_backend() {
        // Cas-OFFinder's CLI only accepts C / G / A — `X` would be
        // a typo. Reject up front so the user catches it before a
        // multi-second binary spawn.
        let d = tempdir("cas-offinder-badback");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cas-offinder.search"

[bio.cas_offinder]
input   = "query.in"
output  = "hits.tsv"
backend = "X"
"#,
        )
        .unwrap();
        let err = CasOffinderInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("backend"), "msg: {msg}");
        assert!(msg.contains('X'), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
