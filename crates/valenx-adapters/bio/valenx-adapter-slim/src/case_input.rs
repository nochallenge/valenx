//! `[bio.slim]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "slim.simulate"
//!
//! [bio.slim]
//! script          = "model.slim"
//! seed            = 42                 # optional — passed via -s when present
//! output_basename = "sim"
//! extra_args      = ["-d", "MU=1e-7"]  # optional, defaults to []
//! ```
//!
//! SLiM is a forward-time population-genetics simulator driven by
//! Eidos scripts (its own embedded scripting language, syntactically
//! close to R). The adapter doesn't generate Eidos; the user authors
//! a `.slim` model file referenced from `script`. Output paths are
//! determined by the script itself (typically via `writeFile()` /
//! `treeSeqOutput()` calls), so the adapter only needs to know the
//! basename for `collect()`-time labelling — it doesn't try to
//! predict the exact output filenames the script will write.
//!
//! `seed` maps to `slim -s <N>`, used to make a stochastic run
//! reproducible. Optional — when omitted, SLiM picks its own random
//! seed and prints it on the run banner.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct SlimInput {
    /// Path to the Eidos `.slim` model file (relative to the case
    /// directory, or absolute).
    pub script: PathBuf,
    /// Optional PRNG seed surfaced via `slim -s <N>`. None lets
    /// SLiM pick its own seed (printed on the run banner).
    pub seed: Option<u64>,
    /// Filename stem the user's script uses for outputs. Surfaced
    /// here so `collect()` can label artefacts uniformly even
    /// though SLiM scripts choose their own output paths.
    pub output_basename: String,
    /// Additional CLI arguments appended after the script path.
    /// `-d KEY=VALUE` is the canonical way to inject Eidos
    /// constants from outside the script.
    pub extra_args: Vec<String>,
}

impl SlimInput {
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
            .and_then(|v| v.get("slim"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.slim] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.slim].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.slim].script must not be empty"
            )));
        }

        let seed = match block.get("seed") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.slim].seed must be an integer"))
                })?;
                if raw < 0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.slim].seed must be non-negative, got {raw}"
                    )));
                }
                Some(raw as u64)
            }
            None => None,
        };

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.slim].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.slim].output_basename must not be empty"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.slim].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.slim].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            script: PathBuf::from(script),
            seed,
            output_basename: output_basename.to_string(),
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("slim-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "slim.simulate"

[bio.slim]
script          = "model.slim"
output_basename = "sim"
"#,
        )
        .unwrap();
        let input = SlimInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("model.slim"));
        assert_eq!(input.seed, None);
        assert_eq!(input.output_basename, "sim");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_seed() {
        // Reproducible run — pinning a seed plus an Eidos constant
        // injection via `-d MU=...` exercises both override paths.
        let d = tempdir("slim-seed");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "slim.simulate"

[bio.slim]
script          = "neutral.slim"
seed            = 1234567890
output_basename = "neutral"
extra_args      = ["-d", "MU=1e-7", "-d", "L=1000000"]
"#,
        )
        .unwrap();
        let input = SlimInput::from_case_dir(&d).unwrap();
        assert_eq!(input.seed, Some(1_234_567_890));
        assert_eq!(input.output_basename, "neutral");
        assert_eq!(
            input.extra_args,
            vec![
                "-d".to_string(),
                "MU=1e-7".to_string(),
                "-d".to_string(),
                "L=1000000".to_string(),
            ]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_script() {
        // The Eidos script drives the entire simulation — empty
        // string means SLiM has no model to run. Reject so the
        // failure is fast and obvious.
        let d = tempdir("slim-noscript");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "slim.simulate"

[bio.slim]
script          = ""
output_basename = "sim"
"#,
        )
        .unwrap();
        let err = SlimInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("script"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_basename() {
        // Output basename anchors collect()'s artefact labels;
        // empty string would leave the user with unlabelled
        // artefacts. Reject up front.
        let d = tempdir("slim-nobase");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "slim.simulate"

[bio.slim]
script          = "model.slim"
output_basename = ""
"#,
        )
        .unwrap();
        let err = SlimInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("output_basename"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
