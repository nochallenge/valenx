//! `[bio.beast2]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "beast2.mcmc"
//!
//! [bio.beast2]
//! xml        = "model.xml"
//! seed       = 42                # optional — passed via -seed when present
//! threads    = 1                 # optional, defaults to 1
//! overwrite  = false             # optional, defaults to false
//! extra_args = ["-prefix", "r1_"]  # optional, defaults to []
//! ```
//!
//! BEAST 2's primary input is a BEAUti-generated XML model file
//! describing the phylogenetic prior, partition data, clock model,
//! tree prior, and operator schedule. The adapter doesn't generate
//! the XML; the user authors / generates it (typically through
//! BEAUti) and references it from `xml`.
//!
//! `seed` maps to `beast -seed <N>`, used to make a stochastic MCMC
//! run reproducible. `threads` maps to `beast -threads N` (BEAST 2
//! parallelises tree-likelihood evaluation across threads).
//! `overwrite = true` adds `-overwrite` so an existing output set
//! from a previous run is replaced rather than triggering a fail.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct Beast2Input {
    /// Path to the BEAST 2 XML model file (relative to the case
    /// directory, or absolute).
    pub xml: PathBuf,
    /// Optional PRNG seed surfaced via `beast -seed <N>`. None lets
    /// BEAST pick its own seed (printed on the run banner).
    pub seed: Option<u64>,
    /// Number of threads BEAST should use for tree-likelihood
    /// evaluation. Maps to `-threads N`. Defaults to 1.
    pub threads: u32,
    /// When true, append `-overwrite` so existing output files
    /// from a prior run are replaced. Defaults to false (BEAST's
    /// own default — fail rather than overwrite).
    pub overwrite: bool,
    /// Additional CLI arguments appended to the invocation.
    pub extra_args: Vec<String>,
}

impl Beast2Input {
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
            .and_then(|v| v.get("beast2"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.beast2] section",
                    case_toml.display()
                ))
            })?;

        let xml = block
            .get("xml")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.beast2].xml required")))?;
        if xml.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.beast2].xml must not be empty"
            )));
        }

        let seed = match block.get("seed") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.beast2].seed must be an integer"))
                })?;
                if raw < 0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.beast2].seed must be non-negative, got {raw}"
                    )));
                }
                Some(raw as u64)
            }
            None => None,
        };

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.beast2].threads must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.beast2].threads must be >= 1, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.beast2].threads `{raw}` exceeds u32::MAX"
                    )));
                }
                raw as u32
            }
            None => 1,
        };

        let overwrite = block
            .get("overwrite")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.beast2].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.beast2].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            xml: PathBuf::from(xml),
            seed,
            threads,
            overwrite,
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
        let d = tempdir("beast2-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "beast2.mcmc"

[bio.beast2]
xml = "model.xml"
"#,
        )
        .unwrap();
        let input = Beast2Input::from_case_dir(&d).unwrap();
        assert_eq!(input.xml, PathBuf::from("model.xml"));
        assert_eq!(input.seed, None);
        // Defaults: 1 thread, no overwrite, no extras.
        assert_eq!(input.threads, 1);
        assert!(!input.overwrite);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_seed_and_overwrite() {
        // Reproducible run with overwrite enabled and explicit
        // multi-threading. Mirrors a typical "redo this MCMC chain
        // with the same seed for diagnostics" workflow.
        let d = tempdir("beast2-seed");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "beast2.mcmc"

[bio.beast2]
xml        = "model.xml"
seed       = 1234567890
threads    = 8
overwrite  = true
extra_args = ["-prefix", "r1_"]
"#,
        )
        .unwrap();
        let input = Beast2Input::from_case_dir(&d).unwrap();
        assert_eq!(input.seed, Some(1_234_567_890));
        assert_eq!(input.threads, 8);
        assert!(input.overwrite);
        assert_eq!(
            input.extra_args,
            vec!["-prefix".to_string(), "r1_".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_xml() {
        // The XML is the whole model — empty string means BEAST has
        // no work to do. Reject so the failure is fast and obvious.
        let d = tempdir("beast2-noxml");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "beast2.mcmc"

[bio.beast2]
xml = ""
"#,
        )
        .unwrap();
        let err = Beast2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("xml"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_threads() {
        // BEAST treats `-threads 0` as a hard error; we reject up
        // front so the validation failure is surfaced before BEAST
        // starts.
        let d = tempdir("beast2-zerot");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "beast2.mcmc"

[bio.beast2]
xml     = "model.xml"
threads = 0
"#,
        )
        .unwrap();
        let err = Beast2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("threads"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
