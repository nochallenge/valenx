//! `[bio.samtools]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "samtools.view"
//!
//! [bio.samtools]
//! action     = "view"             # one of "view", "sort", "index", "flagstat"
//! input      = "aligned.sam"      # SAM/BAM/CRAM input
//! output     = "aligned.bam"      # required for view/sort, ignored otherwise
//! threads    = 4                  # optional, defaults to 1
//! extra_args = ["-b"]             # optional, defaults to []
//! ```
//!
//! `action` selects which samtools subcommand the adapter wraps:
//!
//! - `view`     — convert SAM↔BAM (or count records with `-c`).
//!   Requires `output`.
//! - `sort`     — sort by coordinate. Requires `output`.
//! - `index`    — index a sorted BAM (writes `<input>.bai` next to
//!   the input). `output` ignored.
//! - `flagstat` — alignment QC summary. Output is implicitly
//!   `flagstat.txt` in the workdir; `output` ignored.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical samtools action list. Module-public so the UI can
/// surface the supported values without redefining them here.
pub const SUPPORTED_ACTIONS: &[&str] = &["view", "sort", "index", "flagstat"];

/// Returns true if the action requires an explicit `output` path.
/// `view` / `sort` produce a new SAM/BAM file at the given path.
/// `index` writes a sidecar (`.bai`) next to the input. `flagstat`
/// emits its summary to stdout (the adapter redirects it to a fixed
/// file in the workdir).
pub fn action_requires_output(action: &str) -> bool {
    matches!(action, "view" | "sort")
}

#[derive(Clone, Debug, PartialEq)]
pub struct SamtoolsInput {
    pub action: String,
    pub input: PathBuf,
    pub output: Option<PathBuf>,
    pub threads: u32,
    pub extra_args: Vec<String>,
}

impl SamtoolsInput {
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
            .and_then(|v| v.get("samtools"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.samtools] section",
                    case_toml.display()
                ))
            })?;

        let action = block
            .get("action")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.samtools].action required (one of {SUPPORTED_ACTIONS:?})"
                ))
            })?;
        if !SUPPORTED_ACTIONS.contains(&action) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.samtools].action `{action}` not recognised — \
                 expected one of {SUPPORTED_ACTIONS:?}"
            )));
        }

        let input_str = block.get("input").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.samtools].input required (path to SAM/BAM/CRAM file)"
            ))
        })?;
        if input_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.samtools].input must not be empty"
            )));
        }

        let output_raw = block.get("output").and_then(|v| v.as_str());
        if action_requires_output(action) {
            // view / sort require a valid output path.
            let s = output_raw.ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.samtools].output required for action `{action}`"
                ))
            })?;
            if s.is_empty() {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.samtools].output must not be empty"
                )));
            }
        }
        // For index / flagstat we accept an `output` key but ignore
        // it — there's no surface for it in the CLI invocation. The
        // None-when-absent shape is preserved so consumers can
        // distinguish "not set" from "empty string".
        let output = output_raw.filter(|s| !s.is_empty()).map(PathBuf::from);

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.samtools].threads must be an integer"
                    ))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.samtools].threads must be >= 1, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 1,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.samtools].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.samtools].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            action: action.to_string(),
            input: PathBuf::from(input_str),
            output,
            threads,
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_view_minimal() {
        // view with input + output. Defaults: 1 thread, no extras.
        let d = tempdir("samtools");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "samtools.view"

[bio.samtools]
action = "view"
input  = "aligned.sam"
output = "aligned.bam"
"#,
        )
        .unwrap();
        let input = SamtoolsInput::from_case_dir(&d).unwrap();
        assert_eq!(input.action, "view");
        assert_eq!(input.input, PathBuf::from("aligned.sam"));
        assert_eq!(input.output, Some(PathBuf::from("aligned.bam")));
        assert_eq!(input.threads, 1);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn view_requires_output() {
        // view without output must fail at parse time — there's no
        // sensible default and silently writing to stdout would
        // bypass the artifact-collection contract.
        let d = tempdir("samtools");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "samtools.view"

[bio.samtools]
action = "view"
input  = "aligned.sam"
"#,
        )
        .unwrap();
        let err = SamtoolsInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("output"), "msg: {msg}");
        assert!(msg.contains("view"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn flagstat_does_not_require_output() {
        // flagstat emits to stdout; the adapter pins the output
        // filename, so no explicit `output` is needed.
        let d = tempdir("samtools");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "samtools.flagstat"

[bio.samtools]
action = "flagstat"
input  = "aligned.bam"
"#,
        )
        .unwrap();
        let input = SamtoolsInput::from_case_dir(&d).unwrap();
        assert_eq!(input.action, "flagstat");
        assert_eq!(input.output, None);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_action() {
        // `samtools depth` is a real subcommand but isn't one the
        // adapter wraps — must be rejected up front.
        let d = tempdir("samtools");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "samtools.depth"

[bio.samtools]
action = "depth"
input  = "aligned.bam"
"#,
        )
        .unwrap();
        let err = SamtoolsInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("flagstat"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_full_sort_with_threads_and_extras() {
        // sort with explicit threading and `-O cram` extras.
        let d = tempdir("samtools");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "samtools.sort"

[bio.samtools]
action     = "sort"
input      = "aligned.bam"
output     = "sorted.bam"
threads    = 8
extra_args = ["-O", "bam", "-l", "9"]
"#,
        )
        .unwrap();
        let input = SamtoolsInput::from_case_dir(&d).unwrap();
        assert_eq!(input.action, "sort");
        assert_eq!(input.threads, 8);
        assert_eq!(
            input.extra_args,
            vec![
                "-O".to_string(),
                "bam".to_string(),
                "-l".to_string(),
                "9".to_string()
            ]
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
