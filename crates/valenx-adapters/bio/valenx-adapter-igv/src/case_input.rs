//! `[bio.igv]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "igv.index"
//!
//! [bio.igv]
//! action      = "index"           # one of "index", "count", "sort", "tile"
//! input       = "aligned.bam"     # BAM/VCF/SAM input
//! output      = "aligned.tdf"     # required for count/sort/tile, ignored for index
//! window_size = 25                # optional, defaults to 25 (used by `count`)
//! extra_args  = ["-z", "5"]       # optional, defaults to []
//! ```
//!
//! `action` selects which igvtools subcommand the adapter wraps:
//!
//! - `index` — write a `.bai` (BAM) or `.idx` (VCF) sidecar next to
//!   the input. `output` is ignored (igvtools writes the sidecar
//!   next to the input, not in the workdir).
//! - `count` — generate a `.tdf` density file from a BAM. Requires
//!   `output`. Honours `window_size` (the `-w` flag).
//! - `sort`  — coordinate-sort a SAM/BAM/VCF. Requires `output`.
//! - `tile`  — generate a `.tdf` tile from a coverage track.
//!   Requires `output`.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical igvtools action list. Module-public so the UI can
/// surface the supported values without redefining them here.
pub const SUPPORTED_ACTIONS: &[&str] = &["index", "count", "sort", "tile"];

/// Returns true if the action requires an explicit `output` path.
/// `count` / `sort` / `tile` produce a new file at the given path;
/// `index` writes a sidecar next to the input and ignores `output`.
pub fn action_requires_output(action: &str) -> bool {
    matches!(action, "count" | "sort" | "tile")
}

#[derive(Clone, Debug, PartialEq)]
pub struct IgvInput {
    /// One of `index`, `count`, `sort`, `tile`.
    pub action: String,
    /// Input BAM/VCF/SAM path.
    pub input: PathBuf,
    /// Output path for count / sort / tile. None for index.
    pub output: Option<PathBuf>,
    /// Window size for the `count` action. Defaults to 25 (igvtools'
    /// own default). Must be >= 1.
    pub window_size: u32,
    /// Additional CLI arguments appended after the action's
    /// canonical positional / flag set.
    pub extra_args: Vec<String>,
}

impl IgvInput {
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
            .and_then(|v| v.get("igv"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.igv] section",
                    case_toml.display()
                ))
            })?;

        let action = block
            .get("action")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.igv].action required (one of {SUPPORTED_ACTIONS:?})"
                ))
            })?;
        if !SUPPORTED_ACTIONS.contains(&action) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.igv].action `{action}` not recognised — \
                 expected one of {SUPPORTED_ACTIONS:?}"
            )));
        }

        let input_str = block.get("input").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.igv].input required (path to BAM/VCF/SAM file)"
            ))
        })?;
        if input_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.igv].input must not be empty"
            )));
        }

        let output_raw = block.get("output").and_then(|v| v.as_str());
        if action_requires_output(action) {
            // count / sort / tile require a valid output path.
            let s = output_raw.ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.igv].output required for action `{action}`"
                ))
            })?;
            if s.is_empty() {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "[bio.igv].output must not be empty"
                )));
            }
        }
        // For `index` we ignore any provided `output` — igvtools
        // writes the sidecar next to the input. Preserve the
        // None-when-absent shape so consumers can distinguish "not
        // set" from "empty string".
        let output = if action_requires_output(action) {
            output_raw.filter(|s| !s.is_empty()).map(PathBuf::from)
        } else {
            None
        };

        let window_size = match block.get("window_size") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.igv].window_size must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.igv].window_size must be >= 1, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 25,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.igv].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.igv].extra_args entries must be strings"
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
            window_size,
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_index_minimal() {
        // index: input only, no output. window_size defaults to 25
        // even though `index` doesn't use it — we still parse the
        // default through.
        let d = tempdir("igv");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "igv.index"

[bio.igv]
action = "index"
input  = "aligned.bam"
"#,
        )
        .unwrap();
        let input = IgvInput::from_case_dir(&d).unwrap();
        assert_eq!(input.action, "index");
        assert_eq!(input.input, PathBuf::from("aligned.bam"));
        assert_eq!(input.output, None, "index ignores output");
        assert_eq!(input.window_size, 25);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_count_with_window() {
        // count: explicit window_size override + extras.
        let d = tempdir("igv");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "igv.count"

[bio.igv]
action      = "count"
input       = "aligned.bam"
output      = "aligned.tdf"
window_size = 100
extra_args  = ["-z", "5"]
"#,
        )
        .unwrap();
        let input = IgvInput::from_case_dir(&d).unwrap();
        assert_eq!(input.action, "count");
        assert_eq!(input.window_size, 100);
        assert_eq!(input.output, Some(PathBuf::from("aligned.tdf")));
        assert_eq!(input.extra_args, vec!["-z".to_string(), "5".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn sort_requires_output() {
        // sort without output must fail at parse time — there's no
        // sensible default and silently writing back over the input
        // would be unsafe.
        let d = tempdir("igv");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "igv.sort"

[bio.igv]
action = "sort"
input  = "aligned.bam"
"#,
        )
        .unwrap();
        let err = IgvInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("output"), "msg: {msg}");
        assert!(msg.contains("sort"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_action() {
        // `igvtools toTDF` is a real subcommand alias but we
        // wrap a fixed set; reject anything outside the list.
        let d = tempdir("igv");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "igv.totdf"

[bio.igv]
action = "totdf"
input  = "aligned.bam"
output = "aligned.tdf"
"#,
        )
        .unwrap();
        let err = IgvInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("tile"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_window() {
        // window_size must be >= 1; zero / negative are nonsensical
        // for bin-counting and igvtools would error at runtime.
        // Catch it at parse time.
        let d = tempdir("igv");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "igv.count"

[bio.igv]
action      = "count"
input       = "aligned.bam"
output      = "aligned.tdf"
window_size = 0
"#,
        )
        .unwrap();
        let err = IgvInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("window_size"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
