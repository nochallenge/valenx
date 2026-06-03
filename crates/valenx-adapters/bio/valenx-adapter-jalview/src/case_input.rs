//! `[bio.jalview]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "jalview.view"
//!
//! [bio.jalview]
//! jar             = "/opt/jalview/jalview.jar"
//! input           = "alignment.fa"
//! output_basename = "image"
//! output_format   = "png"     # optional, defaults to "png"
//! extra_args      = []         # optional, defaults to []
//! ```
//!
//! Jalview is the Barton group's Java-based multiple sequence
//! alignment viewer. Its headless mode (`-nodisplay`) consumes
//! an alignment file and emits an image (PNG / SVG), HTML
//! rendering, or a re-formatted alignment in any supported
//! format (FASTA / Clustal / etc.).
//!
//! Jalview ships as a Java JAR (`jalview.jar`) — there's no
//! `jalview` launcher on PATH for headless invocation. The user
//! supplies the jar path via `jar`; we probe that `java` itself
//! is installed and invoke `java -jar <jar> -nodisplay ...`
//! from `prepare()`.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct JalviewInput {
    /// Absolute path to the `jalview.jar` distributed by the
    /// Barton group.
    pub jar: PathBuf,
    /// Path to the input alignment (FASTA / Clustal / Stockholm /
    /// any format Jalview can read).
    pub input: PathBuf,
    /// Filename stem Jalview uses to label the output
    /// (`<basename>.<ext>` where `<ext>` is derived from
    /// `output_format`).
    pub output_basename: String,
    /// Output format flag passed to Jalview (`png`, `html`, `svg`,
    /// `fasta`, `clustal`, ...). Defaults to `"png"`.
    pub output_format: String,
    /// Additional CLI arguments appended to the `java -jar` call —
    /// useful for `-colour <scheme>` or `-features <file>` overrides.
    pub extra_args: Vec<String>,
}

impl JalviewInput {
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
            .and_then(|v| v.get("jalview"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.jalview] section",
                    case_toml.display()
                ))
            })?;

        let jar = block.get("jar").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.jalview].jar required (path to jalview.jar)"
            ))
        })?;
        if jar.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.jalview].jar must not be empty"
            )));
        }

        let input = block
            .get("input")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.jalview].input required")))?;
        if input.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.jalview].input must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.jalview].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.jalview].output_basename must not be empty"
            )));
        }

        // `output_format` is optional — defaults to `"png"` (the
        // most common headless use case is batch image export).
        let output_format = match block.get("output_format") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.jalview].output_format must be a string"
                    ))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.jalview].output_format must not be empty"
                    )));
                }
                s.to_string()
            }
            None => "png".to_string(),
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.jalview].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.jalview].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            jar: PathBuf::from(jar),
            input: PathBuf::from(input),
            output_basename: output_basename.to_string(),
            output_format,
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_with_default_format() {
        // No `output_format` → defaults to `"png"`. This is the
        // most common headless use case (batch image export of
        // alignments).
        let d = tempdir("jalview-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "jalview.view"

[bio.jalview]
jar             = "/opt/jalview/jalview.jar"
input           = "alignment.fa"
output_basename = "image"
"#,
        )
        .unwrap();
        let input = JalviewInput::from_case_dir(&d).unwrap();
        assert_eq!(input.jar, PathBuf::from("/opt/jalview/jalview.jar"));
        assert_eq!(input.input, PathBuf::from("alignment.fa"));
        assert_eq!(input.output_basename, "image");
        assert_eq!(input.output_format, "png");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_jar() {
        // Without the JAR path, we can't compose a `java -jar`
        // invocation. Reject up front so the user catches the typo
        // at validation time rather than after `java` spins up only
        // to error on a missing jar.
        let d = tempdir("jalview-nojar");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "jalview.view"

[bio.jalview]
jar             = ""
input           = "alignment.fa"
output_basename = "image"
"#,
        )
        .unwrap();
        let err = JalviewInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("jar"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_input() {
        // The alignment input is what Jalview renders — without it
        // there's nothing to view. Reject up front.
        let d = tempdir("jalview-noinput");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "jalview.view"

[bio.jalview]
jar             = "/opt/jalview/jalview.jar"
input           = ""
output_basename = "image"
"#,
        )
        .unwrap();
        let err = JalviewInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("input"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
