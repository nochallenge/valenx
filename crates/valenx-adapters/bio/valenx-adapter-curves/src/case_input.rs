//! `[bio.curves]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "curves.analyze"
//!
//! [bio.curves]
//! input_pdb       = "structure.pdb"
//! output_basename = "analysis"
//! first_residue   = 1
//! last_residue    = 12
//! extra_args      = []           # optional, defaults to []
//! ```
//!
//! Curves+ takes its parameters via stdin (an `&inp ... &end`
//! Fortran-namelist-style block followed by residue-range cards).
//! The adapter authors that block at `prepare()` time and pipes it
//! into `Cur+`'s stdin at `run()` time. The user only surfaces the
//! input PDB, the residue range to analyse, and the basename
//! Curves+ should use for output filenames (`<basename>.lis`,
//! `<basename>.cda`, etc.).
//!
//! `first_residue` / `last_residue` define the inclusive residue
//! range Curves+ should analyse on the input strand; the adapter
//! computes the strand length and writes the appropriate strand /
//! axis cards. `first_residue <= last_residue` is enforced — a
//! reverse range is rejected up front.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct CurvesInput {
    /// Path to the input `.pdb` file (relative to the case
    /// directory, or absolute).
    pub input_pdb: PathBuf,
    /// Filename stem Curves+ uses for outputs
    /// (`<basename>.lis`, `<basename>.cda`, etc.).
    pub output_basename: String,
    /// First (inclusive) residue index in the strand to analyse.
    pub first_residue: u32,
    /// Last (inclusive) residue index in the strand to analyse.
    /// Must be >= `first_residue`.
    pub last_residue: u32,
    /// Additional CLI arguments appended to the `Cur+` invocation.
    pub extra_args: Vec<String>,
}

impl CurvesInput {
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
            .and_then(|v| v.get("curves"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.curves] section",
                    case_toml.display()
                ))
            })?;

        let input_pdb = block
            .get("input_pdb")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.curves].input_pdb required"))
            })?;
        if input_pdb.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.curves].input_pdb must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.curves].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.curves].output_basename must not be empty"
            )));
        }

        let first_residue = parse_residue(block, "first_residue")?;
        let last_residue = parse_residue(block, "last_residue")?;
        if first_residue > last_residue {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.curves].first_residue ({first_residue}) must be \
                 <= last_residue ({last_residue})"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.curves].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.curves].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            input_pdb: PathBuf::from(input_pdb),
            output_basename: output_basename.to_string(),
            first_residue,
            last_residue,
            extra_args,
        })
    }
}

/// Extract a non-negative integer residue index by key. The Curves+
/// stdin schema expects 1-based residue indices, so 0 is allowed
/// in case the user wants to bias-grab from the start of a strand,
/// but negative values are rejected as a typo.
fn parse_residue(block: &toml::Value, key: &str) -> Result<u32, AdapterError> {
    let v = block
        .get(key)
        .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.curves].{key} required")))?;
    let raw = v.as_integer().ok_or_else(|| {
        AdapterError::Other(anyhow::anyhow!("[bio.curves].{key} must be an integer"))
    })?;
    if raw < 0 {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "[bio.curves].{key} must be non-negative, got {raw}"
        )));
    }
    if raw > u32::MAX as i64 {
        return Err(AdapterError::Other(anyhow::anyhow!(
            "[bio.curves].{key} `{raw}` exceeds u32::MAX"
        )));
    }
    Ok(raw as u32)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("curves-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "curves.analyze"

[bio.curves]
input_pdb       = "structure.pdb"
output_basename = "analysis"
first_residue   = 1
last_residue    = 12
"#,
        )
        .unwrap();
        let input = CurvesInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input_pdb, PathBuf::from("structure.pdb"));
        assert_eq!(input.output_basename, "analysis");
        assert_eq!(input.first_residue, 1);
        assert_eq!(input.last_residue, 12);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_first_after_last() {
        // A reverse range (first > last) leaves Curves+ with a
        // negative-length strand. Reject up front so the failure is
        // fast and obvious.
        let d = tempdir("curves-reverse");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "curves.analyze"

[bio.curves]
input_pdb       = "structure.pdb"
output_basename = "analysis"
first_residue   = 20
last_residue    = 5
"#,
        )
        .unwrap();
        let err = CurvesInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("first_residue") && msg.contains("last_residue"),
            "msg: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_input_pdb() {
        // The PDB is the entire input — empty string means Cur+
        // has no structure to work on. Reject up front.
        let d = tempdir("curves-nopdb");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "curves.analyze"

[bio.curves]
input_pdb       = ""
output_basename = "analysis"
first_residue   = 1
last_residue    = 12
"#,
        )
        .unwrap();
        let err = CurvesInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("input_pdb"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_basename() {
        // Output basename anchors collect()'s artefact labels;
        // empty string would leave the user with unlabelled
        // artefacts. Reject up front.
        let d = tempdir("curves-nobase");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "curves.analyze"

[bio.curves]
input_pdb       = "structure.pdb"
output_basename = ""
first_residue   = 1
last_residue    = 12
"#,
        )
        .unwrap();
        let err = CurvesInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("output_basename"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
