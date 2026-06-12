//! `[bio.psi4]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "psi4.compute"
//!
//! [bio.psi4]
//! input      = "input.dat"
//! output     = "output.dat"
//! threads    = 4              # optional, defaults to 1
//! memory     = "4 gb"         # optional, defaults to "1 gb"
//! extra_args = ["-l", "/opt/psi4"]   # optional, defaults to []
//! ```
//!
//! Psi4 takes a Psithon input file (essentially a Python script with
//! Psi4 directives) and writes a human-readable text output to a
//! second file. `threads` maps to `psi4 -n <N>`, `memory` maps to
//! `psi4 -m <amount>` — the canonical examples in the manual use
//! `1 gb` / `2 gb` / `500 mb` style strings, so we accept either case
//! suffix and a small whole-number prefix without resorting to a
//! regex dependency.

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Cheap regex-free check: `<digits>` followed by optional whitespace
/// and either `mb` / `gb` (case-insensitive). Pinned here because
/// adding a regex crate for one validation would be wasteful, and
/// Psi4 itself only documents these two unit suffixes.
pub fn is_valid_memory(s: &str) -> bool {
    let s = s.trim();
    if s.is_empty() {
        return false;
    }
    // Walk leading digits; require at least one.
    let digit_count = s.chars().take_while(|c| c.is_ascii_digit()).count();
    if digit_count == 0 {
        return false;
    }
    let rest = s[digit_count..].trim_start();
    let suffix = rest.to_ascii_lowercase();
    matches!(suffix.as_str(), "mb" | "gb")
}

#[derive(Clone, Debug, PartialEq)]
pub struct Psi4Input {
    pub input: PathBuf,
    pub output: PathBuf,
    pub threads: u32,
    pub memory: String,
    pub extra_args: Vec<String>,
}

impl Psi4Input {
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
            .and_then(|v| v.get("psi4"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.psi4] section",
                    case_toml.display()
                ))
            })?;

        let input_str = block.get("input").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.psi4].input required (path to Psithon input file)"
            ))
        })?;
        if input_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.psi4].input must not be empty"
            )));
        }

        let output_str = block
            .get("output")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.psi4].output required (path for the Psi4 output file)"
                ))
            })?;
        if output_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.psi4].output must not be empty"
            )));
        }

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.psi4].threads must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.psi4].threads must be >= 1, got {raw}"
                    )));
                }
                if raw > u32::MAX as i64 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.psi4].threads `{raw}` exceeds u32::MAX"
                    )));
                }
                raw as u32
            }
            None => 1,
        };

        let memory = match block.get("memory") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.psi4].memory must be a string like \"4 gb\""
                    ))
                })?;
                if !is_valid_memory(s) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.psi4].memory `{s}` is not a valid Psi4 memory \
                         spec — expected `<N> mb` or `<N> gb` (e.g. \"4 gb\")"
                    )));
                }
                s.to_string()
            }
            None => "1 gb".to_string(),
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.psi4].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.psi4].extra_args entries must be strings"
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
            threads,
            memory,
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
        let d = tempdir("psi4");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "psi4.compute"

[bio.psi4]
input  = "input.dat"
output = "output.dat"
"#,
        )
        .unwrap();
        let input = Psi4Input::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("input.dat"));
        assert_eq!(input.output, PathBuf::from("output.dat"));
        assert_eq!(input.threads, 1);
        assert_eq!(input.memory, "1 gb");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_overrides() {
        // Multi-threaded run with a bumped memory ceiling and a Psi4
        // command-line knob (`-l` points at a non-default Psi4 data
        // directory — common in HPC installs).
        let d = tempdir("psi4");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "psi4.compute"

[bio.psi4]
input      = "h2o.in"
output     = "h2o.out"
threads    = 8
memory     = "16 gb"
extra_args = ["-l", "/opt/psi4/data"]
"#,
        )
        .unwrap();
        let input = Psi4Input::from_case_dir(&d).unwrap();
        assert_eq!(input.threads, 8);
        assert_eq!(input.memory, "16 gb");
        assert_eq!(
            input.extra_args,
            vec!["-l".to_string(), "/opt/psi4/data".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_threads() {
        let d = tempdir("psi4");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "psi4.compute"

[bio.psi4]
input   = "in.dat"
output  = "out.dat"
threads = 0
"#,
        )
        .unwrap();
        let err = Psi4Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("threads"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_invalid_memory() {
        // "two gigs" is not parseable by Psi4. Reject up front so the
        // user sees the failure at validation time rather than after a
        // long subprocess startup.
        let d = tempdir("psi4");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "psi4.compute"

[bio.psi4]
input  = "in.dat"
output = "out.dat"
memory = "two gigs"
"#,
        )
        .unwrap();
        let err = Psi4Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("memory"), "msg: {msg}");
        assert!(msg.contains("two gigs"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn is_valid_memory_helper() {
        // Sanity-check the regex-free helper directly.
        assert!(is_valid_memory("1 gb"));
        assert!(is_valid_memory("16gb"));
        assert!(is_valid_memory("500 MB"));
        assert!(is_valid_memory("4 GB"));
        assert!(!is_valid_memory(""));
        assert!(!is_valid_memory("gb"));
        assert!(!is_valid_memory("1 kb"));
        assert!(!is_valid_memory("two gigs"));
        assert!(!is_valid_memory("1.5 gb"));
    }
}
