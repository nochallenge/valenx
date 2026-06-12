//! `[bio.alphafold2]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "alphafold2.predict"
//!
//! [bio.alphafold2]
//! run_script         = "/path/to/alphafold/run_alphafold.py"
//! python             = "python3"
//! query_fasta        = "query.fasta"
//! data_dir           = "/data/alphafold/databases"
//! max_template_date  = "2022-01-01"
//! model_preset       = "monomer"        # "monomer" | "monomer_ptm" | "multimer"
//! ```
//!
//! AF2's `run_alphafold.py` has no `--num_recycles` CLI flag — recycle
//! count is baked into the model config. The schema deliberately
//! omits it to avoid the misleading impression that bumping a TOML
//! key changes solver behaviour.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct AlphaFold2Input {
    pub run_script: PathBuf,
    pub python: String,
    pub query_fasta: PathBuf,
    pub data_dir: PathBuf,
    /// ISO-8601 date (YYYY-MM-DD) used as the AF2 template cutoff.
    /// Validated to match the date shape on parse — AF2's flag will
    /// reject anything else but this lets us fail fast.
    pub max_template_date: String,
    /// One of `"monomer"` / `"monomer_ptm"` / `"multimer"`.
    pub model_preset: String,
}

/// Set of valid `model_preset` values for AlphaFold 2.
pub const MODEL_PRESETS: &[&str] = &["monomer", "monomer_ptm", "multimer"];

/// Validate an ISO-8601 date (`YYYY-MM-DD`). Strict 10-char check —
/// digits + dashes in the right slots. The validator is small enough
/// that a `regex` dep would be overkill.
pub fn is_iso_date(s: &str) -> bool {
    if s.len() != 10 {
        return false;
    }
    let bytes = s.as_bytes();
    // Positions 4 and 7 are dashes.
    if bytes[4] != b'-' || bytes[7] != b'-' {
        return false;
    }
    // Everything else must be a digit.
    for (i, b) in bytes.iter().enumerate() {
        if i == 4 || i == 7 {
            continue;
        }
        if !b.is_ascii_digit() {
            return false;
        }
    }
    true
}

impl AlphaFold2Input {
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
            .and_then(|v| v.get("alphafold2"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.alphafold2] section",
                    case_toml.display()
                ))
            })?;
        let run_script = block
            .get("run_script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.alphafold2].run_script required"))
            })?;
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        let query_fasta = block
            .get("query_fasta")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.alphafold2].query_fasta required"))
            })?;
        let data_dir = block
            .get("data_dir")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.alphafold2].data_dir required"))
            })?;
        let max_template_date = block
            .get("max_template_date")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.alphafold2].max_template_date required"
                ))
            })?;
        if !is_iso_date(max_template_date) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.alphafold2].max_template_date `{max_template_date}` \
                 must match YYYY-MM-DD"
            )));
        }
        let model_preset = block
            .get("model_preset")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.alphafold2].model_preset required"))
            })?;
        if !MODEL_PRESETS.contains(&model_preset) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.alphafold2].model_preset `{model_preset}` not recognised; \
                 must be one of {MODEL_PRESETS:?}"
            )));
        }
        Ok(Self {
            run_script: PathBuf::from(run_script),
            python,
            query_fasta: PathBuf::from(query_fasta),
            data_dir: PathBuf::from(data_dir),
            max_template_date: max_template_date.to_string(),
            model_preset: model_preset.to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_case() {
        let d = tempdir("alphafold2-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "alphafold2.predict"

[bio.alphafold2]
run_script        = "/opt/alphafold/run_alphafold.py"
query_fasta       = "query.fasta"
data_dir          = "/data/alphafold"
max_template_date = "2022-01-01"
model_preset      = "monomer"
"#,
        )
        .unwrap();
        let input = AlphaFold2Input::from_case_dir(&d).unwrap();
        assert_eq!(
            input.run_script,
            PathBuf::from("/opt/alphafold/run_alphafold.py")
        );
        assert_eq!(input.query_fasta, PathBuf::from("query.fasta"));
        assert_eq!(input.data_dir, PathBuf::from("/data/alphafold"));
        assert_eq!(input.max_template_date, "2022-01-01");
        assert_eq!(input.model_preset, "monomer");
        // Defaults.
        assert_eq!(input.python, "python3");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_preset() {
        let d = tempdir("alphafold2-badpreset");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "alphafold2.predict"

[bio.alphafold2]
run_script        = "/opt/alphafold/run_alphafold.py"
query_fasta       = "q.fasta"
data_dir          = "/d"
max_template_date = "2022-01-01"
model_preset      = "deluxe_monomer"
"#,
        )
        .unwrap();
        let err = AlphaFold2Input::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("model_preset"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_bad_date() {
        let d = tempdir("alphafold2-baddate");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "alphafold2.predict"

[bio.alphafold2]
run_script        = "/opt/alphafold/run_alphafold.py"
query_fasta       = "q.fasta"
data_dir          = "/d"
max_template_date = "yesterday"
model_preset      = "monomer"
"#,
        )
        .unwrap();
        let err = AlphaFold2Input::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("max_template_date"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("alphafold2-nosec");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = AlphaFold2Input::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.alphafold2]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn iso_date_validator_accepts_canonical_form() {
        assert!(is_iso_date("2022-01-01"));
        assert!(is_iso_date("1999-12-31"));
        assert!(is_iso_date("0000-00-00")); // Shape-only check.
    }

    #[test]
    fn iso_date_validator_rejects_malformed_strings() {
        assert!(!is_iso_date("2022-1-1"));
        assert!(!is_iso_date("2022/01/01"));
        assert!(!is_iso_date("yesterday"));
        assert!(!is_iso_date("2022-01-01T"));
        assert!(!is_iso_date(""));
    }

    /// Round-17 M2 RED→GREEN (representative): a `case.toml` larger than
    /// the project-loader cap (`MAX_PROJECT_FILE_BYTES = 1 MiB`) is
    /// refused at the bounded-read boundary instead of being slurped
    /// into memory.
    ///
    /// Pre-fix every adapter's `from_case_dir` called the unbounded
    /// `std::fs::read_to_string` variant on the case-TOML path, so a
    /// hostile multi-GB file under a derived workdir would allocate
    /// the full payload before the TOML parser saw it. The sweep
    /// replaced the raw read with
    /// `valenx_core::io_caps::read_capped_to_string` parameterised on
    /// `MAX_PROJECT_FILE_BYTES`, which rejects oversize files with an
    /// `InvalidData` IO error — this test is the representative for
    /// the 100+ adapters touched by the sweep.
    #[test]
    fn from_case_dir_rejects_oversize_case_toml() {
        use valenx_core::project::loader::MAX_PROJECT_FILE_BYTES;
        let d = tempdir("alphafold2-oversize");
        // Build a payload one byte larger than the cap. The TOML
        // content doesn't have to parse — the bounded read fires
        // BEFORE the parser gets a chance.
        let oversize = (MAX_PROJECT_FILE_BYTES as usize) + 1024;
        // Fill with valid-looking TOML padding so a future move of
        // the cap check past the parse step would still hit a parse
        // error, not a hand-rolled garbage rejection.
        let mut payload = String::with_capacity(oversize + 256);
        payload.push_str("[case]\nphysics = \"bio\"\nsolver = \"alphafold2.predict\"\n");
        payload.push_str("[bio.alphafold2]\n");
        payload.push_str("run_script        = \"/opt/x\"\n");
        payload.push_str("query_fasta       = \"q.fasta\"\n");
        payload.push_str("data_dir          = \"/d\"\n");
        payload.push_str("max_template_date = \"2022-01-01\"\n");
        payload.push_str("model_preset      = \"monomer\"\n");
        // Pad to ≥ oversize bytes with a long comment.
        payload.push_str("# ");
        while payload.len() < oversize {
            payload.push('x');
        }
        std::fs::write(d.join("case.toml"), &payload).unwrap();
        let err = AlphaFold2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        // The cap helper emits `"file ... exceeds N-byte cap"`. We
        // assert on the cap-shape so the test pins the read-capped
        // path, not a downstream parse error.
        assert!(
            msg.contains("exceeds") && msg.contains("cap"),
            "expected oversize-cap error, got: {msg}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
