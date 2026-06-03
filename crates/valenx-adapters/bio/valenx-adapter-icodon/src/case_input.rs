//! `[bio.icodon]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "icodon.predict"
//!
//! [bio.icodon]
//! script           = "predict.R"
//! rscript          = "Rscript"           # optional, defaults to "Rscript"
//! # input_fasta    = "transcripts.fa"    # optional, omit to skip
//! output_basename  = "stability"
//! ```
//!
//! iCodon (Vejnar lab) is an R package that predicts codon-level
//! mRNA stability from sequence. The adapter itself doesn't generate
//! R; the user supplies a `predict.R` that does
//! `library(iCodon)` and the actual prediction work. We just spawn
//! `Rscript <predict.R>` after staging the script (and any optional
//! input FASTA) into the workdir.
//!
//! `input_fasta` is optional: omit it for scripts that fetch their
//! own data, or supply a path to a FASTA file the script reads via
//! the staged filename.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct IcodonInput {
    /// Path to the user-authored R driver script (relative to the
    /// case directory, or absolute). Must end in `.R`
    /// (case-insensitive).
    pub script: PathBuf,
    /// Rscript binary name / path. Defaults to `Rscript` so the
    /// adapter walks PATH; can be pinned to an absolute path for
    /// users with multiple R installs.
    pub rscript: String,
    /// Optional path to an input FASTA the script reads as its
    /// starting point. `None` means the script fetches /
    /// synthesises its own data.
    pub input_fasta: Option<PathBuf>,
    /// Filename stem for outputs. The script writes
    /// `<basename>*.csv` / `<basename>*.tsv` (stability tables),
    /// `<basename>*.rds` (R objects), and `<basename>*.png` (plots)
    /// into the workdir.
    pub output_basename: String,
}

impl IcodonInput {
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
            .and_then(|v| v.get("icodon"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.icodon] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.icodon].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.icodon].script must not be empty"
            )));
        }
        // Enforce a `.R` extension (case-insensitive). The R
        // ecosystem doesn't have a hard rule but `Rscript` itself
        // expects `.R` / `.r`; flagging this up front saves a
        // confusing runtime error from R.
        let ext_ok = std::path::Path::new(script)
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("R"))
            .unwrap_or(false);
        if !ext_ok {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.icodon].script `{script}` must end in `.R`"
            )));
        }

        let rscript = block
            .get("rscript")
            .and_then(|v| v.as_str())
            .unwrap_or("Rscript")
            .to_string();

        let input_fasta = match block.get("input_fasta").and_then(|v| v.as_str()) {
            Some(s) if !s.is_empty() => Some(PathBuf::from(s)),
            _ => None,
        };

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.icodon].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.icodon].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            rscript,
            input_fasta,
            output_basename: output_basename.to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("icodon-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "icodon.predict"

[bio.icodon]
script          = "predict.R"
output_basename = "stability"
"#,
        )
        .unwrap();
        let input = IcodonInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("predict.R"));
        assert_eq!(input.rscript, "Rscript");
        // No input_fasta — script fetches / synthesises its own.
        assert_eq!(input.input_fasta, None);
        assert_eq!(input.output_basename, "stability");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_input_fasta() {
        // Pinned Rscript path + an input FASTA the script reads.
        let d = tempdir("icodon-fasta");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "icodon.predict"

[bio.icodon]
script          = "predict.R"
rscript         = "/opt/R/4.3.0/bin/Rscript"
input_fasta     = "transcripts.fa"
output_basename = "stability"
"#,
        )
        .unwrap();
        let input = IcodonInput::from_case_dir(&d).unwrap();
        assert_eq!(input.rscript, "/opt/R/4.3.0/bin/Rscript");
        assert_eq!(input.input_fasta, Some(PathBuf::from("transcripts.fa")));
        assert_eq!(input.output_basename, "stability");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_non_r_extension() {
        // Wrong extension is the most common typo (`.r` works,
        // `.py` / `.txt` do not). Catch it at parse time so the
        // user gets a clear error before Rscript is invoked.
        let d = tempdir("icodon-badext");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "icodon.predict"

[bio.icodon]
script          = "predict.py"
output_basename = "stability"
"#,
        )
        .unwrap();
        let err = IcodonInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains(".R"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
