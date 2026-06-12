//! `[bio.dnachisel]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "dnachisel.optimize"
//!
//! [bio.dnachisel]
//! script           = "optimize.py"
//! python           = "python3"           # optional, defaults to python3
//! # input_fasta    = "protein.fa"         # optional, omit when generating from scratch
//! output_basename  = "optimized"
//! ```
//!
//! DNA Chisel is the Edinburgh Genome Foundry's Python codon-
//! optimization / sequence-design library — given a target protein or
//! a degenerate template, it solves the constraint-satisfaction
//! problem of choosing codons that satisfy GC-content windows, avoid
//! restriction sites / homopolymers, match a host's codon-usage
//! table, etc. The adapter itself doesn't generate Python; the user
//! authors an `optimize.py` that does
//! `from dnachisel import DnaOptimizationProblem, ...` (or similar)
//! and the actual constraint setup. We just spawn `python <script>`
//! after staging the script (and any optional `.fa` / `.fasta` input)
//! into the workdir.
//!
//! `input_fasta` is optional: omit it for scripts that generate
//! sequences from scratch (e.g. from a protein string baked into the
//! script), or supply a path to an existing `.fa` / `.fasta` file the
//! script reads as its template.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct DnaChiselInput {
    /// Path to the user-authored Python driver script (relative to
    /// the case directory). Must end in `.py` (case-insensitive).
    pub script: PathBuf,
    /// Python interpreter binary name / path. Defaults to `python3`
    /// so the adapter walks PATH; can be pinned to an absolute path
    /// for users with multiple Python installs / venvs.
    pub python: String,
    /// Optional path to an input FASTA the script reads as a starting
    /// template. `None` means the script generates / synthesises its
    /// own input sequences (e.g. from a protein string in the
    /// script).
    pub input_fasta: Option<PathBuf>,
    /// Filename stem for outputs. The script writes
    /// `<basename>*.fasta` (FASTA), `<basename>*.gb` / `.genbank`
    /// (GenBank with annotated constraints), `<basename>*.json`
    /// (constraint reports), and `<basename>*.png` (optional
    /// objective-curve plots) into the workdir.
    pub output_basename: String,
}

impl DnaChiselInput {
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
            .and_then(|v| v.get("dnachisel"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.dnachisel] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.dnachisel].script required"))
            })?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.dnachisel].script must not be empty"
            )));
        }
        // Enforce a `.py` extension (case-insensitive). Python
        // tolerates other extensions but `import dnachisel` workflows
        // are conventionally `.py`; flagging this up front saves a
        // confusing runtime error from the interpreter.
        let ext_ok = std::path::Path::new(script)
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("py"))
            .unwrap_or(false);
        if !ext_ok {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.dnachisel].script `{script}` must end in `.py`"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let input_fasta = match block.get("input_fasta").and_then(|v| v.as_str()) {
            Some(s) if !s.is_empty() => Some(PathBuf::from(s)),
            _ => None,
        };

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.dnachisel].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.dnachisel].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
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
        let d = tempdir("dnachisel-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "dnachisel.optimize"

[bio.dnachisel]
script          = "optimize.py"
output_basename = "optimized"
"#,
        )
        .unwrap();
        let input = DnaChiselInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("optimize.py"));
        assert_eq!(input.python, "python3");
        // No input_fasta — script generates its own sequences (e.g.
        // from a protein string baked in).
        assert_eq!(input.input_fasta, None);
        assert_eq!(input.output_basename, "optimized");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_input_fasta() {
        // Pinned conda interpreter + an existing FASTA template the
        // script optimizes (e.g. a protein -> codon-optimized DNA
        // pipeline).
        let d = tempdir("dnachisel-input");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "dnachisel.optimize"

[bio.dnachisel]
script          = "optimize.py"
python          = "/opt/conda/envs/codon/bin/python"
input_fasta     = "protein.fa"
output_basename = "spike-optimized"
"#,
        )
        .unwrap();
        let input = DnaChiselInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/codon/bin/python");
        assert_eq!(input.input_fasta, Some(PathBuf::from("protein.fa")));
        assert_eq!(input.output_basename, "spike-optimized");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_non_py_extension() {
        // Wrong extension is the most common typo (`.fasta`, `.fa`
        // from a copy-paste off the input field); catch it at parse
        // time so the user gets a clear error before Python is
        // invoked.
        let d = tempdir("dnachisel-badext");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "dnachisel.optimize"

[bio.dnachisel]
script          = "optimize.fa"
output_basename = "optimized"
"#,
        )
        .unwrap();
        let err = DnaChiselInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains(".py"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
