//! `[bio.pydna]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "pydna.design"
//!
//! [bio.pydna]
//! script           = "design.py"
//! python           = "python3"           # optional, defaults to python3
//! # input_genbank  = "template.gb"        # optional, omit when designing from scratch
//! output_basename  = "design"
//! ```
//!
//! pydna is Bjorn Johansson's Python plasmid / clone-design library —
//! the de-facto choice for in-silico molecular cloning (PCR primer
//! design, restriction-enzyme digests, Gibson / Golden-Gate assembly,
//! homologous-recombination simulation). The adapter itself doesn't
//! generate Python; the user authors a `design.py` that does
//! `from pydna.dseqrecord import Dseqrecord` (or similar) and the
//! actual cloning logic. We just spawn `python <script>` after
//! staging the script (and any optional `.gb` input) into the
//! workdir.
//!
//! `input_genbank` is optional: omit it for scripts that fetch
//! sequences from GenBank / synthesise them from primers, or supply
//! a path to an existing `.gb` / `.genbank` file the script reads as
//! its template.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct PydnaInput {
    /// Path to the user-authored Python driver script (relative to
    /// the case directory, or absolute). Must end in `.py`
    /// (case-insensitive).
    pub script: PathBuf,
    /// Python interpreter binary name / path. Defaults to `python3`
    /// so the adapter walks PATH; can be pinned to an absolute path
    /// for users with multiple Python installs / venvs.
    pub python: String,
    /// Optional path to an input GenBank `.gb` / `.genbank` the
    /// script reads as a starting template. `None` means the script
    /// fetches / synthesises its own input sequences.
    pub input_genbank: Option<PathBuf>,
    /// Filename stem for outputs. The script writes
    /// `<basename>*.gb` / `.genbank` (GenBank), `<basename>*.fasta`
    /// (FASTA), and `<basename>*.csv` (analysis tables) into the
    /// workdir.
    pub output_basename: String,
}

impl PydnaInput {
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
            .and_then(|v| v.get("pydna"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.pydna] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.pydna].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.pydna].script must not be empty"
            )));
        }
        // Enforce a `.py` extension (case-insensitive). Python
        // tolerates other extensions but `import pydna` workflows are
        // conventionally `.py`; flagging this up front saves a
        // confusing runtime error from the interpreter.
        let ext_ok = std::path::Path::new(script)
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("py"))
            .unwrap_or(false);
        if !ext_ok {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.pydna].script `{script}` must end in `.py`"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let input_genbank = match block.get("input_genbank").and_then(|v| v.as_str()) {
            Some(s) if !s.is_empty() => Some(PathBuf::from(s)),
            _ => None,
        };

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.pydna].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.pydna].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_genbank,
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
        let d = tempdir("pydna-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pydna.design"

[bio.pydna]
script          = "design.py"
output_basename = "design"
"#,
        )
        .unwrap();
        let input = PydnaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("design.py"));
        assert_eq!(input.python, "python3");
        // No input_genbank — script fetches / synthesises its own.
        assert_eq!(input.input_genbank, None);
        assert_eq!(input.output_basename, "design");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_input_genbank() {
        // Pinned conda interpreter + an existing GenBank template the
        // script extends with new restriction sites.
        let d = tempdir("pydna-input");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pydna.design"

[bio.pydna]
script          = "extend.py"
python          = "/opt/conda/envs/cloning/bin/python"
input_genbank   = "template.gb"
output_basename = "extended"
"#,
        )
        .unwrap();
        let input = PydnaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/cloning/bin/python");
        assert_eq!(input.input_genbank, Some(PathBuf::from("template.gb")));
        assert_eq!(input.output_basename, "extended");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_non_py_extension() {
        // Wrong extension is the most common typo (`.fasta`, `.gb`
        // from a copy-paste off the input field); catch it at parse
        // time so the user gets a clear error before Python is
        // invoked.
        let d = tempdir("pydna-badext");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "pydna.design"

[bio.pydna]
script          = "design.gb"
output_basename = "design"
"#,
        )
        .unwrap();
        let err = PydnaInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains(".py"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
