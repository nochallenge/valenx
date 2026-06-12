//! `[bio.chopchop]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "chopchop.design"
//!
//! [bio.chopchop]
//! script          = "chopchop.py"
//! python          = "python3"             # optional, defaults to python3
//! target          = "target.fa"
//! genome          = "hg38"
//! cas_variant     = "Cas9"                # Cas9 | Cas12a | Cas13 | TALEN
//! pam             = "NGG"
//! output_basename = "guides"
//! ```
//!
//! CHOPCHOP is a Python-driven CRISPR guide-RNA design tool. The user
//! authors `chopchop.py` (or whatever filename) referenced from
//! `script`; the adapter stages the script + target FASTA into the
//! workdir and emits a flat `valenx_params.json` so the script can
//! read the parsed knobs without re-parsing `case.toml`.
//!
//! `cas_variant` is restricted to the four nuclease families CHOPCHOP
//! actively supports — `Cas9`, `Cas12a`, `Cas13`, and `TALEN` (the
//! TALEN-design path predates CRISPR but lives in the same code).
//! `pam` is the protospacer-adjacent motif the design step keys off
//! (`NGG` for SpCas9, `TTTV` for Cas12a, etc.) — left as a free-form
//! string so non-canonical variants stay expressible.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct ChopchopInput {
    pub script: PathBuf,
    pub python: String,
    pub target: PathBuf,
    pub genome: String,
    pub cas_variant: String,
    pub pam: String,
    pub output_basename: String,
}

/// Recognised `cas_variant` values. The four nuclease families
/// CHOPCHOP's design pipeline supports today.
pub const CAS_VARIANTS: &[&str] = &["Cas9", "Cas12a", "Cas13", "TALEN"];

impl ChopchopInput {
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
            .and_then(|v| v.get("chopchop"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.chopchop] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.chopchop].script required"))
            })?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.chopchop].script must not be empty"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let target = block
            .get("target")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.chopchop].target required"))
            })?;
        if target.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.chopchop].target must not be empty"
            )));
        }

        let genome = block
            .get("genome")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.chopchop].genome required"))
            })?;
        if genome.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.chopchop].genome must not be empty"
            )));
        }

        let cas_variant = block
            .get("cas_variant")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.chopchop].cas_variant required"))
            })?;
        if cas_variant.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.chopchop].cas_variant must not be empty"
            )));
        }
        if !CAS_VARIANTS.contains(&cas_variant) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.chopchop].cas_variant `{cas_variant}` not recognised; \
                 must be one of {CAS_VARIANTS:?}"
            )));
        }

        let pam = block
            .get("pam")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.chopchop].pam required")))?;
        if pam.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.chopchop].pam must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.chopchop].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.chopchop].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            target: PathBuf::from(target),
            genome: genome.to_string(),
            cas_variant: cas_variant.to_string(),
            pam: pam.to_string(),
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
        let d = tempdir("chopchop-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "chopchop.design"

[bio.chopchop]
script          = "chopchop.py"
target          = "target.fa"
genome          = "hg38"
cas_variant     = "Cas9"
pam             = "NGG"
output_basename = "guides"
"#,
        )
        .unwrap();
        let input = ChopchopInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("chopchop.py"));
        assert_eq!(input.target, PathBuf::from("target.fa"));
        assert_eq!(input.genome, "hg38");
        assert_eq!(input.cas_variant, "Cas9");
        assert_eq!(input.pam, "NGG");
        assert_eq!(input.output_basename, "guides");
        // Default Python interpreter when not explicitly pinned.
        assert_eq!(input.python, "python3");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_cas12a() {
        // Cas12a (formerly Cpf1) is the second-most-common nuclease
        // CHOPCHOP supports. Its canonical PAM is `TTTV` (T-rich,
        // immediately 5' of the protospacer) — distinct from SpCas9's
        // `NGG`. Pinning the alternate PAM here exercises the
        // free-form `pam` field alongside the variant whitelist.
        let d = tempdir("chopchop-cas12a");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "chopchop.design"

[bio.chopchop]
script          = "chopchop.py"
python          = "/opt/conda/envs/chopchop/bin/python"
target          = "exon3.fa"
genome          = "mm10"
cas_variant     = "Cas12a"
pam             = "TTTV"
output_basename = "cas12a_guides"
"#,
        )
        .unwrap();
        let input = ChopchopInput::from_case_dir(&d).unwrap();
        assert_eq!(input.cas_variant, "Cas12a");
        assert_eq!(input.pam, "TTTV");
        assert_eq!(input.genome, "mm10");
        assert_eq!(input.python, "/opt/conda/envs/chopchop/bin/python");
        assert_eq!(input.output_basename, "cas12a_guides");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_cas_variant() {
        // CHOPCHOP doesn't ship a `Cas14` design path. Reject up
        // front so the user sees the failure at validation time
        // rather than after a script-side traceback.
        let d = tempdir("chopchop-badvar");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "chopchop.design"

[bio.chopchop]
script          = "chopchop.py"
target          = "target.fa"
genome          = "hg38"
cas_variant     = "Cas14"
pam             = "NGG"
output_basename = "guides"
"#,
        )
        .unwrap();
        let err = ChopchopInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("cas_variant"), "msg: {msg}");
        assert!(msg.contains("Cas14"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_pam() {
        // PAM is the protospacer-adjacent motif the design pass keys
        // off; an empty string means CHOPCHOP would scan every
        // position in the target. Reject explicitly so the user
        // catches the typo before a multi-minute script start.
        let d = tempdir("chopchop-nopam");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "chopchop.design"

[bio.chopchop]
script          = "chopchop.py"
target          = "target.fa"
genome          = "hg38"
cas_variant     = "Cas9"
pam             = ""
output_basename = "guides"
"#,
        )
        .unwrap();
        let err = ChopchopInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("pam"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
