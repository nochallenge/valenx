//! `[bio.crispor]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "crispor.design"
//!
//! [bio.crispor]
//! script          = "crispor.py"
//! python          = "python3"             # optional, defaults to python3
//! target          = "target.fa"
//! genome          = "hg38"
//! pam             = "NGG"
//! batch_id        = "batch-2026-04-30"    # optional, default null
//! output_basename = "guides"
//! ```
//!
//! CRISPOR is the second-generation CRISPR guide-RNA design tool from
//! the Tefor Infrastructure / UCSC Genome Browser team. Compared with
//! CHOPCHOP, CRISPOR's distinguishing feature is its rigorous
//! off-target prediction pass — it scores potential off-targets against
//! a reference genome with the CFD scoring model and surfaces a
//! per-guide MIT specificity score. The Python entry point is the same
//! script (`crispor.py`) used by the public web service.
//!
//! `batch_id` is an optional CRISPOR-side identifier the user can
//! attach to a run for cross-referencing CRISPOR's own batch history;
//! it's emitted into `valenx_params.json` as either a JSON string or
//! `null` so the user script doesn't need to special-case its absence.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct CrisporInput {
    pub script: PathBuf,
    pub python: String,
    pub target: PathBuf,
    pub genome: String,
    pub pam: String,
    pub batch_id: Option<String>,
    pub output_basename: String,
}

impl CrisporInput {
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
            .and_then(|v| v.get("crispor"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.crispor] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.crispor].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.crispor].script must not be empty"
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
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.crispor].target required")))?;
        if target.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.crispor].target must not be empty"
            )));
        }

        let genome = block
            .get("genome")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.crispor].genome required")))?;
        if genome.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.crispor].genome must not be empty"
            )));
        }

        let pam = block
            .get("pam")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.crispor].pam required")))?;
        if pam.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.crispor].pam must not be empty"
            )));
        }

        // batch_id is optional. When the field is present we still
        // require it to be a string (TOML doesn't have nullable
        // strings, so an empty string is the closest analogue —
        // accept it but treat as absent so downstream code can rely
        // on `Option::is_some` meaning "user provided a real id").
        let batch_id = match block.get("batch_id") {
            None => None,
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.crispor].batch_id must be a string"))
                })?;
                if s.is_empty() {
                    None
                } else {
                    Some(s.to_string())
                }
            }
        };

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.crispor].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.crispor].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            target: PathBuf::from(target),
            genome: genome.to_string(),
            pam: pam.to_string(),
            batch_id,
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
        let d = tempdir("crispor-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "crispor.design"

[bio.crispor]
script          = "crispor.py"
target          = "target.fa"
genome          = "hg38"
pam             = "NGG"
output_basename = "guides"
"#,
        )
        .unwrap();
        let input = CrisporInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("crispor.py"));
        assert_eq!(input.target, PathBuf::from("target.fa"));
        assert_eq!(input.genome, "hg38");
        assert_eq!(input.pam, "NGG");
        assert_eq!(input.output_basename, "guides");
        // batch_id is optional and defaults to None when omitted.
        assert!(input.batch_id.is_none());
        // Default Python interpreter when not explicitly pinned.
        assert_eq!(input.python, "python3");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_batch_id() {
        // batch_id is a CRISPOR-side identifier the user can attach
        // for cross-referencing CRISPOR's batch history. When present
        // it surfaces verbatim through the parsed input.
        let d = tempdir("crispor-batch");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "crispor.design"

[bio.crispor]
script          = "crispor.py"
python          = "/opt/conda/envs/crispor/bin/python"
target          = "exon3.fa"
genome          = "mm10"
pam             = "NGG"
batch_id        = "batch-2026-04-30"
output_basename = "mm10_guides"
"#,
        )
        .unwrap();
        let input = CrisporInput::from_case_dir(&d).unwrap();
        assert_eq!(input.batch_id.as_deref(), Some("batch-2026-04-30"));
        assert_eq!(input.genome, "mm10");
        assert_eq!(input.python, "/opt/conda/envs/crispor/bin/python");
        assert_eq!(input.output_basename, "mm10_guides");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_pam() {
        // PAM is the protospacer-adjacent motif the design pass keys
        // off; an empty string is meaningless. Reject explicitly so
        // the user catches the typo before a multi-minute script start.
        let d = tempdir("crispor-nopam");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "crispor.design"

[bio.crispor]
script          = "crispor.py"
target          = "target.fa"
genome          = "hg38"
pam             = ""
output_basename = "guides"
"#,
        )
        .unwrap();
        let err = CrisporInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("pam"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_genome() {
        // CRISPOR's off-target pass needs a concrete reference
        // genome assembly — an empty string would silently disable
        // the off-target check. Reject up front.
        let d = tempdir("crispor-nogenome");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "crispor.design"

[bio.crispor]
script          = "crispor.py"
target          = "target.fa"
genome          = ""
pam             = "NGG"
output_basename = "guides"
"#,
        )
        .unwrap();
        let err = CrisporInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("genome"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
