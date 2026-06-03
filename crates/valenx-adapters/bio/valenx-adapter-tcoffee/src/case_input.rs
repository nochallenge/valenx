//! `[bio.tcoffee]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "tcoffee.align"
//!
//! [bio.tcoffee]
//! input           = "input.fa"
//! output_basename = "alignment"
//! outfmt          = "clustalw"           # optional, defaults to "clustalw"
//! mode            = "regular"            # optional, omitted = t_coffee's default
//! extra_args      = ["-n_core=4"]        # optional, defaults to []
//! ```
//!
//! `outfmt` selects T-Coffee's `-output=` flag. T-Coffee accepts a long
//! list of format names — `clustalw`, `fasta_aln`, `phylip`, `msf`,
//! `score_html`, `score_ascii`, etc. We don't gate on the list; T-Coffee
//! itself errors clearly on unknown values, and gating would lock the
//! adapter against future format additions.
//!
//! `mode` is optional. When omitted, T-Coffee picks its built-in default
//! (`regular`, the classic progressive alignment). Setting it routes
//! through one of T-Coffee's higher-level pipelines like `expresso`
//! (structure-aware) or `psicoffee` (PSI-BLAST profile-driven).

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct TCoffeeInput {
    pub input: PathBuf,
    pub output_basename: String,
    pub outfmt: String,
    pub mode: Option<String>,
    pub extra_args: Vec<String>,
}

impl TCoffeeInput {
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
            .and_then(|v| v.get("tcoffee"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.tcoffee] section",
                    case_toml.display()
                ))
            })?;

        let input_str = block.get("input").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.tcoffee].input required (path to multi-FASTA input)"
            ))
        })?;
        if input_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.tcoffee].input must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.tcoffee].output_basename required \
                     (basename for the alignment output, e.g. \"alignment\")"
                ))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.tcoffee].output_basename must not be empty"
            )));
        }

        let outfmt = match block.get("outfmt") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.tcoffee].outfmt must be a string"))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.tcoffee].outfmt must not be empty"
                    )));
                }
                s.to_string()
            }
            None => "clustalw".to_string(),
        };

        let mode = match block.get("mode") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.tcoffee].mode must be a string"))
                })?;
                if s.is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.tcoffee].mode must not be empty when set; \
                         omit the key entirely to use t_coffee's default"
                    )));
                }
                Some(s.to_string())
            }
            None => None,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.tcoffee].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.tcoffee].extra_args entries must be strings"
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
            output_basename: output_basename.to_string(),
            outfmt,
            mode,
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_case() {
        // Only the two required keys: input, output_basename. Defaults:
        // outfmt = "clustalw", mode = None, no extras.
        let d = tempdir("tcoffee");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "tcoffee.align"

[bio.tcoffee]
input           = "input.fa"
output_basename = "alignment"
"#,
        )
        .unwrap();
        let input = TCoffeeInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("input.fa"));
        assert_eq!(input.output_basename, "alignment");
        assert_eq!(input.outfmt, "clustalw");
        assert_eq!(input.mode, None);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_full_case_with_overrides() {
        // FASTA-style output, expresso mode, an n_core extra. Confirms
        // every optional key threads through cleanly — including `mode`,
        // which becomes `Some(_)` rather than the default `None`.
        let d = tempdir("tcoffee");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "tcoffee.align"

[bio.tcoffee]
input           = "input.fa"
output_basename = "aligned"
outfmt          = "fasta_aln"
mode            = "expresso"
extra_args      = ["-n_core=4", "-quiet"]
"#,
        )
        .unwrap();
        let input = TCoffeeInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("input.fa"));
        assert_eq!(input.output_basename, "aligned");
        assert_eq!(input.outfmt, "fasta_aln");
        assert_eq!(input.mode, Some("expresso".to_string()));
        assert_eq!(
            input.extra_args,
            vec!["-n_core=4".to_string(), "-quiet".to_string()]
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("tcoffee");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = TCoffeeInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.tcoffee]"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
