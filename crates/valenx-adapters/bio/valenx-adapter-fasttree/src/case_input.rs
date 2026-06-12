//! `[bio.fasttree]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "fasttree.tree"
//!
//! [bio.fasttree]
//! alignment  = "aln.fa"
//! output     = "tree.nwk"
//! seq_type   = "nt"               # "nt" | "aa"
//! use_gtr    = true               # nucleotide-only; ignored for aa
//! gamma      = false              # optional, defaults to false
//! extra_args = ["-spr", "4"]      # optional, defaults to []
//! ```
//!
//! `seq_type` selects the alphabet:
//!
//! - `nt` — nucleotide alignment; FastTree's default model is JC
//!   (Jukes-Cantor). Setting `use_gtr = true` switches to the
//!   General Time Reversible model (more parameters, more accurate).
//! - `aa` — amino-acid alignment; FastTree's default model is JTT;
//!   `use_gtr` is meaningless here (FastTree has no AA-GTR mode)
//!   and is ignored when emitting the command line.
//!
//! `gamma` enables a discrete gamma rate-heterogeneity model
//! (`-gamma`) — slower but more biologically defensible.
//!
//! FastTree writes the Newick tree to **stdout** with no `-o` flag,
//! so the adapter's `run()` must redirect stdout to a file. The
//! `output` field names that file (resolved against the workdir).

use std::path::PathBuf;
use valenx_core::AdapterError;

/// Canonical FastTree alphabet list. Module-public so the adapter
/// can surface the supported values to the UI.
pub const SUPPORTED_SEQ_TYPES: &[&str] = &["nt", "aa"];

#[derive(Clone, Debug, PartialEq)]
pub struct FastTreeInput {
    pub alignment: PathBuf,
    pub output: PathBuf,
    pub seq_type: String,
    pub use_gtr: bool,
    pub gamma: bool,
    pub extra_args: Vec<String>,
}

impl FastTreeInput {
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
            .and_then(|v| v.get("fasttree"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.fasttree] section",
                    case_toml.display()
                ))
            })?;

        let alignment_str = block
            .get("alignment")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.fasttree].alignment required (path to multi-FASTA alignment)"
                ))
            })?;
        if alignment_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.fasttree].alignment must not be empty"
            )));
        }

        let output_str = block
            .get("output")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.fasttree].output required (path for FastTree's Newick output, \
                     resolved against the workdir)"
                ))
            })?;
        if output_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.fasttree].output must not be empty"
            )));
        }

        let seq_type = match block.get("seq_type") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.fasttree].seq_type must be a string"))
                })?;
                if !SUPPORTED_SEQ_TYPES.contains(&s) {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.fasttree].seq_type `{s}` not recognised — \
                         expected one of {SUPPORTED_SEQ_TYPES:?}"
                    )));
                }
                s.to_string()
            }
            None => "nt".to_string(),
        };

        let use_gtr = block
            .get("use_gtr")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let gamma = block
            .get("gamma")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.fasttree].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.fasttree].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            alignment: PathBuf::from(alignment_str),
            output: PathBuf::from(output_str),
            seq_type,
            use_gtr,
            gamma,
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_nt_minimal() {
        // Minimal nucleotide config: alignment + output + nt seq_type.
        // Other fields take defaults: no GTR, no gamma, no extras.
        let d = tempdir("fasttree");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "fasttree.tree"

[bio.fasttree]
alignment = "aln.fa"
output    = "tree.nwk"
seq_type  = "nt"
"#,
        )
        .unwrap();
        let input = FastTreeInput::from_case_dir(&d).unwrap();
        assert_eq!(input.alignment, PathBuf::from("aln.fa"));
        assert_eq!(input.output, PathBuf::from("tree.nwk"));
        assert_eq!(input.seq_type, "nt");
        assert!(!input.use_gtr);
        assert!(!input.gamma);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_aa_with_gamma() {
        // Amino-acid alignment with gamma rate-heterogeneity. We
        // accept use_gtr = true here even though FastTree ignores
        // it for AA — validation isn't this adapter's job, the CLI
        // composer in lib.rs is responsible for skipping the flag.
        let d = tempdir("fasttree");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "fasttree.tree"

[bio.fasttree]
alignment = "aln.fa"
output    = "tree.nwk"
seq_type  = "aa"
gamma     = true
"#,
        )
        .unwrap();
        let input = FastTreeInput::from_case_dir(&d).unwrap();
        assert_eq!(input.seq_type, "aa");
        assert!(input.gamma);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_gtr_off() {
        // Nucleotide alignment with use_gtr explicitly set to false.
        // This is the JC default, but exercising the path proves the
        // bool parser doesn't accidentally promote `false` to true.
        let d = tempdir("fasttree");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "fasttree.tree"

[bio.fasttree]
alignment  = "aln.fa"
output     = "tree.nwk"
seq_type   = "nt"
use_gtr    = false
extra_args = ["-spr", "4"]
"#,
        )
        .unwrap();
        let input = FastTreeInput::from_case_dir(&d).unwrap();
        assert_eq!(input.seq_type, "nt");
        assert!(!input.use_gtr);
        assert_eq!(input.extra_args, vec!["-spr".to_string(), "4".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_unknown_seq_type() {
        // `dna` is the obvious wrong-name; reject so the user fixes
        // the typo instead of silently getting wrong-alphabet output.
        let d = tempdir("fasttree");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "fasttree.tree"

[bio.fasttree]
alignment = "aln.fa"
output    = "tree.nwk"
seq_type  = "dna"
"#,
        )
        .unwrap();
        let err = FastTreeInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not recognised"), "msg: {msg}");
        assert!(msg.contains("nt"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
