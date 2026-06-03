//! `[bio.openbabel]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "openbabel.convert"
//!
//! [bio.openbabel]
//! input         = "ligand.mol2"
//! output        = "ligand.sdf"
//! input_format  = "mol2"        # optional, override sniff-by-extension
//! output_format = "sdf"         # optional, override sniff-by-extension
//! gen_3d        = false         # optional, defaults to false (--gen3D)
//! add_hydrogens = false         # optional, defaults to false (-h)
//! extra_args    = ["--canonical"]
//! ```
//!
//! Open Babel sniffs format from extension by default; the explicit
//! `input_format` / `output_format` keys are escape hatches for files
//! whose extensions don't match Open Babel's table or for stripped
//! filenames. `gen_3d = true` lifts a 2D structure into 3D coordinates;
//! `add_hydrogens = true` adds explicit hydrogens before write.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct OpenBabelInput {
    pub input: PathBuf,
    pub output: PathBuf,
    pub input_format: Option<String>,
    pub output_format: Option<String>,
    pub gen_3d: bool,
    pub add_hydrogens: bool,
    pub extra_args: Vec<String>,
}

impl OpenBabelInput {
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
            .and_then(|v| v.get("openbabel"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.openbabel] section",
                    case_toml.display()
                ))
            })?;

        let input_str = block.get("input").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!("[bio.openbabel].input required"))
        })?;
        if input_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.openbabel].input must not be empty"
            )));
        }

        let output_str = block
            .get("output")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.openbabel].output required"))
            })?;
        if output_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.openbabel].output must not be empty"
            )));
        }

        let input_format = block
            .get("input_format")
            .and_then(|v| v.as_str())
            .map(|s| s.to_string());
        let output_format = block
            .get("output_format")
            .and_then(|v| v.as_str())
            .map(|s| s.to_string());

        let gen_3d = block
            .get("gen_3d")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);
        let add_hydrogens = block
            .get("add_hydrogens")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.openbabel].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.openbabel].extra_args entries must be strings"
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
            input_format,
            output_format,
            gen_3d,
            add_hydrogens,
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
        let d = tempdir("openbabel");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "openbabel.convert"

[bio.openbabel]
input  = "ligand.mol2"
output = "ligand.sdf"
"#,
        )
        .unwrap();
        let input = OpenBabelInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("ligand.mol2"));
        assert_eq!(input.output, PathBuf::from("ligand.sdf"));
        assert!(input.input_format.is_none());
        assert!(input.output_format.is_none());
        assert!(!input.gen_3d);
        assert!(!input.add_hydrogens);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_explicit_formats() {
        // Files with stripped extensions or non-standard naming need
        // the format hints (Open Babel can't sniff them).
        let d = tempdir("openbabel");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "openbabel.convert"

[bio.openbabel]
input         = "ligand.dat"
output        = "ligand.out"
input_format  = "mol2"
output_format = "sdf"
extra_args    = ["--canonical"]
"#,
        )
        .unwrap();
        let input = OpenBabelInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input_format, Some("mol2".to_string()));
        assert_eq!(input.output_format, Some("sdf".to_string()));
        assert_eq!(input.extra_args, vec!["--canonical".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_3d_and_hydrogens() {
        // The classic "lift a SMILES into a docking-ready 3D ligand"
        // pipeline — generate 3D coords + add explicit hydrogens.
        let d = tempdir("openbabel");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "openbabel.convert"

[bio.openbabel]
input         = "smiles.smi"
output        = "ligand.pdb"
gen_3d        = true
add_hydrogens = true
"#,
        )
        .unwrap();
        let input = OpenBabelInput::from_case_dir(&d).unwrap();
        assert!(input.gen_3d);
        assert!(input.add_hydrogens);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_input() {
        // Empty `input` would let the bare `obabel` invocation drop
        // through and produce a confusing tool-side error.
        let d = tempdir("openbabel");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "openbabel.convert"

[bio.openbabel]
input  = ""
output = "ligand.sdf"
"#,
        )
        .unwrap();
        let err = OpenBabelInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.openbabel].input"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
