//! `[bio.seurat]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "seurat.analyze"
//!
//! [bio.seurat]
//! script           = "analysis.R"
//! rscript          = "Rscript"           # optional, defaults to "Rscript"
//! # input_data     = "matrix.h5"         # optional, omit to skip
//! output_basename  = "analysis"
//! ```
//!
//! Seurat is the dominant R-based toolkit for single-cell RNA-seq:
//! QC, normalisation, dimensionality reduction (PCA / UMAP),
//! clustering, differential expression, integration. The adapter
//! itself doesn't generate R; the user supplies an `analysis.R`
//! that does `library(Seurat)` and the actual data work. We just
//! spawn `Rscript <analysis.R>` after staging the script (and any
//! optional input matrix) into the workdir.
//!
//! `input_data` is optional: omit it for scripts that fetch their
//! own data (e.g. `SeuratData::InstallData("pbmc3k")`), or supply
//! a path to a matrix file (`.h5` / `.mtx` / `.rds`) that the
//! script reads via the staged filename.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct SeuratInput {
    /// Path to the user-authored R driver script (relative to the
    /// case directory, or absolute). Must end in `.R`
    /// (case-insensitive).
    pub script: PathBuf,
    /// Rscript binary name / path. Defaults to `Rscript` so the
    /// adapter walks PATH; can be pinned to an absolute path for
    /// users with multiple R installs.
    pub rscript: String,
    /// Optional path to an input matrix the script reads as its
    /// starting point (`.h5` / `.mtx` / `.rds` are typical). `None`
    /// means the script fetches / synthesises its own data.
    pub input_data: Option<PathBuf>,
    /// Filename stem for outputs. The script writes `<basename>*.rds`
    /// (Seurat objects), `<basename>*.csv` (tables), and
    /// `<basename>*.png` (plots) into the workdir.
    pub output_basename: String,
}

impl SeuratInput {
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
            .and_then(|v| v.get("seurat"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.seurat] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.seurat].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.seurat].script must not be empty"
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
                "[bio.seurat].script `{script}` must end in `.R`"
            )));
        }

        let rscript = block
            .get("rscript")
            .and_then(|v| v.as_str())
            .unwrap_or("Rscript")
            .to_string();

        let input_data = match block.get("input_data").and_then(|v| v.as_str()) {
            Some(s) if !s.is_empty() => Some(PathBuf::from(s)),
            _ => None,
        };

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.seurat].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.seurat].output_basename must not be empty"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            rscript,
            input_data,
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
        let d = tempdir("seurat-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "seurat.analyze"

[bio.seurat]
script          = "analysis.R"
output_basename = "analysis"
"#,
        )
        .unwrap();
        let input = SeuratInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("analysis.R"));
        assert_eq!(input.rscript, "Rscript");
        // No input_data — script fetches / synthesises its own.
        assert_eq!(input.input_data, None);
        assert_eq!(input.output_basename, "analysis");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_input_data() {
        // Pinned Rscript path + an input matrix the script reads.
        let d = tempdir("seurat-input");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "seurat.analyze"

[bio.seurat]
script          = "qc.R"
rscript         = "/opt/R/4.3.0/bin/Rscript"
input_data      = "matrix.h5"
output_basename = "qc"
"#,
        )
        .unwrap();
        let input = SeuratInput::from_case_dir(&d).unwrap();
        assert_eq!(input.rscript, "/opt/R/4.3.0/bin/Rscript");
        assert_eq!(input.input_data, Some(PathBuf::from("matrix.h5")));
        assert_eq!(input.output_basename, "qc");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_non_r_extension() {
        // Wrong extension is the most common typo (`.r` works,
        // `.py` / `.txt` do not). Catch it at parse time so the
        // user gets a clear error before Rscript is invoked.
        let d = tempdir("seurat-badext");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "seurat.analyze"

[bio.seurat]
script          = "analysis.py"
output_basename = "analysis"
"#,
        )
        .unwrap();
        let err = SeuratInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains(".R"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
