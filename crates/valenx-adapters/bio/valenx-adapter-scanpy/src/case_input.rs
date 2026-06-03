//! `[bio.scanpy]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "scanpy.analyse"
//!
//! [bio.scanpy]
//! script        = "analyse.py"
//! python        = "python3"           # optional, default python3
//! input_h5ad    = "raw_counts.h5ad"
//! output_h5ad   = "annotated.h5ad"
//! n_top_genes   = 2000                # optional, default 2000
//! n_pcs         = 50                  # optional, default 50
//! n_neighbors   = 15                  # optional, default 15
//! resolution    = 1.0                 # optional, default 1.0
//! ```

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct ScanpyInput {
    pub script: PathBuf,
    pub python: String,
    pub input_h5ad: PathBuf,
    pub output_h5ad: String,
    /// Number of highly-variable genes to retain via
    /// `sc.pp.highly_variable_genes`. Scanpy's standard recipe
    /// defaults to 2000.
    pub n_top_genes: u32,
    /// Number of principal components to compute via `sc.tl.pca`.
    /// 50 is the canonical default for single-cell workflows.
    pub n_pcs: u32,
    /// k for `sc.pp.neighbors` — the kNN graph used by UMAP / Leiden.
    pub n_neighbors: u32,
    /// Leiden / Louvain clustering resolution. Higher resolution
    /// yields more clusters; 1.0 is the standard starting point.
    pub resolution: f64,
}

const DEFAULT_N_TOP_GENES: u32 = 2000;
const DEFAULT_N_PCS: u32 = 50;
const DEFAULT_N_NEIGHBORS: u32 = 15;
const DEFAULT_RESOLUTION: f64 = 1.0;

impl ScanpyInput {
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
            .and_then(|v| v.get("scanpy"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.scanpy] section",
                    case_toml.display()
                ))
            })?;

        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .map(|s| s.trim().to_string())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.scanpy].script required")))?;
        if script.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scanpy].script must be non-empty"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let input_h5ad = block
            .get("input_h5ad")
            .and_then(|v| v.as_str())
            .map(|s| s.trim().to_string())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.scanpy].input_h5ad required"))
            })?;
        if input_h5ad.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scanpy].input_h5ad must be non-empty"
            )));
        }

        let output_h5ad = block
            .get("output_h5ad")
            .and_then(|v| v.as_str())
            .map(|s| s.trim().to_string())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.scanpy].output_h5ad required"))
            })?;
        if output_h5ad.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scanpy].output_h5ad must be non-empty"
            )));
        }

        let n_top_genes = block
            .get("n_top_genes")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_N_TOP_GENES);
        if n_top_genes < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scanpy].n_top_genes must be >= 1, got {n_top_genes}"
            )));
        }

        let n_pcs = block
            .get("n_pcs")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_N_PCS);
        if n_pcs < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scanpy].n_pcs must be >= 1, got {n_pcs}"
            )));
        }

        let n_neighbors = block
            .get("n_neighbors")
            .and_then(|v| v.as_integer())
            .map(|n| n as u32)
            .unwrap_or(DEFAULT_N_NEIGHBORS);
        if n_neighbors < 1 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scanpy].n_neighbors must be >= 1, got {n_neighbors}"
            )));
        }

        let resolution = block
            .get("resolution")
            .and_then(|v| v.as_float())
            .or_else(|| {
                // TOML allows integer literals for floats — accept
                // `resolution = 1` as `1.0`.
                block
                    .get("resolution")
                    .and_then(|v| v.as_integer())
                    .map(|i| i as f64)
            })
            .unwrap_or(DEFAULT_RESOLUTION);
        if !(resolution > 0.0 && resolution.is_finite()) {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.scanpy].resolution must be > 0 and finite, got {resolution}"
            )));
        }

        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_h5ad: PathBuf::from(input_h5ad),
            output_h5ad,
            n_top_genes,
            n_pcs,
            n_neighbors,
            resolution,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("scanpy-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "scanpy.analyse"

[bio.scanpy]
script      = "analyse.py"
input_h5ad  = "raw.h5ad"
output_h5ad = "annotated.h5ad"
"#,
        )
        .unwrap();
        let input = ScanpyInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("analyse.py"));
        assert_eq!(input.python, "python3");
        assert_eq!(input.input_h5ad, PathBuf::from("raw.h5ad"));
        assert_eq!(input.output_h5ad, "annotated.h5ad");
        // Defaults pinned to Scanpy's canonical recipe.
        assert_eq!(input.n_top_genes, 2000);
        assert_eq!(input.n_pcs, 50);
        assert_eq!(input.n_neighbors, 15);
        assert!((input.resolution - 1.0).abs() < f64::EPSILON);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_overrides() {
        let d = tempdir("scanpy-overrides");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "scanpy.analyse"

[bio.scanpy]
script       = "analyse.py"
python       = "/opt/conda/envs/sc/bin/python"
input_h5ad   = "raw.h5ad"
output_h5ad  = "annotated.h5ad"
n_top_genes  = 4000
n_pcs        = 30
n_neighbors  = 20
resolution   = 0.5
"#,
        )
        .unwrap();
        let input = ScanpyInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/sc/bin/python");
        assert_eq!(input.n_top_genes, 4000);
        assert_eq!(input.n_pcs, 30);
        assert_eq!(input.n_neighbors, 20);
        assert!((input.resolution - 0.5).abs() < f64::EPSILON);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_pcs() {
        let d = tempdir("scanpy-zeropcs");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "scanpy.analyse"

[bio.scanpy]
script      = "analyse.py"
input_h5ad  = "raw.h5ad"
output_h5ad = "annotated.h5ad"
n_pcs       = 0
"#,
        )
        .unwrap();
        let err = ScanpyInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("n_pcs"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_resolution() {
        let d = tempdir("scanpy-zerores");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "scanpy.analyse"

[bio.scanpy]
script      = "analyse.py"
input_h5ad  = "raw.h5ad"
output_h5ad = "annotated.h5ad"
resolution  = 0.0
"#,
        )
        .unwrap();
        let err = ScanpyInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("resolution"));
        let _ = std::fs::remove_dir_all(&d);
    }
}
