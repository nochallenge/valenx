//! JSON manifest for sweep datasets — the top-level descriptor written
//! alongside per-sample `.npy` files in an exported sweep dataset.
//!
//! Schema per RFC 0012 §"Output layout."

use std::path::Path;

use crate::ExportError;

/// Top-level manifest written alongside a sweep export's per-sample
/// `.npy` files. Schema per RFC 0012 §"Output layout."
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct ExportManifest {
    /// Bumped per breaking schema change. Loaders should refuse
    /// versions they don't know.
    pub valenx_export_version: String,
    /// Source the dataset was built from (sweep workdir path or
    /// project name; free-form for now).
    pub sweep_source: String,
    /// Number of samples in the dataset.
    pub sample_count: usize,
    /// Per-input declaration. Order matters — loaders concatenate
    /// inputs in declaration order.
    pub inputs: ExportSchemaSection,
    /// Per-output declaration.
    pub outputs: ExportSchemaSection,
    /// Train/val/test split, if applied. None = single-split dataset.
    #[serde(default)]
    pub split: Option<DatasetSplit>,
    /// Provenance bag — git commit, Valenx version, adapter versions,
    /// per-case source-toml hashes. Free-form for the v0 scaffold;
    /// schema firms up when the manifest writer wires into the
    /// real sweep pipeline.
    #[serde(default)]
    pub provenance: serde_json::Value,
}

/// One inputs/outputs declaration: a list of array schemas.
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct ExportSchemaSection {
    pub schema: Vec<ExportArraySchema>,
}

/// Per-array metadata. The actual numpy arrays live in
/// `inputs/sample_NNNN.npz` (or `.npy` per-array, depending on the
/// future `--format` flag).
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct ExportArraySchema {
    /// Array key inside the sample's npz/npy file.
    pub name: String,
    /// Tensor shape per sample. `[1]` for a scalar metric, `[3]`
    /// for a velocity vector, `[N]` for a per-node field.
    pub shape: Vec<usize>,
    /// Short SI symbol or `"1"` for dimensionless.
    pub units: String,
    /// numpy dtype string ("f32", "f64", "i32", etc.).
    pub dtype: String,
    /// Optional location for fields ("OnNode" / "OnCell"); omitted
    /// for scalars.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub location: Option<String>,
    /// Optional mesh id when this array is defined on a common mesh
    /// (geometry sweeps that interpolate to a baseline). Omitted
    /// for scalars.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mesh_id: Option<String>,
}

/// Train/val/test fractions. Sums must be ≤ 1.0; remainder is
/// dropped as a "holdout" pool.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct DatasetSplit {
    pub train: f32,
    pub val: f32,
    pub test: f32,
}

impl ExportManifest {
    /// Convenience constructor stamping the current export version.
    pub fn new(sweep_source: impl Into<String>, sample_count: usize) -> Self {
        Self {
            valenx_export_version: "0.1".to_string(),
            sweep_source: sweep_source.into(),
            sample_count,
            inputs: ExportSchemaSection::default(),
            outputs: ExportSchemaSection::default(),
            split: None,
            provenance: serde_json::json!({}),
        }
    }
}

/// Write an [`ExportManifest`] to `path` as pretty-printed JSON.
pub fn write_manifest(manifest: &ExportManifest, path: &Path) -> Result<(), ExportError> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| ExportError::Io {
            path: parent.to_path_buf(),
            source: e,
        })?;
    }
    let text = serde_json::to_string_pretty(manifest).map_err(|e| ExportError::Io {
        path: path.to_path_buf(),
        source: std::io::Error::other(e.to_string()),
    })?;
    valenx_core::io_caps::atomic_write_str(path, &text).map_err(|e| ExportError::Io {
        path: path.to_path_buf(),
        source: e,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn manifest_round_trips_through_json() {
        let mut manifest = ExportManifest::new("test-sweep", 200);
        manifest.inputs.schema.push(ExportArraySchema {
            name: "aoa_deg".into(),
            shape: vec![1],
            units: "deg".into(),
            dtype: "f32".into(),
            location: None,
            mesh_id: None,
        });
        manifest.outputs.schema.push(ExportArraySchema {
            name: "pressure".into(),
            shape: vec![12345],
            units: "Pa".into(),
            dtype: "f32".into(),
            location: Some("OnNode".into()),
            mesh_id: Some("common-airfoil".into()),
        });
        manifest.split = Some(DatasetSplit {
            train: 0.7,
            val: 0.15,
            test: 0.15,
        });
        manifest.provenance = serde_json::json!({
            "valenx_version": "0.1.0-alpha.1",
            "git_commit": "deadbeef",
        });

        let json = serde_json::to_string(&manifest).expect("serialize");
        let parsed: ExportManifest = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(parsed.sample_count, 200);
        assert_eq!(parsed.inputs.schema.len(), 1);
        assert_eq!(parsed.outputs.schema.len(), 1);
        assert_eq!(parsed.outputs.schema[0].location.as_deref(), Some("OnNode"));
        assert_eq!(parsed.split.as_ref().map(|s| s.train), Some(0.7));
    }

    #[test]
    fn manifest_omits_array_optional_fields_when_unset() {
        // ExportArraySchema's location / mesh_id are Option +
        // #[serde(skip_serializing_if = "Option::is_none")] so they
        // don't appear in the JSON when unset. Loaders shouldn't
        // have to special-case nulls.
        let mut manifest = ExportManifest::new("simple", 10);
        manifest.inputs.schema.push(ExportArraySchema {
            name: "x".into(),
            shape: vec![1],
            units: "1".into(),
            dtype: "f32".into(),
            location: None,
            mesh_id: None,
        });
        let json = serde_json::to_string(&manifest).expect("serialize");
        assert!(!json.contains("\"location\""), "got: {json}");
        assert!(!json.contains("\"mesh_id\""), "got: {json}");
    }

    #[test]
    fn write_manifest_creates_a_file_loaders_can_consume() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-export-manifest-{}.json",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let manifest = ExportManifest::new("smoke", 3);
        write_manifest(&manifest, &tmp).expect("write");
        let text = std::fs::read_to_string(&tmp).expect("read");
        // Round-trip the file back through serde to confirm shape.
        let parsed: ExportManifest = serde_json::from_str(&text).expect("parse");
        assert_eq!(parsed.sample_count, 3);
        assert_eq!(parsed.valenx_export_version, "0.1");
        let _ = std::fs::remove_file(&tmp);
    }
}
