//! Sweep-dataset batch exporter — packs N samples' inputs/outputs
//! into the canonical (inputs.npy + outputs.npy + sample_ids.json +
//! manifest.json) layout per RFC 0012.

use std::path::Path;

use thiserror::Error;

use valenx_fields::Results;

use crate::manifest::{
    write_manifest, DatasetSplit, ExportArraySchema, ExportManifest, ExportSchemaSection,
};
use crate::npy::write_npy_f64_nd;
use crate::ExportError;

/// One sample in a sweep export — a derived case's input parameter
/// values + its collected scalar outputs.
///
/// Inputs are pre-extracted to (name, f64) pairs because the caller
/// owns the policy decision of which substitutions count as numeric
/// (string substitutions like `"kEpsilon"` need either one-hot
/// encoding or to be tracked separately as categorical features —
/// neither is in scope for the v0 exporter).
pub struct Sample<'a> {
    /// Stable id within the sweep (matches `DerivedCase::id`).
    pub id: String,
    /// Numeric inputs that vary across the sweep.
    pub inputs: Vec<(String, f64)>,
    /// Adapter-collected results for this run.
    pub outputs: &'a Results,
}

/// Configuration for [`export_sweep_dataset`]. Captures which output
/// scalars to extract per sample, optional split fractions, and the
/// provenance bag stamped into the manifest.
#[derive(Clone, Debug)]
pub struct DatasetExportConfig {
    /// Scalar names to extract from each sample's `Results.scalars`.
    /// Order matters — the dataset's per-sample output vector packs
    /// values in this order. Missing scalars produce a structured
    /// error; the export stops rather than silently emitting NaNs.
    pub output_names: Vec<String>,
    /// Train/val/test split applied to the sample list. Stamped into
    /// the manifest as informational metadata; the exporter does NOT
    /// shuffle or partition the data — that's the loader's job.
    pub split: Option<DatasetSplit>,
    /// Provenance bag stamped into the manifest. Loaders surface this
    /// as the dataset's "where did this come from?" record.
    pub provenance: serde_json::Value,
}

impl DatasetExportConfig {
    /// Builder convenience: a config with a list of output scalar
    /// names and no split / provenance.
    pub fn from_output_names<I, S>(names: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        Self {
            output_names: names.into_iter().map(Into::into).collect(),
            split: None,
            provenance: serde_json::json!({}),
        }
    }
}

/// Errors that bubble up from [`export_sweep_dataset`] — distinct
/// from [`ExportError`] because the dataset exporter has higher-level
/// failure modes (missing scalar, ragged input vectors) on top of the
/// raw filesystem errors.
#[derive(Debug, Error)]
pub enum DatasetExportError {
    #[error("filesystem: {0}")]
    Io(#[from] ExportError),
    #[error("sample `{sample_id}` is missing scalar `{name}` declared in output_names")]
    MissingScalar { sample_id: String, name: String },
    #[error("sample `{sample_id}` has {actual} inputs, but the first sample declared {expected}")]
    RaggedInputs {
        sample_id: String,
        actual: usize,
        expected: usize,
    },
    #[error(
        "sample `{sample_id}` input #{index} is named `{actual}`, but the first sample named it `{expected}`"
    )]
    InputNameMismatch {
        sample_id: String,
        index: usize,
        expected: String,
        actual: String,
    },
}

/// Write a sweep's samples as an ML-ready dataset.
///
/// Layout per RFC 0012 §"Output layout" (Stacked variant — N samples
/// concatenated into one inputs.npy + one outputs.npy):
///
/// ```text
/// out_dir/
/// ├── inputs.npy         # f64, shape (n_samples, n_inputs)
/// ├── outputs.npy        # f64, shape (n_samples, n_outputs)
/// ├── sample_ids.json    # ["sweep-000", "sweep-001", …]
/// └── manifest.json      # ExportManifest
/// ```
///
/// All samples must share the same input schema (same names in the
/// same order). Scalar outputs are pulled from each sample's
/// `Results.scalars` by the names listed in `config.output_names`.
///
/// Returns the manifest written to disk so the caller can stash it
/// for downstream use without re-reading the file.
pub fn export_sweep_dataset(
    samples: &[Sample<'_>],
    config: &DatasetExportConfig,
    out_dir: &Path,
    sweep_source: &str,
) -> Result<ExportManifest, DatasetExportError> {
    std::fs::create_dir_all(out_dir).map_err(|e| {
        DatasetExportError::Io(ExportError::Io {
            path: out_dir.to_path_buf(),
            source: e,
        })
    })?;

    // Empty sample list = empty manifest, no .npy files. Loader can
    // detect zero-sample datasets and skip rather than choking on a
    // shape-(0, N) array.
    if samples.is_empty() {
        let manifest = ExportManifest {
            valenx_export_version: "0.1".into(),
            sweep_source: sweep_source.into(),
            sample_count: 0,
            inputs: ExportSchemaSection::default(),
            outputs: ExportSchemaSection::default(),
            split: config.split.clone(),
            provenance: config.provenance.clone(),
        };
        write_manifest(&manifest, &out_dir.join("manifest.json"))
            .map_err(DatasetExportError::Io)?;
        return Ok(manifest);
    }

    // Lock in the input schema from the first sample.
    let first = &samples[0];
    let n_inputs = first.inputs.len();
    let input_names: Vec<String> = first.inputs.iter().map(|(n, _)| n.clone()).collect();

    // Sanity-check the rest of the samples agree on the schema. We
    // catch ragged sweeps here rather than producing silently-wrong
    // arrays.
    for s in samples.iter().skip(1) {
        if s.inputs.len() != n_inputs {
            return Err(DatasetExportError::RaggedInputs {
                sample_id: s.id.clone(),
                actual: s.inputs.len(),
                expected: n_inputs,
            });
        }
        for (idx, ((expected, _), (actual, _))) in
            first.inputs.iter().zip(s.inputs.iter()).enumerate()
        {
            if expected != actual {
                return Err(DatasetExportError::InputNameMismatch {
                    sample_id: s.id.clone(),
                    index: idx,
                    expected: expected.clone(),
                    actual: actual.clone(),
                });
            }
        }
    }

    // Pack inputs into one (n_samples, n_inputs) f64 array.
    let n_samples = samples.len();
    let mut inputs_flat: Vec<f64> = Vec::with_capacity(n_samples * n_inputs);
    for s in samples {
        for (_, v) in &s.inputs {
            inputs_flat.push(*v);
        }
    }
    write_npy_f64_nd(
        &out_dir.join("inputs.npy"),
        &inputs_flat,
        &[n_samples, n_inputs],
    )
    .map_err(DatasetExportError::Io)?;

    // Pack outputs into one (n_samples, n_outputs) f64 array. Look
    // up scalars by name; missing scalars are a hard error so the
    // user knows the dataset is incomplete before training starts.
    let n_outputs = config.output_names.len();
    let mut outputs_flat: Vec<f64> = Vec::with_capacity(n_samples * n_outputs);
    let mut output_units: Vec<String> = Vec::with_capacity(n_outputs);
    for name in &config.output_names {
        let units = first
            .outputs
            .scalars
            .get(name)
            .map(|r| r.units.display.unwrap_or("1").to_string())
            .unwrap_or_else(|| "1".to_string());
        output_units.push(units);
    }
    for s in samples {
        for name in &config.output_names {
            let record =
                s.outputs
                    .scalars
                    .get(name)
                    .ok_or_else(|| DatasetExportError::MissingScalar {
                        sample_id: s.id.clone(),
                        name: name.clone(),
                    })?;
            outputs_flat.push(record.value);
        }
    }
    write_npy_f64_nd(
        &out_dir.join("outputs.npy"),
        &outputs_flat,
        &[n_samples, n_outputs],
    )
    .map_err(DatasetExportError::Io)?;

    // Sidecar: sample ids in declaration order so the loader can map
    // a row index back to the sweep run that produced it.
    let ids: Vec<String> = samples.iter().map(|s| s.id.clone()).collect();
    let ids_json = serde_json::to_string_pretty(&ids).map_err(|e| {
        DatasetExportError::Io(ExportError::Io {
            path: out_dir.join("sample_ids.json"),
            source: std::io::Error::other(e.to_string()),
        })
    })?;
    valenx_core::io_caps::atomic_write_str(&out_dir.join("sample_ids.json"), &ids_json).map_err(
        |e| {
            DatasetExportError::Io(ExportError::Io {
                path: out_dir.join("sample_ids.json"),
                source: e,
            })
        },
    )?;

    // Build + write the top-level manifest.
    let manifest = ExportManifest {
        valenx_export_version: "0.1".into(),
        sweep_source: sweep_source.into(),
        sample_count: n_samples,
        inputs: ExportSchemaSection {
            schema: input_names
                .iter()
                .map(|name| ExportArraySchema {
                    name: name.clone(),
                    shape: vec![1],
                    units: "1".into(),
                    dtype: "f64".into(),
                    location: None,
                    mesh_id: None,
                })
                .collect(),
        },
        outputs: ExportSchemaSection {
            schema: config
                .output_names
                .iter()
                .zip(output_units.iter())
                .map(|(name, units)| ExportArraySchema {
                    name: name.clone(),
                    shape: vec![1],
                    units: units.clone(),
                    dtype: "f64".into(),
                    location: None,
                    mesh_id: None,
                })
                .collect(),
        },
        split: config.split.clone(),
        provenance: config.provenance.clone(),
    };
    write_manifest(&manifest, &out_dir.join("manifest.json")).map_err(DatasetExportError::Io)?;
    Ok(manifest)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_fields::scalar::ScalarSource;
    use valenx_fields::units::DIMENSIONLESS;
    use valenx_fields::{ScalarRecord, TimeKey};

    fn results_with_two_named_scalars(drag: f64, lift: f64) -> Results {
        use valenx_fields::provenance::Sha256Hex;
        let prov = valenx_fields::Provenance {
            adapter: "test".into(),
            adapter_version: "0".into(),
            tool: "Test".into(),
            tool_version: "0".into(),
            case_hash: Sha256Hex::new(""),
            mesh_hash: Sha256Hex::new(""),
            input_hash: Sha256Hex::new(""),
            tools_lock_hash: Sha256Hex::new(""),
            run_id: "00000000-0000-0000-0000-000000000000".into(),
            wall_time_seconds: 0.0,
            completed_at: "1970-01-01T00:00:00Z".into(),
            ancestors: Vec::new(),
        };
        let mut r = Results::empty("test", prov);
        r.scalars.insert(ScalarRecord {
            name: "drag".into(),
            value: drag,
            units: DIMENSIONLESS,
            time: TimeKey::Steady,
            source: ScalarSource::Extracted,
            description: None,
        });
        r.scalars.insert(ScalarRecord {
            name: "lift".into(),
            value: lift,
            units: DIMENSIONLESS,
            time: TimeKey::Steady,
            source: ScalarSource::Extracted,
            description: None,
        });
        r
    }

    fn fresh_out_dir(label: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "valenx-dataset-{label}-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }

    #[test]
    fn export_sweep_dataset_writes_inputs_outputs_ids_and_manifest() {
        let r0 = results_with_two_named_scalars(1.0, 5.0);
        let r1 = results_with_two_named_scalars(2.0, 6.0);
        let r2 = results_with_two_named_scalars(3.0, 7.0);
        let samples = vec![
            Sample {
                id: "sweep-000".into(),
                inputs: vec![("aoa".into(), 0.0), ("re".into(), 1e6)],
                outputs: &r0,
            },
            Sample {
                id: "sweep-001".into(),
                inputs: vec![("aoa".into(), 5.0), ("re".into(), 1e6)],
                outputs: &r1,
            },
            Sample {
                id: "sweep-002".into(),
                inputs: vec![("aoa".into(), 10.0), ("re".into(), 1e6)],
                outputs: &r2,
            },
        ];
        let cfg = DatasetExportConfig::from_output_names(["drag", "lift"]);
        let out = fresh_out_dir("happy");
        let manifest = export_sweep_dataset(&samples, &cfg, &out, "test-sweep").expect("export ok");

        assert_eq!(manifest.sample_count, 3);
        assert_eq!(manifest.inputs.schema.len(), 2);
        assert_eq!(manifest.outputs.schema.len(), 2);
        assert_eq!(manifest.inputs.schema[0].name, "aoa");
        assert_eq!(manifest.outputs.schema[1].name, "lift");

        // All four files landed.
        assert!(out.join("inputs.npy").is_file());
        assert!(out.join("outputs.npy").is_file());
        assert!(out.join("sample_ids.json").is_file());
        assert!(out.join("manifest.json").is_file());

        // sample_ids.json round-trips back to the same list.
        let ids: Vec<String> =
            serde_json::from_str(&std::fs::read_to_string(out.join("sample_ids.json")).unwrap())
                .unwrap();
        assert_eq!(ids, vec!["sweep-000", "sweep-001", "sweep-002"]);

        // outputs.npy has the right shape header and the lift values
        // pack in the right order.
        let bytes = std::fs::read(out.join("outputs.npy")).unwrap();
        let header_len = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
        let header = std::str::from_utf8(&bytes[10..10 + header_len]).unwrap();
        assert!(header.contains("'shape': (3, 2)"));
        let data_start = 6 + 2 + 2 + header_len;
        // Row 1, col 1 = sample[1].lift = 6.0.
        let v = f64::from_le_bytes(bytes[data_start + 24..data_start + 32].try_into().unwrap());
        assert_eq!(v, 6.0);

        let _ = std::fs::remove_dir_all(&out);
    }

    #[test]
    fn export_sweep_dataset_errors_when_a_sample_is_missing_a_scalar() {
        // Sample 1's Results doesn't contain "lift" — exporter must
        // fail loudly rather than silently emitting a NaN.
        let r0 = results_with_two_named_scalars(1.0, 5.0);
        let mut r1 = Results::empty(
            "test",
            valenx_fields::Provenance {
                adapter: "test".into(),
                adapter_version: "0".into(),
                tool: "Test".into(),
                tool_version: "0".into(),
                case_hash: valenx_fields::provenance::Sha256Hex::new(""),
                mesh_hash: valenx_fields::provenance::Sha256Hex::new(""),
                input_hash: valenx_fields::provenance::Sha256Hex::new(""),
                tools_lock_hash: valenx_fields::provenance::Sha256Hex::new(""),
                run_id: "00000000-0000-0000-0000-000000000000".into(),
                wall_time_seconds: 0.0,
                completed_at: "1970-01-01T00:00:00Z".into(),
                ancestors: Vec::new(),
            },
        );
        r1.scalars.insert(ScalarRecord {
            name: "drag".into(),
            value: 2.0,
            units: DIMENSIONLESS,
            time: TimeKey::Steady,
            source: ScalarSource::Extracted,
            description: None,
        });
        let samples = vec![
            Sample {
                id: "sweep-000".into(),
                inputs: vec![("aoa".into(), 0.0)],
                outputs: &r0,
            },
            Sample {
                id: "sweep-001".into(),
                inputs: vec![("aoa".into(), 5.0)],
                outputs: &r1,
            },
        ];
        let cfg = DatasetExportConfig::from_output_names(["drag", "lift"]);
        let out = fresh_out_dir("missing");
        let err = export_sweep_dataset(&samples, &cfg, &out, "test-sweep")
            .expect_err("missing-scalar must error");
        match err {
            DatasetExportError::MissingScalar { sample_id, name } => {
                assert_eq!(sample_id, "sweep-001");
                assert_eq!(name, "lift");
            }
            other => panic!("wrong error variant: {other:?}"),
        }
        let _ = std::fs::remove_dir_all(&out);
    }

    #[test]
    fn export_sweep_dataset_rejects_ragged_inputs() {
        let r0 = results_with_two_named_scalars(1.0, 5.0);
        let r1 = results_with_two_named_scalars(2.0, 6.0);
        let samples = vec![
            Sample {
                id: "sweep-000".into(),
                inputs: vec![("aoa".into(), 0.0), ("re".into(), 1e6)],
                outputs: &r0,
            },
            Sample {
                id: "sweep-001".into(),
                inputs: vec![("aoa".into(), 5.0)], // missing `re`
                outputs: &r1,
            },
        ];
        let cfg = DatasetExportConfig::from_output_names(["drag"]);
        let out = fresh_out_dir("ragged");
        let err =
            export_sweep_dataset(&samples, &cfg, &out, "ragged").expect_err("ragged must error");
        match err {
            DatasetExportError::RaggedInputs {
                sample_id,
                actual,
                expected,
            } => {
                assert_eq!(sample_id, "sweep-001");
                assert_eq!(actual, 1);
                assert_eq!(expected, 2);
            }
            other => panic!("wrong error variant: {other:?}"),
        }
        let _ = std::fs::remove_dir_all(&out);
    }

    #[test]
    fn export_sweep_dataset_handles_zero_samples_with_an_empty_manifest() {
        let cfg = DatasetExportConfig::from_output_names(["drag"]);
        let out = fresh_out_dir("empty");
        let manifest =
            export_sweep_dataset(&[], &cfg, &out, "empty-sweep").expect("empty must succeed");
        assert_eq!(manifest.sample_count, 0);
        assert!(out.join("manifest.json").is_file());
        // No data arrays for an empty dataset — loaders should detect
        // sample_count == 0 and skip rather than reading a stub file.
        assert!(!out.join("inputs.npy").exists());
        assert!(!out.join("outputs.npy").exists());
        let _ = std::fs::remove_dir_all(&out);
    }
}
