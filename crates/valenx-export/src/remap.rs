//! Field-output remap helper for geometry-varying sweeps.
//!
//! Bridges to `valenx_fields::interp::nearest_neighbour_remap` so the
//! sweep-dataset exporter can pack per-sample mesh fields into one
//! aligned tensor over a shared reference mesh.

use thiserror::Error;

/// One sample's mesh + the field arrays the export pipeline wants
/// remapped onto a shared reference. Bundles the per-sample inputs
/// + sample id so the output rows stay aligned to the sample list.
pub struct FieldSample<'a> {
    pub id: String,
    pub points: &'a [nalgebra::Vector3<f64>],
    /// Each named scalar field defined on `points`. Vector / tensor
    /// outputs are deferred to follow-up commits — the v0 helper
    /// handles OnNode scalars only (the most common case for
    /// surrogate-model inputs).
    pub fields: Vec<&'a valenx_fields::Field>,
}

/// Remap every sample's listed field onto a shared reference mesh
/// via [`valenx_fields::interp::nearest_neighbour_remap`]. Returns
/// one remapped Field per (sample, field name) pair — flattened so
/// downstream packing into npy can iterate in (n_samples,
/// n_field_names, n_reference_nodes) order.
///
/// All samples must declare fields in the same order; inconsistent
/// orderings produce a structured error rather than silently
/// shuffling columns.
pub fn remap_sample_fields(
    samples: &[FieldSample<'_>],
    reference_points: &[nalgebra::Vector3<f64>],
) -> Result<Vec<Vec<valenx_fields::Field>>, FieldRemapError> {
    if samples.is_empty() {
        return Ok(Vec::new());
    }
    let first = &samples[0];
    let n_fields = first.fields.len();
    let names: Vec<String> = first.fields.iter().map(|f| f.name.clone()).collect();

    // Sanity-check the rest of the samples match the field schema.
    for s in samples.iter().skip(1) {
        if s.fields.len() != n_fields {
            return Err(FieldRemapError::FieldCountMismatch {
                sample_id: s.id.clone(),
                expected: n_fields,
                got: s.fields.len(),
            });
        }
        for (idx, (expected, actual)) in names.iter().zip(s.fields.iter()).enumerate() {
            if expected != &actual.name {
                return Err(FieldRemapError::FieldNameMismatch {
                    sample_id: s.id.clone(),
                    index: idx,
                    expected: expected.clone(),
                    actual: actual.name.clone(),
                });
            }
        }
    }

    let mut out: Vec<Vec<valenx_fields::Field>> = Vec::with_capacity(samples.len());
    for s in samples {
        let mut row: Vec<valenx_fields::Field> = Vec::with_capacity(n_fields);
        for f in &s.fields {
            let sample = valenx_fields::interp::SampleField {
                points: s.points,
                field: f,
            };
            let remapped =
                valenx_fields::interp::nearest_neighbour_remap(reference_points, &sample).map_err(
                    |e| FieldRemapError::Interp {
                        sample_id: s.id.clone(),
                        source: e,
                    },
                )?;
            row.push(remapped);
        }
        out.push(row);
    }
    Ok(out)
}

/// Errors raised by the field-remap helpers.
#[derive(Debug, Error)]
pub enum FieldRemapError {
    /// A later sample has a different number of fields than the first.
    #[error("sample `{sample_id}` declared {got} fields; first sample declared {expected}")]
    FieldCountMismatch {
        /// Identifier of the offending sample.
        sample_id: String,
        /// Field count seen on the first (reference) sample.
        expected: usize,
        /// Field count seen on the offending sample.
        got: usize,
    },
    /// A later sample names one of its fields differently from the
    /// first sample at the same index.
    #[error(
        "sample `{sample_id}` field #{index} is named `{actual}`, but the first sample named it `{expected}`"
    )]
    FieldNameMismatch {
        /// Identifier of the offending sample.
        sample_id: String,
        /// 0-based field index.
        index: usize,
        /// Field name on the reference sample.
        expected: String,
        /// Field name on the offending sample.
        actual: String,
    },
    /// Interpolation of one sample's field failed (see the wrapped
    /// [`valenx_fields::interp::InterpError`] for the precise reason).
    #[error("sample `{sample_id}` interpolation failed: {source}")]
    Interp {
        /// Identifier of the offending sample.
        sample_id: String,
        /// Underlying interpolation error.
        #[source]
        source: valenx_fields::interp::InterpError,
    },
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_fields::TimeKey;

    fn make_scalar_field(name: &str, data: Vec<f64>) -> valenx_fields::Field {
        valenx_fields::Field {
            name: name.into(),
            kind: valenx_fields::FieldKind::Scalar,
            location: valenx_fields::Location::OnNode,
            region: valenx_fields::RegionRef("default".into()),
            units: valenx_fields::units::DIMENSIONLESS,
            time: TimeKey::Steady,
            data,
            range: None,
        }
    }

    #[test]
    fn remap_sample_fields_returns_one_row_per_sample() {
        use nalgebra::Vector3;
        let pts0 = vec![Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0)];
        let pts1 = vec![Vector3::new(0.5, 0.0, 0.0), Vector3::new(1.5, 0.0, 0.0)];
        let f0 = make_scalar_field("p", vec![10.0, 20.0]);
        let f1 = make_scalar_field("p", vec![100.0, 200.0]);
        let samples = vec![
            FieldSample {
                id: "s0".into(),
                points: &pts0,
                fields: vec![&f0],
            },
            FieldSample {
                id: "s1".into(),
                points: &pts1,
                fields: vec![&f1],
            },
        ];
        let ref_pts = vec![Vector3::new(0.4, 0.0, 0.0), Vector3::new(1.1, 0.0, 0.0)];
        let out = remap_sample_fields(&samples, &ref_pts).expect("remap");
        assert_eq!(out.len(), 2);
        // Sample 0: pts0=[0,1]; ref 0.4 nearest pts0[0]=0 (10.0);
        // 1.1 nearest pts0[1]=1 (20.0).
        assert_eq!(out[0][0].data, vec![10.0, 20.0]);
        // Sample 1: pts1=[0.5,1.5]; ref 0.4 nearest pts1[0]=0.5
        // (100.0); 1.1 nearest pts1[1]=1.5 (200.0).
        assert_eq!(out[1][0].data, vec![100.0, 200.0]);
    }

    #[test]
    fn remap_sample_fields_rejects_field_count_mismatch() {
        use nalgebra::Vector3;
        let pts = vec![Vector3::new(0.0, 0.0, 0.0)];
        let f0 = make_scalar_field("p", vec![1.0]);
        let f1a = make_scalar_field("p", vec![2.0]);
        let f1b = make_scalar_field("T", vec![300.0]);
        let samples = vec![
            FieldSample {
                id: "s0".into(),
                points: &pts,
                fields: vec![&f0],
            },
            FieldSample {
                id: "s1".into(),
                points: &pts,
                fields: vec![&f1a, &f1b],
            },
        ];
        let ref_pts = vec![Vector3::new(0.0, 0.0, 0.0)];
        let err = remap_sample_fields(&samples, &ref_pts).unwrap_err();
        match err {
            FieldRemapError::FieldCountMismatch {
                sample_id,
                expected,
                got,
            } => {
                assert_eq!(sample_id, "s1");
                assert_eq!(expected, 1);
                assert_eq!(got, 2);
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn remap_sample_fields_rejects_field_name_mismatch() {
        use nalgebra::Vector3;
        let pts = vec![Vector3::new(0.0, 0.0, 0.0)];
        let f0 = make_scalar_field("p", vec![1.0]);
        let f1 = make_scalar_field("T", vec![2.0]);
        let samples = vec![
            FieldSample {
                id: "s0".into(),
                points: &pts,
                fields: vec![&f0],
            },
            FieldSample {
                id: "s1".into(),
                points: &pts,
                fields: vec![&f1],
            },
        ];
        let ref_pts = vec![Vector3::new(0.0, 0.0, 0.0)];
        let err = remap_sample_fields(&samples, &ref_pts).unwrap_err();
        match err {
            FieldRemapError::FieldNameMismatch {
                sample_id,
                expected,
                actual,
                ..
            } => {
                assert_eq!(sample_id, "s1");
                assert_eq!(expected, "p");
                assert_eq!(actual, "T");
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn remap_sample_fields_handles_empty_input() {
        use nalgebra::Vector3;
        let ref_pts = vec![Vector3::new(0.0, 0.0, 0.0)];
        let out = remap_sample_fields(&[], &ref_pts).expect("empty -> empty");
        assert!(out.is_empty());
    }

    #[test]
    fn remap_sample_fields_propagates_interp_errors() {
        use nalgebra::Vector3;
        // Empty reference triggers InterpError::EmptyReference; the
        // helper must propagate it under a FieldRemapError::Interp.
        let pts = vec![Vector3::new(0.0, 0.0, 0.0)];
        let f0 = make_scalar_field("p", vec![1.0]);
        let samples = vec![FieldSample {
            id: "s0".into(),
            points: &pts,
            fields: vec![&f0],
        }];
        let err = remap_sample_fields(&samples, &[]).unwrap_err();
        match err {
            FieldRemapError::Interp { sample_id, .. } => assert_eq!(sample_id, "s0"),
            other => panic!("wrong error: {other:?}"),
        }
    }
}
