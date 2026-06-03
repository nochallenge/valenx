//! CSV writer for `Results.scalars` catalogs.
//!
//! Pulled out of `lib.rs` so the CSV-specific helpers (`format_csv_row`,
//! `csv_quote`, `format_time_key`) live next to the public writer that
//! uses them.

use std::fmt::Write as _;
use std::io::Write as _;
use std::path::Path;

use valenx_fields::{Results, ScalarRecord, TimeKey};

use crate::ExportError;

/// Write the `Results.scalars` catalog as a CSV file. Columns:
///
/// 1. `name` — record name as the adapter wrote it.
/// 2. `value` — `f64` formatted via `{:.10e}` so precision survives
///    the round-trip into Excel / pandas without rounding.
/// 3. `units` — short symbol from `Units::display`, or empty.
/// 4. `time_key` — encoded as `steady` / `iter:<N>` / `time:<seconds>`
///    so a downstream loader can reconstruct the variant.
/// 5. `source` — `Computed` / `Extracted` / `UserDefined`.
pub fn write_scalars_csv(results: &Results, path: &Path) -> Result<(), ExportError> {
    // R30: stream the catalog through the crash-safe *streaming* atomic
    // writer (sidecar → fsync → rename). `atomic_write_streaming` hands
    // the closure a `&mut BufWriter<File>`, so a large catalog is never
    // buffered whole in memory — yet a torn/concurrent write can't leave
    // a half-written CSV a downstream loader would mis-parse.
    valenx_core::io_caps::atomic_write_streaming(path, |f| {
        writeln!(f, "name,value,units,time_key,source")?;
        for name in results.scalars.names() {
            for record in results.scalars.all(name) {
                writeln!(f, "{}", format_csv_row(record))?;
            }
        }
        Ok(())
    })
    .map_err(|e| ExportError::Io {
        path: path.to_path_buf(),
        source: e,
    })
}

/// Render one ScalarRecord as a CSV row. Pulled out of the writer
/// for testability — the main writer just emits one of these per
/// record.
pub fn format_csv_row(record: &ScalarRecord) -> String {
    let units_label = record.units.display.unwrap_or("");
    let time_key = format_time_key(record.time);
    let source = format!("{:?}", record.source);
    // Quote any value that contains commas / quotes / newlines.
    let name = csv_quote(&record.name);
    let units = csv_quote(units_label);
    let mut row = String::with_capacity(96);
    let _ = write!(
        &mut row,
        "{name},{value:.10e},{units},{time_key},{source}",
        value = record.value,
    );
    row
}

fn format_time_key(time: TimeKey) -> String {
    match time {
        TimeKey::Steady => "steady".to_string(),
        TimeKey::Iteration(n) => format!("iter:{n}"),
        TimeKey::Time { value, .. } => format!("time:{value}"),
    }
}

fn csv_quote(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') {
        let escaped = s.replace('"', "\"\"");
        format!("\"{escaped}\"")
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_fields::scalar::ScalarSource;
    use valenx_fields::units::{Units, DIMENSIONLESS, SECOND};

    fn pa_units() -> Units {
        Units::new([-1, 1, -2, 0, 0, 0, 0], 1.0, Some("Pa"))
    }

    fn results_with_a_few_scalars() -> Results {
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
            name: "drag_coefficient".into(),
            value: 0.123456,
            units: DIMENSIONLESS,
            time: TimeKey::Steady,
            source: ScalarSource::Extracted,
            description: None,
        });
        r.scalars.insert(ScalarRecord {
            name: "p_inlet".into(),
            value: 101325.0,
            units: pa_units(),
            time: TimeKey::Iteration(500),
            source: ScalarSource::Extracted,
            description: None,
        });
        r.scalars.insert(ScalarRecord {
            name: "T_at_t1ms".into(),
            value: 293.15,
            units: DIMENSIONLESS,
            time: TimeKey::Time {
                value: 0.001,
                units: SECOND,
            },
            source: ScalarSource::Computed,
            description: None,
        });
        r
    }

    #[test]
    fn csv_writer_produces_header_and_rows() {
        let results = results_with_a_few_scalars();
        let tmp = std::env::temp_dir().join(format!(
            "valenx-export-{}.csv",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        write_scalars_csv(&results, &tmp).expect("write");
        let text = std::fs::read_to_string(&tmp).expect("read");
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "name,value,units,time_key,source");
        // Three records → three data rows.
        assert_eq!(lines.len(), 4);
        assert!(text.contains("drag_coefficient"));
        assert!(text.contains("Pa"));
        assert!(text.contains("iter:500"));
        assert!(text.contains("time:0.001"));
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn time_key_encoding_round_trip_handles_each_variant() {
        assert_eq!(format_time_key(TimeKey::Steady), "steady");
        assert_eq!(format_time_key(TimeKey::Iteration(500)), "iter:500");
        let t = TimeKey::Time {
            value: 0.005,
            units: SECOND,
        };
        let s = format_time_key(t);
        assert!(s.starts_with("time:0.005"), "got {s}");
    }

    #[test]
    fn csv_quote_escapes_commas_and_quotes() {
        assert_eq!(csv_quote("plain"), "plain");
        assert_eq!(csv_quote("has,comma"), "\"has,comma\"");
        assert_eq!(csv_quote("has \"quotes\""), "\"has \"\"quotes\"\"\"");
        assert_eq!(csv_quote(""), "");
    }

    #[test]
    fn format_csv_row_uses_scientific_notation_for_value() {
        let r = ScalarRecord {
            name: "p".into(),
            value: 101325.0,
            units: pa_units(),
            time: TimeKey::Steady,
            source: ScalarSource::Extracted,
            description: None,
        };
        let row = format_csv_row(&r);
        assert!(row.starts_with("p,1.0132500000e5"));
        assert!(row.ends_with(",Pa,steady,Extracted"));
    }
}
