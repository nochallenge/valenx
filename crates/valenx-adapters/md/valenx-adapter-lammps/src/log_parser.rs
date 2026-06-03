//! Minimal `log.lammps` thermo-block parser.
//!
//! LAMMPS prints a header block + per-step thermo values to its
//! stdout, which the run captures as `log.lammps`. The block we care
//! about is bracketed by:
//!
//! ```text
//! Per MPI rank memory allocation (min/avg/max) = 3.124 | 3.124 | 3.124 Mbytes
//!    Step          Temp          PotEng         KinEng         …
//!         0   0              -6.7723705      0              …
//!       100   1.4985         …
//!       200   1.4972         …
//! Loop time of 0.123456 on 1 procs for 1000 steps with 4000 atoms
//! ```
//!
//! This module finds that block, reads the header to learn column
//! names, and emits one [`ThermoRow`] per data line — each row a
//! `step` plus a `BTreeMap<column_name, value>`.
//!
//! Multiple thermo blocks (one per `run` command in the deck) are
//! concatenated into a single time series.

use std::collections::BTreeMap;

/// One row of LAMMPS thermo output: a step and its named values.
#[derive(Clone, Debug)]
pub struct ThermoRow {
    pub step: u64,
    pub values: BTreeMap<String, f64>,
}

/// Per-column time series collected from every thermo block in a
/// `log.lammps`. `columns` preserves header-line order so callers
/// (chart axes, table widgets) can render in the same order LAMMPS
/// printed them — `BTreeMap` inside each `ThermoRow` would sort
/// alphabetically, which doesn't match what users expect.
#[derive(Clone, Debug, Default)]
pub struct ThermoSeries {
    pub rows: Vec<ThermoRow>,
    pub columns: Vec<String>,
}

impl ThermoSeries {
    /// Header column names in first-block first-encounter order.
    pub fn column_names(&self) -> Vec<&str> {
        self.columns.iter().map(|s| s.as_str()).collect()
    }
}

/// Parse a full `log.lammps` text into a [`ThermoSeries`]. Quietly
/// skips text that doesn't fit the thermo-block shape; LAMMPS logs
/// have a lot of header noise we don't care about.
pub fn parse_log(text: &str) -> ThermoSeries {
    let mut series = ThermoSeries::default();
    let mut lines = text.lines();
    while let Some(line) = lines.next() {
        if !line.contains("Per MPI rank memory allocation") {
            continue;
        }
        // Next non-empty line is the header (`Step Temp PotEng …`).
        let header_line = match lines.find(|l| !l.trim().is_empty()) {
            Some(h) => h,
            None => break,
        };
        let columns: Vec<String> = header_line
            .split_ascii_whitespace()
            .map(|s| s.to_string())
            .collect();
        if columns.is_empty() {
            continue;
        }
        // Record header order on the first block we see — subsequent
        // blocks usually repeat the same header but we don't override
        // an existing one.
        if series.columns.is_empty() {
            series.columns = columns.clone();
        }
        // The "Step" column is canonical — find its index so we know
        // where to read the step number from.
        let step_idx = columns.iter().position(|c| c.eq_ignore_ascii_case("Step"));

        // Read data rows until "Loop time" or another marker.
        for data_line in lines.by_ref() {
            let trimmed = data_line.trim();
            if trimmed.starts_with("Loop time")
                || trimmed.starts_with("Per MPI rank")
                || trimmed.starts_with("ERROR")
            {
                break;
            }
            if trimmed.is_empty() {
                continue;
            }
            let toks: Vec<&str> = trimmed.split_ascii_whitespace().collect();
            if toks.len() != columns.len() {
                // Malformed / partial line — skip rather than fail.
                continue;
            }
            // Parse every column as f64; if any fails we skip the
            // row (this filters out string columns we don't expect).
            let mut values: BTreeMap<String, f64> = BTreeMap::new();
            let mut ok = true;
            for (col, tok) in columns.iter().zip(toks.iter()) {
                match tok.parse::<f64>() {
                    Ok(v) => {
                        values.insert(col.clone(), v);
                    }
                    Err(_) => {
                        ok = false;
                        break;
                    }
                }
            }
            if !ok {
                continue;
            }
            let step = match step_idx
                .and_then(|i| toks.get(i))
                .and_then(|t| t.parse::<u64>().ok())
            {
                Some(s) => s,
                None => continue,
            };
            series.rows.push(ThermoRow { step, values });
        }
    }
    series
}

/// Convert a [`ThermoSeries`] into one canonical
/// [`valenx_fields::ScalarRecord`] per (column, step) pair. Each
/// record gets the step encoded as `TimeKey::Iteration(step)`.
///
/// Units are filled from a small lookup of LAMMPS thermo column names
/// (Temp → K, PotEng/KinEng/TotEng → kcal/mol or whatever the user's
/// `units` block declared, treated as dimensionless for now since
/// LAMMPS units are configurable and we don't track them yet).
/// Pressure → Pa, Volume → m³ in real units; in `units lj` they're
/// dimensionless. Marking everything dimensionless is the honest
/// answer until we wire the LAMMPS `units` setting through.
pub fn to_canonical_scalars(series: &ThermoSeries) -> Vec<valenx_fields::ScalarRecord> {
    let mut out: Vec<valenx_fields::ScalarRecord> = Vec::new();
    for row in &series.rows {
        for (name, value) in &row.values {
            // Skip the step column — it's the time key, not a value.
            if name.eq_ignore_ascii_case("Step") {
                continue;
            }
            out.push(valenx_fields::ScalarRecord {
                name: name.clone(),
                value: *value,
                units: valenx_fields::units::DIMENSIONLESS,
                time: valenx_fields::TimeKey::Iteration(row.step),
                source: valenx_fields::scalar::ScalarSource::Extracted,
                description: None,
            });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE_LOG: &str = r#"
LAMMPS (29 Aug 2024)
Reading data file...
  4000 atoms
  ...

Per MPI rank memory allocation (min/avg/max) = 3.124 | 3.124 | 3.124 Mbytes
   Step          Temp          PotEng         KinEng         TotEng         Press           Volume
         0   0              -6.7723705      0              -6.7723705     -2.6975796      4
       100   1.4985          -5.2123          1.4953          -3.7170     -0.5123          4
       200   1.4972          -5.2056          1.4937          -3.7119     -0.4987          4
Loop time of 0.123456 on 1 procs for 200 steps with 4000 atoms
"#;

    #[test]
    fn parses_one_thermo_block() {
        let series = parse_log(SAMPLE_LOG);
        assert_eq!(series.rows.len(), 3);
        assert_eq!(series.rows[0].step, 0);
        assert_eq!(series.rows[1].step, 100);
        assert_eq!(series.rows[2].step, 200);
        // Column extraction: Temp at step 100 = 1.4985.
        let t = series.rows[1].values.get("Temp").unwrap();
        assert!((t - 1.4985).abs() < 1e-12);
        // Volume column present.
        assert!(series.rows[0].values.contains_key("Volume"));
    }

    #[test]
    fn column_names_preserves_first_seen_order() {
        let series = parse_log(SAMPLE_LOG);
        let cols = series.column_names();
        // First three columns from the header order.
        assert_eq!(cols[0], "Step");
        assert_eq!(cols[1], "Temp");
        assert_eq!(cols[2], "PotEng");
    }

    #[test]
    fn to_canonical_scalars_emits_one_per_column_per_step() {
        let series = parse_log(SAMPLE_LOG);
        let scalars = to_canonical_scalars(&series);
        // 3 rows × 6 non-Step columns = 18 records.
        assert_eq!(scalars.len(), 18);
        // Time keys are step-indexed.
        let temp_records: Vec<_> = scalars.iter().filter(|s| s.name == "Temp").collect();
        assert_eq!(temp_records.len(), 3);
        assert_eq!(temp_records[0].time, valenx_fields::TimeKey::Iteration(0));
        assert_eq!(temp_records[1].time, valenx_fields::TimeKey::Iteration(100));
        assert_eq!(temp_records[2].time, valenx_fields::TimeKey::Iteration(200));
    }

    #[test]
    fn handles_multiple_thermo_blocks() {
        // Two `run` commands → two thermo blocks → concatenated rows.
        let text = r#"
Per MPI rank memory allocation (min/avg/max) = X
   Step          Temp
         0   0
       100   1.5
Loop time of 0.1 on 1 procs

Per MPI rank memory allocation (min/avg/max) = X
   Step          Temp
       200   1.6
       300   1.7
Loop time of 0.1 on 1 procs
"#;
        let series = parse_log(text);
        assert_eq!(series.rows.len(), 4);
        assert_eq!(series.rows[3].step, 300);
    }

    #[test]
    fn empty_log_returns_empty_series() {
        let series = parse_log("");
        assert!(series.rows.is_empty());
        assert!(series.column_names().is_empty());
    }

    #[test]
    fn malformed_data_rows_are_skipped_not_fatal() {
        let text = r#"
Per MPI rank memory allocation (min/avg/max) = X
   Step          Temp          PotEng
         0   0              -6.7
       100   not-a-number    -5.0
       200   1.5             -4.5
Loop time of 0.1 on 1 procs
"#;
        let series = parse_log(text);
        // Two valid rows — the malformed middle row got skipped.
        assert_eq!(series.rows.len(), 2);
        assert_eq!(series.rows[0].step, 0);
        assert_eq!(series.rows[1].step, 200);
    }
}
