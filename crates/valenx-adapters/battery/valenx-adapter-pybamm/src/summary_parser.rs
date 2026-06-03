//! Parse PyBaMM's `summary.json`.

use std::path::Path;

use serde::{Deserialize, Serialize};

use valenx_core::AdapterError;

/// Hard cap on summary.json bytes read into memory.
///
/// Round-19 M3: pre-fix `parse_file` slurped the file via
/// `fs::read_to_string` unbounded — a poisoned workdir with a
/// multi-GB `summary.json` would OOM the renderer before serde_json
/// saw the first token. 1 MiB is generous (real pybamm summary.json
/// files are typically a few KB to ~10 KB — `samples` + scalar
/// fields, never time-series) while staying well below a hostile
/// `cat /dev/zero > summary.json` DoS.
pub const MAX_SUMMARY_FILE_BYTES: usize = 1024 * 1024;

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct PyBammSummary {
    #[serde(default)]
    pub valenx_adapter: String,
    #[serde(default)]
    pub model: String,
    #[serde(default)]
    pub parameter_set: String,
    #[serde(default)]
    pub samples: u64,
    #[serde(default)]
    pub duration_s: f64,
    #[serde(default)]
    pub voltage_start_v: f64,
    #[serde(default)]
    pub voltage_end_v: f64,
    #[serde(default)]
    pub initial_soc: f64,
}

pub fn parse_file(path: &Path) -> Result<PyBammSummary, AdapterError> {
    // Round-19 M3: replace the unbounded `fs::read_to_string` with the
    // shared cap helper so a poisoned summary.json can't slurp into
    // memory before serde_json sees it. 1 MiB is generous — real
    // pybamm summary files are <10 KiB in practice.
    let text = valenx_core::io_caps::read_capped_to_string(path, MAX_SUMMARY_FILE_BYTES)?;
    serde_json::from_str::<PyBammSummary>(&text).map_err(|e| AdapterError::ParseOutput {
        file: path.to_path_buf(),
        reason: format!("pybamm summary.json: {e}"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_full_summary() {
        let json = r#"{
            "valenx_adapter": "pybamm",
            "model": "DFN",
            "parameter_set": "Chen2020",
            "samples": 120,
            "duration_s": 3600.0,
            "voltage_start_v": 4.19,
            "voltage_end_v": 2.8,
            "initial_soc": 1.0
        }"#;
        let s: PyBammSummary = serde_json::from_str(json).unwrap();
        assert_eq!(s.samples, 120);
        assert_eq!(s.model, "DFN");
        assert!((s.voltage_start_v - 4.19).abs() < 1e-6);
    }

    #[test]
    fn tolerates_partial() {
        let json = r#"{ "valenx_adapter": "pybamm" }"#;
        let s: PyBammSummary = serde_json::from_str(json).unwrap();
        assert_eq!(s.valenx_adapter, "pybamm");
        assert_eq!(s.samples, 0);
    }

    /// Round-19 M3 RED→GREEN: a `summary.json` whose advertised size
    /// exceeds `MAX_SUMMARY_FILE_BYTES` must produce an error rather
    /// than slurp the file into memory. Pre-fix `parse_file` used
    /// `fs::read_to_string` unbounded — a poisoned 50 MiB summary
    /// would allocate 50 MiB of RAM before serde_json saw it.
    ///
    /// Sparse-file trick: `set_len(cap+1)` advertises a size past
    /// the cap without writing the bytes to disk; the cap helper
    /// rejects on stat alone, before any read happens. (NTFS / ext4
    /// / APFS all support sparse files, so the test stays in the
    /// millisecond range.)
    #[test]
    fn parse_file_rejects_oversize_summary() {
        let p = std::env::temp_dir().join(format!(
            "valenx-pybamm-summary-toolarge-{}.json",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let f = std::fs::File::create(&p).unwrap();
        f.set_len((MAX_SUMMARY_FILE_BYTES + 1) as u64).unwrap();
        drop(f);
        let err = parse_file(&p).expect_err("oversize summary must error");
        let _ = std::fs::remove_file(&p);
        // The cap helper surfaces this as `AdapterError::Io` with
        // InvalidData kind — both variants are acceptable; assert
        // we got *some* error rather than silently succeeding.
        match err {
            AdapterError::Io(io_err) => {
                assert_eq!(
                    io_err.kind(),
                    std::io::ErrorKind::InvalidData,
                    "expected InvalidData kind for oversize cap, got: {io_err:?}"
                );
            }
            other => panic!("wrong error variant: {other:?}"),
        }
    }
}
