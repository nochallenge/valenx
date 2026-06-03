//! Parse `summary.json` emitted by the MuJoCo Python script.

use std::path::Path;

use serde::{Deserialize, Serialize};

use valenx_core::AdapterError;

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct MuJoCoSummary {
    #[serde(default)]
    pub valenx_adapter: String,
    #[serde(default)]
    pub model: String,
    #[serde(default)]
    pub duration_s: f64,
    #[serde(default)]
    pub timestep_s: f64,
    #[serde(default)]
    pub step_count: u64,
    #[serde(default)]
    pub nq: u32,
    #[serde(default)]
    pub nv: u32,
    #[serde(default)]
    pub nu: u32,
}

pub fn parse_file(path: &Path) -> Result<MuJoCoSummary, AdapterError> {
    // Round-23 named finding: bound the summary read at
    // MAX_MUJOCO_SUMMARY_BYTES (1 MiB) — MuJoCo summaries record
    // only model metadata (always < 10 KiB) so 1 MiB is generous
    // while refusing poisoned/runaway summaries.
    let text = valenx_core::io_caps::read_capped_to_string(
        path,
        valenx_core::io_caps::MAX_MUJOCO_SUMMARY_BYTES as usize,
    )?;
    serde_json::from_str::<MuJoCoSummary>(&text).map_err(|e| AdapterError::ParseOutput {
        file: path.to_path_buf(),
        reason: format!("mujoco summary.json: {e}"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_summary() {
        let json = r#"{
            "valenx_adapter": "mujoco",
            "model": "pendulum.xml",
            "duration_s": 5.0,
            "timestep_s": 0.002,
            "step_count": 2500,
            "nq": 1, "nv": 1, "nu": 1
        }"#;
        let s: MuJoCoSummary = serde_json::from_str(json).unwrap();
        assert_eq!(s.step_count, 2500);
        assert_eq!(s.nq, 1);
        assert!((s.timestep_s - 0.002).abs() < 1e-9);
    }

    #[test]
    fn partial_still_parses() {
        let s: MuJoCoSummary = serde_json::from_str(r#"{"valenx_adapter":"mujoco"}"#).unwrap();
        assert_eq!(s.step_count, 0);
    }

    /// Round-23 named finding RED→GREEN: `parse_file` rejects an
    /// over-cap summary at the read-cap layer.
    #[test]
    fn parse_file_rejects_oversize() {
        let tmp = std::env::temp_dir().join("valenx_mujoco_summary_oversize.json");
        let oversize_bytes = (valenx_core::io_caps::MAX_MUJOCO_SUMMARY_BYTES as usize) + 1024;
        std::fs::write(&tmp, vec![b'x'; oversize_bytes]).unwrap();
        let err = parse_file(&tmp).expect_err("must reject oversize");
        let msg = format!("{err}");
        assert!(
            msg.contains("exceeds") || msg.contains("cap"),
            "expected cap-exceeded msg, got: {msg}"
        );
        let _ = std::fs::remove_file(&tmp);
    }
}
