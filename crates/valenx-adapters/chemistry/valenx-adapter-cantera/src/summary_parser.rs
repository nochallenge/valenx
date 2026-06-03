//! Parse Cantera's `summary.json` — the output of the generated
//! Python script — into a typed [`CanteraSummary`].

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

use valenx_core::AdapterError;

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct CanteraSummary {
    #[serde(default)]
    pub valenx_adapter: String,
    #[serde(default)]
    pub analysis: String,
    #[serde(default)]
    pub mechanism: String,
    #[serde(default)]
    pub initial: Option<ThermoFrame>,
    #[serde(default)]
    pub final_: Option<ThermoFrame>,
    #[serde(default)]
    pub mole_fractions: BTreeMap<String, f64>,
    #[serde(default)]
    pub mean_molecular_weight: Option<f64>,
    #[serde(default)]
    pub species_count_kept: Option<u64>,
    #[serde(default)]
    pub species_count_total: Option<u64>,
}

/// The JSON keyword `final` is reserved in Python but not in JSON —
/// serde lets us rename the field to avoid shadowing Rust's `final`
/// keyword (which is reserved-for-future-use).
impl CanteraSummary {
    /// Deserialisation hook that maps JSON's `final` key onto the
    /// renamed `final_` field. serde_json handles this via the
    /// custom `Deserialize` impl we attach below.
    pub fn parse_str(text: &str) -> Result<Self, serde_json::Error> {
        // Use an untyped Value + manual rename so we don't need a
        // custom Deserialize for every field.
        let mut v: serde_json::Value = serde_json::from_str(text)?;
        if let Some(obj) = v.as_object_mut() {
            if let Some(final_value) = obj.remove("final") {
                obj.insert("final_".to_string(), final_value);
            }
        }
        serde_json::from_value(v)
    }
}

#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize)]
pub struct ThermoFrame {
    #[serde(rename = "T")]
    pub temperature_k: f64,
    #[serde(rename = "P")]
    pub pressure_pa: f64,
    #[serde(rename = "H")]
    pub enthalpy_mass: f64,
    #[serde(rename = "S")]
    pub entropy_mass: f64,
    #[serde(rename = "rho")]
    pub density: f64,
}

pub fn parse_file(path: &Path) -> Result<CanteraSummary, AdapterError> {
    // Round-23 named finding: bound the summary read at
    // MAX_CANTERA_SUMMARY_BYTES (1 MiB) — Cantera summaries are
    // tiny (adapter ID + analysis tag + handful of thermo frames)
    // so 1 MiB is generous while refusing poisoned/runaway summaries.
    let text = valenx_core::io_caps::read_capped_to_string(
        path,
        valenx_core::io_caps::MAX_CANTERA_SUMMARY_BYTES as usize,
    )?;
    CanteraSummary::parse_str(&text).map_err(|e| AdapterError::ParseOutput {
        file: path.to_path_buf(),
        reason: format!("cantera summary.json: {e}"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_full_summary() {
        let json = r#"{
          "valenx_adapter": "cantera",
          "analysis": "HP",
          "mechanism": "gri30.yaml",
          "initial": {"T": 300.0, "P": 101325.0, "H": -72000.0, "S": 6900.0, "rho": 1.177},
          "final":   {"T": 2226.0, "P": 101325.0, "H": -72000.0, "S": 8900.0, "rho": 0.160},
          "mole_fractions": {"CO2": 0.0982, "H2O": 0.1865, "N2": 0.7108},
          "mean_molecular_weight": 28.5,
          "species_count_kept": 3,
          "species_count_total": 53
        }"#;
        let s = CanteraSummary::parse_str(json).unwrap();
        assert_eq!(s.analysis, "HP");
        assert_eq!(s.mole_fractions.len(), 3);
        assert!((s.mole_fractions["CO2"] - 0.0982).abs() < 1e-6);
        assert!((s.final_.unwrap().temperature_k - 2226.0).abs() < 0.1);
        assert_eq!(s.species_count_total, Some(53));
    }

    #[test]
    fn tolerates_partial_summary() {
        let json = r#"{ "valenx_adapter": "cantera" }"#;
        let s = CanteraSummary::parse_str(json).unwrap();
        assert!(s.mole_fractions.is_empty());
        assert!(s.initial.is_none());
    }

    #[test]
    fn malformed_json_errors() {
        let err = CanteraSummary::parse_str("nope").unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("expected") || msg.contains("invalid"));
    }

    /// Round-23 named finding RED→GREEN: `parse_file` rejects an
    /// over-cap summary at the read-cap layer rather than slurping
    /// it into memory. We write a 2 MiB file with 1 MiB cap to
    /// trigger the bound.
    #[test]
    fn parse_file_rejects_oversize() {
        let tmp = std::env::temp_dir().join("valenx_cantera_summary_oversize.json");
        let oversize_bytes = (valenx_core::io_caps::MAX_CANTERA_SUMMARY_BYTES as usize) + 1024;
        std::fs::write(&tmp, vec![b'x'; oversize_bytes]).unwrap();
        let err = parse_file(&tmp).expect_err("must reject oversize");
        // The Io variant from valenx_core::io_caps::read_capped_to_string
        // wraps an InvalidData std::io::Error; AdapterError's #[from]
        // From<std::io::Error> impl turns it into AdapterError::Io.
        let msg = format!("{err}");
        assert!(
            msg.contains("exceeds") || msg.contains("cap"),
            "expected cap-exceeded msg, got: {msg}"
        );
        let _ = std::fs::remove_file(&tmp);
    }
}
