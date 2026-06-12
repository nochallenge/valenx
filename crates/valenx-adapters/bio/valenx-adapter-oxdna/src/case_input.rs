//! `[bio.oxdna]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "oxdna.batch"
//!
//! [bio.oxdna]
//! input    = "input.dat"           # required
//! topology = "duplex.top"          # optional, default reads from input.dat
//! ```
//!
//! oxDNA reads everything (initial conf, topology, force-field
//! flags, integration parameters) from a single `input.dat` file. We
//! still surface the optional `topology` key because oxDNA resolves
//! the `.top` path relative to its working directory, and staging
//! the topology explicitly into the workdir avoids surprises when
//! the user's case directory differs from the adapter's workdir.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct OxDnaInput {
    pub input: PathBuf,
    /// Optional explicit topology path. When set, the adapter
    /// stages this file alongside `input.dat`. When unset, oxDNA
    /// resolves the topology path itself from the `topology = ...`
    /// line inside `input.dat`.
    pub topology: Option<PathBuf>,
}

impl OxDnaInput {
    pub fn from_case_dir(case_dir: &std::path::Path) -> Result<Self, AdapterError> {
        let case_toml = case_dir.join("case.toml");
        let text = valenx_core::io_caps::read_capped_to_string(
            &case_toml,
            valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
        )
        .map_err(|e| AdapterError::Other(anyhow::anyhow!("read {}: {e}", case_toml.display())))?;
        let parsed: toml::Value = toml::from_str(&text).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("parse {}: {e}", case_toml.display()))
        })?;
        let block = parsed
            .get("bio")
            .and_then(|v| v.get("oxdna"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.oxdna] section",
                    case_toml.display()
                ))
            })?;
        let input = block
            .get("input")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.oxdna].input required")))?;
        let topology = block
            .get("topology")
            .and_then(|v| v.as_str())
            .map(PathBuf::from);
        Ok(Self {
            input: PathBuf::from(input),
            topology,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal_case() {
        let d = tempdir("oxdna");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "oxdna.batch"

[bio.oxdna]
input = "input.dat"
"#,
        )
        .unwrap();
        let input = OxDnaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("input.dat"));
        // Topology defaults to None — oxDNA reads it from input.dat.
        assert!(input.topology.is_none());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_missing_section() {
        let d = tempdir("oxdna");
        std::fs::write(
            d.join("case.toml"),
            "[case]\nphysics=\"bio\"\nsolver=\"x\"\n",
        )
        .unwrap();
        let err = OxDnaInput::from_case_dir(&d).unwrap_err();
        assert!(format!("{err}").contains("[bio.oxdna]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn honours_explicit_topology_override() {
        // A duplex case with an explicit topology — verify both
        // input + topology round-trip.
        let d = tempdir("oxdna");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "oxdna.batch"

[bio.oxdna]
input    = "input.dat"
topology = "duplex.top"
"#,
        )
        .unwrap();
        let input = OxDnaInput::from_case_dir(&d).unwrap();
        assert_eq!(input.input, PathBuf::from("input.dat"));
        assert_eq!(input.topology, Some(PathBuf::from("duplex.top")));
        let _ = std::fs::remove_dir_all(&d);
    }
}
