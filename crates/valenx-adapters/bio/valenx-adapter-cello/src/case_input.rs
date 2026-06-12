//! `[bio.cello]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "cello.compile"
//!
//! [bio.cello]
//! jar               = "/opt/cello/cello.jar"
//! verilog           = "circuit.v"
//! user_constraints  = "Eco1C1G1T1.UCF.json"
//! input_sensors     = "Eco1C1G1T1.input.json"
//! output_devices    = "Eco1C1G1T1.output.json"
//! output_basename   = "circuit_out"
//! extra_args        = ["-logFilename", "cello.log"]   # optional, defaults to []
//! ```
//!
//! Cello is the canonical genetic-circuit DNA compiler — it
//! consumes a Verilog netlist describing the desired logic
//! function, plus a triplet of JSON constraint files (a user
//! constraint file pinning the chassis / library, an input sensor
//! file pinning the input promoters, an output device file pinning
//! the reporter), and emits a fully assembled DNA construct that
//! implements the logic in a living cell. Outputs include a
//! human-readable text report, a circuit diagram PNG, and a
//! Graphviz `.dot` netlist of the synthesized circuit.
//!
//! Like j5, Cello v2 is JAR-distributed — there's no `cello`
//! launcher binary on PATH. The user supplies the jar path via
//! `jar`; we probe `java` itself.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct CelloInput {
    /// Absolute path to the Cello v2 jar.
    pub jar: PathBuf,
    /// Path to the Verilog netlist describing the target logic
    /// function. Cello reads it via `-inputNetlist <verilog>`.
    pub verilog: PathBuf,
    /// Path to the user constraint file (JSON) pinning the chassis
    /// / library of available gates.
    pub user_constraints: PathBuf,
    /// Path to the input sensor file (JSON) pinning the input
    /// promoters Cello can wire to the circuit's primary inputs.
    pub input_sensors: PathBuf,
    /// Path to the output device file (JSON) pinning the reporter
    /// genes Cello can wire to the circuit's primary outputs.
    pub output_devices: PathBuf,
    /// Filename stem for outputs. Cello writes results into a
    /// directory it creates under the workdir, prefixed by this
    /// basename. Anchors `collect()`'s artefact filter.
    pub output_basename: String,
    /// Additional CLI arguments appended to the `java -jar` call —
    /// useful for `-logFilename`, `-options`, or pinning a custom
    /// log directory.
    pub extra_args: Vec<String>,
}

impl CelloInput {
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
            .and_then(|v| v.get("cello"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.cello] section",
                    case_toml.display()
                ))
            })?;

        let jar = block.get("jar").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.cello].jar required (path to cello.jar)"
            ))
        })?;
        if jar.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cello].jar must not be empty"
            )));
        }

        let verilog = block
            .get("verilog")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.cello].verilog required (path to .v netlist)"
                ))
            })?;
        if verilog.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cello].verilog must not be empty"
            )));
        }

        let user_constraints = block
            .get("user_constraints")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.cello].user_constraints required (path to UCF JSON)"
                ))
            })?;
        if user_constraints.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cello].user_constraints must not be empty"
            )));
        }

        let input_sensors = block
            .get("input_sensors")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.cello].input_sensors required (path to input sensor JSON)"
                ))
            })?;
        if input_sensors.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cello].input_sensors must not be empty"
            )));
        }

        let output_devices = block
            .get("output_devices")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.cello].output_devices required (path to output device JSON)"
                ))
            })?;
        if output_devices.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cello].output_devices must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!("[bio.cello].output_basename required"))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.cello].output_basename must not be empty"
            )));
        }

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.cello].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.cello].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            jar: PathBuf::from(jar),
            verilog: PathBuf::from(verilog),
            user_constraints: PathBuf::from(user_constraints),
            input_sensors: PathBuf::from(input_sensors),
            output_devices: PathBuf::from(output_devices),
            output_basename: output_basename.to_string(),
            extra_args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("cello-min");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cello.compile"

[bio.cello]
jar               = "/opt/cello/cello.jar"
verilog           = "circuit.v"
user_constraints  = "Eco1C1G1T1.UCF.json"
input_sensors     = "Eco1C1G1T1.input.json"
output_devices    = "Eco1C1G1T1.output.json"
output_basename   = "circuit_out"
"#,
        )
        .unwrap();
        let input = CelloInput::from_case_dir(&d).unwrap();
        assert_eq!(input.jar, PathBuf::from("/opt/cello/cello.jar"));
        assert_eq!(input.verilog, PathBuf::from("circuit.v"));
        assert_eq!(input.user_constraints, PathBuf::from("Eco1C1G1T1.UCF.json"));
        assert_eq!(input.input_sensors, PathBuf::from("Eco1C1G1T1.input.json"));
        assert_eq!(
            input.output_devices,
            PathBuf::from("Eco1C1G1T1.output.json")
        );
        assert_eq!(input.output_basename, "circuit_out");
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_verilog() {
        // Without a Verilog netlist Cello has no logic to compile.
        // Reject up front.
        let d = tempdir("cello-nover");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cello.compile"

[bio.cello]
jar               = "/opt/cello/cello.jar"
verilog           = ""
user_constraints  = "ucf.json"
input_sensors     = "in.json"
output_devices    = "out.json"
output_basename   = "circuit_out"
"#,
        )
        .unwrap();
        let err = CelloInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("verilog"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_user_constraints() {
        // The UCF pins the chassis and the available gate library;
        // empty UCF would crash Cello before any work is done.
        let d = tempdir("cello-noucf");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cello.compile"

[bio.cello]
jar               = "/opt/cello/cello.jar"
verilog           = "circuit.v"
user_constraints  = ""
input_sensors     = "in.json"
output_devices    = "out.json"
output_basename   = "circuit_out"
"#,
        )
        .unwrap();
        let err = CelloInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("user_constraints"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_basename() {
        // Output basename anchors collect()'s artefact filter; empty
        // string would surface every file in the workdir, including
        // the user's input JSONs and Verilog. Reject up front.
        let d = tempdir("cello-nobase");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "cello.compile"

[bio.cello]
jar               = "/opt/cello/cello.jar"
verilog           = "circuit.v"
user_constraints  = "ucf.json"
input_sensors     = "in.json"
output_devices    = "out.json"
output_basename   = ""
"#,
        )
        .unwrap();
        let err = CelloInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("output_basename"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
