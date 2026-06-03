//! `[bio.nupack]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "nupack.analyze"
//!
//! [bio.nupack]
//! script          = "design.py"
//! python          = "python3"          # optional, defaults to "python3"
//! input           = "complex.fa"       # optional — staged next to script
//! output_basename = "design"
//! temperature     = 37.0               # optional, defaults to 37.0 (Celsius)
//! sodium          = 1.0                # optional, defaults to 1.0 M
//! ```
//!
//! NUPACK's modern Python API has no canonical CLI — the user
//! authors a script that calls into `nupack.tubes`, `nupack.complex`,
//! `nupack.thermodynamics`, etc. The adapter stages the script
//! (and optional input file) into the workdir, drops a flat
//! `valenx_params.json` with the parsed knobs (input filename,
//! output_basename, temperature, sodium), and invokes
//! `python <script>`. Scripts read `valenx_params.json` and pass the
//! values through to NUPACK themselves — same convention as
//! RFdiffusion / DeepChem / ESMFold.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct NupackInput {
    pub script: PathBuf,
    pub python: String,
    pub input: Option<PathBuf>,
    pub output_basename: String,
    pub temperature: f64,
    pub sodium: f64,
}

impl NupackInput {
    pub fn from_case_dir(case_dir: &std::path::Path) -> Result<Self, AdapterError> {
        let case_toml = case_dir.join("case.toml");
        let text = valenx_core::io_caps::read_capped_to_string(
            &case_toml,
            valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
        )
        .map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("read {}: {e}", case_toml.display()))
        })?;
        let parsed: toml::Value = toml::from_str(&text).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("parse {}: {e}", case_toml.display()))
        })?;
        let block = parsed
            .get("bio")
            .and_then(|v| v.get("nupack"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.nupack] section",
                    case_toml.display()
                ))
            })?;

        let script_str = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[bio.nupack].script required")))?;
        if script_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.nupack].script must not be empty"
            )));
        }

        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();

        let input = block
            .get("input")
            .and_then(|v| v.as_str())
            .filter(|s| !s.is_empty())
            .map(PathBuf::from);

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.nupack].output_basename required \
                     (filename stem for NUPACK outputs)"
                ))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.nupack].output_basename must not be empty"
            )));
        }

        let temperature = match block.get("temperature") {
            Some(v) => {
                let raw = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.nupack].temperature must be a number"
                        ))
                    })?;
                if !raw.is_finite() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.nupack].temperature must be finite, got {raw}"
                    )));
                }
                raw
            }
            None => 37.0,
        };

        let sodium = match block.get("sodium") {
            Some(v) => {
                let raw = v
                    .as_float()
                    .or_else(|| v.as_integer().map(|i| i as f64))
                    .ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!("[bio.nupack].sodium must be a number"))
                    })?;
                if !raw.is_finite() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.nupack].sodium must be finite, got {raw}"
                    )));
                }
                if raw <= 0.0 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.nupack].sodium must be > 0.0 (Molar), got {raw}"
                    )));
                }
                raw
            }
            None => 1.0,
        };

        Ok(Self {
            script: PathBuf::from(script_str),
            python,
            input,
            output_basename: output_basename.to_string(),
            temperature,
            sodium,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        // Minimum-config case: just script + output_basename.
        // Defaults must match documented values.
        let d = tempdir("nupack");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "nupack.analyze"

[bio.nupack]
script          = "design.py"
output_basename = "design"
"#,
        )
        .unwrap();
        let input = NupackInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("design.py"));
        assert_eq!(input.python, "python3");
        assert_eq!(input.input, None);
        assert_eq!(input.output_basename, "design");
        assert_eq!(input.temperature, 37.0);
        assert_eq!(input.sodium, 1.0);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_input_and_overrides() {
        // Optional input file plus all numeric knobs explicitly set.
        let d = tempdir("nupack");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "nupack.analyze"

[bio.nupack]
script          = "design.py"
python          = "/opt/conda/envs/nupack/bin/python"
input           = "complex.fa"
output_basename = "tube"
temperature     = 25.0
sodium          = 0.15
"#,
        )
        .unwrap();
        let input = NupackInput::from_case_dir(&d).unwrap();
        assert_eq!(input.python, "/opt/conda/envs/nupack/bin/python");
        assert_eq!(input.input, Some(PathBuf::from("complex.fa")));
        assert_eq!(input.output_basename, "tube");
        assert_eq!(input.temperature, 25.0);
        assert!((input.sodium - 0.15).abs() < 1e-9);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_basename() {
        // Empty basename would mean every NUPACK output collides at
        // the workdir root — reject up front.
        let d = tempdir("nupack");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "nupack.analyze"

[bio.nupack]
script          = "design.py"
output_basename = ""
"#,
        )
        .unwrap();
        let err = NupackInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("output_basename"), "msg: {msg}");
        assert!(msg.contains("must not be empty"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_sodium() {
        // Sodium concentration must be strictly positive — 0 M isn't
        // physical for NUPACK's salt-correction model.
        let d = tempdir("nupack");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "nupack.analyze"

[bio.nupack]
script          = "design.py"
output_basename = "design"
sodium          = 0.0
"#,
        )
        .unwrap();
        let err = NupackInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("sodium"), "msg: {msg}");
        assert!(msg.contains("> 0.0"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
