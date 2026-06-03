//! `[cad.occt]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "geometry"
//! solver  = "occt.run"
//!
//! [cad.occt]
//! script           = "model.py"
//! python           = "python3"          # optional, defaults to python3
//! # input_geometry = "input.step"       # optional .step / .stp / .iges / .igs / .brep
//! output_basename  = "model"
//! ```
//!
//! Modeled after the OpenMM adapter's `case_input.rs`: a flat
//! `toml::Value` parse so we don't pull serde-derive in for a tiny
//! struct.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct OcctInput {
    /// User-supplied `.py` script that imports `OCC.Core` to drive
    /// pythonocc-core. Resolved against the case directory and
    /// staged into the workdir before the run.
    pub script: PathBuf,
    /// Python interpreter — defaults to `python3`. Can be a bare
    /// name (resolved via PATH) or an absolute / relative path to a
    /// specific interpreter binary.
    pub python: String,
    /// Optional pre-existing geometry the script can load. If set,
    /// the file is staged alongside the script and surfaced to the
    /// script via `valenx_params.json`.
    pub input_geometry: Option<PathBuf>,
    /// Filename stem the user's script writes its outputs under
    /// (e.g. "model" → "model.step", "model.stl"). `collect()`
    /// walks the workdir for files that start with this stem.
    pub output_basename: String,
}

impl OcctInput {
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
            .get("cad")
            .and_then(|v| v.get("occt"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [cad.occt] section",
                    case_toml.display()
                ))
            })?;
        let script = block
            .get("script")
            .and_then(|v| v.as_str())
            .ok_or_else(|| AdapterError::Other(anyhow::anyhow!("[cad.occt].script required")))?;
        let python = block
            .get("python")
            .and_then(|v| v.as_str())
            .unwrap_or("python3")
            .to_string();
        let input_geometry = block
            .get("input_geometry")
            .and_then(|v| v.as_str())
            .map(PathBuf::from);
        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .unwrap_or("model")
            .to_string();
        Ok(Self {
            script: PathBuf::from(script),
            python,
            input_geometry,
            output_basename,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    #[test]
    fn parses_minimal() {
        let d = tempdir("occt");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "geometry"
solver  = "occt.run"

[cad.occt]
script          = "model.py"
output_basename = "model"
"#,
        )
        .unwrap();
        let input = OcctInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("model.py"));
        assert_eq!(input.python, "python3");
        assert_eq!(input.input_geometry, None);
        assert_eq!(input.output_basename, "model");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_input_geometry() {
        let d = tempdir("occt");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "geometry"
solver  = "occt.run"

[cad.occt]
script          = "boolean.py"
python          = "/opt/conda/envs/occt/bin/python"
input_geometry  = "bracket.step"
output_basename = "result"
"#,
        )
        .unwrap();
        let input = OcctInput::from_case_dir(&d).unwrap();
        assert_eq!(input.script, PathBuf::from("boolean.py"));
        assert_eq!(input.python, "/opt/conda/envs/occt/bin/python");
        assert_eq!(input.input_geometry, Some(PathBuf::from("bracket.step")));
        assert_eq!(input.output_basename, "result");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_non_py_extension() {
        // The case_input parser itself doesn't enforce the `.py`
        // extension — that check lives in `prepare()` so the user
        // gets a structured `InvalidCase` error pointing at
        // case.toml. Here we exercise that the parser still happily
        // accepts the path; the contract is that prepare() rejects.
        let d = tempdir("occt");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "geometry"
solver  = "occt.run"

[cad.occt]
script          = "model.sh"
output_basename = "model"
"#,
        )
        .unwrap();
        let input = OcctInput::from_case_dir(&d).expect("parser is permissive");
        assert_eq!(input.script, PathBuf::from("model.sh"));
        // The extension is preserved verbatim — verifying this lets
        // us assert in lib.rs that prepare() rejects on the
        // extension check rather than on a parse error.
        assert_eq!(
            input
                .script
                .extension()
                .and_then(|e| e.to_str())
                .map(str::to_ascii_lowercase),
            Some("sh".into())
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
