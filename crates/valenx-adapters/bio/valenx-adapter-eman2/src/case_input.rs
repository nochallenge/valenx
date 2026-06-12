//! `[bio.eman2]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "eman2.refine"
//!
//! [bio.eman2]
//! particles         = "particles.lst"
//! model             = "initial_model.hdf"
//! output_basename   = "refine_01"
//! target_resolution = 8.0       # Angstrom
//! symmetry          = "c1"      # optional, defaults to "c1"
//! threads           = 8         # optional, defaults to 1
//! extra_args        = ["--speed=5"]   # optional, defaults to []
//! ```
//!
//! `e2refine_easy.py` is EMAN2's high-level "easy" 3D refinement
//! driver. It expects a particle list (`.lst` / `.hdf` stack), an
//! initial 3D model (`.hdf`), an output basename that EMAN2 turns into
//! a `<basename>_NN/` results directory, a target resolution to
//! refine toward, and a point-group symmetry.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct Eman2Input {
    pub particles: PathBuf,
    pub model: PathBuf,
    pub output_basename: String,
    pub target_resolution: f64,
    pub symmetry: String,
    pub threads: u32,
    pub extra_args: Vec<String>,
}

impl Eman2Input {
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
            .and_then(|v| v.get("eman2"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.eman2] section",
                    case_toml.display()
                ))
            })?;

        let particles_str = block
            .get("particles")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.eman2].particles required (path to particle list / stack)"
                ))
            })?;
        if particles_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.eman2].particles must not be empty"
            )));
        }

        let model_str = block.get("model").and_then(|v| v.as_str()).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "[bio.eman2].model required (path to initial 3D model HDF)"
            ))
        })?;
        if model_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.eman2].model must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.eman2].output_basename required (run prefix EMAN2 turns into <basename>_NN)"
                ))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.eman2].output_basename must not be empty"
            )));
        }

        let target_resolution = block
            .get("target_resolution")
            .and_then(|v| v.as_float().or_else(|| v.as_integer().map(|i| i as f64)))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.eman2].target_resolution required (Angstrom, > 0)"
                ))
            })?;
        if !target_resolution.is_finite() || target_resolution <= 0.0 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.eman2].target_resolution must be a finite positive number, got {target_resolution}"
            )));
        }

        let symmetry = match block.get("symmetry") {
            Some(v) => {
                let s = v.as_str().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.eman2].symmetry must be a string (e.g. \"c1\", \"d7\", \"icos\")"
                    ))
                })?;
                if s.trim().is_empty() {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.eman2].symmetry must not be empty (use \"c1\" for no symmetry)"
                    )));
                }
                s.to_string()
            }
            None => "c1".to_string(),
        };

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.eman2].threads must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.eman2].threads must be >= 1, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 1,
        };

        let extra_args = match block.get("extra_args") {
            Some(arr) => {
                let arr = arr.as_array().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.eman2].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.eman2].extra_args entries must be strings"
                        ))
                    })?;
                    out.push(s.to_string());
                }
                out
            }
            None => Vec::new(),
        };

        Ok(Self {
            particles: PathBuf::from(particles_str),
            model: PathBuf::from(model_str),
            output_basename: output_basename.to_string(),
            target_resolution,
            symmetry,
            threads,
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
        // Bare EMAN2 refinement: particles + model + output prefix +
        // target resolution. Defaults: c1 symmetry, 1 thread, no
        // extras.
        let d = tempdir("eman2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "eman2.refine"

[bio.eman2]
particles         = "particles.lst"
model             = "initial_model.hdf"
output_basename   = "refine_01"
target_resolution = 8.0
"#,
        )
        .unwrap();
        let input = Eman2Input::from_case_dir(&d).unwrap();
        assert_eq!(input.particles, PathBuf::from("particles.lst"));
        assert_eq!(input.model, PathBuf::from("initial_model.hdf"));
        assert_eq!(input.output_basename, "refine_01");
        assert!((input.target_resolution - 8.0).abs() < 1e-9);
        assert_eq!(input.symmetry, "c1");
        assert_eq!(input.threads, 1);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_symmetry_and_threads() {
        // Icosahedral particle (typical of viral capsids) with high
        // thread count and a speed-tuning extra.
        let d = tempdir("eman2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "eman2.refine"

[bio.eman2]
particles         = "particles.lst"
model             = "initial_model.hdf"
output_basename   = "refine_01"
target_resolution = 4.5
symmetry          = "icos"
threads           = 16
extra_args        = ["--speed=5"]
"#,
        )
        .unwrap();
        let input = Eman2Input::from_case_dir(&d).unwrap();
        assert_eq!(input.symmetry, "icos");
        assert_eq!(input.threads, 16);
        assert_eq!(input.extra_args, vec!["--speed=5".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_empty_symmetry() {
        // An empty symmetry string would be silently accepted by EMAN2
        // and produce undefined behaviour; reject up front. Whitespace
        // also counts as empty after trim.
        let d = tempdir("eman2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "eman2.refine"

[bio.eman2]
particles         = "particles.lst"
model             = "initial_model.hdf"
output_basename   = "refine_01"
target_resolution = 8.0
symmetry          = "   "
"#,
        )
        .unwrap();
        let err = Eman2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("symmetry"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_resolution() {
        // 0 Angstrom target resolution is meaningless; reject up
        // front.
        let d = tempdir("eman2");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "eman2.refine"

[bio.eman2]
particles         = "particles.lst"
model             = "initial_model.hdf"
output_basename   = "refine_01"
target_resolution = 0.0
"#,
        )
        .unwrap();
        let err = Eman2Input::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("target_resolution"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
