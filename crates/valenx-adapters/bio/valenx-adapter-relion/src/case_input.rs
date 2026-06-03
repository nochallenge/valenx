//! `[bio.relion]` case-input parsing. Schema:
//!
//! ```toml
//! [case]
//! physics = "bio"
//! solver  = "relion.refine"
//!
//! [bio.relion]
//! particles       = "particles.star"
//! reference       = "ref.mrc"
//! output_basename = "Refine3D/run"
//! angpix          = 1.06         # pixel size in Angstrom
//! mpi_procs       = 4            # optional, defaults to 1
//! threads         = 8            # optional, defaults to 1
//! extra_args      = ["--auto_refine"]   # optional, defaults to []
//! ```
//!
//! RELION's `relion_refine` is the workhorse 3D refinement entry-point.
//! In single-process mode the binary is invoked directly; with
//! `mpi_procs > 1` the adapter switches to `mpirun -n <N>
//! relion_refine_mpi`.
//!
//! `output_basename` is the RELION-style prefix that the binary
//! prepends to every artifact name (`<basename>_class001.mrc`,
//! `<basename>_data.star`, `<basename>_model.star`, …). Pinning it
//! here lets `collect()` walk the workdir for files matching the
//! prefix without having to scan and guess.

use std::path::PathBuf;
use valenx_core::AdapterError;

#[derive(Clone, Debug, PartialEq)]
pub struct RelionInput {
    pub particles: PathBuf,
    pub reference: PathBuf,
    pub output_basename: String,
    pub angpix: f64,
    pub mpi_procs: u32,
    pub threads: u32,
    pub extra_args: Vec<String>,
}

impl RelionInput {
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
            .and_then(|v| v.get("relion"))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "{} missing [bio.relion] section",
                    case_toml.display()
                ))
            })?;

        let particles_str = block
            .get("particles")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.relion].particles required (path to particle stack STAR file)"
                ))
            })?;
        if particles_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.relion].particles must not be empty"
            )));
        }

        let reference_str = block
            .get("reference")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.relion].reference required (path to initial reference MRC)"
                ))
            })?;
        if reference_str.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.relion].reference must not be empty"
            )));
        }

        let output_basename = block
            .get("output_basename")
            .and_then(|v| v.as_str())
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.relion].output_basename required (run prefix RELION prepends to outputs)"
                ))
            })?;
        if output_basename.is_empty() {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.relion].output_basename must not be empty"
            )));
        }

        let angpix = block
            .get("angpix")
            .and_then(|v| v.as_float().or_else(|| v.as_integer().map(|i| i as f64)))
            .ok_or_else(|| {
                AdapterError::Other(anyhow::anyhow!(
                    "[bio.relion].angpix required (pixel size in Angstrom, > 0)"
                ))
            })?;
        if !angpix.is_finite() || angpix <= 0.0 {
            return Err(AdapterError::Other(anyhow::anyhow!(
                "[bio.relion].angpix must be a finite positive number, got {angpix}"
            )));
        }

        let mpi_procs = match block.get("mpi_procs") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!(
                        "[bio.relion].mpi_procs must be an integer"
                    ))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.relion].mpi_procs must be >= 1, got {raw}"
                    )));
                }
                raw as u32
            }
            None => 1,
        };

        let threads = match block.get("threads") {
            Some(v) => {
                let raw = v.as_integer().ok_or_else(|| {
                    AdapterError::Other(anyhow::anyhow!("[bio.relion].threads must be an integer"))
                })?;
                if raw < 1 {
                    return Err(AdapterError::Other(anyhow::anyhow!(
                        "[bio.relion].threads must be >= 1, got {raw}"
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
                        "[bio.relion].extra_args must be an array of strings"
                    ))
                })?;
                let mut out = Vec::with_capacity(arr.len());
                for entry in arr {
                    let s = entry.as_str().ok_or_else(|| {
                        AdapterError::Other(anyhow::anyhow!(
                            "[bio.relion].extra_args entries must be strings"
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
            reference: PathBuf::from(reference_str),
            output_basename: output_basename.to_string(),
            angpix,
            mpi_procs,
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
        // Minimum-viable RELION run: particles + reference + output
        // prefix + pixel size. Defaults: 1 MPI proc, 1 thread, no
        // extras.
        let d = tempdir("relion");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "relion.refine"

[bio.relion]
particles       = "particles.star"
reference       = "ref.mrc"
output_basename = "Refine3D/run"
angpix          = 1.06
"#,
        )
        .unwrap();
        let input = RelionInput::from_case_dir(&d).unwrap();
        assert_eq!(input.particles, PathBuf::from("particles.star"));
        assert_eq!(input.reference, PathBuf::from("ref.mrc"));
        assert_eq!(input.output_basename, "Refine3D/run");
        assert!((input.angpix - 1.06).abs() < 1e-9);
        assert_eq!(input.mpi_procs, 1);
        assert_eq!(input.threads, 1);
        assert!(input.extra_args.is_empty());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn parses_with_mpi_and_threads() {
        // Full HPC config: 8 MPI ranks, 4 threads each (= 32-way
        // parallelism), with the canonical `--auto_refine` 3D
        // auto-refine extra.
        let d = tempdir("relion");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "relion.refine"

[bio.relion]
particles       = "particles.star"
reference       = "ref.mrc"
output_basename = "Refine3D/run"
angpix          = 1.06
mpi_procs       = 8
threads         = 4
extra_args      = ["--auto_refine"]
"#,
        )
        .unwrap();
        let input = RelionInput::from_case_dir(&d).unwrap();
        assert_eq!(input.mpi_procs, 8);
        assert_eq!(input.threads, 4);
        assert_eq!(input.extra_args, vec!["--auto_refine".to_string()]);
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_angpix() {
        // 0 pixel size is meaningless and would produce nonsense
        // metric reconstructions; reject up front.
        let d = tempdir("relion");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "relion.refine"

[bio.relion]
particles       = "particles.star"
reference       = "ref.mrc"
output_basename = "Refine3D/run"
angpix          = 0.0
"#,
        )
        .unwrap();
        let err = RelionInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("angpix"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_mpi() {
        // 0 MPI procs is undefined; we floor at 1 (single-process
        // mode invokes `relion_refine` directly) and reject below.
        let d = tempdir("relion");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "relion.refine"

[bio.relion]
particles       = "particles.star"
reference       = "ref.mrc"
output_basename = "Refine3D/run"
angpix          = 1.06
mpi_procs       = 0
"#,
        )
        .unwrap();
        let err = RelionInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("mpi_procs"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_threads() {
        // 0 threads collapses to ambiguous behaviour in relion's
        // OpenMP layer; reject up front.
        let d = tempdir("relion");
        std::fs::write(
            d.join("case.toml"),
            r#"[case]
physics = "bio"
solver  = "relion.refine"

[bio.relion]
particles       = "particles.star"
reference       = "ref.mrc"
output_basename = "Refine3D/run"
angpix          = 1.06
threads         = 0
"#,
        )
        .unwrap();
        let err = RelionInput::from_case_dir(&d).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("threads"), "msg: {msg}");
        let _ = std::fs::remove_dir_all(&d);
    }
}
