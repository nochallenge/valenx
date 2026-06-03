//! `[cfd.su2]` case-input parsing for the SU2 adapter.

use std::path::PathBuf;

use serde::Deserialize;

use valenx_core::AdapterError;

/// Parsed `[cfd.su2]` block.
#[derive(Clone, Debug, PartialEq)]
pub struct Su2Input {
    /// SU2 .cfg file (relative to case dir). The cfg references the
    /// mesh internally via `MESH_FILENAME=`.
    pub config: PathBuf,
    /// Optional mesh file. When set, staged into the workdir
    /// alongside the cfg so the cfg's relative `MESH_FILENAME=`
    /// resolves. Common extensions: .su2, .cgns, .msh.
    pub mesh: Option<PathBuf>,
    /// Optional OpenMP threads — passed via OMP_NUM_THREADS env var.
    /// `None` lets SU2 / OpenMP pick the system default.
    pub n_threads: Option<u32>,
}

#[derive(Deserialize)]
struct CaseToml {
    case: Option<CaseHeader>,
    cfd: Option<CfdTable>,
}

#[derive(Deserialize)]
struct CaseHeader {
    #[serde(default)]
    physics: String,
}

#[derive(Deserialize)]
struct CfdTable {
    su2: Option<Su2Toml>,
}

#[derive(Deserialize)]
struct Su2Toml {
    config: String,
    #[serde(default)]
    mesh: Option<String>,
    #[serde(default)]
    n_threads: Option<u32>,
}

impl Su2Input {
    pub fn from_case_dir(case_dir: &std::path::Path) -> Result<Self, AdapterError> {
        let toml_path = case_dir.join("case.toml");
        // Round-18 H1 (R17 sweep gap): cap the case.toml read at the
        // shared `MAX_PROJECT_FILE_BYTES` so a hostile multi-GB file
        // can't be slurped into memory before serde sees it.
        let text = valenx_core::io_caps::read_capped_to_string(
            &toml_path,
            valenx_core::project::loader::MAX_PROJECT_FILE_BYTES as usize,
        )
        .map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("read {}: {e}", toml_path.display()))
        })?;
        let parsed: CaseToml = toml::from_str(&text).map_err(|e| {
            AdapterError::Other(anyhow::anyhow!("parse {}: {e}", toml_path.display()))
        })?;
        if let Some(ref hdr) = parsed.case {
            if !hdr.physics.is_empty()
                && !matches!(hdr.physics.as_str(), "cfd" | "compressible" | "aero")
            {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "case physics is `{}` — SU2 handles cfd / compressible / aero",
                    hdr.physics
                )));
            }
        }
        let block = parsed.cfd.and_then(|c| c.su2).ok_or_else(|| {
            AdapterError::Other(anyhow::anyhow!(
                "{} has no [cfd.su2] section — add `config = \"...\"`",
                toml_path.display()
            ))
        })?;
        if let Some(n) = block.n_threads {
            if n == 0 {
                return Err(AdapterError::Other(anyhow::anyhow!(
                    "n_threads must be > 0 (use omit-or-1 for serial)"
                )));
            }
        }
        Ok(Su2Input {
            config: PathBuf::from(block.config),
            mesh: block.mesh.map(PathBuf::from),
            n_threads: block.n_threads,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_test_utils::tempdir;

    fn write_case_toml(dir: &std::path::Path, content: &str) {
        std::fs::write(dir.join("case.toml"), content).unwrap();
    }
    #[test]
    fn parses_minimal_config() {
        let d = tempdir("su2-min");
        write_case_toml(
            &d,
            r#"
[case]
physics = "cfd"

[cfd.su2]
config = "wing.cfg"
"#,
        );
        let input = Su2Input::from_case_dir(&d).expect("parse");
        assert_eq!(input.config, PathBuf::from("wing.cfg"));
        assert!(input.mesh.is_none());
        assert!(input.n_threads.is_none());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn picks_up_mesh_and_threads() {
        let d = tempdir("su2-full");
        write_case_toml(
            &d,
            r#"
[case]
physics = "compressible"

[cfd.su2]
config = "naca0012.cfg"
mesh = "naca0012.su2"
n_threads = 4
"#,
        );
        let input = Su2Input::from_case_dir(&d).expect("parse");
        assert_eq!(input.config, PathBuf::from("naca0012.cfg"));
        assert_eq!(input.mesh, Some(PathBuf::from("naca0012.su2")));
        assert_eq!(input.n_threads, Some(4));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn missing_section_actionable() {
        let d = tempdir("su2-missing");
        write_case_toml(&d, "[case]\nphysics = \"cfd\"\n");
        let r = Su2Input::from_case_dir(&d);
        assert!(format!("{}", r.unwrap_err()).contains("[cfd.su2]"));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_zero_threads() {
        let d = tempdir("su2-zero");
        write_case_toml(&d, "[cfd.su2]\nconfig = \"x.cfg\"\nn_threads = 0\n");
        let r = Su2Input::from_case_dir(&d);
        assert!(r.is_err());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn rejects_wrong_physics() {
        let d = tempdir("su2-wrong");
        write_case_toml(
            &d,
            "[case]\nphysics = \"meshing\"\n[cfd.su2]\nconfig = \"x.cfg\"\n",
        );
        let r = Su2Input::from_case_dir(&d);
        assert!(format!("{}", r.unwrap_err()).contains("meshing"));
        let _ = std::fs::remove_dir_all(&d);
    }

    /// Round-18 H1 RED→GREEN: a case.toml above
    /// `MAX_PROJECT_FILE_BYTES` is rejected before serde sees the
    /// bytes. Sparse-file trick so we don't write 5 MiB of zeros to
    /// disk — `Seek + write_all` makes `metadata.len()` report
    /// past the cap while consuming ~0 disk blocks.
    #[test]
    fn rejects_oversize_case_toml() {
        use std::io::{Seek, SeekFrom, Write};
        let d = tempdir("su2-oversize");
        let toml_path = d.join("case.toml");
        let cap = valenx_core::project::loader::MAX_PROJECT_FILE_BYTES;
        {
            let mut f = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&toml_path)
                .unwrap();
            // 5 MiB worth of "size" via sparse-file seek; cap is 1 MiB.
            f.seek(SeekFrom::Start(5 * 1024 * 1024)).unwrap();
            f.write_all(b"x").unwrap();
        }
        let err = Su2Input::from_case_dir(&d).expect_err("must reject oversize");
        let msg = format!("{err}");
        assert!(
            msg.contains("cap") || msg.contains("exceeds"),
            "expected size-cap error, got: {msg}"
        );
        assert!(
            cap < 5 * 1024 * 1024,
            "test assumption: cap ({cap}) < 5 MiB"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}
