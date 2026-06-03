//! 3MF reader and writer — **deferred to v1.5**.
//!
//! 3MF (3D Manufacturing Format) is a ZIP archive containing
//! `3D/3dmodel.model` XML plus a few metadata files. Implementing it
//! end-to-end requires a `zip` crate (or hand-rolled DEFLATE), which
//! is not currently in `Cargo.toml`'s workspace dependencies. Rather
//! than vendor a hand-rolled zip writer that no one has time to
//! audit, the v1 ship surfaces a clear error and lets users round-
//! trip via OBJ/PLY/STL instead.
//!
//! v1.5 plan: add `zip = "0.6"` to workspace deps, implement
//! `read_path` + `write_path` using `zip::ZipArchive` + a minimal
//! XML writer (`<resources>` + `<object id="1" type="model">` +
//! `<mesh>` with `<vertices>` and `<triangles>`).

use std::io;
use std::path::Path;

use thiserror::Error;

use crate::mesh::Mesh;

/// 3MF errors (currently always `Unsupported`).
#[derive(Debug, Error)]
pub enum ThreeMfError {
    /// IO problem (kept for future-compat).
    #[error(transparent)]
    Io(#[from] io::Error),
    /// 3MF support is not yet implemented.
    #[error(
        "3MF support is not yet implemented — v1 ships OBJ/PLY/STL. \
         3MF will land in v1.5 once a `zip` crate is added to \
         workspace dependencies. See docs/superpowers/plans/\
         2026-05-16-mesh-expansion-phase7.md (Task 25)."
    )]
    Unsupported,
}

/// Read a 3MF file. Always returns `ThreeMfError::Unsupported` in v1.
pub fn read_path(_path: impl AsRef<Path>) -> Result<Mesh, ThreeMfError> {
    Err(ThreeMfError::Unsupported)
}

/// Write a mesh as 3MF. Always returns `ThreeMfError::Unsupported` in v1.
pub fn write_path(_mesh: &Mesh, _path: impl AsRef<Path>) -> Result<(), ThreeMfError> {
    Err(ThreeMfError::Unsupported)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_returns_unsupported() {
        let err = read_path("nonexistent.3mf").unwrap_err();
        assert!(matches!(err, ThreeMfError::Unsupported));
    }

    #[test]
    fn write_returns_unsupported() {
        let m = Mesh::new("empty");
        let err = write_path(&m, "out.3mf").unwrap_err();
        assert!(matches!(err, ThreeMfError::Unsupported));
    }
}
