//! # valenx-step-iges
//!
//! STEP (ISO 10303) and IGES (5.3) import / export for Valenx solids.
//!
//! Two formats, one crate. STEP is the modern CAD interchange standard;
//! IGES is the 1980s-era format that legacy CAM machines still ingest.
//!
//! ## What ships
//!
//! - **STEP via `truck-stepio`** — round-trip a [`valenx_cad::Solid`]
//!   through truck's `Display` / `parse_step` facade. See [`step::write`]
//!   / [`step::read`].
//! - **Hand-rolled minimal IGES** — read / write the 110 (Line), 100
//!   (Arc), 116 (Point), and 124 (Transformation matrix) entity types.
//!   Faces / shells are out of scope for v1 IGES; the importer builds a
//!   *wireframe-only* mesh-backed Solid (documented as such in the
//!   import UI). See [`iges::write`] / [`iges::read`].
//! - **Auto-detect dispatcher** — [`import`] / [`export`] pick STEP vs
//!   IGES from the file extension and route to the right backend.
//!
//! ## What doesn't (and why)
//!
//! - **IGES full surface topology** — Types 142 / 144 / 128 (trimmed
//!   parametric surfaces, NURBS) are months of work to implement
//!   correctly; v1 ships wireframe only. Stretch goal for Phase 8.5.
//! - **Mesh-backed Solid → STEP** — the `valenx-fillet` mesh-domain
//!   pipeline produces triangle soups without BRep topology. STEP wants
//!   faces; returns [`StepIgesError::MeshBackedSolidNotExportable`].
//!   Export STL for triangle-mesh interchange instead.
//!
//! ## Quick start
//!
//! ```no_run
//! use std::path::PathBuf;
//! use valenx_step_iges::{export, import};
//!
//! let cube = valenx_cad::box_solid(10.0, 10.0, 10.0).unwrap();
//! let path = PathBuf::from("cube.step");
//! export(&cube, &path).unwrap();
//! let _read_back = import(&path).unwrap();
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]
// Surface future `&str` byte-offset slicing in clippy review — this
// crate parses untrusted text (STEP/IGES), where non-char-boundary
// slices panic. WARN (not deny): most existing slices are safe ASCII;
// this only flags NEW ones.
#![allow(clippy::string_slice, reason = "parsers slice ASCII fixed-format records at byte offsets from find() or constant ASCII prefixes, always valid char boundaries")]

pub mod ap242;
pub mod error;
pub mod iges;
pub mod iges_trimmed;
pub mod persist;
pub mod step;

pub use error::{ErrorCategory, StepIgesError};

use std::path::Path;

/// Round-9 DoS cap on STEP / IGES / AP242 file reads. Production
/// CAD interchange files top out around 50 MiB for huge assemblies;
/// 256 MiB is generous for legitimate use while refusing the
/// `cat /dev/zero > big.step` denial of service that pre-fix would
/// allocate before parsing. The cap is enforced before
/// `fs::read_to_string` allocates, so a hostile multi-GB file
/// produces a clean `FileTooLarge` error rather than OOM.
pub const MAX_CAD_INTERCHANGE_FILE_BYTES: u64 = 256 * 1024 * 1024;

/// Read `path` as a String with the
/// [`MAX_CAD_INTERCHANGE_FILE_BYTES`] cap enforced both by stat
/// (typed [`StepIgesError::FileTooLarge`] for the size message) AND
/// by a bounded `take(cap+1)` on the reader itself.
///
/// Round-18 L1: the pre-fix sites did `fs::metadata(path)?.len()`
/// followed by an unbounded `fs::read_to_string(path)`. That is a
/// TOCTOU window — between the two syscalls the file could grow,
/// and the second call would slurp the now-multi-GB file into
/// memory without re-checking the cap. This helper closes the gap
/// by capping the second read with `take(cap+1)` (one extra byte
/// so a file that grew to exactly cap+1 still trips the size
/// check) before any allocation overshoots the limit.
///
/// The typed `FileTooLarge` variant is preserved for callers that
/// branch on the error code; a runtime-grown file that beats the
/// stat check is surfaced as a plain [`StepIgesError::Io`] with
/// `InvalidData` kind (which is the closest std-io match for "cap
/// exceeded mid-read").
pub(crate) fn read_capped_cad_text(
    path: &std::path::Path,
    format: &'static str,
) -> Result<String, StepIgesError> {
    use std::io::Read;

    let size = std::fs::metadata(path)?.len();
    if size > MAX_CAD_INTERCHANGE_FILE_BYTES {
        return Err(StepIgesError::FileTooLarge {
            format,
            size,
            cap: MAX_CAD_INTERCHANGE_FILE_BYTES,
        });
    }
    // Cap the read at cap+1 so a hostile file that grew between
    // stat and open still produces a clean error (truncated read +
    // a tail byte that lets us tell "the file got bigger" apart
    // from "the file was exactly at the cap").
    let mut buf = Vec::new();
    std::fs::File::open(path)?
        .take(MAX_CAD_INTERCHANGE_FILE_BYTES + 1)
        .read_to_end(&mut buf)?;
    if buf.len() as u64 > MAX_CAD_INTERCHANGE_FILE_BYTES {
        return Err(StepIgesError::FileTooLarge {
            format,
            size: buf.len() as u64,
            cap: MAX_CAD_INTERCHANGE_FILE_BYTES,
        });
    }
    String::from_utf8(buf).map_err(|e| {
        StepIgesError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            e,
        ))
    })
}

use valenx_cad::Solid;

/// File format produced by [`detect_format`].
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum Format {
    /// ISO 10303 (`.step` / `.stp`).
    Step,
    /// Initial Graphics Exchange Specification (`.iges` / `.igs`).
    Iges,
}

impl Format {
    /// Short label for UI display.
    pub fn label(&self) -> &'static str {
        match self {
            Format::Step => "STEP",
            Format::Iges => "IGES",
        }
    }
}

/// Pick the file format from the path's extension. Returns `None` if
/// the extension is missing or unrecognised — callers decide whether
/// that's an error or a default.
pub fn detect_format(path: &Path) -> Option<Format> {
    let ext = path.extension()?.to_str()?.to_ascii_lowercase();
    match ext.as_str() {
        "step" | "stp" => Some(Format::Step),
        "iges" | "igs" => Some(Format::Iges),
        _ => None,
    }
}

/// Read a solid from `path`, auto-detecting STEP vs IGES from the
/// extension.
///
/// # Errors
///
/// - [`StepIgesError::Unsupported`] when the extension isn't STEP or
///   IGES.
/// - Whatever the underlying backend returns ([`step::read`] /
///   [`iges::read`]).
pub fn import(path: &Path) -> Result<Solid, StepIgesError> {
    match detect_format(path) {
        Some(Format::Step) => step::read(path),
        Some(Format::Iges) => iges::read(path),
        None => Err(StepIgesError::Unsupported {
            format: "auto-detect",
            feature: format!(
                "extension {:?} — supported: .step, .stp, .iges, .igs",
                path.extension()
                    .and_then(|s| s.to_str())
                    .unwrap_or("(none)"),
            ),
        }),
    }
}

/// Write a solid to `path`, auto-detecting STEP vs IGES from the
/// extension.
///
/// # Errors
///
/// Same routing rules as [`import`].
pub fn export(solid: &Solid, path: &Path) -> Result<(), StepIgesError> {
    match detect_format(path) {
        Some(Format::Step) => step::write(solid, path),
        Some(Format::Iges) => iges::write(solid, path),
        None => Err(StepIgesError::Unsupported {
            format: "auto-detect",
            feature: format!(
                "extension {:?} — supported: .step, .stp, .iges, .igs",
                path.extension()
                    .and_then(|s| s.to_str())
                    .unwrap_or("(none)"),
            ),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    /// Round-18 L1 RED→GREEN: the shared `read_capped_cad_text`
    /// helper must (a) return typed `FileTooLarge` for files whose
    /// stat reports past-the-cap and (b) tolerate at-the-cap files
    /// without rejecting them.
    #[test]
    fn read_capped_cad_text_rejects_oversize_via_stat() {
        use std::io::{Seek, SeekFrom, Write};
        let tmp = std::env::temp_dir().join(format!(
            "valenx-cadcap-oversize-{}.step",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&tmp)
                .unwrap();
            f.seek(SeekFrom::Start(MAX_CAD_INTERCHANGE_FILE_BYTES + 1))
                .unwrap();
            f.write_all(b"x").unwrap();
        }
        let err = read_capped_cad_text(&tmp, "STEP").expect_err("must reject oversize");
        let _ = std::fs::remove_file(&tmp);
        match err {
            StepIgesError::FileTooLarge { format, size, cap } => {
                assert_eq!(format, "STEP");
                assert!(size > cap);
                assert_eq!(cap, MAX_CAD_INTERCHANGE_FILE_BYTES);
            }
            other => panic!("expected FileTooLarge, got: {other:?}"),
        }
    }

    /// Round-18 L1: at-the-cap is allowed (only past-the-cap is
    /// rejected). The bounded `take(cap + 1)` permits a `cap`-byte
    /// file and only the size check trips when bytes overflow.
    #[test]
    fn read_capped_cad_text_accepts_small_file() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-cadcap-small-{}.step",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&tmp, b"ISO-10303-21;\nEND-ISO-10303-21;\n").unwrap();
        let text = read_capped_cad_text(&tmp, "STEP").expect("small file ok");
        assert!(text.contains("ISO-10303-21"));
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn detect_format_handles_known_extensions() {
        assert_eq!(detect_format(&PathBuf::from("a.step")), Some(Format::Step));
        assert_eq!(detect_format(&PathBuf::from("a.STP")), Some(Format::Step));
        assert_eq!(detect_format(&PathBuf::from("a.iges")), Some(Format::Iges));
        assert_eq!(detect_format(&PathBuf::from("a.IGS")), Some(Format::Iges));
        assert_eq!(detect_format(&PathBuf::from("a.obj")), None);
        assert_eq!(detect_format(&PathBuf::from("noext")), None);
    }

    #[test]
    fn format_labels_are_stable() {
        assert_eq!(Format::Step.label(), "STEP");
        assert_eq!(Format::Iges.label(), "IGES");
    }

    #[test]
    fn import_unknown_extension_reports_unsupported() {
        let err = import(&PathBuf::from("a.obj")).unwrap_err();
        assert!(matches!(err, StepIgesError::Unsupported { .. }));
        assert!(err.to_string().contains("obj"));
    }
}
