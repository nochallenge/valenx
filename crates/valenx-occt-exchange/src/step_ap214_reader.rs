//! Phase 104 — STEP AP214 (ISO 10303-214) reader.
//!
//! ## What OCCT does
//!
//! `STEPControl_Reader` with the AP214-aware `STEPCAFControl_Reader`
//! sibling resolves the geometry path the same way AP203 does, then
//! also pulls in the AP214-specific authority block + per-face
//! colour into an `XCAFDoc_DocumentTool`. The returned shape carries
//! the colour attribution alongside the BRep, which downstream tools
//! consume via `XCAFDoc_ColorTool::GetColor`.
//!
//! ## v1 status
//!
//! **Honest implementation** (Phase 104.5). The geometry is imported
//! by the proven [`fn@crate::step_ap203_reader`] backend
//! (truck-stepio reads AP214 geometry transparently — the schemas
//! differ only in PDM metadata, not geometry). The AP214-specific
//! attributes are then recovered with a text-level scan:
//!
//! - per-solid colour from `COLOUR_RGB` / `COLOUR_RGB` entities, via
//!   the existing [`valenx_step_iges::ap242::parse_metadata`];
//! - PDM authority strings from `PERSON_AND_ORGANIZATION`,
//!   `APPROVAL_PERSON_ORGANIZATION`, `APPROVAL_STATUS`, and
//!   `SECURITY_CLASSIFICATION_LEVEL` entities.
//!
//! The scan is keyword-driven (the same approach `parse_metadata`
//! uses) — robust to the AP214 dialect variations different CAD
//! packages emit, and it round-trips the colours written by
//! [`fn@crate::step_ap214_writer`].

use std::path::Path;

use valenx_cad::Solid;
use valenx_step_iges::ap242::Ap242Color;

use crate::error::OcctExchangeError;

/// Bundle of geometry + per-face colour + PDM authority recovered
/// from an AP214 file.
#[derive(Clone, Debug)]
pub struct Ap214Import {
    /// The imported geometry, as for [`crate::step_ap203_reader()`].
    pub solid: Solid,
    /// Per-face / per-solid colour attribution recovered from the
    /// AP214 colour entities. Empty when the file carried no colour.
    pub colors: Vec<Ap242Color>,
    /// Free-form PDM authority strings recovered from
    /// `PERSON_AND_ORGANIZATION` / `APPROVAL_*`. Empty when none.
    pub authority: Vec<String>,
}

/// Read a STEP AP214 file from `path`, returning geometry plus the
/// AP214 colour + authority attributes.
///
/// # Errors
///
/// - [`OcctExchangeError::BadInput`] if `path` is not `.step`/`.stp`.
/// - [`OcctExchangeError::Parse`] if the geometry is malformed.
/// - [`OcctExchangeError::Backend`] for other backend failures.
/// - [`OcctExchangeError::Io`] for filesystem failures.
pub fn step_ap214_reader(path: &Path) -> Result<Ap214Import, OcctExchangeError> {
    // Geometry via the proven AP203 backend (handles AP214 geometry).
    let solid = crate::step_ap203_reader(path)?;

    // Re-read the raw text to scan AP214 attributes. step_ap203_reader
    // already validated the extension and that the file parses.
    //
    // Round-21 M1: cap the re-read at MAX_CAD_INTERCHANGE_FILE_BYTES
    // (256 MiB — sister to the step_ap203_reader cap). Pre-fix this
    // was a bare `fs::read_to_string`; the AP203 reader already
    // capped its own read but a file that grew between the two
    // reads would have slurped unbounded here.
    let text = valenx_core::io_caps::read_capped_to_string(
        path,
        valenx_step_iges::MAX_CAD_INTERCHANGE_FILE_BYTES as usize,
    )?;
    let metadata = valenx_step_iges::ap242::parse_metadata(&text);
    let authority = scan_authority(&text);

    Ok(Ap214Import {
        solid,
        colors: metadata.colors,
        authority,
    })
}

/// Scan STEP text for AP214 PDM authority entities. Each recovered
/// string is prefixed with its entity kind so the caller can tell a
/// person from an approval-status.
fn scan_authority(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    for line in text.lines() {
        // Strip the leading `#id =` so keyword matching is positional.
        let body = line
            .split_once('=')
            .map(|(_, b)| b.trim_start())
            .unwrap_or(line);
        let kind = if body.starts_with("PERSON_AND_ORGANIZATION") {
            Some("person")
        } else if body.starts_with("APPROVAL_PERSON_ORGANIZATION") {
            Some("approval")
        } else if body.starts_with("APPROVAL_STATUS") {
            Some("approval-status")
        } else if body.starts_with("SECURITY_CLASSIFICATION_LEVEL") {
            Some("security")
        } else if body.starts_with("PERSON(") {
            Some("person")
        } else if body.starts_with("ORGANIZATION(") {
            Some("organization")
        } else {
            None
        };
        if let Some(kind) = kind {
            let label = first_quoted(body).unwrap_or_default();
            out.push(format!("{kind}: {label}"));
        }
    }
    out
}

/// Extract the first `'quoted'` substring from a STEP entity body.
fn first_quoted(s: &str) -> Option<String> {
    let start = s.find('\'')?;
    let end = s[start + 1..].find('\'')?;
    Some(s[start + 1..start + 1 + end].to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scan_authority_recovers_person_and_approval() {
        let text = "#10 = PERSON('p1', 'Ada', 'Lovelace', $, $, $);\n\
                    #11 = ORGANIZATION('o1', 'Analytical Engines', $);\n\
                    #12 = APPROVAL_STATUS('approved');\n\
                    #13 = SECURITY_CLASSIFICATION_LEVEL('confidential');\n\
                    #14 = CARTESIAN_POINT('', (0.0,0.0,0.0));\n";
        let auth = scan_authority(text);
        assert!(auth.iter().any(|a| a.starts_with("person:")));
        assert!(auth.iter().any(|a| a.starts_with("organization:")));
        assert!(auth.iter().any(|a| a == "approval-status: approved"));
        assert!(auth.iter().any(|a| a == "security: confidential"));
        // The geometry entity must not be picked up.
        assert!(!auth.iter().any(|a| a.contains("CARTESIAN")));
    }

    #[test]
    fn scan_authority_empty_when_geometry_only() {
        let text = "#1 = CARTESIAN_POINT('', (0.0,0.0,0.0));\n\
                    #2 = DIRECTION('', (0.0,0.0,1.0));\n";
        assert!(scan_authority(text).is_empty());
    }

    #[test]
    fn first_quoted_extracts_label() {
        assert_eq!(
            first_quoted("PERSON('p1', 'Ada', $);"),
            Some("p1".to_string())
        );
        assert_eq!(first_quoted("NO_QUOTES_HERE"), None);
    }

    #[test]
    fn reader_rejects_non_step_extension() {
        let err = step_ap214_reader(std::path::Path::new("a.obj")).unwrap_err();
        assert_eq!(err.code(), "occt_exchange.bad_input");
    }

    /// Round-21 M1 RED→GREEN: an oversize `.step` file is refused
    /// at the second-read attribute scan via
    /// `valenx_core::io_caps::read_capped_to_string` rather than
    /// being slurped into memory. (Allocating 257 MiB on disk in CI
    /// would be slow; verify the helper's behaviour against a small
    /// scratch file with a 1 KiB cap — the production code uses
    /// MAX_CAD_INTERCHANGE_FILE_BYTES = 256 MiB.)
    #[test]
    fn oversize_attribute_rescan_returns_invalid_data() {
        use std::io::Write;
        let tmp = std::env::temp_dir().join("valenx_r21_ap214_oversize.step");
        let mut f = std::fs::File::create(&tmp).unwrap();
        f.write_all(&vec![b'X'; 4096]).unwrap();
        drop(f);
        let err = valenx_core::io_caps::read_capped_to_string(&tmp, 1024).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        let _ = std::fs::remove_file(&tmp);
        // Sanity: the production cap is the 256 MiB CAD-interchange
        // sister cap, ensuring we use the shared shape rather than a
        // local one.
        assert_eq!(
            valenx_step_iges::MAX_CAD_INTERCHANGE_FILE_BYTES,
            256u64 * 1024 * 1024
        );
    }
}
