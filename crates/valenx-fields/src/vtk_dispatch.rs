//! Auto-dispatcher between the two VTK readers in this crate.
//!
//! VTK files come in two on-disk flavours:
//!
//! - **`.vtu` (XML)** — modern unstructured-grid wrapper, parsed by
//!   [`crate::vtu::parse_ascii`]. ASCII-DataArray-only today;
//!   appended-binary lands later.
//! - **`.vtk` legacy binary** — the older `# vtk DataFile Version`
//!   format mixing ASCII headers with raw big-endian binary blocks.
//!   Parsed by [`crate::vtk_legacy::parse_binary`].
//!
//! The two formats are unambiguously distinguishable by their first
//! few bytes, so this module sniffs the buffer and routes to the
//! right reader. Frees callers from having to branch on the file
//! extension (which is unreliable — solvers sometimes write XML to
//! `.vtk` and vice versa) or pre-load the bytes twice.

use thiserror::Error;

/// Which on-disk VTK flavour a buffer is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VtkFormat {
    /// XML wrapper — `.vtu` typically. Use [`crate::vtu::parse_ascii`].
    VtuXml,
    /// Legacy binary `# vtk DataFile Version` files.
    /// Use [`crate::vtk_legacy::parse_binary`].
    VtkLegacyBinary,
    /// Legacy ASCII `# vtk DataFile Version` files. We recognise
    /// the format but don't have a reader for it yet — callers get
    /// [`DispatchError::UnsupportedLegacyAscii`].
    VtkLegacyAscii,
}

/// Errors raised by the dispatcher's parse / detect helpers.
#[derive(Debug, Error)]
pub enum DispatchError {
    /// Input buffer is shorter than the format-magic prefix we need to
    /// classify it.
    #[error("buffer too short to identify VTK format ({0} bytes)")]
    TooShort(usize),
    /// Buffer didn't match any known VTK format magic.
    #[error("buffer doesn't look like any known VTK format")]
    Unrecognised,
    /// File was recognised as legacy ASCII VTK, which the dispatcher
    /// does not yet handle.
    #[error(
        "VTK legacy ASCII format isn't implemented yet — convert to VTK legacy binary or .vtu"
    )]
    UnsupportedLegacyAscii,
    /// `.vtu` input was not valid UTF-8.
    #[error("vtu requires valid UTF-8: {0}")]
    NotUtf8(String),
    /// Forwarding for [`crate::vtu::ParseError`].
    #[error("vtu parse: {0}")]
    Vtu(#[from] crate::vtu::ParseError),
    /// Forwarding for [`crate::vtk_legacy::ParseError`].
    #[error("vtk legacy parse: {0}")]
    VtkLegacy(#[from] crate::vtk_legacy::ParseError),
}

/// Inspect the first N bytes of a buffer and decide which VTK reader
/// to use. Returns `None` on buffers shorter than the magic prefix.
///
/// Detection rules (ordered):
/// 1. Starts with `# vtk DataFile Version` → VTK legacy. Walk to
///    the third line: `BINARY` → `VtkLegacyBinary`, `ASCII` →
///    `VtkLegacyAscii`. Anything else also returns `VtkLegacyBinary`
///    so the legacy parser can produce a more specific error.
/// 2. Starts with `<?xml` or contains `<VTKFile` in the first 256
///    bytes → `VtuXml`.
/// 3. Otherwise: `None`.
pub fn sniff(bytes: &[u8]) -> Option<VtkFormat> {
    const LEGACY_MAGIC: &[u8] = b"# vtk DataFile Version";
    if bytes.starts_with(LEGACY_MAGIC) {
        // Walk the first three lines; line 3 is the format keyword.
        let mut newlines = 0usize;
        let mut line3_start = 0usize;
        for (i, &b) in bytes.iter().enumerate() {
            if b == b'\n' {
                newlines += 1;
                if newlines == 2 {
                    line3_start = i + 1;
                }
                if newlines == 3 {
                    let line = &bytes[line3_start..i];
                    let trimmed = trim_ascii(line);
                    if trimmed == b"BINARY" {
                        return Some(VtkFormat::VtkLegacyBinary);
                    } else if trimmed == b"ASCII" {
                        return Some(VtkFormat::VtkLegacyAscii);
                    }
                    // Unknown format keyword — treat as binary so the
                    // legacy parser surfaces a structured error.
                    return Some(VtkFormat::VtkLegacyBinary);
                }
            }
        }
        // Truncated header; default to binary so the parser produces
        // a structured truncation error rather than us silently
        // returning None.
        return Some(VtkFormat::VtkLegacyBinary);
    }
    // VTU XML detection.
    let head_end = bytes.len().min(256);
    let head = &bytes[..head_end];
    if head.starts_with(b"<?xml") {
        return Some(VtkFormat::VtuXml);
    }
    if window_contains(head, b"<VTKFile") {
        return Some(VtkFormat::VtuXml);
    }
    None
}

/// High-level loader: sniff the format, hand off to the right
/// reader, return the canonical (Mesh, Fields) pair. One-line
/// upgrade for callers that don't care which format the file is.
///
/// The `mesh_id` is forwarded to the underlying canonical converter.
pub fn load_canonical(
    bytes: &[u8],
    mesh_id: impl Into<String>,
) -> Result<(valenx_mesh::Mesh, Vec<crate::Field>), DispatchError> {
    if bytes.len() < 4 {
        return Err(DispatchError::TooShort(bytes.len()));
    }
    match sniff(bytes) {
        Some(VtkFormat::VtuXml) => {
            // Detect appended-binary first because that path can't
            // round-trip through &str (the binary tail breaks UTF-8).
            // The substring check is cheap relative to a full parse.
            if window_contains(bytes, b"<AppendedData") {
                let data = crate::vtu::parse_appended_raw(bytes)?;
                return Ok(data.to_canonical(mesh_id));
            }
            let text =
                std::str::from_utf8(bytes).map_err(|e| DispatchError::NotUtf8(e.to_string()))?;
            let data = crate::vtu::parse_ascii(text)?;
            Ok(data.to_canonical(mesh_id))
        }
        Some(VtkFormat::VtkLegacyBinary) => {
            let data = crate::vtk_legacy::parse_binary(bytes)?;
            Ok(data.to_canonical(mesh_id))
        }
        Some(VtkFormat::VtkLegacyAscii) => Err(DispatchError::UnsupportedLegacyAscii),
        None => Err(DispatchError::Unrecognised),
    }
}

/// Small helper: substring search inside a byte slice. Faster than
/// pulling in `memchr` for one user.
fn window_contains(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() || needle.len() > haystack.len() {
        return false;
    }
    haystack.windows(needle.len()).any(|w| w == needle)
}

/// Strip leading + trailing ASCII whitespace from a byte slice. The
/// stdlib's `trim_ascii` exists on `[u8]` only on a recent toolchain;
/// we inline our own to keep the MSRV pin where it is.
fn trim_ascii(bytes: &[u8]) -> &[u8] {
    let mut start = 0usize;
    let mut end = bytes.len();
    while start < end && bytes[start].is_ascii_whitespace() {
        start += 1;
    }
    while end > start && bytes[end - 1].is_ascii_whitespace() {
        end -= 1;
    }
    &bytes[start..end]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sniff_recognises_legacy_binary_header() {
        let bytes = b"# vtk DataFile Version 3.0\ntitle\nBINARY\nDATASET UNSTRUCTURED_GRID\n";
        assert_eq!(sniff(bytes), Some(VtkFormat::VtkLegacyBinary));
    }

    #[test]
    fn sniff_recognises_legacy_ascii_header() {
        let bytes = b"# vtk DataFile Version 3.0\ntitle\nASCII\nDATASET UNSTRUCTURED_GRID\n";
        assert_eq!(sniff(bytes), Some(VtkFormat::VtkLegacyAscii));
    }

    #[test]
    fn sniff_recognises_vtu_xml_with_xml_declaration() {
        let bytes = b"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" />";
        assert_eq!(sniff(bytes), Some(VtkFormat::VtuXml));
    }

    #[test]
    fn sniff_recognises_vtu_xml_without_xml_declaration() {
        // Some writers omit the <?xml prologue. The fallback substring
        // search for <VTKFile must catch it.
        let bytes = b"<!-- comment -->\n<VTKFile type=\"UnstructuredGrid\" />";
        assert_eq!(sniff(bytes), Some(VtkFormat::VtuXml));
    }

    #[test]
    fn sniff_returns_none_on_arbitrary_garbage() {
        assert_eq!(sniff(b"not vtk at all"), None);
        assert_eq!(sniff(b""), None);
        assert_eq!(sniff(b"\x00\x00\x00\x00"), None);
    }

    #[test]
    fn sniff_handles_truncated_legacy_header() {
        // Header missing the BINARY/ASCII line — still classified as
        // legacy so the legacy parser can produce a Truncated error.
        let bytes = b"# vtk DataFile Version 3.0\ntitle\n";
        assert_eq!(sniff(bytes), Some(VtkFormat::VtkLegacyBinary));
    }

    #[test]
    fn sniff_rejects_legacy_with_unknown_format_keyword() {
        // Falls through to "treat as binary" so the parser surfaces
        // a structured error — better than silently returning None
        // and confusing the caller.
        let bytes = b"# vtk DataFile Version 3.0\ntitle\nUNKNOWN\nDATASET UNSTRUCTURED_GRID\n";
        assert_eq!(sniff(bytes), Some(VtkFormat::VtkLegacyBinary));
    }

    #[test]
    fn sniff_xml_check_only_scans_first_256_bytes() {
        // <VTKFile far beyond the head window -> miss. Confirms the
        // scan is bounded.
        let mut bytes = vec![b' '; 1024];
        bytes.extend_from_slice(b"<VTKFile />");
        assert_eq!(sniff(&bytes), None);
    }

    #[test]
    fn load_canonical_reads_a_minimal_legacy_binary_tet() {
        // Build the same minimal-tet payload the vtk_legacy module
        // tests with, but route through load_canonical to verify the
        // dispatch happens.
        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        buf.extend_from_slice(b"dispatch test\n");
        buf.extend_from_slice(b"BINARY\n");
        buf.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        buf.extend_from_slice(b"POINTS 4 float\n");
        for v in [
            0.0_f32, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
        ] {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        buf.push(b'\n');
        buf.extend_from_slice(b"CELLS 1 5\n");
        for v in [4u32, 0, 1, 2, 3] {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        buf.push(b'\n');
        buf.extend_from_slice(b"CELL_TYPES 1\n");
        buf.extend_from_slice(&10u32.to_be_bytes());
        buf.push(b'\n');
        let (mesh, fields) = load_canonical(&buf, "dispatch-tet").expect("load");
        assert_eq!(mesh.id, "dispatch-tet");
        assert_eq!(mesh.nodes.len(), 4);
        assert!(fields.is_empty());
    }

    #[test]
    fn load_canonical_rejects_legacy_ascii_today() {
        let bytes = b"# vtk DataFile Version 3.0\ntitle\nASCII\nDATASET UNSTRUCTURED_GRID\n";
        let err = load_canonical(bytes, "x").expect_err("must fail");
        assert!(matches!(err, DispatchError::UnsupportedLegacyAscii));
    }

    #[test]
    fn load_canonical_too_short_returns_too_short() {
        let bytes = b"";
        match load_canonical(bytes, "x") {
            Err(DispatchError::TooShort(n)) => assert_eq!(n, 0),
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn load_canonical_unrecognised_returns_unrecognised() {
        let bytes = b"random binary garbage that isn't vtk";
        match load_canonical(bytes, "x") {
            Err(DispatchError::Unrecognised) => {}
            other => panic!("wrong error: {other:?}"),
        }
    }
}
