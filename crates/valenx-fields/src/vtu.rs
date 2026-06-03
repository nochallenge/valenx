//! Minimal `.vtu` (VTK XML UnstructuredGrid) ASCII parser.
//!
//! VTK's XML format has an enormous surface area — appended-binary
//! data, parallel decomposition (PVTU), compressed binary, all sorts
//! of optional metadata. This parser handles the slice OpenFOAM /
//! Elmer / CalculiX actually emit when they write ASCII output:
//!
//! - `<VTKFile type="UnstructuredGrid" …>`
//! - one `<Piece>` with point/cell counts
//! - `<Points>` block with one `<DataArray>` of float coordinates
//! - `<Cells>` block with `connectivity` / `offsets` / `types` arrays
//! - `<PointData>` and/or `<CellData>` with named scalar / vector
//!   `<DataArray>` entries
//!
//! Binary appended data, base64 encoding, and PVTU all return
//! `ParseError::Unsupported` rather than partial / wrong output —
//! the contract is "ASCII in, real data out, nothing else".
//!
//! No XML library dependency. The format is rigid enough that
//! `find` / `strip_prefix` on the relevant element shapes is enough,
//! and avoiding the dep keeps `valenx-fields`' license perimeter
//! tight.

use std::collections::BTreeMap;

use thiserror::Error;

/// Round-11 hardening (R11-4): hard upper bound on `NumberOfPoints`
/// and `NumberOfCells` declared in the `<Piece>` element. A malicious
/// VTU with `NumberOfPoints="18446744073709551615"` would force the
/// loop that allocates the points / connectivity / offsets buffers to
/// ask the allocator for hundreds of exabytes via the (legitimate)
/// `n_points * 3` and `Vec::with_capacity(n_points)` patterns. 256 M
/// is far above any realistic mesh anyone would ship in a `.vtu`
/// (production solvers serialize meshes well below this — billion-cell
/// runs split into parallel PVTU pieces).
pub const MAX_VTU_POINTS: usize = 256 * 1024 * 1024;

/// Cell types we recognise. Codes match VTK's UnstructuredGrid
/// convention; see `vtkCellType.h` upstream. Anything outside this
/// set parses as `VtuCellType::Unknown(code)` so callers can decide
/// whether to skip or hard-fail.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VtuCellType {
    Vertex,   // 1
    Line,     // 3
    Triangle, // 5
    Quad,     // 9
    Tet,      // 10
    Hex,      // 12
    Wedge,    // 13 (also called Prism)
    Pyramid,  // 14
    Tri6,     // 22 (quadratic triangle)
    Tet10,    // 24 (quadratic tet)
    Hex20,    // 25 (quadratic hex)
    Unknown(u8),
}

impl VtuCellType {
    /// Decode a VTK cell-type integer code into the matching enum
    /// variant. Unknown codes fall through to [`VtuCellType::Unknown`].
    pub fn from_code(code: u8) -> Self {
        match code {
            1 => Self::Vertex,
            3 => Self::Line,
            5 => Self::Triangle,
            9 => Self::Quad,
            10 => Self::Tet,
            12 => Self::Hex,
            13 => Self::Wedge,
            14 => Self::Pyramid,
            22 => Self::Tri6,
            24 => Self::Tet10,
            25 => Self::Hex20,
            other => Self::Unknown(other),
        }
    }

    /// Return the VTK numeric code corresponding to this variant.
    /// [`VtuCellType::Unknown`] round-trips through its inner byte.
    pub fn code(self) -> u8 {
        match self {
            Self::Vertex => 1,
            Self::Line => 3,
            Self::Triangle => 5,
            Self::Quad => 9,
            Self::Tet => 10,
            Self::Hex => 12,
            Self::Wedge => 13,
            Self::Pyramid => 14,
            Self::Tri6 => 22,
            Self::Tet10 => 24,
            Self::Hex20 => 25,
            Self::Unknown(c) => c,
        }
    }
}

/// Mesh data extracted from a `.vtu` file. Layout matches the on-disk
/// VTK convention: a flat connectivity array indexed by per-cell
/// offsets, plus a parallel `cell_types` array.
#[derive(Clone, Debug, Default)]
pub struct VtuMesh {
    /// Point coordinates as `[x, y, z]` triples. VTK is always 3D
    /// even for 2D meshes (the unused axis is zero).
    pub points: Vec<[f64; 3]>,
    /// Flat connectivity buffer; cell `i` reads
    /// `connectivity[offsets[i-1] .. offsets[i]]` (with
    /// `offsets[-1] = 0`).
    pub connectivity: Vec<u32>,
    /// End-index in `connectivity` for each cell (cumulative).
    pub offsets: Vec<u32>,
    /// VTK cell-type codes, parallel to `offsets`.
    pub cell_types: Vec<VtuCellType>,
}

impl VtuMesh {
    /// How many cells the mesh holds.
    pub fn cell_count(&self) -> usize {
        self.offsets.len()
    }

    /// Iterator over `(cell_type, slice_of_node_ids)` for each cell.
    /// Node ids are 0-based (VTK's native indexing).
    ///
    /// # Panics
    ///
    /// Round-5 hardening: this iterator skips malformed cells rather
    /// than panicking. Specifically, when:
    /// - `cell_types.len() != offsets.len()` (mismatched arrays),
    /// - `offsets` is not monotonically non-decreasing, or
    /// - any offset exceeds `connectivity.len()`,
    ///
    /// the offending cell is dropped from the iterator. Strict callers
    /// can invoke [`Self::validate`] up front for an explicit error.
    pub fn cells(&self) -> impl Iterator<Item = (VtuCellType, &[u32])> + '_ {
        let conn_len = self.connectivity.len();
        let types_len = self.cell_types.len();
        self.offsets
            .iter()
            .enumerate()
            .filter_map(move |(i, &end)| {
                if i >= types_len {
                    return None;
                }
                let start = if i == 0 {
                    0
                } else {
                    self.offsets[i - 1] as usize
                };
                let end = end as usize;
                if end > conn_len || start > end {
                    return None;
                }
                Some((self.cell_types[i], &self.connectivity[start..end]))
            })
    }

    /// Validate the mesh's structural invariants. Round-5: returns
    /// `Err` when offsets are non-monotonic, when `cell_types.len()`
    /// disagrees with `offsets.len()`, or when an offset exceeds the
    /// connectivity buffer. Callers that need a hard failure mode
    /// (rather than the [`Self::cells`] iterator's filter-malformed
    /// behaviour) call this before iterating.
    pub fn validate(&self) -> Result<(), ParseError> {
        if self.cell_types.len() != self.offsets.len() {
            return Err(ParseError::CountMismatch {
                what: "Cells/types vs Cells/offsets",
                expected: self.offsets.len(),
                got: self.cell_types.len(),
            });
        }
        let mut prev = 0u64;
        for (i, &end) in self.offsets.iter().enumerate() {
            let end_u64 = end as u64;
            if end_u64 < prev {
                return Err(ParseError::BadAttribute(format!(
                    "Cells/offsets must be monotonically non-decreasing; \
                     offsets[{i}]={end} < previous {prev}"
                )));
            }
            if end_u64 > self.connectivity.len() as u64 {
                return Err(ParseError::CountMismatch {
                    what: "Cells/offsets entry",
                    expected: self.connectivity.len(),
                    got: end as usize,
                });
            }
            prev = end_u64;
        }
        Ok(())
    }
}

/// One named field from a `.vtu` `<PointData>` or `<CellData>` block.
#[derive(Clone, Debug)]
pub struct VtuField {
    pub name: String,
    /// 1 for scalar, 3 for vector, 9 for tensor (VTK convention).
    pub components: usize,
    /// Flat `[v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, …]` buffer.
    pub data: Vec<f64>,
}

impl VtuField {
    /// How many sample points the field carries (data length / components).
    pub fn samples(&self) -> usize {
        if self.components == 0 {
            0
        } else {
            self.data.len() / self.components
        }
    }
}

/// Top-level result of parsing a `.vtu` file.
#[derive(Clone, Debug, Default)]
pub struct VtuData {
    pub mesh: VtuMesh,
    /// Fields stored at mesh nodes.
    pub point_fields: Vec<VtuField>,
    /// Fields stored per cell.
    pub cell_fields: Vec<VtuField>,
}

/// Errors a `.vtu` parser can emit. We keep these structured so the
/// UI can render them as actionable diagnostics rather than just a
/// raw "parse failed" string.
#[derive(Debug, Error, PartialEq)]
pub enum ParseError {
    #[error("missing required element: {0}")]
    MissingElement(&'static str),
    #[error("missing required attribute `{attr}` on `{element}`")]
    MissingAttribute {
        element: &'static str,
        attr: &'static str,
    },
    #[error("unexpected attribute value: {0}")]
    BadAttribute(String),
    #[error("malformed numeric data in `{0}`: {1}")]
    BadNumeric(&'static str, String),
    #[error("count mismatch in `{what}`: expected {expected}, got {got}")]
    CountMismatch {
        what: &'static str,
        expected: usize,
        got: usize,
    },
    #[error(
        "unsupported VTU feature: {0} \
         (valenx-fields parses ASCII UnstructuredGrid only — appended-binary, \
          base64, compressed, or PVTU data needs a real VTK library)"
    )]
    Unsupported(&'static str),
}

/// Parse the ASCII shape of a `.vtu` file. Returns `Unsupported` for
/// the appended / compressed / parallel variants — call
/// [`parse_appended_raw`] for the appended-binary path.
pub fn parse_ascii(text: &str) -> Result<VtuData, ParseError> {
    // Reject formats we don't handle up-front so users get a clear
    // error instead of a corrupted parse.
    if text.contains("<AppendedData") {
        return Err(ParseError::Unsupported(
            "AppendedData section (use parse_appended_raw)",
        ));
    }
    if text.contains("compressor=") {
        return Err(ParseError::Unsupported("compressed payload"));
    }
    if text.contains("PUnstructuredGrid") {
        return Err(ParseError::Unsupported("parallel PVTU file"));
    }

    // Anchor on the UnstructuredGrid element.
    if !text.contains("<UnstructuredGrid") {
        return Err(ParseError::MissingElement("UnstructuredGrid"));
    }
    let piece_start = text
        .find("<Piece")
        .ok_or(ParseError::MissingElement("Piece"))?;
    let piece_end = text[piece_start..]
        .find("</Piece>")
        .ok_or(ParseError::MissingElement("Piece"))?
        + piece_start;
    let piece = &text[piece_start..piece_end];

    let n_points = checked_count(piece, "NumberOfPoints")?;
    let n_cells = checked_count(piece, "NumberOfCells")?;

    // -------- Points --------
    let points_section =
        slice_block(piece, "<Points>", "</Points>").ok_or(ParseError::MissingElement("Points"))?;
    let points_arr = first_data_array(points_section)?;
    let points_components = points_arr.components.max(1);
    if points_components != 3 {
        return Err(ParseError::BadAttribute(format!(
            "Points DataArray must have NumberOfComponents=3; got {points_components}"
        )));
    }
    let coords = parse_floats(points_arr.body, "Points")?;
    // Round-11: `n_points * 3` was unchecked pre-fix; the cap above
    // already keeps n_points within MAX_VTU_POINTS, so `* 3` cannot
    // overflow usize on any reasonable platform.
    let expected_coords = n_points
        .checked_mul(3)
        .ok_or(ParseError::BadAttribute(
            "Points: n_points * 3 overflows usize".to_string(),
        ))?;
    if coords.len() != expected_coords {
        return Err(ParseError::CountMismatch {
            what: "Points",
            expected: expected_coords,
            got: coords.len(),
        });
    }
    let mut points: Vec<[f64; 3]> = Vec::with_capacity(n_points);
    for chunk in coords.chunks_exact(3) {
        points.push([chunk[0], chunk[1], chunk[2]]);
    }

    // -------- Cells --------
    let cells_section =
        slice_block(piece, "<Cells>", "</Cells>").ok_or(ParseError::MissingElement("Cells"))?;
    let arrays = all_data_arrays(cells_section)?;
    let connectivity = arrays
        .get("connectivity")
        .ok_or(ParseError::MissingElement("Cells/connectivity"))?;
    let offsets = arrays
        .get("offsets")
        .ok_or(ParseError::MissingElement("Cells/offsets"))?;
    let types = arrays
        .get("types")
        .ok_or(ParseError::MissingElement("Cells/types"))?;
    let connectivity = parse_u32s(connectivity.body, "Cells/connectivity")?;
    let offsets = parse_u32s(offsets.body, "Cells/offsets")?;
    let types_raw = parse_u32s(types.body, "Cells/types")?;
    if offsets.len() != n_cells {
        return Err(ParseError::CountMismatch {
            what: "Cells/offsets",
            expected: n_cells,
            got: offsets.len(),
        });
    }
    if types_raw.len() != n_cells {
        return Err(ParseError::CountMismatch {
            what: "Cells/types",
            expected: n_cells,
            got: types_raw.len(),
        });
    }
    let cell_types: Vec<VtuCellType> = types_raw
        .iter()
        .map(|&c| VtuCellType::from_code((c & 0xff) as u8))
        .collect();
    let mesh = VtuMesh {
        points,
        connectivity,
        offsets,
        cell_types,
    };

    // -------- PointData (optional) --------
    let point_fields = parse_data_section(piece, "<PointData", "</PointData>", n_points)?;
    // -------- CellData (optional) --------
    let cell_fields = parse_data_section(piece, "<CellData", "</CellData>", n_cells)?;

    Ok(VtuData {
        mesh,
        point_fields,
        cell_fields,
    })
}

// ---------------------------------------------------------------------------
// Appended-raw binary parser
// ---------------------------------------------------------------------------

/// Parse a `.vtu` file in **appended-raw** binary format. The XML
/// prefix declares each `<DataArray format="appended" offset="N"/>`;
/// the raw data lives after the `<AppendedData encoding="raw">_`
/// marker as a sequence of `[u32 size header][little-endian payload]`
/// blocks indexed by `offset`.
///
/// Scope (matches what OpenFOAM 11+, Elmer, ParaView write by
/// default for binary export):
/// - `header_type="UInt32"` only — the v0 size header. v1.1+ also
///   uses `UInt64`; out of scope.
/// - DataArray types: `Float32`, `Float64`, `Int32`, `UInt32`,
///   `UInt8`, `Int8`. Anything else returns
///   [`ParseError::Unsupported`].
/// - No `compressor=` (raw uncompressed only). Compressed appended
///   payloads need zlib / lz4 deps; out of scope for v0.
///
/// On format detection failures (missing AppendedData, no `_`
/// separator) returns [`ParseError::Unsupported`] with a specific
/// message so the caller can route a fallback parser.
pub fn parse_appended_raw(bytes: &[u8]) -> Result<VtuData, ParseError> {
    // Locate the AppendedData section + the `_` data-start marker.
    let needle: &[u8] = b"<AppendedData";
    let appended_open =
        window_index(bytes, needle).ok_or(ParseError::MissingElement("AppendedData"))?;
    let after_open = appended_open + needle.len();
    let tag_close = bytes[after_open..]
        .iter()
        .position(|&b| b == b'>')
        .ok_or(ParseError::Unsupported("malformed AppendedData open tag"))?
        + after_open;
    // The text up to `<AppendedData` is the XML prefix the regular
    // parser handles. Everything after `_` is the raw payload.
    let xml_text = std::str::from_utf8(&bytes[..appended_open])
        .map_err(|e| ParseError::BadAttribute(format!("vtu xml prefix not UTF-8: {e}")))?;
    if xml_text.contains("compressor=") {
        return Err(ParseError::Unsupported(
            "compressed appended payload (zlib / lz4 not supported)",
        ));
    }
    // Validate the encoding attribute on the AppendedData tag itself.
    let appended_attrs = std::str::from_utf8(&bytes[after_open..tag_close]).unwrap_or("");
    let encoding = parse_attribute_str(appended_attrs, "encoding").unwrap_or("raw");
    if encoding != "raw" {
        return Err(ParseError::Unsupported(
            "only AppendedData encoding=\"raw\" is supported (base64 deferred)",
        ));
    }
    // Skip whitespace + find the `_` separator that marks the start
    // of the raw payload per the VTK spec.
    let underscore_pos = (tag_close + 1..bytes.len())
        .find(|&i| bytes[i] == b'_')
        .ok_or(ParseError::Unsupported(
            "missing `_` separator after <AppendedData ...>",
        ))?;
    let raw = &bytes[underscore_pos + 1..];

    // Walk the XML prefix the same way parse_ascii does, but each
    // DataArray's data comes from the raw section instead of the
    // text body.
    if !xml_text.contains("<UnstructuredGrid") {
        return Err(ParseError::MissingElement("UnstructuredGrid"));
    }
    let piece_start = xml_text
        .find("<Piece")
        .ok_or(ParseError::MissingElement("Piece"))?;
    let piece_end = xml_text[piece_start..]
        .find("</Piece>")
        .ok_or(ParseError::MissingElement("Piece"))?
        + piece_start;
    let piece = &xml_text[piece_start..piece_end];

    let n_points = checked_count(piece, "NumberOfPoints")?;
    let n_cells = checked_count(piece, "NumberOfCells")?;

    // Points
    let points_section =
        slice_block(piece, "<Points>", "</Points>").ok_or(ParseError::MissingElement("Points"))?;
    let points_attrs = first_data_array_attrs(points_section)?;
    let points_components = points_attrs.components.max(1);
    if points_components != 3 {
        return Err(ParseError::BadAttribute(format!(
            "Points DataArray must have NumberOfComponents=3; got {points_components}"
        )));
    }
    let coords = read_appended_floats(raw, &points_attrs, "Points")?;
    // Round-11: same `n_points * 3` overflow guard as parse_ascii.
    let expected_coords = n_points
        .checked_mul(3)
        .ok_or(ParseError::BadAttribute(
            "Points: n_points * 3 overflows usize".to_string(),
        ))?;
    if coords.len() != expected_coords {
        return Err(ParseError::CountMismatch {
            what: "Points",
            expected: expected_coords,
            got: coords.len(),
        });
    }
    let mut points: Vec<[f64; 3]> = Vec::with_capacity(n_points);
    for chunk in coords.chunks_exact(3) {
        points.push([chunk[0], chunk[1], chunk[2]]);
    }

    // Cells
    let cells_section =
        slice_block(piece, "<Cells>", "</Cells>").ok_or(ParseError::MissingElement("Cells"))?;
    let cells_arrays = all_data_array_attrs(cells_section)?;
    let connectivity_attrs = cells_arrays
        .get("connectivity")
        .ok_or(ParseError::MissingElement("Cells/connectivity"))?;
    let offsets_attrs = cells_arrays
        .get("offsets")
        .ok_or(ParseError::MissingElement("Cells/offsets"))?;
    let types_attrs = cells_arrays
        .get("types")
        .ok_or(ParseError::MissingElement("Cells/types"))?;
    let connectivity = read_appended_u32s(raw, connectivity_attrs, "Cells/connectivity")?;
    let offsets = read_appended_u32s(raw, offsets_attrs, "Cells/offsets")?;
    let types_raw = read_appended_u32s(raw, types_attrs, "Cells/types")?;
    if offsets.len() != n_cells {
        return Err(ParseError::CountMismatch {
            what: "Cells/offsets",
            expected: n_cells,
            got: offsets.len(),
        });
    }
    let cell_types: Vec<VtuCellType> = types_raw
        .iter()
        .map(|&c| VtuCellType::from_code((c & 0xff) as u8))
        .collect();
    let mesh = VtuMesh {
        points,
        connectivity,
        offsets,
        cell_types,
    };

    // PointData (optional)
    let point_fields =
        parse_appended_data_section(piece, raw, "<PointData", "</PointData>", n_points)?;
    let cell_fields = parse_appended_data_section(piece, raw, "<CellData", "</CellData>", n_cells)?;

    Ok(VtuData {
        mesh,
        point_fields,
        cell_fields,
    })
}

/// Attributes-only view of a DataArray — same as DataArrayBody but
/// without the inline body (which doesn't exist for appended).
#[derive(Clone, Debug)]
struct AppendedDataArrayAttrs {
    name: String,
    components: usize,
    type_name: String,
    offset: usize,
}

/// First DataArray's attributes inside a `<Points>` block.
fn first_data_array_attrs(section: &str) -> Result<AppendedDataArrayAttrs, ParseError> {
    iter_data_array_attrs(section)
        .next()
        .ok_or(ParseError::MissingElement("DataArray"))?
}

fn all_data_array_attrs(
    section: &str,
) -> Result<BTreeMap<String, AppendedDataArrayAttrs>, ParseError> {
    let mut out = BTreeMap::new();
    for arr in iter_data_array_attrs(section) {
        let arr = arr?;
        out.insert(arr.name.clone(), arr);
    }
    Ok(out)
}

fn iter_data_array_attrs(
    section: &str,
) -> impl Iterator<Item = Result<AppendedDataArrayAttrs, ParseError>> + '_ {
    let mut cursor = 0;
    std::iter::from_fn(move || {
        let rest = &section[cursor..];
        let open = rest.find("<DataArray")?;
        let after_open = open + "<DataArray".len();
        let close_tag = rest[after_open..].find('>')?;
        let attrs_end = after_open + close_tag;
        let attrs_str = &rest[after_open..attrs_end];
        // Self-closed `<DataArray ... />` is the appended variant;
        // body-bearing `<DataArray ...>...</DataArray>` is ASCII.
        // We accept either shape here — the format attribute is
        // what classifies them.
        let is_self_closed = rest[..attrs_end].ends_with('/')
            || rest.as_bytes().get(attrs_end) == Some(&b'>')
                && rest.as_bytes().get(attrs_end - 1) == Some(&b'/');
        let advance = if is_self_closed {
            attrs_end + 1
        } else {
            // Skip past the body so the next iteration starts after
            // `</DataArray>`.
            let body_start = attrs_end + 1;
            let body_end = rest[body_start..].find("</DataArray>")?;
            body_start + body_end + "</DataArray>".len()
        };
        cursor += advance;

        let name = match parse_attribute_str(attrs_str, "Name") {
            Some(n) => n.to_string(),
            None => "".to_string(),
        };
        let components = parse_attribute_u64(attrs_str, "NumberOfComponents").unwrap_or(1) as usize;
        let format = parse_attribute_str(attrs_str, "format").unwrap_or("ascii");
        if format != "appended" {
            return Some(Err(ParseError::Unsupported(
                "DataArray format must be `appended` in parse_appended_raw",
            )));
        }
        let type_name = match parse_attribute_str(attrs_str, "type") {
            Some(t) => t.to_string(),
            None => {
                return Some(Err(ParseError::MissingAttribute {
                    element: "DataArray",
                    attr: "type",
                }));
            }
        };
        let offset = match parse_attribute_u64(attrs_str, "offset") {
            Some(o) => o as usize,
            None => {
                return Some(Err(ParseError::MissingAttribute {
                    element: "DataArray",
                    attr: "offset",
                }));
            }
        };
        Some(Ok(AppendedDataArrayAttrs {
            name,
            components,
            type_name,
            offset,
        }))
    })
}

/// Read an appended DataArray's bytes as f64 (promoting f32 if
/// declared as Float32). Used for Points + scalar/vector fields.
fn read_appended_floats(
    raw: &[u8],
    arr: &AppendedDataArrayAttrs,
    what: &'static str,
) -> Result<Vec<f64>, ParseError> {
    let (block, elem_size) = read_appended_block(raw, arr, what)?;
    let mut out = Vec::with_capacity(block.len() / elem_size);
    match arr.type_name.as_str() {
        "Float32" => {
            for chunk in block.chunks_exact(4) {
                out.push(f32::from_le_bytes(chunk.try_into().unwrap()) as f64);
            }
        }
        "Float64" => {
            for chunk in block.chunks_exact(8) {
                out.push(f64::from_le_bytes(chunk.try_into().unwrap()));
            }
        }
        other => {
            return Err(ParseError::Unsupported(match other {
                "Int32" => "Int32 in a float-shaped DataArray",
                "UInt32" => "UInt32 in a float-shaped DataArray",
                _ => "non-float DataArray type for floats",
            }));
        }
    }
    Ok(out)
}

/// Read an appended DataArray's bytes as u32 (promoting smaller
/// integer types). Used for connectivity / offsets / types.
fn read_appended_u32s(
    raw: &[u8],
    arr: &AppendedDataArrayAttrs,
    what: &'static str,
) -> Result<Vec<u32>, ParseError> {
    let (block, _) = read_appended_block(raw, arr, what)?;
    let mut out = Vec::new();
    match arr.type_name.as_str() {
        "Int32" => {
            for chunk in block.chunks_exact(4) {
                out.push(i32::from_le_bytes(chunk.try_into().unwrap()) as u32);
            }
        }
        "UInt32" => {
            for chunk in block.chunks_exact(4) {
                out.push(u32::from_le_bytes(chunk.try_into().unwrap()));
            }
        }
        "Int64" => {
            for chunk in block.chunks_exact(8) {
                out.push(i64::from_le_bytes(chunk.try_into().unwrap()) as u32);
            }
        }
        "UInt64" => {
            for chunk in block.chunks_exact(8) {
                out.push(u64::from_le_bytes(chunk.try_into().unwrap()) as u32);
            }
        }
        "UInt8" => {
            // VTK writes cell types as UInt8. Promote to u32 so the
            // downstream cell-type lookup table can stay one shape.
            for &b in block {
                out.push(b as u32);
            }
        }
        "Int8" => {
            for &b in block {
                out.push(b as i8 as u32);
            }
        }
        _ => {
            return Err(ParseError::Unsupported("non-integer DataArray type"));
        }
    }
    Ok(out)
}

/// Common: read the size header at `offset` and return the data
/// block + the element-size of the DataArray's declared `type`.
///
/// Round-11 hardening (R11-3): pre-fix this routine did `arr.offset +
/// 4` and `block_start + size` without overflow checks. A malicious VTU
/// could declare `offset=usize::MAX - 2` so the addition wrapped to a
/// small positive number that *did* fit inside `raw.len()`, letting
/// the slice indexer return whatever bytes the underlying allocator
/// happened to have past the buffer. `checked_add` makes every offset
/// arithmetic explicit.
fn read_appended_block<'r>(
    raw: &'r [u8],
    arr: &AppendedDataArrayAttrs,
    what: &'static str,
) -> Result<(&'r [u8], usize), ParseError> {
    let block_start = arr.offset.checked_add(4).ok_or(ParseError::BadAttribute(
        format!("{what}: offset+4 overflows usize (offset={})", arr.offset),
    ))?;
    if block_start > raw.len() {
        return Err(ParseError::CountMismatch {
            what,
            expected: block_start,
            got: raw.len(),
        });
    }
    let size = u32::from_le_bytes(raw[arr.offset..block_start].try_into().unwrap()) as usize;
    let block_end = block_start
        .checked_add(size)
        .ok_or(ParseError::BadAttribute(format!(
            "{what}: block_start+size overflows usize \
             (block_start={block_start}, size={size})"
        )))?;
    if block_end > raw.len() {
        return Err(ParseError::CountMismatch {
            what,
            expected: block_end,
            got: raw.len(),
        });
    }
    let block = &raw[block_start..block_end];
    let elem_size = match arr.type_name.as_str() {
        "Float32" | "Int32" | "UInt32" => 4,
        "Float64" | "Int64" | "UInt64" => 8,
        "Int8" | "UInt8" => 1,
        _ => return Err(ParseError::Unsupported("unknown DataArray element type")),
    };
    Ok((block, elem_size))
}

/// Appended variant of `parse_data_section`. Walks every DataArray
/// in the PointData / CellData block and reads its bytes from `raw`
/// at the declared `offset`.
fn parse_appended_data_section(
    piece: &str,
    raw: &[u8],
    open_tag: &str,
    close_tag: &str,
    expected_samples: usize,
) -> Result<Vec<VtuField>, ParseError> {
    let Some(start_idx) = piece.find(open_tag) else {
        return Ok(Vec::new());
    };
    let close_of_open = piece[start_idx..].find('>').map(|i| i + start_idx);
    if let Some(close_of_open) = close_of_open {
        if piece[start_idx..=close_of_open].ends_with("/>") {
            return Ok(Vec::new());
        }
    }
    let section_end = piece[start_idx..]
        .find(close_tag)
        .ok_or(ParseError::MissingElement("PointData/CellData close"))?
        + start_idx;
    let section = &piece[start_idx..section_end];
    let arrays = iter_data_array_attrs(section);
    let mut out: Vec<VtuField> = Vec::new();
    for arr in arrays {
        let arr = arr?;
        let data = read_appended_floats(raw, &arr, "PointData/CellData")?;
        let expected = expected_samples * arr.components.max(1);
        if data.len() != expected {
            return Err(ParseError::CountMismatch {
                what: "PointData/CellData entry",
                expected,
                got: data.len(),
            });
        }
        out.push(VtuField {
            name: arr.name.clone(),
            components: arr.components.max(1),
            data,
        });
    }
    Ok(out)
}

/// Substring search inside a byte slice. Mirrors the helper in
/// vtk_dispatch — both are tiny enough to keep inlined per-crate.
fn window_index(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || needle.len() > haystack.len() {
        return None;
    }
    haystack.windows(needle.len()).position(|w| w == needle)
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// One `<DataArray>` element's parsed shape.
struct DataArrayBody<'a> {
    name: String,
    components: usize,
    format: String,
    body: &'a str,
}

/// Find the FIRST `<DataArray>` inside a section. Used for
/// `<Points>` which has only one (the coordinates).
fn first_data_array(section: &str) -> Result<DataArrayBody<'_>, ParseError> {
    iter_data_arrays(section)
        .next()
        .ok_or(ParseError::MissingElement("DataArray"))?
}

/// Collect every `<DataArray>` keyed by `Name=` attribute. Used for
/// `<Cells>` where we need to find connectivity / offsets / types.
fn all_data_arrays(section: &str) -> Result<BTreeMap<String, DataArrayBody<'_>>, ParseError> {
    let mut out = BTreeMap::new();
    for arr in iter_data_arrays(section) {
        let arr = arr?;
        out.insert(arr.name.clone(), arr);
    }
    Ok(out)
}

/// Iterate every `<DataArray>` element in a section. Each yields a
/// `Result` because attribute parsing can fail.
fn iter_data_arrays(section: &str) -> impl Iterator<Item = Result<DataArrayBody<'_>, ParseError>> {
    let mut cursor = 0;
    std::iter::from_fn(move || {
        let rest = &section[cursor..];
        let open = rest.find("<DataArray")?;
        let after_open = open + "<DataArray".len();
        let close_tag = rest[after_open..].find('>')?;
        let attrs_end = after_open + close_tag;
        let attrs_str = &rest[after_open..attrs_end];

        let body_start = attrs_end + 1;
        let body_end = rest[body_start..].find("</DataArray>")?;
        let body = &rest[body_start..body_start + body_end];
        cursor += body_start + body_end + "</DataArray>".len();

        let name = match parse_attribute_str(attrs_str, "Name") {
            Some(n) => n.to_string(),
            None => "".to_string(), // Points has no Name; allow blank
        };
        let components = parse_attribute_u64(attrs_str, "NumberOfComponents").unwrap_or(1) as usize;
        let format = parse_attribute_str(attrs_str, "format")
            .unwrap_or("ascii")
            .to_string();
        if format != "ascii" {
            return Some(Err(ParseError::Unsupported("non-ascii DataArray format")));
        }
        Some(Ok(DataArrayBody {
            name,
            components,
            format,
            body,
        }))
    })
}

/// Parse a `<PointData …>` or `<CellData …>` section if present.
/// Either tag may appear self-closed with no children, in which
/// case we return an empty Vec.
fn parse_data_section(
    piece: &str,
    open_tag: &str,
    close_tag: &str,
    expected_samples: usize,
) -> Result<Vec<VtuField>, ParseError> {
    let Some(start_idx) = piece.find(open_tag) else {
        return Ok(Vec::new());
    };
    // Self-closed `<PointData/>` or `<CellData/>` — no payload.
    let close_of_open = piece[start_idx..].find('>').map(|i| i + start_idx);
    if let Some(close_of_open) = close_of_open {
        if piece[start_idx..=close_of_open].ends_with("/>") {
            return Ok(Vec::new());
        }
    }
    let Some(end_idx) = piece[start_idx..].find(close_tag).map(|i| i + start_idx) else {
        return Ok(Vec::new());
    };
    let section = &piece[start_idx..end_idx];

    let mut fields: Vec<VtuField> = Vec::new();
    for arr in iter_data_arrays(section) {
        let arr = arr?;
        let data = parse_floats(arr.body, "DataArray body")?;
        let expected_len = expected_samples * arr.components.max(1);
        if data.len() != expected_len {
            return Err(ParseError::CountMismatch {
                what: "DataArray body",
                expected: expected_len,
                got: data.len(),
            });
        }
        fields.push(VtuField {
            name: arr.name,
            components: arr.components.max(1),
            data,
        });
        // Touch `format` so the field stays alive in error logs and
        // Clippy doesn't warn. We already validated it == "ascii"
        // inside `iter_data_arrays`.
        let _ = arr.format;
    }
    Ok(fields)
}

/// Find a substring bracketed by `open` and `close` within `text`.
/// Returns the inner slice (without the brackets). `None` if either
/// marker is absent.
fn slice_block<'a>(text: &'a str, open: &str, close: &str) -> Option<&'a str> {
    let start = text.find(open)? + open.len();
    let end = text[start..].find(close)? + start;
    Some(&text[start..end])
}

/// Pull `attr="value"` out of a tag's attribute string. Returns the
/// value without quotes.
fn parse_attribute_str<'a>(attrs: &'a str, attr: &str) -> Option<&'a str> {
    let needle_double = format!("{attr}=\"");
    let needle_single = format!("{attr}='");
    if let Some(start) = attrs.find(&needle_double) {
        let value_start = start + needle_double.len();
        let value_end = attrs[value_start..].find('"')? + value_start;
        return Some(&attrs[value_start..value_end]);
    }
    if let Some(start) = attrs.find(&needle_single) {
        let value_start = start + needle_single.len();
        let value_end = attrs[value_start..].find('\'')? + value_start;
        return Some(&attrs[value_start..value_end]);
    }
    None
}

fn parse_attribute_u64(attrs: &str, attr: &str) -> Option<u64> {
    parse_attribute_str(attrs, attr).and_then(|s| s.trim().parse::<u64>().ok())
}

/// Round-11 hardening (R11-4): read a count attribute (`NumberOfPoints` /
/// `NumberOfCells`) from a `<Piece>` element and cap it at
/// [`MAX_VTU_POINTS`]. Pre-fix a malicious VTU declaring
/// `NumberOfPoints="18446744073709551615"` would let the downstream
/// `n_points * 3` and `Vec::with_capacity(n_points)` ask the
/// allocator for impossible amounts of memory. The shared helper
/// keeps both `parse_ascii` and `parse_appended_raw` honest.
fn checked_count(piece: &str, attr: &'static str) -> Result<usize, ParseError> {
    let raw = parse_attribute_u64(piece, attr).ok_or(ParseError::MissingAttribute {
        element: "Piece",
        attr,
    })?;
    if raw > MAX_VTU_POINTS as u64 {
        return Err(ParseError::BadAttribute(format!(
            "{attr}={raw} exceeds MAX_VTU_POINTS={MAX_VTU_POINTS} cap (DoS protection)"
        )));
    }
    Ok(raw as usize)
}

/// Parse whitespace-separated f64 values out of a DataArray body.
fn parse_floats(body: &str, what: &'static str) -> Result<Vec<f64>, ParseError> {
    let mut out: Vec<f64> = Vec::new();
    for tok in body.split_ascii_whitespace() {
        let v: f64 = tok
            .parse()
            .map_err(|_| ParseError::BadNumeric(what, tok.to_string()))?;
        out.push(v);
    }
    Ok(out)
}

/// Parse whitespace-separated u32 values out of a DataArray body.
fn parse_u32s(body: &str, what: &'static str) -> Result<Vec<u32>, ParseError> {
    let mut out: Vec<u32> = Vec::new();
    for tok in body.split_ascii_whitespace() {
        let v: u32 = tok
            .parse()
            .map_err(|_| ParseError::BadNumeric(what, tok.to_string()))?;
        out.push(v);
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Canonical conversion — VtuData -> valenx_mesh::Mesh + valenx_fields::Field
// ---------------------------------------------------------------------------

/// Map a VTK cell type to the canonical `valenx_mesh::ElementType`.
///
/// Returns `None` for cell types valenx-mesh doesn't represent
/// (Vertex, Pyramid before it has a region, anything in
/// `VtuCellType::Unknown`). Callers typically skip those cells via
/// the per-block group-by below.
pub fn vtu_to_element_type(ct: VtuCellType) -> Option<valenx_mesh::ElementType> {
    use valenx_mesh::ElementType as ET;
    match ct {
        VtuCellType::Line => Some(ET::Line2),
        VtuCellType::Triangle => Some(ET::Tri3),
        VtuCellType::Quad => Some(ET::Quad4),
        VtuCellType::Tet => Some(ET::Tet4),
        VtuCellType::Hex => Some(ET::Hex8),
        VtuCellType::Wedge => Some(ET::Prism6),
        VtuCellType::Pyramid => Some(ET::Pyr5),
        VtuCellType::Tri6 => Some(ET::Tri6),
        VtuCellType::Tet10 => Some(ET::Tet10),
        VtuCellType::Hex20 => Some(ET::Hex20),
        // Vertex (1) has no meaningful canonical mapping — points are
        // already in the Mesh::nodes array. Unknown codes likewise fall
        // through silently rather than fail-on-import.
        VtuCellType::Vertex | VtuCellType::Unknown(_) => None,
    }
}

impl VtuData {
    /// Convert this `VtuData` into canonical valenx types: a
    /// [`valenx_mesh::Mesh`] (one [`valenx_mesh::ElementBlock`] per
    /// distinct cell type) plus a list of [`crate::Field`]s for the
    /// point + cell data.
    ///
    /// `mesh_id` is the identifier the new `Mesh` carries — typically
    /// the source case name or workdir basename so users can tell
    /// which run produced the mesh.
    ///
    /// Cells with no canonical mapping ([`VtuCellType::Vertex`],
    /// [`VtuCellType::Unknown`]) are silently skipped — they're rare
    /// in real solver output and forcing an error here would block
    /// every result-load on weird-but-harmless data.
    pub fn to_canonical(
        &self,
        mesh_id: impl Into<String>,
    ) -> (valenx_mesh::Mesh, Vec<crate::Field>) {
        use std::collections::BTreeMap;
        use valenx_mesh::{ElementBlock, Mesh};

        let mut mesh = Mesh::new(mesh_id);
        // Coordinates: f64 [x,y,z] -> Vector3<f64>.
        mesh.nodes = self
            .mesh
            .points
            .iter()
            .map(|p| nalgebra::Vector3::new(p[0], p[1], p[2]))
            .collect();

        // Group cells by canonical element type. We can't pre-allocate
        // the connectivity vectors because we don't know the per-type
        // count up-front; the BTreeMap ordering is fine because typical
        // meshes have only a few distinct types (often just one).
        let mut blocks: BTreeMap<u8, ElementBlock> = BTreeMap::new();
        for (ct, nodes) in self.mesh.cells() {
            let Some(canonical) = vtu_to_element_type(ct) else {
                // Skip Vertex / Unknown cells — they don't have a
                // canonical representation, and dropping them keeps
                // the mesh consistent (count-of-cells matches sum-of-
                // block-counts).
                continue;
            };
            let entry = blocks
                .entry(canonical as u8)
                .or_insert_with(|| ElementBlock::new(canonical));
            entry.connectivity.extend_from_slice(nodes);
        }
        mesh.element_blocks = blocks.into_values().collect();
        // Populate counts AND quality scalars in one pass so consumers
        // reading mesh.stats see AR / skewness / orthogonality without
        // needing a follow-up quality_report call.
        let _ = mesh.recompute_quality_stats();

        // Convert every PointData / CellData entry into a canonical
        // Field. VTU doesn't carry units or a time-step index, so we
        // default to dimensionless + steady; downstream consumers can
        // re-stamp these as needed.
        let mut fields: Vec<crate::Field> = Vec::new();
        for vf in &self.point_fields {
            fields.push(vtu_to_field(vf, crate::Location::OnNode));
        }
        for vf in &self.cell_fields {
            fields.push(vtu_to_field(vf, crate::Location::OnCell));
        }

        (mesh, fields)
    }
}

/// Convert one `VtuField` into a canonical [`crate::Field`]. Picks a
/// `FieldKind` from the component count: 1 → Scalar, 3 → Vector,
/// 9 → Tensor 3×3, anything else falls back to a flat 1-D vector
/// of `n` components so the data isn't silently lost.
fn vtu_to_field(vf: &VtuField, location: crate::Location) -> crate::Field {
    let kind = match vf.components {
        1 => crate::FieldKind::Scalar,
        3 => crate::FieldKind::Vector { dim: 3 },
        9 => crate::FieldKind::Tensor { rows: 3, cols: 3 },
        n => crate::FieldKind::Vector { dim: n as u8 },
    };
    let range = field_range(&vf.data);
    crate::Field {
        name: vf.name.clone(),
        kind,
        location,
        // VTU has no notion of named regions — the whole dataset is
        // implicitly one region. Downstream consumers can split this
        // further if they parse mesh boundary metadata separately.
        region: crate::RegionRef("default".to_string()),
        // VTU files don't carry units. Mark the field dimensionless;
        // adapters that know the physics (e.g. CFD fields are
        // velocity / pressure / temperature) can re-stamp these via
        // a follow-up pass before serving them to the report layer.
        units: crate::units::DIMENSIONLESS,
        // Steady by default — VTU encodes a single timestep per file.
        // Multi-timestep series come in as a sequence of separate
        // .vtu files indexed externally (e.g. PVD), which would
        // re-stamp these to TimeKey::Time { ... } at load time.
        time: crate::TimeKey::Steady,
        data: vf.data.clone(),
        range,
    }
}

/// Min/max across every f64 in the buffer. `None` for empty data.
fn field_range(data: &[f64]) -> Option<(f64, f64)> {
    let mut iter = data.iter();
    let &first = iter.next()?;
    let mut min = first;
    let mut max = first;
    for &v in iter {
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
    }
    Some((min, max))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal one-tet mesh with a scalar point field. Matches the
    /// shape OpenFOAM's foamToVTK writes for a tiny test case.
    fn one_tet_with_pressure() -> &'static str {
        r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="4" NumberOfCells="1">
      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
          0 0 0
          1 0 0
          0 1 0
          0 0 1
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">0 1 2 3</DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">4</DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">10</DataArray>
      </Cells>
      <PointData Scalars="p">
        <DataArray type="Float32" Name="p" NumberOfComponents="1" format="ascii">
          1.0 2.0 3.0 4.0
        </DataArray>
        <DataArray type="Float32" Name="U" NumberOfComponents="3" format="ascii">
          1 0 0  2 0 0  3 0 0  4 0 0
        </DataArray>
      </PointData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>"#
    }

    #[test]
    fn parses_one_tet_with_scalar_and_vector_fields() {
        let data = parse_ascii(one_tet_with_pressure()).expect("parse");
        // Mesh
        assert_eq!(data.mesh.points.len(), 4);
        assert_eq!(data.mesh.points[0], [0.0, 0.0, 0.0]);
        assert_eq!(data.mesh.points[3], [0.0, 0.0, 1.0]);
        assert_eq!(data.mesh.cell_count(), 1);
        assert_eq!(data.mesh.cell_types[0], VtuCellType::Tet);
        // Iterate cells: one tet with nodes [0,1,2,3].
        let cells: Vec<_> = data.mesh.cells().collect();
        assert_eq!(cells.len(), 1);
        assert_eq!(cells[0].0, VtuCellType::Tet);
        assert_eq!(cells[0].1, &[0u32, 1, 2, 3]);
        // Point data — both scalar p and vector U.
        assert_eq!(data.point_fields.len(), 2);
        let p = data
            .point_fields
            .iter()
            .find(|f| f.name == "p")
            .expect("p field");
        assert_eq!(p.components, 1);
        assert_eq!(p.samples(), 4);
        assert_eq!(p.data, &[1.0, 2.0, 3.0, 4.0]);
        let u = data
            .point_fields
            .iter()
            .find(|f| f.name == "U")
            .expect("U field");
        assert_eq!(u.components, 3);
        assert_eq!(u.samples(), 4);
        assert_eq!(
            u.data,
            &[1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0, 4.0, 0.0, 0.0]
        );
        assert!(data.cell_fields.is_empty());
    }

    #[test]
    fn rejects_appended_data_block_explicitly() {
        let text = r#"<VTKFile type="UnstructuredGrid">
<UnstructuredGrid><Piece NumberOfPoints="0" NumberOfCells="0"></Piece></UnstructuredGrid>
<AppendedData encoding="raw">_garbage</AppendedData>
</VTKFile>"#;
        let err = parse_ascii(text).unwrap_err();
        assert!(matches!(err, ParseError::Unsupported(_)), "got {err:?}");
    }

    #[test]
    fn rejects_compressed_payload() {
        let text = r#"<VTKFile type="UnstructuredGrid" compressor="vtkZLibDataCompressor">
</VTKFile>"#;
        let err = parse_ascii(text).unwrap_err();
        assert!(matches!(err, ParseError::Unsupported(_)));
    }

    #[test]
    fn rejects_parallel_pvtu() {
        let text = r#"<VTKFile type="PUnstructuredGrid"></VTKFile>"#;
        let err = parse_ascii(text).unwrap_err();
        assert!(matches!(err, ParseError::Unsupported(_)));
    }

    #[test]
    fn rejects_missing_unstructuredgrid() {
        let text = r#"<VTKFile type="StructuredGrid"></VTKFile>"#;
        let err = parse_ascii(text).unwrap_err();
        assert_eq!(err, ParseError::MissingElement("UnstructuredGrid"));
    }

    #[test]
    fn rejects_count_mismatch_in_points() {
        // Says NumberOfPoints=4 but only writes 3 coordinates' worth.
        let text = r#"<VTKFile type="UnstructuredGrid">
<UnstructuredGrid><Piece NumberOfPoints="4" NumberOfCells="0">
  <Points><DataArray type="Float32" NumberOfComponents="3" format="ascii">
    0 0 0  1 0 0  0 1 0
  </DataArray></Points>
  <Cells>
    <DataArray Name="connectivity" format="ascii"></DataArray>
    <DataArray Name="offsets" format="ascii"></DataArray>
    <DataArray Name="types" format="ascii"></DataArray>
  </Cells>
</Piece></UnstructuredGrid></VTKFile>"#;
        let err = parse_ascii(text).unwrap_err();
        match err {
            ParseError::CountMismatch {
                what,
                expected,
                got,
            } => {
                assert_eq!(what, "Points");
                assert_eq!(expected, 12);
                assert_eq!(got, 9);
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn cell_type_codes_round_trip() {
        for code in [1u8, 3, 5, 9, 10, 12, 13, 14, 22, 24, 25] {
            let ct = VtuCellType::from_code(code);
            assert_ne!(ct, VtuCellType::Unknown(code));
            assert_eq!(ct.code(), code);
        }
        // Unknown code preserves itself.
        let ct = VtuCellType::from_code(99);
        assert_eq!(ct, VtuCellType::Unknown(99));
        assert_eq!(ct.code(), 99);
    }

    #[test]
    fn empty_pointdata_section_is_ok() {
        // A self-closed <PointData/> block is valid VTU and should
        // produce an empty point_fields vec, not an error.
        let text = r#"<VTKFile type="UnstructuredGrid">
<UnstructuredGrid><Piece NumberOfPoints="2" NumberOfCells="1">
  <Points><DataArray type="Float32" NumberOfComponents="3" format="ascii">
    0 0 0  1 0 0
  </DataArray></Points>
  <Cells>
    <DataArray Name="connectivity" format="ascii">0 1</DataArray>
    <DataArray Name="offsets" format="ascii">2</DataArray>
    <DataArray Name="types" format="ascii">3</DataArray>
  </Cells>
  <PointData/>
</Piece></UnstructuredGrid></VTKFile>"#;
        let data = parse_ascii(text).expect("parse");
        assert_eq!(data.mesh.points.len(), 2);
        assert_eq!(data.mesh.cell_count(), 1);
        assert_eq!(data.mesh.cell_types[0], VtuCellType::Line);
        assert!(data.point_fields.is_empty());
    }

    #[test]
    fn to_canonical_one_tet_produces_real_mesh_and_fields() {
        let data = parse_ascii(one_tet_with_pressure()).expect("parse");
        let (mesh, fields) = data.to_canonical("smoke");
        assert_eq!(mesh.id, "smoke");
        assert_eq!(mesh.nodes.len(), 4);
        // First node = origin.
        assert_eq!(mesh.nodes[0].x, 0.0);
        assert_eq!(mesh.nodes[3].z, 1.0);
        // One element block holding the single tet.
        assert_eq!(mesh.element_blocks.len(), 1);
        assert_eq!(
            mesh.element_blocks[0].element_type,
            valenx_mesh::ElementType::Tet4
        );
        assert_eq!(mesh.element_blocks[0].connectivity, vec![0, 1, 2, 3]);
        assert_eq!(mesh.stats.node_count, 4);
        assert_eq!(mesh.stats.element_count, 1);

        // Two point fields: scalar p and vector U.
        assert_eq!(fields.len(), 2);
        let p = fields.iter().find(|f| f.name == "p").expect("p field");
        assert_eq!(p.kind, crate::FieldKind::Scalar);
        assert_eq!(p.location, crate::Location::OnNode);
        assert_eq!(p.region.0, "default");
        assert_eq!(p.data, &[1.0, 2.0, 3.0, 4.0]);
        // Range cached: scalar p ∈ [1, 4].
        assert_eq!(p.range, Some((1.0, 4.0)));

        let u = fields.iter().find(|f| f.name == "U").expect("U field");
        assert_eq!(u.kind, crate::FieldKind::Vector { dim: 3 });
        assert_eq!(u.location, crate::Location::OnNode);
        // U.x ranges from 1 to 4; cached range is min/max across ALL
        // components, not per-component, so y/z zeros pull min to 0.
        assert_eq!(u.range, Some((0.0, 4.0)));
    }

    #[test]
    fn to_canonical_skips_vertex_and_unknown_cells() {
        // Hand-built case: one tet (cell type 10 = Tet) plus one
        // vertex (cell type 1) which has no canonical mapping. The
        // vertex should be dropped silently and the tet should land
        // in its own block.
        let text = r#"<VTKFile type="UnstructuredGrid">
<UnstructuredGrid><Piece NumberOfPoints="5" NumberOfCells="2">
  <Points><DataArray type="Float32" NumberOfComponents="3" format="ascii">
    0 0 0  1 0 0  0 1 0  0 0 1  2 0 0
  </DataArray></Points>
  <Cells>
    <DataArray Name="connectivity" format="ascii">0 1 2 3 4</DataArray>
    <DataArray Name="offsets" format="ascii">4 5</DataArray>
    <DataArray Name="types" format="ascii">10 1</DataArray>
  </Cells>
</Piece></UnstructuredGrid></VTKFile>"#;
        let data = parse_ascii(text).expect("parse");
        let (mesh, fields) = data.to_canonical("with-vertex");
        assert!(fields.is_empty());
        // One block (the tet); vertex got dropped.
        assert_eq!(mesh.element_blocks.len(), 1);
        assert_eq!(
            mesh.element_blocks[0].element_type,
            valenx_mesh::ElementType::Tet4
        );
        assert_eq!(mesh.stats.element_count, 1);
        // Nodes are kept as-is — node 4 (the unused vertex point)
        // is still in the array because dropping it would invalidate
        // the connectivity indices of the kept tet.
        assert_eq!(mesh.nodes.len(), 5);
    }

    #[test]
    fn to_canonical_groups_by_cell_type() {
        // Two tets + one hex → expect two ElementBlocks.
        let text = r#"<VTKFile type="UnstructuredGrid">
<UnstructuredGrid><Piece NumberOfPoints="13" NumberOfCells="3">
  <Points><DataArray type="Float32" NumberOfComponents="3" format="ascii">
    0 0 0  1 0 0  0 1 0  0 0 1
    0 0 2  1 0 2  0 1 2
    0 0 3  1 0 3  1 1 3  0 1 3  0 0 4  1 1 4
  </DataArray></Points>
  <Cells>
    <DataArray Name="connectivity" format="ascii">
      0 1 2 3
      4 5 6 0
      7 8 9 10 11 12 0 1
    </DataArray>
    <DataArray Name="offsets" format="ascii">4 8 16</DataArray>
    <DataArray Name="types" format="ascii">10 10 12</DataArray>
  </Cells>
</Piece></UnstructuredGrid></VTKFile>"#;
        let data = parse_ascii(text).expect("parse");
        let (mesh, _) = data.to_canonical("mixed");
        assert_eq!(mesh.element_blocks.len(), 2);
        // Tets go in one block of length 8 (2 cells × 4 nodes).
        let tet_block = mesh
            .element_blocks
            .iter()
            .find(|b| b.element_type == valenx_mesh::ElementType::Tet4)
            .expect("tet block");
        assert_eq!(tet_block.connectivity.len(), 8);
        // Hex is one cell × 8 nodes.
        let hex_block = mesh
            .element_blocks
            .iter()
            .find(|b| b.element_type == valenx_mesh::ElementType::Hex8)
            .expect("hex block");
        assert_eq!(hex_block.connectivity.len(), 8);
        assert_eq!(mesh.stats.element_count, 3);
    }

    #[test]
    fn vtu_to_element_type_covers_supported_set() {
        use valenx_mesh::ElementType as ET;
        assert_eq!(vtu_to_element_type(VtuCellType::Tet), Some(ET::Tet4));
        assert_eq!(vtu_to_element_type(VtuCellType::Hex), Some(ET::Hex8));
        assert_eq!(vtu_to_element_type(VtuCellType::Wedge), Some(ET::Prism6));
        assert_eq!(vtu_to_element_type(VtuCellType::Tri6), Some(ET::Tri6));
        assert_eq!(vtu_to_element_type(VtuCellType::Tet10), Some(ET::Tet10));
        assert_eq!(vtu_to_element_type(VtuCellType::Hex20), Some(ET::Hex20));
        // Vertex + Unknown have no canonical mapping.
        assert_eq!(vtu_to_element_type(VtuCellType::Vertex), None);
        assert_eq!(vtu_to_element_type(VtuCellType::Unknown(99)), None);
    }

    #[test]
    fn field_range_handles_empty_and_constant_data() {
        assert_eq!(field_range(&[]), None);
        assert_eq!(field_range(&[3.0]), Some((3.0, 3.0)));
        assert_eq!(field_range(&[5.0, 5.0, 5.0]), Some((5.0, 5.0)));
        assert_eq!(field_range(&[-1.0, 0.0, 2.5, -3.0]), Some((-3.0, 2.5)));
    }

    /// Round-5 RED→GREEN: `VtuMesh::cells` used to panic with slice
    /// OOB when `offsets` was non-monotonic (an offset earlier in
    /// the array was larger than a later one). The fix is the new
    /// filter-malformed iterator + `VtuMesh::validate` for callers
    /// that want a hard error.
    #[test]
    fn cells_rejects_non_monotonic_offsets() {
        let mesh = VtuMesh {
            points: vec![[0.0; 3]; 4],
            connectivity: vec![0u32, 1, 2, 3],
            // Non-monotonic: 4 (legitimate end of one tet) then 2
            // (rewinds past previous end) — used to panic in `cells()`.
            offsets: vec![4u32, 2u32],
            cell_types: vec![VtuCellType::Tet, VtuCellType::Line],
        };
        // Iterator must not panic — malformed entries are dropped.
        let collected: Vec<_> = mesh.cells().collect();
        // The first cell (offsets[0]=4) is well-formed; the second
        // (offsets[1]=2 < previous 4) is malformed and dropped.
        assert_eq!(collected.len(), 1);
        // Strict `validate()` flags it explicitly.
        let err = mesh.validate().expect_err("validate must reject");
        let msg = format!("{err}");
        assert!(
            msg.contains("monotonic") || msg.contains("non-decreasing"),
            "msg: {msg}"
        );
    }

    /// Round-5: `cell_types.len() != offsets.len()` is another
    /// hostile shape the pre-fix iterator panicked on.
    #[test]
    fn cells_rejects_cell_types_length_mismatch() {
        let mesh = VtuMesh {
            points: vec![[0.0; 3]; 4],
            connectivity: vec![0u32, 1, 2, 3],
            offsets: vec![4u32],
            // Two cell-type entries for one offset — malformed.
            cell_types: vec![VtuCellType::Tet, VtuCellType::Line],
        };
        // Iterator drops the second entry (no parallel offset).
        let collected: Vec<_> = mesh.cells().collect();
        assert_eq!(collected.len(), 1);
        // validate() flags it explicitly.
        let err = mesh.validate().expect_err("validate must reject");
        match err {
            ParseError::CountMismatch { what, .. } => {
                assert!(what.contains("Cells/types"), "what: {what}");
            }
            other => panic!("expected CountMismatch, got {other:?}"),
        }
    }

    /// Round-5: an offset that overruns `connectivity.len()` used to
    /// slice past the buffer and panic.
    #[test]
    fn cells_rejects_offset_past_connectivity_end() {
        let mesh = VtuMesh {
            points: vec![[0.0; 3]; 4],
            connectivity: vec![0u32, 1, 2, 3], // 4 entries
            offsets: vec![10u32],               // claims 10 — past end
            cell_types: vec![VtuCellType::Tet],
        };
        let collected: Vec<_> = mesh.cells().collect();
        assert!(collected.is_empty(), "cell with overrun offset must be dropped");
        let err = mesh.validate().expect_err("validate must reject");
        let msg = format!("{err}");
        assert!(msg.contains("Cells/offsets"), "msg: {msg}");
    }

    #[test]
    fn cell_data_field_parses() {
        let text = r#"<VTKFile type="UnstructuredGrid">
<UnstructuredGrid><Piece NumberOfPoints="3" NumberOfCells="1">
  <Points><DataArray type="Float32" NumberOfComponents="3" format="ascii">
    0 0 0  1 0 0  0 1 0
  </DataArray></Points>
  <Cells>
    <DataArray Name="connectivity" format="ascii">0 1 2</DataArray>
    <DataArray Name="offsets" format="ascii">3</DataArray>
    <DataArray Name="types" format="ascii">5</DataArray>
  </Cells>
  <CellData Scalars="cellId">
    <DataArray type="Int32" Name="cellId" NumberOfComponents="1" format="ascii">42</DataArray>
  </CellData>
</Piece></UnstructuredGrid></VTKFile>"#;
        let data = parse_ascii(text).expect("parse");
        assert_eq!(data.cell_fields.len(), 1);
        assert_eq!(data.cell_fields[0].name, "cellId");
        assert_eq!(data.cell_fields[0].data, &[42.0]);
        assert_eq!(data.mesh.cell_types[0], VtuCellType::Triangle);
    }

    // -----------------------------------------------------------------
    // Appended-raw binary parser
    // -----------------------------------------------------------------

    /// Synthesize a minimal VTU appended-raw file: 1 tetrahedron
    /// + 1 scalar PointData "T". Used as the happy-path fixture.
    fn synth_appended_tet() -> Vec<u8> {
        // Build the binary tail first so we know each block's offset.
        // VTU spec: each block = u32 size header (LE) + payload.
        let mut binary: Vec<u8> = Vec::new();

        // Block 0: Points — 4 nodes × 3 components × 4 bytes = 48
        let points_offset = binary.len();
        let points: [f32; 12] = [
            0.0, 0.0, 0.0, // p0
            1.0, 0.0, 0.0, // p1
            0.0, 1.0, 0.0, // p2
            0.0, 0.0, 1.0, // p3
        ];
        binary.extend_from_slice(&((4 * 3 * 4) as u32).to_le_bytes());
        for v in points {
            binary.extend_from_slice(&v.to_le_bytes());
        }

        // Block 1: connectivity — 4 i32 values = 16 bytes
        let conn_offset = binary.len();
        binary.extend_from_slice(&((4 * 4) as u32).to_le_bytes());
        for v in [0i32, 1, 2, 3] {
            binary.extend_from_slice(&v.to_le_bytes());
        }

        // Block 2: offsets — 1 i32 value = 4 bytes
        let off_offset = binary.len();
        binary.extend_from_slice(&(4u32).to_le_bytes());
        binary.extend_from_slice(&4i32.to_le_bytes());

        // Block 3: types — 1 u8 value = 1 byte
        let types_offset = binary.len();
        binary.extend_from_slice(&(1u32).to_le_bytes());
        binary.push(10u8); // VTK_TETRA

        // Block 4: T scalar — 4 f32 values = 16 bytes
        let t_offset = binary.len();
        binary.extend_from_slice(&((4 * 4) as u32).to_le_bytes());
        for v in [273.15_f32, 280.0, 290.0, 300.0] {
            binary.extend_from_slice(&v.to_le_bytes());
        }

        let xml = format!(
            r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32">
<UnstructuredGrid>
<Piece NumberOfPoints="4" NumberOfCells="1">
<Points>
  <DataArray type="Float32" NumberOfComponents="3" format="appended" offset="{points_offset}"/>
</Points>
<Cells>
  <DataArray type="Int32" Name="connectivity" format="appended" offset="{conn_offset}"/>
  <DataArray type="Int32" Name="offsets" format="appended" offset="{off_offset}"/>
  <DataArray type="UInt8" Name="types" format="appended" offset="{types_offset}"/>
</Cells>
<PointData Scalars="T">
  <DataArray type="Float32" Name="T" NumberOfComponents="1" format="appended" offset="{t_offset}"/>
</PointData>
</Piece>
</UnstructuredGrid>
<AppendedData encoding="raw">
_"#
        );
        let mut out: Vec<u8> = xml.into_bytes();
        out.extend_from_slice(&binary);
        out.extend_from_slice(b"\n</AppendedData>\n</VTKFile>\n");
        out
    }

    #[test]
    fn parse_appended_raw_handles_a_minimal_tet() {
        let bytes = synth_appended_tet();
        let data = parse_appended_raw(&bytes).expect("parse");
        assert_eq!(data.mesh.points.len(), 4);
        assert_eq!(data.mesh.points[1], [1.0, 0.0, 0.0]);
        assert_eq!(data.mesh.cell_types, vec![VtuCellType::Tet]);
        assert_eq!(data.mesh.connectivity, vec![0, 1, 2, 3]);
        assert_eq!(data.mesh.offsets, vec![4]);
        assert_eq!(data.point_fields.len(), 1);
        assert_eq!(data.point_fields[0].name, "T");
        assert!((data.point_fields[0].data[0] - 273.15).abs() < 1e-3);
    }

    #[test]
    fn parse_appended_raw_rejects_compressed_payloads() {
        // Inject a `compressor="vtkZLibDataCompressor"` attribute on
        // VTKFile so the early-bail check trips.
        let mut bytes = synth_appended_tet();
        // Find VTKFile open tag and inject the compressor attribute.
        // Simpler: rebuild a tiny stub that triggers the path.
        let _ = &mut bytes;
        let stub = b"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" compressor=\"vtkZLibDataCompressor\">\n<AppendedData encoding=\"raw\">_</AppendedData></VTKFile>\n";
        let err = parse_appended_raw(stub).unwrap_err();
        assert!(matches!(err, ParseError::Unsupported(_)), "got {err:?}");
    }

    #[test]
    fn parse_appended_raw_rejects_base64_encoding() {
        let bytes =
            b"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\">\n<AppendedData encoding=\"base64\">_garbage</AppendedData>\n</VTKFile>\n";
        let err = parse_appended_raw(bytes).unwrap_err();
        assert!(matches!(err, ParseError::Unsupported(_)), "got {err:?}");
    }

    #[test]
    fn parse_appended_raw_rejects_missing_underscore_separator() {
        // <AppendedData> tag present but no `_` after it.
        let bytes = b"<VTKFile type=\"UnstructuredGrid\"><UnstructuredGrid><Piece NumberOfPoints=\"0\" NumberOfCells=\"0\"></Piece></UnstructuredGrid><AppendedData encoding=\"raw\">  no underscore here  </AppendedData></VTKFile>";
        // We strip the `_` so the parser must hit the missing-separator branch.
        // To make this real: include `_` would let it parse and return a different error.
        // The fixture above genuinely lacks `_`, so the parser's
        // (after_open..bytes.len()).find(|i| bytes[i] == b'_') returns None.
        // BUT — `_` appears in nothing else here, and the chained byte
        // search covers the whole tail; verify the error variant.
        let err = parse_appended_raw(bytes).unwrap_err();
        // Either "missing `_` separator" or another structural error
        // — both are Unsupported / MissingElement.
        assert!(
            matches!(
                err,
                ParseError::Unsupported(_) | ParseError::MissingElement(_)
            ),
            "got {err:?}"
        );
    }

    #[test]
    fn parse_appended_raw_yields_canonical_mesh_with_one_block() {
        // Round-trip: parse + canonicalise yields a Tet4 block.
        let bytes = synth_appended_tet();
        let data = parse_appended_raw(&bytes).expect("parse");
        let (mesh, fields) = data.to_canonical("dispatch-appended");
        assert_eq!(mesh.id, "dispatch-appended");
        assert_eq!(mesh.element_blocks.len(), 1);
        assert_eq!(
            mesh.element_blocks[0].element_type,
            valenx_mesh::ElementType::Tet4
        );
        assert_eq!(fields.len(), 1);
        assert_eq!(fields[0].name, "T");
    }

    /// Round-11 RED→GREEN — R11-4. Pre-fix `parse_ascii` /
    /// `parse_appended_raw` accepted `NumberOfPoints="..."` verbatim
    /// from the `<Piece>` element. A malicious VTU could declare
    /// `NumberOfPoints="18446744073709551615"` (u64::MAX), then the
    /// downstream `Vec::with_capacity(n_points)` / `n_points * 3`
    /// arithmetic would ask the allocator for hundreds of exabytes
    /// before the OOM killer fired. The cap rejects at parse time.
    #[test]
    fn parse_ascii_rejects_number_of_points_above_cap() {
        let hostile = r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="18446744073709551615" NumberOfCells="1">
      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
          0 0 0
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">0</DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">1</DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">10</DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>"#;
        let err = parse_ascii(hostile).expect_err("oversized NumberOfPoints must reject");
        match err {
            ParseError::BadAttribute(msg) => {
                assert!(
                    msg.contains("NumberOfPoints"),
                    "error must name NumberOfPoints; got: {msg}"
                );
                assert!(
                    msg.contains("MAX_VTU_POINTS"),
                    "error must mention the cap; got: {msg}"
                );
            }
            other => panic!("expected ParseError::BadAttribute, got {other:?}"),
        }
    }

    /// Round-11 RED→GREEN — R11-4 sister: NumberOfCells gets the
    /// same cap. The `n_cells` value flows into multiple
    /// `Vec::with_capacity` sites downstream (offsets / types / cell
    /// fields), so the same DoS class applies.
    #[test]
    fn parse_ascii_rejects_number_of_cells_above_cap() {
        let hostile = r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="1" NumberOfCells="18446744073709551615">
      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
          0 0 0
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">0</DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">1</DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">10</DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>"#;
        let err = parse_ascii(hostile).expect_err("oversized NumberOfCells must reject");
        match err {
            ParseError::BadAttribute(msg) => {
                assert!(
                    msg.contains("NumberOfCells"),
                    "error must name NumberOfCells; got: {msg}"
                );
            }
            other => panic!("expected ParseError::BadAttribute, got {other:?}"),
        }
    }

    /// Round-11 RED→GREEN — R11-3. Pre-fix `read_appended_block` did
    /// raw `arr.offset + 4` and `block_start + size` arithmetic. A
    /// malicious appended VTU declaring `offset="18446744073709551612"`
    /// would wrap around inside `arr.offset + 4` and produce a tiny
    /// positive number that *did* index into `raw.len()`, then read
    /// arbitrary bytes outside the intended bounds. The
    /// `checked_add` rejects the overflow at parse time.
    ///
    /// Note: we can't directly synthesise a VTU with such an offset
    /// because the `<DataArray offset="...">` value would never
    /// survive `parse_attribute_u64` followed by `as usize` on a
    /// 32-bit platform — so we test the internal helper directly.
    #[test]
    fn read_appended_block_rejects_offset_overflow() {
        // Build the smallest possible "raw" buffer and a fake
        // DataArray with offset = usize::MAX - 2. The pre-fix code
        // computed `usize::MAX - 2 + 4 = usize::MAX + 2 = 1`
        // (wrapping), then sliced `raw[usize::MAX-2..1]` which
        // panics or returns wrong data. The fix returns
        // BadAttribute.
        let raw = vec![0u8; 64];
        let arr = AppendedDataArrayAttrs {
            name: "evil".to_string(),
            components: 1,
            type_name: "Float32".to_string(),
            offset: usize::MAX - 2,
        };
        let err = read_appended_block(&raw, &arr, "evil")
            .expect_err("offset overflow must be caught by checked_add");
        match err {
            ParseError::BadAttribute(msg) => {
                assert!(
                    msg.contains("overflow"),
                    "error must mention overflow; got: {msg}"
                );
            }
            other => panic!("expected ParseError::BadAttribute, got {other:?}"),
        }
    }
}
