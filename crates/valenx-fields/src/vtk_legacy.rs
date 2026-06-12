//! Minimal VTK legacy **binary** file format parser.
//!
//! VTK has two file formats: the modern XML wrapper (`.vtu`, `.vtp`,
//! etc. — see [`crate::vtu`]) and the older legacy format that ParaView
//! still writes by default with the "VTK legacy binary" exporter. The
//! legacy format mixes ASCII headers with raw big-endian binary data
//! sections.
//!
//! ## Scope
//!
//! Real-world coverage: UNSTRUCTURED_GRID datasets with POINTS,
//! CELLS, CELL_TYPES, plus POINT_DATA / CELL_DATA SCALARS and
//! VECTORS in float / double / int. That covers the slice OpenFOAM,
//! CalculiX, Elmer, and ParaView legacy-export emit when run with
//! their default flags.
//!
//! Out of scope: STRUCTURED_GRID, STRUCTURED_POINTS, RECTILINEAR_GRID,
//! POLYDATA, FIELD-data records, COLOR_SCALARS / TENSORS / NORMALS —
//! callers get [`ParseError::Unsupported`] with the missing keyword.
//!
//! ## On-disk layout
//!
//! ```text
//! # vtk DataFile Version <X>.<Y>      (line 1: magic header)
//! <free-form title>                    (line 2: max 256 chars)
//! BINARY                               (line 3: format keyword)
//! DATASET UNSTRUCTURED_GRID            (line 4: dataset type)
//! POINTS <n> <type>                    (geometry header)
//! <binary block: n * 3 * sizeof(type) bytes, BIG-ENDIAN>
//! CELLS <n_cells> <total_size>
//! <binary block: total_size * 4 bytes (i32 BE), per-cell prefix is
//!  the count followed by the connectivity indices>
//! CELL_TYPES <n_cells>
//! <binary block: n_cells * 4 bytes (i32 BE)>
//! POINT_DATA <n>                       (optional)
//! SCALARS <name> <type> <comps>        (data header)
//! LOOKUP_TABLE <name>                  (data header pre-amble)
//! <binary block: n * comps * sizeof(type) bytes, BIG-ENDIAN>
//! ```
//!
//! ## Endianness
//!
//! VTK legacy binary stores **everything big-endian** regardless of
//! host. Modern x86 hosts byte-swap on read. The parser uses
//! `u32::from_be_bytes` / `f32::from_be_bytes` / `f64::from_be_bytes`
//! everywhere — never `from_ne_bytes`.

use thiserror::Error;

/// One scalar / vector array attached to either points or cells.
#[derive(Clone, Debug, PartialEq)]
pub struct LegacyArray {
    pub name: String,
    /// Number of components per tuple: 1 for SCALARS, 3 for VECTORS.
    pub components: usize,
    /// Number of tuples = `points.len()` for POINT_DATA or
    /// `cells.len()` for CELL_DATA.
    pub tuples: usize,
    /// Flat row-major sample buffer; length = `tuples * components`.
    pub data: Vec<f64>,
}

impl LegacyArray {
    /// Sample count (handy for catalog-side validation).
    pub fn samples(&self) -> usize {
        self.tuples
    }
}

/// One unstructured grid + its attached arrays. Mirrors the VTU
/// extraction shape so the colour-mapping / mesh-conversion paths can
/// share code with the XML parser.
#[derive(Clone, Debug, Default)]
pub struct LegacyData {
    /// Free-form title from line 2 of the file.
    pub title: String,
    /// Point coordinates as `[x, y, z]` triples.
    pub points: Vec<[f64; 3]>,
    /// Per-cell connectivity buffer. Each cell is stored as
    /// `[n_points_in_cell, p0, p1, ..., p(n-1)]` flattened — the same
    /// shape VTK writes on disk. Use [`Self::iter_cells`] to walk it.
    pub cells: Vec<u32>,
    /// VTK cell-type code per cell (see `crate::vtu::VtuCellType` for
    /// the codebook).
    pub cell_types: Vec<u8>,
    /// Per-point arrays.
    pub point_data: Vec<LegacyArray>,
    /// Per-cell arrays.
    pub cell_data: Vec<LegacyArray>,
}

impl LegacyData {
    /// Iterate cells, yielding `(cell_type_code, &[point_indices])`.
    pub fn iter_cells(&self) -> CellIter<'_> {
        CellIter {
            data: &self.cells,
            types: &self.cell_types,
            cursor: 0,
            cell_idx: 0,
        }
    }
}

/// Iterator over `(cell_type, connectivity_slice)` tuples returned by
/// [`LegacyData::iter_cells`].
pub struct CellIter<'a> {
    data: &'a [u32],
    types: &'a [u8],
    cursor: usize,
    cell_idx: usize,
}

impl<'a> Iterator for CellIter<'a> {
    type Item = (u8, &'a [u32]);
    fn next(&mut self) -> Option<Self::Item> {
        if self.cursor >= self.data.len() {
            return None;
        }
        let n = self.data[self.cursor] as usize;
        let start = self.cursor + 1;
        // Round-5 fix: a malformed CELLS block with `n = u32::MAX`
        // (the count prefix at `data[cursor]`) used to overflow
        // `start + n` and / or slice past the end of `data`, panicking
        // with a slice OOB. We now bounds-check the cell-end against
        // `data.len()` and short-circuit cleanly when the implied end
        // would escape the buffer. Same applies to `cell_idx`
        // exceeding `types.len()` — a malformed legacy file can have
        // a CELLS block declare more cells than CELL_TYPES carries.
        let end = match start.checked_add(n) {
            Some(e) if e <= self.data.len() => e,
            _ => {
                // Stop iteration rather than panic. Callers iterating
                // a partially-broken file get only the cells that were
                // safe to decode.
                self.cursor = self.data.len();
                return None;
            }
        };
        if self.cell_idx >= self.types.len() {
            self.cursor = self.data.len();
            return None;
        }
        let cell = &self.data[start..end];
        let ct = self.types[self.cell_idx];
        self.cursor = end;
        self.cell_idx += 1;
        Some((ct, cell))
    }
}

impl LegacyData {
    /// Convert this `LegacyData` into canonical valenx types: a
    /// [`valenx_mesh::Mesh`] (one [`valenx_mesh::ElementBlock`] per
    /// distinct cell type) plus a list of [`crate::Field`]s for the
    /// point + cell data. Mirrors
    /// [`crate::vtu::VtuData::to_canonical`] so the downstream catalog
    /// path is the same regardless of which VTK flavour the file
    /// landed in.
    ///
    /// Cells with no canonical mapping (`Vertex`, `Unknown`) are
    /// silently skipped — same policy as the VTU bridge.
    pub fn to_canonical(
        &self,
        mesh_id: impl Into<String>,
    ) -> (valenx_mesh::Mesh, Vec<crate::Field>) {
        use crate::vtu::{vtu_to_element_type, VtuCellType};
        use std::collections::BTreeMap;
        use valenx_mesh::{ElementBlock, Mesh};

        let mut mesh = Mesh::new(mesh_id);
        mesh.nodes = self
            .points
            .iter()
            .map(|p| nalgebra::Vector3::new(p[0], p[1], p[2]))
            .collect();

        let mut blocks: BTreeMap<u8, ElementBlock> = BTreeMap::new();
        for (ct_code, nodes) in self.iter_cells() {
            let ct = VtuCellType::from_code(ct_code);
            let Some(canonical) = vtu_to_element_type(ct) else {
                continue;
            };
            let entry = blocks
                .entry(canonical as u8)
                .or_insert_with(|| ElementBlock::new(canonical));
            entry.connectivity.extend_from_slice(nodes);
        }
        mesh.element_blocks = blocks.into_values().collect();
        // Populate counts AND quality scalars (see vtu.rs for rationale).
        let _ = mesh.recompute_quality_stats();

        let mut fields: Vec<crate::Field> = Vec::new();
        for arr in &self.point_data {
            fields.push(legacy_to_field(arr, crate::Location::OnNode));
        }
        for arr in &self.cell_data {
            fields.push(legacy_to_field(arr, crate::Location::OnCell));
        }
        (mesh, fields)
    }
}

/// Convert one [`LegacyArray`] into a canonical [`crate::Field`].
/// Same component-count -> FieldKind mapping as the VTU bridge:
/// 1 = Scalar, 3 = Vector, 9 = Tensor 3x3, anything else falls back
/// to a flat n-component vector so the data isn't silently lost.
fn legacy_to_field(arr: &LegacyArray, location: crate::Location) -> crate::Field {
    let kind = match arr.components {
        1 => crate::FieldKind::Scalar,
        3 => crate::FieldKind::Vector { dim: 3 },
        9 => crate::FieldKind::Tensor { rows: 3, cols: 3 },
        n => crate::FieldKind::Vector { dim: n as u8 },
    };
    let range = legacy_field_range(&arr.data);
    crate::Field {
        name: arr.name.clone(),
        kind,
        location,
        // Legacy VTK has no concept of named regions — the whole
        // dataset is one implicit region (matches the VTU bridge).
        region: crate::RegionRef("default".to_string()),
        // Legacy VTK headers don't carry units. Mark dimensionless;
        // adapters that know the physics can re-stamp via a follow-up
        // pass before serving them to the report layer.
        units: crate::units::DIMENSIONLESS,
        // Steady by default — single-timestep per file. PVD-style
        // multi-step series re-stamp these to TimeKey::Time at load.
        time: crate::TimeKey::Steady,
        data: arr.data.clone(),
        range,
    }
}

fn legacy_field_range(data: &[f64]) -> Option<(f64, f64)> {
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

/// Upper bound on the per-section element / point / tuple count we'll
/// honour from a VTK legacy header before refusing the file outright.
///
/// Round-14 H2 (sister to the round-11 `MAX_VTU_POINTS` cap on
/// `valenx-fields::vtu`): pre-fix a hostile or accidentally-truncated
/// header like `POINTS 10000000000 float` would walk straight into
/// `bytes_needed = n * 3 * elem_size`, multiply usize unchecked, and
/// either overflow into a small number (so the cursor's truncation
/// check passes despite the file being a few hundred bytes long) or
/// allocate hundreds of GB before the OOM-killer notices. 256 M
/// matches the VTU cap so the two parsers refuse the same class of
/// pathological input.
pub const MAX_VTK_LEGACY_POINTS: usize = 256_000_000;

/// Errors raised by [`parse_binary`].
#[derive(Debug, Error)]
pub enum ParseError {
    /// First line is not `# vtk DataFile Version …`.
    #[error("not a VTK legacy file: missing `# vtk DataFile Version` header")]
    MissingMagic,
    /// Line 3 declared a non-`BINARY` format (only binary is supported).
    #[error("expected BINARY format declaration on line 3, got `{0}`")]
    NotBinary(String),
    /// A recognised section uses a feature the parser doesn't handle.
    #[error("unsupported VTK legacy feature: {0}")]
    Unsupported(String),
    /// A section header couldn't be parsed (missing count, bad keyword).
    #[error("malformed `{section}` header at byte {byte}: {reason}")]
    MalformedHeader {
        /// Name of the offending section (`POINTS`, `CELLS`, …).
        section: String,
        /// Approximate byte offset of the bad header.
        byte: usize,
        /// Short human-readable explanation.
        reason: String,
    },
    /// A binary data block was shorter than the declared length.
    #[error("truncated binary data section `{section}`: needed {needed} bytes, had {available}")]
    Truncated {
        /// Name of the offending section.
        section: String,
        /// Number of bytes the header declared.
        needed: usize,
        /// Bytes still available in the buffer.
        available: usize,
    },
    /// A `SCALARS` line declared a type the parser doesn't recognise.
    #[error("unsupported scalar type `{0}`")]
    UnsupportedScalarType(String),
    /// A section header declared a count larger than
    /// [`MAX_VTK_LEGACY_POINTS`] (or the count × element-size product
    /// overflowed usize). Refusing this class of input keeps a hostile
    /// or accidentally-truncated header from allocating hundreds of GB.
    ///
    /// Round-14 H2.
    #[error(
        "VTK legacy section `{what}` declared {count} elements (max {max}) — refusing to allocate"
    )]
    TooLarge {
        /// Which section header tripped the cap (`POINTS`, `CELLS`,
        /// `CELL_TYPES`, etc.).
        what: String,
        /// Element / point / tuple count declared in the header.
        count: usize,
        /// Cap that was exceeded.
        max: usize,
    },
}

/// Parse a VTK legacy binary file into [`LegacyData`].
///
/// The byte buffer is the entire file contents (legacy VTK files mix
/// ASCII headers with binary payloads, so we can't stream-parse with
/// `BufRead::lines`). Call sites typically `std::fs::read(path)?` and
/// pass the resulting `Vec<u8>`.
pub fn parse_binary(bytes: &[u8]) -> Result<LegacyData, ParseError> {
    let mut cur = Cursor::new(bytes);
    // Line 1: magic header.
    let line1 = cur.read_ascii_line()?;
    if !line1.starts_with("# vtk DataFile Version") {
        return Err(ParseError::MissingMagic);
    }
    // Line 2: title.
    let title = cur.read_ascii_line()?;
    // Line 3: format declaration.
    let format = cur.read_ascii_line()?;
    if format.trim() != "BINARY" {
        return Err(ParseError::NotBinary(format));
    }
    // Line 4: dataset type.
    let dataset = cur.read_ascii_line()?;
    let trimmed = dataset.trim();
    if trimmed != "DATASET UNSTRUCTURED_GRID" {
        return Err(ParseError::Unsupported(format!(
            "DATASET other than UNSTRUCTURED_GRID: `{trimmed}`"
        )));
    }

    let mut data = LegacyData {
        title,
        ..Default::default()
    };

    // Walk subsequent section headers until EOF.
    while let Some(line) = cur.read_ascii_line_or_none()? {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix("POINTS ") {
            // POINTS <n> <type>
            let mut it = rest.split_ascii_whitespace();
            let n: usize = parse_size("POINTS", cur.pos(), it.next())?;
            let dtype = it.next().ok_or_else(|| ParseError::MalformedHeader {
                section: "POINTS".into(),
                byte: cur.pos(),
                reason: "missing data type".into(),
            })?;
            data.points = read_points_block(&mut cur, n, dtype)?;
        } else if let Some(rest) = line.strip_prefix("CELLS ") {
            // CELLS <n_cells> <total_size>
            let mut it = rest.split_ascii_whitespace();
            let _n_cells: usize = parse_size("CELLS", cur.pos(), it.next())?;
            let total_size: usize = parse_size("CELLS", cur.pos(), it.next())?;
            // Round-14 H2: same cap shape as the POINTS path — a
            // hostile `CELLS 0 999999999999` would otherwise overflow
            // `total_size * 4` into a small number and pass the
            // truncation check.
            if total_size > MAX_VTK_LEGACY_POINTS {
                return Err(ParseError::TooLarge {
                    what: "CELLS".into(),
                    count: total_size,
                    max: MAX_VTK_LEGACY_POINTS,
                });
            }
            let bytes_needed = total_size.checked_mul(4).ok_or(ParseError::TooLarge {
                what: "CELLS".into(),
                count: total_size,
                max: MAX_VTK_LEGACY_POINTS,
            })?;
            let raw = cur.read_exact_bytes("CELLS", bytes_needed)?;
            data.cells = decode_be_u32(raw);
            cur.skip_trailing_newlines();
        } else if let Some(rest) = line.strip_prefix("CELL_TYPES ") {
            let n: usize = parse_size(
                "CELL_TYPES",
                cur.pos(),
                rest.split_ascii_whitespace().next(),
            )?;
            // Round-14 H2: refuse oversized counts. Same shape as the
            // POINTS / CELLS paths.
            if n > MAX_VTK_LEGACY_POINTS {
                return Err(ParseError::TooLarge {
                    what: "CELL_TYPES".into(),
                    count: n,
                    max: MAX_VTK_LEGACY_POINTS,
                });
            }
            let bytes_needed = n.checked_mul(4).ok_or(ParseError::TooLarge {
                what: "CELL_TYPES".into(),
                count: n,
                max: MAX_VTK_LEGACY_POINTS,
            })?;
            let raw = cur.read_exact_bytes("CELL_TYPES", bytes_needed)?;
            data.cell_types = decode_be_u32(raw).into_iter().map(|v| v as u8).collect();
            cur.skip_trailing_newlines();
        } else if let Some(rest) = line.strip_prefix("POINT_DATA ") {
            let n: usize = parse_size(
                "POINT_DATA",
                cur.pos(),
                rest.split_ascii_whitespace().next(),
            )?;
            // Round-15 L1: extend the round-14 H2 cap to the
            // POINT_DATA tuples count. Pre-fix the only check was via
            // `tuples * components` inside `read_typed_block`, which
            // would overflow usize for `n ≥ usize::MAX / 3` and slip
            // a few-GB allocation past the bytes-needed check.
            if n > MAX_VTK_LEGACY_POINTS {
                return Err(ParseError::TooLarge {
                    what: "POINT_DATA".into(),
                    count: n,
                    max: MAX_VTK_LEGACY_POINTS,
                });
            }
            data.point_data = read_data_arrays(&mut cur, n)?;
        } else if let Some(rest) = line.strip_prefix("CELL_DATA ") {
            let n: usize =
                parse_size("CELL_DATA", cur.pos(), rest.split_ascii_whitespace().next())?;
            // Round-15 L1 sister: same cap on the CELL_DATA path.
            if n > MAX_VTK_LEGACY_POINTS {
                return Err(ParseError::TooLarge {
                    what: "CELL_DATA".into(),
                    count: n,
                    max: MAX_VTK_LEGACY_POINTS,
                });
            }
            data.cell_data = read_data_arrays(&mut cur, n)?;
        } else if line.starts_with("FIELD ")
            || line.starts_with("COLOR_SCALARS ")
            || line.starts_with("TENSORS ")
            || line.starts_with("NORMALS ")
            || line.starts_with("TEXTURE_COORDINATES ")
        {
            return Err(ParseError::Unsupported(format!(
                "VTK legacy section `{line}` (out of v0 scope)"
            )));
        } else {
            // Unknown line — bail rather than silently skipping; that
            // way we don't miss data the caller cares about.
            return Err(ParseError::Unsupported(format!(
                "unrecognised section header `{line}`"
            )));
        }
    }
    Ok(data)
}

/// Internal: read the per-point coordinate block. VTK's POINTS section
/// stores `n` triples of `dtype`. We support `float` (f32) and `double`
/// (f64), promoting to f64 for the canonical [`LegacyData`].
fn read_points_block(cur: &mut Cursor, n: usize, dtype: &str) -> Result<Vec<[f64; 3]>, ParseError> {
    let elem_size = match dtype {
        "float" => 4,
        "double" => 8,
        other => {
            return Err(ParseError::UnsupportedScalarType(other.to_string()));
        }
    };
    // Round-14 H2: refuse a header that declares more points than the
    // hard cap, AND refuse a count whose `n * 3 * elem_size` would
    // overflow usize. Pre-fix `n * 3 * elem_size` for `n = u64::MAX`
    // would wrap to a small number, the truncation check would pass
    // against a normal-sized file, and the decoder would then walk
    // off the end allocating up to hundreds of GB.
    if n > MAX_VTK_LEGACY_POINTS {
        return Err(ParseError::TooLarge {
            what: "POINTS".into(),
            count: n,
            max: MAX_VTK_LEGACY_POINTS,
        });
    }
    let bytes_needed = n
        .checked_mul(3)
        .and_then(|m| m.checked_mul(elem_size))
        .ok_or(ParseError::TooLarge {
            what: "POINTS".into(),
            count: n,
            max: MAX_VTK_LEGACY_POINTS,
        })?;
    let raw = cur.read_exact_bytes("POINTS", bytes_needed)?;
    let mut out: Vec<[f64; 3]> = Vec::with_capacity(n);
    for i in 0..n {
        let off = i * 3 * elem_size;
        let mut p = [0f64; 3];
        for (c, slot) in p.iter_mut().enumerate() {
            let s = off + c * elem_size;
            let e = s + elem_size;
            *slot = match dtype {
                "float" => f32::from_be_bytes(raw[s..e].try_into().unwrap()) as f64,
                "double" => f64::from_be_bytes(raw[s..e].try_into().unwrap()),
                _ => unreachable!(),
            };
        }
        out.push(p);
    }
    cur.skip_trailing_newlines();
    Ok(out)
}

/// Internal: read a sequence of SCALARS / VECTORS arrays under a
/// POINT_DATA / CELL_DATA section. Stops on the next top-level
/// section keyword or EOF.
fn read_data_arrays(cur: &mut Cursor, tuples: usize) -> Result<Vec<LegacyArray>, ParseError> {
    let mut out: Vec<LegacyArray> = Vec::new();
    while let Some(line) = cur.peek_ascii_line()? {
        let trimmed = line.trim();
        if let Some(rest) = trimmed.strip_prefix("SCALARS ") {
            cur.consume_peeked_line(); // consume the SCALARS header
                                       // SCALARS <name> <type> [<components>=1]
            let mut it = rest.split_ascii_whitespace();
            let name = it
                .next()
                .ok_or_else(|| ParseError::MalformedHeader {
                    section: "SCALARS".into(),
                    byte: cur.pos(),
                    reason: "missing name".into(),
                })?
                .to_string();
            let dtype = it.next().ok_or_else(|| ParseError::MalformedHeader {
                section: "SCALARS".into(),
                byte: cur.pos(),
                reason: "missing type".into(),
            })?;
            let components: usize = it.next().map(|s| s.parse().unwrap_or(1)).unwrap_or(1);
            // SCALARS is followed by a LOOKUP_TABLE line; consume it.
            let lookup = cur.read_ascii_line()?;
            if !lookup.trim_start().starts_with("LOOKUP_TABLE") {
                return Err(ParseError::MalformedHeader {
                    section: "SCALARS".into(),
                    byte: cur.pos(),
                    reason: format!(
                        "expected LOOKUP_TABLE after SCALARS header, got `{}`",
                        lookup.trim()
                    ),
                });
            }
            // Round-15 L1: `tuples * components` uses checked_mul so a
            // hostile or accidentally-broken file with `tuples` near
            // `usize::MAX / components` can't overflow and slip a
            // small bytes_needed past the cap inside read_typed_block.
            // Pre-fix `tuples = 2^63` with `components = 2` would
            // wrap to 0, allocate nothing, and the caller would see
            // a fictitious empty array.
            let n_elements = tuples.checked_mul(components).ok_or(ParseError::TooLarge {
                what: "SCALARS".into(),
                count: tuples,
                max: MAX_VTK_LEGACY_POINTS,
            })?;
            let data = read_typed_block(cur, &name, n_elements, dtype)?;
            out.push(LegacyArray {
                name,
                components,
                tuples,
                data,
            });
        } else if let Some(rest) = trimmed.strip_prefix("VECTORS ") {
            cur.consume_peeked_line();
            // VECTORS <name> <type>
            let mut it = rest.split_ascii_whitespace();
            let name = it
                .next()
                .ok_or_else(|| ParseError::MalformedHeader {
                    section: "VECTORS".into(),
                    byte: cur.pos(),
                    reason: "missing name".into(),
                })?
                .to_string();
            let dtype = it.next().ok_or_else(|| ParseError::MalformedHeader {
                section: "VECTORS".into(),
                byte: cur.pos(),
                reason: "missing type".into(),
            })?;
            // Round-15 L1: same checked_mul shield as the SCALARS path.
            // For VECTORS the components is fixed at 3.
            let n_elements = tuples.checked_mul(3).ok_or(ParseError::TooLarge {
                what: "VECTORS".into(),
                count: tuples,
                max: MAX_VTK_LEGACY_POINTS,
            })?;
            let data = read_typed_block(cur, &name, n_elements, dtype)?;
            out.push(LegacyArray {
                name,
                components: 3,
                tuples,
                data,
            });
        } else {
            // Not a data-array header — leave the line in the buffer
            // for the outer loop to interpret.
            break;
        }
    }
    Ok(out)
}

fn read_typed_block(
    cur: &mut Cursor,
    section: &str,
    n_elements: usize,
    dtype: &str,
) -> Result<Vec<f64>, ParseError> {
    let elem_size: usize = match dtype {
        "float" => 4,
        "double" => 8,
        "int" => 4,
        "unsigned_int" => 4,
        "short" => 2,
        "unsigned_short" => 2,
        "char" => 1,
        "unsigned_char" => 1,
        other => {
            return Err(ParseError::UnsupportedScalarType(other.to_string()));
        }
    };
    // Round-14 H2: data arrays nested under POINT_DATA / CELL_DATA go
    // through this path; same overflow / cap shape as the POINTS /
    // CELLS sections so a hostile SCALARS / VECTORS header can't
    // overflow `n_elements * elem_size` and slip a few-GB allocation
    // past the truncation check.
    if n_elements > MAX_VTK_LEGACY_POINTS {
        return Err(ParseError::TooLarge {
            what: section.to_string(),
            count: n_elements,
            max: MAX_VTK_LEGACY_POINTS,
        });
    }
    let bytes_needed = n_elements
        .checked_mul(elem_size)
        .ok_or(ParseError::TooLarge {
            what: section.to_string(),
            count: n_elements,
            max: MAX_VTK_LEGACY_POINTS,
        })?;
    let raw = cur.read_exact_bytes(section, bytes_needed)?;
    let mut out: Vec<f64> = Vec::with_capacity(n_elements);
    for i in 0..n_elements {
        let s = i * elem_size;
        let e = s + elem_size;
        let v: f64 = match dtype {
            "float" => f32::from_be_bytes(raw[s..e].try_into().unwrap()) as f64,
            "double" => f64::from_be_bytes(raw[s..e].try_into().unwrap()),
            "int" => i32::from_be_bytes(raw[s..e].try_into().unwrap()) as f64,
            "unsigned_int" => u32::from_be_bytes(raw[s..e].try_into().unwrap()) as f64,
            "short" => i16::from_be_bytes(raw[s..e].try_into().unwrap()) as f64,
            "unsigned_short" => u16::from_be_bytes(raw[s..e].try_into().unwrap()) as f64,
            "char" => raw[s] as i8 as f64,
            "unsigned_char" => raw[s] as f64,
            _ => unreachable!(),
        };
        out.push(v);
    }
    cur.skip_trailing_newlines();
    Ok(out)
}

fn decode_be_u32(raw: &[u8]) -> Vec<u32> {
    let n = raw.len() / 4;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let s = i * 4;
        out.push(u32::from_be_bytes(raw[s..s + 4].try_into().unwrap()));
    }
    out
}

fn parse_size(section: &str, byte: usize, s: Option<&str>) -> Result<usize, ParseError> {
    let s = s.ok_or_else(|| ParseError::MalformedHeader {
        section: section.into(),
        byte,
        reason: "missing count".into(),
    })?;
    s.parse().map_err(|_| ParseError::MalformedHeader {
        section: section.into(),
        byte,
        reason: format!("count `{s}` is not a non-negative integer"),
    })
}

/// A tiny cursor that knows how to read either ASCII text lines or
/// raw binary blocks from the same backing buffer. Pulled out of the
/// parser so the seek bookkeeping is centralised.
struct Cursor<'a> {
    bytes: &'a [u8],
    pos: usize,
    /// One-line look-ahead buffer for the data-array reader, which
    /// needs to peek at the next section header without consuming it.
    peeked: Option<(usize, String)>,
}

impl<'a> Cursor<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            pos: 0,
            peeked: None,
        }
    }

    fn pos(&self) -> usize {
        self.pos
    }

    /// Read up to and including the next `\n`. Returns the line
    /// contents (without the trailing newline). Errors at EOF.
    fn read_ascii_line(&mut self) -> Result<String, ParseError> {
        if let Some((end, line)) = self.peeked.take() {
            self.pos = end;
            return Ok(line);
        }
        match self.read_ascii_line_or_none()? {
            Some(line) => Ok(line),
            None => Err(ParseError::MalformedHeader {
                section: "<eof>".into(),
                byte: self.pos,
                reason: "unexpected EOF reading ASCII line".into(),
            }),
        }
    }

    fn read_ascii_line_or_none(&mut self) -> Result<Option<String>, ParseError> {
        if let Some((end, line)) = self.peeked.take() {
            self.pos = end;
            return Ok(Some(line));
        }
        if self.pos >= self.bytes.len() {
            return Ok(None);
        }
        let start = self.pos;
        let mut end = start;
        while end < self.bytes.len() && self.bytes[end] != b'\n' {
            end += 1;
        }
        let line_bytes = &self.bytes[start..end];
        // Strip any trailing CR for CRLF-encoded files.
        let trimmed = if line_bytes.ends_with(b"\r") {
            &line_bytes[..line_bytes.len() - 1]
        } else {
            line_bytes
        };
        let line = std::str::from_utf8(trimmed)
            .map_err(|e| ParseError::MalformedHeader {
                section: "<header-line>".into(),
                byte: start,
                reason: format!("non-UTF8 ASCII header: {e}"),
            })?
            .to_string();
        self.pos = if end < self.bytes.len() { end + 1 } else { end };
        Ok(Some(line))
    }

    fn peek_ascii_line(&mut self) -> Result<Option<String>, ParseError> {
        if let Some((_, ref line)) = self.peeked {
            return Ok(Some(line.clone()));
        }
        let saved_pos = self.pos;
        let line = self.read_ascii_line_or_none()?;
        let new_pos = self.pos;
        if let Some(ref l) = line {
            self.peeked = Some((new_pos, l.clone()));
            self.pos = saved_pos;
        }
        Ok(line)
    }

    fn consume_peeked_line(&mut self) {
        if let Some((end, _)) = self.peeked.take() {
            self.pos = end;
        }
    }

    fn read_exact_bytes(&mut self, section: &str, n: usize) -> Result<&'a [u8], ParseError> {
        if self.pos + n > self.bytes.len() {
            return Err(ParseError::Truncated {
                section: section.into(),
                needed: n,
                available: self.bytes.len() - self.pos,
            });
        }
        let slice = &self.bytes[self.pos..self.pos + n];
        self.pos += n;
        Ok(slice)
    }

    /// VTK legacy writes a `\n` between binary sections; consume any
    /// trailing newlines so the next ASCII header read picks up the
    /// right line.
    fn skip_trailing_newlines(&mut self) {
        while self.pos < self.bytes.len()
            && (self.bytes[self.pos] == b'\n' || self.bytes[self.pos] == b'\r')
        {
            self.pos += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build the smallest possible legal VTK legacy binary file:
    /// 4 points forming a single tetrahedron, with one POINT_DATA
    /// scalar.
    fn synthesize_tet_with_scalar() -> Vec<u8> {
        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        buf.extend_from_slice(b"valenx test tet\n");
        buf.extend_from_slice(b"BINARY\n");
        buf.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        buf.extend_from_slice(b"POINTS 4 float\n");
        // 4 points * 3 floats = 12 floats, big-endian.
        let coords: [f32; 12] = [
            0.0, 0.0, 0.0, // p0
            1.0, 0.0, 0.0, // p1
            0.0, 1.0, 0.0, // p2
            0.0, 0.0, 1.0, // p3
        ];
        for v in coords {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        buf.push(b'\n');
        // CELLS 1 5: total_size = 1 (count) + 4 (point indices) = 5
        buf.extend_from_slice(b"CELLS 1 5\n");
        for v in [4u32, 0, 1, 2, 3] {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        buf.push(b'\n');
        buf.extend_from_slice(b"CELL_TYPES 1\n");
        buf.extend_from_slice(&10u32.to_be_bytes()); // VTK_TETRA = 10
        buf.push(b'\n');
        // POINT_DATA: one scalar named "T", one component, 4 tuples.
        buf.extend_from_slice(b"POINT_DATA 4\n");
        buf.extend_from_slice(b"SCALARS T float 1\n");
        buf.extend_from_slice(b"LOOKUP_TABLE default\n");
        for v in [273.15_f32, 280.0, 290.0, 300.0] {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        buf.push(b'\n');
        buf
    }

    /// Deterministic xorshift64 PRNG — reproducible, no external rng dep.
    fn xorshift64(state: &mut u64) -> u64 {
        let mut x = *state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        *state = x;
        x
    }

    /// Robustness: `parse_binary` must NEVER panic on arbitrary input — it
    /// may only return `Ok` or a typed `Err`. This feeds truncations,
    /// single-byte corruptions, and deterministic pseudo-random buffers; if
    /// any input panicked, the test would fail via panic propagation. The
    /// DoS caps + `read_exact_bytes` length-gating keep every iteration
    /// small and fast (counts can't drive a large allocation).
    #[test]
    fn parse_binary_never_panics_on_adversarial_input() {
        let valid = synthesize_tet_with_scalar();

        // 1. Every truncated prefix of a valid file.
        for k in 0..=valid.len() {
            let _ = parse_binary(&valid[..k]);
        }

        // 2. Single-byte corruption (bit-flip and NUL) at every position.
        for i in 0..valid.len() {
            let mut flipped = valid.clone();
            flipped[i] ^= 0xFF;
            let _ = parse_binary(&flipped);
            let mut nulled = valid.clone();
            nulled[i] = 0;
            let _ = parse_binary(&nulled);
        }

        // 3. Deterministic pseudo-random buffers of assorted small sizes.
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15; // fixed seed
        for _ in 0..2000 {
            let len = (xorshift64(&mut state) % 256) as usize;
            let buf: Vec<u8> = (0..len)
                .map(|_| (xorshift64(&mut state) & 0xFF) as u8)
                .collect();
            let _ = parse_binary(&buf);
        }

        // 4. Structured fuzz that genuinely reaches the binary section
        //    decoders: a valid header followed by a POINTS (and sometimes a
        //    POINT_DATA SCALARS) section whose declared count, type token, and
        //    randomized — possibly truncated — binary body are all random. So
        //    read_points_block / read_typed_block / read_exact_bytes get driven
        //    with arbitrary bytes, not bailed at the header. Counts stay tiny
        //    (0..=8) so nothing allocates large. Independent fixed seed.
        let mut state2: u64 = 0x2545_F491_4F6C_DD1D;
        let types = ["float", "double", "int", "short", "char", "wobble"];
        for _ in 0..3000 {
            let mut buf: Vec<u8> =
                b"# vtk DataFile Version 3.0\nt\nBINARY\nDATASET UNSTRUCTURED_GRID\n".to_vec();
            let n = xorshift64(&mut state2) % 9;
            let pt = types[(xorshift64(&mut state2) % types.len() as u64) as usize];
            buf.extend_from_slice(format!("POINTS {n} {pt}\n").as_bytes());
            let body = (xorshift64(&mut state2) % 80) as usize;
            buf.extend((0..body).map(|_| (xorshift64(&mut state2) & 0xFF) as u8));
            if xorshift64(&mut state2) & 1 == 0 {
                let m = xorshift64(&mut state2) % 9;
                let st = types[(xorshift64(&mut state2) % types.len() as u64) as usize];
                buf.extend_from_slice(
                    format!("POINT_DATA {m}\nSCALARS s {st} 1\nLOOKUP_TABLE default\n").as_bytes(),
                );
                let sb = (xorshift64(&mut state2) % 80) as usize;
                buf.extend((0..sb).map(|_| (xorshift64(&mut state2) & 0xFF) as u8));
            }
            let _ = parse_binary(&buf);
        }
    }

    #[test]
    fn parse_binary_handles_a_minimal_tet_with_scalars() {
        let bytes = synthesize_tet_with_scalar();
        let data = parse_binary(&bytes).expect("parse");
        assert_eq!(data.title, "valenx test tet");
        assert_eq!(data.points.len(), 4);
        assert_eq!(data.points[1], [1.0, 0.0, 0.0]);
        assert_eq!(data.points[3], [0.0, 0.0, 1.0]);
        assert_eq!(data.cell_types, vec![10]);
        assert_eq!(data.point_data.len(), 1);
        let arr = &data.point_data[0];
        assert_eq!(arr.name, "T");
        assert_eq!(arr.components, 1);
        assert_eq!(arr.tuples, 4);
        assert!((arr.data[0] - 273.15).abs() < 1e-3);
    }

    #[test]
    fn iter_cells_yields_one_tet() {
        let bytes = synthesize_tet_with_scalar();
        let data = parse_binary(&bytes).expect("parse");
        let cells: Vec<(u8, Vec<u32>)> =
            data.iter_cells().map(|(t, ix)| (t, ix.to_vec())).collect();
        assert_eq!(cells.len(), 1);
        assert_eq!(cells[0].0, 10);
        assert_eq!(cells[0].1, vec![0, 1, 2, 3]);
    }

    #[test]
    fn parse_binary_rejects_ascii_format_files() {
        // Minimal header but format = ASCII, not BINARY.
        let mut bytes: Vec<u8> = Vec::new();
        bytes.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        bytes.extend_from_slice(b"some title\n");
        bytes.extend_from_slice(b"ASCII\n");
        bytes.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        let err = parse_binary(&bytes).expect_err("must reject");
        assert!(matches!(err, ParseError::NotBinary(_)));
    }

    #[test]
    fn parse_binary_rejects_missing_magic_header() {
        let bytes = b"BINARY\nDATASET UNSTRUCTURED_GRID\n";
        let err = parse_binary(bytes).expect_err("must reject");
        assert!(matches!(err, ParseError::MissingMagic));
    }

    #[test]
    fn parse_binary_rejects_unsupported_dataset_types() {
        let mut bytes: Vec<u8> = Vec::new();
        bytes.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        bytes.extend_from_slice(b"\n");
        bytes.extend_from_slice(b"BINARY\n");
        bytes.extend_from_slice(b"DATASET STRUCTURED_GRID\n");
        let err = parse_binary(&bytes).expect_err("must reject");
        assert!(matches!(err, ParseError::Unsupported(_)));
    }

    #[test]
    fn parse_binary_rejects_truncated_data_blocks() {
        // POINTS 100 float promises 100 * 3 * 4 = 1200 bytes of float
        // data; provide just one float so the cursor must report
        // truncation.
        let mut bytes: Vec<u8> = Vec::new();
        bytes.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        bytes.extend_from_slice(b"truncated\n");
        bytes.extend_from_slice(b"BINARY\n");
        bytes.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        bytes.extend_from_slice(b"POINTS 100 float\n");
        bytes.extend_from_slice(&1.0_f32.to_be_bytes()); // only 4 bytes
        let err = parse_binary(&bytes).expect_err("must reject");
        match err {
            ParseError::Truncated {
                section,
                needed,
                available,
            } => {
                assert_eq!(section, "POINTS");
                assert_eq!(needed, 1200);
                assert!(available < needed);
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn parse_binary_handles_double_precision_points() {
        let mut bytes: Vec<u8> = Vec::new();
        bytes.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        bytes.extend_from_slice(b"two-doubles\n");
        bytes.extend_from_slice(b"BINARY\n");
        bytes.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        bytes.extend_from_slice(b"POINTS 1 double\n");
        for v in [1.5_f64, 2.5, 3.5] {
            bytes.extend_from_slice(&v.to_be_bytes());
        }
        bytes.push(b'\n');
        bytes.extend_from_slice(b"CELLS 0 0\n");
        bytes.extend_from_slice(b"CELL_TYPES 0\n");
        let data = parse_binary(&bytes).expect("parse");
        assert_eq!(data.points, vec![[1.5, 2.5, 3.5]]);
    }

    #[test]
    fn to_canonical_one_tet_yields_real_mesh_and_one_field() {
        let bytes = synthesize_tet_with_scalar();
        let data = parse_binary(&bytes).expect("parse");
        let (mesh, fields) = data.to_canonical("smoke-tet");
        assert_eq!(mesh.id, "smoke-tet");
        assert_eq!(mesh.nodes.len(), 4);
        // One block of one Tet4.
        assert_eq!(mesh.element_blocks.len(), 1);
        let blk = &mesh.element_blocks[0];
        assert_eq!(blk.element_type, valenx_mesh::ElementType::Tet4);
        // 1 tet × 4 nodes per tet = 4 connectivity indices.
        assert_eq!(blk.connectivity, vec![0, 1, 2, 3]);
        // One scalar field landed.
        assert_eq!(fields.len(), 1);
        assert_eq!(fields[0].name, "T");
        assert!(matches!(fields[0].kind, crate::FieldKind::Scalar));
        assert!(matches!(fields[0].location, crate::Location::OnNode));
        let (lo, hi) = fields[0].range.expect("range");
        assert!((lo - 273.15).abs() < 1e-3);
        assert!((hi - 300.0).abs() < 1e-3);
    }

    #[test]
    fn to_canonical_skips_unknown_cell_types() {
        // Synthesise a file with an unknown cell type code (99).
        // The mesh should come out empty (no canonical mapping)
        // while points / fields still land.
        let mut bytes: Vec<u8> = Vec::new();
        bytes.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        bytes.extend_from_slice(b"unknown\n");
        bytes.extend_from_slice(b"BINARY\n");
        bytes.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        bytes.extend_from_slice(b"POINTS 3 float\n");
        for v in [0.0_f32, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0] {
            bytes.extend_from_slice(&v.to_be_bytes());
        }
        bytes.push(b'\n');
        bytes.extend_from_slice(b"CELLS 1 4\n");
        for v in [3u32, 0, 1, 2] {
            bytes.extend_from_slice(&v.to_be_bytes());
        }
        bytes.push(b'\n');
        bytes.extend_from_slice(b"CELL_TYPES 1\n");
        bytes.extend_from_slice(&99u32.to_be_bytes()); // unknown
        bytes.push(b'\n');
        let data = parse_binary(&bytes).expect("parse");
        let (mesh, fields) = data.to_canonical("unknown-test");
        assert_eq!(mesh.nodes.len(), 3);
        assert!(mesh.element_blocks.is_empty(), "unknown cells skipped");
        assert!(fields.is_empty());
    }

    /// Round-5 RED→GREEN: a malformed VTK legacy CELLS block whose
    /// first u32 count is `u32::MAX` would previously panic the
    /// `CellIter` with a slice OOB (`start + n` overflows AND/OR
    /// slices past `data.len()`). The fix bounds-checks the implied
    /// end and stops iteration cleanly instead.
    #[test]
    fn cells_iter_rejects_oversized_n() {
        // Construct a LegacyData where `cells = [u32::MAX]` — one
        // pretend cell whose declared length is u32::MAX. Previously
        // iterating this panicked; now it short-circuits.
        let data = LegacyData {
            title: "stub".into(),
            points: Vec::new(),
            cells: vec![u32::MAX],
            cell_types: vec![10u8], // VTK_TETRA placeholder
            point_data: Vec::new(),
            cell_data: Vec::new(),
        };
        // Iteration must NOT panic; it must yield zero cells (the
        // malformed entry is dropped) and terminate.
        let collected: Vec<_> = data.iter_cells().collect();
        assert!(
            collected.is_empty(),
            "malformed CELLS block must drop the oversized entry, got {} cells",
            collected.len()
        );
    }

    /// Round-5: same iterator must also not panic when `cell_types`
    /// has fewer entries than CELLS (a class of malformed file where
    /// the two parallel arrays disagree on length).
    #[test]
    fn cells_iter_handles_cell_types_shorter_than_cells() {
        // One cell with one node, but no cell_types entry. The
        // pre-fix iterator panicked on `self.types[self.cell_idx]`.
        let data = LegacyData {
            title: "stub".into(),
            points: Vec::new(),
            cells: vec![1u32, 0u32],
            cell_types: Vec::new(),
            point_data: Vec::new(),
            cell_data: Vec::new(),
        };
        let collected: Vec<_> = data.iter_cells().collect();
        assert!(collected.is_empty());
    }

    /// Round-14 H2 RED→GREEN: a POINTS header declaring 10 billion
    /// points must be refused before the parser tries to allocate
    /// the implied bytes (~120 GB at 12 bytes/point). Pre-fix the
    /// `n * 3 * elem_size` multiplication overflowed usize on 64-bit,
    /// either underallocating or panicking via Vec; the cap fires
    /// cleanly with a TooLarge error and zero allocation.
    #[test]
    fn parse_binary_rejects_oversized_points_count() {
        let mut bytes: Vec<u8> = Vec::new();
        bytes.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        bytes.extend_from_slice(b"oversized\n");
        bytes.extend_from_slice(b"BINARY\n");
        bytes.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        // 10 billion points — well above the 256 M cap.
        bytes.extend_from_slice(b"POINTS 10000000000 float\n");
        // No payload follows on purpose — the cap must fire first.
        let err = parse_binary(&bytes).expect_err("must reject oversized POINTS");
        match err {
            ParseError::TooLarge { what, count, max } => {
                assert_eq!(what, "POINTS");
                assert_eq!(count, 10_000_000_000);
                assert_eq!(max, MAX_VTK_LEGACY_POINTS);
            }
            other => panic!("expected TooLarge, got {other:?}"),
        }
    }

    /// Round-14 H2 RED→GREEN sister: the CELLS path uses the same
    /// guard so a malicious `CELLS 0 18446744073709551615` (which
    /// would overflow `total_size * 4`) is refused before any
    /// allocation.
    #[test]
    fn parse_binary_rejects_oversized_cells_total_size() {
        let mut bytes: Vec<u8> = Vec::new();
        bytes.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        bytes.extend_from_slice(b"oversized cells\n");
        bytes.extend_from_slice(b"BINARY\n");
        bytes.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        bytes.extend_from_slice(b"POINTS 0 float\n");
        // Total connectivity-buffer size well past the cap.
        bytes.extend_from_slice(b"CELLS 1 9999999999\n");
        let err = parse_binary(&bytes).expect_err("must reject oversized CELLS");
        match err {
            ParseError::TooLarge { what, count, max } => {
                assert_eq!(what, "CELLS");
                assert_eq!(count, 9_999_999_999);
                assert_eq!(max, MAX_VTK_LEGACY_POINTS);
            }
            other => panic!("expected TooLarge, got {other:?}"),
        }
    }

    // -----------------------------------------------------------------
    // Round-15 L1 RED→GREEN: extend the round-14 H2 cap to the
    // POINT_DATA / CELL_DATA tuples count. Pre-fix the only check on
    // the POINT_DATA `n` was implicit via the eventual `tuples *
    // components` multiplication in `read_typed_block`, which would
    // overflow usize for `n = u64::MAX`-class values BEFORE the cap
    // was reached. The fix caps `n` at the POINT_DATA / CELL_DATA
    // parse site so we refuse the file before allocating the
    // pathological inner array.
    // -----------------------------------------------------------------

    #[test]
    fn parse_binary_rejects_oversized_point_data_count() {
        let mut bytes: Vec<u8> = Vec::new();
        bytes.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        bytes.extend_from_slice(b"oversized point_data\n");
        bytes.extend_from_slice(b"BINARY\n");
        bytes.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        bytes.extend_from_slice(b"POINTS 0 float\n");
        bytes.extend_from_slice(b"CELLS 0 0\n");
        bytes.extend_from_slice(b"CELL_TYPES 0\n");
        // 10 billion tuples — well above the 256 M cap.
        bytes.extend_from_slice(b"POINT_DATA 10000000000\n");
        bytes.extend_from_slice(b"SCALARS T float 1\n");
        bytes.extend_from_slice(b"LOOKUP_TABLE default\n");
        // No payload follows — the cap must fire before
        // read_typed_block is reached.
        let err = parse_binary(&bytes).expect_err("must reject oversized POINT_DATA");
        match err {
            ParseError::TooLarge { what, count, max } => {
                assert_eq!(what, "POINT_DATA");
                assert_eq!(count, 10_000_000_000);
                assert_eq!(max, MAX_VTK_LEGACY_POINTS);
            }
            other => panic!("expected TooLarge for POINT_DATA, got {other:?}"),
        }
    }

    #[test]
    fn parse_binary_rejects_oversized_cell_data_count() {
        let mut bytes: Vec<u8> = Vec::new();
        bytes.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        bytes.extend_from_slice(b"oversized cell_data\n");
        bytes.extend_from_slice(b"BINARY\n");
        bytes.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        bytes.extend_from_slice(b"POINTS 0 float\n");
        bytes.extend_from_slice(b"CELLS 0 0\n");
        bytes.extend_from_slice(b"CELL_TYPES 0\n");
        // 9 billion tuples on the CELL_DATA path — sister cap.
        bytes.extend_from_slice(b"CELL_DATA 9999999999\n");
        bytes.extend_from_slice(b"VECTORS U float\n");
        let err = parse_binary(&bytes).expect_err("must reject oversized CELL_DATA");
        match err {
            ParseError::TooLarge { what, count, max } => {
                assert_eq!(what, "CELL_DATA");
                assert_eq!(count, 9_999_999_999);
                assert_eq!(max, MAX_VTK_LEGACY_POINTS);
            }
            other => panic!("expected TooLarge for CELL_DATA, got {other:?}"),
        }
    }

    /// Round-15 L1 sister: the SCALARS / VECTORS `tuples * components`
    /// multiplication must use `checked_mul` so `tuples` near
    /// `usize::MAX / 2` with `components = 3` doesn't overflow into a
    /// small value before the cap inside `read_typed_block` is
    /// reached. With the POINT_DATA cap already at 256 M, the
    /// multiplication path is now unreachable via legit headers, but
    /// the checked_mul defends against future cap changes and remains
    /// the right hardening.
    #[test]
    fn parse_binary_rejects_tuples_times_components_overflow() {
        // We can't actually express `tuples > 256M * components > 1`
        // through a single POINT_DATA header without first tripping
        // the new POINT_DATA cap — which is the point. Pin that the
        // cap fires for the tuples-side input regardless of the
        // downstream multiplication.
        let mut bytes: Vec<u8> = Vec::new();
        bytes.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        bytes.extend_from_slice(b"overflow check\n");
        bytes.extend_from_slice(b"BINARY\n");
        bytes.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        bytes.extend_from_slice(b"POINTS 0 float\n");
        bytes.extend_from_slice(b"CELLS 0 0\n");
        bytes.extend_from_slice(b"CELL_TYPES 0\n");
        // Just over the cap — must be rejected.
        bytes.extend_from_slice(b"POINT_DATA 300000000\n");
        let err = parse_binary(&bytes).expect_err("must reject");
        assert!(matches!(err, ParseError::TooLarge { .. }));
    }

    #[test]
    fn parse_binary_handles_vector_arrays_in_point_data() {
        let mut bytes: Vec<u8> = Vec::new();
        bytes.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        bytes.extend_from_slice(b"vec test\n");
        bytes.extend_from_slice(b"BINARY\n");
        bytes.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        bytes.extend_from_slice(b"POINTS 2 float\n");
        for v in [0.0_f32, 0.0, 0.0, 1.0, 0.0, 0.0] {
            bytes.extend_from_slice(&v.to_be_bytes());
        }
        bytes.push(b'\n');
        bytes.extend_from_slice(b"CELLS 0 0\n");
        bytes.extend_from_slice(b"CELL_TYPES 0\n");
        bytes.extend_from_slice(b"POINT_DATA 2\n");
        bytes.extend_from_slice(b"VECTORS U float\n");
        // 2 points * 3 components = 6 floats
        for v in [10.0_f32, 0.0, 0.0, 20.0, 0.0, 0.0] {
            bytes.extend_from_slice(&v.to_be_bytes());
        }
        bytes.push(b'\n');
        let data = parse_binary(&bytes).expect("parse");
        assert_eq!(data.point_data.len(), 1);
        let v = &data.point_data[0];
        assert_eq!(v.name, "U");
        assert_eq!(v.components, 3);
        assert_eq!(v.tuples, 2);
        assert_eq!(v.data, vec![10.0, 0.0, 0.0, 20.0, 0.0, 0.0]);
    }
}
