//! Stanford PLY reader and writer.
//!
//! Scope: positions (`element vertex N` with `property float x/y/z`)
//! and triangular / fan-triangulated polygon faces (`element face M`
//! with `property list uchar int vertex_indices`).
//!
//! ## Format coverage
//!
//! - **ASCII PLY** (`format ascii 1.0`) — read + write.
//! - **Binary PLY** (`format binary_little_endian 1.0` /
//!   `binary_big_endian 1.0`) — read (Phase 26.5). The body decoder
//!   honours every PLY scalar type (`char`/`int8` … `double`/
//!   `float64`, plus the `uint8` … aliases) and variable-length
//!   list properties with an arbitrary count type. The writer still
//!   emits ASCII; a binary writer is a follow-up.
//!
//! Custom properties beyond `x/y/z` on vertices (e.g. colors,
//! normals) are tolerated: the value is parsed and discarded so the
//! per-element record stays correctly aligned.

use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

use nalgebra::Vector3;
use thiserror::Error;

use crate::element::{ElementBlock, ElementType};
use crate::mesh::Mesh;

/// Upper bound on the per-element list count we will trust from a PLY
/// file. PLY list-property headers store the count as a `uchar`, `ushort`,
/// `int`, etc.; a hostile or corrupted file can claim a per-face vertex
/// count of `i32::MAX`, which a naive `Vec::with_capacity(n)` translates
/// into a ~16 GB allocation. Legitimate mesh data never approaches this:
/// even a 1k-gon NURBS tessellation tops out around 10k indices. Cap at
/// 1M for slack but reject anything larger.
pub const MAX_PLY_LIST_LEN: usize = 1_000_000;

/// R29 A: hard cap on the total bytes [`read_path`] will pull off
/// disk for a single PLY file. Mirrors
/// `valenx_core::io_caps::MAX_PLY_ASCII_BYTES` (1 GiB) — duplicated
/// inline rather than imported because valenx-core sits *upstream* of
/// this crate in the dependency graph (`valenx-core → valenx-fields →
/// valenx-mesh`), so valenx-mesh cannot add a `valenx-core` dependency
/// without forming a cycle. See the same rationale on
/// `MAX_OBJ_LINE_BYTES` in `format/obj.rs`. 1 GiB comfortably covers
/// dense real-world ASCII/binary meshes while refusing the
/// read-the-whole-file-into-RAM DoS shape that an unbounded
/// `read_to_end` exposes.
pub const MAX_PLY_FILE_BYTES: u64 = 1024 * 1024 * 1024;

/// All PLY-related errors.
#[derive(Debug, Error)]
pub enum PlyError {
    /// IO problem.
    #[error(transparent)]
    Io(#[from] io::Error),
    /// Malformed PLY content.
    #[error("malformed PLY: {0}")]
    Malformed(String),
    /// Per-element list count exceeds [`MAX_PLY_LIST_LEN`] — almost
    /// certainly a hostile or corrupted file. Refused before allocation
    /// so an attacker cannot trigger a multi-GB `Vec` reserve.
    #[error("PLY list length {count} exceeds maximum {max}")]
    ListTooLarge {
        /// Count claimed in the file header.
        count: usize,
        /// Cap enforced by this implementation.
        max: usize,
    },
    /// A face vertex-index list referenced a value the canonical
    /// connectivity layout (`u32`) cannot represent — negative,
    /// NaN, ±inf, or larger than the file's vertex count. Pre-fix the
    /// reader would do `value as u32` and silently saturate negatives
    /// to 0 and `> u32::MAX` to `u32::MAX`, producing wrong
    /// connectivity that pointed at random vertices.
    #[error("PLY vertex index {value} is invalid (negative / non-finite / out of range)")]
    InvalidIndex {
        /// The offending raw value read from the file.
        value: f64,
    },
    /// The file on disk is larger than the byte cap we are willing to
    /// read into memory ([`MAX_PLY_FILE_BYTES`]). Refused before the
    /// `read_to_end` so a hostile multi-GB file cannot OOM the import.
    #[error("PLY file size {size} bytes exceeds maximum {cap} bytes")]
    TooLarge {
        /// File length reported by the filesystem.
        size: u64,
        /// Cap enforced by this reader.
        cap: u64,
    },
}

/// PLY scalar type — covers every name and alias in the spec.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ScalarType {
    /// 8-bit signed.
    I8,
    /// 8-bit unsigned.
    U8,
    /// 16-bit signed.
    I16,
    /// 16-bit unsigned.
    U16,
    /// 32-bit signed.
    I32,
    /// 32-bit unsigned.
    U32,
    /// 32-bit IEEE float.
    F32,
    /// 64-bit IEEE float.
    F64,
}

impl ScalarType {
    /// Parse a PLY type token (`float`, `float32`, `uchar`, `uint8`, …).
    fn parse(tok: &str) -> Option<ScalarType> {
        Some(match tok {
            "char" | "int8" => ScalarType::I8,
            "uchar" | "uint8" => ScalarType::U8,
            "short" | "int16" => ScalarType::I16,
            "ushort" | "uint16" => ScalarType::U16,
            "int" | "int32" => ScalarType::I32,
            "uint" | "uint32" => ScalarType::U32,
            "float" | "float32" => ScalarType::F32,
            "double" | "float64" => ScalarType::F64,
            _ => return None,
        })
    }

    /// Byte width of the scalar.
    fn width(self) -> usize {
        match self {
            ScalarType::I8 | ScalarType::U8 => 1,
            ScalarType::I16 | ScalarType::U16 => 2,
            ScalarType::I32 | ScalarType::U32 | ScalarType::F32 => 4,
            ScalarType::F64 => 8,
        }
    }
}

/// A single property declaration.
#[derive(Debug, Clone)]
struct PropertySpec {
    /// Property name (e.g. `x`, `red`, `vertex_indices`).
    name: String,
    /// For a scalar property: its type. For a list property: the
    /// element type.
    value_type: ScalarType,
    /// `Some(count_type)` for a `property list` declaration, `None`
    /// for a plain scalar property.
    list_count_type: Option<ScalarType>,
}

#[derive(Debug, Clone)]
struct ElementSpec {
    name: String,
    count: usize,
    properties: Vec<PropertySpec>,
}

/// Wire format declared by the PLY header.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PlyFormat {
    /// Whitespace-delimited text body.
    Ascii,
    /// Little-endian binary body.
    BinaryLe,
    /// Big-endian binary body.
    BinaryBe,
}

/// Read a PLY file from `path` into a canonical [`Mesh`].
///
/// Both ASCII and binary (`binary_little_endian` /
/// `binary_big_endian`) PLY are supported. The mesh `id` is set to
/// the file stem.
pub fn read_path(path: impl AsRef<Path>) -> Result<Mesh, PlyError> {
    read_path_with_cap(path, MAX_PLY_FILE_BYTES)
}

/// Cap-aware core of [`read_path`]. Stats the file first and refuses
/// anything larger than `cap` with [`PlyError::TooLarge`], then bounds
/// the actual read with [`Read::take`] so a file that grows between the
/// `stat` and the `read_to_end` still cannot pull more than `cap` bytes
/// into RAM. Extracted so the cap can be exercised by a unit test with a
/// tiny `cap` (a real 1 GiB fixture is impractical to materialise).
fn read_path_with_cap(path: impl AsRef<Path>, cap: u64) -> Result<Mesh, PlyError> {
    let path_ref = path.as_ref();
    let mut file = File::open(path_ref)?;
    let id = path_ref
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("ply")
        .to_string();
    let size = file.metadata()?.len();
    if size > cap {
        return Err(PlyError::TooLarge { size, cap });
    }
    let mut bytes = Vec::new();
    // Defensive second line of bounding: even if the file grew after the
    // stat above, `take` guarantees we never read past `cap` bytes.
    Read::take(&mut file, cap).read_to_end(&mut bytes)?;
    read_bytes(id, &bytes)
}

/// Read an ASCII PLY from an in-memory string. The string must
/// contain the full header + body; line endings are platform-tolerant.
///
/// Binary PLY cannot round-trip through a `&str` (the body is not
/// UTF-8) — use [`read_bytes`] or [`read_path`] for binary files.
pub fn read_str(id: String, text: &str) -> Result<Mesh, PlyError> {
    read_bytes(id, text.as_bytes())
}

/// Parse the header text of a PLY file. Returns the format, the
/// element/property declarations, and the byte offset where the body
/// starts (the byte just past the `end_header` line terminator).
fn parse_header(bytes: &[u8]) -> Result<(PlyFormat, Vec<ElementSpec>, usize), PlyError> {
    // Locate `end_header` followed by a newline. The header is always
    // ASCII even for binary files, so we can scan bytewise.
    let marker = b"end_header";
    let mut header_end = None;
    let mut i = 0;
    while i + marker.len() <= bytes.len() {
        if &bytes[i..i + marker.len()] == marker {
            // Advance past the marker and the line terminator
            // (`\n` or `\r\n`).
            let mut j = i + marker.len();
            if j < bytes.len() && bytes[j] == b'\r' {
                j += 1;
            }
            if j < bytes.len() && bytes[j] == b'\n' {
                j += 1;
            }
            header_end = Some(j);
            break;
        }
        i += 1;
    }
    let body_start = header_end.ok_or_else(|| PlyError::Malformed("missing end_header".into()))?;
    let header_text = std::str::from_utf8(&bytes[..body_start])
        .map_err(|_| PlyError::Malformed("PLY header is not valid UTF-8".into()))?;

    let mut lines = header_text.lines();
    let first = lines
        .next()
        .ok_or_else(|| PlyError::Malformed("empty file".into()))?;
    if first.trim() != "ply" {
        return Err(PlyError::Malformed(format!(
            "expected 'ply' magic, got {first:?}"
        )));
    }
    let mut format: Option<PlyFormat> = None;
    let mut elements: Vec<ElementSpec> = Vec::new();
    for line in lines {
        let line = line.trim();
        if line.is_empty() || line.starts_with("comment") || line.starts_with("obj_info") {
            continue;
        }
        if let Some(rest) = line.strip_prefix("format ") {
            let mut iter = rest.split_whitespace();
            let kind = iter.next().unwrap_or("");
            format = Some(match kind {
                "ascii" => PlyFormat::Ascii,
                "binary_little_endian" => PlyFormat::BinaryLe,
                "binary_big_endian" => PlyFormat::BinaryBe,
                other => return Err(PlyError::Malformed(format!("unknown PLY format {other:?}"))),
            });
        } else if let Some(rest) = line.strip_prefix("element ") {
            let mut iter = rest.split_whitespace();
            let name = iter.next().unwrap_or("").to_string();
            let count: usize = iter
                .next()
                .and_then(|s| s.parse().ok())
                .ok_or_else(|| PlyError::Malformed(format!("bad element line: {line}")))?;
            elements.push(ElementSpec {
                name,
                count,
                properties: Vec::new(),
            });
        } else if let Some(rest) = line.strip_prefix("property ") {
            let elem = elements
                .last_mut()
                .ok_or_else(|| PlyError::Malformed(format!("property without element: {line}")))?;
            let mut iter = rest.split_whitespace();
            let first = iter.next().unwrap_or("");
            if first == "list" {
                // property list <count_type> <index_type> <name>
                let count_tok = iter.next().unwrap_or("");
                let elem_tok = iter.next().unwrap_or("");
                let name = iter.next().unwrap_or("").to_string();
                let count_type = ScalarType::parse(count_tok).ok_or_else(|| {
                    PlyError::Malformed(format!("bad list count type: {count_tok:?}"))
                })?;
                let value_type = ScalarType::parse(elem_tok).ok_or_else(|| {
                    PlyError::Malformed(format!("bad list element type: {elem_tok:?}"))
                })?;
                elem.properties.push(PropertySpec {
                    name,
                    value_type,
                    list_count_type: Some(count_type),
                });
            } else {
                // property <type> <name>
                let name = iter.next().unwrap_or("").to_string();
                let value_type = ScalarType::parse(first)
                    .ok_or_else(|| PlyError::Malformed(format!("bad property type: {first:?}")))?;
                elem.properties.push(PropertySpec {
                    name,
                    value_type,
                    list_count_type: None,
                });
            }
        }
        // (`end_header` is consumed by the byte scan above.)
    }
    let format = format.ok_or_else(|| PlyError::Malformed("missing format declaration".into()))?;
    Ok((format, elements, body_start))
}

/// Read a PLY file from a raw byte buffer (header + body). Dispatches
/// to the ASCII or binary body decoder based on the `format` line.
pub fn read_bytes(id: String, bytes: &[u8]) -> Result<Mesh, PlyError> {
    let (format, elements, body_start) = parse_header(bytes)?;
    match format {
        PlyFormat::Ascii => {
            let body = std::str::from_utf8(&bytes[body_start..])
                .map_err(|_| PlyError::Malformed("ASCII PLY body is not UTF-8".into()))?;
            read_ascii_body(id, &elements, body)
        }
        PlyFormat::BinaryLe => read_binary_body(id, &elements, &bytes[body_start..], false),
        PlyFormat::BinaryBe => read_binary_body(id, &elements, &bytes[body_start..], true),
    }
}

/// Decode the ASCII body of a PLY file given the parsed header.
fn read_ascii_body(id: String, elements: &[ElementSpec], body: &str) -> Result<Mesh, PlyError> {
    let mut lines = body.lines();
    let mut mesh = Mesh::new(id);
    let mut tri_conn: Vec<u32> = Vec::new();
    // Body: for each element, read `count` lines, each having N values
    // matching the property declarations.
    for elem in elements {
        match elem.name.as_str() {
            "vertex" => {
                let ix = elem.properties.iter().position(|p| p.name == "x");
                let iy = elem.properties.iter().position(|p| p.name == "y");
                let iz = elem.properties.iter().position(|p| p.name == "z");
                let (ix, iy, iz) = match (ix, iy, iz) {
                    (Some(a), Some(b), Some(c)) => (a, b, c),
                    _ => {
                        return Err(PlyError::Malformed(
                            "vertex element missing x/y/z property".into(),
                        ))
                    }
                };
                for _ in 0..elem.count {
                    let line = lines.next().ok_or_else(|| {
                        PlyError::Malformed("ran out of lines reading vertices".into())
                    })?;
                    let toks: Vec<&str> = line.split_whitespace().collect();
                    if toks.len() < elem.properties.len() {
                        return Err(PlyError::Malformed(format!(
                            "vertex needs {} tokens, got {}: {line:?}",
                            elem.properties.len(),
                            toks.len()
                        )));
                    }
                    let x: f64 = toks[ix]
                        .parse()
                        .map_err(|_| PlyError::Malformed(format!("bad x: {:?}", toks[ix])))?;
                    let y: f64 = toks[iy]
                        .parse()
                        .map_err(|_| PlyError::Malformed(format!("bad y: {:?}", toks[iy])))?;
                    let z: f64 = toks[iz]
                        .parse()
                        .map_err(|_| PlyError::Malformed(format!("bad z: {:?}", toks[iz])))?;
                    mesh.nodes.push(Vector3::new(x, y, z));
                }
            }
            "face" => {
                // Single list property `vertex_indices`. We support the
                // canonical layout: <n> v0 v1 ... v(n-1).
                //
                // R34 L1: PLY declares `element vertex N` before
                // `element face M`, and we walk elements in declaration
                // order, so `mesh.nodes` is fully built here. Reject any
                // face index `>= vertex_count` — a face citing vertex
                // 999 in a 3-vertex file otherwise parses into
                // connectivity and panics in the shared consumers
                // (`decimate` / `boolean`) that raw-index
                // `positions[tri[k]]`.
                let vertex_count = mesh.nodes.len() as i64;
                for _ in 0..elem.count {
                    let line = lines.next().ok_or_else(|| {
                        PlyError::Malformed("ran out of lines reading faces".into())
                    })?;
                    let toks: Vec<&str> = line.split_whitespace().collect();
                    if toks.is_empty() {
                        return Err(PlyError::Malformed("empty face line".into()));
                    }
                    let n: usize = toks[0].parse().map_err(|_| {
                        PlyError::Malformed(format!("bad face count: {:?}", toks[0]))
                    })?;
                    if n > MAX_PLY_LIST_LEN {
                        return Err(PlyError::ListTooLarge {
                            count: n,
                            max: MAX_PLY_LIST_LEN,
                        });
                    }
                    if toks.len() < 1 + n {
                        return Err(PlyError::Malformed(format!(
                            "face line has count {n} but only {} indices",
                            toks.len() - 1
                        )));
                    }
                    let mut idxs: Vec<u32> = Vec::with_capacity(n);
                    for k in 0..n {
                        let raw = toks[1 + k];
                        // Round-3 fix: parse as i64 first so we can
                        // distinguish "out-of-range u32" (e.g. negative,
                        // which is a real spec-allowed-but-corrupt index)
                        // from "junk that doesn't even look like an
                        // integer". Negatives become InvalidIndex,
                        // garbage stays Malformed.
                        match raw.parse::<i64>() {
                            Ok(signed) => {
                                // R34 L1: the upper bound is the vertex
                                // count, not just u32::MAX — an index
                                // within u32 but past the mesh's vertices
                                // is still out of range and would panic
                                // downstream.
                                if signed < 0 || signed >= vertex_count {
                                    return Err(PlyError::InvalidIndex {
                                        value: signed as f64,
                                    });
                                }
                                idxs.push(signed as u32);
                            }
                            Err(_) => {
                                return Err(PlyError::Malformed(format!("bad index: {raw:?}")));
                            }
                        }
                    }
                    // Fan triangulate.
                    for k in 1..(n.saturating_sub(1)) {
                        tri_conn.extend_from_slice(&[idxs[0], idxs[k], idxs[k + 1]]);
                    }
                }
            }
            _ => {
                // Unknown element: skip `count` lines.
                for _ in 0..elem.count {
                    if lines.next().is_none() {
                        return Err(PlyError::Malformed(format!(
                            "ran out of lines skipping element {}",
                            elem.name
                        )));
                    }
                }
            }
        }
    }
    if !tri_conn.is_empty() {
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = tri_conn;
        mesh.element_blocks.push(blk);
    }
    mesh.recompute_stats();
    Ok(mesh)
}

/// Cursor over the binary body — tracks a byte offset and decodes
/// typed scalars with the file's endianness.
struct BinCursor<'a> {
    buf: &'a [u8],
    pos: usize,
    big_endian: bool,
}

impl<'a> BinCursor<'a> {
    fn new(buf: &'a [u8], big_endian: bool) -> Self {
        Self {
            buf,
            pos: 0,
            big_endian,
        }
    }

    /// Take the next `n` bytes, erroring on a short body.
    fn take(&mut self, n: usize) -> Result<&'a [u8], PlyError> {
        if self.pos + n > self.buf.len() {
            return Err(PlyError::Malformed(format!(
                "binary PLY body truncated: needed {n} more bytes at offset {}",
                self.pos
            )));
        }
        let slice = &self.buf[self.pos..self.pos + n];
        self.pos += n;
        Ok(slice)
    }

    /// Decode one scalar of `ty`, widening it to `f64` (integers and
    /// floats alike — the caller narrows as needed).
    fn scalar(&mut self, ty: ScalarType) -> Result<f64, PlyError> {
        let raw = self.take(ty.width())?;
        Ok(match ty {
            ScalarType::I8 => raw[0] as i8 as f64,
            ScalarType::U8 => raw[0] as f64,
            ScalarType::I16 => {
                let v = if self.big_endian {
                    i16::from_be_bytes([raw[0], raw[1]])
                } else {
                    i16::from_le_bytes([raw[0], raw[1]])
                };
                v as f64
            }
            ScalarType::U16 => {
                let v = if self.big_endian {
                    u16::from_be_bytes([raw[0], raw[1]])
                } else {
                    u16::from_le_bytes([raw[0], raw[1]])
                };
                v as f64
            }
            ScalarType::I32 => {
                let a = [raw[0], raw[1], raw[2], raw[3]];
                let v = if self.big_endian {
                    i32::from_be_bytes(a)
                } else {
                    i32::from_le_bytes(a)
                };
                v as f64
            }
            ScalarType::U32 => {
                let a = [raw[0], raw[1], raw[2], raw[3]];
                let v = if self.big_endian {
                    u32::from_be_bytes(a)
                } else {
                    u32::from_le_bytes(a)
                };
                v as f64
            }
            ScalarType::F32 => {
                let a = [raw[0], raw[1], raw[2], raw[3]];
                let v = if self.big_endian {
                    f32::from_be_bytes(a)
                } else {
                    f32::from_le_bytes(a)
                };
                v as f64
            }
            ScalarType::F64 => {
                let a = [
                    raw[0], raw[1], raw[2], raw[3], raw[4], raw[5], raw[6], raw[7],
                ];
                if self.big_endian {
                    f64::from_be_bytes(a)
                } else {
                    f64::from_le_bytes(a)
                }
            }
        })
    }
}

/// Decode the binary body of a PLY file given the parsed header.
///
/// Walks elements in declaration order. For the `vertex` element it
/// extracts x/y/z; for the `face` element it extracts the
/// `vertex_indices` list and fan-triangulates. Every other property
/// (and every other element) is consumed byte-for-byte so the cursor
/// stays aligned.
/// Decode the binary body of a PLY file.
///
/// List-property counts are validated against [`MAX_PLY_LIST_LEN`] before
/// any `Vec::with_capacity` call. A file claiming `count = i32::MAX`
/// returns [`PlyError::ListTooLarge`] rather than attempting to reserve
/// gigabytes of memory.
fn read_binary_body(
    id: String,
    elements: &[ElementSpec],
    body: &[u8],
    big_endian: bool,
) -> Result<Mesh, PlyError> {
    let mut cur = BinCursor::new(body, big_endian);
    let mut mesh = Mesh::new(id);
    let mut tri_conn: Vec<u32> = Vec::new();

    for elem in elements {
        // Resolve x/y/z property indices once for the vertex element.
        let xyz = if elem.name == "vertex" {
            let ix = elem.properties.iter().position(|p| p.name == "x");
            let iy = elem.properties.iter().position(|p| p.name == "y");
            let iz = elem.properties.iter().position(|p| p.name == "z");
            match (ix, iy, iz) {
                (Some(a), Some(b), Some(c)) => Some((a, b, c)),
                _ => {
                    return Err(PlyError::Malformed(
                        "vertex element missing x/y/z property".into(),
                    ))
                }
            }
        } else {
            None
        };
        let vi_index = if elem.name == "face" {
            elem.properties
                .iter()
                .position(|p| p.list_count_type.is_some())
        } else {
            None
        };

        for _ in 0..elem.count {
            // Decode each property of this record, capturing the ones
            // we care about.
            let mut scalars: Vec<f64> = Vec::with_capacity(elem.properties.len());
            let mut lists: Vec<Vec<f64>> = vec![Vec::new(); elem.properties.len()];
            for (pi, prop) in elem.properties.iter().enumerate() {
                match prop.list_count_type {
                    None => {
                        scalars.push(cur.scalar(prop.value_type)?);
                    }
                    Some(count_ty) => {
                        // List: a count followed by `count` elements.
                        scalars.push(f64::NAN); // placeholder to keep indices aligned
                        let n = cur.scalar(count_ty)? as usize;
                        if n > MAX_PLY_LIST_LEN {
                            return Err(PlyError::ListTooLarge {
                                count: n,
                                max: MAX_PLY_LIST_LEN,
                            });
                        }
                        let mut items = Vec::with_capacity(n);
                        for _ in 0..n {
                            items.push(cur.scalar(prop.value_type)?);
                        }
                        lists[pi] = items;
                    }
                }
            }

            if let Some((ix, iy, iz)) = xyz {
                mesh.nodes
                    .push(Vector3::new(scalars[ix], scalars[iy], scalars[iz]));
            }
            if let Some(vi) = vi_index {
                let idxs = &lists[vi];
                if idxs.len() >= 3 {
                    // Round-3 fix: validate each index before the `as
                    // u32` cast. Negative values would otherwise silently
                    // saturate to 0 (which points at a real vertex —
                    // wrong but never out-of-bounds), and `> u32::MAX`
                    // would saturate to `u32::MAX`. Either produces
                    // garbage connectivity. NaN-as-u32 is also 0 by
                    // Rust's saturating-cast contract.
                    //
                    // R34 L1: the `vertex` element precedes `face` in
                    // declaration order, so `mesh.nodes` is fully built
                    // — pass the vertex count so an index within u32 but
                    // past the mesh's vertices is rejected too (else it
                    // panics in the shared `decimate`/`boolean`
                    // consumers that raw-index `positions[tri[k]]`).
                    let vertex_count = mesh.nodes.len();
                    let i0 = checked_index(idxs[0], vertex_count)?;
                    for k in 1..idxs.len() - 1 {
                        tri_conn.extend_from_slice(&[
                            i0,
                            checked_index(idxs[k], vertex_count)?,
                            checked_index(idxs[k + 1], vertex_count)?,
                        ]);
                    }
                }
            }
        }
    }

    if !tri_conn.is_empty() {
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = tri_conn;
        mesh.element_blocks.push(blk);
    }
    mesh.recompute_stats();
    Ok(mesh)
}

/// Validate a face-vertex index read as an `f64` (the canonical
/// scalar in this reader) before casting to `u32`. Rejects negatives,
/// NaN, ±inf, values > `u32::MAX`, and (R34 L1) any index `>=
/// vertex_count` — an in-`u32` value past the mesh's vertices is still
/// out of range and would panic the shared `decimate`/`boolean`
/// consumers that raw-index `positions[tri[k]]`. Without this, `value
/// as u32` silently saturates negatives to 0 (which points at a real
/// vertex — wrong but never out-of-bounds), so downstream never notices.
fn checked_index(value: f64, vertex_count: usize) -> Result<u32, PlyError> {
    if !value.is_finite() || value < 0.0 || value > u32::MAX as f64 {
        return Err(PlyError::InvalidIndex { value });
    }
    // Finite, non-negative, in u32 range — now bound by the vertex
    // array length so the index actually addresses a real vertex.
    if value >= vertex_count as f64 {
        return Err(PlyError::InvalidIndex { value });
    }
    // Safe to cast — finite, non-negative, in range.
    Ok(value as u32)
}

/// Write the Tri3 surface of `mesh` to `path` as ASCII PLY.
pub fn write_path(mesh: &Mesh, path: impl AsRef<Path>) -> Result<(), PlyError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    write_to(mesh, &mut writer)?;
    writer.flush()?;
    Ok(())
}

/// Write Tri3 blocks of `mesh` to a `Write` sink as ASCII PLY.
pub fn write_to<W: Write>(mesh: &Mesh, mut sink: W) -> io::Result<()> {
    let tri_count: usize = mesh
        .element_blocks
        .iter()
        .filter(|b| b.element_type == ElementType::Tri3)
        .map(|b| b.connectivity.len() / 3)
        .sum();
    writeln!(sink, "ply")?;
    writeln!(sink, "format ascii 1.0")?;
    writeln!(sink, "comment generated by valenx-mesh")?;
    writeln!(sink, "element vertex {}", mesh.nodes.len())?;
    writeln!(sink, "property float x")?;
    writeln!(sink, "property float y")?;
    writeln!(sink, "property float z")?;
    writeln!(sink, "element face {tri_count}")?;
    writeln!(sink, "property list uchar int vertex_indices")?;
    writeln!(sink, "end_header")?;
    for n in &mesh.nodes {
        writeln!(sink, "{} {} {}", fmt_f(n.x), fmt_f(n.y), fmt_f(n.z))?;
    }
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            writeln!(sink, "3 {} {} {}", tri[0], tri[1], tri[2])?;
        }
    }
    Ok(())
}

fn fmt_f(x: f64) -> String {
    let s = format!("{x:.10}");
    let trimmed = s.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() || trimmed == "-" {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_minimal_triangle() {
        let ply = "\
ply
format ascii 1.0
element vertex 3
property float x
property float y
property float z
element face 1
property list uchar int vertex_indices
end_header
0 0 0
1 0 0
0 1 0
3 0 1 2
";
        let m = read_str("tri".into(), ply).unwrap();
        assert_eq!(m.nodes.len(), 3);
        assert_eq!(m.element_blocks[0].connectivity, vec![0, 1, 2]);
    }

    #[test]
    fn fan_triangulates_quad() {
        let ply = "\
ply
format ascii 1.0
element vertex 4
property float x
property float y
property float z
element face 1
property list uchar int vertex_indices
end_header
0 0 0
1 0 0
1 1 0
0 1 0
4 0 1 2 3
";
        let m = read_str("quad".into(), ply).unwrap();
        assert_eq!(m.element_blocks[0].connectivity, vec![0, 1, 2, 0, 2, 3]);
    }

    #[test]
    fn extra_vertex_properties_are_tolerated() {
        // Adds a red/green/blue triple after xyz — must still parse.
        let ply = "\
ply
format ascii 1.0
element vertex 3
property float x
property float y
property float z
property uchar red
property uchar green
property uchar blue
element face 1
property list uchar int vertex_indices
end_header
0 0 0 255 0 0
1 0 0 0 255 0
0 1 0 0 0 255
3 0 1 2
";
        let m = read_str("rgb".into(), ply).unwrap();
        assert_eq!(m.nodes.len(), 3);
        assert!((m.nodes[0] - Vector3::new(0.0, 0.0, 0.0)).norm() < 1e-9);
    }

    #[test]
    fn unknown_format_is_rejected() {
        let ply = "\
ply
format banana 1.0
end_header
";
        let err = read_str("bad".into(), ply).unwrap_err();
        assert!(matches!(err, PlyError::Malformed(_)));
    }

    #[test]
    fn binary_le_triangle_round_trips() {
        // Hand-assemble a binary_little_endian PLY: 3 float32 vertices,
        // one uchar-count + int32 face. Phase 26.5.
        let header = b"ply\n\
format binary_little_endian 1.0\n\
element vertex 3\n\
property float x\n\
property float y\n\
property float z\n\
element face 1\n\
property list uchar int vertex_indices\n\
end_header\n";
        let mut bytes = header.to_vec();
        // Vertices (x, y, z) as little-endian f32.
        for &(x, y, z) in &[
            (0.0_f32, 0.0_f32, 0.0_f32),
            (2.0, 0.0, 0.0),
            (0.0, 3.0, 0.0),
        ] {
            bytes.extend_from_slice(&x.to_le_bytes());
            bytes.extend_from_slice(&y.to_le_bytes());
            bytes.extend_from_slice(&z.to_le_bytes());
        }
        // Face: count=3 (uchar), then three int32 indices.
        bytes.push(3u8);
        for i in [0i32, 1, 2] {
            bytes.extend_from_slice(&i.to_le_bytes());
        }
        let m = read_bytes("binle".into(), &bytes).unwrap();
        assert_eq!(m.nodes.len(), 3);
        assert!((m.nodes[1] - Vector3::new(2.0, 0.0, 0.0)).norm() < 1e-6);
        assert!((m.nodes[2] - Vector3::new(0.0, 3.0, 0.0)).norm() < 1e-6);
        assert_eq!(m.element_blocks[0].connectivity, vec![0, 1, 2]);
    }

    #[test]
    fn binary_be_quad_fan_triangulates() {
        // binary_big_endian, 4 double vertices, one quad face — must
        // fan-triangulate to two triangles.
        let header = b"ply\n\
format binary_big_endian 1.0\n\
element vertex 4\n\
property double x\n\
property double y\n\
property double z\n\
element face 1\n\
property list uchar int vertex_indices\n\
end_header\n";
        let mut bytes = header.to_vec();
        for &(x, y, z) in &[
            (0.0_f64, 0.0_f64, 0.0_f64),
            (1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
        ] {
            bytes.extend_from_slice(&x.to_be_bytes());
            bytes.extend_from_slice(&y.to_be_bytes());
            bytes.extend_from_slice(&z.to_be_bytes());
        }
        bytes.push(4u8);
        for i in [0i32, 1, 2, 3] {
            bytes.extend_from_slice(&i.to_be_bytes());
        }
        let m = read_bytes("binbe".into(), &bytes).unwrap();
        assert_eq!(m.nodes.len(), 4);
        assert_eq!(m.element_blocks[0].connectivity, vec![0, 1, 2, 0, 2, 3]);
    }

    #[test]
    fn binary_extra_properties_keep_record_aligned() {
        // A binary vertex record with x/y/z floats plus a uchar
        // colour triple — the colour bytes must be consumed so the
        // following vertex still decodes correctly.
        let header = b"ply\n\
format binary_little_endian 1.0\n\
element vertex 2\n\
property float x\n\
property float y\n\
property float z\n\
property uchar red\n\
property uchar green\n\
property uchar blue\n\
end_header\n";
        let mut bytes = header.to_vec();
        for &(x, y, z, r, g, b) in &[
            (1.0_f32, 2.0_f32, 3.0_f32, 10u8, 20u8, 30u8),
            (4.0, 5.0, 6.0, 40, 50, 60),
        ] {
            bytes.extend_from_slice(&x.to_le_bytes());
            bytes.extend_from_slice(&y.to_le_bytes());
            bytes.extend_from_slice(&z.to_le_bytes());
            bytes.push(r);
            bytes.push(g);
            bytes.push(b);
        }
        let m = read_bytes("binrgb".into(), &bytes).unwrap();
        assert_eq!(m.nodes.len(), 2);
        assert!((m.nodes[0] - Vector3::new(1.0, 2.0, 3.0)).norm() < 1e-6);
        assert!((m.nodes[1] - Vector3::new(4.0, 5.0, 6.0)).norm() < 1e-6);
    }

    #[test]
    fn truncated_binary_body_errors() {
        let header = b"ply\n\
format binary_little_endian 1.0\n\
element vertex 2\n\
property float x\n\
property float y\n\
property float z\n\
end_header\n";
        let mut bytes = header.to_vec();
        // Only one vertex's worth of data — the second read must fail.
        bytes.extend_from_slice(&0.0_f32.to_le_bytes());
        bytes.extend_from_slice(&0.0_f32.to_le_bytes());
        bytes.extend_from_slice(&0.0_f32.to_le_bytes());
        let err = read_bytes("trunc".into(), &bytes).unwrap_err();
        assert!(matches!(err, PlyError::Malformed(_)));
    }

    #[test]
    fn missing_ply_magic_is_rejected() {
        let err = read_str("bad".into(), "nope\n").unwrap_err();
        assert!(matches!(err, PlyError::Malformed(_)));
    }

    #[test]
    fn ascii_oversized_face_count_is_rejected() {
        // ASCII PLY claiming a face count of 5 million — exceeds
        // MAX_PLY_LIST_LEN, must be refused before allocation. The
        // trailing data is deliberately truncated; the cap fires first.
        let ply = "\
ply
format ascii 1.0
element vertex 3
property float x
property float y
property float z
element face 1
property list uchar int vertex_indices
end_header
0 0 0
1 0 0
0 1 0
5000000 0 1 2
";
        let err = read_str("oversize".into(), ply).unwrap_err();
        assert!(
            matches!(
                err,
                PlyError::ListTooLarge {
                    count: 5_000_000,
                    max: MAX_PLY_LIST_LEN
                }
            ),
            "expected ListTooLarge, got {err:?}"
        );
    }

    #[test]
    fn binary_oversized_list_count_is_rejected() {
        // binary_little_endian PLY with a list count = u32::MAX disguised
        // as int — must reject before reserving the dropper allocation.
        let header = b"ply\n\
format binary_little_endian 1.0\n\
element vertex 3\n\
property float x\n\
property float y\n\
property float z\n\
element face 1\n\
property list int int vertex_indices\n\
end_header\n";
        let mut bytes = header.to_vec();
        for _ in 0..3 {
            bytes.extend_from_slice(&0.0_f32.to_le_bytes());
            bytes.extend_from_slice(&0.0_f32.to_le_bytes());
            bytes.extend_from_slice(&0.0_f32.to_le_bytes());
        }
        // Hostile face: claims i32::MAX index entries follow.
        bytes.extend_from_slice(&i32::MAX.to_le_bytes());
        let err = read_bytes("hostile".into(), &bytes).unwrap_err();
        assert!(
            matches!(err, PlyError::ListTooLarge { .. }),
            "expected ListTooLarge, got {err:?}"
        );
    }

    #[test]
    fn write_then_read_round_trip() {
        let mut m = Mesh::new("rt");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = vec![0, 1, 2];
        m.element_blocks = vec![blk];

        let mut buf: Vec<u8> = Vec::new();
        write_to(&m, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let m2 = read_str("rt".into(), &text).unwrap();
        assert_eq!(m2.nodes.len(), 3);
        assert_eq!(m2.element_blocks[0].connectivity, vec![0, 1, 2]);
    }

    /// Round-3 fix: a PLY face index list with a signed-int property
    /// type can carry negative values. Pre-fix the reader did `value
    /// as u32`, which saturates negatives to 0 and produces silent
    /// connectivity corruption (the face would still parse but point
    /// at vertex 0 instead of the intended one). Confirm the reader
    /// now refuses with InvalidIndex.
    #[test]
    fn rejects_negative_face_index() {
        // 4 vertices, 1 face — but the face references index -1 (which
        // in a signed-int list is "spec-legal but semantically wrong";
        // attackers craft this to poison the connectivity table).
        let ply = "\
ply
format ascii 1.0
element vertex 4
property float x
property float y
property float z
element face 1
property list uchar int vertex_indices
end_header
0 0 0
1 0 0
0 1 0
0 0 1
3 0 1 -1
";
        let err =
            read_str("neg-idx".into(), ply).expect_err("negative face index must be rejected");
        match err {
            PlyError::InvalidIndex { value } => {
                assert_eq!(value, -1.0);
            }
            other => panic!("expected InvalidIndex, got {other:?}"),
        }
    }

    /// R34 L1 (RED→GREEN, ASCII): a face citing a vertex index past
    /// the declared `element vertex` count passes the negative /
    /// non-finite / `> u32::MAX` checks but is still out of range for
    /// the mesh. Pre-fix it parsed into connectivity, then panicked in
    /// the shared consumers (`decimate` / `boolean`) which raw-index
    /// `positions[tri[k]]`. PLY declares vertices before faces, so the
    /// vertex count is known at face-parse time — reject there.
    #[test]
    fn ascii_rejects_face_index_past_vertex_count() {
        let ply = "\
ply
format ascii 1.0
element vertex 3
property float x
property float y
property float z
element face 1
property list uchar int vertex_indices
end_header
0 0 0
1 0 0
0 1 0
3 0 1 999
";
        let err = read_str("oob".into(), ply).expect_err("index past vertex count must error");
        match err {
            PlyError::InvalidIndex { value } => assert_eq!(value, 999.0),
            other => panic!("expected InvalidIndex, got {other:?}"),
        }
    }

    /// R34 L1 (RED→GREEN, binary): same hazard through the
    /// `binary_little_endian` `checked_index` path. 3 vertices, one
    /// face whose third index is 999.
    #[test]
    fn binary_rejects_face_index_past_vertex_count() {
        let header = b"ply\n\
format binary_little_endian 1.0\n\
element vertex 3\n\
property float x\n\
property float y\n\
property float z\n\
element face 1\n\
property list uchar int vertex_indices\n\
end_header\n";
        let mut bytes = header.to_vec();
        for _ in 0..3 {
            bytes.extend_from_slice(&0.0_f32.to_le_bytes());
            bytes.extend_from_slice(&0.0_f32.to_le_bytes());
            bytes.extend_from_slice(&0.0_f32.to_le_bytes());
        }
        bytes.push(3u8);
        for i in [0i32, 1, 999] {
            bytes.extend_from_slice(&i.to_le_bytes());
        }
        let err = read_bytes("oob-bin".into(), &bytes).expect_err("binary index past count");
        match err {
            PlyError::InvalidIndex { value } => assert_eq!(value, 999.0),
            other => panic!("expected InvalidIndex, got {other:?}"),
        }
    }

    #[test]
    fn checked_index_rejects_negative_nan_inf_overflow() {
        // A generous vertex count so these assertions exercise only the
        // finite / u32-range guard, not the R34 L1 vertex-count bound.
        let big = u32::MAX as usize + 1;
        assert!(matches!(
            checked_index(-1.0, big),
            Err(PlyError::InvalidIndex { .. })
        ));
        assert!(matches!(
            checked_index(f64::NAN, big),
            Err(PlyError::InvalidIndex { .. })
        ));
        assert!(matches!(
            checked_index(f64::INFINITY, big),
            Err(PlyError::InvalidIndex { .. })
        ));
        // u32::MAX + 1 overflows u32 but is finite as f64.
        assert!(matches!(
            checked_index((u32::MAX as f64) + 1.0, big),
            Err(PlyError::InvalidIndex { .. })
        ));
        // Valid edges (index strictly less than the vertex count).
        assert_eq!(checked_index(0.0, big).unwrap(), 0);
        assert_eq!(checked_index(u32::MAX as f64, big).unwrap(), u32::MAX);
        assert_eq!(checked_index(42.0, 43).unwrap(), 42);
    }

    /// R34 L1: the vertex-count bound — an index equal to or past the
    /// vertex count is rejected even though it is a valid finite u32.
    #[test]
    fn checked_index_rejects_index_at_or_past_vertex_count() {
        // 3 vertices → valid indices are 0, 1, 2.
        assert_eq!(checked_index(2.0, 3).unwrap(), 2);
        assert!(matches!(
            checked_index(3.0, 3),
            Err(PlyError::InvalidIndex { value }) if value == 3.0
        ));
        assert!(matches!(
            checked_index(999.0, 3),
            Err(PlyError::InvalidIndex { value }) if value == 999.0
        ));
        // Empty vertex array → every index is out of range.
        assert!(matches!(
            checked_index(0.0, 0),
            Err(PlyError::InvalidIndex { .. })
        ));
    }

    /// RED→GREEN (R29 A): `read_path_with_cap` refuses a file larger
    /// than its byte cap with `PlyError::TooLarge` *before* slurping it
    /// into RAM. Pre-fix `read_path` did an unbounded `read_to_end` with
    /// no size check, so a multi-GB file would OOM the import. We test
    /// the helper with a deliberately tiny 16-byte cap against a ~100-byte
    /// PLY rather than materialising a 1 GiB fixture.
    #[test]
    fn read_path_with_cap_rejects_oversize_file() {
        use std::io::Write;
        let tmp = std::env::temp_dir().join(format!(
            "valenx_mesh_ply_cap_red_{}.ply",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::File::create(&tmp).unwrap();
            // A minimal, well-formed single-vertex PLY (~90 bytes) — well
            // past the 16-byte cap, but a clean parse under the real cap.
            f.write_all(
                b"ply\nformat ascii 1.0\nelement vertex 1\nproperty float x\nproperty float y\nproperty float z\nend_header\n0 0 0\n",
            )
            .unwrap();
        }
        let res = read_path_with_cap(&tmp, 16);
        // Sanity: with the real (huge) cap the same file parses fine, so
        // the rejection is purely the cap, not a malformed fixture.
        let ok = read_path_with_cap(&tmp, MAX_PLY_FILE_BYTES);
        let _ = std::fs::remove_file(&tmp);
        match res {
            Err(PlyError::TooLarge { size, cap }) => {
                assert_eq!(cap, 16);
                assert!(size > 16, "file should be larger than the 16-byte cap");
            }
            other => panic!("expected PlyError::TooLarge, got {other:?}"),
        }
        assert!(ok.is_ok(), "file under the real cap must still parse");
    }
}
