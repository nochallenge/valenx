//! STL loader — the format the viewport uses to show imported
//! geometry before the BRep kernel is wired up.
//!
//! STL comes in two flavours:
//!
//! - **ASCII STL** — the `solid …` / `facet normal … / outer loop …`
//!   text form. Still common for tiny test fixtures and as a
//!   human-debuggable export target.
//! - **Binary STL** — 80-byte header, 4-byte little-endian triangle
//!   count, then `N * 50` bytes of `{normal:f32×3, v0:f32×3, v1:f32×3,
//!   v2:f32×3, attribute:u16}`.
//!
//! Both forms are normalized into a common [`TriangleMesh`] with
//! per-face normals and flat vertex lists. The loader auto-detects
//! which variant a file is by sniffing the first bytes.
//!
//! No `unsafe`, no `bytemuck` — the byte juggling is explicit and
//! goes through `u32::from_le_bytes` / `f32::from_le_bytes`.

use std::fs;
use std::io;
use std::path::Path;

use thiserror::Error;

/// A single triangle: three corner positions plus a face normal. All
/// values are `[x, y, z]` tuples in the STL file's native coordinate
/// system (mesh units undeclared — STL does not store them).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StlTriangle {
    pub normal: [f32; 3],
    pub vertices: [[f32; 3]; 3],
}

impl StlTriangle {
    /// Compute the face normal from the triangle's winding order
    /// (right-hand rule). Falls back to `[0, 0, 1]` if the triangle is
    /// degenerate.
    pub fn computed_normal(&self) -> [f32; 3] {
        let v0 = self.vertices[0];
        let v1 = self.vertices[1];
        let v2 = self.vertices[2];
        let a = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
        let b = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
        let n = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if len < 1e-20 {
            [0.0, 0.0, 1.0]
        } else {
            [n[0] / len, n[1] / len, n[2] / len]
        }
    }
}

/// The flavour the loader sniffed out of a file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StlFormat {
    Ascii,
    Binary,
}

/// A loaded STL file as a triangle soup. The viewport turns this
/// into GPU buffers; `valenx-mesh` can also promote it to a proper
/// `Mesh` with boundary groups inferred from connectivity.
#[derive(Clone, Debug, Default)]
pub struct TriangleMesh {
    pub format: Option<StlFormat>,
    pub name: Option<String>,
    pub triangles: Vec<StlTriangle>,
}

impl TriangleMesh {
    /// New, empty triangle mesh (no triangles, no format hint).
    pub fn new() -> Self {
        Self::default()
    }

    /// Total triangle count.
    pub fn triangle_count(&self) -> usize {
        self.triangles.len()
    }

    /// World-space axis-aligned bounding box over all vertex
    /// positions. `None` for an empty mesh.
    pub fn bounding_box(&self) -> Option<([f32; 3], [f32; 3])> {
        let first = self.triangles.first()?.vertices[0];
        let mut min = first;
        let mut max = first;
        for tri in &self.triangles {
            for v in &tri.vertices {
                for i in 0..3 {
                    if v[i] < min[i] {
                        min[i] = v[i];
                    }
                    if v[i] > max[i] {
                        max[i] = v[i];
                    }
                }
            }
        }
        Some((min, max))
    }
}

/// Errors the loader may emit.
#[derive(Debug, Error)]
pub enum StlError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),

    #[error("file is empty or too short to be an STL")]
    TooShort,

    #[error("binary STL triangle count {declared} mismatches body length {actual_tris} triangles")]
    BinaryTriangleCountMismatch { declared: u32, actual_tris: usize },

    #[error("ASCII STL parse error at line {line}: {reason}")]
    AsciiParse { line: usize, reason: String },

    /// Round-9 DoS hardening: the file's metadata reports a size
    /// larger than [`MAX_STL_FILE_BYTES`]. Pre-fix, `fs::read` would
    /// happily allocate a multi-GB buffer before any sanity check
    /// fired. 512 MiB is generous (legit STL exports of large
    /// product meshes routinely top 100 MiB) while still refusing
    /// `cat /dev/zero > big.stl` style denial of service.
    #[error("STL file too large: {size} bytes > {cap} cap (DoS guard)")]
    FileTooLarge { size: u64, cap: u64 },
}

/// Round-9 DoS cap for STL file reads. STL is the most common viz
/// payload and can legitimately reach hundreds of MiB for dense
/// product meshes; the cap is set high enough that real files load
/// while still refusing the multi-GB pathological case.
pub const MAX_STL_FILE_BYTES: u64 = 512 * 1024 * 1024;

/// Load an STL file by path — auto-detects ASCII vs binary.
pub fn load<P: AsRef<Path>>(path: P) -> Result<TriangleMesh, StlError> {
    // Round-9 DoS hardening: file-size cap before allocating.
    let size = fs::metadata(path.as_ref())?.len();
    if size > MAX_STL_FILE_BYTES {
        return Err(StlError::FileTooLarge {
            size,
            cap: MAX_STL_FILE_BYTES,
        });
    }
    let bytes = fs::read(path.as_ref())?;
    load_from_bytes(&bytes)
}

/// Load an STL from an in-memory byte slice.
pub fn load_from_bytes(bytes: &[u8]) -> Result<TriangleMesh, StlError> {
    if bytes.len() < 15 {
        return Err(StlError::TooShort);
    }
    if looks_like_ascii(bytes) {
        parse_ascii(bytes)
    } else {
        parse_binary(bytes)
    }
}

/// ASCII detection heuristic: the first 5 non-whitespace bytes match
/// `solid`, AND the body *also* contains the ASCII-only `facet normal`
/// token. Using both guards avoids the common trap where binary STL
/// files start with a 5-byte ASCII `solid` inside their 80-byte
/// header.
fn looks_like_ascii(bytes: &[u8]) -> bool {
    let start_ascii = bytes
        .iter()
        .copied()
        .skip_while(|b| b.is_ascii_whitespace())
        .take(5)
        .eq(b"solid".iter().copied());
    if !start_ascii {
        return false;
    }
    // Only scan the first ~2 KiB for the token to keep this cheap on
    // huge files.
    let scan_end = bytes.len().min(2048);
    bytes[..scan_end].windows(12).any(|w| w == b"facet normal")
}

fn parse_ascii(bytes: &[u8]) -> Result<TriangleMesh, StlError> {
    let text = std::str::from_utf8(bytes).map_err(|e| StlError::AsciiParse {
        line: 0,
        reason: format!("file is not UTF-8: {e}"),
    })?;

    let mut mesh = TriangleMesh::new();
    mesh.format = Some(StlFormat::Ascii);
    let mut current: Option<StlTriangle> = None;
    let mut vertex_idx: usize = 0;

    for (line_no, raw_line) in text.lines().enumerate() {
        let line = raw_line.trim();
        let line_no = line_no + 1; // 1-indexed for humans
        let mut it = line.split_whitespace();
        let first = match it.next() {
            Some(s) => s,
            None => continue,
        };

        match first {
            "solid" => {
                mesh.name = it.next().map(|s| s.to_string());
            }
            "facet" => {
                // facet normal nx ny nz
                let kw = it.next();
                if kw != Some("normal") {
                    return Err(StlError::AsciiParse {
                        line: line_no,
                        reason: "expected `facet normal`".into(),
                    });
                }
                let normal = parse_vec3(&mut it, line_no)?;
                current = Some(StlTriangle {
                    normal,
                    vertices: [[0.0; 3]; 3],
                });
                vertex_idx = 0;
            }
            "vertex" => {
                let v = parse_vec3(&mut it, line_no)?;
                let tri = current.as_mut().ok_or_else(|| StlError::AsciiParse {
                    line: line_no,
                    reason: "vertex outside a facet".into(),
                })?;
                if vertex_idx > 2 {
                    return Err(StlError::AsciiParse {
                        line: line_no,
                        reason: "more than 3 vertices in facet".into(),
                    });
                }
                tri.vertices[vertex_idx] = v;
                vertex_idx += 1;
            }
            "endfacet" => {
                let tri = current.take().ok_or_else(|| StlError::AsciiParse {
                    line: line_no,
                    reason: "endfacet without facet".into(),
                })?;
                if vertex_idx != 3 {
                    return Err(StlError::AsciiParse {
                        line: line_no,
                        reason: format!("facet has {vertex_idx} vertices, expected 3"),
                    });
                }
                mesh.triangles.push(tri);
            }
            "outer" | "endloop" | "endsolid" => { /* scaffolding */ }
            other => {
                // Spurious lines in well-formed files are rare; we
                // skip quietly rather than fail so third-party
                // exporters that add extras (comments, metadata) don't
                // break import.
                tracing::debug!(target: "valenx-viz", line = line_no, token = other, "unrecognised STL token, skipping");
            }
        }
    }

    if current.is_some() {
        return Err(StlError::AsciiParse {
            line: text.lines().count(),
            reason: "unterminated facet at end of file".into(),
        });
    }
    Ok(mesh)
}

fn parse_vec3<'a>(
    it: &mut impl Iterator<Item = &'a str>,
    line_no: usize,
) -> Result<[f32; 3], StlError> {
    let parse_one = |s: Option<&str>| -> Result<f32, StlError> {
        s.ok_or_else(|| StlError::AsciiParse {
            line: line_no,
            reason: "missing vec3 component".into(),
        })
        .and_then(|v| {
            v.parse::<f32>().map_err(|e| StlError::AsciiParse {
                line: line_no,
                reason: format!("bad float {v:?}: {e}"),
            })
        })
    };
    Ok([
        parse_one(it.next())?,
        parse_one(it.next())?,
        parse_one(it.next())?,
    ])
}

fn parse_binary(bytes: &[u8]) -> Result<TriangleMesh, StlError> {
    // 80-byte header + 4-byte triangle count + 50 bytes per triangle.
    if bytes.len() < 84 {
        return Err(StlError::TooShort);
    }
    let declared = u32::from_le_bytes(bytes[80..84].try_into().unwrap());
    let body = &bytes[84..];
    let per_triangle = 12 + 12 * 3 + 2; // normal + 3 vertices + attribute word
    let actual_tris = body.len() / per_triangle;
    if actual_tris != declared as usize {
        return Err(StlError::BinaryTriangleCountMismatch {
            declared,
            actual_tris,
        });
    }

    let mut mesh = TriangleMesh::new();
    mesh.format = Some(StlFormat::Binary);
    mesh.triangles.reserve_exact(actual_tris);
    for i in 0..actual_tris {
        let base = i * per_triangle;
        let t = StlTriangle {
            normal: read_vec3(&body[base..base + 12]),
            vertices: [
                read_vec3(&body[base + 12..base + 24]),
                read_vec3(&body[base + 24..base + 36]),
                read_vec3(&body[base + 36..base + 48]),
            ],
        };
        mesh.triangles.push(t);
        // attribute byte count at [base + 48..base + 50] ignored.
    }
    Ok(mesh)
}

fn read_vec3(slice: &[u8]) -> [f32; 3] {
    [
        f32::from_le_bytes(slice[0..4].try_into().unwrap()),
        f32::from_le_bytes(slice[4..8].try_into().unwrap()),
        f32::from_le_bytes(slice[8..12].try_into().unwrap()),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    const ASCII_CUBE: &str = "solid cube\n\
facet normal 0 0 -1\n  outer loop\n    vertex 0 0 0\n    vertex 1 0 0\n    vertex 1 1 0\n  endloop\nendfacet\n\
facet normal 0 0 -1\n  outer loop\n    vertex 0 0 0\n    vertex 1 1 0\n    vertex 0 1 0\n  endloop\nendfacet\n\
endsolid cube\n";

    #[test]
    fn parses_ascii_cube() {
        let mesh = load_from_bytes(ASCII_CUBE.as_bytes()).expect("parse");
        assert_eq!(mesh.format, Some(StlFormat::Ascii));
        assert_eq!(mesh.name.as_deref(), Some("cube"));
        assert_eq!(mesh.triangles.len(), 2);
        assert_eq!(mesh.triangles[0].normal, [0.0, 0.0, -1.0]);
        assert_eq!(mesh.triangles[0].vertices[2], [1.0, 1.0, 0.0]);
    }

    #[test]
    fn bounding_box_matches() {
        let mesh = load_from_bytes(ASCII_CUBE.as_bytes()).expect("parse");
        let (min, max) = mesh.bounding_box().expect("non-empty");
        assert_eq!(min, [0.0, 0.0, 0.0]);
        assert_eq!(max, [1.0, 1.0, 0.0]);
    }

    #[test]
    fn computed_normal_is_right_handed() {
        let t = StlTriangle {
            normal: [0.0; 3],
            vertices: [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        };
        let n = t.computed_normal();
        // Right-hand rule: vertices wound CCW in XY plane → normal +Z.
        assert!((n[0]).abs() < 1e-6);
        assert!((n[1]).abs() < 1e-6);
        assert!((n[2] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn binary_stl_roundtrip() {
        // Hand-construct a 2-triangle binary STL.
        let mut bytes: Vec<u8> = Vec::new();
        bytes.extend_from_slice(&[0u8; 80]); // header
        bytes.extend_from_slice(&2u32.to_le_bytes()); // triangle count

        let tri = |n: [f32; 3], v: [[f32; 3]; 3]| {
            let mut b = Vec::<u8>::new();
            for c in n {
                b.extend_from_slice(&c.to_le_bytes());
            }
            for p in v {
                for c in p {
                    b.extend_from_slice(&c.to_le_bytes());
                }
            }
            b.extend_from_slice(&0u16.to_le_bytes());
            b
        };
        bytes.extend(tri(
            [0.0, 0.0, 1.0],
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        ));
        bytes.extend(tri(
            [0.0, 0.0, 1.0],
            [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        ));

        let mesh = load_from_bytes(&bytes).expect("parse binary");
        assert_eq!(mesh.format, Some(StlFormat::Binary));
        assert_eq!(mesh.triangle_count(), 2);
        assert_eq!(mesh.triangles[0].vertices[1], [1.0, 0.0, 0.0]);
    }

    #[test]
    fn rejects_truncated_binary() {
        // header + triangle count say 2, but body only holds 0.
        let mut bytes = vec![0u8; 80];
        bytes.extend_from_slice(&2u32.to_le_bytes());
        let err = load_from_bytes(&bytes).unwrap_err();
        assert!(matches!(
            err,
            StlError::BinaryTriangleCountMismatch {
                declared: 2,
                actual_tris: 0
            }
        ));
    }

    #[test]
    fn rejects_too_short() {
        let err = load_from_bytes(b"hi").unwrap_err();
        assert!(matches!(err, StlError::TooShort));
    }

    /// Round-9 RED→GREEN: `load` must refuse a file larger than
    /// `MAX_STL_FILE_BYTES` before allocating. Use a sparse file
    /// so the test doesn't actually write 512 MiB.
    #[test]
    fn load_rejects_file_above_max_stl_file_bytes() {
        use std::io::{Seek, SeekFrom, Write};
        let tmp = std::env::temp_dir().join(format!(
            "valenx-stl-toobig-{}.stl",
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
            f.seek(SeekFrom::Start(super::MAX_STL_FILE_BYTES + 1))
                .unwrap();
            f.write_all(b"x").unwrap();
        }
        let err = load(&tmp).expect_err("must reject oversize STL");
        match err {
            StlError::FileTooLarge { size, cap } => {
                assert!(size > cap, "size={size} cap={cap}");
                assert_eq!(cap, super::MAX_STL_FILE_BYTES);
            }
            other => panic!("expected FileTooLarge, got: {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }
}
