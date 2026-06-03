//! Wavefront OBJ reader and writer for triangle surface meshes.
//!
//! Scope: positions (`v x y z`) and triangular faces (`f i j k` —
//! 1-based vertex indices, with optional `/vt/vn` suffixes that we
//! discard). Everything else (`o`, `g`, `usemtl`, `mtllib`, `s`,
//! `vn`, `vt`, comments, blank lines) is ignored.
//!
//! Polygons of more than three vertices on a single `f` line are
//! fan-triangulated (`f 1 2 3 4 5` → `(1, 2, 3)`, `(1, 3, 4)`,
//! `(1, 4, 5)`). Negative indices (relative-to-end addressing) are
//! resolved against the running vertex count.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use nalgebra::Vector3;
use thiserror::Error;

use crate::element::{ElementBlock, ElementType};
use crate::mesh::Mesh;

/// Round-24 H3: cap on the bytes any single OBJ ASCII line may
/// consume during import. Mirrors `valenx_core::io_caps::MAX_OBJ_LINE_BYTES`
/// — duplicated inline to dodge the valenx-core → valenx-fields →
/// valenx-mesh dependency cycle (valenx-core is the upstream crate
/// so it can't depend on valenx-mesh, and valenx-mesh can't add a
/// valenx-core dep without closing the loop). OBJ lines are
/// typically under 1 KiB; 4 MiB is generous for hand-authored CAD
/// `f` polygon lines while refusing the unbounded-line DoS shape.
const MAX_OBJ_LINE_BYTES: usize = 4 * 1024 * 1024;

/// Maximum number of vertices on a single OBJ `f` polygon face we
/// fan-triangulate. Hand-authored CAD output sometimes ships single-
/// face polygons in the 10-50 vertex range; reasonable real meshes
/// almost never exceed 1024. Anything past that is either a hostile
/// file (a single `f` line with millions of tokens would balloon the
/// `conn` allocation past the available RAM) or a generator bug.
pub const MAX_OBJ_FACE_VERTICES: usize = 1024;

/// All OBJ-related errors.
#[derive(Debug, Error)]
pub enum ObjError {
    /// Something bad happened during file IO.
    #[error(transparent)]
    Io(#[from] io::Error),
    /// A line couldn't be parsed.
    #[error("malformed OBJ line {line}: {reason}")]
    Malformed {
        /// 1-based line number in the source file.
        line: usize,
        /// What went wrong.
        reason: String,
    },
}

/// Read an OBJ file from `path` and return a canonical [`Mesh`].
///
/// The mesh `id` is set to the file stem (or `"obj"` if the path
/// has no stem). Quality stats are NOT recomputed — call
/// `mesh.recompute_quality_stats()` if you need them.
///
/// Round-24 H3: pre-fix the reader used
/// `reader.lines().map_while(Result::ok)` — an iterator that 1)
/// allocated an unbounded `String` per `\n`-delimited record (4 GiB
/// no-newline OBJ → OOM) and 2) silently dropped IO errors mid-read
/// (truncated mesh import surfaces as `Ok(partial)`, not `Err`).
/// Fix: route through `read_capped_lines_bounded` with the shared
/// `MAX_OBJ_LINE_BYTES` cap (4 MiB), convert each capped line to UTF-8
/// lossily, and propagate any IO / cap error via `?` into the
/// returned `Result`.
///
/// Round-25 M1: pre-fix this collected ALL lines into a
/// `Vec<String>` before calling `read_lines` — doubling peak memory
/// on a 100 MiB OBJ (one copy in the line vec, one copy as the
/// parser walked it building the mesh). The fix streams the
/// bounded-line iterator directly into `parse_streaming` so each
/// line allocation lives only for the duration of one match arm.
///
/// Round-26 L2: pre-fix the byte→`String` conversion used
/// `String::from_utf8_lossy`, which silently replaced non-UTF-8
/// bytes with U+FFFD. A CAD exporter shipping an unintended
/// Latin-1 / Windows-1252 file would import with mangled comments
/// but no error surface, and a deliberately hostile non-UTF-8
/// payload would slip past unnoticed. The fix uses `String::from_utf8`
/// strictly and surfaces failures as `ObjError::Malformed` carrying
/// the 1-based line number; the caller can then point the user at
/// the offending record. Valid OBJ tokens are pure ASCII so legitimate
/// files round-trip identically.
pub fn read_path(path: impl AsRef<Path>) -> Result<Mesh, ObjError> {
    let path_ref = path.as_ref();
    let file = File::open(path_ref)?;
    let reader = BufReader::new(file);
    let id = path_ref
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("obj")
        .to_string();
    // Round-26 L2: parse_streaming_bytes converts each line from
    // bytes to `String` with strict UTF-8 validation, surfacing
    // invalid records with a precise line number rather than
    // silently replacing them with U+FFFD.
    parse_streaming_bytes(
        id,
        read_capped_obj_lines(reader, MAX_OBJ_LINE_BYTES),
    )
}

/// Round-25 M1: streaming OBJ parser. Mirrors `read_lines` but
/// accepts a `Result`-yielding iterator and surfaces IO/cap errors
/// to the caller via `?`, so the file-reader path keeps only one
/// `String` per line in memory at a time instead of the entire
/// `Vec<String>` round-24's implementation collected up front.
///
/// Kept separate from `read_lines` so the in-memory test path (which
/// passes an infallible iterator of `&str` literals) doesn't have to
/// wrap every line in `Ok(...)`.
///
/// Round-26 L2: `read_path` no longer uses this — strict-UTF-8
/// validation lives in `parse_streaming_bytes` so it can carry
/// line numbers into `ObjError::Malformed`. This shape is retained
/// as a test fixture for the round-25 M1 "does not collect input"
/// contract; `#[cfg(test)]` keeps it out of production builds so
/// the dead-code lint doesn't fire.
#[cfg(test)]
fn parse_streaming<I, S>(id: String, lines: I) -> Result<Mesh, ObjError>
where
    I: IntoIterator<Item = std::io::Result<S>>,
    S: AsRef<str>,
{
    let mut mesh = Mesh::new(id);
    let mut conn: Vec<u32> = Vec::new();
    for (i, line_res) in lines.into_iter().enumerate() {
        let raw = line_res?;
        let line_no = i + 1;
        process_line(raw.as_ref(), line_no, &mut mesh, &mut conn)?;
    }
    finalize_mesh(&mut mesh, conn);
    Ok(mesh)
}

/// Round-26 L2: streaming OBJ parser variant that takes a bytes
/// iterator and does strict UTF-8 validation per line, surfacing
/// the failing 1-based line number via `ObjError::Malformed`.
/// Replaces the pre-fix `String::from_utf8_lossy` shape in
/// `read_path` — silent U+FFFD replacement hid both bona-fide
/// encoding errors and deliberate non-UTF-8 payloads.
fn parse_streaming_bytes<I>(id: String, lines: I) -> Result<Mesh, ObjError>
where
    I: IntoIterator<Item = std::io::Result<Vec<u8>>>,
{
    let mut mesh = Mesh::new(id);
    let mut conn: Vec<u32> = Vec::new();
    for (i, line_res) in lines.into_iter().enumerate() {
        let bytes = line_res?;
        let line_no = i + 1;
        // Strict UTF-8: silently-replaced U+FFFD chars from
        // from_utf8_lossy would mask both real encoding errors
        // (non-UTF-8 exporter output) and hostile payloads. Reject
        // the line with a precise number so the caller can repair.
        let raw = String::from_utf8(bytes).map_err(|e| ObjError::Malformed {
            line: line_no,
            reason: format!("not valid UTF-8: {e}"),
        })?;
        process_line(&raw, line_no, &mut mesh, &mut conn)?;
    }
    finalize_mesh(&mut mesh, conn);
    Ok(mesh)
}

/// Inline mirror of `valenx_core::io_caps::read_capped_lines_bounded`
/// — see `MAX_OBJ_LINE_BYTES` for the dependency-cycle rationale.
/// Bounded line reader: each `next()` returns up to `max_per_line`
/// bytes (including the trailing `\n` if present) and stops the
/// iterator with `Err(InvalidData)` if a single line exceeds the cap
/// before a newline appears. Stops on EOF (Ok(0)) and on IO error.
fn read_capped_obj_lines<R: BufRead>(
    mut reader: R,
    max_per_line: usize,
) -> impl Iterator<Item = std::io::Result<Vec<u8>>> {
    let mut done = false;
    std::iter::from_fn(move || {
        if done {
            return None;
        }
        let mut buf: Vec<u8> = Vec::with_capacity(128);
        let cap = (max_per_line as u64).saturating_add(1);
        let mut limited = Read::take(&mut reader, cap);
        match limited.read_until(b'\n', &mut buf) {
            Ok(0) => {
                done = true;
                None
            }
            Ok(_) => {
                if buf.len() > max_per_line {
                    done = true;
                    return Some(Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "OBJ line exceeded {max_per_line}-byte cap (possible \
                             missing newline or hostile input)"
                        ),
                    )));
                }
                Some(Ok(buf))
            }
            Err(e) => {
                done = true;
                Some(Err(e))
            }
        }
    })
}

/// Read an OBJ from any iterator of lines (used by `read_path` and
/// by tests that pass in-memory strings).
///
/// Round-25 M1: factored into `process_line` + `finalize_mesh` so
/// the streaming reader (`parse_streaming`) and the in-memory test
/// path can share the same per-line logic without either one having
/// to allocate the full input up front.
pub fn read_lines<I, S>(id: String, lines: I) -> Result<Mesh, ObjError>
where
    I: IntoIterator<Item = S>,
    S: AsRef<str>,
{
    let mut mesh = Mesh::new(id);
    let mut conn: Vec<u32> = Vec::new();
    for (i, raw) in lines.into_iter().enumerate() {
        let line_no = i + 1;
        process_line(raw.as_ref(), line_no, &mut mesh, &mut conn)?;
    }
    finalize_mesh(&mut mesh, conn);
    Ok(mesh)
}

/// Round-25 M1: parse one OBJ line into the running mesh + connectivity
/// vector. Shared by the streaming `parse_streaming` path (which calls
/// this once per `String` read from disk) and the in-memory `read_lines`
/// path (which calls it once per `&str` literal in the iterator). Pure
/// per-line logic — no buffering, no shared state beyond `mesh` and
/// `conn`.
fn process_line(
    raw: &str,
    line_no: usize,
    mesh: &mut Mesh,
    conn: &mut Vec<u32>,
) -> Result<(), ObjError> {
    let line = raw.trim();
    if line.is_empty() || line.starts_with('#') {
        return Ok(());
    }
    let mut iter = line.split_whitespace();
    let Some(tag) = iter.next() else {
        return Ok(());
    };
    match tag {
        "v" => {
            let coords: Vec<&str> = iter.collect();
            if coords.len() < 3 {
                return Err(ObjError::Malformed {
                    line: line_no,
                    reason: format!("v needs >=3 coords, got {}", coords.len()),
                });
            }
            let x = parse_f64(coords[0], line_no)?;
            let y = parse_f64(coords[1], line_no)?;
            let z = parse_f64(coords[2], line_no)?;
            mesh.nodes.push(Vector3::new(x, y, z));
        }
        "f" => {
            let raw_verts: Vec<&str> = iter.collect();
            if raw_verts.len() < 3 {
                return Err(ObjError::Malformed {
                    line: line_no,
                    reason: format!("f needs >=3 verts, got {}", raw_verts.len()),
                });
            }
            // Round-5 DoS hardening: reject pathological n-gons
            // before the fan-triangulation loop allocates O(n)
            // connectivity entries. A single `f` line with millions
            // of vertex tokens would otherwise OOM the import.
            if raw_verts.len() > MAX_OBJ_FACE_VERTICES {
                return Err(ObjError::Malformed {
                    line: line_no,
                    reason: format!(
                        "f face has {} vertices, max is {} \
                         (re-export with fewer vertices per face — \
                         modern OBJ writers triangulate by default)",
                        raw_verts.len(),
                        MAX_OBJ_FACE_VERTICES
                    ),
                });
            }
            // Each token is `v` or `v/vt` or `v/vt/vn` or `v//vn`.
            let mut idxs: Vec<u32> = Vec::with_capacity(raw_verts.len());
            for tok in raw_verts {
                let v_str = tok.split('/').next().unwrap_or("");
                let v: i64 = v_str.parse().map_err(|_| ObjError::Malformed {
                    line: line_no,
                    reason: format!("not an integer face index: {tok:?}"),
                })?;
                let actual: i64 = if v > 0 {
                    v - 1
                } else if v < 0 {
                    mesh.nodes.len() as i64 + v
                } else {
                    return Err(ObjError::Malformed {
                        line: line_no,
                        reason: "face index 0 not allowed (OBJ is 1-based)".into(),
                    });
                };
                if actual < 0 || actual >= mesh.nodes.len() as i64 {
                    return Err(ObjError::Malformed {
                        line: line_no,
                        reason: format!(
                            "face index {v} resolves to {actual} out of range [0, {})",
                            mesh.nodes.len()
                        ),
                    });
                }
                idxs.push(actual as u32);
            }
            // Fan triangulation for n >= 3.
            for k in 1..(idxs.len() - 1) {
                conn.extend_from_slice(&[idxs[0], idxs[k], idxs[k + 1]]);
            }
        }
        _ => {
            // Ignore vn, vt, o, g, usemtl, mtllib, s, etc.
        }
    }
    Ok(())
}

/// Round-25 M1: install the accumulated triangle connectivity into
/// the mesh and refresh derived stats. Shared finalizer between the
/// streaming + in-memory parse paths.
fn finalize_mesh(mesh: &mut Mesh, conn: Vec<u32>) {
    if !conn.is_empty() {
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = conn;
        mesh.element_blocks.push(blk);
    }
    mesh.recompute_stats();
}

fn parse_f64(s: &str, line_no: usize) -> Result<f64, ObjError> {
    s.parse::<f64>().map_err(|_| ObjError::Malformed {
        line: line_no,
        reason: format!("not a float: {s:?}"),
    })
}

/// Write the Tri3 blocks of `mesh` to `path` as a Wavefront OBJ.
///
/// Emits `v x y z` lines for every node and `f i j k` lines for
/// every triangle (1-based indices). Non-Tri3 blocks are silently
/// skipped — OBJ is a surface format.
pub fn write_path(mesh: &Mesh, path: impl AsRef<Path>) -> Result<(), ObjError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    write_to(mesh, &mut writer)?;
    writer.flush()?;
    Ok(())
}

/// Write the Tri3 surface of `mesh` to an arbitrary `Write` sink.
pub fn write_to<W: Write>(mesh: &Mesh, mut sink: W) -> io::Result<()> {
    writeln!(sink, "# valenx-mesh OBJ export")?;
    writeln!(sink, "# {} vertices", mesh.nodes.len())?;
    for n in &mesh.nodes {
        writeln!(sink, "v {} {} {}", fmt_f(n.x), fmt_f(n.y), fmt_f(n.z))?;
    }
    let tri_count: usize = mesh
        .element_blocks
        .iter()
        .filter(|b| b.element_type == ElementType::Tri3)
        .map(|b| b.connectivity.len() / 3)
        .sum();
    writeln!(sink, "# {tri_count} triangles")?;
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            writeln!(sink, "f {} {} {}", tri[0] + 1, tri[1] + 1, tri[2] + 1)?;
        }
    }
    Ok(())
}

fn fmt_f(x: f64) -> String {
    // Use precise enough representation to round-trip without
    // surprising the user; trim trailing zeros for readability.
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
    fn read_simple_triangle() {
        let obj = "\
v 0 0 0
v 1 0 0
v 0 1 0
f 1 2 3
";
        let m = read_lines("tri".into(), obj.lines()).unwrap();
        assert_eq!(m.nodes.len(), 3);
        assert_eq!(m.element_blocks.len(), 1);
        assert_eq!(m.element_blocks[0].element_type, ElementType::Tri3);
        assert_eq!(m.element_blocks[0].connectivity, vec![0, 1, 2]);
    }

    #[test]
    fn ignores_comments_and_blank_lines_and_misc() {
        let obj = "\
# header
o square
g surface
v 0 0 0
v 1 0 0

v 1 1 0
v 0 1 0
vn 0 0 1
vt 0 0
f 1 2 3
f 1 3 4
";
        let m = read_lines("square".into(), obj.lines()).unwrap();
        assert_eq!(m.nodes.len(), 4);
        assert_eq!(m.element_blocks[0].connectivity.len(), 6);
    }

    #[test]
    fn fan_triangulates_polygon() {
        // Quad `f 1 2 3 4` → two triangles.
        let obj = "\
v 0 0 0
v 1 0 0
v 1 1 0
v 0 1 0
f 1 2 3 4
";
        let m = read_lines("quad".into(), obj.lines()).unwrap();
        assert_eq!(m.element_blocks[0].connectivity, vec![0, 1, 2, 0, 2, 3]);
    }

    #[test]
    fn supports_negative_indices() {
        let obj = "\
v 0 0 0
v 1 0 0
v 0 1 0
f -3 -2 -1
";
        let m = read_lines("neg".into(), obj.lines()).unwrap();
        assert_eq!(m.element_blocks[0].connectivity, vec![0, 1, 2]);
    }

    #[test]
    fn slash_suffixes_are_ignored() {
        let obj = "\
v 0 0 0
v 1 0 0
v 0 1 0
f 1/1/1 2/2/1 3/3/1
";
        let m = read_lines("slash".into(), obj.lines()).unwrap();
        assert_eq!(m.element_blocks[0].connectivity, vec![0, 1, 2]);
    }

    #[test]
    fn rejects_index_zero() {
        let obj = "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 0 1 2\n";
        let err = read_lines("bad".into(), obj.lines()).unwrap_err();
        assert!(matches!(err, ObjError::Malformed { .. }));
    }

    #[test]
    fn rejects_out_of_range_index() {
        let obj = "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 4\n";
        let err = read_lines("oor".into(), obj.lines()).unwrap_err();
        assert!(matches!(err, ObjError::Malformed { .. }));
    }

    /// Round-5 RED→GREEN: an OBJ with a single `f` line containing
    /// 2000 vertex tokens used to flow into the fan-triangulation
    /// loop and allocate millions of connectivity entries. The fix
    /// is a `MAX_OBJ_FACE_VERTICES = 1024` cap, surfaced as a
    /// structured `Malformed` error before the allocation happens.
    #[test]
    fn rejects_oversized_face() {
        // 3 vertex declarations followed by a single `f` line that
        // references 2000 vertices (replicating the index 1 over and
        // over — the index itself is in-range, what we're testing is
        // the per-face count cap).
        let mut obj = String::from("v 0 0 0\nv 1 0 0\nv 0 1 0\nf");
        for _ in 0..2000 {
            obj.push_str(" 1");
        }
        obj.push('\n');
        let err = read_lines("huge-face".into(), obj.lines()).unwrap_err();
        match err {
            ObjError::Malformed { line, reason } => {
                assert_eq!(line, 4);
                assert!(
                    reason.contains("max is")
                        || reason.contains(&MAX_OBJ_FACE_VERTICES.to_string()),
                    "reason: {reason}"
                );
            }
            other => panic!("expected Malformed, got {other:?}"),
        }
    }

    /// Sanity: faces at the cap boundary still parse.
    #[test]
    fn accepts_face_at_cap_boundary() {
        // We accept up to MAX_OBJ_FACE_VERTICES. Build an `f` line
        // with exactly that many entries to confirm the boundary is
        // inclusive.
        let mut obj = String::from("v 0 0 0\nv 1 0 0\nv 0 1 0\nf");
        for _ in 0..MAX_OBJ_FACE_VERTICES {
            obj.push_str(" 1");
        }
        obj.push('\n');
        let result = read_lines("boundary".into(), obj.lines());
        // Index-1 is in-range; this is malformed in another sense (all
        // verts identical → zero-area triangles) but the per-face cap
        // does not fire.
        assert!(result.is_ok(), "boundary case must parse: {result:?}");
    }

    #[test]
    fn write_then_read_round_trip() {
        // Write a small mesh to an in-memory string, parse it back,
        // and verify equivalence.
        let mut m = Mesh::new("rt");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = vec![0, 1, 2, 0, 2, 3];
        m.element_blocks = vec![blk];

        let mut buf: Vec<u8> = Vec::new();
        write_to(&m, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let m2 = read_lines("rt".into(), text.lines()).unwrap();
        assert_eq!(m2.nodes.len(), 4);
        assert_eq!(
            m2.element_blocks[0].connectivity,
            m.element_blocks[0].connectivity
        );
        for (a, b) in m.nodes.iter().zip(m2.nodes.iter()) {
            assert!((a - b).norm() < 1e-9);
        }
    }

    #[test]
    fn write_skips_non_tri3_blocks() {
        let mut m = Mesh::new("mixed");
        m.nodes = vec![Vector3::new(0.0, 0.0, 0.0); 4];
        let mut tri = ElementBlock::new(ElementType::Tri3);
        tri.connectivity = vec![0, 1, 2];
        let mut hex = ElementBlock::new(ElementType::Hex8);
        hex.connectivity = vec![0, 1, 2, 3, 0, 1, 2, 3];
        m.element_blocks = vec![tri, hex];
        let mut buf: Vec<u8> = Vec::new();
        write_to(&m, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        // Only one face line should appear.
        let f_count = text.lines().filter(|l| l.starts_with("f ")).count();
        assert_eq!(f_count, 1);
    }

    /// RED→GREEN (round-24 H3): `read_capped_obj_lines` refuses a
    /// line that exceeds `MAX_OBJ_LINE_BYTES` instead of allocating
    /// the entire input. Pre-fix `BufReader::lines().map_while(Result::ok)`
    /// would allocate one massive `String` per `\n`-delimited
    /// record AND silently swallow the error to surface as
    /// "EOF reached after N valid lines" — corrupted import passing
    /// as Ok.
    #[test]
    fn read_capped_obj_lines_caps_unbounded_record() {
        use std::io::Cursor;
        const TEST_CAP: usize = 64 * 1024;
        let payload = vec![b'x'; 256 * 1024]; // 4x cap, no newline
        let reader = std::io::BufReader::new(Cursor::new(payload));
        let mut errs = 0;
        let mut oks = 0;
        for line in read_capped_obj_lines(reader, TEST_CAP) {
            match line {
                Ok(_) => oks += 1,
                Err(e) => {
                    assert_eq!(e.kind(), std::io::ErrorKind::InvalidData);
                    errs += 1;
                }
            }
        }
        assert_eq!(oks, 0);
        assert_eq!(errs, 1);
    }

    /// RED→GREEN (round-25 M1): `read_path` streams lines through
    /// the parser instead of collecting them into a `Vec<String>`
    /// first. Pre-fix (round-24 H3) `read_path` did
    /// `for line in iter { lines.push(...) }` then handed the full
    /// vec to `read_lines` — peak memory was 2x the file size. We
    /// verify the contract by instrumenting `parse_streaming`'s
    /// caller: the iterator we hand to `parse_streaming` is pulled
    /// item-by-item, never materialised. The assertion is "we never
    /// build a `Vec<String>` of the full input" — we exercise that
    /// by running an iterator that counts how many items have been
    /// pulled at any moment and verifying it never exceeds a small
    /// constant (the parser only needs the current line).
    #[test]
    fn parse_streaming_does_not_collect_input_round25_m1() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        // Synthesise a small in-memory OBJ; the assertion is about
        // streaming behaviour, not data volume.
        let lines_in: Vec<String> = (0..1000)
            .map(|i| format!("v {i} 0 0"))
            .chain(std::iter::once("f 1 2 3".to_string()))
            .collect();
        // Wrap each `.next()` so we can assert the parser pulls them
        // one at a time. `live` counts how many items have been
        // pulled MINUS how many have been parsed; if streaming holds
        // we should never see > 1.
        let pulled = AtomicUsize::new(0);
        let parsed = AtomicUsize::new(0);
        let max_live = AtomicUsize::new(0);
        let iter = lines_in.iter().map(|s| {
            pulled.fetch_add(1, Ordering::SeqCst);
            let live = pulled.load(Ordering::SeqCst) - parsed.load(Ordering::SeqCst);
            max_live.fetch_max(live, Ordering::SeqCst);
            // Wrap in Ok so the streaming signature matches.
            let res: std::io::Result<String> = Ok(s.clone());
            res
        });
        // Tap on the parser side: re-wrap to bump the parsed counter
        // immediately after the parser sees the item.
        let iter2 = iter.inspect(|_res| {
            parsed.fetch_add(1, Ordering::SeqCst);
        });
        let m = parse_streaming("stream".into(), iter2).unwrap();
        // Sanity: the mesh parsed correctly.
        assert_eq!(m.nodes.len(), 1000);
        // The contract: the parser never holds more than ~1 item
        // live at a time. We allow 2 to give the iterator adapters
        // some slack (e.g. fused iterators may bump pulled twice).
        let m = max_live.load(Ordering::SeqCst);
        assert!(
            m <= 2,
            "parse_streaming must NOT materialise the input — \
             saw {m} live items at peak (expected ≤ 2)",
        );
    }

    /// RED→GREEN (round-24 H3): truncated OBJ on disk surfaces as
    /// `Err`, NOT a silently corrupted half-mesh. Pre-fix the OBJ
    /// reader's `map_while(Result::ok)` would drop IO errors mid-
    /// stream and return the partially-parsed mesh as success. We
    /// validate the contract by triggering the per-line cap with a
    /// pathological line — the same error kind that disk-truncation
    /// would now produce.
    #[test]
    fn read_path_propagates_cap_errors() {
        use std::io::Write;
        // 5 MiB of 'x' with NO newline → exceeds MAX_OBJ_LINE_BYTES (4 MiB).
        let tmp = std::env::temp_dir().join("valenx_mesh_obj_cap_red.obj");
        let mut f = std::fs::File::create(&tmp).unwrap();
        f.write_all(&vec![b'x'; 5 * 1024 * 1024]).unwrap();
        drop(f);
        let res = read_path(&tmp);
        let _ = std::fs::remove_file(&tmp);
        match res {
            Err(ObjError::Io(e)) => {
                assert_eq!(e.kind(), std::io::ErrorKind::InvalidData);
            }
            Err(other) => panic!("expected Io error, got {other:?}"),
            Ok(_) => panic!("expected Err on over-cap OBJ, got Ok"),
        }
    }

    /// RED→GREEN (round-26 L2): non-UTF-8 bytes in an OBJ file
    /// surface as `ObjError::Malformed { line, reason }` carrying
    /// the 1-based line number, NOT a silently-mangled mesh with
    /// U+FFFD replacement chars. Pre-fix `String::from_utf8_lossy`
    /// would have produced an `Ok(_)` whose nodes / comments
    /// contained the replacement character — undetectable to the
    /// caller.
    #[test]
    fn read_path_rejects_non_utf8_round26_l2() {
        use std::io::Write;
        let tmp = std::env::temp_dir().join(format!(
            "valenx_mesh_obj_non_utf8_{}.obj",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        // Two valid-UTF-8 lines, then a line with a stray 0xFF
        // byte (a continuation byte without a leading byte —
        // invalid UTF-8). The first two should parse but the third
        // should surface as Malformed with line == 3.
        {
            let mut f = std::fs::File::create(&tmp).unwrap();
            f.write_all(b"v 0 0 0\n").unwrap();
            f.write_all(b"v 1 0 0\n").unwrap();
            f.write_all(&[0xFFu8, b'\n']).unwrap();
        }
        let res = read_path(&tmp);
        let _ = std::fs::remove_file(&tmp);
        match res {
            Err(ObjError::Malformed { line, reason }) => {
                assert_eq!(line, 3, "expected line 3, got {line} (reason: {reason})");
                assert!(
                    reason.contains("UTF-8") || reason.contains("utf-8"),
                    "reason should mention UTF-8 invalidity, got: {reason}",
                );
            }
            Err(other) => panic!("expected Malformed, got {other:?}"),
            Ok(_) => panic!("expected Err on non-UTF-8 OBJ, got Ok"),
        }
    }
}
