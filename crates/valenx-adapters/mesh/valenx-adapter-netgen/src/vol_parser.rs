//! A narrow parser for Netgen's `.vol` ASCII format — enough to
//! promote a Netgen mesh output into a canonical
//! [`valenx_mesh::Mesh`] for the linear element types Valenx
//! consumes today.
//!
//! ## Format outline (the bits we honour)
//!
//! ```text
//! mesh3d
//! dimension
//! 3
//!
//! geomtype
//! 0
//!
//! surfaceelements
//! <N>
//!  surfnr bcnr domin domout np p1 p2 p3 [p4]      # np = 3 (Tri3) or 4 (Quad4)
//!  ... (N lines)
//!
//! volumeelements
//! <N>
//!  matnr np p1 p2 p3 p4 [p5..p8]                  # np = 4..8
//!  ... (N lines)
//!
//! points
//! <N>
//!  x y z
//!  ... (N lines)
//! ```
//!
//! Sections may appear in any order. Indices in element blocks are
//! 1-based in the file; we rewrite to a dense 0-based scheme so the
//! mesh round-trips cleanly through the rest of the workspace.
//!
//! Sections we don't yet consume (`identifications`, `materials`,
//! `bcnames`, `points_3d_periodic`, …) are tolerated silently — we
//! recognise the section header, read the count, and skip the listed
//! number of lines before resuming.
//!
//! The parser is intentionally first-party: pulling a Netgen-specific
//! crate for I/O would widen the dependency footprint for a stable,
//! well-documented text format.

use std::path::Path;

use nalgebra::Vector3;
use thiserror::Error;

use valenx_mesh::{ElementBlock, ElementType, Mesh};

/// Errors the parser may report.
#[derive(Debug, Error)]
pub enum VolError {
    /// I/O failed reading the file.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// The .vol text is structurally invalid.
    #[error("malformed .vol at line {line}: {reason}")]
    Malformed { line: usize, reason: String },

    /// File starts with `mesh2d`/`mesh3d` etc. that we don't support.
    #[error("unsupported .vol header `{header}`; parser handles `mesh3d` only")]
    UnsupportedHeader { header: String },

    /// A volume / surface element block uses a node count we don't
    /// have a canonical [`ElementType`] for (e.g. `np = 10` for
    /// quadratic tetrahedra).
    #[error("unsupported element node count np={np} in section `{section}`")]
    UnsupportedElement { section: &'static str, np: u32 },
}

/// Parse a `.vol` file into a canonical [`Mesh`].
///
/// Round-23 named finding: the read is bounded at
/// [`valenx_core::io_caps::MAX_VOL_FILE_BYTES`] (4 GiB) — sister
/// to the gmsh cap; Netgen `.vol` ASCII is denser per element than
/// gmsh but high-resolution adaptive meshes can still cross GiB.
pub fn parse_file(path: &Path, mesh_id: &str) -> Result<Mesh, VolError> {
    let text = valenx_core::io_caps::read_capped_to_string(
        path,
        valenx_core::io_caps::MAX_VOL_FILE_BYTES as usize,
    )?;
    parse_str(&text, mesh_id)
}

/// Parse an in-memory `.vol` buffer.
pub fn parse_str(text: &str, mesh_id: &str) -> Result<Mesh, VolError> {
    let mut mesh = Mesh::new(mesh_id);

    // We collect surface elements only when the file declares
    // `dimension 2` — for 3D meshes the surface block describes
    // boundaries, which we'll wire into Mesh::boundaries in a future
    // commit; for 2D meshes it IS the element set.
    let mut dimension: Option<u32> = None;

    let mut iter = text.lines().enumerate().peekable();

    while let Some((line_no_zero, raw_line)) = iter.next() {
        let line_no = line_no_zero + 1;
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        match line {
            "mesh3d" => {
                // Format header — must appear before any data.
            }
            "mesh2d" => {
                return Err(VolError::UnsupportedHeader {
                    header: line.to_string(),
                });
            }
            "dimension" => {
                let (l, val) = read_uint_line(&mut iter)?;
                if val != 2 && val != 3 {
                    return Err(VolError::Malformed {
                        line: l,
                        reason: format!("unsupported dimension {val}; expected 2 or 3"),
                    });
                }
                dimension = Some(val as u32);
            }
            "geomtype" => {
                // Geometry-source tag (0=CSG, 11=OCC, 12=Mesh, 13=STL).
                // We don't act on it; just consume the integer that
                // follows so it doesn't trip the section dispatcher.
                let _ = read_uint_line(&mut iter)?;
            }
            "points" => {
                parse_points(&mut iter, &mut mesh, line_no)?;
            }
            "volumeelements" => {
                parse_volume_elements(&mut iter, &mut mesh, line_no)?;
            }
            "surfaceelements" | "surfaceelementsgi" => {
                if dimension == Some(2) {
                    parse_surface_elements_as_2d(&mut iter, &mut mesh, line_no)?;
                } else {
                    skip_section(&mut iter, line_no)?;
                }
            }
            "edgesegments" | "edgesegmentsgi" | "edgesegmentsgi2" => {
                skip_section(&mut iter, line_no)?;
            }
            "identifications"
            | "identificationtypes"
            | "materials"
            | "bcnames"
            | "cd2names"
            | "singular_points"
            | "singular_edge_left"
            | "singular_edge_right"
            | "singular_face_inside"
            | "singular_face_outside" => {
                skip_section(&mut iter, line_no)?;
            }
            "endmesh" => {
                break;
            }
            _ => {
                // Unknown header — try to skip the count + listed
                // lines so subsequent recognised sections still parse.
                // If the next line isn't an integer, surface a clear
                // error rather than silently skipping arbitrary text.
                let next_is_count = iter
                    .peek()
                    .map(|(_, t)| {
                        let trimmed = t.trim();
                        !trimmed.is_empty()
                            && !trimmed.starts_with('#')
                            && trimmed.parse::<u64>().is_ok()
                    })
                    .unwrap_or(false);
                if next_is_count {
                    skip_section(&mut iter, line_no)?;
                }
                // Otherwise treat as stray text (some legacy .vol
                // files include comment-like lines without a leading
                // `#`); do nothing.
            }
        }
    }

    // R34 H1: post-parse connectivity bounds check. `.vol` section
    // order is free (a `volumeelements` block may precede `points`),
    // so a per-line check has no final node count to compare against.
    // Now that every section is consumed and `mesh.nodes` is final,
    // reject any element node index that points past the node array —
    // mirrors the OBJ reject-at-parse guard (`obj.rs`) and seals the
    // hazard at the source: the bad index would otherwise reach the
    // shared mesh consumers (`valenx_mesh::decimate` / `boolean`),
    // which raw-index `positions[tri[k]]` and panic.
    let node_count = mesh.nodes.len() as u32;
    for block in &mesh.element_blocks {
        for &idx in &block.connectivity {
            if idx >= node_count {
                return Err(VolError::Malformed {
                    line: 0,
                    reason: format!(
                        "element references node index {} (1-based {}) but the mesh has only {} points",
                        idx,
                        idx as u64 + 1,
                        node_count
                    ),
                });
            }
        }
    }

    let _ = mesh.recompute_quality_stats();
    Ok(mesh)
}

/// Read a line that must contain a single non-negative integer.
/// Skips blank / comment lines on the way.
fn read_uint_line(
    iter: &mut std::iter::Peekable<std::iter::Enumerate<std::str::Lines<'_>>>,
) -> Result<(usize, u64), VolError> {
    for (line_no_zero, raw) in iter.by_ref() {
        let line_no = line_no_zero + 1;
        let trimmed = raw.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let val = trimmed.parse::<u64>().map_err(|_| VolError::Malformed {
            line: line_no,
            reason: format!("expected integer count, got `{trimmed}`"),
        })?;
        return Ok((line_no, val));
    }
    Err(VolError::Malformed {
        line: 0,
        reason: "expected integer count, got EOF".into(),
    })
}

/// Parse a `points <N>` block: `N` lines of `x y z`.
fn parse_points(
    iter: &mut std::iter::Peekable<std::iter::Enumerate<std::str::Lines<'_>>>,
    mesh: &mut Mesh,
    section_line: usize,
) -> Result<(), VolError> {
    let (_, n) = read_uint_line(iter)?;
    let mut read = 0u64;
    while read < n {
        let Some((line_no_zero, raw)) = iter.next() else {
            return Err(VolError::Malformed {
                line: section_line,
                reason: format!(
                    "points block declared {n} entries but only {read} present before EOF"
                ),
            });
        };
        let line_no = line_no_zero + 1;
        let trimmed = raw.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut parts = trimmed.split_ascii_whitespace();
        let x = parse_f64(&mut parts, line_no, "x")?;
        let y = parse_f64(&mut parts, line_no, "y")?;
        let z = parse_f64(&mut parts, line_no, "z")?;
        mesh.nodes.push(Vector3::new(x, y, z));
        read += 1;
    }
    Ok(())
}

fn parse_f64(
    parts: &mut std::str::SplitAsciiWhitespace<'_>,
    line: usize,
    field: &str,
) -> Result<f64, VolError> {
    let Some(tok) = parts.next() else {
        return Err(VolError::Malformed {
            line,
            reason: format!("expected `{field}`, got nothing"),
        });
    };
    tok.parse::<f64>().map_err(|_| VolError::Malformed {
        line,
        reason: format!("expected `{field}` as f64, got `{tok}`"),
    })
}

fn parse_u32(
    parts: &mut std::str::SplitAsciiWhitespace<'_>,
    line: usize,
    field: &str,
) -> Result<u32, VolError> {
    let Some(tok) = parts.next() else {
        return Err(VolError::Malformed {
            line,
            reason: format!("expected `{field}`, got nothing"),
        });
    };
    tok.parse::<u32>().map_err(|_| VolError::Malformed {
        line,
        reason: format!("expected `{field}` as u32, got `{tok}`"),
    })
}

/// Parse a `volumeelements <N>` block.
///
/// Each line is: `matnr np p1 p2 ... p_np`. `p*` are 1-based node
/// indices in the file; we rewrite to 0-based on insert.
fn parse_volume_elements(
    iter: &mut std::iter::Peekable<std::iter::Enumerate<std::str::Lines<'_>>>,
    mesh: &mut Mesh,
    section_line: usize,
) -> Result<(), VolError> {
    let (_, n) = read_uint_line(iter)?;
    let mut by_type: std::collections::BTreeMap<u32, ElementBlock> =
        std::collections::BTreeMap::new();
    let mut read = 0u64;
    while read < n {
        let Some((line_no_zero, raw)) = iter.next() else {
            return Err(VolError::Malformed {
                line: section_line,
                reason: format!(
                    "volumeelements declared {n} entries but only {read} present before EOF"
                ),
            });
        };
        let line_no = line_no_zero + 1;
        let trimmed = raw.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut parts = trimmed.split_ascii_whitespace();
        let _matnr = parse_u32(&mut parts, line_no, "matnr")?;
        let np = parse_u32(&mut parts, line_no, "np")?;
        let element_type = match np {
            4 => ElementType::Tet4,
            5 => ElementType::Pyr5,
            6 => ElementType::Prism6,
            8 => ElementType::Hex8,
            other => {
                return Err(VolError::UnsupportedElement {
                    section: "volumeelements",
                    np: other,
                });
            }
        };
        let block = by_type
            .entry(np)
            .or_insert_with(|| ElementBlock::new(element_type));
        for k in 0..np {
            let idx = parse_u32(&mut parts, line_no, &format!("p{}", k + 1))?;
            if idx == 0 {
                return Err(VolError::Malformed {
                    line: line_no,
                    reason: "0-based node index in .vol (Netgen indices are 1-based)".into(),
                });
            }
            block.connectivity.push(idx - 1);
        }
        read += 1;
    }
    for (_, block) in by_type {
        mesh.element_blocks.push(block);
    }
    Ok(())
}

/// Parse a `surfaceelements <N>` block, treating each entry as a 2D
/// element (Tri3/Quad4). Used only when the file declares
/// `dimension 2`; for 3D files the surface block describes boundaries
/// and is skipped (we may wire it into [`Mesh::boundaries`] in a
/// follow-up).
///
/// Each line is: `surfnr bcnr domin domout np p1..p_np`. We discard
/// the metadata and only read the connectivity.
fn parse_surface_elements_as_2d(
    iter: &mut std::iter::Peekable<std::iter::Enumerate<std::str::Lines<'_>>>,
    mesh: &mut Mesh,
    section_line: usize,
) -> Result<(), VolError> {
    let (_, n) = read_uint_line(iter)?;
    let mut by_type: std::collections::BTreeMap<u32, ElementBlock> =
        std::collections::BTreeMap::new();
    let mut read = 0u64;
    while read < n {
        let Some((line_no_zero, raw)) = iter.next() else {
            return Err(VolError::Malformed {
                line: section_line,
                reason: format!(
                    "surfaceelements declared {n} entries but only {read} present before EOF"
                ),
            });
        };
        let line_no = line_no_zero + 1;
        let trimmed = raw.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut parts = trimmed.split_ascii_whitespace();
        let _surfnr = parse_u32(&mut parts, line_no, "surfnr")?;
        let _bcnr = parse_u32(&mut parts, line_no, "bcnr")?;
        let _domin = parse_u32(&mut parts, line_no, "domin")?;
        let _domout = parse_u32(&mut parts, line_no, "domout")?;
        let np = parse_u32(&mut parts, line_no, "np")?;
        let element_type = match np {
            3 => ElementType::Tri3,
            4 => ElementType::Quad4,
            other => {
                return Err(VolError::UnsupportedElement {
                    section: "surfaceelements",
                    np: other,
                });
            }
        };
        let block = by_type
            .entry(np)
            .or_insert_with(|| ElementBlock::new(element_type));
        for k in 0..np {
            let idx = parse_u32(&mut parts, line_no, &format!("p{}", k + 1))?;
            if idx == 0 {
                return Err(VolError::Malformed {
                    line: line_no,
                    reason: "0-based node index in .vol (Netgen indices are 1-based)".into(),
                });
            }
            block.connectivity.push(idx - 1);
        }
        read += 1;
    }
    for (_, block) in by_type {
        mesh.element_blocks.push(block);
    }
    Ok(())
}

/// Skip a count-prefixed section we don't care about. Reads the next
/// integer N, then drops the next N non-blank, non-comment lines.
fn skip_section(
    iter: &mut std::iter::Peekable<std::iter::Enumerate<std::str::Lines<'_>>>,
    section_line: usize,
) -> Result<(), VolError> {
    let (_, n) = read_uint_line(iter)?;
    let mut dropped = 0u64;
    while dropped < n {
        let Some((_, raw)) = iter.next() else {
            return Err(VolError::Malformed {
                line: section_line,
                reason: format!(
                    "section declared {n} entries but only {dropped} present before EOF"
                ),
            });
        };
        let trimmed = raw.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        dropped += 1;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One-tet mesh — the simplest Netgen output that exercises the
    /// `points` + `volumeelements` parse paths.
    #[test]
    fn parses_single_tet_mesh() {
        let text = "\
mesh3d

dimension
3

geomtype
0

surfaceelements
0

volumeelements
1
1 4 1 2 3 4

# X         Y         Z
points
4
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0

endmesh
";
        let mesh = parse_str(text, "tet1").expect("parse single tet");
        assert_eq!(mesh.nodes.len(), 4);
        assert_eq!(mesh.element_blocks.len(), 1);
        let block = &mesh.element_blocks[0];
        assert!(matches!(block.element_type, ElementType::Tet4));
        // 1-based -> 0-based on import.
        assert_eq!(block.connectivity, vec![0, 1, 2, 3]);
        assert_eq!(mesh.stats.node_count, 4);
        assert_eq!(mesh.stats.element_count, 1);
    }

    /// Mixed Tet4 + Hex8 in the same volumeelements section — tests
    /// that we open separate ElementBlocks per element type.
    #[test]
    fn separates_tet_and_hex_into_distinct_blocks() {
        // 12 nodes: 4 for the tet (corner), 8 for the hex.
        let text = "\
mesh3d
dimension
3

points
12
0 0 0
1 0 0
0 1 0
0 0 1
2 0 0
3 0 0
3 1 0
2 1 0
2 0 1
3 0 1
3 1 1
2 1 1

volumeelements
2
1 4 1 2 3 4
1 8 5 6 7 8 9 10 11 12

endmesh
";
        let mesh = parse_str(text, "mixed").expect("parse mixed");
        assert_eq!(mesh.nodes.len(), 12);
        assert_eq!(mesh.element_blocks.len(), 2);
        // BTreeMap order: Tet4 (np=4) before Hex8 (np=8).
        assert!(matches!(
            mesh.element_blocks[0].element_type,
            ElementType::Tet4
        ));
        assert!(matches!(
            mesh.element_blocks[1].element_type,
            ElementType::Hex8
        ));
        assert_eq!(mesh.element_blocks[0].connectivity, vec![0, 1, 2, 3]);
        assert_eq!(
            mesh.element_blocks[1].connectivity,
            vec![4, 5, 6, 7, 8, 9, 10, 11]
        );
    }

    /// 2D meshes use surfaceelements as their element set.
    #[test]
    fn parses_2d_mesh_via_surface_elements() {
        let text = "\
mesh3d
dimension
2

points
4
0 0 0
1 0 0
1 1 0
0 1 0

surfaceelements
2
1 0 1 0 3 1 2 3
1 0 1 0 3 1 3 4

endmesh
";
        let mesh = parse_str(text, "square").expect("parse 2d");
        assert_eq!(mesh.nodes.len(), 4);
        assert_eq!(mesh.element_blocks.len(), 1);
        assert!(matches!(
            mesh.element_blocks[0].element_type,
            ElementType::Tri3
        ));
        assert_eq!(mesh.element_blocks[0].connectivity, vec![0, 1, 2, 0, 2, 3]);
    }

    /// 3D meshes leave surfaceelements unparsed (it's a boundary
    /// description, not a volume element set).
    #[test]
    fn three_d_mesh_skips_surface_elements_block() {
        let text = "\
mesh3d
dimension
3

surfaceelements
2
1 0 1 0 3 1 2 3
1 0 1 0 3 1 3 4

points
4
0 0 0
1 0 0
0 1 0
0 0 1

volumeelements
1
1 4 1 2 3 4

endmesh
";
        let mesh = parse_str(text, "skip3d").expect("parse 3d");
        assert_eq!(mesh.element_blocks.len(), 1);
        assert!(matches!(
            mesh.element_blocks[0].element_type,
            ElementType::Tet4
        ));
    }

    /// Identifications / materials / bcnames blocks should be skipped
    /// without poisoning subsequent sections.
    #[test]
    fn tolerates_unknown_count_prefixed_sections() {
        let text = "\
mesh3d
dimension
3

materials
1
1 default

bcnames
2
1 wall
2 inlet

points
1
0 0 0

volumeelements
0

endmesh
";
        let mesh = parse_str(text, "tolerant").expect("parse tolerant");
        assert_eq!(mesh.nodes.len(), 1);
        assert_eq!(mesh.element_blocks.len(), 0);
    }

    /// Quadratic / unsupported element types surface a structured
    /// error so the caller can act on it (rather than silently
    /// producing a malformed mesh).
    #[test]
    fn rejects_quadratic_tet_with_structured_error() {
        let text = "\
mesh3d
dimension
3

volumeelements
1
1 10 1 2 3 4 5 6 7 8 9 10

points
10
0 0 0
1 0 0
0 1 0
0 0 1
0.5 0 0
0.5 0.5 0
0 0.5 0
0 0 0.5
0.5 0 0.5
0 0.5 0.5

endmesh
";
        let err = parse_str(text, "quad-tet").expect_err("must reject np=10");
        assert!(matches!(
            err,
            VolError::UnsupportedElement {
                section: "volumeelements",
                np: 10
            }
        ));
    }

    /// R34 H1 (RED→GREEN): a volume element whose node index points
    /// past the `points` count is a hostile / corrupt `.vol`. Section
    /// order in `.vol` is free, so this can only be caught with a
    /// post-parse pass once `mesh.nodes` is final. Pre-fix the parser
    /// stored `idx - 1` blindly; the bad index then reached the shared
    /// mesh consumers (`valenx_mesh::decimate` / `boolean`) which
    /// raw-index `positions[tri[k]]` → panic. We must reject at parse.
    #[test]
    fn rejects_volume_element_index_past_node_count() {
        // 4 points (valid indices 1..=4), but the tet cites node 99.
        let text = "\
mesh3d
dimension
3

points
4
0 0 0
1 0 0
0 1 0
0 0 1

volumeelements
1
1 4 1 2 3 99

endmesh
";
        let err = parse_str(text, "oob-vol").expect_err("index past node count must error");
        match err {
            VolError::Malformed { reason, .. } => {
                assert!(
                    reason.contains("node index") && reason.contains("99"),
                    "expected out-of-range node index in reason; got {reason}"
                );
            }
            other => panic!("expected Malformed; got {other:?}"),
        }
    }

    /// R34 H1 (RED→GREEN): the same hazard via a 2D surface element.
    /// Also guards the free-section-order case: the offending element
    /// appears *before* the `points` block, so a per-line check would
    /// have nothing to compare against — only the post-parse pass
    /// catches it.
    #[test]
    fn rejects_surface_element_index_past_node_count_sections_reordered() {
        // surfaceelements precede points; the second tri cites node 7
        // but only 4 points are declared further down.
        let text = "\
mesh3d
dimension
2

surfaceelements
2
1 0 1 0 3 1 2 3
1 0 1 0 3 1 3 7

points
4
0 0 0
1 0 0
1 1 0
0 1 0

endmesh
";
        let err =
            parse_str(text, "oob-surf").expect_err("surface index past node count must error");
        match err {
            VolError::Malformed { reason, .. } => {
                assert!(
                    reason.contains("node index") && reason.contains("7"),
                    "expected out-of-range node index in reason; got {reason}"
                );
            }
            other => panic!("expected Malformed; got {other:?}"),
        }
    }

    /// 0-based indices in a Netgen file are a fatal corruption — the
    /// spec mandates 1-based.
    #[test]
    fn flags_zero_index_as_malformed() {
        let text = "\
mesh3d
dimension
3

points
4
0 0 0
1 0 0
0 1 0
0 0 1

volumeelements
1
1 4 0 1 2 3

endmesh
";
        let err = parse_str(text, "bad").expect_err("0 index must error");
        match err {
            VolError::Malformed { reason, .. } => {
                assert!(
                    reason.contains("0-based"),
                    "expected `0-based` in reason; got {reason}"
                );
            }
            other => panic!("expected Malformed; got {other:?}"),
        }
    }
}
