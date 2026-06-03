//! Minimal ASCII `.frd` (CalculiX result format) parser.
//!
//! CalculiX's full `.frd` format is large — a sparse mix of fixed-
//! column ASCII headers and node/element/result blocks, optionally
//! switchable to a binary variant. This parser handles the slice
//! that's actually useful for Valenx's results-rendering pipeline:
//!
//! - Top-level node coordinate block (`2C` … `-1` rows … `-3`).
//! - Per-step nodal field result blocks (`1PSTEP` then a `100C`
//!   group with `-4 <FIELD>` and one `-1 <node> <value …>` per
//!   node, terminated by `-3`).
//!
//! Element connectivity blocks (`3C`) are recognised and skipped —
//! the canonical mesh comes from the upstream gmsh adapter or from
//! the user's own mesh import; reconstructing it from the .frd would
//! risk inconsistency. Element-level results (stress at integration
//! points) are also skipped — surface rendering of element results
//! needs the boundary-extraction step that hasn't landed yet.
//!
//! Binary `.frd` variants return [`ParseError::Unsupported`]. The
//! contract is "ASCII in, real Fields out, nothing else."

use std::collections::HashMap;

use thiserror::Error;

/// One nodal field block from an ASCII `.frd`.
#[derive(Clone, Debug)]
pub struct FrdField {
    /// Field name as written in the `-4` line (e.g. "DISP", "STRESS",
    /// "NDTEMP", "PE"). Trimmed of leading/trailing whitespace.
    pub name: String,
    /// CalculiX step number (1-indexed). Encoded into the canonical
    /// `TimeKey::Iteration` by callers when the field is converted.
    pub step: u32,
    /// Number of components per node — derived from the `-4` line's
    /// declared component count.
    pub components: usize,
    /// Flat per-node values: `[v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, …]`
    /// for a 3-component field, etc. Indexed by `node_id - 1`.
    pub data: Vec<f64>,
}

impl FrdField {
    /// Number of nodes this field covers (`data.len() / components`,
    /// or `0` for the degenerate `components == 0` case).
    pub fn samples(&self) -> usize {
        if self.components == 0 {
            0
        } else {
            self.data.len() / self.components
        }
    }
}

/// Result of parsing a `.frd` file.
#[derive(Clone, Debug, Default)]
pub struct FrdData {
    /// Per-node `[x, y, z]` coordinates. Indexed by `node_id - 1`
    /// (CalculiX uses 1-based node IDs; we 0-index after parse).
    pub points: Vec<[f64; 3]>,
    /// Nodal fields, in the order the parser encountered them.
    pub fields: Vec<FrdField>,
}

/// Maximum number of nodes a single FRD field block can declare,
/// derived from the maximum `-1` row node id. Round-5 hardening:
/// the unbounded `vec![0.0; n_nodes * field_components]` allocation
/// (and the unbounded `data.points.resize(...)` in the coordinate
/// section) lets a malformed FRD with `node_id = u32::MAX` allocate
/// gigabytes before the parse fails downstream. 1M nodes covers
/// every legitimate engineering mesh and is roughly 80 MB at 8
/// f64-per-row.
const MAX_FRD_LIST_LEN: usize = 1_000_000;

/// Errors raised by [`parse_ascii`].
#[derive(Debug, Error, PartialEq)]
pub enum ParseError {
    /// The input file looked binary (high NUL-byte density), and the
    /// parser only handles the ASCII variant.
    #[error(
        "binary .frd is not supported — re-run CalculiX with the default \
         ASCII output (omit `*FILE FORMAT, BINARY`) or convert with ccx2paraview"
    )]
    Unsupported,
    /// A numeric token (coordinate, field value, etc.) failed to parse.
    #[error("malformed numeric data in {context}: {token}")]
    BadNumeric {
        /// Section of the file the bad token appeared in (e.g. `"node"`,
        /// `"field"`).
        context: &'static str,
        /// The offending text.
        token: String,
    },
    /// A field referenced a node id outside the range covered by the
    /// parsed coordinate table.
    #[error("missing or out-of-range node id {id} (max parsed: {max})")]
    BadNodeId {
        /// The bad node id.
        id: u32,
        /// Highest node id that was parsed from the coordinate table.
        max: u32,
    },
    /// A list (node block, field rows) declared a length above the
    /// safe cap. Round-5 DoS hardening — the canonical hostile shape
    /// is `node_id = i32::MAX` which would force a multi-GB allocation
    /// before the rest of the parse fails.
    #[error(
        "FRD list length {len} exceeds the safe cap {cap} \
         (likely a malformed or hostile file)"
    )]
    ListTooLarge {
        /// Length the file declared.
        len: usize,
        /// Cap (== `MAX_FRD_LIST_LEN`).
        cap: usize,
    },
}

/// Parse the ASCII shape of a `.frd` file.
pub fn parse_ascii(text: &str) -> Result<FrdData, ParseError> {
    // Heuristic binary check: real ASCII .frd files are pure printable
    // ASCII + newlines. A high non-printable density suggests binary.
    let nul_count = text.bytes().filter(|&b| b == 0).count();
    if nul_count > 0 {
        return Err(ParseError::Unsupported);
    }

    let mut data = FrdData::default();
    let mut node_ids_seen: HashMap<u32, usize> = HashMap::new();
    let mut current_step: u32 = 0;

    let mut lines = text.lines().peekable();
    while let Some(line) = lines.next() {
        // Section markers — these are fixed-column lines starting with
        // a small integer in the first 5 chars. Whitespace-tolerant.
        let trimmed = line.trim_start();
        if trimmed.starts_with("2C") {
            // Node coordinate block. Read `-1` lines until `-3`.
            for inner in lines.by_ref() {
                let it = inner.trim_start();
                if it.starts_with("-3") {
                    break;
                }
                if let Some(rest) = it.strip_prefix("-1") {
                    let mut iter = rest.split_ascii_whitespace();
                    let id_token = match iter.next() {
                        Some(t) => t,
                        None => continue,
                    };
                    let id: u32 = id_token.parse().map_err(|_| ParseError::BadNumeric {
                        context: "Node id",
                        token: id_token.to_string(),
                    })?;
                    let xs: Vec<f64> = iter
                        .take(3)
                        .map(|t| {
                            t.parse::<f64>().map_err(|_| ParseError::BadNumeric {
                                context: "Node coordinate",
                                token: t.to_string(),
                            })
                        })
                        .collect::<Result<_, _>>()?;
                    if xs.len() != 3 {
                        return Err(ParseError::BadNumeric {
                            context: "Node coordinate triple",
                            token: rest.trim().to_string(),
                        });
                    }
                    // CalculiX node ids are 1-indexed and may not be
                    // contiguous; we track them in a HashMap and grow
                    // `points` to fit. Out-of-band ids are tolerated.
                    let zero_idx = (id as usize).saturating_sub(1);
                    // Round-5 DoS hardening: cap the implied resize
                    // BEFORE it happens. A malformed FRD with
                    // `node_id = u32::MAX` would otherwise allocate
                    // ~96 GB (= u32::MAX * sizeof([f64;3])) before
                    // the rest of the parse rejected the input.
                    if zero_idx >= MAX_FRD_LIST_LEN {
                        return Err(ParseError::ListTooLarge {
                            len: zero_idx + 1,
                            cap: MAX_FRD_LIST_LEN,
                        });
                    }
                    if zero_idx >= data.points.len() {
                        data.points.resize(zero_idx + 1, [0.0; 3]);
                    }
                    data.points[zero_idx] = [xs[0], xs[1], xs[2]];
                    node_ids_seen.insert(id, zero_idx);
                }
            }
            continue;
        }
        if trimmed.starts_with("3C") {
            // Element connectivity block — skip everything until `-3`.
            // We don't reconstruct the mesh from .frd; the canonical
            // mesh comes from the meshing adapter upstream.
            for inner in lines.by_ref() {
                if inner.trim_start().starts_with("-3") {
                    break;
                }
            }
            continue;
        }
        if trimmed.starts_with("1PSTEP") {
            // Step header — extract the step number.
            let mut iter = trimmed.split_ascii_whitespace();
            iter.next(); // "1PSTEP"
            if let Some(s) = iter.next() {
                if let Ok(n) = s.parse::<u32>() {
                    current_step = n;
                }
            }
            continue;
        }
        if trimmed.starts_with("100C") || trimmed.starts_with("100CL") {
            // Field group header: the next `-4` line names the field
            // and declares its component count.
            let mut field_name: Option<String> = None;
            let mut field_components: usize = 1;
            let mut field_data: Vec<(u32, Vec<f64>)> = Vec::new();
            for inner in lines.by_ref() {
                let it = inner.trim_start();
                if it.starts_with("-3") {
                    // End of this field's data.
                    break;
                }
                if let Some(rest) = it.strip_prefix("-4") {
                    let mut iter = rest.split_ascii_whitespace();
                    field_name = iter.next().map(|s| s.to_string());
                    if let Some(c) = iter.next() {
                        if let Ok(n) = c.parse::<usize>() {
                            // CCX `-4 NAME components type` — components
                            // for DISP is 4 (3 components + magnitude
                            // marker), but we only emit the first 3
                            // when we see `-5 D1/D2/D3`. For simplicity
                            // we treat `components - 1` as the real
                            // count when components > 1, except for
                            // single-scalar blocks (NDTEMP etc.).
                            field_components = match n {
                                4 => 3,
                                7 => 6,
                                other => other.max(1),
                            };
                        }
                    }
                    continue;
                }
                if it.starts_with("-5") {
                    // Per-component header — we don't need the names
                    // for anything yet; the data rows are positional.
                    continue;
                }
                if let Some(rest) = it.strip_prefix("-1") {
                    let mut iter = rest.split_ascii_whitespace();
                    let id_token = match iter.next() {
                        Some(t) => t,
                        None => continue,
                    };
                    let id: u32 = id_token.parse().map_err(|_| ParseError::BadNumeric {
                        context: "Field row node id",
                        token: id_token.to_string(),
                    })?;
                    let values: Vec<f64> = iter
                        .take(field_components)
                        .map(|t| {
                            t.parse::<f64>().map_err(|_| ParseError::BadNumeric {
                                context: "Field value",
                                token: t.to_string(),
                            })
                        })
                        .collect::<Result<_, _>>()?;
                    field_data.push((id, values));
                }
            }
            // Materialise the field into the FrdData if we saw both
            // a name and data.
            if let Some(name) = field_name {
                if !field_data.is_empty() {
                    let max_id = field_data.iter().map(|(id, _)| *id).max().unwrap_or(0);
                    let n_nodes = data.points.len().max(max_id as usize);
                    // Round-5 DoS hardening: cap the field allocation
                    // size BEFORE the vec! macro reserves it. A
                    // hostile `-1 4294967295 ...` row would otherwise
                    // force a multi-GB allocation. Cap matches the
                    // same MAX_FRD_LIST_LEN used in the coordinate
                    // block above.
                    if n_nodes > MAX_FRD_LIST_LEN {
                        return Err(ParseError::ListTooLarge {
                            len: n_nodes,
                            cap: MAX_FRD_LIST_LEN,
                        });
                    }
                    let total = n_nodes
                        .checked_mul(field_components)
                        .ok_or(ParseError::ListTooLarge {
                            len: n_nodes,
                            cap: MAX_FRD_LIST_LEN,
                        })?;
                    let mut flat: Vec<f64> = vec![0.0; total];
                    for (id, values) in &field_data {
                        let zero_idx = (*id as usize).saturating_sub(1);
                        if zero_idx >= n_nodes {
                            return Err(ParseError::BadNodeId {
                                id: *id,
                                max: n_nodes as u32,
                            });
                        }
                        let start = zero_idx * field_components;
                        for (i, v) in values.iter().take(field_components).enumerate() {
                            flat[start + i] = *v;
                        }
                    }
                    data.fields.push(FrdField {
                        name: name.trim().to_string(),
                        step: current_step.max(1),
                        components: field_components,
                        data: flat,
                    });
                }
            }
            continue;
        }
        // Anything else (header lines `1C`, etc.) — skip.
    }
    let _ = node_ids_seen;
    Ok(data)
}

/// Convert this [`FrdData`] into canonical [`crate::Field`] entries.
/// Each `FrdField` becomes one `Field` with `Location::OnNode`,
/// dimensionless units (CalculiX's .frd doesn't carry units), and a
/// `TimeKey::Iteration(step)` keyed off the field's step number.
pub fn to_canonical_fields(data: &FrdData) -> Vec<crate::Field> {
    data.fields
        .iter()
        .map(|f| {
            let kind = match f.components {
                1 => crate::FieldKind::Scalar,
                3 => crate::FieldKind::Vector { dim: 3 },
                6 => crate::FieldKind::Tensor { rows: 3, cols: 3 },
                n => crate::FieldKind::Vector { dim: n as u8 },
            };
            let range = field_range(&f.data);
            crate::Field {
                name: f.name.clone(),
                kind,
                location: crate::Location::OnNode,
                region: crate::RegionRef("default".to_string()),
                units: crate::units::DIMENSIONLESS,
                time: crate::TimeKey::Iteration(f.step as u64),
                data: f.data.clone(),
                range,
            }
        })
        .collect()
}

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

    /// Hand-crafted minimal ASCII .frd: 4 nodes + one DISP block at
    /// step 1. Trimmed of the long fixed-column padding CCX writes
    /// since our parser is whitespace-tolerant.
    const ONE_DISP_FRD: &str = r#"
    1C  UTILITY_PROGRAMS
    2C
 -1   1  0.00000E+00 0.00000E+00 0.00000E+00
 -1   2  1.00000E+00 0.00000E+00 0.00000E+00
 -1   3  0.00000E+00 1.00000E+00 0.00000E+00
 -1   4  0.00000E+00 0.00000E+00 1.00000E+00
 -3
    3C
 -1   1   1  0  STEEL
 -2   1   2   3   4
 -3
    1PSTEP    1    1    1
  100CL  102  0.00000E+00       4    1    1
 -4  DISP    4    1
 -5  D1    1    2    1    0
 -5  D2    1    2    2    0
 -5  D3    1    2    3    0
 -1   1  0.00000E+00 0.00000E+00 0.00000E+00
 -1   2  1.23000E-04 0.00000E+00 0.00000E+00
 -1   3  0.00000E+00 4.56000E-04 0.00000E+00
 -1   4  0.00000E+00 0.00000E+00 7.89000E-04
 -3
"#;

    #[test]
    fn parses_minimal_disp_frd() {
        let data = parse_ascii(ONE_DISP_FRD).expect("parse");
        assert_eq!(data.points.len(), 4);
        // Last node = (0, 0, 1).
        assert_eq!(data.points[3], [0.0, 0.0, 1.0]);
        assert_eq!(data.fields.len(), 1);
        let disp = &data.fields[0];
        assert_eq!(disp.name, "DISP");
        assert_eq!(disp.step, 1);
        // CCX `-4 DISP 4` declares 4 components but only 3 are real
        // (D1/D2/D3); the parser strips the trailing magnitude
        // marker and reports 3.
        assert_eq!(disp.components, 3);
        assert_eq!(disp.samples(), 4);
        // Node 4 displacement = (0, 0, 7.89e-4).
        assert!((disp.data[3 * 3] - 0.0).abs() < 1e-12);
        assert!((disp.data[3 * 3 + 1] - 0.0).abs() < 1e-12);
        assert!((disp.data[3 * 3 + 2] - 7.89e-4).abs() < 1e-12);
    }

    #[test]
    fn rejects_binary_frd_with_clear_error() {
        // Binary .frd contains NUL bytes interspersed with text.
        let bytes = b"    1C  UTILITY_PROGRAMS\n\0\0\0BINARYDATA\0\0\n".to_vec();
        let text = String::from_utf8_lossy(&bytes).to_string();
        let err = parse_ascii(&text).unwrap_err();
        assert_eq!(err, ParseError::Unsupported);
    }

    #[test]
    fn to_canonical_disp_becomes_vector_field() {
        let data = parse_ascii(ONE_DISP_FRD).expect("parse");
        let fields = to_canonical_fields(&data);
        assert_eq!(fields.len(), 1);
        let disp = &fields[0];
        assert_eq!(disp.name, "DISP");
        assert_eq!(disp.kind, crate::FieldKind::Vector { dim: 3 });
        assert_eq!(disp.location, crate::Location::OnNode);
        assert_eq!(disp.time, crate::TimeKey::Iteration(1));
        // Cached range: includes the 7.89e-4 max and the 0.0 min.
        assert!(disp.range.is_some());
        let (min, max) = disp.range.unwrap();
        assert!((min - 0.0).abs() < 1e-12);
        assert!((max - 7.89e-4).abs() < 1e-12);
    }

    #[test]
    fn handles_temperature_scalar_field() {
        // Single-component nodal field (NDTEMP).
        let text = r#"
    1C
    2C
 -1   1  0 0 0
 -1   2  1 0 0
 -3
    1PSTEP    1    1    1
  100CL  102  0.00000E+00       1    1    1
 -4  NDTEMP    1    1
 -5  T    1    1    0    0
 -1   1  300.0
 -1   2  350.0
 -3
"#;
        let data = parse_ascii(text).expect("parse");
        assert_eq!(data.fields.len(), 1);
        let t = &data.fields[0];
        assert_eq!(t.name, "NDTEMP");
        assert_eq!(t.components, 1);
        assert_eq!(t.data, vec![300.0, 350.0]);

        let canonical = to_canonical_fields(&data);
        assert_eq!(canonical[0].kind, crate::FieldKind::Scalar);
    }

    /// Round-5 RED→GREEN: a malformed `.frd` with an oversized node
    /// id (here i32::MAX, the canonical hostile shape) used to force
    /// the parser to allocate ~96 GB (= max_id * sizeof([f64;3])) in
    /// the coordinate block before any later check rejected the
    /// input. The fix is a `MAX_FRD_LIST_LEN` cap on the implied
    /// resize, surfacing as `ParseError::ListTooLarge`.
    #[test]
    fn rejects_oversized_node_count() {
        // Hand-craft an FRD with one node at id 2147483647 (i32::MAX,
        // which is what CalculiX-style external code can legitimately
        // synthesise from a corrupted MPI rank index). Without the cap,
        // the parser would attempt `data.points.resize(2147483647 + 1, ...)`
        // before failing.
        let text = "    1C\n    2C\n -1   2147483647  0 0 0\n -3\n";
        let err = parse_ascii(text).expect_err("must reject oversized node id");
        match err {
            ParseError::ListTooLarge { len, cap } => {
                assert!(
                    len >= 2_147_483_647,
                    "len {len} should be >= i32::MAX, got {len}"
                );
                assert_eq!(cap, MAX_FRD_LIST_LEN);
            }
            other => panic!("expected ListTooLarge, got {other:?}"),
        }
    }

    /// Round-5 sister test: the same cap fires when the FIELD block
    /// has an oversized node id (the second allocation site in the
    /// parser).
    #[test]
    fn rejects_oversized_field_row_node_id() {
        // Build an FRD where the coordinate block has a tiny single
        // node but the field block references an oversized id. This
        // exercises the second allocation site (the `vec![0.0; n_nodes
        // * field_components]` in the field-block path).
        let text = "    1C\n    2C\n -1   1  0 0 0\n -3\n    1PSTEP    1    1    1\n  \
                    100CL  102  0  1  1  1\n -4  NDTEMP    1    1\n -5  T    1    1    0    0\n \
                    -1   2147483647  300.0\n -3\n";
        let err = parse_ascii(text).expect_err("must reject oversized field node id");
        match err {
            ParseError::ListTooLarge { len, cap } => {
                assert!(
                    len >= 2_147_483_647,
                    "len {len} should be >= i32::MAX, got {len}"
                );
                assert_eq!(cap, MAX_FRD_LIST_LEN);
            }
            other => panic!("expected ListTooLarge, got {other:?}"),
        }
    }

    #[test]
    fn multi_step_frd_keeps_separate_fields_per_step() {
        // Same DISP field at two steps — should produce two FrdFields.
        let text = r#"
    1C
    2C
 -1   1  0 0 0
 -1   2  1 0 0
 -3
    1PSTEP    1    1    1
  100CL  102  0  4  1  1
 -4  DISP    4    1
 -5  D1    1    2    1    0
 -5  D2    1    2    2    0
 -5  D3    1    2    3    0
 -1   1  0 0 0
 -1   2  1e-4 0 0
 -3
    1PSTEP    2    1    1
  100CL  102  0  4  1  1
 -4  DISP    4    1
 -5  D1    1    2    1    0
 -5  D2    1    2    2    0
 -5  D3    1    2    3    0
 -1   1  0 0 0
 -1   2  2e-4 0 0
 -3
"#;
        let data = parse_ascii(text).expect("parse");
        assert_eq!(data.fields.len(), 2);
        assert_eq!(data.fields[0].step, 1);
        assert_eq!(data.fields[1].step, 2);
        // Step 2's node 2 displacement x = 2e-4.
        assert!((data.fields[1].data[3] - 2e-4).abs() < 1e-12);

        let canonical = to_canonical_fields(&data);
        assert_eq!(canonical[0].time, crate::TimeKey::Iteration(1));
        assert_eq!(canonical[1].time, crate::TimeKey::Iteration(2));
    }
}
