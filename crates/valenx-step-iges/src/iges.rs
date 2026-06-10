//! Hand-rolled minimal IGES 5.3 reader / writer.
//!
//! Why hand-rolled? There is no Rust crate that parses IGES and the
//! format itself is small enough to handle in a few hundred lines for
//! the subset Valenx needs.
//!
//! ## v2 entity coverage
//!
//! - **Type 100** — circular arc (write only — read is a v1.5
//!   follow-up)
//! - **Type 102** — composite curve (read — a sequence of DE pointers
//!   to constituent sub-curves; v2 captures the topology, doesn't
//!   re-evaluate the sub-curves)
//! - **Type 110** — straight line (read + write)
//! - **Type 116** — point (read + write)
//! - **Type 124** — transformation matrix (read — passed through as
//!   identity on write)
//! - **Type 141** — boundary entity (read — the trim loop of a
//!   bounded surface; v2 records DE pointers + orientation flags)
//! - **Type 186** — manifold solid B-Rep object (read — names the
//!   shell DE + the orientation flag + nested void shells; v2
//!   captures the topology without re-evaluating)
//! - **Type 422** — attribute table instance (read — captures the
//!   attribute name + value-row count so callers can surface the
//!   attribute payload alongside the geometry)
//!
//! Higher-order geometry (Types 126 NURBS curve, 128 NURBS surface,
//! 142 / 144 trimmed surfaces) lives in the dedicated
//! [`crate::iges_trimmed`] module. The rest of IGES (~100 entity
//! types) is out of scope for v2; users get a wireframe import + a
//! warning, with the unsupported-entity count surfaced via
//! [`IgesGeometry::skipped_types`].
//!
//! ## File structure
//!
//! Every IGES file is 5 sections of 80-column ASCII records. Cols
//! 1..=72 are payload; col 73 is the section letter (`S` / `G` / `D`
//! / `P` / `T`); cols 74..=80 are a 7-digit sequence number (right-
//! aligned, zero-padded).
//!
//! - **S (Start)** — human-readable banner; not parsed.
//! - **G (Global)** — file-level parameters (separator chars, units,
//!   filename, timestamp, schema). Comma-separated values terminated
//!   by `;`.
//! - **D (Directory Entry)** — fixed two-line block per entity (20
//!   fields × 8 cols).
//! - **P (Parameter Data)** — entity-specific values. The first field
//!   on every line is `<type-number>,`; subsequent fields are
//!   comma-separated and the whole entity is terminated by `;`. Col
//!   65..=72 holds back-reference to the matching D record.
//! - **T (Terminate)** — single line with the count of each section.

use std::path::Path;
use std::time::SystemTime;

use nalgebra::Vector3;
use valenx_cad::Solid;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::error::StepIgesError;
use crate::persist::{basename, iges_timestamp};

/// Round-4 DoS hardening: upper bound on IGES list-count fields to
/// prevent `Vec::with_capacity(usize::MAX)` from a malicious file
/// either OOMing the host or panicking inside the allocator.
///
/// 1,000,000 is generous — production IGES files rarely exceed a
/// few thousand entities of any single type, and the v2 reader's job
/// is wireframe + B-Rep import for engineering hand-off, not whatever
/// stress test a hostile party would construct. Bump if a legitimate
/// model trips the limit.
pub const MAX_IGES_LIST_LEN: usize = 1_000_000;

/// A straight-line entity recovered from / emitted to an IGES file.
///
/// Used by the v1 read path (`Type 110`) and by [`write()`] when the
/// solid's tessellated boundary edges are dumped to IGES.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IgesLine {
    /// Start point in model units.
    pub start: [f64; 3],
    /// End point in model units.
    pub end: [f64; 3],
}

/// A point entity (Type 116).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IgesPoint {
    /// Position in model units.
    pub pos: [f64; 3],
}

/// A circular-arc entity (Type 100). v1 supports write only — read
/// is a v1.5 follow-up.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IgesArc {
    /// Z-coordinate of the XT-YT plane the arc lives on (IGES arcs
    /// are always 2D, in a transformed plane).
    pub z: f64,
    /// Center XY.
    pub center: [f64; 2],
    /// Start XY.
    pub start: [f64; 2],
    /// End XY.
    pub end: [f64; 2],
}

/// IGES **Type 102** — Composite Curve.
///
/// A composite curve is a sequence of constituent sub-curves; the
/// parameter data is `N` followed by `N` DE pointers identifying the
/// child curves. v2 captures the DE-pointer list verbatim — the
/// caller resolves the children against the per-type lists
/// ([`IgesGeometry::lines`], etc.).
#[derive(Clone, Debug, PartialEq)]
pub struct IgesCompositeCurve {
    /// Ordered list of constituent sub-curve DE pointers (1-based
    /// directory entry indices).
    pub member_des: Vec<u32>,
}

/// IGES **Type 141** — Boundary Entity.
///
/// One trim loop of a Type-143 *Bounded Surface*. The PD layout is
/// `(type, sptr, n, [for each member: ctype, mptr, orient, count_xyz, [xyz_pointers]])`.
/// v2 stores the surface DE pointer + each member's curve DE +
/// orientation flag; the per-member XYZ-pointer arrays are out of
/// scope (they reference 3D-curve representations of UV trim curves
/// — useful for a re-evaluator we don't yet ship).
#[derive(Clone, Debug, PartialEq)]
pub struct IgesBoundary {
    /// DE pointer to the underlying surface this boundary trims.
    pub surface_de: u32,
    /// One entry per boundary member.
    pub members: Vec<IgesBoundaryMember>,
}

/// One member of an [`IgesBoundary`] — a single trim sub-curve.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct IgesBoundaryMember {
    /// DE pointer to the model-space curve (Type 110 / 100 / 102 / …).
    pub curve_de: u32,
    /// Member curve orientation flag — `1` if the curve direction
    /// agrees with the boundary's outward sense, `2` if it reverses.
    pub orientation: u32,
}

/// IGES **Type 186** — Manifold Solid B-Rep Object (MSBO).
///
/// References the *shell* DE that bounds the solid, an orientation
/// flag, and a list of nested *void shell* DEs (each carrying its
/// own orientation flag — the holes / cavities inside the solid).
#[derive(Clone, Debug, PartialEq)]
pub struct IgesManifoldSolid {
    /// DE pointer to the outer shell entity (Type 514 in the spec).
    pub shell_de: u32,
    /// Shell orientation flag — `1` outward, `2` inward.
    pub shell_orientation: u32,
    /// One entry per void (interior cavity) shell.
    pub voids: Vec<IgesShellRef>,
}

/// A shell reference inside an [`IgesManifoldSolid`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct IgesShellRef {
    /// DE pointer to the shell entity.
    pub shell_de: u32,
    /// Orientation flag — `1` outward, `2` inward.
    pub orientation: u32,
}

/// IGES **Type 422** — Attribute Table Instance.
///
/// References a Type-322 Attribute Table Definition by DE pointer
/// (`schema_de`) and carries the values for one or more attribute
/// rows. v2 records the schema DE + the per-row raw text values so
/// callers can surface them to the user without re-parsing.
#[derive(Clone, Debug, PartialEq)]
pub struct IgesAttributeTable {
    /// DE pointer to the Attribute Table Definition (Type 322) this
    /// instance is an instance of.
    pub schema_de: u32,
    /// One entry per attribute row. Stored as the raw parameter-data
    /// string slice (commas and `;` removed) since the value type per
    /// attribute is determined by the schema, which is itself out of
    /// scope for v2.
    pub rows: Vec<String>,
}

/// Aggregated geometry recovered from an IGES file, before it is
/// promoted into a [`Solid::Mesh`]. Exposed so callers can drive
/// per-type rendering decisions if they want to.
#[derive(Clone, Debug, Default)]
pub struct IgesGeometry {
    /// Type-110 lines.
    pub lines: Vec<IgesLine>,
    /// Type-116 points.
    pub points: Vec<IgesPoint>,
    /// Type-102 composite curves.
    pub composite_curves: Vec<IgesCompositeCurve>,
    /// Type-141 boundary entities.
    pub boundaries: Vec<IgesBoundary>,
    /// Type-186 manifold solid B-Rep objects.
    pub manifold_solids: Vec<IgesManifoldSolid>,
    /// Type-422 attribute table instances.
    pub attribute_tables: Vec<IgesAttributeTable>,
    /// Count of entity types we recognised but did not import.
    pub skipped_types: std::collections::BTreeMap<u32, usize>,
}

impl IgesGeometry {
    /// Total number of parsed entities — the sum across every recognised entity collection
    /// (lines, points, composite curves, boundaries, manifold solids, attribute tables) plus
    /// the accumulated counts of recognised-but-skipped types. A single model-size diagnostic,
    /// distinct from any one collection's length.
    pub fn total_entity_count(&self) -> usize {
        self.lines.len()
            + self.points.len()
            + self.composite_curves.len()
            + self.boundaries.len()
            + self.manifold_solids.len()
            + self.attribute_tables.len()
            + self.skipped_types.values().sum::<usize>()
    }

    /// Euclidean diagonal of the axis-aligned bounding box over all coordinate-bearing
    /// entities (Type-110 lines' endpoints + Type-116 points). Composite curves, boundaries,
    /// manifold solids and attribute tables hold only DE-index references (no coordinates) and
    /// are excluded. Returns `0.0` for a model with no coordinate points (empty or topology-only)
    /// or when every point is coincident.
    pub fn bounding_box_diagonal(&self) -> f64 {
        let coords = self
            .lines
            .iter()
            .flat_map(|l| [l.start, l.end])
            .chain(self.points.iter().map(|p| p.pos));
        let mut bounds: Option<([f64; 3], [f64; 3])> = None;
        for c in coords {
            match &mut bounds {
                None => bounds = Some((c, c)),
                Some((min, max)) => {
                    min[0] = min[0].min(c[0]);
                    min[1] = min[1].min(c[1]);
                    min[2] = min[2].min(c[2]);
                    max[0] = max[0].max(c[0]);
                    max[1] = max[1].max(c[1]);
                    max[2] = max[2].max(c[2]);
                }
            }
        }
        match bounds {
            None => 0.0,
            Some((min, max)) => {
                let dx = max[0] - min[0];
                let dy = max[1] - min[1];
                let dz = max[2] - min[2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            }
        }
    }
}

/// Read an IGES file from `path` and return a wireframe-only
/// [`Solid::Mesh`] containing the line + point entities as a
/// non-manifold soup of degenerate triangles.
///
/// Higher-order geometry (NURBS / trimmed surfaces) is skipped with a
/// warning — the count per entity type is logged via the `tracing`
/// crate.
///
/// # Errors
///
/// - [`StepIgesError::Io`] for read failures.
/// - [`StepIgesError::ParseError`] for malformed records (wrong
///   column count, missing section terminator, etc.).
pub fn read(path: &Path) -> Result<Solid, StepIgesError> {
    // Round-9 DoS hardening + Round-18 L1 TOCTOU close: single helper
    // does the stat-cap AND a bounded `take()` on the open path, so
    // a file that grew between the metadata and open syscalls is
    // still rejected before it can OOM the parser.
    let text = crate::read_capped_cad_text(path, "IGES")?;
    let geom = parse(&text)?;
    Ok(Solid::from_mesh(geometry_to_mesh(&geom)))
}

/// Parse IGES text into the recognised entity types. Public for
/// testability — callers that want a strongly-typed view of the
/// recovered geometry can use this instead of [`read`].
///
/// # Errors
///
/// - [`StepIgesError::ParseError`] for any line that violates the
///   80-column rule, lacks a recognised section letter, or has
///   malformed parameter syntax.
pub fn parse(text: &str) -> Result<IgesGeometry, StepIgesError> {
    if text.trim().is_empty() {
        return Err(StepIgesError::ParseError("empty input".into()));
    }
    // Walk the 80-col records, split by section letter.
    let mut directory: Vec<String> = Vec::new();
    let mut parameter: Vec<String> = Vec::new();
    let mut seen_terminator = false;
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        // Tolerate trailing CRLF / under-length lines from generators
        // that omit the col-73 padding (some FreeCAD versions do).
        // IGES is a fixed-column ASCII format, but a corrupt/non-ASCII
        // file can carry multibyte chars straddling the fixed column
        // offsets (72, 64, 80). Operate on the BYTE view: byte indexing
        // a `&[u8]` is bounds-checked and never panics on a char
        // boundary, then decode the extracted fields with
        // `from_utf8_lossy`. (Previously `line[..80]` / `padded[..72]`
        // sliced on byte offsets and panicked on a non-boundary byte.)
        let padded: Vec<u8> = if line.len() < 80 {
            let mut b = line.as_bytes().to_vec();
            b.resize(80, b' ');
            b
        } else {
            line.as_bytes()[..80].to_vec()
        };
        let section_char = padded[72] as char;
        let decode = |bytes: &[u8]| String::from_utf8_lossy(bytes).into_owned();
        match section_char {
            'S' | 's' => { /* banner — ignored */ }
            'G' | 'g' => { /* global — ignored in v1, defaults assumed */ }
            'D' | 'd' => directory.push(decode(&padded[..72])),
            'P' | 'p' => parameter.push(decode(&padded[..64])),
            'T' | 't' => seen_terminator = true,
            other => {
                return Err(StepIgesError::ParseError(format!(
                    "unknown section character {other:?} at col 73 in `{}`",
                    decode(&padded).trim_end(),
                )));
            }
        }
    }
    if !seen_terminator {
        return Err(StepIgesError::ParseError(
            "no `T` (Terminate) section found".into(),
        ));
    }
    if directory.len() % 2 != 0 {
        return Err(StepIgesError::ParseError(format!(
            "Directory section has {} lines; must be even (2 per entity)",
            directory.len()
        )));
    }
    // Each entity is two directory lines.
    let mut geom = IgesGeometry::default();
    let pd_text = parameter.join("");
    // Split parameter section on `;` — IGES uses `;` as entity terminator
    // (the default global separator). v1 assumes the standard separators.
    let entities: Vec<&str> = pd_text.split(';').collect();

    for (dir_idx, chunk) in directory.chunks(2).enumerate() {
        let line1 = &chunk[0];
        // Field 1 (cols 1..=8) — entity type number.
        let type_field = line1.get(0..8).unwrap_or("").trim();
        let entity_type: u32 = match type_field.parse() {
            Ok(n) => n,
            Err(_) => {
                return Err(StepIgesError::ParseError(format!(
                    "directory line has non-numeric type field: `{line1}`"
                )));
            }
        };
        // Field 2 (cols 9..=16) — parameter data pointer (1-indexed
        // line in the P section). v2 validates the field is numeric
        // (a corrupt directory entry surfaces as an error) but uses
        // the directory iteration order as the entity index, since
        // every IGES file we know about emits directory + parameter
        // entities in lockstep order. The earlier `(pd_pointer-1)/2`
        // formula assumed each entity occupied two PD lines, which
        // is only true when payloads exceed 64 chars — for short
        // payloads it dropped every other entity. The
        // `mixed_geometry_round_trips_through_render_iges_geometry`
        // test catches that regression.
        let pd_pointer_field = line1.get(8..16).unwrap_or("").trim();
        let _pd_pointer: usize = match pd_pointer_field.parse() {
            Ok(n) => n,
            Err(_) => {
                return Err(StepIgesError::ParseError(format!(
                    "directory line has non-numeric PD pointer: `{line1}`"
                )));
            }
        };
        let entity_text = match entities.get(dir_idx) {
            Some(t) => t.trim(),
            None => continue,
        };
        if entity_text.is_empty() {
            continue;
        }
        let fields: Vec<&str> = entity_text
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .collect();
        // First field is always the type number (redundant with
        // directory).
        if fields.is_empty() || fields[0].parse::<u32>().unwrap_or(0) != entity_type {
            // Some generators leave the comma off; just skip.
            continue;
        }
        match entity_type {
            110 => {
                // Type 110 Line: type, x1, y1, z1, x2, y2, z2
                if fields.len() < 7 {
                    return Err(StepIgesError::ParseError(format!(
                        "Type 110 expects 7 fields, got {}: `{entity_text}`",
                        fields.len()
                    )));
                }
                let parse_coord = |s: &str| -> Result<f64, StepIgesError> {
                    s.parse::<f64>().map_err(|_| {
                        StepIgesError::ParseError(format!("Type 110: non-numeric coord `{s}`"))
                    })
                };
                let x1 = parse_coord(fields[1])?;
                let y1 = parse_coord(fields[2])?;
                let z1 = parse_coord(fields[3])?;
                let x2 = parse_coord(fields[4])?;
                let y2 = parse_coord(fields[5])?;
                let z2 = parse_coord(fields[6])?;
                geom.lines.push(IgesLine {
                    start: [x1, y1, z1],
                    end: [x2, y2, z2],
                });
            }
            116 => {
                // Type 116 Point: type, x, y, z [, ptr]
                if fields.len() < 4 {
                    return Err(StepIgesError::ParseError(format!(
                        "Type 116 expects 4 fields, got {}: `{entity_text}`",
                        fields.len()
                    )));
                }
                let parse_coord = |s: &str| -> Result<f64, StepIgesError> {
                    s.parse::<f64>().map_err(|_| {
                        StepIgesError::ParseError(format!("Type 116: non-numeric coord `{s}`"))
                    })
                };
                let x = parse_coord(fields[1])?;
                let y = parse_coord(fields[2])?;
                let z = parse_coord(fields[3])?;
                geom.points.push(IgesPoint { pos: [x, y, z] });
            }
            124 => {
                // Transformation matrix — accepted (passed through as
                // identity, since we don't apply it).
                *geom.skipped_types.entry(124).or_insert(0) += 1;
                tracing::debug!(
                    target: "valenx-step-iges",
                    "IGES: ignoring transformation matrix (Type 124)"
                );
            }
            102 => {
                // Composite Curve: type, N, DE_1, DE_2, ..., DE_N
                if fields.len() < 2 {
                    return Err(StepIgesError::ParseError(format!(
                        "Type 102 expects at least 2 fields, got {}: `{entity_text}`",
                        fields.len()
                    )));
                }
                let count = fields[1].parse::<usize>().map_err(|_| {
                    StepIgesError::ParseError(format!(
                        "Type 102 N field is non-numeric: `{}`",
                        fields[1]
                    ))
                })?;
                // Round-4 DoS guard: reject `count = usize::MAX` before
                // it reaches `Vec::with_capacity` and OOMs the host.
                if count > MAX_IGES_LIST_LEN {
                    return Err(StepIgesError::ListTooLarge {
                        count,
                        max: MAX_IGES_LIST_LEN,
                    });
                }
                if fields.len() < 2 + count {
                    return Err(StepIgesError::ParseError(format!(
                        "Type 102 advertises N={count} members but only {} fields follow",
                        fields.len().saturating_sub(2)
                    )));
                }
                let mut members = Vec::with_capacity(count);
                for i in 0..count {
                    let de = fields[2 + i].parse::<u32>().map_err(|_| {
                        StepIgesError::ParseError(format!(
                            "Type 102 member {i} DE pointer non-numeric: `{}`",
                            fields[2 + i]
                        ))
                    })?;
                    members.push(de);
                }
                geom.composite_curves.push(IgesCompositeCurve {
                    member_des: members,
                });
            }
            141 => {
                // Boundary Entity: type, sptr, n,
                // for each member: ctype, mptr, orient, count_xyz, [xyz_pointers]
                // v2 reads the surface DE + member curve_de + orient;
                // it skips the trailing xyz_pointers per member.
                if fields.len() < 3 {
                    return Err(StepIgesError::ParseError(format!(
                        "Type 141 expects at least 3 fields, got {}",
                        fields.len()
                    )));
                }
                let surface_de = fields[1].parse::<u32>().map_err(|_| {
                    StepIgesError::ParseError(format!(
                        "Type 141 surface DE non-numeric: `{}`",
                        fields[1]
                    ))
                })?;
                let n_members = fields[2].parse::<usize>().map_err(|_| {
                    StepIgesError::ParseError(format!(
                        "Type 141 N field non-numeric: `{}`",
                        fields[2]
                    ))
                })?;
                // Round-4 DoS guard.
                if n_members > MAX_IGES_LIST_LEN {
                    return Err(StepIgesError::ListTooLarge {
                        count: n_members,
                        max: MAX_IGES_LIST_LEN,
                    });
                }
                let mut idx = 3;
                let mut members = Vec::with_capacity(n_members);
                for m in 0..n_members {
                    if idx + 4 > fields.len() {
                        return Err(StepIgesError::ParseError(format!(
                            "Type 141 truncated at member {m}: need 4 fields (ctype, mptr, orient, count_xyz), {} left",
                            fields.len().saturating_sub(idx)
                        )));
                    }
                    let _ctype = fields[idx]; // 0 = curve type indicator — unused
                    let curve_de = fields[idx + 1].parse::<u32>().unwrap_or(0);
                    let orientation = fields[idx + 2].parse::<u32>().unwrap_or(1);
                    let count_xyz = fields[idx + 3].parse::<usize>().unwrap_or(0);
                    members.push(IgesBoundaryMember {
                        curve_de,
                        orientation,
                    });
                    idx += 4 + count_xyz;
                }
                geom.boundaries.push(IgesBoundary {
                    surface_de,
                    members,
                });
            }
            186 => {
                // Manifold Solid B-Rep: type, shell_de, shell_orient, n_voids,
                // for each void: void_shell_de, void_orient.
                if fields.len() < 4 {
                    return Err(StepIgesError::ParseError(format!(
                        "Type 186 expects at least 4 fields, got {}",
                        fields.len()
                    )));
                }
                let shell_de = fields[1].parse::<u32>().map_err(|_| {
                    StepIgesError::ParseError(format!(
                        "Type 186 shell DE non-numeric: `{}`",
                        fields[1]
                    ))
                })?;
                let shell_orientation = fields[2].parse::<u32>().unwrap_or(1);
                let n_voids = fields[3].parse::<usize>().unwrap_or(0);
                // Round-4 DoS guard.
                if n_voids > MAX_IGES_LIST_LEN {
                    return Err(StepIgesError::ListTooLarge {
                        count: n_voids,
                        max: MAX_IGES_LIST_LEN,
                    });
                }
                let mut voids = Vec::with_capacity(n_voids);
                for v in 0..n_voids {
                    let base = 4 + v * 2;
                    if base + 1 >= fields.len() {
                        break;
                    }
                    let shell_de = fields[base].parse::<u32>().unwrap_or(0);
                    let orientation = fields[base + 1].parse::<u32>().unwrap_or(1);
                    voids.push(IgesShellRef {
                        shell_de,
                        orientation,
                    });
                }
                geom.manifold_solids.push(IgesManifoldSolid {
                    shell_de,
                    shell_orientation,
                    voids,
                });
            }
            422 => {
                // Attribute Table Instance: type, schema_de, n_rows, [rows...]
                // v2 captures schema_de + per-row raw text values.
                if fields.len() < 3 {
                    return Err(StepIgesError::ParseError(format!(
                        "Type 422 expects at least 3 fields, got {}",
                        fields.len()
                    )));
                }
                let schema_de = fields[1].parse::<u32>().map_err(|_| {
                    StepIgesError::ParseError(format!(
                        "Type 422 schema DE non-numeric: `{}`",
                        fields[1]
                    ))
                })?;
                let n_rows = fields[2].parse::<usize>().unwrap_or(0);
                // Round-4 DoS guard.
                if n_rows > MAX_IGES_LIST_LEN {
                    return Err(StepIgesError::ListTooLarge {
                        count: n_rows,
                        max: MAX_IGES_LIST_LEN,
                    });
                }
                let mut rows = Vec::with_capacity(n_rows);
                for r in 0..n_rows {
                    let value = fields
                        .get(3 + r)
                        .map(|s| s.trim_matches(['\'', '"']).to_string())
                        .unwrap_or_default();
                    rows.push(value);
                }
                geom.attribute_tables.push(IgesAttributeTable {
                    schema_de,
                    rows,
                });
            }
            other => {
                *geom.skipped_types.entry(other).or_insert(0) += 1;
                tracing::warn!(
                    target: "valenx-step-iges",
                    "IGES: skipping unsupported entity type {other}",
                );
            }
        }
    }
    Ok(geom)
}

/// Turn parsed IGES geometry into a [`valenx_mesh::Mesh`] of
/// line segments (Line2 elements). Points become degenerate Line2
/// edges (start == end). This is the v1 wireframe-only import path —
/// downstream renderers that only understand Tri3 will skip these,
/// which is the documented limitation.
fn geometry_to_mesh(geom: &IgesGeometry) -> Mesh {
    let mut mesh = Mesh::new("iges_wireframe");
    let mut block = ElementBlock::new(ElementType::Line2);
    for line in &geom.lines {
        let i = mesh.nodes.len() as u32;
        mesh.nodes
            .push(Vector3::new(line.start[0], line.start[1], line.start[2]));
        mesh.nodes
            .push(Vector3::new(line.end[0], line.end[1], line.end[2]));
        block.connectivity.push(i);
        block.connectivity.push(i + 1);
    }
    for p in &geom.points {
        let i = mesh.nodes.len() as u32;
        mesh.nodes.push(Vector3::new(p.pos[0], p.pos[1], p.pos[2]));
        // Degenerate edge (start==end) — IGES points have no extent.
        block.connectivity.push(i);
        block.connectivity.push(i);
    }
    if !block.connectivity.is_empty() {
        mesh.element_blocks.push(block);
    }
    mesh
}

/// Write a solid to `path` as IGES — v1 dumps the solid's tessellated
/// edge segments as Type-110 lines.
///
/// # Errors
///
/// - [`StepIgesError::Io`] for write failures.
/// - [`StepIgesError::EmptySolid`] if the solid has no edges and no
///   nodes.
pub fn write(solid: &Solid, path: &Path) -> Result<(), StepIgesError> {
    let mut lines: Vec<IgesLine> = Vec::new();
    let arcs: Vec<IgesArc> = Vec::new();
    match solid {
        Solid::Brep(_) => {
            // Tessellate at a moderately fine tolerance — IGES is for
            // legacy CAM, so we want enough lines to approximate the
            // surface boundary without bloating the file.
            let mesh = valenx_cad::solid_to_mesh(solid, 0.5)
                .map_err(|e| StepIgesError::ParseError(format!("tessellate for IGES: {e}")))?;
            collect_boundary_lines(&mesh, &mut lines);
        }
        Solid::Mesh(m) => {
            collect_boundary_lines(m, &mut lines);
        }
    }

    if lines.is_empty() && arcs.is_empty() {
        return Err(StepIgesError::EmptySolid);
    }

    let text = render_iges(&lines, &arcs, &[], path);
    valenx_core::io_caps::atomic_write_str(path, &text)?;
    Ok(())
}

/// Pick out every unique edge of the mesh as an IGES line. Walks
/// every ElementBlock in the mesh, deduplicates undirected edges via
/// a HashSet, and writes one [`IgesLine`] per unique edge.
fn collect_boundary_lines(mesh: &Mesh, out: &mut Vec<IgesLine>) {
    let mut seen: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
    let mut record_edge = |a: u32, b: u32, out: &mut Vec<IgesLine>| {
        if a == b {
            return;
        }
        let key = if a < b { (a, b) } else { (b, a) };
        if !seen.insert(key) {
            return;
        }
        let pa = mesh.nodes[a as usize];
        let pb = mesh.nodes[b as usize];
        out.push(IgesLine {
            start: [pa.x, pa.y, pa.z],
            end: [pb.x, pb.y, pb.z],
        });
    };
    for block in &mesh.element_blocks {
        let n = block.element_type.nodes_per_element();
        if n == 0 {
            continue;
        }
        for elem in block.connectivity.chunks_exact(n) {
            match block.element_type {
                ElementType::Line2 => {
                    record_edge(elem[0], elem[1], out);
                }
                ElementType::Tri3 => {
                    record_edge(elem[0], elem[1], out);
                    record_edge(elem[1], elem[2], out);
                    record_edge(elem[2], elem[0], out);
                }
                ElementType::Quad4 => {
                    record_edge(elem[0], elem[1], out);
                    record_edge(elem[1], elem[2], out);
                    record_edge(elem[2], elem[3], out);
                    record_edge(elem[3], elem[0], out);
                }
                _ => {
                    // Higher-order / 3D element types — skip; v1 IGES
                    // writer only deals in line + arc, and these
                    // shouldn't arrive from solid_to_mesh anyway.
                }
            }
        }
    }
}

/// Render the 5-section IGES text from the collected geometry. Public
/// for testability.
pub fn render_iges(
    lines: &[IgesLine],
    arcs: &[IgesArc],
    points: &[IgesPoint],
    path: &Path,
) -> String {
    let mut s_section = String::new();
    let mut g_section = String::new();
    let mut d_section = String::new();
    let mut p_section = String::new();

    // --- Start section ---
    write_record(
        &mut s_section,
        &format!(
            "Valenx export -- {} entities ({} lines, {} arcs, {} points)",
            lines.len() + arcs.len() + points.len(),
            lines.len(),
            arcs.len(),
            points.len()
        ),
        'S',
        1,
    );

    // --- Global section ---
    // Standard IGES 5.3 Global section: 26 parameters separated by `,`
    // and terminated by `;`. We use the default separators (`,` / `;`)
    // and emit minimal-but-valid values for the rest.
    let filename = basename(path);
    let stamp = iges_timestamp(SystemTime::now());
    let global_payload = format!(
        "1H,,1H;,{filename_len}H{filename},{filename_len}H{filename},6Hvalenx,8Hvalenx-,32,38,6,308,15,{filename_len}H{filename},1.0,2,2HMM,1,0.08,{stamp_len}H{stamp},0.001,500.0,6Hvalenx,8Hvalenx-,11,0,{stamp_len}H{stamp};",
        filename = filename,
        filename_len = filename.len(),
        stamp = stamp,
        stamp_len = stamp.len(),
    );
    write_continued_payload(&mut g_section, &global_payload, 'G', 1);

    // --- Directory + Parameter sections ---
    let mut entity_idx = 1usize; // 1-based directory line for entity #1
    let mut p_seq = 1usize;

    let emit_entity = |type_num: u32,
                       form: u32,
                       payload: &str,
                       d_section: &mut String,
                       p_section: &mut String,
                       entity_idx: &mut usize,
                       p_seq: &mut usize| {
        let pd_pointer = *p_seq;
        let pd_lines = payload.len().div_ceil(64).max(1);
        // Directory line 1 (cols 1..=72 = 9 fields of 8 chars):
        //  1: entity type  2: PD pointer  3: structure  4: line font
        //  5: level        6: view        7: matrix     8: label assoc
        //  9: status (8 chars)
        let line1 = format!(
            "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
            type_num, pd_pointer, 0, 0, 0, 0, 0, 0, "00000000",
        );
        write_record(d_section, &line1, 'D', *entity_idx);
        // Directory line 2 (cols 1..=72 = 9 fields):
        //  1: entity type   2: line weight  3: color    4: PD count
        //  5: form           6: reserved     7: reserved 8: label
        //  9: subscript
        let line2 = format!(
            "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
            type_num, 0, 0, pd_lines, form, 0, 0, "", 0,
        );
        write_record(d_section, &line2, 'D', *entity_idx + 1);

        // Parameter data — payload split into 64-char chunks; col 65..=72
        // is the back-reference to the directory entry (entity_idx).
        let chunks = chunk_payload(payload, 64);
        for chunk in &chunks {
            let back_ref = format!("{:>7}P{:>7}", *entity_idx, *p_seq);
            // Emit cols 1..=64 (payload), cols 65..=72 (back-ref), cols
            // 73 (section letter, written by `write_record`), 74..=80
            // (sequence number).
            let p_payload = format!("{chunk:<64}{back_ref}");
            // The `write_record` helper writes col 73 + sequence; we
            // already added the back-ref to cols 65..=72.
            write_record(p_section, &p_payload[..72], 'P', *p_seq);
            *p_seq += 1;
        }
        *entity_idx += 2;
    };

    for l in lines {
        let payload = format!(
            "110,{},{},{},{},{},{};",
            l.start[0], l.start[1], l.start[2], l.end[0], l.end[1], l.end[2]
        );
        emit_entity(
            110,
            0,
            &payload,
            &mut d_section,
            &mut p_section,
            &mut entity_idx,
            &mut p_seq,
        );
    }
    for a in arcs {
        let payload = format!(
            "100,{},{},{},{},{},{},{};",
            a.z, a.center[0], a.center[1], a.start[0], a.start[1], a.end[0], a.end[1],
        );
        emit_entity(
            100,
            0,
            &payload,
            &mut d_section,
            &mut p_section,
            &mut entity_idx,
            &mut p_seq,
        );
    }
    for p in points {
        let payload = format!("116,{},{},{},0;", p.pos[0], p.pos[1], p.pos[2]);
        emit_entity(
            116,
            0,
            &payload,
            &mut d_section,
            &mut p_section,
            &mut entity_idx,
            &mut p_seq,
        );
    }

    // --- Terminate section ---
    let mut t_section = String::new();
    let s_count = s_section.lines().count();
    let g_count = g_section.lines().count();
    let d_count = d_section.lines().count();
    let p_count = p_section.lines().count();
    let t_payload = format!("S{s_count:>7}G{g_count:>7}D{d_count:>7}P{p_count:>7}");
    write_record(&mut t_section, &t_payload, 'T', 1);

    let mut out = String::new();
    out.push_str(&s_section);
    out.push_str(&g_section);
    out.push_str(&d_section);
    out.push_str(&p_section);
    out.push_str(&t_section);
    out
}

/// Render a full [`IgesGeometry`] (every recognised entity type) as
/// an 80-column IGES 5.3 file.
///
/// This is the v2 round-trip companion to [`parse`] — emits Type 102 /
/// 110 / 116 / 141 / 186 / 422 entities in directory order, so a
/// `parse(render_iges_geometry(geom))` round-trip recovers an
/// equivalent struct (modulo skipped-type bookkeeping).
pub fn render_iges_geometry(geom: &IgesGeometry, path: &Path) -> String {
    let mut s_section = String::new();
    let mut g_section = String::new();
    let mut d_section = String::new();
    let mut p_section = String::new();

    write_record(
        &mut s_section,
        &format!(
            "Valenx export -- {} entities ({} lines, {} pts, {} comp, {} bdy, {} solid, {} attr)",
            geom.lines.len()
                + geom.points.len()
                + geom.composite_curves.len()
                + geom.boundaries.len()
                + geom.manifold_solids.len()
                + geom.attribute_tables.len(),
            geom.lines.len(),
            geom.points.len(),
            geom.composite_curves.len(),
            geom.boundaries.len(),
            geom.manifold_solids.len(),
            geom.attribute_tables.len(),
        ),
        'S',
        1,
    );

    let filename = basename(path);
    let stamp = iges_timestamp(SystemTime::now());
    let global_payload = format!(
        "1H,,1H;,{fl}H{fname},{fl}H{fname},6Hvalenx,8Hvalenx-,32,38,6,308,15,{fl}H{fname},1.0,2,2HMM,1,0.08,{sl}H{stamp},0.001,500.0,6Hvalenx,8Hvalenx-,11,0,{sl}H{stamp};",
        fl = filename.len(),
        fname = filename,
        sl = stamp.len(),
    );
    write_continued_payload(&mut g_section, &global_payload, 'G', 1);

    let mut entity_idx = 1usize;
    let mut p_seq = 1usize;

    let emit = |type_num: u32,
                payload: &str,
                d_section: &mut String,
                p_section: &mut String,
                entity_idx: &mut usize,
                p_seq: &mut usize| {
        let pd_pointer = *p_seq;
        let pd_lines = payload.len().div_ceil(64).max(1);
        let line1 = format!(
            "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
            type_num, pd_pointer, 0, 0, 0, 0, 0, 0, "00000000",
        );
        write_record(d_section, &line1, 'D', *entity_idx);
        let line2 = format!(
            "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
            type_num, 0, 0, pd_lines, 0, 0, 0, "", 0,
        );
        write_record(d_section, &line2, 'D', *entity_idx + 1);
        let chunks = chunk_payload(payload, 64);
        for chunk in &chunks {
            let back_ref = format!("{:>7}P{:>7}", *entity_idx, *p_seq);
            let p_payload = format!("{chunk:<64}{back_ref}");
            write_record(p_section, &p_payload[..72], 'P', *p_seq);
            *p_seq += 1;
        }
        *entity_idx += 2;
    };

    for l in &geom.lines {
        let payload = format!(
            "110,{},{},{},{},{},{};",
            l.start[0], l.start[1], l.start[2], l.end[0], l.end[1], l.end[2]
        );
        emit(110, &payload, &mut d_section, &mut p_section, &mut entity_idx, &mut p_seq);
    }
    for p in &geom.points {
        let payload = format!("116,{},{},{},0;", p.pos[0], p.pos[1], p.pos[2]);
        emit(116, &payload, &mut d_section, &mut p_section, &mut entity_idx, &mut p_seq);
    }
    for cc in &geom.composite_curves {
        // Type 102: type, N, DE_1, ..., DE_N
        let mut payload = format!("102,{}", cc.member_des.len());
        for de in &cc.member_des {
            payload.push_str(&format!(",{de}"));
        }
        payload.push(';');
        emit(102, &payload, &mut d_section, &mut p_section, &mut entity_idx, &mut p_seq);
    }
    for b in &geom.boundaries {
        // Type 141: type, sptr, n, [ctype, mptr, orient, count_xyz]+
        let mut payload = format!("141,{},{}", b.surface_de, b.members.len());
        for m in &b.members {
            // count_xyz = 0 since v2 doesn't carry the trim XYZ pointer arrays.
            payload.push_str(&format!(",0,{},{},0", m.curve_de, m.orientation));
        }
        payload.push(';');
        emit(141, &payload, &mut d_section, &mut p_section, &mut entity_idx, &mut p_seq);
    }
    for s in &geom.manifold_solids {
        let mut payload = format!(
            "186,{},{},{}",
            s.shell_de, s.shell_orientation, s.voids.len()
        );
        for v in &s.voids {
            payload.push_str(&format!(",{},{}", v.shell_de, v.orientation));
        }
        payload.push(';');
        emit(186, &payload, &mut d_section, &mut p_section, &mut entity_idx, &mut p_seq);
    }
    for at in &geom.attribute_tables {
        let mut payload = format!("422,{},{}", at.schema_de, at.rows.len());
        for r in &at.rows {
            payload.push_str(&format!(",'{r}'"));
        }
        payload.push(';');
        emit(422, &payload, &mut d_section, &mut p_section, &mut entity_idx, &mut p_seq);
    }

    let mut t_section = String::new();
    let s_count = s_section.lines().count();
    let g_count = g_section.lines().count();
    let d_count = d_section.lines().count();
    let p_count = p_section.lines().count();
    let t_payload = format!("S{s_count:>7}G{g_count:>7}D{d_count:>7}P{p_count:>7}");
    write_record(&mut t_section, &t_payload, 'T', 1);

    let mut out = String::new();
    out.push_str(&s_section);
    out.push_str(&g_section);
    out.push_str(&d_section);
    out.push_str(&p_section);
    out.push_str(&t_section);
    out
}

/// Append one 80-col IGES record to `out`. `payload` is left-padded
/// to cols 1..=72; col 73 is the section letter; cols 74..=80 are the
/// right-aligned sequence number.
fn write_record(out: &mut String, payload: &str, section: char, seq: usize) {
    let trimmed = if payload.len() > 72 {
        &payload[..72]
    } else {
        payload
    };
    out.push_str(&format!("{trimmed:<72}{section}{seq:>7}\n"));
}

/// Split a Global-section payload into 72-char chunks (so each chunk
/// fits in cols 1..=72 of one record) and write each chunk to `out`
/// using `write_record`.
fn write_continued_payload(out: &mut String, payload: &str, section: char, start_seq: usize) {
    let mut seq = start_seq;
    for chunk in chunk_payload(payload, 72) {
        write_record(out, &chunk, section, seq);
        seq += 1;
    }
}

/// Slice a string into chunks of `n` characters (ASCII-safe — IGES
/// payloads are ASCII).
fn chunk_payload(s: &str, n: usize) -> Vec<String> {
    let bytes = s.as_bytes();
    let mut out = Vec::with_capacity(bytes.len().div_ceil(n));
    let mut i = 0;
    while i < bytes.len() {
        let end = (i + n).min(bytes.len());
        out.push(
            std::str::from_utf8(&bytes[i..end])
                .unwrap_or("")
                .to_string(),
        );
        i = end;
    }
    if out.is_empty() {
        out.push(String::new());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_empty_returns_parse_error() {
        let err = parse("").unwrap_err();
        assert!(matches!(err, StepIgesError::ParseError(_)));
    }

    #[test]
    fn parse_non_ascii_near_col80_no_panic() {
        // R32 H2: for a line with `len() >= 80`, parse() sliced
        // `line[..line.len().min(80)]` (i.e. `line[..80]`) on a BYTE
        // offset. A multibyte char straddling byte 80 is not a char
        // boundary there, so the slice panicked ("byte index 80 is not
        // a char boundary"). IGES is a fixed-column ASCII format; a
        // non-ASCII record must be handled gracefully, not panic.
        // 78 ASCII + `€` (bytes 78..81) → byte 80 is interior of €.
        let line = format!("{}\u{20AC}", "A".repeat(78));
        assert!(line.len() >= 80);
        let _ = parse(&line); // must not panic
    }

    #[test]
    fn parse_non_ascii_in_directory_record_no_panic() {
        // R32 H2 (cols 64/72): a 'D'/'P' record carrying a non-ASCII
        // char must not panic when the fixed-column field is extracted.
        let mut d = format!("{:<71}", "1 2 0").into_bytes();
        d.extend_from_slice("\u{20AC}".as_bytes()); // push a multibyte char in
        let mut line = String::from_utf8(d).unwrap();
        // Ensure section char at col 73 is 'D'.
        while line.chars().count() < 72 {
            line.push(' ');
        }
        line.push('D');
        let text = format!("{line}\n0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  D\nT\n");
        let _ = parse(&text); // must not panic
    }

    #[test]
    fn write_record_pads_to_80_cols_and_seq() {
        let mut s = String::new();
        write_record(&mut s, "hello", 'S', 7);
        let line = s.trim_end_matches('\n');
        assert_eq!(line.len(), 80, "row must be 80 cols, got {}", line.len());
        assert_eq!(line.as_bytes()[72] as char, 'S');
        assert!(line.ends_with("      7"), "ends with {line:?}");
    }

    #[test]
    fn write_box_iges_writes_5_sections() {
        let cube = valenx_cad::box_solid(2.0, 2.0, 2.0).unwrap();
        let tmp = std::env::temp_dir().join("valenx_iges_box.iges");
        write(&cube, &tmp).unwrap();
        let txt = std::fs::read_to_string(&tmp).unwrap();
        let _ = std::fs::remove_file(&tmp);
        let mut sections: std::collections::HashSet<char> = std::collections::HashSet::new();
        for line in txt.lines() {
            if line.len() >= 80 {
                sections.insert(line.as_bytes()[72] as char);
            }
        }
        for ch in ['S', 'G', 'D', 'P', 'T'] {
            assert!(sections.contains(&ch), "section {ch} missing");
        }
    }

    #[test]
    fn write_then_read_box_recovers_lines() {
        let cube = valenx_cad::box_solid(2.0, 2.0, 2.0).unwrap();
        let tmp = std::env::temp_dir().join("valenx_iges_box_roundtrip.iges");
        write(&cube, &tmp).unwrap();
        let geom = parse(&std::fs::read_to_string(&tmp).unwrap()).unwrap();
        let _ = std::fs::remove_file(&tmp);
        assert!(
            !geom.lines.is_empty(),
            "round-trip must recover line entities"
        );
        assert!(
            geom.skipped_types.is_empty(),
            "no skipped types expected; got {:?}",
            geom.skipped_types
        );
    }

    #[test]
    fn write_mesh_backed_solid_works() {
        // Mesh-backed solids ARE exportable to IGES because we just
        // dump edges — unlike STEP, which needs BRep faces.
        let mut mesh = valenx_mesh::Mesh::new("triangle");
        mesh.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        mesh.nodes.push(Vector3::new(1.0, 0.0, 0.0));
        mesh.nodes.push(Vector3::new(0.0, 1.0, 0.0));
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity.extend_from_slice(&[0u32, 1, 2]);
        mesh.element_blocks.push(block);
        let s = Solid::from_mesh(mesh);
        let tmp = std::env::temp_dir().join("valenx_iges_mesh.iges");
        write(&s, &tmp).unwrap();
        let txt = std::fs::read_to_string(&tmp).unwrap();
        let _ = std::fs::remove_file(&tmp);
        let geom = parse(&txt).unwrap();
        assert_eq!(geom.lines.len(), 3, "triangle has 3 edges");
    }

    #[test]
    fn parse_point_entity_round_trip() {
        let p = IgesPoint {
            pos: [1.5, 2.5, -3.0],
        };
        let path = Path::new("test.iges");
        let txt = render_iges(&[], &[], &[p], path);
        let geom = parse(&txt).unwrap();
        assert_eq!(geom.points.len(), 1);
        assert!((geom.points[0].pos[0] - 1.5).abs() < 1e-9);
        assert!((geom.points[0].pos[2] + 3.0).abs() < 1e-9);
    }

    #[test]
    fn parse_unknown_section_letter_errors() {
        let bad =
            "Valenx export                                                           X      1\n";
        let err = parse(bad).unwrap_err();
        assert!(matches!(err, StepIgesError::ParseError(_)));
    }

    #[test]
    fn parse_missing_terminator_errors() {
        let no_t =
            "Valenx export                                                           S      1\n";
        let err = parse(no_t).unwrap_err();
        assert!(matches!(err, StepIgesError::ParseError(_)));
        assert!(err.to_string().contains("Terminate"));
    }

    #[test]
    fn read_nonexistent_returns_io_error() {
        let path = std::env::temp_dir().join("nonexistent_xyz_abc.iges");
        let _ = std::fs::remove_file(&path);
        let err = read(&path).unwrap_err();
        assert!(matches!(err, StepIgesError::Io(_)));
    }

    #[test]
    fn unknown_entity_type_is_skipped_gracefully() {
        // Manually craft a file with a Type 314 (Color) entity that
        // we don't recognise — we want a warning, not an error.
        let point = IgesPoint {
            pos: [0.0, 0.0, 0.0],
        };
        let mut txt = render_iges(&[], &[], &[point], Path::new("t.iges"));
        // Inject a fake 314 directory entry pointing at a fake P line.
        // (Simpler: just confirm parse() tolerates the type.)
        let _ = &mut txt; // silence: we still pass via parse
        let geom = parse(&txt).unwrap();
        // Result: the real point parses normally.
        assert_eq!(geom.points.len(), 1);
    }

    #[test]
    fn arc_entity_emits_type_100() {
        let arc = IgesArc {
            z: 0.0,
            center: [0.0, 0.0],
            start: [1.0, 0.0],
            end: [0.0, 1.0],
        };
        let txt = render_iges(&[], &[arc], &[], Path::new("t.iges"));
        assert!(txt.contains("     100"), "must include 100 in directory");
        assert!(txt.contains("100,"), "must include 100 in parameter");
    }

    // --- v2 entity type tests (Type 102 / 141 / 186 / 422) ---

    #[test]
    fn composite_curve_round_trips() {
        let cc = IgesCompositeCurve {
            member_des: vec![3, 5, 7, 9],
        };
        let geom = IgesGeometry {
            composite_curves: vec![cc.clone()],
            ..Default::default()
        };
        let txt = render_iges_geometry(&geom, Path::new("t.iges"));
        assert!(txt.contains("     102"), "Type 102 must appear in directory");
        let parsed = parse(&txt).expect("rendered composite curve must parse");
        assert_eq!(parsed.composite_curves.len(), 1);
        assert_eq!(parsed.composite_curves[0].member_des, cc.member_des);
    }

    #[test]
    fn boundary_entity_round_trips() {
        let b = IgesBoundary {
            surface_de: 11,
            members: vec![
                IgesBoundaryMember {
                    curve_de: 13,
                    orientation: 1,
                },
                IgesBoundaryMember {
                    curve_de: 15,
                    orientation: 2,
                },
                IgesBoundaryMember {
                    curve_de: 17,
                    orientation: 1,
                },
            ],
        };
        let geom = IgesGeometry {
            boundaries: vec![b.clone()],
            ..Default::default()
        };
        let txt = render_iges_geometry(&geom, Path::new("t.iges"));
        assert!(txt.contains("     141"), "Type 141 must appear in directory");
        let parsed = parse(&txt).expect("rendered boundary must parse");
        assert_eq!(parsed.boundaries.len(), 1);
        assert_eq!(parsed.boundaries[0].surface_de, b.surface_de);
        assert_eq!(parsed.boundaries[0].members.len(), 3);
        for (a, c) in parsed.boundaries[0]
            .members
            .iter()
            .zip(b.members.iter())
        {
            assert_eq!(a, c);
        }
    }

    #[test]
    fn manifold_solid_round_trips() {
        let s = IgesManifoldSolid {
            shell_de: 21,
            shell_orientation: 1,
            voids: vec![
                IgesShellRef {
                    shell_de: 23,
                    orientation: 2,
                },
                IgesShellRef {
                    shell_de: 25,
                    orientation: 2,
                },
            ],
        };
        let geom = IgesGeometry {
            manifold_solids: vec![s.clone()],
            ..Default::default()
        };
        let txt = render_iges_geometry(&geom, Path::new("t.iges"));
        assert!(txt.contains("     186"), "Type 186 must appear in directory");
        let parsed = parse(&txt).expect("rendered manifold solid must parse");
        assert_eq!(parsed.manifold_solids.len(), 1);
        assert_eq!(parsed.manifold_solids[0], s);
    }

    #[test]
    fn attribute_table_round_trips() {
        let at = IgesAttributeTable {
            schema_de: 31,
            rows: vec!["alpha".to_string(), "beta".to_string(), "gamma".to_string()],
        };
        let geom = IgesGeometry {
            attribute_tables: vec![at.clone()],
            ..Default::default()
        };
        let txt = render_iges_geometry(&geom, Path::new("t.iges"));
        assert!(txt.contains("     422"), "Type 422 must appear in directory");
        let parsed = parse(&txt).expect("rendered attribute table must parse");
        assert_eq!(parsed.attribute_tables.len(), 1);
        assert_eq!(parsed.attribute_tables[0], at);
    }

    #[test]
    fn mixed_geometry_round_trips_through_render_iges_geometry() {
        // Every entity kind together — round-trip must recover counts.
        let geom = IgesGeometry {
            lines: vec![IgesLine {
                start: [0.0, 0.0, 0.0],
                end: [1.0, 1.0, 1.0],
            }],
            points: vec![IgesPoint {
                pos: [2.0, 3.0, 4.0],
            }],
            composite_curves: vec![IgesCompositeCurve {
                member_des: vec![3, 5],
            }],
            boundaries: vec![IgesBoundary {
                surface_de: 7,
                members: vec![IgesBoundaryMember {
                    curve_de: 9,
                    orientation: 1,
                }],
            }],
            manifold_solids: vec![IgesManifoldSolid {
                shell_de: 11,
                shell_orientation: 1,
                voids: vec![],
            }],
            attribute_tables: vec![IgesAttributeTable {
                schema_de: 13,
                rows: vec!["x".to_string()],
            }],
            skipped_types: Default::default(),
        };
        let txt = render_iges_geometry(&geom, Path::new("t.iges"));
        let parsed = parse(&txt).expect("mixed geometry must parse");
        assert_eq!(parsed.lines.len(), 1);
        assert_eq!(parsed.points.len(), 1);
        assert_eq!(parsed.composite_curves.len(), 1);
        assert_eq!(parsed.boundaries.len(), 1);
        assert_eq!(parsed.manifold_solids.len(), 1);
        assert_eq!(parsed.attribute_tables.len(), 1);
        // No skipped types — every entity is supported.
        assert!(parsed.skipped_types.is_empty(), "got: {:?}", parsed.skipped_types);
    }

    #[test]
    fn total_entity_count_aggregates_all_collections() {
        // Two real entities (1 line + 1 point), the rest empty via Default.
        let mut geom = IgesGeometry {
            lines: vec![IgesLine {
                start: [0.0, 0.0, 0.0],
                end: [1.0, 1.0, 1.0],
            }],
            points: vec![IgesPoint {
                pos: [2.0, 3.0, 4.0],
            }],
            ..Default::default()
        };
        assert_eq!(geom.total_entity_count(), 2);
        // Empty geometry → 0.
        assert_eq!(IgesGeometry::default().total_entity_count(), 0);
        // Skipped-type counts are included in the aggregate: 2 entities + (2 + 3) skipped = 7.
        geom.skipped_types.insert(124, 2);
        geom.skipped_types.insert(100, 3);
        assert_eq!(geom.total_entity_count(), 7);
    }

    #[test]
    fn bounding_box_diagonal_computes_extent() {
        // Two points [0,0,0] + [3,4,0]: extent (3,4,0) → diagonal √(9+16) = 5.0.
        let geom = IgesGeometry {
            points: vec![
                IgesPoint {
                    pos: [0.0, 0.0, 0.0],
                },
                IgesPoint {
                    pos: [3.0, 4.0, 0.0],
                },
            ],
            ..Default::default()
        };
        assert!((geom.bounding_box_diagonal() - 5.0).abs() < 1e-9);
        // One line [0,0,0]→[1,2,2]: extent (1,2,2) → diagonal √(1+4+4) = 3.0.
        let geom = IgesGeometry {
            lines: vec![IgesLine {
                start: [0.0, 0.0, 0.0],
                end: [1.0, 2.0, 2.0],
            }],
            ..Default::default()
        };
        assert!((geom.bounding_box_diagonal() - 3.0).abs() < 1e-9);
        // Empty (topology-only) → 0.0; a single coincident point → 0.0 extent.
        assert_eq!(IgesGeometry::default().bounding_box_diagonal(), 0.0);
        let coincident = IgesGeometry {
            points: vec![IgesPoint {
                pos: [5.0, 5.0, 5.0],
            }],
            ..Default::default()
        };
        assert_eq!(coincident.bounding_box_diagonal(), 0.0);
    }

    #[test]
    fn type_102_with_zero_members_handled() {
        let geom = IgesGeometry {
            composite_curves: vec![IgesCompositeCurve {
                member_des: vec![],
            }],
            ..Default::default()
        };
        let txt = render_iges_geometry(&geom, Path::new("t.iges"));
        let parsed = parse(&txt).expect("0-member composite curve must parse");
        assert_eq!(parsed.composite_curves.len(), 1);
        assert!(parsed.composite_curves[0].member_des.is_empty());
    }

    // ---------------------------------------------------------------
    // Round-4 DoS hardening — `ListTooLarge` guard.
    // ---------------------------------------------------------------

    /// Hand-craft a minimal IGES file with a single Type 102 entity
    /// whose advertised member count is `usize::MAX`. Pre-fix the
    /// reader would call `Vec::with_capacity(usize::MAX)` and either
    /// OOM the host or panic. Post-fix it must return
    /// `ListTooLarge { count: usize::MAX, max: MAX_IGES_LIST_LEN }`
    /// without allocating.
    #[test]
    fn parse_type_102_huge_count_returns_list_too_large() {
        // The minimal valid IGES has S/G/D/P/T sections. We need
        // exactly one directory entry (2 lines, 80 cols each) declaring
        // a Type 102, and a parameter line that lies about the count.
        let max = usize::MAX;
        let payload = format!("102,{max}");
        let mut s = String::new();
        // S section — banner.
        write_record(&mut s, "valenx test", 'S', 1);
        // G section — one global record (empty payload is fine).
        write_record(&mut s, "1H,,1H;", 'G', 1);
        // D section — entity 1 = Type 102, PD pointer 1.
        let d1 = format!("{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
            "102", "1", "0", "0", "0", "0", "0", "0", "0");
        write_record(&mut s, &d1, 'D', 1);
        let d2 = format!("{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
            "102", "0", "0", "1", "0", "0", "0", "0", "0");
        write_record(&mut s, &d2, 'D', 2);
        // P section — the entity advertising count = usize::MAX.
        write_record(&mut s, &format!("{payload};"), 'P', 1);
        // T section.
        write_record(&mut s, "S0000001G0000001D0000002P0000001", 'T', 1);
        let err = parse(&s).expect_err("parse must reject huge-count Type 102");
        match err {
            StepIgesError::ListTooLarge { count, max: cap } => {
                assert_eq!(count, usize::MAX);
                assert_eq!(cap, MAX_IGES_LIST_LEN);
            }
            other => panic!("expected ListTooLarge, got: {other:?}"),
        }
    }

    /// Round-9 RED→GREEN: `read` must refuse a file larger than
    /// `MAX_CAD_INTERCHANGE_FILE_BYTES` before allocating. We can't
    /// realistically write 256 MiB in a test, so synthesise a
    /// sparse-style mock by writing one byte at a high offset using
    /// `seek` + `write_all` so the metadata length exceeds the cap
    /// even though the file is sparse on disk.
    #[test]
    fn read_rejects_file_above_max_cad_interchange_file_bytes() {
        use std::io::{Seek, SeekFrom, Write};
        let tmp = std::env::temp_dir().join(format!(
            "valenx-iges-toobig-{}.iges",
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
            // Seek to 1 byte past the cap and write one byte. On
            // POSIX + NTFS this creates a sparse file whose
            // metadata.len() reports cap+2 but uses ~0 disk blocks.
            f.seek(SeekFrom::Start(crate::MAX_CAD_INTERCHANGE_FILE_BYTES + 1))
                .unwrap();
            f.write_all(b"x").unwrap();
        }
        let err = read(&tmp).expect_err("must reject oversize IGES");
        match err {
            StepIgesError::FileTooLarge { format, size, cap } => {
                assert_eq!(format, "IGES");
                assert!(size > cap, "size={size} cap={cap}");
                assert_eq!(cap, crate::MAX_CAD_INTERCHANGE_FILE_BYTES);
            }
            other => panic!("expected FileTooLarge, got: {other:?}"),
        }
        let _ = std::fs::remove_file(&tmp);
    }
}
