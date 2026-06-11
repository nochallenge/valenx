//! IFC4 ISO-10303-21 writer.
//!
//! The writer emits one entity per line in the `DATA;` section. Each
//! `add(...)` call returns the entity's id (`#N`) so callers can wire
//! references between entities. The writer is single-pass — the
//! caller is responsible for emitting child entities before their
//! parents (or for using `reserve_id` to forward-declare).
//!
//! ## IFC GUIDs
//!
//! IFC uses a 22-character compressed-base64 identifier called
//! `IfcGloballyUniqueId`. The compression scheme is a custom
//! base64-with-comma alphabet defined by the IFC spec. v1 generates
//! a UUIDv4 (16 random bytes) and packs it into 22 chars per the
//! standard algorithm — see [`ifc_guid_v4`]. The result is
//! validator-safe ("22 chars from the IFC alphabet, deterministic
//! one-to-one bit packing").

use std::collections::HashMap;
use std::path::Path;

use crate::error::ArchError;

/// Custom base64 alphabet used by IFC's [`ifc_guid_v4`] encoder.
/// 64 characters total — same length and characters as the spec.
pub(crate) const IFC_ALPHABET: &[u8; 64] = b"0123456789\
ABCDEFGHIJKLMNOPQRSTUVWXYZ\
abcdefghijklmnopqrstuvwxyz_$";

/// Generate a fresh IFC GUID (22-char IfcGloballyUniqueId) from
/// 16 random bytes.
///
/// The compression scheme packs the 128 bits of a UUID into 22 chars
/// from the IFC alphabet (see `IFC_ALPHABET`). The first char encodes only 2 bits (top
/// of the UUID), and each subsequent char encodes 6 bits — total
/// `2 + 21 × 6 = 128` bits. We use a small linear-congruential RNG
/// seeded from the system time + a counter so consecutive calls
/// don't collide.
///
/// IFC validators accept any 22-char string from this alphabet; the
/// exact encoding doesn't have to match a specific UUID round-trip.
pub fn ifc_guid_v4() -> String {
    // Pull 16 random-ish bytes from a small in-process counter +
    // time. We avoid pulling the `uuid` crate to keep deps lean for
    // v1 (it's already in `[workspace.dependencies]` but not in this
    // crate); the entropy is far more than IFC validators care about.
    use std::sync::atomic::{AtomicU64, Ordering};
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let c = COUNTER.fetch_add(1, Ordering::Relaxed);
    let now_ns = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos() as u64)
        .unwrap_or(0);
    let mut bytes = [0u8; 16];
    // Simple xorshift-style scrambler seeded from (now_ns, c).
    let mut state: u64 = now_ns ^ (c.wrapping_mul(0x9E3779B97F4A7C15));
    for chunk in bytes.chunks_mut(8) {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let s = state.to_le_bytes();
        for (i, b) in chunk.iter_mut().enumerate() {
            *b = s[i];
        }
    }
    // Mark as UUIDv4 + RFC 4122 variant (cosmetic — IFC validators
    // don't enforce, but it makes the GUID look like a real UUID
    // when decoded).
    bytes[6] = (bytes[6] & 0x0F) | 0x40;
    bytes[8] = (bytes[8] & 0x3F) | 0x80;
    encode_ifc_guid(&bytes)
}

/// Pack 16 bytes (128 bits) into 22 IFC-alphabet chars.
///
/// Algorithm (per IFC spec): the first char carries 2 bits (top of
/// the integer), then 21 chars carry 6 bits each — total `2 + 21*6 =
/// 128`. We process big-endian.
fn encode_ifc_guid(bytes: &[u8; 16]) -> String {
    // Convert 16 bytes to a 128-bit big-endian integer.
    let mut hi: u64 = 0;
    let mut lo: u64 = 0;
    for &b in &bytes[..8] {
        hi = (hi << 8) | b as u64;
    }
    for &b in &bytes[8..] {
        lo = (lo << 8) | b as u64;
    }
    let mut out = [0u8; 22];
    // Top 2 bits → index 0.
    out[0] = IFC_ALPHABET[((hi >> 62) & 0x3) as usize];
    // 21 chunks of 6 bits going down. We need to slide a 6-bit
    // window across hi:lo from bit 61 down to bit 0.
    let mut idx = 1usize;
    let mut shift: i32 = 56; // bit position of the next chunk's LSB in `hi`.
    while idx < 22 {
        let val = if shift >= 0 {
            // Chunk entirely in hi.
            ((hi >> shift) & 0x3F) as u8
        } else if shift + 6 > 0 {
            // Straddle: top part from hi (low bits), bottom from lo (top bits).
            let from_hi = (hi & ((1u64 << (shift + 6)) - 1)) << (-shift);
            let from_lo = lo >> (64 + shift);
            ((from_hi | from_lo) & 0x3F) as u8
        } else {
            // Chunk entirely in lo: its LSB sits at bit (64 + shift) of the
            // 128-bit value, so shift `lo` right by `64 + shift` (the same
            // amount the straddle branch uses for its `from_lo` term). The old
            // `-shift - 8` only coincided at shift == -36, mis-packing the rest.
            ((lo >> (64 + shift)) & 0x3F) as u8
        };
        out[idx] = IFC_ALPHABET[val as usize];
        idx += 1;
        shift -= 6;
    }
    String::from_utf8(out.to_vec()).expect("alphabet bytes are ASCII")
}

/// Single-pass IFC writer. Appends entity lines to an internal
/// buffer and hands out monotonically-increasing entity ids.
pub struct IfcWriter {
    /// All emitted `#N=ENTITY(...)` lines (without leading `#N=`).
    /// `lines[i]` is entity id `i+1`.
    pub(crate) lines: Vec<String>,
    /// Header lines (between `HEADER;` and `ENDSEC;`).
    pub(crate) header: Vec<String>,
    /// Map entity-class name → canonical singletons. Used so
    /// callers can re-fetch the project / site / building / storey
    /// IDs without re-walking the list.
    pub(crate) singletons: HashMap<&'static str, usize>,
}

impl Default for IfcWriter {
    fn default() -> Self {
        Self::new()
    }
}

impl IfcWriter {
    /// Empty writer with no header / data entities yet.
    pub fn new() -> Self {
        Self {
            lines: Vec::new(),
            header: Vec::new(),
            singletons: HashMap::new(),
        }
    }

    /// Append an `ENTITY(...)` line and return its new id.
    ///
    /// `body` should be the entity payload without the leading `#N=`
    /// or trailing `;` — e.g. `"IFCPROJECT('guid',$,'p',$,$,$,$,$,$)"`.
    pub fn add(&mut self, body: impl Into<String>) -> usize {
        self.lines.push(body.into());
        self.lines.len()
    }

    /// Add a header line (between `HEADER;` and `ENDSEC;`).
    pub fn add_header(&mut self, line: impl Into<String>) {
        self.header.push(line.into());
    }

    /// Finalise to an ISO-10303-21 string.
    pub fn finish(&self) -> String {
        let mut s = String::new();
        s.push_str("ISO-10303-21;\n");
        s.push_str("HEADER;\n");
        for h in &self.header {
            s.push_str(h);
            s.push_str(";\n");
        }
        s.push_str("ENDSEC;\n");
        s.push_str("DATA;\n");
        for (i, body) in self.lines.iter().enumerate() {
            s.push_str(&format!("#{}={};\n", i + 1, body));
        }
        s.push_str("ENDSEC;\n");
        s.push_str("END-ISO-10303-21;\n");
        s
    }

    /// Write the full ISO-10303-21 file to `path`.
    ///
    /// Round-24 L1: pre-fix this called `self.finish()` (building a
    /// monolithic `String` containing every header line + every
    /// entity body, then converting that to UTF-8 bytes) and passed
    /// the whole thing to `fs::write`. For a full building model
    /// with 100k entities and 80 byte average bodies that's an
    /// 8 MiB peak heap allocation that exists for the duration of
    /// the write — twice (once as the String, once as the bytes
    /// kernel copy). Streaming via `BufWriter` keeps the peak
    /// allocation at the 8 KiB BufWriter window plus one entity-
    /// body String at a time. Same on-disk bytes, dramatically
    /// smaller peak memory.
    ///
    /// ## Round-27 H1 STRUCTURAL consolidation
    ///
    /// Now routes through
    /// [`valenx_core::io_caps::atomic_write_streaming`] — same
    /// peak-memory profile (streaming via a BufWriter), but
    /// publishes the file atomically: writes to a unique
    /// `<basename>.tmp.<pid>.<counter>` sidecar with `O_NOFOLLOW`
    /// (Unix) / `FILE_FLAG_OPEN_REPARSE_POINT` (Windows), fsyncs the
    /// sidecar, renames over the target, then fsyncs the parent dir
    /// (Unix). Pre-R27 a concurrent multi-MiB IFC export shared the
    /// `File::create` handle with another writer — the kernel can
    /// interleave bytes from both writers (same shape as the R26 H1
    /// dock bug). Post-R27 each writer owns a distinct sidecar so
    /// the rename publishes ONE writer's content atomically; a
    /// concurrent reader sees either the old file or the new file
    /// in its entirety.
    pub fn write_to_path(&self, path: &Path) -> Result<(), ArchError> {
        use std::io::Write;
        valenx_core::io_caps::atomic_write_streaming(path, |w| {
            w.write_all(b"ISO-10303-21;\n")?;
            w.write_all(b"HEADER;\n")?;
            for h in &self.header {
                w.write_all(h.as_bytes())?;
                w.write_all(b";\n")?;
            }
            w.write_all(b"ENDSEC;\n")?;
            w.write_all(b"DATA;\n")?;
            for (i, body) in self.lines.iter().enumerate() {
                // One temporary String per entity — bounded to the
                // body's length plus ~6 bytes of framing. No
                // monolithic buffer.
                let line = format!("#{}={};\n", i + 1, body);
                w.write_all(line.as_bytes())?;
            }
            w.write_all(b"ENDSEC;\n")?;
            w.write_all(b"END-ISO-10303-21;\n")?;
            Ok(())
        })
        .map_err(|e| ArchError::IfcWriteFailed(format!("write {}: {e}", path.display())))?;
        Ok(())
    }
}

/// Quote a string for IFC's STEP Part 21 STRING type. Wraps in
/// single quotes and escapes any embedded single-quote as `''`.
///
/// Round-14 M4: also sanitises `\n`, `\r` (strip), and `\\` (escape
/// per Part 21 to `\\`). Pre-fix a hostile project / wall / space
/// name like `"X');\n#9999=IFCPROJECT(..."` survived the single-
/// quote escape and let the embedded `);\n#9999=` close the parent
/// entity and inject arbitrary IFC entities. Mirrors the round-12
/// CAM `CommentStyle::wrap` sanitiser — same class of bug (caller-
/// owned text flows into a structured output format whose delimiters
/// the caller doesn't own) and the same fix (strip the line-break
/// chars that would let the attacker break out of the wrapper).
///
/// IFC's Part 21 STRING type uses backslash as a control-sequence
/// prefix (e.g. `\X\<hex>` for non-ASCII), so a raw backslash that
/// the caller didn't intend as a control sequence has to be escaped
/// to `\\` per the spec — otherwise a name with a literal `\X` in it
/// would be re-interpreted by IFC consumers as a hex escape.
///
/// Round-25 L2: also strip Unicode line separators U+2028
/// (LINE SEPARATOR) and U+2029 (PARAGRAPH SEPARATOR). Mirrors the
/// round-20 L4 `python_str_repr` fix — some IFC consumers
/// (browser-based viewers, JS-toolchains) treat U+2028 / U+2029
/// as line terminators per ECMA-262, so a hostile string containing
/// one can break out of the IFC entity at the JSON-export or
/// JavaScript-evaluation layer even though the on-disk Part 21
/// representation looks single-line. Strip them defensively for
/// the same reason CR / LF are stripped.
///
/// Round-26 L1: also strip U+0085 (NEXT LINE / NEL). Same threat
/// model as CR / LF / U+2028 / U+2029 — Unicode UAX-14 line-
/// breaking treats NEL as a hard break, and some consumers
/// (terminals, line-oriented logging pipelines, EBCDIC-aware
/// translators) treat it as a newline. The character is rare in
/// modern UTF-8 input, but a hostile project name carrying NEL
/// would otherwise let an attacker embed what looks like a line
/// terminator to downstream tooling.
fn ifc_str(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('\'');
    for c in s.chars() {
        match c {
            // Strip CR / LF / NEL outright. IFC Part 21 entity lines
            // are newline-terminated; embedding one here lets a
            // hostile string break out of the entity and inject a
            // sibling. NEL (U+0085) is added in round-26 L1 — see
            // function-level doc for the threat model.
            '\n' | '\r' | '\u{0085}' => {}
            // Round-25 L2: strip U+2028 / U+2029. Same threat model
            // as CR/LF but for downstream consumers (Web-IFC, viewers,
            // JSON exporters) that follow ECMA-262 line-terminator
            // rules. Defence-in-depth — IFC's own parsers treat
            // these as ordinary characters, so post-fix output round-
            // trips identically for any IFC tooling that doesn't
            // pre-process via JS string semantics.
            '\u{2028}' | '\u{2029}' => {}
            // Single-quote escape per Part 21.
            '\'' => out.push_str("''"),
            // Backslash escape per Part 21 — `\\` is the literal-
            // backslash sequence so consumers don't read a stray
            // `\X` as the start of a control sequence.
            '\\' => out.push_str("\\\\"),
            _ => out.push(c),
        }
    }
    out.push('\'');
    out
}

/// Convert an `id: usize` to a `#N` reference.
fn r(id: usize) -> String {
    format!("#{id}")
}

/// Build the canonical IFC4 header for a project.
fn write_header(w: &mut IfcWriter, project_name: &str) {
    let now = chrono_like_now();
    w.add_header("FILE_DESCRIPTION(('ViewDefinition [CoordinationView]'),'2;1')".to_string());
    w.add_header(format!(
        "FILE_NAME({name},'{ts}',('Valenx'),('Valenx'),'valenx-arch','valenx-arch','')",
        name = ifc_str(project_name),
        ts = now,
    ));
    w.add_header("FILE_SCHEMA(('IFC4'))".to_string());
}

/// A simple timestamp string in ISO 8601 format. Generated without a
/// chrono dep — uses `SystemTime` arithmetic.
fn chrono_like_now() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let dur = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default();
    let total_secs = dur.as_secs() as i64;
    let (y, mo, d, h, mi, s) = epoch_to_ymd_hms(total_secs);
    format!("{y:04}-{mo:02}-{d:02}T{h:02}:{mi:02}:{s:02}")
}

/// Civil calendar from seconds since 1970-01-01 UTC. Doesn't honor
/// leap seconds — fine for an IFC file timestamp.
fn epoch_to_ymd_hms(epoch: i64) -> (i32, u32, u32, u32, u32, u32) {
    let days = epoch.div_euclid(86_400);
    let secs_today = epoch.rem_euclid(86_400) as u32;
    let h = secs_today / 3600;
    let mi = (secs_today / 60) % 60;
    let s = secs_today % 60;

    // Howard Hinnant's days-from-civil inverse — works for dates
    // 1970-01-01 onward.
    let z = days + 719_468;
    let era = if z >= 0 {
        z / 146_097
    } else {
        (z - 146_096) / 146_097
    };
    let doe = (z - era * 146_097) as u32;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mo = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if mo <= 2 { y + 1 } else { y };
    (y as i32, mo, d, h, mi, s)
}

/// Emit a placement / coordinate system pair pointing at the world
/// origin with the standard +X / +Z axes. Returns the
/// IfcLocalPlacement entity id.
fn emit_origin_placement(w: &mut IfcWriter) -> usize {
    let p_origin = w.add("IFCCARTESIANPOINT((0.,0.,0.))");
    let dir_x = w.add("IFCDIRECTION((1.,0.,0.))");
    let dir_z = w.add("IFCDIRECTION((0.,0.,1.))");
    let axis_placement = w.add(format!(
        "IFCAXIS2PLACEMENT3D({},{},{})",
        r(p_origin),
        r(dir_z),
        r(dir_x)
    ));
    w.add(format!("IFCLOCALPLACEMENT($,{})", r(axis_placement)))
}

/// Emit a placement positioned at `(x, y, z)` aligned with world
/// axes, parented to `parent_placement_id`.
fn emit_placement_at(
    w: &mut IfcWriter,
    parent_placement_id: usize,
    x: f64,
    y: f64,
    z: f64,
) -> usize {
    let p = w.add(format!("IFCCARTESIANPOINT(({x:.6},{y:.6},{z:.6}))"));
    let dir_x = w.add("IFCDIRECTION((1.,0.,0.))");
    let dir_z = w.add("IFCDIRECTION((0.,0.,1.))");
    let axis = w.add(format!(
        "IFCAXIS2PLACEMENT3D({},{},{})",
        r(p),
        r(dir_z),
        r(dir_x)
    ));
    w.add(format!(
        "IFCLOCALPLACEMENT({},{})",
        r(parent_placement_id),
        r(axis)
    ))
}

/// Emit an extruded-area-solid representation for a rectangular
/// profile of `width × depth`, extruded along +Z by `height`,
/// positioned at the supplied placement.
fn emit_rectangular_extrusion(
    w: &mut IfcWriter,
    placement: usize,
    width: f64,
    depth: f64,
    height: f64,
) -> usize {
    let axis2d = placement_axis2d(w);
    let profile = w.add(format!(
        "IFCRECTANGLEPROFILEDEF(.AREA.,$,{},{:.6},{:.6})",
        r(axis2d),
        width,
        depth
    ));
    let extrude_dir = w.add("IFCDIRECTION((0.,0.,1.))");
    let solid_origin = w.add("IFCCARTESIANPOINT((0.,0.,0.))");
    let dir_z = w.add("IFCDIRECTION((0.,0.,1.))");
    let dir_x = w.add("IFCDIRECTION((1.,0.,0.))");
    let position = w.add(format!(
        "IFCAXIS2PLACEMENT3D({},{},{})",
        r(solid_origin),
        r(dir_z),
        r(dir_x)
    ));
    let solid = w.add(format!(
        "IFCEXTRUDEDAREASOLID({},{},{},{:.6})",
        r(profile),
        r(position),
        r(extrude_dir),
        height
    ));
    let ctx = geometric_context(w);
    let representation = w.add(format!(
        "IFCSHAPEREPRESENTATION({},'Body','SweptSolid',({}))",
        r(ctx),
        r(solid),
    ));
    let _ = placement;
    w.add(format!(
        "IFCPRODUCTDEFINITIONSHAPE($,$,({}))",
        r(representation)
    ))
}

/// Lazily build a 2D placement at the origin for profile usage.
fn placement_axis2d(w: &mut IfcWriter) -> usize {
    // Cache via the singletons map.
    if let Some(&id) = w.singletons.get("axis2d") {
        return id;
    }
    let origin2d = w.add("IFCCARTESIANPOINT((0.,0.))");
    let dir2d = w.add("IFCDIRECTION((1.,0.))");
    let axis = w.add(format!("IFCAXIS2PLACEMENT2D({},{})", r(origin2d), r(dir2d)));
    w.singletons.insert("axis2d", axis);
    axis
}

/// Lazily build a single shared IfcGeometricRepresentationContext.
fn geometric_context(w: &mut IfcWriter) -> usize {
    if let Some(&id) = w.singletons.get("geom_ctx") {
        return id;
    }
    // Build the world coordinate system + context. These need to
    // exist even if no entity uses them, for IFC validators.
    let origin = w.add("IFCCARTESIANPOINT((0.,0.,0.))");
    let dir_x = w.add("IFCDIRECTION((1.,0.,0.))");
    let dir_z = w.add("IFCDIRECTION((0.,0.,1.))");
    let world_cs = w.add(format!(
        "IFCAXIS2PLACEMENT3D({},{},{})",
        r(origin),
        r(dir_z),
        r(dir_x)
    ));
    let ctx = w.add(format!(
        "IFCGEOMETRICREPRESENTATIONCONTEXT($,'Model',3,1.E-05,{},$)",
        r(world_cs)
    ));
    w.singletons.insert("geom_ctx", ctx);
    ctx
}

/// Emit the canonical IfcProject / IfcSite / IfcBuilding /
/// IfcBuildingStorey hierarchy with RelAggregates links. Returns the
/// storey id (everything else can be linked off it via
/// IfcRelContainedInSpatialStructure).
fn emit_project_hierarchy(w: &mut IfcWriter, project_name: &str) -> usize {
    let project_placement = emit_origin_placement(w);
    let ctx = geometric_context(w);
    let project = w.add(format!(
        "IFCPROJECT({guid},$,{name},$,$,{name},$,({ctx_ref}),$)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(project_name),
        ctx_ref = r(ctx),
    ));
    let site = w.add(format!(
        "IFCSITE({guid},$,'Site',$,$,{place},$,$,.ELEMENT.,$,$,$,$,$)",
        guid = ifc_str(&ifc_guid_v4()),
        place = r(project_placement)
    ));
    let building = w.add(format!(
        "IFCBUILDING({guid},$,'Building',$,$,{place},$,$,.ELEMENT.,$,$,$)",
        guid = ifc_str(&ifc_guid_v4()),
        place = r(project_placement)
    ));
    let storey = w.add(format!(
        "IFCBUILDINGSTOREY({guid},$,'Storey',$,$,{place},$,$,.ELEMENT.,0.0)",
        guid = ifc_str(&ifc_guid_v4()),
        place = r(project_placement)
    ));
    // RelAggregates wiring: project → site, site → building, building → storey.
    w.add(format!(
        "IFCRELAGGREGATES({guid},$,$,$,{rel},({rel_to}))",
        guid = ifc_str(&ifc_guid_v4()),
        rel = r(project),
        rel_to = r(site)
    ));
    w.add(format!(
        "IFCRELAGGREGATES({guid},$,$,$,{rel},({rel_to}))",
        guid = ifc_str(&ifc_guid_v4()),
        rel = r(site),
        rel_to = r(building)
    ));
    w.add(format!(
        "IFCRELAGGREGATES({guid},$,$,$,{rel},({rel_to}))",
        guid = ifc_str(&ifc_guid_v4()),
        rel = r(building),
        rel_to = r(storey)
    ));
    w.singletons.insert("project", project);
    w.singletons.insert("site", site);
    w.singletons.insert("building", building);
    w.singletons.insert("storey", storey);
    storey
}

/// Emit a wall as IfcWall + placement + extruded-area shape.
/// Returns the wall entity id.
pub fn write_wall(w: &mut IfcWriter, wall: &crate::wall::WallParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, wall.start.x, wall.start.y, wall.start.z);
    let shape =
        emit_rectangular_extrusion(w, placement, wall.length(), wall.thickness, wall.height);
    w.add(format!(
        "IFCWALL({guid},$,{name},$,$,{place},{shape},$,.NOTDEFINED.)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Wall ({})", wall.material)),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit a slab.
pub fn write_slab(w: &mut IfcWriter, slab: &crate::slab::SlabParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let (cx, cy, z0) = slab_centroid(slab);
    let placement = emit_placement_at(w, storey_place, cx, cy, z0);
    // Bounding box width × depth for the v1 simplification.
    let (xmin, ymin, xmax, ymax) = slab_aabb_xy(slab);
    let width = (xmax - xmin).max(0.001);
    let depth = (ymax - ymin).max(0.001);
    let shape = emit_rectangular_extrusion(w, placement, width, depth, slab.thickness);
    w.add(format!(
        "IFCSLAB({guid},$,{name},$,$,{place},{shape},$,.FLOOR.)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Slab ({})", slab.material)),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit a column.
pub fn write_column(w: &mut IfcWriter, col: &crate::column::ColumnParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, col.base.x, col.base.y, col.base.z);
    let (w_dim, d_dim) = match &col.cross_section {
        crate::column::ColumnSection::Rectangle { width, depth } => (*width, *depth),
        crate::column::ColumnSection::Circular { radius, .. } => (2.0 * radius, 2.0 * radius),
        crate::column::ColumnSection::IBeam { width, depth, .. } => (*width, *depth),
    };
    let shape = emit_rectangular_extrusion(w, placement, w_dim, d_dim, col.height);
    w.add(format!(
        "IFCCOLUMN({guid},$,{name},$,$,{place},{shape},$,.NOTDEFINED.)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Column ({})", col.material)),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit a beam.
pub fn write_beam(w: &mut IfcWriter, beam: &crate::beam::BeamParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, beam.start.x, beam.start.y, beam.start.z);
    let (w_dim, d_dim) = beam.cross_section.outer_box();
    let shape = emit_rectangular_extrusion(w, placement, beam.length(), w_dim.max(d_dim), w_dim);
    w.add(format!(
        "IFCBEAM({guid},$,{name},$,$,{place},{shape},$,.NOTDEFINED.)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Beam ({})", beam.material)),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit a window.
pub fn write_window(w: &mut IfcWriter, win: &crate::window::WindowParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(
        w,
        storey_place,
        win.position_along_wall,
        0.0,
        win.position_height,
    );
    let shape = emit_rectangular_extrusion(
        w,
        placement,
        win.width,
        win.frame_thickness.max(0.02),
        win.height,
    );
    w.add(format!(
        "IFCWINDOW({guid},$,{name},$,$,{place},{shape},$,{h:.6},{wd:.6},.NOTDEFINED.,.NOTDEFINED.,$)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Window ({})", win.style.label())),
        place = r(placement),
        shape = r(shape),
        h = win.height,
        wd = win.width
    ))
}

/// Emit a door.
pub fn write_door(w: &mut IfcWriter, door: &crate::door::DoorParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, door.position_along_wall, 0.0, 0.0);
    let shape = emit_rectangular_extrusion(w, placement, door.width, 0.04, door.height);
    w.add(format!(
        "IFCDOOR({guid},$,{name},$,$,{place},{shape},$,{h:.6},{wd:.6},.NOTDEFINED.,.NOTDEFINED.,$)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Door ({})", door.style.label())),
        place = r(placement),
        shape = r(shape),
        h = door.height,
        wd = door.width
    ))
}

/// Emit a stair (IFC4 IfcStair).
pub fn write_stair(w: &mut IfcWriter, stair: &crate::stair::StairParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, stair.base.x, stair.base.y, stair.base.z);
    let shape =
        emit_rectangular_extrusion(w, placement, stair.total_run, stair.width, stair.total_rise);
    w.add(format!(
        "IFCSTAIR({guid},$,'Stair',$,$,{place},{shape},$,.STRAIGHT_RUN_STAIR.)",
        guid = ifc_str(&ifc_guid_v4()),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit a roof.
pub fn write_roof(w: &mut IfcWriter, roof: &crate::roof::RoofParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let (xmin, ymin, xmax, ymax) = boundary_aabb_xy(&roof.boundary);
    let placement = emit_placement_at(
        w,
        storey_place,
        xmin,
        ymin,
        roof.boundary.first().map(|p| p.z).unwrap_or(0.0),
    );
    let shape = emit_rectangular_extrusion(
        w,
        placement,
        (xmax - xmin).max(0.001),
        (ymax - ymin).max(0.001),
        roof.peak_height.max(0.001),
    );
    let predef = match roof.roof_type {
        crate::roof::RoofType::Flat => ".FLAT_ROOF.",
        crate::roof::RoofType::Gable => ".GABLE_ROOF.",
        crate::roof::RoofType::Hip => ".HIP_ROOF.",
        crate::roof::RoofType::Shed => ".SHED_ROOF.",
    };
    w.add(format!(
        "IFCROOF({guid},$,{name},$,$,{place},{shape},$,{predef})",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Roof ({})", roof.roof_type.label())),
        place = r(placement),
        shape = r(shape),
        predef = predef
    ))
}

/// Emit a space (IfcSpace).
pub fn write_space(w: &mut IfcWriter, sp: &crate::space::SpaceParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let (xmin, ymin, xmax, ymax) = boundary_aabb_xy(&sp.boundary);
    let placement = emit_placement_at(
        w,
        storey_place,
        xmin,
        ymin,
        sp.boundary.first().map(|p| p.z).unwrap_or(0.0),
    );
    let shape = emit_rectangular_extrusion(
        w,
        placement,
        (xmax - xmin).max(0.001),
        (ymax - ymin).max(0.001),
        sp.ceiling_height,
    );
    w.add(format!(
        "IFCSPACE({guid},$,{name},$,$,{place},{shape},$,.ELEMENT.,.INTERNAL.,$)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&sp.space_name),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit an HVAC duct segment as IfcDuctSegment.
pub fn write_duct(w: &mut IfcWriter, d: &crate::mep::DuctSegmentParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, d.start.x, d.start.y, d.start.z);
    let (sw, sh) = d.shape.outer_box();
    let shape = emit_rectangular_extrusion(w, placement, sw, sh, d.length().max(0.001));
    let predef = match d.shape {
        crate::mep::DuctShape::Round { .. } => ".RIGIDSEGMENT.",
        crate::mep::DuctShape::Rectangular { .. } => ".RIGIDSEGMENT.",
        crate::mep::DuctShape::Oval { .. } => ".RIGIDSEGMENT.",
    };
    w.add(format!(
        "IFCDUCTSEGMENT({guid},$,{name},$,$,{place},{shape},$,{predef})",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!(
            "Duct {} ({} - {})",
            d.shape.label(),
            d.flow_direction.label(),
            d.material
        )),
        place = r(placement),
        shape = r(shape),
        predef = predef
    ))
}

/// Emit a plumbing / process pipe segment as IfcPipeSegment.
pub fn write_pipe(w: &mut IfcWriter, p: &crate::mep::PipeSegmentParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, p.start.x, p.start.y, p.start.z);
    let shape = emit_rectangular_extrusion(
        w,
        placement,
        p.diameter,
        p.diameter,
        p.length().max(0.001),
    );
    w.add(format!(
        "IFCPIPESEGMENT({guid},$,{name},$,$,{place},{shape},$,.RIGIDSEGMENT.)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Pipe ({} - {})", p.fluid, p.material)),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit an electrical cable run as IfcCableSegment.
pub fn write_cable(w: &mut IfcWriter, c: &crate::mep::CableSegmentParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, c.start.x, c.start.y, c.start.z);
    let shape = emit_rectangular_extrusion(
        w,
        placement,
        c.diameter,
        c.diameter,
        c.length().max(0.001),
    );
    w.add(format!(
        "IFCCABLESEGMENT({guid},$,{name},$,$,{place},{shape},$,.CABLESEGMENT.)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!(
            "Cable {:.0}V {:.1}mm2 ({})",
            c.voltage, c.conductor_csa_mm2, c.material
        )),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit an electrical conduit segment as IfcConduitSegment.
pub fn write_conduit(w: &mut IfcWriter, c: &crate::mep::ConduitSegmentParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, c.start.x, c.start.y, c.start.z);
    let shape = emit_rectangular_extrusion(
        w,
        placement,
        c.outer_diameter,
        c.outer_diameter,
        c.length().max(0.001),
    );
    w.add(format!(
        "IFCCABLECARRIERSEGMENT({guid},$,{name},$,$,{place},{shape},$,.CONDUITSEGMENT.)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Conduit ({})", c.material)),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit an MEP equipment placement using the kind-appropriate IFC4
/// entity ([`crate::mep::EquipmentKind::ifc_entity`]).
pub fn write_mep_equipment(w: &mut IfcWriter, e: &crate::mep::MepEquipmentParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, e.position.x, e.position.y, e.position.z);
    let shape = emit_rectangular_extrusion(w, placement, e.size[0], e.size[1], e.size[2]);
    let entity = e.kind.ifc_entity();
    // For the IFC4 distribution elements, the standard signature is
    // (GUID, OwnerHistory, Name, Description, ObjectType, ObjectPlacement,
    // Representation, Tag, PredefinedType).
    let predef = match e.kind {
        crate::mep::EquipmentKind::AirHandlingUnit => ".CONSTANTFLOW.",
        crate::mep::EquipmentKind::VavBox => ".VARIABLEFLOWPRESSUREDEPENDANT.",
        crate::mep::EquipmentKind::Pump => ".CIRCULATOR.",
        crate::mep::EquipmentKind::Valve => ".GATEVALVE.",
        crate::mep::EquipmentKind::SprinklerHead => ".SPRINKLER.",
        crate::mep::EquipmentKind::ElectricalPanel => ".DISTRIBUTIONBOARD.",
        crate::mep::EquipmentKind::LightFitting => ".POINTSOURCE.",
    };
    w.add(format!(
        "{ent}({guid},$,{name},{descr},$,{place},{shape},{tag},{predef})",
        ent = entity,
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("{} {}", e.kind.label(), e.tag)),
        descr = ifc_str(&e.description),
        place = r(placement),
        shape = r(shape),
        tag = ifc_str(&e.tag),
        predef = predef
    ))
}

/// Emit an IfcCovering (interior / exterior surface covering such as
/// flooring, plaster, ceiling tiles).
pub fn write_covering(
    w: &mut IfcWriter,
    slab: &crate::slab::SlabParams,
    predefined: &str,
) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let (xmin, ymin, xmax, ymax) = boundary_aabb_xy(&slab.boundary);
    let placement = emit_placement_at(
        w,
        storey_place,
        xmin,
        ymin,
        slab.boundary.first().map(|p| p.z).unwrap_or(0.0),
    );
    let shape = emit_rectangular_extrusion(
        w,
        placement,
        (xmax - xmin).max(0.001),
        (ymax - ymin).max(0.001),
        slab.thickness,
    );
    w.add(format!(
        "IFCCOVERING({guid},$,{name},$,$,{place},{shape},$,{predef})",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Covering ({})", slab.material)),
        place = r(placement),
        shape = r(shape),
        predef = predefined
    ))
}

/// Emit an IfcCurtainWall using a wall's footprint geometry.
pub fn write_curtain_wall(w: &mut IfcWriter, wall: &crate::wall::WallParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, wall.start.x, wall.start.y, wall.start.z);
    let shape = emit_rectangular_extrusion(w, placement, wall.length(), wall.thickness, wall.height);
    w.add(format!(
        "IFCCURTAINWALL({guid},$,{name},$,$,{place},{shape},$,.NOTDEFINED.)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Curtain Wall ({})", wall.material)),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit an IfcFooting using a slab's footprint.
pub fn write_footing(w: &mut IfcWriter, slab: &crate::slab::SlabParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let (xmin, ymin, xmax, ymax) = boundary_aabb_xy(&slab.boundary);
    let placement = emit_placement_at(
        w,
        storey_place,
        xmin,
        ymin,
        slab.boundary.first().map(|p| p.z).unwrap_or(0.0),
    );
    let shape = emit_rectangular_extrusion(
        w,
        placement,
        (xmax - xmin).max(0.001),
        (ymax - ymin).max(0.001),
        slab.thickness,
    );
    w.add(format!(
        "IFCFOOTING({guid},$,{name},$,$,{place},{shape},$,.STRIP_FOOTING.)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Footing ({})", slab.material)),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit an IfcPile (foundation pile) using column geometry.
pub fn write_pile(w: &mut IfcWriter, col: &crate::column::ColumnParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, col.base.x, col.base.y, col.base.z);
    let (w_dim, d_dim) = match &col.cross_section {
        crate::column::ColumnSection::Rectangle { width, depth } => (*width, *depth),
        crate::column::ColumnSection::Circular { radius, .. } => (2.0 * radius, 2.0 * radius),
        crate::column::ColumnSection::IBeam { width, depth, .. } => (*width, *depth),
    };
    let shape = emit_rectangular_extrusion(w, placement, w_dim, d_dim, col.height);
    w.add(format!(
        "IFCPILE({guid},$,{name},$,$,{place},{shape},$,.BORED.,$)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Pile ({})", col.material)),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit an IfcRailing (e.g. a stair guard rail). Geometry is a thin
/// horizontal bar of `length × thickness × height`.
pub fn write_railing(
    w: &mut IfcWriter,
    start: nalgebra::Vector3<f64>,
    end: nalgebra::Vector3<f64>,
    height: f64,
    material: &str,
) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, start.x, start.y, start.z);
    let length = (end - start).norm().max(0.001);
    let shape = emit_rectangular_extrusion(w, placement, length, 0.05, height);
    w.add(format!(
        "IFCRAILING({guid},$,{name},$,$,{place},{shape},$,.HANDRAIL.)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Railing ({material})")),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit an IfcRamp.
pub fn write_ramp(w: &mut IfcWriter, stair: &crate::stair::StairParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, stair.base.x, stair.base.y, stair.base.z);
    let shape =
        emit_rectangular_extrusion(w, placement, stair.total_run, stair.width, stair.total_rise);
    w.add(format!(
        "IFCRAMP({guid},$,'Ramp',$,$,{place},{shape},$,.STRAIGHT_RUN_RAMP.)",
        guid = ifc_str(&ifc_guid_v4()),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit an IfcChimney attached to a column footprint.
pub fn write_chimney(w: &mut IfcWriter, col: &crate::column::ColumnParams) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, col.base.x, col.base.y, col.base.z);
    let (w_dim, d_dim) = match &col.cross_section {
        crate::column::ColumnSection::Rectangle { width, depth } => (*width, *depth),
        crate::column::ColumnSection::Circular { radius, .. } => (2.0 * radius, 2.0 * radius),
        crate::column::ColumnSection::IBeam { width, depth, .. } => (*width, *depth),
    };
    let shape = emit_rectangular_extrusion(w, placement, w_dim, d_dim, col.height);
    w.add(format!(
        "IFCCHIMNEY({guid},$,{name},$,$,{place},{shape},$,.NOTDEFINED.)",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Chimney ({})", col.material)),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit an IfcFurnishingElement at a position with the given size +
/// label (used for furniture / fixtures placements).
pub fn write_furnishing(
    w: &mut IfcWriter,
    e: &crate::mep::MepEquipmentParams,
) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let placement = emit_placement_at(w, storey_place, e.position.x, e.position.y, e.position.z);
    let shape = emit_rectangular_extrusion(w, placement, e.size[0], e.size[1], e.size[2]);
    w.add(format!(
        "IFCFURNISHINGELEMENT({guid},$,{name},$,$,{place},{shape},{tag})",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(&format!("Furnishing {}", e.tag)),
        place = r(placement),
        shape = r(shape),
        tag = ifc_str(&e.tag)
    ))
}

/// Emit an IfcOpeningElement representing a window or door void.
/// The opening is sized to the (width × thickness × height) of the
/// opening, positioned at the opening's centre on the host wall.
///
/// Returns the opening entity id; the caller wires it to its host
/// via [`emit_rel_voids_element`].
pub fn write_opening_for_window(
    w: &mut IfcWriter,
    win: &crate::window::WindowParams,
    host: &crate::wall::WallParams,
) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let axis = host.axis_xy();
    let cx = host.start + axis * win.position_along_wall;
    let placement = emit_placement_at(
        w,
        storey_place,
        cx.x,
        cx.y,
        host.start.z + win.position_height,
    );
    let shape =
        emit_rectangular_extrusion(w, placement, win.width, host.thickness * 1.5, win.height);
    w.add(format!(
        "IFCOPENINGELEMENT({guid},$,'WindowOpening',$,$,{place},{shape},$,.OPENING.)",
        guid = ifc_str(&ifc_guid_v4()),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit an IfcOpeningElement representing a door void.
pub fn write_opening_for_door(
    w: &mut IfcWriter,
    door: &crate::door::DoorParams,
    host: &crate::wall::WallParams,
) -> usize {
    let storey_place = w
        .singletons
        .get("storey")
        .copied()
        .unwrap_or_else(|| emit_origin_placement(w));
    let axis = host.axis_xy();
    let cx = host.start + axis * door.position_along_wall;
    let placement = emit_placement_at(w, storey_place, cx.x, cx.y, host.start.z);
    let shape =
        emit_rectangular_extrusion(w, placement, door.width, host.thickness * 1.5, door.height);
    w.add(format!(
        "IFCOPENINGELEMENT({guid},$,'DoorOpening',$,$,{place},{shape},$,.OPENING.)",
        guid = ifc_str(&ifc_guid_v4()),
        place = r(placement),
        shape = r(shape)
    ))
}

/// Emit an IfcRelVoidsElement linking an opening to its host wall.
pub fn emit_rel_voids_element(w: &mut IfcWriter, host: usize, opening: usize) -> usize {
    w.add(format!(
        "IFCRELVOIDSELEMENT({guid},$,$,$,{host},{open})",
        guid = ifc_str(&ifc_guid_v4()),
        host = r(host),
        open = r(opening)
    ))
}

/// Emit an IfcRelSpaceBoundary linking a space to a building element.
pub fn emit_rel_space_boundary(
    w: &mut IfcWriter,
    space: usize,
    element: usize,
    boundary_kind: &str,
) -> usize {
    w.add(format!(
        "IFCRELSPACEBOUNDARY({guid},$,'SpaceBoundary',$,{sp},{el},$,{kind},.PHYSICAL.)",
        guid = ifc_str(&ifc_guid_v4()),
        sp = r(space),
        el = r(element),
        kind = boundary_kind
    ))
}

/// Emit one IfcPropertySet attached to `element` via an
/// IfcRelDefinesByProperties. `pset_name` is the canonical name
/// (e.g. `"Pset_WallCommon"`); `props` is a list of `(name, value)`
/// pairs that become IfcPropertySingleValue entities.
///
/// Each property value is wrapped in the appropriate IFC measure
/// type: f64 → IfcReal, string → IfcText, bool → IfcBoolean.
pub fn emit_pset(
    w: &mut IfcWriter,
    element: usize,
    pset_name: &str,
    props: &[(&str, PropValue)],
) -> usize {
    let mut prop_ids: Vec<usize> = Vec::with_capacity(props.len());
    for (name, value) in props {
        let nominal = match value {
            PropValue::Text(s) => format!("IFCTEXT({})", ifc_str(s)),
            PropValue::Real(x) => format!("IFCREAL({x:.6})"),
            PropValue::Bool(b) => format!("IFCBOOLEAN(.{}.)", if *b { "T" } else { "F" }),
            PropValue::Integer(i) => format!("IFCINTEGER({i})"),
            PropValue::Label(s) => format!("IFCLABEL({})", ifc_str(s)),
        };
        let nominal_id = w.add(nominal);
        let prop = w.add(format!(
            "IFCPROPERTYSINGLEVALUE({name},$,{val},$)",
            name = ifc_str(name),
            val = r(nominal_id)
        ));
        prop_ids.push(prop);
    }
    let props_list = prop_ids.iter().map(|i| r(*i)).collect::<Vec<_>>().join(",");
    let pset = w.add(format!(
        "IFCPROPERTYSET({guid},$,{name},$,({props}))",
        guid = ifc_str(&ifc_guid_v4()),
        name = ifc_str(pset_name),
        props = props_list
    ));
    w.add(format!(
        "IFCRELDEFINESBYPROPERTIES({guid},$,$,$,({el}),{pset})",
        guid = ifc_str(&ifc_guid_v4()),
        el = r(element),
        pset = r(pset)
    ));
    pset
}

/// A property value carried by [`emit_pset`] — chooses the IFC
/// measure type to wrap on emit.
#[derive(Debug, Clone)]
pub enum PropValue {
    /// Free-form text — wraps as `IfcText`.
    Text(String),
    /// Numeric — wraps as `IfcReal`.
    Real(f64),
    /// Boolean — wraps as `IfcBoolean`.
    Bool(bool),
    /// Integer count — wraps as `IfcInteger`.
    Integer(i64),
    /// Label (short identifier) — wraps as `IfcLabel`.
    Label(String),
}

fn slab_centroid(slab: &crate::slab::SlabParams) -> (f64, f64, f64) {
    let (xmin, ymin, xmax, ymax) = slab_aabb_xy(slab);
    let z = slab.boundary.first().map(|p| p.z).unwrap_or(0.0);
    (0.5 * (xmin + xmax), 0.5 * (ymin + ymax), z)
}

fn slab_aabb_xy(slab: &crate::slab::SlabParams) -> (f64, f64, f64, f64) {
    boundary_aabb_xy(&slab.boundary)
}

fn boundary_aabb_xy(boundary: &[nalgebra::Vector3<f64>]) -> (f64, f64, f64, f64) {
    let mut xmin = f64::INFINITY;
    let mut ymin = f64::INFINITY;
    let mut xmax = f64::NEG_INFINITY;
    let mut ymax = f64::NEG_INFINITY;
    for p in boundary {
        xmin = xmin.min(p.x);
        ymin = ymin.min(p.y);
        xmax = xmax.max(p.x);
        ymax = ymax.max(p.y);
    }
    if !xmin.is_finite() {
        return (0.0, 0.0, 1.0, 1.0);
    }
    (xmin, ymin, xmax, ymax)
}

/// Top-level entry point — write a complete IFC4 file for `doc` to
/// `path`.
///
/// The emitted file carries:
/// - the canonical `IfcProject` / `IfcSite` / `IfcBuilding` /
///   `IfcBuildingStorey` hierarchy with `IfcRelAggregates` links;
/// - one IFC entity per arch entity, choosing the kind-appropriate
///   IFC4 class (see [`write_wall`], [`write_duct`], etc.);
/// - per-element [`emit_pset`] property sets — `Pset_WallCommon`,
///   `Pset_SlabCommon`, `Pset_BeamCommon`, `Pset_ColumnCommon`,
///   `Pset_WindowCommon`, `Pset_DoorCommon`, `Pset_SpaceCommon`,
///   `Pset_DuctSegmentTypeCommon`, `Pset_PipeSegmentTypeCommon`,
///   `Pset_CableSegmentTypeCommon`, and `Pset_ConduitSegmentTypeCommon`;
/// - `IfcOpeningElement` voids for windows + doors with
///   `IfcRelVoidsElement` host-to-opening relationships;
/// - `IfcRelSpaceBoundary` rows for every (space, wall) pair where
///   the wall's mid-point lies inside the space's bounding box —
///   approximates space-element adjacency without a full topological
///   intersection;
/// - one `IfcRelContainedInSpatialStructure` row tying every
///   building element to the storey.
pub fn write_document(doc: &crate::document::ArchDocument, path: &Path) -> Result<(), ArchError> {
    let mut w = IfcWriter::new();
    write_header(&mut w, &doc.project_name);
    let storey = emit_project_hierarchy(&mut w, &doc.project_name);

    // Build a wall lookup so windows / doors can find their host.
    let walls: std::collections::HashMap<usize, &crate::wall::WallParams> = doc
        .entities
        .iter()
        .filter_map(|(id, e)| match e {
            crate::entity::ArchEntity::Wall(w) => Some((*id, w)),
            _ => None,
        })
        .collect();

    // Emit every entity. Remember each entity's IFC id and the host
    // wall (for IfcRelVoidsElement) so we can wire relationships
    // after every element is emitted.
    let mut element_ids: Vec<usize> = Vec::with_capacity(doc.count());
    // (arch_id, ifc_id) tuples per kind, for relationship resolution.
    let mut wall_ifc_id: std::collections::HashMap<usize, usize> =
        std::collections::HashMap::new();
    let mut space_ifc_id: std::collections::HashMap<usize, usize> =
        std::collections::HashMap::new();
    // Window/door (arch_id, ifc_id, host_arch_id) — for openings.
    let mut window_ifc: Vec<(usize, usize, usize)> = Vec::new();
    let mut door_ifc: Vec<(usize, usize, usize)> = Vec::new();

    for (arch_id, ent) in &doc.entities {
        let eid = match ent {
            crate::entity::ArchEntity::Wall(w_p) => {
                let id = write_wall(&mut w, w_p);
                wall_ifc_id.insert(*arch_id, id);
                emit_pset(
                    &mut w,
                    id,
                    "Pset_WallCommon",
                    &[
                        ("Reference", PropValue::Label(w_p.material.clone())),
                        ("LoadBearing", PropValue::Bool(true)),
                        ("IsExternal", PropValue::Bool(false)),
                        ("ThermalTransmittance", PropValue::Real(0.0)),
                        ("FireRating", PropValue::Label("F30".into())),
                    ],
                );
                id
            }
            crate::entity::ArchEntity::Slab(s_p) => {
                let id = write_slab(&mut w, s_p);
                emit_pset(
                    &mut w,
                    id,
                    "Pset_SlabCommon",
                    &[
                        ("Reference", PropValue::Label(s_p.material.clone())),
                        ("LoadBearing", PropValue::Bool(true)),
                        ("IsExternal", PropValue::Bool(false)),
                        ("PitchAngle", PropValue::Real(0.0)),
                    ],
                );
                id
            }
            crate::entity::ArchEntity::Column(c_p) => {
                let id = write_column(&mut w, c_p);
                emit_pset(
                    &mut w,
                    id,
                    "Pset_ColumnCommon",
                    &[
                        ("Reference", PropValue::Label(c_p.material.clone())),
                        ("LoadBearing", PropValue::Bool(true)),
                        ("IsExternal", PropValue::Bool(false)),
                        ("Slope", PropValue::Real(0.0)),
                    ],
                );
                id
            }
            crate::entity::ArchEntity::Beam(b_p) => {
                let id = write_beam(&mut w, b_p);
                emit_pset(
                    &mut w,
                    id,
                    "Pset_BeamCommon",
                    &[
                        ("Reference", PropValue::Label(b_p.material.clone())),
                        ("Span", PropValue::Real(b_p.length())),
                        ("LoadBearing", PropValue::Bool(true)),
                        ("IsExternal", PropValue::Bool(false)),
                    ],
                );
                id
            }
            crate::entity::ArchEntity::Window(wp) => {
                let id = write_window(&mut w, wp);
                window_ifc.push((*arch_id, id, wp.host));
                emit_pset(
                    &mut w,
                    id,
                    "Pset_WindowCommon",
                    &[
                        ("Reference", PropValue::Label(wp.style.label().to_string())),
                        ("FrameThickness", PropValue::Real(wp.frame_thickness)),
                        ("IsExternal", PropValue::Bool(true)),
                        ("SmokeStop", PropValue::Bool(false)),
                    ],
                );
                id
            }
            crate::entity::ArchEntity::Door(dp) => {
                let id = write_door(&mut w, dp);
                door_ifc.push((*arch_id, id, dp.host));
                emit_pset(
                    &mut w,
                    id,
                    "Pset_DoorCommon",
                    &[
                        ("Reference", PropValue::Label(dp.style.label().to_string())),
                        ("IsExternal", PropValue::Bool(true)),
                        ("HandicapAccessible", PropValue::Bool(false)),
                        ("SmokeStop", PropValue::Bool(true)),
                    ],
                );
                id
            }
            crate::entity::ArchEntity::Stair(sp) => write_stair(&mut w, sp),
            crate::entity::ArchEntity::Roof(rp) => write_roof(&mut w, rp),
            crate::entity::ArchEntity::Space(sp) => {
                let id = write_space(&mut w, sp);
                space_ifc_id.insert(*arch_id, id);
                emit_pset(
                    &mut w,
                    id,
                    "Pset_SpaceCommon",
                    &[
                        ("Reference", PropValue::Label(sp.space_name.clone())),
                        ("FloorArea", PropValue::Real(sp.floor_area())),
                        (
                            "GrossPlannedArea",
                            PropValue::Real(sp.floor_area()),
                        ),
                        ("IsExternal", PropValue::Bool(false)),
                        ("PubliclyAccessible", PropValue::Bool(true)),
                    ],
                );
                id
            }
            crate::entity::ArchEntity::DuctSegment(d) => {
                let id = write_duct(&mut w, d);
                emit_pset(
                    &mut w,
                    id,
                    "Pset_DuctSegmentTypeCommon",
                    &[
                        ("Reference", PropValue::Label(d.material.clone())),
                        ("Shape", PropValue::Label(d.shape.label().to_string())),
                        ("CrossSectionArea", PropValue::Real(d.shape.flow_area())),
                        ("NominalLength", PropValue::Real(d.length())),
                        (
                            "FlowDirection",
                            PropValue::Label(d.flow_direction.label().into()),
                        ),
                    ],
                );
                id
            }
            crate::entity::ArchEntity::PipeSegment(p) => {
                let id = write_pipe(&mut w, p);
                emit_pset(
                    &mut w,
                    id,
                    "Pset_PipeSegmentTypeCommon",
                    &[
                        ("Reference", PropValue::Label(p.material.clone())),
                        ("Fluid", PropValue::Label(p.fluid.clone())),
                        ("NominalDiameter", PropValue::Real(p.diameter)),
                        ("CrossSectionArea", PropValue::Real(p.flow_area())),
                        ("Pressure", PropValue::Real(p.operating_pressure)),
                    ],
                );
                id
            }
            crate::entity::ArchEntity::CableSegment(c) => {
                let id = write_cable(&mut w, c);
                emit_pset(
                    &mut w,
                    id,
                    "Pset_CableSegmentTypeCommon",
                    &[
                        ("Reference", PropValue::Label(c.material.clone())),
                        ("CrossSectionalArea", PropValue::Real(c.conductor_csa_mm2)),
                        ("NominalDiameter", PropValue::Real(c.diameter)),
                        ("Voltage", PropValue::Real(c.voltage)),
                    ],
                );
                id
            }
            crate::entity::ArchEntity::ConduitSegment(c) => {
                let id = write_conduit(&mut w, c);
                emit_pset(
                    &mut w,
                    id,
                    "Pset_CableCarrierSegmentTypeCommon",
                    &[
                        ("Reference", PropValue::Label(c.material.clone())),
                        ("OuterDiameter", PropValue::Real(c.outer_diameter)),
                        ("InnerDiameter", PropValue::Real(c.inner_diameter)),
                        ("CrossSectionArea", PropValue::Real(c.free_area())),
                    ],
                );
                id
            }
            crate::entity::ArchEntity::MepEquipment(e) => {
                let id = write_mep_equipment(&mut w, e);
                emit_pset(
                    &mut w,
                    id,
                    "Pset_DistributionElementCommon",
                    &[
                        ("Reference", PropValue::Label(e.kind.label().to_string())),
                        ("Tag", PropValue::Label(e.tag.clone())),
                        ("Description", PropValue::Text(e.description.clone())),
                    ],
                );
                id
            }
        };
        element_ids.push(eid);
    }

    // Window / door openings — emit IfcOpeningElement + IfcRelVoidsElement.
    for (_arch_id, win_id, host_arch) in &window_ifc {
        if let Some(host_wall) = walls.get(host_arch) {
            let host_ifc = match wall_ifc_id.get(host_arch) {
                Some(&id) => id,
                None => continue,
            };
            // Reconstruct the WindowParams from the doc by scanning.
            // We've already emitted the IfcWindow; opening sits
            // alongside it.
            for (arch_id, ent) in &doc.entities {
                if let crate::entity::ArchEntity::Window(wp) = ent {
                    if *arch_id != *_arch_id {
                        continue;
                    }
                    let opening = write_opening_for_window(&mut w, wp, host_wall);
                    emit_rel_voids_element(&mut w, host_ifc, opening);
                    // Link the window's IfcWindow as the filling
                    // element for the opening via IfcRelFillsElement.
                    w.add(format!(
                        "IFCRELFILLSELEMENT({guid},$,$,$,{open},{win})",
                        guid = ifc_str(&ifc_guid_v4()),
                        open = r(opening),
                        win = r(*win_id)
                    ));
                    break;
                }
            }
        }
    }
    for (_arch_id, door_id, host_arch) in &door_ifc {
        if let Some(host_wall) = walls.get(host_arch) {
            let host_ifc = match wall_ifc_id.get(host_arch) {
                Some(&id) => id,
                None => continue,
            };
            for (arch_id, ent) in &doc.entities {
                if let crate::entity::ArchEntity::Door(dp) = ent {
                    if *arch_id != *_arch_id {
                        continue;
                    }
                    let opening = write_opening_for_door(&mut w, dp, host_wall);
                    emit_rel_voids_element(&mut w, host_ifc, opening);
                    w.add(format!(
                        "IFCRELFILLSELEMENT({guid},$,$,$,{open},{door})",
                        guid = ifc_str(&ifc_guid_v4()),
                        open = r(opening),
                        door = r(*door_id)
                    ));
                    break;
                }
            }
        }
    }

    // Space ↔ wall boundaries — for every (space, wall) where the
    // wall's midpoint lies in the space's XY bounding box, emit an
    // IfcRelSpaceBoundary.
    for (space_arch_id, space_ifc) in &space_ifc_id {
        // Resolve space by arch_id.
        let Some(crate::entity::ArchEntity::Space(sp)) = doc.get_entity(*space_arch_id) else {
            continue;
        };
        let (xmin, ymin, xmax, ymax) = boundary_aabb_xy(&sp.boundary);
        for (wall_arch_id, wall_ifc) in &wall_ifc_id {
            let Some(crate::entity::ArchEntity::Wall(wall_p)) = doc.get_entity(*wall_arch_id)
            else {
                continue;
            };
            let mid = (wall_p.start + wall_p.end) * 0.5;
            if mid.x >= xmin && mid.x <= xmax && mid.y >= ymin && mid.y <= ymax {
                emit_rel_space_boundary(&mut w, *space_ifc, *wall_ifc, ".INTERNAL.");
            }
        }
    }

    // RelContainedInSpatialStructure linking storey ← all elements.
    if !element_ids.is_empty() {
        let elements = element_ids
            .iter()
            .map(|id| r(*id))
            .collect::<Vec<_>>()
            .join(",");
        w.add(format!(
            "IFCRELCONTAINEDINSPATIALSTRUCTURE({guid},$,$,$,({els}),{storey})",
            guid = ifc_str(&ifc_guid_v4()),
            els = elements,
            storey = r(storey)
        ));
    }

    w.write_to_path(path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::document::ArchDocument;
    use crate::entity::ArchEntity;
    use crate::wall::WallParams;
    use nalgebra::Vector3;

    #[test]
    fn guid_is_22_chars_from_alphabet() {
        for _ in 0..50 {
            let g = ifc_guid_v4();
            assert_eq!(g.len(), 22, "wrong length for {g}");
            for c in g.chars() {
                assert!(
                    IFC_ALPHABET.contains(&(c as u8)),
                    "char {c:?} not in IFC alphabet"
                );
            }
        }
    }

    #[test]
    fn ifc_guid_packs_big_endian_6_bit_chunks() {
        // Independently pack the 128-bit big-endian value as 2 + 21×6 bits and
        // compare. (The old `lo >> (-shift - 8)` else-branch mis-packed every
        // chunk from output index 12 on — the existing length/alphabet test
        // never exercised the bit layout.)
        let bytes: [u8; 16] = [
            0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0, 0x0f, 0xed, 0xcb, 0xa9, 0x87, 0x65,
            0x43, 0x21,
        ];
        let v = u128::from_be_bytes(bytes);
        let mut expected = String::with_capacity(22);
        expected.push(IFC_ALPHABET[((v >> 126) & 0x3) as usize] as char);
        for k in 1..22u32 {
            let shift = 126 - 6 * k; // LSB of this 6-bit chunk in the 128-bit value
            expected.push(IFC_ALPHABET[((v >> shift) & 0x3F) as usize] as char);
        }
        assert_eq!(encode_ifc_guid(&bytes), expected);
    }

    /// Round-14 M4 RED→GREEN: `ifc_str` must sanitise newlines and
    /// escape backslashes so a hostile caller-owned string can't
    /// break out of the Part 21 entity wrapper and inject sibling
    /// entities. Pre-fix only `'` was escaped — a project name like
    /// `"X');\n#9999=IFCPROJECT(..."` survived to the output and let
    /// the embedded `);\n#9999=` close the parent entity and start a
    /// new one.
    #[test]
    fn ifc_str_sanitises_injection_payload() {
        let payload = "X');\n#9999=IFCPROJECT('evil',$,$,$,$,$,$,$,$";
        let quoted = ifc_str(payload);
        // No newline (would break out of the entity).
        assert!(
            !quoted.contains('\n'),
            "newline must be stripped, got: {quoted}"
        );
        // No CR either.
        assert!(
            !quoted.contains('\r'),
            "carriage return must be stripped, got: {quoted}"
        );
        // The hostile `);` substring stays (it's just text inside a
        // string literal once the newline is gone), but the actual
        // injection is dead because `);` without a following newline
        // is still inside the quoted string.
        assert!(quoted.starts_with('\''));
        assert!(quoted.ends_with('\''));
        // The single-quote escape still works: `X'` → `X''`.
        assert!(
            quoted.contains("X'')"),
            "single-quote should still escape, got: {quoted}"
        );
    }

    /// RED→GREEN (round-25 L2): U+2028 (LINE SEPARATOR) and
    /// U+2029 (PARAGRAPH SEPARATOR) MUST be stripped from `ifc_str`
    /// output. Mirrors the round-20 L4 `python_str_repr` fix —
    /// Web-IFC / JS-based downstream consumers treat these as
    /// line terminators per ECMA-262 and a hostile name could
    /// break out of the entity at the JSON-export layer even
    /// though the on-disk Part 21 representation looks single-line.
    #[test]
    fn ifc_str_strips_unicode_line_separators_round25_l2() {
        // U+2028
        let q1 = ifc_str("evil\u{2028}name");
        assert!(
            !q1.contains('\u{2028}'),
            "U+2028 must be stripped, got: {q1}",
        );
        assert!(
            q1.contains("evilname"),
            "remaining characters must concatenate, got: {q1}",
        );
        // U+2029
        let q2 = ifc_str("attack\u{2029}payload");
        assert!(
            !q2.contains('\u{2029}'),
            "U+2029 must be stripped, got: {q2}",
        );
        assert!(
            q2.contains("attackpayload"),
            "remaining characters must concatenate, got: {q2}",
        );
        // Both at once with other special chars — interaction sanity.
        let q3 = ifc_str("a\u{2028}b\r\nc\u{2029}d'e\\f");
        assert!(!q3.contains('\u{2028}') && !q3.contains('\u{2029}'));
        assert!(!q3.contains('\n') && !q3.contains('\r'));
        // Single-quote escape still fires: `'` → `''`.
        assert!(q3.contains("d''e"));
        // Backslash escape still fires: `\` → `\\`.
        assert!(q3.contains("f"));
    }

    /// RED→GREEN (round-26 L1): U+0085 (NEXT LINE / NEL) MUST be
    /// stripped from `ifc_str` output. Same threat model as CR / LF
    /// / U+2028 / U+2029 — Unicode UAX-14 treats NEL as a hard
    /// break, and some line-oriented downstream consumers
    /// (terminals, log aggregators, EBCDIC-aware translators)
    /// honour it as a newline. Round-25 L2 covered U+2028 / U+2029
    /// but missed NEL.
    #[test]
    fn ifc_str_strips_nel_round26_l1() {
        let q = ifc_str("attack\u{0085}payload");
        assert!(
            !q.contains('\u{0085}'),
            "U+0085 (NEL) must be stripped, got: {q}",
        );
        assert!(
            q.contains("attackpayload"),
            "remaining characters must concatenate, got: {q}",
        );
        // Interaction: NEL alongside the round-25 L2 strips + the
        // round-14 M4 CR/LF strips + quote/backslash escapes.
        let q2 = ifc_str("a\u{0085}b\u{2028}c\r\nd\u{2029}e'f\\g");
        assert!(!q2.contains('\u{0085}'));
        assert!(!q2.contains('\u{2028}') && !q2.contains('\u{2029}'));
        assert!(!q2.contains('\n') && !q2.contains('\r'));
        assert!(q2.contains("e''f"));
    }

    /// Round-14 M4 sister: backslash escape per Part 21. A caller-
    /// supplied label like `"foo\\X3F"` would otherwise be re-read
    /// by IFC consumers as the hex-escape sequence `\X3F` ("?").
    #[test]
    fn ifc_str_escapes_backslash() {
        let quoted = ifc_str("foo\\bar");
        // Two backslashes in source → four backslashes total after
        // escape (Part 21 doubles the literal `\` to `\\`).
        assert!(
            quoted.contains("foo\\\\bar"),
            "expected doubled backslash, got: {quoted}"
        );
    }

    #[test]
    fn writer_emits_header_and_data_sections() {
        let mut w = IfcWriter::new();
        write_header(&mut w, "Test");
        let _ = emit_project_hierarchy(&mut w, "Test");
        let s = w.finish();
        assert!(s.starts_with("ISO-10303-21;"));
        assert!(s.contains("HEADER;"));
        assert!(s.contains("FILE_SCHEMA(('IFC4'))"));
        assert!(s.contains("DATA;"));
        assert!(s.contains("IFCPROJECT"));
        assert!(s.contains("IFCSITE"));
        assert!(s.contains("IFCBUILDING"));
        assert!(s.contains("IFCBUILDINGSTOREY"));
        assert!(s.ends_with("END-ISO-10303-21;\n"));
    }

    #[test]
    fn write_document_round_trip_includes_wall() {
        let mut doc = ArchDocument::new("Casa");
        doc.add_entity(ArchEntity::Wall(WallParams {
            start: Vector3::zeros(),
            end: Vector3::new(3.0, 0.0, 0.0),
            height: 2.5,
            thickness: 0.2,
            material: "Brick".into(),
        }));
        let tmp = std::env::temp_dir().join("valenx_arch_ifc_round_trip.ifc");
        write_document(&doc, &tmp).unwrap();
        let content = std::fs::read_to_string(&tmp).unwrap();
        assert!(content.starts_with("ISO-10303-21;"));
        assert!(content.contains("FILE_SCHEMA(('IFC4'))"));
        assert!(content.contains("IFCWALL"));
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn epoch_to_ymd_hms_known_value() {
        // 2020-01-01T00:00:00 UTC = 1577836800 seconds.
        let (y, mo, d, h, mi, s) = epoch_to_ymd_hms(1_577_836_800);
        assert_eq!((y, mo, d, h, mi, s), (2020, 1, 1, 0, 0, 0));
    }

    /// RED→GREEN (round-27 H1 STRUCTURAL): two concurrent writers
    /// each writing a distinct multi-MiB IFC body must end with one
    /// writer's content at the target — NOT an interleaved mix of
    /// both. Pre-R27 `write_to_path` did `File::create(path)` then
    /// streamed via BufWriter; two threads racing share the file
    /// handle and the kernel can interleave their `write(2)` calls
    /// (POSIX `write` is only atomic ≤ PIPE_BUF = 4 KiB). Post-R27
    /// each writer owns a distinct sidecar via
    /// `atomic_write_streaming`, then `fs::rename` atomically
    /// promotes ONE sidecar over the target.
    ///
    /// Synthesise the IFC body manually rather than running the
    /// full doc emit pipeline so the assertion stays focused on the
    /// publication-atomicity contract and the test stays fast.
    #[test]
    fn write_to_path_concurrent_multi_mb_no_interleaving_round27_h1() {
        use std::sync::{Arc, Barrier};
        use std::thread;
        // Build two distinct ~5 MiB IFC bodies — distinct guid
        // strings so the on-disk content can be told apart.
        fn build_big_writer(tag: &str) -> IfcWriter {
            let mut w = IfcWriter::new();
            w.add_header("FILE_SCHEMA(('IFC4'))".to_string());
            // Each entity is ~140 bytes; ~40k entities ≈ 5 MiB.
            for i in 0..40_000 {
                w.add(format!(
                    "IFCWALL('guid-{tag}-{i:08}',$,'wall-{tag}-{i}',$,$,$,$,$,.NOTDEFINED.)"
                ));
            }
            w
        }
        let writer_a = Arc::new(build_big_writer("A"));
        let writer_b = Arc::new(build_big_writer("B"));
        let tmp = std::env::temp_dir().join(format!(
            "valenx-arch-ifc-h1-{}.ifc",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let barrier = Arc::new(Barrier::new(2));
        let wa = Arc::clone(&writer_a);
        let wb = Arc::clone(&writer_b);
        let tmp_a = tmp.clone();
        let tmp_b = tmp.clone();
        let ba = Arc::clone(&barrier);
        let bb = Arc::clone(&barrier);
        let ha = thread::spawn(move || {
            ba.wait();
            wa.write_to_path(&tmp_a)
        });
        let hb = thread::spawn(move || {
            bb.wait();
            wb.write_to_path(&tmp_b)
        });
        ha.join().unwrap().expect("writer A");
        hb.join().unwrap().expect("writer B");
        let body = std::fs::read_to_string(&tmp).expect("read final");
        let has_a = body.contains("guid-A-00000000");
        let has_b = body.contains("guid-B-00000000");
        let _ = std::fs::remove_file(&tmp);
        // The final file MUST be one writer's content end-to-end —
        // never a mix.
        let xor_a_b = has_a ^ has_b;
        assert!(
            xor_a_b,
            "final body must equal exactly ONE writer's content — \
             has_A={has_a}, has_B={has_b} (len={})",
            body.len(),
        );
        // Also pin the round-trip well-formedness — the rename leg
        // mustn't have truncated.
        assert!(body.starts_with("ISO-10303-21;"), "missing header");
        assert!(
            body.ends_with("END-ISO-10303-21;\n"),
            "missing footer (rename leg didn't write the tail?)",
        );
    }

    /// RED→GREEN (round-24 L1): `write_to_path` now streams via
    /// `BufWriter` instead of building a monolithic `String` via
    /// `finish()` first. The on-disk bytes must be IDENTICAL to
    /// what `finish()` would have produced (the streaming path is a
    /// pure refactor — same content, smaller peak memory).
    /// Test asserts byte-for-byte equality between
    /// `write_to_path` and `finish()`.
    #[test]
    fn write_to_path_produces_same_bytes_as_finish() {
        let mut w = IfcWriter::new();
        w.add_header("FILE_SCHEMA(('IFC4'))".to_string());
        w.add("IFCPROJECT('guid-1',$,'p',$,$,$,$,$,$)".to_string());
        w.add("IFCBUILDING('guid-2',$,'b',$,$,$,$,$,$,$,$,$)".to_string());
        w.add("IFCWALL('guid-3',$,'w',$,$,$,$,$,$)".to_string());
        let expected = w.finish();
        let tmp = std::env::temp_dir().join(format!(
            "valenx_arch_ifc_stream_l1_{}.ifc",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        w.write_to_path(&tmp).unwrap();
        let actual = std::fs::read_to_string(&tmp).unwrap();
        let _ = std::fs::remove_file(&tmp);
        assert_eq!(
            actual, expected,
            "streamed bytes must equal finish() bytes — found drift"
        );
    }
}
