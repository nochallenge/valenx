//! ParaView Data (`.pvd`) time-series wrapper.
//!
//! `.pvd` is the XML manifest ParaView writes for multi-timestep
//! datasets. Each entry points at one `.vtu` (or `.vtp`, `.vtk`,
//! etc.) snapshot plus the time value to associate with it. Loaders
//! walk the entries in order and present a time-keyed view of the
//! field catalog.
//!
//! Layout:
//!
//! ```xml
//! <?xml version="1.0"?>
//! <VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
//!   <Collection>
//!     <DataSet timestep="0.0" group="" part="0" file="step_0.vtu"/>
//!     <DataSet timestep="0.5" group="" part="0" file="step_1.vtu"/>
//!     <DataSet timestep="1.0" group="" part="0" file="step_2.vtu"/>
//!   </Collection>
//! </VTKFile>
//! ```
//!
//! This module owns the writer + reader for the manifest only;
//! resolving each referenced file's bytes is the caller's job
//! (typically through [`crate::vtk_dispatch::load_canonical`]).

use std::io;
use std::path::{Path, PathBuf};

use thiserror::Error;

/// One time-step entry in a PVD collection.
#[derive(Clone, Debug, PartialEq)]
pub struct PvdEntry {
    /// Time value associated with this snapshot. Units are case-
    /// dependent (seconds for transient CFD, dimensionless iteration
    /// counters for steady cases that emit pseudo-time).
    pub timestep: f64,
    /// Path to the snapshot file. Stored verbatim from the
    /// manifest's `file=` attribute; can be relative (resolved
    /// against the PVD's parent dir) or absolute.
    pub file: PathBuf,
    /// Optional `group=` attribute. ParaView uses this to bin
    /// snapshots into named groups; defaults to empty when omitted.
    pub group: String,
    /// Optional `part=` attribute for parallel decompositions. 0
    /// for the common single-rank case.
    pub part: u32,
}

/// Top-level PVD collection — a list of entries plus the byte order
/// declared on `VTKFile`.
#[derive(Clone, Debug, PartialEq)]
pub struct PvdCollection {
    pub byte_order: String,
    pub entries: Vec<PvdEntry>,
}

impl PvdCollection {
    /// Convenience: a fresh collection with `LittleEndian` byte
    /// order (matches every modern x86_64 / ARM64 host).
    pub fn new() -> Self {
        Self {
            byte_order: "LittleEndian".into(),
            entries: Vec::new(),
        }
    }

    /// Push an entry. Returns the new entry count for fluent use.
    pub fn add(&mut self, entry: PvdEntry) -> usize {
        self.entries.push(entry);
        self.entries.len()
    }

    /// Sort entries by `timestep` ascending. ParaView expects this
    /// ordering for its time-slider; out-of-order entries trigger
    /// "ParaView reordered the dataset" warnings.
    pub fn sort_by_time(&mut self) {
        self.entries.sort_by(|a, b| {
            a.timestep
                .partial_cmp(&b.timestep)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Render the collection as a PVD manifest string.
    pub fn to_xml(&self) -> String {
        let mut s = String::with_capacity(256 + self.entries.len() * 80);
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str(&format!(
            "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"{}\">\n",
            self.byte_order
        ));
        s.push_str("  <Collection>\n");
        for e in &self.entries {
            // Forward-slash file paths so the manifest works on
            // Linux + macOS hosts even when written from Windows.
            let file = e.file.to_string_lossy().replace('\\', "/");
            s.push_str(&format!(
                "    <DataSet timestep=\"{ts}\" group=\"{group}\" part=\"{part}\" file=\"{file}\"/>\n",
                ts = e.timestep,
                group = xml_escape_attr(&e.group),
                part = e.part,
                file = xml_escape_attr(&file),
            ));
        }
        s.push_str("  </Collection>\n");
        s.push_str("</VTKFile>\n");
        s
    }
}

impl Default for PvdCollection {
    fn default() -> Self {
        Self::new()
    }
}

/// Write a PVD manifest to disk. Creates parent directories as
/// needed.
pub fn write_pvd(collection: &PvdCollection, path: &Path) -> io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    atomic_write_local(path, collection.to_xml().as_bytes())
}

/// Crash-safe write: serialise to a process-unique sidecar, fsync it,
/// then atomically `rename` it over the destination. A reader (e.g.
/// ParaView) therefore never observes a half-written manifest, and an
/// interrupted write can't leave a truncated `.pvd` at the real path.
///
/// This is a *local* re-implementation of
/// `valenx_core::io_caps::atomic_write_bytes`. valenx-fields cannot
/// depend on valenx-core: the workspace dependency chain runs
/// `valenx-core → valenx-fields → valenx-mesh`, so valenx-core is
/// *upstream* of this crate and adding a `valenx-core` dependency here
/// would close a cycle Cargo rejects outright (the same constraint
/// documented on `MAX_OBJ_LINE_BYTES` / `MAX_PLY_FILE_BYTES` in
/// valenx-mesh). Keep this small helper in sync with the canonical one.
fn atomic_write_local(path: &Path, contents: &[u8]) -> io::Result<()> {
    use std::io::Write as _;

    let file_name = path.file_name().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "atomic_write_local: target path has no file name component",
        )
    })?;
    // Process-unique sidecar name so two concurrent writers to sibling
    // paths in the same directory cannot collide on the temp file.
    let mut tmp_name = file_name.to_os_string();
    tmp_name.push(format!(".tmp.{}", std::process::id()));
    let tmp_path: PathBuf = match path.parent() {
        Some(parent) => parent.join(&tmp_name),
        None => PathBuf::from(&tmp_name),
    };

    // Scope the file handle so it is closed before the rename.
    {
        let mut f = std::fs::File::create(&tmp_path)?;
        if let Err(e) = f.write_all(contents).and_then(|()| f.sync_all()) {
            // Best-effort cleanup of the orphan sidecar on any write/sync
            // failure so we don't litter the directory.
            let _ = std::fs::remove_file(&tmp_path);
            return Err(e);
        }
    }
    if let Err(e) = std::fs::rename(&tmp_path, path) {
        let _ = std::fs::remove_file(&tmp_path);
        return Err(e);
    }
    Ok(())
}

/// Errors raised by [`parse_pvd`] and [`resolve_entry_path`].
#[derive(Debug, Error)]
pub enum PvdError {
    /// The XML did not start with `<VTKFile type="Collection">`.
    #[error("missing required `<VTKFile type=\"Collection\">` root tag")]
    MissingRoot,
    /// No `<Collection>` child element under the root.
    #[error("missing required `<Collection>` element")]
    MissingCollection,
    /// A `<DataSet>` element was malformed (missing required attribute,
    /// unparseable timestamp, etc.).
    #[error("malformed DataSet entry at position {index}: {reason}")]
    BadEntry {
        /// 0-based index of the offending `<DataSet>` in source order.
        index: usize,
        /// Short human-readable explanation.
        reason: String,
    },
    /// Round-19 M2: a `<DataSet file="...">` entry resolved to a path
    /// that escapes the PVD manifest's parent directory. A hostile
    /// manifest could otherwise read arbitrary host files via
    /// `<DataSet file="../../etc/passwd"/>`; resolved paths must stay
    /// rooted under the PVD's containing directory.
    #[error("DataSet file `{file}` escapes the PVD manifest's parent directory")]
    PathEscape {
        /// The offending `file=` value verbatim from the manifest.
        file: String,
    },
}

/// Parse a PVD manifest string. Returns the entries in declaration
/// order — call [`PvdCollection::sort_by_time`] if the source isn't
/// already time-ordered.
pub fn parse_pvd(text: &str) -> Result<PvdCollection, PvdError> {
    if !text.contains("<VTKFile") || !text.contains("type=\"Collection\"") {
        return Err(PvdError::MissingRoot);
    }
    let byte_order = parse_attr(text, "byte_order").unwrap_or_else(|| "LittleEndian".into());
    let coll_start = text
        .find("<Collection>")
        .ok_or(PvdError::MissingCollection)?;
    let coll_end = text[coll_start..]
        .find("</Collection>")
        .ok_or(PvdError::MissingCollection)?
        + coll_start;
    let body = &text[coll_start + "<Collection>".len()..coll_end];

    let mut entries: Vec<PvdEntry> = Vec::new();
    let mut cursor = 0usize;
    let mut idx = 0usize;
    loop {
        let rest = &body[cursor..];
        let Some(open) = rest.find("<DataSet") else {
            break;
        };
        let after_open = open + "<DataSet".len();
        let Some(close) = rest[after_open..].find('>') else {
            return Err(PvdError::BadEntry {
                index: idx,
                reason: "missing `>` after DataSet open tag".into(),
            });
        };
        let attrs_end = after_open + close;
        let attrs_str = &rest[after_open..attrs_end];
        // Self-closed `<DataSet ... />` (the conventional form) or
        // body-bearing — both work. Skip past either to advance.
        cursor += attrs_end + 1;

        let timestep_str = parse_attr(attrs_str, "timestep").ok_or_else(|| PvdError::BadEntry {
            index: idx,
            reason: "missing `timestep=` attribute".into(),
        })?;
        let timestep = timestep_str
            .parse::<f64>()
            .map_err(|e| PvdError::BadEntry {
                index: idx,
                reason: format!("timestep `{timestep_str}` not a number: {e}"),
            })?;
        let file = parse_attr(attrs_str, "file").ok_or_else(|| PvdError::BadEntry {
            index: idx,
            reason: "missing `file=` attribute".into(),
        })?;
        let group = parse_attr(attrs_str, "group").unwrap_or_default();
        let part = parse_attr(attrs_str, "part")
            .and_then(|s| s.parse::<u32>().ok())
            .unwrap_or(0);

        entries.push(PvdEntry {
            timestep,
            file: PathBuf::from(file),
            group,
            part,
        });
        idx += 1;
    }
    Ok(PvdCollection {
        byte_order,
        entries,
    })
}

/// Resolve a [`PvdEntry::file`] against the PVD manifest's parent
/// directory. Mirrors how `valenx-plugin::loader::resolve_payload`
/// handles relative manifests — relative paths become
/// `<pvd_dir>/<file>`, absolute paths pass through unchanged.
///
/// Round-19 M2 hardening: relative entries are canonicalised and
/// prefix-checked against the canonicalised PVD parent directory.
/// A hostile manifest containing `<DataSet file="../../etc/passwd"/>`
/// (or any other `..`-bearing path that resolves above the PVD's
/// parent) is rejected with [`PvdError::PathEscape`]. Without the
/// check, a poisoned PVD shared between operators could exfiltrate
/// arbitrary host files at load time.
///
/// Absolute paths pass through unchanged — they're the caller's
/// problem to validate; the legitimate use case is a PVD that
/// references a sibling directory's pre-staged snapshots via an
/// explicit absolute path.
///
/// When the canonicalisation fails (file doesn't exist yet — fine
/// during write-then-read flows), we fall back to a lexical
/// `..`-component scan so the check still catches the obvious
/// traversal shape before any IO happens.
pub fn resolve_entry_path(pvd_path: &Path, entry: &PvdEntry) -> Result<PathBuf, PvdError> {
    if entry.file.is_absolute() {
        return Ok(entry.file.clone());
    }
    let parent = pvd_path.parent().unwrap_or(Path::new("."));
    let joined = parent.join(&entry.file);
    // Best-effort canonicalisation. If both sides canonicalise, the
    // prefix-check is authoritative. If the PVD parent doesn't
    // canonicalise (very rare — the manifest itself was read from
    // disk, so its parent exists) we still run the lexical
    // `..`-component scan below as defence-in-depth.
    let parent_canon = std::fs::canonicalize(parent);
    let joined_canon = std::fs::canonicalize(&joined);
    if let (Ok(pc), Ok(jc)) = (parent_canon, joined_canon) {
        if !jc.starts_with(&pc) {
            return Err(PvdError::PathEscape {
                file: entry.file.to_string_lossy().into_owned(),
            });
        }
        return Ok(jc);
    }
    // Fallback: lexical scan for `..` traversal. A `..` component
    // anywhere in the relative path is the canonical escape vector;
    // refuse it without depending on the filesystem state.
    if entry
        .file
        .components()
        .any(|c| matches!(c, std::path::Component::ParentDir))
    {
        return Err(PvdError::PathEscape {
            file: entry.file.to_string_lossy().into_owned(),
        });
    }
    Ok(joined)
}

fn parse_attr(attrs: &str, name: &str) -> Option<String> {
    // Look for `name="value"` substring — sufficient for the
    // tightly-constrained attribute shapes a PVD file uses.
    let needle = format!("{name}=\"");
    let start = attrs.find(&needle)? + needle.len();
    let end = attrs[start..].find('"')? + start;
    Some(attrs[start..end].to_string())
}

fn xml_escape_attr(s: &str) -> String {
    // Minimal escape: `&`, `<`, `"`. PVD attributes don't use `>`
    // or `'` in practice (filenames + integer IDs).
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('"', "&quot;")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_collection() -> PvdCollection {
        let mut coll = PvdCollection::new();
        coll.add(PvdEntry {
            timestep: 0.0,
            file: PathBuf::from("step_0.vtu"),
            group: String::new(),
            part: 0,
        });
        coll.add(PvdEntry {
            timestep: 0.5,
            file: PathBuf::from("step_1.vtu"),
            group: String::new(),
            part: 0,
        });
        coll.add(PvdEntry {
            timestep: 1.0,
            file: PathBuf::from("step_2.vtu"),
            group: String::new(),
            part: 0,
        });
        coll
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

    /// Robustness: `parse_pvd` must NEVER panic on adversarial input — only
    /// return `Ok` or a typed `PvdError`. Feeds truncations, single-byte
    /// corruptions, random buffers, and burst-corruptions of a valid PVD; a
    /// panic in any would fail the test via propagation. parse_pvd terminates
    /// for all inputs (the DataSet cursor advances >= 9 bytes per match) and
    /// has no count-driven allocation, so every iteration is small and fast.
    #[test]
    fn parse_pvd_never_panics_on_adversarial_input() {
        let valid = sample_collection().to_xml();
        let bytes = valid.as_bytes();

        // 1. Every truncated prefix of a valid PVD manifest.
        for k in 0..=bytes.len() {
            let _ = parse_pvd(&String::from_utf8_lossy(&bytes[..k]));
        }
        // 2. Single-byte corruption at every position.
        for i in 0..bytes.len() {
            let mut m = bytes.to_vec();
            m[i] ^= 0xFF;
            let _ = parse_pvd(&String::from_utf8_lossy(&m));
        }
        // 3. Deterministic pseudo-random buffers (mostly the MissingRoot guard).
        let mut state: u64 = 0x6A09_E667_F3BC_C908; // fixed seed
        for _ in 0..2000 {
            let len = (xorshift64(&mut state) % 256) as usize;
            let buf: Vec<u8> = (0..len).map(|_| (xorshift64(&mut state) & 0xFF) as u8).collect();
            let _ = parse_pvd(&String::from_utf8_lossy(&buf));
        }
        // 4. Burst corruption of the valid manifest: the <VTKFile…Collection>
        //    structure survives a few byte changes, so the parser proceeds
        //    into the DataSet attribute loop (timestep/file/part parsing) with
        //    hostile values — the deep-coverage phase.
        for _ in 0..2000 {
            let mut m = bytes.to_vec();
            let flips = 1 + xorshift64(&mut state) % 6;
            for _ in 0..flips {
                let pos = (xorshift64(&mut state) as usize) % m.len();
                m[pos] = (xorshift64(&mut state) & 0xFF) as u8;
            }
            let _ = parse_pvd(&String::from_utf8_lossy(&m));
        }
    }

    #[test]
    fn to_xml_emits_well_formed_pvd_root() {
        let xml = sample_collection().to_xml();
        assert!(xml.starts_with("<?xml version=\"1.0\"?>"), "got:\n{xml}");
        assert!(xml.contains("<VTKFile type=\"Collection\""));
        assert!(xml.contains("<Collection>"));
        assert!(xml.contains("</Collection>"));
        assert!(xml.contains("</VTKFile>"));
    }

    #[test]
    fn to_xml_emits_one_dataset_per_entry() {
        let xml = sample_collection().to_xml();
        let count = xml.matches("<DataSet ").count();
        assert_eq!(count, 3);
    }

    #[test]
    fn to_xml_normalises_backslashes_to_forward_slash() {
        // PVD manifests are cross-platform — Windows backslashes must
        // be translated to forward slashes so the same file works on
        // Linux + macOS hosts.
        let mut coll = PvdCollection::new();
        coll.add(PvdEntry {
            timestep: 0.0,
            file: PathBuf::from("snapshots\\step_0.vtu"),
            group: String::new(),
            part: 0,
        });
        let xml = coll.to_xml();
        assert!(xml.contains("file=\"snapshots/step_0.vtu\""), "got:\n{xml}");
        assert!(!xml.contains("\\"), "backslash leaked: {xml}");
    }

    #[test]
    fn round_trip_preserves_entries() {
        let original = sample_collection();
        let xml = original.to_xml();
        let parsed = parse_pvd(&xml).expect("parse");
        assert_eq!(parsed, original);
    }

    #[test]
    fn parse_rejects_input_without_collection_root() {
        let text = "<?xml version=\"1.0\"?>\n<NotVTK/>\n";
        let err = parse_pvd(text).unwrap_err();
        assert!(matches!(err, PvdError::MissingRoot));
    }

    #[test]
    fn parse_rejects_dataset_missing_timestep() {
        let text = r#"<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>
    <DataSet group="" part="0" file="step_0.vtu"/>
  </Collection>
</VTKFile>"#;
        let err = parse_pvd(text).unwrap_err();
        match err {
            PvdError::BadEntry { index, reason } => {
                assert_eq!(index, 0);
                assert!(reason.contains("timestep"));
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn parse_handles_optional_group_and_part() {
        // Entries that omit `group=` and `part=` should default to
        // empty / 0 rather than failing.
        let text = r#"<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>
    <DataSet timestep="0.0" file="step_0.vtu"/>
  </Collection>
</VTKFile>"#;
        let coll = parse_pvd(text).expect("parse");
        assert_eq!(coll.entries.len(), 1);
        assert_eq!(coll.entries[0].group, "");
        assert_eq!(coll.entries[0].part, 0);
    }

    #[test]
    fn sort_by_time_orders_entries_ascending() {
        let mut coll = PvdCollection::new();
        coll.add(PvdEntry {
            timestep: 1.0,
            file: PathBuf::from("a.vtu"),
            group: String::new(),
            part: 0,
        });
        coll.add(PvdEntry {
            timestep: 0.5,
            file: PathBuf::from("b.vtu"),
            group: String::new(),
            part: 0,
        });
        coll.add(PvdEntry {
            timestep: 0.0,
            file: PathBuf::from("c.vtu"),
            group: String::new(),
            part: 0,
        });
        coll.sort_by_time();
        assert_eq!(coll.entries[0].file, PathBuf::from("c.vtu"));
        assert_eq!(coll.entries[1].file, PathBuf::from("b.vtu"));
        assert_eq!(coll.entries[2].file, PathBuf::from("a.vtu"));
    }

    #[test]
    fn write_pvd_round_trips_through_disk() {
        let path = std::env::temp_dir().join(format!(
            "valenx-pvd-roundtrip-{}.pvd",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let coll = sample_collection();
        write_pvd(&coll, &path).expect("write");
        let text = std::fs::read_to_string(&path).expect("read");
        let parsed = parse_pvd(&text).expect("parse");
        assert_eq!(parsed, coll);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn resolve_entry_path_handles_relative_and_absolute() {
        // Use a non-existent parent so canonicalize fails on both
        // sides; the fallback lexical-scan path then kicks in and
        // returns the plain join (no `..` to refuse).
        let pvd_path = PathBuf::from("/runs/case-A/transient.pvd");
        let rel = PvdEntry {
            timestep: 0.0,
            file: PathBuf::from("step_0.vtu"),
            group: String::new(),
            part: 0,
        };
        assert_eq!(
            resolve_entry_path(&pvd_path, &rel).expect("resolve"),
            PathBuf::from("/runs/case-A/step_0.vtu"),
        );
        let abs_path = if cfg!(windows) {
            PathBuf::from(r"C:\absolute\step_0.vtu")
        } else {
            PathBuf::from("/absolute/step_0.vtu")
        };
        let abs = PvdEntry {
            timestep: 0.0,
            file: abs_path.clone(),
            group: String::new(),
            part: 0,
        };
        assert_eq!(resolve_entry_path(&pvd_path, &abs).expect("resolve"), abs_path);
    }

    /// Round-19 M2 RED→GREEN: a PVD whose `<DataSet>` entry contains
    /// `file="../../etc/passwd"` (or any other `..` traversal that
    /// resolves above the manifest's parent directory) must be
    /// rejected at resolve time. Pre-fix `resolve_entry_path`
    /// returned the joined `pvd_dir/../../etc/passwd` verbatim and
    /// any subsequent `fs::read_to_string` would happily follow the
    /// traversal out of the case bundle.
    #[test]
    fn resolve_entry_path_rejects_dotdot_traversal() {
        // Use a real on-disk parent so the canonicalize-then-prefix
        // branch is exercised (the lexical fallback is also covered
        // by the second test below).
        let parent_dir = std::env::temp_dir().join(format!(
            "valenx-pvd-traversal-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&parent_dir).unwrap();
        let pvd_path = parent_dir.join("transient.pvd");
        std::fs::write(&pvd_path, b"<placeholder/>").unwrap();
        let traversal = PvdEntry {
            timestep: 0.0,
            file: PathBuf::from("../../etc/passwd"),
            group: String::new(),
            part: 0,
        };
        let err = resolve_entry_path(&pvd_path, &traversal)
            .expect_err("traversal must be refused");
        match err {
            PvdError::PathEscape { file } => {
                assert!(file.contains(".."), "got: {file:?}");
            }
            other => panic!("wrong error variant: {other:?}"),
        }
        let _ = std::fs::remove_dir_all(&parent_dir);
    }

    /// Defence-in-depth: even when the PVD parent doesn't exist on
    /// disk (so canonicalize fails), the lexical `..` scan must
    /// still catch the obvious traversal shape. Hits the fallback
    /// branch deterministically by pointing at a synthetic
    /// /nonexistent prefix.
    #[test]
    fn resolve_entry_path_lexical_fallback_rejects_dotdot() {
        let pvd_path = PathBuf::from("/nonexistent-pvd-fallback/transient.pvd");
        let traversal = PvdEntry {
            timestep: 0.0,
            file: PathBuf::from("../../../etc/passwd"),
            group: String::new(),
            part: 0,
        };
        let err = resolve_entry_path(&pvd_path, &traversal)
            .expect_err("lexical fallback must refuse `..`");
        assert!(matches!(err, PvdError::PathEscape { .. }), "got: {err:?}");
    }

    #[test]
    fn xml_escape_attr_handles_special_chars() {
        assert_eq!(xml_escape_attr("plain"), "plain");
        assert_eq!(xml_escape_attr("a&b"), "a&amp;b");
        assert_eq!(xml_escape_attr("a<b"), "a&lt;b");
        assert_eq!(xml_escape_attr("a\"b"), "a&quot;b");
    }
}
