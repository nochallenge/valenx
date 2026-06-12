//! Mesh discovery + loading from a solver workdir. The auto-loader
//! routes through here after a run finishes — pick the most-recent
//! VTU / VTK / PVD snapshot, parse it via `valenx_fields::vtk_dispatch`,
//! and hand the canonical [`Mesh`] back to the viewport host. Also
//! the mesh-bounding-box utility the viewport uses to frame a fresh
//! load.

use std::path::PathBuf;

use valenx_mesh::Mesh;

/// Axis-aligned bounding box over a canonical mesh's node
/// coordinates. Returns `None` for an empty mesh so the viewport
/// can fall back to its default camera.
pub(crate) fn mesh_bounding_box(mesh: &Mesh) -> Option<([f32; 3], [f32; 3])> {
    let first = mesh.nodes.first()?;
    let mut min = [first.x as f32, first.y as f32, first.z as f32];
    let mut max = min;
    for n in &mesh.nodes {
        let v = [n.x as f32, n.y as f32, n.z as f32];
        for i in 0..3 {
            if v[i] < min[i] {
                min[i] = v[i];
            }
            if v[i] > max[i] {
                max[i] = v[i];
            }
        }
    }
    Some((min, max))
}

/// Pick the most-recent snapshot the auto-loader should drop into the
/// viewport. Prefers a `.pvd` time-series manifest when one is present
/// (the manifest is the solver's curated record of which snapshot is
/// "latest"), falling back to a lexicographic walk of `.vtu`/`.vtk`
/// files via [`latest_vtk_in_workdir`].
///
/// PVD selection: the entry with the highest `timestep` wins; on ties
/// the last-declared entry wins (stable behaviour for transient
/// solvers that write the manifest in-order).
pub(crate) fn latest_snapshot_in_workdir(workdir: &std::path::Path) -> Option<PathBuf> {
    if let Some(pvd) = first_pvd_in_workdir(workdir) {
        if let Some(snap) = latest_entry_from_pvd(&pvd) {
            return Some(snap);
        }
        // Manifest present but unparseable — log and fall through to
        // the raw .vtu/.vtk walk so the user still sees *something*
        // instead of a blank viewport.
        tracing::warn!(
            target: "valenx",
            ?pvd,
            "PVD present but couldn't extract latest entry; falling back to .vtu/.vtk walk"
        );
    }
    latest_vtk_in_workdir(workdir)
}

/// First `.pvd` file found under the workdir (lexicographic). Most
/// solvers emit one manifest per case; the few that emit several
/// (e.g. one per output group) are not yet a priority — picking the
/// first is a useful behaviour even for those because it gives the
/// user some snapshot rather than none.
pub(crate) fn first_pvd_in_workdir(workdir: &std::path::Path) -> Option<PathBuf> {
    let mut found: Vec<PathBuf> = Vec::new();
    let mut stack: Vec<PathBuf> = vec![workdir.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let Ok(entries) = std::fs::read_dir(&dir) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if let Ok(ft) = entry.file_type() {
                if ft.is_dir() {
                    stack.push(path);
                    continue;
                }
            }
            let is_pvd = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case("pvd"))
                .unwrap_or(false);
            if is_pvd {
                found.push(path);
            }
        }
    }
    found.sort();
    found.into_iter().next()
}

/// Parse a `.pvd` manifest and return the resolved path of the entry
/// with the highest `timestep`. `None` when the file is unreadable,
/// the XML is malformed, the manifest exceeds the size cap, or it is
/// empty.
///
/// Round-18 M2: the read is capped at
/// [`valenx_core::io_caps::MAX_PVD_FILE_BYTES`] so a hostile
/// multi-GB manifest can't slurp into memory before the XML parser
/// sees any of it.
pub(crate) fn latest_entry_from_pvd(pvd_path: &std::path::Path) -> Option<PathBuf> {
    let text = valenx_core::io_caps::read_capped_to_string(
        pvd_path,
        valenx_core::io_caps::MAX_PVD_FILE_BYTES,
    )
    .ok()?;
    let coll = valenx_fields::pvd::parse_pvd(&text).ok()?;
    if coll.entries.is_empty() {
        return None;
    }
    // Find the max timestep, breaking ties on later index. partial_cmp
    // returns None for NaN — treat those as "smaller than anything"
    // so a single bad timestep doesn't poison the choice.
    let mut best_idx: usize = 0;
    let mut best_ts = coll.entries[0].timestep;
    for (i, entry) in coll.entries.iter().enumerate().skip(1) {
        let beats = entry
            .timestep
            .partial_cmp(&best_ts)
            .map(|o| !matches!(o, std::cmp::Ordering::Less))
            .unwrap_or(false);
        if beats {
            best_idx = i;
            best_ts = entry.timestep;
        }
    }
    // Round-19 M2: resolve_entry_path now returns `Result` to refuse
    // hostile `<DataSet file="../../etc/passwd"/>` traversal. The
    // mesh-loader's contract is "best-effort, return None on failure",
    // so we drop the path-escape error onto the floor exactly like
    // every other failure mode upstream.
    valenx_fields::pvd::resolve_entry_path(pvd_path, &coll.entries[best_idx]).ok()
}

/// Walk a workdir for `.vtu` (XML) or `.vtk` (legacy binary) files
/// and return the lexicographically latest path. For OpenFOAM's
/// `<case>_<N>.vtu` naming this maps to the highest-time-step
/// snapshot; for Code Aster's `RESU.vtk` it returns the unique file.
/// `None` if neither extension matches.
pub(crate) fn latest_vtk_in_workdir(workdir: &std::path::Path) -> Option<PathBuf> {
    let mut found: Vec<PathBuf> = Vec::new();
    let mut stack: Vec<PathBuf> = vec![workdir.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let Ok(entries) = std::fs::read_dir(&dir) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if let Ok(ft) = entry.file_type() {
                if ft.is_dir() {
                    stack.push(path);
                    continue;
                }
            }
            let ext_match = path
                .extension()
                .and_then(|s| s.to_str())
                .map(|s| {
                    let s = s.to_ascii_lowercase();
                    s == "vtu" || s == "vtk"
                })
                .unwrap_or(false);
            if ext_match {
                found.push(path);
            }
        }
    }
    found.sort();
    found.into_iter().next_back()
}

/// Read + parse a `.vtu` / `.vtk` file and convert its mesh half
/// into a canonical [`Mesh`]. Doesn't touch the field data — that
/// path is owned by the adapter's `collect()` pipeline (which
/// populates `Results.fields`); this function exists so the
/// post-run hook can populate the viewport without going through
/// the catalog.
///
/// Routes through `valenx_fields::vtk_dispatch` so both VTU XML
/// and VTK legacy binary files load via the same entry point.
///
/// Round-22 M1 (R20 L1 sister gap): pre-fix the read was a bare
/// `std::fs::read(path)` — a corrupted or adversarial workdir with
/// a multi-GB `.vtk` file would have slurped into memory before
/// `vtk_dispatch` saw the magic bytes. The CFD adapters
/// (SU2 / OpenFOAM) and the Elmer FEA adapter got this cap in
/// R20 L1 / R21 M3; the app-side auto-loader (which runs straight
/// after a run finishes, on the same workdir) shared the gap until
/// now. Capped at [`MAX_VTK_FILE_BYTES`] (4 GiB) to match the
/// adapter-side caps.
pub(crate) fn load_mesh_from_vtk(path: &std::path::Path) -> Result<Mesh, String> {
    let bytes =
        valenx_core::io_caps::read_capped_to_bytes(path, valenx_core::io_caps::MAX_VTK_FILE_BYTES)
            .map_err(|e| format!("read {}: {e}", path.display()))?;
    let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("vtk");
    let (mesh, _fields) = valenx_fields::vtk_dispatch::load_canonical(&bytes, stem)
        .map_err(|e| format!("parse {}: {e}", path.display()))?;
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn latest_vtk_in_workdir_picks_highest_timestep() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-app-vtu-walk-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let vtk_dir = tmp.join("VTK");
        std::fs::create_dir_all(&vtk_dir).unwrap();
        // OpenFOAM's typical naming: case_<N>.vtu. Lexicographic sort
        // puts the highest N last when N is zero-padded; the helper
        // doesn't pad, so we test with hand-picked values that sort
        // correctly anyway.
        std::fs::write(vtk_dir.join("case_100.vtu"), "x").unwrap();
        std::fs::write(vtk_dir.join("case_500.vtu"), "x").unwrap();
        std::fs::write(vtk_dir.join("case_50.vtu"), "x").unwrap();
        // Sneak in a non-vtu — should be ignored.
        std::fs::write(vtk_dir.join("notes.txt"), "noise").unwrap();
        let latest = latest_vtk_in_workdir(&tmp).expect("found a vtu");
        // Lexicographic last of {case_100, case_50, case_500} is
        // case_500 — exactly the OpenFOAM convention.
        assert!(
            latest
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .ends_with("case_500.vtu"),
            "got {}",
            latest.display()
        );
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn latest_vtk_in_workdir_returns_none_when_empty() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-app-vtu-walk-empty-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        assert!(latest_vtk_in_workdir(&tmp).is_none());
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn load_mesh_from_vtk_round_trips_a_one_tet() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-app-vtu-load-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let path = tmp.join("cavity_500.vtu");
        std::fs::write(
            &path,
            r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid">
<UnstructuredGrid><Piece NumberOfPoints="4" NumberOfCells="1">
  <Points><DataArray type="Float32" NumberOfComponents="3" format="ascii">
    0 0 0  1 0 0  0 1 0  0 0 1
  </DataArray></Points>
  <Cells>
    <DataArray Name="connectivity" format="ascii">0 1 2 3</DataArray>
    <DataArray Name="offsets" format="ascii">4</DataArray>
    <DataArray Name="types" format="ascii">10</DataArray>
  </Cells>
</Piece></UnstructuredGrid></VTKFile>"#,
        )
        .unwrap();
        let mesh = load_mesh_from_vtk(&path).expect("load mesh");
        assert_eq!(mesh.nodes.len(), 4);
        assert_eq!(mesh.element_blocks.len(), 1);
        assert_eq!(
            mesh.element_blocks[0].element_type,
            valenx_mesh::ElementType::Tet4
        );
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn load_mesh_from_vtk_handles_legacy_binary_format() {
        // The auto-load path needs to accept .vtk legacy binary
        // alongside .vtu XML so SU2 / Code Aster outputs land in
        // the viewport without manual intervention.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-vtk-legacy-load-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let path = tmp.join("RESU.vtk");
        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        buf.extend_from_slice(b"smoke legacy load\n");
        buf.extend_from_slice(b"BINARY\n");
        buf.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
        buf.extend_from_slice(b"POINTS 4 float\n");
        for v in [
            0.0_f32, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
        ] {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        buf.push(b'\n');
        buf.extend_from_slice(b"CELLS 1 5\n");
        for v in [4u32, 0, 1, 2, 3] {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        buf.push(b'\n');
        buf.extend_from_slice(b"CELL_TYPES 1\n");
        buf.extend_from_slice(&10u32.to_be_bytes());
        buf.push(b'\n');
        std::fs::write(&path, buf).unwrap();
        let mesh = load_mesh_from_vtk(&path).expect("load legacy");
        assert_eq!(mesh.nodes.len(), 4);
        assert_eq!(
            mesh.element_blocks[0].element_type,
            valenx_mesh::ElementType::Tet4
        );
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn latest_vtk_in_workdir_picks_vtk_over_an_older_vtu() {
        // Mixed-format workdirs (e.g. an older OpenFOAM run + a
        // newer Code Aster legacy export) should sort
        // lexicographically across BOTH extensions, not bias one
        // way. This test confirms .vtk and .vtu are both visible to
        // the picker.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-vtk-mixed-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::fs::write(tmp.join("a_older.vtu"), b"<placeholder/>").unwrap();
        std::fs::write(tmp.join("z_newer.vtk"), b"placeholder").unwrap();
        let latest = latest_vtk_in_workdir(&tmp).expect("found");
        assert_eq!(latest.file_name().unwrap(), "z_newer.vtk");
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn latest_snapshot_in_workdir_prefers_pvd_over_raw_vtu() {
        // Workdir contains both a `.pvd` manifest pointing at one VTU
        // file and a separate "newer" VTU (lex order). The PVD-aware
        // walker should pick the manifest's latest entry, not the
        // unrelated VTU.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-app-pvd-prefer-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        // Create a VTU "newer" by name to confirm the PVD wins on
        // priority alone, not on lex sort.
        std::fs::write(tmp.join("zzz_newer.vtu"), "<placeholder/>").unwrap();
        // PVD with two timesteps; the second (timestep=1.0) is what
        // should be picked.
        std::fs::write(tmp.join("step_a.vtu"), "<placeholder/>").unwrap();
        std::fs::write(tmp.join("step_b.vtu"), "<placeholder/>").unwrap();
        let pvd = r#"<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>
    <DataSet timestep="0.0" file="step_a.vtu"/>
    <DataSet timestep="1.0" file="step_b.vtu"/>
  </Collection>
</VTKFile>"#;
        std::fs::write(tmp.join("transient.pvd"), pvd).unwrap();
        let snap = latest_snapshot_in_workdir(&tmp).expect("found");
        assert_eq!(
            snap.file_name().unwrap(),
            "step_b.vtu",
            "got {}",
            snap.display()
        );
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn latest_snapshot_in_workdir_falls_back_to_vtk_walk_when_no_pvd() {
        // No `.pvd` present -> behaviour matches `latest_vtk_in_workdir`.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-app-pvd-fallback-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::fs::write(tmp.join("case_010.vtu"), "x").unwrap();
        std::fs::write(tmp.join("case_020.vtu"), "x").unwrap();
        let snap = latest_snapshot_in_workdir(&tmp).expect("found");
        assert_eq!(snap.file_name().unwrap(), "case_020.vtu");
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn latest_snapshot_in_workdir_returns_none_for_truly_empty() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-app-pvd-empty-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        assert!(latest_snapshot_in_workdir(&tmp).is_none());
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn latest_snapshot_falls_back_when_pvd_is_unparseable() {
        // PVD malformed -> the helper should warn + still return a
        // VTU/VTK from the fallback walk so the user isn't left with
        // a blank viewport.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-app-pvd-bad-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        std::fs::write(tmp.join("garbage.pvd"), "this is not xml").unwrap();
        std::fs::write(tmp.join("snap.vtu"), "x").unwrap();
        let snap = latest_snapshot_in_workdir(&tmp).expect("fallback");
        assert_eq!(snap.file_name().unwrap(), "snap.vtu");
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn latest_entry_from_pvd_picks_highest_timestep_even_when_unsorted() {
        // PVD entries declared out of order -> max-timestep wins.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-app-pvd-unsorted-{}.pvd",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let pvd = r#"<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>
    <DataSet timestep="2.5" file="mid.vtu"/>
    <DataSet timestep="10.0" file="late.vtu"/>
    <DataSet timestep="0.5" file="early.vtu"/>
  </Collection>
</VTKFile>"#;
        std::fs::write(&tmp, pvd).unwrap();
        let latest = latest_entry_from_pvd(&tmp).expect("parsed");
        assert_eq!(latest.file_name().unwrap(), "late.vtu");
        let _ = std::fs::remove_file(&tmp);
    }

    /// Round-18 M2 RED→GREEN: `latest_entry_from_pvd` must refuse a
    /// `.pvd` manifest larger than `MAX_PVD_FILE_BYTES` and return
    /// `None` (rather than slurping the multi-GB file into memory
    /// before the XML parser sees any of it).
    #[test]
    fn latest_entry_from_pvd_rejects_oversize_manifest() {
        use std::io::{Seek, SeekFrom, Write};
        let tmp = std::env::temp_dir().join(format!(
            "valenx-app-pvd-oversize-{}.pvd",
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
            // 1 byte past the PVD cap → metadata.len() reports oversize.
            f.seek(SeekFrom::Start(
                valenx_core::io_caps::MAX_PVD_FILE_BYTES as u64 + 1,
            ))
            .unwrap();
            f.write_all(b"x").unwrap();
        }
        let snap = latest_entry_from_pvd(&tmp);
        let _ = std::fs::remove_file(&tmp);
        assert!(snap.is_none(), "expected None for oversize manifest");
    }

    #[test]
    fn first_pvd_in_workdir_finds_one_in_subdir() {
        // PVD nested under a per-step subdirectory should still be
        // discoverable — the recursive walk handles that.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-app-pvd-nested-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let nested = tmp.join("results").join("vtk");
        std::fs::create_dir_all(&nested).unwrap();
        let pvd_path = nested.join("output.pvd");
        std::fs::write(
            &pvd_path,
            "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\"></VTKFile>",
        )
        .unwrap();
        let found = first_pvd_in_workdir(&tmp).expect("found");
        assert_eq!(found, pvd_path);
        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// Round-22 M1 RED→GREEN (R20 L1 sister gap closed): a `.vtk`
    /// file larger than `MAX_VTK_FILE_BYTES` (4 GiB) must be
    /// rejected as an IO error WITHOUT being slurped into memory.
    /// Pre-fix `load_mesh_from_vtk` did a bare `std::fs::read` and
    /// would have allocated the full file size before `vtk_dispatch`
    /// saw the magic bytes.
    ///
    /// Uses `set_len` to create a sparse over-cap file without
    /// writing 5 GiB of zeros on every CI run.
    #[test]
    fn load_mesh_from_vtk_rejects_oversize_file() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-app-vtu-r22m1-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let path = tmp.join("oversize.vtu");
        // Past the 4 GiB MAX_VTK_FILE_BYTES cap.
        let f = std::fs::File::create(&path).unwrap();
        f.set_len(valenx_core::io_caps::MAX_VTK_FILE_BYTES + 1)
            .unwrap();
        drop(f);
        let err = load_mesh_from_vtk(&path)
            .expect_err("round-22 M1: a 4 GiB+ vtk must be rejected as IO error before reading");
        assert!(
            err.starts_with("read "),
            "must surface as a read error, not parse: {err}"
        );
        let _ = std::fs::remove_dir_all(&tmp);
    }
}
