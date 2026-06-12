//! STEP (ISO 10303 AP203/AP242) import / export via `truck-stepio`.
//!
//! ## Overview
//!
//! `truck-stepio` provides a `Display`-trait-based writer
//! (`truck_stepio::out::CompleteStepDisplay`) and an `ruststep`-based
//! reader (`truck_stepio::in::Table::from_data_section` +
//! `Table::to_compressed_shell`). We thin-wrap both so the rest of
//! Valenx talks in terms of [`valenx_cad::Solid`].
//!
//! ## Round-trip semantics
//!
//! - Export uses [`truck_modeling::Solid::compress`] to turn a runtime
//!   [`truck_modeling::Solid`] into a [`truck_topology::compress::CompressedSolid`]
//!   that the writer wants, then formats with the standard schema
//!   header.
//! - Import parses the file with `ruststep`, picks the first
//!   `SHELL_BASED_SURFACE_MODEL` / `ADVANCED_BREP_SHAPE_REPRESENTATION`
//!   shell, converts it to a `CompressedShell`, and wraps it back into
//!   a `truck_modeling::Solid` via the standard `Solid::new`
//!   constructor.
//!
//! Multiple shells per file → we keep the first and emit a warning;
//! v1 doesn't merge them automatically.

use std::path::Path;

use valenx_cad::Solid;

use crate::error::StepIgesError;

/// Write `solid` to `path` as ISO 10303 AP203/AP214 STEP text.
///
/// # Errors
///
/// - [`StepIgesError::MeshBackedSolidNotExportable`] if `solid` is a
///   mesh-backed solid (no BRep topology).
/// - [`StepIgesError::EmptySolid`] if the underlying truck solid has
///   no boundary shells.
/// - [`StepIgesError::Io`] for filesystem failures.
/// - [`StepIgesError::Unsupported`] when the crate was built without
///   the `step` feature.
pub fn write(solid: &Solid, path: &Path) -> Result<(), StepIgesError> {
    #[cfg(not(feature = "step"))]
    {
        let _ = (solid, path);
        Err(StepIgesError::Unsupported {
            format: "STEP",
            feature: "truck-stepio dep gated by `step` feature".to_string(),
        })
    }
    #[cfg(feature = "step")]
    {
        let inner = match solid {
            Solid::Brep(b) => b,
            Solid::Mesh(_) => {
                return Err(StepIgesError::MeshBackedSolidNotExportable {
                    format: "STEP",
                    reason: "STEP needs BRep faces; mesh-backed solids have \
                             only triangles. Apply the fillet/chamfer last \
                             or export STL instead."
                        .to_string(),
                });
            }
        };
        if inner.boundaries().is_empty() {
            return Err(StepIgesError::EmptySolid);
        }
        let text = solid_to_step_text(inner, path);
        valenx_core::io_caps::atomic_write_str(path, &text)?;
        Ok(())
    }
}

#[cfg(feature = "step")]
fn solid_to_step_text(solid: &truck_modeling::Solid, path: &Path) -> String {
    use truck_stepio::out;

    let compressed = solid.compress();
    let header = out::StepHeaderDescriptor {
        file_name: crate::persist::basename(path),
        time_stamp: chrono_like_iso_now(),
        authors: vec!["valenx".to_string()],
        organization: vec!["valenx".to_string()],
        organization_system: "valenx-step-iges".to_string(),
        authorization: String::new(),
    };
    out::CompleteStepDisplay::new(out::StepModel::from(&compressed), header).to_string()
}

/// truck-stepio's `Default` impl pulls in `chrono::Utc::now()` for the
/// timestamp. We rebuild our own so we don't link to chrono directly
/// and stay testable.
#[cfg(feature = "step")]
fn chrono_like_iso_now() -> String {
    use std::time::SystemTime;
    let t = crate::persist::iges_timestamp(SystemTime::now());
    // STEP wants `YYYY-MM-DDTHH:MM:SS`; transform `YYYYMMDD.HHMMSS`.
    let bytes = t.as_bytes();
    if bytes.len() == 15 {
        format!(
            "{}-{}-{}T{}:{}:{}",
            std::str::from_utf8(&bytes[0..4]).unwrap_or("0000"),
            std::str::from_utf8(&bytes[4..6]).unwrap_or("00"),
            std::str::from_utf8(&bytes[6..8]).unwrap_or("00"),
            std::str::from_utf8(&bytes[9..11]).unwrap_or("00"),
            std::str::from_utf8(&bytes[11..13]).unwrap_or("00"),
            std::str::from_utf8(&bytes[13..15]).unwrap_or("00"),
        )
    } else {
        t
    }
}

/// Read a STEP file from `path` and return the first solid (or shell
/// promoted to a solid) it contains.
///
/// v1 picks the first shell when the file contains many; future work
/// could surface a per-shell list.
///
/// # Errors
///
/// - [`StepIgesError::Io`] for read failures.
/// - [`StepIgesError::ParseError`] for malformed STEP text or files
///   with no shells we can convert.
/// - [`StepIgesError::Unsupported`] when the crate was built without
///   the `step` feature.
pub fn read(path: &Path) -> Result<Solid, StepIgesError> {
    #[cfg(not(feature = "step"))]
    {
        let _ = path;
        Err(StepIgesError::Unsupported {
            format: "STEP",
            feature: "truck-stepio dep gated by `step` feature".to_string(),
        })
    }
    #[cfg(feature = "step")]
    {
        // Round-19 M1: migrate from stat-then-read to the shared
        // `read_capped_cad_text` helper that round-18 L1 introduced
        // for the IGES / IGES-trimmed / AP242 readers. The previous
        // stat-then-read shape was TOCTOU-vulnerable — a file that
        // raced past the cap between `metadata()` and `read_to_string()`
        // would bypass the cap entirely; the shared helper caps the
        // second read with `take(cap+1)` so a runtime-grown file gets
        // rejected at the read step too.
        let text = crate::read_capped_cad_text(path, "STEP")?;
        parse_step_text(&text)
    }
}

/// STEP import strategy
/// =====================
///
/// truck-stepio's `Table::to_compressed_shell` returns a
/// `CompressedShell<Point3, Curve3D, Surface>` where `Curve3D` /
/// `Surface` are stepio's own enum types — **not** the same Curve /
/// Surface that `truck_modeling::Solid` uses. There is no built-in
/// conversion, so v1 takes the pragmatic path: tessellate the imported
/// shell directly via truck-meshalgo's `robust_triangulation`, then
/// wrap the resulting polygon mesh as a [`Solid::Mesh`]. This loses
/// BRep topology (so re-export to STEP will fail with
/// [`StepIgesError::MeshBackedSolidNotExportable`]), but it covers the
/// 80% case: "import a STEP file from SolidWorks → see it in the
/// viewport → 3D-print or mesh-process it". True BRep round-trip is
/// Phase 8.5 work.
#[cfg(feature = "step")]
fn parse_step_text(text: &str) -> Result<Solid, StepIgesError> {
    use truck_meshalgo::tessellation::*;
    use truck_stepio::r#in::Table;

    if text.trim().is_empty() {
        return Err(StepIgesError::ParseError("empty input".into()));
    }
    let exchange = ruststep::parser::parse(text)
        .map_err(|e| StepIgesError::ParseError(format!("ruststep: {e}")))?;
    let data = exchange
        .data
        .first()
        .ok_or_else(|| StepIgesError::ParseError("no DATA section in STEP file".into()))?;
    let table = Table::from_data_section(data);
    if table.shell.is_empty() {
        return Err(StepIgesError::ParseError(
            "no shells (no `ADVANCED_BREP_SHAPE_REPRESENTATION` / \
             `MANIFOLD_SOLID_BREP` / `SHELL_BASED_SURFACE_MODEL`) in DATA"
                .into(),
        ));
    }
    let shell_count = table.shell.len();
    let first = table
        .shell
        .values()
        .next()
        .expect("checked non-empty above");
    if shell_count > 1 {
        tracing::warn!(
            target: "valenx-step-iges",
            "STEP file has {shell_count} shells; importing only the first",
        );
    }
    let cshell = table
        .to_compressed_shell(first)
        .map_err(|e| StepIgesError::ParseError(format!("to_compressed_shell: {e}")))?;
    // Tessellate at a v1 default tolerance — fine enough to keep
    // small features visible without bloating the viewport mesh.
    let poly_shell = cshell.robust_triangulation(0.5);
    let valenx_mesh = polyshell_to_valenx_mesh(&poly_shell);
    Ok(Solid::from_mesh(valenx_mesh))
}

/// Concatenate every per-face PolygonMesh in the triangulated shell
/// into a single valenx [`Mesh`] of Tri3 elements.
#[cfg(feature = "step")]
fn polyshell_to_valenx_mesh(
    poly_shell: &truck_topology::compress::CompressedShell<
        truck_modeling::Point3,
        truck_meshalgo::prelude::PolylineCurve<truck_modeling::Point3>,
        Option<truck_meshalgo::prelude::PolygonMesh>,
    >,
) -> valenx_mesh::Mesh {
    use valenx_mesh::{ElementBlock, ElementType, Mesh};

    let mut mesh = Mesh::new("step_import");
    let mut block = ElementBlock::new(ElementType::Tri3);

    for cface in &poly_shell.faces {
        let Some(poly) = &cface.surface else { continue };
        let base = mesh.nodes.len() as u32;
        for v in poly.positions() {
            mesh.nodes.push(nalgebra::Vector3::new(v[0], v[1], v[2]));
        }
        // Walk every Tri3 face (truck PolygonMesh has tri_faces() +
        // quad_faces() + other_faces() — we triangulate quads on the
        // fly into two Tri3s).
        for tri in poly.tri_faces().iter() {
            block.connectivity.extend_from_slice(&[
                base + tri[0].pos as u32,
                base + tri[1].pos as u32,
                base + tri[2].pos as u32,
            ]);
        }
        for quad in poly.quad_faces().iter() {
            block.connectivity.extend_from_slice(&[
                base + quad[0].pos as u32,
                base + quad[1].pos as u32,
                base + quad[2].pos as u32,
            ]);
            block.connectivity.extend_from_slice(&[
                base + quad[0].pos as u32,
                base + quad[2].pos as u32,
                base + quad[3].pos as u32,
            ]);
        }
    }
    if !block.connectivity.is_empty() {
        mesh.element_blocks.push(block);
    }
    mesh
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[cfg(feature = "step")]
    fn read_empty_file_returns_parse_error() {
        let tmp = std::env::temp_dir().join("valenx_step_empty.step");
        std::fs::write(&tmp, "").unwrap();
        let err = read(&tmp).unwrap_err();
        let _ = std::fs::remove_file(&tmp);
        assert!(matches!(err, StepIgesError::ParseError(_)));
    }

    #[test]
    #[cfg(feature = "step")]
    fn read_nonexistent_file_returns_io_error() {
        let path = std::env::temp_dir().join("valenx_step_nonexistent_xyz.step");
        let _ = std::fs::remove_file(&path);
        let err = read(&path).unwrap_err();
        assert!(matches!(err, StepIgesError::Io(_)), "got {err:?}");
    }

    #[test]
    #[cfg(feature = "step")]
    fn read_malformed_returns_parse_error() {
        let tmp = std::env::temp_dir().join("valenx_step_malformed.step");
        std::fs::write(&tmp, "this is not STEP text!!\nrandom bytes\n").unwrap();
        let err = read(&tmp).unwrap_err();
        let _ = std::fs::remove_file(&tmp);
        assert!(matches!(err, StepIgesError::ParseError(_)));
    }

    /// Round-19 M1 RED→GREEN: a .step file whose size exceeds
    /// `MAX_CAD_INTERCHANGE_FILE_BYTES` must produce `FileTooLarge`
    /// instead of slurping the whole payload into memory before the
    /// parser sees it. The fix migrated the `step::read` path off the
    /// stat-then-`fs::read_to_string` shape (TOCTOU-vulnerable to a
    /// file that grew between stat and read) onto the shared
    /// `read_capped_cad_text` helper round-18 L1 introduced — the
    /// same helper sister readers (iges, iges_trimmed, ap242) use.
    ///
    /// Sparse-file trick: `set_len(cap+1)` advertises a length past
    /// the cap without actually writing gigabytes to disk. NTFS /
    /// ext4 / APFS all report the set length to stat without
    /// allocating data blocks.
    #[test]
    #[cfg(feature = "step")]
    fn read_rejects_oversize_via_helper() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx_step_toolarge-{}.step",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let f = std::fs::File::create(&tmp).unwrap();
        f.set_len(crate::MAX_CAD_INTERCHANGE_FILE_BYTES + 1)
            .unwrap();
        drop(f);
        let err = read(&tmp).expect_err("oversize .step must error");
        let _ = std::fs::remove_file(&tmp);
        assert!(
            matches!(err, StepIgesError::FileTooLarge { format: "STEP", .. }),
            "wrong error variant: {err:?}"
        );
    }

    #[test]
    #[cfg(feature = "step")]
    fn write_mesh_backed_solid_is_rejected() {
        let mesh = valenx_mesh::Mesh::new("rejected");
        let s = valenx_cad::Solid::from_mesh(mesh);
        let tmp = std::env::temp_dir().join("valenx_step_meshbacked.step");
        let err = write(&s, &tmp).unwrap_err();
        let _ = std::fs::remove_file(&tmp);
        assert!(
            matches!(err, StepIgesError::MeshBackedSolidNotExportable { .. }),
            "got {err:?}",
        );
    }

    #[test]
    #[cfg(feature = "step")]
    fn write_box_round_trip_preserves_bounding_box() {
        // STEP export keeps BRep topology; STEP import returns a
        // mesh-backed solid (see `parse_step_text` for the rationale).
        // The bbox assertion confirms the round-trip preserves the
        // model's overall dimensions.
        let cube = valenx_cad::box_solid(10.0, 10.0, 10.0).unwrap();
        let tmp = std::env::temp_dir().join("valenx_step_roundtrip_box.step");
        write(&cube, &tmp).expect("write");
        let read_back = read(&tmp).expect("read");
        let _ = std::fs::remove_file(&tmp);

        let mesh_a = valenx_cad::solid_to_mesh(&cube, 1.0).unwrap();
        let mesh_b = valenx_cad::solid_to_mesh(&read_back, 1.0).unwrap();
        assert!(!mesh_a.nodes.is_empty(), "cube tessellates");
        assert!(!mesh_b.nodes.is_empty(), "read-back cube tessellates");

        let bbox = |m: &valenx_mesh::Mesh| {
            let mut mn = [f64::INFINITY; 3];
            let mut mx = [f64::NEG_INFINITY; 3];
            for n in &m.nodes {
                for i in 0..3 {
                    mn[i] = mn[i].min(n[i]);
                    mx[i] = mx[i].max(n[i]);
                }
            }
            [mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]]
        };
        let dims_a = bbox(&mesh_a);
        let dims_b = bbox(&mesh_b);
        for i in 0..3 {
            assert!(
                (dims_a[i] - dims_b[i]).abs() < 1.0,
                "bbox dim {i} mismatch: a={} b={}",
                dims_a[i],
                dims_b[i],
            );
        }
    }

    #[test]
    #[cfg(feature = "step")]
    fn write_cylinder_round_trip_exercises_curved_surface() {
        // A box has only planar faces; a cylinder exercises the curved-
        // surface path through truck-stepio's exporter + the parser's
        // CYLINDRICAL_SURFACE handler on read-back. Equivalent to the
        // Pad-replay test from the plan but without dragging in the
        // feature-tree crate (which has step-iges as a dev-dep, so
        // step-iges can't dev-dep on it without forming a cycle).
        let cyl = valenx_cad::cylinder(5.0, 10.0).expect("cylinder builds");
        let tmp = std::env::temp_dir().join("valenx_step_cylinder_roundtrip.step");
        write(&cyl, &tmp).expect("write");
        let read_back = read(&tmp).expect("read");
        let _ = std::fs::remove_file(&tmp);

        let mesh_a = valenx_cad::solid_to_mesh(&cyl, 0.5).unwrap();
        let mesh_b = valenx_cad::solid_to_mesh(&read_back, 0.5).unwrap();
        assert!(!mesh_a.nodes.is_empty());
        assert!(!mesh_b.nodes.is_empty());
    }
}
