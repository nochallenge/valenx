//! Slicer-bundle export.
//!
//! v1 defers actual G-code generation to the user's existing slicer
//! (Cura, PrusaSlicer, etc.) and emits a directory of STL files plus
//! a manifest the slicer can consume. The slicer is invoked
//! externally — Valenx never spawns it (lockdown).

use std::fs;
use std::path::Path;

use valenx_mesh::write_stl_binary;

use crate::error::PrintBedError;
use crate::printer::{Part, Printer};

/// Caller-supplied slicer hints written to the manifest. Format is
/// intentionally opaque so any slicer's settings dialect fits.
#[derive(Clone, Debug, Default)]
pub struct SlicerSettings {
    /// Free-form lines (e.g. `"layer_height: 0.2"`, `"print_speed: 50"`).
    pub lines: Vec<String>,
}

/// Write every part's mesh to `path/<name>.stl` plus a `manifest.txt`
/// listing the parts, the printer info, and the slicer settings.
pub fn export_layout(
    parts: &[Part],
    printer: &Printer,
    slicer_settings: &SlicerSettings,
    path: impl AsRef<Path>,
) -> Result<(), PrintBedError> {
    let path = path.as_ref();
    fs::create_dir_all(path).map_err(|e| PrintBedError::Io {
        path: path.display().to_string(),
        reason: e.to_string(),
    })?;
    let mut manifest = String::new();
    manifest.push_str(&format!(
        "# Valenx Print Bed Layout\nbed_size: {:?}\nbed_type: {:?}\nbed_material: {:?}\n",
        printer.bed_size, printer.bed_type, printer.bed_material
    ));
    manifest.push_str("\n## Parts\n");
    for p in parts {
        // Round-10 H2 fix: `Part::name` is an unvalidated `String`
        // that flowed straight into `path.join(format!("{}.stl", ...))`.
        // A hostile / mistyped name like `"../etc/passwd"` wrote the
        // STL outside the export directory. Validate it as a
        // basename — no `/`, no `\`, no `..`, no absolute paths —
        // before the join, using the same canonical helper every
        // adapter goes through.
        valenx_core::adapter_helpers::validate_output_basename(&p.name, "Part.name").map_err(
            |e| PrintBedError::BadParameter {
                name: "Part.name",
                reason: format!("{e}"),
            },
        )?;
        let stl_path = path.join(format!("{}.stl", p.name));
        // valenx-mesh's STL writer takes &Mesh + Path.
        write_stl_binary(&p.mesh, &stl_path).map_err(|e| PrintBedError::StlWrite {
            path: stl_path.display().to_string(),
            reason: e.to_string(),
        })?;
        manifest.push_str(&format!(
            "- name: {}\n  stl: {}\n  bed_position: [{:.3}, {:.3}]\n",
            p.name,
            stl_path.display(),
            p.bed_position[0],
            p.bed_position[1],
        ));
    }
    if !slicer_settings.lines.is_empty() {
        manifest.push_str("\n## Slicer settings\n");
        for line in &slicer_settings.lines {
            manifest.push_str(&format!("{line}\n"));
        }
    }
    let manifest_path = path.join("manifest.txt");
    valenx_core::io_caps::atomic_write_str(&manifest_path, &manifest).map_err(|e| {
        PrintBedError::Io {
            path: manifest_path.display().to_string(),
            reason: e.to_string(),
        }
    })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::printer::{BedMaterial, BedType};
    use nalgebra::Vector3;
    use valenx_mesh::element::{ElementBlock, ElementType};
    use valenx_mesh::Mesh;

    fn tmp_outdir(tag: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "valenx_print_bed_{}_{}_{}",
            tag,
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_nanos())
                .unwrap_or(0)
        ));
        let _ = fs::remove_dir_all(&dir);
        dir
    }

    fn triangle_part(name: &str) -> Part {
        let mut m = Mesh::new(name);
        m.nodes.push(Vector3::zeros());
        m.nodes.push(Vector3::new(10.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(0.0, 10.0, 0.0));
        let mut b = ElementBlock::new(ElementType::Tri3);
        b.connectivity.extend_from_slice(&[0, 1, 2]);
        m.element_blocks.push(b);
        m.recompute_stats();
        Part::new(name, m)
    }

    #[test]
    fn export_creates_stl_and_manifest() {
        let dir = tmp_outdir("export");
        let printer = Printer::new((220.0, 220.0, 250.0), BedType::Heated, BedMaterial::Pei);
        let parts = vec![triangle_part("a"), triangle_part("b")];
        let settings = SlicerSettings {
            lines: vec!["layer_height: 0.2".into()],
        };
        export_layout(&parts, &printer, &settings, &dir).unwrap();
        assert!(dir.join("a.stl").exists());
        assert!(dir.join("b.stl").exists());
        assert!(dir.join("manifest.txt").exists());
    }

    /// Round-10 H2 RED→GREEN: `Part::name` is unvalidated `String`
    /// that pre-fix flowed straight into `path.join(format!("{}.stl",
    /// p.name))`. A hostile name `"../etc/passwd"` wrote the STL
    /// outside `path`. The export now rejects the part before any
    /// STL byte is written.
    #[test]
    fn export_rejects_part_name_path_traversal() {
        let dir = tmp_outdir("traversal");
        let printer = Printer::new((220.0, 220.0, 250.0), BedType::Heated, BedMaterial::Pei);
        let parts = vec![triangle_part("../etc/passwd")];
        let settings = SlicerSettings::default();
        let err = export_layout(&parts, &printer, &settings, &dir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("Part.name"),
            "expected Part.name in error, got: {msg}"
        );
        // The STL must NOT have leaked outside the export dir.
        assert!(!dir.parent().unwrap().join("etc").exists());
    }

    #[test]
    fn export_rejects_part_name_with_absolute_path() {
        let dir = tmp_outdir("abs");
        let printer = Printer::new((220.0, 220.0, 250.0), BedType::Heated, BedMaterial::Pei);
        // Use a platform-flavoured absolute path: `/foo` is treated
        // as drive-root on Windows by Path::is_absolute → not quite,
        // but the basename validator also rejects `\` separators.
        let evil_name = if cfg!(windows) { "C:\\evil" } else { "/evil" };
        let parts = vec![triangle_part(evil_name)];
        let settings = SlicerSettings::default();
        let err = export_layout(&parts, &printer, &settings, &dir).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("Part.name"),
            "expected Part.name in error, got: {msg}"
        );
    }
}
