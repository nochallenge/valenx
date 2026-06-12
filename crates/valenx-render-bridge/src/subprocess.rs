//! Phase 30.5 — real subprocess adapters for the external renderers.
//!
//! [`crate::emit`] produces the scene-file *text* for Cycles / LuxRender /
//! POV-Ray but stops there — the caller had to launch the renderer
//! itself. This module closes that gap with genuine subprocess
//! adapters for **Cycles** and **LuxCoreRender**:
//!
//! - Cycles: [`emit::scene_file`](crate::emit::scene_file) already
//!   writes a Cycles XML scene; [`run_cycles`] writes it to a temp
//!   file and invokes the standalone `cycles` renderer on it.
//! - LuxCoreRender: [`emit_luxcore_scn`] + [`emit_luxcore_cfg`]
//!   serialise the job to LuxCore's modern **`.scn` (scene) + `.cfg`
//!   (render config)** SDL pair (distinct from the legacy `.lxs`
//!   `emit` produces for the classic LuxRender). [`run_luxcore`]
//!   writes both files and invokes `luxcoreconsole`.
//!
//! Every adapter is honest about a missing renderer: if the
//! executable is not on `PATH` it returns
//! [`RenderError::ToolNotAvailable`] naming the tool — never a silent
//! no-op, never a half-written output.
//!
//! ## Testing note
//!
//! The crate's tests cover the **scene-file serialisation** and the
//! **argument-vector construction** only. No test launches a
//! renderer — [`run_cycles`] / [`run_luxcore`] spawn a subprocess and
//! are exercised interactively.

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::engine::RenderEngine;
use crate::error::RenderError;
use crate::scene::RenderJob;

/// Outcome of a successful render dispatch — the path the renderer
/// wrote its output image to, plus the scene file(s) the adapter
/// generated (kept so a caller can inspect / archive them).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RenderArtifacts {
    /// The rendered image path (`job.output_path`).
    pub image_path: PathBuf,
    /// Every scene / config file the adapter wrote, in the order it
    /// wrote them.
    pub scene_files: Vec<PathBuf>,
}

// ===========================================================================
// LuxCoreRender `.scn` + `.cfg` serialisation
// ===========================================================================

/// Serialise `job` to a LuxCoreRender **scene** (`.scn`) SDL string.
///
/// LuxCore's modern format splits the description in two: the `.scn`
/// holds the camera, materials, and geometry (`scene.camera.*`,
/// `scene.materials.*`, `scene.objects.*`); the companion `.cfg`
/// ([`emit_luxcore_cfg`]) holds the render settings. This is the
/// format `luxcoreconsole` / the LuxCore SDK consume — *not* the
/// legacy `.lxs` that [`crate::emit::scene_file`] emits for classic
/// LuxRender.
///
/// `geometry_dir` is the directory the per-object `.ply` mesh files
/// are written into (LuxCore references external geometry by path);
/// the returned `.scn` text points `scene.objects.*.ply` at
/// `<geometry_dir>/<object>.ply`.
///
/// # Errors
///
/// [`RenderError::EmptyScene`] / [`RenderError::BadParameter`] via
/// [`RenderJob::validate`].
pub fn emit_luxcore_scn(job: &RenderJob, geometry_dir: &Path) -> Result<String, RenderError> {
    use std::fmt::Write;
    job.validate()?;
    let mut s = String::new();
    let _ = writeln!(s, "# LuxCoreRender scene — emitted by valenx-render-bridge");
    let _ = writeln!(s, "# Title: {}", job.title);

    // --- Camera. ---
    let cam = &job.camera;
    let _ = writeln!(
        s,
        "scene.camera.lookat.orig = {:.6} {:.6} {:.6}",
        cam.position.x, cam.position.y, cam.position.z
    );
    let _ = writeln!(
        s,
        "scene.camera.lookat.target = {:.6} {:.6} {:.6}",
        cam.target.x, cam.target.y, cam.target.z
    );
    let _ = writeln!(
        s,
        "scene.camera.up = {:.6} {:.6} {:.6}",
        cam.up.x, cam.up.y, cam.up.z
    );
    // LuxCore's field-of-view is the *horizontal* fov in degrees.
    let fov_h = 2.0 * ((cam.fov_v_rad * 0.5).tan() * cam.aspect()).atan();
    let _ = writeln!(s, "scene.camera.fieldofview = {:.6}", fov_h.to_degrees());

    // --- Materials → LuxCore `matte` / `metal` / `glass`. ---
    for (id, mat) in &job.materials {
        let lux_id = sanitize(id);
        if mat.metallic > 0.5 {
            // Metallic → LuxCore `metal2` with the specular tint.
            let _ = writeln!(s, "scene.materials.{lux_id}.type = metal2");
            let _ = writeln!(
                s,
                "scene.materials.{lux_id}.n = {:.4} {:.4} {:.4}",
                mat.specular_color[0], mat.specular_color[1], mat.specular_color[2]
            );
            let _ = writeln!(
                s,
                "scene.materials.{lux_id}.uroughness = {:.4}",
                mat.roughness.max(1e-3)
            );
            let _ = writeln!(
                s,
                "scene.materials.{lux_id}.vroughness = {:.4}",
                mat.roughness.max(1e-3)
            );
        } else if mat.roughness < 1e-3 && mat.ior > 1.05 {
            // Smooth dielectric → LuxCore `glass`.
            let _ = writeln!(s, "scene.materials.{lux_id}.type = glass");
            let _ = writeln!(
                s,
                "scene.materials.{lux_id}.kt = {:.4} {:.4} {:.4}",
                mat.diffuse_color[0], mat.diffuse_color[1], mat.diffuse_color[2]
            );
            let _ = writeln!(s, "scene.materials.{lux_id}.interiorior = {:.4}", mat.ior);
        } else {
            // Everything else → matte (Lambertian) diffuse.
            let _ = writeln!(s, "scene.materials.{lux_id}.type = matte");
            let _ = writeln!(
                s,
                "scene.materials.{lux_id}.kd = {:.4} {:.4} {:.4}",
                mat.diffuse_color[0], mat.diffuse_color[1], mat.diffuse_color[2]
            );
        }
        // Emissive surfaces become area lights.
        if mat.emissive.iter().any(|&c| c > 0.0) {
            let _ = writeln!(
                s,
                "scene.materials.{lux_id}.emission = {:.4} {:.4} {:.4}",
                mat.emissive[0], mat.emissive[1], mat.emissive[2]
            );
        }
    }

    // --- Objects: one per mesh, geometry by external `.ply` path. ---
    for m in &job.meshes {
        let obj_id = sanitize(&m.name);
        let ply_path = geometry_dir.join(format!("{obj_id}.ply"));
        let _ = writeln!(
            s,
            "scene.objects.{obj_id}.material = {}",
            sanitize(&m.material_id)
        );
        let _ = writeln!(s, "scene.objects.{obj_id}.ply = {}", ply_path.display());
    }

    // --- Environment / sky. LuxCore needs at least one light. ---
    if let Some(env) = &job.environment {
        let _ = writeln!(s, "scene.lights.env.type = infinite");
        let _ = writeln!(s, "scene.lights.env.file = {}", env.hdr_path.display());
        let _ = writeln!(
            s,
            "scene.lights.env.gain = {0:.4} {0:.4} {0:.4}",
            env.intensity
        );
    } else {
        // A neutral constant-colour sky so the scene is not black.
        let _ = writeln!(s, "scene.lights.sky.type = constantinfinite");
        let _ = writeln!(s, "scene.lights.sky.color = 0.8 0.85 0.95");
    }
    Ok(s)
}

/// Serialise the LuxCoreRender **render config** (`.cfg`) for `job`.
///
/// The `.cfg` holds the engine settings — render engine
/// (`PATHCPU`), film resolution, halt condition — and points at the
/// `.scn` via `scene.file`. `scn_path` is the path the matching
/// [`emit_luxcore_scn`] output will be written to.
///
/// # Errors
///
/// [`RenderError::EmptyScene`] / [`RenderError::BadParameter`] via
/// [`RenderJob::validate`].
pub fn emit_luxcore_cfg(job: &RenderJob, scn_path: &Path) -> Result<String, RenderError> {
    use std::fmt::Write;
    job.validate()?;
    let mut s = String::new();
    let _ = writeln!(
        s,
        "# LuxCoreRender config — emitted by valenx-render-bridge"
    );
    let _ = writeln!(s, "scene.file = {}", scn_path.display());
    let _ = writeln!(s, "renderengine.type = PATHCPU");
    let cam = &job.camera;
    let _ = writeln!(s, "film.width = {}", cam.image_width);
    let _ = writeln!(s, "film.height = {}", cam.image_height);
    // Halt after a fixed sample budget so the render terminates.
    let _ = writeln!(s, "batch.haltspp = 128");
    let _ = writeln!(s, "film.outputs.0.type = RGBA_IMAGEPIPELINE");
    let _ = writeln!(s, "film.outputs.0.filename = {}", job.output_path.display());
    Ok(s)
}

// ===========================================================================
// Subprocess command construction
// ===========================================================================

/// Build the `cycles` standalone-renderer argument vector for a scene
/// XML at `scene_xml` writing its image to `output`.
///
/// Returns the argument list **excluding** the program name. The
/// standalone Cycles renderer takes the scene XML positionally plus
/// `--output` for the image and `--background` for headless runs.
///
/// Exposed (and unit-tested) separately from [`run_cycles`] so the
/// command is verifiable without launching a subprocess.
pub fn cycles_command(scene_xml: &Path, output: &Path) -> Vec<String> {
    vec![
        "--background".into(), // headless, no GUI window
        "--output".into(),
        output.display().to_string(),
        "--samples".into(),
        "128".into(),
        scene_xml.display().to_string(), // scene file, positional
    ]
}

/// Build the `luxcoreconsole` argument vector for a render config at
/// `cfg_path`.
///
/// Returns the argument list **excluding** the program name.
/// `luxcoreconsole` (the headless LuxCore CLI) takes the `.cfg`
/// positionally; the `.cfg` itself references the `.scn` and the
/// output path, so no other flags are needed for a batch render.
pub fn luxcore_command(cfg_path: &Path) -> Vec<String> {
    vec![cfg_path.display().to_string()]
}

// ===========================================================================
// Subprocess adapters
// ===========================================================================

/// Render `job` with the standalone **Cycles** renderer.
///
/// Writes the Cycles XML scene (via [`crate::emit::scene_file`]) into
/// `work_dir`, then invokes the `cycles` executable on it. The
/// rendered image lands at `job.output_path`.
///
/// `job.engine` must be [`RenderEngine::Cycles`].
///
/// # Errors
///
/// - [`RenderError::BadParameter`] if `job.engine` is not Cycles.
/// - [`RenderError::ToolNotAvailable`] if `cycles` is not on `PATH`.
/// - [`RenderError::RendererFailed`] if Cycles exits non-zero.
/// - [`RenderError::Io`] for a file-write or spawn failure.
///
/// Launches a subprocess — not exercised by the crate tests.
pub fn run_cycles(job: &RenderJob, work_dir: &Path) -> Result<RenderArtifacts, RenderError> {
    if job.engine != RenderEngine::Cycles {
        return Err(RenderError::BadParameter {
            name: "engine",
            reason: format!(
                "run_cycles requires RenderEngine::Cycles, got {:?}",
                job.engine
            ),
        });
    }
    let scene_text = crate::emit::scene_file(job)?;
    let scene_path = work_dir.join("valenx_cycles_scene.xml");
    // R29 H1: route through the canonical tmp+fsync+rename helper rather
    // than a bare std::fs::write so a crash mid-write can't leave a
    // truncated scene file the renderer then chokes on.
    valenx_core::io_caps::atomic_write_str(&scene_path, &scene_text).map_err(RenderError::Io)?;

    let args = cycles_command(&scene_path, &job.output_path);
    let status = spawn_renderer("cycles", &args)?;
    if !status.success() {
        return Err(RenderError::RendererFailed {
            tool: "cycles",
            detail: status.to_string(),
        });
    }
    Ok(RenderArtifacts {
        image_path: job.output_path.clone(),
        scene_files: vec![scene_path],
    })
}

/// Render `job` with **LuxCoreRender** (`luxcoreconsole`).
///
/// Writes the `.scn` scene + `.cfg` config pair (via
/// [`emit_luxcore_scn`] / [`emit_luxcore_cfg`]) plus the per-object
/// `.ply` geometry files into `work_dir`, then invokes
/// `luxcoreconsole` on the `.cfg`. The rendered image lands at
/// `job.output_path`.
///
/// `job.engine` must be [`RenderEngine::LuxRender`] (Valenx's single
/// Lux engine variant covers LuxCoreRender).
///
/// # Errors
///
/// - [`RenderError::BadParameter`] if `job.engine` is not LuxRender.
/// - [`RenderError::ToolNotAvailable`] if `luxcoreconsole` is not on
///   `PATH`.
/// - [`RenderError::RendererFailed`] if LuxCore exits non-zero.
/// - [`RenderError::Io`] for a file-write or spawn failure.
///
/// Launches a subprocess — not exercised by the crate tests.
pub fn run_luxcore(job: &RenderJob, work_dir: &Path) -> Result<RenderArtifacts, RenderError> {
    if job.engine != RenderEngine::LuxRender {
        return Err(RenderError::BadParameter {
            name: "engine",
            reason: format!(
                "run_luxcore requires RenderEngine::LuxRender, got {:?}",
                job.engine
            ),
        });
    }
    let scn_path = work_dir.join("valenx_luxcore_scene.scn");
    let cfg_path = work_dir.join("valenx_luxcore_render.cfg");

    let scn_text = emit_luxcore_scn(job, work_dir)?;
    let cfg_text = emit_luxcore_cfg(job, &scn_path)?;
    // R29 H1: canonical atomic_write (tmp+fsync+rename) for the scene /
    // config pair and each .ply, replacing bare std::fs::write so a
    // crash can't leave LuxCore a half-written input.
    valenx_core::io_caps::atomic_write_str(&scn_path, &scn_text).map_err(RenderError::Io)?;
    valenx_core::io_caps::atomic_write_str(&cfg_path, &cfg_text).map_err(RenderError::Io)?;

    // LuxCore references geometry by external `.ply` — write one per
    // mesh next to the scene file.
    let mut written = vec![scn_path.clone(), cfg_path.clone()];
    for m in &job.meshes {
        let ply_path = work_dir.join(format!("{}.ply", sanitize(&m.name)));
        valenx_core::io_caps::atomic_write_str(&ply_path, &ascii_ply(&m.mesh))
            .map_err(RenderError::Io)?;
        written.push(ply_path);
    }

    let args = luxcore_command(&cfg_path);
    let status = spawn_renderer("luxcoreconsole", &args)?;
    if !status.success() {
        return Err(RenderError::RendererFailed {
            tool: "luxcoreconsole",
            detail: status.to_string(),
        });
    }
    Ok(RenderArtifacts {
        image_path: job.output_path.clone(),
        scene_files: written,
    })
}

/// Spawn `tool` with `args`, mapping a missing executable onto
/// [`RenderError::ToolNotAvailable`].
fn spawn_renderer(
    tool: &'static str,
    args: &[String],
) -> Result<std::process::ExitStatus, RenderError> {
    match Command::new(tool)
        .args(args)
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
    {
        Ok(status) => Ok(status),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            Err(RenderError::ToolNotAvailable { tool })
        }
        Err(e) => Err(RenderError::Io(e)),
    }
}

/// Probe whether the renderer for `engine` is installed on `PATH`.
///
/// Runs the renderer's `--help` (or `--version`) with all I/O
/// suppressed; returns `true` only if the process spawns. Lets a UI
/// grey out an engine ahead of time. Launches a subprocess — not
/// called by the crate tests.
pub fn renderer_available(engine: RenderEngine) -> bool {
    let tool = match engine {
        RenderEngine::Cycles => "cycles",
        RenderEngine::LuxRender => "luxcoreconsole",
        RenderEngine::PovRay => "povray",
        RenderEngine::Native => return true, // in-process, always "available"
    };
    Command::new(tool)
        .arg("--help")
        .stdin(std::process::Stdio::null())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .is_ok()
}

// ===========================================================================
// Helpers
// ===========================================================================

/// Replace characters LuxCore's SDL identifier grammar dislikes with
/// `_` so a user-supplied mesh / material name is always a valid id.
fn sanitize(name: &str) -> String {
    if name.is_empty() {
        return "unnamed".to_string();
    }
    name.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect()
}

/// Serialise a [`valenx_mesh::Mesh`] to a minimal ASCII PLY — the
/// geometry container LuxCore loads per object. Vertices + `Tri3`
/// faces only.
fn ascii_ply(mesh: &valenx_mesh::Mesh) -> String {
    use std::fmt::Write;
    let mut tris: Vec<[u32; 3]> = Vec::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::element::ElementType::Tri3 {
            continue;
        }
        for c in block.connectivity.chunks_exact(3) {
            tris.push([c[0], c[1], c[2]]);
        }
    }
    let mut s = String::new();
    let _ = writeln!(s, "ply");
    let _ = writeln!(s, "format ascii 1.0");
    let _ = writeln!(s, "comment emitted by valenx-render-bridge");
    let _ = writeln!(s, "element vertex {}", mesh.nodes.len());
    let _ = writeln!(s, "property float x");
    let _ = writeln!(s, "property float y");
    let _ = writeln!(s, "property float z");
    let _ = writeln!(s, "element face {}", tris.len());
    let _ = writeln!(s, "property list uchar int vertex_indices");
    let _ = writeln!(s, "end_header");
    for v in &mesh.nodes {
        let _ = writeln!(s, "{:.6} {:.6} {:.6}", v.x, v.y, v.z);
    }
    for t in &tris {
        let _ = writeln!(s, "3 {} {} {}", t[0], t[1], t[2]);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::camera::Camera;
    use crate::material::Material;
    use nalgebra::Vector3;
    use std::path::PathBuf;

    fn unit_tri_mesh() -> valenx_mesh::Mesh {
        let mut m = valenx_mesh::Mesh::new("t");
        m.nodes = vec![
            Vector3::zeros(),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut block =
            valenx_mesh::element::ElementBlock::new(valenx_mesh::element::ElementType::Tri3);
        block.connectivity = vec![0, 1, 2];
        m.element_blocks.push(block);
        m
    }

    fn job(engine: RenderEngine) -> RenderJob {
        let mut job = RenderJob::new();
        job.engine = engine;
        job.camera = Camera::default();
        job.output_path = PathBuf::from("out.png");
        let id = job.add_material(Material::matte("paint", [0.5, 0.6, 0.7]));
        job.add_mesh("widget", unit_tri_mesh(), id);
        job
    }

    #[test]
    fn luxcore_scn_has_camera_materials_objects() {
        let scn = emit_luxcore_scn(&job(RenderEngine::LuxRender), &PathBuf::from("/tmp")).unwrap();
        assert!(scn.contains("scene.camera.lookat.orig"));
        assert!(scn.contains("scene.camera.fieldofview"));
        // The matte material block.
        assert!(scn.contains("scene.materials.paint.type = matte"));
        assert!(scn.contains("scene.materials.paint.kd"));
        // The object references the material + an external .ply.
        assert!(scn.contains("scene.objects.widget.material = paint"));
        assert!(scn.contains("scene.objects.widget.ply"));
        // No environment → a constant sky so the render is not black.
        assert!(scn.contains("constantinfinite"));
    }

    #[test]
    fn luxcore_scn_maps_metal_to_metal2() {
        let mut j = job(RenderEngine::LuxRender);
        let id = j.add_material(Material::polished_metal("chrome", [0.9, 0.9, 0.9]));
        j.add_mesh("ring", unit_tri_mesh(), id);
        let scn = emit_luxcore_scn(&j, &PathBuf::from("/tmp")).unwrap();
        assert!(scn.contains("scene.materials.chrome.type = metal2"));
    }

    #[test]
    fn luxcore_scn_maps_glass_to_glass() {
        let mut j = job(RenderEngine::LuxRender);
        let id = j.add_material(Material::glass("pane"));
        j.add_mesh("window", unit_tri_mesh(), id);
        let scn = emit_luxcore_scn(&j, &PathBuf::from("/tmp")).unwrap();
        assert!(scn.contains("scene.materials.pane.type = glass"));
        assert!(scn.contains("scene.materials.pane.interiorior"));
    }

    #[test]
    fn luxcore_cfg_points_at_the_scn_and_output() {
        let j = job(RenderEngine::LuxRender);
        let cfg = emit_luxcore_cfg(&j, &PathBuf::from("/tmp/s.scn")).unwrap();
        assert!(cfg.contains("scene.file = "));
        assert!(cfg.contains("s.scn"));
        assert!(cfg.contains("renderengine.type = PATHCPU"));
        assert!(cfg.contains("film.width = 1920"));
        assert!(cfg.contains("film.height = 1080"));
        assert!(cfg.contains("out.png"));
    }

    #[test]
    fn luxcore_scn_rejects_empty_scene() {
        let empty = RenderJob::new();
        assert!(emit_luxcore_scn(&empty, &PathBuf::from("/tmp")).is_err());
    }

    #[test]
    fn cycles_command_is_headless_with_output_and_scene() {
        let cmd = cycles_command(&PathBuf::from("scene.xml"), &PathBuf::from("out.png"));
        assert!(cmd.iter().any(|a| a == "--background"));
        assert!(cmd.iter().any(|a| a == "--output"));
        // The scene file appears as a (positional) argument.
        assert!(cmd.iter().any(|a| a.ends_with("scene.xml")));
        assert!(cmd.iter().any(|a| a.ends_with("out.png")));
    }

    #[test]
    fn luxcore_command_passes_the_cfg() {
        let cmd = luxcore_command(&PathBuf::from("render.cfg"));
        assert_eq!(cmd.len(), 1);
        assert!(cmd[0].ends_with("render.cfg"));
    }

    #[test]
    fn run_cycles_rejects_wrong_engine() {
        // Must fail on the engine check BEFORE any subprocess work.
        let j = job(RenderEngine::LuxRender);
        let err = run_cycles(&j, &PathBuf::from(".")).unwrap_err();
        assert_eq!(err.code(), "render.bad_parameter");
    }

    #[test]
    fn run_luxcore_rejects_wrong_engine() {
        let j = job(RenderEngine::Cycles);
        let err = run_luxcore(&j, &PathBuf::from(".")).unwrap_err();
        assert_eq!(err.code(), "render.bad_parameter");
    }

    #[test]
    fn tool_not_available_error_codes() {
        // The error variant codes the subprocess adapters raise.
        let e = RenderError::ToolNotAvailable { tool: "cycles" };
        assert_eq!(e.code(), "render.tool_not_available");
        let e = RenderError::RendererFailed {
            tool: "luxcoreconsole",
            detail: "exit 1".into(),
        };
        assert_eq!(e.code(), "render.renderer_failed");
    }

    #[test]
    fn ascii_ply_round_trips_a_triangle() {
        let ply = ascii_ply(&unit_tri_mesh());
        assert!(ply.starts_with("ply\n"));
        assert!(ply.contains("element vertex 3"));
        assert!(ply.contains("element face 1"));
        assert!(ply.contains("3 0 1 2"));
    }

    #[test]
    fn sanitize_replaces_sdl_unfriendly_chars() {
        assert_eq!(sanitize("my mesh.001"), "my_mesh_001");
        assert_eq!(sanitize(""), "unnamed");
        assert_eq!(sanitize("ok_name"), "ok_name");
    }

    // NOTE: run_cycles / run_luxcore / renderer_available against a
    // real renderer are intentionally NOT tested — they spawn a
    // subprocess. UI/subprocess-coupled — run interactively only.
}
