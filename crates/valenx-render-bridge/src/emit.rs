//! Scene-file emitters per render engine.
//!
//! v1 emits a minimal-but-valid file for LuxRender / Cycles /
//! POV-Ray and refuses Native (callers handle Native via the wgpu
//! viewport screenshot path).

use std::fmt::Write;

use crate::engine::RenderEngine;
use crate::error::RenderError;
use crate::scene::{RenderJob, SceneMesh};

/// Emit the scene-file text for `job.engine`. The caller writes the
/// returned string to disk and dispatches the renderer subprocess.
pub fn scene_file(job: &RenderJob) -> Result<String, RenderError> {
    job.validate()?;
    match job.engine {
        RenderEngine::LuxRender => Ok(emit_lux(job)),
        RenderEngine::Cycles => Ok(emit_cycles(job)),
        RenderEngine::PovRay => Ok(emit_povray(job)),
        RenderEngine::Native => Err(RenderError::EngineNotImplemented {
            engine: "Native",
            reason: "Native engine renders via the desktop shell's wgpu viewport \
                     screenshot path — no scene file applies"
                .into(),
        }),
    }
}

/// Emit a LuxCoreRender / LuxRender `.lxs` scene.
///
/// v1 emits the global block + camera + film + a one-line
/// `WorldBegin` … `WorldEnd` shell with per-mesh `Shape "trianglemesh"`
/// entries. No texture maps; flat material maps to `Matte` BSDF.
fn emit_lux(job: &RenderJob) -> String {
    let mut s = String::new();
    let _ = writeln!(
        s,
        "# LuxCoreRender scene — emitted by valenx-render-bridge v1"
    );
    let _ = writeln!(s, "# Title: {}", job.title);
    let cam = &job.camera;
    let _ = writeln!(s, "Film \"fleximage\"");
    let _ = writeln!(s, "  \"integer xresolution\" [{}]", cam.image_width);
    let _ = writeln!(s, "  \"integer yresolution\" [{}]", cam.image_height);
    let _ = writeln!(
        s,
        "  \"string filename\" [\"{}\"]",
        job.output_path.display()
    );
    let _ = writeln!(
        s,
        "LookAt {:.6} {:.6} {:.6}  {:.6} {:.6} {:.6}  {:.6} {:.6} {:.6}",
        cam.position.x,
        cam.position.y,
        cam.position.z,
        cam.target.x,
        cam.target.y,
        cam.target.z,
        cam.up.x,
        cam.up.y,
        cam.up.z,
    );
    let _ = writeln!(
        s,
        "Camera \"perspective\" \"float fov\" [{:.6}]",
        cam.fov_v_rad.to_degrees()
    );
    let _ = writeln!(s, "WorldBegin");
    for (id, mat) in &job.materials {
        let _ = writeln!(
            s,
            "MakeNamedMaterial \"{}\" \"string type\" [\"matte\"] \"color Kd\" [{:.4} {:.4} {:.4}]",
            id, mat.diffuse_color[0], mat.diffuse_color[1], mat.diffuse_color[2]
        );
    }
    for m in &job.meshes {
        let _ = writeln!(s, "AttributeBegin");
        let _ = writeln!(s, "NamedMaterial \"{}\"", m.material_id);
        emit_lux_mesh(&mut s, m);
        let _ = writeln!(s, "AttributeEnd");
    }
    let _ = writeln!(s, "WorldEnd");
    s
}

fn emit_lux_mesh(s: &mut String, m: &SceneMesh) {
    let _ = writeln!(s, "# mesh: {}", m.name);
    let _ = write!(s, "Shape \"trianglemesh\" \"point P\" [");
    for v in &m.mesh.nodes {
        let _ = write!(s, " {:.6} {:.6} {:.6}", v.x, v.y, v.z);
    }
    let _ = writeln!(s, " ]");
    let _ = write!(s, "  \"integer triindices\" [");
    for block in &m.mesh.element_blocks {
        if block.element_type != valenx_mesh::element::ElementType::Tri3 {
            continue;
        }
        for chunk in block.connectivity.chunks_exact(3) {
            let _ = write!(s, " {} {} {}", chunk[0], chunk[1], chunk[2]);
        }
    }
    let _ = writeln!(s, " ]");
}

/// Emit a Cycles XML scene. v1 emits the bare-minimum shader-graph
/// (`diffuse_bsdf`) per material plus a `mesh` node per object.
fn emit_cycles(job: &RenderJob) -> String {
    let mut s = String::new();
    let _ = writeln!(
        s,
        "<!-- Cycles XML scene — emitted by valenx-render-bridge v1 -->"
    );
    let _ = writeln!(s, "<cycles>");
    let cam = &job.camera;
    let _ = writeln!(
        s,
        "  <camera width=\"{}\" height=\"{}\" fov=\"{:.6}\" />",
        cam.image_width, cam.image_height, cam.fov_v_rad
    );
    let _ = writeln!(
        s,
        "  <transform matrix=\"1 0 0 {:.6} 0 1 0 {:.6} 0 0 1 {:.6} 0 0 0 1\">",
        cam.position.x, cam.position.y, cam.position.z
    );
    let _ = writeln!(s, "    <camera />");
    let _ = writeln!(s, "  </transform>");
    for (id, mat) in &job.materials {
        let _ = writeln!(s, "  <shader name=\"{id}\">");
        let _ = writeln!(
            s,
            "    <diffuse_bsdf name=\"bsdf\" color=\"{:.4} {:.4} {:.4}\" />",
            mat.diffuse_color[0], mat.diffuse_color[1], mat.diffuse_color[2]
        );
        let _ = writeln!(
            s,
            "    <connect from=\"bsdf bsdf\" to=\"output surface\" />"
        );
        let _ = writeln!(s, "  </shader>");
    }
    for m in &job.meshes {
        let _ = writeln!(s, "  <state shader=\"{}\">", m.material_id);
        let mut verts = String::new();
        for v in &m.mesh.nodes {
            let _ = write!(verts, "{:.6} {:.6} {:.6}  ", v.x, v.y, v.z);
        }
        let mut tris = String::new();
        for block in &m.mesh.element_blocks {
            if block.element_type != valenx_mesh::element::ElementType::Tri3 {
                continue;
            }
            for chunk in block.connectivity.chunks_exact(3) {
                let _ = write!(tris, "{} {} {}  ", chunk[0], chunk[1], chunk[2]);
            }
        }
        let _ = writeln!(
            s,
            "    <mesh P=\"{}\" verts=\"{}\" nverts=\"3\" />",
            verts.trim(),
            tris.trim()
        );
        let _ = writeln!(s, "  </state>");
    }
    let _ = writeln!(s, "</cycles>");
    s
}

/// Emit a POV-Ray SDL `.pov` scene. v1 emits the camera + a point
/// light + per-mesh `mesh2 { vertex_vectors face_indices texture }`
/// blocks with a plain `pigment { rgb }` texture per material.
fn emit_povray(job: &RenderJob) -> String {
    let mut s = String::new();
    let _ = writeln!(s, "// POV-Ray scene — emitted by valenx-render-bridge v1");
    let _ = writeln!(s, "// Title: {}", job.title);
    let _ = writeln!(s, "global_settings {{ assumed_gamma 1.0 }}");
    let cam = &job.camera;
    let _ = writeln!(s, "camera {{");
    let _ = writeln!(
        s,
        "  location <{:.6}, {:.6}, {:.6}>",
        cam.position.x, cam.position.y, cam.position.z
    );
    let _ = writeln!(
        s,
        "  look_at <{:.6}, {:.6}, {:.6}>",
        cam.target.x, cam.target.y, cam.target.z
    );
    let _ = writeln!(
        s,
        "  up <{:.6}, {:.6}, {:.6}>",
        cam.up.x, cam.up.y, cam.up.z
    );
    let _ = writeln!(s, "  angle {:.6}", cam.fov_v_rad.to_degrees());
    let _ = writeln!(s, "}}");
    // Default white point light if user didn't supply any.
    if job.lights.is_empty() {
        let _ = writeln!(
            s,
            "light_source {{ <10, 10, 10> color rgb <1.0, 1.0, 1.0> }}"
        );
    } else {
        for light in &job.lights {
            if let crate::light::Light::Point {
                position, color, ..
            } = light
            {
                let _ = writeln!(
                    s,
                    "light_source {{ <{:.6}, {:.6}, {:.6}> color rgb <{:.4}, {:.4}, {:.4}> }}",
                    position.x, position.y, position.z, color[0], color[1], color[2]
                );
            }
        }
    }
    for m in &job.meshes {
        let mat = job.materials.get(&m.material_id).expect("validated");
        let _ = writeln!(s, "mesh2 {{");
        let _ = write!(s, "  vertex_vectors {{ {}", m.mesh.nodes.len());
        for v in &m.mesh.nodes {
            let _ = write!(s, ", <{:.6}, {:.6}, {:.6}>", v.x, v.y, v.z);
        }
        let _ = writeln!(s, " }}");
        // Triangle indices
        let mut tris: Vec<[usize; 3]> = Vec::new();
        for block in &m.mesh.element_blocks {
            if block.element_type != valenx_mesh::element::ElementType::Tri3 {
                continue;
            }
            for chunk in block.connectivity.chunks_exact(3) {
                tris.push([chunk[0] as usize, chunk[1] as usize, chunk[2] as usize]);
            }
        }
        let _ = write!(s, "  face_indices {{ {}", tris.len());
        for t in &tris {
            let _ = write!(s, ", <{}, {}, {}>", t[0], t[1], t[2]);
        }
        let _ = writeln!(s, " }}");
        let _ = writeln!(
            s,
            "  texture {{ pigment {{ rgb <{:.4}, {:.4}, {:.4}> }} }}",
            mat.diffuse_color[0], mat.diffuse_color[1], mat.diffuse_color[2]
        );
        let _ = writeln!(s, "}}");
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::camera::Camera;
    use crate::engine::RenderEngine;
    use crate::material::Material;
    use crate::scene::RenderJob;
    use nalgebra::Vector3;

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

    fn one_mesh_job(engine: RenderEngine) -> RenderJob {
        let mut job = RenderJob::new();
        job.engine = engine;
        job.camera = Camera::default();
        let id = job.add_material(Material::matte("paint", [0.5, 0.5, 0.5]));
        job.add_mesh("triangle", unit_tri_mesh(), id);
        job
    }

    #[test]
    fn empty_scene_errors() {
        let job = RenderJob::new();
        assert!(scene_file(&job).is_err());
    }

    #[test]
    fn lux_emits_with_film() {
        let s = scene_file(&one_mesh_job(RenderEngine::LuxRender)).unwrap();
        assert!(s.contains("Film \"fleximage\""));
        assert!(s.contains("WorldBegin"));
        assert!(s.contains("MakeNamedMaterial"));
    }

    #[test]
    fn cycles_emits_xml() {
        let s = scene_file(&one_mesh_job(RenderEngine::Cycles)).unwrap();
        assert!(s.contains("<cycles>"));
        assert!(s.contains("diffuse_bsdf"));
    }

    #[test]
    fn povray_emits_mesh2() {
        let s = scene_file(&one_mesh_job(RenderEngine::PovRay)).unwrap();
        assert!(s.contains("mesh2"));
        assert!(s.contains("vertex_vectors"));
    }

    #[test]
    fn native_engine_returns_not_implemented() {
        assert!(matches!(
            scene_file(&one_mesh_job(RenderEngine::Native)).unwrap_err(),
            RenderError::EngineNotImplemented {
                engine: "Native",
                ..
            }
        ));
    }

    #[test]
    fn missing_material_errors() {
        let mut job = RenderJob::new();
        job.engine = RenderEngine::LuxRender;
        // Add a mesh that references a missing material id.
        job.add_mesh("tri", unit_tri_mesh(), "ghost".to_string());
        assert!(matches!(
            scene_file(&job).unwrap_err(),
            RenderError::BadParameter { .. }
        ));
    }
}
