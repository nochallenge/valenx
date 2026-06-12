//! Scene description — meshes + camera + lights + materials +
//! engine + output target.

use std::collections::HashMap;
use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use crate::camera::Camera;
use crate::engine::RenderEngine;
use crate::light::Light;
use crate::material::{Material, MaterialId};

/// One mesh in the scene, tagged with the material it should render
/// with.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SceneMesh {
    /// Display name.
    pub name: String,
    /// The triangulated mesh (Phase 7 tessellation output).
    pub mesh: valenx_mesh::Mesh,
    /// Material id — must exist in the parent
    /// [`RenderJob::materials`].
    pub material_id: MaterialId,
}

/// A scene reference to an HDR environment map for image-based
/// lighting.
///
/// This is the **serialisable** form stored on a [`RenderJob`]: just
/// the file path and the two exposure / orientation knobs. The actual
/// pixel data is loaded at render time into a
/// [`crate::environment::EnvironmentMap`] (which is intentionally not
/// serialised — a multi-megabyte float buffer does not belong in a
/// RON scene file).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EnvironmentRef {
    /// Path to the Radiance `.hdr` equirectangular environment map.
    pub hdr_path: PathBuf,
    /// Exposure multiplier applied to every environment sample.
    pub intensity: f32,
    /// Yaw rotation of the environment about world `+Y`, in radians.
    pub yaw: f32,
}

impl EnvironmentRef {
    /// A reference to `hdr_path` with neutral exposure (intensity 1,
    /// no yaw).
    pub fn new(hdr_path: impl Into<PathBuf>) -> Self {
        Self {
            hdr_path: hdr_path.into(),
            intensity: 1.0,
            yaw: 0.0,
        }
    }

    /// Load the referenced `.hdr` file into an
    /// [`crate::environment::EnvironmentMap`], applying this ref's
    /// `intensity` and `yaw`.
    ///
    /// # Errors
    ///
    /// [`crate::error::RenderError::Io`] if the file cannot be read
    /// (including size-cap exceeded — see
    /// [`valenx_core::io_caps::MAX_HDR_FILE_BYTES`]), or
    /// [`crate::error::RenderError::BadParameter`] if it is not a
    /// valid Radiance HDR.
    ///
    /// Round-20 M2: pre-fix this was a bare `std::fs::read(&self.hdr_path)`
    /// — a serialised `.ron` scene with a hostile or stale `hdr_path`
    /// would slurp a multi-GB file into memory before the HDR parser
    /// saw the magic bytes. The cap is 256 MiB (generous for an
    /// 8K equirectangular RGBE map; refuses the
    /// `cat /dev/zero > big.hdr` DoS).
    pub fn load(&self) -> Result<crate::environment::EnvironmentMap, crate::error::RenderError> {
        let bytes = valenx_core::io_caps::read_capped_to_bytes(
            &self.hdr_path,
            valenx_core::io_caps::MAX_HDR_FILE_BYTES,
        )?;
        let mut map = crate::environment::EnvironmentMap::from_radiance_hdr(&bytes)?;
        map.intensity = self.intensity;
        map.yaw = self.yaw;
        Ok(map)
    }
}

/// A full render job — everything the emitter needs to write a scene
/// file.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RenderJob {
    /// Optional display title.
    pub title: String,
    /// Meshes to render.
    pub meshes: Vec<SceneMesh>,
    /// Material library keyed by id. Every `SceneMesh::material_id`
    /// must resolve here.
    pub materials: HashMap<MaterialId, Material>,
    /// Camera.
    pub camera: Camera,
    /// Lights.
    pub lights: Vec<Light>,
    /// Optional HDR environment map for image-based lighting. `None`
    /// leaves the scene lit only by the explicit [`Light`] list.
    #[serde(default)]
    pub environment: Option<EnvironmentRef>,
    /// Output image path (or scene-file path for engines that need a
    /// separate dispatch step).
    pub output_path: PathBuf,
    /// Engine choice.
    pub engine: RenderEngine,
}

impl Default for RenderJob {
    fn default() -> Self {
        Self {
            title: "Untitled".into(),
            meshes: Vec::new(),
            materials: HashMap::new(),
            camera: Camera::default(),
            lights: Vec::new(),
            environment: None,
            output_path: PathBuf::from("render.png"),
            engine: RenderEngine::Native,
        }
    }
}

impl RenderJob {
    /// Empty job with a sensible default camera and the Native
    /// engine.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a mesh to the scene.
    pub fn add_mesh(
        &mut self,
        name: impl Into<String>,
        mesh: valenx_mesh::Mesh,
        material_id: impl Into<MaterialId>,
    ) {
        self.meshes.push(SceneMesh {
            name: name.into(),
            mesh,
            material_id: material_id.into(),
        });
    }

    /// Insert / replace a material in the library.
    pub fn add_material(&mut self, material: Material) -> MaterialId {
        let id = material.name.clone();
        self.materials.insert(id.clone(), material);
        id
    }

    /// Append a light.
    pub fn add_light(&mut self, light: Light) {
        self.lights.push(light);
    }

    /// Set (or clear, with `None`) the HDR environment map used for
    /// image-based lighting.
    pub fn set_environment(&mut self, environment: Option<EnvironmentRef>) {
        self.environment = environment;
    }

    /// Verify every mesh's `material_id` resolves in `materials`.
    pub fn validate(&self) -> Result<(), crate::error::RenderError> {
        if self.meshes.is_empty() {
            return Err(crate::error::RenderError::EmptyScene);
        }
        for m in &self.meshes {
            if !self.materials.contains_key(&m.material_id) {
                return Err(crate::error::RenderError::BadParameter {
                    name: "material_id",
                    reason: format!(
                        "mesh `{}` references missing material `{}`",
                        m.name, m.material_id
                    ),
                });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn environment_ref_defaults_to_neutral_exposure() {
        let env = EnvironmentRef::new("studio.hdr");
        assert_eq!(env.intensity, 1.0);
        assert_eq!(env.yaw, 0.0);
        assert_eq!(env.hdr_path, PathBuf::from("studio.hdr"));
    }

    #[test]
    fn render_job_environment_round_trips_through_ron() {
        let mut job = RenderJob::new();
        job.set_environment(Some(EnvironmentRef {
            hdr_path: PathBuf::from("sky.hdr"),
            intensity: 2.5,
            yaw: 1.0,
        }));
        let text = ron::to_string(&job).unwrap();
        let back: RenderJob = ron::from_str(&text).unwrap();
        assert_eq!(back.environment, job.environment);
    }

    #[test]
    fn render_job_without_environment_still_deserialises() {
        // The `#[serde(default)]` on `environment` means an older
        // scene file with no environment field still loads.
        let mut job = RenderJob::new();
        job.set_environment(None);
        let text = ron::to_string(&job).unwrap();
        let back: RenderJob = ron::from_str(&text).unwrap();
        assert!(back.environment.is_none());
    }

    /// Round-20 M2 RED→GREEN: a `.hdr` file larger than
    /// `MAX_HDR_FILE_BYTES` (256 MiB) must be rejected as an IO
    /// error BEFORE `from_radiance_hdr` allocates a 500 MiB byte
    /// vec. Pre-fix the loader did a bare `std::fs::read(&path)` —
    /// a hostile multi-GB file would have slurped into memory before
    /// the HDR magic-byte check happened.
    #[test]
    fn environment_ref_load_rejects_oversize_hdr_file() {
        use std::io::Write;
        let tmp = std::env::temp_dir().join(format!(
            "valenx_r20m2_oversize_{}.hdr",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        // 257 MiB — just past the 256 MiB cap. Use set_len + a single
        // byte write so the test doesn't write 500 MiB on every CI.
        let mut f = std::fs::File::create(&tmp).unwrap();
        f.set_len(257 * 1024 * 1024).unwrap();
        f.write_all(b"x").unwrap();
        drop(f);
        let env = EnvironmentRef::new(&tmp);
        let err = env
            .load()
            .expect_err("round-20 M2: 257 MiB HDR must be rejected as IO error");
        assert_eq!(
            err.code(),
            "render.io",
            "an oversize HDR file is an IO error (size-cap exceeded)"
        );
        let _ = std::fs::remove_file(&tmp);
    }
}
