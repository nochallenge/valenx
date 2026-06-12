//! # valenx-render-bridge
//!
//! Render bridge — describe a scene (geometry + camera + lights +
//! materials) and emit it as the scene-format file expected by an
//! external renderer (LuxRender / Cycles / POV-Ray). The actual
//! subprocess launch is the caller's responsibility — this crate
//! stays pure-data so it compiles on every workspace target without
//! pulling in a process / tokio dep.
//!
//! The FreeCAD `Render` community workbench equivalent.
//!
//! Phase 30 of the FreeCAD-parity roadmap.
//!
//! # Pipeline
//!
//! 1. Build a [`RenderJob`] by combining a [`Camera`], one or more
//!    [`Light`]s, a material library (`HashMap<MaterialId, Material>`),
//!    and the meshes to render.
//! 2. Pick a [`RenderEngine`].
//! 3. Call [`emit::scene_file`] to produce the text payload the
//!    renderer expects.
//! 4. Write the payload to disk + launch the external renderer
//!    (subprocess launch is owned by the desktop shell).
//!
//! The [`RenderEngine::Native`] variant short-circuits the
//! scene-file step — it documents that the desktop shell should
//! capture the active wgpu viewport as a PNG instead of dispatching
//! to an external renderer.
//!
//! # Real subprocess adapters (Phase 30.5)
//!
//! [`subprocess`] closes the loop for the external path tracers:
//! [`subprocess::run_cycles`] dispatches the Cycles XML scene to the
//! standalone `cycles` renderer, and [`subprocess::run_luxcore`]
//! serialises the job to LuxCoreRender's modern `.scn` + `.cfg` SDL
//! pair ([`subprocess::emit_luxcore_scn`] /
//! [`subprocess::emit_luxcore_cfg`]) and invokes `luxcoreconsole`. A
//! renderer that is not on `PATH` yields a clear
//! [`error::RenderError::ToolNotAvailable`] — never a silent no-op.
//!
//! # Environment lighting (Phase 30 IBL)
//!
//! [`environment`] adds HDR image-based lighting: [`EnvironmentMap`]
//! decodes a Radiance `.hdr` file and computes the diffuse-irradiance
//! convolution; [`scene::EnvironmentRef`] is the serialisable scene
//! reference (file path + intensity + yaw) carried on [`RenderJob`].
//!
//! # Real-time PBR shading (Phases 30.6–30.7)
//!
//! [`pbr`] adds the **Cook-Torrance microfacet BRDF** — the GGX
//! normal distribution, the Smith geometry term, Fresnel-Schlick —
//! and a forward-shading evaluator ([`pbr::shade_surface`]) that
//! lights a surface point from analytic lights plus the IBL
//! environment term. This is the real-time PBR path (the closed-form
//! BRDF a GPU fragment shader runs per pixel), *not* a path tracer.
//!
//! The **split-sum specular IBL** (Phase 30.7) completes the pipeline:
//! [`environment::EnvironmentMap::prefilter_specular`] convolves the
//! environment into a roughness-indexed mip chain
//! ([`PrefilteredEnvironment`]) by GGX importance sampling, and
//! [`pbr::compute_brdf_lut`] precomputes the environment-BRDF
//! scale/bias table ([`BrdfLut`]); [`pbr::specular_ibl`] reconstructs
//! the roughness-aware specular ambient term from the two factors.
//!
//! The CPU PBR + IBL library is complete and test-covered. The
//! **WGSL forward shader** that ports the BRDF to the GPU ships in
//! [`wgsl_pbr`] (see below) — GPU-unverified, as that module documents.
//!
//! # Real-time GI — irradiance-volume light probes
//!
//! [`irradiance_volume`] adds **baked real-time global illumination**:
//! a 3-D grid of spherical-harmonic light probes
//! ([`IrradianceVolume`]). Each probe's incident radiance is gathered
//! from a caller-supplied scene-radiance closure (a path-tracer
//! ray-gather) and stored as `L1` / `L2` SH; a
//! [`IrradianceVolume::sample_irradiance`] lookup trilinearly blends
//! the surrounding probes for the indirect bounced light at any point.
//!
//! # WGSL PBR forward shader
//!
//! [`wgsl_pbr`] is the **GPU port** of the [`pbr`] BRDF —
//! [`wgsl_pbr::PBR_FORWARD_WGSL`] is the complete Cook-Torrance forward
//! shader (the identical GGX/Smith/Fresnel maths) plus `#[repr(C)]`
//! uniform-block layouts a `wgpu` host casts straight into device
//! buffers. **The WGSL is GPU-unverified** — it compiles and `cargo
//! check`s and its BRDF is cross-checked term-by-term against the CPU
//! reference, but it has not been run on hardware (see the module's
//! honest-scope note). The CPU library remains the verified path.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod camera;
pub mod emit;
pub mod engine;
pub mod environment;
pub mod error;
pub mod irradiance_volume;
pub mod light;
pub mod material;
pub mod pbr;
pub mod persist;
pub mod scene;
pub mod subprocess;
pub mod wgsl_pbr;

pub use camera::Camera;
pub use emit::scene_file;
pub use engine::RenderEngine;
pub use environment::{EnvironmentMap, IrradianceMap, PrefilteredEnvironment, PrefilteredLevel};
pub use error::{ErrorCategory, RenderError};
pub use irradiance_volume::{fibonacci_sphere, sh_basis, IrradianceVolume, LightProbe, ShOrder};
pub use light::Light;
pub use material::{Material, MaterialId};
pub use pbr::{
    ambient_ibl, brdf_direct, compute_brdf_lut, incident_light, shade_surface, specular_ibl,
    BrdfLut, IncidentLight, ShadedColor, SurfacePoint,
};
pub use persist::RenderFile;
pub use scene::{EnvironmentRef, RenderJob, SceneMesh};
pub use subprocess::{
    cycles_command, emit_luxcore_cfg, emit_luxcore_scn, luxcore_command, renderer_available,
    run_cycles, run_luxcore, RenderArtifacts,
};
pub use wgsl_pbr::{
    reference_brdf_rgb, PbrFrameUniform, PbrLightUniform, PbrMaterialUniform, ShL2Uniform,
    MAX_LIGHTS, PBR_FORWARD_WGSL,
};
