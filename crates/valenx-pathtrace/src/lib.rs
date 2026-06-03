//! # valenx-pathtrace
//!
//! A CPU **Monte-Carlo path tracer** with real global illumination —
//! the offline-rendering counterpart to the real-time Cook-Torrance
//! forward shader in `valenx-render-bridge::pbr`.
//!
//! ## What this is
//!
//! A genuine, unbiased path tracer. It solves the rendering equation by
//! tracing light-transport paths and Monte-Carlo-integrating the
//! result, so it captures the effects a forward shader cannot —
//! multi-bounce indirect light, colour bleeding, soft shadows from area
//! lights, ambient occlusion — all as emergent consequences of the
//! integral, not special-cased.
//!
//! The pipeline:
//!
//! 1. **Acceleration** — a [`bvh::Bvh`] (binned surface-area-heuristic
//!    bounding-volume hierarchy) over the scene triangles, so a ray
//!    cast is `O(log triangles)`.
//! 2. **Intersection** — the Möller-Trumbore ray-triangle test
//!    ([`geometry::Triangle::intersect`]).
//! 3. **Sampling** — a seedable PCG32 RNG ([`sampling::Rng`]) feeding
//!    cosine-weighted hemisphere importance sampling
//!    ([`sampling::cosine_hemisphere`]) and GGX half-vector sampling.
//! 4. **Shading** — the Cook-Torrance microfacet BRDF + a Lambert
//!    diffuse lobe, reusing the GGX / Smith / Fresnel functions from
//!    [`valenx_render_bridge::pbr`].
//! 5. **Integration** — [`tracer::render`] runs the path integral per
//!    pixel with **next-event estimation** (direct connections to
//!    sampled emitters), **Russian-roulette** path termination, and
//!    **HDR environment lighting** (an escaped ray samples the
//!    [`valenx_render_bridge::environment::EnvironmentMap`] — the same
//!    Radiance `.hdr` decoder the IBL path uses).
//! 6. **Output** — an [`framebuffer::HdrFramebuffer`] of raw radiance,
//!    tone-mapped (ACES filmic + sRGB) to a displayable
//!    [`framebuffer::LdrImage`].
//!
//! ## Honest scope — a real v1, not Cycles
//!
//! Every algorithm here is the genuine article and the renderer
//! converges, without bias, to the correct image — the white-furnace
//! test and the analytic direct-lighting test in [`tracer`] verify
//! that. It is deliberately a **v1**: it does *not* match a
//! decade-matured production renderer's *speed* or *feature breadth*.
//!
//! ## Graduated extensions
//!
//! Several of the originally-deferred follow-ups have since shipped as
//! real, test-covered modules — Cycles-class building blocks:
//!
//! - [`mis`] — **multiple importance sampling**. [`mis::render_mis`]
//!   combines next-event (light) sampling and BSDF sampling under the
//!   power heuristic — unbiased (it converges to the same radiance as
//!   the NEE [`tracer::render`] path) with markedly lower variance on
//!   a glossy surface lit by an area light.
//! - [`dielectric`] — a **transmission / refraction BSDF**: the
//!   dielectric Fresnel equations, Snell refraction, total internal
//!   reflection, and a GGX rough-dielectric ("frosted glass") variant.
//! - [`denoise`] — the **edge-avoiding à-trous wavelet denoiser**
//!   (Dammertz et al. 2010), a classical, non-ML denoiser guided by
//!   albedo / normal / depth feature buffers.
//! - [`volume`] — **volumetric rendering**: ray-marched
//!   emission-absorption + single scattering + the Henyey-Greenstein
//!   phase function ([`volume::render_volume`]).
//! - [`light_tree`] — a **hierarchical light importance tree**
//!   (Conty-Estevez & Kulla 2018 — Cycles, PBRT v4): emitter sampling
//!   that descends a power × geometric-importance binary hierarchy in
//!   place of a uniform light pick. The variance falls dramatically
//!   on many-light scenes; the [`tracer`] NEE path and the [`mis`]
//!   light-sampling path now both go through it.
//! - [`bdpt`] — **bidirectional path tracing** (Veach 1997). Traces a
//!   camera and a light subpath and connects every prefix-pair of
//!   vertices, weighted across all strategies under the MIS power
//!   heuristic. Captures specular-diffuse-specular (caustic) paths
//!   that unidirectional + NEE cannot deliver in finite samples.
//!   Opt-in render mode — [`bdpt::render_bdpt`] is a peer to
//!   [`tracer::render`] / [`mis::render_mis`].
//! - [`sss`] — **subsurface scattering** via a random-walk BSSRDF
//!   (skin, marble, wax). A subsurface ray enters the surface,
//!   diffuses through the medium with Henyey-Greenstein phase
//!   sampling and per-channel Beer-Lambert extinction, and exits at a
//!   nearby surface point with the correct multiple-scattered exitance.
//!
//! What still deliberately waits for documented follow-ups:
//!
//! - **No spectral rendering, no Metropolis-style mutation chains,
//!   no adaptive sampling, no photon mapping.**
//! - **Single-threaded.** The render loop is structured so it *could*
//!   parallelise (each pixel's RNG is independently seeded and the
//!   scene is read-only), but it does not spawn threads — that keeps
//!   the crate dependency-free.
//!
//! None of those omissions affects the correctness of what ships; each
//! is an additive, well-understood extension.
//!
//! ## Example
//!
//! ```
//! use valenx_pathtrace::{
//!     PtCamera, PtMaterial, RenderParams, SceneBuilder, Vec3, render,
//! };
//! use valenx_pathtrace::math::vec3;
//!
//! // A camera looking at the origin.
//! let camera = PtCamera::look_at(
//!     vec3(0.0, 0.0, 4.0),
//!     Vec3::ZERO,
//!     vec3(0.0, 1.0, 0.0),
//!     50f32.to_radians(),
//!     32,
//!     32,
//! );
//! let mut builder = SceneBuilder::new(camera);
//! let white = builder.add_material(PtMaterial::diffuse([0.8, 0.8, 0.8]));
//! builder.add_quad(
//!     vec3(-5.0, -5.0, 0.0),
//!     vec3(5.0, -5.0, 0.0),
//!     vec3(5.0, 5.0, 0.0),
//!     vec3(-5.0, 5.0, 0.0),
//!     white,
//! );
//! let scene = builder.build();
//!
//! let params = RenderParams { samples_per_pixel: 8, ..RenderParams::default() };
//! let framebuffer = render(&scene, &params).expect("render small framebuffer");
//! let image = framebuffer.to_ldr(params.exposure);
//! assert_eq!(image.width, 32);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bdpt;
pub mod bvh;
pub mod denoise;
pub mod dielectric;
pub mod error;
pub mod framebuffer;
pub mod geometry;
pub mod light_tree;
pub mod math;
pub mod mis;
pub mod sampling;
pub mod scene;
pub mod sss;
pub mod tracer;
pub mod volume;

pub use bdpt::{render_bdpt, BdptParams};
pub use bvh::Bvh;
pub use denoise::{denoise_atrous, AtrousParams, GuideBuffers};
pub use dielectric::{
    fresnel_dielectric, refract, sample_rough_dielectric, sample_smooth_dielectric,
    DielectricEvent, DielectricSample,
};
pub use error::PathTraceError;
pub use framebuffer::{tonemap_pixel, HdrFramebuffer, LdrImage};
pub use geometry::{Aabb, Hit, Ray, Triangle};
pub use light_tree::{LightSample, LightTree};
pub use math::{ortho_basis, vec3, Vec3};
pub use mis::{balance_heuristic, bsdf_pdf, power_heuristic, render_mis};
pub use sampling::{cosine_hemisphere, cosine_hemisphere_pdf, Rng};
pub use scene::{PtCamera, PtMaterial, Scene, SceneBuilder, Subsurface};
pub use sss::{random_walk_sss, RandomWalkResult};
pub use tracer::{render, trace_single_ray, RenderParams};
pub use volume::{
    henyey_greenstein, march_ray, render_volume, transmittance, DensityGrid, Medium,
    VolumeBox, VolumeLight, VolumeParams, VolumeResult,
};
