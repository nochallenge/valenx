//! wgpu **Cook-Torrance PBR forward render pass** — the GPU-side
//! companion to `valenx-render-bridge`'s shading library.
//!
//! # What this is
//!
//! `valenx-render-bridge::pbr` is the CPU-side physically-based
//! shading library, and `valenx_render_bridge::wgsl_pbr` ports its
//! Cook-Torrance BRDF to a WGSL fragment shader
//! ([`valenx_render_bridge::PBR_FORWARD_WGSL`]). This module is the
//! **`wgpu` render-pass plumbing** that drives that shader: it builds
//! the pipeline, the bind-group layout, and the uniform / vertex
//! buffers, and encodes a forward pass that shades geometry against
//!
//! - analytic lights (point / directional) — the Cook-Torrance BRDF;
//! - the IBL ambient term — a diffuse + Fresnel-specular environment;
//! - the **irradiance-volume GI** term — an `L2` spherical-harmonic
//!   probe ([`valenx_render_bridge::irradiance_volume`]);
//! - the material's emissive colour.
//!
//! # Self-contained — does NOT touch the live viewport loop
//!
//! [`PbrForwardPass`] is a **standalone module**: it owns its own
//! pipeline, its own offscreen colour + depth targets, and its own
//! buffers. It is deliberately *not* wired into the app's existing
//! `wgpu_renderer::WgpuRenderer` viewport loop — that flat-shaded loop
//! keeps running untouched. A caller that wants the PBR pass
//! constructs a [`PbrForwardPass`] from the same
//! `egui_wgpu::RenderState`, renders into its offscreen target, and
//! displays the returned [`egui::TextureId`]; the two render paths are
//! independent and cannot destabilise each other.
//!
//! # HONEST REQUIREMENT — this code is GPU-unverified
//!
//! **This module is GPU shader plumbing. It compiles and `cargo
//! check`s cleanly, the WGSL it loads is cross-checked term-by-term
//! against the CPU BRDF reference (see
//! `valenx_render_bridge::wgsl_pbr`'s tests), and the uniform layouts
//! are size-asserted — but the render pass has NOT been executed on a
//! GPU.** Producing a correct on-screen image requires running the
//! app on real graphics hardware, which the automated test suite
//! does not do. Treat this pass as *written-correct-against-the-CPU-
//! reference* but *not pixel-validated*. The verified shading path is
//! the CPU library in `valenx-render-bridge`; this is its faithful
//! GPU wiring, pending a hardware run.

use std::sync::Arc;

use bytemuck::{Pod, Zeroable};
use eframe::egui;
use eframe::egui_wgpu;
use eframe::wgpu::{self, util::DeviceExt};

use valenx_render_bridge::wgsl_pbr::{
    PbrFrameUniform, PbrLightUniform, PbrMaterialUniform, ShL2Uniform, MAX_LIGHTS, PBR_FORWARD_WGSL,
};
use valenx_render_bridge::{IrradianceVolume, Material};

/// Vertex layout the PBR pipeline consumes — position + normal, the
/// same `[f32; 3] × 2` layout the existing viewport pipeline
/// (`wgpu_renderer::Vertex`) uses, so geometry buffers are
/// interchangeable.
#[repr(C)]
#[derive(Copy, Clone, Debug, Default, Pod, Zeroable)]
pub struct PbrVertex {
    /// World-space vertex position.
    pub position: [f32; 3],
    /// World-space vertex normal.
    pub normal: [f32; 3],
}

const PBR_VERTEX_ATTRS: [wgpu::VertexAttribute; 2] = wgpu::vertex_attr_array![
    0 => Float32x3, // position
    1 => Float32x3, // normal
];

fn pbr_vertex_layout() -> wgpu::VertexBufferLayout<'static> {
    wgpu::VertexBufferLayout {
        array_stride: std::mem::size_of::<PbrVertex>() as wgpu::BufferAddress,
        step_mode: wgpu::VertexStepMode::Vertex,
        attributes: &PBR_VERTEX_ATTRS,
    }
}

/// The depth format of the PBR pass's own depth target.
const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;

/// The fixed-size GPU light array — `MAX_LIGHTS` slots, padded with
/// inactive lights beyond the active count. `#[repr(C)]` so it casts
/// straight into a uniform buffer matching the WGSL
/// `array<LightUniform, MAX_LIGHTS>`.
#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct LightArray {
    lights: [PbrLightUniform; MAX_LIGHTS],
}

impl LightArray {
    /// Pack a slice of lights into the fixed array, padding the unused
    /// tail with [`PbrLightUniform::inactive`].
    fn from_slice(lights: &[PbrLightUniform]) -> LightArray {
        let mut arr = [PbrLightUniform::inactive(); MAX_LIGHTS];
        for (slot, light) in arr.iter_mut().zip(lights.iter()) {
            *slot = *light;
        }
        LightArray { lights: arr }
    }
}

/// One frame's worth of shading inputs for [`PbrForwardPass::render`].
///
/// This is the CPU description the pass uploads to its uniform buffers
/// each frame — keeping it a plain struct lets a caller assemble it
/// without touching `wgpu` types.
pub struct PbrFrame<'a> {
    /// Model-view-projection matrix, column-major.
    pub view_proj: [[f32; 4]; 4],
    /// World-space camera position.
    pub camera_pos: [f32; 3],
    /// Ambient-IBL environment radiance (linear RGB).
    pub env_color: [f32; 3],
    /// Analytic lights — at most [`MAX_LIGHTS`] are used; extras are
    /// ignored.
    pub lights: &'a [PbrLightUniform],
    /// The surface material.
    pub material: &'a Material,
    /// Ambient-occlusion multiplier for the ambient term (`1.0` = no
    /// occlusion).
    pub ambient_occlusion: f32,
    /// The `L2` spherical-harmonic GI probe for the indirect bounced
    /// light — typically one probe sampled from an
    /// [`IrradianceVolume`]. Use [`ShL2Uniform::zero`] for no GI.
    pub gi_probe: ShL2Uniform,
    /// Triangle geometry — position + normal per vertex, three
    /// vertices per triangle.
    pub vertices: &'a [PbrVertex],
}

/// A self-contained Cook-Torrance PBR forward render pass.
///
/// Construct one with [`PbrForwardPass::new`] from the app's
/// `egui_wgpu::RenderState`; drive it per frame with
/// [`PbrForwardPass::render`]. It owns its pipeline and offscreen
/// targets and never touches the viewport's flat-shaded loop.
pub struct PbrForwardPass {
    device: Arc<wgpu::Device>,
    queue: Arc<wgpu::Queue>,
    target_format: wgpu::TextureFormat,

    pipeline: wgpu::RenderPipeline,
    // Kept alive for the pipeline; the bind group references it.
    #[allow(dead_code)]
    bind_group_layout: wgpu::BindGroupLayout,
    frame_buffer: wgpu::Buffer,
    material_buffer: wgpu::Buffer,
    light_buffer: wgpu::Buffer,
    gi_buffer: wgpu::Buffer,
    bind_group: wgpu::BindGroup,

    vertex_buffer: wgpu::Buffer,
    vertex_capacity: u64,
    vertex_count: u32,

    offscreen: Option<PbrOffscreen>,
}

/// The PBR pass's offscreen colour + depth render targets.
struct PbrOffscreen {
    size: [u32; 2],
    #[allow(dead_code)]
    color_texture: wgpu::Texture,
    color_view: wgpu::TextureView,
    #[allow(dead_code)]
    depth_texture: wgpu::Texture,
    depth_view: wgpu::TextureView,
    egui_id: egui::TextureId,
}

impl PbrForwardPass {
    /// Build the PBR pipeline and its initial (empty) buffers from the
    /// app's `egui_wgpu::RenderState`.
    ///
    /// The offscreen colour / depth targets are created lazily on the
    /// first [`PbrForwardPass::render`] — the viewport pixel size is
    /// not known at construction.
    pub fn new(render_state: &egui_wgpu::RenderState) -> Self {
        let device = render_state.device.clone();
        let queue = render_state.queue.clone();
        let target_format = render_state.target_format;

        // The WGSL Cook-Torrance forward shader, straight from the
        // render-bridge crate.
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("valenx.pbr_forward.shader"),
            source: wgpu::ShaderSource::Wgsl(PBR_FORWARD_WGSL.into()),
        });

        // Bind group 0: frame / material / light-array / GI-probe
        // uniforms — bindings 0..3, matching the WGSL declarations.
        let uniform_entry = |binding: u32| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::VERTEX_FRAGMENT,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Uniform,
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("valenx.pbr_forward.bind_group_layout"),
            entries: &[
                uniform_entry(0),
                uniform_entry(1),
                uniform_entry(2),
                uniform_entry(3),
            ],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("valenx.pbr_forward.pipeline_layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("valenx.pbr_forward.pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: "vs_main",
                compilation_options: Default::default(),
                buffers: &[pbr_vertex_layout()],
            },
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: Some(wgpu::Face::Back),
                polygon_mode: wgpu::PolygonMode::Fill,
                unclipped_depth: false,
                conservative: false,
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: DEPTH_FORMAT,
                depth_write_enabled: true,
                depth_compare: wgpu::CompareFunction::Less,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: "fs_main",
                compilation_options: Default::default(),
                targets: &[Some(wgpu::ColorTargetState {
                    format: target_format,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            multiview: None,
        });

        // Uniform buffers, initialised to zeroed contents.
        let frame_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("valenx.pbr_forward.frame"),
            contents: bytemuck::bytes_of(&PbrFrameUniform::new(
                [[0.0; 4]; 4],
                [0.0; 3],
                [0.0; 3],
                0,
            )),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let material_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("valenx.pbr_forward.material"),
            contents: bytemuck::bytes_of(&PbrMaterialUniform::from_material(
                &Material::default(),
                1.0,
            )),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let light_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("valenx.pbr_forward.lights"),
            contents: bytemuck::bytes_of(&LightArray::from_slice(&[])),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let gi_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("valenx.pbr_forward.gi_probe"),
            contents: bytemuck::bytes_of(&ShL2Uniform::zero()),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("valenx.pbr_forward.bind_group"),
            layout: &bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: frame_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: material_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: light_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: gi_buffer.as_entire_binding(),
                },
            ],
        });

        let initial_capacity: u64 = 4096;
        let vertex_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("valenx.pbr_forward.vertices"),
            size: initial_capacity * std::mem::size_of::<PbrVertex>() as u64,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        Self {
            device,
            queue,
            target_format,
            pipeline,
            bind_group_layout,
            frame_buffer,
            material_buffer,
            light_buffer,
            gi_buffer,
            bind_group,
            vertex_buffer,
            vertex_capacity: initial_capacity,
            vertex_count: 0,
            offscreen: None,
        }
    }

    /// Create / resize the PBR pass's offscreen colour + depth targets
    /// to `size` if the viewport changed.
    fn ensure_offscreen_for(&mut self, renderer: &mut egui_wgpu::Renderer, size: [u32; 2]) {
        let [w, h] = size;
        if w == 0 || h == 0 {
            return;
        }
        if let Some(os) = &self.offscreen {
            if os.size == size {
                return;
            }
        }

        let color_texture = self.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("valenx.pbr_forward.color"),
            size: wgpu::Extent3d {
                width: w,
                height: h,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: self.target_format,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });
        let color_view = color_texture.create_view(&wgpu::TextureViewDescriptor::default());

        let depth_texture = self.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("valenx.pbr_forward.depth"),
            size: wgpu::Extent3d {
                width: w,
                height: h,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: DEPTH_FORMAT,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        });
        let depth_view = depth_texture.create_view(&wgpu::TextureViewDescriptor::default());

        let egui_id = match &self.offscreen {
            Some(prev) => {
                renderer.update_egui_texture_from_wgpu_texture(
                    &self.device,
                    &color_view,
                    wgpu::FilterMode::Linear,
                    prev.egui_id,
                );
                prev.egui_id
            }
            None => renderer.register_native_texture(
                &self.device,
                &color_view,
                wgpu::FilterMode::Linear,
            ),
        };

        self.offscreen = Some(PbrOffscreen {
            size,
            color_texture,
            color_view,
            depth_texture,
            depth_view,
            egui_id,
        });
    }

    /// Grow / write the vertex buffer with this frame's geometry.
    fn upload_vertices(&mut self, vertices: &[PbrVertex]) {
        self.vertex_count = vertices.len() as u32;
        if vertices.is_empty() {
            return;
        }
        let needed_capacity = vertices.len() as u64;
        if needed_capacity > self.vertex_capacity {
            let new_capacity = needed_capacity.next_power_of_two().max(4096);
            self.vertex_buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("valenx.pbr_forward.vertices"),
                size: new_capacity * std::mem::size_of::<PbrVertex>() as u64,
                usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            self.vertex_capacity = new_capacity;
        }
        self.queue
            .write_buffer(&self.vertex_buffer, 0, bytemuck::cast_slice(vertices));
    }

    /// Render `frame` into the PBR pass's offscreen target and return
    /// the [`egui::TextureId`] a caller displays via `egui::Image`.
    ///
    /// Returns `None` if the viewport size is zero. The geometry is
    /// drawn with the Cook-Torrance forward shader; the uniforms
    /// (frame / material / lights / GI probe) are uploaded fresh each
    /// call.
    ///
    /// Like the whole module, this is **GPU-unverified** — it encodes
    /// a correct-by-construction pass but has not been run on
    /// hardware.
    pub fn render(
        &mut self,
        renderer: &mut egui_wgpu::Renderer,
        size: [u32; 2],
        frame: &PbrFrame<'_>,
    ) -> Option<egui::TextureId> {
        if size[0] == 0 || size[1] == 0 {
            return None;
        }
        self.ensure_offscreen_for(renderer, size);
        self.upload_vertices(frame.vertices);

        // Upload the per-frame uniforms.
        let frame_uniform = PbrFrameUniform::new(
            frame.view_proj,
            frame.camera_pos,
            frame.env_color,
            frame.lights.len().min(MAX_LIGHTS) as u32,
        );
        self.queue
            .write_buffer(&self.frame_buffer, 0, bytemuck::bytes_of(&frame_uniform));
        let material_uniform =
            PbrMaterialUniform::from_material(frame.material, frame.ambient_occlusion);
        self.queue.write_buffer(
            &self.material_buffer,
            0,
            bytemuck::bytes_of(&material_uniform),
        );
        let light_array = LightArray::from_slice(frame.lights);
        self.queue
            .write_buffer(&self.light_buffer, 0, bytemuck::bytes_of(&light_array));
        self.queue
            .write_buffer(&self.gi_buffer, 0, bytemuck::bytes_of(&frame.gi_probe));

        let os = self.offscreen.as_ref()?;

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("valenx.pbr_forward.encoder"),
            });
        {
            let mut rpass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("valenx.pbr_forward.pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &os.color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.05,
                            g: 0.05,
                            b: 0.06,
                            a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &os.depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            if self.vertex_count > 0 {
                rpass.set_pipeline(&self.pipeline);
                rpass.set_bind_group(0, &self.bind_group, &[]);
                rpass.set_vertex_buffer(0, self.vertex_buffer.slice(..));
                rpass.draw(0..self.vertex_count, 0..1);
            }
        }
        self.queue.submit(Some(encoder.finish()));
        Some(os.egui_id)
    }
}

/// Sample one [`ShL2Uniform`] GI probe from an [`IrradianceVolume`] at
/// a world position — a convenience for a caller that bakes an
/// irradiance volume and wants the SH probe to feed [`PbrFrame`].
///
/// The volume's trilinearly-blended coefficients at `position` would
/// ideally be packed directly; this v1 samples the volume's *nearest*
/// stored probe (the volume keeps per-probe SH, not a blended-coeff
/// accessor), which is exact at a probe and a close approximation
/// between them. A blended-coefficient accessor on `IrradianceVolume`
/// is a documented follow-up.
pub fn gi_probe_at(volume: &IrradianceVolume, position: [f32; 3]) -> ShL2Uniform {
    // Map the position to the nearest grid index on each axis.
    let nearest = |k: usize| -> usize {
        let span = volume.max[k] - volume.min[k];
        let n = match k {
            0 => volume.dims.0,
            1 => volume.dims.1,
            _ => volume.dims.2,
        };
        if span.abs() < 1e-12 {
            0
        } else {
            let t = ((position[k] - volume.min[k]) / span).clamp(0.0, 1.0);
            (t * (n - 1) as f32).round() as usize
        }
    };
    let ix = nearest(0).min(volume.dims.0 - 1);
    let iy = nearest(1).min(volume.dims.1 - 1);
    let iz = nearest(2).min(volume.dims.2 - 1);
    let idx = ix + iy * volume.dims.0 + iz * volume.dims.0 * volume.dims.1;
    let probe = &volume.probes[idx];
    ShL2Uniform::from_coeffs(&probe.coeffs)
}

/// Build PBR vertex data from a `valenx_viz::TriangleMesh` — three
/// flat-shaded vertices per triangle, each carrying the face normal.
///
/// The same convention as `wgpu_renderer::triangles_to_vertices`, so a
/// caller can feed either render path from one mesh.
pub fn triangles_to_pbr_vertices(mesh: &valenx_viz::TriangleMesh) -> Vec<PbrVertex> {
    let mut out = Vec::with_capacity(mesh.triangles.len() * 3);
    for tri in &mesh.triangles {
        let n = tri.computed_normal();
        for v in &tri.vertices {
            out.push(PbrVertex {
                position: *v,
                normal: n,
            });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_render_bridge::ShOrder;

    /// The PBR vertex is `#[repr(C)]` and exactly 24 bytes — the same
    /// layout as the viewport pipeline's vertex, so geometry buffers
    /// are interchangeable between the two render paths.
    #[test]
    fn pbr_vertex_is_24_bytes() {
        assert_eq!(std::mem::size_of::<PbrVertex>(), 24);
    }

    /// The fixed GPU light array is the WGSL-matching size — one
    /// `LightUniform` (32 bytes) per `MAX_LIGHTS` slot.
    #[test]
    fn light_array_has_the_wgsl_size() {
        assert_eq!(
            std::mem::size_of::<LightArray>(),
            32 * MAX_LIGHTS,
            "light array must match array<LightUniform, MAX_LIGHTS>"
        );
    }

    /// `LightArray::from_slice` pads the unused tail with inactive
    /// lights and never overruns the fixed array.
    #[test]
    fn light_array_pads_and_clamps() {
        // Fewer lights than slots — the tail is inactive.
        let one = LightArray::from_slice(&[PbrLightUniform::point(
            [1.0, 2.0, 3.0],
            [1.0, 1.0, 1.0],
            10.0,
        )]);
        assert_eq!(one.lights[0].direction_or_pos, [1.0, 2.0, 3.0, 1.0]);
        assert_eq!(one.lights[MAX_LIGHTS - 1].color_intensity, [0.0; 4]);
        // More lights than slots — only the first MAX_LIGHTS are kept,
        // no panic.
        let many =
            vec![PbrLightUniform::directional([0.0, -1.0, 0.0], [1.0; 3], 1.0); MAX_LIGHTS + 5];
        let arr = LightArray::from_slice(&many);
        assert_eq!(arr.lights.len(), MAX_LIGHTS);
    }

    /// `triangles_to_pbr_vertices` emits three vertices per triangle,
    /// each sharing the face normal.
    #[test]
    fn triangles_to_pbr_vertices_three_per_triangle() {
        use valenx_viz::{StlFormat, StlTriangle, TriangleMesh};
        let mesh = TriangleMesh {
            format: Some(StlFormat::Ascii),
            name: None,
            triangles: vec![StlTriangle {
                normal: [0.0, 0.0, 1.0],
                vertices: [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            }],
        };
        let verts = triangles_to_pbr_vertices(&mesh);
        assert_eq!(verts.len(), 3);
        assert_eq!(verts[0].normal, verts[1].normal);
        assert_eq!(verts[1].normal, verts[2].normal);
    }

    /// `gi_probe_at` returns the nearest probe's coefficients — exact
    /// when the query lands on a probe.
    #[test]
    fn gi_probe_at_reads_back_a_probe() {
        let mut vol = IrradianceVolume::new([0.0; 3], [1.0; 3], (2, 2, 2), ShOrder::L2).unwrap();
        // Bake a directional-ish scene so the probes are non-trivial.
        vol.bake(64, |_o, d| [d[0].max(0.0), 0.2, 0.4]);
        // A query at the (0,0,0) corner should return that probe.
        let probe = gi_probe_at(&vol, [0.0, 0.0, 0.0]);
        let expected = &vol.probes[0];
        for (i, c) in probe.coeffs.iter().enumerate() {
            assert!(
                (c[0] - expected.coeffs[i][0]).abs() < 1e-5,
                "gi_probe_at coefficient {i} should match the corner probe"
            );
        }
    }
}
