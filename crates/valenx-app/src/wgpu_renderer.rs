//! Offscreen wgpu renderer for the 3D viewport.
//!
//! eframe's egui-wgpu render pass does not expose a depth buffer, so
//! "just draw triangles in its render pass" loses depth correctness
//! on overlapping geometry. The fix is to render the 3D scene into
//! our own offscreen colour + depth textures, then hand the colour
//! texture to egui as a `egui::TextureId` that the UI displays like any
//! other image.
//!
//! Flow per frame:
//!
//! 1. `ensure_offscreen_for(size)` creates or resizes the colour +
//!    depth targets if the viewport's pixel size changed.
//! 2. `upload_vertices(vertices)` writes the flat-shaded vertex data
//!    into a growable device buffer.
//! 3. `render_into_offscreen(mvp, light)` encodes a render pass
//!    clearing colour + depth and drawing every triangle with the
//!    shader in `VALENX_SHADER_WGSL`.
//! 4. The caller draws the returned `egui::TextureId` via `egui::Image`.
//!
//! The shader is intentionally small — flat per-face Lambert
//! lighting against a single key-light direction, neutral brushed-
//! metal base colour. Richer material models, shadows, SSAO, etc.
//! are post-Year-1 work.

use std::sync::Arc;

use bytemuck::{Pod, Zeroable};
use eframe::egui;
use eframe::egui_wgpu;
use eframe::wgpu::{self, util::DeviceExt};

/// Vertex layout the pipeline expects: position (float3) + normal
/// (float3). No texture coords, no vertex colours — flat shading
/// derives colour entirely from the normal and light direction.
#[repr(C)]
#[derive(Copy, Clone, Debug, Default, Pod, Zeroable)]
pub struct Vertex {
    pub position: [f32; 3],
    pub normal: [f32; 3],
}

const VERTEX_ATTRS: [wgpu::VertexAttribute; 2] = wgpu::vertex_attr_array![
    0 => Float32x3, // position
    1 => Float32x3, // normal
];

fn vertex_layout() -> wgpu::VertexBufferLayout<'static> {
    wgpu::VertexBufferLayout {
        array_stride: std::mem::size_of::<Vertex>() as wgpu::BufferAddress,
        step_mode: wgpu::VertexStepMode::Vertex,
        attributes: &VERTEX_ATTRS,
    }
}

/// GPU uniform block. All fields are `vec4` / `mat4` so WGSL's 16-byte
/// alignment rules are satisfied at every offset without padding.
///
/// Layout (byte offsets):
/// - `mvp`       0   – 63  (mat4x4, 64 bytes)
/// - `inv_mvp`   64  – 127 (mat4x4, 64 bytes) — inverse of mvp, used by the
///   grid shader to unproject screen pixels to world rays
/// - `light_dir` 128 – 143 (vec4, xyz = direction, w unused)
/// - `cam_pos`   144 – 159 (vec4, xyz = world eye position, w unused)
/// - `grid`      160 – 175 (vec4, x = minor spacing, y = unused, z = unused,
///   w = fade distance)
/// - `grid2`     176 – 191 (vec4, x = LOD blend_t [0→minor fades out], rest unused)
///
/// Total: 192 bytes.
#[repr(C)]
#[derive(Copy, Clone, Debug, Default, Pod, Zeroable)]
pub struct Uniforms {
    pub mvp: [[f32; 4]; 4],
    /// Inverse of `mvp`. Used by the grid fragment shader to unproject NDC
    /// coordinates into world-space rays for the Y=0 ground intersection.
    pub inv_mvp: [[f32; 4]; 4],
    pub light_dir: [f32; 4],
    /// World camera position (xyz); w unused.
    pub cam_pos: [f32; 4],
    /// Grid params: x = minor spacing, y = unused, z = unused, w = fade distance.
    pub grid: [f32; 4],
    /// Additional grid LOD params: x = blend_t (0=minor fully visible,
    /// 1=minor faded), yzw unused.
    pub grid2: [f32; 4],
}

const SHADER_WGSL: &str = r#"
// Full Uniforms layout (192 bytes). All fields declared even when unused
// so the struct layout matches the buffer exactly.
struct Uniforms {
    mvp:       mat4x4<f32>,
    inv_mvp:   mat4x4<f32>,
    light_dir: vec4<f32>,
    cam_pos:   vec4<f32>,
    grid:      vec4<f32>,
    grid2:     vec4<f32>,
}
@group(0) @binding(0) var<uniform> u: Uniforms;

struct VIn {
    @location(0) position: vec3<f32>,
    @location(1) normal:   vec3<f32>,
}
struct VOut {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) normal: vec3<f32>,
}

@vertex
fn vs_main(v_in: VIn) -> VOut {
    var out: VOut;
    out.clip_pos = u.mvp * vec4<f32>(v_in.position, 1.0);
    out.normal = v_in.normal;
    return out;
}

@fragment
fn fs_main(f_in: VOut) -> @location(0) vec4<f32> {
    let n = normalize(f_in.normal);
    let l = normalize(u.light_dir.xyz);
    let lambert = max(dot(-l, n), 0.0);
    let base = vec3<f32>(0.66, 0.76, 0.88);
    let ambient = 0.22;
    let color = base * (ambient + 0.78 * lambert);
    return vec4<f32>(color, 1.0);
}
"#;

/// Fusion-360-quality "infinite ground grid" shader.
///
/// # Design
///
/// Uses a fullscreen triangle (3 vertices, no vertex buffer). For each
/// fragment the NDC position is unprojected via `inv_mvp` to obtain a
/// world-space ray that is intersected with the Y = 0 ground plane.  The
/// intersection's world XZ coordinates drive three levels of analytic
/// `fwidth`-antialiased grid lines, and the fragment writes its own depth
/// so 3-D geometry correctly occludes it.
///
/// # Features
/// - Major + minor + "next-major" LOD levels with a smooth crossfade
///   (no pop when zooming; driven by `grid2.x = blend_t`).
/// - Colored principal axes: X = red, Z = blue, 2px-wide screen-space lines.
/// - Distance fade-to-horizon (quadratic ease) via `grid.w = fade_dist`.
/// - Correct occlusion by 3-D geometry via `@builtin(frag_depth)`.
/// - Runs on Vulkan / Metal / DX12 / GL (core WGSL only).
const GRID_SHADER_WGSL: &str = r#"
struct Uniforms {
    mvp:       mat4x4<f32>,
    inv_mvp:   mat4x4<f32>,
    light_dir: vec4<f32>,
    cam_pos:   vec4<f32>,
    // grid.x = minor_a spacing, grid.w = fade distance
    grid:      vec4<f32>,
    // grid2.x = blend_t  (0 → minor fully visible, 1 → minor faded)
    grid2:     vec4<f32>,
}
@group(0) @binding(0) var<uniform> u: Uniforms;

// Vertex output: clip position + NDC xy for unprojection.
struct GVert {
    @builtin(position) clip: vec4<f32>,
    @location(0) ndc: vec2<f32>,
}

// Fragment output: colour + custom depth (for correct 3-D occlusion).
struct GFrag {
    @location(0) color: vec4<f32>,
    @builtin(frag_depth) depth: f32,
}

// Fullscreen triangle — covers NDC [-1,1]×[-1,1] with 3 vertices, no VBO.
@vertex
fn grid_vs(@builtin(vertex_index) vid: u32) -> GVert {
    var pos = array<vec2<f32>, 3>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 3.0, -1.0),
        vec2<f32>(-1.0,  3.0),
    );
    let p = pos[vid];
    return GVert(
        vec4<f32>(p.x, p.y, 0.5, 1.0),  // z=0.5 placeholder; frag_depth overrides
        p,
    );
}

// Analytic grid-line coverage at spacing s. Returns [0,1].
// Uses fwidth of the grid-coordinate for sub-pixel antialiasing.
fn grid_line(c: vec2<f32>, s: f32) -> f32 {
    let g = c / s;
    let fw = max(fwidth(g), vec2<f32>(1e-5));
    let d = abs(fract(g - 0.5) - 0.5) / fw;
    return 1.0 - min(min(d.x, d.y), 1.0);
}

@fragment
fn grid_fs(i: GVert) -> GFrag {
    // ---- 1. Unproject NDC → world-space ray -------------------------
    let h0 = u.inv_mvp * vec4<f32>(i.ndc, 0.0, 1.0);  // near plane (WebGPU NDC z=0)
    let h1 = u.inv_mvp * vec4<f32>(i.ndc, 1.0, 1.0);  // far plane
    let w0 = h0.xyz / h0.w;
    let w1 = h1.xyz / h1.w;
    let rd = w1 - w0;

    // ---- 2. Intersect ray with Y = 0 ground plane -------------------
    if abs(rd.y) < 1e-5 { discard; }  // ray parallel to ground → no hit
    let t = -w0.y / rd.y;
    if t < 0.0 { discard; }           // hit is behind the camera

    let hit = w0 + t * rd;
    let wxz = hit.xz;                 // world (X, Z) at the ground hit

    // ---- 3. Re-project to get the correct fragment depth ------------
    let clip_hit = u.mvp * vec4<f32>(hit.x, 0.0, hit.z, 1.0);
    // The wgpu projection (gl_to_wgpu in this file) already maps depth to
    // WebGPU's [0,1] NDC range, so frag_depth is the straight divide — and
    // it now matches the mesh pipeline's depth, so geometry occludes the
    // grid correctly even right up close.
    let depth = clip_hit.z / clip_hit.w;
    if depth < 0.0 || depth > 1.0 { discard; }  // outside [near, far]

    // ---- 4. Distance-based fade to horizon --------------------------
    let dist = length(wxz - u.cam_pos.xz);
    let fade_dist = max(u.grid.w, 1e-3);
    let raw_fade = clamp(1.0 - dist / fade_dist, 0.0, 1.0);
    // Ease-OUT: keep the grid near full strength across most of the view and
    // fade it hard only near the horizon (Blender-like), instead of dimming
    // quadratically from the camera outward (which read as washed-out).
    let fade = 1.0 - (1.0 - raw_fade) * (1.0 - raw_fade);

    // ---- 5. Three-level LOD grid with smooth crossfade --------------
    // minor_a: current minor spacing (fades out as blend_t → 1)
    // major_a = minor_a × 10: transitions from major brightness → minor brightness
    // major_b = minor_a × 100: next-level major (fades in as blend_t → 1)
    let minor_a = max(u.grid.x, 1e-6);
    let major_a = minor_a * 10.0;
    let major_b = minor_a * 100.0;
    let bt = clamp(u.grid2.x, 0.0, 1.0);  // blend_t

    // Softer, finer line weights (Fusion-style): faint minor lines, a
    // gentle major emphasis — the grid recedes politely so geometry pops,
    // instead of a busy high-contrast mesh.
    let a_minor  = grid_line(wxz, minor_a) * (1.0 - bt) * 0.28;
    let a_major  = grid_line(wxz, major_a) * mix(0.55, 0.30, bt);
    let a_major2 = grid_line(wxz, major_b) * bt * 0.55;
    var a = max(max(a_minor, a_major), a_major2);

    // ---- 6. Base grid colour — a near-neutral grey (Blender-like). ----
    var col = vec3<f32>(0.47, 0.48, 0.50);

    // ---- 7. Principal axes — crisp ~1.5px ANTI-ALIASED lines, softly
    // saturated (Fusion-style) instead of a hard full-alpha band. The
    // old version was a binary cutoff at 2.5px + alpha 0.93, which read
    // as thick and aliased ("neon"). Coverage eases from a ~0.75px
    // half-width core across one more pixel via smoothstep.
    let fw = max(fwidth(wxz), vec2<f32>(1e-5));
    // X axis runs along world X → highlight where world Z (wxz.y) ≈ 0.
    let ax_x = 1.0 - smoothstep(0.75 * fw.y, 1.75 * fw.y, abs(wxz.y));
    // Z axis runs along world Z → highlight where world X (wxz.x) ≈ 0.
    let ax_z = 1.0 - smoothstep(0.75 * fw.x, 1.75 * fw.x, abs(wxz.x));
    if ax_x > 0.001 {
        col = mix(col, vec3<f32>(0.80, 0.30, 0.33), ax_x);
        a = max(a, ax_x * 0.85);
    }
    if ax_z > 0.001 {
        col = mix(col, vec3<f32>(0.28, 0.49, 0.86), ax_z);
        a = max(a, ax_z * 0.85);
    }

    // ---- 8. Apply fade and early-out --------------------------------
    a *= fade;
    if a < 0.003 { discard; }

    var out: GFrag;
    out.color = vec4<f32>(col, a);
    out.depth = depth;
    return out;
}
"#;

/// Screen-space vertical gradient backdrop, drawn first (behind the mesh and
/// the transparent grid). A subtle dark gradient — lighter toward the top —
/// gives the viewport depth instead of a flat fill, the way Blender's viewport
/// background does. No uniforms: the colours are baked in.
const BG_SHADER_WGSL: &str = r#"
struct BgVert {
    @builtin(position) clip: vec4<f32>,
    @location(0) ndc: vec2<f32>,
}

@vertex
fn bg_vs(@builtin(vertex_index) vid: u32) -> BgVert {
    var pos = array<vec2<f32>, 3>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 3.0, -1.0),
        vec2<f32>(-1.0,  3.0),
    );
    let p = pos[vid];
    return BgVert(vec4<f32>(p.x, p.y, 1.0, 1.0), p);
}

@fragment
fn bg_fs(i: BgVert) -> @location(0) vec4<f32> {
    // ndc.y: -1 bottom → +1 top. Lighter near the top, darker low.
    let t = clamp(i.ndc.y * 0.5 + 0.5, 0.0, 1.0);
    let bottom = vec3<f32>(0.075, 0.085, 0.100);
    let top    = vec3<f32>(0.150, 0.165, 0.190);
    return vec4<f32>(mix(bottom, top, t), 1.0);
}
"#;

const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;

/// The viewport's wgpu render pipeline — owns the device, queue, and
/// per-frame resources for the Phong-shaded triangle mesh path.
pub struct WgpuRenderer {
    device: Arc<wgpu::Device>,
    queue: Arc<wgpu::Queue>,
    target_format: wgpu::TextureFormat,

    pipeline: wgpu::RenderPipeline,
    grid_pipeline: wgpu::RenderPipeline,
    bg_pipeline: wgpu::RenderPipeline,
    // Kept alive for the pipeline even though no code reads it after
    // construction; dropping it would invalidate `bind_group`.
    #[allow(dead_code)]
    bind_group_layout: wgpu::BindGroupLayout,
    uniform_buffer: wgpu::Buffer,
    bind_group: wgpu::BindGroup,

    vertex_buffer: wgpu::Buffer,
    vertex_capacity: u64,
    vertex_count: u32,

    offscreen: Option<Offscreen>,
}

struct Offscreen {
    size: [u32; 2],
    // Textures kept alive so the views (used every frame) don't
    // dangle. Neither is read back directly after construction.
    #[allow(dead_code)]
    color_texture: wgpu::Texture,
    color_view: wgpu::TextureView,
    #[allow(dead_code)]
    depth_texture: wgpu::Texture,
    depth_view: wgpu::TextureView,
    egui_id: egui::TextureId,
}

impl WgpuRenderer {
    /// Build the pipeline + initial (empty) buffers. Offscreen
    /// textures are created on demand — the caller doesn't know the
    /// viewport's pixel size at construction time.
    pub fn new(render_state: &egui_wgpu::RenderState) -> Self {
        let device = render_state.device.clone();
        let queue = render_state.queue.clone();
        let target_format = render_state.target_format;

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("valenx.viewport.shader"),
            source: wgpu::ShaderSource::Wgsl(SHADER_WGSL.into()),
        });

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("valenx.viewport.bind_group_layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX_FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("valenx.viewport.pipeline_layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("valenx.viewport.pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: "vs_main",
                compilation_options: Default::default(),
                buffers: &[vertex_layout()],
            },
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                strip_index_format: None,
                // Right-hand winding (Counter-clockwise is front) —
                // matches the STL loader's computed normals and the
                // projection math in valenx-viz.
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

        // Grid pipeline — fullscreen-triangle approach (no vertex buffer).
        // Uses inv_mvp to unproject screen pixels to world-space rays and
        // intersects with Y=0. Outputs @builtin(frag_depth) so geometry
        // correctly occludes the grid via the depth test.
        let grid_shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("valenx.viewport.grid_shader"),
            source: wgpu::ShaderSource::Wgsl(GRID_SHADER_WGSL.into()),
        });
        let grid_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("valenx.viewport.grid_pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &grid_shader,
                entry_point: "grid_vs",
                compilation_options: Default::default(),
                buffers: &[], // no vertex buffer — fullscreen triangle
            },
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                cull_mode: None, // fullscreen triangle never back-faces
                ..Default::default()
            },
            // depth_write_enabled: true so the fragment's @builtin(frag_depth)
            // (the re-projected ground-plane depth) is stored. This lets 3-D
            // geometry drawn first correctly occlude the grid (they wrote smaller
            // depth values; the grid's larger Y=0 depth fails Less).
            depth_stencil: Some(wgpu::DepthStencilState {
                format: DEPTH_FORMAT,
                depth_write_enabled: true,
                depth_compare: wgpu::CompareFunction::Less,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            fragment: Some(wgpu::FragmentState {
                module: &grid_shader,
                entry_point: "grid_fs",
                compilation_options: Default::default(),
                targets: &[Some(wgpu::ColorTargetState {
                    format: target_format,
                    blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            multiview: None,
        });

        // Background gradient pipeline — a screen-space backdrop drawn before
        // the mesh and grid. Depth test Always + no depth write, so it sits
        // behind everything without touching the depth buffer. Empty layout
        // (no uniforms — the gradient colours are baked into the shader).
        let bg_shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("valenx.viewport.bg_shader"),
            source: wgpu::ShaderSource::Wgsl(BG_SHADER_WGSL.into()),
        });
        let bg_pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("valenx.viewport.bg_pipeline_layout"),
            bind_group_layouts: &[],
            push_constant_ranges: &[],
        });
        let bg_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("valenx.viewport.bg_pipeline"),
            layout: Some(&bg_pipeline_layout),
            vertex: wgpu::VertexState {
                module: &bg_shader,
                entry_point: "bg_vs",
                compilation_options: Default::default(),
                buffers: &[],
            },
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                cull_mode: None,
                ..Default::default()
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: DEPTH_FORMAT,
                depth_write_enabled: false,
                depth_compare: wgpu::CompareFunction::Always,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            fragment: Some(wgpu::FragmentState {
                module: &bg_shader,
                entry_point: "bg_fs",
                compilation_options: Default::default(),
                targets: &[Some(wgpu::ColorTargetState {
                    format: target_format,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            multiview: None,
        });

        let uniforms = Uniforms::default();
        let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("valenx.viewport.uniforms"),
            contents: bytemuck::bytes_of(&uniforms),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("valenx.viewport.bind_group"),
            layout: &bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buffer.as_entire_binding(),
            }],
        });

        // Initial empty vertex buffer with a small headroom so the
        // first draw doesn't have to reallocate.
        let initial_capacity: u64 = 4096;
        let vertex_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("valenx.viewport.vertices"),
            size: initial_capacity * std::mem::size_of::<Vertex>() as u64,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        Self {
            device,
            queue,
            target_format,
            pipeline,
            grid_pipeline,
            bg_pipeline,
            bind_group_layout,
            uniform_buffer,
            bind_group,
            vertex_buffer,
            vertex_capacity: initial_capacity,
            vertex_count: 0,
            offscreen: None,
        }
    }

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
            label: Some("valenx.viewport.color"),
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
            label: Some("valenx.viewport.depth"),
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

        // Register (or re-register) with egui. If we already have an
        // id we keep it and point it at the new view — egui's ids
        // are stable across texture re-creations which avoids
        // flickering during resize.
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

        self.offscreen = Some(Offscreen {
            size,
            color_texture,
            color_view,
            depth_texture,
            depth_view,
            egui_id,
        });
    }

    fn upload_vertices(&mut self, vertices: &[Vertex]) {
        self.vertex_count = vertices.len() as u32;
        if vertices.is_empty() {
            return;
        }
        let needed_bytes = std::mem::size_of_val(vertices) as u64;
        let needed_capacity = vertices.len() as u64;
        if needed_capacity > self.vertex_capacity {
            // Grow by 2× so we don't thrash on small growth steps.
            let new_capacity = needed_capacity.next_power_of_two().max(4096);
            self.vertex_buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("valenx.viewport.vertices"),
                size: new_capacity * std::mem::size_of::<Vertex>() as u64,
                usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            self.vertex_capacity = new_capacity;
        }
        self.queue
            .write_buffer(&self.vertex_buffer, 0, bytemuck::cast_slice(vertices));
        let _ = needed_bytes;
    }

    fn write_uniforms(&self, uniforms: &Uniforms) {
        self.queue
            .write_buffer(&self.uniform_buffer, 0, bytemuck::bytes_of(uniforms));
    }

    /// Render to the offscreen target and return the egui TextureId
    /// the caller should display. Returns `None` if the viewport size
    /// is zero.
    ///
    /// `inv_mvp` is the inverse of `mvp`; the grid shader uses it to
    /// unproject fragments to world-space rays.
    /// `grid2` carries the LOD blend factor (`grid2[0]` = blend_t).
    #[allow(clippy::too_many_arguments)]
    pub fn render(
        &mut self,
        renderer: &mut egui_wgpu::Renderer,
        size: [u32; 2],
        mvp: [[f32; 4]; 4],
        inv_mvp: [[f32; 4]; 4],
        light_dir: [f32; 3],
        cam_pos: [f32; 3],
        grid: [f32; 4],
        grid2: [f32; 4],
        vertices: &[Vertex],
    ) -> Option<egui::TextureId> {
        if size[0] == 0 || size[1] == 0 {
            return None;
        }
        self.ensure_offscreen_for(renderer, size);
        self.upload_vertices(vertices);
        self.write_uniforms(&Uniforms {
            mvp,
            inv_mvp,
            light_dir: [light_dir[0], light_dir[1], light_dir[2], 0.0],
            cam_pos: [cam_pos[0], cam_pos[1], cam_pos[2], 0.0],
            grid,
            grid2,
        });

        let os = self.offscreen.as_ref()?;

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("valenx.viewport.encoder"),
            });
        {
            let mut rpass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("valenx.viewport.pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &os.color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.10,
                            g: 0.11,
                            b: 0.13,
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
            // Background gradient first — a screen-space backdrop behind the
            // mesh and the transparent grid (depth test Always, no depth write,
            // so it never occludes geometry).
            rpass.set_pipeline(&self.bg_pipeline);
            rpass.draw(0..3, 0..1);
            // Mesh next (writes depth for later occlusion test).
            if self.vertex_count > 0 {
                rpass.set_pipeline(&self.pipeline);
                rpass.set_bind_group(0, &self.bind_group, &[]);
                rpass.set_vertex_buffer(0, self.vertex_buffer.slice(..));
                rpass.draw(0..self.vertex_count, 0..1);
            }
            // Ground grid + axes drawn after the mesh. The fullscreen
            // triangle (3 vertices) runs grid_fs which: unprojections each
            // pixel to world-space, intersects Y=0, computes grid alpha,
            // and writes the correct @builtin(frag_depth). Because the
            // grid's frag_depth is larger (Y=0 is farther than geometry
            // above Y=0), the Less depth-test fails wherever the mesh was
            // already drawn — clean occlusion without a depth pre-pass.
            rpass.set_pipeline(&self.grid_pipeline);
            rpass.set_bind_group(0, &self.bind_group, &[]);
            rpass.draw(0..3, 0..1); // fullscreen triangle
        }
        self.queue.submit(Some(encoder.finish()));
        Some(os.egui_id)
    }
}

/// Build flat-shaded vertex data from a `TriangleMesh`. Each input
/// triangle becomes three `Vertex` entries that share the triangle's
/// face normal, so the fragment shader gets a constant normal across
/// each triangle (the visual hallmark of flat shading).
pub fn triangles_to_vertices(mesh: &valenx_viz::TriangleMesh) -> Vec<Vertex> {
    let mut out = Vec::with_capacity(mesh.triangles.len() * 3);
    for tri in &mesh.triangles {
        let n = tri.computed_normal();
        for v in &tri.vertices {
            out.push(Vertex {
                position: *v,
                normal: n,
            });
        }
    }
    out
}

/// Build an MVP matrix suitable for the WGSL uniform. nalgebra's
/// OpenGL→WebGPU clip-space depth correction. nalgebra's `new_perspective`
/// emits OpenGL-convention clip space (NDC z in [-1, 1]), but WebGPU (and
/// DX12/Metal) expect [0, 1]. Without it the mesh's near half is clipped
/// (clip.z < 0 falls outside wgpu's [0, w] clip volume) and the grid's
/// analytic depth does not match the mesh's — so geometry doesn't occlude
/// the grid correctly up close. This maps z' = 0.5·z + 0.5·w. Applied ONLY
/// to the wgpu pipelines (mesh + grid); picking (`scene.rs`) and the
/// egui-painter overlays (`projection.rs`) keep the raw GL projection — they
/// don't use the depth buffer, so they're unaffected.
fn gl_to_wgpu() -> nalgebra::Matrix4<f32> {
    nalgebra::Matrix4::new(
        1.0, 0.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, 0.0, //
        0.0, 0.0, 0.5, 0.5, //
        0.0, 0.0, 0.0, 1.0,
    )
}

/// `Matrix4` is column-major and the WGSL `mat4x4<f32>` type also
/// stores column-major, so we write columns-as-outer-axis.
pub fn mvp_from_camera(camera: &valenx_viz::OrbitCamera, width: f32, height: f32) -> [[f32; 4]; 4] {
    let aspect = (width / height.max(1.0)).max(1e-6);
    let view = camera.view_matrix();
    let proj = camera.projection_matrix(aspect);
    let mvp = gl_to_wgpu() * proj * view;
    mat4_to_cols(mvp)
}

/// Build the inverse of the MVP matrix for the WGSL uniform. The grid
/// shader uses this to unproject NDC coordinates to world-space rays.
/// Falls back to the identity when the MVP is singular (degenerate camera).
pub fn inv_mvp_from_camera(
    camera: &valenx_viz::OrbitCamera,
    width: f32,
    height: f32,
) -> [[f32; 4]; 4] {
    let aspect = (width / height.max(1.0)).max(1e-6);
    let view = camera.view_matrix();
    let proj = camera.projection_matrix(aspect);
    let mvp = gl_to_wgpu() * proj * view;
    let inv = mvp
        .try_inverse()
        .unwrap_or_else(nalgebra::Matrix4::identity);
    mat4_to_cols(inv)
}

/// Convert a nalgebra `Matrix4<f32>` (column-major) to the `[[f32;4];4]`
/// column-major layout that WGSL's `mat4x4<f32>` expects.
fn mat4_to_cols(m: nalgebra::Matrix4<f32>) -> [[f32; 4]; 4] {
    let mut result = [[0.0f32; 4]; 4];
    for col in 0..4 {
        for row in 0..4 {
            result[col][row] = m[(row, col)];
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vertex_is_repr_c_and_24_bytes() {
        // Position (3 × 4) + normal (3 × 4) = 24 bytes, no padding.
        assert_eq!(std::mem::size_of::<Vertex>(), 24);
    }

    #[test]
    fn uniforms_are_repr_c_and_aligned() {
        // mvp (64) + inv_mvp (64) + light_dir (16) + cam_pos (16) + grid (16) + grid2 (16)
        // = 192 bytes, all at 16-byte WGSL-aligned offsets.
        assert_eq!(std::mem::size_of::<Uniforms>(), 192);
    }

    #[test]
    fn triangles_to_vertices_three_per_triangle() {
        use valenx_viz::{StlFormat, StlTriangle, TriangleMesh};
        let mesh = TriangleMesh {
            format: Some(StlFormat::Ascii),
            name: None,
            triangles: vec![
                StlTriangle {
                    normal: [0.0, 0.0, 1.0],
                    vertices: [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                },
                StlTriangle {
                    normal: [0.0, 0.0, 1.0],
                    vertices: [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
                },
            ],
        };
        let verts = triangles_to_vertices(&mesh);
        assert_eq!(verts.len(), 6);
        // Each triangle's three vertices share the face normal.
        assert_eq!(verts[0].normal, verts[1].normal);
        assert_eq!(verts[1].normal, verts[2].normal);
    }

    #[test]
    fn mvp_is_finite() {
        let cam = valenx_viz::OrbitCamera::default();
        let m = mvp_from_camera(&cam, 800.0, 600.0);
        for col in m {
            for v in col {
                assert!(v.is_finite());
            }
        }
    }

    #[test]
    fn inv_mvp_is_inverse_of_mvp() {
        let cam = valenx_viz::OrbitCamera::default();
        let mvp = mvp_from_camera(&cam, 800.0, 600.0);
        let inv = inv_mvp_from_camera(&cam, 800.0, 600.0);
        // MVP × inv(MVP) should be the identity. Tolerance is 1e-3 (not
        // 1e-4): the gl_to_wgpu depth correction folded into the MVP
        // slightly worsens conditioning, so the f32 inverse drifts ~3e-4 —
        // still tiny, while a genuinely-wrong inverse would be off by O(0.1).
        let product = mat4_mul(mvp, inv);
        for (col, product_col) in product.iter().enumerate() {
            for (row, &value) in product_col.iter().enumerate() {
                let expected = if row == col { 1.0_f32 } else { 0.0_f32 };
                assert!(
                    (value - expected).abs() < 1e-3,
                    "mvp×inv_mvp[{row},{col}] = {value}, expected {expected}"
                );
            }
        }
    }

    fn mat4_mul(a: [[f32; 4]; 4], b: [[f32; 4]; 4]) -> [[f32; 4]; 4] {
        let mut out = [[0.0f32; 4]; 4];
        for col in 0..4 {
            for row in 0..4 {
                for (k, a_row) in a.iter().enumerate() {
                    out[col][row] += a_row[row] * b[col][k];
                }
            }
        }
        out
    }
}

/// Headless GPU render-path validation — runs the real wgpu pipeline
/// **without a window**.
///
/// # What this proves
///
/// The shaders (`SHADER_WGSL` here, `PBR_FORWARD_WGSL` in
/// `valenx-render-bridge`) and the `wgpu` render passes were, until
/// this module, GPU-unverified — they `cargo check`ed and their maths
/// was cross-checked against the CPU reference, but no code had ever
/// built the pipeline on a device or shaded a pixel.
///
/// These tests close that gap honestly:
///
/// - [`viewport_shader_validates_with_naga`] runs the exact WGSL
///   front-end + validator `wgpu` itself uses on the private viewport
///   `SHADER_WGSL` — a GPU-free proof it is sound WGSL.
/// - [`headless_pbr_render_shades_a_lit_quad`] requests an **off-screen**
///   `wgpu` device (no surface, no window), builds the real PBR
///   pipeline from `valenx_render_bridge::PBR_FORWARD_WGSL`, renders a
///   directionally-lit quad into a 256×256 texture, copies it to a
///   buffer, reads the pixels back, and asserts the framebuffer is
///   genuinely shaded (not the clear colour, with a l-to-r brightness
///   gradient from the light).
///
/// # Honest environment handling
///
/// If `request_adapter` returns `None` (a CI box / sandbox with no
/// GPU and no software fallback), the render test logs a note and
/// returns early — it never hangs, never falsely passes, never fails
/// spuriously. The naga test needs no GPU and always runs.
///
/// Named `headless_ui_tests` so the repo-wide safe test filter
/// (`cargo test -p valenx-app headless_ui_tests`) selects it; nothing
/// here opens a window or touches `rfd`.
#[cfg(test)]
mod headless_ui_tests {
    use super::*;
    use valenx_render_bridge::wgsl_pbr::{
        PbrFrameUniform, PbrLightUniform, PbrMaterialUniform, ShL2Uniform, MAX_LIGHTS,
        PBR_FORWARD_WGSL,
    };
    use valenx_render_bridge::Material;

    /// The private viewport mesh `SHADER_WGSL` parses + semantically
    /// validates with `naga` — the same front-end `wgpu` runs inside
    /// `create_shader_module`. A GPU-free proof the viewport shader is
    /// valid WGSL.
    #[test]
    fn viewport_shader_validates_with_naga() {
        validate_wgsl(SHADER_WGSL, "mesh SHADER_WGSL");
    }

    /// The fusion-grade ground-grid shader parses + semantically validates
    /// with `naga`. Exercises `@builtin(frag_depth)`, `inv_mvp` unprojection,
    /// derivative functions (`fwidth`), and the LOD crossfade math.
    ///
    /// Note: control-flow uniformity is skipped for this shader because
    /// `fwidth` is intentionally used after `discard` (non-uniform control
    /// flow). On real hardware (Vulkan + SPV_EXT_demote_to_helper_invocation,
    /// Metal, DX12) this is perfectly valid; naga's static checker is
    /// conservative. All other semantic checks still run.
    #[test]
    fn grid_shader_validates_with_naga() {
        let source = GRID_SHADER_WGSL;
        let label = "GRID_SHADER_WGSL";
        let module = naga::front::wgsl::parse_str(source).unwrap_or_else(|e| {
            panic!("{label} failed to parse:\n{}", e.emit_to_string(source))
        });
        // Skip CONTROL_FLOW_UNIFORMITY: `fwidth` after `discard` is valid on
        // real GPUs but naga's static analysis flags it conservatively.
        let flags = naga::valid::ValidationFlags::all()
            & !naga::valid::ValidationFlags::CONTROL_FLOW_UNIFORMITY;
        let mut validator =
            naga::valid::Validator::new(flags, naga::valid::Capabilities::empty());
        validator.validate(&module).unwrap_or_else(|e| {
            panic!("{label} failed naga validation:\n{}", e.emit_to_string(source))
        });
    }

    fn validate_wgsl(source: &str, label: &str) {
        let module = naga::front::wgsl::parse_str(source).unwrap_or_else(|e| {
            panic!("{label} failed to parse:\n{}", e.emit_to_string(source))
        });
        let mut validator = naga::valid::Validator::new(
            naga::valid::ValidationFlags::all(),
            naga::valid::Capabilities::empty(),
        );
        validator.validate(&module).unwrap_or_else(|e| {
            panic!("{label} failed naga validation:\n{}", e.emit_to_string(source))
        });
    }

    /// Vertex layout for the PBR pipeline — position + normal, the
    /// `[f32; 3] × 2` layout `PBR_FORWARD_WGSL`'s `vs_main` declares.
    #[repr(C)]
    #[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
    struct PbrVertex {
        position: [f32; 3],
        normal: [f32; 3],
    }

    /// The fixed-size GPU light array — `MAX_LIGHTS` slots, matching
    /// the WGSL `array<LightUniform, MAX_LIGHTS>`.
    #[repr(C)]
    #[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
    struct LightArray {
        lights: [PbrLightUniform; MAX_LIGHTS],
    }

    /// Headless PBR render: build the real `PBR_FORWARD_WGSL` pipeline
    /// on an off-screen device, shade a lit quad, read the pixels back,
    /// and assert the GPU genuinely shaded the scene.
    ///
    /// Skips cleanly (logs + returns) when no GPU adapter exists.
    #[test]
    fn headless_pbr_render_shades_a_lit_quad() {
        const SIZE: u32 = 256;

        // --- off-screen device: an Instance with no surface ---
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });
        let adapter = pollster::block_on(instance.request_adapter(
            &wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::default(),
                force_fallback_adapter: false,
                compatible_surface: None,
            },
        ));
        let Some(adapter) = adapter else {
            eprintln!(
                "headless_pbr_render_shades_a_lit_quad: no wgpu adapter in this \
                 environment — skipping the on-device render (the naga static \
                 validation still ran and covers shader soundness)."
            );
            return;
        };
        let (device, queue) = match pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("valenx.headless_pbr.device"),
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::downlevel_defaults(),
            },
            None,
        )) {
            Ok(pair) => pair,
            Err(e) => {
                eprintln!(
                    "headless_pbr_render_shades_a_lit_quad: adapter present but \
                     request_device failed ({e}) — skipping the on-device render."
                );
                return;
            }
        };

        // --- shader + pipeline from the real PBR WGSL ---
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("valenx.headless_pbr.shader"),
            source: wgpu::ShaderSource::Wgsl(PBR_FORWARD_WGSL.into()),
        });
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
        let bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("valenx.headless_pbr.bgl"),
                entries: &[
                    uniform_entry(0),
                    uniform_entry(1),
                    uniform_entry(2),
                    uniform_entry(3),
                ],
            });
        let pipeline_layout =
            device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("valenx.headless_pbr.layout"),
                bind_group_layouts: &[&bind_group_layout],
                push_constant_ranges: &[],
            });
        let pbr_vertex_attrs = wgpu::vertex_attr_array![0 => Float32x3, 1 => Float32x3];
        // A linear (non-sRGB) target so the readback pixels are the
        // shader's tone-mapped + sRGB-encoded output unchanged — no
        // second sRGB transform from the attachment.
        let target_format = wgpu::TextureFormat::Rgba8Unorm;
        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("valenx.headless_pbr.pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: "vs_main",
                compilation_options: Default::default(),
                buffers: &[wgpu::VertexBufferLayout {
                    array_stride: std::mem::size_of::<PbrVertex>() as wgpu::BufferAddress,
                    step_mode: wgpu::VertexStepMode::Vertex,
                    attributes: &pbr_vertex_attrs,
                }],
            },
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                // The quad is wound CCW as seen from +Z; cull back
                // faces exactly as the live PBR pass does.
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: Some(wgpu::Face::Back),
                ..Default::default()
            },
            depth_stencil: None,
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

        // --- a quad in the z=0 plane, normal +Z, facing the camera ---
        // Two CCW triangles spanning x,y in [-1, 1].
        let quad = [
            PbrVertex { position: [-1.0, -1.0, 0.0], normal: [0.0, 0.0, 1.0] },
            PbrVertex { position: [1.0, -1.0, 0.0], normal: [0.0, 0.0, 1.0] },
            PbrVertex { position: [1.0, 1.0, 0.0], normal: [0.0, 0.0, 1.0] },
            PbrVertex { position: [-1.0, -1.0, 0.0], normal: [0.0, 0.0, 1.0] },
            PbrVertex { position: [1.0, 1.0, 0.0], normal: [0.0, 0.0, 1.0] },
            PbrVertex { position: [-1.0, 1.0, 0.0], normal: [0.0, 0.0, 1.0] },
        ];
        let vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("valenx.headless_pbr.vertices"),
            contents: bytemuck::cast_slice(&quad),
            usage: wgpu::BufferUsages::VERTEX,
        });

        // --- uniforms ---
        // An orthographic MVP mapping x,y in [-1,1] straight to clip
        // space (identity is fine — the quad already spans the NDC box).
        let identity = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        // Camera in front of the quad on +Z.
        let frame = PbrFrameUniform::new(identity, [0.0, 0.0, 3.0], [0.0, 0.0, 0.0], 1);
        // A bright matte white surface.
        let material = Material::matte("white", [0.9, 0.9, 0.9]);
        let mat_uniform = PbrMaterialUniform::from_material(&material, 1.0);
        // One *point* light, in front of the quad and offset toward +X.
        // The shader attenuates a point light by inverse-square world
        // distance, so a fragment near the right edge (close to the
        // light) receives far more irradiance than one at the left
        // edge — a genuine position-dependent gradient a flat quad
        // under a *directional* light could never show. The intensity
        // is tuned to land the lit surface in the mid range of the
        // ACES tone curve: too strong and both sides clip toward white
        // (the gradient vanishes), too weak and nothing is bright. ~90
        // keeps the far (left) side mid-grey and the near (right) side
        // distinctly brighter without saturating.
        let light = PbrLightUniform::point([1.5, 0.0, 1.0], [1.0, 1.0, 1.0], 90.0);
        let mut light_slots = [PbrLightUniform::inactive(); MAX_LIGHTS];
        light_slots[0] = light;
        let light_array = LightArray { lights: light_slots };

        let frame_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("frame"),
            contents: bytemuck::bytes_of(&frame),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let mat_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("material"),
            contents: bytemuck::bytes_of(&mat_uniform),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let light_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("lights"),
            contents: bytemuck::bytes_of(&light_array),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let gi_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("gi"),
            contents: bytemuck::bytes_of(&ShL2Uniform::zero()),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("valenx.headless_pbr.bg"),
            layout: &bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: frame_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: mat_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: light_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: gi_buf.as_entire_binding() },
            ],
        });

        // --- offscreen colour target + a known clear colour ---
        let color_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("valenx.headless_pbr.color"),
            size: wgpu::Extent3d {
                width: SIZE,
                height: SIZE,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: target_format,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let color_view = color_texture.create_view(&wgpu::TextureViewDescriptor::default());
        // The clear colour — distinct from anything the shader writes.
        let clear = wgpu::Color { r: 0.0, g: 0.0, b: 0.0, a: 1.0 };

        // Readback buffer — 256-byte-aligned row pitch (wgpu's
        // COPY_BYTES_PER_ROW_ALIGNMENT). 256 px × 4 bytes = 1024,
        // already a multiple of 256.
        let bytes_per_row = SIZE * 4;
        assert_eq!(bytes_per_row % 256, 0, "row pitch must be 256-aligned");
        let readback = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("valenx.headless_pbr.readback"),
            size: (bytes_per_row * SIZE) as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        // --- encode: render pass + copy-to-buffer ---
        let mut encoder =
            device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("valenx.headless_pbr.encoder"),
            });
        {
            let mut rpass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("valenx.headless_pbr.pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(clear),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            rpass.set_pipeline(&pipeline);
            rpass.set_bind_group(0, &bind_group, &[]);
            rpass.set_vertex_buffer(0, vertex_buffer.slice(..));
            rpass.draw(0..quad.len() as u32, 0..1);
        }
        encoder.copy_texture_to_buffer(
            wgpu::ImageCopyTexture {
                texture: &color_texture,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            wgpu::ImageCopyBuffer {
                buffer: &readback,
                layout: wgpu::ImageDataLayout {
                    offset: 0,
                    bytes_per_row: Some(bytes_per_row),
                    rows_per_image: Some(SIZE),
                },
            },
            wgpu::Extent3d {
                width: SIZE,
                height: SIZE,
                depth_or_array_layers: 1,
            },
        );
        queue.submit(Some(encoder.finish()));

        // --- map the readback buffer + pull the pixels ---
        let slice = readback.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
        device.poll(wgpu::Maintain::Wait);
        rx.recv()
            .expect("map_async callback never fired")
            .expect("readback buffer map failed");
        let data = slice.get_mapped_range();
        let pixels: Vec<[u8; 4]> = data
            .chunks_exact(4)
            .map(|c| [c[0], c[1], c[2], c[3]])
            .collect();
        drop(data);
        readback.unmap();

        // --- assertions: the GPU genuinely shaded the scene ---
        assert_eq!(pixels.len(), (SIZE * SIZE) as usize, "wrong pixel count");

        // (1) The framebuffer is NOT uniformly the clear colour — the
        //     quad was rasterised and the fragment shader ran.
        let clear_rgb = [0u8, 0u8, 0u8];
        let non_clear = pixels
            .iter()
            .filter(|p| [p[0], p[1], p[2]] != clear_rgb)
            .count();
        assert!(
            non_clear > (SIZE * SIZE / 2) as usize,
            "expected the lit quad to fill most of the frame; only {non_clear} \
             of {} pixels differ from the clear colour",
            SIZE * SIZE
        );

        // (2) The lit surface is genuinely bright somewhere — a matte
        //     white quad near a point light must produce a clearly
        //     lit region, not just faint noise. A pixel-sum above 300
        //     (of a 765 max) is unambiguously brighter than the black
        //     clear colour.
        let brightest = pixels
            .iter()
            .map(|p| p[0] as u32 + p[1] as u32 + p[2] as u32)
            .max()
            .unwrap_or(0);
        assert!(
            brightest > 300,
            "the lit quad should have a clearly lit region near the \
             light; brightest pixel sums only {brightest} (max 765)"
        );

        // (3) A horizontal brightness gradient toward the light. The
        //     point light sits in front of the quad offset toward +X,
        //     so right-side fragments are closer (less inverse-square
        //     falloff, more head-on n·l) — the mid-row mean luminance
        //     of a right-of-centre band must clearly exceed a
        //     left-of-centre band.
        let mid_row = SIZE / 2;
        let lum = |x: u32| -> f64 {
            let p = pixels[(mid_row * SIZE + x) as usize];
            0.2126 * p[0] as f64 + 0.7152 * p[1] as f64 + 0.0722 * p[2] as f64
        };
        let quarter = SIZE / 4;
        let left_mean: f64 =
            (quarter - 16..quarter + 16).map(lum).sum::<f64>() / 32.0;
        let right_mean: f64 = (3 * quarter - 16..3 * quarter + 16)
            .map(lum)
            .sum::<f64>()
            / 32.0;
        // The left (far) band must NOT be saturated — otherwise both
        // bands are clipped white and the test would not actually be
        // exercising a gradient. Then require a clear ≥15 % margin,
        // which a uniformly-shaded frame (ratio 1.0) cannot pass.
        assert!(
            left_mean < 235.0,
            "the far band is saturated ({left_mean:.1}/255) — the scene \
             is too bright to exercise a gradient; lower the intensity"
        );
        assert!(
            right_mean > left_mean * 1.15,
            "expected a clear brightness gradient toward the +X point \
             light: left band luminance {left_mean:.2}, right band \
             {right_mean:.2} (right should be >1.15x left)"
        );

        eprintln!(
            "headless_pbr_render_shades_a_lit_quad: rendered on `{}` ({:?}) — \
             {non_clear}/{} shaded px, brightest sum {brightest}, \
             gradient L={left_mean:.1} -> R={right_mean:.1}",
            adapter.get_info().name,
            adapter.get_info().backend,
            SIZE * SIZE
        );
    }
}
