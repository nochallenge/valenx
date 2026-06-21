//! The central viewport panel.
//!
//! Two render styles ship today, selectable through the View menu or
//! the command palette:
//!
//! - **Shaded** (default): flat-shaded filled triangles with
//!   back-face culling and a dot-product light model, painted via
//!   `egui::Mesh` + painter's-algorithm depth sort.
//! - **Wireframe**: three line segments per triangle, no culling.
//!
//! Loaded canonical meshes built from surface elements (Tri3 / Quad4 —
//! the procedural CAD models such as springs, rockets and engines)
//! render through the same shaded path as STLs. Volumetric meshes, and
//! any mesh carrying a field overlay, keep the element-wireframe view.
//!
//! Both paths share `valenx-viz::projection`, so when the real
//! `wgpu` render pass arrives it only has to replace the rasteriser;
//! the camera and projection math are already test-covered.
//!
//! Mouse interaction:
//!
//! - Left-drag or middle-drag — orbit.
//! - Shift + drag — pan (Phase 1 tail; no-op today).
//! - Scroll — zoom.
//! - Double-click — reframe on the mesh bounding box.

use eframe::egui;
use eframe::egui_wgpu;
use valenx_mesh::{ElementType, Mesh};
use valenx_viz::{
    project_point, project_triangle, OrbitCamera, ScreenPoint, StlTriangle, TriangleMesh,
};

use crate::wgpu_renderer::{
    inv_mvp_from_camera, mvp_from_camera, triangles_to_vertices, WgpuRenderer,
};
use crate::{LoadedMesh, LoadedStl};

/// Which renderer state the viewport draws in.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ShadingMode {
    /// Per-vertex lighting on filled triangles (default).
    #[default]
    Shaded,
    /// Edge-only line rendering, ignoring face normals.
    Wireframe,
}

/// Borrowed bundle of everything the viewport needs to render one
/// frame — handed to the viewport render entry point each frame so
/// the app's own state stays cheap to clone.
pub struct ViewportState<'a> {
    /// Mutable orbit camera (the panel may dolly / orbit it).
    pub camera: &'a mut OrbitCamera,
    /// Currently-loaded STL, if any.
    pub stl: Option<&'a LoadedStl>,
    /// Currently-loaded canonical mesh, if any.
    pub mesh: Option<&'a LoadedMesh>,
    /// Shading mode to draw the geometry in.
    pub shading: ShadingMode,
    /// Optional wgpu render context (the renderer falls back to
    /// software path when this is `None`).
    pub wgpu: Option<WgpuCtx<'a>>,
    /// Optional per-node scalar overlay. When `Some`, the mesh
    /// wireframe colour-codes each edge by the average of its two
    /// endpoints' field values via the `cool_to_warm` ramp. `None`
    /// keeps the default flat-blue wireframe colour.
    pub field_overlay: Option<&'a valenx_fields::Field>,
    /// When `Some((point, normal))`, draw the cross-section where
    /// that plane slices the loaded geometry. Driven by the Mesh
    /// Toolbox's "Show cut overlay" checkbox.
    pub cut_overlay: Option<([f64; 3], [f64; 3])>,
    /// When `Some(...)`, draw the active sketch as a 2D line / circle
    /// / arc overlay on top of the geometry. Driven by the Sketcher
    /// panel's "Show sketch overlay" checkbox.
    pub sketch_overlay: Option<crate::sketch_overlay::SketchOverlayState<'a>>,
    /// When `Some(...)`, draw the active draft document (2D
    /// construction lines, dimensions, text) on top of the geometry.
    /// Driven by the Draft panel's "Show draft overlay" checkbox.
    pub draft_overlay: Option<crate::draft_overlay::DraftOverlayState<'a>>,
    /// When `Some(...)`, draw the last-generated CAM toolpath as
    /// coloured polylines (rapid/cut/plunge) on top of the geometry.
    /// Driven by the CAM panel's "Toggle Simulate Overlay" button.
    pub cam_overlay: Option<&'a valenx_cam::Toolpath>,
    /// Whether the cursor snaps to the ground grid: snaps the live
    /// cursor coordinate to the nearest grid node and draws a marker
    /// there (Fusion-style). Toggled from the View menu (defaults on).
    pub snap_to_grid: bool,
}

/// Per-frame wgpu bundle. `None` means eframe was built without the
/// wgpu backend (or the render state wasn't available at startup);
/// the viewport falls back to the painter-only shaded path.
pub struct WgpuCtx<'a> {
    pub renderer: &'a mut WgpuRenderer,
    pub render_state: &'a egui_wgpu::RenderState,
    pub pixels_per_point: f32,
}

/// Render the viewport into the current `ui`'s available rect.
pub fn show(ui: &mut egui::Ui, mut state: ViewportState<'_>) {
    let available = ui.available_size();
    let (rect, response) = ui.allocate_exact_size(
        egui::vec2(available.x.max(32.0), available.y.max(32.0)),
        egui::Sense::click_and_drag(),
    );

    // Orbit on primary or middle drag.
    if response.dragged_by(egui::PointerButton::Primary)
        || response.dragged_by(egui::PointerButton::Middle)
    {
        let delta = response.drag_delta();
        state.camera.orbit(delta.x * 0.5, -delta.y * 0.5);
    }

    if response.hovered() {
        let scroll = ui.input(|i| i.raw_scroll_delta.y);
        if scroll.abs() > f32::EPSILON {
            state.camera.zoom(scroll * 0.01);
        }
    }

    if response.double_clicked() {
        // Prefer the canonical mesh bounds if one's loaded; STL
        // fallback otherwise.
        if let Some(m) = state.mesh {
            if let Some((min, max)) = mesh_aabb(&m.mesh) {
                state.camera.frame_bounds(min, max);
            }
        } else if let Some(stl) = state.stl {
            if let Some((min, max)) = stl.mesh.bounding_box() {
                state.camera.frame_bounds(min, max);
            }
        }
    }

    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_gray(42));
    painter.rect_stroke(
        rect,
        0.0,
        egui::Stroke::new(1.0, egui::Color32::from_gray(28)),
    );

    // Header strip — prefer the canonical mesh when loaded, else
    // the STL, else a hint.
    let header_text = match (state.mesh, state.stl) {
        (Some(m), _) => format!(
            "Viewport — {} ({} nodes · {} elements) · drag to orbit · scroll to zoom · F to frame",
            m.path.file_name().unwrap_or_default().to_string_lossy(),
            m.mesh.stats.node_count,
            m.mesh.stats.element_count,
        ),
        (None, Some(stl)) => format!(
            "Viewport — {} ({} triangles) · {} · drag to orbit · scroll to zoom · F to frame",
            stl.path.file_name().unwrap_or_default().to_string_lossy(),
            stl.mesh.triangle_count(),
            match state.shading {
                ShadingMode::Shaded => "shaded",
                ShadingMode::Wireframe => "wireframe",
            }
        ),
        (None, None) => "Viewport — no geometry loaded".to_string(),
    };
    painter.text(
        rect.min + egui::vec2(8.0, 4.0),
        egui::Align2::LEFT_TOP,
        header_text,
        egui::FontId::proportional(12.0),
        egui::Color32::from_gray(190),
    );

    // Triangles: draw the STL first, then overlay the mesh (if any)
    // so element edges sit on top of the surface shell. Empty case
    // shows a help hint.
    let anything_loaded = state.stl.is_some() || state.mesh.is_some();

    // GPU scene background: the ground grid + world axes render ALWAYS
    // when the wgpu backend is available — even with no geometry — into a
    // depth-buffered offscreen pass. In Shaded mode with an STL the same
    // pass also draws the shaded mesh, so geometry correctly occludes the
    // grid. Empty slice => grid only.
    // The shaded wgpu pass draws the ground grid + axes ALWAYS, plus an
    // optional shaded surface: the STL when one's loaded (Shaded mode),
    // otherwise a canonical-mesh surface (Shaded mode, no field overlay)
    // so loaded meshes — springs, rockets, engines — render as solids
    // instead of bare wireframe. The mesh surface is owned so its
    // vertices outlive the draw call.
    let mesh_surface = match (state.mesh, state.shading, state.field_overlay) {
        (Some(m), ShadingMode::Shaded, None) if state.stl.is_none() => {
            mesh_to_triangle_surface(&m.mesh)
        }
        _ => None,
    };
    let wgpu_drew = if let Some(ctx) = state.wgpu.as_mut() {
        if let Some(stl) = state
            .stl
            .filter(|_| matches!(state.shading, ShadingMode::Shaded))
        {
            draw_shaded_wgpu(&painter, rect, state.camera, stl, ctx)
        } else if let Some(surf) = mesh_surface.as_ref() {
            let vertices = triangles_to_vertices(surf);
            render_wgpu_scene(&painter, rect, state.camera, ctx, &vertices)
        } else {
            render_wgpu_scene(&painter, rect, state.camera, ctx, &[])
        }
    } else {
        false
    };
    // True when the GPU drew the mesh (not STL) surface as a shaded
    // solid: the mesh block below then skips its wireframe so the solid
    // reads clean.
    let mesh_shaded = mesh_surface.is_some() && wgpu_drew;

    if let Some(stl) = state.stl {
        match state.shading {
            ShadingMode::Shaded => {
                // wgpu drew the shaded mesh above; otherwise fall back to
                // the painter-only Lambert (eframe without a wgpu backend).
                if !wgpu_drew {
                    draw_shaded(ui, &painter, rect, state.camera, stl);
                }
            }
            ShadingMode::Wireframe => draw_wireframe(&painter, rect, state.camera, stl),
        }
        draw_bbox(&painter, rect, state.camera, stl);
    }
    if let Some(mesh) = state.mesh {
        // Filled-triangle field rendering: when an OnNode scalar
        // overlay is present, paint surface triangles with per-vertex
        // colours from the colormap (Gouraud-style smooth shading via
        // egui's per-vertex tinting). Wireframe lands on top for edge
        // definition.
        if let Some(field) = state.field_overlay {
            draw_mesh_filled_field(&painter, rect, state.camera, mesh, field);
            draw_mesh_wireframe(&painter, rect, state.camera, mesh, state.field_overlay);
        } else if !mesh_shaded {
            // Wireframe mode, no wgpu backend, or a mesh with no
            // extractable surface (volumetric elements): fall back to the
            // element wireframe. When the shaded solid drew, leave it clean.
            draw_mesh_wireframe(&painter, rect, state.camera, mesh, None);
        }
    }
    // Cut-plane cross-section overlay. Drawn after the geometry so
    // the bright slice line sits on top of both the shaded shell and
    // the wireframe. Only runs when the Mesh Toolbox's "Show cut
    // overlay" checkbox is ticked — zero cost otherwise.
    if let Some((point, normal)) = state.cut_overlay {
        use nalgebra::Vector3;
        let pt = Vector3::new(point[0], point[1], point[2]);
        let nm = Vector3::new(normal[0], normal[1], normal[2]);
        // Compute the cross-section against whatever geometry is
        // loaded — prefer the canonical mesh (the cut math is native
        // to it), fall back to the STL triangle soup.
        let segments: Vec<_> = if let Some(m) = state.mesh {
            valenx_mesh::cut::intersect_plane(&m.mesh, pt, nm)
        } else if let Some(stl) = state.stl {
            // STL path: feed the loose triangles straight through the
            // triangle-soup intersection variant — no transient Mesh
            // allocation per frame.
            let tris: Vec<[[f64; 3]; 3]> = stl
                .mesh
                .triangles
                .iter()
                .map(|t| {
                    [
                        [
                            t.vertices[0][0] as f64,
                            t.vertices[0][1] as f64,
                            t.vertices[0][2] as f64,
                        ],
                        [
                            t.vertices[1][0] as f64,
                            t.vertices[1][1] as f64,
                            t.vertices[1][2] as f64,
                        ],
                        [
                            t.vertices[2][0] as f64,
                            t.vertices[2][1] as f64,
                            t.vertices[2][2] as f64,
                        ],
                    ]
                })
                .collect();
            valenx_mesh::cut::intersect_plane_triangles(&tris, pt, nm)
        } else {
            Vec::new()
        };
        if !segments.is_empty() {
            draw_cut_overlay(&painter, rect, state.camera, &segments);
        }
    }

    // Sketch overlay (Phase 1H). Drawn after the cut overlay so the
    // sketch lines sit on top of everything — sketches are an
    // interactive editing surface and need to read clearly. Disabled
    // by default; the Sketcher panel's "Show sketch overlay" checkbox
    // forwards the state here when on. Reuses `state.camera` — the
    // overlay borrows it immutably, so this section runs after every
    // camera-mutating interaction above.
    if let Some(overlay) = state.sketch_overlay {
        crate::sketch_overlay::draw(&painter, rect, state.camera, overlay);
    }

    // Draft overlay (Phase 4E). Same layering rationale as the sketch
    // overlay — sits on top of the cut overlay and the underlying
    // geometry so 2D construction lines and dimensions read clearly.
    if let Some(overlay) = state.draft_overlay {
        crate::draft_overlay::draw(&painter, rect, state.camera, overlay);
    }

    // CAM toolpath overlay (Phase 10E). Drawn last so the coloured
    // polylines (rapid/cut/plunge) sit on top of everything else.
    if let Some(tp) = state.cam_overlay {
        crate::cam_overlay::draw(&painter, rect, state.camera, tp);
    }

    // Colour-bar legend in the bottom-right when an overlay is
    // active. Without this users see colours but can't read the
    // field name or the value range — the colours are decorative
    // rather than informative.
    if let Some(field) = state.field_overlay {
        draw_field_legend(&painter, rect, field);
    }
    if !anything_loaded {
        // In-project empty state. The no-project / no-mesh state shows
        // the welcome landing page instead — by the time `viewport::show`
        // runs we know the user has at least a project loaded but
        // hasn't picked geometry yet. Keep the hint quiet and pointed
        // at what they probably want next: run a case (which generates
        // mesh + results), or drop an STL to inspect.
        // The grid + axes already read as a "3D workspace ready for
        // geometry", so keep this to one quiet hint pinned to the
        // bottom-centre — out of the way of the axes that cross at the
        // viewport centre (a big centred label there looked unfinished).
        painter.text(
            egui::pos2(rect.center().x, rect.bottom() - 18.0),
            egui::Align2::CENTER_BOTTOM,
            "Run a case · drop a .stl · or insert a primitive from the Part panel",
            egui::FontId::proportional(12.5),
            egui::Color32::from_gray(115),
        );
    }

    // Camera readout
    painter.text(
        rect.left_bottom() + egui::vec2(8.0, -22.0),
        egui::Align2::LEFT_BOTTOM,
        format!(
            "camera · az {:>6.1}° · el {:>6.1}° · d {:>7.3}",
            state.camera.azimuth_deg, state.camera.elevation_deg, state.camera.distance
        ),
        egui::FontId::monospace(11.0),
        egui::Color32::from_gray(150),
    );

    // Orientation gizmo + live cursor-coordinate HUD, painted on top of
    // everything. A click on a gizmo axis snaps the camera to that view.
    let grid_spacing = valenx_viz::grid_lod_params(state.camera.distance).0;
    let gizmo_hover = crate::scene_overlay::draw_overlay(
        &painter,
        rect,
        state.camera,
        response.hover_pos(),
        grid_spacing,
        state.snap_to_grid,
    );
    if response.clicked() {
        if let Some(face) = gizmo_hover {
            state.camera.set_view(valenx_viz::gizmo_view_for_face(face));
        }
    }
}

// ---------------------------------------------------------------------------
// Shaded via wgpu (real depth-buffered 3D)
// ---------------------------------------------------------------------------

fn draw_shaded_wgpu(
    painter: &egui::Painter,
    rect: egui::Rect,
    camera: &OrbitCamera,
    stl: &LoadedStl,
    ctx: &mut WgpuCtx<'_>,
) -> bool {
    let vertices = triangles_to_vertices(&stl.mesh);
    render_wgpu_scene(painter, rect, camera, ctx, &vertices)
}

/// Render the offscreen wgpu scene — the ground grid + axes **always**,
/// plus the optional shaded mesh `vertices` — and blit it into `rect`.
/// Pass an empty slice to draw just the grid (e.g. the empty viewport).
fn render_wgpu_scene(
    painter: &egui::Painter,
    rect: egui::Rect,
    camera: &OrbitCamera,
    ctx: &mut WgpuCtx<'_>,
    vertices: &[crate::wgpu_renderer::Vertex],
) -> bool {
    let w_logical = rect.width();
    let h_logical = rect.height();
    if w_logical < 1.0 || h_logical < 1.0 {
        return false;
    }
    let ppp = ctx.pixels_per_point.max(0.5);
    let size_px = [(w_logical * ppp) as u32, (h_logical * ppp) as u32];

    let mvp = mvp_from_camera(camera, w_logical, h_logical);
    let inv_mvp = inv_mvp_from_camera(camera, w_logical, h_logical);
    let light_dir = [-0.3f32, -0.8, -0.5];
    let eye = camera.eye();
    let cam_pos = [eye.x, eye.y, eye.z];
    let (minor, blend_t) = valenx_viz::grid_lod_params(camera.distance);
    // x = minor spacing, yz unused, w = fade distance.
    let grid = [minor, 0.0, 0.0, camera.distance * 14.0];
    // x = LOD blend_t (0 = minor fully visible, 1 = minor faded out).
    let grid2 = [blend_t, 0.0, 0.0, 0.0];

    let mut render_guard = ctx.render_state.renderer.write();
    let texture_id = match ctx.renderer.render(
        &mut render_guard,
        size_px,
        mvp,
        inv_mvp,
        light_dir,
        cam_pos,
        grid,
        grid2,
        vertices,
    ) {
        Some(id) => id,
        None => return false,
    };
    drop(render_guard);

    // Blit the offscreen texture into the viewport rect via egui's
    // 2D image shape. uv = full texture.
    let uv = egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0));
    painter.image(texture_id, rect, uv, egui::Color32::WHITE);
    true
}

/// Build a renderable triangle **surface** from a canonical mesh's
/// surface elements — Tri3 directly, Quad4 split into two triangles —
/// with per-face normals, so loaded meshes can be drawn as shaded solids
/// through the same wgpu pass as STLs. Volumetric element types (Tet4,
/// Hex8, ...) are skipped: those keep the wireframe path. `None` when no
/// surface triangle could be produced.
///
/// `pub(crate)` so the dockable `workspace:<n>` tile can build the same
/// shaded surface for its own per-tile 3-D view (see
/// [`crate::dock_layout`]'s `render_workspace_body`).
pub(crate) fn mesh_to_triangle_surface(mesh: &Mesh) -> Option<TriangleMesh> {
    let nodes = &mesh.nodes;
    let tri = |a: u32, b: u32, c: u32| -> Option<StlTriangle> {
        let p = nodes.get(a as usize)?;
        let q = nodes.get(b as usize)?;
        let r = nodes.get(c as usize)?;
        let v = |x: &nalgebra::Vector3<f64>| [x.x as f32, x.y as f32, x.z as f32];
        let mut t = StlTriangle {
            normal: [0.0, 0.0, 1.0],
            vertices: [v(p), v(q), v(r)],
        };
        t.normal = t.computed_normal();
        Some(t)
    };
    let mut tris: Vec<StlTriangle> = Vec::new();
    for block in &mesh.element_blocks {
        let conn = &block.connectivity;
        match block.element_type {
            ElementType::Tri3 => {
                for f in conn.chunks_exact(3) {
                    tris.extend(tri(f[0], f[1], f[2]));
                }
            }
            ElementType::Quad4 => {
                for q in conn.chunks_exact(4) {
                    tris.extend(tri(q[0], q[1], q[2]));
                    tris.extend(tri(q[0], q[2], q[3]));
                }
            }
            _ => {}
        }
    }
    if tris.is_empty() {
        None
    } else {
        Some(TriangleMesh {
            format: None,
            name: None,
            triangles: tris,
        })
    }
}

// ---------------------------------------------------------------------------
// Shaded (painter fallback — filled triangles + per-face light)
// ---------------------------------------------------------------------------

fn draw_shaded(
    ui: &egui::Ui,
    painter: &egui::Painter,
    rect: egui::Rect,
    camera: &OrbitCamera,
    stl: &LoadedStl,
) {
    let width = rect.width();
    let height = rect.height();
    let origin = rect.min;
    let eye = camera.eye();

    // Over-the-shoulder key light + soft fill. Normalised once.
    let light = normalize([-0.3f32, -0.8, -0.5]);

    #[derive(Clone, Copy)]
    struct Entry {
        depth: f32,
        pts: [ScreenPoint; 3],
        color: egui::Color32,
    }

    let mut entries: Vec<Entry> = Vec::with_capacity(stl.mesh.triangle_count());
    for tri in &stl.mesh.triangles {
        // Compute the normal from the triangle's winding rather than
        // trusting the file — many exporters ship zeros.
        let normal = tri.computed_normal();

        // Back-face cull against the view direction. For an orbit
        // camera we take the average centroid's direction from the
        // eye as the view vector.
        let centroid = tri_centroid(tri);
        let view = normalize([
            centroid[0] - eye.x,
            centroid[1] - eye.y,
            centroid[2] - eye.z,
        ]);
        if dot(view, normal) > 0.0 {
            continue;
        }

        if let Some(pts) = project_triangle(camera, width, height, &tri.vertices) {
            let depth = (pts[0].depth + pts[1].depth + pts[2].depth) / 3.0;
            let lambert = (-dot(light, normal)).clamp(0.0, 1.0);
            let intensity = 0.22 + 0.78 * lambert;
            let color = shade_color(intensity);
            entries.push(Entry { depth, pts, color });
        }
    }

    // Painter's-algorithm: far first.
    entries.sort_by(|a, b| {
        b.depth
            .partial_cmp(&a.depth)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Build one Mesh per ~1024 triangles — egui's tessellator handles
    // batches fine, but we cap so a pathological input doesn't OOM.
    let mut mesh = egui::Mesh::default();
    let mut count = 0usize;
    for e in &entries {
        push_triangle(&mut mesh, origin, e.pts, e.color);
        count += 1;
        if count.is_multiple_of(1024) {
            painter.add(egui::Shape::mesh(std::mem::take(&mut mesh)));
        }
    }
    if !mesh.is_empty() {
        painter.add(egui::Shape::mesh(mesh));
    }

    // Keep ui around for future interactivity hooks; currently only
    // used for input reads at the caller. Binding it avoids a
    // `dead_code`-style warning when this function stays narrow.
    let _ = ui;
}

fn push_triangle(
    mesh: &mut egui::Mesh,
    origin: egui::Pos2,
    pts: [ScreenPoint; 3],
    color: egui::Color32,
) {
    let base = mesh.vertices.len() as u32;
    for p in &pts {
        mesh.vertices.push(egui::epaint::Vertex {
            pos: origin + egui::vec2(p.x, p.y),
            uv: egui::epaint::WHITE_UV,
            color,
        });
    }
    mesh.add_triangle(base, base + 1, base + 2);
}

fn tri_centroid(tri: &StlTriangle) -> [f32; 3] {
    [
        (tri.vertices[0][0] + tri.vertices[1][0] + tri.vertices[2][0]) / 3.0,
        (tri.vertices[0][1] + tri.vertices[1][1] + tri.vertices[2][1]) / 3.0,
        (tri.vertices[0][2] + tri.vertices[1][2] + tri.vertices[2][2]) / 3.0,
    ]
}

fn normalize(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-9 {
        [0.0, 0.0, 1.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

fn dot(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn shade_color(intensity: f32) -> egui::Color32 {
    // Neutral brushed-metal base so CAD models read naturally.
    // Intensity 0..1 scales across a tonal band rather than to pure
    // black — pure black on a dark viewport looks like a hole.
    let base = [170.0f32, 195.0, 225.0];
    let r = (base[0] * intensity).clamp(28.0, 255.0);
    let g = (base[1] * intensity).clamp(28.0, 255.0);
    let b = (base[2] * intensity).clamp(28.0, 255.0);
    egui::Color32::from_rgb(r as u8, g as u8, b as u8)
}

// ---------------------------------------------------------------------------
// Wireframe
// ---------------------------------------------------------------------------

fn draw_wireframe(
    painter: &egui::Painter,
    rect: egui::Rect,
    camera: &OrbitCamera,
    stl: &LoadedStl,
) {
    let width = rect.width();
    let height = rect.height();
    let origin = rect.min;
    let stroke = egui::Stroke::new(0.8, egui::Color32::from_rgb(150, 190, 230));

    let mut projected: Vec<(f32, [ScreenPoint; 3])> = Vec::with_capacity(stl.mesh.triangle_count());
    for tri in &stl.mesh.triangles {
        if let Some(pts) = project_triangle(camera, width, height, &tri.vertices) {
            let avg_depth = (pts[0].depth + pts[1].depth + pts[2].depth) / 3.0;
            projected.push((avg_depth, pts));
        }
    }
    projected.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

    for (_d, pts) in &projected {
        let p0 = origin + egui::vec2(pts[0].x, pts[0].y);
        let p1 = origin + egui::vec2(pts[1].x, pts[1].y);
        let p2 = origin + egui::vec2(pts[2].x, pts[2].y);
        painter.line_segment([p0, p1], stroke);
        painter.line_segment([p1, p2], stroke);
        painter.line_segment([p2, p0], stroke);
    }
}

fn draw_bbox(painter: &egui::Painter, rect: egui::Rect, camera: &OrbitCamera, stl: &LoadedStl) {
    let (min, max) = match stl.mesh.bounding_box() {
        Some(bb) => bb,
        None => return,
    };
    let stroke = egui::Stroke::new(1.0, egui::Color32::from_rgb(220, 180, 120));
    let origin = rect.min;
    let w = rect.width();
    let h = rect.height();

    let corners = [
        [min[0], min[1], min[2]],
        [max[0], min[1], min[2]],
        [max[0], max[1], min[2]],
        [min[0], max[1], min[2]],
        [min[0], min[1], max[2]],
        [max[0], min[1], max[2]],
        [max[0], max[1], max[2]],
        [min[0], max[1], max[2]],
    ];
    let mut proj: [Option<ScreenPoint>; 8] = Default::default();
    for (i, c) in corners.iter().enumerate() {
        proj[i] = project_point(camera, w, h, *c);
    }
    const EDGES: [(usize, usize); 12] = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ];
    for (a, b) in EDGES {
        if let (Some(pa), Some(pb)) = (proj[a], proj[b]) {
            let p0 = origin + egui::vec2(pa.x, pa.y);
            let p1 = origin + egui::vec2(pb.x, pb.y);
            painter.line_segment([p0, p1], stroke);
        }
    }
}

/// Draw the cross-section edges where the cut plane intersects the
/// loaded geometry. Bright cyan so it reads clearly over both the
/// shaded shell and the wireframe. Mirrors `draw_bbox`'s
/// project-then-line_segment pattern.
fn draw_cut_overlay(
    painter: &egui::Painter,
    rect: egui::Rect,
    camera: &OrbitCamera,
    segments: &[valenx_mesh::cut::LineSegment],
) {
    let stroke = egui::Stroke::new(2.0, egui::Color32::from_rgb(80, 230, 230));
    let origin = rect.min;
    let w = rect.width();
    let h = rect.height();
    for seg in segments {
        let a = [seg.a.x as f32, seg.a.y as f32, seg.a.z as f32];
        let b = [seg.b.x as f32, seg.b.y as f32, seg.b.z as f32];
        if let (Some(pa), Some(pb)) = (
            project_point(camera, w, h, a),
            project_point(camera, w, h, b),
        ) {
            painter.line_segment(
                [
                    origin + egui::vec2(pa.x, pa.y),
                    origin + egui::vec2(pb.x, pb.y),
                ],
                stroke,
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Canonical-mesh wireframe overlay
// ---------------------------------------------------------------------------

/// Draw every element's face edges as line segments, projected into
/// the viewport rect via the shared `OrbitCamera` math. Intended as
/// a "what did the mesher actually produce" preview — pairs nicely
/// with the shaded STL shell underneath, but also usable standalone.
fn draw_mesh_wireframe(
    painter: &egui::Painter,
    rect: egui::Rect,
    camera: &OrbitCamera,
    loaded: &LoadedMesh,
    field_overlay: Option<&valenx_fields::Field>,
) {
    if loaded.mesh.nodes.is_empty() {
        return;
    }
    let w = rect.width();
    let h = rect.height();
    let origin = rect.min;
    // Default stroke when no overlay is requested. With overlay we
    // build a fresh per-edge stroke each iteration so we can paint
    // the colour from the colormap.
    let default_stroke = egui::Stroke::new(0.8, egui::Color32::from_rgb(180, 210, 240));

    // Pre-project all nodes once.
    let projected: Vec<Option<ScreenPoint>> = loaded
        .mesh
        .nodes
        .iter()
        .map(|n| project_point(camera, w, h, [n.x as f32, n.y as f32, n.z as f32]))
        .collect();

    // Resolve the field overlay into a per-node-index value lookup +
    // cached range. We only support scalar OnNode fields here; vector
    // / tensor / per-cell would need surface rendering to make sense.
    let scalar_lookup = field_overlay
        .filter(|f| {
            matches!(f.kind, valenx_fields::FieldKind::Scalar)
                && matches!(f.location, valenx_fields::Location::OnNode)
                && f.data.len() == loaded.mesh.nodes.len()
        })
        .map(|f| {
            let (min, max) = f.range.unwrap_or_else(|| field_min_max(&f.data));
            (f.data.as_slice(), min, max)
        });

    let line = |a: usize, b: usize| -> Option<(egui::Pos2, egui::Pos2)> {
        let pa = projected.get(a)?.as_ref()?;
        let pb = projected.get(b)?.as_ref()?;
        Some((
            origin + egui::vec2(pa.x, pa.y),
            origin + egui::vec2(pb.x, pb.y),
        ))
    };

    for block in &loaded.mesh.element_blocks {
        let npe = block.element_type.nodes_per_element();
        if npe == 0 {
            continue;
        }
        let element_count = block.connectivity.len() / npe;
        for i in 0..element_count {
            let start = i * npe;
            let nodes = &block.connectivity[start..start + npe];
            for &(a_rel, b_rel) in edges_for(block.element_type) {
                let a = nodes[a_rel] as usize;
                let b = nodes[b_rel] as usize;
                if let Some((pa, pb)) = line(a, b) {
                    let stroke = match &scalar_lookup {
                        Some((data, min, max)) => {
                            // Average the two endpoints' field values
                            // for the edge colour. Slightly thicker
                            // line than the default so the overlay
                            // reads as "the data layer" rather than
                            // bleeding into the geometry layer.
                            let va = data[a];
                            let vb = data[b];
                            let avg = 0.5 * (va + vb);
                            let [r, g, b_] =
                                valenx_fields::colormap::cool_to_warm_in_range(avg, *min, *max);
                            egui::Stroke::new(1.2, egui::Color32::from_rgb(r, g, b_))
                        }
                        None => default_stroke,
                    };
                    painter.line_segment([pa, pb], stroke);
                }
            }
        }
    }
}

/// Filled-triangle field overlay. Walks every element block, expands
/// each element into surface triangles via [`triangles_for`], projects
/// the vertices into screen space, samples the colormap from the
/// scalar field (per-vertex for OnNode, per-cell flat for OnCell),
/// and paints them via `egui::Mesh`.
///
/// Supported field kinds:
/// - **OnNode + Scalar** with `data.len() == nodes.len()` —
///   Gouraud-style smooth shading via per-vertex colours.
/// - **OnCell + Scalar** with `data.len() == sum(cells across all
///   blocks)` — flat-fill: every triangle of an element gets the
///   same colour from the cell's value.
///
/// Vector / tensor fields fall through silently. The wireframe
/// overlay still colour-codes their edges (per-edge averaging) so
/// the user gets some visual feedback.
///
/// Painted before the wireframe so the wireframe lands on top for
/// edge definition.
fn draw_mesh_filled_field(
    painter: &egui::Painter,
    rect: egui::Rect,
    camera: &OrbitCamera,
    loaded: &LoadedMesh,
    field: &valenx_fields::Field,
) {
    if loaded.mesh.nodes.is_empty() {
        return;
    }
    if !matches!(field.kind, valenx_fields::FieldKind::Scalar) {
        return;
    }
    // Validate the data length matches the location BEFORE doing any
    // projection work so a wrong-shaped field doesn't half-render.
    let location_ok = match field.location {
        valenx_fields::Location::OnNode => field.data.len() == loaded.mesh.nodes.len(),
        valenx_fields::Location::OnCell => field.data.len() == total_cell_count(&loaded.mesh),
        _ => false,
    };
    if !location_ok {
        return;
    }

    let w = rect.width();
    let h = rect.height();
    let origin = rect.min;
    let (min_v, max_v) = field.range.unwrap_or_else(|| field_min_max(&field.data));

    let projected: Vec<Option<ScreenPoint>> = loaded
        .mesh
        .nodes
        .iter()
        .map(|n| project_point(camera, w, h, [n.x as f32, n.y as f32, n.z as f32]))
        .collect();

    let value_to_color = |v: f64| -> egui::Color32 {
        let [r, g, b] = valenx_fields::colormap::cool_to_warm_in_range(v, min_v, max_v);
        // Slight alpha so the wireframe overlay reads as "the data
        // layer" without completely hiding the shaded geometry below.
        egui::Color32::from_rgba_unmultiplied(r, g, b, 200)
    };

    // For OnCell we need a running cell index across element blocks
    // so the field.data lookup matches the canonical iteration order
    // (block 0 cells first, then block 1, etc.).
    let mut global_cell_offset = 0usize;
    for block in &loaded.mesh.element_blocks {
        let npe = block.element_type.nodes_per_element();
        if npe == 0 {
            continue;
        }
        let triangles = triangles_for(block.element_type);
        if triangles.is_empty() {
            global_cell_offset += block.connectivity.len() / npe;
            continue;
        }
        let element_count = block.connectivity.len() / npe;
        let mut egui_mesh = egui::Mesh::default();
        for i in 0..element_count {
            let start = i * npe;
            let nodes = &block.connectivity[start..start + npe];
            // OnCell: one colour per element, applied to every
            // triangle of that element. OnNode: per-vertex lookup
            // inside the inner loop.
            let cell_color = match field.location {
                valenx_fields::Location::OnCell => {
                    Some(value_to_color(field.data[global_cell_offset + i]))
                }
                _ => None,
            };
            for &(a_rel, b_rel, c_rel) in triangles {
                let a = nodes[a_rel] as usize;
                let b = nodes[b_rel] as usize;
                let c = nodes[c_rel] as usize;
                let pa = match projected.get(a).and_then(|p| p.as_ref()) {
                    Some(p) => p,
                    None => continue,
                };
                let pb = match projected.get(b).and_then(|p| p.as_ref()) {
                    Some(p) => p,
                    None => continue,
                };
                let pc = match projected.get(c).and_then(|p| p.as_ref()) {
                    Some(p) => p,
                    None => continue,
                };
                let (ca, cb, cc) = match cell_color {
                    Some(c) => (c, c, c),
                    None => (
                        value_to_color(field.data[a]),
                        value_to_color(field.data[b]),
                        value_to_color(field.data[c]),
                    ),
                };
                let v0 = egui_mesh.vertices.len() as u32;
                egui_mesh.vertices.push(egui::epaint::Vertex {
                    pos: origin + egui::vec2(pa.x, pa.y),
                    uv: egui::epaint::WHITE_UV,
                    color: ca,
                });
                egui_mesh.vertices.push(egui::epaint::Vertex {
                    pos: origin + egui::vec2(pb.x, pb.y),
                    uv: egui::epaint::WHITE_UV,
                    color: cb,
                });
                egui_mesh.vertices.push(egui::epaint::Vertex {
                    pos: origin + egui::vec2(pc.x, pc.y),
                    uv: egui::epaint::WHITE_UV,
                    color: cc,
                });
                egui_mesh.indices.push(v0);
                egui_mesh.indices.push(v0 + 1);
                egui_mesh.indices.push(v0 + 2);
            }
        }
        if !egui_mesh.indices.is_empty() {
            painter.add(egui::Shape::Mesh(egui_mesh));
        }
        global_cell_offset += element_count;
    }
}

/// Sum of cells across every element block in a canonical mesh.
/// Used by the per-cell field overlay to validate that the field's
/// data length matches the cell count + to compute the global cell
/// offset for each block.
fn total_cell_count(mesh: &Mesh) -> usize {
    mesh.element_blocks
        .iter()
        .map(|b| {
            let npe = b.element_type.nodes_per_element();
            if npe == 0 {
                0
            } else {
                b.connectivity.len() / npe
            }
        })
        .sum()
}

/// Per-element-type triangulation table. Each tuple is a triangle's
/// three local-vertex indices (matching the connectivity layout
/// `valenx_mesh::ElementBlock` uses). Surface elements (Tri / Quad)
/// just produce their own face triangles; volume elements emit one
/// triangle per surface face (Tet -> 4, Pyr -> 6, Prism -> 8, Hex
/// -> 12).
///
/// 1-D elements (`Line2`) produce no triangles — return an empty
/// slice so the filled overlay quietly skips them.
fn triangles_for(et: ElementType) -> &'static [(usize, usize, usize)] {
    match et {
        ElementType::Line2 => &[],
        ElementType::Tri3 | ElementType::Tri6 => &[(0, 1, 2)],
        ElementType::Quad4 => &[(0, 1, 2), (0, 2, 3)],
        ElementType::Tet4 | ElementType::Tet10 => &[
            (0, 1, 2), // bottom face
            (0, 1, 3), // side 1
            (1, 2, 3), // side 2
            (2, 0, 3), // side 3
        ],
        ElementType::Pyr5 => &[
            // Quad base split into two triangles.
            (0, 1, 2),
            (0, 2, 3),
            // Four triangular sides up to the apex (vertex 4).
            (0, 1, 4),
            (1, 2, 4),
            (2, 3, 4),
            (3, 0, 4),
        ],
        ElementType::Prism6 => &[
            // Two triangular caps.
            (0, 1, 2),
            (3, 4, 5),
            // Three quadrilateral sides, each split into two triangles.
            (0, 1, 4),
            (0, 4, 3),
            (1, 2, 5),
            (1, 5, 4),
            (2, 0, 3),
            (2, 3, 5),
        ],
        ElementType::Hex8 | ElementType::Hex20 => &[
            // 6 quad faces × 2 triangles each = 12 triangles.
            // Bottom (0-1-2-3)
            (0, 1, 2),
            (0, 2, 3),
            // Top (4-5-6-7)
            (4, 5, 6),
            (4, 6, 7),
            // Front (0-1-5-4)
            (0, 1, 5),
            (0, 5, 4),
            // Right (1-2-6-5)
            (1, 2, 6),
            (1, 6, 5),
            // Back (2-3-7-6)
            (2, 3, 7),
            (2, 7, 6),
            // Left (3-0-4-7)
            (3, 0, 4),
            (3, 4, 7),
        ],
    }
}

/// Paint a small colour-bar legend in the bottom-right of the viewport
/// rect. Shows the field name on top, the cool-to-warm ramp as a
/// vertical strip, and min/max value labels alongside.
///
/// Pure painter work — no input, no state, just a few hundred pixels
/// in the corner. Looks like:
///
/// ```text
/// ┌──────────┐
/// │ p        │   ← field name
/// │ ▓ 1.234e2│   ← max
/// │ ▓        │
/// │ ▓        │
/// │ ▒        │
/// │ ▒        │
/// │ ░ -3.5e1 │   ← min
/// └──────────┘
/// ```
///
/// The strip itself is a stack of 32 thin horizontal rectangles, one
/// per ramp sample. That's enough resolution for the eye to read it
/// as a smooth gradient without measurable cost on the painter.
fn draw_field_legend(painter: &egui::Painter, rect: egui::Rect, field: &valenx_fields::Field) {
    // Only meaningful for scalar fields with a known range. Vector /
    // tensor / unranged fields fall through silently — the wireframe
    // overlay above already filtered to scalar OnNode, but we re-
    // check defensively so this helper stays standalone.
    if !matches!(field.kind, valenx_fields::FieldKind::Scalar) {
        return;
    }
    let (min, max) = match field.range {
        Some(r) => r,
        None => return,
    };

    // Layout: 14 px wide strip, 140 px tall, 8 px padding from the
    // viewport edges. The text labels sit to the strip's right.
    let strip_w = 14.0_f32;
    let strip_h = 140.0_f32;
    let pad = 8.0_f32;
    let label_room = 90.0_f32; // horizontal room for "1.234e+02" etc.

    let strip_right = rect.right() - pad - label_room;
    let strip_left = strip_right - strip_w;
    let strip_top = rect.bottom() - pad - strip_h;
    let strip_bottom = rect.bottom() - pad;

    // Backing card so the legend stays readable over busy meshes.
    let card = egui::Rect::from_min_max(
        egui::pos2(strip_left - 6.0, strip_top - 22.0),
        egui::pos2(rect.right() - pad + 2.0, strip_bottom + 6.0),
    );
    painter.rect_filled(
        card,
        4.0,
        egui::Color32::from_rgba_premultiplied(20, 22, 26, 200),
    );

    // Field name above the strip, with the timestep label on the
    // line below if it's anything other than steady. Two-line header
    // keeps both bits readable; transient runs need the timestep
    // visible to know which snapshot they're seeing.
    painter.text(
        egui::pos2(strip_left, strip_top - 18.0),
        egui::Align2::LEFT_TOP,
        &field.name,
        egui::FontId::monospace(11.0),
        egui::Color32::from_gray(220),
    );
    if !matches!(field.time, valenx_fields::TimeKey::Steady) {
        painter.text(
            egui::pos2(strip_left, strip_top - 6.0),
            egui::Align2::LEFT_TOP,
            crate::format_time_key(field.time),
            egui::FontId::monospace(9.0),
            egui::Color32::from_gray(170),
        );
    }

    // The gradient strip itself — 32 stacked slices, top = max
    // (warm), bottom = min (cool), matching scientific-plot
    // convention.
    const SLICES: usize = 32;
    let slice_h = strip_h / SLICES as f32;
    for i in 0..SLICES {
        // i=0 → top (t=1, warm); i=SLICES-1 → bottom (t=0, cool).
        let t = 1.0 - (i as f32 + 0.5) / SLICES as f32;
        let [r, g, b] = valenx_fields::colormap::cool_to_warm(t);
        let y0 = strip_top + i as f32 * slice_h;
        let slice = egui::Rect::from_min_max(
            egui::pos2(strip_left, y0),
            egui::pos2(strip_right, y0 + slice_h + 0.5),
        );
        painter.rect_filled(slice, 0.0, egui::Color32::from_rgb(r, g, b));
    }

    // Outline the strip so the colour boundaries don't bleed into the
    // dark card background.
    painter.rect_stroke(
        egui::Rect::from_min_max(
            egui::pos2(strip_left, strip_top),
            egui::pos2(strip_right, strip_bottom),
        ),
        0.0,
        egui::Stroke::new(1.0, egui::Color32::from_gray(80)),
    );

    // Min / max labels.
    painter.text(
        egui::pos2(strip_right + 6.0, strip_top),
        egui::Align2::LEFT_TOP,
        format_field_value(max),
        egui::FontId::monospace(10.0),
        egui::Color32::from_gray(220),
    );
    painter.text(
        egui::pos2(strip_right + 6.0, strip_bottom),
        egui::Align2::LEFT_BOTTOM,
        format_field_value(min),
        egui::FontId::monospace(10.0),
        egui::Color32::from_gray(220),
    );
}

/// Format a single scalar value for the colour-bar legend. Picks
/// scientific notation outside `[1e-3, 1e6)` and trimmed-decimal
/// otherwise — wide enough range to cover most CFD pressures and
/// thermal temperatures without label overflow.
fn format_field_value(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    let abs = v.abs();
    if !(1e-3..1e6).contains(&abs) {
        format!("{v:.3e}")
    } else {
        let s = format!("{v:.4}");
        let s = s.trim_end_matches('0').trim_end_matches('.').to_string();
        if s.is_empty() {
            "0".to_string()
        } else {
            s
        }
    }
}

/// Min/max sweep over a scalar field buffer. Used as a fallback when
/// the field's cached `range` somehow wasn't set at conversion time.
fn field_min_max(data: &[f64]) -> (f64, f64) {
    if data.is_empty() {
        return (0.0, 1.0);
    }
    let mut min = data[0];
    let mut max = data[0];
    for &v in &data[1..] {
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
    }
    (min, max)
}

/// Per-element-type edge list (index pairs into the element's
/// `connectivity` slice). Only covers the linear types the parser
/// produces; quadratic variants fall through the linear edge set.
fn edges_for(et: ElementType) -> &'static [(usize, usize)] {
    match et {
        ElementType::Line2 => &[(0, 1)],
        ElementType::Tri3 | ElementType::Tri6 => &[(0, 1), (1, 2), (2, 0)],
        ElementType::Quad4 => &[(0, 1), (1, 2), (2, 3), (3, 0)],
        ElementType::Tet4 | ElementType::Tet10 => &[(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)],
        ElementType::Pyr5 => &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (0, 4),
            (1, 4),
            (2, 4),
            (3, 4),
        ],
        ElementType::Prism6 => &[
            (0, 1),
            (1, 2),
            (2, 0),
            (3, 4),
            (4, 5),
            (5, 3),
            (0, 3),
            (1, 4),
            (2, 5),
        ],
        ElementType::Hex8 | ElementType::Hex20 => &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ],
    }
}

/// Axis-aligned bounding box over a canonical mesh's node
/// coordinates. Mirrors `lib::mesh_bounding_box` locally so the
/// viewport module doesn't need to import from the parent.
pub(crate) fn mesh_aabb(mesh: &Mesh) -> Option<([f32; 3], [f32; 3])> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shading_mode_default_is_shaded() {
        assert_eq!(ShadingMode::default(), ShadingMode::Shaded);
    }

    #[test]
    fn triangles_for_tet_emits_four_faces() {
        // A linear tet has 4 triangular faces — anything fewer means
        // we're hiding part of the surface, anything more means we're
        // double-counting and overdrawing.
        assert_eq!(triangles_for(ElementType::Tet4).len(), 4);
        assert_eq!(triangles_for(ElementType::Tet10).len(), 4);
    }

    #[test]
    fn triangles_for_hex_emits_twelve_faces() {
        // 6 quad faces × 2 triangles each = 12.
        assert_eq!(triangles_for(ElementType::Hex8).len(), 12);
        assert_eq!(triangles_for(ElementType::Hex20).len(), 12);
    }

    #[test]
    fn triangles_for_pyramid_and_prism_match_face_topology() {
        // Pyramid: square base (2 triangles) + 4 triangular sides = 6.
        assert_eq!(triangles_for(ElementType::Pyr5).len(), 6);
        // Prism: 2 triangular caps + 3 quad sides (2 triangles each) = 8.
        assert_eq!(triangles_for(ElementType::Prism6).len(), 8);
    }

    #[test]
    fn triangles_for_surface_elements_emits_one_or_two_triangles() {
        assert_eq!(triangles_for(ElementType::Tri3).len(), 1);
        assert_eq!(triangles_for(ElementType::Tri6).len(), 1);
        assert_eq!(triangles_for(ElementType::Quad4).len(), 2);
    }

    #[test]
    fn triangles_for_line_emits_no_triangles() {
        // 1-D elements have no surface — the filled overlay must
        // skip them rather than panic.
        assert_eq!(triangles_for(ElementType::Line2).len(), 0);
    }

    #[test]
    fn total_cell_count_sums_blocks() {
        use nalgebra::Vector3;
        use valenx_mesh::{ElementBlock, ElementType, Mesh};
        let mut m = Mesh::new("smoke");
        m.nodes = vec![Vector3::zeros(); 12];
        // 2 tets in block 0 (4 nodes each = 8 connectivity entries)
        let mut tets = ElementBlock::new(ElementType::Tet4);
        tets.connectivity = vec![0; 8];
        // 3 hexes in block 1 (8 nodes each = 24 connectivity entries)
        let mut hexes = ElementBlock::new(ElementType::Hex8);
        hexes.connectivity = vec![0; 24];
        m.element_blocks = vec![tets, hexes];
        assert_eq!(total_cell_count(&m), 5);
    }

    #[test]
    fn total_cell_count_handles_unknown_element_types() {
        use valenx_mesh::{ElementBlock, ElementType, Mesh};
        // Line2 has 2 nodes per element; should still count as 1 cell.
        let mut m = Mesh::new("lines");
        m.nodes = vec![nalgebra::Vector3::zeros(); 4];
        let mut lines = ElementBlock::new(ElementType::Line2);
        lines.connectivity = vec![0, 1, 2, 3];
        m.element_blocks = vec![lines];
        // 4 connectivity entries / 2 npe = 2 line cells
        assert_eq!(total_cell_count(&m), 2);
    }

    #[test]
    fn total_cell_count_for_empty_mesh_is_zero() {
        use valenx_mesh::Mesh;
        let m = Mesh::new("empty");
        assert_eq!(total_cell_count(&m), 0);
    }

    #[test]
    fn mesh_surface_from_tris_keeps_every_triangle() {
        use nalgebra::Vector3;
        use valenx_mesh::{ElementBlock, ElementType, Mesh};
        let mut m = Mesh::new("two-tris");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut tris = ElementBlock::new(ElementType::Tri3);
        tris.connectivity = vec![0, 1, 2, 0, 2, 3];
        m.element_blocks = vec![tris];
        let surf = mesh_to_triangle_surface(&m).expect("tri mesh yields a surface");
        assert_eq!(surf.triangles.len(), 2);
        // Flat square in the z = 0 plane → unit +/-z face normals.
        for t in &surf.triangles {
            assert!(
                (t.normal[2].abs() - 1.0).abs() < 1e-5,
                "expected an axial normal, got {:?}",
                t.normal
            );
        }
    }

    #[test]
    fn mesh_surface_splits_each_quad_into_two_triangles() {
        use nalgebra::Vector3;
        use valenx_mesh::{ElementBlock, ElementType, Mesh};
        let mut m = Mesh::new("one-quad");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut quads = ElementBlock::new(ElementType::Quad4);
        quads.connectivity = vec![0, 1, 2, 3];
        m.element_blocks = vec![quads];
        let surf = mesh_to_triangle_surface(&m).expect("quad mesh yields a surface");
        assert_eq!(surf.triangles.len(), 2);
    }

    #[test]
    fn mesh_surface_is_none_for_volumetric_and_empty() {
        use valenx_mesh::{ElementBlock, ElementType, Mesh};
        // Volumetric elements emit no direct surface here → wireframe path.
        let mut m = Mesh::new("tets");
        m.nodes = vec![nalgebra::Vector3::zeros(); 4];
        let mut tets = ElementBlock::new(ElementType::Tet4);
        tets.connectivity = vec![0, 1, 2, 3];
        m.element_blocks = vec![tets];
        assert!(mesh_to_triangle_surface(&m).is_none());
        // Empty mesh → None.
        assert!(mesh_to_triangle_surface(&Mesh::new("empty")).is_none());
    }

    #[test]
    fn triangles_for_indices_stay_within_element_node_count() {
        // Local-index sanity: every triangle's three indices must be
        // < the element type's nodes_per_element. Catches off-by-one
        // typos in the topology table.
        for et in [
            ElementType::Tri3,
            ElementType::Quad4,
            ElementType::Tet4,
            ElementType::Pyr5,
            ElementType::Prism6,
            ElementType::Hex8,
        ] {
            let npe = et.nodes_per_element();
            for &(a, b, c) in triangles_for(et) {
                assert!(a < npe, "{et:?} triangle vertex a={a} >= npe={npe}");
                assert!(b < npe, "{et:?} triangle vertex b={b} >= npe={npe}");
                assert!(c < npe, "{et:?} triangle vertex c={c} >= npe={npe}");
            }
        }
    }

    #[test]
    fn format_field_value_picks_readable_shape() {
        // Whole numbers stay short and integer-looking.
        assert_eq!(format_field_value(0.0), "0");
        assert_eq!(format_field_value(101325.0), "101325");
        assert_eq!(format_field_value(-1.0), "-1");
        // Small fractionals keep up to 4 decimals, trim trailing zeros.
        assert_eq!(format_field_value(1.5), "1.5");
        assert_eq!(format_field_value(0.001234), "0.0012");
        // Below 1e-3 falls into scientific notation so labels stay
        // narrow even for tiny gradients (e.g. residuals or species
        // mass fractions in chemistry).
        assert_eq!(format_field_value(1e-5), "1.000e-5");
        // Above 1e6 also goes scientific so the label doesn't
        // overflow the legend's label_room budget.
        assert_eq!(format_field_value(2.5e8), "2.500e8");
    }

    #[test]
    fn shade_color_within_bounds() {
        let c = shade_color(0.0);
        // Clamp floor should keep the fully-in-shadow colour above
        // pure black so faces never disappear.
        assert!(c.r() >= 28);
        let bright = shade_color(1.0);
        assert!(bright.r() > c.r());
    }

    #[test]
    fn normalize_of_zero_is_safe() {
        let n = normalize([0.0, 0.0, 0.0]);
        assert_eq!(n, [0.0, 0.0, 1.0]);
    }

    #[test]
    fn mesh_aabb_over_unit_tet() {
        use nalgebra::Vector3;
        use valenx_mesh::Mesh;
        let mut m = Mesh::new("t");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let (min, max) = super::mesh_aabb(&m).expect("non-empty");
        assert_eq!(min, [0.0, 0.0, 0.0]);
        assert_eq!(max, [1.0, 1.0, 1.0]);
    }
}
